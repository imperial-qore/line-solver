#!/usr/bin/env python3
"""Dramatiq actor and dispatcher behind ``run-tests.sh --dramatiq``.

``--parallel`` runs a wave's phases side by side ON ONE HOST, in tmux sessions.
``--dramatiq`` runs the same waves ACROSS THE PICARD CLUSTER: every phase of a
wave is one Dramatiq message, one ephemeral worker per host drains its own
queue, and the driver waits for the wave exactly as it does for the tmux group
-- by polling the per-phase ``.rc`` file each phase leaves behind.

Three things are deliberate here.

* THE BROKER IS SHARED WITH THE SOLVERNN DATAGEN CLUSTER, THE QUEUES ARE NOT.
  Redis on picard00:6379 is the same server, but the datagen cluster lives in
  db 0 on ``dramatiq:default`` and this runs in db 3 (``LINE_DRAMATIQ_BROKER``)
  on queues named after the run id, so a test run neither sees nor delays a
  datagen task. The workers are this run's own -- started over ssh when the run
  starts, killed when it ends -- and never the deployed ``dramatiq-worker``
  containers, whose image holds SolverNN, not a LINE checkout or its toolchain.

* ONE QUEUE PER HOST, BECAUSE THE ASSIGNMENT IS THE DRIVER'S TO MAKE. A single
  shared queue would let whichever worker happened to poll first take the work,
  which is neither round-robin nor capability-aware. Instead the driver picks
  the host -- round-robin over those that can actually run the phase -- and
  enqueues on that host's own queue, which only that host's worker subscribes
  to. ``LINE_DRAMATIQ_QUEUES`` carries the full set so both sides declare the
  same queues.

* NO RESULT BACKEND. The actor writes the phase's output, its host, its elapsed
  time and its exit status to files under the shared run directory, in the same
  ``phase<id>-<slug>.{txt,secs,rc,host}`` layout ``--parallel`` uses, and writes
  the ``.rc`` LAST so its appearance is the completion signal. The driver's
  waiting and log-splicing code is therefore the same for both modes.

Entry points (all take the broker URL, run id and queue set from the
environment):

    python -m dramatiq dramatiq_phases --processes 1 --threads 1 --queues <q>
    python dramatiq_phases.py enqueue <payload.json>
    python dramatiq_phases.py purge
    python dramatiq_phases.py purge-all
"""

import json
import os
import socket
import subprocess
import sys
import time

import dramatiq
from dramatiq.brokers.redis import RedisBroker

# db 3, not the datagen cluster's db 0. Overridable for a different broker host.
DEFAULT_BROKER = "redis://picard00:6379/3"

# A phase may legitimately run for hours (the parity sub-suites have no timeout
# of their own), so the actor's time limit is a runaway guard, not a schedule.
DEFAULT_TIME_LIMIT_MS = 7 * 24 * 60 * 60 * 1000

RUN_ID = os.environ.get("LINE_DRAMATIQ_RUN", "")
BROKER_URL = os.environ.get("LINE_DRAMATIQ_BROKER", DEFAULT_BROKER)
TIME_LIMIT_MS = int(os.environ.get("LINE_DRAMATIQ_TIME_LIMIT_MS", DEFAULT_TIME_LIMIT_MS))

# Every per-host queue of this run, so the enqueuing side and each worker agree
# on what exists. A worker still CONSUMES only the one queue named on its
# command line; declaring the rest costs nothing and keeps the module symmetric.
QUEUES = [q for q in os.environ.get("LINE_DRAMATIQ_QUEUES", "").split(",") if q]

# The actor's own queue name has to be fixed at decoration time and is never
# sent to: every message is redirected to a host queue before it is enqueued.
DISPATCH_QUEUE = "linetests_%s_dispatch" % (RUN_ID or "unset")


def _connect():
    if not RUN_ID:
        sys.exit("LINE_DRAMATIQ_RUN is unset: refusing to share queues between runs.")
    broker = RedisBroker(url=BROKER_URL)
    dramatiq.set_broker(broker)
    broker.declare_queue(DISPATCH_QUEUE)
    for queue in QUEUES:
        broker.declare_queue(queue)
    return broker


broker = _connect()


@dramatiq.actor(queue_name=DISPATCH_QUEUE, max_retries=0, time_limit=TIME_LIMIT_MS)
def run_phase(phase_id, label, command, stem, env):
    """Run one phase of the suite on this worker's host.

    ``command`` is the very string the sequential driver would hand to ``bash``,
    already pointing at the staged tree; ``stem`` is the shared-filesystem path
    prefix the four result files hang off; ``env`` holds the few LINE_* settings
    the phase needs on top of the worker's own login environment.
    """
    host = socket.gethostname()
    started = time.time()
    with open(stem + ".host", "w") as handle:
        handle.write("%s %s\n" % (host, time.strftime("%Y-%m-%dT%H:%M:%S%z")))

    environ = dict(os.environ)
    environ.update({str(k): str(v) for k, v in (env or {}).items()})

    with open(stem + ".txt", "wb", buffering=0) as log:
        log.write(("=== phase %s (%s) on %s at %s\n\n"
                   % (phase_id, label, host, time.strftime("%Y-%m-%dT%H:%M:%S%z"))).encode())
        log.flush()
        try:
            rc = subprocess.call(["bash", "-c", command], stdout=log,
                                 stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                 env=environ)
        except BaseException as exc:            # the worker must still report an outcome
            log.write(("\n=== phase %s aborted on %s: %r\n" % (phase_id, host, exc)).encode())
            rc = 1
            raise
        finally:
            elapsed = int(time.time() - started)
            # .secs before .rc: the driver treats the .rc file as the completion
            # signal and reads the elapsed time straight after seeing it.
            with open(stem + ".secs", "w") as handle:
                handle.write("%d\n" % elapsed)
            with open(stem + ".rc", "w") as handle:
                handle.write("%d\n" % rc)


def _enqueue(payload_path):
    """Enqueue one wave.

    Each phase names the queue it goes on, which the driver has already chosen:
    the host is picked round-robin among those whose probe reported the
    capability the phase needs.
    """
    with open(payload_path) as handle:
        phases = json.load(handle)
    for phase in phases:
        queue = phase["queue"]
        if queue not in QUEUES:
            sys.exit("phase %s names queue %r, which this run never declared"
                     % (phase.get("id"), queue))
        message = run_phase.message(phase["id"], phase["label"], phase["command"],
                                    phase["stem"], phase.get("env", {}))
        broker.enqueue(message.copy(queue_name=queue))
        print("enqueued phase %s (%s) -> %s" % (phase["id"], phase["label"],
                                                phase.get("host", queue)))


def _purge():
    """Delete this run's queues, so a killed worker cannot resurrect a phase."""
    for queue in [DISPATCH_QUEUE] + QUEUES:
        try:
            broker.flush(queue)
        except Exception as exc:                # a queue that never existed is not an error
            print("purge %s: %s" % (queue, exc))


def _purge_all():
    """Delete EVERY run's queues, not just this one's: the ``--reset`` path.

    A reset has already killed the drivers and the workers of every run in
    flight, so what is left here is a message nobody will ever ack sitting on a
    queue nobody will drain -- and a worker that came back would pick it up and
    run a phase against a tree the reset is about to unstage.

    THE QUEUE NAMES ARE READ BACK OUT OF REDIS rather than taken from
    LINE_DRAMATIQ_QUEUES, because the whole point of the sweep is the runs this
    process was never told about. The key layout is the redis broker's own
    (dispatch.lua): ``<ns>:<queue>`` is the message-id list and ``.msgs``,
    ``.XQ``, ``.XQ.msgs``, ``.DQ`` and ``.DQ.msgs`` hang off it, so the base name
    is whatever is left once one of those suffixes is stripped.
    """
    namespace = broker.namespace
    prefix = ("%s:" % namespace).encode()
    # Longest first: ".XQ.msgs" also ends in ".msgs", and stripping the shorter
    # one would leave ".XQ" behind as a queue of its own.
    suffixes = (".XQ.msgs", ".DQ.msgs", ".XQ", ".DQ", ".msgs")
    names = set()
    for key in broker.client.scan_iter(match="%s:linetests_*" % namespace, count=500):
        name = key[len(prefix):].decode()
        for suffix in suffixes:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        names.add(name)

    for name in sorted(names):
        try:
            broker.flush(name)
            print("purged queue %s" % name)
        except Exception as exc:                # a queue mid-deletion is not an error
            print("purge %s: %s" % (name, exc))

    # ``purge`` in the lua drops only the ack set of the worker asking, so the
    # ack sets of the workers this reset just killed outlive their queues.
    acks = [key for key in broker.client.scan_iter(
        match="%s:__acks__.*linetests_*" % namespace, count=500)]
    if acks:
        broker.client.delete(*acks)
        print("dropped %d stale ack set(s)" % len(acks))

    if not names and not acks:
        print("no linetests_* queue on %s" % BROKER_URL)
    # Left alone deliberately: a worker's heartbeat expires on its own
    # (heartbeat_timeout, 60s by default) and dramatiq's maintenance pass drops
    # it. Reported so a reset that found workers says what is still winding down.
    # OUR OWN GOES FIRST, because every dispatch registers one -- the purges above
    # included -- and counting it would report this sweep to itself as a worker.
    heartbeats = "%s:__heartbeats__" % namespace
    broker.client.zrem(heartbeats, broker.broker_id)
    stale = broker.client.zcard(heartbeats)
    if stale:
        print("%d worker heartbeat(s) still registered; they expire on their own" % stale)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "enqueue":
        _enqueue(sys.argv[2])
    elif sys.argv[1] == "purge":
        _purge()
    elif sys.argv[1] == "purge-all":
        _purge_all()
    else:
        sys.exit("unknown subcommand: %s" % sys.argv[1])
