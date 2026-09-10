"""Execution backends of the reversed-rate fixed point.

WHAT MAKES THIS SAFE IS THE DECOMPOSITION, NOT THE SCHEDULING. Agent k's
generator is

    Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c

so an agent reads the rest of the model only through the scalar reversed rates
x, and it writes only its own slot of the sweep's output. The sweep is Jacobi --
every x_a is read off the PREVIOUS sweep's stationary vectors and only then do
the agents re-solve -- so the agent order is immaterial. A parallel or
distributed sweep therefore produces the SAME iterates as the serial one.

How far that survives floating point depends on WHO runs the agent solve:

* ``threads`` is bit-identical to ``serial``. Same code, same process, only the
  order of evaluation differs, and the order is what does not matter.
* ``cluster`` is bit-identical only when the worker runs the same implementation
  as the coordinator. The wire itself is exact -- JSON round-trips a double
  without loss -- but the stationary vector comes back from the WORKER's solve,
  so a Python coordinator driving a Java ag-worker agrees to a few ulp rather
  than bit for bit (measured: 1.1e-16 on the M/M/1 tandem). That is the ordinary
  cross-codebase difference the parity harness already tolerates, not a protocol
  defect, but it means a cluster run must not be used to regenerate a seeded
  golden that a same-language run will later read.

Three backends, differing only in who evaluates an agent:

* ``serial``  -- the caller's own loop; represented by ``None``, so the default
  path carries no scheduling machinery at all.
* ``threads`` -- a thread pool, one task per agent, barrier per sweep. Threads
  rather than processes because the cost of an agent is ``ctmc_solve``, i.e. a
  LAPACK call that RELEASES THE GIL; a process pool would have to pickle every
  agent's matrices per sweep and would lose more than it gained.
* ``cluster`` -- agents partitioned over ag-worker processes, which may be on
  other machines. The static half of each agent is shipped once and only x
  crosses the wire per sweep.

A cluster worker that is missing, slow or broken is NOT fatal: its agents are
solved on the coordinator through the same agent path, so the answer is the
run's answer either way and only the wall clock changes.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import json
import socket
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np

EXEC_SERIAL = 'serial'
#: Fan the agents of a sweep out over a local thread pool.
#:
#: Named 'parallel', with 'para' accepted as an alias -- the same pair SolverSSA's
#: replica analyzer answers to, so one spelling convention covers both solvers. It
#: was called 'threads' until 2026-08-19; that name is no longer accepted, and the
#: alias is normalised in ``create`` so nothing downstream sees two spellings.
EXEC_PARALLEL = 'parallel'
EXEC_PARA = 'para'
EXEC_CLUSTER = 'cluster'
_MODES = (EXEC_SERIAL, EXEC_PARALLEL, EXEC_CLUSTER)


def create(config, method):
    """Resolve the backend named by ``config``; None means the serial loop."""
    if not config:
        return None
    mode = str(config.get('exec', EXEC_SERIAL)).lower()
    if mode == EXEC_PARA:
        mode = EXEC_PARALLEL
    if mode not in _MODES:
        # 'threads' was this backend's name until 2026-08-19. Name the rename
        # rather than reporting a backend that still exists as unknown.
        if mode == 'threads':
            raise ValueError("The 'threads' execution backend was renamed to 'parallel' "
                             "(alias 'para').")
        raise ValueError("Unknown AG execution backend %r. Use 'serial', 'parallel' "
                         "(alias 'para') or 'cluster'." % mode)
    if mode == EXEC_SERIAL:
        return None
    if mode == EXEC_CLUSTER and method == 'inapinf':
        # The remote worker implements the FINITE agent solve. 'inapinf' replaces
        # it with the matrix-geometric treatment of an open agent -- Neuts' R
        # matrix and the scalar-tail detection that precedes it -- which the
        # worker does not carry, and answering with the finite solve instead
        # would silently change the method.
        raise ValueError("The 'cluster' execution backend does not carry the 'inapinf' "
                         "agent solve (the matrix-geometric tail of an open agent runs on "
                         "the coordinator only). Use exec 'serial' or 'parallel' with "
                         "'inapinf', or method 'inap'/'inapplus' with 'cluster'.")
    if mode == EXEC_PARALLEL:
        return ParallelBackend(int(config.get('nworkers', 0) or 0))
    endpoints = list(config.get('endpoints') or [])
    if not endpoints:
        raise ValueError("The 'cluster' execution backend needs worker endpoints: set "
                         "options.config['endpoints'] to a list of 'host:port' strings, each "
                         "one an ag-worker started with "
                         "'java -cp jline.jar jline.cli.AgWorker -p <port>'.")
    return ClusterBackend(endpoints, float(config.get('worker_timeout', 30.0)))


class ParallelBackend:
    """Fan the agents of a sweep out over a thread pool."""

    def __init__(self, nworkers=0):
        self.nworkers = nworkers or None
        self._pool = None

    def sweep(self, x, model, agent_solve, agent_generator, solve_component):
        if self._pool is None:
            self._pool = ThreadPoolExecutor(max_workers=self.nworkers)
        n = model.num_processes
        futures = [self._pool.submit(agent_solve, k, x, model) for k in range(n)]
        Q_list = [None] * n
        pi_list = [None] * n
        for k, f in enumerate(futures):
            # An agent solve that raises is a defect in the agent, not in the
            # scheduling, so it must surface as itself rather than as a
            # half-filled sweep that fails later somewhere unrelated.
            Q_list[k], pi_list[k] = f.result()
        return pi_list, Q_list

    def sweep_qbd(self, x, model, is_open_proc, agent_solve_qbd):
        """The 'inapinf' twin of sweep: the same fan-out over the qbd agent."""
        if self._pool is None:
            self._pool = ThreadPoolExecutor(max_workers=self.nworkers)
        n = model.num_processes
        futures = [self._pool.submit(agent_solve_qbd, k, x, model, is_open_proc)
                   for k in range(n)]
        Q_list = [None] * n
        pi_list = [None] * n
        rho_proc = np.zeros(n)
        is_geom_proc = np.zeros(n, dtype=bool)
        geom_data = [None] * n
        for k, f in enumerate(futures):
            Q_list[k], pi_list[k], rho_proc[k], is_geom_proc[k], geom_data[k] = f.result()
        return pi_list, Q_list, rho_proc, is_geom_proc, geom_data

    def close(self):
        if self._pool is not None:
            self._pool.shutdown(wait=False)
            self._pool = None


class ClusterBackend:
    """Partition the agents over remote ag-worker processes.

    Speaks the same newline-delimited JSON protocol as ``jline.solvers.ag.AgWire``:
    one ``assign`` per solve carrying the static half of each owned agent, then one
    ``sweep`` per iteration carrying the reversed rates.
    """

    def __init__(self, endpoints, timeout=30.0):
        self.endpoints = list(endpoints)
        self.timeout = timeout
        self._conns = [None] * len(self.endpoints)
        self._owns = None
        self._live = [True] * len(self.endpoints)
        self._assigned = False

    def _partition(self, num_processes):
        # Round-robin in agent index order, computed before any connection is
        # attempted so a dead worker does not shift the others' agents. A pure
        # function of (agent count, worker count), so a rerun assigns the same
        # agents to the same workers.
        nw = len(self.endpoints)
        self._owns = [[] for _ in range(nw)]
        for k in range(num_processes):
            self._owns[k % nw].append(k)

    def _connect(self, w):
        host, _, port = self.endpoints[w].rpartition(':')
        sock = socket.create_connection((host, int(port)), timeout=self.timeout)
        sock.settimeout(self.timeout)
        self._conns[w] = (sock, sock.makefile('rw', encoding='utf-8', newline='\n'))

    def _call(self, w, msg):
        sock, f = self._conns[w]
        f.write(json.dumps(msg) + '\n')
        f.flush()
        line = f.readline()
        if not line:
            raise IOError('worker closed the connection')
        reply = json.loads(line)
        if 'error' in reply:
            raise IOError(reply['error'])
        return reply

    @staticmethod
    def _triplets(M):
        """Non-zero entries of M as [row, col, value], 0-based."""
        if M is None:
            return []
        rows, cols = np.nonzero(M)
        return [[int(i), int(j), float(M[i, j])] for i, j in zip(rows, cols)]

    def _agent_payload(self, k, model):
        num_actions = model.num_actions
        N = model.N
        AP = model.AP
        R = model.R
        n = int(N[k])
        L = R.get((num_actions, k), np.zeros((n, n)))
        passive = []
        active = []
        for c in range(num_actions):
            if AP[c, 1] == k:
                passive.append({'c': int(c),
                                'M': self._triplets(R.get((c, 1), np.zeros((n, n))))})
            elif AP[c, 0] == k:
                active.append({'c': int(c),
                               'M': self._triplets(R.get((c, 0), np.zeros((n, n))))})
        return {
            'k': int(k),
            'n': n,
            'mph': int(model.mph[k]),
            'nlev': int(model.nlev[k]),
            'level': [int(v) for v in np.asarray(model.level[k]).ravel()],
            'L': self._triplets(L),
            'passive': passive,
            'active': active,
        }

    def sweep(self, x, model, agent_solve, agent_generator, solve_component):
        n = model.num_processes
        if self._owns is None:
            self._partition(n)

        if not self._assigned:
            for w in range(len(self.endpoints)):
                if not self._owns[w]:
                    continue
                try:
                    self._connect(w)
                    payload = [self._agent_payload(k, model) for k in self._owns[w]]
                    reply = self._call(w, {'op': 'assign', 'agents': payload})
                    if reply.get('op') != 'assigned':
                        raise IOError('worker did not acknowledge the assignment')
                except Exception as e:      # noqa: BLE001 - any failure is non-fatal
                    warnings.warn('AG worker %s is unreachable (%s); its %d agent(s) run on '
                                  'the coordinator instead.'
                                  % (self.endpoints[w], e, len(self._owns[w])))
                    self._live[w] = False
            self._assigned = True

        # The generator is rebuilt here in any case: the metrics stage reads it,
        # and rebuilding is cheaper than transporting it.
        Q_list = [agent_generator(k, x, model) for k in range(n)]
        pi_list = [None] * n

        for w in range(len(self.endpoints)):
            if not self._live[w] or not self._owns[w]:
                continue
            try:
                reply = self._call(w, {'op': 'sweep',
                                       'x': [float(v) for v in np.asarray(x).ravel()]})
                if reply.get('op') != 'swept':
                    raise IOError('unexpected reply %r' % reply.get('op'))
                for entry in reply['agents']:
                    k = int(entry['k'])
                    if k not in self._owns[w]:
                        raise IOError('worker answered for agent %d, which it was never '
                                      'assigned' % k)
                    pi_list[k] = np.asarray(entry['pi'], dtype=float)
            except Exception as e:          # noqa: BLE001 - any failure is non-fatal
                warnings.warn('AG worker %s failed mid-sweep (%s); its %d agent(s) are solved '
                              'on the coordinator for the rest of the run.'
                              % (self.endpoints[w], e, len(self._owns[w])))
                self._live[w] = False

        # Whatever no worker answered for, solve here.
        for k in range(n):
            if pi_list[k] is None:
                pi_list[k] = solve_component(Q_list[k], int(model.mph[k]),
                                             int(model.nlev[k]), model.level[k])
        return pi_list, Q_list

    def close(self):
        for w, conn in enumerate(self._conns):
            if conn is None:
                continue
            sock, f = conn
            try:
                f.write(json.dumps({'op': 'bye'}) + '\n')
                f.flush()
            except Exception:               # noqa: BLE001 - already gone
                pass
            try:
                sock.close()
            except Exception:               # noqa: BLE001
                pass
            self._conns[w] = None
