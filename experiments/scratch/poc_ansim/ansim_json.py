"""Model builders for the ansim priority comparison, in LINE's portable JSON
format, driven through the C++ `line-cli` binary.

This is the lock-free path: it needs no MATLAB, so it can run while the host
MATLAB lock is contended. It covers every integer-population variant. The
fractional-field 'mn' mode has no JSON route -- the schema types a class
population as an integer -- so that one stays in MATLAB.
"""
import json, subprocess, tempfile, os, time

LINE_CLI = "/data/gcasale/line-dev.git/common/line-cli"

# ---------------------------------------------------------------- distributions

def svc(mean, law):
    """Service law at the priority queue, at fixed mean."""
    if law == "exp":
        return {"type": "Exp", "params": {"lambda": 1.0 / mean}}
    if law == "hypoexp":                      # SCV = 0.5 -> Erlang-2
        return {"type": "Erlang", "params": {"lambda": 2.0 / mean, "k": 2}}
    if law == "hyperexp":                     # SCV = 4, balanced-means fit
        scv = 4.0
        # LINE's HyperExp.fitMeanAndSCV: balanced means, p from the SCV
        p = 0.5 * (1 + ((scv - 1) / (scv + 1)) ** 0.5)
        l1 = 2 * p / mean
        l2 = 2 * (1 - p) / mean
        return {"type": "HyperExp", "params": {"p": [p, 1 - p], "lambda": [l1, l2]}}
    raise ValueError(law)


def dps_weights(cfg):
    mx = max(cfg["prio"])
    return [cfg["Wdps"] ** (mx - p) for p in cfg["prio"]]


SCHED = {"psprio": "PSPRIO", "dps": "DPS", "hol": "HOL",
         "prs": "FCFSPRPRIO", "fcfs": "FCFS"}

# ---------------------------------------------------------------- base model

def build_model(cfg, variant):
    K = len(cfg["N"])
    w = dps_weights(cfg)
    prio = [0] * K if variant in ("dps", "fcfs") else cfg["prio"]
    names = ["Class%d" % (r + 1) for r in range(K)]

    think = {"name": "Think", "type": "Delay", "service": {}}
    q = {"name": "Q", "type": "Queue", "scheduling": SCHED[variant], "service": {}}
    for r, nm in enumerate(names):
        think["service"][nm] = {"type": "Exp", "params": {"lambda": 1.0 / cfg["Z"][r]}}
        q["service"][nm] = svc(cfg["D"][r], cfg["law"])
    if variant == "dps":
        q["schedParams"] = {nm: w[r] for r, nm in enumerate(names)}

    classes = [{"name": nm, "type": "Closed", "population": int(cfg["N"][r]),
                "refNode": "Think", "priority": int(prio[r])}
               for r, nm in enumerate(names)]

    return _net("prio_" + variant, [think, q], classes,
                {nm: [("Think", "Q"), ("Q", "Think")] for nm in names})


# ---------------------------------------------------------------- ansim surrogate

def build_surrogate(cfg, r, field, mode):
    """One tagged class-r job cycling over a frozen self-looping field.

    mode 'ctmcP' pins every other class; mode 'live' keeps the classes strictly
    higher in priority than the tagged one dynamic, because under preemption a
    pinned high-priority job never releases the server and starves the tagged
    job outright.
    """
    K = len(cfg["N"])
    live = [cfg["prio"][s] < cfg["prio"][r] for s in range(K)] if mode == "live" \
        else [False] * K

    pop = [round_preserving(field[s], cfg["N"][s] - (1 if s == r else 0))
           for s in range(K)]

    think = {"name": "Think", "type": "Delay", "service": {}}
    q = {"name": "Q", "type": "Queue", "scheduling": "PSPRIO", "service": {}}
    classes, routes = [], {}

    for s in range(K):
        if live[s]:
            nm = "L%d" % s
            classes.append({"name": nm, "type": "Closed",
                            "population": int(cfg["N"][s]), "refNode": "Think",
                            "priority": int(cfg["prio"][s])})
            think["service"][nm] = {"type": "Exp", "params": {"lambda": 1.0 / cfg["Z"][s]}}
            q["service"][nm] = svc(cfg["D"][s], cfg["law"])
            routes[nm] = [("Think", "Q"), ("Q", "Think")]
            continue
        for i, node in enumerate(("Think", "Q")):
            nm = "F%d%d" % (i, s)
            classes.append({"name": nm, "type": "SelfLooping",
                            "population": int(pop[s][i]), "refNode": node,
                            "priority": int(cfg["prio"][s])})
            # A pinned class never leaves its station: it self-loops there and is
            # Disabled at the other one. Routing it round the serial chain would
            # send it to a station where it has no service.
            if i == 0:
                think["service"][nm] = {"type": "Exp", "params": {"lambda": 1.0 / cfg["Z"][s]}}
                q["service"][nm] = {"type": "Disabled"}
            else:
                think["service"][nm] = {"type": "Disabled"}
                q["service"][nm] = svc(cfg["D"][s], cfg["law"])
            routes[nm] = [(node, node)]

    classes.append({"name": "Moving", "type": "Closed", "population": 1,
                    "refNode": "Think", "priority": int(cfg["prio"][r])})
    think["service"]["Moving"] = {"type": "Exp", "params": {"lambda": 1.0 / cfg["Z"][r]}}
    q["service"]["Moving"] = svc(cfg["D"][r], cfg["law"])
    routes["Moving"] = [("Think", "Q"), ("Q", "Think")]

    return _net("ansim_surrogate", [think, q], classes, routes)


def round_preserving(x, total):
    """Largest-remainder rounding of x to integers summing exactly to total."""
    x = [max(v, 0.0) for v in x]
    n = [0] * len(x)
    if total <= 0:
        return n
    s = sum(x)
    if s <= 0:
        n[0] = total
        return n
    x = [v * total / s for v in x]
    n = [int(v // 1) for v in x]
    short = total - sum(n)
    order = sorted(range(len(x)), key=lambda i: x[i] - n[i], reverse=True)
    for i in order[:short]:
        n[i] += 1
    return n


def _net(name, nodes, classes, routes):
    matrix = {}
    for nm, arcs in routes.items():
        matrix["%s,%s" % (nm, nm)] = {a: {b: 1.0} for a, b in arcs}
    return {"format": "line-model", "version": "1.0",
            "model": {"type": "Network", "name": name, "nodes": nodes,
                      "classes": classes,
                      "routing": {"type": "matrix", "matrix": matrix}}}


# ---------------------------------------------------------------- solving

def solve(model, solver, method=None, samples=None, seed=None, timeout=5400):
    """Return (QLen as {(station,class): value}, seconds, error-or-None)."""
    cmd = [LINE_CLI, "-s", solver, "-a", "avg", "-o", "json", "-v", "silent"]
    if method:
        cmd += ["--method", method]
    if samples:
        # --seed is accepted only by the simulation solvers, so it rides with
        # --samples and never goes to an analytical solver.
        cmd += ["--samples", str(int(samples))]
        if seed is not None:
            cmd += ["--seed", str(int(seed))]
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(model, fh)
        path = fh.name
    cmd += ["-f", path]
    t0 = time.time()
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        os.unlink(path)
        return None, time.time() - t0, "timeout"
    os.unlink(path)
    dt = time.time() - t0
    blob = None
    for line in out.stdout.splitlines():
        line = line.strip()
        if line.startswith("{"):
            blob = line
    if blob is None:
        err = (out.stderr or out.stdout).strip().replace("\n", " ")
        return None, dt, err[:200] or "no output"
    try:
        d = json.loads(blob)["avg"]
    except Exception as e:
        return None, dt, "parse: %s" % e
    q = {}
    for st, cl, v in zip(d["Station"], d["JobClass"], d["QLen"]):
        q[(st, cl)] = v
    return q, dt, None
