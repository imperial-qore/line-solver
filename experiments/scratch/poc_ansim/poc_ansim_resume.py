"""Finish the C++ (line-cli) sweep, skipping whatever is already on disk.

Competitor results land in poc_ansim_prio_cpp.json and ansim results in
ansim_only.json; both are re-read on start, so a run that is interrupted picks
up where it stopped instead of repeating ~10s-per-solve work.
"""
import sys, json, os, time
sys.path.insert(0, "/data/gcasale/line-dev.git/experiments/scratch/poc_ansim")
import ansim_json as A
from poc_ansim_prio import ansim, qmat, l1, BASE, LAWS, NS, SEED, LDES_SAMPLES

HERE = "/data/gcasale/line-dev.git/experiments/scratch/poc_ansim"
COMP = os.path.join(HERE, "poc_ansim_prio_cpp.json")
ANSIM = os.path.join(HERE, "ansim_only.json")

JOBS = [
    # Approximations and simulation first: they are seconds each and are the
    # numbers the comparison actually turns on.
    ("FLDmn",   "dps",    "fluid", "minnormal", None),
    ("MVAdps",  "dps",    "mva",   None,        None),
    ("MVAcl",   "hol",    "mva",   None,        None),
    ("MVAprs",  "prs",    "mva",   None,        None),
    ("MVAfcfs", "fcfs",   "mva",   None,        None),
    ("LDES1e4", "psprio", "ldes",  None,        LDES_SAMPLES[0]),
    ("LDES1e5", "psprio", "ldes",  None,        LDES_SAMPLES[1]),
    ("SSA1e5",  "psprio", "ssa",   None,        100000),
    # Exact-CTMC floors last. HOL and FCFSPRPRIO track queue ORDER, so their
    # state spaces blow up well before PSPRIO's does: both are out of reach at
    # N=8 and are capped so they fail fast instead of eating the run.
    ("dpsgap",  "dps",    "ctmc",  None,        None),
    ("prsgap",  "prs",    "ctmc",  None,        None),
    ("fcfsgap", "fcfs",   "ctmc",  None,        None),
    ("holgap",  "hol",    "ctmc",  None,        None),
]


def load(path):
    if os.path.exists(path):
        try:
            return json.load(open(path))
        except Exception:
            pass
    return []


def find(rows, law, N):
    for r in rows:
        if r["law"] == law and r["N"] == N:
            return r
    return None


def main():
    only_law = sys.argv[1] if len(sys.argv) > 1 else None
    only_N = int(sys.argv[2]) if len(sys.argv) > 2 else None
    comp, ans = load(COMP), load(ANSIM)
    for law, scv in LAWS:
        for N in NS:
            if only_law and law != only_law:
                continue
            if only_N and N != only_N:
                continue
            cfg = dict(BASE, law=law, N=[N, N])
            classes = ["Class1", "Class2"]
            print("\n=== %s SCV=%.1f N=%d ===" % (law, scv, N), flush=True)

            crow = find(comp, law, N)
            if crow is None:
                crow = {"law": law, "scv": scv, "N": N}
                comp.append(crow)
            arow = find(ans, law, N)
            if arow is None:
                arow = {"law": law, "scv": scv, "N": N}
                ans.append(arow)

            ref = crow.get("ref") or arow.get("ref")
            if ref is None:
                q, dt, e = A.solve(A.build_model(cfg, "psprio"), "ctmc")
                if e:
                    print("  REFERENCE FAILED:", e[:120], flush=True)
                    continue
                ref = qmat(q, classes)
                crow["ref"] = arow["ref"] = ref
            print("  ref %s" % [round(v, 4) for v in ref], flush=True)

            for name, variant, solver, method, samples in JOBS:
                key = "err_" + name
                if key in crow:
                    have = crow[key]
                    print("  %-9s skip (%s)" % (name, "have %.4f" % have if have
                                                is not None else "recorded failure"), flush=True)
                    continue
                tmo = (240 if N >= 8 else 900) if solver == "ctmc" else 600
                q, dt, e = A.solve(A.build_model(cfg, variant), solver,
                                   method=method, samples=samples,
                                   seed=SEED if samples else None, timeout=tmo)
                if e:
                    crow[key] = None
                    print("  %-9s FAILED %s" % (name, e[:90]), flush=True)
                else:
                    crow[key] = l1(qmat(q, classes), ref)
                    crow["t_" + name] = dt
                    print("  %-9s L1=%8.4f (%.1fs)" % (name, crow[key], dt), flush=True)
                json.dump(comp, open(COMP, "w"), indent=1)

            for mode in ("ctmcP", "live"):
                key = "err_ansim_" + mode
                if key in arow:
                    have = arow[key]
                    print("  ansim-%-6s skip (%s)" % (mode, "have %.4f" % have if have
                                                      is not None else "recorded failure"), flush=True)
                    continue
                flat, dt, it, e = ansim(cfg, mode)
                if e:
                    arow[key] = None
                    print("  ansim-%-6s FAILED %s" % (mode, str(e)[:90]), flush=True)
                else:
                    arow[key] = l1(flat, ref)
                    arow["t_ansim_" + mode] = dt
                    arow["it_ansim_" + mode] = it
                    arow["q_ansim_" + mode] = flat
                    print("  ansim-%-6s L1=%8.4f (%.1fs, %d it) %s"
                          % (mode, arow[key], dt, it, [round(x, 4) for x in flat]), flush=True)
                json.dump(ans, open(ANSIM, "w"), indent=1)

    print("\nRESUME DONE", flush=True)


if __name__ == "__main__":
    main()
