"""Priority-station comparison for the ansim surrogate, via the C++ line-cli.

Model: Delay + Queue1 under PSPRIO (only the most urgent non-empty priority
group is served, PS inside the group), D=(2,3), Z=(4,6), priorities (0,1),
N1=N2=N, service law at the queue swept over exp / hyper-exp (SCV=4) /
hypo-exp (SCV=0.5) at fixed mean. Reference: exact PSPRIO CTMC. Error is L1
over the whole station-by-class queue-length matrix.

Every competitor is priced against the exact CTMC of the discipline it is
actually forced to solve, so its approximation error is separated from the cost
of the relaxation itself:
  dps-gap   DPS-W CTMC vs PSPRIO CTMC   floor for FLD / MVA-DPS
  hol-gap   HOL CTMC   vs PSPRIO CTMC   floor for MVA-HOL
  prs-gap   PRS CTMC   vs PSPRIO CTMC   floor for MVA-PRS
  fcfs-gap  FCFS CTMC  vs PSPRIO CTMC   floor for Marie (priority-blind)
"""
import sys, json, csv, time
sys.path.insert(0, "/data/gcasale/line-dev.git/experiments/scratch/poc_ansim")
import ansim_json as A

STATIONS = ["Think", "Q"]

BASE = dict(D=[2, 3], Z=[4, 6], prio=[0, 1], Wdps=100)
LAWS = [("exp", 1.0), ("hyperexp", 4.0), ("hypoexp", 0.5)]
NS = [2, 4, 8]
LDES_SAMPLES = [10000, 100000]
SEED = 23000
MAXIT, TOL = 40, 1e-6


def qmat(q, classes):
    """Queue-length matrix as a flat list, station-major over the given classes."""
    return [q.get((st, cl), 0.0) for st in STATIONS for cl in classes]


def l1(a, b):
    return sum(abs(x - y) for x, y in zip(a, b))


def ansim(cfg, mode):
    """ansim v2 fixed point. Returns (flat M*K queue lengths, seconds, iters, err)."""
    K = len(cfg["N"])
    # start with the whole population thinking
    QN = [[float(cfg["N"][s]), 0.0] for s in range(K)]   # QN[s] = [think, queue]
    t0 = time.time()
    for it in range(1, MAXIT + 1):
        prev = [row[:] for row in QN]
        for r in range(K):
            field = [row[:] for row in QN]
            if cfg["N"][r] > 0:
                f = (cfg["N"][r] - 1.0) / cfg["N"][r]
                field[r] = [v * f for v in QN[r]]
            m = A.build_surrogate(cfg, r, field, mode)
            q, _, err = A.solve(m, "ctmc")
            if err:
                return None, time.time() - t0, it, err
            # A row absent from the table is a zero: under a frozen field a
            # starved tagged job never returns to the delay station at all.
            occ = [q.get(("Think", "Moving"), 0.0), q.get(("Q", "Moving"), 0.0)]
            s = sum(occ)
            if s <= 0:
                return None, time.time() - t0, it, "tagged occupancy vanished"
            QN[r] = [cfg["N"][r] * v / s for v in occ]
        delta = sum(abs(QN[s][i] - prev[s][i]) for s in range(K) for i in range(2))
        if delta < TOL:
            break
    flat = [QN[s][i] for i in range(2) for s in range(K)]   # station-major
    return flat, time.time() - t0, it, None


def main():
    rows = []
    for law, scv in LAWS:
        for N in NS:
            cfg = dict(BASE, law=law, N=[N, N])
            classes = ["Class1", "Class2"]
            print("\n=== %s SCV=%.1f N=%d ===" % (law, scv, N), flush=True)
            row = {"law": law, "scv": scv, "N": N}

            qref, tref, err = A.solve(A.build_model(cfg, "psprio"), "ctmc")
            if err:
                print("  REFERENCE FAILED:", err, flush=True)
                rows.append(row)
                continue
            ref = qmat(qref, classes)
            print("  reference %s (%.1fs)" % ([round(v, 4) for v in ref], tref), flush=True)

            jobs = [
                ("dpsgap",   "dps",    "ctmc",  None,  None),
                ("holgap",   "hol",    "ctmc",  None,  None),
                ("prsgap",   "prs",    "ctmc",  None,  None),
                ("fcfsgap",  "fcfs",   "ctmc",  None,  None),
                ("FLDmn",    "dps",    "fluid", "minnormal", None),
                ("MVAdps",   "dps",    "mva",   None,  None),
                ("MVAcl",    "hol",    "mva",   None,  None),
                ("MVAprs",   "prs",    "mva",   None,  None),
                ("MVAfcfs",  "fcfs",   "mva",   None,  None),
                ("LDES1e4",  "psprio", "ldes",  None,  LDES_SAMPLES[0]),
                ("LDES1e5",  "psprio", "ldes",  None,  LDES_SAMPLES[1]),
                ("SSA1e5",   "psprio", "ssa",   None,  100000),
            ]
            for name, variant, solver, method, samples in jobs:
                q, dt, e = A.solve(A.build_model(cfg, variant), solver,
                                   method=method, samples=samples,
                                   seed=SEED if samples else None)
                if e:
                    row["err_" + name] = None
                    print("  %-9s FAILED %s" % (name, e[:100]), flush=True)
                else:
                    v = l1(qmat(q, classes), ref)
                    row["err_" + name] = v
                    row["t_" + name] = dt
                    print("  %-9s L1=%8.4f (%.1fs)" % (name, v, dt), flush=True)

            for mode in ("ctmcP", "live"):
                flat, dt, it, e = ansim(cfg, mode)
                name = "ansim_" + mode
                if e:
                    row["err_" + name] = None
                    print("  %-9s FAILED %s" % (name, str(e)[:100]), flush=True)
                else:
                    v = l1(flat, ref)
                    row["err_" + name] = v
                    row["t_" + name] = dt
                    row["it_" + name] = it
                    print("  %-9s L1=%8.4f (%.1fs, %d it) %s"
                          % (name, v, dt, it, [round(x, 4) for x in flat]), flush=True)

            rows.append(row)
            json.dump(rows, open("/data/gcasale/line-dev.git/experiments/scratch/"
                                 "poc_ansim/poc_ansim_prio_cpp.json", "w"), indent=1)

    keys = ["law", "scv", "N"] + sorted({k for r in rows for k in r if k != "law"
                                         and k != "scv" and k != "N"})
    with open("/data/gcasale/line-dev.git/experiments/scratch/poc_ansim/"
              "poc_ansim_prio_cpp.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
