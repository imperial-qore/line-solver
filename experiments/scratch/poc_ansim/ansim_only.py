"""ansim fixed point only, re-run after the self-loop routing fix.
Runs alongside the competitor sweep; results merge on (law, N)."""
import sys, json, time
sys.path.insert(0, "/data/gcasale/line-dev.git/experiments/scratch/poc_ansim")
import ansim_json as A
from poc_ansim_prio import ansim, qmat, l1, BASE, LAWS, NS

out = []
for law, scv in LAWS:
    for N in NS:
        cfg = dict(BASE, law=law, N=[N, N])
        qref, _, err = A.solve(A.build_model(cfg, "psprio"), "ctmc")
        if err:
            print(f"{law} N={N} REF FAIL {err[:90]}", flush=True); continue
        ref = qmat(qref, ["Class1", "Class2"])
        rec = {"law": law, "scv": scv, "N": N, "ref": ref}
        print(f"\n=== {law} SCV={scv} N={N} ===  ref={[round(v,4) for v in ref]}", flush=True)
        for mode in ("ctmcP", "live"):
            flat, dt, it, e = ansim(cfg, mode)
            if e:
                rec["err_ansim_" + mode] = None
                print(f"  ansim-{mode:6s} FAILED {str(e)[:90]}", flush=True)
            else:
                v = l1(flat, ref)
                rec["err_ansim_" + mode] = v
                rec["t_ansim_" + mode] = dt
                rec["it_ansim_" + mode] = it
                rec["q_ansim_" + mode] = flat
                print(f"  ansim-{mode:6s} L1={v:8.4f} ({dt:.1f}s, {it} it) {[round(x,4) for x in flat]}", flush=True)
        out.append(rec)
        json.dump(out, open("/data/gcasale/line-dev.git/experiments/scratch/poc_ansim/ansim_only.json","w"), indent=1)
print("\nANSIM DONE", flush=True)
