#!/usr/bin/env sage-python
"""
Self-test for the LINE symbolic service, run inside the image.

Exercises every handler in server.py without going through HTTP, so a failure
points at the algebra rather than at the transport. The oracle is the two
station closed model that the symbolic CTMC tests in all three codebases use:
Delay Exp(1) and FCFS Queue Exp(2) with N = 3 jobs.

    docker run --rm -v "$PWD/sage:/srv:ro" imperialqore/line-sage-rest:latest \
        sage -python /srv/selftest.py

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import server  # noqa: E402

FAILURES = []


def call(fn, payload):
    """Calls a handler the way the HTTP layer does, turning a refusal into a
    response rather than an exception."""
    try:
        return fn(payload)
    except server.LineSymError as err:
        return {"status": "error", "code": err.code, "message": err.message}


def ok(r):
    """Fails loudly with the server's message instead of a KeyError later."""
    if r.get("status") != "ok":
        raise SystemExit("unexpected error [%s]: %s" % (r.get("code"), r.get("message")))
    return r


def check(name, got, want):
    ok = got == want
    if not ok:
        FAILURES.append("%s: got %r, want %r" % (name, got, want))
    print("%-34s %s" % (name, "ok" if ok else "FAIL  got=%r want=%r" % (got, want)))


def check_true(name, cond, detail=""):
    if not cond:
        FAILURES.append("%s: %s" % (name, detail))
    print("%-34s %s" % (name, "ok" if cond else "FAIL  " + detail))


# --- two state chain -------------------------------------------------------
r = call(server.handle_ctmc_solve, {"symbols": ["x1", "x2"],
                              "Q": [["-x1", "x1"], ["x2", "-x2"]]})
check("2-state status", r["status"], "ok")
check_true("2-state pi", r["pi"] == ["x2/(x1 + x2)", "x1/(x1 + x2)"], str(r.get("pi")))
check("2-state nConnComp", r["nConnComp"], 1)

# --- the cross-codebase oracle --------------------------------------------
# Delay Exp(1) + FCFS Queue Exp(2), N = 3. States are ordered by the number of
# jobs at the queue; the generator below is the one getSymbolicGenerator
# produces after normalizing each event filtration by its minimum rate.
Q = [["-3*x1", "3*x1", "0", "0"],
     ["x2", "-x2 - 2*x1", "2*x1", "0"],
     ["0", "x2", "-x2 - x1", "x1"],
     ["0", "0", "x2", "-x2"]]
r = call(server.handle_ctmc_solve, {"symbols": ["x1", "x2"], "Q": Q})
check("oracle status", r["status"], "ok")
e = ok(call(server.handle_eval, {"exprs": r["pi"], "values": {"x1": "1", "x2": "2"}}))
want = [0.21052631578947367, 0.3157894736842105, 0.3157894736842105, 0.15789473684210528]
got = e["values"]
check_true("oracle pi at x1=1,x2=2",
           all(abs(a - b) < 1e-12 for a, b in zip(got, want)), str(got))
check_true("oracle pi sums to 1", abs(sum(got) - 1.0) < 1e-12, str(sum(got)))
check("oracle exact pi[0]", e["exact"][0], "4/19")
# Every entry must print over the same denominator: a rational constant is a
# unit in the fraction field, so without the primitive normal form the same
# probability comes back scaled differently from one state to the next.
dens = [p.split("/", 1)[1] for p in r["pi"]]
check_true("oracle shares one denominator", len(set(dens)) == 1, str(dens))
check("oracle pi[3]", r["pi"][3],
      "6*x1^3/(6*x1^3 + 6*x1^2*x2 + 3*x1*x2^2 + x2^3)")
check("oracle den", r["den"], "6*x1^3 + 6*x1^2*x2 + 3*x1*x2^2 + x2^3")
check("oracle num", r["num"], ["x2^3", "3*x1*x2^2", "6*x1^2*x2", "6*x1^3"])

# --- exact decimal reading -------------------------------------------------
# 0.1 must become 1/10, not the binary double nearest to it.
r = call(server.handle_ctmc_solve, {"symbols": ["x1"],
                              "Q": [["-0.1*x1", "0.1*x1"], ["0.3*x1", "-0.3*x1"]]})
check("decimal is exact", r["pi"], ["3/4", "1/4"])

# --- reducible generator ---------------------------------------------------
# Two disjoint two-state chains: each is solved on its own and the pair is
# renormalized, which is the initial-vector-free MATLAB behaviour.
r = call(server.handle_ctmc_solve, {"symbols": ["x1"],
                              "Q": [["-x1", "x1", "0", "0"],
                                    ["x1", "-x1", "0", "0"],
                                    ["0", "0", "-x1", "x1"],
                                    ["0", "0", "x1", "-x1"]]})
check("reducible nConnComp", r["nConnComp"], 2)
check("reducible pi", r["pi"], ["1/4", "1/4", "1/4", "1/4"])

# --- absorbing generator is refused, not fabricated ------------------------
r = call(server.handle_ctmc_solve, {"symbols": ["x1"], "Q": [["-x1", "x1"], ["0", "0"]]})
check("absorbing is an error", r["status"], "error")
check("absorbing error code", r.get("code"), "absorbing")

# --- sensitivity -----------------------------------------------------------
# d/dx1 of x1/(x1+x2) is x2/(x1+x2)^2, which at x1=1, x2=2 is 2/9.
r = call(server.handle_ctmc_sensitivity, {"symbols": ["x1", "x2"],
                                    "Q": [["-x1", "x1"], ["x2", "-x2"]],
                                    "theta": "x1",
                                    "reward": ["0", "1"]})
check("sensitivity status", r["status"], "ok")
e = ok(call(server.handle_eval, {"exprs": [r["S"]], "values": {"x1": "1", "x2": "2"}}))
check_true("dE[r]/dx1 at (1,2)", abs(e["values"][0] - 2.0 / 9.0) < 1e-12, str(e["values"]))

# --- simplify, diff, fluid -------------------------------------------------
r = call(server.handle_simplify, {"exprs": ["(x1^2 - 1)/(x1 - 1)"], "form": "cancel"})
check("cancel", r["results"], ["x1 + 1"])
r = call(server.handle_simplify, {"exprs": ["x1^2 - 1"], "form": "factor"})
check("factor", r["results"], ["(x1 + 1)*(x1 - 1)"])
r = call(server.handle_diff, {"exprs": ["x1^3"], "var": "x1", "order": 2})
check("diff order 2", r["results"], ["6*x1"])
r = call(server.handle_fluid_odes, {"rhs": ["-a*x + b*y", "a*x - b*y"], "vars": ["x", "y"],
                              "want": ["jacobian", "latex", "equilibria"]})
check("jacobian", r["jacobian"], [["-a", "b"], ["a", "-b"]])
check_true("fluid latex", r["latex"][0].find("\\") >= 0 or "x" in r["latex"][0],
           str(r.get("latex")))

# --- fractional exponents, which int() would silently truncate -------------
# (1 + t^8)^(1/8) with the exponent truncated to 0 collapses to 1, and the
# p-norm smoothed fluid drift then comes back wrong by a thousandth with
# nothing reported.
r = ok(call(server.handle_eval, {"exprs": ["x/(1 + (x/4)^8)^(1/8)"],
                                 "values": {"x": "0.9"}}))
check_true("fractional exponent", abs(r["values"][0] - 0.899999260702) < 1e-9,
           str(r["values"]))
r = call(server.handle_ctmc_solve, {"symbols": ["x1"],
                                    "Q": [["-x1^(1/2)", "x1^(1/2)"],
                                          ["x1", "-x1"]]})
check("fractional exponent refused in a field", r["status"], "error")

# --- min() in the drift, which no polynomial ring can hold -----------------
r = call(server.handle_fluid_odes, {"rhs": ["-min(x, 1)"], "vars": ["x"], "want": ["latex"]})
check("min in drift", r["status"], "ok")

# --- the parser refuses what it cannot evaluate safely ---------------------
r = call(server.handle_simplify, {"exprs": ["__import__('os').system('id')"], "form": "cancel"})
check("no eval escape", r["status"], "error")
r = call(server.handle_ctmc_solve, {"symbols": ["x1"], "Q": [["-x1", "x9"], ["x1", "-x1"]]})
check("undeclared symbol refused", r["status"], "error")

# --- readiness -------------------------------------------------------------
r = call(server.handle_ready, {})
check("ready", r["ready"], True)

print("")
if FAILURES:
    print("%d FAILED" % len(FAILURES))
    for f in FAILURES:
        print("  " + f)
    sys.exit(1)
print("all checks passed")
