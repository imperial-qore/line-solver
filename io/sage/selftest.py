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

# --- /ctmc/measures: a weighted sum of pi ----------------------------------
# The premise of the endpoint is that every mean SolverCTMC reports is a
# NUMERIC linear functional of pi, so the oracles here are the closed forms of
# those functionals, never a normal form read back out of the engine.
r = ok(call(server.handle_ctmc_measures,
            {"symbols": ["x1", "x2"], "Q": [["-x1", "x1"], ["x2", "-x2"]],
             "weights": [["0", "1"], ["1", "0"], ["0", "0"]],
             "names": ["P1", "P0", "Zero"],
             "ratios": [{"name": "odds", "num": 0, "den": 1},
                        {"name": "bad", "num": 0, "den": 2}]}))
check("measures P1", r["measures"][0]["expr"], "x1/(x1 + x2)")
check("measures P0", r["measures"][1]["expr"], "x2/(x1 + x2)")
# The ratio must be REDUCED: _primitive clears integer content but a common
# polynomial factor is the field's job, so the split has to canonicalize first.
check("measures ratio reduced", r["ratios"][0]["expr"], "x1/x2")
check("measures ratio num", r["ratios"][0]["num"], "x1")
# A denominator that is identically zero has no value. Answering 0, as the
# numeric arms' isnan sweep does, would be a lie about an exact result.
check("zero denominator undefined", r["ratios"][1]["expr"], "undefined")
check_true("zero denominator reason", "identically zero" in r["ratios"][1].get("reason", ""),
           str(r["ratios"][1]))

# Mean queue length on the cross-codebase oracle. pi = [4, 6, 6, 3]/19 at
# x1 = 1, x2 = 2, and the states are ordered by the number of jobs at the
# queue, so QN = (0*4 + 1*6 + 2*6 + 3*3)/19 = 27/19.
r = ok(call(server.handle_ctmc_measures,
            {"symbols": ["x1", "x2"], "Q": Q,
             "weights": [["0", "1", "2", "3"]], "names": ["QN"]}))
e = ok(call(server.handle_eval, {"exprs": [r["measures"][0]["expr"]],
                                 "values": {"x1": "1", "x2": "2"}}))
check_true("oracle QN", abs(e["values"][0] - 27.0 / 19.0) < 1e-12, str(e["values"]))

# --- /ctmc/passage: transform and moments ----------------------------------
# One transient state absorbing at rate x1 is Exp(x1): L = x1/(x1 + s),
# E[T] = 1/x1, E[T^2] = 2/x1^2. Analytic, not read back from the engine.
r = ok(call(server.handle_ctmc_passage,
            {"symbols": ["x1"], "S": [["-x1"]], "s0": ["x1"], "alpha": ["1"],
             "svar": "s", "want": ["lst", "moments"], "nmax": 2}))
check("passage lst", r["lst"], "x1/(x1 + s)")
check("passage E[T]", r["moments"][0], "1/x1")
check("passage E[T^2]", r["moments"][1], "2/x1^2")

# First passage 0 -> 3 on the symbolic M/M/1/3 birth-death chain. At x1 = 1,
# x2 = 2 an independent numeric solve of (-S)m1 = 1 and (-S)m2 = 2*m1 gives
# 11 and 228.
Spass = [["-x1", "x1", "0"], ["x2", "-x1-x2", "x1"], ["0", "x2", "-x1-x2"]]
r = ok(call(server.handle_ctmc_passage,
            {"symbols": ["x1", "x2"], "S": Spass, "s0": ["0", "0", "x1"],
             "alpha": ["1", "0", "0"], "want": ["moments"], "nmax": 2}))
e = ok(call(server.handle_eval, {"exprs": r["moments"], "values": {"x1": "1", "x2": "2"}}))
check_true("passage bd moments", abs(e["values"][0] - 11.0) < 1e-12
           and abs(e["values"][1] - 228.0) < 1e-12, str(e["values"]))

# A state that cannot reach the target has an infinite passage time and makes
# the sub-generator singular there. Refuse by name; do not let the solve fail
# with a generic message, and never return a finite number for it.
r = call(server.handle_ctmc_passage,
         {"symbols": ["x2"], "S": [["0", "0"], ["x2", "-x2"]], "s0": ["0", "x2"],
          "alpha": ["1", "0"], "want": ["moments"]})
check("passage unreachable refused", r["status"], "error")
check("passage unreachable code", r.get("code"), "absorbing")

# Every start already inside the target: T is identically 0. Well defined, not
# an error, and the numeric arm agrees.
r = ok(call(server.handle_ctmc_passage,
            {"symbols": ["x1"], "S": [["-x1"]], "s0": ["x1"], "alpha": ["0"],
             "atom": "1", "want": ["moments"], "nmax": 2}))
check("passage atom moments", r["moments"], ["0", "0"])

# The transform symbol must not also be a rate symbol.
r = call(server.handle_ctmc_passage,
         {"symbols": ["s"], "S": [["-s"]], "s0": ["s"], "alpha": ["1"],
          "svar": "s", "want": ["lst"]})
check("passage svar collision refused", r["status"], "error")

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
