/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The --api invocation path: the JSON conversion policy, the dispatch table
 * and every refusal it owes a caller.
 *
 * Most of this exercises line::reg::api_invoke directly, which is the entry
 * point main() and a future pybind11 binding both call, so the test covers the
 * code that ships rather than a re-implementation of it. The end-to-end cases
 * go through the built binary, because the exit code and the stderr message on
 * a refusal are part of the contract a subprocess caller depends on and they
 * exist nowhere else.
 *
 * The numeric fixtures are MATLAB-verified values, not values read back out of
 * this port.
 */
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>

#include "doctest.h"
#include "line/reg/api_dispatch.h"
#include "line/solvers/wrappers/lqns/lqns_probe.h"
#include "line/util/method_type.h"

using line::reg::api_invoke;
using line::reg::Json;
using line::reg::parse_arith;

namespace {

// MATLAB-verified: pfqn_ca(L=[0.6 0.4], N=[2 1], Z=[1 0.5]).
constexpr double CA_LG = 0.610851937833;
// MATLAB-verified: pfqn_comomrm_ms(same, m=1, S=2).
constexpr double COMOMRM_MS_LG = 0.17227122094;
// MATLAB-verified: pfqn_conv(L=[0.6 0.4; 0.3 0.7], N=[2 1], Z=[1 0.5]).
constexpr double CONV_LG = 1.46093790412;

Json ca_args() {
    return Json::parse(R"({"L": [[0.6, 0.4]], "N": [2, 1], "Z": [1, 0.5]})");
}

/** Run the CLI binary, capturing stdout and the exit status. */
struct Run {
    int status = 0;
    std::string out;
};

Run run_cli(const std::string& args) {
    const std::string cmd = std::string(LINE_MP_CLI_BINARY) + " " + args + " 2>&1";
    Run r;
    FILE* p = popen(cmd.c_str(), "r");
    REQUIRE(p != nullptr);
    char buf[512];
    while (std::fgets(buf, sizeof(buf), p) != nullptr) r.out += buf;
    const int rc = pclose(p);
    r.status = WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
    return r;
}

std::string write_temp(const std::string& name, const std::string& text) {
    const std::string path = std::string("/tmp/line_mp_test_") + name;
    std::ofstream out(path.c_str());
    out << text;
    out.close();
    return path;
}

}  // namespace

TEST_CASE("--arith parses the three accepted forms and refuses the rest") {
    CHECK(parse_arith("double").str() == "double");
    CHECK(parse_arith("exact").str() == "exact");
    CHECK(parse_arith("real:50").str() == "real:50");
    CHECK(parse_arith("real:100").str() == "real:100");
    CHECK(parse_arith("real:200").str() == "real:200");

    // Between tiers, round UP: a caller never receives less precision than asked.
    CHECK(parse_arith("real:60").str() == "real:100");
    CHECK(parse_arith("real:1").str() == "real:50");

    CHECK_THROWS_AS(parse_arith("float"), line::InputError);
    CHECK_THROWS_AS(parse_arith("real:"), line::InputError);
    CHECK_THROWS_AS(parse_arith("real:abc"), line::InputError);
    CHECK_THROWS_AS(parse_arith("exact:50"), line::InputError);
    // Wider than any instantiated tier: refused, not silently downgraded.
    CHECK_THROWS_AS(parse_arith("real:400"), line::UnsupportedError);
}

TEST_CASE("pfqn_ca over --api at double matches the MATLAB fixture") {
    const Json r = api_invoke("pfqn_ca", "double", ca_args());
    CHECK(r["function"] == "pfqn_ca");
    CHECK(r["arith"] == "double");
    CHECK(r["results"]["lG"].get<double>() == doctest::Approx(CA_LG).epsilon(1e-11));
    // At double, a scalar crosses as a bare number.
    CHECK(r["results"]["G"].is_number());
    CHECK(r["results"]["G"].get<double>() == doctest::Approx(std::exp(CA_LG)).epsilon(1e-11));
}

TEST_CASE("pfqn_ca at exact returns num/den whose ratio is the double") {
    const Json r = api_invoke("pfqn_ca", "exact", ca_args());
    CHECK(r["arith"] == "exact");
    const Json& G = r["results"]["G"];
    REQUIRE(G.is_object());
    const std::string num = G["num"].get<std::string>();
    const std::string den = G["den"].get<std::string>();
    const double ratio = std::stod(num) / std::stod(den);
    CHECK(ratio == doctest::Approx(G["double"].get<double>()).epsilon(1e-15));
    CHECK(ratio == doctest::Approx(std::exp(CA_LG)).epsilon(1e-11));

    // The decimal-literal policy: 0.6 in the JSON is the rational 3/5, so the
    // constant is 921/500 exactly, as pinned in test_pfqn_comom_family.cpp.
    // Reading 0.6 as the nearest dyadic rational instead would give the same
    // double and a meaningless fraction.
    CHECK(num == "921");
    CHECK(den == "500");

    // lG is a C++ double whatever the arithmetic: it is finite where G is not
    // representable, so it crosses as a bare number and is never num/den.
    CHECK(r["results"]["lG"].is_number());
    CHECK(r["results"]["lG"].get<double>() == doctest::Approx(CA_LG).epsilon(1e-11));
}

TEST_CASE("pfqn_ca at real:50 carries a decimal string alongside the double") {
    const Json r = api_invoke("pfqn_ca", "real:50", ca_args());
    CHECK(r["arith"] == "real:50");
    const Json& G = r["results"]["G"];
    REQUIRE(G.is_object());
    CHECK(G["double"].get<double>() == doctest::Approx(std::exp(CA_LG)).epsilon(1e-11));
    const std::string dec = G["dec"].get<std::string>();
    CHECK(dec.size() > 20);  // 50 digits, not a double reprinted
    CHECK(std::stod(dec) == doctest::Approx(std::exp(CA_LG)).epsilon(1e-11));
}

TEST_CASE("pfqn_conv over --api matches the MATLAB fixture") {
    const Json args =
        Json::parse(R"({"L": [[0.6, 0.4], [0.3, 0.7]], "N": [2, 1], "Z": [1, 0.5]})");
    CHECK(api_invoke("pfqn_conv", "double", args)["results"]["lG"].get<double>() ==
          doctest::Approx(CONV_LG).epsilon(1e-11));
    CHECK(api_invoke("pfqn_conv", "exact", args)["results"]["lG"].get<double>() ==
          doctest::Approx(CONV_LG).epsilon(1e-11));
}

TEST_CASE("pfqn_comomrm_ms over --api matches the MATLAB fixture") {
    const Json args =
        Json::parse(R"({"L": [[0.6, 0.4]], "N": [2, 1], "Z": [1, 0.5], "m": 1, "S": 2})");
    const Json r = api_invoke("pfqn_comomrm_ms", "exact", args);
    CHECK(r["results"]["lG"].get<double>() == doctest::Approx(COMOMRM_MS_LG).epsilon(1e-10));
    // The queue-length marginal is a probability vector.
    const Json& prob = r["results"]["prob"];
    REQUIRE(prob.is_array());
    double total = 0.0;
    for (const Json& p : prob) total += p["double"].get<double>();
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("pfqn_mva over --api returns every output of the MATLAB signature") {
    const Json r = api_invoke("pfqn_mva", "double", ca_args());
    const Json& res = r["results"];
    CHECK(res["XN"].is_array());
    CHECK(res["QN"].is_array());
    CHECK(res["UN"].is_array());
    CHECK(res["CN"].is_array());
    CHECK(res["lG"].get<double>() == doctest::Approx(CA_LG).epsilon(1e-11));
    // Population conservation: sum_i Q_i + X Z = N, per class.
    REQUIRE(res["QN"].size() == 1);
    for (std::size_t r_i = 0; r_i < 2; ++r_i) {
        const double Z[2] = {1.0, 0.5};
        const double N[2] = {2.0, 1.0};
        const double q = res["QN"][0][r_i].get<double>();
        const double x = res["XN"][r_i].get<double>();
        CHECK(q + x * Z[r_i] == doctest::Approx(N[r_i]).epsilon(1e-10));
    }
}

TEST_CASE("ctmc_solve and dtmc_solve over --api") {
    // Two-state chain: rates 1 and 3, stationary law (3/4, 1/4).
    const Json q = Json::parse(R"({"Q": [[-1, 1], [3, -3]]})");
    const Json r = api_invoke("ctmc_solve", "exact", q);
    REQUIRE(r["results"]["p"].size() == 2);
    CHECK(r["results"]["p"][0]["num"].get<std::string>() == "3");
    CHECK(r["results"]["p"][0]["den"].get<std::string>() == "4");
    CHECK(r["results"]["p"][1]["num"].get<std::string>() == "1");
    CHECK(r["results"]["p"][1]["den"].get<std::string>() == "4");

    const Json p = Json::parse(R"({"P": [[0.5, 0.5], [0.25, 0.75]]})");
    // pi P = pi with pi = (1/3, 2/3): 0.5/3 + 0.25*2/3 = 1/3.
    const Json d = api_invoke("dtmc_solve", "exact", p);
    REQUIRE(d["results"]["PROB"].size() == 2);
    CHECK(d["results"]["PROB"][0]["num"].get<std::string>() == "1");
    CHECK(d["results"]["PROB"][0]["den"].get<std::string>() == "3");

    const Json g = api_invoke("ctmc_makeinfgen", "exact",
                              Json::parse(R"({"Q": [[0, 1], [3, 0]]})"));
    CHECK(g["results"]["Q"][0][0]["num"].get<std::string>() == "-1");
}

TEST_CASE("a function that is registered but not wired is refused by name") {
    // ctmc_stochcomp is ported, tested and in the registry, but has no dispatch
    // entry: it must say so, not return an empty result.
    REQUIRE(line::find_api("ctmc_stochcomp") != nullptr);
    REQUIRE(!line::reg::api_is_exposed("ctmc_stochcomp"));
    try {
        api_invoke("ctmc_stochcomp", "double", Json::parse(R"({"Q": [[-1, 1], [3, -3]]})"));
        FAIL("ctmc_stochcomp must be refused");
    } catch (const line::UnsupportedError& e) {
        const std::string msg = e.what();
        CHECK(msg.find("ctmc_stochcomp") != std::string::npos);
        CHECK(msg.find("not yet exposed") != std::string::npos);
    }
}

TEST_CASE("a function that is not ported at all is refused as such") {
    REQUIRE(line::find_api("pfqn_mom") == nullptr);
    try {
        api_invoke("pfqn_mom", "double", ca_args());
        FAIL("pfqn_mom must be refused");
    } catch (const line::UnsupportedError& e) {
        const std::string msg = e.what();
        CHECK(msg.find("pfqn_mom") != std::string::npos);
        CHECK(msg.find("not ported") != std::string::npos);
    }
}

TEST_CASE("an arithmetic the function does not support is refused by name") {
    // pfqn_marie needs logs and is registered for double and real only.
    const line::ApiEntry* e = line::find_api("pfqn_marie");
    REQUIRE(e != nullptr);
    REQUIRE(!line::api_supports(*e, line::Arith::Exact));
    try {
        api_invoke("pfqn_marie", "exact", Json::object());
        FAIL("pfqn_marie at exact must be refused");
    } catch (const line::UnsupportedError& err) {
        const std::string msg = err.what();
        CHECK(msg.find("pfqn_marie") != std::string::npos);
        CHECK(msg.find("does not support") != std::string::npos);
        // The supported modes come from the registry entry, not a hardcoded list.
        CHECK(msg.find("double") != std::string::npos);
        CHECK(msg.find("real") != std::string::npos);
    }
}

TEST_CASE("argument errors are refused rather than defaulted") {
    // An unknown key would otherwise take its default silently.
    CHECK_THROWS_AS(api_invoke("pfqn_ca", "double",
                               Json::parse(R"({"L": [[0.6]], "N": [1], "ZZ": [1]})")),
                    line::InputError);
    // A missing required argument.
    CHECK_THROWS_AS(api_invoke("pfqn_ca", "double", Json::parse(R"({"L": [[0.6]]})")),
                    line::InputError);
    // A non-integral population.
    CHECK_THROWS_AS(api_invoke("pfqn_ca", "double", Json::parse(R"({"L": [[0.6]], "N": [1.5]})")),
                    line::InputError);
    // A ragged matrix.
    CHECK_THROWS_AS(
        api_invoke("pfqn_ca", "double", Json::parse(R"({"L": [[0.6, 0.4], [0.3]], "N": [1, 1]})")),
        line::InputError);
    // An argument the JSON boundary cannot carry is named, not ignored.
    CHECK_THROWS_AS(api_invoke("pfqn_conv", "double",
                               Json::parse(R"({"L": [[0.6]], "N": [1], "cdscaling": [1]})")),
                    line::UnsupportedError);
    // The argument payload has to be an object keyed by parameter names.
    CHECK_THROWS_AS(api_invoke("pfqn_ca", "double", Json::parse("[1, 2, 3]")), line::InputError);
}

TEST_CASE("numeric strings are taken literally at every arithmetic") {
    // "3/5" and 0.6 are the same value under the decimal policy.
    const Json byString =
        api_invoke("pfqn_ca", "exact",
                   Json::parse(R"({"L": [["3/5", "2/5"]], "N": [2, 1], "Z": [1, "1/2"]})"));
    CHECK(byString["results"]["G"]["num"].get<std::string>() == "921");
    CHECK(byString["results"]["G"]["den"].get<std::string>() == "500");

    // A value with no short decimal form: 1/3 is exact as a string.
    const Json third = api_invoke("pfqn_ca", "exact",
                                 Json::parse(R"({"L": [["1/3"]], "N": [1]})"));
    CHECK(third["results"]["G"]["num"].get<std::string>() == "1");
    CHECK(third["results"]["G"]["den"].get<std::string>() == "3");
}

TEST_CASE("a decimal below one is decimal, not octal") {
    // Regression: the big-integer constructor reads a leading 0 as an OCTAL
    // prefix, so the digit string "025" assembled from 0.25 arrived as 21.
    // Integers and decimals whose digits are all below 8 were unaffected, which
    // is why pfqn_ca on 0.6 looked right while dtmc_solve on 0.25 did not.
    // ctmc_solve of P - I for P = [[0.5, 0.5], [0.25, 0.75]] is (1/3, 2/3).
    const Json r = api_invoke("ctmc_solve", "exact",
                              Json::parse(R"({"Q": [[-0.5, 0.5], [0.25, -0.25]]})"));
    CHECK(r["results"]["p"][0]["num"].get<std::string>() == "1");
    CHECK(r["results"]["p"][0]["den"].get<std::string>() == "3");
    CHECK(r["results"]["p"][1]["num"].get<std::string>() == "2");
    CHECK(r["results"]["p"][1]["den"].get<std::string>() == "3");

    // Digits 8 and 9 after a leading zero are not even valid octal, so this
    // would have thrown rather than answered: 0.09 = 9/100.
    const Json g = api_invoke("ctmc_makeinfgen", "exact",
                              Json::parse(R"({"Q": [[0, 0.09], [0.008, 0]]})"));
    CHECK(g["results"]["Q"][0][1]["num"].get<std::string>() == "9");
    CHECK(g["results"]["Q"][0][1]["den"].get<std::string>() == "100");
    CHECK(g["results"]["Q"][1][0]["num"].get<std::string>() == "1");
    CHECK(g["results"]["Q"][1][0]["den"].get<std::string>() == "125");
}

TEST_CASE("the CLI reports results and refusals with the right exit status") {
    const std::string args = write_temp("ca.json", R"({"L": [[0.6, 0.4]], "N": [2, 1],
                                                       "Z": [1, 0.5]})");

    const Run ok = run_cli("--api pfqn_ca --args " + args + " -o json");
    CHECK(ok.status == 0);
    const Json parsed = Json::parse(ok.out);
    CHECK(parsed["results"]["lG"].get<double>() == doctest::Approx(CA_LG).epsilon(1e-11));

    const Run readable = run_cli("--api pfqn_ca --args " + args + " --arith exact");
    CHECK(readable.status == 0);
    CHECK(readable.out.find("921 / 500") != std::string::npos);

    // Arguments on standard input when --args is absent.
    const Run piped = run_cli("--api pfqn_ca -o json < " + args);
    CHECK(piped.status == 0);
    CHECK(Json::parse(piped.out)["results"]["lG"].get<double>() ==
          doctest::Approx(CA_LG).epsilon(1e-11));

    // Every refusal is nonzero and names what went wrong.
    const Run unwired = run_cli("--api ctmc_stochcomp --args " + args);
    CHECK(unwired.status != 0);
    CHECK(unwired.out.find("ctmc_stochcomp") != std::string::npos);

    const Run badArith = run_cli("--api pfqn_ca --args " + args + " --arith quad");
    CHECK(badArith.status != 0);
    CHECK(badArith.out.find("double, exact, real:<digits>") != std::string::npos);

    const std::string broken = write_temp("broken.json", R"({"L": [[0.6, 0.4], "N")");
    const Run malformed = run_cli("--api pfqn_ca --args " + broken);
    CHECK(malformed.status != 0);
    CHECK(malformed.out.find("malformed") != std::string::npos);

    const Run missing = run_cli("--api pfqn_ca --args /tmp/line_mp_test_does_not_exist.json");
    CHECK(missing.status != 0);
    CHECK(missing.out.find("cannot open") != std::string::npos);

    std::remove(args.c_str());
    std::remove(broken.c_str());
}

// 2026-07-28: `-s ssa` ignored `--method`. solve_model_ssa called
// solver_ssa_nrm_analyzer directly, so every method ran the NRM. SolverSSA.m:55
// lists 'serial', 'para' and 'parallel' as valid methods and
// solver_ssa_analyzer.m:136-176 dispatches them, so the CLI routes through
// ssa::solver_ssa and reports the engine that actually ran.
TEST_CASE("-s ssa honours --method and reports the engine that actually ran") {
    const std::string model = write_temp("ssa_cqn.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "ssa_cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})");
    const std::string base = "-s ssa -f " + model + " --samples 2000 --seed 4242";

    const Run def = run_cli(base);
    CHECK(def.status == 0);
    CHECK(def.out.find("method=nrm") != std::string::npos);

    // The method that used to be swallowed: a different estimator entirely.
    const Run serial = run_cli(base + " --method serial");
    CHECK(serial.status == 0);
    CHECK(serial.out.find("method=serial") != std::string::npos);

    const Run alias = run_cli(base + " --method ssa");
    CHECK(alias.status == 0);
    CHECK(alias.out.find("method=serial") != std::string::npos);

    // The banner must report what RAN, not what was asked for: this model is
    // NRM-eligible, so the dispatcher answers 'parallel' with the NRM and the
    // CLI has to say 'nrm'.
    const Run par = run_cli(base + " --method parallel");
    CHECK(par.status == 0);
    CHECK(par.out.find("method=nrm") != std::string::npos);

    const Run bad = run_cli(base + " --method gillespie");
    CHECK(bad.status != 0);
    CHECK(bad.out.find("gillespie") != std::string::npos);

    std::remove(model.c_str());
}

// 2026-07-28: pinning what the SSA path had lost. Every solver's runAnalyzer
// validates its method against listValidMethods before analysing
// (matlab/src/solvers/Solver.m checkOptions, and each solver's own
// listValidMethods, e.g. matlab/src/solvers/MVA/@SolverMVA/SolverMVA.m:43), so
// an unknown name is an error in the reference and must not be answered here
// either; and
// the reference reports `self.result.method`, the algorithm that RAN, which is
// why MATLAB prints 'default/comom' rather than 'default'. Both were unpinned
// on every path but SSA, and a regression to an engine entry would restore the
// exact defect fixed above without failing anything.
TEST_CASE("-s ag takes --max-states and refuses the knobs AgOptions does not have") {
    const std::string mm1 = write_temp("ag_knobs_mm1.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "mm1",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 2.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Sink", "type": "Sink"}],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue": 1.0}, "Queue": {"Sink": 1.0}}}}}})");

    // --max-states is the truncation level of an OPEN agent's queue-length
    // dimension, so it is honoured where there is an agent to truncate...
    const Run ok = run_cli("-s ag --max-states 40 " + mm1);
    CHECK(ok.status == 0);
    CHECK(ok.out.find("SolverAG") != std::string::npos);

    // ...and refused where there is not, rather than accepted and dropped: a
    // caller who believes a truncation applied when nothing truncated has been
    // told the wrong thing about their own answer.
    const Run elsewhere = run_cli("-s mva --max-states 40 " + mm1);
    CHECK(elsewhere.status != 0);
    CHECK(elsewhere.out.find("max-states") != std::string::npos);

    // AgOptions carries tol and iter_max and NO iter_tol -- the same shape as
    // MamOptions -- so the flag is named rather than quietly ignored.
    const Run itol = run_cli("-s ag --iter_tol 1e-6 " + mm1);
    CHECK(itol.status != 0);
    CHECK(itol.out.find("iter_tol") != std::string::npos);

    // A truncation level is a state COUNT; zero and negatives are input errors.
    const Run zero = run_cli("-s ag --max-states 0 " + mm1);
    CHECK(zero.status != 0);
}

TEST_CASE("every CLI solver path validates --method and reports the one that ran") {
    const std::string cqn = write_temp("sweep_cqn.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})");
    const std::string mm1 = write_temp("sweep_mm1.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "mm1",
        "nodes": [
          {"name": "Source", "type": "Source",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 2.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "FCFS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Sink", "type": "Sink"}],
        "classes": [{"name": "Class1", "type": "Open"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Source": {"Queue": 1.0}, "Queue": {"Sink": 1.0}}}}}})");

    // An unknown method is refused with a nonzero status that names it, on
    // every path the CLI owns. `-s ssa` is covered by the case above.
    const std::string paths[][2] = {{"-s mva", cqn},         {"-s nc", cqn},
                                    {"-s mam", mm1},         {"-s ba", cqn},
                                    {"-s ctmc", cqn},        {"-s fluid", cqn},
                                    {"-s fluid -a odes", cqn}, {"-s mva -a prob", cqn},
                                    {"-s nc -a prob", cqn}};
    for (const std::string* p : paths) {
        const Run bad = run_cli(p[0] + " -f " + p[1] + " --method nosuchmethod");
        INFO("path ", p[0]);
        CHECK(bad.status != 0);
        CHECK(bad.out.find("nosuchmethod") != std::string::npos);
    }

    // The banner names the algorithm the dispatcher RESOLVED, never the knob:
    // a banner echoing k.method would read `method=default` (or an empty field)
    // on all five of these. SolverCTMC is excluded because `default` is the
    // genuine name of its only algorithm, so there the two coincide.
    const std::string resolved[][2] = {{"-s mva", cqn},   {"-s nc", cqn},    {"-s mam", mm1},
                                       {"-s ba", cqn},    {"-s fluid", cqn}};
    for (const std::string* p : resolved) {
        const Run def = run_cli(p[0] + " -f " + p[1]);
        INFO("path ", p[0]);
        CHECK(def.status == 0);
        CHECK(def.out.find("method=") != std::string::npos);
        CHECK(def.out.find("method=default ") == std::string::npos);
        CHECK(def.out.find("method= ") == std::string::npos);
    }

    std::remove(cqn.c_str());
    std::remove(mm1.c_str());
}

// The layered path was a second executable (line-ln) until it was folded into
// this binary. What that merge can silently break is not the LN arithmetic --
// test_ln_ofbiz.cpp covers that directly -- but the ROUTING to it: an input
// format that stops being recognised, or a Network solver method name that starts
// being accepted for a LayeredNetwork and answers with numbers from the wrong
// engine. Both are invisible to every other test, so they are pinned here.
//
// THE `auto` ARM IS MACHINE-DEPENDENT, as it is in the reference: LINE ships no
// lqns binary and `chooseAvgSolverHeur.m` gates its LQNS candidate on
// `SolverLQNS.isAvailable()`, so a bare `.lqnx` path reaches SolverLN only
// where lqns is absent. Every case below that means "the LN engine" therefore
// names `-s ln` explicitly; the cases that are about the `auto` arm itself are
// asserted against the probe.
TEST_CASE("the merged CLI routes a layered model to SolverLN") {
    const std::string ofbiz =
        std::string(LINE_MP_REPO_ROOT) + "/matlab/examples/basic/layeredModel/lqn_ofbiz.xml";
    const bool lqns_here = line::lqns::lqns_is_available();

    // A .xml path with no -i at all: the extension alone selects the layered
    // reader, which is how line-ln was invoked and what bench_ln_ofbiz.sh does.
    const Run sniffed = run_cli(ofbiz);
    CHECK(sniffed.status == 0);
    if (lqns_here) {
        CHECK(sniffed.out.find("SolverAUTO selected lqns") != std::string::npos);
        CHECK(sniffed.out.find("SolverLQNS(") != std::string::npos);
    } else {
        CHECK(sniffed.out.find("SolverLN(SolverMVA)") != std::string::npos);
        CHECK(sniffed.out.find("layers=14") != std::string::npos);
        CHECK(sniffed.out.find("converged=1") != std::string::npos);
    }

    // -s ln names the LN engine whatever is installed.
    const Run ln = run_cli("-s ln " + ofbiz);
    CHECK(ln.status == 0);
    CHECK(ln.out.find("SolverLN(SolverMVA)") != std::string::npos);
    CHECK(ln.out.find("layers=14") != std::string::npos);
    CHECK(ln.out.find("converged=1") != std::string::npos);

    // -i lqnx and -s ln name the same path explicitly.
    const Run named = run_cli("-i lqnx -s ln -f " + ofbiz);
    CHECK(named.status == 0);
    CHECK(named.out.find("SolverLN(SolverMVA)") != std::string::npos);

    // -o layers dumps the ensemble instead of the AvgTable.
    const Run layers = run_cli("-o layers -s ln " + ofbiz);
    CHECK(layers.status == 0);
    CHECK(layers.out.find("LAYER 1 ") != std::string::npos);
    CHECK(layers.out.find("STATION ") != std::string::npos);
    CHECK(layers.out.find("SolverLN(SolverMVA)") == std::string::npos);
    // ...and lqns has no ensemble to dump, so it says so rather than answering.
    if (lqns_here) {
        const Run layers_lqns = run_cli("-o layers -s lqns " + ofbiz);
        CHECK(layers_lqns.status != 0);
        CHECK(layers_lqns.out.find("-o layers") != std::string::npos);
    }

    // ln.comom asks for NC LAYERS, which the port gained on 2026-07-31. The bare
    // `ln` token deliberately keeps its MVA layers rather than following the
    // Java CLI, so the two tokens must resolve to DIFFERENT engines and the
    // banner has to say which one ran.
    const Run comom = run_cli("-s ln.comom -f " + ofbiz);
    CHECK(comom.status == 0);
    CHECK(comom.out.find("SolverLN(SolverNC)") != std::string::npos);

    // `-s lqns` runs the binary where there is one and is refused BY NAME where
    // there is not -- never served by the MVA layers under another name.
    {
        const Run bad = run_cli("-s lqns " + ofbiz);
        if (lqns_here) {
            CHECK(bad.status == 0);
            CHECK(bad.out.find("SolverLQNS(") != std::string::npos);
        } else {
            CHECK(bad.status != 0);
            CHECK(bad.out.find("lqns") != std::string::npos);
        }
        CHECK(bad.out.find("SolverLN(") == std::string::npos);
    }

    // A Network solver method name is not silently applied to a LayeredNetwork.
    const Run wrong = run_cli("-s mva " + ofbiz);
    CHECK(wrong.status != 0);

    // A Network-only knob is refused rather than dropped, the same discipline
    // the json path applies to every knob its chosen solver lacks -- and on
    // EITHER layered engine, so which one `auto` picked cannot decide whether
    // the knob is honoured or quietly ignored.
    const Run knob = run_cli("--cutoff 4 -s ln " + ofbiz);
    CHECK(knob.status != 0);
    CHECK(knob.out.find("cutoff") != std::string::npos);
    const Run knob_auto = run_cli("--cutoff 4 " + ofbiz);
    CHECK(knob_auto.status != 0);
    CHECK(knob_auto.out.find("cutoff") != std::string::npos);

    // ...and symmetrically, a layered knob does not reach a Network solve.
    const std::string cqn = write_temp("lqnmerge_cqn.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})");
    const Run layerknob = run_cli("-s mva -f " + cqn + " --layer-solver fluid");
    CHECK(layerknob.status != 0);
    CHECK(layerknob.out.find("layer-solver") != std::string::npos);
    std::remove(cqn.c_str());
}

// 2026-07-28: nothing referenced util::method_type, so the banner's type= field
// was unpinned on every path. A method with no method name of its own silently takes
// its solver's default label, which turns an honest 'approximate,randomized'
// into 'exact,deterministic' or the reverse without failing anything. The
// reference classifies the RUNTIME-RESOLVED method for this reason
// (matlab/src/solvers/Solver.m:78-89, which names 'default/imci' as the case),
// on the axis matlab/src/solvers/Solver.m:92 isStochasticMethod defines;
// matlab/src/solvers/NC/@SolverNC/SolverNC.m:174 fixes NC's stochastic set as
// {mci, imci, ls, sampling, is}, which is what the NC rows below assert. The
// two-axis label itself is matlab/src/io/line_method_type.m, ported here as
// util::method_type alongside jline/solvers/MethodType.java and the
// _METHOD_TYPES registry of python/line_solver/solvers/base.py:41-45. Those
// three arrived in 5f07ef92f, which is on master and NOT on worktree-mp-cpp,
// so they resolve from master and not from a checkout of this branch.
TEST_CASE("method_type classifies every method the CLI can dispatch") {
    using line::util::method_type;

    // The discriminating cases: each of these differs from what its solver's
    // own default label would give, so a missing token shows up here.
    CHECK(method_type("SolverNC", "comom") == "exact,deterministic");
    CHECK(method_type("SolverNC", "ca") == "exact,deterministic");
    CHECK(method_type("SolverMVA", "exact") == "exact,deterministic");
    CHECK(method_type("SolverNC", "mci") == "approximate,randomized");
    CHECK(method_type("SolverNC", "imci") == "approximate,randomized");
    CHECK(method_type("SolverNC", "ls") == "approximate,randomized");
    CHECK(method_type("SolverNC", "is") == "approximate,randomized");
    CHECK(method_type("SolverNC", "sampling") == "approximate,randomized");
    // The mixed limited-load-dependent route evaluates the product form itself
    // (Bruell-Balbo-Afshari effective capacity), so it is exact and not an
    // expansion -- it sat in the expansion block until 2026-08-29.
    CHECK(method_type("SolverNC", "ncldmx") == "exact,deterministic");
    CHECK(method_type("SolverNC", "default/ncldmx") == "exact,deterministic");
    // Deterministic NC methods must NOT be swept in with the five above.
    CHECK(method_type("SolverNC", "le") == "approximate,deterministic");
    CHECK(method_type("SolverNC", "propfair") == "approximate,deterministic");

    // The banner form: a resolved default carries its prefix, and the reference
    // classifies on the resolved tail rather than on the word 'default'.
    CHECK(method_type("SolverNC", "default/comom") == "exact,deterministic");
    CHECK(method_type("SolverNC", "default/imci") == "approximate,randomized");
    CHECK(method_type("SolverMAM", "default/dec.source") == "approximate,deterministic");

    // Every SSA method is a sample path, in this port as in
    // matlab/src/solvers/SSA/@SolverSSA/SolverSSA.m:58-62.
    CHECK(method_type("SolverSSA", "nrm") == "approximate,randomized");
    CHECK(method_type("SolverSSA", "serial") == "approximate,randomized");
    CHECK(method_type("SolverSSA", "para") == "approximate,randomized");
    CHECK(method_type("SolverSSA", "parallel") == "approximate,randomized");

    // A bound reports the side it lies on; a family with no sided variant does not.
    CHECK(method_type("SolverBA", "gb.upper") == "upper bound,deterministic");
    CHECK(method_type("SolverBA", "gb.lower") == "lower bound,deterministic");
    CHECK(method_type("SolverBA", "lr") == "bound,deterministic");

    // Reached by the head-of-name lookup: only 'dec' is a token, and neither
    // 'source.mmap' nor 'mmap' is, so the name resolves on its first component.
    CHECK(method_type("SolverMAM", "dec.source.mmap") == "approximate,deterministic");
    CHECK(method_type("SolverMVA", "amva.qd") == "approximate,deterministic");

    CHECK(method_type("SolverCTMC", "gpu") == "exact,deterministic");
    CHECK(method_type("SolverFLD", "matrix") == "approximate,deterministic");
    CHECK(method_type("SolverFLD", "rmf") == "approximate,deterministic");

    // 'default' is the one name that legitimately falls to the per-solver
    // label, because no algorithm has been resolved yet when it is printed.
    CHECK(method_type("SolverCTMC", "default") == "exact,deterministic");
    CHECK(method_type("SolverSSA", "default") == "approximate,randomized");
    CHECK(method_type("SolverBA", "default") == "bound,deterministic");
}

// 2026-07-30: `-o json` was ACCEPTED AND IGNORED by `-s ssa` and `-s fluid`.
// Both arms hand-rolled their own table printer -- their solution types are
// SsaSolution/FluidSolution rather than mva::AvgResult<T> -- and neither
// consulted g_json_output, so a caller that asked for JSON got the readable
// table with no brace in it. The Python lang='cpp' bridge hit exactly that as
// "no JSON object found in solver output". Accepting a flag and not applying it
// is the silent-acceptance defect this CLI refuses everywhere else, so the
// rendering now goes through one emitter (emit_avg_table) for every arm and this
// case pins that every -a avg path answers -o json with an object.
TEST_CASE("-o json is honoured by every solver arm, not only the AvgResult ones") {
    const std::string cqn = write_temp("json_cqn.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})");

    // A value carrying more digits than %12.6g could print is what says the JSON
    // path was taken and not the readable one reformatted. Every arm below has at
    // least one such value on this model EXCEPT `-s ba`, whose bounds here are
    // short decimals, so the check is asked of it only where it is meaningful.
    struct Local {
        static bool has_long_fraction(const std::string& s) {
            for (std::size_t i = 0; i + 1 < s.size(); ++i) {
                if (s[i] != '.') continue;
                std::size_t n = 0;
                while (i + 1 + n < s.size() && std::isdigit(static_cast<unsigned char>(s[i + 1 + n])))
                    ++n;
                if (n >= 9) return true;
            }
            return false;
        }
    };

    const std::string args[] = {"-s mva", "-s nc", "-s ba", "-s ctmc",
                                "-s fluid --method matrix",
                                "-s ssa --samples 2000 --seed 4242"};
    for (const std::string& a : args) {
        const Run r = run_cli(a + " -f " + cqn + " -a avg -o json");
        INFO("path ", a);
        CHECK(r.status == 0);
        // The banner precedes the object on stdout (it carries no brace), so the
        // contract is that an object is THERE and carries the table's keys.
        const std::size_t brace = r.out.find('{');
        REQUIRE(brace != std::string::npos);
        const std::string json = r.out.substr(brace);
        CHECK(json.find("\"avg\"") != std::string::npos);
        CHECK(json.find("\"AvgTable\"") != std::string::npos);
        CHECK(json.find("\"QLen\"") != std::string::npos);
        CHECK(json.find("\"ResidT\"") != std::string::npos);
        CHECK(json.find("\"arith\"") != std::string::npos);
        CHECK(json.find("\"method\"") != std::string::npos);
        if (a != "-s ba") CHECK(Local::has_long_fraction(json));
    }

    std::remove(cqn.c_str());
}

/**
 * The model the non-average arms are exercised on, written once.
 *
 * Delay -> PS Queue with one closed class of 2 jobs: a three-state chain, small
 * enough that a generator, a state space and a stationary law can all be checked
 * by value rather than only by shape.
 */
namespace {

std::string write_nonavg_model(const std::string& name, bool with_reward) {
    std::string rewards;
    if (with_reward)
        rewards = R"(, "rewards": [{"name": "qlen_q", "type": "QLen", "node": "Queue"}])";
    return write_temp(name, std::string(R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}})") + rewards + "}}");
}

}  // namespace

TEST_CASE("-o json answers every non-average analysis under a key named after its -a") {
    const std::string cqn = write_nonavg_model("nonavg_cqn.json", false);

    // The pair is (invocation, the key its answer must arrive under). A payload
    // under any other key would be an answer to a question the caller did not
    // ask, which is exactly what reading the object positionally would hide.
    struct Case {
        const char* args;
        const char* key;
    };
    const Case cases[] = {
        {"-s mva -a prob", "prob"},
        {"-s mva -a marg", "marg"},
        {"-s mva -a normconst", "normconst"},
        {"-s nc -a prob", "prob"},
        {"-s nc -a normconst", "normconst"},
        {"-s ctmc -a gen", "gen"},
        {"-s ctmc -a states", "states"},
        {"-s ctmc -a sens", "sens"},
        {"-s ctmc -a cdf", "cdf"},
        {"-s ctmc -a tranprob --tspan 0:2", "tranprob"},
        {"-s ctmc -a sample --samples 25 --seed 4242", "sample"},
        {"-s fluid --method matrix -a odes", "odes"},
    };
    for (const Case& c : cases) {
        const Run r = run_cli(std::string(c.args) + " -f " + cqn + " -o json");
        INFO("path ", c.args);
        CHECK(r.status == 0);
        const std::size_t brace = r.out.find('{');
        REQUIRE(brace != std::string::npos);
        const Json j = Json::parse(r.out.substr(brace));
        CHECK(j.contains(c.key));
        CHECK(j.contains("arith"));
        // `-a reward` and `-a cdf` return no solved chain, so they carry no
        // resolved method and the key is ABSENT rather than echoing "default"
        // back; every other arm reports the method that ran.
        if (std::string(c.key) != "cdf") CHECK(j.contains("method"));
    }

    std::remove(cqn.c_str());
}

/**
 * `-a marg` and `-a normconst`, the two @@SolverMVA probability getters that
 * had no CLI flag, against MATLAB on the Delay -> PS queue with 2 closed jobs.
 *
 * The figures are `SolverMVA(model).getProbMarg(i,1)` and
 * `getProbNormConstAggr()` run on that model in MATLAB R2025a: the Schmidt
 * binomial fitted to Q = (1.2, 0.8) over N = 2, and log(5/4). CHECKED BY VALUE
 * and not only by shape, because a marginal that sums to one is not thereby the
 * right marginal -- any binomial does.
 */
TEST_CASE("-a marg and -a normconst match MATLAB on the two-job closed model") {
    const std::string cqn = write_nonavg_model("nonavg_marg.json", false);

    const Run m = run_cli("-s mva -a marg -f " + cqn + " -o json");
    REQUIRE(m.status == 0);
    // The parsed object is NAMED, not bound through a reference to a
    // subobject of a temporary: that reference dangles the moment the full
    // expression ends, and reads back as an empty array rather than as an error.
    const Json mj = Json::parse(m.out.substr(m.out.find('{')));
    const Json& curves = mj.at("marg").at("marginal");
    REQUIRE(curves.size() == 2);  // two stations, one class
    const double ref[2][3] = {{0.16, 0.48, 0.36}, {0.36, 0.48, 0.16}};
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("station ", i);
        CHECK(curves[i].at("station").get<std::size_t>() == i);
        CHECK(curves[i].at("jobclass").get<std::size_t>() == 0);
        const std::vector<double> P = curves[i].at("P").get<std::vector<double> >();
        const std::vector<long> jobs = curves[i].at("Jobs").get<std::vector<long> >();
        REQUIRE(P.size() == 3);  // n = 0, 1, 2 over the class population
        for (std::size_t n = 0; n < 3; ++n) {
            CHECK(jobs[n] == static_cast<long>(n));
            CHECK(P[n] == doctest::Approx(ref[i][n]).epsilon(1e-9));
        }
    }

    // --marg-states selects FROM that curve, so the same numbers must come back
    // under the n the caller named and in that order.
    const Run sel = run_cli("-s mva -a marg --node 2 --class 1 --marg-states 2,0 -f " + cqn +
                            " -o json");
    REQUIRE(sel.status == 0);
    const Json sj = Json::parse(sel.out.substr(sel.out.find('{')));
    const Json& one = sj.at("marg").at("marginal");
    REQUIRE(one.size() == 1);
    const std::vector<long> jobs_sel = one[0].at("Jobs").get<std::vector<long> >();
    REQUIRE(jobs_sel.size() == 2);
    CHECK(jobs_sel[0] == 2);
    CHECK(jobs_sel[1] == 0);
    const std::vector<double> Psel = one[0].at("P").get<std::vector<double> >();
    CHECK(Psel[0] == doctest::Approx(0.16).epsilon(1e-9));
    CHECK(Psel[1] == doctest::Approx(0.36).epsilon(1e-9));

    // THE TWO SOLVERS MUST AGREE ON log G, and do: it is the same constant of
    // the same product-form model, however each of them reaches it.
    for (const char* solver : {"-s mva", "-s nc"}) {
        const Run g = run_cli(std::string(solver) + " -a normconst -f " + cqn + " -o json");
        INFO("solver ", solver);
        REQUIRE(g.status == 0);
        const Json j = Json::parse(g.out.substr(g.out.find('{')));
        CHECK(j.at("normconst").at("logNormConstAggr").get<double>() ==
              doctest::Approx(std::log(1.25)).epsilon(1e-9));
    }

    // The state list is the reference's own bound: a state above the class
    // population is an error there and must be one here.
    const Run over = run_cli("-s mva -a marg --marg-states 5 -f " + cqn);
    CHECK(over.status != 0);
    CHECK(over.out.find("exceeds the maximum population") != std::string::npos);

    // --class and --marg-states belong to ONE analysis; anywhere else they are
    // refused rather than dropped, which is what would make them look honoured.
    const Run stray = run_cli("-s mva -a prob --class 1 -f " + cqn);
    CHECK(stray.status != 0);
    CHECK(stray.out.find("getProbMarg") != std::string::npos);

    std::remove(cqn.c_str());
}

TEST_CASE("the -a gen payload rebuilds Q and its filtration exactly") {
    const std::string cqn = write_nonavg_model("nonavg_gen.json", false);
    const Run r = run_cli("-s ctmc -a gen -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    const Json j = Json::parse(r.out.substr(r.out.find('{')));
    const Json& g = j.at("gen");
    CHECK(g.at("indexBase").get<int>() == 0);
    const std::size_t n = g.at("size").get<std::size_t>();
    CHECK(n == 3);
    CHECK(g.at("states").get<std::size_t>() == n);
    const std::vector<std::size_t> from = g.at("From").get<std::vector<std::size_t> >();
    const std::vector<std::size_t> to = g.at("To").get<std::vector<std::size_t> >();
    const std::vector<double> rate = g.at("Rate").get<std::vector<double> >();
    REQUIRE(from.size() == to.size());
    REQUIRE(from.size() == rate.size());
    CHECK(from.size() == g.at("nnz").get<std::size_t>());
    // EVERY ROW OF A GENERATOR SUMS TO ZERO, and rebuilding it from the triplets
    // is the only way to see that the sparse form is the matrix and not a
    // selection of it.
    std::vector<double> rowsum(n, 0.0);
    for (std::size_t k = 0; k < rate.size(); ++k) {
        REQUIRE(from[k] < n);
        REQUIRE(to[k] < n);
        rowsum[from[k]] += rate[k];
    }
    for (std::size_t i = 0; i < n; ++i) CHECK(rowsum[i] == doctest::Approx(0.0).epsilon(1e-12));

    // sum_a filt[a] is the OFF-DIAGONAL part of Q: the filtration is half of what
    // getInfGen returns, and the halves have to add up.
    std::vector<std::vector<double> > acc(n, std::vector<double>(n, 0.0));
    for (const Json& e : g.at("sync")) {
        const std::vector<std::size_t> ff = e.at("From").get<std::vector<std::size_t> >();
        const std::vector<std::size_t> ft = e.at("To").get<std::vector<std::size_t> >();
        const std::vector<double> fr = e.at("Rate").get<std::vector<double> >();
        CHECK(ff.size() == e.at("nnz").get<std::size_t>());
        for (std::size_t k = 0; k < fr.size(); ++k) acc[ff[k]][ft[k]] += fr[k];
    }
    for (std::size_t k = 0; k < rate.size(); ++k)
        if (from[k] != to[k])
            CHECK(acc[from[k]][to[k]] == doctest::Approx(rate[k]).epsilon(1e-12));

    std::remove(cqn.c_str());
}

TEST_CASE("the -a states payload pairs the space with the law that indexes it") {
    const std::string cqn = write_nonavg_model("nonavg_states.json", false);
    const Run r = run_cli("-s ctmc -a states -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    // The parsed object is bound to a VALUE first: binding a reference straight
    // to `.at(...)` of the temporary leaves it dangling, which doctest reports as
    // "cannot use at() with null" rather than as the lifetime bug it is.
    const Json parsed = Json::parse(r.out.substr(r.out.find('{')));
    const Json& s = parsed.at("states");
    const std::vector<double> pi = s.at("pi").get<std::vector<double> >();
    CHECK(s.at("space").size() == pi.size());
    CHECK(s.at("spaceAggr").size() == pi.size());
    // The payload key is "states" and the state COUNT is also "states"; nesting
    // the meta inside the payload is what keeps the second from replacing the
    // first, which it did when both sat at envelope level.
    CHECK(s.at("states").get<std::size_t>() == pi.size());
    double total = 0.0;
    for (double p : pi) total += p;
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));
    std::size_t widths = 0;
    for (const Json& w : s.at("NodeWidths")) widths += w.get<std::size_t>();
    CHECK(s.at("space").at(0).size() == widths);

    std::remove(cqn.c_str());
}

TEST_CASE("the -a sample payload carries the space its indices point into") {
    const std::string cqn = write_nonavg_model("nonavg_sample.json", false);
    const Run r = run_cli("-s ctmc -a sample --samples 25 --seed 4242 -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    const Json parsed = Json::parse(r.out.substr(r.out.find('{')));
    const Json& p = parsed.at("sample");
    CHECK(p.at("seed").get<unsigned long>() == 4242UL);
    CHECK(p.at("events").get<std::size_t>() == 25);
    const std::vector<std::size_t> st = p.at("state").get<std::vector<std::size_t> >();
    const Json& space = p.at("space");
    REQUIRE(!st.empty());
    // A trajectory of indices into a chain the caller cannot see is not a
    // trajectory: the space travels with it, and the population is conserved on
    // every row of it because the model is closed.
    for (std::size_t i = 0; i < st.size(); ++i) {
        REQUIRE(st[i] < space.size());
        double jobs = 0.0;
        for (const Json& v : space.at(st[i])) jobs += v.get<double>();
        CHECK(jobs == doctest::Approx(2.0).epsilon(1e-12));
    }
    // The same seed twice is the same trace; that is what makes it a measurement.
    const Run again = run_cli("-s ctmc -a sample --samples 25 --seed 4242 -f " + cqn + " -o json");
    CHECK(again.out == r.out);

    std::remove(cqn.c_str());
}

TEST_CASE("a declared reward crosses the wire and is evaluated on the chain") {
    const std::string cqn = write_nonavg_model("nonavg_reward.json", true);
    const Run r = run_cli("-s ctmc -a reward -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    const Json j = Json::parse(r.out.substr(r.out.find('{')));
    const Json& rw = j.at("reward");
    CHECK(rw.at("Reward").at(0).get<std::string>() == "qlen_q");
    // The QLen template at the PS queue is the queue's mean queue length, so the
    // reward's expectation IS the AvgTable's QLen for that station -- the reward
    // reader is checked against a number the solver already reports, not against a
    // recomputation of it here.
    const Run avg = run_cli("-s ctmc -a avg -f " + cqn + " -o json");
    REQUIRE(avg.status == 0);
    const Json a = Json::parse(avg.out.substr(avg.out.find('{')));
    const Json& tbl = a.at("avg");
    double qlen_queue = 0.0;
    for (std::size_t i = 0; i < tbl.at("Station").size(); ++i)
        if (tbl.at("Station").at(i).get<std::string>() == "Queue")
            qlen_queue += tbl.at("QLen").at(i).get<double>();
    CHECK(rw.at("E").at(0).get<double>() == doctest::Approx(qlen_queue).epsilon(1e-9));
    // No resolved method on this arm: it returns expectations and no solved chain.
    CHECK(!j.contains("method"));

    std::remove(cqn.c_str());
}

TEST_CASE("a reward the reader cannot reproduce is refused by name") {
    const std::string bad = write_temp("nonavg_badreward.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}},
        "rewards": [{"name": "mystery", "type": "Custom", "node": "Queue"}]}})");
    const Run r = run_cli("-s ctmc -a reward -f " + bad + " -o json");
    CHECK(r.status != 0);
    CHECK(r.out.find("Custom") != std::string::npos);
    CHECK(r.out.find("QLen, Util and Blocking") != std::string::npos);
    std::remove(bad.c_str());
}

TEST_CASE("--notation reaches the ODE export and is refused everywhere else") {
    const std::string cqn = write_nonavg_model("nonavg_odes.json", false);
    const Run scalar =
        run_cli("-s fluid --method matrix -a odes --notation scalar -f " + cqn + " -o json");
    const Run matrix =
        run_cli("-s fluid --method matrix -a odes --notation matrix -f " + cqn + " -o json");
    REQUIRE(scalar.status == 0);
    REQUIRE(matrix.status == 0);
    const Json js = Json::parse(scalar.out.substr(scalar.out.find('{')));
    const Json jm = Json::parse(matrix.out.substr(matrix.out.find('{')));
    CHECK(js.at("odes").at("notation").get<std::string>() == "scalar");
    CHECK(jm.at("odes").at("notation").get<std::string>() == "matrix");
    // Two documents of one drift: accepting the flag and writing the scalar form
    // for both would be the silent acceptance every other knob here refuses.
    CHECK(js.at("odes").at("latex").get<std::string>() !=
          jm.at("odes").at("latex").get<std::string>());

    const Run wrong = run_cli("-s mva -a avg --notation matrix -f " + cqn);
    CHECK(wrong.status != 0);
    CHECK(wrong.out.find("--notation") != std::string::npos);
    const Run unknown =
        run_cli("-s fluid --method matrix -a odes --notation runes -f " + cqn);
    CHECK(unknown.status != 0);
    CHECK(unknown.out.find("runes") != std::string::npos);

    std::remove(cqn.c_str());
}

TEST_CASE("JobClass::completes defaults to true, as in the other three codebases") {
    // The response-time law needs a completing class in the tagged chain to
    // absorb at; with the default at false, `-a cdf` refused on every ordinary
    // model read from JSON, since nothing on the wire carries the flag.
    const std::string cqn = write_nonavg_model("nonavg_completes.json", false);
    const Run r = run_cli("-s ctmc -a cdf -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    const Json parsed = Json::parse(r.out.substr(r.out.find('{')));
    const Json& c = parsed.at("cdf");
    REQUIRE(!c.at("respt").empty());
    const Json& curve = c.at("respt").at(0);
    const std::vector<double> F = curve.at("F").get<std::vector<double> >();
    REQUIRE(F.size() > 1);
    // A CDF is non-decreasing and reaches 1; a law built from an all-zero
    // filtration -- what an empty completing set would have produced -- does not.
    for (std::size_t i = 1; i < F.size(); ++i) CHECK(F[i] >= F[i - 1] - 1e-12);
    CHECK(F.back() == doctest::Approx(1.0).epsilon(1e-3));

    std::remove(cqn.c_str());
}

TEST_CASE("-s uq expands a Prior and reports the design it averaged over") {
    // M/M/1 at lambda = 1 whose service rate is 4 with weight 1/4 and 2 with
    // weight 3/4: the alternatives are rho = 1/4 and rho = 1/2, whose queue
    // lengths are 1/3 and 1, so the expectation is 0.25/3 + 0.75 = 0.8333.
    const std::string mm1 = write_temp("uq_mm1.json", R"({
      "format":"line-model","version":"1.0","model":{
       "type":"Network","name":"uqmm1",
       "nodes":[
        {"name":"Source 1","type":"Source","service":{
           "Class1":{"type":"Exp","params":{"lambda":1.0}}}},
        {"name":"Queue 1","type":"Queue","scheduling":"FCFS","service":{
           "Class1":{"type":"Prior",
                     "distributions":[{"type":"Exp","params":{"lambda":4.0}},
                                      {"type":"Exp","params":{"lambda":2.0}}],
                     "probabilities":[0.25,0.75]}}},
        {"name":"Sink 1","type":"Sink"}],
       "classes":[{"name":"Class1","type":"Open"}],
       "routing":{"type":"matrix","matrix":{
         "Class1,Class1":{"Source 1":{"Queue 1":1.0},"Queue 1":{"Sink 1":1.0}}}}
      }})");

    const Run avg = run_cli("-s uq --uq-solver mva -o json -f " + mm1);
    REQUIRE(avg.status == 0);
    const Json a = Json::parse(avg.out.substr(avg.out.find('{'))).at("avg");
    const std::vector<double> Q = a.at("QLen").get<std::vector<double> >();
    REQUIRE(Q.size() == 2);
    CHECK(Q[1] == doctest::Approx(0.25 / 3.0 + 0.75).epsilon(1e-9));

    SUBCASE("-a posterior carries the points, the weights and what was substituted") {
        const Run post = run_cli("-s uq --uq-solver mva -a posterior -o json -f " + mm1);
        REQUIRE(post.status == 0);
        const Json p = Json::parse(post.out.substr(post.out.find('{'))).at("posterior");
        CHECK(p.at("design") == "quadrature");
        CHECK(p.at("stage") == "mva");
        CHECK(p.at("priors").size() == 1);
        CHECK(p.at("priors").at(0).at("kind") == "service");
        // The substituted MEANS, 1/4 and 1/2, are the provenance of the two rows.
        const std::vector<std::vector<double> > sub =
            p.at("substitutedMean").get<std::vector<std::vector<double> > >();
        REQUIRE(sub.size() == 2);
        CHECK(sub[0][0] == doctest::Approx(0.25));
        CHECK(sub[1][0] == doctest::Approx(0.5));
        const std::vector<double> w = p.at("Weight").get<std::vector<double> >();
        const std::vector<double> pq = p.at("QLen").get<std::vector<double> >();
        double acc = 0.0;
        for (std::size_t i = 0; i < w.size(); ++i) acc += w[i] * pq[i];
        CHECK(acc == doctest::Approx(0.25 / 3.0 + 0.75).epsilon(1e-9));
    }
    SUBCASE("the stage solver is required, and named when it is not ported") {
        const Run missing = run_cli("-s uq -f " + mm1);
        CHECK(missing.status != 0);
        CHECK(missing.out.find("--uq-solver") != std::string::npos);
        const Run unported = run_cli("-s uq --uq-solver ldes -f " + mm1);
        CHECK(unported.status != 0);
        CHECK(unported.out.find("ldes") != std::string::npos);
    }
    SUBCASE("--uq-solver is refused by the solvers that solve one model") {
        const Run wrong = run_cli("-s mva --uq-solver nc -f " + mm1);
        CHECK(wrong.status != 0);
        CHECK(wrong.out.find("--uq-solver") != std::string::npos);
    }
    SUBCASE("every other solver refuses the Prior rather than averaging it away") {
        const Run plain = run_cli("-s mva -f " + mm1);
        CHECK(plain.status != 0);
        CHECK(plain.out.find("Prior") != std::string::npos);
    }
    SUBCASE("-a interval on an open model falls back and says it is not an enclosure") {
        const Run iv = run_cli("-s uq --uq-solver mva -a interval -o json -f " + mm1);
        REQUIRE(iv.status == 0);
        const Json p = Json::parse(iv.out.substr(iv.out.find('{'))).at("interval");
        CHECK(p.at("exact") == false);
        CHECK(p.at("intervalMethod") == "sampled");
        // The Prior sits on a SERVICE process here, so the first condition to
        // fire is the closedness of the class, not the arrival test.
        CHECK(std::string(p.at("why")).find("not closed") != std::string::npos);
        // The sampled range spans the two alternatives: rho = 1/4 and 1/2.
        const std::vector<double> lo = p.at("QLen_lo").get<std::vector<double> >();
        const std::vector<double> up = p.at("QLen_up").get<std::vector<double> >();
        CHECK(lo.back() == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
        CHECK(up.back() == doctest::Approx(1.0).epsilon(1e-9));
    }
    SUBCASE("the Monte Carlo design draws options.samples points of equal weight") {
        const Run mc = run_cli("-s uq --uq-solver mva --method montecarlo --samples 8 --seed 7 "
                               "-a posterior -o json -f " + mm1);
        REQUIRE(mc.status == 0);
        const Json p = Json::parse(mc.out.substr(mc.out.find('{'))).at("posterior");
        CHECK(p.at("design") == "montecarlo");
        CHECK(p.at("substitutedMean").size() == 8);
        const std::vector<double> w = p.at("Weight").get<std::vector<double> >();
        for (double x : w) CHECK(x == doctest::Approx(0.125));
    }

    std::remove(mm1.c_str());
}

TEST_CASE("-s uq -a interval takes the exact hull on a qualifying closed model") {
    // Single class, closed, one delay and one PS queue whose service is
    // uncertain between mean 0.25 and mean 0.5: the monotonicity theorems apply,
    // so the interval is ATTAINED at the two demand endpoints and no design
    // point is solved at all.
    const std::string cqn = write_temp("uq_closed.json", R"({
      "format":"line-model","version":"1.0","model":{
       "type":"Network","name":"uqclosed",
       "nodes":[
        {"name":"Think","type":"Delay","service":{
           "C1":{"type":"Exp","params":{"lambda":1.0}}}},
        {"name":"Queue 1","type":"Queue","scheduling":"PS","service":{
           "C1":{"type":"Prior",
                 "distributions":[{"type":"Exp","params":{"lambda":2.0}},
                                  {"type":"Exp","params":{"lambda":4.0}}],
                 "probabilities":[0.5,0.5]}}}],
       "classes":[{"name":"C1","type":"Closed","population":2,"refNode":"Think"}],
       "routing":{"type":"matrix","matrix":{
         "C1,C1":{"Think":{"Queue 1":1.0},"Queue 1":{"Think":1.0}}}}
      }})");

    const Run iv = run_cli("-s uq --uq-solver mva -a interval -o json -f " + cqn);
    REQUIRE(iv.status == 0);
    const Json p = Json::parse(iv.out.substr(iv.out.find('{'))).at("interval");
    CHECK(p.at("exact") == true);
    CHECK(p.at("intervalMethod") == "mvainterval");
    CHECK_FALSE(p.contains("why"));
    // The queue's own MVA answers at D = 0.25 and D = 0.5 with Z = 1, N = 2.
    const std::vector<double> lo = p.at("QLen_lo").get<std::vector<double> >();
    const std::vector<double> up = p.at("QLen_up").get<std::vector<double> >();
    CHECK(lo.back() == doctest::Approx(6.0 / 13.0).epsilon(1e-9));
    CHECK(up.back() == doctest::Approx(0.8).epsilon(1e-9));
    const std::vector<double> X = p.at("X").get<std::vector<double> >();
    CHECK(X[0] == doctest::Approx(1.2).epsilon(1e-9));
    CHECK(X[1] == doctest::Approx(20.0 / 13.0).epsilon(1e-9));

    SUBCASE("and the expectation lies inside the hull it brackets") {
        const Run avg = run_cli("-s uq --uq-solver mva -o json -f " + cqn);
        REQUIRE(avg.status == 0);
        const Json a = Json::parse(avg.out.substr(avg.out.find('{'))).at("avg");
        const std::vector<double> Q = a.at("QLen").get<std::vector<double> >();
        CHECK(Q.back() >= lo.back() - 1e-9);
        CHECK(Q.back() <= up.back() + 1e-9);
    }

    std::remove(cqn.c_str());
}

// ---------------------------------------------------------------------------
// -a node: getAvgNodeTable, the analysis in the NODE index space
// ---------------------------------------------------------------------------

TEST_CASE("-a node reports the nodes the station table cannot have a row for") {
    // Delay -> Router -> Queue -> Delay. The Router is the point of the fixture:
    // it is not a station, so `-a avg` cannot report it at all, yet every job in
    // the closed chain passes through it once per cycle.
    const std::string cqn = write_temp("nodetable_router.json", R"({
      "format": "line-model", "version": "1.0",
      "model": {"type": "Network", "name": "cqn_router",
        "nodes": [
          {"name": "Delay", "type": "Delay",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 1.0}}}},
          {"name": "Fan", "type": "Router"},
          {"name": "Queue", "type": "Queue", "scheduling": "PS",
           "service": {"Class1": {"type": "Exp", "fit": {"method": "fitMean", "mean": 0.5}}}}],
        "classes": [{"name": "Class1", "type": "Closed", "population": 2, "refNode": "Delay"}],
        "routing": {"type": "matrix", "matrix": {"Class1,Class1": {
          "Delay": {"Fan": 1.0}, "Fan": {"Queue": 1.0}, "Queue": {"Delay": 1.0}}}}}})");

    // MATLAB-verified exact MVA on Z = 1, D = 0.5, N = 2: R = 0.5(1 + 1/3),
    // X = 2/(1 + 2/3) = 1.2, Q_queue = 0.8 and so Q_delay = 1.2.
    const double X = 1.2, QQ = 0.8, QD = 1.2, UQ = 0.6;

    const Run r = run_cli("-s mva -a node -f " + cqn + " -o json");
    REQUIRE(r.status == 0);
    const Json j = Json::parse(r.out.substr(r.out.find('{')));
    const Json& p = j.at("node");
    CHECK(p.at("type").get<std::string>() == "AvgNodeTable");
    const std::vector<std::string> node = p.at("Node").get<std::vector<std::string> >();
    const std::vector<double> Q = p.at("QLen").get<std::vector<double> >();
    const std::vector<double> U = p.at("Util").get<std::vector<double> >();
    const std::vector<double> A = p.at("ArvR").get<std::vector<double> >();
    const std::vector<double> T = p.at("Tput").get<std::vector<double> >();

    std::map<std::string, std::size_t> row;
    for (std::size_t i = 0; i < node.size(); ++i) row[node[i]] = i;
    REQUIRE(row.count("Fan") == 1);  // the whole reason the analysis exists
    REQUIRE(row.count("Delay") == 1);
    REQUIRE(row.count("Queue") == 1);

    // The station rows are the AvgTable's, scattered to their node indices.
    CHECK(Q[row["Queue"]] == doctest::Approx(QQ).epsilon(1e-9));
    CHECK(U[row["Queue"]] == doctest::Approx(UQ).epsilon(1e-9));
    CHECK(Q[row["Delay"]] == doctest::Approx(QD).epsilon(1e-9));
    CHECK(T[row["Queue"]] == doctest::Approx(X).epsilon(1e-9));

    // A ROUTER HOLDS NO JOBS AND SERVES NONE, so its queue length and
    // utilization are zero while its flow is the whole chain's: that pair is
    // the difference between this table and the AvgTable, not a rounding of it.
    CHECK(Q[row["Fan"]] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(U[row["Fan"]] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(T[row["Fan"]] == doctest::Approx(X).epsilon(1e-9));
    CHECK(A[row["Fan"]] == doctest::Approx(X).epsilon(1e-9));

    SUBCASE("and the station rows agree with -a avg column for column") {
        const Run avg = run_cli("-s mva -f " + cqn + " -o json");
        REQUIRE(avg.status == 0);
        const Json a = Json::parse(avg.out.substr(avg.out.find('{'))).at("avg");
        const std::vector<std::string> st = a.at("Station").get<std::vector<std::string> >();
        const std::vector<double> aq = a.at("QLen").get<std::vector<double> >();
        const std::vector<double> at = a.at("Tput").get<std::vector<double> >();
        for (std::size_t i = 0; i < st.size(); ++i) {
            INFO("station ", st[i]);
            REQUIRE(row.count(st[i]) == 1);
            CHECK(Q[row[st[i]]] == doctest::Approx(aq[i]).epsilon(1e-12));
            CHECK(T[row[st[i]]] == doctest::Approx(at[i]).epsilon(1e-12));
        }
    }

    SUBCASE("every solver that returns an AvgResult answers it, and agrees") {
        for (const char* s : {"-s nc", "-s ctmc", "-s ba"}) {
            const Run o = run_cli(std::string(s) + " -a node -f " + cqn + " -o json");
            INFO("solver ", s);
            REQUIRE(o.status == 0);
            const Json oj = Json::parse(o.out.substr(o.out.find('{'))).at("node");
            CHECK(oj.at("type").get<std::string>() == "AvgNodeTable");
            const std::vector<std::string> n2 = oj.at("Node").get<std::vector<std::string> >();
            CHECK(std::find(n2.begin(), n2.end(), std::string("Fan")) != n2.end());
        }
        // The two exact engines must give the Router the same flow; BA is a
        // bound and is deliberately not held to it.
        const Run nc = run_cli("-s nc -a node -f " + cqn + " -o json");
        const Json nj = Json::parse(nc.out.substr(nc.out.find('{'))).at("node");
        const std::vector<std::string> n2 = nj.at("Node").get<std::vector<std::string> >();
        const std::vector<double> t2 = nj.at("Tput").get<std::vector<double> >();
        for (std::size_t i = 0; i < n2.size(); ++i)
            if (n2[i] == "Fan") CHECK(t2[i] == doctest::Approx(X).epsilon(1e-9));
    }

    SUBCASE("a solver whose runner returns no AvgResult refuses it by name") {
        // 'qns' reaches the -a node gate and is turned away by it. `ssa` and
        // `fluid` were HERE until their runners started returning the station
        // AvgResult; the arm's own message lists them now, and the subcase
        // below asserts they answer rather than refuse. `ldes`, `jmt`, `uq` and
        // `env` are refused EARLIER, by their own analysis ladders, so they do
        // not exercise this gate.
        const Run o = run_cli("-s qns -a node -f " + cqn);
        CHECK(o.status != 0);
        CHECK(o.out.find("-a node reports the per-node table") != std::string::npos);
    }

    SUBCASE("the sample-path and fluid engines answer it too") {
        // Both runners return the station AvgResult the node table is built
        // from, so `-a node` is served rather than refused. The Router is the
        // row `-a avg` cannot have, which is the whole point of the analysis.
        for (const char* s : {"-s ssa", "-s fluid"}) {
            const Run o = run_cli(std::string(s) + " -a node -f " + cqn + " -o json");
            INFO("solver ", s);
            REQUIRE(o.status == 0);
            const Json oj = Json::parse(o.out.substr(o.out.find('{'))).at("node");
            CHECK(oj.at("type").get<std::string>() == "AvgNodeTable");
            const std::vector<std::string> n2 = oj.at("Node").get<std::vector<std::string> >();
            CHECK(std::find(n2.begin(), n2.end(), std::string("Fan")) != n2.end());
        }
    }

    SUBCASE("the exact arithmetic backend serves it too") {
        // -s mva and not -s nc: the normalizing-constant analyzer forms
        // X = exp(lG(N-1) - lG(N)) and refuses exact for its own reason, which
        // is a property of that solver and not of this analysis.
        const Run e = run_cli("-s mva -a node --arith exact -f " + cqn + " -o json");
        REQUIRE(e.status == 0);
        CHECK(Json::parse(e.out.substr(e.out.find('{'))).at("arith").get<std::string>() ==
              "exact");
    }

    std::remove(cqn.c_str());
}


TEST_CASE("a model-argument api function is registered and refused for the right reason") {
    // lqn_boxbounds and the sn predicates are ported and belong in --list-api,
    // which is the port's coverage manifest; none of them can cross --api,
    // because each takes a struct rather than matrices. The refusal must say
    // THAT, and not that the caller should wait for a later version.
    CHECK(line::find_api("lqn_boxbounds") != nullptr);
    CHECK(line::find_api("sn_has_product_form") != nullptr);
    CHECK_FALSE(line::reg::api_is_exposed("lqn_boxbounds"));
    for (const char* n : {"lqn_boxbounds", "sn_has_product_form"}) {
        INFO("function ", n);
        try {
            api_invoke(n, "double", Json::object());
            FAIL("expected a refusal");
        } catch (const line::UnsupportedError& e) {
            const std::string msg = e.what();
            CHECK(msg.find("not matrices") != std::string::npos);
            CHECK(msg.find("-s <solver> -a <analysis>") != std::string::npos);
        }
    }
    // The registry's own domain strings stay discoverable: lsn_max_multiplicity
    // is filed under "lqn" and not under an "lsn" domain of one member.
    CHECK(line::find_api("lsn_max_multiplicity")->domain == "lqn");
    CHECK(line::find_api("lqn_boxbounds")->domain == "lqn");
}

TEST_CASE("-i pnml reaches the PNML reader and the solver that can answer") {
    // The document is the one all four writers emit for the three-token,
    // two-place cycle, byte for byte. Reading THAT rather than a hand-written
    // net is what makes this a test of the interchange and not of a fixture:
    // the same bytes are what MATLAB, the JAR and python were run on, and all
    // four returned QLen 1.148571... at P1 and 1.851428... at P2. The digits
    // asserted below are this CLI's own six significant figures; the other
    // three print five, so the shared prefix is what the codebases agree on.
    const std::string doc =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">\n"
        "  <net id=\"cyclicspn\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">\n"
        "    <name><text>cyclicspn</text></name>\n"
        "    <page id=\"page0\">\n"
        "      <place id=\"P1\">\n"
        "        <name><text>P1</text></name>\n"
        "        <initialMarking><text>3</text></initialMarking>\n"
        "      </place>\n"
        "      <place id=\"P2\">\n"
        "        <name><text>P2</text></name>\n"
        "        <initialMarking><text>0</text></initialMarking>\n"
        "      </place>\n"
        "      <transition id=\"T1\">\n"
        "        <name><text>T1</text></name>\n"
        "        <toolspecific tool=\"LINE\" version=\"3.0\">\n"
        "          <mode transition=\"T1\" name=\"Mode1\" timing=\"timed\" servers=\"1\""
        " priority=\"1\" weight=\"1\">\n"
        "            <distribution name=\"Exp\">\n"
        "              <parameter name=\"lambda\" value=\"2\"/>\n"
        "            </distribution>\n"
        "          </mode>\n"
        "        </toolspecific>\n"
        "      </transition>\n"
        "      <transition id=\"T2\">\n"
        "        <name><text>T2</text></name>\n"
        "        <toolspecific tool=\"LINE\" version=\"3.0\">\n"
        "          <mode transition=\"T2\" name=\"Mode2\" timing=\"timed\" servers=\"1\""
        " priority=\"1\" weight=\"1\">\n"
        "            <distribution name=\"Exp\">\n"
        "              <parameter name=\"lambda\" value=\"1.5\"/>\n"
        "            </distribution>\n"
        "          </mode>\n"
        "        </toolspecific>\n"
        "      </transition>\n"
        "      <arc id=\"a1\" source=\"P1\" target=\"T1\">\n"
        "        <inscription><text>1</text></inscription>\n"
        "      </arc>\n"
        "      <arc id=\"a2\" source=\"T1\" target=\"P2\">\n"
        "        <inscription><text>1</text></inscription>\n"
        "      </arc>\n"
        "      <arc id=\"a3\" source=\"T2\" target=\"P1\">\n"
        "        <inscription><text>1</text></inscription>\n"
        "      </arc>\n"
        "      <arc id=\"a4\" source=\"P2\" target=\"T2\">\n"
        "        <inscription><text>1</text></inscription>\n"
        "      </arc>\n"
        "    </page>\n"
        "  </net>\n"
        "</pnml>\n";
    const std::string path = write_temp("cyclicspn.pnml", doc);

    SUBCASE("an explicit -i pnml solves it") {
        const Run r = run_cli("-f " + path + " -i pnml -s ctmc -a avg");
        CHECK(r.status == 0);
        CHECK(r.out.find("P1") != std::string::npos);
        CHECK(r.out.find("P2") != std::string::npos);
        CHECK(r.out.find("1.14857") != std::string::npos);
        CHECK(r.out.find("1.85143") != std::string::npos);
    }

    SUBCASE("a bare .pnml path is enough, as a .lqnx path is") {
        const Run r = run_cli("-f " + path + " -s ctmc -a avg");
        CHECK(r.status == 0);
        CHECK(r.out.find("1.14857") != std::string::npos);
    }

    SUBCASE("an unknown input token is refused and names the accepted ones") {
        const Run r = run_cli("-f " + path + " -i pnmlx -s ctmc -a avg");
        CHECK(r.status != 0);
        CHECK(r.out.find("pnml") != std::string::npos);
    }

    std::remove(path.c_str());
}

/**
 * qsys_mapphc over --api. The MAP/PH/c solve reduces to M/M/c when the arrival
 * MAP is Poisson and the PH service has one phase, so the closed form pins the
 * dispatch: lambda = mu = 1, c = 2 gives Erlang-C P(Wq > 0) = 1/3, Lq = 1/3 and
 * L = 4/3, and the per-server utilization is 1/2.
 */
TEST_CASE("qsys_mapphc over --api reduces to M/M/2 and returns its waiting-time block") {
    const Json args = Json::parse(
        R"({"D0": [[-1.0]], "D1": [[1.0]], "alpha": [1.0], "S": [[-1.0]], "c": 2})");
    const Json r = api_invoke("qsys_mapphc", "double", args);
    CHECK(r["function"] == "qsys_mapphc");
    const Json& res = r["results"];
    CHECK(res["meanQueueLength"].get<double>() == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    CHECK(res["meanWaitingTime"].get<double>() == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(res["meanSojournTime"].get<double>() == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    CHECK(res["utilization"].get<double>() == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(res["probWait"].get<double>() == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(res["phaseCount"].get<std::size_t>() == 1u);
    REQUIRE(res["waitingTimeMoments"].size() == 3);
    CHECK(res["waitingTimeMoments"][0].get<double>() == doctest::Approx(1.0 / 3.0).epsilon(1e-9));

    SUBCASE("the three tuning arguments are named together or not at all") {
        Json partial = args;
        partial["num_w_moms"] = 2;
        CHECK_THROWS_AS(api_invoke("qsys_mapphc", "double", partial), line::InputError);
    }

    SUBCASE("all three named selects the tuned overload") {
        Json tuned = args;
        tuned["dist_size"] = 20;
        tuned["num_w_moms"] = 2;
        tuned["w_points"] = Json::parse("[0.5, 1.0]");
        const Json t = api_invoke("qsys_mapphc", "double", tuned);
        CHECK(t["results"]["queueLengthDist"].size() <= 20);
        REQUIRE(t["results"]["waitingTimeMoments"].size() == 2);
        REQUIRE(t["results"]["waitingTimeCCDF"].size() == 2);
        // M/M/2 waiting time: P(Wq > t) = P(Wq > 0) exp(-(c mu - lambda) t).
        CHECK(t["results"]["waitingTimeCCDF"][0].get<double>() ==
              doctest::Approx(std::exp(-0.5) / 3.0).epsilon(1e-9));
        CHECK(t["results"]["waitingTimePoints"][1].get<double>() ==
              doctest::Approx(1.0).epsilon(1e-12));
    }
}
