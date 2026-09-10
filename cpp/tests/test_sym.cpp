/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `api/sym`: the `SageRestEngine` client, the `sym_resolve` ladder, and the
 * three solver getters that delegate to the backend.
 *
 * TWO HALVES, AND ONLY ONE OF THEM CAN ALWAYS RUN. The encoders, the URL
 * normalization, the decimal printer and every refusal in the ladder are this
 * port's own code and are checked unconditionally, WITHOUT TOUCHING THE
 * NETWORK: `none`, an `https://` URL and a port nothing listens on are all
 * decided locally. The rest needs a live line-sage-rest, which is a 3 GB image
 * LINE does not ship, so those cases are SKIPPED BY NAME -- never silently
 * passed, since a green suite that exercised no computer algebra is exactly the
 * failure mode that matters here.
 *
 * HOW TO RUN THE LIVE HALF. Point `LINE_SAGE_URL` (or `LINE_SYM_TEST_URL`, which
 * wins and exists so a suite run can name a service without changing what the
 * SOLVERS resolve) at a service, or leave one listening on 8085 or 8080:
 *
 *   docker run -d --rm -p 8085:8080 imperialqore/line-sage-rest:latest
 *   LINE_SAGE_URL=http://localhost:8085 ./line_mp_tests -ts=api/sym
 *
 * THE GATE DELIBERATELY STOPS AT STEP 3 OF THE LADDER. A bare `sym_resolve`
 * would go on to START A CONTAINER from a locally present image, and a suite
 * that pulls a multi-gigabyte image into the middle of an unrelated run is not
 * one anybody can predict the cost of. That step has its own case, behind
 * `LINE_SYM_TEST_DOCKER=1`.
 *
 * THE ORACLES ARE ANALYTIC, never a normal form read back out of the engine.
 * Printed expressions are not comparable across engine versions -- the header of
 * `sym_engine.h` says so -- so every check here substitutes numbers and compares
 * those: pi of a two-state chain against x2/(x1+x2), the symbolic stationary law
 * against the one `solver_ctmc_analyzer` computes with no backend at all, and
 * the symbolic sensitivity against the finite-difference branch.
 */

#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sym/sym_engines.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_sens.h"
#include "line/solvers/ctmc/solver_ctmc_symbolic.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Sets or clears an environment variable, restoring the previous state after. */
class EnvGuard {
public:
    EnvGuard(const std::string& name, const std::string& value) : name_(name) { save(); set(value); }
    explicit EnvGuard(const std::string& name) : name_(name) {
        save();
        ::unsetenv(name_.c_str());
    }
    ~EnvGuard() {
        if (had_)
            ::setenv(name_.c_str(), prev_.c_str(), 1);
        else
            ::unsetenv(name_.c_str());
    }

private:
    void save() {
        const char* prev = std::getenv(name_.c_str());
        had_ = prev != nullptr;
        if (had_) prev_ = prev;
    }
    void set(const std::string& value) { ::setenv(name_.c_str(), value.c_str(), 1); }

    std::string name_, prev_;
    bool had_ = false;
};

/** Name of the variable a suite run uses to point these cases at a service. */
const char* const TEST_URL_ENV = "LINE_SYM_TEST_URL";

/**
 * A line-sage-rest this suite may use, or a null pointer.
 *
 * Steps 1 to 3 of the ladder and no further: an explicit URL from either
 * variable, else a service already listening on a conventional port, verified
 * through `/api/v1/info` exactly as `sym_resolve` verifies it. Resolved once,
 * since an unreachable probe costs a connect per case otherwise.
 */
std::shared_ptr<sym::SymEngine> live_engine() {
    static bool resolved = false;
    static std::shared_ptr<sym::SymEngine> engine;
    if (resolved) return engine;
    resolved = true;

    std::string url;
    const char* named = std::getenv(TEST_URL_ENV);
    if (named == nullptr || util::trim(named).empty()) named = std::getenv(sym::SYM_URL_ENV);
    if (named != nullptr) url = util::trim(named);

    if (!url.empty()) {
        if (url.compare(0, 7, "http://") != 0) return engine;  // no TLS in this port's client
        std::shared_ptr<sym::SageRestEngine> e = std::make_shared<sym::SageRestEngine>(url);
        if (sym::detail::is_sage_service(*e)) engine = e;
        return engine;
    }

    const std::vector<int> ports = sym::sym_probe_ports();
    for (std::size_t i = 0; i < ports.size(); ++i) {
        std::shared_ptr<sym::SageRestEngine> e = std::make_shared<sym::SageRestEngine>(
            "http://localhost:" + std::to_string(ports[i]));
        if (sym::detail::is_sage_service(*e)) {
            engine = e;
            return engine;
        }
    }
    return engine;
}

/** The base URL of the live engine, for the cases that re-resolve it themselves. */
std::string live_url() {
    const std::shared_ptr<sym::SymEngine> e = live_engine();
    const sym::SageRestEngine* rest = dynamic_cast<const sym::SageRestEngine*>(e.get());
    return rest != nullptr ? rest->getBaseUrl() : std::string();
}

/** True when the live half must be skipped, having said so. */
bool no_service() {
    if (live_engine()) return false;
    // std::string, not the pointer: doctest streams a bare `const char*` as an
    // address, which would turn the one line telling a reader how to enable
    // these cases into two hex numbers.
    MESSAGE("no line-sage-rest reachable: skipping the cases that need one (set "
            << std::string(TEST_URL_ENV) << " or " << std::string(sym::SYM_URL_ENV)
            << ", or run one on port 8085)");
    return true;
}

/** The generator of a two-state chain, x1 up and x2 down. */
std::vector<std::vector<std::string>> two_state_Q() {
    std::vector<std::vector<std::string>> Q(2, std::vector<std::string>(2));
    Q[0][0] = "-x1";
    Q[0][1] = "x1";
    Q[1][0] = "x2";
    Q[1][1] = "-x2";
    return Q;
}

/** The symbols of `two_state_Q`. */
std::vector<std::string> two_state_symbols() {
    std::vector<std::string> s;
    s.push_back("x1");
    s.push_back("x2");
    return s;
}

/** x1 = 2, x2 = 3, the point every live case substitutes. */
std::map<std::string, double> two_state_point() {
    std::map<std::string, double> at;
    at["x1"] = 2.0;
    at["x2"] = 3.0;
    return at;
}

/** Source -> FCFS Queue (capacity K) -> Sink, one open class. */
qn::Network<double> mm1k(double lambda, double mu, int K) {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, D::exp_rate(lambda));
    m.set_service(q, c, D::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// The half that needs nothing but this port
// ---------------------------------------------------------------------------

TEST_CASE("api/sym: the base URL is normalized once and an empty one is refused") {
    // A trailing slash would make every request path a double slash, which the
    // service's router does not match; it is stripped at construction rather
    // than at each call site, so a caller cannot get it wrong.
    CHECK(sym::SageRestEngine("http://localhost:8085/").getBaseUrl() == "http://localhost:8085");
    CHECK(sym::SageRestEngine("  http://localhost:8085///  ").getBaseUrl() ==
          "http://localhost:8085");
    CHECK(sym::SageRestEngine("http://localhost:8085").name() == "sage");
    CHECK_THROWS_AS(sym::SageRestEngine(""), InputError);
    CHECK_THROWS_AS(sym::SageRestEngine("   "), InputError);

    // The timeout travels in the request body as well as bounding the socket
    // read, so it is engine state and not a per-call argument.
    sym::SageRestEngine e("http://localhost:8085");
    CHECK(e.getTimeoutSeconds() == sym::SageRestEngine::DEFAULT_TIMEOUT_SECONDS);
    CHECK(e.setTimeoutSeconds(7).getTimeoutSeconds() == 7);
}

TEST_CASE("api/sym: sym_resolve decides 'none' and 'https' without touching the network") {
    // 'none' is honoured as an ANSWER and not read as "not given": a caller that
    // turned the backend off must never have a container started under it.
    CHECK(sym::sym_resolve("none") == nullptr);
    CHECK(sym::sym_resolve("off") == nullptr);
    CHECK(sym::sym_resolve("NONE") == nullptr);

    // https is REFUSED BY NAME rather than downgraded or reported as "no
    // backend": this port's HTTP client has no TLS, and reporting absence for a
    // service that is up and merely unreachable sends the caller looking in the
    // wrong place.
    CHECK_THROWS_AS(sym::sym_resolve("https://sage.example.org"), UnsupportedError);

    // The same URL through the environment is IGNORED, not raised on: the
    // variable is ambient and may well have been set for another codebase's
    // client, so it must not turn every solve into an error.
    EnvGuard env(sym::SYM_URL_ENV, "https://sage.example.org");
    CHECK_NOTHROW(sym::sym_resolve("none"));
}

TEST_CASE("api/sym: a URL nothing answers resolves to no backend, not to an error") {
    // Port 1 is reserved and nothing binds it, so this is a refused connect and
    // not a hang. `isAvailable` swallows it, and the ladder reports absence --
    // which is what a caller that can fall back to native algebra needs.
    const sym::SageRestEngine dead("http://127.0.0.1:1");
    CHECK_FALSE(dead.isAvailable());
    CHECK(sym::sym_resolve("http://127.0.0.1:1") == nullptr);

    // The refusal a solver hands the user names the image and the variable, so
    // an absent backend is actionable without reading the source.
    qn::Network<double> m = mm1k(0.5, 1.0, 2);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcSymbolicOptions off;
    off.backend = "none";
    try {
        ctmc::ctmc_symbolic_solution(sn, ctmc::CtmcOptions(), off);
        FAIL("ctmc_symbolic_solution must refuse without a backend");
    } catch (const sym::SymEngineError& err) {
        const std::string what = err.what();
        CHECK(what.find(sym::SYM_DOCKER_IMAGE) != std::string::npos);
        CHECK(what.find(sym::SYM_URL_ENV) != std::string::npos);
    }
}

TEST_CASE("api/sym: a coefficient crosses the wire as its shortest round-tripping text") {
    // A COEFFICIENT'S TEXT IS ITS VALUE: the service reads the literal as an
    // exact rational, so 0.1 must arrive as "0.1" and become 1/10, not as the
    // 17-digit expansion of the nearest double. This is also why the printer is
    // shortest-round-tripping and not "%.17g".
    CHECK(sym::detail::decimal_string(0.1) == "0.1");
    CHECK(sym::detail::decimal_string(1.0) == "1");
    CHECK(sym::detail::decimal_string(-2.5) == "-2.5");
    CHECK(std::strtod(sym::detail::decimal_string(1.0 / 3.0).c_str(), nullptr) == 1.0 / 3.0);
    CHECK(std::strtod(sym::detail::decimal_string(1e-8).c_str(), nullptr) == 1e-8);
}

TEST_CASE("api/sym: the encoders refuse a ragged Q and drop the inactive symbols") {
    // A ragged Q is caught HERE and not by the service: a 400 carrying a Python
    // traceback is not a diagnosis of the caller's matrix.
    std::vector<std::vector<std::string>> ragged(2);
    ragged[0].push_back("-x1");
    ragged[0].push_back("x1");
    ragged[1].push_back("x2");
    CHECK_THROWS_AS(sym::detail::to_json_matrix(ragged), InputError);
    CHECK_THROWS_AS(sym::detail::to_json_matrix(std::vector<std::vector<std::string>>()),
                    InputError);

    // An empty ENTRY is a structural zero and is sent as "0"; an empty SYMBOL
    // marks an event with no positive rate and is DROPPED, since the server
    // would otherwise be asked to declare a variable named "".
    std::vector<std::vector<std::string>> holes(1, std::vector<std::string>(1, ""));
    CHECK(sym::detail::to_json_matrix(holes)[0][0].get<std::string>() == "0");
    std::vector<std::string> syms;
    syms.push_back("x1");
    syms.push_back("");
    syms.push_back("x3");
    CHECK(sym::detail::to_json_array(syms).size() == 2);
}

// ---------------------------------------------------------------------------
// The half that needs a live line-sage-rest
// ---------------------------------------------------------------------------

TEST_CASE("api/sym [live]: the service identifies itself and the ladder reaches it") {
    if (no_service()) return;
    const std::string url = live_url();

    // Step 1, an explicit URL. The engine names the ALGEBRA it fronts, not the
    // transport, so a caller can branch on "sage" without knowing about REST.
    const std::shared_ptr<sym::SymEngine> byUrl = sym::sym_resolve(url);
    REQUIRE(byUrl != nullptr);
    CHECK(byUrl->name() == "sage");
    CHECK(byUrl->isAvailable());

    // IDENTITY IS READ, NOT ASSUMED. Every imperialqore line-*-rest listens on
    // 8080 by convention, so the ladder's port probe accepts a service only when
    // /api/v1/info carries sage_version -- a health probe alone would happily
    // take the LQNS one and fail on the first symbolic request.
    sym::SageRestEngine* rest = dynamic_cast<sym::SageRestEngine*>(byUrl.get());
    REQUIRE(rest != nullptr);
    const sym::detail::Json info = rest->info();
    CHECK(info.contains("sage_version"));
    CHECK(sym::detail::is_sage_service(*rest));

    // Step 2, the environment. It is consulted only when the caller named no URL
    // of its own, so an explicit 'none' still wins over an ambient variable.
    {
        EnvGuard env(sym::SYM_URL_ENV, url);
        const std::shared_ptr<sym::SymEngine> byEnv = sym::sym_resolve("auto");
        REQUIRE(byEnv != nullptr);
        CHECK(byEnv->name() == "sage");
        CHECK(sym::sym_resolve("none") == nullptr);
    }
}

TEST_CASE("api/sym [live]: solveCTMC returns the exact stationary law of a two-state chain") {
    if (no_service()) return;
    const std::shared_ptr<sym::SymEngine> engine = live_engine();

    const sym::CtmcSolution sol = engine->solveCTMC(two_state_Q(), two_state_symbols());
    REQUIRE(sol.pi.size() == 2);
    CHECK(sol.nConnComp == 1);

    // pi = (x2, x1)/(x1 + x2). Checked by SUBSTITUTION, not by string: at
    // x1 = 2, x2 = 3 the birth-death balance gives 3/5 and 2/5 exactly, and the
    // service returns exact rationals, so this is an equality and not a
    // tolerance.
    const std::vector<double> pi = engine->eval(sol.pi, two_state_point());
    REQUIRE(pi.size() == 2);
    CHECK(pi[0] == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(pi[1] == doctest::Approx(0.4).epsilon(1e-12));

    // num and den are the same vector over ONE common denominator, which is what
    // a caller differentiating by hand needs; they must agree with pi entry by
    // entry rather than being an independently normalized second answer.
    REQUIRE(sol.num.size() == 2);
    std::vector<std::string> ratio;
    ratio.push_back("(" + sol.num[0] + ")/(" + sol.den + ")");
    ratio.push_back("(" + sol.num[1] + ")/(" + sol.den + ")");
    const std::vector<double> byParts = engine->eval(ratio, two_state_point());
    REQUIRE(byParts.size() == 2);
    CHECK(byParts[0] == doctest::Approx(pi[0]).epsilon(1e-12));
    CHECK(byParts[1] == doctest::Approx(pi[1]).epsilon(1e-12));
}

TEST_CASE("api/sym [live]: ctmcSensitivity returns the scaled and unscaled sensitivity") {
    if (no_service()) return;
    const std::shared_ptr<sym::SymEngine> engine = live_engine();

    // Reward 1 on the up state, so E[r] = x1/(x1 + x2); then dE/dx1 =
    // x2/(x1 + x2)^2 = 3/25 and the scaled one is (x1/E) dE/dx1 = x2/(x1 + x2) =
    // 3/5. Both are closed forms of the fixture, not values read back out.
    std::vector<std::string> reward;
    reward.push_back("0");
    reward.push_back("1");
    const sym::SymSensitivity s =
        engine->ctmcSensitivity(two_state_Q(), two_state_symbols(), "x1", reward);
    CHECK(s.hasReward);
    REQUIRE(s.dpi.size() == 2);

    std::vector<std::string> exprs;
    exprs.push_back(s.S);
    exprs.push_back(s.SS);
    exprs.push_back(s.Er);
    exprs.push_back(s.dpi[0]);
    exprs.push_back(s.dpi[1]);
    const std::vector<double> v = engine->eval(exprs, two_state_point());
    REQUIRE(v.size() == 5);
    CHECK(v[0] == doctest::Approx(0.12).epsilon(1e-12));
    CHECK(v[1] == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(v[2] == doctest::Approx(0.4).epsilon(1e-12));
    // pi sums to one for every x1, so its derivative sums to zero: the
    // normalization is an identity in the symbols and survives differentiation.
    CHECK(v[3] + v[4] == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
}

TEST_CASE("api/sym [live]: simplify, diff and eval are exact, and a free symbol is NaN") {
    if (no_service()) return;
    const std::shared_ptr<sym::SymEngine> engine = live_engine();
    const std::map<std::string, double> at = two_state_point();

    // cancel is a RATIONAL FUNCTION operation: the removable singularity at
    // x1 = x2 must go, leaving something that evaluates to 5 at (2,3).
    std::vector<std::string> ratio;
    ratio.push_back("(x1^2 - x2^2)/(x1 - x2)");
    const std::vector<std::string> cancelled = engine->simplify(ratio, "cancel");
    REQUIRE(cancelled.size() == 1);
    CHECK(engine->eval(cancelled, at)[0] == doctest::Approx(5.0).epsilon(1e-12));

    // The order argument is a REPEATED derivative and not a multi-index: the
    // second derivative of x1^3 x2 in x1 is 6 x1 x2, which is 36 at (2,3).
    std::vector<std::string> cubic;
    cubic.push_back("x1^3*x2");
    const std::vector<std::string> d2 = engine->diff(cubic, "x1", 2);
    REQUIRE(d2.size() == 1);
    CHECK(engine->eval(d2, at)[0] == doctest::Approx(36.0).epsilon(1e-12));
    CHECK(engine->eval(engine->diff(cubic, "x1", 1), at)[0] ==
          doctest::Approx(36.0).epsilon(1e-12));  // 3 x1^2 x2 = 36 at (2,3) too

    // The decimal is read as 1/10 and not as the nearest double, so the sum is
    // the exact 2.1 and not 2.0999999999999996.
    std::vector<std::string> shifted;
    shifted.push_back("x1 + 0.1");
    CHECK(engine->eval(shifted, at)[0] == doctest::Approx(2.1).epsilon(1e-15));

    // AN UNASSIGNED SYMBOL IS NaN, NOT ZERO. A zero would read as a legitimate
    // value and quietly turn a partial substitution into a wrong number.
    std::vector<std::string> free;
    free.push_back("x1 + x3");
    CHECK(std::isnan(engine->eval(free, at)[0]));
}

TEST_CASE("api/sym [live]: fluidODEs returns the Jacobian, the LaTeX and the equilibria") {
    if (no_service()) return;
    const std::shared_ptr<sym::SymEngine> engine = live_engine();

    // dx1/dt = -x1 x2, dx2/dt = x1 x2 - x2. J = [[-x2, -x1], [x2, x1 - 1]],
    // which at (2,3) is [[-3,-2],[3,1]].
    std::vector<std::string> rhs;
    rhs.push_back("-x1*x2");
    rhs.push_back("x1*x2 - x2");
    std::vector<std::string> vars;
    vars.push_back("x1");
    vars.push_back("x2");
    std::vector<std::string> want;
    want.push_back("jacobian");
    want.push_back("latex");
    want.push_back("equilibria");
    const sym::FluidODEs odes = engine->fluidODEs(rhs, vars, want);

    REQUIRE(odes.hasJacobian);
    REQUIRE(odes.jacobian.size() == 2);
    REQUIRE(odes.jacobian[0].size() == 2);
    std::vector<std::string> flat;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) flat.push_back(odes.jacobian[i][j]);
    const std::vector<double> J = engine->eval(flat, two_state_point());
    REQUIRE(J.size() == 4);
    CHECK(J[0] == doctest::Approx(-3.0).epsilon(1e-12));
    CHECK(J[1] == doctest::Approx(-2.0).epsilon(1e-12));
    CHECK(J[2] == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(J[3] == doctest::Approx(1.0).epsilon(1e-12));

    CHECK(odes.hasLatex);
    CHECK(odes.latex.size() == 2);

    // f(x) = 0 here is the whole line x2 = 0, so the solution carries a FREE
    // PARAMETER: `hasEquilibria` distinguishes "asked and answered" from "never
    // asked", which an empty list alone could not.
    CHECK(odes.hasEquilibria);
    REQUIRE(!odes.equilibria.empty());
    CHECK(odes.equilibria[0].count("x1") == 1);
    CHECK(odes.equilibria[0].count("x2") == 1);

    // WANT IS A REQUEST AND NOT A HINT: what was not asked for is absent, so a
    // caller never pays for an equilibria solve it did not want.
    std::vector<std::string> jacOnly;
    jacOnly.push_back("jacobian");
    const sym::FluidODEs bare = engine->fluidODEs(rhs, vars, jacOnly);
    CHECK(bare.hasJacobian);
    CHECK_FALSE(bare.hasEquilibria);
    CHECK_FALSE(bare.hasLatex);
}

TEST_CASE("api/sym [live]: a rejected request comes back as SymEngineError") {
    if (no_service()) return;
    const std::shared_ptr<sym::SymEngine> engine = live_engine();

    // The service parses with an explicit AST walker and never with eval, so a
    // malformed expression is a NAMED refusal carrying the offending text -- not
    // a 500, and not a silently dropped entry that would shorten the result list.
    std::vector<std::string> bad;
    bad.push_back("x1 +* x2");
    CHECK_THROWS_AS(engine->simplify(bad, "cancel"), sym::SymEngineError);

    // A non-square Q is caught client side; an UNKNOWN NORMAL FORM is not
    // something this port can enumerate, so the service decides it and the
    // failure still arrives as the same exception type.
    std::vector<std::string> ok;
    ok.push_back("x1 + x2");
    CHECK_THROWS_AS(engine->simplify(ok, "no-such-form"), sym::SymEngineError);
}

TEST_CASE("ctmc symbolic [live]: getSymbolicSolution at the nominal rates is the numeric pi") {
    if (no_service()) return;

    // M/M/1/2 with lambda = 0.5 and mu = 1: pi is proportional to (1, 1/2, 1/4),
    // i.e. (4, 2, 1)/7. THE ORACLE IS THE NUMERIC ANALYZER, which shares no code
    // with the symbolic path -- it solves pi Q = 0 in doubles here, while the
    // expression was solved over a rational function field on the service.
    qn::Network<double> m = mm1k(0.5, 1.0, 2);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    ctmc::CtmcSymbolicOptions symopt;
    symopt.backend = live_url();

    const ctmc::CtmcSymbolicSolution<double> s = ctmc::ctmc_symbolic_solution(sn, opt, symopt);
    CHECK(s.engine == "sage");
    REQUIRE(s.pi.size() == s.space.size());

    // Substituting each event's nominal rate is what makes the two comparable:
    // the filtration was divided by that rate when the symbol was introduced.
    const ctmc::CtmcSymbolicGenerator<double> g = ctmc::ctmc_symbolic_generator(sn, opt);
    std::map<std::string, double> at;
    for (std::size_t e = 0; e < g.symbols.size(); ++e)
        if (!g.symbols[e].empty()) at[g.symbols[e]] = g.rate0[e];
    const std::vector<double> symbolic = live_engine()->eval(s.pi, at);

    const ctmc::CtmcSolution<double> numeric = ctmc::solver_ctmc_analyzer(sn, opt);
    REQUIRE(symbolic.size() == numeric.pi.size());
    for (std::size_t i = 0; i < symbolic.size(); ++i) {
        CAPTURE(i);
        CHECK(symbolic[i] == doctest::Approx(numeric.pi[i]).epsilon(1e-10));
    }
}

TEST_CASE("ctmc sens [live]: method 'symbolic' agrees with method 'fd'") {
    if (no_service()) return;

    // The two branches share nothing but the model: 'fd' differences two
    // independently solved generators, while 'symbolic' differentiates the
    // rational function the service returned and applies the chain rule over a
    // differenced RATE MAP. Agreement is therefore evidence about both.
    qn::Network<double> m = mm1k(0.5, 1.0, 2);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    ctmc::CtmcSymbolicOptions symopt;
    symopt.backend = live_url();

    ctmc::CtmcSensParam<double> mu;
    mu.name = "mu";
    mu.value = 1.0;
    mu.set = [](qn::NetworkStruct<double>& x, double v) {
        x.service[1][0] = D::exp_rate(v);
        x.refresh_rates();
    };

    // Reward 1 on every state with a job present, so E[r] is the utilization and
    // dE/dmu is negative: a faster server empties the queue.
    const ctmc::CtmcSolution<double> base = ctmc::solver_ctmc_analyzer(sn, opt);
    std::vector<double> reward(base.chain.space.size(), 1.0);
    reward[0] = 0.0;

    const ctmc::CtmcSens<double> sy =
        ctmc::solver_ctmc_sensitivity(sn, opt, mu, reward, "symbolic", symopt);
    const ctmc::CtmcSens<double> fd = ctmc::solver_ctmc_sensitivity(sn, opt, mu, reward, "fd");

    REQUIRE(sy.dpi.size() == fd.dpi.size());
    for (std::size_t i = 0; i < sy.dpi.size(); ++i) {
        CAPTURE(i);
        CHECK(sy.dpi[i] == doctest::Approx(fd.dpi[i]).epsilon(1e-6));
    }
    CHECK(sy.S == doctest::Approx(fd.S).epsilon(1e-6));
    CHECK(sy.SS == doctest::Approx(fd.SS).epsilon(1e-6));
    CHECK(sy.S < 0.0);

    // pi sums to one for every mu, so dpi sums to zero -- the residual of the
    // normalizing equation, read back out of the symbolic branch.
    double total = 0.0;
    for (std::size_t i = 0; i < sy.dpi.size(); ++i) total += sy.dpi[i];
    CHECK(total == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
}

// ---------------------------------------------------------------------------
// Step 4 of the ladder, which starts a container and is therefore opt-in
// ---------------------------------------------------------------------------

TEST_CASE("api/sym [docker]: sym_resolve starts the service itself and stops it again") {
    // OPT-IN BY NAME. This case starts a 3 GB image, so it must never run as a
    // side effect of an unrelated suite; and a pull is refused here in any case,
    // since the resolution below passes "auto" rather than the "sage" keyword.
    if (std::getenv("LINE_SYM_TEST_DOCKER") == nullptr) {
        MESSAGE("LINE_SYM_TEST_DOCKER is not set: skipping the container case");
        return;
    }
    if (!io::docker_daemon_available()) {
        MESSAGE("no Docker daemon: skipping the container case");
        return;
    }
    const std::string image = sym::sym_find_image();
    if (image.empty()) {
        MESSAGE("no local " << std::string(sym::SYM_DOCKER_IMAGE)
                            << ": skipping the container case (docker pull it to run this)");
        return;
    }
    // A service already listening would be taken at step 3, so there would be no
    // container to start OR to stop, and the assertions below would be about
    // somebody else's process.
    if (live_engine()) {
        MESSAGE("a line-sage-rest is already reachable: skipping the container case, which "
                "asserts on a container this process owns");
        return;
    }

    EnvGuard clearNamed(TEST_URL_ENV);
    EnvGuard clearUrl(sym::SYM_URL_ENV);
    const std::shared_ptr<sym::SymEngine> engine = sym::sym_resolve("auto");
    REQUIRE(engine != nullptr);
    CHECK(engine->name() == "sage");
    CHECK(engine->isAvailable());

    // The port is EPHEMERAL, so a second process, or a hand-started service on
    // the conventional port, does not collide with this one.
    const sym::SageRestEngine* rest = dynamic_cast<const sym::SageRestEngine*>(engine.get());
    REQUIRE(rest != nullptr);
    CHECK(rest->getBaseUrl().compare(0, 17, "http://localhost:") == 0);
    CHECK(rest->getBaseUrl() != "http://localhost:8080");

    std::vector<std::string> ex;
    ex.push_back("(x^2 - 1)/(x - 1)");
    const std::vector<std::string> cancelled = engine->simplify(ex, "cancel");
    REQUIRE(cancelled.size() == 1);
    std::map<std::string, double> at;
    at["x"] = 4.0;
    CHECK(engine->eval(cancelled, at)[0] == doctest::Approx(5.0).epsilon(1e-12));

    // The same call the atexit handler makes. It is exercised HERE because
    // atexit does NOT run on a signal, so a killed suite would otherwise leave
    // the container behind and nothing would have checked that stopping works.
    sym::sym_stop_container();
    CHECK_FALSE(engine->isAvailable());
}
