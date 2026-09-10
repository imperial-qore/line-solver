/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * End-to-end tests of `solver_nc_lossn_analyzer`: a Network with a Finite
 * Capacity Region in, an NC metric table out, with nothing hand-assembled in
 * between.
 *
 * THE ERLANG FIXED POINT IS AN APPROXIMATION -- it assumes the rows of the
 * admission rule block independently -- so nothing here compares it to an exact
 * solver on a region whose rows share a class. It is compared to a closed form
 * only in the two cases where the independence assumption is VACUOUS and the
 * fixed point provably collapses to Erlang's loss formula:
 *
 *   1. one row, one class: the reduced load carries the factor (1-E)^1 that the
 *      division by (1-E) cancels, so E = B(nu, C) exactly;
 *   2. per-class rows only: each row sees exactly one class, so the same
 *      cancellation applies row by row and Loss_r = B(nu_r, C_r) exactly.
 *
 * Erlang B itself is computed in the test by its rational recursion
 * B_k = nu B_{k-1} / (k + nu B_{k-1}), B_0 = 1, which is an independent
 * evaluation and not a number read back out of the implementation. The one
 * literal is the MATLAB-printed value of B(4,5), cited where it appears.
 *
 * Everything else is an IDENTITY of the assembly -- Little's law at the
 * infinite server, the complementarity of carried and offered load, monotonicity
 * in the offered load, and the vanishing of blocking as the capacity grows --
 * or a refusal that must arrive by name.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/nc/solver_nc_lossn.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/util/error.h"

namespace qn = line::qn;
namespace nc = line::nc;
using line::Rational;
using line::lang::DropStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Erlang B by its rational recursion, exact for an integer capacity. */
double erlangB(double nu, int C) {
    double b = 1.0;
    for (int k = 1; k <= C; ++k) b = nu * b / (k + nu * b);
    return b;
}

/**
 * Source -> Delay -> Sink with `K` open classes, the Delay alone inside one
 * DROP region. `classCap[r]` and `globalCap` use the -1 sentinel for unbounded.
 */
qn::Network<double> lossnet(const std::vector<double>& lambda, const std::vector<double>& mu,
                            const std::vector<double>& classCap, double globalCap) {
    qn::Network<double> m("lossnet");
    const std::size_t src = m.add_source("Src");
    const std::size_t d = m.add_delay("Region");
    const std::size_t snk = m.add_sink("Sink");
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 0; r < lambda.size(); ++r) {
        const std::size_t c = m.add_open_class("C" + std::to_string(r + 1));
        m.set_arrival(src, c, Dist::exp_rate(lambda[r]));
        m.set_service(d, c, Dist::exp_rate(mu[r]));
        P.set(c, c, src, d, 1.0);
        P.set(c, c, d, snk, 1.0);
    }
    m.link(P);
    m.add_region({d}, classCap, globalCap,
                 std::vector<DropStrategy>(lambda.size(), DropStrategy::DROP));
    return m;
}

}  // namespace

TEST_CASE("nclossn: a single-row region reproduces Erlang B exactly") {
    // nu = lambda/mu = 4, capacity 5. B(4,5) as MATLAB prints it.
    const double lambda = 2.0, mu = 0.5;
    const int C = 5;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, static_cast<double>(C));
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(nc::nc_is_lossn_model(sn));

    nc::NcSolverOptions opt;
    opt.method = "erlangfp";
    const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(sn, opt);

    const double nu = lambda / mu;
    CHECK(r.nu[0] == doctest::Approx(nu).epsilon(1e-12));
    CHECK(erlangB(nu, C) == doctest::Approx(0.199066874027994).epsilon(1e-12));  // MATLAB
    CHECK(r.Loss[0] == doctest::Approx(erlangB(nu, C)).epsilon(1e-7));

    // One row only: the global job cap, right-hand side C.
    REQUIRE(r.A.rows() == 1);
    CHECK(r.A(0, 0) == doctest::Approx(1.0));
    CHECK(r.Cvec[0] == doctest::Approx(static_cast<double>(C)));

    // Little's law at the infinite server, and the source emitting the
    // post-drop rate.
    const std::size_t delay = 2 - 1;  // Source is station 1, Delay station 2
    const double Xc = lambda * (1.0 - r.Loss[0]);
    CHECK(r.sol.sol.X[0] == doctest::Approx(Xc).epsilon(1e-12));
    CHECK(r.sol.sol.Q(delay, 0) == doctest::Approx(nu * (1.0 - r.Loss[0])).epsilon(1e-9));
    CHECK(r.sol.sol.Q(delay, 0) == doctest::Approx(Xc / mu).epsilon(1e-9));
    CHECK(r.sol.sol.U(delay, 0) == doctest::Approx(r.sol.sol.Q(delay, 0)).epsilon(1e-12));
    CHECK(r.sol.sol.R(delay, 0) == doctest::Approx(1.0 / mu).epsilon(1e-12));
    CHECK(r.sol.sol.Tp(delay, 0) == doctest::Approx(Xc).epsilon(1e-12));
    CHECK(r.sol.sol.Tp(0, 0) == doctest::Approx(Xc).epsilon(1e-12));
    CHECK(r.sol.actualmethod == "lossn.erlangfp");
    // The fixed point carries no normalizing constant, and reports none.
    CHECK(std::isnan(r.sol.sol.lG));
}

TEST_CASE("nclossn: per-class caps give per-class Erlang B and no coupling") {
    // Two classes, no global cap, so each class owns exactly one row of A and
    // the fixed point decouples into two independent Erlang systems.
    const std::vector<double> lambda{3.0, 3.0}, mu{1.0, 1.0};
    qn::Network<double> m = lossnet(lambda, mu, {2.0, 12.0}, -1.0);
    nc::NcSolverOptions opt;
    opt.method = "erlangfp";
    const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(m.get_struct(), opt);

    REQUIRE(r.A.rows() == 2);
    CHECK(r.Loss[0] == doctest::Approx(erlangB(3.0, 2)).epsilon(1e-7));
    CHECK(r.Loss[1] == doctest::Approx(erlangB(3.0, 12)).epsilon(1e-7));
    CHECK(r.Loss[0] > r.Loss[1]);  // the tighter cap blocks more
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(r.Loss[c] >= 0.0);
        CHECK(r.Loss[c] <= 1.0);
        CHECK(r.sol.sol.X[c] == doctest::Approx(lambda[c] * (1.0 - r.Loss[c])).epsilon(1e-12));
    }
}

TEST_CASE("nclossn: blocking rises with the offered load") {
    double prev = -1.0;
    for (double lambda : {0.5, 2.0, 6.0, 20.0}) {
        qn::Network<double> m = lossnet({lambda}, {1.0}, {-1.0}, 4.0);
        nc::NcSolverOptions opt;
        opt.method = "erlangfp";
        const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(m.get_struct(), opt);
        CHECK(r.Loss[0] >= 0.0);
        CHECK(r.Loss[0] <= 1.0);
        CHECK(r.Loss[0] > prev);
        prev = r.Loss[0];
    }
}

TEST_CASE("nclossn: capacity far above the offered load blocks nothing") {
    // Not the unbounded region, which has no admission rule at all and is
    // refused below; this is the limit approached from inside the model.
    const double lambda = 1.0, mu = 1.0;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, 100.0);
    nc::NcSolverOptions opt;
    opt.method = "erlangfp";
    const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(m.get_struct(), opt);
    CHECK(r.Loss[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.sol.sol.X[0] == doctest::Approx(lambda).epsilon(1e-12));
    CHECK(r.sol.sol.Q(1, 0) == doctest::Approx(lambda / mu).epsilon(1e-9));
}

TEST_CASE("nclossn: a region with no bounded row is refused, not solved") {
    qn::Network<double> m = lossnet({1.0}, {1.0}, {-1.0}, -1.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(nc::nc_is_lossn_model(sn));  // it is shaped like one
    nc::NcSolverOptions opt;
    opt.method = "erlangfp";
    CHECK_THROWS_AS(nc::solver_nc_lossn_analyzer(sn, opt), line::UnsupportedError);
}

TEST_CASE("nclossn: the exact transform is the default on an integral rule") {
    // nu = 4, capacity 5, one row: the transform must return Erlang B itself,
    // not the fixed point's approximation to it, and it must be reached without
    // being asked for -- 'default' resolves to it whenever the rule is integral,
    // as in the reference.
    const double lambda = 2.0, mu = 0.5;
    const int C = 5;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, static_cast<double>(C));
    const qn::NetworkStruct<double>& sn = m.get_struct();

    for (const std::string& meth :
         {std::string("exact"), std::string("ms"), std::string("default")}) {
        nc::NcSolverOptions opt;
        opt.method = meth;
        const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(sn, opt);
        INFO("method ", meth);
        CHECK(r.sol.actualmethod == "lossn.exact");
        CHECK(r.Loss[0] == doctest::Approx(erlangB(4.0, C)).epsilon(1e-13));
        CHECK(r.sol.sol.iter == 1);  // direct, not iterative
        CHECK(r.E.empty());          // no per-row blocking outside the fixed point
        // g(C) = sum_{n=0}^{C} nu^n/n! for the single-row Erlang system, which
        // the fixed point cannot report at all.
        double g = 0.0, t = 1.0;
        for (int n = 0; n <= C; ++n) {
            if (n > 0) t *= 4.0 / n;
            g += t;
        }
        CHECK(r.sol.sol.lG == doctest::Approx(std::log(g)).epsilon(1e-12));
        // Little's law at the infinite server still holds on this path.
        CHECK(r.sol.sol.Q(1, 0) == doctest::Approx(r.sol.sol.X[0] / mu).epsilon(1e-12));
    }
}

TEST_CASE("nclossn: the exact transform and the fixed point differ on coupled rows") {
    // Two classes under a global job cap AND a per-class cap, so both rows see
    // class 1 and the independence assumption behind the Erlang fixed point is
    // false. The point of the test is that the two methods are NOT
    // interchangeable: the transform is exact and the approximation is off by a
    // visible margin, which is why 'exact' must never silently downgrade.
    const std::vector<double> lambda{3.0, 2.0}, mu{1.0, 1.0};
    qn::Network<double> m = lossnet(lambda, mu, {3.0, -1.0}, 5.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    nc::NcSolverOptions ex;
    ex.method = "exact";
    const nc::NcLossnSolution<double> e = nc::solver_nc_lossn_analyzer(sn, ex);
    nc::NcSolverOptions fp;
    fp.method = "erlangfp";
    const nc::NcLossnSolution<double> a = nc::solver_nc_lossn_analyzer(sn, fp);

    REQUIRE(e.A.rows() == 2);
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(e.Loss[r] > 0.0);
        CHECK(e.Loss[r] < 1.0);
    }
    // The tighter-capped class blocks more under both, but not by the same
    // amount; a tolerance that made these agree would be hiding the difference.
    CHECK(e.Loss[0] > e.Loss[1]);
    CHECK(std::fabs(e.Loss[0] - a.Loss[0]) > 1e-3);

    // The sampler is a third, independent algorithm and must bracket the exact
    // value it does not share a line of code with.
    nc::NcSolverOptions mc;
    mc.method = "mci";
    mc.samples = 200000;
    mc.seed = 7;
    const nc::NcLossnSolution<double> s = nc::solver_nc_lossn_analyzer(sn, mc);
    for (std::size_t r = 0; r < 2; ++r) {
        INFO("class ", r, ": exact ", e.Loss[r], " mci ", s.Loss[r]);
        CHECK(s.Loss[r] == doctest::Approx(e.Loss[r]).epsilon(2e-2));
    }
    CHECK(s.sol.sol.lG == doctest::Approx(e.sol.sol.lG).epsilon(1e-2));
}

TEST_CASE("nclossn: an explicit exact request on a fractional rule is refused") {
    // A memory budget with a fractional class size. 'default' downgrades to
    // erlangfp (which then refuses for its own reason), but an EXPLICIT 'exact'
    // must reach lossn_manjunath and be refused there: the residue argument counts
    // whole units of capacity, and answering with an approximation under the
    // name of an exact method is the failure this guards.
    qn::Network<double> m("lossnet-frac-exact");
    const std::size_t src = m.add_source("Src");
    const std::size_t d = m.add_delay("Region");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(1.0));
    m.set_service(d, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, d, 1.0);
    P.set(d, snk, 1.0);
    m.link(P);
    m.add_region({d}, {-1.0}, -1.0, {DropStrategy::DROP}, {}, {2.5}, 10.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    nc::NcSolverOptions ex;
    ex.method = "exact";
    CHECK_THROWS_AS(nc::solver_nc_lossn_analyzer(sn, ex), line::InputError);
}

TEST_CASE("nclossn: mci agrees with Erlang B and carries a normalizing constant") {
    const double lambda = 4.0, mu = 1.0;
    const int C = 5;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, static_cast<double>(C));
    nc::NcSolverOptions opt;
    opt.method = "mci";
    opt.samples = 50000;
    opt.seed = 2026;
    const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(m.get_struct(), opt);

    CHECK(r.Loss[0] == doctest::Approx(erlangB(4.0, C)).epsilon(5e-2));
    // g(C) = sum_{n=0}^{C} nu^n / n! for the single-link Erlang system.
    double g = 0.0, fact = 1.0;
    for (int n = 0; n <= C; ++n) {
        if (n > 0) fact *= n;
        g += std::pow(4.0, n) / fact;
    }
    CHECK(r.sol.sol.lG == doctest::Approx(std::log(g)).epsilon(1e-2));
    CHECK(r.sol.actualmethod == "lossn.mci");
    CHECK(r.sol.sol.X[0] == doctest::Approx(lambda * (1.0 - r.Loss[0])).epsilon(1e-12));
    CHECK(r.E.empty());  // no link blocking vector outside the fixed point
}

TEST_CASE("nclossn: a fractional admission rule is refused under erlangfp") {
    // A memory budget with fractional class sizes makes one row of A and its
    // right-hand side fractional; the ported Erlang B takes an integer capacity
    // and integer circuit requirements, so it must refuse rather than truncate.
    qn::Network<double> m("lossnet-frac");
    const std::size_t src = m.add_source("Src");
    const std::size_t d = m.add_delay("Region");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(1.0));
    m.set_service(d, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, d, 1.0);
    P.set(d, snk, 1.0);
    m.link(P);
    m.add_region({d}, {-1.0}, -1.0, {DropStrategy::DROP}, {}, {2.5}, 10.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    nc::NcSolverOptions ef;
    ef.method = "erlangfp";
    CHECK_THROWS_AS(nc::solver_nc_lossn_analyzer(sn, ef), line::UnsupportedError);

    // The sampler compares in real arithmetic and accepts the same region.
    nc::NcSolverOptions mc;
    mc.method = "mci";
    mc.samples = 20000;
    const nc::NcLossnSolution<double> r = nc::solver_nc_lossn_analyzer(sn, mc);
    CHECK(r.Loss[0] >= 0.0);
    CHECK(r.Loss[0] <= 1.0);
    CHECK(r.Cvec[0] == doctest::Approx(10.0));
    CHECK(r.A(0, 0) == doctest::Approx(2.5));
}

TEST_CASE("nclossn: the exact transform is reached through solver_nc_run_analyzer, not only directly") {
    // Everything above calls the analyzer directly, which is how the wiring gap
    // of 2026-07-28 stayed invisible. This goes through the dispatch, so the
    // interception, the residual region gates and the analyzer are all on the
    // path -- and it is a `default` request, the one the reference resolves to
    // the transform.
    const double lambda = 2.0, mu = 0.5;
    const int C = 5;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, static_cast<double>(C));
    nc::NcSolverOptions opt;  // method = "default"
    const line::mva::AvgResult<double> r = nc::solver_nc_run_analyzer(m.get_struct(), opt);

    // The runner's 'default/<algorithm>' convention, so the algorithm that
    // produced the numbers survives a default request.
    CHECK(r.method == "default");
    CHECK(r.actualmethod == "default/lossn.exact");
    const double nu = lambda / mu;
    const double loss = erlangB(nu, C);
    CHECK(r.XN[0] == doctest::Approx(lambda * (1.0 - loss)).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(nu * (1.0 - loss)).epsilon(1e-12));
    CHECK(r.RN(1, 0) == doctest::Approx(1.0 / mu).epsilon(1e-12));
}

TEST_CASE("nclossn: 'ms' reaches the exact transform through solver_nc_run_analyzer") {
    // 'ms' names the same lossn_manjunath transform as 'exact' (solver_nc_lossn.h's
    // own method name alias) and was already accepted by solver_nc_lossn_analyzer
    // directly, but list_valid_methods omitted it, so a caller going through
    // the ordinary entry point solver_nc_run_analyzer got "method 'ms' is unsupported
    // by this solver" instead of an answer. Pinning it through the same public
    // path the 'default' wiring test above uses, not only the analyzer.
    const double lambda = 2.0, mu = 0.5;
    const int C = 5;
    qn::Network<double> m = lossnet({lambda}, {mu}, {-1.0}, static_cast<double>(C));
    nc::NcSolverOptions opt;
    opt.method = "ms";
    const line::mva::AvgResult<double> r = nc::solver_nc_run_analyzer(m.get_struct(), opt);

    CHECK(r.method == "ms");
    CHECK(r.actualmethod == "lossn.exact");
    const double nu = lambda / mu;
    const double loss = erlangB(nu, C);
    CHECK(r.XN[0] == doctest::Approx(lambda * (1.0 - loss)).epsilon(1e-12));
    CHECK(r.QN(1, 0) == doctest::Approx(nu * (1.0 - loss)).epsilon(1e-12));
}

TEST_CASE("nclossn: the exact transform answers under exact arithmetic") {
    // The analyzer used to be gated whole on has_transcendental, so a loss
    // network was refused outright under --arith exact. Only the two inexact
    // methods are gated now: the transform is rational throughout and the
    // metrics are ratios, so Loss is an exact rational here.
    qn::Network<Rational> m("lossnet-exact");
    const std::size_t src = m.add_source("Src");
    const std::size_t d = m.add_delay("Region");
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, line::lang::Distrib<Rational>::exp_rate(Rational(2)));
    m.set_service(d, c, line::lang::Distrib<Rational>::exp_rate(Rational(1, 2)));
    qn::RoutingMatrix<Rational> P;
    P.set(src, d, Rational(1));
    P.set(d, snk, Rational(1));
    m.link(P);
    m.add_region({d}, {-1.0}, 5.0, {DropStrategy::DROP});
    const qn::NetworkStruct<Rational>& sn = m.get_struct();
    REQUIRE(nc::nc_is_lossn_model(sn));

    nc::NcSolverOptions ex;
    ex.method = "exact";
    const nc::NcLossnSolution<Rational> r = nc::solver_nc_lossn_analyzer(sn, ex);

    // Erlang B(4, 5) as an exact rational, by the same recursion the double
    // tests use: no tolerance is involved on either side of this comparison.
    Rational b(1);
    for (int k = 1; k <= 5; ++k) b = Rational(4) * b / (Rational(k) + Rational(4) * b);
    CHECK(r.Loss[0] == b);
    CHECK(r.nu[0] == Rational(4));
    CHECK(r.sol.actualmethod == "lossn.exact");
    // and it is the same number the double instantiation reports.
    CHECK(static_cast<double>(r.Loss[0]) == doctest::Approx(erlangB(4.0, 5)).epsilon(1e-13));

    // The two inexact methods refuse by name rather than pretending.
    for (const std::string& meth : {std::string("erlangfp"), std::string("mci")}) {
        nc::NcSolverOptions bad;
        bad.method = meth;
        INFO("method ", meth);
        CHECK_THROWS_AS(nc::solver_nc_lossn_analyzer(sn, bad), line::UnsupportedError);
    }
}

TEST_CASE("nclossn: the shape test rejects a model that is not a loss network") {
    // No region at all.
    qn::Network<double> plain("plain");
    const std::size_t src = plain.add_source("Src");
    const std::size_t q = plain.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t snk = plain.add_sink("Sink");
    const std::size_t c = plain.add_open_class("C1");
    plain.set_arrival(src, c, Dist::exp_rate(0.5));
    plain.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    plain.link(P);
    CHECK_FALSE(nc::nc_is_lossn_model(plain.get_struct()));

    // A region that holds jobs back instead of discarding them is a blocking
    // network, not a loss network, and must not be routed here.
    qn::Network<double> hold("holdnet");
    const std::size_t s2 = hold.add_source("Src");
    const std::size_t d2 = hold.add_delay("Region");
    const std::size_t k2 = hold.add_sink("Sink");
    const std::size_t c2 = hold.add_open_class("C1");
    hold.set_arrival(s2, c2, Dist::exp_rate(1.0));
    hold.set_service(d2, c2, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> P2;
    P2.set(s2, d2, 1.0);
    P2.set(d2, k2, 1.0);
    hold.link(P2);
    hold.add_region({d2}, {-1.0}, 4.0, {DropStrategy::WAITQ});
    CHECK_FALSE(nc::nc_is_lossn_model(hold.get_struct()));
}
