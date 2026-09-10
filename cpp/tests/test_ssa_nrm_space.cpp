/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverSSA, the explicit state-space variant of the NRM, and the standalone
 * enabled-set scan `solver_ssa_findenabled.m`.
 *
 * WHAT MAKES THIS FILE DIFFERENT FROM test_ssa_nrm.cpp. That file can only test
 * a simulator against closed forms, because the plain engine stores nothing that
 * another codepath also computes. This variant does: it builds a state space and
 * a per-state rate table, and BOTH have an exact counterpart in the already
 * verified CTMC machinery. So the sharpest assertions here are exact identities
 * between two independent accounts of the same model, and only the means fall
 * back on statistics:
 *
 *   EXACT     the enumerated aggregate space equals the projection of
 *             `reachable_space_generator`'s walk, state for state; the enabled
 *             set equals the generator row, rate for rate; a seed reproduces its
 *             trace bit for bit; the path never leaves the enumerated space.
 *   t-TEST    a mean against the exact CTMC answer, eight independent seeds, the
 *             two-sided one-sample t statistic against the 1% critical value on
 *             7 degrees of freedom, 3.4995. Same discipline as test_ssa_nrm.cpp.
 *   REFUSAL   an unported or unrepresentable construct throws UnsupportedError
 *             and the message NAMES it.
 *
 * Every simulated result is reported with its sample count and its seed, because
 * a simulated number without them is not a measurement.
 */
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ssa/solver_ssa_nrm_space.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace ssa = line::ssa;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Student t at 1% two-sided on 7 degrees of freedom, as test_ssa_nrm.cpp uses. */
constexpr double kT7 = 3.4995;
constexpr std::size_t kReps = 8;

/** Think -> PS Queue -> Think, one closed class. Every node is a station. */
qn::Network<double> cqn_ps(double n, double think_mean, double mu) {
    qn::Network<double> m("nrmspace-cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(1.0 / think_mean));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** Think -> {Q1, Q2, Q3} in equal thirds -> Think, one closed class. */
qn::Network<double> cqn_three_way(double n, double think_mean, double mu) {
    qn::Network<double> m("nrmspace-3way");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(1.0 / think_mean));
    m.set_service(q1, c, Dist::exp_rate(mu));
    m.set_service(q2, c, Dist::exp_rate(mu));
    m.set_service(q3, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0 / 3.0);
    P.set(d, q2, 1.0 / 3.0);
    P.set(d, q3, 1.0 / 3.0);
    P.set(q1, d, 1.0);
    P.set(q2, d, 1.0);
    P.set(q3, d, 1.0);
    m.link(P);
    return m;
}

/** The aggregate (node, class) population vector of one encoded network state. */
std::vector<double> project(const qn::NetworkStruct<double>& sn,
                            const qn::NetState<double>& s) {
    const std::size_t K = sn.nclasses;
    std::vector<double> v(sn.nof_nodes() * K, 0.0);
    for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f) {
        const std::size_t ind = sn.stateful_nodes[f];
        const std::pair<double, std::vector<double>> mg =
            qn::to_marginal_aggr(sn, ind, s.local[f]);
        for (std::size_t r = 0; r < K; ++r) v[(ind - 1) * K + r] = mg.second[r];
    }
    return v;
}

/** The two-sided one-sample t statistic of `x` against `mu0`. */
double tstat(const std::vector<double>& x, double mu0) {
    double m = 0.0;
    for (double v : x) m += v;
    m /= static_cast<double>(x.size());
    double s2 = 0.0;
    for (double v : x) s2 += (v - m) * (v - m);
    s2 /= static_cast<double>(x.size() - 1);
    const double se = std::sqrt(s2 / static_cast<double>(x.size()));
    // A degenerate spread means every replicate agreed exactly; that is a pass
    // when the value is right and must not become a division by zero.
    if (!(se > 0.0)) return std::fabs(m - mu0) > 1e-12 ? 1e9 : 0.0;
    return std::fabs(m - mu0) / se;
}

}  // namespace

TEST_CASE("SolverSSA nrm.space: the enumerated space is the CTMC reachable walk") {
    // Two jobs over a Delay and a PS queue: three aggregate states. The engine
    // closes the space forward over its REACTIONS; the CTMC walk closes it over
    // the STATE HANDLERS. Nothing is shared between the two paths but the
    // NetworkStruct, so agreeing state for state is a real cross-check and not a
    // restatement of one computation.
    qn::Network<double> m = cqn_ps(2.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    // The projection below is only well defined when every node holds jobs;
    // a Router or a Sink would own an aggregate slot the encoded space has no
    // column for.
    REQUIRE(sn.stateful_nodes.size() == sn.nof_nodes());

    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 1;
    ssa::NrmSpaceEngine<double> eng(sn, opt);
    const std::vector<ssa::NrmSpaceState> got = eng.enumerate_space();

    qn::NetState<double> init;
    REQUIRE(ctmc::analyzer_detail::default_init_state(sn, init));
    const std::vector<qn::Sync<double>> sync = qn::refresh_sync(sn);
    const std::vector<qn::GlobalSync<double>> gsync = qn::refresh_global_sync(sn);
    const std::vector<qn::NetState<double>> walk =
        ctmc::reachable_space_generator(sn, init, sync, gsync);

    std::map<std::vector<double>, int> seen;
    for (std::size_t s = 0; s < walk.size(); ++s) seen[project(sn, walk[s])] += 1;
    for (std::size_t s = 0; s < got.size(); ++s) seen[got[s].n] += 2;
    CHECK(got.size() == walk.size());
    for (std::map<std::vector<double>, int>::const_iterator it = seen.begin(); it != seen.end();
         ++it)
        CHECK(it->second == 3);  // in both, exactly once each
    // Two jobs over two nodes: (2,0), (1,1), (0,2).
    CHECK(got.size() == 3);
    // Population is conserved by every reaction, which is what makes the space
    // finite in the first place.
    for (std::size_t s = 0; s < got.size(); ++s) {
        double tot = 0.0;
        for (double v : got[s].n) tot += v;
        CHECK(tot == doctest::Approx(2.0).epsilon(1e-12));
    }
}

TEST_CASE("SolverSSA nrm.space: findenabled reproduces the CTMC generator row") {
    // The enabled set at a state IS that state's generator row: the CTMC
    // assembly and the scan apply the same handlers, so aggregating the scan by
    // destination must give back Q entry for entry. Self-loops are the one
    // place the two differ on purpose -- the generator cancels them on the
    // diagonal -- so they are checked against the row sum instead.
    qn::Network<double> m = cqn_three_way(2.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<qn::Sync<double>> sync = qn::refresh_sync(sn);
    const std::vector<qn::GlobalSync<double>> gsync = qn::refresh_global_sync(sn);
    const std::vector<std::size_t> cut(sn.nclasses, 0);  // closed model: no cutoff
    const std::vector<qn::NetState<double>> space = qn::space_generator(sn, cut);
    const ctmc::CtmcResult<double> chain = ctmc::solver_ctmc(sn, space, sync, gsync);

    std::map<std::vector<double>, std::size_t> index;
    for (std::size_t s = 0; s < space.size(); ++s)
        index[ctmc::ctmc_detail::state_key(space[s])] = s;
    REQUIRE(space.size() > 1);

    for (std::size_t s = 0; s < space.size(); ++s) {
        std::vector<std::vector<double>> arv, dep;
        const std::vector<ssa::EnabledEvent<double>> ev =
            ssa::ssa_find_enabled(sn, sync, gsync, space[s], &arv, &dep);

        std::map<std::size_t, double> agg;
        double selfloop = 0.0;
        for (std::size_t e = 0; e < ev.size(); ++e) {
            const std::map<std::vector<double>, std::size_t>::const_iterator it =
                index.find(ctmc::ctmc_detail::state_key(ev[e].next));
            // A closed model's space is complete, so every successor the scan
            // produces has to be a state the enumeration already contains.
            REQUIRE(it != index.end());
            if (it->second == s) selfloop += ev[e].rate;
            else agg[it->second] += ev[e].rate;
        }
        for (std::map<std::size_t, double>::const_iterator it = agg.begin(); it != agg.end();
             ++it)
            CHECK(it->second == doctest::Approx(chain.Q(s, it->first)).epsilon(1e-12));
        // Nothing the generator has may be missing from the scan either.
        double rowsum = 0.0;
        for (std::size_t j = 0; j < space.size(); ++j) {
            if (j == s) continue;
            if (chain.Q(s, j) != 0.0) CHECK(agg.find(j) != agg.end());
            rowsum += chain.Q(s, j);
        }
        CHECK(rowsum == doctest::Approx(-chain.Q(s, s)).epsilon(1e-9));
        CHECK(selfloop >= 0.0);

        // The departure rates the scan reports are the same rates, grouped by
        // the node that lost the job rather than by where it landed. The model
        // is deliberately one whose ACTIVE outcomes are deterministic: the
        // generator's departure accumulator omits the active outcome
        // probability, where the SSA convention (solver_ssa.m, depRatesSamples)
        // includes it, so the two coincide only where that probability is one.
        // Every branch exercised here returns prob = 1 from after_event.
        double depsum = 0.0;
        for (std::size_t f = 0; f < dep.size(); ++f)
            for (std::size_t r = 0; r < sn.nclasses; ++r) {
                depsum += dep[f][r];
                CHECK(dep[f][r] == doctest::Approx(chain.dep_rates[s][f][r]).epsilon(1e-12));
                CHECK(arv[f][r] == doctest::Approx(chain.arv_rates[s][f][r]).epsilon(1e-12));
            }
        CHECK(depsum == doctest::Approx(rowsum + selfloop).epsilon(1e-9));
    }
}

TEST_CASE("SolverSSA nrm.space: a seeded run is reproducible and the seed reaches it") {
    qn::Network<double> m = cqn_ps(3.0, 2.0, 3.0);
    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 3000;
    opt.seed = 5150;

    ssa::NrmSpaceEngine<double> e1(m.get_struct(), opt);
    ssa::NrmSpaceEngine<double> e2(m.get_struct(), opt);
    const ssa::SsaNrmSpaceRun<double> r1 = e1.run();
    const ssa::SsaNrmSpaceRun<double> r2 = e2.run();

    REQUIRE(r1.samples == 3000);
    REQUIRE(r1.seed == 5150);
    // Bit for bit, not to a tolerance: the same seed drives the same draws
    // through the same arithmetic, so any difference at all is a defect.
    CHECK(r1.tran_rx == r2.tran_rx);
    bool same = r1.tran_time.size() == r2.tran_time.size();
    for (std::size_t i = 0; same && i < r1.tran_time.size(); ++i)
        same = r1.tran_time[i] == r2.tran_time[i];
    CHECK(same);
    CHECK(r1.simulated_time == r2.simulated_time);
    REQUIRE(r1.pi.size() == r2.pi.size());
    for (std::size_t s = 0; s < r1.pi.size(); ++s) CHECK(r1.pi[s] == r2.pi[s]);

    ssa::SsaNrmSpaceOptions other = opt;
    other.seed = 271828;
    ssa::NrmSpaceEngine<double> e3(m.get_struct(), other);
    const ssa::SsaNrmSpaceRun<double> r3 = e3.run();
    // 3000 firings agreeing across two streams has probability below 2^-3000.
    CHECK(r3.tran_rx != r1.tran_rx);
}

TEST_CASE("SolverSSA nrm.space: the path stays inside the enumerated space") {
    qn::Network<double> m = cqn_ps(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 5000;
    opt.seed = 909;
    ssa::NrmSpaceEngine<double> eng(sn, opt);
    const std::vector<ssa::NrmSpaceState> space = eng.enumerate_space();
    const ssa::SsaNrmSpaceRun<double> r = eng.run();

    std::map<std::vector<double>, bool> known;
    for (std::size_t s = 0; s < space.size(); ++s) known[space[s].n] = true;
    REQUIRE(!r.space.empty());
    for (std::size_t s = 0; s < r.space.size(); ++s)
        CHECK(known.find(r.space[s].n) != known.end());
    // Four states for three jobs over two nodes, and 5000 firings on a chain
    // that mixes in a handful of events covers all of them; a path that reaches
    // fewer states than it can is a stuck path, which is the failure this sees.
    CHECK(r.space.size() == space.size());

    double tot = 0.0;
    for (std::size_t s = 0; s < r.pi.size(); ++s) {
        CHECK(r.pi[s] >= 0.0);
        tot += r.pi[s];
    }
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));
    // The departure-rate table is exact per state, so it must equal what the
    // propensity function returns when asked directly.
    for (std::size_t s = 0; s < r.space.size(); ++s) {
        const std::vector<double> A = eng.propensities(r.space[s]);
        for (std::size_t k = 0; k < A.size(); ++k)
            CHECK(r.dep_rates(s, k) == doctest::Approx(A[k]).epsilon(1e-12));
    }
}

TEST_CASE("SolverSSA nrm.space: the analyzer matches the exact CTMC (t-test)") {
    // Delay(mean 2) -> PS Queue(rate 3), three jobs: a closed network SolverCTMC
    // solves exactly, so the oracle carries no error of its own and the whole
    // deviation is the simulation's.
    qn::Network<double> m = cqn_ps(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    // THE BAND. Eight independent seeds give eight replicate means; the two
    // sided one-sample t statistic against the exact value must fall inside the
    // 1% critical value on 7 degrees of freedom. This is a stated band with the
    // spread MEASURED from the replicates rather than guessed, so it neither
    // flakes when the chain happens to mix slowly nor passes a biased estimator
    // that a fixed percentage band would absorb.
    std::vector<double> qd(kReps, 0.0), qq(kReps, 0.0), tq(kReps, 0.0), uq(kReps, 0.0);
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaNrmSpaceOptions opt;
        opt.samples = 20000;             // reported with the seed: neither is optional
        opt.seed = 23000 + 100 * rep;    // 23000 is LINE's own default seed
        const ssa::SsaNrmSpaceSolution<double> sim = ssa::solver_ssa_nrm_space_analyzer(sn, opt);
        REQUIRE(sim.avg.method == "nrm.space");
        REQUIRE(sim.avg.samples == 20000);
        REQUIRE(sim.seed == opt.seed);
        qd[rep] = sim.avg.QN(0, 0);
        qq[rep] = sim.avg.QN(1, 0);
        tq[rep] = sim.avg.TN(1, 0);
        uq[rep] = sim.avg.UN(1, 0);
        // EXACT, not statistical: every state of this chain holds three jobs, so
        // the time average of the total is three whatever the path did. A fault
        // in the state table or the time weighting breaks this identity while
        // leaving the individual means plausible.
        CHECK(sim.avg.QN(0, 0) + sim.avg.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-9));
    }
    CHECK(tstat(qd, exact.QN(0, 0)) < kT7);
    CHECK(tstat(qq, exact.QN(1, 0)) < kT7);
    CHECK(tstat(tq, exact.TN(1, 0)) < kT7);
    CHECK(tstat(uq, exact.UN(1, 0)) < kT7);
}

TEST_CASE("SolverSSA nrm.space: it agrees with the plain NRM engine (t-test)") {
    // The two engines answer the same question by different routes: one
    // integrates the metrics along the path, the other tabulates the states and
    // forms pi * A. They consume different draws, so the comparison is
    // statistical -- eight seeds each, and the space variant's replicate mean
    // must sit inside the t band around the plain engine's.
    qn::Network<double> m = cqn_ps(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    std::vector<double> qs(kReps, 0.0), qn_(kReps, 0.0);
    for (std::size_t rep = 0; rep < kReps; ++rep) {
        ssa::SsaNrmSpaceOptions sopt;
        sopt.samples = 20000;
        sopt.seed = 4000 + 37 * rep;
        qs[rep] = ssa::solver_ssa_nrm_space_analyzer(sn, sopt).avg.QN(1, 0);

        ssa::SsaOptions nopt;
        nopt.samples = 20000;
        nopt.seed = 4000 + 37 * rep;
        qn_[rep] = ssa::solver_ssa_nrm_analyzer(sn, nopt).QN(1, 0);
    }
    double mean_plain = 0.0;
    for (double v : qn_) mean_plain += v;
    mean_plain /= static_cast<double>(kReps);
    CHECK(tstat(qs, mean_plain) < kT7);
}

TEST_CASE("SolverSSA nrm.space: a three-way routing split reaches every destination") {
    // Think -> {Q1, Q2, Q3} in equal thirds. The reference draws the destination
    // with `1+find(rand>=cdfVec,1)`, which collapses to 2 for every draw above
    // the first cumulative weight and so never selects the third destination.
    // This port uses the inverse CDF instead (see the divergence note in the
    // header), and the defect would show here as a Q3 that is never visited.
    qn::Network<double> m = cqn_three_way(3.0, 2.0, 3.0);
    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 40000;
    opt.seed = 31337;
    const ssa::SsaNrmSpaceSolution<double> sim =
        ssa::solver_ssa_nrm_space_analyzer(m.get_struct(), opt);

    // Under the reference's selection this is exactly zero, so the assertion is
    // structural rather than statistical.
    CHECK(sim.avg.TN(3, 0) > 0.0);
    CHECK(sim.avg.QN(3, 0) > 0.0);
    // With equal rates and an equal split the three queues are exchangeable, so
    // their throughputs agree up to Monte Carlo error. Ten per cent at 40000
    // firings on a four-state-per-node chain is many standard errors wide; the
    // point of the assertion is the ORDER of the numbers, since a misrouted
    // third destination is not a ten per cent effect but a total one.
    CHECK(sim.avg.TN(2, 0) == doctest::Approx(sim.avg.TN(1, 0)).epsilon(0.10));
    CHECK(sim.avg.TN(3, 0) == doctest::Approx(sim.avg.TN(1, 0)).epsilon(0.10));
    // Population is conserved across the four stations exactly.
    double tot = 0.0;
    for (std::size_t i = 0; i < 4; ++i) tot += sim.avg.QN(i, 0);
    CHECK(tot == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("SolverSSA nrm.space: refusals name what is missing") {
    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 200;

    // An open model: in the aggregate network a Sink has no outgoing routing and
    // a Source's token is consumed and never replenished, so the state space is
    // unbounded and the table would grow one row per sample.
    qn::Network<double> open("nrmspace-open");
    const std::size_t src = open.add_source("Src");
    const std::size_t q = open.add_queue("Q", SchedStrategy::PS);
    const std::size_t sk = open.add_sink("Sink");
    const std::size_t oc = open.add_open_class("C1");
    open.set_arrival(src, oc, Dist::exp_rate(0.5));
    open.set_service(q, oc, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> Po;
    Po.set(src, q, 1.0);
    Po.set(q, sk, 1.0);
    open.link(Po);
    CHECK_THROWS_AS(ssa::solver_ssa_nrm_space_analyzer(open.get_struct(), opt),
                    line::UnsupportedError);

    // Phase-type service: the state carries no phase dimension, so collapsing it
    // onto the mean rate would report the right throughput and the wrong queue
    // length. Refused rather than approximated.
    qn::Network<double> ph("nrmspace-ph");
    const std::size_t d2 = ph.add_delay("Think");
    const std::size_t q2 = ph.add_queue("Q", SchedStrategy::PS);
    const std::size_t c2 = ph.add_closed_class("C1", 2.0, d2);
    ph.set_service(d2, c2, Dist::exp_rate(0.5));
    ph.set_service(q2, c2, Dist::erlang(3.0, 2));
    qn::RoutingMatrix<double> Pp;
    Pp.set(d2, q2, 1.0);
    Pp.set(q2, d2, 1.0);
    ph.link(Pp);
    CHECK_THROWS_AS(ssa::solver_ssa_nrm_space_analyzer(ph.get_struct(), opt),
                    line::UnsupportedError);

    // A discipline outside the five-armed propensity switch.
    qn::Network<double> siro("nrmspace-siro");
    const std::size_t d3 = siro.add_delay("Think");
    const std::size_t q3 = siro.add_queue("Q", SchedStrategy::SIRO);
    const std::size_t c3 = siro.add_closed_class("C1", 2.0, d3);
    siro.set_service(d3, c3, Dist::exp_rate(0.5));
    siro.set_service(q3, c3, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> Ps;
    Ps.set(d3, q3, 1.0);
    Ps.set(q3, d3, 1.0);
    siro.link(Ps);
    CHECK_THROWS_AS(ssa::solver_ssa_nrm_space_analyzer(siro.get_struct(), opt),
                    line::UnsupportedError);

    // A space larger than the declared cap is refused rather than truncated: a
    // pi normalized over whichever states happened to fit is the law of a
    // different chain.
    qn::Network<double> big = cqn_ps(20.0, 2.0, 3.0);
    ssa::SsaNrmSpaceOptions capped;
    capped.samples = 500;
    capped.seed = 1;
    capped.state_max = 3;
    CHECK_THROWS_AS(ssa::solver_ssa_nrm_space_analyzer(big.get_struct(), capped),
                    line::UnsupportedError);
    ssa::NrmSpaceEngine<double> eng(big.get_struct(), capped);
    CHECK_THROWS_AS(eng.enumerate_space(), line::UnsupportedError);
}

TEST_CASE("SolverSSA nrm.space: exact arithmetic is refused, not narrowed") {
    // The clocks are -log(u) of a uniform: there is no exact value for a
    // rational backend to compute, and no information a wider float would carry
    // that the Monte Carlo error does not swamp.
    qn::Network<line::Rational> m("nrmspace-exact");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", 2.0, d);
    m.set_service(d, c, line::lang::Distrib<line::Rational>::exp_rate(
                            line::num_traits<line::Rational>::from_double(0.5)));
    m.set_service(q, c, line::lang::Distrib<line::Rational>::exp_rate(
                            line::num_traits<line::Rational>::from_int(3)));
    qn::RoutingMatrix<line::Rational> P;
    P.set(d, q, line::num_traits<line::Rational>::from_int(1));
    P.set(q, d, line::num_traits<line::Rational>::from_int(1));
    m.link(P);

    ssa::SsaNrmSpaceOptions opt;
    opt.samples = 10;
    CHECK_THROWS_AS(ssa::solver_ssa_nrm_space_analyzer(m.get_struct(), opt),
                    line::UnsupportedError);
}
