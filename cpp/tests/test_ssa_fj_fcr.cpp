/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The two constructs the SSA serial engine used to refuse: fork-join and finite
 * capacity regions.
 *
 * BOTH ARE MEASURED AGAINST THE EXACT CHAIN THAT MODELS THE SAME RULE, and the
 * pairing is the point. A fork-join model is simulated on the tag-augmented copy
 * and compared with SolverCTMC, which solves the augmented copy too, so a
 * disagreement is in the sample path and not in the transformation. A DROP
 * region is compared with the CENSORED generator and a WAITQ region with the
 * token-FIFO generator, because those are the two different models the rules
 * name; comparing either against the other would fail for a reason that has
 * nothing to do with simulation.
 *
 * WHERE AN EXACT IDENTITY EXISTS IT IS ASSERTED FIRST, because no amount of
 * Monte Carlo error can perturb one and it detects a bias the statistical bands
 * cannot. A DROP region on a closed model conserves the population and caps
 * every visited state; a WAITQ model conserves it only once the parked jobs are
 * added back. A fork-join model has NO such identity -- one parent becomes B
 * siblings between the fork and the join, so the summed queue lengths exceed the
 * declared population by the siblings in flight -- and the test says so instead
 * of asserting the closed-network identity that holds everywhere else.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/ssa/solver_ssa_serial.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace ssa = line::ssa;
using line::lang::DropStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> Fork -> {Q1, Q2} -> Join -> Delay, one closed class of `n` jobs. */
qn::Network<double> fj_cqn(double n) {
    qn::Network<double> m("ssa-fj-cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t f = m.add_fork("Fork");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, f, 1.0);
    P.set(f, q1, 1.0);
    P.set(f, q2, 1.0);
    P.set(q1, j, 1.0);
    P.set(q2, j, 1.0);
    P.set(j, d, 1.0);
    m.link(P);
    return m;
}

/** Think -> Q1 -> Q2 -> Think, one closed class; nodes 1, 2, 3 in that order. */
qn::Network<double> cycle3(double njobs) {
    qn::Network<double> m("ssa-fcr-cycle3");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q1, c, Dist::exp_rate(2.0));
    m.set_service(q2, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q1, 1.0);
    P.set(q1, q2, 1.0);
    P.set(q2, d, 1.0);
    m.link(P);
    return m;
}

double sum_all(const line::Matrix<double>& M) {
    double s = 0;
    for (std::size_t i = 0; i < M.rows(); ++i)
        for (std::size_t r = 0; r < M.cols(); ++r) s += M(i, r);
    return s;
}

double sum_all(const std::vector<double>& v) {
    double s = 0;
    for (std::size_t i = 0; i < v.size(); ++i) s += v[i];
    return s;
}

}  // namespace

TEST_CASE("SolverSSA serial: a fork-join model fires the fork and folds the siblings back") {
    const double N = 2.0;
    qn::Network<double> m = fj_cqn(N);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;  // reported with the seed: neither number is optional
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

    // THE TABLE IS IN THE CALLER'S CLASS SPACE, the path in the augmented one.
    REQUIRE(sim.avg.QN.cols() == sn.nclasses);
    REQUIRE(sim.fjclassmap.size() > sn.nclasses);
    REQUIRE(sim.avg.samples == 200000);
    REQUIRE(sim.seed == 23000);

    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    // THE POPULATION IS NOT CONSERVED IN A FORK-JOIN MODEL and the test says so
    // rather than asserting the closed-network identity that holds everywhere
    // else: one parent task becomes B siblings between the fork and the join, so
    // the summed queue lengths exceed N by the mean number of siblings in
    // flight. What must hold is agreement with the chain that models the same
    // augmentation, which is the comparison below.
    CHECK(sum_all(sim.avg.QN) > N);
    // THE BAND: at 2e5 firings on a chain this small the relative standard error
    // of a time-averaged mean is well under one per cent, so 10 per cent is
    // roughly ten standard errors and a violation is a bias rather than an
    // unlucky path.
    const double band = 0.10;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        INFO("station ", sn.stations[i].name);
        CHECK(sim.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(band));
        CHECK(sim.avg.TN(i, 0) == doctest::Approx(exact.TN(i, 0)).epsilon(band));
    }
    // The branches are genuinely loaded: a fork that never fired would leave
    // both queues empty and the population identity above would still hold.
    double branch = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == "Q1" || sn.stations[i].name == "Q2")
            branch += sim.avg.QN(i, 0);
    CHECK(branch > 1e-3);
}

TEST_CASE("SolverSSA: the dispatcher answers a fork-join model instead of refusing it") {
    // `default` is the NRM, whose eligibility gate excludes fork-join in every
    // codebase; the reference falls back to the serial engine there and so does
    // this port. The method the ANSWER carries must say which engine ran.
    qn::Network<double> m = fj_cqn(2.0);
    ssa::SsaOptions opt;
    opt.samples = 20000;
    opt.seed = 23000;
    const ssa::SsaSolution r = ssa::solver_ssa(m.get_struct(), opt);
    CHECK(r.method == "serial");
    CHECK(r.samples == 20000);
    // The branches carry load, so the fork fired: a dispatcher that answered
    // with an engine which never fires a fork would report them empty.
    CHECK(sum_all(r.QN) > 2.0);
    CHECK(r.XN[0] > 0.0);
}

TEST_CASE("SolverSSA serial: a DROP region censors the path where SolverCTMC censors the chain") {
    const double N = 3.0, cap = 1.0;
    const std::size_t Q1 = 2, Q2 = 3;
    qn::Network<double> m = cycle3(N);
    m.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                 std::vector<DropStrategy>{DropStrategy::DROP});
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

    // EXACT: under censoring no job is destroyed on a closed model, so the
    // population is intact, and NO state the path visits may put more than `cap`
    // jobs inside the region. The second is the check that the censoring
    // actually happened: a run that ignored the region satisfies the first.
    CHECK(sum_all(sim.avg.QN) == doctest::Approx(N).epsilon(1e-9));
    const std::size_t K = sn.nclasses;
    const std::size_t i1 = sn.nodes[Q1 - 1].station, i2 = sn.nodes[Q2 - 1].station;
    bool capped = true;
    for (std::size_t s = 0; s < sim.run.ssq.rows(); ++s)
        if (sim.run.ssq(s, (i1 - 1) * K) + sim.run.ssq(s, (i2 - 1) * K) > cap + 1e-9) capped = false;
    CHECK(capped);
    CHECK(sum_all(sim.parked) == doctest::Approx(0.0).epsilon(1e-12));

    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());
    const double band = 0.10;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        INFO("station ", sn.stations[i].name);
        CHECK(sim.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(band));
        CHECK(sim.avg.TN(i, 0) == doctest::Approx(exact.TN(i, 0)).epsilon(band));
    }
}

TEST_CASE("SolverSSA serial: a WAITQ region parks jobs and the population balances") {
    const double N = 3.0, cap = 1.0;
    const std::size_t Q1 = 2, Q2 = 3;
    qn::Network<double> m = cycle3(N);
    m.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                 std::vector<DropStrategy>{DropStrategy::WAITQ});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(ctmc::ctmc_has_waitq_region(sn));

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim = ssa::solver_ssa_serial_analyzer(sn, opt);

    // EXACT, per visited state: a parked job is in no station, so the station
    // counts plus the FIFO contents must account for every job in every state.
    // A cascade that duplicated a job on release, or lost one on parking, breaks
    // this at the first such state.
    REQUIRE(sim.run.buf.size() == sim.run.space.size());
    const std::size_t K = sn.nclasses;
    const std::size_t i1 = sn.nodes[Q1 - 1].station, i2 = sn.nodes[Q2 - 1].station;
    bool balanced = true, capped = true, any_parked = false;
    for (std::size_t s = 0; s < sim.run.ssq.rows(); ++s) {
        double at_stations = 0;
        for (std::size_t c = 0; c < sim.run.ssq.cols(); ++c) at_stations += sim.run.ssq(s, c);
        double park = 0;
        for (std::size_t f = 0; f < sim.run.buf[s].size(); ++f) park += sim.run.buf[s][f].size();
        if (park > 0) any_parked = true;
        if (std::fabs(at_stations + park - N) > 1e-9) balanced = false;
        if (sim.run.ssq(s, (i1 - 1) * K) + sim.run.ssq(s, (i2 - 1) * K) > cap + 1e-9)
            capped = false;
    }
    CHECK(balanced);
    CHECK(capped);
    CHECK(any_parked);  // the region is binding, so the identities are not vacuous

    // The same balance on the means, which is what a caller reads.
    CHECK(sum_all(sim.avg.QN) + sum_all(sim.parked) == doctest::Approx(N).epsilon(1e-9));
    CHECK(sum_all(sim.avg.QN) < N - 1e-6);

    // Against the exact WAITQ chain, which is a DIFFERENT model from the DROP
    // one: the simulator must reproduce the rule it was given.
    const ctmc::WaitqSolution<double> exact =
        ctmc::solver_ctmc_waitq_analyzer(sn, ctmc::CtmcOptions());
    const double band = 0.10;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        INFO("station ", sn.stations[i].name);
        CHECK(sim.avg.QN(i, 0) == doctest::Approx(exact.sol.avg.QN(i, 0)).epsilon(band));
        CHECK(sim.avg.TN(i, 0) == doctest::Approx(exact.sol.avg.TN(i, 0)).epsilon(band));
    }
    CHECK(sum_all(sim.parked) == doctest::Approx(sum_all(exact.parked)).epsilon(band));
}

TEST_CASE("SolverSSA serial: WAITQ and DROP disagree on the same capped model") {
    // The discriminator: two models differing only in the rule. If the engine
    // served a WAITQ region as DROP -- the silent failure this machinery exists
    // to prevent -- the two runs would agree to within Monte Carlo error.
    const double N = 3.0, cap = 1.0;
    const std::size_t Q1 = 2, Q2 = 3;
    qn::Network<double> mdrop = cycle3(N);
    mdrop.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                     std::vector<DropStrategy>{DropStrategy::DROP});
    qn::Network<double> mwait = cycle3(N);
    mwait.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                     std::vector<DropStrategy>{DropStrategy::WAITQ});

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sd = ssa::solver_ssa_serial_analyzer(mdrop.get_struct(), opt);
    const ssa::SsaSerialSolution<double> sw = ssa::solver_ssa_serial_analyzer(mwait.get_struct(), opt);

    CHECK(sum_all(sd.avg.QN) == doctest::Approx(N).epsilon(1e-9));
    CHECK(sum_all(sw.avg.QN) < sum_all(sd.avg.QN) - 1e-3);
    CHECK(sum_all(sw.parked) > 1e-3);
}

/*
 * The remaining constructs the DISPATCHER now reaches through the serial engine.
 *
 * Each one is a construct the NRM gate refuses by name, so `default` can only
 * answer it by falling back -- which is what these assert, together with an
 * invariant sharp enough to catch a fallback that ran but simulated the wrong
 * model. They are the evidence behind `ssa_feature_set` declaring them: a
 * featset entry is a CLAIM, and AUTO chooses SolverSSA on the strength of it.
 */

TEST_CASE("SolverSSA: a cache model reaches the serial engine and its reads split") {
    // Source -> Cache(2 items, 1 slot, LRU) -> Sink, one read class switching
    // into Hit or Miss. Every read leaves as exactly one of the two, so the
    // realized hit ratio is a probability and the two departure streams divide
    // the read rate between them.
    const double lambda = 0.5;
    qn::Network<double> m("ssa-cache");
    const std::size_t src = m.add_source("Source");
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t mis = m.add_open_class("Miss");

    qn::CacheParam<double> cp;
    cp.nitems = 2;
    cp.itemcap.push_back(1);
    cp.replacestrat = line::lang::ReplacementStrategy::LRU;
    cp.pread.assign(3, std::vector<double>());
    cp.pread[rd - 1] = std::vector<double>{0.5, 0.5};
    cp.hitclass.assign(3, 0);
    cp.missclass.assign(3, 0);
    cp.hitclass[rd - 1] = hit;
    cp.missclass[rd - 1] = mis;
    const std::size_t ca = m.add_cache("C", cp);

    m.set_arrival(src, rd, Dist::exp_rate(lambda));
    qn::RoutingMatrix<double> P;
    P.set(rd, rd, src, ca, 1.0);
    P.set(hit, hit, ca, sk, 1.0);
    P.set(mis, mis, ca, sk, 1.0);
    m.link(P);

    ssa::SsaSerialOptions opt;
    opt.samples = 50000;
    opt.seed = 23000;
    opt.cutoff = 1.0;
    const ssa::SsaSerialSolution<double> sim =
        ssa::solver_ssa_serial_analyzer(m.get_struct(), opt);

    REQUIRE(sim.cache.size() == 1);
    CHECK(sim.cache[0].node == ca);
    const double h = sim.cache[0].hitprob[rd - 1], mi = sim.cache[0].missprob[rd - 1];
    CHECK(h + mi == doctest::Approx(1.0).epsilon(1e-9));
    // With one slot and two equally read items under LRU the request that just
    // arrived is a hit exactly when it repeats the previous one, so the ratio is
    // strictly inside (0,1) and the cache is neither always cold nor always warm.
    CHECK(h > 0.0);
    CHECK(h < 1.0);
}

TEST_CASE("SolverSSA: a stochastic Petri net fires its transitions on the sample path") {
    // P1 -T1-> P2 -T2-> P1 with one token: the two-state marked graph whose
    // stationary law is pi(P1) = r2/(r1+r2). The firing is a GLOBAL
    // synchronization, atomic across both places, which no ordinary sync can
    // express -- so a run that reproduces the balance has fired them properly.
    const double r1 = 2.0, r2 = 3.0;
    qn::Network<double> m("ssa-marked-graph");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, Dist::exp_rate(1.0));
    m.set_service(p2, c, Dist::exp_rate(1.0));
    qn::TransitionParam<double> t1;
    t1.nmodes = 1;
    t1.modenames.push_back("fire");
    t1.enabling.assign(1, line::Matrix<double>(4, 1, 0.0));
    t1.inhibiting.assign(1, line::Matrix<double>(4, 1, std::numeric_limits<double>::infinity()));
    t1.firing.assign(1, line::Matrix<double>(4, 1, 0.0));
    t1.enabling[0](p1 - 1, 0) = 1.0;
    t1.firing[0](p2 - 1, 0) = 1.0;
    t1.nmodeservers.push_back(1.0);
    t1.firingphases.push_back(1);
    t1.timing.push_back(line::lang::TimingStrategy::TIMED);
    t1.fireweight.push_back(1.0);
    t1.firingproc.push_back(Dist::exp_rate(r1));
    qn::TransitionParam<double> t2 = t1;
    t2.enabling[0] = line::Matrix<double>(4, 1, 0.0);
    t2.firing[0] = line::Matrix<double>(4, 1, 0.0);
    t2.enabling[0](p2 - 1, 0) = 1.0;
    t2.firing[0](p1 - 1, 0) = 1.0;
    t2.firingproc[0] = Dist::exp_rate(r2);
    m.add_transition("T1", t1);
    m.add_transition("T2", t2);

    ssa::SsaSerialOptions opt;
    opt.samples = 200000;
    opt.seed = 23000;
    const ssa::SsaSerialSolution<double> sim =
        ssa::solver_ssa_serial_analyzer(m.get_struct(), opt);

    // EXACT: the token is never created or destroyed by a firing.
    CHECK(sum_all(sim.avg.QN) == doctest::Approx(1.0).epsilon(1e-9));
    // The marking balance, inside a Monte Carlo band: P1 holds the token a
    // fraction r2/(r1+r2) of the time.
    CHECK(sim.avg.QN(0, 0) == doctest::Approx(r2 / (r1 + r2)).epsilon(0.05));
}

TEST_CASE("SolverSSA: PAS and the priority shares reach the engine and conserve population") {
    // One closed class over Think -> Q under each discipline the NRM gate
    // refuses. The population identity is exact and holds under any correct
    // sample path; a discipline whose rate law was silently absent would either
    // throw or freeze the queue at zero, and the second check catches that.
    const SchedStrategy scheds[3] = {SchedStrategy::PAS, SchedStrategy::OI,
                                     SchedStrategy::PSPRIO};
    for (std::size_t t = 0; t < 3; ++t) {
        qn::Network<double> m("ssa-sched");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Q", scheds[t]);
        const std::size_t c = m.add_closed_class("C1", 2.0, d);
        m.set_service(d, c, Dist::exp_rate(1.0));
        m.set_service(q, c, Dist::exp_rate(2.0));
        qn::RoutingMatrix<double> P;
        P.set(d, q, 1.0);
        P.set(q, d, 1.0);
        m.link(P);
        // PAS and OI take their rate law as mu(c) over the ORDERED job list
        // rather than from `sn.rates`, so the station has none until it is set:
        // mu(c) = 2|c| here, which is the rate an infinite server would deliver.
        if (scheds[t] == SchedStrategy::PAS || scheds[t] == SchedStrategy::OI)
            m.set_pas(q, [](const std::vector<std::size_t>& lst) {
                return 2.0 * static_cast<double>(lst.size());
            });

        ssa::SsaOptions opt;
        opt.samples = 20000;
        opt.seed = 23000;
        INFO("discipline index ", t);
        const ssa::SsaSolution r = ssa::solver_ssa(m.get_struct(), opt);
        // The NRM gate refuses PAS and OI by name -- their rate laws are genuine
        // sub-engines that are not ported -- and ACCEPTS PSPRIO, which joined the
        // whitelist on 2026-08-15 to match `solver_ssa_analyzer_nrm.m`. Asserting
        // one verdict for all three is what left this red once that landed.
        const bool nrmRefuses =
            scheds[t] == SchedStrategy::PAS || scheds[t] == SchedStrategy::OI;
        CHECK(r.method == (nrmRefuses ? "serial" : "nrm"));
        CHECK(sum_all(r.QN) == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(r.QN(1, 0) > 1e-6);
        CHECK(r.TN(1, 0) > 0.0);
        // PAS and OI with mu(c) = 2|c| are an INFINITE SERVER of rate 2, so the
        // model is Think(mean 1) -> IS(mean 0.5) and the exact mean queue length
        // is N * (0.5 / 1.5) = 2/3. The ORACLE IS ANALYTIC AND NOT SolverCTMC
        // here: the CTMC enumerates a PAS station's states through
        // `from_marginal_node`, which has no ordered-list arm, so its space is
        // in the [buffer | server] encoding while its handlers read it as a
        // list -- it reports 0.5 on this model. That gap is recorded in
        // `_kb/06-solver-catalog.md`; it is upstream of SolverSSA, which walks
        // the encoding the handlers themselves produce.
        if (scheds[t] == SchedStrategy::PAS || scheds[t] == SchedStrategy::OI)
            CHECK(r.QN(1, 0) == doctest::Approx(2.0 / 3.0).epsilon(0.10));
    }
}
