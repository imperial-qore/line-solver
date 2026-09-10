/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The replicated SolverSSA analyzer (`solver_ssa_analyzer_parallel.m`) and the
 * per-replica event memo (`EventCache.m`).
 *
 * THE TWO THINGS THESE TESTS EXIST TO CATCH.
 *
 *   REPLICATION SILENTLY DEGENERATING INTO A LONGER RUN. If the replicas were
 *   pooled into one sample path, or seeded identically, or if the reported
 *   error bar were computed from the total firing count rather than from the
 *   spread between the R replica means, the point estimate would still look
 *   right and the error bar would be too small. So the tests assert the
 *   combination rule ITSELF: `avg` is the arithmetic mean of the R replica
 *   tables and `*_sem` is their sample standard deviation over sqrt(R), both
 *   recomputed here from `replica` and compared exactly, and the R replicas are
 *   checked to be genuinely different streams.
 *
 *   A CACHE HIT DIFFERING FROM RECOMPUTATION. That is the memo's entire
 *   correctness condition, so it is tested directly and exhaustively: over
 *   every reachable state, every stateful node, every class, every event and
 *   both settings of the no-promote flag, the cached answer is compared field
 *   by field against the free `qn::after_event`. The miss count is checked too,
 *   against the number of distinct queries counted independently in the test,
 *   so a key that MERGES two different queries fails even when the two happen
 *   to have the same value at the states visited.
 *
 * Every simulated number below is reported with its replica count, its
 * per-replica sample count and its base seed, because a simulated number
 * without them is not a measurement. Nothing simulated is compared against an
 * exact answer at a tight tolerance; the one such comparison uses a band stated
 * and justified at the assertion.
 */
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/state_events.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ssa/solver_ssa_parallel.h"
#include "line/solvers/ssa/ssa_event_cache.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace ssa = line::ssa;
using line::lang::EventType;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Think -> FCFS Queue -> Think, one closed class of `n` jobs. */
qn::Network<double> par_cqn(double n, double think_mean, double mu) {
    qn::Network<double> m("ssa-parallel-cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, Dist::exp_rate(1.0 / think_mean));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** True when two event outcomes agree in every field, bit for bit. */
bool same_outcome(const qn::EventOutcome<double>& a, const qn::EventOutcome<double>& b) {
    return a.space == b.space && a.rate == b.rate && a.prob == b.prob;
}

/** The three events a queueing station can be asked about here. */
const EventType kEvents[3] = {EventType::ARV, EventType::DEP, EventType::PHASE};

}  // namespace

TEST_CASE("SSA event cache: a hit returns exactly what recomputation returns") {
    qn::Network<double> m = par_cqn(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ssa::SsaReachability<double> reach = ssa::solver_ssa_reachability(sn);
    REQUIRE(reach.space.size() > 1);

    ssa::SsaEventCache<double> cache = ssa::SsaEventCache<double>::create(true, sn);
    REQUIRE(cache.enabled());

    // The distinct queries, counted with a key built HERE and not by the header,
    // so the two must agree for the miss count below to come out right.
    std::set<std::tuple<std::size_t, int, std::size_t, bool, std::vector<double> > > distinct;
    std::size_t queries = 0;
    bool all_match = true;

    for (std::size_t s = 0; s < reach.space.size(); ++s)
        for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f) {
            const std::size_t ind = sn.stateful_nodes[f];
            const std::vector<double>& row = reach.space[s].local[f];
            for (std::size_t cls = 1; cls <= sn.nclasses; ++cls)
                for (std::size_t e = 0; e < 3; ++e)
                    for (int np = 0; np < 2; ++np) {
                        const bool no_promote = np == 1;
                        const qn::EventOutcome<double> truth =
                            qn::after_event(sn, ind, row, kEvents[e], cls, no_promote);
                        // The FIRST call populates, the SECOND is served from
                        // the memo; both must equal recomputation, which is what
                        // makes a hit unobservable in the value.
                        const qn::EventOutcome<double> first =
                            cache.after_event(sn, ind, row, kEvents[e], cls, no_promote);
                        const qn::EventOutcome<double> second =
                            cache.after_event(sn, ind, row, kEvents[e], cls, no_promote);
                        all_match = all_match && same_outcome(truth, first) &&
                                    same_outcome(truth, second);
                        distinct.insert(std::make_tuple(ind, static_cast<int>(kEvents[e]), cls,
                                                        no_promote, row));
                        queries += 2;
                    }
        }

    CHECK(all_match);
    // Every query was made twice, so the memo must have missed exactly once per
    // DISTINCT query and hit on everything else. An over-merging key would
    // undershoot the miss count; a key carrying something it should not would
    // overshoot it.
    CHECK(cache.misses() == distinct.size());
    CHECK(cache.size() == distinct.size());
    CHECK(cache.hits() == queries - distinct.size());
    CHECK(cache.hits() > 0);
}

TEST_CASE("SSA event cache: the no-promote flag is part of the key") {
    // `noPromote` distinguishes the departure half of an immediate-feedback
    // self-loop from an ordinary departure at the same station in the same
    // state. Whether the two answers differ is a property of the discipline;
    // that they are stored SEPARATELY is a property of the key, and it is the
    // key that is tested here -- so the oracle is the miss count, which holds
    // even where the two values coincide.
    qn::Network<double> m = par_cqn(3.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ssa::SsaReachability<double> reach = ssa::solver_ssa_reachability(sn);
    REQUIRE(!reach.space.empty());

    ssa::SsaEventCache<double> cache = ssa::SsaEventCache<double>::create(true, sn);
    const std::size_t ind = sn.stateful_nodes.back();
    const std::vector<double>& row = reach.space.back().local[sn.stateful_nodes.size() - 1];

    cache.after_event(sn, ind, row, EventType::DEP, 1, false);
    CHECK(cache.misses() == 1);
    cache.after_event(sn, ind, row, EventType::DEP, 1, true);
    CHECK(cache.misses() == 2);
    CHECK(cache.hits() == 0);
    // ... and each is then served, from its own entry.
    CHECK(same_outcome(cache.after_event(sn, ind, row, EventType::DEP, 1, false),
                       qn::after_event(sn, ind, row, EventType::DEP, 1, false)));
    CHECK(same_outcome(cache.after_event(sn, ind, row, EventType::DEP, 1, true),
                       qn::after_event(sn, ind, row, EventType::DEP, 1, true)));
    CHECK(cache.hits() == 2);
    CHECK(cache.size() == 2);
}

TEST_CASE("SSA event cache: a disabled cache computes and stores nothing") {
    // `EventCache.create(false, sn)` returns [] in the reference and afterEvent
    // then takes the uncached path. Here the object survives but memoizes
    // nothing, so the caller keeps ONE code path and the two configurations
    // cannot drift apart.
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ssa::SsaReachability<double> reach = ssa::solver_ssa_reachability(sn);
    ssa::SsaEventCache<double> off = ssa::SsaEventCache<double>::create(false, sn);
    CHECK(!off.enabled());

    const std::size_t ind = sn.stateful_nodes[0];
    const std::vector<double>& row = reach.space[0].local[0];
    for (int i = 0; i < 3; ++i)
        CHECK(same_outcome(off.after_event(sn, ind, row, EventType::DEP, 1),
                           qn::after_event(sn, ind, row, EventType::DEP, 1)));
    CHECK(off.size() == 0);
    CHECK(off.hits() == 0);
    CHECK(off.misses() == 0);
}

TEST_CASE("SSA event cache: a memo is bound to the struct it was built against") {
    // The memoized rows are a function of one struct's rates and capacities, so
    // reusing the memo across structs would answer about another model. The
    // reference's create() takes sn and ignores it; here the binding is checked.
    qn::Network<double> a = par_cqn(2.0, 2.0, 3.0);
    qn::Network<double> b = par_cqn(2.0, 2.0, 7.0);
    const qn::NetworkStruct<double>& sa = a.get_struct();
    const qn::NetworkStruct<double>& sb = b.get_struct();
    const ssa::SsaReachability<double> reach = ssa::solver_ssa_reachability(sa);

    ssa::SsaEventCache<double> cache = ssa::SsaEventCache<double>::create(true, sa);
    const std::vector<double>& row = reach.space[0].local[0];
    CHECK_THROWS_AS(cache.after_event(sb, sb.stateful_nodes[0], row, EventType::DEP, 1),
                    line::InputError);
}

TEST_CASE("SolverSSA parallel: the same base seed reproduces every replica") {
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    ssa::SsaParallelOptions opt;
    opt.samples = 4000;
    opt.seed = 777;
    opt.nreplicas = 4;

    const ssa::SsaParallelSolution<double> r1 = ssa::solver_ssa_parallel_analyzer(m.get_struct(), opt);
    const ssa::SsaParallelSolution<double> r2 = ssa::solver_ssa_parallel_analyzer(m.get_struct(), opt);

    REQUIRE(r1.nreplicas == 4);
    REQUIRE(r1.replica.size() == 4);
    REQUIRE(r1.base_seed == 777);
    // Bit for bit, not to a tolerance: every replica is a pure function of its
    // seed and its budget, so any difference at all is a defect.
    bool identical = true;
    for (std::size_t r = 0; r < 4; ++r) {
        identical = identical && r1.seed[r] == r2.seed[r];
        for (std::size_t i = 0; i < r1.replica[r].QN.rows(); ++i)
            identical = identical && r1.replica[r].QN(i, 0) == r2.replica[r].QN(i, 0);
    }
    CHECK(identical);
    CHECK(r1.avg.QN(1, 0) == r2.avg.QN(1, 0));
    CHECK(r1.QN_sem(1, 0) == r2.QN_sem(1, 0));

    // A different base seed is a different set of R streams.
    ssa::SsaParallelOptions other = opt;
    other.seed = 31337;
    const ssa::SsaParallelSolution<double> r3 =
        ssa::solver_ssa_parallel_analyzer(m.get_struct(), other);
    CHECK(r3.base_seed == 31337);
    CHECK(r3.avg.QN(1, 0) != r1.avg.QN(1, 0));
}

TEST_CASE("SolverSSA parallel: the budget is split across replicas, not multiplied") {
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    ssa::SsaParallelOptions opt;
    opt.samples = 4001;  // deliberately not divisible by R
    opt.seed = 4242;
    opt.nreplicas = 4;
    const ssa::SsaParallelSolution<double> d = ssa::solver_ssa_parallel_analyzer(m.get_struct(), opt);

    // ceil(4001/4) = 1001 firings EACH, so 4004 in total: the request is
    // rounded up, never down, and the total is the budget rather than R times
    // the budget. Reading `avg.samples` as a precision is the mistake the SEM
    // fields exist to prevent, and it is 4004 here for four estimates.
    CHECK(d.samples_requested == 4001);
    CHECK(d.samples_per_replica == 1001);
    CHECK(d.avg.samples == 4004);
    CHECK(d.avg.method == "parallel");
    bool budgets_ok = true, seeds_ok = true;
    for (std::size_t r = 0; r < 4; ++r) {
        budgets_ok = budgets_ok && d.replica[r].samples == 1001;
        seeds_ok = seeds_ok && d.seed[r] == 4242 + r;
        budgets_ok = budgets_ok && d.replica[r].method == "serial";
    }
    CHECK(budgets_ok);
    CHECK(seeds_ok);
}

TEST_CASE("SolverSSA parallel: replica r is fixed by base_seed+r and its own budget") {
    // The reference's worker-count invariance, in the form that survives having
    // no workers: replica r depends on (base_seed + r, per-replica budget) and
    // on nothing else, so two calls that give replica r the same pair must
    // produce the same replica -- even when the two calls run a different
    // NUMBER of replicas. This is what makes the answer a function of
    // (seed, samples, R) alone, which the earlier spmd implementation was not.
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    ssa::SsaParallelOptions few;
    few.samples = 4000;
    few.seed = 900;
    few.nreplicas = 4;  // 1000 firings each
    ssa::SsaParallelOptions many = few;
    many.samples = 8000;
    many.nreplicas = 8;  // also 1000 firings each

    const ssa::SsaParallelSolution<double> a = ssa::solver_ssa_parallel_analyzer(m.get_struct(), few);
    const ssa::SsaParallelSolution<double> b = ssa::solver_ssa_parallel_analyzer(m.get_struct(), many);
    REQUIRE(a.samples_per_replica == 1000);
    REQUIRE(b.samples_per_replica == 1000);

    bool shared_prefix = true;
    for (std::size_t r = 0; r < 4; ++r) {
        shared_prefix = shared_prefix && a.seed[r] == b.seed[r];
        for (std::size_t i = 0; i < a.replica[r].QN.rows(); ++i)
            shared_prefix = shared_prefix && a.replica[r].QN(i, 0) == b.replica[r].QN(i, 0);
    }
    CHECK(shared_prefix);
    // The MEANS differ, because b averages four more replicas -- which is the
    // point: more replicas change the estimate's precision, not what replica r
    // is.
    CHECK(a.avg.QN(1, 0) != b.avg.QN(1, 0));
}

TEST_CASE("SolverSSA parallel: the estimate is a replica mean carrying its own spread") {
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    ssa::SsaParallelOptions opt;
    opt.samples = 16000;
    opt.seed = 20260728;
    opt.nreplicas = 8;
    const ssa::SsaParallelSolution<double> d = ssa::solver_ssa_parallel_analyzer(m.get_struct(), opt);
    REQUIRE(d.replica.size() == 8);

    // Recomputed here from the R replica tables, which is the definition of the
    // combination rule. A pooled sample path would give a similar mean and an
    // error bar smaller by the path's autocorrelation factor; the exact match
    // below is what rules that out.
    for (std::size_t i = 0; i < d.avg.QN.rows(); ++i) {
        double mu = 0.0;
        for (std::size_t r = 0; r < 8; ++r) mu += d.replica[r].QN(i, 0);
        mu /= 8.0;
        double ss = 0.0;
        for (std::size_t r = 0; r < 8; ++r)
            ss += (d.replica[r].QN(i, 0) - mu) * (d.replica[r].QN(i, 0) - mu);
        const double sem = std::sqrt(ss / 7.0) / std::sqrt(8.0);
        CHECK(d.avg.QN(i, 0) == doctest::Approx(mu).epsilon(1e-13));
        CHECK(d.QN_sem(i, 0) == doctest::Approx(sem).epsilon(1e-13));
    }

    // The R replicas must be R DIFFERENT streams. If they were all seeded alike
    // the spread would be exactly zero and the reported SEM would claim an
    // exact answer from 8 copies of one run.
    CHECK(d.QN_sem(1, 0) > 0.0);
    CHECK(d.TN_sem(1, 0) > 0.0);
    CHECK(d.XN_sem[0] > 0.0);

    // A single replica has no observable spread, and NaN is the honest report:
    // zero would assert the estimate is exact.
    ssa::SsaParallelOptions one = opt;
    one.nreplicas = 1;
    const ssa::SsaParallelSolution<double> s = ssa::solver_ssa_parallel_analyzer(m.get_struct(), one);
    CHECK(s.replica.size() == 1);
    CHECK(s.samples_per_replica == 16000);
    CHECK(std::isnan(s.QN_sem(1, 0)));
    CHECK(s.avg.QN(1, 0) == s.replica[0].QN(1, 0));
}

TEST_CASE("SolverSSA parallel: the replica mean sits inside a stated band of the exact chain") {
    // 16 replicas x 8000 firings from base seed 5150, with the leading 10% of
    // each path discarded. The discard matters: the path starts with every job
    // at the reference station, and without it the O(1/n) transient bias would
    // be the same size as the Monte Carlo error at this budget, so a failure
    // would not distinguish the two.
    const double n = 2.0;
    qn::Network<double> m = par_cqn(n, 2.0, 3.0);
    ssa::SsaParallelOptions opt;
    opt.nreplicas = 16;
    opt.samples = 16 * 8000;
    opt.seed = 5150;
    opt.warmupfrac = 0.1;
    const ssa::SsaParallelSolution<double> d = ssa::solver_ssa_parallel_analyzer(m.get_struct(), opt);
    REQUIRE(d.samples_per_replica == 8000);

    ctmc::CtmcOptions copt;
    const line::mva::AvgResult<double> ex = ctmc::solver_ctmc_run_analyzer(m.get_struct(), copt);

    // EXACT AND STRUCTURAL FIRST. Every sampled state of a closed network has
    // its queue lengths summing to the population, so every replica's time
    // average does, and so does their mean. No amount of Monte Carlo error can
    // perturb this, which makes it the sharpest bias detector available.
    double total = 0.0;
    for (std::size_t i = 0; i < d.avg.QN.rows(); ++i) total += d.avg.QN(i, 0);
    CHECK(total == doctest::Approx(n).epsilon(1e-12));

    // THEN THE STATISTICAL BAND. Six standard errors of the replica mean, with
    // the standard error ESTIMATED from the same 16 replicas (15 degrees of
    // freedom). The two-sided t tail beyond 6 sigma at 15 dof is below 3e-5, so
    // a failure here is a bias and not a rare draw; the run is seeded, so the
    // outcome is in fact deterministic and the band only has to be wide enough
    // to have been safe before the seed was chosen. The 1e-9 floor keeps a
    // metric whose replicas happen to agree exactly from demanding an exact
    // match against the chain.
    for (std::size_t i = 0; i < d.avg.QN.rows(); ++i) {
        const double band = 6.0 * d.QN_sem(i, 0) + 1e-9;
        CHECK(std::abs(d.avg.QN(i, 0) - ex.QN(i, 0)) < band);
    }
    const double tband = 6.0 * d.TN_sem(1, 0) + 1e-9;
    CHECK(std::abs(d.avg.TN(1, 0) - ex.TN(1, 0)) < tband);
    const double xband = 6.0 * d.XN_sem[0] + 1e-9;
    CHECK(std::abs(d.avg.XN[0] - ex.XN[0]) < xband);
}

TEST_CASE("SolverSSA parallel: the entry refuses by name what it does not replicate") {
    qn::Network<double> m = par_cqn(2.0, 2.0, 3.0);
    ssa::SsaParallelOptions opt;
    opt.samples = 200;
    opt.nreplicas = 2;

    opt.method = "parallel";
    CHECK_NOTHROW(ssa::solver_ssa_parallel(m.get_struct(), opt));
    opt.method = "para";
    CHECK_NOTHROW(ssa::solver_ssa_parallel(m.get_struct(), opt));

    // One replica is the same estimator at a sqrt(R)-wider error bar, so the
    // entry says so rather than answering the wrong question quietly.
    opt.method = "serial";
    CHECK_THROWS_AS(ssa::solver_ssa_parallel(m.get_struct(), opt), line::UnsupportedError);
    opt.method = "nrm";
    CHECK_THROWS_AS(ssa::solver_ssa_parallel(m.get_struct(), opt), line::UnsupportedError);

    // The memo is ported and verified above, but no engine consults it yet, so
    // asking for it is refused rather than accepted and ignored.
    opt.method = "parallel";
    opt.eventcache = true;
    CHECK_THROWS_AS(ssa::solver_ssa_parallel(m.get_struct(), opt), line::UnsupportedError);
}
