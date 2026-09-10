/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The WAITQ finite capacity region: the augmented-state generator and the
 * feature gates around it.
 *
 * WHY THESE ORACLES. WAITQ has no closed form and no reference table here, so
 * every check is an invariant the construction must satisfy however it computes:
 *
 *   - INERTNESS. A region whose caps bound nothing refuses nothing, so no token
 *     is ever parked and the WAITQ generator must reproduce the region-less
 *     solution. It is the one check that pins the whole walk against a path that
 *     shares no line of code with it: the default analyzer enumerates a lattice
 *     and censors it by connectivity, this walks the reachable space.
 *   - CONSERVATION IS THE DISCRIMINATOR BETWEEN THE TWO RULES. Under DROP a
 *     refused job is destroyed and the chain is censored, so the station queue
 *     lengths still sum to the population. Under WAITQ the job is PARKED, so
 *     they sum to strictly less and only station jobs plus parked jobs recover
 *     the population. That difference is exactly what makes WAITQ a different
 *     mechanism rather than a variant of DROP, and it is checked per STATE and
 *     not only on the means, which no accidental cancellation could satisfy.
 *   - EVERY GATE MUST STOP THE SOLVE. A combination served as though it were
 *     supported returns plausible numbers for a model nobody described, so the
 *     refusal is the behaviour and its absence is the defect.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_fcr.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::DropStrategy;
using line::lang::RoutingStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/**
 * Think -> Q1 -> Q2 -> Think, one closed class; nodes 1, 2, 3 in that order.
 *
 * Two queues so a region can cap a SUM: no per-station capacity expresses "at
 * most one job between these two", and it is the sum that makes the refusal
 * depend on state the refused job's own destination does not carry.
 */
qn::Network<double> waitq_cycle3(double njobs) {
    qn::Network<double> m("waitq_cycle3");
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

/** Total jobs held at stations, over every station and class. */
double station_jobs(const line::Matrix<double>& QN) {
    double s = 0;
    for (std::size_t i = 0; i < QN.rows(); ++i)
        for (std::size_t r = 0; r < QN.cols(); ++r) s += QN(i, r);
    return s;
}

/** Total mean parked jobs, over every class. */
double parked_jobs(const std::vector<double>& parked) {
    double s = 0;
    for (std::size_t r = 0; r < parked.size(); ++r) s += parked[r];
    return s;
}

}  // namespace

TEST_CASE("ctmc waitq: an unbounded WAITQ region reproduces the region-less solution") {
    // The caps are all the -1 sentinel, so nothing is ever refused and no token
    // is ever parked. The augmented generator must then agree with the default
    // one on every mean -- and it reaches that agreement by a completely
    // different route, walking the reachable space where the analyzer enumerates
    // the lattice and keeps a connected component of it.
    const double N = 3.0;
    const std::size_t THINK = 1, Q1 = 2, Q2 = 3;

    qn::Network<double> plain = waitq_cycle3(N);
    qn::Network<double> regd = waitq_cycle3(N);
    regd.add_region(std::vector<std::size_t>{THINK, Q1, Q2}, std::vector<double>{-1.0}, -1.0,
                    std::vector<DropStrategy>{DropStrategy::WAITQ});

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(plain.get_struct(), opt);
    const qn::NetworkStruct<double>& sn1 = regd.get_struct();
    REQUIRE(ctmc::ctmc_has_waitq_region(sn1));
    const ctmc::WaitqSolution<double> d1 = ctmc::solver_ctmc_waitq_analyzer(sn1, opt);

    for (std::size_t i = 0; i < sn1.nstations; ++i) {
        CHECK(d1.sol.avg.QN(i, 0) == doctest::Approx(d0.avg.QN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.UN(i, 0) == doctest::Approx(d0.avg.UN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.TN(i, 0) == doctest::Approx(d0.avg.TN(i, 0)).epsilon(1e-9));
        CHECK(d1.sol.avg.RN(i, 0) == doctest::Approx(d0.avg.RN(i, 0)).epsilon(1e-9));
    }
    CHECK(d1.sol.avg.XN[0] == doctest::Approx(d0.avg.XN[0]).epsilon(1e-9));

    // Not one token was parked, which is what "the region bounds nothing" means
    // stated on the states rather than on the means. A walk that parked a job it
    // should have admitted would still agree on some aggregate by luck; it could
    // not leave every buffer empty.
    bool all_empty = true;
    for (std::size_t s = 0; s < d1.buf.size(); ++s)
        for (std::size_t f = 0; f < d1.buf[s].size(); ++f)
            if (!d1.buf[s][f].empty()) all_empty = false;
    CHECK(all_empty);
    CHECK(parked_jobs(d1.parked) == doctest::Approx(0.0).epsilon(1e-12));

    double mass = 0;
    for (std::size_t s = 0; s < d1.sol.pi.size(); ++s) mass += d1.sol.pi[s];
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("ctmc waitq: a cap above the population is inert for the same reason") {
    // A cap that no state can reach exercises the admission test on every
    // transition and must still never fire. It separates "the caps are read" from
    // "the caps are compared correctly": an off-by-one in the comparison shows up
    // here and nowhere in the sentinel case above.
    const double N = 2.0;
    const std::size_t Q1 = 2, Q2 = 3;

    qn::Network<double> plain = waitq_cycle3(N);
    qn::Network<double> regd = waitq_cycle3(N);
    regd.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, N,
                    std::vector<DropStrategy>{DropStrategy::WAITQ});

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(plain.get_struct(), opt);
    const qn::NetworkStruct<double>& sn1 = regd.get_struct();
    const ctmc::WaitqSolution<double> d1 = ctmc::solver_ctmc_waitq_analyzer(sn1, opt);

    for (std::size_t i = 0; i < sn1.nstations; ++i)
        CHECK(d1.sol.avg.QN(i, 0) == doctest::Approx(d0.avg.QN(i, 0)).epsilon(1e-9));
    CHECK(parked_jobs(d1.parked) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(station_jobs(d1.sol.avg.QN) == doctest::Approx(N).epsilon(1e-9));
}

TEST_CASE("ctmc waitq: a parked job is in no station and the population still balances") {
    // Three jobs, at most one of them between Q1 and Q2. A job refused entry has
    // ALREADY LEFT Think, so it is in no station at all: the station queue
    // lengths sum to strictly less than three, and only adding the FIFO back
    // recovers the population.
    const double N = 3.0, cap = 1.0;
    const std::size_t Q1 = 2, Q2 = 3;

    qn::Network<double> m = waitq_cycle3(N);
    m.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                 std::vector<DropStrategy>{DropStrategy::WAITQ});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    const ctmc::WaitqSolution<double> d = ctmc::solver_ctmc_waitq_analyzer(sn, opt);

    const std::size_t K = sn.nclasses;
    const std::size_t i1 = sn.nodes[Q1 - 1].station, i2 = sn.nodes[Q2 - 1].station;

    // PER STATE, not only on the means: every augmented state accounts for all
    // three jobs, and none of them puts more than `cap` inside the region. A walk
    // that duplicated a job on release, or lost one on parking, breaks the first;
    // one that admitted a job it should have refused breaks the second.
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, d.sol.chain.space);
    REQUIRE(A.rows() == d.buf.size());
    bool balanced = true, capped = true, any_parked = false;
    for (std::size_t s = 0; s < A.rows(); ++s) {
        double at_stations = 0;
        for (std::size_t c = 0; c < A.cols(); ++c) at_stations += A(s, c);
        double park = 0;
        for (std::size_t f = 0; f < d.buf[s].size(); ++f) park += d.buf[s][f].size();
        if (park > 0) any_parked = true;
        if (std::fabs(at_stations + park - N) > 1e-9) balanced = false;
        if (A(s, (i1 - 1) * K) + A(s, (i2 - 1) * K) > cap + 1e-9) capped = false;
    }
    CHECK(balanced);
    CHECK(capped);
    // The region is genuinely binding, so the invariants above are not being
    // satisfied vacuously by a chain in which nothing ever parks.
    CHECK(any_parked);

    // The same balance on the means, which is what a caller reads.
    CHECK(station_jobs(d.sol.avg.QN) + parked_jobs(d.parked) == doctest::Approx(N).epsilon(1e-9));
    CHECK(station_jobs(d.sol.avg.QN) < N - 1e-6);
    CHECK(parked_jobs(d.parked) > 1e-6);
}

TEST_CASE("ctmc waitq: WAITQ parks where DROP destroys, on the same capped model") {
    // The sharpest available discriminator between the two rules, run on models
    // that differ ONLY in the rule. DROP censors the chain, so its stations still
    // hold every job; WAITQ moves the refused job out of every station into the
    // region FIFO, so they do not. Both must still report a normalized law.
    const double N = 3.0, cap = 1.0;
    const std::size_t Q1 = 2, Q2 = 3;

    qn::Network<double> mdrop = waitq_cycle3(N);
    mdrop.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                     std::vector<DropStrategy>{DropStrategy::DROP});
    qn::Network<double> mwait = waitq_cycle3(N);
    mwait.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                     std::vector<DropStrategy>{DropStrategy::WAITQ});

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> dd = ctmc::solver_ctmc_analyzer(mdrop.get_struct(), opt);
    const ctmc::WaitqSolution<double> dw =
        ctmc::solver_ctmc_waitq_analyzer(mwait.get_struct(), opt);

    // A DROP region on a closed model censors rather than destroys, so its
    // population is intact and no accounting is needed.
    CHECK(station_jobs(dd.avg.QN) == doctest::Approx(N).epsilon(1e-9));
    // WAITQ holds jobs outside every station, so the same sum is short by
    // exactly what the FIFO holds.
    CHECK(station_jobs(dw.sol.avg.QN) + parked_jobs(dw.parked) == doctest::Approx(N).epsilon(1e-9));
    CHECK(station_jobs(dw.sol.avg.QN) < station_jobs(dd.avg.QN) - 1e-6);

    // The two rules genuinely disagree about the model, which is the point: a
    // WAITQ region served as DROP would report the DROP numbers.
    bool differs = false;
    for (std::size_t i = 0; i < mwait.get_struct().nstations; ++i)
        if (std::fabs(dw.sol.avg.QN(i, 0) - dd.avg.QN(i, 0)) > 1e-6) differs = true;
    CHECK(differs);

    double mass = 0;
    for (std::size_t s = 0; s < dw.sol.pi.size(); ++s) mass += dw.sol.pi[s];
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("ctmc waitq: a per-class cap binds its own class and no other") {
    // Two independent closed classes over Think -> Q1 -> Think, with the region
    // capping class A at one job and leaving class B unbounded. A second A job
    // must park; a B job must be admitted regardless, INCLUDING while A is at its
    // cap. That is rule 2 of the WAITQ semantics -- a fresh arrival satisfying
    // the constraints overtakes a queue stuck on a different one -- observed
    // through its consequence rather than by reading the FIFO order out.
    qn::Network<double> m("waitq_twoclass");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t a = m.add_closed_class("A", 2.0, d);
    const std::size_t b = m.add_closed_class("B", 1.0, d);
    m.set_service(d, a, Dist::exp_rate(0.5));
    m.set_service(q1, a, Dist::exp_rate(2.0));
    m.set_service(d, b, Dist::exp_rate(0.5));
    m.set_service(q1, b, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(a, a, d, q1, 1.0);
    P.set(a, a, q1, d, 1.0);
    P.set(b, b, d, q1, 1.0);
    P.set(b, b, q1, d, 1.0);
    m.link(P);
    m.add_region(std::vector<std::size_t>{q1}, std::vector<double>{1.0, -1.0}, -1.0,
                 std::vector<DropStrategy>{DropStrategy::WAITQ, DropStrategy::WAITQ});

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::WaitqSolution<double> dw =
        ctmc::solver_ctmc_waitq_analyzer(sn, ctmc::CtmcOptions());

    const std::size_t K = sn.nclasses;
    const std::size_t i1 = sn.nodes[q1 - 1].station;
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, dw.sol.chain.space);
    bool acapped = true, b_overtakes = false, balanced = true;
    for (std::size_t s = 0; s < A.rows(); ++s) {
        const double na = A(s, (i1 - 1) * K), nb = A(s, (i1 - 1) * K + 1);
        if (na > 1.0 + 1e-9) acapped = false;
        if (na > 0.5 && nb > 0.5) b_overtakes = true;
        double at_stations = 0;
        for (std::size_t c = 0; c < A.cols(); ++c) at_stations += A(s, c);
        double park = 0;
        for (std::size_t f = 0; f < dw.buf[s].size(); ++f) park += dw.buf[s][f].size();
        if (std::fabs(at_stations + park - 3.0) > 1e-9) balanced = false;
    }
    CHECK(acapped);
    CHECK(b_overtakes);
    CHECK(balanced);
    // Only class A can ever be refused here, so only class A can be parked.
    CHECK(dw.parked[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(dw.parked[0] > 1e-6);
    CHECK(station_jobs(dw.sol.avg.QN) + parked_jobs(dw.parked) ==
          doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("ctmc waitq: the rule gate distinguishes WAITQ from DROP") {
    const std::size_t Q1 = 2, Q2 = 3;
    qn::Network<double> drop = waitq_cycle3(2.0);
    drop.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, 1.0,
                    std::vector<DropStrategy>{DropStrategy::DROP});
    CHECK_FALSE(ctmc::ctmc_has_waitq_region(drop.get_struct()));

    qn::Network<double> plain = waitq_cycle3(2.0);
    CHECK_FALSE(ctmc::ctmc_has_waitq_region(plain.get_struct()));

    // The builder's DEFAULT for a class the region did not name is WAITQ, so a
    // region declared with no rule at all must be recognized as one.
    qn::Network<double> deflt = waitq_cycle3(2.0);
    deflt.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, 1.0);
    CHECK(ctmc::ctmc_has_waitq_region(deflt.get_struct()));
    CHECK_NOTHROW(ctmc::solver_ctmc_waitq_analyzer(deflt.get_struct(), ctmc::CtmcOptions()));
}

TEST_CASE("ctmc waitq: every unsupported combination is refused by name") {
    // Fixtures are built by mutating a COPY of a solved struct: each of these
    // declares a feature whose full builder setup would add nothing to what is
    // being checked, which is that the gate fires at all.
    const std::size_t Q1 = 2, Q2 = 3;
    qn::Network<double> m = waitq_cycle3(2.0);
    m.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, 1.0,
                 std::vector<DropStrategy>{DropStrategy::WAITQ});
    const qn::NetworkStruct<double> good = m.get_struct();
    REQUIRE_NOTHROW(ctmc::ctmc_check_waitq_support(good));

    {
        // A stochastic Petri net transition: a firing is atomic across its arcs
        // and has no single destination for a refused token to park against.
        qn::NetworkStruct<double> bad = good;
        bad.transparam[1] = qn::TransitionParam<double>();
        CHECK_THROWS_AS(ctmc::ctmc_check_waitq_support(bad), line::UnsupportedError);
    }
    {
        qn::NetworkStruct<double> bad = good;
        bad.fj.push_back(std::make_pair(std::size_t(1), std::size_t(2)));
        CHECK_THROWS_AS(ctmc::ctmc_check_waitq_support(bad), line::UnsupportedError);
    }
    {
        // State-dependent routing: the destination a token parks against is
        // fixed when the job is refused and would have to be re-decided later.
        qn::NetworkStruct<double> bad = good;
        bad.nodes[0].routing.assign(bad.nclasses, RoutingStrategy::RROBIN);
        CHECK_THROWS_AS(ctmc::ctmc_check_waitq_support(bad), line::UnsupportedError);
    }
    {
        // True blocking at a station is ORTHOGONAL to the region rule and is
        // refused on its own terms, not folded into WAITQ.
        qn::NetworkStruct<double> bad = good;
        REQUIRE(!bad.droprule.empty());
        bad.droprule[0][0] = DropStrategy::BAS;
        CHECK_THROWS_AS(ctmc::ctmc_check_waitq_support(bad), line::UnsupportedError);
    }

    // A region rule this generator does not implement either, reached through
    // the builder because the rule is what the region declares.
    const DropStrategy blocking[] = {DropStrategy::BAS, DropStrategy::BBS, DropStrategy::RSRD};
    for (std::size_t i = 0; i < sizeof(blocking) / sizeof(blocking[0]); ++i) {
        qn::Network<double> bm = waitq_cycle3(2.0);
        bm.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, 1.0,
                      std::vector<DropStrategy>{blocking[i]});
        CHECK_THROWS_AS(ctmc::ctmc_check_waitq_support(bm.get_struct()), line::UnsupportedError);
    }
}

TEST_CASE("ctmc waitq: an initial state the region cannot hold is rejected") {
    // Every job starts at its reference station, so putting the reference
    // station inside a region capped below the population makes the model
    // unstartable. Reporting a stationary law for it would be reporting one for
    // a chain with no initial state.
    const std::size_t THINK = 1;
    qn::Network<double> m = waitq_cycle3(3.0);
    m.add_region(std::vector<std::size_t>{THINK}, std::vector<double>{-1.0}, 1.0,
                 std::vector<DropStrategy>{DropStrategy::WAITQ});
    CHECK_THROWS_AS(ctmc::solver_ctmc_waitq_analyzer(m.get_struct(), ctmc::CtmcOptions()),
                    line::InputError);
}
