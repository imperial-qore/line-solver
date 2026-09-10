/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Finite capacity regions and parametric sensitivity, both driven through the
 * whole analyzer rather than against a hand-built chain.
 *
 * WHY THE ORACLES HERE ARE IDENTITIES AND NOT TABLES. Neither header has a
 * closed form of its own: a region is a predicate on the state space and a
 * sensitivity is a derivative of whatever the model already reports, so a
 * number copied out of either implementation would only restate it. What can be
 * checked instead is what each must satisfy no matter how it is computed.
 *
 * For the region that is INERTNESS and MEMBERSHIP: a region whose caps bound
 * nothing must leave the chain and every reported mean bit-for-bit where they
 * were, a region that does bound must leave a space every one of whose states
 * satisfies the cap, and a rule that is not DROP must stop the solve instead of
 * being served as DROP -- the one failure mode that produces plausible numbers.
 *
 * For the sensitivity it is AGREEMENT WITH A RESOLVE. dE[r]/dtheta is computed
 * here by differentiating pi Q = 0, one linear solve on the matrix the
 * stationary solve already has; the same derivative can be had by solving the
 * model twice at theta +- H and differencing the reported mean. The two share
 * the analyzer and nothing else -- no line of `ctmc_sens`, of the dQ assembly
 * or of the reward accumulation is on the second path -- so their agreement is
 * evidence about the sensitivity code and not about the model.
 */
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_fcr.h"
#include "line/solvers/ctmc/solver_ctmc_sens.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::DropStrategy;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/**
 * Delay -> FCFS Queue -> Delay, one closed class.
 *
 * Nodes are 1 = Think and 2 = Q, in creation order, which the callers name
 * below rather than rediscover.
 */
qn::Network<double> closed_dq(double njobs, double zrate, double mu) {
    qn::Network<double> m("cqn2");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(zrate));
    m.set_service(q, c, Dist::exp_rate(mu));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * Think -> Q1 -> Q2 -> Think, one closed class; nodes 1, 2, 3 in that order.
 *
 * Two queues because a region is only distinguishable from a per-station
 * capacity when it caps a SUM: no `set_capacity` on Q1 and Q2 separately can
 * express "at most two jobs between them".
 */
qn::Network<double> closed_cycle3(double njobs) {
    qn::Network<double> m("cqn3");
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

/** The aggregate job counts of one state, in the `(ist-1)*K + k` order. */
std::vector<double> aggr_row(const line::Matrix<double>& A, std::size_t s) {
    std::vector<double> nir(A.cols(), 0.0);
    for (std::size_t c = 0; c < A.cols(); ++c) nir[c] = A(s, c);
    return nir;
}

}  // namespace

TEST_CASE("ctmc fcr: a region that bounds nothing leaves the chain untouched") {
    // The strongest available check on the filter, because it needs no oracle
    // of its own: the two models differ ONLY by a region whose caps are all the
    // -1 sentinel, so every state is admissible and the filter must be the
    // identity. A filter that mis-decodes the aggregate row, or that reads the
    // sentinel as a bound of minus one job, fails here before it can fail
    // anywhere subtler.
    const double N = 3.0, z = 0.5, mu = 3.0;
    const std::size_t THINK = 1, Q = 2;

    qn::Network<double> plain = closed_dq(N, z, mu);
    qn::Network<double> regd = closed_dq(N, z, mu);
    regd.add_region(std::vector<std::size_t>{THINK, Q}, std::vector<double>{-1.0}, -1.0,
                    std::vector<DropStrategy>{DropStrategy::DROP});

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(plain.get_struct(), opt);
    const ctmc::CtmcSolution<double> d1 = ctmc::solver_ctmc_analyzer(regd.get_struct(), opt);

    REQUIRE(d1.chain.space.size() == d0.chain.space.size());
    const qn::NetworkStruct<double>& sn1 = regd.get_struct();
    const line::mva::AvgResult<double> r0 = ctmc::solver_ctmc_run_analyzer(plain.get_struct(), opt);
    const line::mva::AvgResult<double> r1 = ctmc::solver_ctmc_run_analyzer(sn1, opt);
    for (std::size_t i = 0; i < sn1.nstations; ++i) {
        CHECK(r1.QN(i, 0) == doctest::Approx(r0.QN(i, 0)).epsilon(1e-12));
        CHECK(r1.UN(i, 0) == doctest::Approx(r0.UN(i, 0)).epsilon(1e-12));
        CHECK(r1.TN(i, 0) == doctest::Approx(r0.TN(i, 0)).epsilon(1e-12));
        CHECK(r1.RN(i, 0) == doctest::Approx(r0.RN(i, 0)).epsilon(1e-12));
    }
    // The region still declares MEMBERSHIP even with no cap, and membership is
    // what the utilization estimator consults to decide whether arrivals can be
    // lost here. A closed class loses none, so the estimator must stay the
    // ordinary one and the queue's utilization must still be the busy
    // probability -- one minus the mass of the empty state.
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn1, d1.chain.space);
    const std::size_t istq = sn1.nodes[Q - 1].station;
    double idle = 0;
    for (std::size_t s = 0; s < d1.chain.space.size(); ++s)
        if (A(s, (istq - 1) * sn1.nclasses) < 0.5) idle += d1.pi[s];
    CHECK(r1.UN(istq - 1, 0) == doctest::Approx(1.0 - idle).epsilon(1e-9));
}

TEST_CASE("ctmc fcr: a region capping two stations censors the state space") {
    // Four jobs in a cycle, at most two of them between Q1 and Q2. The cap is
    // below the population on purpose: a cap at or above N would be satisfied
    // by every state and would not distinguish a filter from a no-op.
    const double N = 4.0, cap = 2.0;
    const std::size_t Q1 = 2, Q2 = 3;

    qn::Network<double> plain = closed_cycle3(N);
    qn::Network<double> regd = closed_cycle3(N);
    regd.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, cap,
                    std::vector<DropStrategy>{DropStrategy::DROP});

    const qn::NetworkStruct<double>& sn0 = plain.get_struct();
    const qn::NetworkStruct<double>& sn1 = regd.get_struct();
    const std::size_t K = sn1.nclasses;
    const std::size_t i1 = sn1.nodes[Q1 - 1].station, i2 = sn1.nodes[Q2 - 1].station;

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(sn0, opt);
    const ctmc::CtmcSolution<double> d1 = ctmc::solver_ctmc_analyzer(sn1, opt);
    CHECK(d1.chain.space.size() < d0.chain.space.size());

    // Every surviving state satisfies the cap. This is the property the region
    // means, stated on the states themselves rather than on any mean derived
    // from them.
    const line::Matrix<double> A1 = ctmc::ctmc_state_space_aggr(sn1, d1.chain.space);
    for (std::size_t s = 0; s < d1.chain.space.size(); ++s) {
        const double inregion = A1(s, (i1 - 1) * K) + A1(s, (i2 - 1) * K);
        CHECK(inregion <= cap + 1e-9);
        // A DROP region on a CLOSED model censors rather than destroys: the
        // forbidden states are unreachable, the transitions into them are
        // deleted and the job stays where it was. The population is therefore
        // still conserved, which a filter that leaked jobs would break.
        double tot = 0;
        for (std::size_t c = 0; c < A1.cols(); ++c) tot += A1(s, c);
        CHECK(tot == doctest::Approx(N).epsilon(1e-12));
    }

    // The surviving COUNT is the number of unfiltered states the predicate
    // admits, which reaches the same answer through `ctmc_region_admissible`
    // applied one state at a time to the space of the region-less model. It
    // pins the filter to its own predicate without hard-coding a size.
    const line::Matrix<double> A0 = ctmc::ctmc_state_space_aggr(sn0, d0.chain.space);
    std::size_t admissible = 0;
    for (std::size_t s = 0; s < d0.chain.space.size(); ++s)
        if (ctmc::ctmc_region_admissible(sn1, aggr_row(A0, s))) ++admissible;
    CHECK(admissible == d1.chain.space.size());

    // And the same count from combinatorics: the states of a single-class
    // closed cycle are the compositions of N over the three stations, of which
    // the region keeps those whose last two parts sum to at most the cap.
    std::size_t all = 0, kept = 0;
    for (int a = 0; a <= static_cast<int>(N); ++a)
        for (int b = 0; a + b <= static_cast<int>(N); ++b) {
            ++all;
            if (b + (static_cast<int>(N) - a - b) <= static_cast<int>(cap)) ++kept;
        }
    CHECK(all == d0.chain.space.size());
    CHECK(kept == d1.chain.space.size());

    // The censored chain is still a chain: its stationary law is normalized.
    double mass = 0;
    for (std::size_t s = 0; s < d1.pi.size(); ++s) mass += d1.pi[s];
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("ctmc fcr: every region rule other than DROP is refused") {
    // WAITQ parks the refused job in a per-region buffer outside every station,
    // and BAS, BBS and RSRD hold it at the sender; all four are extra state
    // this generator does not carry. Serving them as DROP would return a
    // stationary law of a chain the model does not describe, so the refusal is
    // the correct behaviour and its ABSENCE is what these checks guard. Only
    // the refusal is asserted, never the wording.
    const std::size_t THINK = 1, Q = 2;

    // WAITQ is reached on the ordinary solve path, not only by calling the
    // rule check directly, because it is the builder's DEFAULT for any class
    // the region did not name -- passing no rule at all must not slip through.
    qn::Network<double> waitq = closed_dq(2.0, 0.5, 3.0);
    waitq.add_region(std::vector<std::size_t>{THINK, Q}, std::vector<double>{-1.0}, 1.0);
    CHECK_THROWS_AS(ctmc::solver_ctmc_run_analyzer(waitq.get_struct(), ctmc::CtmcOptions()),
                    line::UnsupportedError);

    const DropStrategy refused[] = {DropStrategy::WAITQ, DropStrategy::BAS, DropStrategy::BBS,
                                    DropStrategy::RSRD, DropStrategy::RETRIAL};
    for (std::size_t i = 0; i < sizeof(refused) / sizeof(refused[0]); ++i) {
        qn::Network<double> m = closed_dq(2.0, 0.5, 3.0);
        m.add_region(std::vector<std::size_t>{THINK, Q}, std::vector<double>{-1.0}, 1.0,
                     std::vector<DropStrategy>{refused[i]});
        CHECK_THROWS_AS(ctmc::ctmc_check_region_rules(m.get_struct()), line::UnsupportedError);
    }

    // DROP is the one rule that passes, so the loop above is testing the gate
    // and not merely the presence of a region.
    qn::Network<double> ok = closed_dq(2.0, 0.5, 3.0);
    ok.add_region(std::vector<std::size_t>{THINK, Q}, std::vector<double>{-1.0}, 1.0,
                  std::vector<DropStrategy>{DropStrategy::DROP});
    CHECK_NOTHROW(ctmc::ctmc_check_region_rules(ok.get_struct()));
}

TEST_CASE("ctmc fcr: the drop mask marks the members of a DROP region only") {
    // The mask decides which stations lose the offered-rate utilization
    // estimator, so it has to be MEMBERSHIP OF A DROP REGION and not either
    // half of that on its own: a member of a WAITQ region loses nothing, and a
    // non-member of a DROP region is unaffected by it.
    const std::size_t Q1 = 2, Q2 = 3;

    qn::Network<double> m = closed_cycle3(3.0);
    m.add_region(std::vector<std::size_t>{Q1}, std::vector<double>{-1.0}, 2.0,
                 std::vector<DropStrategy>{DropStrategy::DROP});
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t i1 = sn.nodes[Q1 - 1].station;
    const std::vector<bool> mask = ctmc::ctmc_in_drop_region(sn);
    REQUIRE(mask.size() == sn.nstations);
    for (std::size_t i = 0; i < mask.size(); ++i)
        CHECK(static_cast<bool>(mask[i]) == (i == i1 - 1));

    // A model with no region at all marks nothing, which is the baseline every
    // model without a region relies on.
    qn::Network<double> bare = closed_cycle3(3.0);
    const std::vector<bool> none = ctmc::ctmc_in_drop_region(bare.get_struct());
    for (std::size_t i = 0; i < none.size(); ++i) CHECK(static_cast<bool>(none[i]) == false);

    // A WAITQ region has members but drops nothing. The solve refuses this
    // model, so the mask is only reachable directly; it must still report the
    // truth, since a mask that keyed on membership alone would call these two
    // stations lossy.
    qn::Network<double> wq = closed_cycle3(3.0);
    wq.add_region(std::vector<std::size_t>{Q1, Q2}, std::vector<double>{-1.0}, 2.0,
                  std::vector<DropStrategy>{DropStrategy::WAITQ});
    const std::vector<bool> wmask = ctmc::ctmc_in_drop_region(wq.get_struct());
    for (std::size_t i = 0; i < wmask.size(); ++i) CHECK(static_cast<bool>(wmask[i]) == false);
}

TEST_CASE("ctmc sens: dE[Q]/dmu agrees with a central difference of the solved mean") {
    // Delay -> FCFS Queue -> Delay, three jobs, and the reward is the queue's
    // own length, so the sensitivity is the derivative of a number the AvgTable
    // already reports. That is what makes the resolve-and-difference check
    // possible at all, and the two paths meet only at `solver_ctmc_analyzer`.
    const double N = 3.0, z = 0.5, mu = 3.0;
    const std::size_t Q = 2;

    qn::Network<double> m = closed_dq(N, z, mu);
    const std::size_t istq = m.get_struct().nodes[Q - 1].station;
    const std::size_t K = m.get_struct().nclasses;
    // The reward is declared on the model, as a caller would, and the per-state
    // vector below is that same function evaluated on the aggregate rows: one
    // definition, so the sensitivity and the mean it is differentiated from
    // cannot drift apart.
    m.set_reward("qlen", [istq, K](const std::vector<double>& n) { return n[(istq - 1) * K]; });
    const qn::NetworkStruct<double> sn = m.get_struct();
    REQUIRE(sn.reward.size() == 1);

    // Perturbing theta rewrites the station's service distribution and
    // refreshes the derived rate matrix. Only `refresh_rates` is needed: the
    // routing, the chains and the capacities of a closed model are functions of
    // the topology and the populations, none of which a rate can move.
    const std::function<void(qn::NetworkStruct<double>&, double)> set_mu =
        [istq](qn::NetworkStruct<double>& s, double v) {
            s.service[istq - 1][0] = Dist::exp_rate(v);
            s.refresh_rates();
        };

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> base = ctmc::solver_ctmc_analyzer(sn, opt);
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, base.chain.space);
    std::vector<double> reward(base.chain.space.size(), 0.0);
    for (std::size_t s = 0; s < reward.size(); ++s) reward[s] = sn.reward[0].fn(aggr_row(A, s));

    ctmc::CtmcSensParam<double> p;
    p.name = "mu";
    p.value = mu;
    p.set = set_mu;
    const ctmc::CtmcSens<double> sens = ctmc::solver_ctmc_sensitivity(sn, opt, p, reward);

    // pi sums to one for EVERY theta, so its derivative sums to zero. The
    // normalization is the last equation of the sensitivity solve, and this is
    // the residual of that equation read back out.
    double dtot = 0;
    for (std::size_t s = 0; s < sens.dpi.size(); ++s) dtot += sens.dpi[s];
    CHECK(dtot == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));

    // The oracle. H is far larger than the internal step for dQ, so the two
    // truncation errors are unrelated: this compares the derivative of the
    // stationary law against the difference of two independently solved models.
    const double H = 1e-3;
    qn::NetworkStruct<double> up = sn, dn = sn;
    set_mu(up, mu + H);
    set_mu(dn, mu - H);
    const double qup = ctmc::solver_ctmc_run_analyzer(up, opt).QN(istq - 1, 0);
    const double qdn = ctmc::solver_ctmc_run_analyzer(dn, opt).QN(istq - 1, 0);
    const double fd = (qup - qdn) / (2.0 * H);
    CHECK(sens.S == doctest::Approx(fd).epsilon(1e-4));

    // A faster server holds fewer jobs, at any population and any think time,
    // so the sign is not a matter of the discretization.
    CHECK(sens.S < 0.0);

    // The scaled form is the elasticity, and E[r] is the queue length the
    // AvgTable reports; checking the definition also checks that the reward
    // accumulated over pi is the same mean the analyzer assembled by its own
    // route.
    const double er = ctmc::solver_ctmc_run_analyzer(sn, opt).QN(istq - 1, 0);
    REQUIRE(sens.scaled_valid);
    CHECK(sens.SS == doctest::Approx(mu / er * sens.S).epsilon(1e-9));
}

TEST_CASE("ctmc sens: the symbolic method needs its backend and a stray one is rejected") {
    // 'symbolic' IS delivered by this port now (sens_detail::symbolic_sensitivity),
    // and like the reference it is served by the line-sage-rest backend. Whether
    // that backend resolves is a property of the machine, so what is pinned here
    // is the refusal when it is switched OFF: the request must fail by name
    // rather than fall back to the difference quotient, which would report a
    // different method's answer under the symbolic one's. An unlisted method
    // string is a caller error instead, hence the different exception.
    const std::size_t Q = 2;
    qn::Network<double> m = closed_dq(2.0, 0.5, 3.0);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const std::size_t istq = sn.nodes[Q - 1].station;

    ctmc::CtmcSensParam<double> p;
    p.name = "mu";
    p.value = 3.0;
    p.set = [istq](qn::NetworkStruct<double>& s, double v) {
        s.service[istq - 1][0] = Dist::exp_rate(v);
        s.refresh_rates();
    };
    const ctmc::CtmcOptions opt;
    ctmc::CtmcSymbolicOptions nobackend;
    nobackend.backend = "none";
    CHECK_THROWS_AS(
        ctmc::solver_ctmc_sensitivity(sn, opt, p, std::vector<double>(), "symbolic", nobackend),
        line::sym::SymEngineError);
    CHECK_THROWS_AS(ctmc::solver_ctmc_sensitivity(sn, opt, p, std::vector<double>(), "fd2"),
                    line::InputError);

    // A parameter with no setter cannot be applied to the model, so there is no
    // perturbed generator to difference; that is an input error and not a
    // silently zero sensitivity.
    ctmc::CtmcSensParam<double> noset;
    noset.name = "mu";
    noset.value = 3.0;
    CHECK_THROWS_AS(ctmc::solver_ctmc_sensitivity(sn, opt, noset), line::InputError);
}

TEST_CASE("ctmc sens: the ranking is the same table sorted by scaled influence") {
    // Two parameters of different units -- a service rate and a think rate --
    // which is exactly why the ranking sorts on the SCALED sensitivity: the
    // unscaled ones are not comparable. The check is that the table carries the
    // same numbers a per-parameter call returns and orders them by descending
    // magnitude; which parameter wins is read off those numbers rather than
    // asserted, so the test states the ordering property and not a result.
    const double N = 3.0, z = 0.5, mu = 3.0;
    const std::size_t THINK = 1, Q = 2;

    qn::Network<double> m = closed_dq(N, z, mu);
    const qn::NetworkStruct<double> sn = m.get_struct();
    const std::size_t K = sn.nclasses;
    const std::size_t istq = sn.nodes[Q - 1].station, istd = sn.nodes[THINK - 1].station;

    ctmc::CtmcSensParam<double> pmu, pz;
    pmu.name = "mu";
    pmu.value = mu;
    pmu.set = [istq](qn::NetworkStruct<double>& s, double v) {
        s.service[istq - 1][0] = Dist::exp_rate(v);
        s.refresh_rates();
    };
    pz.name = "zrate";
    pz.value = z;
    pz.set = [istd](qn::NetworkStruct<double>& s, double v) {
        s.service[istd - 1][0] = Dist::exp_rate(v);
        s.refresh_rates();
    };

    const ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> base = ctmc::solver_ctmc_analyzer(sn, opt);
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, base.chain.space);
    std::vector<double> reward(base.chain.space.size(), 0.0);
    for (std::size_t s = 0; s < reward.size(); ++s) reward[s] = A(s, (istq - 1) * K);

    std::vector<ctmc::CtmcSensParam<double>> params;
    params.push_back(pmu);
    params.push_back(pz);
    const std::vector<ctmc::CtmcSensRank<double>> rows =
        ctmc::solver_ctmc_sensitivity_ranking(sn, opt, params, reward);
    REQUIRE(rows.size() == 2);

    const ctmc::CtmcSens<double> smu = ctmc::solver_ctmc_sensitivity(sn, opt, pmu, reward);
    const ctmc::CtmcSens<double> sz = ctmc::solver_ctmc_sensitivity(sn, opt, pz, reward);
    // Both parameters survive the sort, each carrying its own numbers: a rank
    // table that lost one, or that paired a name with another's sensitivity,
    // would still be sorted.
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const ctmc::CtmcSens<double>& want = rows[i].parameter == "mu" ? smu : sz;
        CHECK(rows[i].S == doctest::Approx(want.S).epsilon(1e-12));
        CHECK(rows[i].SS == doctest::Approx(want.SS).epsilon(1e-12));
        CHECK(rows[i].value == doctest::Approx(rows[i].parameter == "mu" ? mu : z));
        CHECK(rows[i].scaled_valid);
    }
    CHECK(rows[0].parameter != rows[1].parameter);
    CHECK(std::fabs(rows[0].SS) >= std::fabs(rows[1].SS));
    // The order agrees with the magnitudes computed outside the ranking.
    const bool mu_first = std::fabs(smu.SS) >= std::fabs(sz.SS);
    CHECK(rows[0].parameter == (mu_first ? "mu" : "zrate"));

    // The signs are opposite and are retained: a faster server shortens the
    // queue, a shorter think time lengthens it. A ranking that sorted on the
    // signed value rather than the magnitude would order these two the other
    // way round.
    CHECK(smu.SS < 0.0);
    CHECK(sz.SS > 0.0);
}
