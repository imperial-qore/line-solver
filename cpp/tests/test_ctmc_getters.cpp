/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The `@SolverCTMC` accessor surface: `getGenerator` / `getInfGen`,
 * `getStateSpace` / `getStateSpaceAggr`, the `getTranProb*` family and the two
 * symbolic refusals.
 *
 * AN ACCESSOR HAS NO NUMBERS OF ITS OWN, so every oracle here is an identity
 * against something that was computed by a different route. The generator is
 * checked against the definition of a generator and against its own filtration;
 * the state space accessors against the chain the analyzer solved; the
 * transient family against the STATIONARY probability family at a large t,
 * which is the one comparison that exercises the labelling and the occupancy
 * together. Reading a value back out of the accessor under test would prove
 * only that it is deterministic.
 *
 * WHY sum_a filt[a] IS THE LOAD-BEARING CHECK. The generator sums every
 * synchronization's contribution into one entry, so a filtration that lost or
 * double counted an event still produces a perfectly valid-looking Q. The only
 * thing that catches it is re-summing the per-event matrices and demanding the
 * off-diagonal back, which is what `getGenerator`'s caller does when it splits
 * the chain on one event to build a response-time MAP.
 */
#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_symbolic.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
namespace sym = line::sym;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> FCFS Queue (capacity K) -> Sink, one open class. */
qn::Network<double> acc_mm1k(double lambda, double mu, int K) {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** Delay -> FCFS Queue -> Delay, one closed class. */
qn::Network<double> acc_cqn(double njobs) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(0.5));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/**
 * The same topology with a matrix-exponential service: D0 carries a NEGATIVE
 * off-diagonal, so the stationary vector is a signed measure and no per-state
 * or transient probability exists. Patched in place because no builder entry
 * point produces an ME.
 *
 * The patched struct is RETURNED rather than re-fetched by the caller, since
 * `get_struct` re-runs the whole refresh chain on every call and the patch is
 * not something the builder can reproduce.
 */
qn::NetworkStruct<double>& make_me(qn::Network<double>& m) {
    qn::NetworkStruct<double>& sn = const_cast<qn::NetworkStruct<double>&>(m.get_struct());
    line::Matrix<double> D0(2, 2, 0.0), D1(2, 2, 0.0);
    D0(0, 0) = -2.0;
    D0(0, 1) = -0.5;
    D0(1, 1) = -3.0;
    D1(0, 1) = 2.5;
    D1(1, 0) = 3.0;
    sn.service[1][0].D0 = D0;
    sn.service[1][0].D1 = D1;
    sn.service[1][0].disabled = false;
    return sn;
}

/** True when two label rows agree, which is how a marginal groups states. */
bool same_row(const line::Matrix<double>& A, std::size_t i, std::size_t j) {
    for (std::size_t c = 0; c < A.cols(); ++c)
        if (A(i, c) != A(j, c)) return false;
    return true;
}

}  // namespace

TEST_CASE("ctmc getters: the event filtration re-sums to the off-diagonal of the generator") {
    const int K = 4;
    qn::Network<double> m = acc_mm1k(0.6, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const ctmc::CtmcGenerator<double> g = ctmc::ctmc_get_generator(sn, opt);

    REQUIRE(g.Q.rows() == g.space.size());
    REQUIRE(g.Q.rows() == g.Q.cols());
    // The accessor forces `keep_filtration` on, so an empty filtration here
    // would mean the flag never reached the assembly.
    REQUIRE(!g.filt.empty());
    // One matrix per synchronization: the pairing is the only thing that says
    // WHICH event a filtration matrix belongs to.
    CHECK(g.filt.size() == g.sync.size());

    for (std::size_t s = 0; s < g.Q.rows(); ++s) {
        double row = 0;
        for (std::size_t t = 0; t < g.Q.cols(); ++t) {
            if (s != t) CHECK(g.Q(s, t) >= 0.0);
            row += g.Q(s, t);
        }
        // A generator's rows sum to zero by construction; a rate that was added
        // without being taken off the diagonal breaks this and nothing else.
        CHECK(row == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        CHECK(g.Q(s, s) <= 0.0);
    }

    // The identity the filtration exists for. It holds OFF the diagonal only:
    // a self-loop is a real firing that `make_infgen` cancels in Q, so the
    // diagonals of the two sides are answering different questions.
    for (std::size_t s = 0; s < g.Q.rows(); ++s)
        for (std::size_t t = 0; t < g.Q.cols(); ++t) {
            if (s == t) continue;
            double acc = 0;
            for (std::size_t a = 0; a < g.filt.size(); ++a) acc += g.filt[a](s, t);
            CHECK(acc == doctest::Approx(g.Q(s, t)).epsilon(1e-12));
        }

    // getInfGen is a pure alias in the reference, so it must not be a second
    // implementation that can drift from the first.
    const ctmc::CtmcGenerator<double> gi = ctmc::ctmc_get_infgen(sn, opt);
    REQUIRE(gi.Q.rows() == g.Q.rows());
    for (std::size_t s = 0; s < g.Q.rows(); ++s)
        for (std::size_t t = 0; t < g.Q.cols(); ++t) CHECK(gi.Q(s, t) == g.Q(s, t));
}

TEST_CASE("ctmc getters: the generator accessor refuses a solve that kept no filtration") {
    qn::Network<double> m = acc_mm1k(0.5, 1.0, 3);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = 3;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    REQUIRE(d.chain.filt.empty());

    // The filtration is not recoverable from Q, so answering with an empty one
    // would hand back a chain that looks eventless rather than one that was not
    // asked to record its events.
    CHECK_THROWS_AS(ctmc::ctmc_get_generator(sn, d), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_infgen(sn, d), line::InputError);

    ctmc::CtmcOptions kept = opt;
    kept.keep_filtration = true;
    const ctmc::CtmcSolution<double> dk = ctmc::solver_ctmc_analyzer(sn, kept);
    const ctmc::CtmcGenerator<double> g = ctmc::ctmc_get_generator(sn, dk);
    const ctmc::CtmcGenerator<double> gi = ctmc::ctmc_get_infgen(sn, dk);
    CHECK(g.filt.size() == gi.filt.size());
    // The accessor reports the chain that was SOLVED, so its rows are the rows
    // the stationary vector is indexed by; a caller may cross-index the two.
    CHECK(g.space.size() == dk.pi.size());
}

TEST_CASE("ctmc getters: the state space accessors reproduce the solved chain") {
    qn::Network<double> m = acc_cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    const ctmc::CtmcStateSpace<double> ss = ctmc::ctmc_get_state_space(sn, d);

    REQUIRE(ss.space.size() == d.chain.space.size());
    REQUIRE(ss.flat.rows() == d.chain.space.size());
    std::size_t total = 0;
    for (std::size_t f = 0; f < ss.node_width.size(); ++f) total += ss.node_width[f];
    // The widths ARE the cut points the reference needs to split the flat matrix
    // back into per-node blocks, so they must tile it exactly.
    CHECK(total == ss.flat.cols());
    CHECK(ss.node_width.size() == sn.stateful_nodes.size());

    // The flat form is the concatenation of the blocks in stateful-node order.
    // Comparing against the chain rather than against the accessor's own split
    // is what makes this an oracle instead of a restatement.
    for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
        std::size_t c = 0;
        for (std::size_t f = 0; f < d.chain.space[s].local.size(); ++f) {
            REQUIRE(d.chain.space[s].local[f].size() == ss.node_width[f]);
            for (std::size_t j = 0; j < d.chain.space[s].local[f].size(); ++j)
                CHECK(ss.flat(s, c++) == d.chain.space[s].local[f][j]);
        }
    }

    const line::Matrix<double> A = ctmc::ctmc_get_state_space_aggr(sn, d);
    const line::Matrix<double> A2 = ctmc::ctmc_state_space_aggr(sn, d.chain.space);
    REQUIRE(A.rows() == A2.rows());
    REQUIRE(A.cols() == A2.cols());
    for (std::size_t s = 0; s < A.rows(); ++s)
        for (std::size_t c = 0; c < A.cols(); ++c) CHECK(A(s, c) == A2(s, c));

    // The aggregate weighted by the stationary vector IS the mean queue length,
    // which ties the column block order `(ist-1)*K + k` to a number the analyzer
    // produced without ever forming this matrix.
    double q = 0;
    for (std::size_t s = 0; s < A.rows(); ++s) q += d.pi[s] * A(s, 1);
    CHECK(q == doctest::Approx(d.avg.QN(1, 0)).epsilon(1e-9));

    // The forms that solve for themselves must land on the same chain: the
    // analyzer is deterministic, so a difference would mean the accessor and
    // the solve disagree about which component was kept.
    const ctmc::CtmcStateSpace<double> ss2 = ctmc::ctmc_get_state_space(sn, opt);
    const line::Matrix<double> A3 = ctmc::ctmc_get_state_space_aggr(sn, opt);
    CHECK(ss2.flat.rows() == ss.flat.rows());
    CHECK(ss2.flat.cols() == ss.flat.cols());
    CHECK(A3.rows() == A.rows());
    CHECK(A3.cols() == A.cols());
}

TEST_CASE("ctmc getters: the transient family converges to the stationary family") {
    qn::Network<double> m = acc_cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    const ctmc::CtmcSolution<double> st = ctmc::solver_ctmc_analyzer(sn, opt);
    // One integration serves all four queries: they differ only in how they
    // LABEL the same occupancy, which is the property under test.
    const ctmc::CtmcTransient<double> tr =
        ctmc::solver_ctmc_transient_analyzer(sn, opt, 0.0, 200.0);
    const std::size_t ist = 2, ind = sn.node_of_station(ist);
    const std::size_t isf = sn.stateful_of_station(ist);
    REQUIRE(isf != 0);

    const ctmc::CtmcTranProb<double> ps = ctmc::ctmc_get_tran_prob_sys(sn, tr);
    const ctmc::CtmcTranProb<double> pa = ctmc::ctmc_get_tran_prob_sys_aggr(sn, tr);
    const ctmc::CtmcTranProb<double> pn = ctmc::ctmc_get_tran_prob(sn, tr, ind);
    const ctmc::CtmcTranProb<double> pna = ctmc::ctmc_get_tran_prob_aggr(sn, tr, ind);

    REQUIRE(!tr.t.empty());
    const std::size_t n = tr.chain.chain.space.size();
    const std::size_t last = tr.t.size() - 1;
    REQUIRE(n > 0);
    REQUIRE(ps.pit.rows() == tr.t.size());
    REQUIRE(ps.pit.cols() == n);
    // Every query reports the SAME occupancy; only the labels differ. The four
    // sharing one pi(t) is what makes them four views and not four solves.
    CHECK(pa.pit.rows() == ps.pit.rows());
    CHECK(pn.t.size() == ps.t.size());
    CHECK(pna.labels.rows() == n);

    // pi(t) is a distribution at every time, not only in the limit.
    for (std::size_t i = 0; i < ps.pit.rows(); ++i) {
        double mass = 0;
        for (std::size_t s = 0; s < n; ++s) mass += ps.pit(i, s);
        CHECK(mass == doctest::Approx(1.0).epsilon(1e-6));
    }

    // The system labels are the whole network state, the per-node labels are
    // one block of it: both are read off the chain the transient solved.
    REQUIRE(pn.labels.cols() == tr.chain.chain.space[0].local[isf - 1].size());
    for (std::size_t s = 0; s < n; ++s)
        for (std::size_t j = 0; j < pn.labels.cols(); ++j)
            CHECK(pn.labels(s, j) == tr.chain.chain.space[s].local[isf - 1][j]);

    // The per-node aggregate decodes the same marginal the station-block
    // aggregate holds, by a different route: one goes through the node, the
    // other through the station column block.
    const line::Matrix<double> A = ctmc::ctmc_state_space_aggr(sn, tr.chain.chain.space);
    for (std::size_t s = 0; s < n; ++s)
        CHECK(pna.labels(s, 0) == doctest::Approx(A(s, (ist - 1) * sn.nclasses)).epsilon(1e-12));

    // By t = 200 on a chain with rates of order 1 the transient has died out, so
    // each transient query must reproduce its stationary counterpart:
    // getTranProbSys -> getProbSys, and the aggregates likewise once the states
    // sharing a label are summed.
    for (std::size_t s = 0; s < n; ++s) {
        const double joint = ctmc::solver_ctmc_joint(sn, st, tr.chain.chain.space[s]);
        CHECK(ps.pit(last, s) == doctest::Approx(joint).epsilon(1e-4));
    }
    for (std::size_t s = 0; s < n; ++s) {
        double grp = 0;
        for (std::size_t s2 = 0; s2 < n; ++s2)
            if (same_row(pa.labels, s, s2)) grp += pa.pit(last, s2);
        CHECK(grp == doctest::Approx(ctmc::solver_ctmc_jointaggr(sn, st,
                                                                 tr.chain.chain.space[s]))
                         .epsilon(1e-4));
    }
    for (std::size_t s = 0; s < n; ++s) {
        double grp = 0;
        for (std::size_t s2 = 0; s2 < n; ++s2)
            if (same_row(pna.labels, s, s2)) grp += pna.pit(last, s2);
        const std::vector<double> marg =
            ctmc::solver_ctmc_margaggr(sn, st, tr.chain.chain.space[s]);
        CHECK(grp == doctest::Approx(marg[ist - 1]).epsilon(1e-4));
    }

    // The labels are also what the analyzer integrates its own means against,
    // so summing them under pi(t) has to give back QNt at that instant. This is
    // the check that the column block order was not permuted.
    for (std::size_t i = 0; i < tr.t.size(); i += 5) {
        double qn = 0;
        for (std::size_t s = 0; s < n; ++s)
            qn += pa.pit(i, s) * pa.labels(s, (ist - 1) * sn.nclasses);
        CHECK(qn == doctest::Approx(tr.QNt[ist - 1][0][i]).epsilon(1e-9));
    }
}

TEST_CASE("ctmc getters: a transient query without a finite horizon is refused") {
    qn::Network<double> m = acc_cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    const std::size_t ind = sn.node_of_station(2);
    const double inf = std::numeric_limits<double>::infinity();

    // pi(t) on [0, Inf) is the stationary vector, which the getProb* family
    // already answers, so an unbounded horizon is a caller error rather than a
    // long integration.
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sn, opt, ind, 0.0, inf), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_aggr(sn, opt, ind, 0.0, inf), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys(sn, opt, 0.0, inf), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys_aggr(sn, opt, 0.0, inf), line::InputError);
    // An empty or reversed span has no trajectory to report either.
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys(sn, opt, 5.0, 5.0), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys(sn, opt, 5.0, 1.0), line::InputError);
}

TEST_CASE("ctmc getters: a stateless or out-of-range node is refused") {
    qn::Network<double> m = acc_cqn(2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;
    const ctmc::CtmcTransient<double> tr =
        ctmc::solver_ctmc_transient_analyzer(sn, opt, 0.0, 1.0);

    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sn, tr, 0), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sn, tr, sn.nodes.size() + 1), line::InputError);

    // A Sink holds no block of the network state, so there is nothing to slice
    // out for it: an empty label matrix would read as a node with no states
    // rather than as a query that cannot be asked.
    qn::Network<double> mo = acc_mm1k(0.5, 1.0, 3);
    const qn::NetworkStruct<double>& sno = mo.get_struct();
    ctmc::CtmcOptions oo;
    oo.cutoff = 3;
    const ctmc::CtmcTransient<double> tro =
        ctmc::solver_ctmc_transient_analyzer(sno, oo, 0.0, 1.0);
    REQUIRE(sno.sinkNode != 0);
    REQUIRE(sno.stateful_index(sno.sinkNode) == 0);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sno, tro, sno.sinkNode), line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_aggr(sno, tro, sno.sinkNode), line::InputError);

    // A trajectory integrated from ANOTHER model. The overloads that take a
    // ready CtmcTransient cannot verify the pairing in general, but a struct
    // with more stateful nodes than the trajectory has blocks would be read
    // past the end of a state, so that case is refused rather than indexed.
    qn::Network<double> wide("wide");
    const std::size_t dl = wide.add_delay("Think");
    const std::size_t q1 = wide.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = wide.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t cw = wide.add_closed_class("C1", 1.0, dl);
    wide.set_service(dl, cw, Dist::exp_rate(1.0));
    wide.set_service(q1, cw, Dist::exp_rate(2.0));
    wide.set_service(q2, cw, Dist::exp_rate(2.0));
    qn::RoutingMatrix<double> Pw;
    Pw.set(dl, q1, 1.0);
    Pw.set(q1, q2, 1.0);
    Pw.set(q2, dl, 1.0);
    wide.link(Pw);
    const qn::NetworkStruct<double>& snw = wide.get_struct();
    REQUIRE(snw.stateful_nodes.size() > sn.stateful_nodes.size());
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(snw, tr, snw.node_of_station(3)),
                    line::InputError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_aggr(snw, tr, snw.node_of_station(3)),
                    line::InputError);
}

TEST_CASE("ctmc getters: a Source node aggregates to zero, not Inf") {
    const int K = 3;
    qn::Network<double> m = acc_mm1k(0.5, 1.0, K);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const std::size_t src = sn.node_of_station(1);
    const ctmc::CtmcTranProb<double> p =
        ctmc::ctmc_get_tran_prob_aggr(sn, opt, src, 0.0, 0.5);

    // `to_marginal` reports an infinite reservoir at an EXT station. That is a
    // statement about the encoding, and an Inf in a label matrix would poison
    // every aggregate a caller builds by summing it under pi(t).
    REQUIRE(p.labels.cols() == sn.nclasses);
    for (std::size_t s = 0; s < p.labels.rows(); ++s) {
        CHECK(std::isfinite(p.labels(s, 0)));
        CHECK(p.labels(s, 0) == 0.0);
    }

    // The system aggregate must agree with it: the two decode the same station
    // by different routes, so a Source treated as a reservoir in one and as
    // empty in the other would make the per-node and system views inconsistent.
    const ctmc::CtmcTranProb<double> ps =
        ctmc::ctmc_get_tran_prob_sys_aggr(sn, opt, 0.0, 0.5);
    REQUIRE(ps.labels.rows() == p.labels.rows());
    for (std::size_t s = 0; s < ps.labels.rows(); ++s) CHECK(ps.labels(s, 0) == 0.0);

    // The unaggregated per-node view is the ENCODING row, not the marginal, so
    // it keeps whatever the Source's phase block holds and is not zeroed.
    const ctmc::CtmcTranProb<double> praw = ctmc::ctmc_get_tran_prob(sn, opt, src, 0.0, 0.5);
    CHECK(praw.labels.rows() == p.labels.rows());
    CHECK(praw.labels.cols() > 0);
}

TEST_CASE("ctmc getters: the ME gate stops every transient query") {
    qn::Network<double> healthy = acc_mm1k(0.5, 1.0, 3);
    ctmc::CtmcOptions opt;
    opt.cutoff = 3;
    const ctmc::CtmcTransient<double> tr =
        ctmc::solver_ctmc_transient_analyzer(healthy.get_struct(), opt, 0.0, 1.0);

    qn::Network<double> m = acc_mm1k(0.5, 1.0, 3);
    const qn::NetworkStruct<double>& sn = make_me(m);
    const std::size_t ind = sn.node_of_station(2);

    // Uniformization on a signed generator diverges, so a transient answer does
    // not exist here at all. The gate is a property of the STRUCT, so it must
    // fire before any state is indexed: the trajectory handed in below comes
    // from the healthy model of the same shape and is never reached.
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sn, tr, ind), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_aggr(sn, tr, ind), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys(sn, tr), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys_aggr(sn, tr), line::UnsupportedError);

    // The gate must also come BEFORE the timespan check on the solving forms,
    // or an ME model with a bad horizon would be reported as a horizon problem.
    const double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob(sn, opt, ind, 0.0, inf), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_aggr(sn, opt, ind, 0.0, inf),
                    line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys(sn, opt, 0.0, inf), line::UnsupportedError);
    CHECK_THROWS_AS(ctmc::ctmc_get_tran_prob_sys_aggr(sn, opt, 0.0, inf),
                    line::UnsupportedError);
}

TEST_CASE("ctmc symbolic: the generator is assembled with no backend at all") {
    qn::Network<double> m = acc_mm1k(0.5, 1.0, 2);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;

    // The generator is LINEAR in the event symbols -- one numeric filtration per
    // synchronization, scaled by x1..xE -- so it needs no computer algebra. Only
    // the stationary solve does; see the case below.
    const ctmc::CtmcSymbolicGenerator<double> g = ctmc::ctmc_symbolic_generator(sn, opt);
    const ctmc::CtmcGenerator<double> num = ctmc::ctmc_get_generator(sn, opt);
    REQUIRE(g.Q.size() == g.space.size());
    REQUIRE(g.space.size() == num.space.size());
    // `symbols`, `filt`, `terms` and `rate0` are all indexed by synchronization,
    // and an event with no positive rate KEEPS its slot rather than renumbering
    // every later symbol.
    CHECK(g.symbols.size() == g.sync.size());
    CHECK(g.rate0.size() == g.sync.size());
    CHECK(g.terms.size() == g.sync.size());
    CHECK(!g.active_symbols().empty());

    // Substituting each event's minimum positive rate must reproduce the numeric
    // generator exactly: the filtration was divided by that rate, so the nominal
    // value of x_e is it. The header names this the check to run before trusting
    // a printed normal form.
    const line::Matrix<double> Qx = ctmc::ctmc_symbolic_eval_infgen(g, g.rate0);
    REQUIRE(Qx.rows() == num.Q.rows());
    REQUIRE(Qx.cols() == num.Q.cols());
    for (std::size_t i = 0; i < Qx.rows(); ++i)
        for (std::size_t j = 0; j < Qx.cols(); ++j) {
            CAPTURE(i);
            CAPTURE(j);
            CHECK(std::fabs(Qx(i, j) - num.Q(i, j)) <= 1e-12);
        }
}

TEST_CASE("ctmc symbolic: the stationary solve is refused when no backend answers") {
    qn::Network<double> m = acc_mm1k(0.5, 1.0, 2);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt;

    // Solving pi Q = 0 symbolically needs cancellation in a rational function
    // field, which is delegated to whatever `api/sym` resolves. `none` keeps the
    // resolution LOCAL, so this case never touches the network and cannot pass
    // or fail on whether a line-sage-rest happens to be up on this box.
    ctmc::CtmcSymbolicOptions symopt;
    symopt.backend = "none";
    CHECK_THROWS_AS(ctmc::ctmc_symbolic_solution(sn, opt, symopt), sym::SymEngineError);
}
