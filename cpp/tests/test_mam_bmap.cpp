/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The batch-arrival and batch-service level blocks, and the two finite-capacity
 * helpers: `solver_mam_bmap_map_1`, `solver_mam_map_bmap_1`,
 * `mam_detect_mmck` and `mam_truncate_renorm`.
 *
 * THE ORACLES ARE IDENTITIES, because the blocks are algebra and algebra can be
 * checked against itself in ways that do not depend on the code that produced
 * it:
 *
 *  1. DEGENERACY. A BMAP whose batch size is always one is a MAP, so the
 *     M/G/1-type blocks must coincide with `qbd_mapmap1_blocks`'s QBD blocks
 *     BIT FOR BIT, not approximately. The one block that must NOT coincide is
 *     the level-zero local one, and the difference must be exactly the service
 *     completion the reference folds back there.
 *  2. RATE PRESERVATION. Splitting a MAP's D1 by a batch pmf leaves the EPOCH
 *     rate alone and multiplies the CUSTOMER rate by the mean batch size. Both
 *     halves are asserted, because reporting the epoch rate as the arrival rate
 *     is exactly the mistake that turns a BMAP back into a MAP.
 *  3. CONSERVATION UNDER CLIPPING. Every boundary level of the GI/M/1-type
 *     chain must still have zero row sums once the oversized batches have been
 *     clipped onto "empty the queue". That identity is the entire justification
 *     for the clipping convention, and it fails immediately if the tail is
 *     dropped instead of lumped.
 *  4. `mam_truncate_renorm` on an M/M/1 must return the M/M/1/K law. That is a
 *     closed form AND a claim the reference's own docstring makes ("For an
 *     M/M/1 input the renormalized distribution coincides exactly with the
 *     M/M/1/K marginal"), so it pins the convention as well as the arithmetic.
 *  5. `mam_detect_mmck` must accept a hand-built M/M/c/K and reject each of the
 *     three conditions its docstring lists, one at a time.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mam/solver_mam_bmap.h"

using namespace line;
using Dd = lang::Distrib<double>;
using lang::SchedStrategy;

namespace {

template <class T>
Matrix<T> bmat(std::size_t r, std::size_t c, const std::vector<double>& v) {
    Matrix<T> m(r, c, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < r; ++i)
        for (std::size_t j = 0; j < c; ++j) m(i, j) = num_traits<T>::from_double(v[i * c + j]);
    return m;
}

/** Largest absolute entrywise difference; 0 means bit-for-bit equal. */
double maxdiff(const Matrix<double>& A, const Matrix<double>& B) {
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    double d = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            d = std::max(d, std::fabs(A(i, j) - B(i, j)));
    return d;
}

/** Largest absolute row sum, which is zero exactly for a generator. */
double worst_rowsum(const Matrix<double>& A) {
    double w = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j);
        w = std::max(w, std::fabs(s));
    }
    return w;
}

Matrix<double> plus(const Matrix<double>& A, const Matrix<double>& B) {
    Matrix<double> C = A;
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < C.cols(); ++j) C(i, j) += B(i, j);
    return C;
}

/** An MMPP2-shaped arrival MAP: two phases with different arrival rates. */
mam::Map<double> two_phase_arrival() {
    mam::Map<double> a;
    a.D0 = bmat<double>(2, 2, {-1.4, 0.2, 0.3, -0.8});
    a.D1 = bmat<double>(2, 2, {1.2, 0.0, 0.0, 0.5});
    return a;
}

/** A two-phase service MAP. */
mam::Map<double> two_phase_service() {
    mam::Map<double> s;
    s.D0 = bmat<double>(2, 2, {-3.0, 0.5, 0.1, -2.0});
    s.D1 = bmat<double>(2, 2, {2.0, 0.5, 1.4, 0.5});
    return s;
}

/** D1 split by a batch pmf: D_k = D1 p_k, which leaves the epoch rate alone. */
std::vector<Matrix<double>> split_by_pmf(const mam::Map<double>& m,
                                         const std::vector<double>& pmf) {
    std::vector<Matrix<double>> D;
    D.push_back(m.D0);
    for (double p : pmf) {
        Matrix<double> Dk = m.D1;
        for (std::size_t i = 0; i < Dk.rows(); ++i)
            for (std::size_t j = 0; j < Dk.cols(); ++j) Dk(i, j) *= p;
        D.push_back(Dk);
    }
    return D;
}

/** Source -> FCFS Queue -> Sink, two open classes, one shared service law. */
qn::Network<double> mmck_model(const std::string& name, const Dd& s1, const Dd& s2,
                               double servers, double cap) {
    qn::Network<double> m(name);
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dd::exp_rate(0.3));
    m.set_arrival(src, c2, Dd::exp_rate(0.4));
    m.set_service(q, c1, s1);
    m.set_service(q, c2, s2);
    if (servers != 1.0) m.set_number_of_servers(q, servers);
    if (cap > 0.0) m.set_capacity(q, cap);
    qn::RoutingMatrix<double> P;
    for (std::size_t r : {c1, c2}) {
        P.set(r, r, src, q, 1.0);
        P.set(r, r, q, k, 1.0);
    }
    m.link(P);
    return m;
}

/** Marked Poisson of the given total rate, the arrival shape M/M/c/K needs. */
mam::Mmap<double> poisson_mmap(double rate) {
    mam::Mmap<double> a;
    a.D0 = bmat<double>(1, 1, {-rate});
    a.D1 = bmat<double>(1, 1, {rate});
    a.Dc.assign(1, a.D1);
    return a;
}

template <class F>
void bmap_refuses(F f, const std::string& needle) {
    try {
        f();
        FAIL("expected a refusal naming: ", needle);
    } catch (const UnsupportedError& e) {
        CHECK(std::string(e.what()).find(needle) != std::string::npos);
    }
}

template <class F>
void bmap_rejects_input(F f, const std::string& needle) {
    try {
        f();
        FAIL("expected an input error naming: ", needle);
    } catch (const InputError& e) {
        CHECK(std::string(e.what()).find(needle) != std::string::npos);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// BMAP/MAP/1, the M/G/1-type blocks
// ---------------------------------------------------------------------------

TEST_CASE("mam bmap: a unit-batch BMAP gives back the MAP/MAP/1 QBD blocks") {
    const mam::Map<double> arr = two_phase_arrival();
    const mam::Map<double> svc = two_phase_service();
    // Every arrival epoch releases exactly one job, so the BMAP IS the MAP.
    const std::vector<Matrix<double>> D = split_by_pmf(arr, {1.0});
    const mam::BmapMap1Blocks<double> b = mam::solver_mam_bmap_map_1_blocks(D, svc);
    const mam::QbdMapMap1Blocks<double> q = mam::qbd_mapmap1_blocks(arr, svc);

    REQUIRE(b.K == 1);
    REQUIRE(b.ma == 2);
    REQUIRE(b.ms == 2);
    REQUIRE(b.m == 4);
    REQUIRE(b.Aup.size() == 1);
    // Bit for bit: both sides are the same Kronecker expressions, so anything
    // other than an exact match means the layout or the phase order differs.
    CHECK(maxdiff(b.A0, q.B) == 0.0);
    CHECK(maxdiff(b.A1, q.L) == 0.0);
    CHECK(maxdiff(b.Aup[0], q.F) == 0.0);

    // The level-zero block is the one place the two chains genuinely differ.
    // solver_mam_bmap_map_1.m folds the service completion back onto it, so the
    // server's phase process keeps running while the queue is empty; qbd_mapmap1
    // stops it. The identity is B0 = A1 + A0, to rounding rather than bit for
    // bit: the reference sums S0 + S1 BEFORE the Kronecker product, so the two
    // sides associate the same three terms differently.
    CHECK(maxdiff(b.B0, plus(b.A1, b.A0)) < 1e-14);
    CHECK(maxdiff(b.B0, q.Lbar) > 0.1);
    // The up-blocks are shared between the boundary and the repeating part.
    CHECK(maxdiff(b.Bup[0], b.Aup[0]) == 0.0);
}

TEST_CASE("mam bmap: the customer rate is the epoch rate times the mean batch") {
    const mam::Map<double> arr = two_phase_arrival();
    const mam::Map<double> svc = two_phase_service();
    const double epochRate = mam::map_lambda(arr);

    // Splitting D1 by a pmf cannot change D0 + sum_k D_k, so the phase process
    // and hence the epoch rate are untouched; only the customer count changes.
    const std::vector<double> pmf = {0.5, 0.3, 0.2};
    double meanBatch = 0.0;
    for (std::size_t k = 0; k < pmf.size(); ++k)
        meanBatch += static_cast<double>(k + 1) * pmf[k];
    CHECK(meanBatch == doctest::Approx(1.7));

    const std::vector<Matrix<double>> D = split_by_pmf(arr, pmf);
    const mam::BmapMap1Blocks<double> b = mam::solver_mam_bmap_map_1_blocks(D, svc);
    REQUIRE(b.K == 3);
    CHECK(b.lambda == doctest::Approx(epochRate * meanBatch).epsilon(1e-12));
    // The service side is an ordinary MAP, so its rate is its own lambda.
    CHECK(b.mu == doctest::Approx(mam::map_lambda(svc)).epsilon(1e-12));
    CHECK(b.rho == doctest::Approx(b.lambda / b.mu).epsilon(1e-12));

    // The degenerate pmf recovers the MAP's own rate, which is the statement
    // that a unit batch adds nothing.
    const mam::BmapMap1Blocks<double> b1 =
        mam::solver_mam_bmap_map_1_blocks(split_by_pmf(arr, {1.0}), svc);
    CHECK(b1.lambda == doctest::Approx(epochRate).epsilon(1e-12));
}

TEST_CASE("mam bmap: the M/G/1-type blocks form a generator at both levels") {
    const mam::Map<double> svc = two_phase_service();
    const std::vector<Matrix<double>> D = split_by_pmf(two_phase_arrival(), {0.5, 0.3, 0.2});
    const mam::BmapMap1Blocks<double> b = mam::solver_mam_bmap_map_1_blocks(D, svc);

    // Repeating level: down, local and every up-block together.
    Matrix<double> rep = plus(b.A0, b.A1);
    for (const Matrix<double>& U : b.Aup) rep = plus(rep, U);
    CHECK(worst_rowsum(rep) < 1e-12);

    // Level zero: the local block already carries the folded-back completion,
    // so it must close WITHOUT A0 being added a second time.
    Matrix<double> lvl0 = b.B0;
    for (const Matrix<double>& U : b.Bup) lvl0 = plus(lvl0, U);
    CHECK(worst_rowsum(lvl0) < 1e-12);
}

// ---------------------------------------------------------------------------
// MAP/BMAP/1, the GI/M/1-type blocks and the clipping convention
// ---------------------------------------------------------------------------

TEST_CASE("mam bmap: a unit-batch service BMAP gives back the MAP/MAP/1 blocks") {
    const mam::Map<double> arr = two_phase_arrival();
    const mam::Map<double> svc = two_phase_service();
    const std::vector<Matrix<double>> D = split_by_pmf(svc, {1.0});
    const mam::MapBmap1Blocks<double> b = mam::solver_mam_map_bmap_1_blocks(arr, D);
    const mam::QbdMapMap1Blocks<double> q = mam::qbd_mapmap1_blocks(arr, svc);

    REQUIRE(b.K == 1);
    CHECK(maxdiff(b.A0, q.F) == 0.0);
    CHECK(maxdiff(b.A1, q.L) == 0.0);
    CHECK(maxdiff(b.Adown[0], q.B) == 0.0);
    // With a single batch size the only boundary transition to level 0 is the
    // ordinary service completion.
    REQUIRE(b.Bto0.size() == 1);
    CHECK(maxdiff(b.Bto0[0], b.Adown[0]) == 0.0);
    // Level zero folds the service back, exactly as the arrival side does.
    CHECK(maxdiff(b.B1, plus(b.A1, b.Adown[0])) == 0.0);
}

TEST_CASE("mam bmap: the batch-service boundary lumps the tail and loses nothing") {
    const mam::Map<double> arr = two_phase_arrival();
    // A service BMAP that can clear up to three jobs at once.
    const std::vector<Matrix<double>> D = split_by_pmf(two_phase_service(), {0.5, 0.3, 0.2});
    const mam::MapBmap1Blocks<double> b = mam::solver_mam_map_bmap_1_blocks(arr, D);
    REQUIRE(b.K == 3);
    REQUIRE(b.Bto0.size() == 3);

    // The repeating level closes.
    Matrix<double> rep = plus(b.A0, b.A1);
    for (const Matrix<double>& Dn : b.Adown) rep = plus(rep, Dn);
    CHECK(worst_rowsum(rep) < 1e-12);

    // EVERY boundary level closes too, which is what "lumped, not dropped"
    // means: from level j the batches of size k < j go to level j-k and every
    // batch of size k >= j is clipped onto level 0, so the total outflow is
    // still the repeating one. Dropping the tail instead would leave these row
    // sums strictly negative.
    for (std::size_t j = 1; j <= b.K; ++j) {
        Matrix<double> lvl = plus(b.A0, b.A1);
        for (std::size_t k = 1; k < j; ++k) lvl = plus(lvl, b.Adown[k - 1]);
        lvl = plus(lvl, b.Bto0[j - 1]);
        CHECK(worst_rowsum(lvl) < 1e-12);
    }
    // Level 0 closes with the arrival block alone, the service having been
    // folded into B1.
    CHECK(worst_rowsum(plus(b.B1, b.A0)) < 1e-12);

    // The clipped blocks are nested: peeling level j off Bto0 leaves Bto0 for
    // level j+1, so no batch is counted twice and none is lost.
    for (std::size_t j = 1; j < b.K; ++j) {
        Matrix<double> diff = b.Bto0[j - 1];
        for (std::size_t i = 0; i < diff.rows(); ++i)
            for (std::size_t c = 0; c < diff.cols(); ++c)
                diff(i, c) -= b.Bto0[j](i, c) + b.Adown[j - 1](i, c);
        CHECK(worst_rowsum(diff) < 1e-12);
        CHECK(maxdiff(diff, Matrix<double>(diff.rows(), diff.cols(), 0.0)) < 1e-12);
    }
    // From level 1 every batch empties the queue, so Bto0[0] carries the whole
    // service mass; that is the extreme case of the lumping.
    Matrix<double> all = b.Adown[0];
    for (std::size_t k = 1; k < b.K; ++k) all = plus(all, b.Adown[k]);
    CHECK(maxdiff(b.Bto0[0], all) == 0.0);
    // From level K only the largest batch is clipped, so nothing is added.
    CHECK(maxdiff(b.Bto0[b.K - 1], b.Adown[b.K - 1]) == 0.0);
}

TEST_CASE("mam bmap: the batch service rate counts customers, not epochs") {
    const mam::Map<double> arr = two_phase_arrival();
    const mam::Map<double> svc = two_phase_service();
    const std::vector<double> pmf = {0.2, 0.5, 0.3};
    double meanBatch = 0.0;
    for (std::size_t k = 0; k < pmf.size(); ++k)
        meanBatch += static_cast<double>(k + 1) * pmf[k];
    CHECK(meanBatch == doctest::Approx(2.1));

    const mam::MapBmap1Blocks<double> b =
        mam::solver_mam_map_bmap_1_blocks(arr, split_by_pmf(svc, pmf));
    CHECK(b.lambda == doctest::Approx(mam::map_lambda(arr)).epsilon(1e-12));
    // Batch service serves 2.1 customers per completion epoch, so the queue is
    // stable at a load that a single-service MAP of the same epoch rate is not.
    CHECK(b.mu == doctest::Approx(mam::map_lambda(svc) * meanBatch).epsilon(1e-12));
    CHECK(b.rho == doctest::Approx(b.lambda / b.mu).epsilon(1e-12));
    CHECK(b.stable);

    const mam::MapBmap1Blocks<double> b1 =
        mam::solver_mam_map_bmap_1_blocks(arr, split_by_pmf(svc, {1.0}));
    CHECK(b1.rho == doctest::Approx(b.rho * meanBatch).epsilon(1e-12));
}

// ---------------------------------------------------------------------------
// mam_truncate_renorm
// ---------------------------------------------------------------------------

TEST_CASE("mam bmap: truncate_renorm on an M/M/1 is exactly the M/M/1/K law") {
    // The reference's docstring claims this outright: "For an M/M/1 input the
    // renormalized distribution coincides exactly with the M/M/1/K marginal."
    // p_n = rho^n / sum_{j=0..K} rho^j, and the blocking probability is p_K.
    const double lambda = 0.6, mu = 1.0;
    const std::size_t capK = 5;
    const mam::Mmap<double> arv = poisson_mmap(lambda);
    std::vector<mam::PhService<double>> svc(1);
    svc[0].sigma = {1.0};
    svc[0].S = bmat<double>(1, 1, {-mu});

    const mam::basic_detail::TruncRenorm<double> r = mam::mam_truncate_renorm(arv, svc, capK);
    REQUIRE(r.p.size() == capK + 1);

    const double rho = lambda / mu;
    double norm = 0.0;
    for (std::size_t n = 0; n <= capK; ++n) norm += std::pow(rho, static_cast<double>(n));
    double mass = 0.0, meanN = 0.0;
    for (std::size_t n = 0; n <= capK; ++n) {
        const double pn = std::pow(rho, static_cast<double>(n)) / norm;
        CHECK(r.p[n] == doctest::Approx(pn).epsilon(1e-9));
        mass += r.p[n];
        meanN += static_cast<double>(n) * pn;
    }
    // Renormalization is the point of the routine, so the mass is an identity.
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.meanQ == doctest::Approx(meanN).epsilon(1e-9));
    CHECK(r.lossProb == doctest::Approx(r.p[capK]).epsilon(1e-12));

    // THE TAIL IS REDISTRIBUTED, NOT PILED ON. Lumping, the convention the
    // batch boundary uses, would put the whole geometric tail rho^capK on the
    // top level; renormalizing instead spreads it over every level, so the
    // reported blocking probability is strictly smaller.
    const double lumped = std::pow(rho, static_cast<double>(capK));
    CHECK(r.lossProb < lumped);
    CHECK(r.lossProb == doctest::Approx(lumped / norm).epsilon(1e-9));
}

TEST_CASE("mam bmap: truncate_renorm keeps the mean inside the buffer") {
    // A heavily loaded queue whose infinite-buffer mean far exceeds the buffer:
    // the clip to [0, capK] is what stops the routine reporting more jobs than
    // the station can hold.
    const mam::Mmap<double> arv = poisson_mmap(0.95);
    std::vector<mam::PhService<double>> svc(1);
    svc[0].sigma = {1.0};
    svc[0].S = bmat<double>(1, 1, {-1.0});
    for (std::size_t capK : {1u, 2u, 6u}) {
        const mam::basic_detail::TruncRenorm<double> r = mam::mam_truncate_renorm(arv, svc, capK);
        CHECK(r.meanQ >= 0.0);
        CHECK(r.meanQ <= static_cast<double>(capK) + 1e-12);
        double mass = 0.0;
        for (double p : r.p) mass += p;
        CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    }
}

// ---------------------------------------------------------------------------
// mam_detect_mmck
// ---------------------------------------------------------------------------

TEST_CASE("mam bmap: detect_mmck accepts an M/M/c/K and rejects each departure") {
    // The docstring's three conditions: single-phase arrivals, Exp service in
    // every active class, and one shared rate. Each is broken on its own below.
    qn::Network<double> ok = mmck_model("mmckA", Dd::exp_rate(2.0), Dd::exp_rate(2.0), 2.0, 6.0);
    const mam::MmckDetection<double> d = mam::mam_detect_mmck(ok.get_struct(), 2, poisson_mmap(0.7));
    CHECK(d.isMmck);
    CHECK(d.muRate == doctest::Approx(2.0).epsilon(1e-12));

    // 1. A two-phase arrival stream is not Poisson, so the birth-death chain
    //    would answer a different arrival process.
    mam::Mmap<double> mmpp;
    mmpp.D0 = bmat<double>(2, 2, {-1.4, 0.2, 0.3, -0.8});
    mmpp.D1 = bmat<double>(2, 2, {1.2, 0.0, 0.0, 0.5});
    mmpp.Dc.assign(1, mmpp.D1);
    CHECK_FALSE(mam::mam_detect_mmck(ok.get_struct(), 2, mmpp).isMmck);

    // 2. A non-exponential service in ANY active class.
    qn::Network<double> erl =
        mmck_model("mmckB", Dd::exp_rate(2.0), Dd::erlang_fit(0.5, 0.5), 2.0, 6.0);
    CHECK_FALSE(mam::mam_detect_mmck(erl.get_struct(), 2, poisson_mmap(0.7)).isMmck);

    // 3. Per-class rates that differ leave the aggregate service
    //    non-exponential even though each class is exponential.
    qn::Network<double> mixed =
        mmck_model("mmckC", Dd::exp_rate(2.0), Dd::exp_rate(3.0), 2.0, 6.0);
    CHECK_FALSE(mam::mam_detect_mmck(mixed.get_struct(), 2, poisson_mmap(0.7)).isMmck);

    // Rates equal to within the reference's 1e-9 relative slack still pass, so
    // the test is a tolerance and not an exact bit comparison.
    qn::Network<double> near =
        mmck_model("mmckD", Dd::exp_rate(2.0), Dd::exp_rate(2.0 + 1e-11), 2.0, 6.0);
    CHECK(mam::mam_detect_mmck(near.get_struct(), 2, poisson_mmap(0.7)).isMmck);

    // The capacity plays no part in the detection: it is the caller that knows
    // the buffer, and a station with none is still an M/M/c.
    qn::Network<double> nocap =
        mmck_model("mmckE", Dd::exp_rate(2.0), Dd::exp_rate(2.0), 2.0, -1.0);
    CHECK(mam::mam_detect_mmck(nocap.get_struct(), 2, poisson_mmap(0.7)).isMmck);
}

// ---------------------------------------------------------------------------
// Refusals and input errors
// ---------------------------------------------------------------------------

TEST_CASE("mam bmap: the M/G/1-type ETAQA measures reproduce their closed forms") {
    // A unit-batch BMAP with Poisson arrivals and exponential service is an
    // M/M/1, whose level is geometric: E[N] = rho/(1-rho) and the higher
    // moments follow, so all three are exact numbers and not a golden.
    std::vector<Matrix<double>> D;
    D.push_back(bmat<double>(1, 1, {-0.6}));
    D.push_back(bmat<double>(1, 1, {0.6}));
    mam::Map<double> svc;
    svc.D0 = bmat<double>(1, 1, {-1.0});
    svc.D1 = bmat<double>(1, 1, {1.0});

    const mam::BmapQueueResult<double> r = mam::solver_mam_bmap_map_1(D, svc);
    const double rho = 0.6;
    CHECK(r.QN == doctest::Approx(rho / (1 - rho)).epsilon(1e-12));
    CHECK(r.UN == doctest::Approx(rho).epsilon(1e-12));
    CHECK(r.TN == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(r.RN == doctest::Approx(2.5).epsilon(1e-12));
    REQUIRE(r.qlenMoments.size() == 3);
    CHECK(r.qlenMoments[1] ==
          doctest::Approx(rho * (1 + rho) / ((1 - rho) * (1 - rho))).epsilon(1e-10));
    CHECK(r.qlenMoments[2] == doctest::Approx(rho * (1 + 4 * rho + rho * rho) /
                                              std::pow(1 - rho, 3.0))
                                  .epsilon(1e-10));
    // A rank-one A0 sends MG1_EG home with G = 1 before cyclic reduction runs.
    CHECK(r.fund(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    // The aggregates are (1-rho, rho(1-rho), rho^2).
    CHECK(r.piAgg(0, 0) == doctest::Approx(0.4).epsilon(1e-10));
    CHECK(r.piAgg(0, 1) == doctest::Approx(0.24).epsilon(1e-10));
    CHECK(r.piAgg(0, 2) == doctest::Approx(0.36).epsilon(1e-10));

    // A GENUINE batch: sizes 1 and 2 with equal probability, batch rate 0.3, so
    // lambda = 0.45. M[X]/M/1 has E[N] = rho/(1-rho) + rho E[X(X-1)]/(2 E[X](1-rho)),
    // which the unit-batch case cannot distinguish and a MAP surrogate would miss.
    std::vector<Matrix<double>> Dx;
    Dx.push_back(bmat<double>(1, 1, {-0.3}));
    Dx.push_back(bmat<double>(1, 1, {0.15}));
    Dx.push_back(bmat<double>(1, 1, {0.15}));
    const mam::BmapQueueResult<double> rx = mam::solver_mam_bmap_map_1(Dx, svc);
    const double rhox = 0.45, EX = 1.5, EXX1 = 1.0;
    CHECK(rx.QN == doctest::Approx(rhox / (1 - rhox) + rhox * (EXX1 / EX) / (2 * (1 - rhox)))
                       .epsilon(1e-10));
    CHECK(rx.TN == doctest::Approx(0.45).epsilon(1e-12));
}

TEST_CASE("mam bmap: the M/G/1-type ETAQA measures match MATLAB on a MAP input") {
    // Two-phase arrivals split half/half over batch sizes 1 and 2, two-phase
    // service: four phases per level, three A blocks above A0. Values from
    // solver_mam_bmap_map_1.m, which is MAMSolver's MG1_*_ETAQA verbatim.
    const std::vector<Matrix<double>> D = split_by_pmf(two_phase_arrival(), {0.5, 0.5});
    const mam::BmapQueueResult<double> r = mam::solver_mam_bmap_map_1(D, two_phase_service());

    CHECK(r.QN == doctest::Approx(2.496123595502).epsilon(1e-10));
    CHECK(r.UN == doctest::Approx(0.610619469027).epsilon(1e-10));
    CHECK(r.RN == doctest::Approx(1.808785214132).epsilon(1e-10));
    CHECK(r.TN == doctest::Approx(1.38).epsilon(1e-12));
    REQUIRE(r.qlenMoments.size() == 3);
    CHECK(r.qlenMoments[1] == doctest::Approx(17.835504350850).epsilon(1e-9));
    CHECK(r.qlenMoments[2] == doctest::Approx(191.796380027147).epsilon(1e-9));

    // The aggregates sum to one, which is the ETAQA normalization, and the
    // first four of them are pi_0 blockwise.
    double mass = 0.0;
    for (std::size_t j = 0; j < r.piAgg.cols(); ++j) mass += r.piAgg(0, j);
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.piAgg(0, 0) == doctest::Approx(0.118941068513).epsilon(1e-9));
    CHECK(r.piAgg(0, 4) == doctest::Approx(0.052543965796).epsilon(1e-9));
    CHECK(r.piAgg(0, 8) == doctest::Approx(0.188514965691).epsilon(1e-9));

    // G is stochastic: the chain is positive recurrent, so a level is left
    // downwards with probability one.
    CHECK(worst_rowsum(r.fund) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.fund(0, 0) == doctest::Approx(0.686039381261).epsilon(1e-9));
    CHECK(r.fund(3, 3) == doctest::Approx(0.217998974669).epsilon(1e-9));
}

TEST_CASE("mam bmap: the GI/M/1-type ETAQA measures match MATLAB") {
    // Unit-batch service first: MAP/BMAP/1 degenerates to the same M/M/1, and
    // R is the scalar rho.
    mam::Map<double> arr;
    arr.D0 = bmat<double>(1, 1, {-0.6});
    arr.D1 = bmat<double>(1, 1, {0.6});
    std::vector<Matrix<double>> D;
    D.push_back(bmat<double>(1, 1, {-1.0}));
    D.push_back(bmat<double>(1, 1, {1.0}));
    const mam::BmapQueueResult<double> r = mam::solver_mam_map_bmap_1(arr, D);
    CHECK(r.QN == doctest::Approx(1.5).epsilon(1e-10));
    CHECK(r.RN == doctest::Approx(2.5).epsilon(1e-10));
    CHECK(r.fund(0, 0) == doctest::Approx(0.6).epsilon(1e-10));
    CHECK(r.piAgg(0, 0) == doctest::Approx(0.4).epsilon(1e-10));
    CHECK(r.piAgg(0, 2) == doctest::Approx(0.36).epsilon(1e-10));

    // A genuine batch service, single phase: the level still falls by one or
    // two per completion, so R is no longer rho.
    mam::Map<double> arr2;
    arr2.D0 = bmat<double>(1, 1, {-0.9});
    arr2.D1 = bmat<double>(1, 1, {0.9});
    std::vector<Matrix<double>> D2;
    D2.push_back(bmat<double>(1, 1, {-1.0}));
    D2.push_back(bmat<double>(1, 1, {0.5}));
    D2.push_back(bmat<double>(1, 1, {0.5}));
    const mam::BmapQueueResult<double> r2 = mam::solver_mam_map_bmap_1(arr2, D2);
    CHECK(r2.QN == doctest::Approx(2.0611000442).epsilon(1e-9));
    CHECK(r2.UN == doctest::Approx(0.6).epsilon(1e-12));
    CHECK(r2.fund(0, 0) == doctest::Approx(0.6733200531).epsilon(1e-9));

    // Two-phase arrivals against a two-phase batch service, the four-phase
    // case. THE MEAN QUEUE LENGTH IS NEGATIVE HERE, and it is negative in
    // MATLAB too, to twelve digits: GIM1_qlen_ETAQA initializes its
    // accumulator with the SCALAR A(3) instead of the third block (see
    // lib/smc/etaqa.h, defect 3), which corrupts the last column of the linear
    // system whenever m > 1. The port reproduces the reference rather than
    // repairing it silently, and this test pins that agreement; pi and R,
    // which the defect does not touch, are the meaningful outputs.
    const std::vector<Matrix<double>> Dsvc = split_by_pmf(two_phase_service(), {0.5, 0.5});
    const mam::BmapQueueResult<double> r3 = mam::solver_mam_map_bmap_1(two_phase_arrival(), Dsvc);
    CHECK(r3.QN == doctest::Approx(-3.434197005219).epsilon(1e-9));
    CHECK(r3.UN == doctest::Approx(0.271386430678).epsilon(1e-10));
    CHECK(r3.TN == doctest::Approx(0.92).epsilon(1e-12));

    double mass = 0.0;
    for (std::size_t j = 0; j < r3.piAgg.cols(); ++j) mass += r3.piAgg(0, j);
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r3.piAgg(0, 0) == doctest::Approx(0.227094812985).epsilon(1e-9));
    CHECK(r3.piAgg(0, 4) == doctest::Approx(0.082609049829).epsilon(1e-9));
    CHECK(r3.piAgg(0, 8) == doctest::Approx(0.050296137186).epsilon(1e-9));
    CHECK(r3.fund(0, 0) == doctest::Approx(0.319301430437).epsilon(1e-9));
    CHECK(r3.fund(3, 3) == doctest::Approx(0.187111846865).epsilon(1e-9));
}

TEST_CASE("mam bmap: the ETAQA measures refuse at any arithmetic but double") {
    // The blocks are exact at Rational; the ETAQA solve on top of them is
    // LAPACK, FFT and SVD, so it says so instead of down-converting.
    mam::Map<Rational> arr;
    arr.D0 = bmat<Rational>(1, 1, {-0.6});
    arr.D1 = bmat<Rational>(1, 1, {0.6});
    std::vector<Matrix<Rational>> D;
    D.push_back(bmat<Rational>(1, 1, {-1.0}));
    D.push_back(bmat<Rational>(1, 1, {1.0}));
    bmap_refuses([&] { mam::solver_mam_map_bmap_1(arr, D); }, "--arith double");
    bmap_refuses([&] { mam::solver_mam_bmap_map_1(D, arr); }, "double precision only");
    // and it names the entry point that DOES answer at every arithmetic
    bmap_refuses([&] { mam::solver_mam_bmap_map_1(D, arr); }, "solver_mam_bmap_map_1_blocks");
}

TEST_CASE("mam bmap: a malformed input is reported before the ETAQA refusal") {
    const mam::Map<double> svc = two_phase_service();
    // A ragged BMAP: the validation must fire first, so the caller learns the
    // real problem instead of the missing third-party solver.
    std::vector<Matrix<double>> ragged;
    ragged.push_back(bmat<double>(2, 2, {-1.0, 0.0, 0.0, -1.0}));
    ragged.push_back(bmat<double>(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1}));
    bmap_rejects_input([&] { mam::solver_mam_bmap_map_1(ragged, svc); }, "must be 2x2");

    std::vector<Matrix<double>> lone;
    lone.push_back(bmat<double>(1, 1, {-1.0}));
    bmap_rejects_input([&] { mam::solver_mam_bmap_map_1(lone, svc); }, "at least D0 and D1");

    // The GI/M/1-type side additionally validates that both processes are
    // generators, which is the reference's own 1e-10 row-sum test.
    mam::Map<double> bad;
    bad.D0 = bmat<double>(1, 1, {-1.0});
    bad.D1 = bmat<double>(1, 1, {0.5});  // rows do not sum to zero
    const std::vector<Matrix<double>> good = split_by_pmf(two_phase_service(), {1.0});
    bmap_rejects_input([&] { mam::solver_mam_map_bmap_1(bad, good); },
                       "C0 + C1 must have zero row sums");

    const mam::Map<double> arr = two_phase_arrival();
    std::vector<Matrix<double>> badsvc;
    badsvc.push_back(bmat<double>(1, 1, {-1.0}));
    badsvc.push_back(bmat<double>(1, 1, {0.25}));
    bmap_rejects_input([&] { mam::solver_mam_map_bmap_1(arr, badsvc); },
                       "must have zero row sums");

    qn::Network<double> two = mmck_model("mmckF", Dd::exp_rate(2.0), Dd::exp_rate(2.0), 1.0, -1.0);
    const qn::NetworkStruct<double>& L = two.get_struct();
    bmap_rejects_input([&] { mam::mam_detect_mmck(L, 99, poisson_mmap(0.5)); }, "out of range");
}

TEST_CASE("mam bmap: the blocks are exact at Rational, the truncation is not") {
    // The block assembly is Kronecker products and one stationary solve, so it
    // carries no tolerance and instantiates at exact arithmetic.
    mam::Map<Rational> arr;
    arr.D0 = bmat<Rational>(1, 1, {-1.0});
    arr.D1 = bmat<Rational>(1, 1, {1.0});
    mam::Map<Rational> svc;
    svc.D0 = bmat<Rational>(1, 1, {-2.0});
    svc.D1 = bmat<Rational>(1, 1, {2.0});
    std::vector<Matrix<Rational>> D;
    D.push_back(svc.D0);
    D.push_back(bmat<Rational>(1, 1, {1.0}));
    D.push_back(bmat<Rational>(1, 1, {1.0}));  // half the epochs clear two jobs

    const mam::MapBmap1Blocks<Rational> b = mam::solver_mam_map_bmap_1_blocks(arr, D);
    // lambda = 1, and the service clears 1*(1) + 2*(1) = 3 customers per unit
    // time, both exactly representable.
    CHECK(b.lambda == num_traits<Rational>::from_int(1));
    CHECK(b.mu == num_traits<Rational>::from_int(3));
    CHECK(b.rho == Rational(num_traits<Rational>::from_int(1) / num_traits<Rational>::from_int(3)));

    // The finite-buffer marginal is not exact at any arithmetic, and says so.
    const mam::Mmap<Rational> arv = [] {
        mam::Mmap<Rational> a;
        a.D0 = bmat<Rational>(1, 1, {-1.0});
        a.D1 = bmat<Rational>(1, 1, {1.0});
        a.Dc.assign(1, a.D1);
        return a;
    }();
    std::vector<mam::PhService<Rational>> ph(1);
    ph[0].sigma = {num_traits<Rational>::from_int(1)};
    ph[0].S = bmat<Rational>(1, 1, {-2.0});
    bmap_refuses([&] { mam::mam_truncate_renorm(arv, ph, 4); }, "--arith double or --arith real");
}
