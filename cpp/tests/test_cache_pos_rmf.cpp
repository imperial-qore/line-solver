/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `cache_miss_fifo_rmf` and `cache_miss_sfifo_rmf`, the position-resolved mean
 * fields, and the general access graph of `cache_miss_rmf`.
 *
 * THE ORACLES ARE MATLAB 3.0.7, TO 1e-9. These are deterministic fixed points of
 * a declared drift, not approximations of a quantity computable another way, so
 * the reference's own numbers are the only oracle that distinguishes a correct
 * port from a plausible one -- an identity like "the miss probabilities are in
 * [0,1]" is satisfied by a drift with a sign error in it. The values below were
 * produced by
 *
 *   lam = [1.0 0.6 0.4 0.3 0.2 0.1];  m = [2 1];
 *   [M,MU,MI,pi0] = cache_miss_fifo_rmf([], m, lam);
 *
 * and its sfifo / accost siblings, printed at %.10f.
 *
 * WHAT EACH CASE SEPARATES. FIFO(m) and strict FIFO(m) share every term of their
 * drift except where a demoted tail is reinserted, so a port that confused them
 * would pass any structural check; pinning both to distinct reference numbers is
 * what catches it. The access-graph cases likewise pin the general drift, whose
 * only structural difference from the linear one is that admission and promotion
 * are item-weighted.
 */

#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_miss_pos_rmf.h"

using namespace line;

namespace {

Matrix<double> pos_lambda() {
    const double w[6] = {1.0, 0.6, 0.4, 0.3, 0.2, 0.1};
    Matrix<double> lam(1, 6, 0.0);
    for (std::size_t k = 0; k < 6; ++k) lam(0, k) = w[k];
    return lam;
}

std::vector<int> pos_caps() {
    std::vector<int> m;
    m.push_back(2);
    m.push_back(1);
    return m;
}

/** A graph that splits admission between the two lists and lets a hit stay. */
std::vector<std::vector<Matrix<double> > > pos_accost() {
    Matrix<double> g(3, 3, 0.0);
    g(0, 1) = 0.5;
    g(0, 2) = 0.5;
    g(1, 1) = 0.3;
    g(1, 2) = 0.7;
    g(2, 2) = 1.0;
    return std::vector<std::vector<Matrix<double> > >(1, std::vector<Matrix<double> >(6, g));
}

void check_pi0(const std::vector<double>& got, const double (&want)[6]) {
    REQUIRE(got.size() == 6);
    for (std::size_t k = 0; k < 6; ++k)
        CHECK(got[k] == doctest::Approx(want[k]).epsilon(1e-9));
}

}  // namespace

TEST_CASE("cache_miss_fifo_rmf: the linear-chain fixed point, pinned to MATLAB") {
    const cache::CacheMissRmfResult<double> r =
        cache::cache_miss_fifo_rmf(std::vector<double>(), pos_caps(), pos_lambda());
    const double want[6] = {0.1903762427, 0.3293031823, 0.4574222220,
                            0.5478684445, 0.6634560154, 0.8115738932};
    check_pi0(r.pi0, want);
    CHECK(r.M == doctest::Approx(0.9491361666).epsilon(1e-9));
    REQUIRE(r.MU.size() == 1);
    CHECK(r.MU[0] == doctest::Approx(0.9491361666).epsilon(1e-9));
}

TEST_CASE("cache_miss_sfifo_rmf: strict reinsertion gives a DIFFERENT fixed point") {
    const cache::CacheMissRmfResult<double> r =
        cache::cache_miss_sfifo_rmf(std::vector<double>(), pos_caps(), pos_lambda());
    const double want[6] = {0.1622347661, 0.3144441223, 0.4577095313,
                            0.5571522680, 0.6800380627, 0.8284212495};
    check_pi0(r.pi0, want);
    CHECK(r.M == doctest::Approx(0.9199804699).epsilon(1e-9));
    // The whole point of the policy: it is not FIFO(m) rewritten.
    const cache::CacheMissRmfResult<double> f =
        cache::cache_miss_fifo_rmf(std::vector<double>(), pos_caps(), pos_lambda());
    CHECK(r.M != doctest::Approx(f.M).epsilon(1e-6));
}

TEST_CASE("cache_miss_rmf: a declared access graph takes the general drift") {
    const cache::CacheMissRmfResult<double> r = cache::cache_miss_rmf(
        std::vector<double>(), pos_caps(), pos_lambda(), 10000.0, pos_accost());
    const double want[6] = {0.2767257221, 0.3727839929, 0.4591136116,
                            0.5234841032, 0.6142328275, 0.7536597427};
    check_pi0(r.pi0, want);
    CHECK(r.M == doctest::Approx(1.0392993332).epsilon(1e-9));
    // The general drift has no 1/N refinement -- it is written for the chain --
    // so the port must REPORT that rather than claim a refined answer.
    CHECK_FALSE(r.refined);
}

TEST_CASE("cache_miss_fifo_rmf / sfifo: the same graph, the same two policies") {
    const cache::CacheMissRmfResult<double> f = cache::cache_miss_fifo_rmf(
        std::vector<double>(), pos_caps(), pos_lambda(), pos_accost());
    const double wf[6] = {0.2692454159, 0.3675809038, 0.4575700394,
                          0.5249003825, 0.6192823018, 0.7614209565};
    check_pi0(f.pi0, wf);
    CHECK(f.M == doctest::Approx(1.0302906448).epsilon(1e-9));

    const cache::CacheMissRmfResult<double> s = cache::cache_miss_sfifo_rmf(
        std::vector<double>(), pos_caps(), pos_lambda(), pos_accost());
    const double ws[6] = {0.2568588474, 0.3599997390, 0.4561478171,
                          0.5279053822, 0.6271265081, 0.7719617063};
    check_pi0(s.pi0, ws);
    CHECK(s.M == doctest::Approx(1.0163109045).epsilon(1e-9));
}

TEST_CASE("cache_miss_rmf: a LINEAR access graph keeps the refined path") {
    // `build_item_graphs` returns empty on the chain, and that emptiness is what
    // selects the 1/N-refined solve. Declaring the chain explicitly must
    // therefore give the SAME answer as declaring nothing -- if it did not, a
    // model that spelled out its default graph would silently change solver.
    Matrix<double> lin(3, 3, 0.0);
    lin(0, 1) = 1.0;
    lin(1, 2) = 1.0;
    lin(2, 2) = 1.0;
    const std::vector<std::vector<Matrix<double> > > accost(
        1, std::vector<Matrix<double> >(6, lin));
    const cache::CacheMissRmfResult<double> a =
        cache::cache_miss_rmf(std::vector<double>(), pos_caps(), pos_lambda());
    const cache::CacheMissRmfResult<double> b = cache::cache_miss_rmf(
        std::vector<double>(), pos_caps(), pos_lambda(), 10000.0, accost);
    CHECK(b.refined == a.refined);
    CHECK(b.M == doctest::Approx(a.M).epsilon(1e-12));
}

TEST_CASE("cache_miss_fifo_rmf: a transient starts cold and ends at the fixed point") {
    // The seeded transient is the reference's TSPAN/X0INIT path. An empty seed
    // takes the routine's own default -- popularity-ordered on the chain -- so
    // the assertion here is the one the OTHER seed makes false: starting from an
    // explicitly EMPTY cache, every item misses at t = 0.
    const std::vector<int> m = pos_caps();
    std::size_t slots = 0;
    for (std::size_t l = 0; l < m.size(); ++l) slots += static_cast<std::size_t>(m[l]);
    const std::vector<double> cold(6 * slots, 0.0);
    const cache::CacheMissRmfResult<double> tr = cache::cache_miss_fifo_rmf_transient(
        std::vector<double>(), m, pos_lambda(), 0.0, 200.0, cold);
    REQUIRE(tr.pi0_t.cols() > 1);
    for (std::size_t k = 0; k < 6; ++k)
        CHECK(tr.pi0_t(k, 0) == doctest::Approx(1.0).epsilon(1e-12));
    // And it converges to the steady fixed point the first case pinned.
    const double want[6] = {0.1903762427, 0.3293031823, 0.4574222220,
                            0.5478684445, 0.6634560154, 0.8115738932};
    for (std::size_t k = 0; k < 6; ++k)
        CHECK(tr.pi0_t(k, tr.pi0_t.cols() - 1) == doctest::Approx(want[k]).epsilon(1e-6));
}
