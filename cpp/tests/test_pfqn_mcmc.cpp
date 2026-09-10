/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Chen-O'Cinneide regularization (pfqn_mcmc).
 *
 * W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm for
 * Closed Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.
 *
 * Two kinds of assertion are made here, and only the first kind is statistical.
 *
 * The throughput checks compare the estimator against the EXACT ratio
 * G(N-e_r)/G(N) from convolution -- pfqn_ca for the single-server models,
 * pfqn_ncld with mu_i(k)=min(k,c_i) for the multiserver one -- at a tolerance
 * the measured error clears with room to spare. They are seeded, so they are
 * deterministic runs of a random algorithm rather than flaky tests.
 *
 * The queue-length check is NOT statistical and holds to machine precision:
 * every state of the regularized chain satisfies sum_i Y(i,r) = N(r), and the
 * reported Q is a weighted average of those states, so with no delay station the
 * columns of Q sum to N exactly whatever the sample path was. That invariant is
 * what catches an indexing or weighting error, which a tolerance on a noisy mean
 * cannot.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_ble.h"
#include "line/api/pfqn/pfqn_mcmc.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/api/pfqn/pfqn_ncld.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

/** Monte Carlo tolerance on the throughput ratios, well above the measured error. */
constexpr double RTOL = 0.02;

/** The 3 single-server stations plus IS station of Example 5.1 of the paper. */
Matrix<double> example51L() {
    const double mu[3] = {0.2, 0.5, 0.8};
    const int sets[4][3] = {{1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
    Matrix<double> L(3, 4, 0.0);
    for (int c = 0; c < 4; ++c)
        for (int i = 0; i < 3; ++i)
            if (sets[c][i]) L(i, c) = 1.0 / mu[i];
    return L;
}

std::vector<double> example51Z() { return std::vector<double>{0.0, 2.0, 2.0, 2.0}; }

Matrix<double> as_row(const std::vector<double>& v) {
    Matrix<double> m(1, v.size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i) m(0, i) = v[i];
    return m;
}

/** X(r) = G(N-e_r)/G(N) by convolution, the quantity pfqn_mcmc estimates. */
std::vector<double> exact_ratios(const Matrix<double>& L, const std::vector<int>& N,
                                 const std::vector<double>& Z) {
    const Matrix<double> Zm = as_row(Z);
    const double lG = pfqn_ca(L, N, Zm).lG;
    std::vector<double> X(N.size(), 0.0);
    for (std::size_t r = 0; r < N.size(); ++r) {
        std::vector<int> Nr = N;
        Nr[r] -= 1;
        X[r] = std::exp(pfqn_ca(L, Nr, Zm).lG - lG);
    }
    return X;
}

}  // namespace

TEST_CASE("pfqn_mcmc: example 5.1 matches the exact ratios") {
    const Matrix<double> L = example51L();
    const std::vector<int> N(4, 3);
    const std::vector<double> Z = example51Z();
    McRng rng(23000UL);
    const McmcResult<double> res =
        pfqn_mcmc(L, N, Z, std::vector<double>(), 200000, 30, 0.1, rng);
    const std::vector<double> exact = exact_ratios(L, N, Z);
    for (std::size_t r = 0; r < 4; ++r) {
        CHECK(res.X[r] == doctest::Approx(exact[r]).epsilon(RTOL));
        // Every interval is a real interval around the estimate, and every standard
        // error is a positive number rather than a NaN from a zero denominator.
        CHECK(res.Xse[r] > 0.0);
        CHECK(res.Xlo[r] < res.X[r]);
        CHECK(res.X[r] < res.Xhi[r]);
    }
    CHECK(res.batches == 30);
    CHECK(res.burnin == res.samples / 10);
}

TEST_CASE("pfqn_mcmc: the estimator is consistent") {
    // Tenfold more work must buy roughly sqrt(10) less error, not a fixed floor. A biased
    // estimator -- a mis-weighted holding time, a warm-up that never ends -- passes the
    // tolerance check above at one sample size and then stops improving. Comparing two
    // sizes separates noise from bias without pinning either run's value.
    const Matrix<double> L = example51L();
    const std::vector<int> N(4, 3);
    const std::vector<double> Z = example51Z();
    const std::vector<double> exact = exact_ratios(L, N, Z);
    double err[2] = {0.0, 0.0};
    const std::size_t sizes[2] = {50000, 500000};
    for (int k = 0; k < 2; ++k) {
        McRng rng(23000UL);
        const McmcResult<double> res =
            pfqn_mcmc(L, N, Z, std::vector<double>(), sizes[k], 30, 0.1, rng);
        for (std::size_t r = 0; r < 4; ++r)
            err[k] = std::max(err[k], std::fabs(res.X[r] - exact[r]) / exact[r]);
    }
    CHECK(err[1] < 0.5 * err[0]);
}

TEST_CASE("pfqn_mcmc: queue lengths conserve the population exactly") {
    // sum_i Y(i,r) = N(r) in every state, so the weighted average inherits it.
    Matrix<double> L(3, 2, 0.0);
    const double d[3][2] = {{0.6, 0.2}, {0.3, 0.5}, {0.1, 0.4}};
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 2; ++r) L(i, r) = d[i][r];
    const std::vector<int> N{4, 3};
    McRng rng(23000UL);
    const McmcResult<double> res = pfqn_mcmc(L, N, std::vector<double>(),
                                             std::vector<double>(), 50000, 30, 0.1, rng);
    for (std::size_t r = 0; r < 2; ++r) {
        double sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i) sum += res.Q(i, r);
        CHECK(sum == doctest::Approx(static_cast<double>(N[r])).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_mcmc: multiserver matches the exact load-dependent constant") {
    // The paper's own selling point (its Tables IV and V). The reference is the
    // load-dependent convolution with mu_i(k) = min(k, c_i), the same product form the
    // regularized chain samples, so any disagreement beyond the Monte Carlo error is an
    // error in the Psi_i(Y_i) = min(s_i, Y_i) rate.
    Matrix<double> L(3, 2, 0.0);
    const double d[3][2] = {{0.6, 0.2}, {0.3, 0.5}, {0.1, 0.4}};
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 2; ++r) L(i, r) = d[i][r];
    const std::vector<int> N{4, 3};
    const std::vector<double> c{2.0, 1.0, 3.0};
    const int Ntot = 7;
    Matrix<double> mu(3, Ntot, 0.0);
    for (int i = 0; i < 3; ++i)
        for (int k = 1; k <= Ntot; ++k) mu(i, k - 1) = std::min<double>(k, c[i]);
    const Matrix<double> Zm(1, 2, 0.0);
    const double lG = pfqn_ncld(L, N, Zm, mu, NcldMethod::Exact, 0.0).lG;
    std::vector<double> exact(2, 0.0);
    for (std::size_t r = 0; r < 2; ++r) {
        std::vector<int> Nr = N;
        Nr[r] -= 1;
        exact[r] = std::exp(pfqn_ncld(L, Nr, Zm, mu, NcldMethod::Exact, 0.0).lG - lG);
    }
    McRng rng(23000UL);
    const McmcResult<double> res =
        pfqn_mcmc(L, N, std::vector<double>(), c, 400000, 30, 0.1, rng);
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(res.X[r] == doctest::Approx(exact[r]).epsilon(RTOL));
        double sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i) sum += res.Q(i, r);
        CHECK(sum == doctest::Approx(static_cast<double>(N[r])).epsilon(1e-12));
    }
}

TEST_CASE("pfqn_mcmc: an open class is refused") {
    // The chain lives on the integer lattice sum_i Y(i,r) = N(r); the negative entry the
    // dispatcher uses to mark an open class has no state space at all.
    Matrix<double> L(2, 2, 0.5);
    const std::vector<int> N{-1, 1};
    McRng rng(23000UL);
    CHECK_THROWS_AS(pfqn_mcmc(L, N, std::vector<double>(), std::vector<double>(), 1000, 30,
                              0.1, rng),
                    line::InputError);
}

TEST_CASE("pfqn_mcmc: a populated class with no demand anywhere is refused") {
    Matrix<double> L(2, 2, 0.0);
    L(0, 0) = 0.6;
    L(1, 0) = 0.3;
    const std::vector<int> N{2, 1};
    McRng rng(23000UL);
    CHECK_THROWS_AS(pfqn_mcmc(L, N, std::vector<double>(), std::vector<double>(), 1000, 30,
                              0.1, rng),
                    line::InputError);
}

TEST_CASE("pfqn_mcmc: an empty network returns zeros without simulating") {
    Matrix<double> L(2, 2, 0.5);
    const std::vector<int> N{0, 0};
    McRng rng(23000UL);
    const McmcResult<double> res = pfqn_mcmc(L, N, std::vector<double>(),
                                             std::vector<double>(), 1000, 30, 0.1, rng);
    for (std::size_t r = 0; r < 2; ++r) CHECK(res.X[r] == 0.0);
    CHECK(res.batches == 0);
    CHECK(res.samples == 0);
    CHECK(res.burnin == 0);
}

TEST_CASE("pfqn_nc: the mcmc arm returns the BLE constant, not a paper result") {
    // The method estimates ratios and never forms G, so pfqn_nc supplies lG from BLE.
    // Pinning it here records that the number is deliberate: it cancels out of every mean
    // value the analyzer reports, and only getProbNormConstAggr reads it.
    const Matrix<double> L = example51L();
    const std::vector<int> N(4, 3);
    const std::vector<double> Z = example51Z();
    NcOptions nopt;
    nopt.samples = 20000;
    nopt.seed = 23000;
    const NcDispatchResult<double> res =
        pfqn_nc(std::vector<double>(4, 0.0), L, N, as_row(Z), NcMethod::Mcmc, 0.0, nopt);
    CHECK(res.method == "mcmc");
    std::vector<double> Nv(4, 3.0);
    CHECK(res.lG == doctest::Approx(pfqn_ble(L, Nv, Z).lG).epsilon(1e-9));
    // and the X/Q channel carried the simulation's mean values out
    CHECK(res.X.size() == 4);
    double xsum = 0.0;
    for (double x : res.X) xsum += x;
    CHECK(xsum > 0.0);
}
