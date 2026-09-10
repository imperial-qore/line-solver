/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Queue-Shift Approximation. The three stress networks of Sect. 5.2 of
 * Schweitzer, Serazzi and Broglia (Tools'98, LNCS 1469) are fully specified in
 * print together with the QSA errors they produce under the paper's own metric
 * eq. (17), so they pin the algorithm rather than this port: any drift in the
 * shift definition or in the extrapolation (15) moves them.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_qsa.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

Matrix<double> mat(const std::vector<std::vector<double>>& a) {
    Matrix<double> m(a.size(), a[0].size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) m(i, j) = a[i][j];
    return m;
}

std::vector<int> ivec(const std::vector<double>& n) {
    std::vector<int> v;
    for (double d : n) v.push_back(static_cast<int>(d));
    return v;
}

Matrix<double> zrow(std::size_t R) { return Matrix<double>(1, R, 0.0); }

/// eq. (17): err(Q) = max |Q_appr - Q_MVA| / K_r, err(U) = max |U_appr - U_MVA|
void eq17(const Matrix<double>& L, const std::vector<double>& N, double& eQ, double& eU) {
    const std::vector<double> Z(N.size(), 0.0);
    const auto ex = pfqn_mva(L, ivec(N), zrow(N.size()));
    const auto q = pfqn_qsa(L, N, Z);
    eQ = 0.0;
    eU = 0.0;
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r) {
            eQ = std::max(eQ, std::fabs(q.QN(i, r) - ex.QN(i, r)) / N[r]);
            eU = std::max(eU, std::fabs(q.UN(i, r) - ex.UN(i, r)));
        }
}

}  // namespace

TEST_CASE("QSA reproduces the published stress case of Table 2") {
    // K = (21,29) on 5 queueing stations, the case with the largest Linearizer
    // error. Published: QSA err(Q) = 0.001812, err(U) = 0.000639.
    const Matrix<double> L =
        mat({{16, 50}, {19, 28}, {0, 42}, {12, 25}, {39, 29}});
    const std::vector<double> N{21.0, 29.0};
    double eQ = 0.0, eU = 0.0;
    eq17(L, N, eQ, eU);
    CHECK(eQ == doctest::Approx(0.001812).epsilon(1e-3));
    CHECK(eU == doctest::Approx(0.000639).epsilon(1e-3));
}

TEST_CASE("QSA reproduces the published stress case of Table 3") {
    // K = (111,89), the case maximizing the Linearizer-to-QSA error ratio.
    // Published: err(Q) = 4.18e-10, err(U) = 5.54e-10, which is that paper's
    // own residual tolerance; here the Newton solve reaches machine precision,
    // consistent with its claim of "errors of zero (at least 6 digits)".
    const Matrix<double> L = mat({{16, 73}, {59, 15}, {6, 26}, {2, 36}, {10, 0}});
    const std::vector<double> N{111.0, 89.0};
    double eQ = 0.0, eU = 0.0;
    eq17(L, N, eQ, eU);
    CHECK(eQ < 4.18e-10);
    CHECK(eU < 5.54e-10);
}

TEST_CASE("QSA reproduces the published three-class stress case of Table 4") {
    // K = (2,2,2) on the Chandy-Neuse Example 2 loadings. Published:
    // err(Q) = 0.000185, err(U) = 0.000415. The tiny population is what makes
    // this one discriminating: the extrapolation (15) has no 1/Ksum to hide in.
    const Matrix<double> L = mat({{10, 1, 1}, {1, 10, 1}, {1, 1, 10}});
    const std::vector<double> N{2.0, 2.0, 2.0};
    double eQ = 0.0, eU = 0.0;
    eq17(L, N, eQ, eU);
    CHECK(eQ == doctest::Approx(0.000185).epsilon(1e-2));
    CHECK(eU == doctest::Approx(0.000415).epsilon(1e-2));
}

TEST_CASE("QSA conserves the population and stays non-negative") {
    const Matrix<double> L = mat({{16, 50}, {19, 28}, {0, 42}, {12, 25}, {39, 29}});
    const std::vector<double> N{21.0, 29.0}, Z{0.0, 0.0};
    const auto q = pfqn_qsa(L, N, Z);
    double sum = 0.0;
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r) {
            CHECK(q.QN(i, r) >= -1e-12);
            sum += q.QN(i, r);
        }
    CHECK(sum == doctest::Approx(50.0).epsilon(1e-12));
    // Per-class conservation, eq. (6).
    for (std::size_t r = 0; r < L.cols(); ++r) {
        double sr = 0.0;
        for (std::size_t i = 0; i < L.rows(); ++i) sr += q.QN(i, r);
        CHECK(sr == doctest::Approx(N[r]).epsilon(1e-10));
    }
}

TEST_CASE("QSA carries think time in the cycle time, not in the queues") {
    const Matrix<double> L = mat({{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}});
    const std::vector<double> N{6.0, 4.0}, Z{2.0, 1.0};
    const auto q = pfqn_qsa(L, N, Z);
    double held = 0.0;
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r) held += q.QN(i, r);
    for (std::size_t r = 0; r < L.cols(); ++r) held += q.XN[r] * Z[r];
    CHECK(held == doctest::Approx(10.0).epsilon(1e-10));
    // The approximation stays close to exact MVA on this model.
    const auto ex = pfqn_mva(L, ivec(N), [&] {
        Matrix<double> m(1, 2, 0.0);
        m(0, 0) = Z[0];
        m(0, 1) = Z[1];
        return m;
    }());
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r)
            CHECK(std::fabs(q.QN(i, r) - ex.QN(i, r)) < 2e-2);
}

TEST_CASE("A delay centre declared through type is served without a queueing term") {
    // Marking station 2 INF must give the same answer as moving its demand into
    // the think time, up to the approximation's own error on the queueing rows.
    const Matrix<double> L = mat({{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}});
    const std::vector<double> N{6.0, 4.0}, Z{0.0, 0.0};
    const std::vector<AmvaSched> type{AmvaSched::PS, AmvaSched::INF, AmvaSched::PS};
    const auto q = pfqn_qsa(L, N, Z, type);
    // The delay row holds exactly X_r * L_ir, which is what (13b) says.
    for (std::size_t r = 0; r < L.cols(); ++r)
        CHECK(q.QN(1, r) == doctest::Approx(q.XN[r] * L(1, r)).epsilon(1e-12));
    // Its residence time is the bare demand: no queueing term.
    for (std::size_t r = 0; r < L.cols(); ++r)
        CHECK(q.RN(1, r) == doctest::Approx(L(1, r)).epsilon(1e-12));
}

TEST_CASE("The two-level variant of eq. (14) is the cruder one") {
    // Sect. 3.1 gives it an error of O(1/Ksum); on Table 2 (Ksum = 50) that is
    // exactly the order observed, and it must not be mistaken for the default.
    const Matrix<double> L = mat({{16, 50}, {19, 28}, {0, 42}, {12, 25}, {39, 29}});
    const std::vector<double> N{21.0, 29.0}, Z{0.0, 0.0};
    const auto ex = pfqn_mva(L, ivec(N), zrow(2));
    const auto q2 = pfqn_qsa(L, N, Z, std::vector<AmvaSched>(), 1e-10, 100, 2);
    const auto q3 = pfqn_qsa(L, N, Z, std::vector<AmvaSched>(), 1e-10, 100, 3);
    double e2 = 0.0, e3 = 0.0;
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r) {
            e2 = std::max(e2, std::fabs(q2.QN(i, r) - ex.QN(i, r)) / N[r]);
            e3 = std::max(e3, std::fabs(q3.QN(i, r) - ex.QN(i, r)) / N[r]);
        }
    CHECK(e3 < e2);
    CHECK(e2 < 0.1);
    CHECK(e2 > 1.0 / 50.0);
}

TEST_CASE("QSA is exact where the model leaves nothing to approximate") {
    // One job per class: every arrival sees the other classes only, so the
    // arrival-instant queue length is exact and so is the shift.
    const Matrix<double> L = mat({{1.0, 2.0}, {3.0, 1.0}, {2.0, 2.0}});
    const std::vector<double> N{1.0, 1.0}, Z{0.0, 0.0};
    const auto ex = pfqn_mva(L, ivec(N), zrow(2));
    const auto q = pfqn_qsa(L, N, Z);
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < L.cols(); ++r)
            CHECK(q.QN(i, r) == doctest::Approx(ex.QN(i, r)).epsilon(1e-10));
}

TEST_CASE("An empty class contributes nothing and does not poison the rest") {
    const Matrix<double> L = mat({{1.0, 2.0}, {3.0, 1.0}, {2.0, 2.0}});
    const std::vector<double> N{3.0, 0.0}, Z{0.0, 0.0};
    const auto q = pfqn_qsa(L, N, Z);
    for (std::size_t i = 0; i < L.rows(); ++i) CHECK(q.QN(i, 1) == 0.0);
    CHECK(q.XN[1] == 0.0);
    double sum = 0.0;
    for (std::size_t i = 0; i < L.rows(); ++i) sum += q.QN(i, 0);
    CHECK(sum == doctest::Approx(3.0).epsilon(1e-10));
}

TEST_CASE("QSA beats Linearizer where the paper says it does") {
    // Sect. 5.2, Table 4: QSA has roughly half the error of Linearizer. The
    // regime matters -- the shift is carried per STATION, so the gain is for
    // class-independent service times at the queueing centres.
    const Matrix<double> L = mat({{10, 1, 1}, {1, 10, 1}, {1, 1, 10}});
    const std::vector<double> N{2.0, 2.0, 2.0};
    double eQ = 0.0, eU = 0.0;
    eq17(L, N, eQ, eU);
    CHECK(eQ < 0.000318);   // the published Linearizer error on the same model
    CHECK(eU < 0.000716);
}
