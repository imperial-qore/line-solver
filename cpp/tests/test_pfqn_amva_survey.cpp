/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The six approximate MVA algorithms of Chapter 2 of H. Wang, "Approximate MVA
 * Algorithms for Solving Queueing Network Models", M.Sc. thesis, University of
 * Toronto, 1997, that LINE lacked until 2026-08-03: Bard LCP, Chow SA, the
 * Hsieh-Lam PAM family, Eager Looping, dSeS-Lavenberg-Muntz Clustering and the
 * dSeS-Muntz Improved Linearizer.
 *
 * The reference values are the MATLAB, JAR and Python ports on the same input;
 * all four agree to the printed digits. Two of the checks pin PROPERTIES rather
 * than numbers, and those are the ones that catch a bad transcription:
 *
 *  - dmlin MUST equal linearizer. IL is a cost reduction, not an
 *    approximation; survey eq. (2.50) drops the Delta^(j) correction and
 *    breaks it.
 *  - looping MUST bracket the exact solution, and lcp MUST be worse than bs,
 *    which is the accuracy ordering of Table 2.1.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_chow.h"
#include "line/api/pfqn/pfqn_clust.h"
#include "line/api/pfqn/pfqn_dmlin.h"
#include "line/api/pfqn/pfqn_lcp.h"
#include "line/api/pfqn/pfqn_linearizer.h"
#include "line/api/pfqn/pfqn_looping.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_pam.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

Matrix<double> demands() {
    Matrix<double> L(3, 2, 0.0);
    L(0, 0) = 0.10; L(0, 1) = 0.30;
    L(1, 0) = 0.20; L(1, 1) = 0.05;
    L(2, 0) = 0.40; L(2, 1) = 0.15;
    return L;
}

const std::vector<double> Nd{5.0, 4.0};
const std::vector<int> Ni{5, 4};
const std::vector<double> Zv{1.0, 2.0};

Matrix<double> zmat() {
    Matrix<double> Z(1, 2, 0.0);
    Z(0, 0) = 1.0;
    Z(0, 1) = 2.0;
    return Z;
}

/// max relative throughput error against exact MVA
double xerr(const std::vector<double>& X, const std::vector<double>& Xex) {
    double e = 0.0;
    for (std::size_t c = 0; c < X.size(); ++c)
        e = std::max(e, std::fabs(X[c] - Xex[c]) / Xex[c]);
    return e;
}

}  // namespace

TEST_CASE("pfqn_lcp: Bard large customer population") {
    const AmvaResult<double> r = pfqn_lcp(demands(), Nd, Zv);
    CHECK(r.XN[0] == doctest::Approx(1.502598).epsilon(1e-6));
    CHECK(r.XN[1] == doctest::Approx(1.188360).epsilon(1e-6));

    // Table 2.1 ranks LCP below the Bard-Schweitzer PE algorithm it seeded,
    // because dropping the (N-1)/N factor can only over-count the queue the
    // arriving customer sees.
    const auto ex = pfqn_mva(demands(), Ni, zmat());
    const AmvaResult<double> bs = pfqn_bs(demands(), Nd, Zv);
    CHECK(xerr(r.XN, ex.XN) > xerr(bs.XN, ex.XN));
}

TEST_CASE("pfqn_chow: second approximation off the LCP solution") {
    const AmvaResult<double> r = pfqn_chow(demands(), Nd, Zv);
    CHECK(r.XN[0] == doctest::Approx(1.721975).epsilon(1e-6));
    CHECK(r.XN[1] == doctest::Approx(1.250086).epsilon(1e-6));

    // Rank 3 against the PE algorithm's rank 4: the theta-correction is what
    // buys the improvement, so a theta of zero would show up here.
    const auto ex = pfqn_mva(demands(), Ni, zmat());
    const AmvaResult<double> bs = pfqn_bs(demands(), Nd, Zv);
    CHECK(xerr(r.XN, ex.XN) < xerr(bs.XN, ex.XN));

    // The backward estimator of eq. (2.16) is a different, admissible answer.
    const AmvaResult<double> b =
        pfqn_chow(demands(), Nd, Zv, std::vector<AmvaSched>(), 1e-6, 1000, Matrix<double>(),
                  ChowVariant::Backward);
    CHECK(b.XN[0] != doctest::Approx(r.XN[0]).epsilon(1e-9));
}

TEST_CASE("pfqn_pam: the three noniterative proportional approximations") {
    const AmvaResult<double> b = pfqn_pam(demands(), Nd, Zv, PamVariant::Basic);
    CHECK(b.XN[0] == doctest::Approx(1.351351).epsilon(1e-6));
    CHECK(b.XN[1] == doctest::Approx(1.024515).epsilon(1e-6));
    CHECK(b.iterations == 1);  // noniterative by construction

    // No centre is overloaded here, so the PAMI capping is inert and PAMB and
    // PAMI must coincide; a capping applied unconditionally would not.
    const AmvaResult<double> i = pfqn_pam(demands(), Nd, Zv, PamVariant::Improved);
    CHECK(i.XN[0] == doctest::Approx(b.XN[0]).epsilon(1e-12));
    CHECK(i.XN[1] == doctest::Approx(b.XN[1]).epsilon(1e-12));

    // PAMT unrolls one MVA step more, which moves the answer.
    const AmvaResult<double> t = pfqn_pam(demands(), Nd, Zv, PamVariant::Two);
    CHECK(t.XN[0] == doctest::Approx(1.672982).epsilon(1e-6));
    CHECK(t.XN[1] == doctest::Approx(1.197762).epsilon(1e-6));
}

TEST_CASE("pfqn_dmlin: the Improved Linearizer IS Linearizer") {
    const LinearizerResult<double> d = pfqn_dmlin(demands(), Ni, zmat());
    const LinearizerResult<double> l = pfqn_linearizer(demands(), Ni, zmat());
    // The defining claim of de Souza e Silva and Muntz (1990): same fixed
    // point, lower cost. Transcribing survey eq. (2.50) literally breaks this.
    for (std::size_t c = 0; c < d.X.size(); ++c)
        CHECK(d.X[c] == doctest::Approx(l.X[c]).epsilon(1e-9));
    for (std::size_t i = 0; i < d.Q.rows(); ++i)
        for (std::size_t c = 0; c < d.Q.cols(); ++c)
            CHECK(d.Q(i, c) == doctest::Approx(l.Q(i, c)).epsilon(1e-9));
    CHECK(d.X[0] == doctest::Approx(1.790568).epsilon(1e-6));
}

TEST_CASE("pfqn_clust: clustering approximation with the automatic decomposition") {
    const AmvaResult<double> r = pfqn_clust(demands(), Nd, Zv);
    CHECK(r.XN[0] == doctest::Approx(1.799264).epsilon(1e-6));
    CHECK(r.XN[1] == doctest::Approx(1.241524).epsilon(1e-6));

    // Linearizer inside is the only setting that buys anything: CA with the PE
    // algorithm inside every subnetwork IS global PE, per the paper.
    const AmvaResult<double> pe =
        pfqn_clust(demands(), Nd, Zv, std::vector<std::vector<std::size_t> >(),
                   std::vector<std::vector<std::size_t> >(),
                   ClustInner::ProportionalEstimation);
    CHECK(pe.XN[0] != doctest::Approx(r.XN[0]).epsilon(1e-9));
}

TEST_CASE("pfqn_looping: the bracket that seeds the multiclass PBH") {
    const LoopingBounds<double> b = pfqn_looping(demands(), Nd, Zv);
    CHECK(b.Xlo[0] == doctest::Approx(1.341502).epsilon(1e-6));
    CHECK(b.Xup[0] == doctest::Approx(2.129780).epsilon(1e-6));

    // What makes it a bound rather than an estimate.
    const auto ex = pfqn_mva(demands(), Ni, zmat());
    for (std::size_t c = 0; c < b.Xlo.size(); ++c) {
        CHECK(b.Xlo[c] <= ex.XN[c] + 1e-9);
        CHECK(b.Xup[c] >= ex.XN[c] - 1e-9);
    }
}
