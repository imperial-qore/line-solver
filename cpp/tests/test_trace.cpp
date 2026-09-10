/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Empirical trace statistics. The oracles are hand computations on a five
 * point trace -- every sample moment, autocovariance and class frequency
 * below can be checked with a pencil -- and the identities the estimators
 * satisfy: the lag-0 autocovariance IS the population variance, the class
 * moments sum to the class-independent moment, the class probabilities and
 * the transition frequencies sum to one, an empirical CDF is nondecreasing
 * and reaches one. All of it is a field computation, so the exact
 * instantiation must reproduce the hand values with no rounding at all; that
 * is the point of running these in Rational rather than only in double.
 */
#include <vector>

#include "doctest.h"
#include "line/api/trace/autocov.h"
#include "line/api/trace/mtrace_backward_moment.h"
#include "line/api/trace/mtrace_count.h"
#include "line/api/trace/mtrace_cov.h"
#include "line/api/trace/mtrace_cross_moment.h"
#include "line/api/trace/mtrace_forward_moment.h"
#include "line/api/trace/mtrace_iat2counts.h"
#include "line/api/trace/mtrace_joint.h"
#include "line/api/trace/mtrace_mean.h"
#include "line/api/trace/mtrace_merge.h"
#include "line/api/trace/mtrace_moment.h"
#include "line/api/trace/mtrace_moment_simple.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/api/trace/mtrace_sigma2.h"
#include "line/api/trace/mtrace_split.h"
#include "line/api/trace/mtrace_summary.h"
#include "line/api/trace/trace_acf.h"
#include "line/api/trace/trace_bicov.h"
#include "line/api/trace/trace_gamma.h"
#include "line/api/trace/trace_iat2bins.h"
#include "line/api/trace/trace_iat2counts.h"
#include "line/api/trace/trace_idc.h"
#include "line/api/trace/trace_idi.h"
#include "line/api/trace/trace_joint.h"
#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_pmf.h"
#include "line/api/trace/trace_scv.h"
#include "line/api/trace/trace_skew.h"
#include "line/api/trace/trace_summary.h"
#include "line/api/trace/trace_var.h"

using line::Rational;
using line::Real50;
using namespace line::trace;

namespace {

/** The five point trace 1,2,3,4,5, exact. */
std::vector<Rational> ramp_q() {
    std::vector<Rational> S;
    for (int k = 1; k <= 5; ++k) S.push_back(Rational(k));
    return S;
}

std::vector<double> ramp_d() {
    std::vector<double> S;
    for (int k = 1; k <= 5; ++k) S.push_back(static_cast<double>(k));
    return S;
}

/** Labels 1,2,1,2,1 on that trace. */
std::vector<int> alt_labels() {
    std::vector<int> A;
    A.push_back(1);
    A.push_back(2);
    A.push_back(1);
    A.push_back(2);
    A.push_back(1);
    return A;
}

}  // namespace

TEST_CASE("sample moments of a five point trace are exact") {
    const std::vector<Rational> S = ramp_q();
    // mean = 15/5 = 3
    CHECK(trace_mean(S) == Rational(3));
    // unbiased variance = (4+1+0+1+4)/4 = 10/4
    CHECK(trace_var(S) == Rational(5, 2));
    // population variance = 10/5
    CHECK(trace_var(S, false) == Rational(2));
    // scv = var/mean^2
    CHECK(trace_scv(S) == Rational(5, 18));
    CHECK(trace_scv(S, false) == Rational(2, 9));
    // a symmetric trace has zero skewness
    CHECK(trace_skew(ramp_d()) == doctest::Approx(0.0));
}

TEST_CASE("autocovariance at lag 0 is the population variance") {
    const std::vector<Rational> S = ramp_q();
    const std::vector<Rational> acv = autocov(S);
    REQUIRE(acv.size() == 4u);
    // The identity that defines the estimator, exactly.
    CHECK(acv[0] == trace_var(S, false));
    // Hand computed: centred trace -2,-1,0,1,2.
    // lag 1: (2+0+0+2)/4 = 1 ; lag 2: (0-1+0)/3 = -1/3
    CHECK(acv[1] == Rational(1));
    CHECK(acv[2] == Rational(-1, 3));

    std::vector<int> lags;
    lags.push_back(1);
    lags.push_back(2);
    const std::vector<Rational> rho = trace_acf(S, lags);
    REQUIRE(rho.size() == 2u);
    CHECK(rho[0] == Rational(1, 2));
    CHECK(rho[1] == Rational(-1, 6));
    // An acf coefficient never leaves [-1,1].
    for (std::size_t i = 0; i < rho.size(); ++i) {
        CHECK(rho[i] <= Rational(1));
        CHECK(rho[i] >= Rational(-1));
    }
    // Lags beyond n-2 are dropped, as in both references.
    std::vector<int> big;
    big.push_back(4);
    CHECK(trace_acf(S, big).empty());
}

TEST_CASE("the lag-1 joint moment reproduces the autocovariance") {
    const std::vector<Rational> S = ramp_q();
    std::vector<int> lag;
    lag.push_back(0);
    lag.push_back(1);
    std::vector<unsigned> order(2, 1u);
    // E[X_t X_{t+1}] over the 4 usable pairs = (2+6+12+20)/4 = 10
    const Rational jm = trace_joint(S, lag, order);
    CHECK(jm == Rational(10));
    // and centring it gives the lag-1 autocovariance exactly
    const std::vector<Rational> acv = autocov(S);
    CHECK(jm - Rational(9) == acv[1]);

    // The bicovariance grid is a set of third order joint moments.
    std::vector<int> grid;
    grid.push_back(1);
    const TraceBicovResult<Rational> bc = trace_bicov(S, grid);
    REQUIRE(bc.bicov.size() == 1u);
    // lags [1,1,1] cumulate to (1,2,3) -> shifted (0,1,2): (1*2*3 + 2*3*4 + 3*4*5)/3
    CHECK(bc.bicov[0] == Rational(6 + 24 + 60, 3));
}

TEST_CASE("counting and binning of a unit rate trace") {
    std::vector<Rational> S(5, Rational(1));
    // window of length 2 from each arrival epoch contains 2 arrivals
    const std::vector<long> C = trace_iat2counts(S, Rational(2));
    REQUIRE(C.size() == 3u);
    for (std::size_t i = 0; i < C.size(); ++i) CHECK(C[i] == 2);

    const TraceBinsResult B = trace_iat2bins(S, Rational(2));
    // bins of width 2 over the span 4: two full bins plus the trailing one
    REQUIRE(B.counts.size() == 3u);
    CHECK(B.counts[0] == 2);
    CHECK(B.counts[1] == 2);
    CHECK(B.counts[2] == 1);
    long total = 0;
    for (std::size_t i = 0; i < B.counts.size(); ++i) total += B.counts[i];
    CHECK(total == static_cast<long>(B.membership.size()));
}

TEST_CASE("index of dispersion for intervals, exactly") {
    const std::vector<Rational> S = ramp_q();
    // k = 2: the reference implementations aggregate S(t)+S(t+1) for
    // t = 1..n-k-1 into a vector of length n-k, leaving a trailing zero.
    // Sk = [3, 5, 0], mean 8/3, unbiased var 19/3, IDI = 2*(19/3)/(64/9)
    const TraceIdiResult<Rational> r = trace_idi(S, std::vector<long>(1, 2));
    REQUIRE(r.idi.size() == 1u);
    CHECK(r.idi[0] == Rational(57, 32));
    CHECK(r.support[0] == 2);
    // Without the spurious zero: Sk = [3,5], mean 4, var 2, IDI = 2*2/16
    const TraceIdiResult<Rational> r2 = trace_idi(S, std::vector<long>(1, 2), 0, true);
    CHECK(r2.idi[0] == Rational(1, 4));
    // A constant trace has no dispersion at all.
    std::vector<Rational> flat(10, Rational(1));
    const TraceIdiResult<Rational> r3 = trace_idi(flat, std::vector<long>(1, 2), 0, true);
    CHECK(r3.idi[0] == Rational(0));
    // trace_idc uses k = ceil(n/30), i.e. 1 here
    CHECK(trace_idc(S) == trace_idi(S, std::vector<long>(1, 1)).idi[0]);
}

TEST_CASE("an empirical pmf is a distribution and its cdf reaches one") {
    std::vector<int> X;
    X.push_back(1);
    X.push_back(1);
    X.push_back(2);
    X.push_back(3);
    X.push_back(3);
    X.push_back(3);
    const TracePmfResult<Rational> p = trace_pmf<Rational>(X);
    REQUIRE(p.values.size() == 3u);
    CHECK(p.values[0] == 1);
    CHECK(p.values[2] == 3);
    CHECK(p.pmf[0] == Rational(1, 3));
    CHECK(p.pmf[1] == Rational(1, 6));
    CHECK(p.pmf[2] == Rational(1, 2));
    Rational cdf(0), prev(0);
    for (std::size_t i = 0; i < p.pmf.size(); ++i) {
        CHECK(p.pmf[i] >= Rational(0));
        cdf += p.pmf[i];
        CHECK(cdf >= prev);  // nondecreasing
        prev = cdf;
    }
    CHECK(cdf == Rational(1));  // and it reaches exactly one
}

TEST_CASE("trace_summary is consistent with its components") {
    std::vector<double> S;
    for (int k = 1; k <= 6; ++k) S.push_back(static_cast<double>(k));
    const TraceSummary<double> s = trace_summary(S);
    CHECK(s.mean == doctest::Approx(3.5));
    CHECK(s.min == doctest::Approx(1.0));
    CHECK(s.max == doctest::Approx(6.0));
    CHECK(s.q50 == doctest::Approx(3.5));
    CHECK(s.iqr == doctest::Approx(s.q75 - s.q25));
    // median absolute deviation about the median of 1..6
    CHECK(s.mad == doctest::Approx(1.5));
    CHECK(s.skew == doctest::Approx(0.0));
    CHECK(s.scv == doctest::Approx(trace_scv(S)));
    CHECK(s.idc == doctest::Approx(trace_idc(S)));
    CHECK(s.idc_scv_ratio == doctest::Approx(s.idc / s.scv));
    REQUIRE(s.acf.size() == 4u);
    const std::vector<double> acv = autocov(S);
    CHECK(s.acf[0] == doctest::Approx(acv[1] / acv[0]));
}

TEST_CASE("trace_gamma fits a geometric acf on its grid") {
    // A deterministic trace with a slow ramp; the fit must land on the grid
    // and reproduce rho0 = (1 - 1/scv)/2 exactly.
    std::vector<Rational> S;
    for (int k = 1; k <= 12; ++k) S.push_back(Rational(k));
    const TraceGammaResult<Rational> g = trace_gamma(S, 1000);
    CHECK(g.rho0 == (Rational(1) - Rational(1) / trace_scv(S, false)) / Rational(2));
    CHECK(g.gamma >= Rational(990, 1000));
    CHECK(g.gamma <= Rational(999, 1000));
    CHECK(g.residuals >= Rational(0));
    // The reported residual is the one of the reported gamma.
    std::vector<int> lags;
    for (int l = 1; l <= 10; ++l) lags.push_back(l);
    const std::vector<Rational> rho = trace_acf(S, lags);
    Rational res(0);
    for (std::size_t i = 0; i < rho.size(); ++i) {
        const Rational d = rho[i] - g.rho0 * line::num_pow_int(g.gamma, static_cast<unsigned>(lags[i]));
        res += d * d;
    }
    CHECK(res == g.residuals);
}

TEST_CASE("class probabilities and transition frequencies sum to one") {
    const std::vector<int> A = alt_labels();
    const std::vector<Rational> pc = mtrace_pc<Rational>(A);
    REQUIRE(pc.size() == 2u);
    CHECK(pc[0] == Rational(3, 5));
    CHECK(pc[1] == Rational(2, 5));
    CHECK(pc[0] + pc[1] == Rational(1));

    const line::Matrix<Rational> sig = mtrace_sigma<Rational>(A);
    REQUIRE(sig.rows() == 2u);
    // the trace strictly alternates, so only the off-diagonal pairs occur
    CHECK(sig(0, 0) == Rational(0));
    CHECK(sig(0, 1) == Rational(1, 2));
    CHECK(sig(1, 0) == Rational(1, 2));
    CHECK(sig(1, 1) == Rational(0));
    Rational tot(0);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) tot += sig(i, j);
    CHECK(tot == Rational(1));

    const line::Matrix<Rational> sig2 = mtrace_sigma2<Rational>(A);
    REQUIRE(sig2.rows() == 2u);
    REQUIRE(sig2.cols() == 4u);
    Rational tot2(0);
    for (std::size_t i = 0; i < sig2.rows(); ++i)
        for (std::size_t j = 0; j < sig2.cols(); ++j) tot2 += sig2(i, j);
    CHECK(tot2 == Rational(1));
    // 1,2,1 and 2,1,2 are the only observed triples, at 1/3 and 2/3
    CHECK(sig2(0, 0 * 2 + 1) == Rational(0));
    CHECK(sig2(0, 1 * 2 + 0) == Rational(2, 3));
    CHECK(sig2(1, 0 * 2 + 1) == Rational(1, 3));
}

TEST_CASE("class moments sum to the class independent moment") {
    const std::vector<Rational> S = ramp_q();
    const std::vector<int> A = alt_labels();
    const std::vector<unsigned> o1(1, 1u);

    // Horvath variables, unnormalized: class 1 gets 1,3,5 and class 2 gets 2,4
    const line::Matrix<Rational> B = mtrace_moment(S, A, o1, false, false);
    REQUIRE(B.rows() == 2u);
    CHECK(B(0, 0) == Rational(9, 5));
    CHECK(B(1, 0) == Rational(6, 5));
    CHECK(B(0, 0) + B(1, 0) == trace_mean(S));  // the defining identity

    // Normalized: each entry becomes the class conditional mean, and the
    // class-probability weighted sum is again the overall mean
    const line::Matrix<Rational> Bn = mtrace_moment(S, A, o1, false, true);
    CHECK(Bn(0, 0) == Rational(3));
    CHECK(Bn(1, 0) == Rational(3));
    const std::vector<Rational> pc = mtrace_pc<Rational>(A);
    CHECK(Bn(0, 0) * pc[0] + Bn(1, 0) * pc[1] == trace_mean(S));

    // Buchholz variables: the interval AFTER a class-c event, over N-1 terms
    const line::Matrix<Rational> F = mtrace_moment(S, A, o1, true, false);
    CHECK(F(0, 0) == Rational(3, 2));  // (2+4)/4
    CHECK(F(1, 0) == Rational(2));     // (3+5)/4
    CHECK(F(0, 0) + F(1, 0) == Rational(7, 2));  // mean of T(2:end)

    // backward/forward are the two branches of mtrace_moment, normalized
    CHECK(mtrace_backward_moment(S, A, o1)(0, 0) == Bn(0, 0));
    CHECK(mtrace_forward_moment(S, A, o1)(1, 0) == mtrace_moment(S, A, o1, true, true)(1, 0));
}

TEST_CASE("cross moments and joint moments of a marked trace") {
    const std::vector<Rational> S = ramp_q();
    const std::vector<int> A = alt_labels();
    const MtraceCrossMomentResult<Rational> cm = mtrace_cross_moment(S, A, 1u);
    // transitions 1->2 carry T = 2 and 4, transitions 2->1 carry T = 3 and 5
    CHECK(cm.count(0, 1) == 2);
    CHECK(cm.count(1, 0) == 2);
    CHECK(cm.count(0, 0) == 0);
    CHECK(cm.mc(0, 1) == Rational(3));
    CHECK(cm.mc(1, 0) == Rational(4));
    // mtrace_moment_simple is the same function under a second name
    const MtraceCrossMomentResult<Rational> ms = mtrace_moment_simple(S, A, 1u);
    CHECK(ms.mc(0, 1) == cm.mc(0, 1));
    CHECK(ms.mc(1, 0) == cm.mc(1, 0));

    // The cross moments, weighted by the transition frequencies, rebuild the
    // mean of T(2:end).
    const line::Matrix<Rational> sig = mtrace_sigma<Rational>(A);
    Rational agg(0);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) agg += sig(i, j) * cm.mc(i, j);
    CHECK(agg == Rational(7, 2));

    // Joint moments, indexed by the raw label: interior events of class 1 are
    // at positions 3 and 5 (1-based), of which only position 3 has a successor
    const MtraceJointResult<Rational> jm = mtrace_joint(S, A, 1u, 1u);
    REQUIRE(jm.jm.size() == 2u);
    CHECK(jm.count[0] == 1);  // only position 3 is interior and of class 1
    CHECK(jm.count[1] == 2);
    CHECK(jm.jm[0] == Rational(2 * 3));               // T2*T3
    CHECK(jm.jm[1] == Rational(1 * 2 + 3 * 4, 2));    // (T1*T2 + T3*T4)/2
}

TEST_CASE("splitting and merging marked traces") {
    const std::vector<Rational> S = ramp_q();
    const std::vector<int> A = alt_labels();
    const MtraceSplitResult<Rational> sp = mtrace_split(S, A);
    REQUIRE(sp.traces.size() == 2u);
    // epochs 1,3,6,10,15; class 1 at 1,6,15 and class 2 at 3,10
    REQUIRE(sp.traces[0].size() == 3u);
    CHECK(sp.traces[0][0] == Rational(1));
    CHECK(sp.traces[0][1] == Rational(5));
    CHECK(sp.traces[0][2] == Rational(9));
    REQUIRE(sp.traces[1].size() == 2u);
    CHECK(sp.traces[1][0] == Rational(3));
    CHECK(sp.traces[1][1] == Rational(7));
    // each per-class trace sums to the epoch of that class's last event
    Rational s0(0);
    for (std::size_t i = 0; i < sp.traces[0].size(); ++i) s0 += sp.traces[0][i];
    CHECK(s0 == Rational(15));

    // Merging two deterministic streams: epochs 2,4 and 3
    std::vector<Rational> t1(2, Rational(2)), t2(1, Rational(3));
    const MtraceMergeResult<Rational> mg = mtrace_merge(t1, t2);
    REQUIRE(mg.times.size() == 3u);
    CHECK(mg.times[0] == Rational(2));
    CHECK(mg.times[1] == Rational(1));
    CHECK(mg.times[2] == Rational(1));
    CHECK(mg.labels[0] == 1);
    CHECK(mg.labels[1] == 2);
    CHECK(mg.labels[2] == 1);
    // the merged trace spans the later of the two streams
    Rational span(0);
    for (std::size_t i = 0; i < mg.times.size(); ++i) span += mg.times[i];
    CHECK(span == Rational(4));
}

TEST_CASE("per-class means, counts and count processes") {
    const std::vector<Rational> S = ramp_q();
    std::vector<int> type;  // 0-indexed types, as the kpctoolbox mtrace_mean wants
    type.push_back(0);
    type.push_back(1);
    type.push_back(0);
    type.push_back(1);
    type.push_back(0);
    const MtraceMeanResult<Rational> m = mtrace_mean(S, 3, type);
    REQUIRE(m.mean.size() == 3u);
    CHECK(m.count[0] == 3);
    CHECK(m.mean[0] == Rational(3));   // (1+3+5)/3
    CHECK(m.mean[1] == Rational(3));   // (2+4)/2
    CHECK(m.count[2] == 0);            // an absent type is undefined, not zero

    // Per-class counting process of a unit rate alternating trace
    std::vector<Rational> U(5, Rational(1));
    const std::vector<int> A = alt_labels();
    const MtraceCountsResult<Rational> cc = mtrace_iat2counts(U, A, Rational(2));
    REQUIRE(cc.counts.rows() == 3u);
    REQUIRE(cc.counts.cols() == 2u);
    for (std::size_t r = 0; r < cc.counts.rows(); ++r)
        CHECK(cc.counts(r, 0) + cc.counts(r, 1) == 2);  // two arrivals per window

    // Fixed resolution counts: every event is counted exactly once
    const MtraceCountResult<Rational> mc = mtrace_count(U, A, Rational(2));
    long tot = 0;
    for (std::size_t r = 0; r < mc.counts.rows(); ++r)
        for (std::size_t c = 0; c < mc.counts.cols(); ++c) tot += mc.counts(r, c);
    CHECK(tot == 4);  // the four events after the first arrival epoch
}

TEST_CASE("class pair covariances are symmetric and nonnegative on the diagonal") {
    const std::vector<Rational> S = ramp_q();
    const std::vector<int> A = alt_labels();
    const std::vector<std::vector<line::Matrix<Rational>>> cov = mtrace_cov(S, A);
    REQUIRE(cov.size() == 2u);
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 2; ++b) {
            CHECK(cov[a][b](0, 1) == cov[a][b](1, 0));
            CHECK(cov[a][b](0, 0) >= Rational(0));
            CHECK(cov[a][b](1, 1) >= Rational(0));
        }
    // The masked class-1 series is 1,0,3,0 with mean 1, so its unbiased
    // variance over the N-1 = 4 observations is (0+1+4+1)/3 = 2.
    CHECK(cov[0][0](0, 0) == Rational(2));
}

TEST_CASE("mtrace_summary assembles the fitting descriptors") {
    std::vector<Rational> S;
    for (int k = 1; k <= 8; ++k) S.push_back(Rational(k));
    std::vector<int> A;
    for (int k = 0; k < 8; ++k) A.push_back(1 + (k % 2));
    const MtraceSummary<Rational> s = mtrace_summary(S, A, 3);
    REQUIRE(s.M.size() == 5u);
    CHECK(s.M[0] == trace_mean(S));
    // raw second moment of 1..8 is 204/8
    CHECK(s.M[1] == Rational(204, 8));
    REQUIRE(s.acf.size() == 3u);
    CHECK(s.Pc[0] + s.Pc[1] == Rational(1));
    CHECK(s.F1.rows() == 2u);
    CHECK(s.B1.rows() == 2u);
    // the backward moments are the Horvath branch of mtrace_moment
    const std::vector<unsigned> o1(1, 1u);
    CHECK(s.B1(0, 0) == mtrace_moment(S, A, o1, false, true)(0, 0));
}

TEST_CASE("double and exact arithmetic agree to 1e-9") {
    const std::vector<double> Sd = ramp_d();
    const std::vector<Rational> Sq = ramp_q();
    CHECK(trace_mean(Sd) == doctest::Approx(static_cast<double>(trace_mean(Sq))).epsilon(1e-9));
    CHECK(trace_var(Sd) == doctest::Approx(static_cast<double>(trace_var(Sq))).epsilon(1e-9));
    CHECK(trace_scv(Sd) == doctest::Approx(static_cast<double>(trace_scv(Sq))).epsilon(1e-9));
    const std::vector<double> ad = autocov(Sd);
    const std::vector<Rational> aq = autocov(Sq);
    REQUIRE(ad.size() == aq.size());
    for (std::size_t i = 0; i < ad.size(); ++i)
        CHECK(ad[i] == doctest::Approx(static_cast<double>(aq[i])).epsilon(1e-9));

    const std::vector<int> A = alt_labels();
    const std::vector<unsigned> o2(1, 2u);
    const line::Matrix<double> Md = mtrace_moment(Sd, A, o2, false, true);
    const line::Matrix<Rational> Mq = mtrace_moment(Sq, A, o2, false, true);
    for (std::size_t i = 0; i < Md.rows(); ++i)
        CHECK(Md(i, 0) == doctest::Approx(static_cast<double>(Mq(i, 0))).epsilon(1e-9));

    // and the 50 digit backend agrees with both
    std::vector<Real50> Sr;
    for (int k = 1; k <= 5; ++k) Sr.push_back(Real50(k));
    CHECK(static_cast<double>(trace_scv(Sr)) ==
          doctest::Approx(static_cast<double>(trace_scv(Sq))).epsilon(1e-9));
    const TraceIdiResult<Real50> ir = trace_idi(Sr, std::vector<long>(1, 2));
    const TraceIdiResult<Rational> iq = trace_idi(Sq, std::vector<long>(1, 2));
    CHECK(static_cast<double>(ir.idi[0]) ==
          doctest::Approx(static_cast<double>(iq.idi[0])).epsilon(1e-9));
}
