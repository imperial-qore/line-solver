/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Fork-join API. The oracles are closed forms known independently of the
 * MATLAB source: the expected maximum of K i.i.d. exponentials is H_K/mu, a
 * fork-join of one branch is an ordinary M/M/1, and a k-of-n quorum at k = n
 * is an AND-join while at k = 1 it is a minimum. Every identity that must hold
 * with no rounding at all (Rmin = Rmax at K = 1, R_2 = R + S_2, the Nelson-
 * Tantawi approximation collapsing onto the exact two-way formula) is asserted
 * as an exact rational equality; the transcendental functions are pinned to
 * MATLAB values computed with matlab -singleCompThread on the same inputs.
 *
 * Two of the MATLAB reference formulas are provably wrong (fj_xmax_hyperexp
 * drops an outer binomial, fj_xmax_erlang at k = 2 drops the m = 0 term and
 * carries a spurious factor 1/2). The port reproduces them exactly, because
 * MATLAB is ground truth for a port, and the tests below both pin the ported
 * value and record how far it sits from the true expected maximum, so that
 * neither the defect nor its eventual repair can pass unnoticed.
 */
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/fj/fj_bounds.h"
#include "line/api/fj/fj_char_max.h"
#include "line/api/fj/fj_gk_bound.h"
#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_order_stat.h"
#include "line/api/fj/fj_quantile.h"
#include "line/api/fj/fj_ordstat_exp.h"
#include "line/api/fj/fj_tail_ordstat.h"
#include "line/api/sn/sn_join_droprate.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/api/fj/fj_quorum_moments.h"
#include "line/api/fj/fj_respt_2way.h"
#include "line/api/fj/fj_respt_nt.h"
#include "line/api/fj/fj_respt_varki.h"
#include "line/api/fj/fj_respt_vm.h"
#include "line/api/fj/fj_rmax.h"
#include "line/api/fj/fj_rmax_erlang.h"
#include "line/api/fj/fj_rmax_evd.h"
#include "line/api/fj/fj_sm_tput.h"
#include "line/api/fj/fj_synch_delay.h"
#include "line/api/fj/fj_xmax_2.h"
#include "line/api/fj/fj_xmax_approx.h"
#include "line/api/fj/fj_xmax_emma.h"
#include "line/api/fj/fj_xmax_erlang.h"
#include "line/api/fj/fj_xmax_exp.h"
#include "line/api/fj/fj_xmax_hyperexp.h"
#include "line/api/fj/fj_xmax_normal.h"
#include "line/api/fj/fj_xmax_pareto.h"

using line::Rational;
using namespace line::fj;

namespace {
Rational Q(long n, long d = 1) { return Rational(n, d); }
}  // namespace

TEST_CASE("harmonic numbers are exact unit fractions") {
    CHECK(fj_harmonic<Rational>(1) == Q(1));
    CHECK(fj_harmonic<Rational>(2) == Q(3, 2));
    CHECK(fj_harmonic<Rational>(3) == Q(11, 6));
    CHECK(fj_harmonic<Rational>(4) == Q(25, 12));
    CHECK(fj_harmonic<Rational>(7) == Q(363, 140));  // MATLAB: 2.5928571428571425
    CHECK(static_cast<double>(fj_harmonic<Rational>(7)) ==
          doctest::Approx(2.5928571428571425).epsilon(1e-15));
    CHECK_THROWS_AS(fj_harmonic<Rational>(0), line::InputError);
}

TEST_CASE("the expected maximum of K identical exponentials is H_K/mu, exactly") {
    // The closed-form oracle: E[max of K iid Exp(mu)] = sum_{i=1..K} 1/(i mu).
    for (unsigned K = 1; K <= 8; ++K) {
        const Rational mu = Q(7, 3);
        Rational oracle(0);
        for (unsigned i = 1; i <= K; ++i) oracle += Q(1) / (Rational(static_cast<long>(i)) * mu);
        CHECK(fj_xmax_exp<Rational>(K, mu) == oracle);
    }
    // Two-branch case, reached independently through the unequal-rate formula.
    CHECK(fj_xmax_2<Rational>(Q(7, 3), Q(7, 3)) == fj_xmax_exp<Rational>(2, Q(7, 3)));
    CHECK(fj_xmax_2<Rational>(Q(2), Q(3)) == Q(19, 30));  // MATLAB: 0.6333333333333333
    // Saturated throughput is the exact reciprocal of the expected maximum.
    CHECK(fj_sm_tput<Rational>(4, Q(2)) * fj_xmax_exp<Rational>(4, Q(2)) == Q(1));
    CHECK(fj_sm_tput<Rational>(4, Q(2)) == Q(24, 25));  // MATLAB: 0.96
}

TEST_CASE("a one-branch fork-join is an ordinary M/M/1, exactly") {
    const Rational lam = Q(1, 4), mu = Q(1);
    const Rational R = Q(1) / (mu - lam);  // 4/3
    // At K = 1 the pessimistic and optimistic bounds coincide with R.
    const FJBoundsResult<Rational> b1 = fj_bounds<Rational>(1, lam, mu);
    CHECK(b1.Rmax == R);
    CHECK(b1.Rmin == R);
    // Both response-time approximations reduce to R as well.
    CHECK(fj_respt_varki<Rational>(1, lam, mu) == R);
    CHECK(fj_respt_vm<Rational>(1, lam, mu) == R);
    CHECK(fj_rmax<Rational>(1, lam, mu) == R);
}

TEST_CASE("bounds bracket the approximations and match MATLAB") {
    const Rational lam = Q(1, 4), mu = Q(1);
    const FJBoundsResult<Rational> b = fj_bounds<Rational>(4, lam, mu);
    CHECK(b.Rmin < b.Rmax);
    CHECK(b.Rmax == Q(25, 9));  // H_4 / (mu (1-rho)) = (25/12)/(3/4)
    CHECK(static_cast<double>(b.Rmax) == doctest::Approx(2.7777777777777772).epsilon(1e-14));
    CHECK(static_cast<double>(b.Rmin) == doctest::Approx(2.5350649350649346).epsilon(1e-14));
    // Rmax is exactly H_K times the M/M/1 response time.
    CHECK(b.Rmax == fj_rmax<Rational>(4, lam, mu));
    // The Varki approximation lies inside the bracket.
    const Rational Rv = fj_respt_varki<Rational>(4, lam, mu);
    CHECK(Rv > b.Rmin);
    CHECK(Rv < b.Rmax);
    CHECK_THROWS_AS(fj_bounds<Rational>(2, Q(2), Q(1)), line::NumericError);
}

TEST_CASE("the two-way identities hold with no rounding at all") {
    const Rational lam = Q(1, 4), mu = Q(1);
    const Rational R = Q(1) / (mu - lam);
    // R_2 = R + S_2: the exact two-way response time is the branch response
    // time plus the synchronization delay.
    CHECK(fj_respt_2way<Rational>(lam, mu) == R + fj_synch_delay<Rational>(lam, mu));
    // Nelson-Tantawi collapses onto the exact formula at K = 2.
    CHECK(fj_respt_nt<Rational>(2, lam, mu) == fj_respt_2way<Rational>(lam, mu));
    CHECK(fj_respt_2way<Rational>(lam, mu) == Q(47, 24));  // MATLAB: 1.9583333333333333
    CHECK(fj_synch_delay<Rational>(lam, mu) == Q(5, 8));   // MATLAB: 0.625
    CHECK_THROWS_AS(fj_respt_nt<Rational>(1, lam, mu), line::InputError);
}

TEST_CASE("the response-time approximations agree with MATLAB to 1e-12") {
    const Rational lam = Q(1, 4), mu = Q(1);
    CHECK(static_cast<double>(fj_respt_nt<Rational>(5, lam, mu)) ==
          doctest::Approx(2.8880471380471375).epsilon(1e-12));
    CHECK(static_cast<double>(fj_respt_varki<Rational>(5, lam, mu)) ==
          doctest::Approx(2.8950178476494264).epsilon(1e-12));
    CHECK(static_cast<double>(fj_respt_vm<Rational>(6, lam, mu)) ==
          doctest::Approx(3.1005256508344758).epsilon(1e-12));
    // Same inputs in double arithmetic: the alternating sum in fj_respt_vm is
    // still well conditioned at K = 6, so the two agree to near machine
    // precision.
    CHECK(fj_respt_vm<double>(6, 0.25, 1.0) ==
          doctest::Approx(static_cast<double>(fj_respt_vm<Rational>(6, lam, mu))).epsilon(1e-12));
    CHECK(fj_respt_varki<double>(5, 0.25, 1.0) ==
          doctest::Approx(static_cast<double>(fj_respt_varki<Rational>(5, lam, mu))).epsilon(1e-12));
}

TEST_CASE("fj_respt_vm: the alternating sum is exact where double loses digits") {
    // A_K alternates with terms growing like C(K,i); by K = 40 the largest term
    // is many orders of magnitude above the result. The exact evaluation is
    // still a well-defined rational and the interpolation stays bounded
    // between the light-traffic and heavy-traffic ends.
    const Rational lam = Q(1, 100), mu = Q(1);
    const Rational R40 = fj_respt_vm<Rational>(40, lam, mu);
    CHECK(R40 > Rational(0));
    // At vanishing load the interpolation tends to H_K / (mu - lambda).
    CHECK(static_cast<double>(R40) ==
          doctest::Approx(static_cast<double>(fj_harmonic<Rational>(40) / (mu - lam))).epsilon(2e-2));
}

TEST_CASE("fj_xmax_hyperexp reproduces MATLAB, which is itself defective") {
    // Port fidelity: these are the numbers MATLAB returns.
    CHECK(static_cast<double>(fj_xmax_hyperexp<Rational>(4, Q(1, 3), Q(1), Q(2))) ==
          doctest::Approx(0.40467372134038804).epsilon(1e-12));
    CHECK(static_cast<double>(fj_xmax_hyperexp<Rational>(5, Q(2, 5), Q(2), Q(2))) ==
          doctest::Approx(0.39166666666666672).epsilon(1e-12));

    // REFERENCE DEFECT, recorded so it cannot be lost. With mu1 = mu2 = mu a
    // hyperexponential IS an exponential, so the answer must be H_K/mu. The
    // MATLAB formula omits the outer binomial C(K,n) from the inclusion-
    // exclusion expansion and instead returns the alternating harmonic sum
    // sum_n (-1)^{n+1}/(n mu). Both values are asserted here: the first is what
    // the port produces (MATLAB fidelity), the second is what it ought to be.
    const Rational mu = Q(2);
    Rational alternating(0);
    for (unsigned n = 1; n <= 5; ++n) {
        const Rational t = Q(1) / (Rational(static_cast<long>(n)) * mu);
        if (n % 2 == 1) alternating += t;
        else alternating -= t;
    }
    CHECK(fj_xmax_hyperexp<Rational>(5, Q(2, 5), mu, mu) == alternating);
    CHECK(fj_xmax_hyperexp<Rational>(5, Q(2, 5), mu, mu) != fj_xmax_exp<Rational>(5, mu));
    // The true expected maximum, from MATLAB's own quadrature of 1 - F^K.
    CHECK(static_cast<double>(fj_xmax_exp<Rational>(5, mu)) ==
          doctest::Approx(1.1416666666666679).epsilon(1e-12));
}

TEST_CASE("fj_xmax_erlang at k = 2 reproduces MATLAB, which is itself defective") {
    CHECK(static_cast<double>(fj_xmax_erlang<Rational>(4, 2, Q(2))) ==
          doctest::Approx(0.36595775462962965).epsilon(1e-12));

    // REFERENCE DEFECT. At K = 1 the expected maximum of a single sample is its
    // mean, which for Erlang-2 with per-phase rate mu is 2/mu. The MATLAB
    // closed form returns 1/(2 mu): its inner sum starts at m = 1, dropping the
    // m = 0 term of the integral of the survival function, and carries an extra
    // factor 1/2. Consistently, at K = 4 and mu = 2 it returns 0.366 where the
    // true expected maximum (MATLAB quadrature of 1 - F^4) is 1.774 -- below
    // the branch mean of 1, which no maximum can be.
    CHECK(fj_xmax_erlang<Rational>(1, 2, Q(1)) == Q(1, 2));
    CHECK(fj_xmax_erlang<Rational>(1, 2, Q(1)) != Q(2));  // the Erlang-2 mean
    CHECK(static_cast<double>(fj_xmax_erlang<Rational>(4, 2, Q(2))) < 1.0);

    // The general-k quadrature branch, which is a separate code path, is sound:
    // at K = 3, k = 3, mu = 1 it agrees with MATLAB's adaptive quadrature and
    // sits above the branch mean k/mu = 3.
    CHECK(fj_xmax_erlang<double>(3, 3, 1.0) == doctest::Approx(4.4956275720164616).epsilon(1e-8));
    CHECK(fj_xmax_erlang<double>(3, 3, 1.0) > 3.0);
    // Exact arithmetic has no quadrature to offer and says so.
    CHECK_THROWS_AS(fj_xmax_erlang<Rational>(3, 3, Q(1)), line::UnsupportedError);
}

TEST_CASE("fj_rmax_erlang: the K = 2 closed form is exact and matches MATLAB") {
    CHECK(fj_rmax_erlang<Rational>(2, 2, Q(1, 4), Q(2)) == Q(55, 32));  // MATLAB: 1.71875
    CHECK(fj_rmax_erlang<double>(2, 2, 0.25, 2.0) == doctest::Approx(1.71875).epsilon(1e-14));
    // The JAR drops the mu_R^{m+n} numerator from this correction and so
    // returns 1.93115 here; MATLAB and this port do not.
    CHECK(fj_rmax_erlang<double>(2, 2, 0.25, 2.0) < 1.8);
    // General K goes through the quadrature branch.
    CHECK(fj_rmax_erlang<double>(3, 2, 0.25, 2.0) ==
          doctest::Approx(2.0081018518518512).epsilon(1e-7));
    CHECK_THROWS_AS(fj_rmax_erlang<Rational>(3, 2, Q(1, 4), Q(2)), line::UnsupportedError);
    CHECK_THROWS_AS(fj_rmax_erlang<Rational>(2, 2, Q(2), Q(2)), line::NumericError);
}

TEST_CASE("fj_xmax_approx: the exponential family is exact, the others are not") {
    const FJXmaxApproxResult<Rational> e = fj_xmax_approx<Rational>(6, Q(1), Q(1));
    CHECK(e.GK == fj_harmonic<Rational>(6) - Q(1));
    CHECK(e.Xmax == Q(1) + e.GK);
    CHECK_THROWS_AS(fj_xmax_approx<Rational>(6, Q(1), Q(1), FJDistType::Evd),
                    line::UnsupportedError);
    const FJXmaxApproxResult<double> v = fj_xmax_approx<double>(6, 1.0, 1.0, FJDistType::Evd);
    CHECK(v.GK == doctest::Approx(1.3970291267372636).epsilon(1e-12));
    CHECK(v.Xmax == doctest::Approx(2.3970291267372636).epsilon(1e-12));
}

TEST_CASE("the transcendental fork-join approximations match MATLAB to 1e-12") {
    const FJGKBoundResult<double> g = fj_gk_bound<double>(7);
    CHECK(g.exponential == doctest::Approx(1.5928571428571425).epsilon(1e-12));
    CHECK(g.uniform == doctest::Approx(1.299038105676658).epsilon(1e-12));
    CHECK(g.evd == doctest::Approx(1.5172199187065736).epsilon(1e-12));
    CHECK(g.upper_bound == doctest::Approx(1.6641005886756874).epsilon(1e-12));
    // David's bound really does dominate the three families at this K.
    CHECK(g.upper_bound > g.exponential);
    CHECK(g.upper_bound > g.uniform);
    CHECK(g.upper_bound > g.evd);

    CHECK(fj_quantile<double>(5, 0.9) == doctest::Approx(3.859805239746545).epsilon(1e-12));
    CHECK(fj_rmax_evd<double>(5, 1.0, 0.5) == doctest::Approx(1.6274367960545366).epsilon(1e-12));
    CHECK(fj_rmax_evd<double>(5, 1.0, 0.5, true) ==
          doctest::Approx(1.4940447213027848).epsilon(1e-12));
    // At K = 1 the extreme-value correction vanishes identically.
    CHECK(fj_rmax_evd<double>(1, 2.0, 0.5) == 2.0);

    CHECK(fj_xmax_emma<double>(5, 1.0) == doctest::Approx(2.242274181360723).epsilon(1e-12));

    const FJXmaxNormalResult<double> n = fj_xmax_normal<double>(10, 0.0, 1.0);
    CHECK(n.Xmax == doctest::Approx(2.2723800666962277).epsilon(1e-12));
    CHECK(n.Vmax == doctest::Approx(0.35718983958614847).epsilon(1e-12));
    // Ordering verified against MATLAB fj_xmax_normal(10,0,1,'arnold'/'corrected'):
    // default 2.272380066696, arnold 2.145966026289, corrected 2.180696323178.
    // The Arnold variant is BELOW the default, not above it.
    CHECK(fj_xmax_normal<double>(10, 0.0, 1.0, FJNormalMethod::Arnold).Xmax <
          n.Xmax);
    CHECK(fj_xmax_normal<double>(10, 0.0, 1.0, FJNormalMethod::Arnold).Xmax ==
          doctest::Approx(2.145966026289).epsilon(1e-11));
    CHECK(fj_xmax_normal<double>(10, 0.0, 1.0, FJNormalMethod::Corrected).Xmax ==
          doctest::Approx(2.180696323178).epsilon(1e-11));
    CHECK(fj_xmax_normal<double>(10, 0.0, 1.0, FJNormalMethod::Corrected).Xmax < n.Xmax);
}

TEST_CASE("fj_char_max bounds the expected maximum and matches MATLAB") {
    const FJCharMaxResult<double> e = fj_char_max<double>(5, 1.0);
    CHECK(e.MK == doctest::Approx(2.2833333333333332).epsilon(1e-12));
    CHECK(e.mK == doctest::Approx(1.6094379124341003).epsilon(1e-12));
    // For the exponential the characteristic maximum IS the expected maximum.
    CHECK(e.MK == doctest::Approx(fj_xmax_exp<double>(5, 1.0)).epsilon(1e-14));

    const FJCharMaxResult<double> r = fj_char_max<double>(5, 3u, 1.0);
    CHECK(r.MK == doctest::Approx(5.7140441230856993).epsilon(1e-10));
    CHECK(r.mK == doctest::Approx(4.2790298601253669).epsilon(1e-10));
    // m_K is by definition the point where the survival function equals 1/K.
    CHECK(detail::erlang_survival<double>(r.mK, 3, 1.0) == doctest::Approx(0.2).epsilon(1e-10));
    // Gravey's M_K is an upper bound on the true expected maximum.
    CHECK(r.MK > fj_xmax_erlang<double>(5, 3, 1.0));
}

TEST_CASE("fj_xmax_pareto: the characteristic maximum is exact, the mean is truncated") {
    const FJXmaxParetoResult<double> p = fj_xmax_pareto<double>(5, 3.0);
    CHECK(p.MK == doctest::Approx(3.1299278400300907).epsilon(1e-12));
    CHECK(p.Xmax == doctest::Approx(2.7282987340142668).epsilon(1e-4));
    // The bound dominates the (truncated) mean, as it must.
    CHECK(p.MK > p.Xmax);
    CHECK_THROWS_AS(fj_xmax_pareto<double>(5, 2.0), line::InputError);
}

TEST_CASE("fj_order_stat: the CDF identities are exact in rational arithmetic") {
    // A rational CDF, so every order-statistic CDF stays in the field.
    const std::function<Rational(const Rational&)> F = [](const Rational& y) {
        return Rational(y / (y + Rational(1)));
    };
    const Rational y = Q(1);
    const Rational Fy = F(y);  // 1/2
    const unsigned K = 3;

    // The maximum and the minimum have their textbook forms.
    CHECK(fj_order_stat<Rational>(y, K, K, F).F_Yk == line::num_pow_int(Fy, K));
    CHECK(fj_order_stat<Rational>(y, K, K, F).F_Yk == Q(1, 8));  // MATLAB: 0.125
    CHECK(fj_order_stat<Rational>(y, 1, K, F).F_Yk ==
          Rational(Q(1) - line::num_pow_int(Rational(Q(1) - Fy), K)));
    CHECK(fj_order_stat<Rational>(y, 2, K, F).F_Yk == Q(1, 2));  // MATLAB: 0.5

    // sum_{k=1..K} F_{Y_k}(y) = K F(y): each of the K samples is counted once
    // in the tail sums. Exact, and false in double at the last bit.
    Rational total(0);
    for (unsigned k = 1; k <= K; ++k) total += fj_order_stat<Rational>(y, k, K, F).F_Yk;
    CHECK(total == Rational(static_cast<long>(K)) * Fy);

    // The order statistics are stochastically ordered.
    CHECK(fj_order_stat<Rational>(y, 1, K, F).F_Yk > fj_order_stat<Rational>(y, 2, K, F).F_Yk);
    CHECK(fj_order_stat<Rational>(y, 2, K, F).F_Yk > fj_order_stat<Rational>(y, 3, K, F).F_Yk);

    // Exact arithmetic has no quadrature, so no mean is offered.
    CHECK(fj_order_stat<Rational>(y, 2, K, F).mean_available == false);
    CHECK_THROWS_AS(fj_order_stat<Rational>(y, 0, K, F), line::InputError);
}

TEST_CASE("fj_order_stat: the expected extremes match their exponential oracles") {
    // E[max of K iid Exp(1)] = H_K, E[min] = 1/K.
    const std::function<double(const double&)> F = [](const double& t) {
        return t <= 0.0 ? 0.0 : 1.0 - std::exp(-t);
    };
    const unsigned K = 3;
    const FJOrderStatResult<double> mx = fj_order_stat<double>(1.0, K, K, F);
    CHECK(mx.mean_available);
    CHECK(mx.E_Yk == doctest::Approx(11.0 / 6.0).epsilon(1e-6));  // H_3
    const FJOrderStatResult<double> mn = fj_order_stat<double>(1.0, 1, K, F);
    CHECK(mn.E_Yk == doctest::Approx(1.0 / 3.0).epsilon(1e-6));
    // The interior order statistic sits between the two.
    const FJOrderStatResult<double> md = fj_order_stat<double>(1.0, 2, K, F);
    CHECK(md.E_Yk > mn.E_Yk);
    CHECK(md.E_Yk < mx.E_Yk);
}

TEST_CASE("fj_ordstat_exp: the exponential order statistics, exactly") {
    // E[X_(1)] = 1 / sum lambda_i, whatever the branch means.
    const std::vector<double> ri{1.0 / 3.0, 1.0 / 4.0, 1.0 / 2.0};
    CHECK(fj_ordstat_exp<double>(ri, 1) == doctest::Approx(1.0 / 9.0).epsilon(1e-12));

    // n i.i.d. Exp(1) branches: E[X_(k)] = sum_{j=n-k+1..n} 1/j.
    const std::vector<double> iid{1.0, 1.0, 1.0, 1.0};
    CHECK(fj_ordstat_exp<double>(iid, 1) == doctest::Approx(0.25).epsilon(1e-12));
    CHECK(fj_ordstat_exp<double>(iid, 2) == doctest::Approx(0.25 + 1.0 / 3.0).epsilon(1e-12));
    CHECK(fj_ordstat_exp<double>(iid, 4) ==
          doctest::Approx(0.25 + 1.0 / 3.0 + 0.5 + 1.0).epsilon(1e-12));

    // At k = n it must reproduce the classical E[max] series term for term: that is
    // what keeps a standard join bit-identical to the pre-quorum fixed point.
    const std::vector<double> het{0.7, 1.3, 2.9, 0.4};
    double expected = 0;
    for (std::size_t mask = 1; mask < (std::size_t(1) << het.size()); ++mask) {
        double s = 0;
        int bits = 0;
        for (std::size_t i = 0; i < het.size(); ++i)
            if (mask & (std::size_t(1) << i)) {
                s += 1.0 / het[i];
                ++bits;
            }
        expected += (bits % 2 == 1 ? 1.0 : -1.0) / s;
    }
    CHECK(fj_ordstat_exp<double>(het, het.size()) == doctest::Approx(expected).epsilon(1e-10));

    // Monotone in k.
    double prev = -1;
    for (std::size_t k = 1; k <= het.size(); ++k) {
        const double m = fj_ordstat_exp<double>(het, k);
        CHECK(m > prev);
        prev = m;
    }

    // A branch of zero mean completes instantly: it counts toward the quorum at once
    // and never delays the join. The reference reaches the same values through 1/Inf.
    const std::vector<double> withzero{0.0, 2.0};
    CHECK(fj_ordstat_exp<double>(withzero, 1) == doctest::Approx(0.0).epsilon(1e-14));
    CHECK(fj_ordstat_exp<double>(withzero, 2) == doctest::Approx(2.0).epsilon(1e-14));

    CHECK(fj_ordstat_exp<double>(std::vector<double>{3.0}, 1) ==
          doctest::Approx(3.0).epsilon(1e-14));
    CHECK(fj_ordstat_exp<double>(std::vector<double>{}, 1) == doctest::Approx(0.0));
    CHECK_THROWS_AS(fj_ordstat_exp<double>(ri, 0), line::InputError);
    CHECK_THROWS_AS(fj_ordstat_exp<double>(ri, 4), line::InputError);

    // The exact path is RATIONAL: no sqrt, no log, so it holds in the field.
    const std::vector<Rational> rq{Q(1, 3), Q(1, 4), Q(1, 2)};
    CHECK(fj_ordstat_exp<Rational>(rq, 1) == Q(1, 9));
    CHECK(static_cast<double>(fj_ordstat_exp<Rational>(rq, 3)) ==
          doctest::Approx(fj_ordstat_exp<double>(ri, 3)).epsilon(1e-12));
}

TEST_CASE("fj_quorum_moments: k = n is the AND-join, k = 1 is the minimum") {
    const std::vector<double> means{1.0, 2.0, 3.0};
    const std::vector<double> vars{1.0, 1.0, 1.0};
    const FJQuorumMomentsResult<double> q2 = fj_quorum_moments<double>(means, vars, 2);
    CHECK(q2.m == doctest::Approx(1.9629629629629628).epsilon(1e-12));  // MATLAB
    CHECK(q2.v == doctest::Approx(0.5727023319615907).epsilon(1e-12));

    const FJQuorumMomentsResult<double> q1of2 =
        fj_quorum_moments<double>(std::vector<double>{1.0, 2.0},
                                  std::vector<double>{0.25, 0.5}, 1);
    CHECK(q1of2.m == doctest::Approx(0.96071628993408054).epsilon(1e-12));  // MATLAB
    CHECK(q1of2.v == doctest::Approx(0.19766714776939592).epsilon(1e-12));

    // The quorum time is monotone in k: a larger quorum waits longer.
    const FJQuorumMomentsResult<double> q1 = fj_quorum_moments<double>(means, vars, 1);
    const FJQuorumMomentsResult<double> q3 = fj_quorum_moments<double>(means, vars, 3);
    CHECK(q1.m < q2.m);
    CHECK(q2.m < q3.m);
    // Identical deterministic branches: every order statistic is that constant.
    const FJQuorumMomentsResult<double> det =
        fj_quorum_moments<double>(std::vector<double>{2.0, 2.0, 2.0},
                                  std::vector<double>{0.0, 0.0, 0.0}, 2);
    CHECK(det.m == doctest::Approx(2.0).epsilon(1e-14));
    CHECK(det.v == doctest::Approx(0.0).epsilon(1e-14));
    // A single branch reproduces its own fitted mean.
    const FJQuorumMomentsResult<double> one =
        fj_quorum_moments<double>(std::vector<double>{3.0}, std::vector<double>{1.0}, 1);
    CHECK(one.m == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(one.v == doctest::Approx(1.0).epsilon(1e-12));

    CHECK_THROWS_AS(fj_quorum_moments<double>(means, vars, 4), line::InputError);
    CHECK_THROWS_AS(fj_quorum_moments<double>(means, std::vector<double>{1.0}, 1),
                    line::InputError);
}

TEST_CASE("fork-join results agree between double and exact arithmetic") {
    const Rational lam = Q(1, 4), mu = Q(1);
    CHECK(fj_harmonic<double>(9) ==
          doctest::Approx(static_cast<double>(fj_harmonic<Rational>(9))).epsilon(1e-9));
    CHECK(fj_bounds<double>(5, 0.25, 1.0).Rmin ==
          doctest::Approx(static_cast<double>(fj_bounds<Rational>(5, lam, mu).Rmin)).epsilon(1e-9));
    CHECK(fj_respt_2way<double>(0.25, 1.0) ==
          doctest::Approx(static_cast<double>(fj_respt_2way<Rational>(lam, mu))).epsilon(1e-9));
    CHECK(fj_respt_nt<double>(6, 0.25, 1.0) ==
          doctest::Approx(static_cast<double>(fj_respt_nt<Rational>(6, lam, mu))).epsilon(1e-9));
    CHECK(fj_synch_delay<double>(0.25, 1.0) ==
          doctest::Approx(static_cast<double>(fj_synch_delay<Rational>(lam, mu))).epsilon(1e-9));
    CHECK(fj_xmax_hyperexp<double>(6, 1.0 / 3.0, 1.0, 2.0) ==
          doctest::Approx(static_cast<double>(fj_xmax_hyperexp<Rational>(6, Q(1, 3), Q(1), Q(2))))
              .epsilon(1e-9));
    CHECK(fj_xmax_erlang<double>(6, 2, 2.0) ==
          doctest::Approx(static_cast<double>(fj_xmax_erlang<Rational>(6, 2, Q(2)))).epsilon(1e-9));
    CHECK(fj_rmax_erlang<double>(2, 3, 0.25, 2.0) ==
          doctest::Approx(static_cast<double>(fj_rmax_erlang<Rational>(2, 3, Q(1, 4), Q(2))))
              .epsilon(1e-9));
}

/**
 * The Join row's loss, which is the one station where LossRate = ArvR - Tput
 * does NOT hold: ArvR counts the SIBLINGS offered (N per parent job) and Tput
 * the PARENT jobs released. Reading the identity there charges (N-1)/N of the
 * offered traffic as lost at EVERY join, standard joins included.
 */
namespace {

using DD = line::lang::Distrib<double>;

/** Delay -> Fork -{Q1,Q2,Q3}- Join(k of 3) -> Delay; k = 0 declares STD. */
line::qn::Network<double> quorum_model(double k) {
    line::qn::Network<double> m("quorum");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q1 = m.add_queue("Queue1", line::lang::SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Queue2", line::lang::SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Queue3", line::lang::SchedStrategy::PS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("class1", 5, d, 0);
    m.set_service(d, c, DD::exp_rate(1.0));
    m.set_service(q1, c, DD::exp_rate(2.0));
    m.set_service(q2, c, DD::exp_rate(2.0));
    m.set_service(q3, c, DD::exp_rate(2.0));
    if (k > 0) m.set_join_strategy(j, line::lang::JoinStrategy::PARTIAL, k);
    line::qn::RoutingMatrix<double> P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, f, q3, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, j, 1.0);
    P.set(c, c, q3, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("sn_join_siblings counts the siblings at the FORK") {
    line::qn::Network<double> m = quorum_model(2.0);
    const line::qn::NetworkStruct<double>& sn = m.get_struct();
    REQUIRE(sn.fj.size() == 1);
    CHECK(sn.join_siblings(sn.fj[0].second) == 3);
    CHECK(sn.quorum_joins().size() == 1);
    // k = n is a full join and not a quorum, and neither is an absent declaration
    line::qn::Network<double> full = quorum_model(3.0);
    CHECK(full.get_struct().quorum_joins().empty());
    line::qn::Network<double> std_join = quorum_model(0.0);
    CHECK(std_join.get_struct().quorum_joins().empty());
}

TEST_CASE("sn_join_droprate: a standard join loses nothing, a quorum loses n-k") {
    const double x = 1.6;
    {
        line::qn::Network<double> m = quorum_model(0.0);
        const line::qn::NetworkStruct<double>& sn = m.get_struct();
        const std::size_t ist = sn.nodes[sn.fj[0].second - 1].station;
        REQUIRE(ist != 0);
        line::Matrix<double> TN(sn.stations.size(), sn.classes.size(), 0.0);
        line::Matrix<double> AN(sn.stations.size(), sn.classes.size(), 0.0);
        TN(ist - 1, 0) = x;
        AN(ist - 1, 0) = 3 * x;  // all three siblings offered, all three consumed
        const line::Matrix<double> d = line::sn::sn_join_droprate(sn, TN, AN);
        CHECK(d(ist - 1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    }
    {
        line::qn::Network<double> m = quorum_model(2.0);
        const line::qn::NetworkStruct<double>& sn = m.get_struct();
        const std::size_t ist = sn.nodes[sn.fj[0].second - 1].station;
        line::Matrix<double> TN(sn.stations.size(), sn.classes.size(), 0.0);
        line::Matrix<double> AN(sn.stations.size(), sn.classes.size(), 0.0);
        TN(ist - 1, 0) = x;
        AN(ist - 1, 0) = 3 * x;  // three offered, two consumed, one discarded
        const line::Matrix<double> d = line::sn::sn_join_droprate(sn, TN, AN);
        CHECK(d(ist - 1, 0) == doctest::Approx(x).epsilon(1e-12));
        // and nowhere else
        for (std::size_t i = 0; i < sn.stations.size(); ++i)
            if (i != ist - 1) CHECK(d(i, 0) == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("sn_join_droprate never reports a negative loss") {
    line::qn::Network<double> m = quorum_model(2.0);
    const line::qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t ist = sn.nodes[sn.fj[0].second - 1].station;
    line::Matrix<double> TN(sn.stations.size(), sn.classes.size(), 0.0);
    line::Matrix<double> AN(sn.stations.size(), sn.classes.size(), 0.0);
    TN(ist - 1, 0) = 5.0;  // inconsistent with AN on purpose
    AN(ist - 1, 0) = 1.0;
    CHECK(line::sn::sn_join_droprate(sn, TN, AN)(ist - 1, 0) >= 0.0);
}

TEST_CASE("a quorum join leaves the chain capacity UNBOUNDED") {
    // The stragglers of an already-fired parent are still in their branches when
    // it forks again, and nothing bounds that backlog, so a branch station holds
    // no more than the class population only under a STANDARD join. Capping it at
    // sum(njobs) made the engine refuse the model outright, dropping a closed job.
    line::qn::Network<double> mq = quorum_model(2.0);
    const line::qn::NetworkStruct<double>& snq = mq.get_struct();
    for (std::size_t i = 0; i < snq.stations.size(); ++i)
        CHECK(std::isinf(snq.classcap[i][0]));
    line::qn::Network<double> ms = quorum_model(0.0);
    const line::qn::NetworkStruct<double>& sns = ms.get_struct();
    for (std::size_t i = 0; i < sns.stations.size(); ++i)
        CHECK(sns.classcap[i][0] == doctest::Approx(5.0));
}

/**
 * `Fork.setTasksPerLink(w)` sends w IDENTICAL tasks down each of B links, so a
 * firing creates w*B siblings. The join synchronises on the order statistic of
 * that many branch times -- each branch REPLICATED w times -- and NOT on w times
 * the order statistic of B, which is w*H_B/mu where the answer is H_(w*B)/mu.
 * The reference values are SolverMVA on the same model in MATLAB, the JAR and
 * native python, which agree with this port to every printed digit.
 */
namespace {

/** Delay -> Fork(w tasks per link) -{nb PS queues}- Join -> Delay, N = 3. */
line::qn::Network<double> tpl_model(std::size_t w, std::size_t nb) {
    line::qn::Network<double> m("FJ-TPL");
    const std::size_t d = m.add_delay("Delay");
    std::vector<std::size_t> qs;
    for (std::size_t b = 0; b < nb; ++b)
        qs.push_back(m.add_queue("Queue" + std::to_string(b + 1), line::lang::SchedStrategy::PS));
    const std::size_t f = m.add_fork("Fork", static_cast<double>(w));
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("class1", 3.0, d, 0);
    m.set_service(d, c, DD::exp_rate(1.0));
    for (std::size_t b = 0; b < nb; ++b) m.set_service(qs[b], c, DD::exp_rate(2.0));
    line::qn::RoutingMatrix<double> P;
    P.set(c, c, d, f, 1.0);
    for (std::size_t b = 0; b < nb; ++b) {
        P.set(c, c, f, qs[b], 1.0);
        P.set(c, c, qs[b], j, 1.0);
    }
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

double tpl_tput(line::qn::Network<double>& m) {
    line::mva::MvaOptions opt;
    opt.method = "default";
    line::Matrix<double> init;
    return line::mva::solver_mva_run_analyzer(m.get_struct(), opt, init).TN(0, 0);
}

}  // namespace

TEST_CASE("tasksPerLink: the join synchronises on w*B siblings") {
    const double ref[2][3] = {{1.167124, 0.702141, 0.507307},
                              {1.058612, 0.662789, 0.487276}};
    for (std::size_t bi = 0; bi < 2; ++bi) {
        const std::size_t nb = bi + 2;
        for (std::size_t wi = 0; wi < 3; ++wi) {
            const std::size_t w = wi + 1;
            line::qn::Network<double> m = tpl_model(w, nb);
            CHECK(tpl_tput(m) == doctest::Approx(ref[bi][wi]).epsilon(1e-5));
        }
    }
}

TEST_CASE("tasksPerLink: a branch station holds w jobs per circulating parent") {
    for (std::size_t nb = 2; nb <= 3; ++nb)
        for (std::size_t w = 1; w <= 3; ++w) {
            line::qn::Network<double> m = tpl_model(w, nb);
            CHECK(m.get_struct().classcap[0][0] ==
                  doctest::Approx(3.0 * static_cast<double>(w)));
        }
}

TEST_CASE("tasksPerLink: replicating is not scaling the order statistic") {
    // w*E[X_(B)] is strictly larger than E[X_(w*B)] on identical branches, so the
    // two cannot both be the synchronisation instant.
    const std::vector<double> ri{0.5, 0.5};
    std::vector<double> rep = ri;
    rep.insert(rep.end(), ri.begin(), ri.end());
    const double scaled = 2.0 * fj_ordstat_exp<double>(ri, 2);
    const double replicated = fj_ordstat_exp<double>(rep, 4);
    CHECK(scaled == doctest::Approx(2.0 * 0.5 * (1.0 + 1.0 / 2.0)).epsilon(1e-12));
    CHECK(replicated ==
          doctest::Approx(0.5 * (1.0 + 1.0 / 2.0 + 1.0 / 3.0 + 1.0 / 4.0)).epsilon(1e-12));
    CHECK(replicated < scaled);
}

/**
 * `fj_tail_ordstat`: the k-of-n (QUORUM) fork-join tail. A quorum join fires on
 * the k-th of n siblings, so the request response time is the k-th ORDER
 * STATISTIC of the branch times and not their maximum; reading the maximum
 * returns the AND-join tail under a quorum's name, the same number for every k.
 * The reference values are agreed across MATLAB, the JAR and native python, and
 * were checked against a 400k-sample Monte Carlo of the SAME fitted GE branches
 * (so the check is of the inversion, not of the fit): every k agreed to within
 * 1.2% at the 99th percentile, which is the sampling error there.
 */
TEST_CASE("fj_tail_ordstat: the k-of-n quorum tail") {
    const std::vector<double> ET(1, 2.0), VT(1, 6.0);
    const double hom[4][2] = {{0.871483, 2.179109},
                              {2.146981, 4.220359},
                              {4.189345, 7.427549},
                              {8.718015, 15.050364}};
    for (std::size_t k = 1; k <= 4; ++k) {
        CHECK(line::fj::fj_tail_ordstat<double>(ET, VT, 4, 0.90, k).xp ==
              doctest::Approx(hom[k - 1][0]).epsilon(1e-5));
        CHECK(line::fj::fj_tail_ordstat<double>(ET, VT, 4, 0.99, k).xp ==
              doctest::Approx(hom[k - 1][1]).epsilon(1e-5));
    }
    const std::vector<double> ETh{1.0, 2.0, 4.0}, VTh{1.0, 8.0, 40.0};
    const double het[3][2] = {{0.980116, 2.397990}, {3.006171, 7.388476}, {12.421746, 29.881069}};
    for (std::size_t k = 1; k <= 3; ++k) {
        CHECK(line::fj::fj_tail_ordstat<double>(ETh, VTh, 1, 0.90, k).xp ==
              doctest::Approx(het[k - 1][0]).epsilon(1e-5));
        CHECK(line::fj::fj_tail_ordstat<double>(ETh, VTh, 1, 0.99, k).xp ==
              doctest::Approx(het[k - 1][1]).epsilon(1e-5));
    }
}

TEST_CASE("fj_tail_ordstat: k = n is fj_tail_forktail EXACTLY") {
    // No existing result may move: a full join must return the ForkTail root
    // itself, not a re-derivation of it.
    const std::vector<double> ET(1, 2.0), VT(1, 6.0);
    const std::vector<double> ETh{1.0, 2.0, 4.0}, VTh{1.0, 8.0, 40.0};
    const std::vector<double> kv(1, 4.0);
    const double ps[4] = {0.5, 0.9, 0.99, 0.999};
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(line::fj::fj_tail_ordstat<double>(ET, VT, 4, ps[i], 4).xp ==
              line::fj::fj_tail_forktail<double>(ET, VT, kv, ps[i]).xp);
        CHECK(line::fj::fj_tail_ordstat<double>(ETh, VTh, 1, ps[i], 3).xp ==
              line::fj::fj_tail_forktail<double>(ETh, VTh, std::vector<double>(), ps[i]).xp);
    }
}

TEST_CASE("fj_tail_ordstat: the percentile grows with the quorum") {
    // Waiting for more siblings can only take longer.
    const std::vector<double> ET(1, 2.0), VT(1, 6.0);
    const std::vector<double> ETh{1.0, 2.0, 4.0}, VTh{1.0, 8.0, 40.0};
    const double ps[3] = {0.5, 0.9, 0.99};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t k = 1; k < 4; ++k)
            CHECK(line::fj::fj_tail_ordstat<double>(ET, VT, 4, ps[i], k).xp <
                  line::fj::fj_tail_ordstat<double>(ET, VT, 4, ps[i], k + 1).xp);
        for (std::size_t k = 1; k < 3; ++k)
            CHECK(line::fj::fj_tail_ordstat<double>(ETh, VTh, 1, ps[i], k).xp <
                  line::fj::fj_tail_ordstat<double>(ETh, VTh, 1, ps[i], k + 1).xp);
    }
}

TEST_CASE("fj_tail_ordstat: an out-of-range quorum or percentile is refused") {
    const std::vector<double> ET(1, 2.0), VT(1, 6.0);
    CHECK_THROWS_AS(line::fj::fj_tail_ordstat<double>(ET, VT, 4, 0.99, 5), line::InputError);
    CHECK_THROWS_AS(line::fj::fj_tail_ordstat<double>(ET, VT, 4, 0.0, 2), line::InputError);
    CHECK_THROWS_AS(line::fj::fj_tail_ordstat<double>(ET, VT, 4, 1.0, 2), line::InputError);
}
