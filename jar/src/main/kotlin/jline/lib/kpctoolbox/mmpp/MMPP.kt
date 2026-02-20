package jline.lib.kpctoolbox.mmpp

import jline.api.mam.*
import jline.api.mam.Mmpp2_fitc.Companion.mmpp2_fitc as api_mmpp2_fitc
import jline.api.mam.Mmpp2_fitc_approx.Companion.mmpp2_fitc_approx as api_mmpp2_fitc_approx
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Markov Modulated Poisson Process (MMPP) functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mmpp/
 * Delegates to jline.api.mam implementations where possible.
 */

/**
 * Fits a 2-state MMPP to match first three moments and G2 parameter.
 *
 * This is the core MMPP2 fitter matching MATLAB's mmpp2_fit3(E1,E2,E3,G2).
 * G2 specifies the ratio of consecutive autocorrelations: rho(i)/rho(i-1).
 *
 * @param E1 First moment (mean)
 * @param E2 Second moment
 * @param E3 Third moment
 * @param G2 Autocorrelation decay ratio
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fit3(E1: Double, E2: Double, E3: Double, G2: Double): MatrixCell {
    return jline.api.mam.mmpp2_fit(E1, E2, E3, G2)
}

/**
 * Fits MMPP2 from mean, SCV, skewness, and index of dispersion for counts.
 *
 * Matches MATLAB: mmpp2_fit1(mean, scv, skew, idc)
 *
 * @param mean Mean inter-arrival time
 * @param scv Squared coefficient of variation
 * @param skew Skewness (-1 for automatic)
 * @param idc Index of dispersion for counts
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fit1(mean: Double, scv: Double, skew: Double, idc: Double): MatrixCell {
    return jline.api.mam.mmpp2_fit1(mean, scv, skew, idc)
}

/**
 * Fits MMPP2 from mean, SCV, skewness, and G2 parameter.
 *
 * Matches MATLAB: mmpp2_fit2(mean, scv, skew, g2)
 *
 * @param mean Mean inter-arrival time
 * @param scv Squared coefficient of variation
 * @param skew Skewness
 * @param g2 Autocorrelation decay ratio
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fit2(mean: Double, scv: Double, skew: Double, g2: Double): MatrixCell {
    if (scv == 1.0) {
        return map_exponential(mean)
    }
    val E1 = mean
    val E2 = (1 + scv) * E1 * E1
    val E3 = -(2 * E1 * E1 * E1 - 3 * E1 * E2 - skew * FastMath.pow(E2 - E1 * E1, 1.5))
    return mmpp2_fit3(E1, E2, E3, g2)
}

/**
 * Fits MMPP2 from mean, SCV, skewness, and lag-1 autocorrelation.
 *
 * Matches MATLAB: mmpp2_fit4(mean, scv, skew, acf1)
 *
 * @param mean Mean inter-arrival time
 * @param scv Squared coefficient of variation
 * @param skew Skewness (-1 for automatic)
 * @param acf1 Lag-1 autocorrelation
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fit4(mean: Double, scv: Double, skew: Double, acf1: Double): MatrixCell {
    val E1 = mean
    val E2 = (1 + scv) * E1 * E1
    val E3 = if (skew == -1.0) {
        -1.0
    } else {
        -(2 * E1 * E1 * E1 - 3 * E1 * E2 - skew * FastMath.pow(E2 - E1 * E1, 1.5))
    }
    val rho0 = (1 - 1 / scv) / 2
    val g2 = acf1 / rho0
    return mmpp2_fit3(E1, E2, E3, g2)
}

/**
 * Fits MMPP2 from counting process statistics using Heffes-Lucantoni method.
 *
 * Matches MATLAB: mmpp2_fitc(mu, bt1, bt2, binf, m3t2, t1, t2)
 *
 * @param mu Arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 Third central moment at scale t2
 * @param t1 First time scale
 * @param t2 Second time scale
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fitc(mu: Double, bt1: Double, bt2: Double, binf: Double,
               m3t2: Double, t1: Double, t2: Double): MatrixCell {
    val result = api_mmpp2_fitc(mu, bt1, bt2, binf, m3t2, t1, t2)
    val MAP = MatrixCell(2)
    MAP[0] = result[0]
    MAP[1] = result[1]
    return MAP
}

/**
 * Fits MMPP2 from counting process statistics using optimization.
 *
 * Matches MATLAB: mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
 *
 * @param a Arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 Third central moment at scale t2
 * @param t1 First time scale
 * @param t2 Second time scale
 * @return Fitted MAP as {D0, D1}
 */
fun mmpp2_fitc_approx(a: Double, bt1: Double, bt2: Double, binf: Double,
                      m3t2: Double, t1: Double, t2: Double): MatrixCell {
    val result = api_mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
    val MAP = MatrixCell(2)
    MAP[0] = result[0]
    MAP[1] = result[1]
    return MAP
}

/**
 * Fits theoretical characteristics of a MAP(n) with a MMPP(2).
 *
 * Matches MATLAB: mmpp2_fitc_theoretical(map, t1, t2, tinf)
 *
 * @param MAP Input MAP as {D0, D1, ...}
 * @param t1 First time scale (default 1)
 * @param t2 Second time scale (default 10)
 * @param tinf Large time scale for asymptotic IDC (default 1e8)
 * @return Fitted MMPP2 as {D0, D1}
 */
fun mmpp2_fitc_theoretical(MAP: MatrixCell, t1: Double = 1.0, t2: Double = 10.0, tinf: Double = 1e8): MatrixCell {
    val a = map_count_mean(MAP, t1) / t1
    val bt1 = map_count_var(MAP, t1) / (a * t1)
    val bt2 = map_count_var(MAP, t2) / (a * t2)
    val binf = map_count_var(MAP, tinf) / (a * tinf)
    val mt2 = map_count_moment(MAP, t2, intArrayOf(1, 2, 3))
    val m3t2 = mt2[2] - 3 * mt2[1] * mt2[0] + 2 * mt2[0] * mt2[0] * mt2[0]

    return mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2)
}

/**
 * Generates a random MMPP with K states.
 *
 * Matches MATLAB: mmpp_rand(K) - generates random D0, D1 (diagonal), normalizes.
 *
 * @param K Number of states
 * @return Random MMPP as {D0, D1}
 */
fun mmpp_rand(K: Int): MatrixCell {
    return jline.api.mam.mmpp_rand(K)
}

/**
 * Generates a random MMPP with 2 states.
 *
 * @return Random MMPP as {D0, D1}
 */
fun mmpp_rand(): MatrixCell {
    return mmpp_rand(2)
}
