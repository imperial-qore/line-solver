/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MANJUNATH_H
#define LINE_API_PFQN_MANJUNATH_H

/**
 * Exact normalizing constant of a closed multiclass product-form network whose
 * state space carries arbitrary linear integer constraints (Manjunath-Sikdar).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_manjunath.m and
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_manjunath.java.
 *
 * This is the queueing-network half of the transform technique of which
 * `lossn_manjunath` is the loss-network half. The two solve the same problem --
 * sum a product form over an irregular integer state space -- from opposite ends
 * of the paper: `lossn_manjunath` implements Section 2.2, a set of '<=' rows over
 * the Poisson terms nu^n/n!, while this routine implements Section 3 together
 * with Section 5.3, a MIXED set of '=', '<=' and '>' rows over the BCMP terms,
 * where the population constraint of a closed network is itself one of the
 * equalities.
 *
 * THE MODEL. M queueing stations (rows of L) and Mz delay stations (rows of Z)
 * serve R closed classes with populations N. With n_i = sum_r n_ir,
 *
 *   p(n) = (1/G) prod_{i queueing} n_i! prod_r L_ir^{n_ir}/n_ir!
 *                prod_{i delay}         prod_r Z_ir^{n_ir}/n_ir!
 *
 * Every state obeys the R population equalities sum_i n_ir = N_r; the caller may
 * impose any number of further rows sum_{i,r} A(j, i + S*r) n_ir {=,<=,>} b(j)
 * with S = M + Mz, i.e. A acts on the (M+Mz) x R occupancy read column by column
 * with the queueing stations first. With no extra rows the result is exactly
 * `pfqn_ca`'s normalizing constant, which is the parity oracle used by the
 * tests.
 *
 * WHY THE GENERATING FUNCTION IS A PRODUCT, AND WHERE THE n_i! GOES. Marking
 * class r by z_r and row j by y_j, and writing u_ir = z_r prod_j y_j^{A(j,i+S r)},
 * the sum over the occupancies of a single QUEUEING station is, by the
 * multinomial theorem,
 *
 *   sum_{n_i.} n_i! prod_r (L_ir u_ir)^{n_ir}/n_ir!
 *     = sum_k (sum_r L_ir u_ir)^k = 1 / (1 - sum_r L_ir u_ir),
 *
 * so the n_i! that couples the classes is exactly what turns the station's
 * factor from an exponential into a geometric one. The paper reaches the same
 * place through the Euler integral n! = int_0^inf e^-t t^n dt (Eqns 16-18),
 * which is that geometric series evaluated; the closed form is used here because
 * there is then no quadrature to discretize -- and, decisively for this port, no
 * transcendental, so the whole computation stays rational.
 *
 * WHY THIS ONE RUNS AT EXACT ARITHMETIC. Every operation on the series is an
 * addition or a multiplication of demands, plus a division by a small integer in
 * the delay convolution, all of which are rational whenever L and Z are. G is
 * therefore exact under Arith::Exact and only `lG` is transcendental, obtained
 * through `num_traits<T>::log_as_double`. The double-only power-of-two rescaling
 * of `pfqn_ca` is applied on the same terms and for the same reason: exact
 * rationals have no exponent range to leave, and scaling them would only inflate
 * their denominators.
 *
 * THE ELIMINATION ORDER IS THE MEMORY BOUND. Variable y_j is created when the
 * first station its row touches is multiplied in and discharged immediately
 * after the last, so peak memory is prod_r (N_r+1) times the product of (b_j+1)
 * over the SIMULTANEOUSLY LIVE rows, not over all rows. The class axes are live
 * throughout, so prod_r (N_r+1) is a floor -- the same lattice `pfqn_ca` walks.
 * The peak is bounded by `PfqnManjunathOptions::max_live_states` and a model
 * above it is refused by name rather than allowed to exhaust the machine.
 *
 * SCOPE. Load-dependent and multiserver stations are NOT covered: their
 * per-station term is not geometric, and the multiclass n_i! coupling used above
 * then breaks. Use `pfqn_gld` or `pfqn_conwayms`. A and b must be integer valued
 * and A nonnegative, since the residue argument counts whole units; a fractional
 * entry is refused rather than rounded.
 *
 * Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 * Evaluation of Product Form Expressions Over Irregular Multidimensional Integer
 * Spaces. Sections 3 and 5.3.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Controls of pfqn_manjunath. The reference has none; the fields are port-local. */
struct PfqnManjunathOptions {
    /**
     * Cap on the number of series coefficients held at once. The default is 2^26
     * coefficients, half a gigabyte at double. Raise it deliberately.
     */
    std::size_t max_live_states = static_cast<std::size_t>(1) << 26;
    /**
     * Also return the per-class decomposition. Requires the ONE configuration in
     * which the truncated product form is the EXACT stationary law: a single
     * queueing station inside the region and a SINGLE DELAY STATION OUTSIDE IT.
     * Anything else is refused by name rather than answered wrongly.
     */
    bool stats = false;
};

/**
 * Result of pfqn_manjunath.
 *
 * The decomposition fields are populated only when `PfqnManjunathOptions::stats`
 * is set. They are `double` in every arithmetic, exactly as `lG` is, because
 * each is a RATIO of normalizing constants taken in the log domain -- which is
 * what keeps the four codebases comparable digit for digit and what removes the
 * range problem from a ratio of two very large constants. `G` itself stays exact
 * under Arith::Exact.
 *
 * `Q + think + blocked == N` exactly: a refused admission is a DELETED
 * transition, so a blocked job never leaves the delay, and because the think
 * time is exponential a held job is indistinguishable from one still thinking.
 * Little's law is what separates the two.
 */
template <class T>
struct PfqnManjunathResult {
    T G;                          ///< normalizing constant in the requested arithmetic
    double lG = 0.0;              ///< log of the constant, always a double
    std::size_t peak_states = 0;  ///< peak live series coefficients, the realised cost
    bool has_stats = false;       ///< whether the decomposition below was computed
    std::vector<double> Q;        ///< mean class r jobs at the queueing station
    std::vector<double> X;        ///< class r cycle throughput
    std::vector<double> U;        ///< class r utilization of the queueing station
    std::vector<double> think;    ///< class r jobs genuinely thinking, X_r * Z_r
    std::vector<double> blocked;  ///< class r jobs held at the delay by the constraint
    std::vector<double> delay;    ///< class r jobs at the delay, think + blocked
};

namespace detail {

/** MATLAB's round: half away from zero, which std::lround also does. */
inline int manjunath_round(double x) { return static_cast<int>(std::lround(x)); }


/** Early exit, carrying an all-zero decomposition when one was asked for. */
template <class T>
PfqnManjunathResult<T> manjunath_empty(const T& G, double lG, std::size_t peak,
                                       const PfqnManjunathOptions& options, std::size_t R) {
    PfqnManjunathResult<T> out;
    out.G = G;
    out.lG = lG;
    out.peak_states = peak;
    if (options.stats) {
        out.has_stats = true;
        out.Q.assign(R, 0.0);
        out.X.assign(R, 0.0);
        out.U.assign(R, 0.0);
        out.think.assign(R, 0.0);
        out.blocked.assign(R, 0.0);
        out.delay.assign(R, 0.0);
    }
    return out;
}

/** Forward declaration: the decomposition calls the entry point recursively. */
template <class T>
void manjunath_stats(PfqnManjunathResult<T>& out, const Matrix<T>& L, const std::vector<int>& N,
                     const Matrix<T>& Z, const std::vector<std::vector<long>>& A,
                     const std::vector<long>& b, const std::string& sense, double lG,
                     std::size_t M, std::size_t Mz, std::size_t S, std::size_t R,
                     const PfqnManjunathOptions& options);

/**
 * Coefficient-domain evaluation of the multiple contour integral of Eqn 9 for a
 * set of '=' and '<=' rows. The series lives in a flat vector indexed
 * column-major, its first R axes the class markers z_r of extent N_r+1
 * throughout and its remaining J axes the row markers y_j, of extent 1 while row
 * j is not live and b_j+1 while it is.
 */
template <class T>
T pfqn_manjunath_series(const Matrix<T>& L, const Matrix<T>& Z, const std::vector<int>& N,
                        const std::vector<std::vector<long>>& A, const std::vector<long>& b,
                        const std::string& sense, const PfqnManjunathOptions& options,
                        std::size_t& peak) {
    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t Mz = Z.empty() ? 0 : Z.rows();
    const std::size_t S = M + Mz;
    const std::size_t R = N.size();
    const std::size_t J = b.size();
    const std::size_t D = R + J;
    const T zero = num_traits<T>::from_int(0);

    // Row j is created at the first station it touches and discharged after the
    // last, so only an induced width of rows is ever live. A row that reached
    // here touches at least one station, the trivial ones having been decided.
    std::vector<std::size_t> first(J, 0), last(J, 0);
    for (std::size_t j = 0; j < J; ++j) {
        bool seen = false;
        for (std::size_t i = 0; i < S; ++i) {
            bool touch = false;
            for (std::size_t r = 0; r < R; ++r)
                if (A[j][i + S * r] != 0) touch = true;
            if (!touch) continue;
            if (!seen) {
                first[j] = i;
                seen = true;
            }
            last[j] = i;
        }
        if (!seen)
            throw NumericError("pfqn_manjunath: a constraint row with no nonzero entry reached "
                               "the series; the rule was not reduced");
    }

    std::vector<std::size_t> dims(D, 1);
    for (std::size_t r = 0; r < R; ++r) dims[r] = static_cast<std::size_t>(N[r]) + 1;
    std::size_t P = 1;
    for (std::size_t d = 0; d < D; ++d) P *= dims[d];
    std::vector<T> ser(P, zero);
    ser[0] = num_traits<T>::from_int(1);
    if (P > peak) peak = P;

    for (std::size_t i = 0; i < S; ++i) {
        for (std::size_t j = 0; j < J; ++j) {
            if (first[j] != i) continue;
            const std::size_t newdim = static_cast<std::size_t>(b[j]) + 1;
            std::size_t pre = 1, post = 1;
            for (std::size_t d = 0; d < R + j; ++d) pre *= dims[d];
            for (std::size_t d = R + j + 1; d < D; ++d) post *= dims[d];
            if (newdim != 0 && pre * post > options.max_live_states / newdim)
                throw UnsupportedError(
                    "pfqn_manjunath: the exact transform would hold more than " +
                    std::to_string(options.max_live_states) +
                    " series coefficients at once. The floor is the product of (N_r+1) over the "
                    "classes, on top of which each simultaneously live constraint row multiplies "
                    "by (b_j+1); raise PfqnManjunathOptions::max_live_states deliberately");
            // The existing content keeps its coefficients and enters at degree
            // zero in the new variable: nothing so far carries a power of it.
            std::vector<T> grown(pre * newdim * post, zero);
            for (std::size_t q = 0; q < post; ++q)
                for (std::size_t p = 0; p < pre; ++p)
                    grown[p + q * pre * newdim] = ser[p + q * pre];
            ser.swap(grown);
            dims[R + j] = newdim;
            if (ser.size() > peak) peak = ser.size();
        }

        std::vector<std::size_t> stride(D, 1);
        for (std::size_t d = 1; d < D; ++d) stride[d] = stride[d - 1] * dims[d - 1];

        // The monomial one class r job at station i contributes: z_r gains one
        // degree and y_j gains A(j, i + S*r). A row not live at this station has
        // a zero entry here by construction of first/last, so a dead axis is
        // never shifted.
        std::vector<std::vector<long>> delta(R, std::vector<long>(D, 0));
        std::vector<std::size_t> off(R, 0);
        std::vector<bool> fits(R, true);
        for (std::size_t r = 0; r < R; ++r) {
            delta[r][r] = 1;
            for (std::size_t j = 0; j < J; ++j) delta[r][R + j] = A[j][i + S * r];
            std::size_t o = 0;
            for (std::size_t d = 0; d < D; ++d) {
                if (delta[r][d] > static_cast<long>(dims[d]) - 1)
                    fits[r] = false;  // a single job already breaks the cut
                o += static_cast<std::size_t>(delta[r][d]) * stride[d];
            }
            off[r] = o;
        }

        if (i < M) {
            // Queueing station: solve (1 - sum_r L_ir u_ir) x = ser in place.
            // Every monomial of the operator raises the total class degree by
            // one, so p - off[r] is always a strictly smaller flat index and a
            // sweep in increasing flat index reads only final coefficients: the
            // sweep IS the solve, not an iterate.
            std::vector<std::size_t> sub(D, 0);
            for (std::size_t p = 0; p < ser.size(); ++p) {
                for (std::size_t r = 0; r < R; ++r) {
                    if (!fits[r] || L(i, r) == zero) continue;
                    bool ok = true;
                    for (std::size_t d = 0; d < D && ok; ++d)
                        if (static_cast<long>(sub[d]) < delta[r][d]) ok = false;
                    if (ok) ser[p] += T(L(i, r) * ser[p - off[r]]);
                }
                for (std::size_t d = 0; d < D; ++d) {
                    if (++sub[d] < dims[d]) break;
                    sub[d] = 0;
                }
            }
        } else {
            // Delay station: no n_i! coupling, so the factor is a product of
            // exponentials, one per class, each convolved in term by term. There
            // is no first-order recurrence to exploit here, which is why the
            // delay costs a factor of the population that the queueing station
            // does not.
            for (std::size_t r = 0; r < R; ++r) {
                if (!fits[r] || Z(i - M, r) == zero || N[r] == 0) continue;
                std::vector<T> nxt(ser);
                std::vector<T> term(ser);
                for (int n = 1; n <= N[r]; ++n) {
                    std::vector<T> shifted(ser.size(), zero);
                    std::vector<std::size_t> sub(D, 0);
                    for (std::size_t p = 0; p < ser.size(); ++p) {
                        bool ok = true;
                        for (std::size_t d = 0; d < D && ok; ++d)
                            if (static_cast<long>(sub[d]) < delta[r][d]) ok = false;
                        if (ok) shifted[p] = term[p - off[r]];
                        for (std::size_t d = 0; d < D; ++d) {
                            if (++sub[d] < dims[d]) break;
                            sub[d] = 0;
                        }
                    }
                    const T c = Z(i - M, r) / num_traits<T>::from_int(n);
                    bool any = false;
                    for (std::size_t p = 0; p < shifted.size(); ++p) {
                        shifted[p] = T(c * shifted[p]);
                        if (!(shifted[p] == zero)) any = true;
                    }
                    term.swap(shifted);
                    if (!any) break;
                    for (std::size_t p = 0; p < nxt.size(); ++p) nxt[p] += term[p];
                }
                ser.swap(nxt);
            }
        }

        for (std::size_t j = 0; j < J; ++j) {
            if (last[j] != i) continue;
            // Discharge marker j. The multiplier (y^{b+1}-1)/(y-1) of a '<=' row
            // turns its residue into the partial sum of the coefficients of
            // degrees 0..b, and the 1/y^{b+1} of an '=' row picks degree b.
            std::size_t pre = 1, post = 1;
            for (std::size_t d = 0; d < R + j; ++d) pre *= dims[d];
            for (std::size_t d = R + j + 1; d < D; ++d) post *= dims[d];
            const std::size_t dj = dims[R + j];
            std::vector<T> out(pre * post, zero);
            if (sense[j] == 'E') {
                const std::size_t rhs = static_cast<std::size_t>(b[j]);
                for (std::size_t q = 0; q < post; ++q)
                    for (std::size_t p = 0; p < pre; ++p)
                        out[p + q * pre] = ser[p + rhs * pre + q * pre * dj];
            } else {
                for (std::size_t q = 0; q < post; ++q)
                    for (std::size_t d = 0; d < dj; ++d)
                        for (std::size_t p = 0; p < pre; ++p)
                            out[p + q * pre] += ser[p + d * pre + q * pre * dj];
            }
            ser.swap(out);
            dims[R + j] = 1;
        }
    }

    std::size_t expect = 1;
    for (std::size_t r = 0; r < R; ++r) expect *= static_cast<std::size_t>(N[r]) + 1;
    if (ser.size() != expect)
        throw NumericError("pfqn_manjunath: a constraint row was never discharged; the "
                           "elimination order is inconsistent");
    // The closed network's own equalities: degree exactly N_r in every class.
    std::size_t flat = 0, st = 1;
    for (std::size_t r = 0; r < R; ++r) {
        flat += static_cast<std::size_t>(N[r]) * st;
        st *= static_cast<std::size_t>(N[r]) + 1;
    }
    return ser[flat];
}

}  // namespace detail

/**
 * @param L      (M x R) service demands at the queueing stations; may be empty
 * @param N      (R) population per class
 * @param Z      (Mz x R) think times at the delay stations; may be empty
 * @param A      (J x (M+Mz)*R) extra constraint coefficients, nonnegative integers
 * @param b      (J) extra constraint right-hand sides, integers
 * @param sense  one character per row, 'E' (=), 'L' (<=) or 'G' (>); empty means all 'L'
 */
template <class T>
PfqnManjunathResult<T> pfqn_manjunath(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z, const Matrix<T>& A,
                                      const std::vector<long>& b, const std::string& sense,
                                      const PfqnManjunathOptions& options = {}) {
    const std::size_t R = N.size();
    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t Mz = Z.empty() ? 0 : Z.rows();
    const std::size_t S = M + Mz;
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_manjunath: L and N disagree on the class count");
    if (!Z.empty() && Z.cols() != R)
        throw InputError("pfqn_manjunath: Z and N disagree on the class count");

    const std::size_t Jin = b.size();
    if (Jin > 0 && (A.rows() != Jin || A.cols() != S * R))
        throw InputError("pfqn_manjunath: A must be J x (M+Mz)*R, acting on the occupancy read "
                         "column by column with the queueing stations first");
    std::string sn = sense;
    if (sn.empty()) sn = std::string(Jin, 'L');
    if (sn.size() != Jin)
        throw InputError("pfqn_manjunath: sense must have one character per row of b");
    for (std::size_t j = 0; j < Jin; ++j)
        if (sn[j] != 'E' && sn[j] != 'L' && sn[j] != 'G')
            throw InputError("pfqn_manjunath: sense must contain only 'E' (=), 'L' (<=) or "
                             "'G' (>)");

    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0)
            return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), 0, options, R);

    // The integer view of the extra rows, refusing a fractional or negative
    // entry by name: the residue argument counts whole units.
    std::vector<std::vector<long>> Ai(Jin, std::vector<long>(S * R, 0));
    for (std::size_t j = 0; j < Jin; ++j)
        for (std::size_t c = 0; c < S * R; ++c) {
            const double a = num_traits<T>::to_double(A(j, c));
            if (a < 0.0 || std::fabs(a - std::round(a)) > 1e-9)
                throw InputError("pfqn_manjunath: A must contain nonnegative integers; the "
                                 "residue argument counts whole units");
            Ai[j][c] = std::lround(a);
        }

    // A row of zeros constrains nothing, so it is decided here rather than
    // carried as a one-coefficient dimension: 0 = b, 0 <= b and 0 > b are each
    // settled by the sign of b alone. Same for a negative right-hand side, which
    // no nonnegative combination can meet ('E', 'L') or can fail to beat ('G').
    std::vector<std::vector<long>> Ak;
    std::vector<long> bk;
    std::string sk;
    for (std::size_t j = 0; j < Jin; ++j) {
        bool trivial = true;
        for (std::size_t c = 0; c < S * R; ++c)
            if (Ai[j][c] != 0) trivial = false;
        if (sn[j] == 'E') {
            if (b[j] < 0 || (trivial && b[j] != 0))
                return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), 0, options, R);
            if (trivial) continue;
        } else if (sn[j] == 'L') {
            if (b[j] < 0)
                return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), 0, options, R);
            if (trivial) continue;
        } else {
            if (b[j] < 0) continue;  // 0 > negative always holds
            if (trivial)
                return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), 0, options, R);
        }
        Ak.push_back(Ai[j]);
        bk.push_back(b[j]);
        sk.push_back(sn[j]);
    }
    const std::size_t J = bk.size();

    if (S == 0) {
        // No station: the only state is empty, admissible when every class is too.
        bool empty = true;
        for (std::size_t r = 0; r < R; ++r)
            if (N[r] != 0) empty = false;
        if (empty) return detail::manjunath_empty<T>(num_traits<T>::from_int(1), 0.0, 1, options, R);
        return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), 0, options, R);
    }

    long Nt = 0;
    for (std::size_t r = 0; r < R; ++r) Nt += N[r];

    // Every monomial that survives the extraction has total degree sum(N) in the
    // demands, so a common power-of-two rescaling moves lG by a known amount and
    // nothing else. Applied only under double, on the same grounds as pfqn_ca:
    // exact rationals have no exponent range to leave.
    int kscale = 0;
    if (std::is_same<T, double>::value && Nt > 0) {
        double lGest = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < M; ++i) {
            double t = 0.0;
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r)
                if (N[r] > 0) {
                    const double lir = num_traits<T>::to_double(L(i, r));
                    if (lir > 0)
                        t += N[r] * std::log(lir);
                    else
                        ok = false;
                }
            if (ok && t > lGest) lGest = t;
        }
        if (Mz > 0) {
            double t = 0.0;
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r)
                if (N[r] > 0) {
                    double zs = 0.0;
                    for (std::size_t k = 0; k < Mz; ++k) zs += num_traits<T>::to_double(Z(k, r));
                    if (zs > 0)
                        t += N[r] * std::log(zs) - std::lgamma(N[r] + 1.0);
                    else
                        ok = false;
                }
            if (ok && t > lGest) lGest = t;
        }
        if (std::isfinite(lGest))
            kscale = detail::manjunath_round(lGest / (static_cast<double>(Nt) * std::log(2.0)));
    }
    Matrix<T> Ls = L;
    Matrix<T> Zs = Z;
    if (kscale != 0) {
        const T c = num_traits<T>::from_int(1) /
                    num_pow_int(num_traits<T>::from_int(2), static_cast<unsigned>(std::abs(kscale)));
        const T f = (kscale > 0) ? c : num_traits<T>::from_int(1) / c;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) Ls(i, r) = T(Ls(i, r) * f);
        for (std::size_t k = 0; k < Mz; ++k)
            for (std::size_t r = 0; r < R; ++r) Zs(k, r) = T(Zs(k, r) * f);
    }

    // A '>' row is the complement of a '<=' row at the same right-hand side,
    // which is how the paper discharges it (Eqn 6). With several such rows the
    // product of the complements expands by inclusion-exclusion, so the series is
    // evaluated once per subset of them, with the subset's rows re-entered as
    // '<=' and the rest dropped. Exact, and the only place the cost is
    // exponential -- in the number of '>' rows, normally zero.
    std::vector<std::size_t> gt;
    for (std::size_t j = 0; j < J; ++j)
        if (sk[j] == 'G') gt.push_back(j);
    const std::size_t K = gt.size();
    T Gs = num_traits<T>::from_int(0);
    std::size_t peak = 0;
    for (std::size_t mask = 0; mask < (static_cast<std::size_t>(1) << K); ++mask) {
        std::vector<bool> on(J, false);
        for (std::size_t j = 0; j < J; ++j) on[j] = (sk[j] != 'G');
        std::size_t bits = 0;
        for (std::size_t t = 0; t < K; ++t)
            if (mask & (static_cast<std::size_t>(1) << t)) {
                on[gt[t]] = true;
                ++bits;
            }
        std::vector<std::vector<long>> As;
        std::vector<long> bs;
        std::string ss;
        for (std::size_t j = 0; j < J; ++j)
            if (on[j]) {
                As.push_back(Ak[j]);
                bs.push_back(bk[j]);
                ss.push_back(sk[j] == 'G' ? 'L' : sk[j]);
            }
        const T g = detail::pfqn_manjunath_series<T>(Ls, Zs, N, As, bs, ss, options, peak);
        if (bits % 2 == 0)
            Gs += g;
        else
            Gs -= g;
    }

    if (!(Gs > num_traits<T>::from_int(0))) {
        // Either the admissible set is empty or the '>' complements cancelled it.
        return detail::manjunath_empty<T>(num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity(), peak, options, R);
    }
    const double lG = num_traits<T>::log_as_double(Gs) +
                      static_cast<double>(Nt) * kscale * std::log(2.0);
    T G = Gs;
    if (kscale != 0) {
        const T two = num_traits<T>::from_int(2);
        const unsigned e = static_cast<unsigned>(std::abs(static_cast<long>(Nt) * kscale));
        const T p = num_pow_int(two, e);
        G = (kscale > 0) ? T(Gs * p) : T(Gs / p);
    }
    PfqnManjunathResult<T> out;
    out.G = G;
    out.lG = lG;
    out.peak_states = peak;
    if (options.stats)
        detail::manjunath_stats<T>(out, L, N, Z, Ak, bk, sk, lG, M, Mz, S, R, options);
    return out;
}

/** Overload without extra constraints: the plain closed-network constant. */
template <class T>
PfqnManjunathResult<T> pfqn_manjunath(const Matrix<T>& L, const std::vector<int>& N,
                                      const Matrix<T>& Z,
                                      const PfqnManjunathOptions& options = {}) {
    return pfqn_manjunath<T>(L, N, Z, Matrix<T>(), std::vector<long>(), std::string(), options);
}

/** Overload without think times or extra constraints. */
template <class T>
PfqnManjunathResult<T> pfqn_manjunath(const Matrix<T>& L, const std::vector<int>& N,
                                      const PfqnManjunathOptions& options = {}) {
    return pfqn_manjunath<T>(L, N, Matrix<T>(), Matrix<T>(), std::vector<long>(), std::string(),
                             options);
}

namespace detail {

/**
 * Per-class decomposition, for the one configuration in which the truncated
 * product form is the exact stationary law: a single queueing station inside the
 * region and a single delay station outside it.
 *
 * WHY THE CONFIGURATION IS NOT A CONVENIENCE. With one queueing station the state
 * is the queue occupancy alone (the delay holds the complement) and every
 * transition moves one job of one class by one unit, so the chain is a
 * multidimensional birth-death process. That process is reversible, and Kelly's
 * truncation theorem then applies verbatim: restricting it to the
 * coordinate-convex set A n <= b and renormalizing gives exactly the truncated
 * product form. Add a second queueing station and the delay -> q1 -> q2 -> delay
 * cycle destroys reversibility; truncation no longer preserves the product form,
 * measured at 131% relative error on the stationary law of a 2-class, N = [2 2]
 * instance. G and lG stay correct in every configuration; only the metrics are
 * withheld.
 *
 * WHERE THE BLOCKED JOBS SIT. Nowhere special: a refused admission is a DELETED
 * transition, so a blocked job never leaves the delay, and because the think time
 * is exponential a held job is indistinguishable from one still thinking. The
 * delay population carries both and Little's law separates them. This is NOT the
 * WAITQ rule of SolverSSA/SolverCTMC/JMT, which moves a refused job out of the
 * delay into a per-region FIFO counted at no station.
 */
template <class T>
void manjunath_stats(PfqnManjunathResult<T>& out, const Matrix<T>& L, const std::vector<int>& N,
                     const Matrix<T>& Z, const std::vector<std::vector<long>>& A,
                     const std::vector<long>& b, const std::string& sense, double lG,
                     std::size_t M, std::size_t Mz, std::size_t S, std::size_t R,
                     const PfqnManjunathOptions& options) {
    if (Mz != 1)
        throw InputError("pfqn_manjunath: the per-class decomposition needs exactly one delay "
                         "station, got " + std::to_string(Mz) +
                         ". Pass Z as a 1xR row of think times");
    if (M != 1)
        throw UnsupportedError(
            "pfqn_manjunath: the per-class decomposition needs exactly one queueing station, got " +
            std::to_string(M) +
            ". With two or more the delay->q1->q2->delay cycle makes the chain irreversible, "
            "Kelly truncation no longer holds, and the truncated product form is not the "
            "stationary law (measured at 131% error). G and lG are still returned and still "
            "correct as a sum over the admissible set");
    const std::size_t J = b.size();
    // The delay must sit OUTSIDE the region: its columns are S*r + (S-1).
    for (std::size_t r = 0; r < R; ++r) {
        const std::size_t dcol = S * r + (S - 1);
        for (std::size_t j = 0; j < J; ++j)
            if (A[j][dcol] != 0)
                throw InputError("pfqn_manjunath: constraint row(s) reference the delay station "
                                 "in class " + std::to_string(r + 1) + " (column " +
                                 std::to_string(dcol) + "). The delay must lie OUTSIDE the finite "
                                 "capacity region, because the decomposition charges every held "
                                 "job to it");
    }

    Matrix<T> Am(J, S * R);
    for (std::size_t j = 0; j < J; ++j)
        for (std::size_t c = 0; c < S * R; ++c) Am(j, c) = num_traits<T>::from_int(A[j][c]);

    PfqnManjunathOptions sub = options;
    sub.stats = false;  // the recursion answers with constants only

    // Ratios of normalizing constants are taken in the LOG domain, so the
    // internal power-of-two rescaling cancels without ever being reconstructed.
    out.Q.assign(R, 0.0);
    out.X.assign(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        const std::size_t qcol = S * r;  // column of (queueing station, class r)

        // Throughput. One class r job removed from the queue leaves a state of
        // population N - e_r whose admission rule is shifted by that job's own
        // requirement column, exactly as the loss network's g(C - A e_r).
        if (N[r] >= 1) {
            std::vector<int> Nr(N);
            Nr[r] -= 1;
            std::vector<long> br(J);
            for (std::size_t j = 0; j < J; ++j) br[j] = b[j] - A[j][qcol];
            const PfqnManjunathResult<T> s = pfqn_manjunath<T>(L, Nr, Z, Am, br, sense, sub);
            if (std::isfinite(s.lG)) out.X[r] = std::exp(s.lG - lG);
        }

        // Mean queue length from the marginal law. An '=' row is discharged by
        // picking a single coefficient, so each call returns the mass of exactly
        // that occupancy.
        for (int k = 1; k <= N[r]; ++k) {
            Matrix<T> Ak(J + 1, S * R);
            for (std::size_t j = 0; j < J; ++j)
                for (std::size_t c = 0; c < S * R; ++c) Ak(j, c) = Am(j, c);
            for (std::size_t c = 0; c < S * R; ++c)
                Ak(J, c) = num_traits<T>::from_int(c == qcol ? 1 : 0);
            std::vector<long> bk(b);
            bk.push_back(k);
            const PfqnManjunathResult<T> s =
                pfqn_manjunath<T>(L, N, Z, Ak, bk, sense + "E", sub);
            if (std::isfinite(s.lG)) out.Q[r] += k * std::exp(s.lG - lG);
        }
    }

    out.U.assign(R, 0.0);
    out.think.assign(R, 0.0);
    out.delay.assign(R, 0.0);
    out.blocked.assign(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        out.U[r] = out.X[r] * num_traits<T>::to_double(L(0, r));   // one server, one visit
        out.think[r] = out.X[r] * num_traits<T>::to_double(Z(0, r));  // Little's law
        out.delay[r] = static_cast<double>(N[r]) - out.Q[r];       // not at the queue
        out.blocked[r] = out.delay[r] - out.think[r];              // held there
    }
    out.has_stats = true;
}

}  // namespace detail

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MANJUNATH_H
