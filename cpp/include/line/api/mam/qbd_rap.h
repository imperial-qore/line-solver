/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_RAP_H
#define LINE_API_MAM_QBD_RAP_H

/**
 * Quasi-birth-death processes with rational arrival process components, and
 * the RAP/RAP/1 queue built on top of them.
 *
 * Templated port of matlab/src/api/mam/qbd_rap.m (including its local
 * qbd_rap_g) and matlab/src/api/mam/qbd_raprap1.m, following
 * N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes with Rational
 * Arrival Process Components", Stochastic Models 26(3), 2010, 309-334. The
 * equilibrium construction is their Theorem 7 and the stability test their
 * Corollary 8.
 *
 * The process is given by its repeating blocks (A0, A1, A2) -- A0 up, A1
 * local, A2 down -- and its boundary blocks (B0, B1). Unlike a Markovian QBD
 * the blocks need not be nonnegative; they are only required to be
 * conservative, (A0 + A1 + A2) e = 0 and (B0 + B1) e = 0. The analysis rests
 * on the prediction-process interpretation of a RAP, which is what lets a QBD
 * argument survive the loss of nonnegativity.
 *
 * Theorem 7, step by step:
 *   1. G solves A0 G^2 + A1 G + A2 = 0.
 *   2. U = A1 + A0 G.
 *   3. R = A0 (-U)^-1.
 *   4. pihat0 (B1 + R A2) = 0 with pihat0 e = 1.
 *   5. pi_0 = K pihat0 with K chosen so pi_0 (I - R)^-1 e = 1.
 *   6. pi_n = pi_0 R^n.
 * Positive recurrence holds iff Sp(R) < 1 and step 4 has a unique solution.
 *
 * COMPUTING G. The blocks are not nonnegative, so logarithmic and cyclic
 * reduction carry no convergence guarantee, and the paper leaves the general
 * case explicitly open (Section 6). Two paths, exactly as in the reference:
 *   - A2 of rank one, A2 = u v: then G = e v / (v e) solves the equation in
 *     closed form. Conservativity gives (A0 + A1) e = -A2 e = -u (v e), and G
 *     is idempotent, so A0 G^2 + A1 G = (A0 + A1) e v/(v e) = -u v = -A2.
 *     This is the case of the paper's own example.
 *   - otherwise, natural functional iteration G <- (-A1)^-1 (A2 + A0 G^2) as a
 *     warm start, then Newton on the Sylvester-form Jacobian
 *     (A0 G + A1) H + A0 H G = -(A0 G^2 + A1 G + A2), solved through its
 *     Kronecker expansion (I (x) (A0 G + A1) + G^T (x) A0) vec(H).
 * An unconverged G is never returned: the residual and the constraint G e = e
 * are both checked and a failure raises NumericError carrying both numbers.
 *
 * DIVERGENCES FROM THE REFERENCE, all in how a quantity is EXTRACTED rather
 * than in what it is, and all tested:
 *   - the right factor v of a rank-one A2 is taken as the row of A2 with the
 *     largest infinity norm instead of the top right singular vector. G
 *     depends on v only through v/(v e), and every nonzero row of a rank-one
 *     matrix is a scalar multiple of v, so the two agree exactly; this keeps
 *     the step inside the templated arithmetic instead of routing it through
 *     a double-precision SVD. The rank test itself still uses the singular
 *     values (util/eig.h), matching the reference's sv(2) <= 1e-10 sv(1).
 *   - the boundary vector of step 4 is obtained from the linear system
 *     x V = 0, sum(x) = 1 (qbd_detail::statvec) rather than from the last
 *     right singular vector of V^T. The solution is unique up to scale
 *     precisely when the reference's own second-smallest-singular-value test
 *     passes, so the accept/reject decision is unchanged, but the vector is
 *     computed at the working precision instead of in double.
 *   - rcond(-U) is replaced by the EXACT reciprocal 1-norm condition number
 *     1 / (||X||_1 ||X^-1||_1). MATLAB's rcond only estimates that quantity.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: the general G is a
 * fixed-point iteration plus Newton driven to a tolerance, and Sp(R) is the
 * modulus of an eigenvalue, which is algebraic and not rational. Sp(R) is
 * computed by converting R to double and calling LAPACK (util/eig.h), so at
 * Real<D> the STABILITY GATE is only double-accurate; every returned quantity
 * is computed at the full working precision. The gate is a comparison against
 * 1 - 1e-12 and any model that close to the null-recurrent boundary has an
 * unbounded queue anyway, which is why the precision loss is confined there.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace rap_detail {

/** Frobenius norm. */
template <class T>
T normfro(const Matrix<T>& A) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j) * A(i, j);
    using std::sqrt;
    return T(sqrt(s));
}

/** Largest absolute entry of a vector, MATLAB's norm(v, inf) for a vector. */
template <class T>
T vecnorminf(const std::vector<T>& v) {
    T best = num_traits<T>::from_int(0);
    for (const T& x : v) {
        const T a = num_abs(T(x));
        if (a > best) best = a;
    }
    return best;
}

/** Exact reciprocal 1-norm condition number, the quantity MATLAB's rcond estimates. */
template <class T>
T rcond1(const Matrix<T>& A) {
    Matrix<T> Ainv;
    try {
        Ainv = inverse(A);
    } catch (const NumericError&) {
        return num_traits<T>::from_int(0);
    }
    const T na = qbd_detail::norm1(A);
    const T ni = qbd_detail::norm1(Ainv);
    const T p = na * ni;
    if (p == num_traits<T>::from_int(0)) return num_traits<T>::from_int(0);
    return T(num_traits<T>::from_int(1) / p);
}

/** Machine epsilon of the working type, the reference's `eps`. */
template <class T>
T working_eps() {
    return num_traits<T>::from_double(2.220446049250313e-16);
}

/** Solves A0 G^2 + A1 G + A2 = 0 for G (the local qbd_rap_g of qbd_rap.m). */
template <class T>
Matrix<T> qbd_rap_g(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2,
                    const T& blockScale) {
    using namespace qbd_detail;
    const std::size_t m = A1.rows();
    const std::vector<T> e = ones<T>(m);
    const T resTol = num_traits<T>::from_double(1e-10) * blockScale;
    const T zero = num_traits<T>::from_int(0);

    // Rank-one A2 admits the closed form G = e v / (v e).
    Matrix<double> A2d(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) A2d(i, j) = num_traits<T>::to_double(A2(i, j));
    const std::vector<double> sv = svd_values(A2d);
    if (sv.size() > 1 && sv[0] > 0.0 && sv[1] <= 1e-10 * sv[0]) {
        // Any nonzero row of a rank-one matrix spans its row space.
        std::size_t best = 0;
        T bestn = zero;
        for (std::size_t i = 0; i < m; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < m; ++j) s += num_abs(T(A2(i, j)));
            if (s > bestn) {
                bestn = s;
                best = i;
            }
        }
        std::vector<T> v(m);
        for (std::size_t j = 0; j < m; ++j) v[j] = A2(best, j);
        T ve = zero;
        for (const T& x : v) ve += x;
        if (num_abs(T(ve)) < num_traits<T>::from_double(1e-12) * vecnorminf(v))
            throw NumericError(
                "qbd_rap: A2 has rank one but its right factor v satisfies v*e = 0, so the "
                "closed form G = e*v/(v*e) is undefined");
        Matrix<T> G(m, m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) = v[j] / ve;
        const T res = normfro(madd(madd(matmul(A0, matmul(G, G)), matmul(A1, G)), A2));
        if (res > resTol)
            throw NumericError(
                "qbd_rap: the rank-one closed form for G leaves a residual "
                "||A0*G^2 + A1*G + A2||_F = " +
                std::to_string(num_traits<T>::to_double(res)) +
                ", above the roundoff level " +
                std::to_string(num_traits<T>::to_double(resTol)));
        return G;
    }

    const Matrix<T> negA1 = mscale(A1, T(num_traits<T>::from_int(-1)));
    if (rcond1(negA1) < working_eps<T>())
        throw NumericError(
            "qbd_rap: the local block A1 is singular, the iteration for G cannot be started");
    const Matrix<T> negA1inv = inverse(negA1);

    Matrix<T> G(m, m, zero);
    for (unsigned it = 0; it < 200; ++it) {
        const Matrix<T> Gnew = matmul(negA1inv, madd(A2, matmul(A0, matmul(G, G))));
        const T nG = normfro(G);
        const T scale = nG > num_traits<T>::from_int(1) ? nG : T(num_traits<T>::from_int(1));
        const bool done = normfro(msub(Gnew, G)) <= num_traits<T>::from_double(1e-14) * scale;
        G = Gnew;
        if (done) break;
    }

    // Newton on F(G) = A0 G^2 + A1 G + A2, through the Kronecker expansion of
    // the Sylvester operator. vec is COLUMN-major, matching MATLAB's res(:).
    for (unsigned it = 0; it < 100; ++it) {
        const Matrix<T> res = madd(madd(matmul(A0, matmul(G, G)), matmul(A1, G)), A2);
        if (normfro(res) <= resTol) break;
        const Matrix<T> M = madd(matmul(A0, G), A1);
        Matrix<T> J(m * m, m * m, zero);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j)
                for (std::size_t pp = 0; pp < m; ++pp)
                    for (std::size_t q = 0; q < m; ++q) {
                        T val = zero;
                        if (j == q) val += M(i, pp);
                        val += G(q, j) * A0(i, pp);
                        J(i + j * m, pp + q * m) = val;
                    }
        if (rcond1(J) < working_eps<T>()) break;
        std::vector<T> rhs(m * m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) rhs[i + j * m] = -res(i, j);
        const std::vector<T> y = solve(J, rhs);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) += y[i + j * m];
    }

    const T res = normfro(madd(madd(matmul(A0, matmul(G, G)), matmul(A1, G)), A2));
    const std::vector<T> Ge = mulvec(G, e);
    T ge = zero;
    for (std::size_t i = 0; i < m; ++i) {
        const T a = num_abs(T(Ge[i] - e[i]));
        if (a > ge) ge = a;
    }
    if (!(res <= resTol) || ge > num_traits<T>::from_double(1e-8))
        throw NumericError(
            "qbd_rap: could not compute the matrix G for this QBD with RAP components: residual "
            "||A0*G^2 + A1*G + A2||_F = " +
            std::to_string(num_traits<T>::to_double(res)) + " against a tolerance of " +
            std::to_string(num_traits<T>::to_double(resTol)) + ", and ||G*e-e||_inf = " +
            std::to_string(num_traits<T>::to_double(ge)) +
            ". The blocks are not nonnegative, so neither the functional iteration nor Newton's "
            "method is guaranteed to converge; the justification of algorithms for G in this "
            "setting is an open problem in Section 6 of Bean and Nielsen (2010). Supply a model "
            "with a rank-one A2, for which G is available in closed form.");
    return G;
}

}  // namespace rap_detail

/** Everything qbd_rap returns. */
template <class T>
struct QbdRapResult {
    std::vector<T> levelProb;  ///< marginal level probabilities, levels 0..numLevels
    T QN;                      ///< exact mean queue length, pi0 R (I-R)^-2 e
    Matrix<T> R, G, U;
    double spr;         ///< Sp(R), from a double eigensolve (see the header note)
    Matrix<T> pqueue;   ///< (numLevels+1) x m, row n holding pi_n
    std::vector<T> pi0;
};

/**
 * Equilibrium analysis of a QBD with RAP components (qbd_rap.m).
 *
 * @param A0,A1,A2 repeating blocks, up / local / down
 * @param B0,B1    boundary up and local blocks at level 0
 * @param numLevels highest level reported
 */
template <class T>
QbdRapResult<T> qbd_rap(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2,
                        const Matrix<T>& B0, const Matrix<T>& B1, std::size_t numLevels) {
    static_assert(num_traits<T>::has_transcendental, "qbd_rap requires transcendental arithmetic");
    using namespace qbd_detail;
    using namespace rap_detail;

    const std::size_t m = A1.rows();
    if (A1.cols() != m || A0.rows() != m || A0.cols() != m || A2.rows() != m || A2.cols() != m ||
        B0.rows() != m || B0.cols() != m || B1.rows() != m || B1.cols() != m)
        throw InputError("qbd_rap: all QBD blocks must be square and of the same order");
    if (m == 0) throw InputError("qbd_rap: empty blocks");

    const std::vector<T> e = ones<T>(m);
    const Matrix<T> I = eye<T>(m);
    T blockScale = num_traits<T>::from_int(1);
    {
        const T c[3] = {normfro(A0), normfro(A1), normfro(A2)};
        for (int k = 0; k < 3; ++k)
            if (c[k] > blockScale) blockScale = c[k];
    }
    const T conservTol = num_traits<T>::from_double(1e-8) * blockScale;

    const std::vector<T> ce = mulvec(madd(madd(A0, A1), A2), e);
    if (vecnorminf(ce) > conservTol)
        throw InputError(
            "qbd_rap: the repeating blocks are not conservative, ||(A0+A1+A2)*e||_inf = " +
            std::to_string(num_traits<T>::to_double(vecnorminf(ce))) +
            ". A QBD with RAP components requires (A0+A1+A2)*e = 0.");
    const std::vector<T> be = mulvec(madd(B0, B1), e);
    if (vecnorminf(be) > conservTol)
        throw InputError(
            "qbd_rap: the boundary blocks are not conservative, ||(B0+B1)*e||_inf = " +
            std::to_string(num_traits<T>::to_double(vecnorminf(be))) +
            ". A QBD with RAP components requires (B0+B1)*e = 0 at level 0.");

    QbdRapResult<T> out;
    out.G = qbd_rap_g(A0, A1, A2, blockScale);
    out.U = madd(A1, matmul(A0, out.G));
    const Matrix<T> negU = mscale(out.U, T(num_traits<T>::from_int(-1)));
    if (rcond1(negU) < working_eps<T>())
        throw NumericError(
            "qbd_rap: the matrix U = A1 + A0*G is singular, R = A0*inv(-U) does not exist");
    out.R = matmul(A0, inverse(negU));

    // Corollary 8(i). See the header note on the precision of this gate.
    Matrix<double> Rd(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) Rd(i, j) = num_traits<T>::to_double(out.R(i, j));
    out.spr = spectral_radius(Rd);
    if (out.spr >= 1.0 - 1e-12)
        throw NumericError("qbd_rap: the process is not positive recurrent, Sp(R) = " +
                           std::to_string(out.spr) +
                           " >= 1 (Corollary 8 of Bean and Nielsen, 2010)");

    // Step 4, with the reference's singular-value diagnostics on V.
    const Matrix<T> V = madd(B1, matmul(out.R, A2));
    Matrix<double> Vd(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) Vd(i, j) = num_traits<T>::to_double(V(i, j));
    const std::vector<double> sv = svd_values(Vd);
    const double nullTol = 1e-8 * (sv[0] > 1.0 ? sv[0] : 1.0);
    if (sv[m - 1] > nullTol)
        throw NumericError(
            "qbd_rap: the boundary equation x*(B1 + R*A2) = 0 has no nontrivial solution "
            "(smallest singular value " +
            std::to_string(sv[m - 1]) + " against tolerance " + std::to_string(nullTol) +
            "), so the process is not positive recurrent (Corollary 8(ii))");
    if (m > 1 && sv[m - 2] <= nullTol)
        throw NumericError(
            "qbd_rap: the boundary equation x*(B1 + R*A2) = 0 has a solution space of dimension "
            "greater than one, the equilibrium vector is not unique");
    const std::vector<T> pihat0 = statvec(V);

    // Step 5.
    const std::vector<T> ImRinv_e = mulvec(inverse(msub(I, out.R)), e);
    T denom = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < m; ++i) denom += pihat0[i] * ImRinv_e[i];
    if (denom == num_traits<T>::from_int(0))
        throw NumericError("qbd_rap: the level-0 vector cannot be normalized");
    out.pi0.resize(m);
    for (std::size_t i = 0; i < m; ++i) out.pi0[i] = pihat0[i] / denom;

    // Consistency of the supplied boundary up-block.
    const std::vector<T> pi0R = vecmul(out.pi0, out.R);
    const std::vector<T> pi0R2 = vecmul(pi0R, out.R);
    std::vector<T> bal = vecmul(out.pi0, B0);
    const std::vector<T> t1 = vecmul(pi0R, A1);
    const std::vector<T> t2 = vecmul(pi0R2, A2);
    for (std::size_t i = 0; i < m; ++i) bal[i] += t1[i] + t2[i];
    const T pnorm = vecnorminf(out.pi0);
    const T balTol =
        conservTol * (pnorm > num_traits<T>::from_int(1) ? pnorm : num_traits<T>::from_int(1));
    if (vecnorminf(bal) > balTol)
        throw NumericError(
            "qbd_rap: the boundary block B0 is inconsistent with the repeating blocks, "
            "||pi0*B0 + pi1*A1 + pi2*A2||_inf = " +
            std::to_string(num_traits<T>::to_double(vecnorminf(bal))) +
            ". The level-0 balance equation of Theorem 7 requires pi0*(B0-A0) = 0.");

    // Step 6.
    out.pqueue = Matrix<T>(numLevels + 1, m);
    std::vector<T> pin = out.pi0;
    out.levelProb.assign(numLevels + 1, num_traits<T>::from_int(0));
    for (std::size_t n = 0; n <= numLevels; ++n) {
        for (std::size_t j = 0; j < m; ++j) {
            out.pqueue(n, j) = pin[j];
            out.levelProb[n] += pin[j];
        }
        pin = vecmul(pin, out.R);
    }

    // Exact mean queue length, sum_n n pi0 R^n e = pi0 R (I-R)^-2 e.
    const Matrix<T> ImRinv = inverse(msub(I, out.R));
    const std::vector<T> w = mulvec(ImRinv, mulvec(ImRinv, e));
    out.QN = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < m; ++i) out.QN += pi0R[i] * w[i];
    return out;
}

/** qbd_rap with the reference defaults, B0 = A0, B1 = A1 and 20 levels. */
template <class T>
QbdRapResult<T> qbd_rap(const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2) {
    return qbd_rap(A0, A1, A2, A0, A1, static_cast<std::size_t>(20));
}

/** Everything qbd_raprap1 returns. */
template <class T>
struct QbdRapRap1Result {
    T XN;              ///< throughput, the arrival rate of the RAP
    T QN;              ///< mean queue length from the TRUNCATED level series
    T UN;              ///< utilization, 1 - P(level 0)
    Matrix<T> pqueue;  ///< level vectors, one row per level kept
    Matrix<T> R, G;
    T eta;                   ///< caudal characteristic, Sp(R)
    Matrix<T> B, L, F;       ///< the QBD blocks
    QbdRapResult<T> core;    ///< the full qbd_rap answer, including the exact QN
};

/**
 * RAP/RAP/1 queue (qbd_raprap1.m).
 *
 * The two RAPs are independent, so the QBD phase space is the product of the
 * two phase spaces with the ARRIVAL phase major, phase index (a-1) ns + s.
 * The Kronecker factors must not be swapped: downstream consumers index
 * pqueue by that convention.
 *
 * The truncation rule of the level series is the one of QBD_pi with
 * MaxNumComp 100: accumulate until the mass reaches 1 - 1e-10, capped at 101
 * level vectors. qbd_rap returns the exact mean queue length in closed form,
 * but QN here is deliberately the TRUNCATED sum, because that is the
 * documented return value and the JAR and Python ports must cut the tail at
 * the same point; core.QN carries the closed form for comparison.
 *
 * @param util if positive, the service RAP is rescaled to mean util/lambda_a
 * @param arrival the arrival RAP
 * @param service_in the service RAP, rescaled to the requested utilization
 */
template <class T>
QbdRapRap1Result<T> qbd_raprap1(const Map<T>& arrival, const Map<T>& service_in, const T& util) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_raprap1 requires transcendental arithmetic");
    const std::size_t na = arrival.order();
    const std::size_t ns = service_in.order();
    Map<T> service = service_in;
    if (util > num_traits<T>::from_int(0))
        service = map_scale(service, T(util / map_lambda(arrival)));

    QbdRapRap1Result<T> out;
    out.F = kron(arrival.D1, eye<T>(ns));
    out.L = qbd_detail::madd(kron(arrival.D0, eye<T>(ns)), kron(eye<T>(na), service.D0));
    out.B = kron(eye<T>(na), service.D1);
    const Matrix<T> B1 = kron(arrival.D0, eye<T>(ns));

    out.core = qbd_rap(out.F, out.L, out.B, out.F, B1, static_cast<std::size_t>(0));
    out.R = out.core.R;
    out.G = out.core.G;
    out.eta = num_traits<T>::from_double(out.core.spr);

    const std::size_t m = na * ns;
    const std::size_t maxNumComp = 100;
    std::vector<std::vector<T>> levels;
    levels.push_back(out.core.pi0);
    T sumpi = num_traits<T>::from_int(0);
    for (const T& v : out.core.pi0) sumpi += v;
    const T target = num_traits<T>::from_int(1) - num_traits<T>::from_double(1e-10);
    while (sumpi < target && levels.size() < 1 + maxNumComp) {
        const std::vector<T> nxt = vecmul(levels.back(), out.R);
        levels.push_back(nxt);
        for (const T& v : nxt) sumpi += v;
    }

    out.pqueue = Matrix<T>(levels.size(), m);
    out.QN = num_traits<T>::from_int(0);
    for (std::size_t n = 0; n < levels.size(); ++n) {
        T lp = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < m; ++j) {
            out.pqueue(n, j) = levels[n][j];
            lp += levels[n][j];
        }
        out.QN += num_traits<T>::from_int(static_cast<long>(n)) * lp;
    }

    T p0 = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < m; ++j) p0 += out.pqueue(0, j);
    out.UN = num_traits<T>::from_int(1) - p0;
    out.XN = map_lambda(arrival);
    return out;
}

/** qbd_raprap1 without rescaling the service process. */
template <class T>
QbdRapRap1Result<T> qbd_raprap1(const Map<T>& arrival, const Map<T>& service) {
    return qbd_raprap1(arrival, service, T(num_traits<T>::from_int(0)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_RAP_H
