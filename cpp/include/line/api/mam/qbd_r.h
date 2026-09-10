/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_R_H
#define LINE_API_MAM_QBD_R_H

/**
 * Quasi-birth-death processes: the rate matrix R, the fundamental matrix G,
 * the caudal characteristic, and the stationary distribution.
 *
 * Templated port of matlab/src/api/mam/qbd_R.m, qbd_R_logred.m, qbd_fundmat.m,
 * and of the boundary solve of matlab/lib/thirdparty/smcsolver/QBD_pi.m and the
 * caudal characteristic of QBD_Caudal.m. Cross-checked against
 * jar/src/main/java/jline/api/mam/Qbd_R.java and Qbd_R_logred.java.
 *
 * Block convention, level-independent QBD in continuous time:
 *
 *     Q = | Lbar  F     0     0    ... |
 *         | B     L     F     0    ... |
 *         | 0     B     L     F    ... |
 *         | ...                        |
 *
 * with B the backward (downward) block, L the local block and F the forward
 * (upward) block. R is the minimal non-negative solution of
 *
 *     F + R L + R^2 B = 0,
 *
 * which is the equation A0 + R A1 + R^2 A2 = 0 with A0 = F, A1 = L, A2 = B.
 * G is the minimal non-negative solution of B + L G + F G^2 = 0.
 *
 * ARITHMETIC. Three routines here compute R or G by a fixed-point iteration
 * driven to a tolerance -- successive substitution, logarithmic reduction and
 * cyclic reduction -- and one (qbd_caudal) brackets a spectral radius. None of
 * them terminates in a finite number of field operations, so the value they
 * return is an approximation no matter how the arithmetic is carried out;
 * running them at exact rational arithmetic would produce a rational number
 * with a denominator doubling at every cyclic-reduction step and still not the
 * exact R. They are therefore gated on num_traits<T>::has_transcendental,
 * which admits double and Real<D> and rejects Rational at compile time.
 *
 * Everything downstream of R is a finite rational computation and is left
 * un-gated: qbd_R_residual, qbd_pi (a null-vector solve plus a geometric tail)
 * and the moment formulas in qbd_mapmap1.h all instantiate at Rational. That
 * split is the point of the port -- given R to whatever accuracy, the boundary
 * probabilities and the queue-length moments carry no additional error.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace qbd_detail {

/** Elementwise A - B. */
template <class T>
Matrix<T> msub(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols()) throw InputError("qbd: shape mismatch");
    Matrix<T> C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) - B(i, j);
    return C;
}

/** Elementwise A + B. */
template <class T>
Matrix<T> madd(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols()) throw InputError("qbd: shape mismatch");
    Matrix<T> C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) + B(i, j);
    return C;
}

/** Elementwise s*A. */
template <class T>
Matrix<T> mscale(const Matrix<T>& A, const T& s) {
    Matrix<T> C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) * s;
    return C;
}

/** MATLAB's norm(X,1): the largest absolute column sum. */
template <class T>
T norm1(const Matrix<T>& A) {
    T best = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < A.cols(); ++j) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < A.rows(); ++i) s += num_abs(T(A(i, j)));
        if (s > best) best = s;
    }
    return best;
}

/** MATLAB's norm(X,inf): the largest absolute row sum. */
template <class T>
T norminf(const Matrix<T>& A) {
    T best = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < A.cols(); ++j) s += num_abs(T(A(i, j)));
        if (s > best) best = s;
    }
    return best;
}

/**
 * Stationary row vector of M: solves x M = 0 with sum(x) = 1 by replacing the
 * last column of M with ones and solving the transposed system, as MATLAB's
 * statvec does. Unlike mc::ctmc_solve this does NOT recompute the diagonal of
 * M: the level-zero block Lbar + R B is already a generator by construction,
 * and repairing its diagonal would silently absorb an error in R instead of
 * letting it show up in the residual.
 */
template <class T>
std::vector<T> statvec(const Matrix<T>& M) {
    const std::size_t m = M.rows();
    if (M.cols() != m) throw InputError("qbd statvec: matrix is not square");
    if (m == 0) throw InputError("qbd statvec: empty matrix");
    if (m == 1) return std::vector<T>{num_traits<T>::from_int(1)};
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> A(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) A(i, j) = (i == m - 1) ? one : M(j, i);
    std::vector<T> b(m, zero);
    b[m - 1] = one;
    return solve(A, b);
}

}  // namespace qbd_detail

/**
 * Residual of the defining equation of R, ||F + R L + R^2 B||_inf.
 *
 * Exact in rational arithmetic, which is what makes it a usable oracle: a
 * residual computed in the same double arithmetic that produced R can be small
 * simply because the error cancels, whereas evaluating the same expression on
 * the exact rationals of B, L, F and on the exact rational lift of R measures
 * only how far R itself is from a solution.
 */
template <class T>
T qbd_R_residual(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F, const Matrix<T>& R) {
    using namespace qbd_detail;
    const Matrix<T> res = madd(madd(F, matmul(R, L)), matmul(matmul(R, R), B));
    return norminf(res);
}

/** Residual of the defining equation of G, ||B + L G + F G^2||_inf. */
template <class T>
T qbd_G_residual(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F, const Matrix<T>& G) {
    using namespace qbd_detail;
    const Matrix<T> res = madd(madd(B, matmul(L, G)), matmul(F, matmul(G, G)));
    return norminf(res);
}

/**
 * R by successive substitutions (qbd_R.m): iterate R <- -(F + R^2 B) L^-1.
 *
 * Linearly convergent, at the rate of the caudal characteristic, so it is slow
 * near saturation; qbd_fundmat is quadratically convergent and is what
 * qbd_mapmap1 uses. Kept because it is the MATLAB entry point of the same name
 * and because its simplicity makes it a useful independent check on the
 * cyclic-reduction result.
 *
 * @param tol stopping tolerance on ||R_k - R_{k+1}||_1 (MATLAB uses 1e-12)
 * @param B backward (level down) block A_2
 * @param L local (within level) block A_1
 * @param F forward (level up) block A_0
 * @param iter_max iteration cap
 */
template <class T>
Matrix<T> qbd_R(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F, unsigned iter_max,
                const T& tol) {
    static_assert(num_traits<T>::has_transcendental, "qbd_R requires transcendental arithmetic");
    using namespace qbd_detail;
    const Matrix<T> Linv = inverse(L);
    const Matrix<T> Fil = matmul(F, Linv);
    const Matrix<T> BiL = matmul(B, Linv);
    const Matrix<T> negFil = mscale(Fil, T(num_traits<T>::from_int(-1)));
    Matrix<T> R = negFil;
    Matrix<T> Rprime = msub(negFil, matmul(matmul(R, R), BiL));
    for (unsigned it = 0; it < iter_max; ++it) {
        R = Rprime;
        Rprime = msub(negFil, matmul(matmul(R, R), BiL));
        if (norm1(msub(R, Rprime)) <= tol) break;
    }
    return Rprime;
}

/** qbd_R with the MATLAB defaults, 100000 iterations and tolerance 1e-12. */
template <class T>
Matrix<T> qbd_R(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F) {
    return qbd_R(B, L, F, 100000u, T(num_traits<T>::from_double(1e-12)));
}

/**
 * R by logarithmic reduction (qbd_R_logred.m).
 *
 * Builds the matrix S of the taboo probabilities of ever going down, doubling
 * the horizon at every step, then recovers R = -F (L + F S)^-1. Quadratically
 * convergent; the stopping test is on how close S e is to e, i.e. on how much
 * mass of the downward passage is still unaccounted for.
 */
template <class T>
Matrix<T> qbd_R_logred(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F, unsigned iter_max,
                       const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_R_logred requires transcendental arithmetic");
    using namespace qbd_detail;
    const std::size_t r = L.rows();
    const Matrix<T> Linv = inverse(L);
    const T minus = num_traits<T>::from_int(-1);
    Matrix<T> iLF = mscale(matmul(Linv, F), minus);
    Matrix<T> iLB = mscale(matmul(Linv, B), minus);
    Matrix<T> T_ = iLF;
    Matrix<T> S = iLB;
    const Matrix<T> I = eye<T>(r);
    const std::vector<T> e = ones<T>(r);
    for (unsigned it = 0; it < iter_max; ++it) {
        const Matrix<T> D = madd(matmul(iLF, iLB), matmul(iLB, iLF));
        const Matrix<T> Minv = inverse(msub(I, D));
        const Matrix<T> iLFn = matmul(Minv, matmul(iLF, iLF));
        const Matrix<T> iLBn = matmul(Minv, matmul(iLB, iLB));
        iLF = iLFn;
        iLB = iLBn;
        S = madd(S, matmul(T_, iLB));
        T_ = matmul(T_, iLF);
        // ||e - S e||_1 over the column vector, i.e. the sum of absolute
        // deficits, exactly as MATLAB's norm(ones - S*ones, 1).
        const std::vector<T> Se = mulvec(S, e);
        T dev = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < r; ++i) dev += num_abs(T(e[i] - Se[i]));
        if (dev <= tol) break;
    }
    const Matrix<T> U = madd(L, matmul(F, S));
    return mscale(matmul(F, inverse(U)), minus);
}

/** qbd_R_logred with the MATLAB defaults, 100000 iterations and tolerance 1e-12. */
template <class T>
Matrix<T> qbd_R_logred(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F) {
    return qbd_R_logred(B, L, F, 100000u, T(num_traits<T>::from_double(1e-12)));
}

/** G and R together, as returned by qbd_fundmat. */
template <class T>
struct QbdFundMat {
    Matrix<T> G;
    Matrix<T> R;
    unsigned iterations = 0;
};

/**
 * G and R by cyclic reduction (qbd_fundmat.m, the Bini-Meini logarithmic
 * reduction on the raw level blocks).
 *
 * The blocks are first uniformized by lambda = max(-diag(L)) into a discrete
 * QBD (Bm, Lm, Fm), G is accumulated over doubling horizons, and R is
 * recovered as R = Fm (I - (Lm + Fm G))^-1. The uniformization does not change
 * either G or R: G is a probability matrix of the embedded jump chain, and the
 * R of the discrete chain solves R = Fm + R Lm + R^2 Bm, which is F + R L +
 * R^2 B = 0 after multiplying through by lambda.
 *
 * Quadratically convergent, so 50 iterations is a generous bound even at
 * utilizations where successive substitution needs millions.
 */
template <class T>
QbdFundMat<T> qbd_fundmat(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F,
                          unsigned iter_max, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_fundmat requires transcendental arithmetic");
    using namespace qbd_detail;
    const std::size_t m = L.rows();
    if (L.cols() != m || B.rows() != m || B.cols() != m || F.rows() != m || F.cols() != m)
        throw InputError("qbd_fundmat: B, L and F must be square and of equal order");
    const Matrix<T> I = eye<T>(m);

    T lamb = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < m; ++i) {
        const T d = -L(i, i);
        if (d > lamb) lamb = d;
    }
    if (lamb <= num_traits<T>::from_int(0))
        throw NumericError("qbd_fundmat: the local block has no negative diagonal entry");
    const T inv_lamb = num_traits<T>::from_int(1) / lamb;
    const Matrix<T> Bm = mscale(B, inv_lamb);
    const Matrix<T> Lm = madd(mscale(L, inv_lamb), I);
    const Matrix<T> Fm = mscale(F, inv_lamb);

    Matrix<T> BF = inverse(msub(I, Lm));
    Matrix<T> BB = matmul(BF, Fm);
    BF = matmul(BF, Bm);
    Matrix<T> G = BF;
    Matrix<T> PI = BB;
    T check = num_traits<T>::from_int(1);
    unsigned numit = 0;
    while (check > tol && numit < iter_max) {
        const Matrix<T> Lstar = madd(matmul(BF, BB), matmul(BB, BF));
        const Matrix<T> Bstar = matmul(BB, BB);
        const Matrix<T> Fstar = matmul(BF, BF);
        const Matrix<T> Minv = inverse(msub(I, Lstar));
        BF = matmul(Minv, Fstar);
        BB = matmul(Minv, Bstar);
        G = madd(G, matmul(PI, BF));
        PI = matmul(PI, BB);
        const T nb = norminf(BB);
        const T nf = norminf(BF);
        check = nb < nf ? nb : nf;
        ++numit;
    }
    QbdFundMat<T> out;
    out.G = G;
    out.R = matmul(Fm, inverse(msub(I, madd(Lm, matmul(Fm, G)))));
    out.iterations = numit;
    return out;
}

/** qbd_fundmat with the MATLAB defaults, 50 iterations and tolerance 1e-14. */
template <class T>
QbdFundMat<T> qbd_fundmat(const Matrix<T>& B, const Matrix<T>& L, const Matrix<T>& F) {
    return qbd_fundmat(B, L, F, 50u, T(num_traits<T>::from_double(1e-14)));
}

/**
 * Caudal characteristic eta = sp(R), the decay rate of the queue-length tail.
 *
 * R is entrywise non-negative, so its spectral radius is its Perron root and
 * is bracketed by the Collatz-Wielandt bounds
 *
 *     min_i (R x)_i / x_i  <=  sp(R)  <=  max_i (R x)_i / x_i
 *
 * for any strictly positive x. Power iteration is run on I + R, which is
 * aperiodic whenever R is irreducible and keeps the iterate strictly positive,
 * and the midpoint of the bracket is returned once it is tighter than tol.
 *
 * This is the same quantity as MATLAB's QBD_Caudal, which brackets it by
 * bisecting on the dominant eigenvalue of A(eta) = B + L eta + F eta^2, and as
 * the JAR's spectralRadiusMapmap1, which takes the largest eigenvalue modulus
 * of R from a full eigendecomposition. The bracket form is preferred here
 * because it needs no eigensolver and reports its own accuracy.
 */
template <class T>
T qbd_caudal(const Matrix<T>& R, unsigned iter_max, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_caudal requires transcendental arithmetic");
    const std::size_t n = R.rows();
    if (R.cols() != n) throw InputError("qbd_caudal: matrix is not square");
    if (n == 0) throw InputError("qbd_caudal: empty matrix");
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> x = ones<T>(n);
    T lo = zero, hi = zero;
    for (unsigned it = 0; it < iter_max; ++it) {
        const std::vector<T> Rx = mulvec(R, x);
        lo = Rx[0] / x[0];
        hi = lo;
        for (std::size_t i = 1; i < n; ++i) {
            const T q = Rx[i] / x[i];
            if (q < lo) lo = q;
            if (q > hi) hi = q;
        }
        if (T(hi - lo) <= tol) break;
        // Advance with I + R: sp(I+R) = 1 + sp(R) and the iterate stays > 0.
        T s = zero;
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += Rx[i];
            s += x[i];
        }
        if (s == zero) throw NumericError("qbd_caudal: iterate collapsed to zero");
        for (std::size_t i = 0; i < n; ++i) x[i] /= s;
    }
    return T((lo + hi) / num_traits<T>::from_int(2));
}

/** qbd_caudal with 10000 iterations and tolerance 1e-14. */
template <class T>
T qbd_caudal(const Matrix<T>& R) {
    return qbd_caudal(R, 10000u, T(num_traits<T>::from_double(1e-14)));
}

/**
 * Stationary distribution of a QBD given R (QBD_pi.m, continuous-time branch,
 * default boundary).
 *
 * The level-zero vector solves pi_0 (Lbar + R B) = 0 normalized by
 * pi_0 (I - R)^-1 e = 1, and pi_{k+1} = pi_k R. Levels are generated until the
 * accumulated mass exceeds 1 - mass_tol or max_levels is reached.
 *
 * Un-gated: given R this is a linear solve, a matrix inverse and a geometric
 * recursion, all finite sequences of field operations. At Rational the level
 * probabilities are the exact ones implied by the R that was supplied.
 *
 * @param B      backward block (QBD_pi's B0)
 * @param Lbar   level-zero local block (QBD_pi's B1)
 * @param R the rate matrix of the QBD
 * @param max_levels level truncation used for the returned distribution
 * @param mass_tol probability mass left above the truncation that is tolerated
 * @return (levels x m) matrix, row k holding pi_k
 */
template <class T>
Matrix<T> qbd_pi(const Matrix<T>& B, const Matrix<T>& Lbar, const Matrix<T>& R,
                 std::size_t max_levels, const T& mass_tol) {
    using namespace qbd_detail;
    const std::size_t m = R.rows();
    if (R.cols() != m) throw InputError("qbd_pi: R is not square");
    if (Lbar.rows() != m || Lbar.cols() != m || B.rows() != m || B.cols() != m)
        throw InputError("qbd_pi: block orders disagree with R");
    if (max_levels == 0) throw InputError("qbd_pi: max_levels must be positive");
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);

    const Matrix<T> ImR = msub(eye<T>(m), R);
    const Matrix<T> ImRinv = inverse(ImR);

    std::vector<T> pi0 = statvec(madd(Lbar, matmul(R, B)));
    // Normalize so that the whole chain has unit mass: sum_k pi_0 R^k e = 1.
    const std::vector<T> t = vecmul(pi0, ImRinv);
    T tot = zero;
    for (const T& v : t) tot += v;
    if (tot == zero) throw NumericError("qbd_pi: degenerate normalization, sp(R) may exceed 1");
    for (T& v : pi0) v /= tot;

    std::vector<std::vector<T>> levels;
    levels.push_back(pi0);
    T mass = zero;
    for (const T& v : pi0) mass += v;
    while (levels.size() < max_levels && T(one - mass) > mass_tol) {
        const std::vector<T> nxt = vecmul(levels.back(), R);
        T add = zero;
        for (const T& v : nxt) add += v;
        levels.push_back(nxt);
        mass += add;
    }
    Matrix<T> out(levels.size(), m);
    for (std::size_t k = 0; k < levels.size(); ++k)
        for (std::size_t j = 0; j < m; ++j) out(k, j) = levels[k][j];
    return out;
}

/** qbd_pi with the MATLAB-side defaults, 20000 levels and mass tolerance 1e-10. */
template <class T>
Matrix<T> qbd_pi(const Matrix<T>& B, const Matrix<T>& Lbar, const Matrix<T>& R) {
    return qbd_pi(B, Lbar, R, static_cast<std::size_t>(20000),
                  T(num_traits<T>::from_double(1e-10)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_R_H
