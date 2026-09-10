/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MMAPGK1_H
#define LINE_API_QSYS_QSYS_MMAPGK1_H

/**
 * The MMAP[K]/G[K]/1 FCFS queue: K customer types with class-dependent GENERAL
 * service, fed by a marked Markovian arrival process.
 *
 * WHY THIS IS NOT MMAPPH1FCFS. That routine needs every type's service to be
 * PHASE TYPE, because it builds a QBD whose phase carries the service phase.
 * Here the service laws are arbitrary and may differ in family across types --
 * deterministic for one, uniform for another -- so no finite phase carries
 * them, and the analysis has to run through transforms instead.
 *
 * THE METHOD, which is He's, theorem for theorem. FCFS makes the actual waiting
 * time of a customer the WORKLOAD it finds on arrival, so everything follows
 * from the joint transform of workload and arrival phase,
 * f(s)_j = E[exp(-s V) 1{phase = j}], which by He's Theorem 4.1 (eq. 4.6)
 * satisfies
 *
 *     f(s) [ s I + D0 + sum_k Dk gk(s) ] = s v0,                          (*)
 *
 * with v0 the idle-phase vector, his y0. v0 needs NO root search: the matrix U
 * solving U = D0 + sum_k Dk Fk(U), Fk(U) = int exp(U t) dFk(t), is his
 * eq. (4.4), the generator of the underlying Markov process obtained by
 * EXCISING the busy periods, and eq. (4.5) with Theorem 4.2 give y0 Q = 0 and
 * y0 e = 1 - rho, i.e. v0 = (1 - rho) pi_U. The same vector is what the
 * analyticity of (*) forces, since for every left eigenpair (w, u) of U one has
 * w [D0 + sum_k Dk gk(-u) - u I] = 0, so the roots of the bracket in the closed
 * right half plane are exactly s = -u over the spectrum of U; the two agree to
 * 2.5e-13, and the stationary route is taken because it needs no complex
 * eigenvector.
 *
 * The per-type actual waiting time is the workload seen by a type-k arrival,
 * biased by that type's own arrival block, his Theorem 5.1 eq. (5.1) summed
 * over the post-arrival phase:
 *
 *     E[exp(-s Wk)] = f(s) Dk e / lambda_k.
 *
 * SCOPE. He allows an arrival to be a BATCH carrying a sequence of types, and
 * his Theorem 5.3 then multiplies the transform by prod_{i<n} f*_{h_i}(s), the
 * service of the customers ahead of the tagged one WITHIN its own batch. This
 * header covers the single-customer-per-arrival case, his Special case 3.3,
 * where that product is empty, which is exactly the MMAP convention LINE
 * carries.
 *
 * MOMENTS WITHOUT INVERSION. Differentiating (*) at s = 0 gives
 * sum_i C(j,i) f_i M_{j-i} = [j = 1] v0. M_0 = D is SINGULAR with right null
 * vector e, so each order fixes f_j only up to a multiple of theta, and that
 * multiple is what the NEXT order's solvability condition supplies. At j = 0
 * the same condition reads theta M_1 e = v0 e, i.e. 1 - rho = 1 - rho, which is
 * the identity that validates the whole setup.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: the U iteration runs
 * to a tolerance, and both the deterministic transform and the inversion need
 * exp. The CDF is the Abate-Whitt Euler sum, which evaluates the transform OFF
 * the real axis, so a small complex layer is carried here rather than in the
 * distribution interface.
 *
 * Reference:
 *   Qi-Ming He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
 *   Queueing Systems 38(4):397-418, 2001.
 */

#include <cmath>
#include <complex>
#include <limits>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mmapgk1, mirroring the MATLAB struct. */
template <class T>
struct MmapGk1Result {
    std::vector<T> lambdas;                 ///< per-type arrival rates
    T arrivalRate;                          ///< sum of the per-type rates
    T utilization;                          ///< rho = sum_k lambda_k E[S_k]
    std::vector<T> idleVector;              ///< v0, summing to 1 - rho
    std::vector<std::vector<T> > waitMoments;  ///< [type][order], E[Wk^j]
    std::vector<T> meanWaitingTime;         ///< per-type E[Wq]
    std::vector<T> meanSojournTime;         ///< per-type E[Wq] + E[S]
    T meanQueueLength;                      ///< E[N] by Little over all types
    std::vector<std::vector<T> > waitCDF;   ///< [type][point], P(Wk <= t)
    std::vector<T> waitPoints;              ///< the requested points
};

namespace mmapgk1detail {

/** Binomial coefficient; declared ahead of the complex block that uses it. */
inline double binom(std::size_t n, std::size_t k);

/** Left null vector of G normalized to sum one. */
template <class T>
std::vector<T> stat_left_null(const Matrix<T>& G) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t n = G.rows();
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(j, i) = G(i, j);
    for (std::size_t j = 0; j < n; ++j) A(n - 1, j) = one;
    std::vector<T> b(n, zero);
    b[n - 1] = one;
    return line::solve(A, b);
}

/** Midpoint nodes and true CDF increments over the support of the law. */
template <class T>
void stieltjes_nodes(const lang::Distrib<T>& law, std::vector<T>& x, std::vector<T>& w) {
    const std::size_t n_grid = 2400;
    const T zero = num_traits<T>::from_int(0);
    const double mean = num_traits<T>::to_double(lang::dist_moment(law, 1u));
    double hi = mean * 60.0;
    const double var = num_traits<T>::to_double((lang::dist_moment(law, 2u) - lang::dist_moment(law, 1u) * lang::dist_moment(law, 1u)));
    if (std::isfinite(var) && var > 0.0) hi = std::max(hi, mean + 12.0 * std::sqrt(var));
    const double step = hi / static_cast<double>(n_grid);
    x.assign(n_grid, zero);
    w.assign(n_grid, zero);
    T mass = zero;
    T prev = lang::dist_cdf(law, num_traits<T>::from_double(0.0));
    for (std::size_t i = 0; i < n_grid; ++i) {
        const T right = num_traits<T>::from_double((i + 1) * step);
        const T cur = lang::dist_cdf(law, right);
        x[i] = num_traits<T>::from_double((i + 0.5) * step);
        w[i] = cur - prev;
        mass += w[i];
        prev = cur;
    }
    if (mass > zero)
        for (std::size_t i = 0; i < n_grid; ++i) w[i] = w[i] / mass;
}

/**
 * Raw moment E[S^j]; dist_moment already evaluates every family exactly.
 *
 * A HEAVY TAIL HAS NO MOMENT of high enough order, and dist_moment says so by
 * THROWING (a Pareto of shape <= j). Here that is not an error but an answer:
 * the moment recursion consumes M_j only from order j onwards, so an infinite
 * M_3 leaves E[Wq] finite and makes E[Wq^2] infinite, which is the truth about
 * such a queue. Letting the throw escape would lose the finite moments too.
 */
template <class T>
T raw_moment(const lang::Distrib<T>& law, std::size_t j) {
    try {
        return lang::dist_moment(law, static_cast<unsigned>(j));
    } catch (const NumericError&) {
        return num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    }
}

/** Kronecker product; C++ has no shared templated kron. */
template <class T>
Matrix<T> gk_kron(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows() * B.rows(), A.cols() * B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            for (std::size_t p = 0; p < B.rows(); ++p)
                for (std::size_t q = 0; q < B.cols(); ++q)
                    C(i * B.rows() + p, j * B.cols() + q) = A(i, j) * B(p, q);
    return C;
}

/** The matrix transform int_0^inf exp(U t) dF(t). */
template <class T>
Matrix<T> matrix_lst(const lang::Distrib<T>& law, const Matrix<T>& U) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t n = U.rows();
    if (lang::process_is_markovian(law.type)) {
        // The density is the SCALAR beta exp(St) s0, so the integral is exact on
        // the Kronecker sum: int exp(Ut) x exp(St) dt = -(U (+) S)^-1.
        const mam::Map<T> rep = lang::dist_to_map(law);
        const Matrix<T>& S = rep.D0;
        const std::size_t ms = S.rows();
        const std::vector<T> beta = mam::map_pie(rep);
        std::vector<T> s0(ms, zero);
        for (std::size_t i = 0; i < ms; ++i) {
            T r = zero;
            for (std::size_t j = 0; j < ms; ++j) r += S(i, j);
            s0[i] = -r;
        }
        Matrix<T> Ims = eye<T>(ms);
        Matrix<T> KS = gk_kron(U, Ims);
        const Matrix<T> In = eye<T>(n);
        const Matrix<T> IkS = gk_kron(In, S);
        for (std::size_t i = 0; i < KS.rows(); ++i)
            for (std::size_t j = 0; j < KS.cols(); ++j) KS(i, j) += IkS(i, j);
        Matrix<T> betaM(1, ms, zero);
        for (std::size_t i = 0; i < ms; ++i) betaM(0, i) = beta[i];
        Matrix<T> s0M(ms, 1, zero);
        for (std::size_t i = 0; i < ms; ++i) s0M(i, 0) = s0[i];
        const Matrix<T> L = gk_kron(In, betaM);
        Matrix<T> R = gk_kron(In, s0M);
        Matrix<T> X = matmul(inverse(KS), R);
        for (std::size_t i = 0; i < X.rows(); ++i)
            for (std::size_t j = 0; j < X.cols(); ++j) X(i, j) = -X(i, j);
        return matmul(L, X);
    }
    if (law.type == lang::ProcessType::DET) {
        Matrix<T> Ud = U;
        const T d = lang::dist_moment(law, 1u);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) Ud(i, j) = Ud(i, j) * d;
        return line::expm(Ud);
    }
    std::vector<T> x, w;
    stieltjes_nodes(law, x, w);
    Matrix<T> F(n, n, zero);
    for (std::size_t i = 0; i < x.size(); ++i) {
        Matrix<T> Ut = U;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) Ut(p, q) = Ut(p, q) * x[i];
        const Matrix<T> E = line::expm(Ut);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) F(p, q) += w[i] * E(p, q);
    }
    (void)one;
    return F;
}

typedef std::complex<double> cdbl;

/** Gaussian elimination with partial pivoting over the complex field. */
inline std::vector<cdbl> complex_solve(std::vector<std::vector<cdbl> > A, std::vector<cdbl> b) {
    const std::size_t n = b.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        double best = std::abs(A[col][col]);
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::abs(A[r][col]) > best) { best = std::abs(A[r][col]); piv = r; }
        }
        if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }
        for (std::size_t r = col + 1; r < n; ++r) {
            const cdbl f = A[r][col] / A[col][col];
            for (std::size_t c = col; c < n; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    std::vector<cdbl> x(n);
    for (std::size_t row = n; row-- > 0;) {
        cdbl s = b[row];
        for (std::size_t c = row + 1; c < n; ++c) s -= A[row][c] * x[c];
        x[row] = s / A[row][row];
    }
    return x;
}

/**
 * Scalar Laplace-Stieltjes transform of a service law at a complex argument.
 *
 * One line, because the LANGUAGE LAYER owns this: dist_lst carries a complex
 * overload beside the real one, with the closed forms, the phase-type solve and
 * the CDF-increment fallback already tiered. A second implementation here would
 * be a second thing to keep true.
 */
template <class T>
cdbl scalar_lst(const lang::Distrib<T>& law, cdbl s) {
    return lang::dist_lst(law, s);
}


/** E[exp(-s Wk)] at a complex argument. */
template <class T>
cdbl wait_lst(const Matrix<T>& D0, const std::vector<Matrix<T> >& Dk,
              const std::vector<lang::Distrib<T> >& svc, const std::vector<T>& v0,
              const std::vector<T>& lambdas, std::size_t k, cdbl s) {
    const std::size_t ma = D0.rows();
    std::vector<std::vector<cdbl> > M(ma, std::vector<cdbl>(ma, cdbl(0.0, 0.0)));
    for (std::size_t i = 0; i < ma; ++i)
        for (std::size_t j = 0; j < ma; ++j) M[i][j] = cdbl(num_traits<T>::to_double(D0(i, j)), 0.0);
    for (std::size_t i = 0; i < ma; ++i) M[i][i] += s;
    for (std::size_t q = 0; q < Dk.size(); ++q) {
        const cdbl g = scalar_lst(svc[q], s);
        for (std::size_t i = 0; i < ma; ++i)
            for (std::size_t j = 0; j < ma; ++j)
                M[i][j] += num_traits<T>::to_double(Dk[q](i, j)) * g;
    }
    // solve f M = s v0, i.e. M' f' = (s v0)'
    std::vector<std::vector<cdbl> > A(ma, std::vector<cdbl>(ma, cdbl(0.0, 0.0)));
    std::vector<cdbl> b(ma, cdbl(0.0, 0.0));
    for (std::size_t i = 0; i < ma; ++i) {
        for (std::size_t j = 0; j < ma; ++j) A[i][j] = M[j][i];
        b[i] = s * num_traits<T>::to_double(v0[i]);
    }
    const std::vector<cdbl> f = complex_solve(A, b);
    cdbl out(0.0, 0.0);
    for (std::size_t i = 0; i < ma; ++i) {
        double rowsum = 0.0;
        for (std::size_t j = 0; j < ma; ++j) rowsum += num_traits<T>::to_double(Dk[k](i, j));
        out += f[i] * rowsum;
    }
    return out / num_traits<T>::to_double(lambdas[k]);
}

/** Abate-Whitt Euler inversion of the type-k waiting time CDF. */
template <class T>
double euler_invert(const Matrix<T>& D0, const std::vector<Matrix<T> >& Dk,
                    const std::vector<lang::Distrib<T> >& svc, const std::vector<T>& v0,
                    const std::vector<T>& lambdas, std::size_t k, double t) {
    if (t <= 0.0) return wait_lst(D0, Dk, svc, v0, lambdas, k, cdbl(1e12, 0.0)).real();
    const double A = 18.4;
    const std::size_t nE = 15, mE = 11;
    const double u = std::exp(A / 2) / t;
    const double x = A / (2 * t);
    std::vector<double> terms(nE + mE + 1, 0.0);
    terms[0] = wait_lst(D0, Dk, svc, v0, lambdas, k, cdbl(x, 0.0)).real() / x / 2.0;
    for (std::size_t j = 1; j <= nE + mE; ++j) {
        const cdbl s(x, M_PI * static_cast<double>(j) / t);
        terms[j] = ((j % 2 == 0) ? 1.0 : -1.0) * (wait_lst(D0, Dk, svc, v0, lambdas, k, s) / s).real();
    }
    std::vector<double> partial(terms.size(), 0.0);
    double run = 0.0;
    for (std::size_t j = 0; j < terms.size(); ++j) { run += terms[j]; partial[j] = run; }
    double F = 0.0;
    for (std::size_t j = 0; j <= mE; ++j) F += binom(mE, j) / std::pow(2.0, static_cast<double>(mE)) * partial[nE + j];
    F = u * F;
    return std::min(std::max(F, 0.0), 1.0);
}


inline double binom(std::size_t n, std::size_t k) {
    double r = 1.0;
    for (std::size_t i = 1; i <= k; ++i) r = r * static_cast<double>(n - k + i) / static_cast<double>(i);
    return std::floor(r + 0.5);
}

}  // namespace mmapgk1detail

/**
 * MMAP[K]/G[K]/1 FCFS, per type.
 *
 * @param MMAP       LINE convention {D0, D1, D^(1), ..., D^(K)}, D1 = sum_k D^(k)
 * @param svc        K service laws, one per marked type; families may differ
 * @param w_points   times at which to evaluate the per-type waiting time CDF
 * @param num_w_moms how many per-type waiting time moments to return
 * @param tol        fixed point tolerance on U
 * @param iter_max   fixed point iteration cap
 */
template <class T>
MmapGk1Result<T> qsys_mmapgk1(const std::vector<Matrix<T> >& MMAP,
                              const std::vector<lang::Distrib<T> >& svc,
                              const std::vector<T>& w_points, std::size_t num_w_moms,
                              double tol, std::size_t iter_max) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mmapgk1 requires transcendental arithmetic");
    using namespace mmapgk1detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (MMAP.size() < 3) throw InputError("qsys_mmapgk1: the MMAP must carry a marked block");
    const std::size_t K = MMAP.size() - 2;
    if (svc.size() != K) throw InputError("qsys_mmapgk1: one service law per marked type");
    const Matrix<T>& D0 = MMAP[0];
    const std::size_t ma = D0.rows();
    std::vector<Matrix<T> > Dk;
    Matrix<T> Dsum = D0;
    for (std::size_t k = 0; k < K; ++k) {
        Dk.push_back(MMAP[k + 2]);
        for (std::size_t i = 0; i < ma; ++i)
            for (std::size_t j = 0; j < ma; ++j) Dsum(i, j) += Dk[k](i, j);
    }

    const std::vector<T> theta = stat_left_null(Dsum);
    std::vector<T> lambdas(K, zero), mean_s(K, zero);
    T rho = zero;
    for (std::size_t k = 0; k < K; ++k) {
        T lam = zero;
        for (std::size_t i = 0; i < ma; ++i)
            for (std::size_t j = 0; j < ma; ++j) lam += theta[i] * Dk[k](i, j);
        lambdas[k] = lam;
        mean_s[k] = lang::dist_moment(svc[k], 1u);
        rho += lam * mean_s[k];
    }
    if (rho >= one) throw InputError("qsys_mmapgk1: load rho must be strictly less than 1");

    Matrix<T> U = D0;
    for (std::size_t it = 0; it < iter_max; ++it) {
        Matrix<T> Unew = D0;
        for (std::size_t k = 0; k < K; ++k) {
            const Matrix<T> Fk = matrix_lst(svc[k], U);
            const Matrix<T> P = matmul(Dk[k], Fk);
            for (std::size_t i = 0; i < ma; ++i)
                for (std::size_t j = 0; j < ma; ++j) Unew(i, j) += P(i, j);
        }
        double diff = 0.0;
        for (std::size_t i = 0; i < ma; ++i)
            for (std::size_t j = 0; j < ma; ++j)
                diff = std::max(diff, std::fabs(num_traits<T>::to_double(Unew(i, j) - U(i, j))));
        U = Unew;
        if (diff <= tol) break;
    }
    // U e = 0 EXACTLY: a property of the fixed point, not of the iterate, and
    // the residue at the tolerance above would move the stationary solve.
    for (std::size_t i = 0; i < ma; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < ma; ++j) rowsum += U(i, j);
        U(i, i) -= rowsum;
    }
    std::vector<T> v0 = stat_left_null(U);
    for (std::size_t i = 0; i < ma; ++i) v0[i] = v0[i] * (one - rho);

    // Moment recursion, see the header comment
    std::vector<Matrix<T> > Mder;
    for (std::size_t j = 0; j <= num_w_moms + 1; ++j) {
        Matrix<T> Mj(ma, ma, zero);
        if (j == 0) {
            Mj = D0;
            for (std::size_t k = 0; k < K; ++k)
                for (std::size_t i = 0; i < ma; ++i)
                    for (std::size_t q = 0; q < ma; ++q) Mj(i, q) += Dk[k](i, q);
        } else {
            for (std::size_t k = 0; k < K; ++k) {
                const T mom = raw_moment(svc[k], j);
                const double sign = (j % 2 == 0) ? 1.0 : -1.0;
                for (std::size_t i = 0; i < ma; ++i)
                    for (std::size_t q = 0; q < ma; ++q)
                        Mj(i, q) += num_traits<T>::from_double(sign) * Dk[k](i, q) * mom;
            }
            if (j == 1)
                for (std::size_t i = 0; i < ma; ++i) Mj(i, i) += one;
        }
        Mder.push_back(Mj);
    }
    const std::vector<T> e(ma, one);
    T denom = zero;
    {
        const std::vector<T> v = mulvec(Mder[1], e);
        for (std::size_t i = 0; i < ma; ++i) denom += theta[i] * v[i];
    }
    std::vector<std::vector<T> > fder;
    fder.push_back(theta);
    // [M_0, e] with f_j^p e = 0 pins the particular solution; solved in the
    // least squares sense, since the system is one equation over-determined.
    Matrix<T> Abase(ma + 1, ma, zero);
    for (std::size_t i = 0; i < ma; ++i)
        for (std::size_t j = 0; j < ma; ++j) Abase(j, i) = Mder[0](i, j);
    for (std::size_t i = 0; i < ma; ++i) Abase(ma, i) = one;
    Matrix<T> AtA(ma, ma, zero);
    for (std::size_t i = 0; i < ma; ++i)
        for (std::size_t j = 0; j < ma; ++j) {
            T s = zero;
            for (std::size_t q = 0; q < ma + 1; ++q) s += Abase(q, i) * Abase(q, j);
            AtA(i, j) = s;
        }
    const Matrix<T> AtAinv = inverse(AtA);
    for (std::size_t j = 1; j <= num_w_moms; ++j) {
        std::vector<T> rhs(ma, zero);
        if (j == 1) for (std::size_t i = 0; i < ma; ++i) rhs[i] = v0[i];
        for (std::size_t i = 0; i < j; ++i) {
            const double c = binom(j, i);
            for (std::size_t q = 0; q < ma; ++q) {
                T acc = zero;
                for (std::size_t p = 0; p < ma; ++p) acc += fder[i][p] * Mder[j - i](p, q);
                rhs[q] -= num_traits<T>::from_double(c) * acc;
            }
        }
        std::vector<T> rhsAug(ma + 1, zero);
        for (std::size_t q = 0; q < ma; ++q) rhsAug[q] = rhs[q];
        std::vector<T> Atb(ma, zero);
        for (std::size_t i = 0; i < ma; ++i) {
            T s = zero;
            for (std::size_t q = 0; q < ma + 1; ++q) s += Abase(q, i) * rhsAug[q];
            Atb[i] = s;
        }
        const std::vector<T> fp = mulvec(AtAinv, Atb);
        T acc2 = zero;
        for (std::size_t i = 0; i < j; ++i) {
            const std::vector<T> v = mulvec(Mder[j + 1 - i], e);
            T inner = zero;
            for (std::size_t p = 0; p < ma; ++p) inner += fder[i][p] * v[p];
            acc2 += num_traits<T>::from_double(binom(j + 1, i)) * inner;
        }
        T fpM1e = zero;
        {
            const std::vector<T> v = mulvec(Mder[1], e);
            for (std::size_t p = 0; p < ma; ++p) fpM1e += fp[p] * v[p];
        }
        const T cfree = (-acc2 / num_traits<T>::from_int(static_cast<long>(j + 1)) - fpM1e) / denom;
        std::vector<T> fj(ma, zero);
        for (std::size_t p = 0; p < ma; ++p) fj[p] = fp[p] + cfree * theta[p];
        fder.push_back(fj);
    }

    MmapGk1Result<T> r;
    r.lambdas = lambdas;
    r.arrivalRate = zero;
    for (std::size_t k = 0; k < K; ++k) r.arrivalRate += lambdas[k];
    r.utilization = rho;
    r.idleVector = v0;
    r.waitMoments.assign(K, std::vector<T>(num_w_moms, zero));
    r.meanWaitingTime.assign(K, zero);
    r.meanSojournTime.assign(K, zero);
    r.meanQueueLength = zero;
    for (std::size_t k = 0; k < K; ++k) {
        for (std::size_t j = 1; j <= num_w_moms; ++j) {
            const std::vector<T> v = mulvec(Dk[k], e);
            T inner = zero;
            for (std::size_t p = 0; p < ma; ++p) inner += fder[j][p] * v[p];
            const double sign = (j % 2 == 0) ? 1.0 : -1.0;
            r.waitMoments[k][j - 1] = num_traits<T>::from_double(sign) * inner / lambdas[k];
        }
        r.meanWaitingTime[k] = num_w_moms >= 1 ? r.waitMoments[k][0] : zero;
        r.meanSojournTime[k] = r.meanWaitingTime[k] + mean_s[k];
        r.meanQueueLength += lambdas[k] * r.meanSojournTime[k];
    }

    r.waitPoints = w_points;
    if (!w_points.empty()) {
        r.waitCDF.assign(K, std::vector<T>(w_points.size(), zero));
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t it = 0; it < w_points.size(); ++it) {
                const double t = num_traits<T>::to_double(w_points[it]);
                r.waitCDF[k][it] = num_traits<T>::from_double(
                    euler_invert(D0, Dk, svc, v0, lambdas, k, t));
            }
        }
    }
    return r;
}


template <class T>
MmapGk1Result<T> qsys_mmapgk1(const std::vector<Matrix<T> >& MMAP,
                              const std::vector<lang::Distrib<T> >& svc) {
    return qsys_mmapgk1(MMAP, svc, std::vector<T>(), static_cast<std::size_t>(3), 1e-12,
                        static_cast<std::size_t>(10000));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MMAPGK1_H
