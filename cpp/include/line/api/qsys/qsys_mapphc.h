/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPPHC_H
#define LINE_API_QSYS_QSYS_MAPPHC_H

/**
 * The MAP/PH/c FCFS queue, solved exactly.
 *
 * THE STATE SPACE, AND WHY IT IS A MULTISET. With c identical servers the
 * server identities carry no information, so the service phases are held as a
 * MULTISET: a configuration is n = (n_1..n_ms) with sum(n) = k servers busy in
 * phase i. There are binomial(ms+k-1,k) of those, the count of Asmussen and
 * Moller (2001), against ms^k for the ordered space -- for ms = 5, c = 6 that
 * is 210 against 15625, which is what makes the exact solve feasible at all.
 * Levels 0..c-1 are the boundary (level = servers busy), levels >= c repeat and
 * carry the queue, so the tail is matrix-geometric in R.
 *
 * THE WAITING TIME, AND WHY IT IS ONE LINEAR ODE. An arrival that finds j
 * customers waiting ahead of it waits for exactly j+1 service completions, so
 * Wq is the (j+1)-st event time of the configuration MAP (Lc, Cdep) started at
 * the arrival-epoch configuration. The level distribution seen by an arrival is
 * matrix-geometric, x_j = pi_c R^j kron(D1,I)/lambda, and folding over j gives
 *
 *     G'(t) = G(t) Lj + R G(t) Cj,   G(0) = (I-R)^-1 kron(D1,I)/lambda,
 *     P(Wq > t) = pi_c G(t) e.
 *
 * That is LINEAR in G, so Wq is matrix-exponential; vectorizing it column-major
 * turns it into a single expm of order (ma*nc)^2. Its Laplace transform obeys
 * the generalized Sylvester equation g(sI-Lj) - R g Cj = G(0), and
 * differentiating that identity leaves the OPERATOR unchanged and only moves
 * the right-hand side, so every moment is one more solve with the same matrix.
 *
 * WHY NOT REUSE qsys_mapmc. That function is MAP/M/c: its phase is the arrival
 * phase alone, because an exponential server has nothing to remember. Here the
 * phase must also carry which service phase each busy server occupies, and the
 * down block is a completion followed by an IMMEDIATE restart at alpha, which
 * has no counterpart in the exponential case. The two agree exactly when
 * ms = 1, and the test file checks that they do.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: logarithmic reduction
 * drives R to a tolerance and never terminates in a finite number of field
 * operations, as in qbd_r.h. Everything consuming R is finite exact matrix
 * algebra, except the CCDF, which needs expm.
 *
 * References:
 *   S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
 *   distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
 *   2001.
 *   D.P. Gaver, P.A. Jacobs, G. Latouche, "Finite birth-and-death models in
 *   randomly changing environments", Adv. Appl. Probab. 16:715-731, 1984.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/ldqbd_mphc.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mapphc, mirroring the MATLAB struct. */
template <class T>
struct MapPhcResult {
    T meanQueueLength;                ///< E[N], number in system
    T meanWaitingTime;                ///< E[Wq], time in queue
    T meanSojournTime;                ///< E[Wq] + E[service]
    T utilization;                    ///< rho = lambda E[service] / c, per server
    std::vector<T> queueLengthDist;   ///< P(N = n), n = 0, 1, ...
    std::vector<T> waitingTimeMoments;///< E[Wq^k], k = 1..num_w_moms
    std::vector<T> waitingTimeCCDF;   ///< P(Wq > t) at the requested points
    std::vector<T> waitingTimePoints; ///< the requested points
    T probWait;                       ///< P(Wq > 0), an arrival finds every server busy
    std::size_t phaseCount;           ///< binomial(ms+c-1,c), the repeating config count
};

namespace mapphcdetail {

/** Elementwise A + B; the qbd_detail twin is not reachable from here. */
template <class T>
Matrix<T> madd2(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw InputError("qsys_mapphc: shape mismatch");
    Matrix<T> C = A;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) += B(i, j);
    return C;
}

/** Kronecker product; C++ has no shared templated kron. */
template <class T>
Matrix<T> mkron(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows() * B.rows(), A.cols() * B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            for (std::size_t p = 0; p < B.rows(); ++p)
                for (std::size_t q = 0; q < B.cols(); ++q)
                    C(i * B.rows() + p, j * B.cols() + q) = A(i, j) * B(p, q);
    return C;
}

/**
 * Compositions of k into ms nonnegative parts, in a fixed order. Shared with the
 * exact M/PH/c LD-QBD blocks, so a configuration index means the same thing in
 * both.
 */
inline std::vector<std::vector<int> > multisets(std::size_t ms, std::size_t k) {
    return mam::ph_multisets(ms, k);
}

inline std::size_t find_cfg(const std::vector<std::vector<int> >& rows, const std::vector<int>& key) {
    for (std::size_t i = 0; i < rows.size(); ++i)
        if (rows[i] == key) return i;
    throw InputError("qsys_mapphc: configuration not found");
}

}  // namespace mapphcdetail

/**
 * MAP/PH/c FCFS, exactly.
 *
 * @param arrival     arrival MAP (D0, D1) of order ma
 * @param alpha       PH service initial vector of order ms
 * @param S           PH service sub-generator of order ms
 * @param c           number of servers, c >= 1
 * @param dist_size   cap on the queue length probabilities materialized
 * @param num_w_moms  how many waiting-time moments to return
 * @param w_points    times at which to evaluate P(Wq > t)
 */
template <class T>
MapPhcResult<T> qsys_mapphc(const mam::Map<T>& arrival, const std::vector<T>& alpha,
                            const Matrix<T>& S, unsigned c, std::size_t dist_size,
                            std::size_t num_w_moms, const std::vector<T>& w_points) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapphc requires transcendental arithmetic");
    using mapphcdetail::find_cfg;
    using mapphcdetail::madd2;
    using mapphcdetail::mkron;
    using mapphcdetail::multisets;

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const Matrix<T>& D0 = arrival.D0;
    const Matrix<T>& D1 = arrival.D1;
    const std::size_t ma = D0.rows();
    const std::size_t ms = S.rows();
    if (D0.cols() != ma || D1.rows() != ma || D1.cols() != ma)
        throw InputError("qsys_mapphc: D0 and D1 must be square and of equal order");
    if (S.cols() != ms || alpha.size() != ms)
        throw InputError("qsys_mapphc: alpha and S must have matching order");
    if (c < 1) throw InputError("qsys_mapphc: c must be a positive integer");
    if (dist_size == 0) throw InputError("qsys_mapphc: dist_size must be positive");

    std::vector<T> s0(ms, zero);
    for (std::size_t i = 0; i < ms; ++i) {
        T r = zero;
        for (std::size_t j = 0; j < ms; ++j) r += S(i, j);
        s0[i] = -r;
    }

    const T lambda = mam::map_lambda(arrival);
    if (lambda <= zero) throw InputError("qsys_mapphc: non-positive arrival rate");
    // E[service] = -alpha S^-1 e
    const std::vector<T> negSinv_e = line::solve(S, ones<T>(ms));
    T mean_service = zero;
    for (std::size_t i = 0; i < ms; ++i) mean_service -= alpha[i] * negSinv_e[i];
    const T ct = num_traits<T>::from_int(static_cast<long>(c));
    const T rho = lambda * mean_service / ct;
    if (rho >= one) throw InputError("qsys_mapphc: load rho must be strictly less than 1");

    std::vector<std::vector<std::vector<int> > > cfg(c + 1);
    for (unsigned k = 0; k <= c; ++k) cfg[k] = multisets(ms, k);

    std::vector<Matrix<T> > Lcfg, Up, Dn;
    for (unsigned k = 0; k <= c; ++k) {
        const std::vector<std::vector<int> >& Ck = cfg[k];
        const std::size_t nk = Ck.size();
        Matrix<T> Lk(nk, nk, zero);
        for (std::size_t row = 0; row < nk; ++row) {
            const std::vector<int>& n = Ck[row];
            for (std::size_t i = 0; i < ms; ++i) {
                if (n[i] == 0) continue;
                const T ni = num_traits<T>::from_int(n[i]);
                for (std::size_t j = 0; j < ms; ++j) {
                    if (j == i) continue;
                    std::vector<int> m = n;
                    --m[i];
                    ++m[j];
                    Lk(row, find_cfg(Ck, m)) += ni * S(i, j);
                }
                Lk(row, row) += ni * S(i, i);
            }
        }
        Lcfg.push_back(Lk);

        if (k < c) {
            const std::vector<std::vector<int> >& Ck1 = cfg[k + 1];
            Matrix<T> Uk(nk, Ck1.size(), zero);
            for (std::size_t row = 0; row < nk; ++row) {
                for (std::size_t j = 0; j < ms; ++j) {
                    std::vector<int> m = Ck[row];
                    ++m[j];
                    Uk(row, find_cfg(Ck1, m)) += alpha[j];
                }
            }
            Up.push_back(Uk);
        } else {
            Up.push_back(Matrix<T>(1, 1, zero));
        }

        if (k > 0) {
            const std::vector<std::vector<int> >& Ckm = cfg[k - 1];
            Matrix<T> Dk(nk, Ckm.size(), zero);
            for (std::size_t row = 0; row < nk; ++row) {
                const std::vector<int>& n = Ck[row];
                for (std::size_t i = 0; i < ms; ++i) {
                    if (n[i] == 0) continue;
                    std::vector<int> m = n;
                    --m[i];
                    Dk(row, find_cfg(Ckm, m)) += num_traits<T>::from_int(n[i]) * s0[i];
                }
            }
            Dn.push_back(Dk);
        } else {
            Dn.push_back(Matrix<T>(1, 1, zero));
        }
    }

    // Completion WITH an immediate restart: the repeating down block
    const std::vector<std::vector<int> >& Cc = cfg[c];
    const std::size_t nc = Cc.size();
    Matrix<T> Cdep(nc, nc, zero);
    for (std::size_t row = 0; row < nc; ++row) {
        const std::vector<int>& n = Cc[row];
        for (std::size_t i = 0; i < ms; ++i) {
            if (n[i] == 0) continue;
            const T ni = num_traits<T>::from_int(n[i]);
            for (std::size_t j = 0; j < ms; ++j) {
                std::vector<int> m = n;
                --m[i];
                ++m[j];
                Cdep(row, find_cfg(Cc, m)) += ni * s0[i] * alpha[j];
            }
        }
    }

    const Matrix<T> Ima = eye<T>(ma);
    const Matrix<T> Inc = eye<T>(nc);
    const Matrix<T> A_up = mkron(D1, Inc);
    const Matrix<T> A_loc = madd2(mkron(D0, Inc), mkron(Ima, Lcfg[c]));
    const Matrix<T> A_dn = mkron(Ima, Cdep);
    const Matrix<T> R = mam::qbd_R_logred(A_dn, A_loc, A_up);

    const std::size_t n_op = nc * ma;
    Matrix<T> ImR(n_op, n_op, zero);
    for (std::size_t i = 0; i < n_op; ++i)
        for (std::size_t j = 0; j < n_op; ++j) ImR(i, j) = (i == j ? one : zero) - R(i, j);
    const Matrix<T> ImRinv = inverse(ImR);
    const std::vector<T> sum_geom = mulvec(ImRinv, ones<T>(n_op));

    // Boundary levels 0..c, with the tail folded into level c through R
    std::vector<std::size_t> sz(c + 1), off(c + 2, 0);
    for (unsigned k = 0; k <= c; ++k) {
        sz[k] = ma * cfg[k].size();
        off[k + 1] = off[k] + sz[k];
    }
    const std::size_t tot = off[c + 1];
    Matrix<T> Q(tot, tot, zero);
    const Matrix<T> RA_dn = matmul(R, A_dn);
    for (unsigned k = 0; k <= c; ++k) {
        const Matrix<T> Ick = eye<T>(cfg[k].size());
        const Matrix<T> diag = (k < c) ? madd2(mkron(D0, Ick), mkron(Ima, Lcfg[k]))
                                       : madd2(A_loc, RA_dn);
        for (std::size_t i = 0; i < sz[k]; ++i)
            for (std::size_t j = 0; j < sz[k]; ++j) Q(off[k] + i, off[k] + j) += diag(i, j);
        if (k < c) {
            const Matrix<T> up = mkron(D1, Up[k]);
            for (std::size_t i = 0; i < up.rows(); ++i)
                for (std::size_t j = 0; j < up.cols(); ++j) Q(off[k] + i, off[k + 1] + j) += up(i, j);
        }
        if (k > 0) {
            const Matrix<T> dn = mkron(Ima, Dn[k]);
            for (std::size_t i = 0; i < dn.rows(); ++i)
                for (std::size_t j = 0; j < dn.cols(); ++j) Q(off[k] + i, off[k - 1] + j) += dn(i, j);
        }
    }

    // pi Q = 0 with the normalization replacing the last equation, transposed
    // into the column-vector solve line::solve expects.
    Matrix<T> M(tot, tot, zero);
    for (std::size_t i = 0; i < tot; ++i)
        for (std::size_t j = 0; j < tot; ++j) M(j, i) = Q(i, j);
    for (std::size_t col = 0; col < tot; ++col) M(tot - 1, col) = zero;
    for (std::size_t i = 0; i < off[c]; ++i) M(tot - 1, i) = one;
    for (std::size_t i = 0; i < sz[c]; ++i) M(tot - 1, off[c] + i) = sum_geom[i];
    std::vector<T> b(tot, zero);
    b[tot - 1] = one;
    const std::vector<T> pi_vec = line::solve(M, b);

    std::vector<T> pi_c(n_op);
    for (std::size_t i = 0; i < n_op; ++i) pi_c[i] = pi_vec[off[c] + i];

    // Queue length distribution
    std::vector<T> ql;
    for (unsigned k = 0; k < c; ++k) {
        T s = zero;
        for (std::size_t i = 0; i < sz[k]; ++i) s += pi_vec[off[k] + i];
        ql.push_back(s);
    }
    std::vector<T> tail = pi_c;
    T acc = zero;
    for (std::size_t i = 0; i < ql.size(); ++i) acc += ql[i];
    {
        T s = zero;
        for (std::size_t i = 0; i < n_op; ++i) s += tail[i];
        ql.push_back(s);
        acc += s;
    }
    while (acc < one - num_traits<T>::from_double(1e-12) && ql.size() < dist_size) {
        std::vector<T> next(n_op, zero);
        for (std::size_t j = 0; j < n_op; ++j)
            for (std::size_t i = 0; i < n_op; ++i) next[j] += tail[i] * R(i, j);
        tail = next;
        T s = zero;
        for (std::size_t i = 0; i < n_op; ++i) s += tail[i];
        ql.push_back(s);
        acc += s;
    }
    // E[N] in CLOSED FORM. dist_size caps the probabilities RETURNED, not the mean:
    // summing the truncated list loses the matrix-geometric tail, which at rho -> 1
    // carries a first-order share of the mass. With pi_{c+j} = pi_c R^j,
    // sum_j (c+j) pi_c R^j e = pi_c [c (I-R)^-1 + R (I-R)^-2] e.
    T meanQL = zero;
    for (unsigned k = 0; k < c; ++k)
        meanQL += num_traits<T>::from_int(static_cast<long>(k)) * ql[k];
    const std::vector<T> u_tail = mulvec(ImRinv, ones<T>(n_op));
    const std::vector<T> r_tail = mulvec(R, mulvec(ImRinv, u_tail));
    const T c_scal = num_traits<T>::from_int(static_cast<long>(c));
    for (std::size_t i = 0; i < n_op; ++i)
        meanQL += pi_c[i] * (c_scal * u_tail[i] + r_tail[i]);

    // Waiting time
    const Matrix<T> Lj = mkron(Ima, Lcfg[c]);
    const Matrix<T> Cj = mkron(Ima, Cdep);
    Matrix<T> G0 = matmul(ImRinv, mkron(D1, Inc));
    for (std::size_t i = 0; i < n_op; ++i)
        for (std::size_t j = 0; j < n_op; ++j) G0(i, j) = G0(i, j) / lambda;
    T probWait = zero;
    {
        const std::vector<T> v = mulvec(G0, ones<T>(n_op));
        for (std::size_t i = 0; i < n_op; ++i) probWait += pi_c[i] * v[i];
    }

    // X(-Lj) - R X Cj = rhs, vectorized column-major
    const std::size_t n2 = n_op * n_op;
    Matrix<T> Kop(n2, n2, zero);
    for (std::size_t a = 0; a < n_op; ++a) {
        for (std::size_t bcol = 0; bcol < n_op; ++bcol) {
            // contribution of X(-Lj): entry (i,a) gets -Lj(bcol,a) X(i,bcol)
            for (std::size_t i = 0; i < n_op; ++i)
                Kop(a * n_op + i, bcol * n_op + i) += -Lj(bcol, a);
            // contribution of -R X Cj: entry (i,a) gets -R(i,p) X(p,bcol) Cj(bcol,a)
            for (std::size_t i = 0; i < n_op; ++i)
                for (std::size_t p = 0; p < n_op; ++p)
                    Kop(a * n_op + i, bcol * n_op + p) += -R(i, p) * Cj(bcol, a);
        }
    }
    std::vector<T> wMoms;
    Matrix<T> gPrev(n_op, n_op, zero);
    for (std::size_t k = 1; k <= num_w_moms; ++k) {
        std::vector<T> rhs(n2, zero);
        for (std::size_t j = 0; j < n_op; ++j)
            for (std::size_t i = 0; i < n_op; ++i)
                rhs[j * n_op + i] = (k == 1)
                        ? G0(i, j)
                        : -num_traits<T>::from_int(static_cast<long>(k - 1)) * gPrev(i, j);
        const std::vector<T> gv = line::solve(Kop, rhs);
        Matrix<T> g(n_op, n_op, zero);
        for (std::size_t j = 0; j < n_op; ++j)
            for (std::size_t i = 0; i < n_op; ++i) g(i, j) = gv[j * n_op + i];
        const std::vector<T> v = mulvec(g, ones<T>(n_op));
        T mk = zero;
        for (std::size_t i = 0; i < n_op; ++i) mk += pi_c[i] * v[i];
        mk = mk * num_traits<T>::from_int(static_cast<long>(k));
        if (k % 2 == 0) mk = -mk;
        wMoms.push_back(mk);
        gPrev = g;
    }

    std::vector<T> wCCDF;
    if (!w_points.empty()) {
        Matrix<T> Kt(n2, n2, zero);
        for (std::size_t a = 0; a < n_op; ++a) {
            for (std::size_t bcol = 0; bcol < n_op; ++bcol) {
                for (std::size_t i = 0; i < n_op; ++i)
                    Kt(a * n_op + i, bcol * n_op + i) += Lj(bcol, a);
                for (std::size_t i = 0; i < n_op; ++i)
                    for (std::size_t p = 0; p < n_op; ++p)
                        Kt(a * n_op + i, bcol * n_op + p) += R(i, p) * Cj(bcol, a);
            }
        }
        std::vector<T> v0(n2, zero);
        for (std::size_t j = 0; j < n_op; ++j)
            for (std::size_t i = 0; i < n_op; ++i) v0[j * n_op + i] = G0(i, j);
        for (std::size_t it = 0; it < w_points.size(); ++it) {
            Matrix<T> Kts(n2, n2, zero);
            for (std::size_t i = 0; i < n2; ++i)
                for (std::size_t j = 0; j < n2; ++j) Kts(i, j) = Kt(i, j) * w_points[it];
            const std::vector<T> vt = mulvec(line::expm(Kts), v0);
            Matrix<T> Gt(n_op, n_op, zero);
            for (std::size_t j = 0; j < n_op; ++j)
                for (std::size_t i = 0; i < n_op; ++i) Gt(i, j) = vt[j * n_op + i];
            const std::vector<T> v = mulvec(Gt, ones<T>(n_op));
            T s = zero;
            for (std::size_t i = 0; i < n_op; ++i) s += pi_c[i] * v[i];
            wCCDF.push_back(s);
        }
    }

    MapPhcResult<T> r;
    r.meanQueueLength = meanQL;
    r.meanWaitingTime = wMoms.empty() ? zero : wMoms[0];
    r.meanSojournTime = r.meanWaitingTime + mean_service;
    r.utilization = rho;
    r.queueLengthDist = ql;
    r.waitingTimeMoments = wMoms;
    r.waitingTimeCCDF = wCCDF;
    r.waitingTimePoints = w_points;
    r.probWait = probWait;
    r.phaseCount = nc;
    return r;
}

template <class T>
MapPhcResult<T> qsys_mapphc(const mam::Map<T>& arrival, const std::vector<T>& alpha,
                            const Matrix<T>& S, unsigned c) {
    return qsys_mapphc(arrival, alpha, S, c, static_cast<std::size_t>(500),
                       static_cast<std::size_t>(3), std::vector<T>());
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPPHC_H
