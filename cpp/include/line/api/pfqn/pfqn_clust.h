/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CLUST_H
#define LINE_API_PFQN_CLUST_H

/**
 * de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_clust.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_clust.java. E. de Souza e Silva,
 * S. S. Lavenberg, R. R. Muntz, "A clustering approximation technique for
 * queueing network models with a large number of chains", IEEE Trans.
 * Computers C-35(5), 1986. The network is covered by subnetworks whose union
 * is the whole network but which need not be disjoint. Every class visiting a
 * subnetwork S is either LOCAL to S, and is then solved inside it, or FOREIGN,
 * and is then seen only through the utilization it leaves behind. Each
 * subnetwork is solved by an ordinary approximate MVA algorithm with two
 * replacements: the complement of S is collapsed into a per-class delay P_c
 * and the foreign classes into a per-centre utilization U_k,
 *
 *   X_c(N) = N_c / (sum_{k in S} R_ck(N) + Z_c + P_c),
 *   Q_k(N) = [sum_{c in LC(S)} R_ck(N) X_c(N) + U_k] / (1 - U_k).
 *
 * Choosing the PE algorithm for every subnetwork reproduces global PE exactly,
 * so the useful setting is Linearizer inside, PE outside: the cost then sits
 * between pfqn_bs and pfqn_linearizer, which is the point of the method.
 *
 * When no decomposition is supplied the criterion of the paper is applied
 * automatically: the cheap PAMB estimate (pfqn_pam) of the centre utilizations
 * is taken, every class is attached to the centre where it loads the most,
 * classes sharing that centre form one cluster, and the subnetwork of a
 * cluster is the set of centres its classes visit. The answer depends on the
 * decomposition, which is why it is an input.
 *
 * The name avoids pfqn_ca, which is the exact convolution algorithm.
 *
 * Arithmetic: field operations only, so each iterate is EXACT in rational
 * arithmetic; both the outer and the inner loop stop on a tolerance.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_pam.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Which algorithm runs inside a subnetwork. */
enum class ClustInner { Linearizer, ProportionalEstimation };

namespace detail {

template <class T>
struct SubnetSolution {
    std::vector<T> X;
    Matrix<T> Q;
};

/** One forward MVA sweep over the local classes, with the foreign share folded in. */
template <class T>
SubnetSolution<T> clust_forward(const Matrix<T>& L, const std::vector<T>& N,
                                const std::vector<T>& Z, const std::vector<T>& Uk,
                                const Matrix<T>& Q,
                                const std::vector<std::vector<std::vector<T> > >& Delta) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    SubnetSolution<T> out;
    out.X.assign(R, zero);
    out.Q = Matrix<T>(M, R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        if (!(N[r] > zero)) continue;
        std::vector<T> W(M, zero);
        T wsum = zero;
        for (std::size_t i = 0; i < M; ++i) {
            // Linearizer estimate of the local queue at N - 1_r
            T qm = zero;
            for (std::size_t s = 0; s < R; ++s) {
                const T ns = (s == r) ? T(N[s] - one) : N[s];
                if (N[s] > zero && ns > zero) qm += ns * T(Q(i, s) / N[s] + Delta[i][s][r]);
            }
            if (qm < zero) qm = zero;
            const T A = T(T(qm + Uk[i]) / T(one - Uk[i]));
            W[i] = L(i, r) * T(one + A);
            wsum += W[i];
        }
        const T den = T(Z[r] + wsum);
        if (den == zero) throw NumericError("pfqn_clust: zero subnetwork cycle time");
        out.X[r] = N[r] / den;
        for (std::size_t i = 0; i < M; ++i) out.Q(i, r) = out.X[r] * W[i];
    }
    return out;
}

template <class T>
Matrix<T> clust_core(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                     const std::vector<T>& Uk, const Matrix<T>& Qin,
                     const std::vector<std::vector<std::vector<T> > >& Delta, double tol,
                     std::size_t maxiter) {
    Matrix<T> Q = Qin;
    for (std::size_t it = 0; it < maxiter; ++it) {
        const Matrix<T> Qold = Q;
        Q = clust_forward(L, N, Z, Uk, Q, Delta).Q;
        double d = 0.0;
        for (std::size_t i = 0; i < Q.rows(); ++i)
            for (std::size_t j = 0; j < Q.cols(); ++j) {
                const double v = std::fabs(num_traits<T>::to_double(T(Q(i, j) - Qold(i, j))));
                if (v > d) d = v;
            }
        if (d < tol) break;
    }
    return Q;
}

/** Approximate MVA restricted to the local classes of one subnetwork. */
template <class T>
SubnetSolution<T> clust_subnet(const Matrix<T>& L, const std::vector<T>& N,
                               const std::vector<T>& Z, const std::vector<T>& Uk,
                               ClustInner inner, double tol, std::size_t maxiter) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    SubnetSolution<T> out;
    out.X.assign(R, zero);
    out.Q = Matrix<T>(M, R, zero);
    if (M == 0 || R == 0) return out;
    const T Mt = num_traits<T>::from_int(static_cast<long>(M));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) out.Q(i, r) = N[r] / Mt;

    if (inner == ClustInner::ProportionalEstimation) {
        for (std::size_t it = 0; it < maxiter; ++it) {
            const Matrix<T> Qold = out.Q;
            for (std::size_t r = 0; r < R; ++r) {
                if (!(N[r] > zero)) continue;
                std::vector<T> W(M, zero);
                T wsum = zero;
                for (std::size_t i = 0; i < M; ++i) {
                    // PE arrival-instant local queue, inflated by the foreign share
                    T loc = zero;
                    for (std::size_t s = 0; s < R; ++s) loc += out.Q(i, s);
                    loc -= out.Q(i, r) / N[r];
                    const T A = T(T(loc + Uk[i]) / T(one - Uk[i]));
                    W[i] = L(i, r) * T(one + A);
                    wsum += W[i];
                }
                const T den = T(Z[r] + wsum);
                if (den == zero) throw NumericError("pfqn_clust: zero subnetwork cycle time");
                out.X[r] = N[r] / den;
                for (std::size_t i = 0; i < M; ++i) out.Q(i, r) = out.X[r] * W[i];
            }
            double d = 0.0;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    const double v =
                        std::fabs(num_traits<T>::to_double(T(out.Q(i, r) - Qold(i, r))));
                    if (v > d) d = v;
                }
            if (d < tol) break;
        }
        return out;
    }

    std::vector<Matrix<T> > Qs(R + 1, out.Q);
    std::vector<std::vector<std::vector<T> > > Delta(
        M, std::vector<std::vector<T> >(R, std::vector<T>(R, zero)));
    for (int pass = 0; pass < 3; ++pass) {
        for (std::size_t s = 0; s <= R; ++s) {
            std::vector<T> Ns(N);
            if (s > 0) Ns[s - 1] = T(Ns[s - 1] - one);
            Qs[s] = clust_core(L, Ns, Z, Uk, Qs[s], Delta, tol, maxiter);
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t s = 1; s <= R; ++s) {
                    const T ns = (r == s - 1) ? T(N[r] - one) : N[r];
                    if (N[r] > zero && ns > zero)
                        Delta[i][r][s - 1] = Qs[s](i, r) / ns - Qs[0](i, r) / N[r];
                    else if (N[r] > zero)
                        Delta[i][r][s - 1] = zero - Qs[0](i, r) / N[r];
                    else
                        Delta[i][r][s - 1] = zero;
                }
    }
    Qs[0] = clust_core(L, N, Z, Uk, Qs[0], Delta, tol, maxiter);
    return clust_forward(L, N, Z, Uk, Qs[0], Delta);
}

}  // namespace detail

/**
 * @param L (M x R) demands, @param N (R) populations, @param Z (R) think times
 * @param subnets      per subnetwork, the 0-based station indices it contains;
 *                     empty for the automatic decomposition described above
 * @param localclasses per subnetwork, the 0-based classes local to it; empty
 *                     for the automatic decomposition
 * @param inner algorithm run inside a subnetwork
 * @param tol convergence tolerance
 * @param maxiter outer-iteration cap
 */
template <class T>
AmvaResult<T> pfqn_clust(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                         const std::vector<std::vector<std::size_t> >& subnets,
                         const std::vector<std::vector<std::size_t> >& localclasses,
                         ClustInner inner = ClustInner::Linearizer, double tol = 1e-6,
                         std::size_t maxiter = 1000) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_clust: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_clust: Z has the wrong length");
    if (subnets.size() != localclasses.size())
        throw InputError("pfqn_clust: subnets and localclasses disagree on the cluster count");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::vector<T> Zv = Z.empty() ? std::vector<T>(R, zero) : Z;

    const AmvaResult<T> seed = pfqn_pam(L, N, Zv, PamVariant::Basic);
    AmvaResult<T> r;
    r.XN = seed.XN;
    r.QN = seed.QN;
    r.UN = Matrix<T>(M, R, zero);
    r.RN = Matrix<T>(M, R, zero);
    if (M == 0) return r;

    std::vector<std::vector<std::size_t> > nets = subnets, locals = localclasses;
    if (nets.empty() || locals.empty()) {
        std::vector<std::size_t> bottleneck(R, 0);
        for (std::size_t s = 0; s < R; ++s) {
            bool any = false;
            T best = zero;
            for (std::size_t i = 0; i < M; ++i) {
                if (!(L(i, s) > zero)) continue;
                const T u = L(i, s) * r.XN[s];
                if (!any || u > best) {
                    best = u;
                    bottleneck[s] = i;
                    any = true;
                }
            }
        }
        std::vector<std::size_t> centres;
        for (std::size_t s = 0; s < R; ++s)
            if (std::find(centres.begin(), centres.end(), bottleneck[s]) == centres.end())
                centres.push_back(bottleneck[s]);
        std::sort(centres.begin(), centres.end());
        nets.clear();
        locals.clear();
        std::vector<bool> covered(M, false);
        for (std::size_t g = 0; g < centres.size(); ++g) {
            std::vector<std::size_t> cls;
            for (std::size_t s = 0; s < R; ++s)
                if (bottleneck[s] == centres[g]) cls.push_back(s);
            std::vector<std::size_t> st;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t a = 0; a < cls.size(); ++a)
                    if (L(i, cls[a]) > zero) {
                        st.push_back(i);
                        break;
                    }
            if (st.empty()) st.push_back(centres[g]);
            for (std::size_t a = 0; a < st.size(); ++a) covered[st[a]] = true;
            nets.push_back(st);
            locals.push_back(cls);
        }
        std::vector<std::size_t> missing;
        for (std::size_t i = 0; i < M; ++i)
            if (!covered[i]) missing.push_back(i);
        if (!missing.empty()) {
            nets.push_back(missing);
            locals.push_back(std::vector<std::size_t>());
        }
    }
    const std::size_t G = nets.size();
    std::vector<int> owner(R, -1);
    for (std::size_t g = 0; g < G; ++g)
        for (std::size_t a = 0; a < locals[g].size(); ++a)
            owner[locals[g][a]] = static_cast<int>(g);

    for (std::size_t it = 1; it <= maxiter; ++it) {
        r.iterations = it;
        const Matrix<T> Qprev = r.QN;
        std::vector<T> Qk(M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t s = 0; s < R; ++s) Qk[i] += r.QN(i, s);

        for (std::size_t g = 0; g < G; ++g) {
            const std::vector<std::size_t>& S = nets[g];
            const std::vector<std::size_t>& LC = locals[g];
            if (LC.empty() || S.empty()) continue;
            std::vector<bool> inS(M, false);
            for (std::size_t a = 0; a < S.size(); ++a) inS[S[a]] = true;
            std::vector<bool> isLocal(R, false);
            for (std::size_t a = 0; a < LC.size(); ++a) isLocal[LC[a]] = true;
            std::vector<bool> isForeign(R, false);
            for (std::size_t s = 0; s < R; ++s) {
                if (isLocal[s]) continue;
                for (std::size_t a = 0; a < S.size(); ++a)
                    if (L(S[a], s) > zero) {
                        isForeign[s] = true;
                        break;
                    }
            }
            // per-class delay in the complement of S
            std::vector<T> Zeff(LC.size(), zero);
            for (std::size_t a = 0; a < LC.size(); ++a) {
                const std::size_t c = LC[a];
                T p = zero;
                if (N[c] > zero)
                    for (std::size_t i = 0; i < M; ++i)
                        if (!inS[i])
                            p += L(i, c) * T(one + Qk[i]) / T(one + L(i, c) * r.XN[c] / N[c]);
                Zeff[a] = T(Zv[c] + p);
            }
            // utilization left in S by the foreign classes
            std::vector<T> Uk(S.size(), zero);
            const T cap = T(one - num_traits<T>::from_int(1) / num_traits<T>::from_int(100000000));
            for (std::size_t b = 0; b < S.size(); ++b) {
                const std::size_t i = S[b];
                T u = zero;
                for (std::size_t s = 0; s < R; ++s)
                    if (isForeign[s] && N[s] > zero)
                        u += L(i, s) * r.XN[s] / T(one + L(i, s) * r.XN[s] / N[s]);
                Uk[b] = (u < cap) ? u : cap;
            }
            Matrix<T> Lsub(S.size(), LC.size(), zero);
            std::vector<T> Nsub(LC.size(), zero);
            for (std::size_t b = 0; b < S.size(); ++b)
                for (std::size_t a = 0; a < LC.size(); ++a) Lsub(b, a) = L(S[b], LC[a]);
            for (std::size_t a = 0; a < LC.size(); ++a) Nsub[a] = N[LC[a]];

            const detail::SubnetSolution<T> sub =
                detail::clust_subnet(Lsub, Nsub, Zeff, Uk, inner, tol, maxiter);
            for (std::size_t a = 0; a < LC.size(); ++a) {
                r.XN[LC[a]] = sub.X[a];
                for (std::size_t b = 0; b < S.size(); ++b) r.QN(S[b], LC[a]) = sub.Q(b, a);
            }
            // the local classes still hold jobs outside S
            for (std::size_t a = 0; a < LC.size(); ++a) {
                const std::size_t c = LC[a];
                if (!(N[c] > zero)) continue;
                for (std::size_t i = 0; i < M; ++i)
                    if (!inS[i])
                        r.QN(i, c) = r.XN[c] * L(i, c) * T(one + Qk[i]) /
                                     T(one + L(i, c) * r.XN[c] / N[c]);
            }
        }
        // a class owned by no subnetwork keeps the seed throughput
        for (std::size_t s = 0; s < R; ++s)
            if (owner[s] < 0)
                for (std::size_t i = 0; i < M; ++i) r.QN(i, s) = r.XN[s] * L(i, s);

        double d = 0.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t s = 0; s < R; ++s) {
                const double v = std::fabs(num_traits<T>::to_double(T(r.QN(i, s) - Qprev(i, s))));
                if (v > d) d = v;
            }
        if (d < tol) {
            r.converged = true;
            break;
        }
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) {
            r.UN(i, s) = r.XN[s] * L(i, s);
            r.RN(i, s) = (N[s] == zero || r.XN[s] == zero) ? zero : T(r.QN(i, s) / r.XN[s]);
        }
    return r;
}

template <class T>
AmvaResult<T> pfqn_clust(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_clust(L, N, Z, std::vector<std::vector<std::size_t> >(),
                      std::vector<std::vector<std::size_t> >());
}

template <class T>
AmvaResult<T> pfqn_clust(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_clust(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CLUST_H
