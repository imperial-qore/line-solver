#pragma once
/**
 * @file mapqn_amva.h
 * @brief Horizontal-cut mean value analysis for a MAP server (SolverMVA method 'amva.mapqn').
 *
 * Closed multiclass network of an exponential infinite-server station (think rate mu_r for
 * class r) and one FCFS single-server station whose class-r service is the MAP (D0_r, D1_r);
 * the MAP of class r moves only while a class-r job is in service and is frozen otherwise,
 * the convention of the CTMC solver. The recursion walks the population lattice n <= N in
 * lexicographic order and solves ONE linear R x R system per point. Its unknowns are the
 * per-phase means Q_r^k = E[n_r 1{k}] over the joint phase k = (k_1..k_R), the busy laws
 * U_r^k = P[serving r, k], the phase law pi_k and the throughputs X_r.
 *
 * Exact relations: the joint phase balance, the class marginals U_r = X_r E[S_r] theta_r
 * and the per-class horizontal cut (generator balance of n_r 1{k}) of Casale-Smirni,
 * "MAP-AMVA: Approximate Mean Value Analysis of Bursty Systems", IEEE/IFIP DSN 2009.
 * Closures: the product busy law theta_r(k_r) prod_{s != r} phi_s(k_s), phi the
 * post-completion law of a frozen MAP, which solves the phase balance identically; the
 * service-age closure of the cross term E[n_r 1{serving s} 1{k}] (class r accumulates at
 * its throughput over the elapsed class-s service, whose mean given the phase is
 * theta_s (-D0_s)^{-1} / theta_s); Little's law resolved by arrival phase with the exact
 * FCFS response of the queue composition seen at n - e_r (the multiclass arrival
 * theorem). K_r = 1 for every class reproduces multiclass FCFS MVA on class means.
 *
 * Port of matlab/src/api/mapqn/mapqn_amva.m; the arithmetic is done in double whatever T,
 * since the recursion interpolates response tables at fractional populations.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mapqn {

/** X: class throughputs; Qq: mean queue lengths at the MAP station (job in service
 *  included); U: busy probability per class, X E[S]; ES: mean service times; pi: joint
 *  phase law at N (class R fastest). */
template <class T>
struct MapqnAmvaResult {
    std::vector<T> X, Qq, U, ES, pi;
};

namespace amva_detail {

inline std::vector<double> solve_linear(std::vector<std::vector<double>> A, std::vector<double> b) {
    const std::size_t n = b.size();
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t i = c + 1; i < n; ++i)
            if (std::fabs(A[i][c]) > std::fabs(A[p][c])) p = i;
        std::swap(A[c], A[p]);
        std::swap(b[c], b[p]);
        const double piv = A[c][c];
        if (piv == 0.0) throw InputError("mapqn_amva: singular linear system");
        for (std::size_t i = c + 1; i < n; ++i) {
            const double f = A[i][c] / piv;
            if (f == 0.0) continue;
            for (std::size_t j = c; j < n; ++j) A[i][j] -= f * A[c][j];
            b[i] -= f * b[c];
        }
    }
    std::vector<double> x(n, 0.0);
    for (std::size_t ii = n; ii-- > 0;) {
        double s = b[ii];
        for (std::size_t j = ii + 1; j < n; ++j) s -= A[ii][j] * x[j];
        x[ii] = s / A[ii][ii];
    }
    return x;
}

inline std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& A) {
    const std::size_t n = A.size();
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
    for (std::size_t c = 0; c < n; ++c) {
        std::vector<double> e(n, 0.0);
        e[c] = 1.0;
        const std::vector<double> col = solve_linear(A, e);
        for (std::size_t i = 0; i < n; ++i) inv[i][c] = col[i];
    }
    return inv;
}

/** theta G = 0, theta 1 = 1: transpose G and replace the last equation by the normalization. */
inline std::vector<double> stationary(const std::vector<std::vector<double>>& G) {
    const std::size_t K = G.size();
    std::vector<std::vector<double>> A(K, std::vector<double>(K, 0.0));
    std::vector<double> rhs(K, 0.0);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) A[i][j] = G[j][i];
    for (std::size_t j = 0; j < K; ++j) A[K - 1][j] = 1.0;
    rhs[K - 1] = 1.0;
    std::vector<double> th = solve_linear(A, rhs);
    double s = 0.0;
    for (double& v : th) { v = std::max(v, 0.0); s += v; }
    for (double& v : th) v /= s;
    return th;
}

}  // namespace amva_detail

template <class T>
MapqnAmvaResult<T> mapqn_amva(const std::vector<T>& mu_in, const std::vector<Matrix<T>>& D0s,
                              const std::vector<Matrix<T>>& D1s, const std::vector<int>& N) {
    using amva_detail::inverse;
    using amva_detail::solve_linear;
    using amva_detail::stationary;
    using Mat = std::vector<std::vector<double>>;
    const std::size_t R = N.size();
    if (mu_in.size() != R || D0s.size() != R || D1s.size() != R)
        throw InputError("mapqn_amva: mu, D0s, D1s and N must all have one entry per class");
    std::vector<double> mu(R);
    std::vector<std::size_t> Ks(R);
    std::vector<Mat> D0(R), D1(R);
    for (std::size_t r = 0; r < R; ++r) {
        mu[r] = num_traits<T>::to_double(mu_in[r]);
        Ks[r] = D0s[r].rows();
        D0[r].assign(Ks[r], std::vector<double>(Ks[r], 0.0));
        D1[r].assign(Ks[r], std::vector<double>(Ks[r], 0.0));
        for (std::size_t i = 0; i < Ks[r]; ++i)
            for (std::size_t j = 0; j < Ks[r]; ++j) {
                D0[r][i][j] = num_traits<T>::to_double(D0s[r](i, j));
                D1[r][i][j] = num_traits<T>::to_double(D1s[r](i, j));
            }
    }
    std::size_t K = 1;
    for (std::size_t r = 0; r < R; ++r) K *= Ks[r];
    std::vector<std::size_t> stride(R, 1);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t t = r + 1; t < R; ++t) stride[r] *= Ks[t];      // class R-1 is the fastest index
    std::vector<std::vector<std::size_t>> krOf(K, std::vector<std::size_t>(R, 0));
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t r = 0; r < R; ++r) krOf[k][r] = (k / stride[r]) % Ks[r];
    int Ntot = 0;
    for (std::size_t r = 0; r < R; ++r) Ntot += N[r];

    std::vector<Mat> G(R), Ainv(R), Tt(R);
    std::vector<std::vector<double>> th(R), phi(R), age(R);
    std::vector<double> ES(R, 0.0), abar(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        const std::size_t Kr = Ks[r];
        G[r].assign(Kr, std::vector<double>(Kr, 0.0));
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) G[r][i][j] = D0[r][i][j] + D1[r][i][j];
        th[r] = stationary(G[r]);
        double rate = 0.0;
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) rate += th[r][i] * D1[r][i][j];
        ES[r] = 1.0 / rate;
        phi[r].assign(Kr, 0.0);                                             // post-completion phase law
        for (std::size_t j = 0; j < Kr; ++j) {
            double a = 0.0;
            for (std::size_t i = 0; i < Kr; ++i) a += th[r][i] * D1[r][i][j];
            phi[r][j] = a * ES[r];
        }
        Mat negD0(Kr, std::vector<double>(Kr, 0.0));
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) negD0[i][j] = -D0[r][i][j];
        const Mat negD0inv = inverse(negD0);
        std::vector<double> s(Kr, 0.0);                                     // mean time to the next completion
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) s[i] += negD0inv[i][j];
        Mat P(Kr, std::vector<double>(Kr, 0.0));                            // embedded phase transition
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) {
                double a = 0.0;
                for (std::size_t m = 0; m < Kr; ++m) a += negD0inv[i][m] * D1[r][m][j];
                P[i][j] = a;
            }
        Tt[r].assign(Ntot + 2, std::vector<double>(Kr, 0.0));               // T[j][k] = e_k'(I+P+..+P^(j-1)) s
        std::vector<double> acc(Kr, 0.0), v = s;
        for (int j = 1; j <= Ntot + 1; ++j) {
            for (std::size_t i = 0; i < Kr; ++i) acc[i] += v[i];
            Tt[r][j] = acc;
            std::vector<double> nv(Kr, 0.0);
            for (std::size_t i = 0; i < Kr; ++i)
                for (std::size_t m = 0; m < Kr; ++m) nv[i] += P[i][m] * v[m];
            v = nv;
        }
        std::vector<double> w(Kr, 0.0);                                     // theta (-D0)^{-1}
        for (std::size_t j = 0; j < Kr; ++j)
            for (std::size_t i = 0; i < Kr; ++i) w[j] += th[r][i] * negD0inv[i][j];
        age[r].assign(Kr, 0.0);
        for (std::size_t j = 0; j < Kr; ++j) { age[r][j] = w[j] / th[r][j]; abar[r] += w[j]; }
        Mat A(Kr, std::vector<double>(Kr, 0.0));
        for (std::size_t i = 0; i < Kr; ++i)
            for (std::size_t j = 0; j < Kr; ++j) A[i][j] = G[r][i][j] - (i == j ? mu[r] : 0.0);
        Ainv[r] = inverse(A);
    }
    // joint-phase shapes: idle law F and class-r busy law u[r]
    std::vector<double> F(K, 1.0);
    std::vector<std::vector<double>> u(R, std::vector<double>(K, 1.0));
    for (std::size_t k = 0; k < K; ++k) {
        for (std::size_t r = 0; r < R; ++r) F[k] *= phi[r][krOf[k][r]];
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t s = 0; s < R; ++s) u[r][k] *= (s == r) ? th[s][krOf[k][s]] : phi[s][krOf[k][s]];
    }
    auto apply_axis = [&](const std::vector<double>& V, const Mat& M, std::size_t r) {
        std::vector<double> out(K, 0.0);
        for (std::size_t k = 0; k < K; ++k) {
            const std::size_t kr = krOf[k][r];
            const std::size_t base = k - kr * stride[r];
            double a = 0.0;
            for (std::size_t h = 0; h < Ks[r]; ++h) a += V[base + h * stride[r]] * M[h][kr];
            out[k] = a;
        }
        return out;
    };
    auto t_at = [&](std::size_t r, double b, std::size_t kr) {
        const Mat& Tr = Tt[r];
        long j0 = static_cast<long>(std::floor(b));
        j0 = std::max<long>(0, std::min<long>(j0, static_cast<long>(Tr.size()) - 2));
        const double f = std::min(std::max(b - static_cast<double>(j0), 0.0), 1.0);
        return (1.0 - f) * Tr[static_cast<std::size_t>(j0)][kr] + f * Tr[static_cast<std::size_t>(j0) + 1][kr];
    };
    // population lattice in lexicographic order: n - e_r always precedes n
    std::vector<std::size_t> lstride(R, 1);
    std::size_t L = 1;
    for (std::size_t rr = R; rr-- > 0;) { lstride[rr] = L; L *= static_cast<std::size_t>(N[rr] + 1); }
    std::vector<std::vector<std::vector<double>>> Qs(L, std::vector<std::vector<double>>(R, std::vector<double>(K, 0.0)));
    std::vector<std::vector<double>> pis(L, std::vector<double>(K, 0.0)), Xs(L, std::vector<double>(R, 0.0));
    pis[0] = F;
    for (std::size_t l = 1; l < L; ++l) {
        std::vector<int> n(R, 0);
        for (std::size_t r = 0; r < R; ++r) n[r] = static_cast<int>((l / lstride[r]) % static_cast<std::size_t>(N[r] + 1));
        std::vector<std::vector<std::vector<double>>> b(R), bN(R);
        std::vector<std::vector<double>> Rk(R, std::vector<double>(K, 0.0));
        for (std::size_t r = 0; r < R; ++r) {
            if (n[r] < 1) continue;
            const std::size_t lp = l - lstride[r];
            b[r].assign(R, std::vector<double>(K, 0.0));
            for (std::size_t t = 0; t < R; ++t)
                for (std::size_t k = 0; k < K; ++k) b[r][t][k] = pis[lp][k] > 0 ? Qs[lp][t][k] / pis[lp][k] : 0.0;
            for (std::size_t k = 0; k < K; ++k) {
                double a = 0.0;
                for (std::size_t t = 0; t < R; ++t) a += t_at(t, b[r][t][k] + (t == r ? 1.0 : 0.0), krOf[k][t]);
                Rk[r][k] = a;
            }
            bN[r].assign(R, std::vector<double>(K, 0.0));
            for (std::size_t t = 0; t < R; ++t) {
                if (t == r || n[t] == 0) continue;
                const double Xt = Xs[lp][t];
                double Qt = 0.0;
                for (std::size_t k = 0; k < K; ++k) Qt += Qs[lp][t][k];
                if (Xt > 0) {
                    const double W = std::max(Qt / Xt - abar[r], 0.0);
                    for (std::size_t k = 0; k < K; ++k)
                        bN[r][t][k] = Xt * std::min(W + age[r][krOf[k][r]], static_cast<double>(n[t]) / Xt);
                }
            }
        }
        // the cut, linear in X: Q_r = c0[r] + sum_s X_s c1[r][s]
        std::vector<std::vector<double>> c0(R);
        std::vector<std::vector<std::vector<double>>> c1(R, std::vector<std::vector<double>>(R));
        for (std::size_t r = 0; r < R; ++r) {
            if (n[r] < 1) continue;
            std::vector<double> v0(K, 0.0);
            for (std::size_t k = 0; k < K; ++k) v0[k] = -mu[r] * n[r] * F[k];
            c0[r] = apply_axis(v0, Ainv[r], r);
            for (std::size_t s = 0; s < R; ++s) {
                std::vector<double> term(K, 0.0);
                for (std::size_t k = 0; k < K; ++k) term[k] = -mu[r] * n[r] * ES[s] * (u[s][k] - F[k]);
                if (s == r) {
                    const std::vector<double> ud = apply_axis(u[r], D1[r], r);
                    for (std::size_t k = 0; k < K; ++k) term[k] += ES[r] * ud[k];
                } else if (n[s] >= 1) {
                    std::vector<double> W(K, 0.0);
                    for (std::size_t k = 0; k < K; ++k) W[k] = ES[s] * u[s][k] * bN[s][r][k];
                    const std::vector<double> a1 = apply_axis(W, G[r], r), a2 = apply_axis(W, G[s], s);
                    for (std::size_t k = 0; k < K; ++k) term[k] += a1[k] - a2[k];
                }
                c1[r][s] = apply_axis(term, Ainv[r], r);
            }
        }
        // Little's law by arrival phase: one R x R solve
        Mat M(R, std::vector<double>(R, 0.0));
        std::vector<double> v(R, 0.0);
        for (std::size_t r = 0; r < R; ++r) M[r][r] = 1.0;
        for (std::size_t r = 0; r < R; ++r) {
            if (n[r] < 1) continue;
            double a = 0.0;
            for (std::size_t k = 0; k < K; ++k) a += (n[r] * F[k] - c0[r][k]) * Rk[r][k];
            v[r] = n[r] - mu[r] * a;
            M[r][r] = 1.0 / mu[r];
            for (std::size_t s = 0; s < R; ++s) {
                if (n[s] < 1) continue;
                double e = 0.0;
                for (std::size_t k = 0; k < K; ++k) e += (n[r] * ES[s] * (u[s][k] - F[k]) - c1[r][s][k]) * Rk[r][k];
                M[r][s] += mu[r] * e;
            }
        }
        const std::vector<double> X = solve_linear(M, v);
        std::vector<double> pi(K, 0.0);
        double psum = 0.0;
        for (std::size_t k = 0; k < K; ++k) {
            double p = F[k];
            for (std::size_t s = 0; s < R; ++s) p += X[s] * ES[s] * (u[s][k] - F[k]);
            pi[k] = std::max(p, 0.0);
            psum += pi[k];
        }
        for (double& p : pi) p /= psum;
        for (std::size_t r = 0; r < R; ++r) {
            if (n[r] < 1) continue;
            std::vector<double> Ur(K, 0.0), Qr(K, 0.0), Wr(K, 0.0);
            double usum = 0.0, wsum = 0.0;
            for (std::size_t k = 0; k < K; ++k) {
                Ur[k] = X[r] * ES[r] * u[r][k];
                Qr[k] = c0[r][k];
                for (std::size_t s = 0; s < R; ++s) Qr[k] += X[s] * c1[r][s][k];
                Wr[k] = std::max(Qr[k] - Ur[k], 0.0);
                usum += Ur[k];
                wsum += Wr[k];
            }
            // project onto Q >= U keeping the flow-balance total n_r - X_r/mu_r
            const double tot = std::max(n[r] - X[r] / mu[r] - usum, 0.0);
            for (std::size_t k = 0; k < K; ++k) Qs[l][r][k] = Ur[k] + (wsum > 0 ? Wr[k] * tot / wsum : 0.0);
        }
        pis[l] = pi;
        Xs[l] = X;
    }
    MapqnAmvaResult<T> out;
    out.X.resize(R); out.Qq.resize(R); out.U.resize(R); out.ES.resize(R); out.pi.resize(K);
    for (std::size_t r = 0; r < R; ++r) {
        double q = 0.0;
        for (std::size_t k = 0; k < K; ++k) q += Qs[L - 1][r][k];
        out.X[r] = num_traits<T>::from_double(Xs[L - 1][r]);
        out.Qq[r] = num_traits<T>::from_double(q);
        out.U[r] = num_traits<T>::from_double(Xs[L - 1][r] * ES[r]);
        out.ES[r] = num_traits<T>::from_double(ES[r]);
    }
    for (std::size_t k = 0; k < K; ++k) out.pi[k] = num_traits<T>::from_double(pis[L - 1][k]);
    return out;
}

}  // namespace mapqn
}  // namespace line
