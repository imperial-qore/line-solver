/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_COMOMRM_H
#define LINE_API_PFQN_COMOMRM_H

/**
 * CoMoM (class-oriented method of moments) for the finite repairman model:
 * one queueing station of multiplicity m plus a delay.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_comomrm.m. This is the routine
 * pfqn_nc dispatches to for the repairman case under both the 'comom' and the
 * 'default' method.
 *
 * The basis h at stage r holds 2r normalizing constants: the r "plus" moments
 * G(n + e_s) and the r constants G(n) and G(n - e_s), s < r. Adding one
 * class-r job applies
 *
 *   h <- ( F1r + F2r / n_r ) h,
 *
 * with F1r a single 1 in the leading position and F2r built from the
 * convolution equation (CE) and the population constraints (PC),
 *
 *   F2r = [ -C^{-1} A12 B2r ; B2r ],   B2r = [ m L_r I , Z_r I ],
 *
 * where C^{-1} is available in closed form for the repairman model (the
 * reference writes it out rather than solving a system, and so does this port,
 * so no linear solver is involved anywhere). Moving from class r-1 to class r
 * expands the basis by interleaving two fresh entries carried over from the
 * PREVIOUS step's basis.
 *
 * SCALING. The reference renormalizes h after every step and folds the
 * discarded factors into a log accumulator; the expansion step then divides
 * the carried-over entries by the last scale factor to put them on the current
 * scale. This port carries the UNSCALED basis and the previous unscaled basis
 * instead, which makes the expansion a plain copy and makes
 *
 *   G(N) = Gremaind * h[R]
 *
 * exact, with h[R] the reference's h(end-(R-1)). The two formulations are
 * algebraically identical: MATLAB's lG = lG0 + log(h(end-(R-1))) + sum(log
 * scale) is the log of exactly that product, since h_scaled * prod(scale) is
 * the unscaled basis by construction.
 *
 * Arithmetic: EXACT-CAPABLE. Every operation is an addition, a multiplication
 * or a division in the field of the inputs. The only transcendentals in the
 * reference are the log-space scale bookkeeping described above and the
 * factln/exp seeding of the zero-think-time head, which is a ratio of
 * factorials and is formed here directly.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_nc_sanitize.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct ComomResult {
    T G;                    ///< normalizing constant
    double lG;              ///< its logarithm
    std::vector<T> basis;   ///< the final unscaled basis h
};

/**
 * @param L (1 x R) demands at the single queueing station
 * @param N (R) populations
 * @param Z (K x R) think times
 * @param m multiplicity of the queueing station
 */
template <class T>
ComomResult<T> pfqn_comomrm(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                            int m) {
    if (!L.empty() && L.rows() != 1)
        throw InputError("pfqn_comomrm: the solver accepts at most a single queueing station");
    if (m < 1) throw InputError("pfqn_comomrm: the multiplicity must be at least one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T mT = num_traits<T>::from_int(m);

    const NcSanitizeResult<T> san = pfqn_nc_sanitize(L, N, Z);
    const std::size_t R = san.N.size();

    ComomResult<T> res;
    if (R == 0) {
        res.G = san.Gremaind;
        res.lG = san.lGremaind;
        res.basis.assign(1, one);
        return res;
    }

    std::vector<T> Lv(R, zero), Zv(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        if (!san.L.empty()) Lv[r] = san.L(0, r);
        for (std::size_t k = 0; k < san.Z.rows(); ++k) Zv[r] += san.Z(k, r);
    }

    // pfqn_nc_sanitize already orders the zero-think-time classes first, so the
    // head is a prefix and the iteration below starts right after it.
    std::size_t nzt = 0;
    while (nzt < R && Zv[nzt] == zero) ++nzt;

    // basis-seeding rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    std::vector<int> nvec(R, 0);
    for (std::size_t z = 0; z < nzt; ++z) nvec[z] = san.N[z];
    const auto headTerm = [&](const std::vector<int>& v, int extra) {
        int tot = 0;
        for (int x : v) tot += x;
        T num = num_factorial<T>(static_cast<unsigned>(tot + extra));
        for (int x : v) num /= num_factorial<T>(static_cast<unsigned>(x));
        return num;
    };

    std::vector<T> h;
    if (nzt > 0) {
        h.assign(2 + 2 * nzt, zero);
        std::size_t k = 0;
        h[k++] = headTerm(nvec, m);  // factln(sum+m+1-1)
        for (std::size_t z = 0; z < nzt; ++z) {
            std::vector<int> v = nvec;
            v[z] -= 1;
            h[k++] = headTerm(v, m);
        }
        h[k++] = headTerm(nvec, m - 1);
        for (std::size_t z = 0; z < nzt; ++z) {
            std::vector<int> v = nvec;
            v[z] -= 1;
            h[k++] = headTerm(v, m - 1);
        }
    } else {
        h.assign(2, one);
    }

    if (nzt == R) {
        // Every class has zero think time: the trivial model is the answer.
        res.basis = h;
        res.G = san.Gremaind * h[h.size() - R - 1];
        res.lG = num_traits<T>::log_as_double(res.G);
        return res;
    }

    // ---- iterate over the remaining classes ---------------------------------
    std::vector<T> hprev = h;
    Matrix<T> F1, F2;  // rebuilt once per class, reused across its jobs
    for (std::size_t r = nzt; r < R; ++r) {
        const std::size_t rr = r + 1;  // the reference's 1-based class index
        for (int Nr = 1; Nr <= san.N[r]; ++Nr) {
            nvec[r] += 1;
            if (Nr == 1) {
                if (rr > nzt + 1) {
                    // basis expansion rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
                    const std::size_t p = rr - 1;
                    std::vector<T> hr(2 * rr, zero);
                    for (std::size_t i = 0; i < p; ++i) hr[i] = h[i];
                    for (std::size_t i = 0; i < p; ++i) hr[rr + i] = h[p + i];
                    hr[p] = hprev[0];
                    hr[2 * rr - 1] = hprev[p];
                    h.swap(hr);
                }
                // A12 (rr x rr)
                Matrix<T> A12(rr, rr, zero);
                A12(0, 0) = -one;
                for (std::size_t s = 0; s + 1 < rr; ++s) {
                    A12(1 + s, 0) = num_traits<T>::from_int(san.N[s]);
                    A12(1 + s, 1 + s) = -Zv[s];
                }
                // B2r (rr x 2rr) = [ m L_r I , Z_r I ]
                Matrix<T> B2r(rr, 2 * rr, zero);
                for (std::size_t i = 0; i < rr; ++i) {
                    B2r(i, i) = mT * Lv[r];
                    B2r(i, rr + i) = Zv[r];
                }
                // iC (rr x rr): closed form of C^{-1} for the repairman model.
                Matrix<T> iC(rr, rr, zero);
                const T minv = -one / mT;
                for (std::size_t j = 0; j < rr; ++j) iC(0, j) = minv;
                for (std::size_t i = 1; i < rr; ++i) iC(i, i) = minv;
                iC(0, 0) = one;

                // F2r = [ -iC * A12 * B2r ; B2r ]
                Matrix<T> W(rr, rr, zero);  // -iC * A12
                for (std::size_t i = 0; i < rr; ++i)
                    for (std::size_t j = 0; j < rr; ++j) {
                        T s = zero;
                        for (std::size_t k = 0; k < rr; ++k) s += iC(i, k) * A12(k, j);
                        W(i, j) = -s;
                    }
                F2 = Matrix<T>(2 * rr, 2 * rr, zero);
                for (std::size_t i = 0; i < rr; ++i)
                    for (std::size_t j = 0; j < 2 * rr; ++j) {
                        T s = zero;
                        for (std::size_t k = 0; k < rr; ++k) s += W(i, k) * B2r(k, j);
                        F2(i, j) = s;
                    }
                for (std::size_t i = 0; i < rr; ++i)
                    for (std::size_t j = 0; j < 2 * rr; ++j) F2(rr + i, j) = B2r(i, j);
                F1 = Matrix<T>(2 * rr, 2 * rr, zero);
                F1(0, 0) = one;
            }
            hprev = h;
            const T inv = one / num_traits<T>::from_int(nvec[r]);
            std::vector<T> hn(2 * rr, zero);
            for (std::size_t i = 0; i < 2 * rr; ++i) {
                T s = zero;
                for (std::size_t j = 0; j < 2 * rr; ++j) s += (F1(i, j) + F2(i, j) * inv) * hprev[j];
                hn[i] = s;
            }
            h.swap(hn);
        }
    }

    res.basis = h;
    res.G = san.Gremaind * h[h.size() - R];
    res.lG = num_traits<T>::log_as_double(res.G);
    return res;
}

/** Overload with the unit multiplicity default. */
template <class T>
ComomResult<T> pfqn_comomrm(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_comomrm(L, N, Z, 1);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_COMOMRM_H
