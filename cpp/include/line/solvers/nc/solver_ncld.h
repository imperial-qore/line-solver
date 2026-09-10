/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NCLD_H
#define LINE_SOLVERS_NC_SOLVER_NCLD_H

/**
 * Port of `solver_ncld.m`: the LOAD-DEPENDENT normalizing-constant analyzer.
 *
 * WHAT MAKES IT A DIFFERENT SOLVER. `solver_nc.h` assumes every station serves
 * at a constant rate, so its constants come from `pfqn_nc`. Here each station
 * carries a rate lattice mu_i(n) and the constants come from `pfqn_ncld`. That
 * is not a refinement of the same formula: the arrival-theorem shortcut used
 * there does not hold, and the queue length is recovered instead from the
 * CONDITIONAL MVA identity
 *
 *   Q_ic = exp(log L_ic + lG(mu^i, N-1_c) - log mu_i(1) - lG(N-1_c)) X_c (1 + CQ_i)
 *
 * where mu^i is the lattice shifted at station i (`pfqn_mushift`) and CQ_i is
 * assembled from the flow-equivalent complement (`pfqn_fnc`). Four constants
 * per (station, chain) rather than one, which is the price of load dependence.
 *
 * WHY MULTISERVER ARRIVES HERE. `@@SolverNC/runAnalyzer` rewrites a genuine
 * c-server station as mu(n) = min(n, c) whenever the model is product-form, and
 * that lattice is EXACT where Seidmann's approximation in `solver_nc.h` is not.
 * The server count is deliberately kept alongside the lattice: utilization is
 * the fraction of the c servers busy, and c cannot be read back from
 * min(1:Nt, c) once the population falls below it.
 *
 * MIXED MODELS take a different route entirely -- the Bruell-Balbo-Ashfari
 * effective-capacity MVA (`pfqn_mvaldmx`) -- because an open chain has no
 * finite population for the constant to be evaluated at.
 *
 * ARITHMETIC. As in `solver_nc.h`, every measure is a difference of logarithms
 * of normalizing constants; a non-transcendental backend is refused by name.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_nonexp_approx.h"
#include "line/api/pfqn/pfqn_fnc.h"
#include "line/api/pfqn/pfqn_mushift.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/api/pfqn/pfqn_ncldmx.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/solver_nc_conv.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

namespace detail {

/**
 * Map `options.method` onto the load-dependent dispatcher.
 *
 * Every name the reference's `compute_norm_const_ld` switches on is now
 * dispatched; `pfqn_ncld` itself refuses the log-domain ones in an exact field.
 * A name that belongs to the load-INDEPENDENT ladder alone (say 'cub') reaches
 * this function only when the model carries a rate lattice, and there is no
 * load-dependent algorithm behind it, so it is refused here by name rather than
 * quietly solved by a different one.
 */
inline pfqn::NcldMethod ncld_pfqn_method(const std::string& method) {
    return pfqn::ncld_method_of(method);
}

/**
 * True when every finite multi-server station already carries mu(n)=min(n,c),
 * i.e. the multiserver is fully described by the load-dependent rates.
 *
 * Port of the local `lld_encodes_multiserver` of the reference.
 */
template <class T>
bool lld_encodes_multiserver(const qn::NetworkStruct<T>& sn) {
    double Ntot = 0.0;
    for (const qn::JobClass& c : sn.classes)
        if (std::isfinite(c.population)) Ntot += c.population;
    if (!std::isfinite(Ntot) || Ntot < 1.0) return false;
    const std::size_t Nt = static_cast<std::size_t>(std::llround(Ntot));
    bool anyLld = false;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (!sn.stations[i].lldscaling.empty()) anyLld = true;
    if (!anyLld) return false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const double c = sn.stations[i].nservers;
        if (!std::isfinite(c) || c <= 1.0) continue;
        const std::vector<T>& lld = sn.stations[i].lldscaling;
        if (lld.size() < Nt) return false;
        for (std::size_t n = 1; n <= Nt; ++n)
            if (num_traits<T>::to_double(lld[n - 1]) != std::min<double>(n, c)) return false;
    }
    return true;
}

/**
 * First column (1-based) of the trailing constant run of a limited
 * load-dependence row, i.e. the level b with mu(n) = mu(b) for every n >= b; 1
 * on a flat row. This is the level `pfqn_ldmx_ec` infers from the row it is
 * handed, so a row cut below it is read as a different, slower station.
 */
template <class T>
std::size_t lld_saturation_level(const Matrix<T>& lldscaling, std::size_t ist) {
    std::size_t b = lldscaling.cols();
    if (b == 0) return 1;
    while (b > 1 && lldscaling(ist, b - 2) == lldscaling(ist, b - 1)) --b;
    return b;
}

}  // namespace detail

/**
 * Port of `solver_ncld.m`.
 *
 * @param sn_in the refreshed struct; its lldscaling is completed in place on a
 *              two-station multiserver model, exactly as the reference does
 * @param opt   solver controls
 */
template <class T>
NcSolution<T> solver_ncld(const qn::NetworkStruct<T>& sn_in, const NcSolverOptions& opt) {
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn_in;
        (void)opt;
        throw UnsupportedError(
            "solver_ncld: the load-dependent normalizing-constant analyzer forms "
            "X = exp(lG(N-1_c) - lG(N)) and needs transcendental arithmetic; this backend has "
            "none");
    } else {
        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        qn::NetworkStruct<T> sn = sn_in;
        const std::size_t M = sn.nstations, K = sn.nclasses, C = sn.nchains;

        std::vector<double> nservers(M, 1.0);
        std::vector<bool> isFCFS(M, false);
        bool anyFCFS = false, anyMulti = false;
        for (std::size_t i = 0; i < M; ++i) {
            nservers[i] = sn.stations[i].nservers;
            isFCFS[i] = sn.stations[i].sched == SchedStrategy::FCFS;
            if (isFCFS[i]) anyFCFS = true;
            if (std::isfinite(nservers[i]) && nservers[i] > 1.0) anyMulti = true;
        }

        double Ntot_d = 0.0;
        bool allClosed = true;
        for (const qn::JobClass& c : sn.classes) {
            if (std::isinf(c.population))
                allClosed = false;
            else
                Ntot_d += c.population;
        }
        const std::size_t Nt = static_cast<std::size_t>(
            std::max<double>(1.0, std::ceil(std::isfinite(Ntot_d) ? Ntot_d : 1.0)));

        // Class-dependent rates beta_{i,r}(n) route to the convolution analyzer
        // BEFORE anything else and unconditionally on the method: every
        // algorithm below reads mu(n) from lldscaling alone and has no way to
        // apply beta, so gating this on 'exact' would silently return the
        // UNSCALED network for every other method.
        for (std::size_t i = 0; i < M; ++i)
            if (static_cast<bool>(sn.stations[i].cdscaling)) return solver_nc_conv(sn, opt);

        bool anyLld = false;
        for (std::size_t i = 0; i < M; ++i)
            if (!sn.stations[i].lldscaling.empty()) anyLld = true;

        if (anyMulti) {
            if (!anyLld && M == 2 && allClosed) {
                // Two stations and no lattice yet: express every station as
                // mu(n) = min(n, c), which is what the reference installs here.
                for (std::size_t i = 0; i < M; ++i) {
                    std::vector<T> lld(Nt, one);
                    for (std::size_t n = 1; n <= Nt; ++n)
                        lld[n - 1] = num_traits<T>::from_double(
                            std::min<double>(static_cast<double>(n), nservers[i]));
                    sn.stations[i].lldscaling = lld;
                }
                anyLld = true;
            } else if (!detail::lld_encodes_multiserver(sn)) {
                throw UnsupportedError(
                    "solver_ncld: the load-dependent solver does not support multi-server "
                    "stations unless they are expressed as limited load dependence mu(n) = "
                    "min(n, c)");
            }
        }

        // The rate lattice, defaulting to a constant rate of one. Its width is the
        // wider of the closed population and the longest row the caller gave
        // setLoadDependence: the mixed route below reads the saturation level off
        // the row itself, and a row cut at Nt is read as a slower station. Past a
        // row's own end the limited load dependence holds its last rate.
        std::size_t Wlld = Nt;
        for (std::size_t i = 0; i < M; ++i)
            Wlld = std::max(Wlld, sn.stations[i].lldscaling.size());
        Matrix<T> lldscaling(M, Wlld, one);
        for (std::size_t i = 0; i < M; ++i) {
            const std::vector<T>& lld = sn.stations[i].lldscaling;
            if (lld.empty()) continue;
            for (std::size_t n = 0; n < Wlld; ++n) lldscaling(i, n) = n < lld.size() ? lld[n] : lld.back();
        }

        mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        Matrix<T> ST = d.ST, ST0 = d.ST;
        const Matrix<T> V = detail::station_visits(sn);
        Matrix<T> SCVnan = sn.scv;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                if (sn.disabled[i][k])
                    SCVnan(i, k) =
                        num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());

        const std::vector<double> Nchain = detail::chain_population(sn);
        std::vector<std::size_t> openChains, closedChains;
        for (std::size_t c = 0; c < C; ++c)
            (std::isinf(Nchain[c]) ? openChains : closedChains).push_back(c);
        std::vector<int> Nnc = detail::nc_population(Nchain);
        std::size_t Ncl = 0;
        for (std::size_t c : closedChains) Ncl += static_cast<std::size_t>(Nnc[c]);

        const pfqn::NcldMethod pmethod = detail::ncld_pfqn_method(opt.method);
        // The estimator controls the log-domain ladder reads; the exact ladder
        // ignores them, so this is passed unconditionally.
        pfqn::NcOptions nopt;
        nopt.samples = opt.samples;
        nopt.seed = opt.seed;
        nopt.tol = opt.tol;
        const T atol = num_traits<T>::from_double(opt.tol);

        std::vector<T> gamma(M, zero), nserv_t(M, one);
        for (std::size_t i = 0; i < M; ++i) nserv_t[i] = num_traits<T>::from_double(nservers[i]);
        std::vector<T> eta(M, one), eta_1(M, zero);
        const int iter_max = anyFCFS ? opt.iter_max : 1;
        int it = 0;
        double lG = 0.0;
        std::string actualmethod = opt.method;
        mva::ClassResults<T> cls;
        std::vector<T> Xchain(C, zero);

        while (it < iter_max) {
            {
                double dev = 0.0;
                for (std::size_t i = 0; i < M; ++i) {
                    const double v = std::fabs(1.0 - num_traits<T>::to_double(eta[i]) /
                                                         num_traits<T>::to_double(eta_1[i]));
                    if (!(v <= dev)) dev = v;
                }
                if (!(dev > opt.iter_tol)) break;
            }
            ++it;
            eta_1 = eta;

            // Chain demands from the CURRENT service times.
            for (std::size_t c = 0; c < C; ++c) {
                const bool open = std::isinf(Nchain[c]);
                const std::size_t rst = sn.classes[sn.inchain[c][0] - 1].refstat;
                for (std::size_t i = 0; i < M; ++i) {
                    T st = zero;
                    for (std::size_t k : sn.inchain[c]) st += ST(i, k - 1) * d.alpha(i, k - 1);
                    d.Lchain(i, c) = T(d.Vchain(i, c) * st);
                    if (open && i + 1 == rst) {
                        // A source row carries 1 / arrival rate, summed over the
                        // classes whose rate is defined.
                        T s = zero;
                        for (std::size_t k : sn.inchain[c])
                            if (!sn.disabled[i][k - 1] &&
                                std::isfinite(num_traits<T>::to_double(ST(i, k - 1))))
                                s += ST(i, k - 1);
                        d.STchain(i, c) = s;
                    } else {
                        d.STchain(i, c) = st;
                    }
                }
            }
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < C; ++c) {
                    if (!std::isfinite(num_traits<T>::to_double(d.STchain(i, c))))
                        d.STchain(i, c) = zero;
                    if (!std::isfinite(num_traits<T>::to_double(d.Lchain(i, c))))
                        d.Lchain(i, c) = zero;
                }

            std::vector<T> lambda(C, zero);
            for (std::size_t c : openChains) {
                const std::size_t rst = sn.classes[sn.inchain[c][0] - 1].refstat;
                if (d.STchain(rst - 1, c) > zero) lambda[c] = T(one / d.STchain(rst - 1, c));
            }

            Matrix<T> L(M, C, zero), Z(M, C, zero), mu(M, Nt, one);
            std::vector<std::size_t> infServers;
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t c = 0; c < C; ++c) L(i, c) = d.Lchain(i, c);
                if (std::isinf(nservers[i])) {
                    infServers.push_back(i);
                    for (std::size_t c = 0; c < C; ++c) Z(i, c) = d.Lchain(i, c);
                    for (std::size_t n = 1; n <= Nt; ++n)
                        mu(i, n - 1) = num_traits<T>::from_double(static_cast<double>(n));
                } else {
                    for (std::size_t n = 0; n < Nt; ++n) mu(i, n) = lldscaling(i, n);
                }
            }

            Matrix<T> Qchain(M, C, zero);
            Xchain.assign(C, zero);

            if (!openChains.empty()) {
                // Mixed limited load-dependent network: the chain-level
                // normalizing constant of Bruell-Balbo-Afshari effective
                // capacity (pfqn_ncldmx), which never enumerates the closed
                // population lattice. The source rows carry only the 1/lambda bookkeeping
                // demand and are excluded from the queueing set; the delay rows
                // fold into the think-time vector.
                std::vector<bool> isSource(M, false), isDelay(M, false);
                for (std::size_t c : openChains)
                    isSource[sn.classes[sn.inchain[c][0] - 1].refstat - 1] = true;
                for (std::size_t i : infServers)
                    if (!isSource[i]) isDelay[i] = true;
                std::vector<std::size_t> queueStations;
                for (std::size_t i = 0; i < M; ++i)
                    if (!isSource[i] && !isDelay[i]) queueStations.push_back(i);
                const std::size_t nq = queueStations.size();

                Matrix<T> Zvec(1, C, zero);
                for (std::size_t i = 0; i < M; ++i)
                    if (isDelay[i])
                        for (std::size_t c = 0; c < C; ++c) Zvec(0, c) += d.Lchain(i, c);

                // pfqn_ldmx_ec reads the limited-load-dependence level b_i off the
                // rate row itself -- the first column equal to the LAST one -- and
                // treats every rate past it as saturated. Cutting the row at the
                // closed population Ncl declares a c-server station saturated at
                // min(n,c) with n < c whenever c exceeds Ncl (and with no closed
                // class at all it flattens the row to mu(1)), so keep every column
                // up to the start of each row's trailing constant run.
                std::size_t ncol = std::max<std::size_t>(1, Ncl);
                for (std::size_t qi = 0; qi < nq; ++qi)
                    ncol = std::max(ncol, detail::lld_saturation_level(lldscaling, queueStations[qi]));
                Matrix<T> Dq(nq, C, zero), muq(nq, ncol, one);
                const std::size_t Wq = lldscaling.cols();
                for (std::size_t qi = 0; qi < nq; ++qi) {
                    for (std::size_t c = 0; c < C; ++c) Dq(qi, c) = d.Lchain(queueStations[qi], c);
                    for (std::size_t n = 0; n < ncol; ++n)
                        muq(qi, n) = lldscaling(queueStations[qi], n < Wq ? n : Wq - 1);
                }
                std::vector<int> Nmx(C, 0);
                for (std::size_t c = 0; c < C; ++c)
                    Nmx[c] = std::isinf(Nchain[c]) ? pfqn::OPEN_CLASS : Nnc[c];
                const pfqn::NcldmxResult<T> mx =
                    pfqn::pfqn_ncldmx(lambda, Dq, Nmx, Zvec, muq, pmethod, atol, nopt);
                lG = mx.lG;
                Xchain = mx.XN;
                for (std::size_t qi = 0; qi < nq; ++qi)
                    for (std::size_t c = 0; c < C; ++c) Qchain(queueStations[qi], c) = mx.QN(qi, c);
                for (std::size_t i = 0; i < M; ++i)
                    if (isDelay[i])
                        for (std::size_t c = 0; c < C; ++c)
                            Qchain(i, c) = T(d.Lchain(i, c) * Xchain[c]);
                actualmethod = "ncldmx";
            } else {
                Matrix<T> Zzero(1, C, zero);
                const pfqn::NcldResult<T> base =
                    pfqn::pfqn_ncld(L, Nnc, Zzero, mu, pmethod, atol, nopt);
                lG = base.lG;
                actualmethod = base.method;

                const bool repairman = (M == 2 && !infServers.empty());
                std::size_t firstDelay = infServers.empty() ? 0 : infServers[0];
                for (std::size_t r = 0; r < C; ++r) {
                    const std::vector<int> Nr = detail::oner(Nnc, r);
                    const double lGr = pfqn::pfqn_ncld(L, Nr, Zzero, mu, pmethod, atol, nopt).lG;
                    Xchain[r] = num_traits<T>::from_double(std::exp(lGr - lG));
                    if (repairman) {
                        const T qd = T(d.Lchain(firstDelay, r) * Xchain[r]);
                        Qchain(firstDelay, r) = qd;
                        for (std::size_t i = 0; i < M; ++i)
                            if (i != firstDelay)
                                Qchain(i, r) = T(num_traits<T>::from_double(Nchain[r]) - qd);
                        continue;
                    }
                    std::size_t nfinite = 0;
                    for (std::size_t i = 0; i < M; ++i)
                        if (std::isfinite(nservers[i])) ++nfinite;
                    for (std::size_t i = 0; i < M; ++i) {
                        if (!(d.Lchain(i, r) > zero)) continue;
                        if (std::isinf(nservers[i])) {
                            Qchain(i, r) = T(d.Lchain(i, r) * Xchain[r]);
                            continue;
                        }
                        if (i + 1 == M && nfinite == 1) {
                            // The only queueing station: give it the balance of
                            // the population rather than a fourth constant.
                            T acc = zero;
                            for (std::size_t i2 : infServers) acc += d.Lchain(i2, r);
                            T q = T(num_traits<T>::from_double(Nchain[r]) - acc * Xchain[r]);
                            for (std::size_t i2 = 0; i2 + 1 < M; ++i2)
                                if (std::isfinite(nservers[i2])) q -= Qchain(i2, r);
                            Qchain(i, r) = q < zero ? zero : q;
                            continue;
                        }
                        // Conditional MVA: the shifted lattice at station i, its
                        // flow-equivalent complement, and the model with station
                        // i removed.
                        const Matrix<T> muhati = pfqn::pfqn_mushift(mu, i);
                        Matrix<T> muhati_row(1, muhati.cols(), zero);
                        for (std::size_t n = 0; n < muhati.cols(); ++n) muhati_row(0, n) = muhati(i, n);
                        const pfqn::FncResult<T> fnc = pfqn::pfqn_fnc(muhati_row);
                        Matrix<T> Lhat(M + 1, C, zero), muhat(M + 1, muhati.cols(), zero);
                        for (std::size_t i2 = 0; i2 < M; ++i2) {
                            for (std::size_t c = 0; c < C; ++c) Lhat(i2, c) = L(i2, c);
                            for (std::size_t n = 0; n < muhati.cols(); ++n)
                                muhat(i2, n) = muhati(i2, n);
                        }
                        for (std::size_t c = 0; c < C; ++c) Lhat(M, c) = L(i, c);
                        for (std::size_t n = 0; n < fnc.mu.cols(); ++n) muhat(M, n) = fnc.mu(0, n);
                        Matrix<T> Lms_i(M - 1, C, zero), mu_i(M - 1, mu.cols(), zero);
                        std::size_t row = 0;
                        for (std::size_t i2 = 0; i2 < M; ++i2) {
                            if (i2 == i) continue;
                            for (std::size_t c = 0; c < C; ++c) Lms_i(row, c) = L(i2, c);
                            for (std::size_t n = 0; n < mu.cols(); ++n) mu_i(row, n) = mu(i2, n);
                            ++row;
                        }
                        const double lGhat_fnci =
                            pfqn::pfqn_ncld(Lhat, Nr, Zzero, muhat, pmethod, atol, nopt).lG;
                        const double lGhatir =
                            pfqn::pfqn_ncld(L, Nr, Zzero, muhati, pmethod, atol, nopt).lG;
                        const double lGr_i =
                            pfqn::pfqn_ncld(Lms_i, Nr, Zzero, mu_i, pmethod, atol, nopt).lG;
                        const double dlGa = lGhat_fnci - lGhatir;
                        const double dlG_i = lGr_i - lGhatir;
                        const T CQ = T(num_traits<T>::from_double(std::exp(dlGa) - 1.0) +
                                       fnc.c[0] * num_traits<T>::from_double(std::exp(dlG_i) - 1.0));
                        const double ldDemand = std::log(num_traits<T>::to_double(L(i, r))) +
                                                lGhatir -
                                                std::log(num_traits<T>::to_double(mu(i, 0))) - lGr;
                        Qchain(i, r) = T(num_traits<T>::from_double(std::exp(ldDemand)) *
                                         Xchain[r] * (one + CQ));
                    }
                }
            }

            Matrix<T> Rchain(M, C, zero), Tchain(M, C, zero);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < C; ++c) {
                    if (Xchain[c] != zero && d.Vchain(i, c) != zero)
                        Rchain(i, c) = T(Qchain(i, c) / Xchain[c] / d.Vchain(i, c));
                    Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
                }
            for (std::size_t i : infServers)
                for (std::size_t c = 0; c < C; ++c)
                    Rchain(i, c) =
                        d.Vchain(i, c) == zero ? zero : T(d.Lchain(i, c) / d.Vchain(i, c));
            for (std::size_t c = 0; c < C; ++c) {
                if (Nnc[c] != 0) continue;
                Xchain[c] = zero;
                for (std::size_t i = 0; i < M; ++i) {
                    Qchain(i, c) = zero;
                    Rchain(i, c) = zero;
                    Tchain(i, c) = zero;
                }
            }
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < C; ++c) {
                    if (!std::isfinite(num_traits<T>::to_double(Qchain(i, c)))) Qchain(i, c) = zero;
                    if (!std::isfinite(num_traits<T>::to_double(Rchain(i, c)))) Rchain(i, c) = zero;
                }

            d.ST = ST;
            cls = mva::sn_deaggregate_chain_results(sn, d, Matrix<T>(), Matrix<T>(), Rchain,
                                                    Tchain, Xchain);
            out.STeff = ST;

            const npfqn::NonexpApproxResult<T> na = npfqn::npfqn_nonexp_approx(
                opt.highvar, isFCFS, sn.rates, ST0, V, SCVnan, cls.Tp, cls.U, gamma, nserv_t);
            ST = na.ST;
            gamma = na.gamma;
            eta = na.eta;
        }

        Matrix<T> Q = cls.Q, U = cls.U, R = cls.R, Tp = cls.Tp;
        std::vector<T> X = cls.X;
        auto absify = [&](Matrix<T>& A) {
            for (std::size_t i = 0; i < A.rows(); ++i)
                for (std::size_t j = 0; j < A.cols(); ++j)
                    if (A(i, j) < zero) A(i, j) = T(-A(i, j));
        };
        absify(Q);
        absify(R);
        absify(U);
        for (T& x : X)
            if (x < zero) x = T(-x);

        // Utilization is re-derived rather than taken from the deaggregation:
        // at a load-dependent station the chain-level U carries the lattice, and
        // what the reference reports is the fraction of the SERVERS busy.
        for (std::size_t i = 0; i < M; ++i) {
            const bool multi = std::isfinite(nservers[i]) && nservers[i] > 1.0;
            if (multi || std::isinf(nservers[i])) {
                const double div = multi ? nservers[i] : 1.0;
                for (std::size_t k = 0; k < K; ++k) {
                    std::size_t c = C;
                    for (std::size_t cc = 0; cc < C; ++cc)
                        if (sn.chains[cc][k]) c = cc;
                    if (c == C) continue;
                    const std::size_t rs = sn.classes[k].refstat;
                    const T vref = sn.visits[c](sn.stateful_of_station(rs) - 1, k);
                    if (vref == zero) continue;
                    const T vi = sn.visits[c](sn.stateful_of_station(i + 1) - 1, k);
                    const bool open = std::isinf(sn.classes[k].population);
                    T rate = zero;
                    if (open) {
                        // the open class's own arrival rate, sn's source rate
                        const std::size_t src = sn.classes[k].refstat;
                        if (!sn.disabled[src - 1][k]) rate = sn.rates(src - 1, k);
                    } else {
                        rate = X[k];
                    }
                    if (!(rate > zero)) continue;
                    U(i, k) = T(rate * vi / vref * ST(i, k) / num_traits<T>::from_double(div));
                }
            } else {
                // `solver_ncld.m:300` divides by `max(lldscaling(ist,:))` over the
                // STATION'S OWN ROW, whose width is whatever the caller gave
                // setLoadDependence -- not over the Nt-column lattice this file
                // builds. Taking the max over the lattice dropped every entry past
                // the population: on a 4-job model with scaling [1 1.6 2 2.2 2.3]
                // the divisor came out 2.2 instead of 2.3 and Util read 0.419042
                // against the reference's 0.400823, with QLen, RespT and Tput all
                // exact. A station with no row of its own keeps the lattice's
                // constant one, which is what the row would hold anyway.
                T mx = zero;
                const std::vector<T>& lldrow = sn.stations[i].lldscaling;
                if (lldrow.empty()) {
                    for (std::size_t n = 0; n < Nt; ++n)
                        if (lldscaling(i, n) > mx) mx = lldscaling(i, n);
                } else {
                    for (std::size_t n = 0; n < lldrow.size(); ++n)
                        if (lldrow[n] > mx) mx = lldrow[n];
                }
                if (mx > zero)
                    for (std::size_t k = 0; k < K; ++k) U(i, k) = T(U(i, k) / mx);
                T s = zero;
                for (std::size_t k = 0; k < K; ++k) s += U(i, k);
                if (num_traits<T>::to_double(s) > 1.0)
                    for (std::size_t k = 0; k < K; ++k) U(i, k) = T(U(i, k) / s);
            }
        }

        auto clear_nonfinite = [&](Matrix<T>& A) {
            for (std::size_t i = 0; i < A.rows(); ++i)
                for (std::size_t j = 0; j < A.cols(); ++j)
                    if (!std::isfinite(num_traits<T>::to_double(A(i, j)))) A(i, j) = zero;
        };
        clear_nonfinite(Q);
        clear_nonfinite(U);
        clear_nonfinite(R);
        for (T& x : X)
            if (!std::isfinite(num_traits<T>::to_double(x))) x = zero;

        for (std::size_t c = 0; c < C; ++c) {
            if (std::isinf(Nchain[c])) continue;
            T qden = zero;
            for (std::size_t k : sn.inchain[c])
                for (std::size_t i = 0; i < M; ++i) qden += Q(i, k - 1);
            const T ratio =
                qden > zero ? T(num_traits<T>::from_double(Nchain[c]) / qden) : zero;
            for (std::size_t k : sn.inchain[c]) {
                X[k - 1] = T(ratio * X[k - 1]);
                for (std::size_t i = 0; i < M; ++i) {
                    Q(i, k - 1) = T(ratio * Q(i, k - 1));
                    Tp(i, k - 1) = T(ratio * Tp(i, k - 1));
                    U(i, k - 1) = T(ratio * U(i, k - 1));
                    R(i, k - 1) = Tp(i, k - 1) == zero ? zero : T(Q(i, k - 1) / Tp(i, k - 1));
                }
            }
        }

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = R;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C.assign(K, zero);
        for (std::size_t k = 0; k < K; ++k) {
            const double njobs = sn.classes[k].population;
            out.sol.C[k] = (std::isfinite(njobs) && X[k] != zero)
                               ? T(num_traits<T>::from_double(njobs) / X[k])
                               : cls.C[k];
        }
        out.sol.lG = lG;
        out.sol.iter = it;
        out.sol.method = actualmethod;
        out.actualmethod = actualmethod;
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NCLD_H
