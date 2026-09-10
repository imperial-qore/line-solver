/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SENS_SOLVER_SENS_TABLE_H
#define LINE_SOLVERS_SENS_SOLVER_SENS_TABLE_H

/**
 * Performance sensitivities with respect to service rates.
 *
 * Port of `matlab/src/solvers/@@NetworkSolver/getSensitivityTable.m`. One row
 * per (station, class) carrying dTput/dRate, dRespT/dRate, dQLen/dRate and
 * dUtil/dRate, produced by one of two branches:
 *
 *   `exact` Analytic differentiation of a product-form recursion. A closed
 *           model goes through `pfqn_sens` (differentiated MVA) evaluated AT
 *           CHAIN LEVEL, an open one through the closed-form BCMP derivatives,
 *           whose stations decouple. Exact to working precision and cheaper
 *           than one extra solve, but defined only where the recursion is:
 *           single-server stations plus delays, and not a mixed model.
 *
 *   `fd`    Forward or central differences on the CALLER'S OWN solver. The
 *           service at (station, class) is rate-scaled by (1 +/- h), the same
 *           solve is re-run and the quotient formed. It costs 1 + M*R solves
 *           (forward) or 2*M*R (central) and is the only branch that applies to
 *           a model with no product form.
 *
 * WHY THE SOLVER IS A CALLBACK. The reference reaches this through a method on
 * NetworkSolver, so `self` names both the model and the engine; there is no
 * solver base class here, and the finite-difference branch needs to re-run
 * THE SAME engine with THE SAME options -- an MVA table produced by perturbing
 * a model that is then solved by CTMC is not a sensitivity of anything. The
 * callback is that engine, already bound to the struct passed in, and it is
 * called after this code has written the perturbed service into it.
 *
 * A SIMULATION SOLVER MUST BE RUN WITH COMMON RANDOM NUMBERS, or the quotient
 * measures Monte Carlo error rather than a derivative. There is no seed to pin
 * here -- the callback owns its options -- so `SensOptions::simulation` only
 * widens the default step to 1e-2, and a caller wiring a stochastic engine in
 * is responsible for handing it a fixed seed.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_sens.h"
#include "line/lang/dist_scale_rate.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/nc/sn_pf_params.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace sens {

/** The name-value contract of getSensitivityTable. */
struct SensOptions {
    std::string method = "auto";     ///< auto | exact | fd
    std::string scheme = "forward";  ///< forward | central
    /** Relative step of the rate perturbation; negative selects the default. */
    double step = -1.0;
    /** True when the callback is a simulator, which widens the default step. */
    bool simulation = false;
};

/** One (station, class) row of the table. */
template <class T>
struct SensRow {
    std::string station, jobclass;
    T dTput, dRespT, dQLen, dUtil;
};

/** What the table carries, plus the branch that produced it. */
template <class T>
struct SensTable {
    std::vector<SensRow<T>> rows;
    std::string method;  ///< "exact" or "fd", the branch actually taken
    /** The analytic Jacobian, set only on the CLOSED exact branch. */
    bool has_jacobian = false;
    pfqn::SensResult<T> jacobian;
};

namespace detail {

/**
 * Scope of the analytic branch: single-server stations (a delay is allowed) and
 * not a mixed open-and-closed model. Class switching IS supported -- the branch
 * aggregates classes into chains before differentiating.
 */
template <class T>
bool sens_exact_in_scope(const qn::NetworkStruct<T>& sn, std::string& why) {
    nc::PfParams<T> p;
    try {
        p = nc::sn_get_product_form_params(sn);
    } catch (const std::exception&) {
        why = "the product-form parameters of this model could not be extracted";
        return false;
    }
    for (double s : p.S)
        if (std::isfinite(s) && s > 1.0) {
            why = "exact sensitivities support single-server stations only";
            return false;
        }
    bool any_open = false, any_closed = false;
    for (double n : p.N) {
        if (std::isinf(n)) any_open = true;
        else any_closed = true;
    }
    if (any_open && any_closed) {
        why = "exact sensitivities do not yet support mixed (open+closed) networks";
        return false;
    }
    return true;
}

/** (nstations x nclasses) mask of the pairs that carry visits. */
template <class T>
std::vector<std::vector<bool>> sens_visit_mask(const qn::NetworkStruct<T>& sn) {
    std::vector<std::vector<bool>> out(sn.nstations, std::vector<bool>(sn.nclasses, false));
    for (std::size_t c = 0; c < sn.nchains; ++c) {
        if (sn.visits[c].empty()) continue;
        for (std::size_t i = 1; i <= sn.nstations; ++i) {
            // Resolved through the node, not stateful_of_station: that helper
            // THROWS on a station with no stateful node, and a mask is not the
            // place to discover one.
            const std::size_t nd = sn.node_of_station(i);
            const std::size_t sf = nd ? sn.stateful_index(nd) : 0;
            if (sf == 0) continue;
            for (std::size_t r = 0; r < sn.nclasses; ++r)
                if (num_traits<T>::to_double(sn.visits[c](sf - 1, r)) > lang::GlobalConstants::Zero)
                    out[i - 1][r] = true;
        }
    }
    return out;
}

}  // namespace detail

/**
 * Build the sensitivity table of `sn` under `solve`.
 *
 * `exact_available` is the reference's `supportsExactSensitivity()`: true for
 * the engines that evaluate the product-form recursion this branch
 * differentiates (MVA and NC), false for every other. Asking for `exact` where
 * it is false is an error rather than a silent downgrade, exactly as in the
 * reference, because the two branches answer to different precision.
 *
 * `solve` is called with `sn` already carrying whatever perturbation this code
 * has written, and must report the metrics of that struct.
 */
template <class T>
SensTable<T> solver_sensitivity_table(qn::NetworkStruct<T>& sn, const SensOptions& opt,
                                      bool exact_available,
                                      const std::function<mva::MvaSolution<T>()>& solve) {
    using lang::Distrib;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (opt.method != "auto" && opt.method != "exact" && opt.method != "fd")
        throw InputError("getSensitivityTable: the method must be 'auto', 'exact' or 'fd'");
    if (opt.scheme != "forward" && opt.scheme != "central")
        throw InputError("getSensitivityTable: the scheme must be 'forward' or 'central'");

    std::string why;
    const bool in_scope = detail::sens_exact_in_scope(sn, why);
    bool use_exact;
    if (opt.method == "exact") {
        if (!exact_available)
            throw UnsupportedError(
                "getSensitivityTable: exact analytic sensitivities differentiate a product-form "
                "recursion and are available on the MVA and NC engines only; use 'fd'");
        if (!in_scope) throw UnsupportedError("getSensitivityTable: " + why);
        use_exact = true;
    } else if (opt.method == "fd") {
        use_exact = false;
    } else {
        use_exact = exact_available && in_scope;
    }

    const nc::PfParams<T> p = nc::sn_get_product_form_params(sn);
    const std::size_t R = sn.nclasses, Mq = p.queue_stations.size(), C = sn.nchains;

    std::vector<std::vector<bool>> mask(Mq, std::vector<bool>(R, false));
    Matrix<T> dT(Mq, R, zero), dR(Mq, R, zero), dQ(Mq, R, zero), dU(Mq, R, zero);
    SensTable<T> out;

    if (use_exact) {
        out.method = "exact";
        bool is_open = false;
        for (double n : p.N)
            if (std::isinf(n)) is_open = true;

        Matrix<T> rates(Mq, R, zero);
        for (std::size_t i = 0; i < Mq; ++i) {
            const std::size_t st = p.queue_stations[i];
            for (std::size_t r = 0; r < R; ++r) {
                rates(i, r) = sn.rates(st - 1, r);
                const double rd = num_traits<T>::to_double(rates(i, r));
                mask[i][r] = std::isfinite(rd) && rd > 0.0 && p.D(i, r) > zero;
            }
        }

        // chainOf(r), and the per-chain aggregates the recursion is evaluated on:
        // a product-form model is solved per CHAIN, and with class switching a
        // class is only a share of one.
        std::vector<std::size_t> chain_of(R, 0);
        Matrix<T> Dc(Mq, C, zero);
        std::vector<T> Zc(C, zero);
        std::vector<int> Nc(C, 0);
        std::vector<T> lambdac(C, zero);
        for (std::size_t c = 0; c < C; ++c) {
            double nsum = 0.0;
            for (std::size_t r = 0; r < R; ++r) {
                if (!sn.chains[c][r]) continue;
                chain_of[r] = c;
                for (std::size_t i = 0; i < Mq; ++i) Dc(i, c) = T(Dc(i, c) + p.D(i, r));
                for (std::size_t z = 0; z < p.Z.rows(); ++z) Zc[c] = T(Zc[c] + p.Z(z, r));
                if (std::isfinite(p.N[r])) nsum += p.N[r];
                lambdac[c] = T(lambdac[c] + p.lambda[r]);
            }
            Nc[c] = static_cast<int>(nsum + 0.5);
        }

        if (!is_open) {
            const pfqn::SensResult<T> s = pfqn::pfqn_sens(Dc, Nc, Zc);
            out.has_jacobian = true;
            out.jacobian = s;
            for (std::size_t i = 0; i < Mq; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    if (!mask[i][r]) continue;
                    const std::size_t c = chain_of[r];
                    if (!(Dc(i, c) > zero)) continue;
                    const T rate = rates(i, r);
                    const T Dir = p.D(i, r);
                    const T visits = T(Dir * rate);   // chain-normalized visit ratio
                    const std::size_t pidx = i * C + c;
                    const T chain = T(-Dir / rate);   // dDc(i,c)/dmu(i,r)
                    const T Xc = s.XN[c];
                    const T Qc = s.QN(i, c);
                    const T dXc = T(s.dX(c, pidx) * chain);
                    const T dQc = T(s.dQ[pidx](i, c) * chain);
                    // Class share of the chain queue here, and its own rate dependence.
                    const T alpha = T(Dir / Dc(i, c));
                    const T dalpha = T(chain * (Dc(i, c) - Dir) / (Dc(i, c) * Dc(i, c)));
                    const T Qir = T(alpha * Qc);
                    const T dQir = T(dalpha * Qc + alpha * dQc);
                    const T Tir = T(Xc * visits);
                    const T dTir = T(dXc * visits);
                    dT(i, r) = dTir;
                    dQ(i, r) = dQir;
                    dU(i, r) = T(dXc * Dir + Xc * chain);
                    // Per-visit response time by Little's law, R = Q/T.
                    if (Tir > zero) dR(i, r) = T((dQir * Tir - Qir * dTir) / (Tir * Tir));
                }
        } else {
            // The stations decouple, so only the own service rate moves the
            // measures at (i, r); the throughput is fixed by the arrival rate.
            Matrix<T> rho(Mq, R, zero);
            std::vector<T> Ui(Mq, zero);
            for (std::size_t i = 0; i < Mq; ++i)
                for (std::size_t r = 0; r < R; ++r) {
                    if (p.D(i, r) > zero) rho(i, r) = T(lambdac[chain_of[r]] * p.D(i, r));
                    Ui[i] = T(Ui[i] + rho(i, r));
                }
            for (std::size_t i = 0; i < Mq; ++i) {
                const T denom = T(one - Ui[i]);
                for (std::size_t r = 0; r < R; ++r) {
                    if (!mask[i][r]) continue;
                    const T rate = rates(i, r);
                    const T svct = T(one / rate);  // per-visit service time
                    const T drho = T(-rho(i, r) / rate);
                    const T dUi = drho;            // own class only
                    const T dsvct = T(-svct / rate);
                    dR(i, r) = T((dsvct * denom + svct * dUi) / (denom * denom));
                    dQ(i, r) = T((drho * denom + rho(i, r) * dUi) / (denom * denom));
                    dU(i, r) = drho;
                    dT(i, r) = zero;  // open throughput = lambda*visits, fixed
                }
            }
        }
    } else {
        out.method = "fd";
        const bool central = opt.scheme == "central";
        double h = opt.step;
        if (!(h > 0.0)) h = opt.simulation ? 1e-2 : 1e-4;
        if (!std::isfinite(h) || h <= 0.0 || h >= 1.0)
            throw InputError("getSensitivityTable: the finite-difference step must be in (0,1)");
        const T hT = num_traits<T>::from_double(h);

        const std::vector<std::vector<bool>> visited = detail::sens_visit_mask(sn);
        for (std::size_t i = 0; i < Mq; ++i) {
            const std::size_t st = p.queue_stations[i];
            for (std::size_t r = 0; r < R; ++r) {
                const double rd = num_traits<T>::to_double(sn.rates(st - 1, r));
                mask[i][r] = std::isfinite(rd) && rd > 0.0 && visited[st - 1][r];
            }
        }

        const mva::MvaSolution<T> base = solve();
        for (std::size_t i = 0; i < Mq; ++i) {
            const std::size_t st = p.queue_stations[i];
            for (std::size_t r = 0; r < R; ++r) {
                if (!mask[i][r]) continue;
                const T rate = sn.rates(st - 1, r);
                const Distrib<T> saved = sn.service[st - 1][r];

                sn.set_service(st, r + 1, lang::dist_scale_rate(saved, T(one + hT)));
                sn.refresh_rates();
                const mva::MvaSolution<T> up = solve();

                mva::MvaSolution<T> down = base;
                T denom = T(rate * hT);
                if (central) {
                    sn.set_service(st, r + 1, lang::dist_scale_rate(saved, T(one - hT)));
                    sn.refresh_rates();
                    down = solve();
                    denom = T(num_traits<T>::from_int(2) * rate * hT);
                }

                dT(i, r) = T((up.Tp(st - 1, r) - down.Tp(st - 1, r)) / denom);
                dR(i, r) = T((up.R(st - 1, r) - down.R(st - 1, r)) / denom);
                dQ(i, r) = T((up.Q(st - 1, r) - down.Q(st - 1, r)) / denom);
                dU(i, r) = T((up.U(st - 1, r) - down.U(st - 1, r)) / denom);

                // Restore before moving on: the sweep perturbs one pair at a
                // time and every later quotient is taken at the base point.
                sn.set_service(st, r + 1, saved);
                sn.refresh_rates();
            }
        }
    }

    for (std::size_t i = 0; i < Mq; ++i) {
        const std::size_t st = p.queue_stations[i];
        const std::size_t nd = sn.node_of_station(st);
        for (std::size_t r = 0; r < R; ++r) {
            if (!mask[i][r]) continue;
            SensRow<T> row;
            row.station = nd > 0 ? sn.nodes[nd - 1].name : sn.stations[st - 1].name;
            row.jobclass = sn.classes[r].name;
            row.dTput = dT(i, r);
            row.dRespT = dR(i, r);
            row.dQLen = dQ(i, r);
            row.dUtil = dU(i, r);
            out.rows.push_back(row);
        }
    }
    return out;
}

}  // namespace sens
}  // namespace line

#endif  // LINE_SOLVERS_SENS_SOLVER_SENS_TABLE_H
