/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_SYMODES_H
#define LINE_SOLVERS_FLUID_FLUID_SYMODES_H

/**
 * A symbolic description of the fluid ODE system: a port of
 * `solver_fluid_symodes.m`, which is what `@@SolverFLD/exportODEs` renders.
 *
 * WHY THE SOLVER CANNOT SIMPLY BE ASKED. The numerical path evaluates the drift
 * at a point; it never holds the drift as an object. Anything that wants to
 * PRINT the system, differentiate it by hand, or hand it to another tool needs
 * the coefficients and the shape of every term, and that is what this builds.
 * It is deliberately a SECOND construction of the same dynamics: it mirrors the
 * numerical code path term by term, so a divergence between the two is a real
 * defect and not a rendering artefact.
 *
 * TWO FORMS, MATCHING THE TWO SOLVER PATHS.
 *
 *   form W (methods default, matrix, pnorm)
 *       dx/dt = W' theta(x) + lambda
 *       theta_s(x) = x_s min(n_i, S_i)/n_i, or the p-norm smoothing of it, and
 *       theta_s = 0 at a Source. This is `solver_fluid_matrix` and is taken
 *       from the SAME assembly the integrator uses, so the two cannot drift.
 *
 *   form J (methods closing, statedep, softmin)
 *       dx/dt = J r(x), r_e(x) = coeff_e * factor_e(x)
 *       with one event per (departure, phase change) and a factor whose shape
 *       is fixed by the scheduling strategy of the station that drives it.
 *
 * THE FACTOR TYPES, which are the whole content of the J form:
 *   lin      x_v                                    (INF, and unhandled policies)
 *   min      x_v min(n_i, S_i)/n_i                  (PS, FCFS under closing)
 *   ext1     1 - sum of the other phases at a Source
 *   dpsmin   x_v min(n_i, S_i)/ntilde_i, the weight w_ir folded into coeff
 *   dpspw    piecewise: x_v below S_i, weighted above it
 *   fcfsw    x_v min(n_i, S_i)/nhat_i, phase weight folded into coeff
 *   fcfsws   the same with softmin in place of min
 *
 * WHAT IS REFUSED BY NAME: `tbi`, `diffusion` and `mfq` have no ODE system of
 * this shape at all -- tbi partitions and re-solves, diffusion adds a noise
 * term, mfq solves a queue analytically -- and `statedep`/`softmin` have no
 * open-model branch, exactly as in the reference.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_matrix.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** The state-dependent factor attached to one event's driving variable. */
struct SymFactor {
    std::string type;                  ///< lin, min, pnorm, ext1, dpsmin, dpspw, fcfsw, fcfsws
    std::size_t station = 0;           ///< 1-based
    std::size_t cls = 0;               ///< 1-based
    std::vector<std::size_t> others;   ///< ext1: the other phases, 1-based state indices
};

/** The symbolic system, in whichever of the two forms the method implies. */
struct FluidSymSystem {
    std::string form;    ///< "W" or "J"
    std::string method;  ///< resolved: default becomes matrix
    std::vector<std::string> station_names, class_names, sched_names;
    std::vector<lang::SchedStrategy> sched;
    std::size_t nstates = 0;
    std::vector<std::size_t> state_station, state_class, state_phase;  ///< 1-based
    std::vector<double> S;  ///< servers per station, infinite already substituted
    std::vector<double> x0;

    // form W
    Matrix<double> W;
    std::vector<double> alambda;
    std::vector<bool> is_source;
    std::vector<bool> is_inf;          ///< (n) state belongs to an INF station
    std::string smoothing;  ///< "min" or "pnorm"
    std::vector<double> pstar;

    // form J
    Matrix<double> J;
    std::vector<double> coeff;
    std::vector<std::size_t> event_var;  ///< 0-based state index driving the event
    std::vector<SymFactor> factor;
    std::size_t nevents = 0;
    Matrix<double> dpsw;
    std::vector<double> fcfs_phase_w;
    double alpha = 0.0;  ///< softmin sharpness, set only for that method
};

namespace detail {

/** The per-strategy factor of an event, and the coefficient folding it implies. */
inline SymFactor sym_factor(const std::string& method, lang::SchedStrategy sched, std::size_t i,
                            std::size_t c, std::size_t ki, double& coeff, const FluidLayout& L,
                            const std::vector<double>& S, const Matrix<double>& dpsw,
                            const std::vector<double>& fcfs_w, std::size_t var) {
    SymFactor f;
    f.station = i + 1;
    f.cls = c + 1;
    if (method == "closing") {
        switch (sched) {
            case lang::SchedStrategy::INF: f.type = "lin"; break;
            case lang::SchedStrategy::EXT:
                if (ki == 0) {
                    f.type = "ext1";
                    for (std::size_t k = 1; k < L.kic[i][c]; ++k)
                        f.others.push_back(L.qidx[i][c] + k + 1);
                } else {
                    f.type = "lin";
                }
                break;
            case lang::SchedStrategy::PS:
            case lang::SchedStrategy::FCFS: f.type = "min"; break;
            case lang::SchedStrategy::DPS:
                // the share w_ir*x/ntilde_i of the capacity min(n_i,S_i), as in
                // fluid_rates_closing_factors: no additive seed, not the full S_i
                f.type = "dpsmin";
                coeff *= dpsw(i, c);
                break;
            default: f.type = "lin"; break;  // no case in ode_rates_closing: rate stays x
        }
        return f;
    }
    // statedep and softmin
    switch (sched) {
        case lang::SchedStrategy::INF: f.type = "lin"; break;
        case lang::SchedStrategy::PS: f.type = "min"; break;
        case lang::SchedStrategy::FCFS:
            f.type = (method == "softmin") ? "fcfsws" : "fcfsw";
            coeff *= fcfs_w[var];
            break;
        case lang::SchedStrategy::DPS: f.type = "dpspw"; break;
        default: f.type = "lin"; break;
    }
    return f;
}

/** True when `statedep`/`softmin` enumerate events out of this station at all. */
inline bool sym_handled_sched(lang::SchedStrategy s) {
    return s == lang::SchedStrategy::INF || s == lang::SchedStrategy::EXT ||
           s == lang::SchedStrategy::PS || s == lang::SchedStrategy::FCFS ||
           s == lang::SchedStrategy::DPS;
}

}  // namespace detail

/**
 * Build the symbolic system of `sn` under `opt.method`.
 *
 * `init_sol` is the initial condition to report; empty takes the solver's
 * default, as `build_x0` in the reference does through `solver_fluid_initsol`.
 */
template <class T>
FluidSymSystem fluid_symodes(const qn::NetworkStruct<T>& sn, const std::string& method_in,
                             double pstar, const std::vector<double>& init_sol) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::string method = method_in;
    if (method.compare(0, 6, "fluid.") == 0) method = method.substr(6);

    FluidSymSystem sys;
    if (method == "default" || method == "matrix" || method == "pnorm") {
        sys.form = "W";
        sys.method = (method == "default") ? "matrix" : method;
    } else if (method == "closing" || method == "statedep" || method == "softmin") {
        sys.form = "J";
        sys.method = method;
    } else {
        throw UnsupportedError(
            "fluid_symodes: the symbolic ODE export is unsupported for method '" + method_in +
            "'; supported are default, matrix, pnorm, closing, statedep and softmin");
    }

    sys.station_names.resize(M);
    sys.sched_names.resize(M);
    sys.sched.resize(M);
    for (std::size_t i = 0; i < M; ++i) {
        sys.station_names[i] = sn.stations[i].name;
        sys.sched[i] = sn.stations[i].sched;
        sys.sched_names[i] = lang::sched_to_text(sn.stations[i].sched);
    }
    sys.class_names.resize(K);
    for (std::size_t r = 0; r < K; ++r) sys.class_names[r] = sn.classes[r].name;

    double closed_pop = 0.0;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isfinite(sn.classes[r].population)) closed_pop += sn.classes[r].population;
    sys.S.assign(M, 1.0);
    for (std::size_t i = 0; i < M; ++i) {
        const double c = sn.stations[i].nservers;
        sys.S[i] = std::isfinite(c) ? c : closed_pop;
    }

    const FluidLayout L = fluid_layout(sn);

    if (sys.form == "W") {
        // Taken from the assembly the integrator itself uses, so the exported
        // system is the system that runs and not a second opinion about it.
        const FluidMatrixSystem ms = fluid_matrix_system(sn, init_sol, pstar);
        sys.nstates = ms.nstates;
        sys.W = ms.W;
        sys.alambda = ms.alambda;
        sys.is_source = ms.is_source;
        sys.is_inf = ms.is_inf;
        sys.x0 = ms.x0;
        sys.smoothing = (pstar > 0.0) ? "pnorm" : "min";
        sys.pstar.assign(M, pstar);
        // Walk the same block enumeration to recover (station, class, phase).
        sys.state_station.reserve(ms.nstates);
        sys.state_class.reserve(ms.nstates);
        sys.state_phase.reserve(ms.nstates);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                const std::size_t p = L.kic[i][r];
                if (p == 0) continue;  // the disabled placeholder is dropped
                for (std::size_t k = 0; k < p; ++k) {
                    sys.state_station.push_back(i + 1);
                    sys.state_class.push_back(r + 1);
                    sys.state_phase.push_back(k + 1);
                }
            }
        if (sys.state_station.size() != sys.nstates)
            throw NumericError("fluid_symodes: the W-form state metadata does not match the drift");
        return sys;
    }

    // ---- form J -----------------------------------------------------------
    const bool sd = (method == "statedep" || method == "softmin");
    if (sd)
        for (std::size_t i = 0; i < M; ++i)
            if (sn.stations[i].sched == lang::SchedStrategy::EXT)
                throw UnsupportedError("fluid_symodes: the '" + method +
                                       "' method does not support open models, so their ODE "
                                       "system cannot be exported; use 'matrix' or 'closing'");

    sys.nstates = L.nstates;
    sys.state_station.assign(L.nstates, 0);
    sys.state_class.assign(L.nstates, 0);
    sys.state_phase.assign(L.nstates, 0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c)
            for (std::size_t k = 0; k < L.kic[i][c]; ++k) {
                sys.state_station[L.qidx[i][c] + k] = i + 1;
                sys.state_class[L.qidx[i][c] + k] = c + 1;
                sys.state_phase[L.qidx[i][c] + k] = k + 1;
            }

    sys.dpsw = Matrix<double>(M, K, 1.0);
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched != lang::SchedStrategy::DPS) continue;
        double tot = 0.0;
        for (std::size_t r = 0; r < K; ++r) {
            const double w = (r < sn.stations[i].schedparam.size())
                                 ? num_traits<T>::to_double(sn.stations[i].schedparam[r])
                                 : 1.0;
            sys.dpsw(i, r) = w;
            tot += w;
        }
        if (tot > 0.0)
            for (std::size_t r = 0; r < K; ++r) sys.dpsw(i, r) /= tot;
    }

    sys.fcfs_phase_w.assign(L.nstates, 0.0);
    if (sd)
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.stations[i].sched != lang::SchedStrategy::FCFS) continue;
            for (std::size_t c = 0; c < K; ++c) {
                if (!L.enabled[i][c]) continue;
                for (std::size_t k = 0; k < L.kic[i][c]; ++k)
                    sys.fcfs_phase_w[L.qidx[i][c] + k] =
                        -1.0 / num_traits<T>::to_double(sn.service[i][c].D0(k, k));
            }
        }

    // The service processes and the station-space routing, as the drift reads them.
    std::vector<std::vector<std::vector<double>>> mu(M, std::vector<std::vector<double>>(K));
    std::vector<std::vector<std::vector<double>>> phi(M, std::vector<std::vector<double>>(K));
    std::vector<std::vector<std::vector<double>>> pie(M, std::vector<std::vector<double>>(K));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) {
                pie[i][r] = std::vector<double>{1.0};
                continue;
            }
            detail::fluid_mu_phi(sn.service[i][r], mu[i][r], phi[i][r]);
            pie[i][r] = detail::fluid_pie(sn.service[i][r]);
        }
    const std::size_t NS = sn.nof_stateful();
    const bool have_rt = sn.rt.rows() == NS * K;
    std::vector<std::size_t> sf(M, 0);
    for (std::size_t i = 0; i < M; ++i) sf[i] = sn.stateful_of_station(i + 1) - 1;
    const auto route = [&](std::size_t i, std::size_t c, std::size_t j, std::size_t l) -> double {
        if (!have_rt) return 0.0;
        return num_traits<T>::to_double(sn.rt(sf[i] * K + c, sf[j] * K + l));
    };

    std::vector<std::vector<double>> cols;  // the J columns, assembled then transposed
    // ---- departures -------------------------------------------------------
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            if (!L.enabled[i][c]) continue;
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t l = 0; l < K; ++l) {
                    if (!(route(i, c, j, l) > 0.0)) continue;
                    for (std::size_t ki = 0; ki < L.kic[i][c]; ++ki)
                        for (std::size_t kj = 0; kj < L.kic[j][l]; ++kj) {
                            if (sd) {
                                // An INF self-loop has no event, and a station
                                // whose policy has no branch contributes none.
                                if (sn.stations[i].sched == lang::SchedStrategy::INF && j == i)
                                    continue;
                                if (!detail::sym_handled_sched(sn.stations[i].sched)) continue;
                            }
                            const double pj = kj < pie[j][l].size() ? pie[j][l][kj] : 0.0;
                            double base = phi[i][c][ki] * mu[i][c][ki] * route(i, c, j, l) * pj;
                            if (!(base > 0.0)) continue;
                            std::vector<double> col(L.nstates, 0.0);
                            col[L.qidx[i][c] + ki] -= 1.0;
                            col[L.qidx[j][l] + kj] += 1.0;
                            const std::size_t var = L.qidx[i][c] + ki;
                            const SymFactor f =
                                detail::sym_factor(method, sn.stations[i].sched, i, c, ki, base, L,
                                                   sys.S, sys.dpsw, sys.fcfs_phase_w, var);
                            cols.push_back(col);
                            sys.coeff.push_back(base);
                            sys.event_var.push_back(var);
                            sys.factor.push_back(f);
                        }
                }
        }
    // ---- phase changes ----------------------------------------------------
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            if (!L.enabled[i][c]) continue;
            if (sd && !detail::sym_handled_sched(sn.stations[i].sched)) continue;
            for (std::size_t ki = 0; ki + 1 < L.kic[i][c]; ++ki)
                for (std::size_t kp = 0; kp < L.kic[i][c]; ++kp) {
                    if (kp == ki) continue;
                    double base = num_traits<T>::to_double(sn.service[i][c].D0(ki, kp));
                    if (!(base > 0.0)) continue;
                    std::vector<double> col(L.nstates, 0.0);
                    col[L.qidx[i][c] + ki] = -1.0;
                    col[L.qidx[i][c] + kp] = 1.0;
                    const std::size_t var = L.qidx[i][c] + ki;
                    const SymFactor f =
                        detail::sym_factor(method, sn.stations[i].sched, i, c, ki, base, L, sys.S,
                                           sys.dpsw, sys.fcfs_phase_w, var);
                    cols.push_back(col);
                    sys.coeff.push_back(base);
                    sys.event_var.push_back(var);
                    sys.factor.push_back(f);
                }
        }

    sys.nevents = cols.size();
    sys.J = Matrix<double>(L.nstates, sys.nevents, 0.0);
    for (std::size_t e = 0; e < sys.nevents; ++e)
        for (std::size_t s = 0; s < L.nstates; ++s) sys.J(s, e) = cols[e][s];
    if (method == "softmin") sys.alpha = 20.0;
    sys.x0 = init_sol;  // the caller supplies the default; see fluid_export_odes
    return sys;
}

/** The drift as expression strings, one per state variable, plus their names. */
struct FluidSymbolicDrift {
    std::vector<std::string> rhs;
    std::vector<std::string> vars;
};

namespace detail {

/** A 1-based index as text, for the `x1 ... xn` variable names. */
inline std::string sym_state_index(std::size_t v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%llu", static_cast<unsigned long long>(v));
    return std::string(buf);
}

/** `num2char`: decimal text a computer-algebra backend reads as exact. */
inline std::string sym_plain_num(double v) {
    char buf[64];
    if (v == std::floor(v) && std::fabs(v) < 1e15) {
        std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>(v));
        return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

/** `joinSum`: the sum of PARTS, with OFFSET added only when it is not zero. */
inline std::string sym_join_sum(const std::vector<std::string>& parts, const std::string& offset) {
    if (parts.empty()) return "(" + offset + ")";
    std::string body;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i) body += " + ";
        body += parts[i];
    }
    return (offset == "0") ? "(" + body + ")" : "(" + offset + " + " + body + ")";
}

/** `stationSum`: the total fluid mass at station i, plus its consumer's offset. */
inline std::string sym_station_sum(const FluidSymSystem& sys, std::size_t station1,
                                  const std::vector<std::string>& vars, const std::string& offset) {
    std::vector<std::string> parts;
    for (std::size_t s = 0; s < sys.nstates; ++s)
        if (sys.state_station[s] == station1) parts.push_back(vars[s]);
    return sym_join_sum(parts, offset);
}


/** `phaseWeightedStationSum`: sum_u w_u x_u with w_u the mean phase time. */
inline std::string sym_phase_weighted_station_sum(const FluidSymSystem& sys, std::size_t station1,
                                                  const std::vector<std::string>& vars,
                                                  const std::string& offset) {
    std::vector<std::string> parts;
    for (std::size_t s = 0; s < sys.nstates; ++s) {
        if (sys.state_station[s] != station1) continue;
        const double w = sys.fcfs_phase_w[s];
        if (w == 0.0) continue;
        parts.push_back("(" + sym_plain_num(w) + ")*" + vars[s]);
    }
    return sym_join_sum(parts, offset);
}

/**
 * `softminExpr`: the smooth minimum in its weighted-average form,
 * (x e^{-a x} + y e^{-a y}) / (e^{-a x} + e^{-a y}).
 *
 * `softmin` itself rewrites this as lo + gap*w/(1+w) only to keep the exponent
 * argument non-positive, an overflow guard that is meaningless symbolically and
 * would introduce the min/max branch this export exists to avoid.
 */
inline std::string sym_softmin_expr(const std::string& x, const std::string& y, double alpha) {
    const std::string a = sym_plain_num(alpha);
    return "((" + x + ")*exp(-(" + a + ")*(" + x + ")) + (" + y + ")*exp(-(" + a + ")*(" + y +
           ")))/(exp(-(" + a + ")*(" + x + ")) + exp(-(" + a + ")*(" + y + ")))";
}

}  // namespace detail

/**
 * Port of `@@SolverFLD/getSymbolicDrift`: the right-hand side of the mean-field ODE
 * system as expression strings, one per state variable.
 *
 * This is the input a computer-algebra backend needs to produce a Jacobian or an
 * equilibrium, and it is the same system `fluid_symodes` describes and
 * `fluid_export_odes` typesets, written out variable by variable instead of in
 * matrix form.
 *
 * ONLY SMOOTH DRIFTS ARE EXPORTED. The default, matrix, closing and statedep
 * methods scale rates by min(n_i, S_i), which is not differentiable at n_i = S_i,
 * so their Jacobian does not exist there; emitting a one-sided derivative would be
 * a silent lie exactly at the regime switch that matters. The p-norm smoothing
 * (`pstar`, methods matrix or pnorm) and the `softmin` method are smooth
 * everywhere, and everything else is refused BY THE FACTOR TYPE that carries the
 * kink rather than by method name, so a new non-smooth branch cannot slip through.
 */
inline FluidSymbolicDrift fluid_symbolic_drift(const FluidSymSystem& sys) {
    FluidSymbolicDrift out;
    out.vars.resize(sys.nstates);
    for (std::size_t s = 0; s < sys.nstates; ++s)
        out.vars[s] = "x" + detail::sym_state_index(s + 1);
    out.rhs.assign(sys.nstates, "0");
    // `eps0`: the FineTol the softmin denominator carries, as the reference writes
    // it into the expression rather than folding it into a number.
    const std::string eps0 = detail::sym_plain_num(lang::GlobalConstants::FineTol);

    if (sys.form == "W") {
        // dx/dt = W' theta(x) + Alambda, with the p-norm smoothing:
        //   theta_s = x_s / (1 + (n_i/S_i)^p_i)^(1/p_i),  theta_s = 0 at a Source.
        if (sys.smoothing != "pnorm")
            throw UnsupportedError(
                "fluid_symbolic_drift: the drift of this method scales rates by min(n_i, S_i), "
                "which is not differentiable at n_i = S_i, so it has no Jacobian there. Set pstar "
                "to use the p-norm smoothing, or use the 'softmin' method");
        std::vector<std::string> theta(sys.nstates, "0");
        for (std::size_t s = 0; s < sys.nstates; ++s) {
            if (sys.is_source[s]) continue;
            const std::size_t i = sys.state_station[s];
            const double S = sys.S[i - 1];
            const double p = sys.pstar[i - 1];
            // An INF station has no min() to smooth, so its theta stays linear;
            // S holds the population there, not infinity
            if (s < sys.is_inf.size() && sys.is_inf[s]) {
                theta[s] = out.vars[s];
                continue;
            }
            if (S <= 0.0 || p <= 0.0) {
                theta[s] = out.vars[s];
                continue;
            }
            const std::string ni = detail::sym_station_sum(sys, i, out.vars, eps0);
            theta[s] = out.vars[s] + "/(1 + (" + ni + "/" + detail::sym_plain_num(S) + ")^" +
                       detail::sym_plain_num(p) + ")^(1/" + detail::sym_plain_num(p) + ")";
        }
        for (std::size_t s = 0; s < sys.nstates; ++s) {
            std::vector<std::string> terms;
            for (std::size_t t = 0; t < sys.nstates; ++t) {
                const double w = sys.W(t, s);
                if (w == 0.0 || theta[t] == "0") continue;
                terms.push_back("(" + detail::sym_plain_num(w) + ")*(" + theta[t] + ")");
            }
            if (sys.alambda[s] != 0.0) terms.push_back(detail::sym_plain_num(sys.alambda[s]));
            if (terms.empty()) continue;
            std::string body;
            for (std::size_t i = 0; i < terms.size(); ++i) {
                if (i) body += " + ";
                body += terms[i];
            }
            out.rhs[s] = body;
        }
        return out;
    }
    if (sys.form != "J")
        throw UnsupportedError("fluid_symbolic_drift: unsupported ODE form '" + sys.form + "'");

    // dx/dt = J r(x), with r_e = coeff(e) * factor_e(x). Only the smooth factor
    // types are exportable: 'min' (PS/FCFS under closing and statedep), 'fcfsw'
    // (statedep FCFS), 'dpsmin' (closing DPS) and 'dpspw' (piecewise DPS) all
    // carry a min or a branch.
    std::vector<std::string> rate(sys.nevents, "0");
    for (std::size_t e = 0; e < sys.nevents; ++e) {
        const SymFactor& f = sys.factor[e];
        const std::string& v = out.vars[sys.event_var[e]];
        std::string factor;
        if (f.type == "lin") {
            factor = v;
        } else if (f.type == "ext1") {
            if (f.others.empty()) {
                factor = "1";
            } else {
                std::string body;
                for (std::size_t k = 0; k < f.others.size(); ++k) {
                    if (k) body += " + ";
                    body += out.vars[f.others[k] - 1];
                }
                factor = "(1 - (" + body + "))";
            }
        } else if (f.type == "fcfsws") {
            const std::string ni = detail::sym_station_sum(sys, f.station, out.vars, "0");
            const std::string nhat =
                detail::sym_phase_weighted_station_sum(sys, f.station, out.vars, eps0);
            factor = v + "*(" +
                     detail::sym_softmin_expr(ni, detail::sym_plain_num(sys.S[f.station - 1]),
                                              sys.alpha) +
                     ")/(" + nhat + ")";
        } else {
            throw UnsupportedError(
                "fluid_symbolic_drift: event " + detail::sym_state_index(e + 1) +
                " scales its rate by the non-smooth factor '" + f.type +
                "', which has no derivative where the regime switches, so the system has no "
                "Jacobian. Use the 'softmin' method, or the p-norm smoothing of the 'matrix' "
                "method");
        }
        rate[e] = "(" + detail::sym_plain_num(sys.coeff[e]) + ")*(" + factor + ")";
    }
    for (std::size_t s = 0; s < sys.nstates; ++s) {
        std::vector<std::string> terms;
        for (std::size_t e = 0; e < sys.nevents; ++e) {
            const double j = sys.J(s, e);
            if (j == 0.0) continue;
            terms.push_back("(" + detail::sym_plain_num(j) + ")*(" + rate[e] + ")");
        }
        if (terms.empty()) continue;
        std::string body;
        for (std::size_t i = 0; i < terms.size(); ++i) {
            if (i) body += " + ";
            body += terms[i];
        }
        out.rhs[s] = body;
    }
    return out;
}

/** What `@@SolverFLD/getJacobian` returns: d f_i / d x_j as expression strings. */
struct FluidSymbolicJacobian {
    std::vector<std::string> vars;                ///< the state variable names
    std::vector<std::string> rhs;                 ///< the drift, one per variable
    std::vector<std::vector<std::string> > J;     ///< J[i][j] = d f_i / d x_j
};

namespace detail {

/** d(n_i)/d(x_m): one when state m belongs to station i, zero otherwise. */
inline bool sym_in_station(const FluidSymSystem& sys, std::size_t station1, std::size_t m) {
    return sys.state_station[m] == station1;
}

/**
 * d/dx_m of the softmin of (n_i, S), by the chain rule through n_i.
 *
 * With a = n_i, b = S constant and D = e^{-alpha a} + e^{-alpha b},
 *
 *   sm(a) = (a e^{-alpha a} + b e^{-alpha b}) / D
 *   sm'(a) = e^{-alpha a} (1 - alpha a + alpha sm(a)) / D
 *
 * which is written out rather than left as a quotient of four exponentials
 * because the form above has one exponential per state and the naive one has
 * four, and a backend that simplifies neither would carry the difference into
 * every entry of the Jacobian.
 */
inline std::string sym_softmin_deriv(const std::string& a, const std::string& b, double alpha) {
    const std::string al = sym_plain_num(alpha);
    const std::string ea = "exp(-(" + al + ")*(" + a + "))";
    const std::string eb = "exp(-(" + al + ")*(" + b + "))";
    const std::string D = "(" + ea + " + " + eb + ")";
    const std::string sm = sym_softmin_expr(a, b, alpha);
    return "(" + ea + "*(1 - (" + al + ")*(" + a + ") + (" + al + ")*(" + sm + ")))/" + D;
}

/** `(u)*(v)`, dropping the term when either side is the literal zero. */
inline std::string sym_mul(const std::string& u, const std::string& v) {
    if (u == "0" || v == "0") return "0";
    if (u == "1") return v;
    if (v == "1") return u;
    return "(" + u + ")*(" + v + ")";
}

/** Sum of the non-zero terms, or the literal zero. */
inline std::string sym_sum(const std::vector<std::string>& parts) {
    std::string body;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (parts[i] == "0") continue;
        if (!body.empty()) body += " + ";
        body += parts[i];
    }
    return body.empty() ? std::string("0") : body;
}

}  // namespace detail

/**
 * Port of `@@SolverFLD/getJacobian`: d f_i / d x_j of the mean-field drift, as
 * expression strings.
 *
 * IT IS DIFFERENTIATED HERE, NOT SENT AWAY. The reference hands the drift
 * strings to the line-sage-rest backend over HTTP and returns what SAGE gives
 * back. This port has no symbolic engine and no HTTP client -- the CTMC port
 * makes the same refusal for its own symbolic getters -- but the drift is not an
 * opaque expression here: `FluidSymSystem` carries it STRUCTURALLY, as a jump
 * matrix times per-event rates whose state dependence is one of a closed set of
 * typed factors. Each of those factors has a derivative that can be written down
 * once, so the Jacobian is produced exactly, by the chain rule over the
 * structure, rather than by parsing the strings the drift printer emitted.
 *
 * The smoothness gate is the SAME one `fluid_symbolic_drift` applies, and for
 * the same reason: `min(n_i, S_i)` has no derivative at n_i = S_i, and a
 * one-sided value there would be a silent lie exactly at the regime switch. The
 * refusal is by FACTOR TYPE, so a new non-smooth branch cannot slip through.
 *
 * `equilibria`, the reference's fourth output, is NOT produced here: solving
 * f(x) = 0 in closed form needs a computer-algebra solver, which differentiating
 * does not. It lives in `fluid_jacobian` (fluid_jacobian.h), which resolves the
 * same line-sage-rest backend the reference uses and asks it; this function stays
 * backend-free so that the Jacobian alone never needs one.
 */
inline FluidSymbolicJacobian fluid_symbolic_jacobian(const FluidSymSystem& sys) {
    const FluidSymbolicDrift d = fluid_symbolic_drift(sys);
    FluidSymbolicJacobian out;
    out.vars = d.vars;
    out.rhs = d.rhs;
    const std::size_t n = sys.nstates;
    out.J.assign(n, std::vector<std::string>(n, "0"));
    const std::string eps0 = detail::sym_plain_num(lang::GlobalConstants::FineTol);

    if (sys.form == "W") {
        // theta_s = x_s / Q_i,  Q_i = (1 + (n_i/S_i)^p_i)^(1/p_i)
        //   d theta_s / d x_m = [m == s]/Q_i - x_s Q_i' / Q_i^2
        //   Q_i' = (1 + (n_i/S_i)^p)^(1/p - 1) (n_i/S_i)^(p-1) / S_i, for m at i
        std::vector<std::vector<std::string> > dtheta(n, std::vector<std::string>(n, "0"));
        for (std::size_t s = 0; s < n; ++s) {
            if (sys.is_source[s]) continue;
            const std::size_t i = sys.state_station[s];
            const double S = sys.S[i - 1], p = sys.pstar[i - 1];
            if ((s < sys.is_inf.size() && sys.is_inf[s]) || S <= 0.0 || p <= 0.0) {
                dtheta[s][s] = "1";  // theta_s = x_s
                continue;
            }
            const std::string ni = detail::sym_station_sum(sys, i, out.vars, eps0);
            const std::string Sn = detail::sym_plain_num(S), pn = detail::sym_plain_num(p);
            const std::string base = "(1 + ((" + ni + ")/" + Sn + ")^" + pn + ")";
            const std::string Q = base + "^(1/" + pn + ")";
            const std::string dQ = base + "^(1/" + pn + " - 1)*((" + ni + ")/" + Sn + ")^(" + pn +
                                   " - 1)/" + Sn;
            for (std::size_t m = 0; m < n; ++m) {
                std::vector<std::string> parts;
                if (m == s) parts.push_back("1/(" + Q + ")");
                if (detail::sym_in_station(sys, i, m))
                    parts.push_back("-(" + out.vars[s] + ")*(" + dQ + ")/(" + Q + ")^2");
                dtheta[s][m] = detail::sym_sum(parts);
            }
        }
        for (std::size_t s = 0; s < n; ++s)
            for (std::size_t m = 0; m < n; ++m) {
                std::vector<std::string> parts;
                for (std::size_t t = 0; t < n; ++t) {
                    const double w = sys.W(t, s);
                    if (w == 0.0 || dtheta[t][m] == "0") continue;
                    parts.push_back(detail::sym_mul(detail::sym_plain_num(w), dtheta[t][m]));
                }
                out.J[s][m] = detail::sym_sum(parts);
            }
        return out;
    }
    if (sys.form != "J")
        throw UnsupportedError("fluid_symbolic_jacobian: unsupported ODE form '" + sys.form + "'");

    // d rate_e / d x_m = coeff_e * d factor_e / d x_m, one row per factor type.
    std::vector<std::vector<std::string> > drate(sys.nevents, std::vector<std::string>(n, "0"));
    for (std::size_t e = 0; e < sys.nevents; ++e) {
        const SymFactor& f = sys.factor[e];
        const std::size_t v = sys.event_var[e];
        const std::string cf = detail::sym_plain_num(sys.coeff[e]);
        if (f.type == "lin") {
            drate[e][v] = cf;
        } else if (f.type == "ext1") {
            for (std::size_t k = 0; k < f.others.size(); ++k)
                drate[e][f.others[k] - 1] = "-(" + cf + ")";
        } else if (f.type == "fcfsws") {
            // factor = x_v * softmin(n_i, S) / nhat_i, with n_i and nhat_i both
            // linear in the states of station i. Product and quotient rules:
            //   d/dx_m = [m == v] sm/nhat
            //          + x_v sm'(n_i) [m at i] / nhat
            //          - x_v sm w_m [m at i] / nhat^2
            const std::string ni = detail::sym_station_sum(sys, f.station, out.vars, "0");
            const std::string nhat =
                detail::sym_phase_weighted_station_sum(sys, f.station, out.vars, eps0);
            const std::string Sn = detail::sym_plain_num(sys.S[f.station - 1]);
            const std::string sm = detail::sym_softmin_expr(ni, Sn, sys.alpha);
            const std::string dsm = detail::sym_softmin_deriv(ni, Sn, sys.alpha);
            for (std::size_t m = 0; m < n; ++m) {
                std::vector<std::string> parts;
                if (m == v) parts.push_back("(" + sm + ")/(" + nhat + ")");
                if (detail::sym_in_station(sys, f.station, m)) {
                    parts.push_back("(" + out.vars[v] + ")*(" + dsm + ")/(" + nhat + ")");
                    const double w = sys.fcfs_phase_w[m];
                    if (w != 0.0)
                        parts.push_back("-(" + out.vars[v] + ")*(" + sm + ")*(" +
                                        detail::sym_plain_num(w) + ")/(" + nhat + ")^2");
                }
                const std::string body = detail::sym_sum(parts);
                drate[e][m] = body == "0" ? std::string("0") : detail::sym_mul(cf, body);
            }
        } else {
            throw UnsupportedError(
                "fluid_symbolic_jacobian: event " + detail::sym_state_index(e + 1) +
                " scales its rate by the non-smooth factor '" + f.type +
                "', which has no derivative where the regime switches, so the system has no "
                "Jacobian. Use the 'softmin' method, or the p-norm smoothing of the 'matrix' "
                "method");
        }
    }

    for (std::size_t s = 0; s < n; ++s)
        for (std::size_t m = 0; m < n; ++m) {
            std::vector<std::string> parts;
            for (std::size_t e = 0; e < sys.nevents; ++e) {
                const double j = sys.J(s, e);
                if (j == 0.0 || drate[e][m] == "0") continue;
                parts.push_back(detail::sym_mul(detail::sym_plain_num(j), drate[e][m]));
            }
            out.J[s][m] = detail::sym_sum(parts);
        }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_SYMODES_H
