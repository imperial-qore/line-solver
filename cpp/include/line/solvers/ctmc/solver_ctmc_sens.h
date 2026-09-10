/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `@@SolverCTMC/getSensitivity` and `getSensitivityRanking`: the
 * parametric sensitivity of a steady-state reward to a scalar model parameter,
 * following Trivedi and Bobbio (2017), Sec. 9.7.
 *
 * WHERE THE ERROR ACTUALLY COMES FROM, which is what the reference's two
 * methods differ about. Differentiating pi Q = 0 gives
 * (dpi/dtheta) Q = -pi (dQ/dtheta), one extra solve with the matrix the
 * stationary solve already factored (`ctmc_sens`, which is EXACT). So the only
 * approximation is dQ/dtheta, and the reference obtains it by central
 * differences on the RATES with the state space held fixed -- legitimate
 * because the state space depends on the topology and the cutoff, never on a
 * rate value. That is O(h^2) on the generator alone, not on the solve.
 *
 * The reference's 'symbolic' method removes even that by solving pi as a
 * rational function of the rate symbols in a computer-algebra backend and
 * differentiating it exactly, leaving only the RATE MAP x_e(theta) differenced.
 * That map is affine in theta in the common cases -- a rate set to theta, or
 * scaled by it -- and a central difference is exact on an affine map, so the
 * whole O(h^2) error of 'fd' disappears. It is served here through `api/sym`,
 * the same line-sage-rest service the reference talks to.
 *
 * WHAT THE SYMBOLIC METHOD REFUSES RATHER THAN APPROXIMATES: a theta that
 * RESHAPES an event's filtration instead of merely scaling it. The chain rule
 * above assumes each event contributes one rate, so a parameter that changes
 * the relative weights inside one filtration breaks the premise; the shape is
 * compared before any round trip and the request is refused by name.
 *
 * ANY PERTURBATION THAT RESIZES THE STATE SPACE IS AN ERROR, not something to
 * paper over: it means theta switched a transition on or off (a rate crossing
 * zero, or an immediate transition appearing), so the two generators describe
 * different chains and their difference is meaningless.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_SENS_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_SENS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include <map>
#include <memory>

#include "line/api/mc/ctmc_sens.h"
#include "line/api/sym/sym_engine.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_symbolic.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/** The scalar parameter a sensitivity is taken with respect to. */
template <class T>
struct CtmcSensParam {
    std::string name;
    double value = 0.0;
    /** Apply theta to a COPY of the struct; the original is never mutated. */
    std::function<void(NetworkStruct<T>&, double)> set;
    /** Central-difference step; <= 0 takes the reference's max(|theta|,1)*1e-6. */
    double step = 0.0;
};

/** What one sensitivity computation returns. */
template <class T>
struct CtmcSens {
    T S = num_traits<T>::from_int(0);   ///< d E[r] / d theta, Eq. (9.79)
    T SS = num_traits<T>::from_int(0);  ///< (theta/E[r]) d E[r] / d theta, Eq. (9.80)
    bool scaled_valid = false;          ///< false when E[r] is zero, MATLAB's NaN
    std::vector<T> dpi, pi;
};

namespace sens_detail {

/**
 * Minimum positive rate of each event filtration, and the filtration normalized
 * by it -- the reference's `eventRates`. An event with no positive rate is
 * inactive and contributes neither.
 */
template <class T>
struct EventRates {
    std::vector<double> rate;
    std::vector<Matrix<T> > shape;
    std::vector<bool> active;
};

template <class T>
EventRates<T> event_rates(const std::vector<Matrix<T> >& F) {
    EventRates<T> out;
    out.rate.assign(F.size(), 0.0);
    out.shape.resize(F.size());
    out.active.assign(F.size(), false);
    for (std::size_t e = 0; e < F.size(); ++e) {
        const symbolic_detail::MinPositive<T> m = symbolic_detail::min_positive(F[e]);
        if (!m.has) continue;
        out.rate[e] = num_traits<T>::to_double(m.value);
        Matrix<T> S(F[e].rows(), F[e].cols(), num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < S.rows(); ++i)
            for (std::size_t j = 0; j < S.cols(); ++j) S(i, j) = T(F[e](i, j) / m.value);
        out.shape[e] = S;
        out.active[e] = true;
    }
    return out;
}

/** True when two normalized filtrations agree entrywise to `tol`. */
template <class T>
bool same_shape(const Matrix<T>& a, const Matrix<T>& b, double tol) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j)
            if (std::fabs(num_traits<T>::to_double(a(i, j)) - num_traits<T>::to_double(b(i, j))) >
                tol)
                return false;
    return true;
}

/**
 * The reference's `symbolicSensitivity`: exact d(pi)/d(x_e) from the backend,
 * combined with a differenced rate map d(x_e)/d(theta) by the chain rule.
 *
 * TWO STEP SIZES ARE IN PLAY AND THEY ARE NOT INTERCHANGEABLE. The rate map is
 * probed at the coarse `1e-3` step because it is expected to be affine, where a
 * coarse step is exact and immune to cancellation. If the second difference says
 * it is NOT affine -- the curvature test below -- the caller's fine step is used
 * instead, trading that exactness for the smaller truncation error a curved map
 * needs.
 */
template <class T>
void symbolic_sensitivity(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                          const CtmcSensParam<T>& param, double theta, double h,
                          const CtmcSymbolicOptions& symopt, std::vector<T>& pi_out,
                          std::vector<T>& dpi_out) {
    const CtmcSymbolicGenerator<T> g = ctmc_symbolic_generator(sn, opt);
    const std::size_t n = g.space.size();
    const std::size_t ne = g.symbols.size();

    std::vector<double> rate0(ne, 0.0);
    for (std::size_t e = 0; e < ne; ++e) rate0[e] = num_traits<T>::to_double(g.rate0[e]);

    const double hrate = std::max(std::fabs(theta), 1.0) * 1e-3;
    NetworkStruct<T> up = sn, dn = sn;
    param.set(up, theta + hrate);
    param.set(dn, theta - hrate);
    const CtmcGenerator<T> gu = ctmc_get_generator(up, opt);
    const CtmcGenerator<T> gd = ctmc_get_generator(dn, opt);
    if (gu.space.size() != n || gd.space.size() != n)
        throw InputError(
            "solver_ctmc_sensitivity: perturbing the parameter changed the state-space size, so "
            "the generators cannot be differenced; this happens when theta switches a transition "
            "on or off (a zero rate, or an immediate transition appearing)");
    if (gu.filt.size() != ne || gd.filt.size() != ne)
        throw InputError(
            "solver_ctmc_sensitivity: perturbing the parameter changed the number of events");

    EventRates<T> ru = event_rates(gu.filt), rd = event_rates(gd.filt);
    for (std::size_t e = 0; e < ne; ++e) {
        if (g.symbols[e].empty()) continue;
        if (!ru.active[e] || !rd.active[e] || !same_shape(ru.shape[e], g.filt[e], 1e-8) ||
            !same_shape(rd.shape[e], g.filt[e], 1e-8))
            throw UnsupportedError(
                "solver_ctmc_sensitivity: perturbing the parameter reshapes the filtration of "
                "event " +
                std::to_string(e + 1) +
                " rather than scaling it, so the generator is not linear in a single rate per "
                "event and the symbolic chain rule does not apply; use method 'fd'");
    }

    double curvature = 0.0, scale = 1.0;
    for (std::size_t e = 0; e < ne; ++e) {
        curvature = std::max(curvature, std::fabs(ru.rate[e] + rd.rate[e] - 2.0 * rate0[e]));
        scale = std::max(scale, std::fabs(rate0[e]));
    }
    double step = hrate;
    if (curvature > 1e-9 * scale) {
        NetworkStruct<T> up2 = sn, dn2 = sn;
        param.set(up2, theta + h);
        param.set(dn2, theta - h);
        ru = event_rates(ctmc_get_generator(up2, opt).filt);
        rd = event_rates(ctmc_get_generator(dn2, opt).filt);
        step = h;
    }
    std::vector<double> drate(ne, 0.0);
    for (std::size_t e = 0; e < ne; ++e) drate[e] = (ru.rate[e] - rd.rate[e]) / (2.0 * step);

    const std::shared_ptr<sym::SymEngine> engine =
        symbolic_detail::require_engine(symopt, "the 'symbolic' sensitivity");
    const std::vector<std::string> pi_expr = engine->solveCTMC(g.Q, g.active_symbols()).pi;
    if (pi_expr.size() != n)
        throw sym::SymEngineError("solver_ctmc_sensitivity: backend '" + engine->name() +
                                  "' returned " + std::to_string(pi_expr.size()) +
                                  " entries for a " + std::to_string(n) + " state chain");

    std::map<std::string, double> assignment;
    for (std::size_t e = 0; e < ne; ++e)
        if (!g.symbols[e].empty()) assignment[g.symbols[e]] = rate0[e];

    const std::vector<double> pi = engine->eval(pi_expr, assignment);
    pi_out.assign(n, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < n && s < pi.size(); ++s)
        pi_out[s] = num_traits<T>::from_double(pi[s]);

    dpi_out.assign(n, num_traits<T>::from_int(0));
    for (std::size_t e = 0; e < ne; ++e) {
        if (g.symbols[e].empty()) continue;
        // An event that does not depend on theta contributes a zero term, and
        // its exact derivative is not worth a round trip.
        if (drate[e] == 0.0) continue;
        const std::vector<double> dvals =
            engine->eval(engine->diff(pi_expr, g.symbols[e], 1), assignment);
        if (dvals.size() != n)
            throw sym::SymEngineError("solver_ctmc_sensitivity: backend '" + engine->name() +
                                      "' returned " + std::to_string(dvals.size()) +
                                      " derivative values for a " + std::to_string(n) +
                                      " state chain");
        for (std::size_t s = 0; s < n; ++s)
            dpi_out[s] += num_traits<T>::from_double(drate[e] * dvals[s]);
    }
}

}  // namespace sens_detail

/**
 * Port of `@@SolverCTMC/getSensitivity`.
 *
 * @param reward reward RATE per state; empty returns dpi and pi with S unset
 * @param method `fd` (the default) or `symbolic`, which needs a backend
 * @param sn the refreshed network struct
 * @param opt CTMC options (state-space cutoff, tolerances, method)
 * @param param the parameter theta being perturbed, and how to set it
 * @param symopt backend selection, read only by `method='symbolic'`
 */
template <class T>
CtmcSens<T> solver_ctmc_sensitivity(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                    const CtmcSensParam<T>& param,
                                    const std::vector<T>& reward = std::vector<T>(),
                                    const std::string& method = "fd",
                                    const CtmcSymbolicOptions& symopt = CtmcSymbolicOptions()) {
    if (method != "fd" && method != "symbolic")
        throw InputError("solver_ctmc_sensitivity: unknown method '" + method +
                         "'; expected 'fd' or 'symbolic'");
    if (!param.set)
        throw InputError("solver_ctmc_sensitivity: the parameter carries no setter, so theta "
                         "cannot be applied to the model");

    const double theta = param.value;
    const double h = param.step > 0 ? param.step : std::max(std::fabs(theta), 1.0) * 1e-6;

    CtmcSens<T> out;
    std::size_t n = 0;
    if (method == "symbolic") {
        sens_detail::symbolic_sensitivity(sn, opt, param, theta, h, symopt, out.pi, out.dpi);
        n = out.pi.size();
    } else {
        const CtmcSolution<T> base = solver_ctmc_analyzer(sn, opt);
        n = base.chain.space.size();

        NetworkStruct<T> up = sn, dn = sn;
        param.set(up, theta + h);
        param.set(dn, theta - h);
        const CtmcSolution<T> su = solver_ctmc_analyzer(up, opt);
        const CtmcSolution<T> sd = solver_ctmc_analyzer(dn, opt);
        if (su.chain.space.size() != n || sd.chain.space.size() != n)
            throw InputError(
                "solver_ctmc_sensitivity: perturbing the parameter changed the state-space size, "
                "so the generators cannot be differenced; this happens when theta switches a "
                "transition on or off (a zero rate, or an immediate transition appearing)");

        Matrix<T> dQ(n, n, num_traits<T>::from_int(0));
        const T twoh = num_traits<T>::from_double(2.0 * h);
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b)
                dQ(a, b) = T((su.chain.Q(a, b) - sd.chain.Q(a, b)) / twoh);

        out.pi = base.pi;
        out.dpi = mc::ctmc_sens(base.chain.Q, dQ, base.pi);
    }
    if (reward.empty()) return out;
    if (reward.size() != n)
        throw InputError("solver_ctmc_sensitivity: the reward must have one entry per state");

    // Eq. (9.83) with dr/dtheta = 0: a reward that itself depends on theta needs
    // the second term, which the reference does not carry either.
    T er = num_traits<T>::from_int(0);
    for (std::size_t s = 0; s < n; ++s) {
        out.S += T(out.dpi[s] * reward[s]);
        er += T(out.pi[s] * reward[s]);
    }
    if (std::fabs(num_traits<T>::to_double(er)) > GlobalConstants::Zero) {
        out.SS = T(num_traits<T>::from_double(theta) / er * out.S);
        out.scaled_valid = true;
    }
    return out;
}

/** One row of the ranking table. */
template <class T>
struct CtmcSensRank {
    std::string parameter;
    double value = 0.0;
    T S = num_traits<T>::from_int(0);
    T SS = num_traits<T>::from_int(0);
    bool scaled_valid = false;
};

/**
 * Port of `@@SolverCTMC/getSensitivityRanking`: rank parameters by influence.
 *
 * The order is by DESCENDING ABSOLUTE SCALED sensitivity, because the scaled
 * form is dimensionless and so is the only one comparable across parameters
 * measured in different units. The SIGN is retained in the table, since it says
 * whether increasing a parameter helps or hurts.
 */
template <class T>
std::vector<CtmcSensRank<T>> solver_ctmc_sensitivity_ranking(
    const NetworkStruct<T>& sn, const CtmcOptions& opt,
    const std::vector<CtmcSensParam<T>>& params, const std::vector<T>& reward) {
    if (reward.empty())
        throw InputError("solver_ctmc_sensitivity_ranking: a reward is required to rank "
                         "parameters");
    std::vector<CtmcSensRank<T>> rows;
    for (std::size_t l = 0; l < params.size(); ++l) {
        const CtmcSens<T> r = solver_ctmc_sensitivity(sn, opt, params[l], reward);
        CtmcSensRank<T> row;
        row.parameter = params[l].name.empty() ? "theta" + std::to_string(l + 1) : params[l].name;
        row.value = params[l].value;
        row.S = r.S;
        row.SS = r.SS;
        row.scaled_valid = r.scaled_valid;
        rows.push_back(row);
    }
    // A parameter whose scaled sensitivity is undefined -- E[r] is zero -- sorts
    // LAST rather than being dropped: it was measured, and its unscaled value
    // is still reported.
    std::stable_sort(rows.begin(), rows.end(),
                     [](const CtmcSensRank<T>& a, const CtmcSensRank<T>& b) {
                         if (a.scaled_valid != b.scaled_valid) return a.scaled_valid;
                         if (!a.scaled_valid) return false;
                         return std::fabs(num_traits<T>::to_double(a.SS)) >
                                std::fabs(num_traits<T>::to_double(b.SS));
                     });
    return rows;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_SENS_H
