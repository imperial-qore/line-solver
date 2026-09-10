/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_LSODA_H
#define LINE_UTIL_LSODA_H

/**
 * LSODA: the LINE-facing wrapper over the vendored solver in
 * `third_party/lsoda.hpp`.
 *
 * WHY LSODA AND NOT THE ROSENBROCK IN `ode.h`. Both integrate a stiff system,
 * but the fluid solver has to AGREE WITH THE OTHER CODEBASES, not merely be
 * accurate. MATLAB drives `ode15s`/`ode23t`, the JAR drives
 * `jline.solvers.fluid.LSODAExt` over imperial-qore/lsoda-java, and native
 * Python drives `line_solver/lib/lsoda.py`; the latter two are ports of the
 * same Heng Li `lsoda.c` that the vendored C++ descends from. Running the same
 * algorithm with the same coefficients is what makes the fluid results line up
 * across languages; a different stiff integrator would be defensible
 * numerically and would still disagree in the last digits everywhere.
 * `ode.h` stays the integrator for everything that has no cross-codebase
 * counterpart to match.
 *
 * DOUBLE ONLY, BY CONSTRUCTION. LSODA selects its own order and step from
 * floating-point error estimates and the machine epsilon of `double` is baked
 * into its coefficients, so there is no meaningful `Real<N>` or `Rational`
 * instantiation of it. Callers templated on `T` must refuse by name for any
 * other backend rather than silently narrowing to `double` -- `solver_fluid`
 * does exactly that.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "line/util/error.h"
#include "lsoda.hpp"

namespace line {

/**
 * The right-hand side dy/dt = f(t, y).
 *
 * `y` and `dydt` are plain 0-based arrays of length neq. The vendored solver
 * keeps its state 1-based internally, as the Fortran original did, and hands
 * the callback pointers already offset past the unused slot, so nothing here
 * has to know about that.
 */
using LsodaRhs = std::function<void(double t, const double* y, double* dydt)>;

/**
 * Integration controls.
 *
 * `max_steps` defaults to the JAR's value, not LSODA's own. Upstream defaults
 * to mxstep = 500, which a stiff fluid model exhausts long before the horizon
 * and then reports istate = -1; `LSODAExt` exists in the JAR precisely to raise
 * it, and this default keeps the two in step.
 */
struct LsodaOptions {
    double rtol = 1e-6;      ///< relative tolerance, applied to every component
    double atol = 1e-6;      ///< absolute tolerance, applied to every component
    /**
     * Per-component tolerances. Empty means "use the scalars above"; otherwise
     * the vector must have one entry per equation and overrides its scalar.
     * A stiff model whose components live on different scales needs this --
     * Robertson is the standard example, with atol 1e-10 on the fast species
     * and 1e-6 on the other two.
     */
    std::vector<double> rtol_vec;
    std::vector<double> atol_vec;
    double h_init = 0.0;     ///< initial step; 0 lets LSODA choose
    double h_min = 0.0;      ///< smallest admissible step; 0 means no bound
    double h_max = 0.0;      ///< largest admissible step; 0 means no bound
    std::size_t max_steps = 10000000;  ///< internal steps between output points
    int max_order_nonstiff = 12;       ///< Adams order cap (mxordn)
    int max_order_stiff = 5;           ///< BDF order cap (mxords)
    /**
     * Start on BDF and never switch to Adams.
     *
     * LSODA starts on Adams and switches only when its own detector says the
     * problem is stiff, and that detector needs `pdlast`, the dominant
     * eigenvalue read off the CORRECTOR's convergence rate. AT A FIXED POINT the
     * corrector converges on the `del <= 100*pnorm*ETA` branch before `pdest` is
     * ever formed, so `pdlast` stays 0, `scaleh`'s stability guard never binds,
     * `h` runs up to `h_max`, and a high-order Adams-Moulton outside its
     * stability region WANDERS around the fixed point instead of settling on
     * it -- which is where a fluid solver spends most of its time. The JAR pins
     * it for exactly this reason (`LSODA.setForceStiff(true)`), MATLAB and
     * native Python carry the same flag, and MATLAB's own stiff slot is
     * `ode15s`, a BDF/NDF method, so pinning is also what matches the reference.
     *
     * NOT THE DEFAULT HERE, and the reason is measured: pinning BDF means a
     * finite-difference Jacobian from the first step, and LSODA sizes that
     * increment as `max(sqrt(eps)*|y_j|, r0/ewt_j)`, which collapses to ~1e-19
     * for a component sitting at exactly zero under a tight atol. Robertson
     * started at (1, 0, 0) with atol (1e-6, 1e-10, 1e-6) diverges that way.
     * A fluid drift is not that problem, but a general-purpose default must not
     * carry the hazard.
     */
    bool force_stiff = false;
    /**
     * A test consulted after every ACCEPTED STEP, ending the integration when it
     * returns true and holding the state for the rest of the grid.
     *
     * WHY THE OPTION EXISTS. `lsoda_integrate` reports at output times only, and
     * a window whose drift has reached a fixed point has nothing left to report:
     * f(y*) = 0 means y is that state for every later t, so the remaining span
     * is known. Without this the step controller grinds -- a stiff controller
     * handed a state it is already at cannot pick a step -- and the window never
     * returns. See `solver_fluid.h` and FLUID_FIXED_POINT_GUARD (MATLAB).
     *
     * Empty by default, and the integration then runs exactly as it always has:
     * the itask = 1 loop below is untouched, so no existing caller changes. Set,
     * it selects the stepwise driver, which reproduces the same output grid one
     * accepted step at a time through `LsodaStepper`.
     */
    std::function<bool(double, const std::vector<double>&)> step_stop;
};

/** Result of an integration, mirroring `OdeSolution` in `ode.h`. */
struct LsodaSolution {
    std::vector<double> t;                   ///< output times, t[0] = t_eval[0]
    std::vector<std::vector<double>> y;      ///< y[i] is the state at t[i]
    std::size_t steps = 0;                   ///< internal steps taken (nst)
    std::size_t f_evals = 0;                 ///< right-hand side evaluations (nfe)
    std::size_t jacobians = 0;               ///< Jacobian evaluations (nje)
    bool success = true;                     ///< false when LSODA returned istate < 0
    int istate = 2;                          ///< the solver's final istate
    /** The method in force at the end: "adams" (nonstiff) or "bdf" (stiff). */
    std::string method = "adams";

    const std::vector<double>& final_state() const {
        if (y.empty()) throw NumericError("LsodaSolution: no state was recorded");
        return y.back();
    }
    double final_time() const {
        if (t.empty()) throw NumericError("LsodaSolution: no state was recorded");
        return t.back();
    }
};

namespace lsoda_detail {

/** Trampoline: the vendored callback is a raw function pointer plus a void*. */
inline void rhs_trampoline(double t, double* y, double* dydt, void* data) {
    (*static_cast<const LsodaRhs*>(data))(t, y, dydt);
}

}  // namespace lsoda_detail

/**
 * Integrate dy/dt = f(t, y) from `t_eval.front()` through every later entry of
 * `t_eval`, returning the state at each.
 *
 * `t_eval` must be non-empty and non-decreasing; its first entry is the initial
 * time and is reported back unchanged with `y0`. This is the `tspan` form the
 * fluid analyzers use, so one call covers both "give me the steady state at
 * T" (two entries) and "give me the transient" (many).
 *
 * A solver failure is reported, not thrown: `success` is false and `istate`
 * carries LSODA's own code, because the fluid iteration treats a step that did
 * not converge as a signal to shorten the horizon rather than as an error.
 */
/**
 * The stepwise driver behind `LsodaOptions::step_stop`; defined below, after
 * `LsodaStepper`, which it drives. Declared here so `lsoda_integrate` can hand
 * off to it without moving either function.
 */
inline LsodaSolution lsoda_integrate_stepwise(const LsodaRhs& f, const std::vector<double>& y0,
                                              const std::vector<double>& t_eval,
                                              const LsodaOptions& opt);

inline LsodaSolution lsoda_integrate(const LsodaRhs& f, const std::vector<double>& y0,
                                     const std::vector<double>& t_eval,
                                     const LsodaOptions& opt = LsodaOptions()) {
    if (t_eval.empty())
        throw InputError("lsoda_integrate: t_eval is empty; it must carry at least the start time");
    if (y0.empty())
        throw InputError("lsoda_integrate: the initial state is empty");
    for (std::size_t i = 1; i < t_eval.size(); ++i)
        if (t_eval[i] < t_eval[i - 1])
            throw InputError("lsoda_integrate: t_eval must be non-decreasing");

    // A stop test needs the trajectory one accepted step at a time, which the
    // itask = 1 loop below does not expose; the stepwise driver does, on the
    // same grid. Only that caller pays for it -- everyone else keeps this loop.
    if (opt.step_stop) return lsoda_integrate_stepwise(f, y0, t_eval, opt);

    const std::size_t neq = y0.size();
    LsodaSolution out;
    out.t.push_back(t_eval[0]);
    out.y.push_back(y0);

    lsoda_impl::LSODA solver;
    // Error weights, 1-based with slot 0 unused. itol follows ODEPACK: 1 both
    // scalar, 2 atol per component, 3 rtol per component, 4 both.
    const bool rvec = !opt.rtol_vec.empty(), avec = !opt.atol_vec.empty();
    if (rvec && opt.rtol_vec.size() != neq)
        throw InputError("lsoda_integrate: rtol_vec has one entry per equation or none");
    if (avec && opt.atol_vec.size() != neq)
        throw InputError("lsoda_integrate: atol_vec has one entry per equation or none");
    std::vector<double> rtol1(neq + 1, opt.rtol), atol1(neq + 1, opt.atol);
    rtol1[0] = 0.0;
    atol1[0] = 0.0;
    for (std::size_t k = 0; k < neq; ++k) {
        if (rvec) rtol1[k + 1] = opt.rtol_vec[k];
        if (avec) atol1[k + 1] = opt.atol_vec[k];
    }
    const int itol = rvec ? (avec ? 4 : 3) : (avec ? 2 : 1);
    solver.set_tolerances(rtol1, atol1, itol);
    solver.set_force_stiff(opt.force_stiff);

    std::vector<double> y = y0;  // 0-based, as lsoda_update expects
    std::vector<double> yout;
    double t = t_eval[0];
    int istate = 1;

    // iworks = {ml, mu, ixpr, mxstep, mxhnil, mxordn, mxords}
    std::array<int, 7> iworks = {{0, 0, 0, static_cast<int>(opt.max_steps), 0,
                                  opt.max_order_nonstiff, opt.max_order_stiff}};
    // rworks = {tcrit, h0, hmax, hmin}; hmax is the step itself, and the solver
    // inverts it into hmxi -- unlike lsoda-java, which is handed the inverse.
    std::array<double, 4> rworks = {{0.0, opt.h_init, opt.h_max, opt.h_min}};
    const int iopt = 1;  // the optional inputs above are in force
    const int jt = 2;    // Jacobian generated internally, full: LSODA's default

    const LsodaRhs* fp = &f;
    for (std::size_t i = 1; i < t_eval.size(); ++i) {
        const double tout = t_eval[i];
        if (tout == t) {  // a repeated output time asks for the state again
            out.t.push_back(t);
            out.y.push_back(y);
            continue;
        }
        yout.assign(neq + 1, 0.0);
        for (std::size_t k = 0; k < neq; ++k) yout[k + 1] = y[k];
        solver.lsoda(lsoda_detail::rhs_trampoline, neq, yout, &t, tout, 1 /*itask*/, &istate,
                     iopt, jt, iworks, rworks, const_cast<LsodaRhs*>(fp));
        for (std::size_t k = 0; k < neq; ++k) y[k] = yout[k + 1];
        out.t.push_back(t);
        out.y.push_back(y);
        if (istate < 0) {  // stop at the first failure; the caller decides what to do
            out.success = false;
            break;
        }
        istate = 2;  // continue the same integration at the next output time
    }

    out.istate = istate;
    out.steps = solver.get_nst();
    out.f_evals = solver.get_nfe();
    out.jacobians = solver.get_nje();
    out.method = solver.get_mused() == 2 ? "bdf" : "adams";
    return out;
}

/**
 * One internal step at a time: ODEPACK's itask = 2.
 *
 * WHY THIS EXISTS. `lsoda_integrate` above reports the state only at the output
 * times it was handed, and the integrator is free to do whatever it likes in
 * between. MATLAB's `odeset('NonNegative')` is not a property of the output
 * grid, it is a rule the STEP CONTROLLER applies to every accepted step -- it
 * charges a negative excursion as error, and it clips the accepted state and
 * resets the divided-difference table. Reproducing that needs the trajectory
 * one accepted step at a time, which is what this exposes; it is the analogue
 * of scipy's `LSODA.step()`, which the native-Python port drives for the same
 * reason.
 *
 * THE RESTART IS THE HISTORY RESET. There is no way to reach into the vendored
 * solver's Nordsieck array and rewrite it, and no need to: constructing a fresh
 * stepper from a modified state is exactly what resetting the difference table
 * accomplishes, because a first call (istate = 1) rebuilds the history from the
 * initial state alone. The caller therefore expresses "clip and reset" by
 * discarding the stepper and building another one.
 *
 * `step()` may carry `t()` PAST `t1`, as itask = 2 always may; `settle_at_end`
 * interpolates back onto `t1` through the same history, which is what itask = 1
 * does on a continuation call.
 */
class LsodaStepper {
  public:
    LsodaStepper(const LsodaRhs& f, const std::vector<double>& y0, double t0, double t1,
                 const LsodaOptions& opt = LsodaOptions())
        : f_(f), y_(y0), t_(t0), t1_(t1), neq_(y0.size()) {
        if (y0.empty()) throw InputError("LsodaStepper: the initial state is empty");
        if (!(t1 > t0))
            throw InputError("LsodaStepper: the final time must exceed the initial time");
        const bool rvec = !opt.rtol_vec.empty(), avec = !opt.atol_vec.empty();
        if (rvec && opt.rtol_vec.size() != neq_)
            throw InputError("LsodaStepper: rtol_vec has one entry per equation or none");
        if (avec && opt.atol_vec.size() != neq_)
            throw InputError("LsodaStepper: atol_vec has one entry per equation or none");
        std::vector<double> rtol1(neq_ + 1, opt.rtol), atol1(neq_ + 1, opt.atol);
        rtol1[0] = 0.0;
        atol1[0] = 0.0;
        for (std::size_t k = 0; k < neq_; ++k) {
            if (rvec) rtol1[k + 1] = opt.rtol_vec[k];
            if (avec) atol1[k + 1] = opt.atol_vec[k];
        }
        solver_.set_tolerances(rtol1, atol1, rvec ? (avec ? 4 : 3) : (avec ? 2 : 1));
        solver_.set_force_stiff(opt.force_stiff);
        iworks_ = {{0, 0, 0, static_cast<int>(opt.max_steps), 0, opt.max_order_nonstiff,
                    opt.max_order_stiff}};
        rworks_ = {{0.0, opt.h_init, opt.h_max, opt.h_min}};
    }

    /** Advance one accepted step. False once the horizon is reached or LSODA gave up. */
    bool step() {
        if (done_ || istate_ < 0) return false;
        call(2, t1_);
        if (istate_ < 0) return false;
        if (t_ >= t1_) done_ = true;
        return true;
    }

    /** Interpolate the state back onto t1 after itask = 2 stepped past it. */
    void settle_at_end() {
        if (istate_ < 0 || !(t_ > t1_)) return;
        call(1, t1_);
    }

    double t() const { return t_; }
    double t_end() const { return t1_; }
    const std::vector<double>& y() const { return y_; }
    int istate() const { return istate_; }
    bool failed() const { return istate_ < 0; }
    bool finished() const { return done_; }
    std::size_t steps() const { return solver_.get_nst(); }
    std::size_t f_evals() const { return solver_.get_nfe(); }

  private:
    void call(int itask, double tout) {
        std::vector<double> yout(neq_ + 1, 0.0);
        for (std::size_t k = 0; k < neq_; ++k) yout[k + 1] = y_[k];
        solver_.lsoda(lsoda_detail::rhs_trampoline, neq_, yout, &t_, tout, itask, &istate_, 1 /*iopt*/,
                      2 /*jt*/, iworks_, rworks_, &f_);
        for (std::size_t k = 0; k < neq_; ++k) y_[k] = yout[k + 1];
        if (istate_ > 0) istate_ = 2;  // continue the same integration on the next call
    }

    LsodaRhs f_;  ///< held by value: a driver builds steppers from temporaries
    lsoda_impl::LSODA solver_;
    std::vector<double> y_;
    double t_ = 0.0;
    double t1_ = 0.0;
    std::size_t neq_ = 0;
    int istate_ = 1;
    bool done_ = false;
    std::array<int, 7> iworks_ = {{0, 0, 0, 0, 0, 12, 5}};
    std::array<double, 4> rworks_ = {{0.0, 0.0, 0.0, 0.0}};
};

/**
 * `lsoda_integrate` on the same output grid, driven one accepted step at a
 * time so `LsodaOptions::step_stop` can end a window early.
 *
 * ONE STEPPER PER OUTPUT INTERVAL, and that is what keeps the grid honest.
 * itask = 2 may carry `t` PAST the interval's end, so each interval is stepped
 * until the stepper reports it is done and then `settle_at_end` interpolates
 * back onto the requested instant through the same history -- exactly what a
 * continuation call with itask = 1 does. The restart between intervals costs a
 * little accuracy against the single itask = 1 call, which is why this driver is
 * NOT the default: only a caller that asked for a stop test pays for it, and the
 * fluid windows that do hand over a two-entry grid, i.e. one interval.
 *
 * WHEN THE TEST FIRES the state is held for the WHOLE remaining grid rather than
 * the trajectory being truncated. The test says the drift is zero to double
 * precision, so y is that state for every later t and the held values are exact,
 * not padded: a caller reading `final_state()` or a transient grid cannot tell
 * this window from one that was stepped to its end, which is the point -- an
 * early return must not read as a failure to any of them.
 */
inline LsodaSolution lsoda_integrate_stepwise(const LsodaRhs& f, const std::vector<double>& y0,
                                              const std::vector<double>& t_eval,
                                              const LsodaOptions& opt) {
    LsodaSolution out;
    out.t.push_back(t_eval[0]);
    out.y.push_back(y0);

    std::vector<double> y = y0;
    double t = t_eval[0];
    std::size_t nst = 0, nfe = 0;
    bool stopped = false;

    for (std::size_t i = 1; i < t_eval.size(); ++i) {
        const double tout = t_eval[i];
        if (stopped || tout == t) {  // settled, or a repeated instant: report again
            out.t.push_back(tout);
            out.y.push_back(y);
            continue;
        }
        LsodaStepper stepper(f, y, t, tout, opt);
        while (stepper.step()) {
            if (opt.step_stop(stepper.t(), stepper.y())) {
                y = stepper.y();
                t = stepper.t();
                stopped = true;
                break;
            }
        }
        nst += stepper.steps();
        nfe += stepper.f_evals();
        if (stepper.failed()) {
            out.istate = stepper.istate();
            out.success = false;
            out.t.push_back(stepper.t());
            out.y.push_back(stepper.y());
            break;
        }
        if (!stopped) {
            stepper.settle_at_end();
            y = stepper.y();
            t = stepper.t_end();
        }
        out.t.push_back(tout);
        out.y.push_back(y);
    }

    out.steps = nst;
    out.f_evals = nfe;
    // The stepper exposes no nje, and a caller reading it off a stopped window
    // would be reading a count this driver never collected: leave it at 0 rather
    // than report a number that is not the Jacobian count.
    out.jacobians = 0;
    return out;
}

/** Convenience form: integrate from t0 to t1 and report only the end state. */
inline std::vector<double> lsoda_final(const LsodaRhs& f, const std::vector<double>& y0, double t0,
                                       double t1, const LsodaOptions& opt = LsodaOptions()) {
    const std::vector<double> span{t0, t1};
    return lsoda_integrate(f, y0, span, opt).final_state();
}

}  // namespace line

#endif  // LINE_UTIL_LSODA_H
