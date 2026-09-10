/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_SOLVER_FLUID_H
#define LINE_SOLVERS_FLUID_SOLVER_FLUID_H

/**
 * SolverFluid: the `closing` method, a port of `solver_fluid.m`,
 * `solver_fluid_iteration.m` and `solver_fluid_closing.m`.
 *
 * WHAT THE SOLVER DOES. The fluid approximation replaces the integer queue
 * lengths of the CTMC with real-valued masses and follows their mean drift.
 * The drift is built in `fluid_odes.h`; this file integrates it and turns the
 * end state into the usual Q/U/R/T/C/X table.
 *
 * WHY THE INTEGRATION IS AN ITERATION RATHER THAN ONE LONG SOLVE. The steady
 * state is the drift's fixed point, and how long it takes to get there is set
 * by the SLOWEST rate in the model. The reference integrates to
 * 10*iter/min(rate) on iteration `iter`, restarting from the previous end
 * state: each pass buys another ten mean events of the slowest transition. A
 * single solve to a guessed horizon either stops short on a stiff model or
 * wastes most of its steps on one that settled early.
 *
 * A NOTE ON THE CONVERGENCE TEST. `movedMassRatio` is the mass moved over ONE
 * window, and for a mode relaxing at rate r it underestimates the distance
 * still to go by (1-exp(-r*window)). On M/M/1 at rho = 0.9 with `minnormal` the
 * fixed point is Q = 7.021524680 and stopping at `iter_tol = 1e-4` lands on
 * 7.014672, out by 0.1%. Both `solver_fluid_iteration.m` and this port
 * therefore ran every one of their `iter_max` passes, which is what made a
 * fluid solve cost a fixed 150 windows however close it started to the answer.
 * The fix is to stop on what the ratio DROPS: summing the geometric tail,
 * ratio*rho/(1-rho) with rho read off the iteration itself, bounds the distance
 * left rather than the distance just travelled, and needs no rate to stand in
 * for the slowest system mode -- when one does, as a bare drift norm must, the
 * stop lands 3% short. `earlystop` (default true, `options.config.fluid_earlystop`)
 * selects it; `iter_tol > 0` remains the caller's own cruder trade. A FINITE
 * `timespan_end` is a transient request and is exempt from both: it integrates
 * to its end time even once the state has settled.
 *
 * THE PRICE OF RUNNING EVERY PASS is this port's own, and it is small: LSODA is
 * restarted once per pass, so at the default `tol = 1e-4` its error accumulates
 * on a state that is already at the fixed point. On Delay(Z=1) -> PS(c=2), N=6,
 * whose `closing` fixed point is exactly (2,4) and which the reference returns
 * to nine digits, this port is right to 1e-8 by pass 8 and 7e-6 by pass 200.
 * `tol = 1e-6` removes it, at the cost of a different trajectory row count.
 *
 * DOUBLE ONLY. The drift is integrated by LSODA, whose coefficients assume
 * `double` (see `util/lsoda.h`), so a non-`double` backend is refused BY NAME
 * rather than silently narrowed.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/npfqn/npfqn_nonexp_approx.h"
#include "line/api/pfqn/pfqn_marie.h"
#include "line/api/sn/sn_nonmarkov_toph.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/solvers/fluid/fluid_conservation_guard.h"
#include "line/solvers/fluid/fluid_nonhyperbolic.h"
#include "line/solvers/fluid/fluid_matrix.h"
#include "line/solvers/fluid/fluid_diffusion.h"
#include "line/solvers/fluid/fluid_aoi.h"
#include "line/solvers/fluid/fluid_mfq.h"
#include "line/solvers/fluid/fluid_mfq_prio.h"
#include "line/solvers/fluid/fluid_tbi.h"
#include "line/solvers/fluid/fluid_mvn_rectangle.h"
#include "line/solvers/fluid/fluid_odes_statedep.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/util/error.h"
#include "line/util/lsoda.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** Controls, defaulting to `SolverOptions('Fluid')` in the reference. */
struct FluidOptions {
    std::string method = "default";
    double tol = 1e-4;            ///< absolute and relative tolerance handed to the integrator
    double iter_tol = 0.0;        ///< >0 stops early when the moved-mass ratio falls below it; 0 runs to iter_max, as the reference does
    bool earlystop = true;        ///< `options.config.fluid_earlystop`: stop on the geometric tail of the window iteration
    std::size_t iter_max = 200;   ///< cap on outer integrations
    /**
     * `options.config.nonmkvorder`: the phase budget `sn_nonmarkov_toph` spends
     * on a non-Markovian service law. The fluid path always takes the PH fit,
     * so this is the Bernstein order.
     */
    std::size_t nonmkv_order = 20;
    double timespan_end = std::numeric_limits<double>::infinity();
    std::vector<double> init_sol;  ///< initial state; empty selects the default below
    /**
     * `options.config.kp_init_sol`: the `kp` method's initial state, in the
     * KO-PENDER layout -- one offset counter walking the stations in order, an
     * arrival-phase block at each EXT station-class and a service-phase block at
     * every other, with no mass returning to the source.
     *
     * NOT `init_sol`, which is laid out for the CLOSING state vector: the two
     * can have the same length on the same model, so sharing one field lets a
     * closing-layout seed be consumed here, silently zeroing the source phase
     * mass and with it the whole network. Empty selects the stationary arrival
     * phase; a wrong-sized seed is refused rather than ignored.
     */
    std::vector<double> kp_init_sol;
    /**
     * `options.config.init_cov`: the `kp` method's initial covariance Sigma(0),
     * dim-by-dim in the same layout as `kp_init_sol`.
     *
     * A caller that carries a DISTRIBUTION across a handoff supplies the second
     * moment beside the mean, so the next stage does not restart from a point
     * mass it never had. Empty keeps the default diag(theta) - theta theta' of
     * the initial arrival phase.
     */
    Matrix<double> init_cov;
    double softmin_alpha = 20.0;   ///< sharpness of the 'softmin' smoothing
    double pstar = 20.0;           ///< exponent of the 'pnorm' smoothing
    /**
     * Opt in to the p-norm under `matrix`/`default` too, which is what
     * `options.config.pstar` does in MATLAB, the JAR and native Python. The
     * exponent above is a default, not a request, so it cannot serve as the
     * flag: leaving it at 20 must keep the hard min() under `matrix`.
     */
    bool pstar_set = false;
    /**
     * `options.config.fork_join`: which fork-join arm the fixed point takes,
     * 'default'/'mmt'/'fjt' or 'ht'. Carried here so that a fluid solve of a
     * fork-join model selects the same transform an MVA or NC solve of it
     * would; see the `has_fork` branch of fluid_runner.h.
     */
    std::string fork_join = "default";
    double timestep = 0.01;        ///< 'diffusion' Euler-Maruyama step
    unsigned long seed = 23000;    ///< 'diffusion' RNG seed
    /**
     * `options.stiff`: integrate the closing family with the explicit stiff arm
     * of `fluid_stiff.h` rather than with LSODA.
     *
     * THE DEFAULT IS FALSE WHERE THE REFERENCE'S IS TRUE, and that is not a
     * downgrade. `options.stiff = true` selects ode15s, a variable-order BDF
     * code; LSODA is a variable-order Adams/BDF code that switches to BDF on
     * its own stiffness detector, so the reference's default arm is the one
     * already taken here. Setting this selects the four-stage Rosenbrock
     * method, which is the family ode23s belongs to -- the reference's OTHER
     * arm -- so the flag names the integrator that is actually different.
     */
    bool stiff = false;
    /**
     * `options.config.hide_immediate`: fold the Immediate-rate transitions into
     * the timed ones by stochastic complementation before integrating. Off in
     * the reference too, which reaches `ode_eliminate_immediate` only when the
     * caller asks for it.
     */
    bool hide_immediate = true;
    /**
     * `options.config.aoi_preemption`: the preemption (bufferless) or
     * replacement (single buffer) probability of the AoI branch of `mfq`.
     * Negative selects the value the scheduling policy implies.
     */
    double aoi_preemption = -1.0;
    /**
     * `options.config.moment_sigma2` and `options.config.moment_cov`: the second
     * moment the drift's non-linear terms are closed with.
     *
     * A CALLER DOES NOT SET THIS. `solver_fluid_moments` does, once per sweep of
     * its outer fixed point, and it is on the options because the mean solve is
     * the ORDINARY closing integration -- the closure has to reach the drift
     * without a second entry point that could drift from the first.
     */
    FluidClosure closure;
    /**
     * `options.config.moment_maxstate`: the largest phase-resolved state the
     * moment-closure methods will build a covariance over. The Lyapunov solve is
     * cubic in it, so this is a refusal threshold and not a tuning knob; it also
     * decides whether `default` resolves to `minnormal` at all.
     */
    std::size_t moment_maxstate = 200;
    /**
     * `options.config.dae_maxstate` and `options.config.dae_maxcov`: the DAE
     * route's own two refusal thresholds, on the simultaneous solve and on the
     * covariance it integrates alongside the mean.
     *
     * ZERO MEANS NOT SET, and that is what makes them options rather than a
     * second copy of the defaults: the route reads them off `FluidDaeOptions`,
     * whose own values a caller may pin directly, and `fluid_dae_options` lets
     * an explicit pin stand wherever the options are silent. A user reaching
     * for the knob writes the option, as in the other three codebases; a test
     * pinning one writes the struct.
     */
    std::size_t dae_maxstate = 0;
    std::size_t dae_maxcov = 0;
    /**
     * `options.config.highvar`: which non-exponential FCFS correction the
     * analyzer's outer refit loop applies, `interp` (the WSC 2020 diffusion
     * interpolation) or `default`/`none`/`hvmva` (no rescaling, so the loop
     * converges after one sweep).
     *
     * THE DEFAULT IS `default`, i.e. NO rescaling, because that is what
     * `SolverOptions.m:127` sets for FLD -- NC is the solver that defaults to
     * `interp`. The refit loop still runs: with no rescaling it refits each FCFS
     * station to a Coxian at its own mean and SCV, which is an identity on a
     * declared Coxian and a two-moment reduction on anything else, and it
     * converges in two sweeps because `eta` is constant. See fluid_nonexp.h.
     */
    std::string highvar = "default";
    /**
     * `options.config.rate_traj = {tgrid, Mmat}`: a caller-supplied per-EVENT
     * multiplier, which is what the coupled LN layer transient injects.
     * `Mmat` must have one row per event of the closing ODE.
     */
    FluidRateMult rate_traj;
    /**
     * `options.config.nhpp_sched`: the (station, class) pairs whose SOURCE
     * carries a non-homogeneous intensity, which the drift is to follow exactly
     * rather than at its time average.
     *
     * IT IS A LIST AND NOT A FLAG, and that is the reference's design. A model
     * can declare an NHPP and still be solved at the nominal -- that is what
     * `solver_fluid.m` does for a steady-state request -- so the schedule enters
     * the drift only when a caller asks for it, which in the reference is
     * `@@SolverFLD/getTranAvg` through `local_detect_nhpp`. `fluid_detect_nhpp`
     * below is that detector; a caller that wants the nominal simply does not
     * call it.
     */
    std::vector<std::pair<std::size_t, std::size_t> > nhpp_sched;  ///< 1-based (station, class)
    /**
     * `options.config.rate_sched`: explicit per-(station, class) rate
     * trajectories, the third source `solver_fluid_ratemult` composes. Used by
     * the coupled LN layer transient to inject time-varying inter-layer demand
     * through the same station-class -> event expansion the NHPP path uses.
     */
    struct RateSched {
        std::size_t station = 0;  ///< 1-based
        std::size_t cls = 0;      ///< 1-based
        std::vector<double> tgrid;
        std::vector<double> rates;
        /** The nominal baked into rate_base; <= 0 selects `Mu{i}{c}(1)`. */
        double nominal = -1.0;
    };
    std::vector<RateSched> rate_sched;
};

/**
 * The second-order results of the moment-closure methods, i.e. what
 * `@@SolverFLD/getMoments` returns.
 *
 * EMPTY FOR EVERY FIRST-ORDER METHOD, which compute no second moment at all --
 * `has_moments` on the solution says which. Reporting zeros instead would be a
 * variance of zero, which is a claim and not an absence.
 */
struct FluidMomentReport {
    Matrix<double> Sigma;            ///< state-level covariance, on range(D)
    Matrix<double> QVar, QStd;       ///< per station and class queue-length variance
    std::vector<double> sigma2;      ///< per-station population variance
    std::vector<double> refinement;  ///< the 1/N correction, `refined` only
    std::size_t outer_iters = 0;
    /// state coordinates of each (station,class): `Sigma` is indexed by SERVICE
    /// PHASE, so reading a per-class population off it needs this map
    std::vector<std::vector<std::vector<std::size_t>>> class_block;
};

/**
 * One point of a transient trajectory: the metrics at time `t`.
 *
 * `getTranAvg` in the reference returns QNt/UNt/TNt as (station x class) cell
 * arrays of time series; this carries the same information sampled at a grid,
 * which is what a caller plots or integrates.
 */
struct FluidTranPoint {
    double t = 0.0;
    Matrix<double> QN, UN, TN;
    /**
     * Per-(station,class) queue-length VARIANCE at this instant, empty where the
     * method carries no second moment along the trajectory. Only the `dae`
     * route fills it, and only below `FluidDaeOptions::maxcov`: the moment
     * closures evaluate their whole transient at the single stationary variance,
     * so a per-point variance would be the same number repeated.
     */
    Matrix<double> QVar;
};

/** What the analyzer returns, in the same shape as the MVA solver's result. */
struct FluidSolution {
    Matrix<double> QN, UN, RN, TN;
    std::vector<double> CN, XN;
    std::vector<double> xvec;   ///< the converged fluid state
    std::size_t iters = 0;
    /**
     * `iter` of `solver_fluid_analyzer.m`: the FCFS non-exponential refit sweeps.
     * Zero when the model has no FCFS station or the method does not refit, which
     * is the reference's own "the loop was never entered".
     */
    std::size_t refit_sweeps = 0;
    std::string method = "closing";
    /**
     * `result.solverSpecific.aoiResults`: set only by the AoI branch of `mfq`,
     * where the age laws, and not QN/RN, are the answer.
     */
    bool has_aoi = false;
    AoiSolution aoi;
    /**
     * `result.solverSpecific.moments`: set only by `minnormal` and `refined`.
     * `closure` is the variance those methods converged to, kept so that a
     * transient asked for afterwards integrates the SAME Gaussian drift the
     * steady-state table was read from rather than the first-order one.
     */
    bool has_moments = false;
    FluidMomentReport moments;
    FluidClosure closure;
};

namespace detail {

/**
 * Port of `solver_fluid_initsol.m`: THE initial condition of every fluid
 * integration, and the only one in this port.
 *
 * IT IS NOT THE `y0` OF `solver_fluid.m`. That vector -- the even spread of a
 * closed population over the stations that serve its class -- is the
 * reference's `ydefault`, reached only when the integrator rejects the real
 * initial point. `solver_fluid_analyzer.m:25-27` fills `options.init_sol` with
 * `solver_fluid_initsol` BEFORE the method switch, so the even spread is never
 * what a solve starts from.
 *
 * WHAT `solver_fluid_initsol` ACTUALLY RETURNS, and why it is this short.
 * It decodes `sn.state`, which for any model that did not call setState is what
 * `Network.initDefault` wrote through `State.fromMarginalAndStarted`. That
 * encoder puts every job it starts in PHASE ONE
 * (`init = spaceClosedSingle(K(r),0); init(1) = si(r)`), never enumerating the
 * phase assignment. Outside a Source the decode then gives phase one
 * `nir(r) - sum_{k>=2} kir(r,k)`, and with every kir(r,k>=2) zero that is the
 * whole per-station population back again. So the round trip through the
 * encoder is an identity, and what is left of `solver_fluid_initsol` is the
 * PLACEMENT `initDefault` computed, written into the phase-one entries.
 *
 * An open class instead holds the unit job pool at its Source, which is
 * `init(1) = 1` in the encoder's EXT branch, and nothing anywhere else.
 *
 * THE PLACEMENT IS `initDefault`'s AND NOT "ALL OF IT AT THE REFERENCE
 * STATION". The reference station takes as much of the population as its
 * `classcap`/`cap` allows and SPILLS the excess onto the remaining stations in
 * ascending order, erroring when it never fits. On a model with no explicit
 * capacity the two coincide -- `refreshCapacity` gives every station the
 * population of the chains that reach it -- but `setCapacity(k)` below the
 * population makes them different vectors, and then the ODE would be started
 * from a point the model is never in.
 */
template <class T>
std::vector<std::vector<double>> fluid_initsol_placement(const qn::NetworkStruct<T>& sn,
                                                         const FluidLayout& L) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> nplace(M, std::vector<double>(K, 0.0));
    std::vector<double> totplace(M, 0.0);
    // An unrefreshed struct carries no capacity table, which is the reference's
    // Inf default and not a zero buffer.
    const bool has_cap = sn.cap.size() == M && sn.classcap.size() == M;

    for (std::size_t r = 0; r < K; ++r) {
        const double pop = sn.classes[r].population;
        if (!std::isfinite(pop)) continue;
        // The reference station first, then every other station in ascending
        // order: `[refist, setdiff(1:M, refist)]`.
        std::vector<std::size_t> order;
        const std::size_t rs = sn.classes[r].refstat;
        if (rs >= 1 && rs <= M) order.push_back(rs - 1);
        for (std::size_t i = 0; i < M; ++i)
            if (order.empty() || i != order[0]) order.push_back(i);

        // A Place takes the WHOLE population and never spills, capacity or not
        // (`initDefault.m:24-28`). It carries no fluid coordinates, so the class
        // then contributes nothing to the state vector -- which is also what the
        // reference's decode does, since a Place has NaN rates and
        // `solver_fluid_initsol.m:29,39` appends a column only for a rated class.
        if (rs >= 1 && rs <= M && sn.stations[rs - 1].nodetype == lang::NodeType::Place) {
            nplace[rs - 1][r] = pop;
            totplace[rs - 1] += pop;
            continue;
        }

        double remaining = pop;
        for (std::size_t oi = 0; oi < order.size() && remaining > 0.0; ++oi) {
            const std::size_t j = order[oi];
            // A Source is the reservoir of the open classes, not a holding place
            // for a closed population; a Place belongs to the Petri net encoding
            // and has no fluid coordinates at all.
            if (sn.stations[j].sched == lang::SchedStrategy::EXT) continue;
            if (sn.stations[j].nodetype == lang::NodeType::Place) continue;
            // A pair with no fluid block would otherwise take its mass through
            // `qidx` into the NEXT block, since a disabled pair's index is the
            // following pair's start.
            if (!L.enabled[j][r]) continue;
            const double ccap = (has_cap && sn.classcap[j].size() > r) ? sn.classcap[j][r] : inf;
            const double scap = has_cap ? sn.cap[j] : inf;
            const double avail = std::min(ccap - nplace[j][r], scap - totplace[j]);
            const double take = std::min(remaining, std::max(0.0, avail));
            nplace[j][r] += take;
            totplace[j] += take;
            remaining -= take;
        }
        if (remaining > 0.0)
            throw InputError("solver_fluid_initsol: cannot place the population of class '" +
                             sn.classes[r].name +
                             "': the total capacity of the stations that serve it is insufficient");
    }
    return nplace;
}

/**
 * The DECLARED initial condition, or false when the model declares none.
 *
 * `solver_fluid_initsol.m` decodes `sn.state{isf}`, which is whatever
 * `setState` or `initFromMarginal` left there and only falls back to
 * `initDefault`'s placement when nothing was set. This port used to recompute
 * that placement unconditionally, so `initFromMarginal([0 0; 4 1])` integrated
 * from the DEFAULT marking instead -- the right answer to a different model, and
 * invisible in a long horizon because every initial condition of a closed model
 * converges to the same stationary point.
 *
 * PHASE ONE IS NOT ASSUMED HERE, unlike in the default placement. A declared
 * state carries `kir(r,k)`, jobs in service in phase k, so the reference writes
 * `nir(r) - sum_{k>=2} kir(r,k)` into phase one and `kir(r,k)` into the rest; an
 * `initFromMarginalAndStarted` state has jobs past phase one and folding them
 * forward would start the integration with a different amount of work in flight.
 *
 * A Source is the EXT branch: it holds no fluid population of its own, and its
 * per-phase entries are the reference's `kir` there too.
 *
 * WHICH ROW OF `statespace` IS THE STATE: THE FIRST ONE CARRYING PRIOR MASS.
 * The pair on the wire is the node's WHOLE local space with a prior over its
 * rows, not a one-row state -- a model saved after `initDefault` sends eight
 * rows for a 3-server FCFS queue, with the prior a point mass on row 0. The
 * reference does not read that pair at all: `solver_fluid_initsol.m` decodes
 * `sn.state{isf}`, the single CURRENT state, and every writer emits that state
 * as row 0 of the space it sends (the invariant `state.h` states for
 * `default_init_state`).
 *
 * A PRIOR OVER SEVERAL ROWS DOES NOT CHANGE THE ANSWER, and must not. The fluid
 * limit is not linear in the initial distribution -- the ODE from the mean of
 * two states is not the mean of the two ODEs -- so there is nothing to average;
 * the reference simply keeps integrating from `sn.state`, which `setStatePrior`
 * does not touch. Declining the mixture and falling back to `initDefault`'s
 * placement instead answered `init_state_fcfs_nonexp`'s Prior 3 with the
 * DEFAULT marking (0.175046 for a reference 0.175821) while its Prior 2, the
 * same state under a point-mass prior, was right.
 *
 * THE FALLBACK IS PER STATION, as `sn_declared_marginal`'s is. A station whose
 * row is absent, undecodable or a mixture keeps `initDefault`'s placement while
 * its neighbours keep their declared rows; an all-or-nothing rule zeroed the
 * whole vector the moment ONE station declared and another did not, which on
 * `init_state_fcfs_nonexp` emptied the network and reported QLen 0.
 */
template <class T>
bool fluid_declared_initsol(const qn::NetworkStruct<T>& sn, const FluidLayout& L,
                            const std::vector<std::vector<double>>& nplace,
                            std::vector<double>& y0) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    bool any = false;
    std::vector<double> out(L.nstates, 0.0);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t ind = sn.node_of_station(i + 1);
        const bool ext = sn.stations[i].sched == lang::SchedStrategy::EXT;
        const typename std::map<std::size_t, Matrix<T>>::const_iterator ss =
            sn.statespace.find(ind);
        const typename std::map<std::size_t, std::vector<T>>::const_iterator sp =
            sn.stateprior.find(ind);
        std::size_t pick = static_cast<std::size_t>(-1);
        if (ss != sn.statespace.end() && sp != sn.stateprior.end() && ss->second.cols() > 0 &&
            ss->second.rows() == sp->second.size()) {
            for (std::size_t r = 0; r < ss->second.rows() && pick == static_cast<std::size_t>(-1);
                 ++r)
                if (num_traits<T>::to_double(sp->second[r]) > 0.0) pick = r;
        }
        qn::Marginal<T> m;
        bool decoded = false;
        if (pick != static_cast<std::size_t>(-1)) {
            std::vector<T> row(ss->second.cols());
            for (std::size_t c = 0; c < ss->second.cols(); ++c) row[c] = ss->second(pick, c);
            std::vector<std::size_t> ph(K, 1), shift(K, 0);
            std::size_t w = 0;
            for (std::size_t r = 0; r < K; ++r) {
                ph[r] = sn.phasessz_of(i + 1, r + 1);
                shift[r] = w;
                w += ph[r];
            }
            try {
                m = qn::to_marginal(sn, i + 1, row, ph, shift, sn.nvars_of(ind));
                decoded = m.nir.size() == K && m.kir.size() == K;
            } catch (const Error&) {
                decoded = false;  // a row the encoding cannot decode is not a state
            }
        }
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) continue;
            if (!decoded) {
                // This station keeps `initDefault`'s placement, all of it in
                // phase one, which is where that encoder starts every job.
                if (!ext && nplace[i][r] > 0.0) out[L.qidx[i][r]] = nplace[i][r];
                if (ext && !std::isfinite(sn.classes[r].population)) out[L.qidx[i][r]] = 1.0;
                continue;
            }
            const std::size_t np = L.kic[i][r];
            for (std::size_t k = 0; k < np; ++k) {
                double v = 0.0;
                if (k < m.kir[r].size()) v = num_traits<T>::to_double(m.kir[r][k]);
                if (k == 0 && !ext) {
                    // Phase one absorbs the waiting buffer: `nir - sum_{k>=2} kir`.
                    double served = 0.0;
                    for (std::size_t j = 1; j < m.kir[r].size(); ++j)
                        served += num_traits<T>::to_double(m.kir[r][j]);
                    v = num_traits<T>::to_double(m.nir[r]) - served;
                }
                // A Source reports nir = +Inf by the EXT sentinel; only its
                // per-phase counts are a quantity, and those are finite.
                if (!std::isfinite(v)) v = 0.0;
                out[L.qidx[i][r] + k] = v;
            }
            any = true;
        }
    }
    if (!any) return false;
    y0.swap(out);
    return true;
}

/** The initial condition itself: the declared state, else the placement above. */
template <class T>
std::vector<double> fluid_default_initsol(const qn::NetworkStruct<T>& sn, const FluidLayout& L) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<double> y0(L.nstates, 0.0);
    const std::vector<std::vector<double>> nplace = fluid_initsol_placement(sn, L);
    if (fluid_declared_initsol(sn, L, nplace, y0)) return y0;
    for (std::size_t r = 0; r < K; ++r) {
        if (std::isfinite(sn.classes[r].population)) {
            // Phase one of a block gets `nir - sum_{k>=2} kir`, and every job
            // `initDefault` starts is in phase one, so that is the whole of the
            // station's share of the population.
            // The enabled guard is load-bearing for a Place: it holds a
            // placement but no fluid block, and a disabled pair's `qidx` is the
            // NEXT pair's start, so writing it would corrupt a neighbour.
            for (std::size_t i = 0; i < M; ++i)
                if (L.enabled[i][r] && nplace[i][r] > 0.0) y0[L.qidx[i][r]] = nplace[i][r];
        } else {
            // An open class holds the unit job pool at its source.
            for (std::size_t i = 0; i < M; ++i)
                if (L.enabled[i][r] && sn.stations[i].sched == lang::SchedStrategy::EXT)
                    y0[L.qidx[i][r]] = 1.0;
        }
    }
    return y0;
}

}  // namespace detail

namespace detail {

/**
 * Snap numerical dust to zero, as `filterMetric` does with
 * `outData(outData < FineTol) = 0`.
 *
 * A station a class never reaches still accumulates a few 1e-15 of mass from
 * the integrator, and printing that as a queue length claims a presence the
 * model does not have. The reference clears it, so a fluid table can be
 * compared with an exact one without every unvisited cell reading as a
 * mismatch.
 */
inline void fluid_snap_fine(Matrix<double>& m) {
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (std::fabs(m(i, j)) < lang::GlobalConstants::FineTol) m(i, j) = 0.0;
}

/**
 * Snap the whole result set, then clear the response time wherever the
 * throughput went with it.
 *
 * RespT is a RATIO of two dusty quantities, so it does not look small even
 * when both of its operands do: 3e-15 over 3e-15 is 1, which would report a
 * unit response time at a station the class never visits. Zeroing it with its
 * throughput is what keeps the row consistent.
 */
inline void fluid_snap_all(Matrix<double>& q, Matrix<double>& u, Matrix<double>& r,
                           Matrix<double>& t) {
    fluid_snap_fine(q);
    fluid_snap_fine(u);
    fluid_snap_fine(t);
    fluid_snap_fine(r);
    for (std::size_t i = 0; i < r.rows(); ++i)
        for (std::size_t j = 0; j < r.cols(); ++j)
            if (t(i, j) == 0.0 && q(i, j) == 0.0) r(i, j) = 0.0;
}

/**
 * Mark the (station, class) pairs the model actually routes a job into, read
 * off the per-chain visit ratios `sn.visits`.
 *
 * A fluid result cannot decide that question from the SIZE of QN or TN. Both
 * carry a decaying remnant of the initial state, spread over pairs the class
 * never reaches, and the remnant is whatever the integrator left behind when it
 * stopped: measured at QN = 1.3e-12 and TN = 1.3e-13 on picard05 for
 * `test_CQN_Cox_CS_7`, i.e. ABOVE `GlobalConstants::Zero`, so a threshold on
 * them divides one remnant by the other and reports the station's own service
 * time, 10.0000086, as a response time. The visit ratios come from the routing
 * solve instead, where an unrouted pair is zero to the last bits (2.7e-17
 * there). `visits` is indexed by STATEFUL node, hence `stateful_of_station`.
 *
 * A struct carrying no visit information decides nothing and every pair is
 * reported visited. Mirrors `fluid_visited_pairs.m`.
 */
template <class T>
std::vector<char> fluid_visited_pairs(const qn::NetworkStruct<T>& sn, std::size_t M, std::size_t K) {
    std::vector<char> visited(M * K, 0);
    bool have = false;
    std::size_t max_cols = 0;
    for (std::size_t c = 0; c < sn.visits.size(); ++c) {
        const Matrix<T>& Vc = sn.visits[c];
        if (Vc.rows() == 0 || Vc.cols() == 0) continue;
        have = true;
        max_cols = std::max(max_cols, Vc.cols());
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t isf = sn.stateful_of_station(i + 1);
            if (isf == 0 || isf > Vc.rows()) {
                for (std::size_t r = 0; r < K; ++r) visited[i * K + r] = 1;
                continue;
            }
            for (std::size_t r = 0; r < K && r < Vc.cols(); ++r)
                if (std::fabs(num_traits<T>::to_double(Vc(isf - 1, r))) >
                    lang::GlobalConstants::Zero)
                    visited[i * K + r] = 1;
        }
    }
    if (!have) {
        std::fill(visited.begin(), visited.end(), static_cast<char>(1));
    } else {
        // A class NO visit matrix reaches is not evidence of a non-visit, only of a
        // struct whose visits were refreshed against fewer classes.
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = max_cols; r < K; ++r) visited[i * K + r] = 1;
    }
    return visited;
}

/**
 * The analyzer-level correction `solver_fluid_analyzer.m:206-262` applies to
 * EVERY method branch, after the switch.
 *
 * The per-method solvers report a utilization read straight off the fluid
 * state: `sum(Xservice/mu)/S`, the server time the drift assigns to the class.
 * That quantity is not a utilization -- it can exceed both 1 and the class's
 * own mean population, because the drift's share is an instantaneous rate and
 * not an occupancy. The reference restates it as the smallest of three
 * quantities that each bound it from above: unity, the mean population per
 * server, and the pre-correction total rescaled by the class's share of
 * TN/rate, the share computed from the TRUE service rates rather than from the
 * drift's approximation of them.
 *
 * Omitting this was worth 111% on `cqn_scheduling_dps`: DPS Queue2/Class2 read
 * 0.2249 (its share of a saturated server) where the reference reports 0.1066
 * (its mean population, which is the binding bound). Both were internally
 * consistent, which is why it survived: the uncorrected value satisfies
 * `U = X E[S]` exactly, and only disagrees with the reference.
 *
 * MATLAB's `min` over a vector SKIPS NaN, so a class whose rate is zero
 * (0/0 in the share) must not poison the minimum; the NaN term is dropped, not
 * propagated.
 */
template <class T>
void fluid_analyzer_correct(const qn::NetworkStruct<T>& sn, const Matrix<double>& Q,
                            Matrix<double>& U, Matrix<double>& R, const Matrix<double>& T_) {
    const std::size_t M = Q.rows(), K = Q.cols();
    // A class the model never routes here has no response time, and QN alone
    // does not say so -- see fluid_visited_pairs.
    const std::vector<char> visited = fluid_visited_pairs(sn, M, K);
    const Matrix<double> U0 = U;
    for (std::size_t i = 0; i < M; ++i) {
        double u0sum = 0.0, share_den = 0.0;
        for (std::size_t r = 0; r < K; ++r) {
            if (!(Q(i, r) > 0.0) || !visited[i * K + r]) continue;
            u0sum += U0(i, r);
            const double rate = num_traits<T>::to_double(sn.rates(i, r));
            if (rate != 0.0) share_den += T_(i, r) / rate;
        }
        // A load-dependent station clears alpha(n) times the nominal work, so the
        // bound that divides by its capacity has to divide by the PEAK scaling:
        // Seff = max(c_i, max_n alpha_i(n)), the same T*S/peak convention
        // `solver_ctmc_avg_from_pi` applies. Without load dependence Seff == c and
        // every expression below is unchanged.
        double c = sn.stations[i].nservers;
        for (std::size_t k = 0; k < sn.stations[i].lldscaling.size(); ++k)
            c = std::max(c, num_traits<T>::to_double(sn.stations[i].lldscaling[k]));
        const bool is_delay = sn.stations[i].sched == lang::SchedStrategy::INF;
        for (std::size_t r = 0; r < K; ++r) {
            if (!(Q(i, r) > 0.0) || !visited[i * K + r]) {
                U(i, r) = 0.0;
                R(i, r) = 0.0;
                continue;
            }
            if (is_delay) {
                U(i, r) = Q(i, r);
                continue;
            }
            double best = 1.0;
            if (std::isfinite(c) && c > 0.0) best = std::min(best, Q(i, r) / c);
            const double rate = num_traits<T>::to_double(sn.rates(i, r));
            if (rate != 0.0 && share_den != 0.0)
                best = std::min(best, u0sum * (T_(i, r) / rate) / share_den);
            U(i, r) = best;
            if (T_(i, r) != 0.0) R(i, r) = Q(i, r) / T_(i, r);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (std::isnan(U(i, r))) U(i, r) = 0.0;
            if (std::isnan(R(i, r))) R(i, r) = 0.0;
        }
}

}  // namespace detail


/**
 * Read Q/U/R/T off ONE fluid state, for the closing family.
 *
 * Factored out because the transient needs exactly this at every point of the
 * trajectory, and a second copy would drift from the steady-state one. `m` is
 * the resolved method name: only `statedep` changes the rules here, and it does
 * so at FCFS stations (see the mean-service-time share below).
 */
template <class T>
void fluid_closing_metrics(const qn::NetworkStruct<T>& sn, const FluidOdeSystem& sys,
                           const std::string& m, const std::vector<double>& xs, Matrix<double>& Q,
                           Matrix<double>& U, Matrix<double>& R, Matrix<double>& T_) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const FluidLayout& L = sys.layout;
    Q = Matrix<double>(M, K, 0.0);
    U = Matrix<double>(M, K, 0.0);
    R = Matrix<double>(M, K, 0.0);
    T_ = Matrix<double>(M, K, 0.0);
    // Queue length is the mass of the (station, class) block.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            double q = 0.0;
            for (std::size_t k = 0; k < L.kic[i][r]; ++k) q += xs[L.qidx[i][r] + k];
            Q(i, r) = q;
        }

    // Throughput, and the per-phase service mass the utilization is read from.
    std::vector<std::vector<std::vector<double>>> xservice(M, std::vector<std::vector<double>>(K));
    for (std::size_t i = 0; i < M; ++i) {
        const lang::SchedStrategy sc = sn.stations[i].sched;
        double xi = 0.0;
        for (std::size_t r = 0; r < K; ++r) xi += Q(i, r);
        double wxi = 0.0;
        for (std::size_t r = 0; r < K; ++r)
            wxi += (r < sn.stations[i].schedparam.size()
                        ? num_traits<T>::to_double(sn.stations[i].schedparam[r])
                        : 1.0) *
                   Q(i, r);
        const double c = sn.stations[i].nservers;
        // The capacity term psi(n) = min(n,c)*alpha(n), and not the bare min: a
        // load-dependent station clears alpha(n) times the nominal work, so the
        // bare min would drop the scaling from Tput and Util while the ODE applied
        // it. With no load dependence `sys.lld[i]` is empty and this IS min(xi,c).
        const double served = fluid_capacity_closure(xi, c, 0.0, sys.lld[i], false).h;

        // `statedep` shares an FCFS server by MEAN SERVICE TIME, not by head
        // count: a class that occupies a server for longer draws a larger share.
        // `solver_fluid_closing.m` applies that weighting for this method only,
        // and additionally overrides TN with sum(Xservice) -- i.e. WITHOUT the
        // completion probability phi that every other branch carries.
        const bool fcfs_statedep = (m == "statedep") && (sc == lang::SchedStrategy::FCFS);
        std::vector<double> wmean(K, 0.0);
        double wni = lang::GlobalConstants::FineTol;
        if (fcfs_statedep) {
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.enabled[i][r]) continue;
                mam::Map<double> mp;
                const std::size_t nn = sn.service[i][r].D0.rows();
                mp.D0 = Matrix<double>(nn, nn, 0.0);
                mp.D1 = Matrix<double>(nn, nn, 0.0);
                for (std::size_t a = 0; a < nn; ++a)
                    for (std::size_t bb = 0; bb < nn; ++bb) {
                        mp.D0(a, bb) = num_traits<T>::to_double(sn.service[i][r].D0(a, bb));
                        mp.D1(a, bb) = num_traits<T>::to_double(sn.service[i][r].D1(a, bb));
                    }
                wmean[r] = mam::map_mean(mp);
                wni += wmean[r] * Q(i, r);
            }
        }

        for (std::size_t r = 0; r < K; ++r) {
            xservice[i][r].assign(L.kic[i][r], 0.0);
            if (!L.enabled[i][r]) continue;
            std::vector<double> mu, phi;
            detail::fluid_mu_phi(sn.service[i][r], mu, phi);
            const std::size_t b = L.qidx[i][r], n = L.kic[i][r];
            double tn = 0.0;
            if (fcfs_statedep) {
                for (std::size_t k = 0; k < n; ++k)
                    xservice[i][r][k] = xs[b + k] * mu[k] * wmean[r] / wni * served;
                double s2 = 0.0;
                for (std::size_t k = 0; k < n; ++k) s2 += xservice[i][r][k];
                T_(i, r) = s2;  // the reference's sum(Xservice) override
                continue;
            }
            for (std::size_t k = 0; k < n; ++k) {
                double mass = xs[b + k];
                switch (sc) {
                    case lang::SchedStrategy::EXT:
                        // The source holds unit mass: phase one carries the rest.
                        if (k == 0) {
                            double rest = 0.0;
                            for (std::size_t p = 1; p < n; ++p) rest += xs[b + p];
                            mass = 1.0 - rest;
                        }
                        tn += mass * mu[k] * phi[k];
                        xservice[i][r][k] = mass * mu[k];
                        break;
                    case lang::SchedStrategy::INF:
                        tn += mass * mu[k] * phi[k];
                        xservice[i][r][k] = mass * mu[k];
                        break;
                    case lang::SchedStrategy::DPS: {
                        const double w = r < sn.stations[i].schedparam.size()
                                             ? num_traits<T>::to_double(sn.stations[i].schedparam[r])
                                             : 1.0;
                        if (wxi > 0.0) {
                            tn += mass * mu[k] * phi[k] * w / wxi * served;
                            xservice[i][r][k] = mass * mu[k] * w / wxi * served;
                        }
                        break;
                    }
                    default:  // PS, FCFS, SIRO and the rest share the servers
                        if (xi > 0.0) {
                            tn += mass * mu[k] * phi[k] / xi * served;
                            xservice[i][r][k] = mass * mu[k] / xi * served;
                        }
                        break;
                }
            }
            T_(i, r) = tn;
        }
    }

    // Utilization: the service mass divided by the phase rate is the time a
    // server spends on it; a delay reports the queue length itself.
    for (std::size_t i = 0; i < M; ++i) {
        const bool is_delay = sn.stations[i].sched == lang::SchedStrategy::INF;
        for (std::size_t r = 0; r < K; ++r) {
            if (!L.enabled[i][r]) continue;
            std::vector<double> mu, phi;
            detail::fluid_mu_phi(sn.service[i][r], mu, phi);
            double u = 0.0;
            for (std::size_t k = 0; k < xservice[i][r].size(); ++k)
                if (xservice[i][r][k] > 0.0 && mu[k] > 0.0) u += xservice[i][r][k] / mu[k];
            // Divided by the PEAK scaling, Seff = max(c, max_n alpha(n)), which is
            // c itself without load dependence.
            double c = sn.stations[i].nservers;
            for (std::size_t k = 0; k < sys.lld[i].size(); ++k) c = std::max(c, sys.lld[i][k]);
            U(i, r) = (is_delay || !std::isfinite(c) || c <= 0.0) ? u : u / c;
        }
    }

    // Response time by Little's law, which the reference also applies here.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            if (T_(i, r) > 0.0) R(i, r) = Q(i, r) / T_(i, r);

    // A Source and a Sink report no queue length, utilization or response
    // time: `getAvgHandles` disables those three metric kinds there, which is
    // the same rule the MVA runner's `filter_metric` applies. The fluid state
    // does carry mass at an EXT source -- the unit job pool the drift needs --
    // and reporting it would show a queue that does not exist. Throughput is
    // kept, since that is the arrival rate.
    for (std::size_t i = 0; i < M; ++i) {
        const qn::NodeType nt = sn.stations[i].nodetype;
        if (nt != qn::NodeType::Source && nt != qn::NodeType::Sink) continue;
        for (std::size_t r = 0; r < K; ++r) {
            Q(i, r) = 0.0;
            U(i, r) = 0.0;
            R(i, r) = 0.0;
        }
    }

    detail::fluid_snap_all(Q, U, R, T_);
}

namespace detail {

// ---------------------------------------------------------------------------
// `solver_fluid_ratemult.m`: the time-varying per-event rate multiplier.
// ---------------------------------------------------------------------------
/**
 * `local_nhpp_steps`: a step-faithful (time, rate) sampling of a
 * piecewise-constant intensity over [t0, thi].
 *
 * Each segment contributes TWO samples, at its start and just before its end,
 * so that the clamped-linear `fluid_interpcols` reproduces a STEP rather than a
 * ramp between segment values. Sampling once per segment would interpolate
 * across the whole segment and integrate an intensity the model never has.
 */
inline void fluid_nhpp_steps(const std::vector<double>& bp, const std::vector<double>& seg_rate,
                             bool cyclic, double t0, double thi, std::vector<double>& seg_t,
                             std::vector<double>& seg_r) {
    seg_t.clear();
    seg_r.clear();
    if (bp.size() < 2 || seg_rate.empty() || !(thi > t0)) return;
    const double period = bp.back() - bp.front();
    std::vector<double> bounds;
    bounds.push_back(t0);
    if (cyclic && period > 0.0) {
        const long kmax = static_cast<long>(std::ceil((thi - t0) / period)) + 2;
        for (long k = -1; k <= kmax; ++k)
            for (std::size_t a = 0; a < bp.size(); ++a)
                bounds.push_back(bp[a] + static_cast<double>(k) * period);
    } else {
        for (std::size_t a = 0; a < bp.size(); ++a) bounds.push_back(bp[a]);
    }
    bounds.push_back(thi);
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::remove_if(bounds.begin(), bounds.end(),
                                [&](double v) { return v < t0 || v > thi; }),
                 bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    if (bounds.size() < 2) return;

    // The rate in force on a segment, read at its MIDPOINT so a boundary never
    // decides which segment is sampled.
    const auto rate_at = [&](double t) -> double {
        double offset = t - bp.front();
        if (cyclic) {
            if (period > 0.0) {
                offset = std::fmod(offset, period);
                if (offset < 0.0) offset += period;
            } else {
                offset = 0.0;
            }
        } else if (offset < 0.0 || offset >= period) {
            return 0.0;  // zero past a non-cyclic horizon, as the reference
        }
        const double pos = bp.front() + offset;
        std::size_t idx = seg_rate.size() - 1;
        for (std::size_t k = 1; k < bp.size(); ++k)
            if (pos < bp[k]) {
                idx = k - 1;
                break;
            }
        return idx < seg_rate.size() ? seg_rate[idx] : 0.0;
    };

    const double neps = std::max(1e-9, 1e-6 * (thi - t0));
    for (std::size_t k = 0; k + 1 < bounds.size(); ++k) {
        const double a = bounds[k], b = bounds[k + 1];
        const double r = rate_at(0.5 * (a + b));
        seg_t.push_back(a);
        seg_r.push_back(r);
        seg_t.push_back(std::max(a + neps, b - neps));
        seg_r.push_back(r);
    }
}

/** `local_merge`: elementwise product of two multipliers on the union grid. */
inline FluidRateMult fluid_ratemult_merge(const FluidRateMult& a, const FluidRateMult& b,
                                          std::size_t nevents) {
    if (a.empty()) return b;
    if (b.empty()) return a;
    std::vector<double> tg = a.tgrid;
    tg.insert(tg.end(), b.tgrid.begin(), b.tgrid.end());
    std::sort(tg.begin(), tg.end());
    tg.erase(std::unique(tg.begin(), tg.end()), tg.end());
    FluidRateMult out;
    out.tgrid = tg;
    out.Mmat = Matrix<double>(nevents, tg.size(), 1.0);
    std::vector<double> ca, cb;
    for (std::size_t j = 0; j < tg.size(); ++j) {
        fluid_interpcols(a.tgrid, a.Mmat, tg[j], ca);
        fluid_interpcols(b.tgrid, b.Mmat, tg[j], cb);
        for (std::size_t e = 0; e < nevents; ++e) {
            const double va = e < ca.size() ? ca[e] : 1.0;
            const double vb = e < cb.size() ? cb[e] : 1.0;
            out.Mmat(e, j) = va * vb;
        }
    }
    return out;
}

/**
 * Whether the OPTIONS make the drift non-autonomous.
 *
 * This is a property of what the caller asked for, not of the model: the same
 * model is autonomous at its nominal and time-varying under `nhpp_sched`. It is
 * what `fluid_minnormal_applicable.m:99` consults to steer `default` away from
 * the moment closures, and what `fluid_moment_terms.m:114` raises on when one
 * was asked for by name.
 */
inline bool fluid_has_time_varying_rates(const FluidOptions& opt) {
    return !opt.rate_traj.empty() || !opt.nhpp_sched.empty() || !opt.rate_sched.empty();
}

/** The events sourced at (station i, class c), both 0-based. */
inline std::vector<std::size_t> fluid_events_of(const FluidOdeSystem& sys, std::size_t i,
                                                std::size_t c) {
    std::vector<std::size_t> rows;
    const std::size_t lo = sys.layout.qidx[i][c];
    const std::size_t hi = lo + sys.layout.kic[i][c];
    for (std::size_t e = 0; e < sys.events.size(); ++e)
        if (sys.events[e].event_idx >= lo && sys.events[e].event_idx < hi) rows.push_back(e);
    return rows;
}

/** One (station, class) trajectory expanded onto the event rows. */
inline FluidRateMult fluid_ratemult_rows(const FluidOdeSystem& sys, std::size_t i, std::size_t c,
                                         const std::vector<double>& seg_t,
                                         const std::vector<double>& seg_r, double nominal) {
    FluidRateMult out;
    if (seg_t.empty() || !(nominal > 0.0)) return out;
    const std::vector<std::size_t> rows = fluid_events_of(sys, i, c);
    if (rows.empty()) return out;
    out.tgrid = seg_t;
    out.Mmat = Matrix<double>(sys.events.size(), seg_t.size(), 1.0);
    for (std::size_t j = 0; j < seg_t.size(); ++j)
        for (std::size_t e : rows) out.Mmat(e, j) = seg_r[j] / nominal;
    return out;
}

/**
 * Port of `solver_fluid_ratemult.m`: compose the three time-varying sources into
 * one per-event multiplier, or return an empty one when none is configured.
 *
 * REFERENCE SCALE MISMATCH, reproduced deliberately. The nominal the multiplier
 * divides by is `Mu{i}{c}(1)`, the FIRST PHASE RATE of the process, while the
 * numerator is `getRateAt(t)`, which for a MAPt is `map_lambda` of the segment,
 * i.e. a STATIONARY ARRIVAL rate. The two coincide for an NHPP, where the
 * process has one phase and the phase rate IS the arrival rate, and that is the
 * case the reference documents and uses. For a multi-phase MAPt they are
 * different quantities and the multiplier is off by their ratio. Ported as
 * written, because parity is the contract; flagged here and in
 * `_kb/06-solver-catalog.md` rather than silently corrected.
 */
template <class T>
FluidRateMult fluid_ratemult(const qn::NetworkStruct<T>& sn, const FluidOdeSystem& sys,
                             const FluidOptions& opt) {
    const std::size_t nev = sys.events.size();
    FluidRateMult out;
    if (!opt.rate_traj.empty()) {
        if (opt.rate_traj.Mmat.rows() != nev)
            throw InputError("solver_fluid_ratemult: rate_traj has " +
                             std::to_string(opt.rate_traj.Mmat.rows()) +
                             " rows but the closing ODE has " + std::to_string(nev) + " events");
        out = opt.rate_traj;
    }

    // The horizon a (possibly cyclic) schedule is expanded over. An unbounded
    // timespan takes a few periods, so a cycle is REPRESENTED rather than
    // clamped after its first segment.
    double t0 = 0.0;
    const double tend = opt.timespan_end;

    FluidRateMult nh;
    for (std::size_t a = 0; a < opt.nhpp_sched.size(); ++a) {
        const std::size_t i = opt.nhpp_sched[a].first - 1, c = opt.nhpp_sched[a].second - 1;
        if (i >= sn.nstations || c >= sn.nclasses) continue;
        if (!sys.layout.enabled[i][c]) continue;
        const lang::Distrib<T>& d = sn.service[i][c];
        if (!d.has_schedule()) continue;
        const std::vector<T> mu = d.mu_vec();
        if (mu.empty()) continue;
        const double nominal = num_traits<T>::to_double(mu[0]);
        if (!(nominal > 0.0)) continue;

        std::vector<double> bp, seg_r;
        for (std::size_t k = 0; k < d.sched_bp.size(); ++k)
            bp.push_back(num_traits<T>::to_double(d.sched_bp[k]));
        // `getRateAt`: the stationary arrival rate of the segment's pair, which
        // for a one-phase process is that phase's rate.
        for (std::size_t k = 0; k < d.sched_D0.size(); ++k) {
            mam::Map<T> m;
            m.D0 = d.sched_D0[k];
            m.D1 = d.sched_D1[k];
            seg_r.push_back(num_traits<T>::to_double(mam::map_lambda(m)));
        }
        const double period = bp.empty() ? 0.0 : bp.back() - bp.front();
        double thi = tend;
        if (!std::isfinite(thi))
            thi = (std::isfinite(period) && period > 0.0) ? t0 + 3.0 * period : t0 + 1.0;
        std::vector<double> seg_t, seg_v;
        fluid_nhpp_steps(bp, seg_r, d.sched_cyclic, t0, thi, seg_t, seg_v);
        nh = fluid_ratemult_merge(nh, fluid_ratemult_rows(sys, i, c, seg_t, seg_v, nominal), nev);
    }

    FluidRateMult rs;
    for (std::size_t a = 0; a < opt.rate_sched.size(); ++a) {
        const FluidOptions::RateSched& e = opt.rate_sched[a];
        const std::size_t i = e.station - 1, c = e.cls - 1;
        if (i >= sn.nstations || c >= sn.nclasses) continue;
        if (!sys.layout.enabled[i][c]) continue;
        if (e.tgrid.size() != e.rates.size() || e.tgrid.empty())
            throw InputError("solver_fluid_ratemult: rate_sched tgrid and rates must be "
                             "non-empty and of equal length");
        double nominal = e.nominal;
        if (!(nominal > 0.0)) {
            const std::vector<T> mu = sn.service[i][c].mu_vec();
            if (mu.empty()) continue;
            nominal = num_traits<T>::to_double(mu[0]);
        }
        if (!(nominal > 0.0)) continue;
        rs = fluid_ratemult_merge(rs, fluid_ratemult_rows(sys, i, c, e.tgrid, e.rates, nominal),
                                  nev);
    }

    out = fluid_ratemult_merge(out, nh, nev);
    out = fluid_ratemult_merge(out, rs, nev);
    return out;
}

/**
 * `slowrate` of the reference: the smallest phase rate among the service
 * processes the layout enables, which is what sets every integration horizon
 * here. Finite rates above `tol` only, falling back to 1 when the model has
 * none, exactly as `solver_fluid.m` does.
 */
template <class T>
double fluid_slow_rate(const qn::NetworkStruct<T>& sn, const FluidLayout& L, double tol) {
    double min_rate = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!L.enabled[i][r]) continue;
            const lang::Distrib<T>& d = sn.service[i][r];
            for (std::size_t k = 0; k < d.D0.rows(); ++k) {
                const double mu = -num_traits<T>::to_double(d.D0(k, k));
                if (mu > tol && std::isfinite(mu)) min_rate = std::min(min_rate, mu);
            }
        }
    return std::isfinite(min_rate) ? min_rate : 1.0;
}

/**
 * The method switch of `solver_fluid_analyzer.m`, without its trailing
 * correction. Kept separate so the correction below runs on EVERY branch, as
 * it does in the reference; folding it into each branch's exit would repeat it
 * six times and let one path drift.
 *
 * Refuses by name anything the port does not cover, rather than returning a
 * number computed by the wrong model.
 */
template <class T>
FluidSolution fluid_dispatch(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "solver_fluid: the fluid solver integrates its drift with LSODA, whose coefficients "
            "assume double precision; rerun with --arith double");
    // `fluid.<name>` is the same method under its qualified spelling.
    std::string m = opt.method;
    if (m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    StateDepKind sd_kind = StateDepKind::StateDep;
    bool statedep_family = false;
    // `solver_fluid_analyzer.m` routes default, matrix AND pnorm to
    // solver_fluid_matrix -- `ode_pnorm.m` is never reached under the name
    // `pnorm`, so pnorm here means "the matrix method with p-norm smoothing".
    // `@@SolverFLD/runAnalyzer.m` resolves the method BEFORE the analyzer sees
    // it, and the resolution depends on the model, not just the name:
    //   a Cache        -> rmf
    //   a DPS station  -> closing   (the matrix method cannot express DPS)
    //   otherwise      -> matrix
    // and `matrix`/`pnorm` asked for EXPLICITLY on a DPS model is an error, not
    // a silent downgrade. Without this gate the C++ ran matrix on a DPS model
    // and redistributed the population differently from every other codebase.
    bool has_dps = false, has_cache = false;
    for (const auto& st : sn.stations)
        if (st.sched == lang::SchedStrategy::DPS) has_dps = true;
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == qn::NodeType::Cache) has_cache = true;
    if ((m == "matrix" || m == "pnorm") && has_dps)
        throw UnsupportedError(
            "solver_fluid: the matrix method does not support DPS scheduling; use method "
            "'closing' (which is what 'default' selects on a DPS model)");
    if (m == "default") {
        if (has_cache)
            throw UnsupportedError(
                "solver_fluid: a Cache model resolves to the 'rmf' fluid method, which is ported "
                "in fluid_cacheqn.h and reached through solver_fluid_run_analyzer (fluid_runner.h); this "
                "function is solver_fluid_analyzer alone and cannot call it without a cyclic "
                "include");
        if (has_dps) m = "closing";
    }
    const bool matrix_family = (m == "default" || m == "matrix" || m == "pnorm");
    if (m == "statedep") {
        statedep_family = true;
        sd_kind = StateDepKind::StateDep;
    } else if (m == "softmin") {
        statedep_family = true;
        sd_kind = StateDepKind::SoftMin;
    } else if (!(matrix_family || m == "closing" || m == "tbi" || m == "diffusion" || m == "mfq")) {
        throw UnsupportedError("solver_fluid: the '" + opt.method +
                               "' fluid method is not solved here; available are 'closing', "
                               "'statedep', 'softmin', 'pnorm', 'matrix', 'tbi', 'diffusion' and "
                               "'mfq', while 'rmf', 'minnormal', 'refined', 'dae' and 'kp' are "
                               "reached through solver_fluid_run_analyzer (fluid_runner.h), which is the "
                               "port of runAnalyzer's resolution");
    }

    const std::size_t M = sn.nstations, K = sn.nclasses;
    FluidOdeSystem sys = fluid_ode_system(sn);
    // The moment closure travels on the options, so it reaches the drift through
    // the SAME system the metrics are read from; `solver_fluid_moments` is the
    // only caller that fills it and an empty one is the first-order drift.
    sys.closure = opt.closure;
    // `solver_fluid_odes.m:127-134`: the time-varying rate multiplier is built
    // once, beside the drift it scales, and an empty one leaves the autonomous
    // closure exactly as it was.
    sys.ratemult = detail::fluid_ratemult(sn, sys, opt);
    const FluidLayout& L = sys.layout;
    if (L.nstates == 0)
        throw InputError("solver_fluid: no station serves any class, so the drift is empty");

    // ---- the slowest rate sets the integration horizon --------------------
    const double min_rate = fluid_slow_rate(sn, L, opt.tol);

    // ---- integrate, restarting from the previous end state ----------------
    std::vector<double> x = opt.init_sol.empty() ? detail::fluid_default_initsol(sn, L) : opt.init_sol;
    if (x.size() != L.nstates)
        throw InputError("solver_fluid: init_sol has " + std::to_string(x.size()) +
                         " entries but the fluid state has " + std::to_string(L.nstates));

    // ---- the exact single Markov-modulated fluid queue --------------------
    if (m == "mfq") {
        // `solver_fluid_analyzer.m` tries the AoI topology FIRST, because it is
        // the more specific one: a capacity-1 or capacity-2 single queue is
        // also a single queue, and the age laws are what that model is for.
        const AoiTopology atop = aoi_is_aoi(sn);
        if (atop.ok) {
            const FluidAoiResult ar = fluid_aoi(sn, atop, opt.aoi_preemption);
            FluidSolution out;
            out.iters = 1;
            out.method = "mfq";
            out.has_aoi = true;
            out.aoi = ar.age;
            out.QN = Matrix<double>(M, K, 0.0);
            out.UN = Matrix<double>(M, K, 0.0);
            out.RN = Matrix<double>(M, K, 0.0);
            out.TN = Matrix<double>(M, K, 0.0);
            out.XN.assign(K, 0.0);
            out.CN.assign(K, 0.0);
            for (std::size_t r = 0; r < K; ++r) {
                out.QN(atop.queue, r) = ar.QN[r];
                out.UN(atop.queue, r) = ar.UN[r];
                out.RN(atop.queue, r) = ar.RN[r];
                out.TN(atop.queue, r) = ar.TN[r];
                out.TN(atop.source, r) = ar.TN[r];
                out.XN[r] = ar.TN[r];
                out.CN[r] = ar.RN[r];
            }
            return out;
        }

        const MfqTopology top = mfq_is_single_queue(sn);
        // Distinct priorities among the open classes send the model to the
        // fluid PRIORITY queue instead of the single fluid-fluid queue.
        bool mixed_prio = false;
        if (top.ok)
            for (std::size_t j = 1; j < top.open_classes.size(); ++j)
                if (sn.classes[top.open_classes[j]].prio != sn.classes[top.open_classes[0]].prio)
                    mixed_prio = true;
        if (mixed_prio) {
            const MfqPrioResult pr = fluid_mfq_prio(sn, top, opt.tol);
            if (pr.fallback) {
                // The reference warns and runs the matrix method instead.
                FluidOptions fb = opt;
                fb.method = "matrix";
                return solver_fluid(sn, fb);
            }
            FluidSolution out;
            out.iters = 1;
            out.method = "mfq";
            out.QN = Matrix<double>(M, K, 0.0);
            out.UN = Matrix<double>(M, K, 0.0);
            out.RN = Matrix<double>(M, K, 0.0);
            out.TN = Matrix<double>(M, K, 0.0);
            out.XN.assign(K, 0.0);
            out.CN.assign(K, 0.0);
            for (std::size_t r = 0; r < K; ++r) {
                out.QN(top.queue, r) = pr.QN[r];
                out.TN(top.queue, r) = pr.TN[r];
                out.TN(top.source, r) = pr.TN[r];
                out.XN[r] = pr.TN[r];
            }
            // The analyzer's post-processing, which is what getAvg reports: a
            // class holding no fluid has neither utilization nor response time,
            // and the utilization of the rest is capped by its own fluid level.
            double ufull = 0.0, tsum = 0.0;
            for (std::size_t r = 0; r < K; ++r)
                if (pr.QN[r] > 0.0) {
                    ufull += pr.UN[r];
                    tsum += pr.TN[r] / num_traits<T>::to_double(sn.rates(top.queue, r));
                }
            const double servers = sn.stations[top.queue].nservers;
            for (std::size_t r = 0; r < K; ++r) {
                if (!(pr.QN[r] > 0.0)) continue;
                const double share =
                    ufull * (pr.TN[r] / num_traits<T>::to_double(sn.rates(top.queue, r))) / tsum;
                out.UN(top.queue, r) = std::min(1.0, std::min(pr.QN[r] / servers, share));
                out.RN(top.queue, r) = pr.QN[r] / pr.TN[r];
                out.CN[r] = out.RN(top.queue, r);
            }
            return out;
        }
        if (top.ok && top.open_classes.size() > 1)
            throw UnsupportedError(
                "fluid mfq: the single fluid-fluid queue analyzes ONE open class, and this model "
                "has several at equal priority; the reference silently reports class 1 only");
        // MFQ IS A SINGLE-QUEUE METHOD AND FALLS BACK, which is what the
        // reference does: solver_fluid_analyzer.m warns "MFQ not applicable:
        // ... Falling back to matrix method" and re-enters solver_fluid_matrix.
        // Refusing instead made 'mfq' -- and therefore its aliases 'butools' and
        // 'aoi' -- reject every multi-station model that MATLAB and native
        // python both answer. This port has no line_warning channel (see
        // mva_dispatch.h), so the substitution is visible in `method` instead,
        // which reports "matrix" exactly as the reference's does.
        if (!top.ok) {
            FluidOptions mopt = opt;
            mopt.method = "matrix";
            return fluid_dispatch(sn, mopt);
        }
        const MfqResult r = fluid_mfq(sn, top, opt.tol);
        FluidSolution out;
        out.iters = 1;  // solved, not iterated
        out.method = "mfq";
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.RN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        out.QN(top.queue, top.cls) = r.QN;
        out.UN(top.queue, top.cls) = r.UN;
        out.RN(top.queue, top.cls) = r.RN;
        out.TN(top.queue, top.cls) = r.TN;
        out.TN(top.source, top.cls) = r.TN;  // the Source reports its arrivals
        out.XN.assign(K, 0.0);
        out.CN.assign(K, 0.0);
        out.XN[top.cls] = r.TN;
        out.CN[top.cls] = r.RN;
        return out;
    }

    // ---- the diffusion approximation: a stochastic trajectory -------------
    if (m == "diffusion") {
        DiffusionOptions dopt;
        dopt.steps = opt.iter_max > 2 ? opt.iter_max : 10000;
        dopt.dt = opt.timestep;
        dopt.seed = opt.seed;
        const DiffusionResult dr = fluid_diffusion(sn, dopt);
        FluidSolution out;
        out.iters = 1;  // a single trajectory
        out.method = "diffusion";
        out.QN = dr.QN;
        out.UN = Matrix<double>(M, K, 0.0);
        out.RN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
            const double c = sn.stations[i].nservers;
            const bool inf_server = !std::isfinite(c);
            for (std::size_t r = 0; r < K; ++r) {
                const double rate = num_traits<T>::to_double(sn.rates(i, r));
                if (rate > 0.0 && std::isfinite(rate)) {
                    // An infinite server clears the whole queue; a single
                    // server clears at most one job's worth at a time.
                    out.TN(i, r) = inf_server ? out.QN(i, r) * rate
                                              : std::min(out.QN(i, r), 1.0) * rate;
                }
                out.UN(i, r) = inf_server ? out.QN(i, r) : std::min(out.QN(i, r) / c, 1.0);
                // TN is zero only to the integrator's accuracy: a class that never visits leaves
                // a ~1e-20 residue in TN too, and a strict > 0 test then divides residue by residue.
                if (out.TN(i, r) > lang::GlobalConstants::Zero)
                    out.RN(i, r) = out.QN(i, r) / out.TN(i, r);
            }
        }
        detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);
        out.XN.assign(K, 0.0);
        out.CN.assign(K, 0.0);
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t rs = sn.classes[r].refstat;
            if (rs >= 1 && rs <= M) out.XN[r] = out.TN(rs - 1, r);
            double q = 0.0;
            for (std::size_t i = 0; i < M; ++i) q += out.QN(i, r);
            if (out.XN[r] > 0.0) out.CN[r] = q / out.XN[r];
        }
        return out;
    }

    LsodaOptions lopt;
    lopt.rtol = opt.tol;
    lopt.atol = opt.tol;

    // ---- the matrix method: one generator, one integration ----------------
    if (matrix_family) {
        // `pnorm` is the matrix method with the smoothing switched on; plain
        // `matrix`/`default` leave pstar at zero, which selects the hard min,
        // unless the caller asked for the smoothing explicitly.
        const double ps = (m == "pnorm" || opt.pstar_set) ? opt.pstar : 0.0;
        const FluidMatrixSystem ms = fluid_matrix_system(sn, x, ps);
        const std::function<void(double, const double*, double*)> mdrift = fluid_matrix_drift(ms);
        const double t1 =
            std::min(opt.timespan_end,
                     10.0 * static_cast<double>(opt.iter_max) / ms.min_rate);
        std::vector<double> xm = fluid_integrate_leg(mdrift, 0.0, t1, ms.x0, lopt);
        for (double& v : xm)
            if (v < 0.0) v = 0.0;

        // DEGENERATE DRIFT: re-integrate with a closed saturation term, do not
        // touch the answer that came back. min(E[n], c) is FLAT above the server
        // count, so a network of saturated stations has a CONTINUUM of fixed
        // points and this method returns whichever one the integrator stopped at
        // -- [9 1] against an exact [5 5] on two identical saturated stations in
        // a closed cycle, and [8 2] with two servers each. The repair is applied
        // to the DRIFT, not to the point: the same trajectory is integrated again
        // with E[min(n, c)] in place of min(E[n], c), which is strictly
        // increasing and so isolates one fixed point.
        //
        // WHY A CLOSURE AND NOT A SMOOTHED min: any smoothing sharp enough to
        // stay faithful to min away from the kink is numerically FLAT far from
        // it. The Boltzmann softmin at alpha = 20 carries a restoring force of
        // exp(-160) at the [9 1] point, and the p-norm trades the two off
        // directly (pstar = 2 recovers [5 5], pstar = 128 gives [8.94 1.06]).
        // The closure escapes the trade-off because its slope comes from the
        // VARIANCE of the marginal rather than from a smoothing width.
        //
        // Only a model that is ACTUALLY degenerate pays for it, so a well-posed
        // model integrates once and is unchanged.
        FluidMatrixSystem msr = ms;
        if (ps <= 0.0 && fluid_matrix_degenerate(ms, mdrift, xm, K)) {
            msr.var_closure = true;
            const std::function<void(double, const double*, double*)> cdrift =
                fluid_matrix_drift(msr);
            std::vector<double> xc = fluid_integrate_leg(cdrift, 0.0, t1, ms.x0, lopt);
            bool finite = xc.size() == xm.size();
            for (std::size_t a = 0; finite && a < xc.size(); ++a)
                if (!std::isfinite(xc[a])) finite = false;
            if (finite) {
                for (double& v : xc)
                    if (v < 0.0) v = 0.0;
                xm = xc;
            } else {
                // A failed repair leaves the unrepaired answer standing rather
                // than turning a wrong number into no number.
                msr.var_closure = false;
            }
        }

        // The same share the drift used, smoothed, closed or neither
        std::vector<double> theta(ms.nstates, 0.0);
        detail::fluid_matrix_theta(msr, xm.data(), theta);

        FluidSolution out;
        out.iters = 1;  // a single integration, unlike the closing iteration
        out.method = (m == "pnorm") ? "pnorm" : "matrix";
        out.xvec = xm;
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.RN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                double q = 0.0, u = 0.0, t = 0.0;
                for (std::size_t a = 0; a < ms.nstates; ++a) {
                    q += ms.sqc(i * K + r, a) * xm[a];
                    u += ms.suc(i * K + r, a) * theta[a];
                    t += ms.stc(i * K + r, a) * theta[a];
                }
                out.QN(i, r) = q;
                // An infinite server reports the queue length itself as its
                // utilization -- there is no capacity to divide by. The SUC map
                // carries 1/S with S substituted by the closed population, so
                // the delay rows have to be restated here, exactly as
                // `solver_fluid.m` does for its UNt.
                out.UN(i, r) =
                    (sn.stations[i].sched == lang::SchedStrategy::INF ||
                     !std::isfinite(sn.stations[i].nservers))
                        ? q
                        : u;
                out.TN(i, r) = t;
                // Little's law, as the reference -- but TN is zero only to the
                // integrator's accuracy, so a class that never visits leaves a
                // residue in q and t alike and a strict > 0 test divides one by
                // the other. See solver_fluid_analyzer.m.
                if (t > lang::GlobalConstants::Zero) out.RN(i, r) = q / t;
            }
        // A Source reports arrivals only. Its states are held at zero with
        // theta = 0, so STC*theta gives it no throughput at all; the reference
        // still shows the arrival rate there, so it is restated from the rates
        // that were injected into the downstream queues.
        for (std::size_t i = 0; i < M; ++i) {
            const qn::NodeType nt = sn.stations[i].nodetype;
            if (nt != qn::NodeType::Source && nt != qn::NodeType::Sink) continue;
            for (std::size_t r = 0; r < K; ++r) {
                out.QN(i, r) = 0.0;
                out.UN(i, r) = 0.0;
                out.RN(i, r) = 0.0;
                if (nt == qn::NodeType::Source && ms.src_arrival.rows() == M)
                    out.TN(i, r) = ms.src_arrival(i, r);
            }
        }
        detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);
        out.XN.assign(K, 0.0);
        out.CN.assign(K, 0.0);
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t rs = sn.classes[r].refstat;
            if (rs >= 1 && rs <= M) out.XN[r] = out.TN(rs - 1, r);
            double q = 0.0;
            for (std::size_t i = 0; i < M; ++i) q += out.QN(i, r);
            if (out.XN[r] > 0.0) out.CN[r] = q / out.XN[r];
        }
        return out;
    }

    // THE PASS LOOP BELOW RESTARTS THE INTEGRATOR ONCE PER PASS, so its local
    // error is paid iter_max times over and lands on a state that is already at
    // the fixed point. At the nominal tol = 1e-4 that accumulated to 3.3e-4 on
    // oqn_basic -- Queue1 QLen 0.10023962 where MATLAB, the JAR and native
    // Python all return 0.10020646880565 -- which reads as a DIFFERENT fluid
    // fixed point and is nothing of the kind: the same run at 1e-6 reproduces
    // the reference to the last digit. So the integrator runs at tol/iter_max
    // while `opt.tol` stays what the caller asked of the FIXED POINT, which is
    // the quantity that tolerance names. The single-integration matrix family
    // above is exempt because nothing is restarted there, and so is
    // `solver_fluid_transient` below, which integrates one grid in one call.
    lopt.rtol = lopt.atol =
        opt.tol / std::max<double>(1.0, static_cast<double>(opt.iter_max));

    // `ode_eliminate_immediate` is applied to the CLOSING drift only, which is
    // the `otherwise` arm of `solver_fluid_odes.m` where the reference applies
    // it: the state-dependent drifts are not a jump/rate system and `tbi`
    // partitions the unreduced one.
    FluidOdeSystem dsys = sys;
    FluidImmediateResult imm_result;
    if (fluid_hide_immediate(sn, opt) && !statedep_family) {
        imm_result = fluid_eliminate_immediate(sn, sys);
        if (imm_result.eliminated) {
            dsys = imm_result.sys;
            // A complemented coordinate receives no transitions at all, so mass left
            // there would sit stranded rather than be integrated. It is PROJECTED, not
            // zeroed: those are jobs, and a cold start puts none there but a warm start
            // from an earlier LN iterate does.
            std::vector<double> xp(x.size(), 0.0);
            for (std::size_t f = 0; f < x.size() && f < imm_result.absorb.rows(); ++f)
                for (std::size_t sidx = 0; sidx < x.size() && sidx < imm_result.absorb.cols();
                     ++sidx)
                    xp[sidx] += x[f] * imm_result.absorb(f, sidx);
            for (std::size_t a = 0; a < x.size(); ++a) x[a] = xp[a];
        }
    }

    const std::function<void(double, const double*, double*)> drift =
        statedep_family ? fluid_drift_statedep(fluid_statedep_system(sn, sd_kind, opt.softmin_alpha,
                                                                     opt.pstar))
                        : fluid_drift(dsys);

    const std::vector<std::vector<std::size_t>> tbi_cells =
        (m == "tbi") ? tbi_partition(sn, TbiOptions().cellsize)
                     : std::vector<std::vector<std::size_t>>();
    // Early stop on the GEOMETRIC TAIL of the window iteration; see the header.
    // The residual cannot go below the integrator's own error, so a request for
    // less than `tol` asks for something unobservable.
    const double drift_tol = std::max(opt.iter_tol, opt.tol);
    const double drift_safety = 0.01;   // headroom, since rho is estimated
    const double min_horizon = 10.0 / min_rate;
    double moved_prev = std::numeric_limits<double>::infinity();
    std::vector<double> rho_hist(3, std::numeric_limits<double>::quiet_NaN());
    int drift_below = 0;
    std::vector<double> drift_buf(x.size(), 0.0);

    double t0 = 0.0;
    std::size_t iter = 0;
    for (; iter < opt.iter_max; ++iter) {
        const double horizon = 10.0 * static_cast<double>(iter + 1) / min_rate;
        const double t1 = std::min(opt.timespan_end, horizon);
        if (!(t1 > t0)) break;
        const std::vector<double> prev = x;
        // A FIXED POINT ENDS THE WINDOW IN CLOSED FORM, and this is what keeps a
        // window that has already converged from becoming a window that never
        // returns. Armed only for an AUTONOMOUS drift, so f(x*) = 0 means
        // x(t) = x* for every later t and the rest of the span is known exactly
        // rather than integrated.
        //
        // THE THRESHOLD IS ROUND-OFF, NOT `drift_tol`. That tolerance (1e-4 by
        // default) says "converged to what the caller asked for", and a state
        // that merely satisfies it is still moving -- cutting the window there
        // was measured to shift results by 1.8e-5 in the Python twin. A
        // normalized residual below GlobalConstants::Zero is the stronger claim
        // that the drift is zero to double precision, which is what makes
        // skipping the remaining span exact instead of approximate.
        //
        // WHY THE WINDOW DOES NOT END ON ITS OWN. A stiff step controller handed
        // a state it is already at cannot pick a step: on the LN layer of
        // test_LQN_13 the Python twin advanced t by 0.011 in 20000 steps from a
        // state with |f| = 1.5e-16, and covered the whole 1000-unit span in 60
        // steps once that state was nudged 1e-6 off the equilibrium. That layer
        // carries an Immediate() coordinate -- an eigenvalue of exactly
        // -GlobalConstants::Immediate = -1e8 that the immediate elimination did
        // not fold out -- so the controller is pinned near 1/1e8 while the
        // window runs to 10*iter/min_rate. 272 windows took 3.0 s between them
        // and the 273rd had not returned after 143 s.
        // A FINITE timespan is a transient request, which must reach its end time
        // rather than stop at the fixed point -- the same gate the geometric-tail
        // test below carries.
        const bool fp_armed = opt.earlystop && !std::isfinite(opt.timespan_end)
                              && !fluid_has_time_varying_rates(opt);
        if (fp_armed) {
            drift(t0, x.data(), drift_buf.data());
            double dn = 0.0, dtot = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                dn += std::fabs(drift_buf[i]);
                dtot += x[i];
            }
            if (dtot > 0.0 && dn / 2.0 / dtot / min_rate < lang::GlobalConstants::Zero) {
                t0 = t1;
                ++iter;
                break;
            }
        }
        // THE TEST ABOVE IS TAKEN ONCE, at the window's first instant. A window
        // that reaches the fixed point AFTER its first step is the same stall
        // entered one step later, and nothing above catches it, so the SAME test
        // rides on the integrator as a per-accepted-step stop. That is where
        // MATLAB's OutputFcn chain and the native-Python step loop take it, so
        // the four codebases stop on one condition at one place. `lsoda.h`
        // reproduces the output grid through `LsodaStepper` when this is set and
        // is untouched when it is not.
        // NOT ON THE TBI ARM. `tbi_advance` integrates one CELL at a time, on a
        // state vector holding only that cell's entries, while `drift` is the
        // WHOLE model's: handing it a cell-local vector reads and writes past
        // the end of both buffers. MATLAB arms the guard in
        // `solver_fluid_iteration.m` and Java in
        // `ClosingAndStateDepMethodsAnalyzer`, neither of which is the tbi arm,
        // so leaving tbi unarmed is also what keeps the four codebases aligned.
        // A stop test for tbi would have to be built from the CELL's drift,
        // inside `tbi_advance`, where the two state spaces agree.
        lopt.step_stop = {};
        if (fp_armed && m != "tbi") {
            const double fp_rate = min_rate;
            lopt.step_stop = [&drift, fp_rate](double tt, const std::vector<double>& yy) {
                if (yy.empty() || !(fp_rate > 0.0)) return false;
                std::vector<double> dy(yy.size(), 0.0);
                drift(tt, yy.data(), dy.data());
                double dn = 0.0, dtot = 0.0;
                for (std::size_t i = 0; i < yy.size(); ++i) {
                    dn += std::fabs(dy[i]);
                    dtot += yy[i];
                }
                return dtot > 0.0 && dn / 2.0 / dtot / fp_rate < lang::GlobalConstants::Zero;
            };
        }
        if (m == "tbi") {
            // Same drift, decomposed over cells; see fluid_tbi.h.
            x = tbi_advance(sys, tbi_cells, x, t0, t1, TbiOptions(), lopt);
        } else if (opt.stiff) {
            // `ode_solve_stiff.m`: same leg, same tolerances, implicit method
            // -- and the same per-restart tightening, since this arm is one leg
            // of the very loop that pays the local error iter_max times.
            FluidStiffOptions sopt;
            sopt.rtol = lopt.rtol;
            sopt.atol = lopt.atol;
            // The stiff arm is the Rosenbrock method, not LSODA, so it carries
            // the stop itself: this is the arm the settling windows run on.
            if (lopt.step_stop) {
                const std::function<bool(double, const std::vector<double>&)> ss_stop =
                    lopt.step_stop;
                sopt.step_stop = [ss_stop](const double& tt, const std::vector<double>& yy) {
                    return ss_stop(tt, yy);
                };
            }
            const OdeSolution<double> ss = fluid_ode_solve_stiff(drift, t0, t1, x, sopt);
            x = ss.final_state();
        } else {
            x = fluid_integrate_leg(drift, t0, t1, x, lopt);
        }
        // The drift conserves mass but the integrator need not, to the last
        // digit; a small negative mass is noise, so it is clamped rather than
        // allowed to feed back as a negative rate.
        for (double& v : x)
            if (v < 0.0) v = 0.0;
        // THAT CLAMP IS ALSO WHERE A DIVERGENCE HIDES. Under the moment closure
        // the drift can leave the simplex, and clamping the result turns a state
        // that is not a solution into one that merely looks like it settled --
        // every later window then integrates from it. The population of a closed
        // class is conserved EXACTLY by the drift, so a deviation is a
        // divergence and nothing else; raise the error the fallback ladder in
        // fluid_runner.h already catches, so 'dae' and then 'matrix'/'closing'
        // answer the model. Gated on gaussian() so the first-order pass, which
        // IS the ladder's own fallback, keeps identical behaviour.
        if (opt.closure.gaussian()) {
            const int bad = fluid_conservation_violation(sn, sys.layout, x);
            if (bad >= 0) {
                throw FluidNonHyperbolicError(
                    "The moment-closure drift left the model: closed chain " +
                    std::to_string(bad) + " moved more than " +
                    std::to_string(static_cast<int>(100 * kFluidConservationTol)) +
                    "% of a population the drift conserves exactly, by t = " +
                    std::to_string(t1) +
                    ", so the excursion is a divergence rather than a solution. "
                    "Falling back to a first-order closure.");
            }
        }
        t0 = t1;

        double moved = 0.0, total = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            moved += std::fabs(x[i] - prev[i]);
            total += prev[i];
        }
        const double ratio = (total > 0.0) ? moved / 2.0 / total : 0.0;
        // A FINITE timespan is a transient request, which must reach its end
        // time rather than stop at the fixed point; see the header.
        //
        // iter_tol = 0, the default, never fires here: the reference runs every
        // one of its iter_max passes and this port now does too.
        if (opt.iter_tol > 0.0 && ratio < opt.iter_tol && !std::isfinite(opt.timespan_end)) {
            ++iter;
            break;
        }
        // THE TERMINATION TEST. `ratio` alone is the mass moved over ONE window
        // and drops the geometric tail behind it; summing that tail,
        // ratio*rho/(1-rho), is what the header says it is missing. rho comes
        // from the iteration itself, so no rate has to stand in for the slowest
        // system mode -- when one does, the stop lands 3% short. The drift, zero
        // AT a fixed point, is an independent second bound. Both must hold on
        // two consecutive windows, past the slowest relaxation time.
        if (opt.earlystop && iter > 0 && !std::isfinite(opt.timespan_end) && t1 >= min_horizon) {
            rho_hist[iter % rho_hist.size()] =
                ratio / std::max(moved_prev, lang::GlobalConstants::Zero);
            double rho = 0.0;
            for (double v : rho_hist)
                if (std::isfinite(v) && v > rho) rho = v;
            drift(t1, x.data(), drift_buf.data());
            double dn = 0.0, dtot = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                dn += std::fabs(drift_buf[i]);
                dtot += x[i];
            }
            const double drift_displ = (dtot > 0.0) ? dn / 2.0 / dtot / min_rate : 0.0;
            // a non-contracting iteration has no tail to sum: it is not converging
            //
            // THE 1e-6 GATE IS NOT AN OVERSIGHT, even though it sits below the
            // integrator's own tol. Relaxing it to "the moved mass reached the
            // integrator floor, so trust the drift residual alone" was TRIED and
            // REVERTED: it stops the M/M/1 rho = 0.9 minnormal solve at 7.018088
            // against the 7.021524680 all four codebases agree on, and it truncates
            // the statedep response-time trajectory to t = 60 instead of its 2000.
            // The accuracy of this loop comes from running the windows, so the stop
            // has to stay conservative. It is also not what makes a solve hang: see
            // _kb/06-solver-catalog.md, where the minnormal closure diverges outright
            // on a bounded multiserver station.
            if (rho < 1.0 && ratio * rho / (1.0 - rho) < drift_safety * drift_tol
                && drift_displ < drift_tol) {
                if (++drift_below >= 2) {
                    ++iter;
                    break;
                }
            } else {
                drift_below = 0;
            }
        }
        moved_prev = ratio;
        if (t1 >= opt.timespan_end) {
            ++iter;
            break;
        }
    }

    // ---- read the metrics off the converged state -------------------------
    FluidSolution out;
    out.iters = iter;
    out.method = (statedep_family || m == "tbi") ? m : std::string("closing");
    out.xvec = x;
    out.QN = Matrix<double>(M, K, 0.0);
    out.UN = Matrix<double>(M, K, 0.0);
    out.RN = Matrix<double>(M, K, 0.0);
    out.TN = Matrix<double>(M, K, 0.0);

    fluid_closing_metrics(sn, sys, m, x, out.QN, out.UN, out.RN, out.TN);

    // THE COMPLETIONS AN ELIMINATED COORDINATE MAKES ARE NOT LOST. The metrics
    // above read throughputs off the STATE, as x_f * mu_f * phi_f summed over
    // phases, and an eliminated coordinate holds no mass there -- so its
    // completions, which are FINITE because mu_f is InfRate, would silently
    // vanish and the station would stop balancing against its neighbours. Their
    // total rate is exactly what `emap` carries: the composed event that replaced
    // the inflow stands for the original completion too.
    if (imm_result.eliminated) {
        const FluidLayout& LL = sys.layout;
        std::vector<std::size_t> cs(LL.nstates, 0), cc(LL.nstates, 0);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t k = 0; k < LL.kic[i][r]; ++k) {
                    cs[LL.qidx[i][r] + k] = i;
                    cc[LL.qidx[i][r] + k] = r;
                }
        std::vector<bool> kept(LL.nstates, false);
        for (std::size_t a = 0; a < imm_result.state_map.size(); ++a)
            kept[imm_result.state_map[a]] = true;
        std::vector<double> rr(x.begin(), x.end());
        fluid_rates_closing(dsys, x.data(), rr);
        for (std::size_t o = 0; o < sys.n_departures && o < sys.events.size(); ++o) {
            const std::size_t c = sys.events[o].event_idx;
            if (c >= LL.nstates || kept[c]) continue;
            double extra = 0.0;
            for (std::size_t e = 0; e < rr.size() && e < imm_result.emap.rows(); ++e)
                extra += imm_result.emap(e, o) * rr[e];
            out.TN(cs[c], cc[c]) += extra;
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (out.TN(i, r) > lang::GlobalConstants::Zero)
                    out.RN(i, r) = out.QN(i, r) / out.TN(i, r);
    }
    detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);

    // System throughput and response time, per chain reference station.
    out.XN.assign(K, 0.0);
    out.CN.assign(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        const std::size_t rs = sn.classes[r].refstat;
        if (rs >= 1 && rs <= M) out.XN[r] = out.TN(rs - 1, r);
        double q = 0.0;
        for (std::size_t i = 0; i < M; ++i) q += out.QN(i, r);
        if (out.XN[r] > 0.0) out.CN[r] = q / out.XN[r];
    }
    return out;
}

// ---------------------------------------------------------------------------
// The FCFS non-exponential refit loop of `solver_fluid_analyzer.m:100-197`.
// ---------------------------------------------------------------------------
/**
 * WHAT IT CORRECTS. The fluid drift of an FCFS station is the drift of a
 * PROCESSOR-SHARING station: a continuous mass has no queueing order to
 * respect, so every class in the buffer is served in proportion to its mass.
 * That is exact for exponential service, where the residual is memoryless and
 * the order does not matter, and wrong for anything else -- the whole point of
 * FCFS is that a long job blocks the ones behind it. The reference therefore
 * does not integrate an FCFS station at its declared service process. It runs
 * the mean-field solve, reads the resulting utilizations back into
 * `npfqn_nonexp_approx`, which rescales each FCFS station's service time by the
 * WSC 2020 diffusion interpolation, refits a COXIAN to that rescaled mean at the
 * station's ORIGINAL SCV, and integrates again. The loop closes on `eta`, the
 * M/G/1 decay rate the interpolation is built from.
 *
 * WHY THE FIT IS A COXIAN AND NOT A RATE CHANGE. `npfqn_nonexp_approx` returns a
 * scaled MEAN only. Writing that mean back as a one-phase exponential would
 * discard the SCV, which is the quantity that made the station non-product-form
 * in the first place; `Coxian.fitMeanAndSCV(1/rate, SCV)` keeps both moments and
 * changes the station's PHASE COUNT, which is why the layout, the initial
 * condition and the drift are all rebuilt inside the loop.
 *
 * THE SCV AND THE RATES ARE THE DECLARED ONES, EVERY SWEEP. `SCV = sn.scv` and
 * `rates0 = sn.rates` are read once, before the loop. Re-reading the SCV from
 * the refitted process would feed the fit its own output -- the Coxian written
 * in sweep n has, by construction, the SCV that was asked for -- so the sequence
 * would freeze at the first fit instead of converging to the interpolation's
 * fixed point. For the same reason the utilization fed back is `TN ./ rates0`
 * and not `TN ./ rates`, and `sn.rates` is never reassigned inside the loop.
 *
 * WHERE THE STATE HANDLING WENT. The reference re-encodes `sn.state` through
 * `State.toMarginal` / `State.fromMarginalAndStarted` whenever the phase count
 * changes, because `solver_fluid_initsol.m` DECODES that state to build the
 * initial condition. `fluid_default_initsol` is the closed form of that round
 * trip (see fluid_closing.h): it writes the `initDefault` placement into the
 * phase-one entries directly and never reads `sn.state`. The re-encoding is
 * therefore an identity here, and rebuilding the initial condition after a phase
 * change is the whole of its observable effect -- which is done.
 */
/**
 * The methods whose FCFS drift is the PS drift, and which the reference
 * therefore refits.
 *
 * The reference's second switch lists `matrix`, `closing`, `tbi`, `minnormal`,
 * `refined` and `dae`, and NOT `statedep`, `softmin`, `pnorm`, `diffusion`,
 * `mfq`, `rmf` or `kp`. The state-dependent family already carries a min()-based
 * capacity term, so its FCFS drift is not the PS one and refitting it would
 * correct a correction; `statedep` is commented in the reference as needing a
 * single iteration, and that comment is the contract. The rest are different
 * solvers, not different closures of the same drift.
 *
 * `dae` refits because it IS the min-normal closure -- same drift, same rate
 * factors, solved simultaneously instead of by substitution -- so the phase
 * count its answer depends on is refitted on the same schedule `minnormal` uses.
 */
inline bool fluid_method_refits_fcfs(const std::string& method) {
    std::string m = method;
    if (m.size() > 6 && m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    return m == "matrix" || m == "closing" || m == "tbi" || m == "minnormal" || m == "refined" ||
           m == "dae";
}

/** `cellsum(sn.visits)` at station level; the reference's `V` argument. */
template <class T>
Matrix<T> fluid_station_visits(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> V(sn.nstations, sn.nclasses, zero);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const std::size_t sf = sn.stateful_of_station(i + 1);
            for (std::size_t k = 0; k < sn.nclasses; ++k)
                V(i, k) = T(V(i, k) + sn.visits[c](sf - 1, k));
        }
    return V;
}

/**
 * MATLAB's max(abs(1 - eta ./ eta_1)) over a vector that MAY contain NaN.
 *
 * `max` SKIPS NaN in MATLAB, and an all-NaN vector makes the comparison
 * `[] > tol` false, i.e. the loop stops. A 0/0 entry -- a station whose decay
 * rate was zero on both sweeps -- must therefore not read as "not converged",
 * which is what a straight `std::max` over NaN would produce.
 */
inline double fluid_eta_gap(const std::vector<double>& eta, const std::vector<double>& eta_1) {
    double best = -std::numeric_limits<double>::infinity();
    bool any = false;
    for (std::size_t i = 0; i < eta.size(); ++i) {
        const double g = std::fabs(1.0 - eta[i] / eta_1[i]);
        if (std::isnan(g)) continue;
        any = true;
        best = std::max(best, g);
    }
    return any ? best : 0.0;
}

/** Elementwise reciprocal of A, with the reference's two sentinels for the
 *  degenerate entries. */
template <class T>
Matrix<T> fluid_reciprocal_guarded(const Matrix<T>& A) {
    const T one = num_traits<T>::from_int(1);
    Matrix<T> B(A.rows(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t r = 0; r < A.cols(); ++r) {
            const double a = num_traits<T>::to_double(A(i, r));
            if (std::isnan(a))
                B(i, r) = num_traits<T>::from_double(lang::GlobalConstants::FineTol);
            else if (a == 0.0 || std::isinf(1.0 / a))
                // A zero rate gives Inf, which the reference replaces by
                // GlobalConstants.Immediate = 1e8 -- a very LARGE service time,
                // not a vanishing one. That reads backwards and is what
                // `solver_fluid_analyzer.m:107-109` does; the pairs it applies to
                // carry no utilization, so the interpolation then leaves them alone.
                B(i, r) = num_traits<T>::from_double(lang::GlobalConstants::Immediate);
            else
                B(i, r) = T(one / A(i, r));
        }
    return B;
}

/**
 * Refit every FCFS station of `sn` to the rescaled service time at its declared
 * SCV, returning true when any station's PHASE COUNT changed.
 *
 * A phase change invalidates the fluid layout, so the caller must rebuild the
 * initial condition rather than carry the previous state vector across.
 */
template <class T>
bool fluid_refit_fcfs_stations(qn::NetworkStruct<T>& sn, const Matrix<T>& rates,
                               const Matrix<T>& SCV,
                               std::vector<std::vector<std::size_t>>& phases) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = sn.nstations, K = sn.nclasses;
    bool changed = false;
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].sched != lang::SchedStrategy::FCFS) continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (!(rates(i, r) > zero) || !(SCV(i, r) > zero)) continue;
            const pfqn::MarieCoxFit<T> cx = pfqn::marie_cox_fit(T(one / rates(i, r)), SCV(i, r));
            // `refresh_rates` is deliberately NOT called: the reference never
            // assigns `sn.rates` inside the loop, and every consumer below the
            // switch (the correction, the next sweep's rho) reads the DECLARED
            // rate. Only the process representation the drift is built from moves.
            sn.service[i][r] = lang::Distrib<T>::coxian(cx.mu, cx.phi);
            if (cx.mu.size() != phases[i][r]) changed = true;
            phases[i][r] = cx.mu.size();
        }
    }
    return changed;
}

/**
 * The loop. `seed` is the first integration, which the caller has already run at
 * the model's declared service processes, and `solve` re-integrates the refitted
 * struct by the SAME method -- UNCORRECTED, because the analyzer applies its
 * correction once, after the loop.
 *
 * Returns the uncorrected table of the final integration, and `iters` is that
 * integration's own pass count.
 *
 * `iters` IS NOT ACCUMULATED ACROSS THE SWEEPS, and that is the reference's
 * split rather than a simplification. `solver_fluid_analyzer.m` returns `iter`,
 * the number of REFIT sweeps, and keeps the summed integration count in
 * `outer_iters`, which it uses for runtime accounting only. This port's `iters`
 * has always meant "passes of the integration the reported table came from", so
 * summing the discarded sweeps into it would change what the field means for
 * every model, refitting or not. The sweep count is reported separately below.
 */
template <class T, class Solve>
FluidSolution fluid_fcfs_nonexp_refit(const qn::NetworkStruct<T>& sn0, const FluidOptions& opt,
                                      const FluidSolution& seed, Solve solve,
                                      qn::NetworkStruct<T>* sn_out = nullptr) {
    const std::size_t M = sn0.nstations, K = sn0.nclasses;
    // `result.solverSpecific.sn` of the reference: the struct the reported
    // solution was integrated on. It is the input one until the loop refits a
    // service process, and a caller that reads the state vector afterwards --
    // the passage time does -- needs the refitted one, whose phase counts the
    // vector is laid out by.
    if (sn_out) *sn_out = sn0;
    if (!fluid_method_refits_fcfs(opt.method)) return seed;
    bool any_fcfs = false;
    for (std::size_t i = 0; i < M; ++i)
        if (sn0.stations[i].sched == lang::SchedStrategy::FCFS) any_fcfs = true;
    if (!any_fcfs) return seed;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const Matrix<T>& rates0 = sn0.rates;
    const Matrix<T>& SCV = sn0.scv;
    const Matrix<T> V = fluid_station_visits(sn0);
    const Matrix<T> ST0 = fluid_reciprocal_guarded(rates0);

    std::vector<bool> isFCFS(M, false);
    std::vector<T> nservers(M, one), gamma(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        isFCFS[i] = sn0.stations[i].sched == lang::SchedStrategy::FCFS;
        nservers[i] = num_traits<T>::from_double(sn0.stations[i].nservers);
    }

    qn::NetworkStruct<T> sn = sn0;
    std::vector<std::vector<std::size_t>> phases(M, std::vector<std::size_t>(K, 0));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) phases[i][r] = sn0.service[i][r].D0.rows();

    FluidSolution cur = seed;
    std::vector<double> eta(M, std::numeric_limits<double>::infinity()), eta_1(M, 0.0);
    std::size_t iter = 0;

    while (fluid_eta_gap(eta, eta_1) > lang::GlobalConstants::CoarseTol && iter <= opt.iter_max) {
        ++iter;
        eta_1 = eta;

        Matrix<T> U(M, K, zero), TN(M, K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                TN(i, r) = num_traits<T>::from_double(cur.TN(i, r));
                if (rates0(i, r) > zero) U(i, r) = T(TN(i, r) / rates0(i, r));
            }

        const npfqn::NonexpApproxResult<T> na = npfqn::npfqn_nonexp_approx(
            opt.highvar, isFCFS, rates0, ST0, V, SCV, TN, U, gamma, nservers);
        gamma = na.gamma;
        for (std::size_t i = 0; i < M; ++i) eta[i] = num_traits<T>::to_double(na.eta[i]);

        const Matrix<T> rates = fluid_reciprocal_guarded(na.ST);
        const bool phase_change = fluid_refit_fcfs_stations(sn, rates, SCV, phases);

        FluidOptions o = opt;
        const std::vector<double> fresh = fluid_default_initsol(sn, fluid_layout(sn));
        o.init_sol = (!phase_change && cur.xvec.size() == fresh.size()) ? cur.xvec : fresh;
        cur = solve(sn, o);
    }

    // The reference re-solves once more from the CLEAN initial condition, so the
    // reported table is the drift of the converged service processes started
    // where the model says the system starts, not where the last sweep left off.
    FluidOptions o = opt;
    o.init_sol = fluid_default_initsol(sn, fluid_layout(sn));
    FluidSolution out = solve(sn, o);
    out.refit_sweeps = iter;
    if (sn_out) *sn_out = sn;
    return out;
}

}  // namespace detail


/**
 * Port of `solver_fluid_analyzer.m`: dispatch on the method, refit the
 * non-exponential FCFS stations the reference refits, then apply the
 * utilization and response-time correction it applies to whatever the branch
 * returned.
 */
template <class T>
FluidSolution solver_fluid(const qn::NetworkStruct<T>& sn_in, const FluidOptions& opt,
                          qn::NetworkStruct<T>* sn_out = nullptr) {
    // `@@SolverFLD/runAnalyzer.m:25` converts the non-Markovian service laws
    // first, and FORCES phfit = 'ph': the ODEs read mu*phi as a flow between
    // phases, and a matrix exponential has no such flow -- its off-diagonal
    // entries are not rates. The default two-moment CME would give a better
    // moment match and a meaningless drift.
    qn::NetworkStruct<T> converted;
    const qn::NetworkStruct<T>* snp = &sn_in;
    if constexpr (num_traits<T>::has_transcendental) {
        if (api::sn_has_nonmarkov(sn_in, false)) {
            converted = sn_in;
            api::NonmarkovOptions no;
            no.order = opt.nonmkv_order;
            no.phfit = api::PhFit::Ph;
            api::sn_nonmarkov_toph(converted, no);
            snp = &converted;
        }
    }
    const qn::NetworkStruct<T>& sn = *snp;

    FluidSolution out = detail::fluid_dispatch(sn, opt);
    out = detail::fluid_fcfs_nonexp_refit(
        sn, opt, out,
        [](const qn::NetworkStruct<T>& s, const FluidOptions& o) { return detail::fluid_dispatch(s, o); },
        sn_out);
    detail::fluid_analyzer_correct(sn, out.QN, out.UN, out.RN, out.TN);
    detail::fluid_snap_all(out.QN, out.UN, out.RN, out.TN);
    return out;
}


/**
 * Port of `local_detect_nhpp` in `@@SolverFLD/getTranAvg.m`: the (station, class)
 * pairs whose SOURCE carries a rate schedule, 1-based.
 *
 * A CALLER OPTS IN, and that is the reference's split rather than a convenience.
 * `getTranAvg` calls this and puts the result on the options; a steady-state
 * request does not, and is answered at the time-averaged nominal. Both are
 * legitimate readings of the same model, so the decision belongs to the entry
 * point and not to the drift builder.
 *
 * The reference tests `ismethod(proc,'getRateSchedule')`, which NHPP, MAPt and
 * PHt all answer; here that is `Distrib::has_schedule()`, the same three.
 */
template <class T>
std::vector<std::pair<std::size_t, std::size_t> > fluid_detect_nhpp(
    const qn::NetworkStruct<T>& sn) {
    std::vector<std::pair<std::size_t, std::size_t> > out;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched != lang::SchedStrategy::EXT) continue;
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            if (!sn.disabled[i][r] && sn.service[i][r].has_schedule())
                out.push_back(std::make_pair(i + 1, r + 1));
    }
    return out;
}

/**
 * Port of `@@SolverFLD/getTranAvg`: the metrics along the trajectory, not just
 * at the fixed point.
 *
 * The reference forces the method to `closing` for a transient (matrix and the
 * smoothed variants are steady-state devices), and so does this. The drift is
 * integrated once over [0, t_end] with the output grid handed to LSODA, and
 * every point is passed through the SAME extraction the steady state uses, so
 * the last point of a long enough run reproduces `solver_fluid` exactly.
 *
 * A transient is only meaningful from a KNOWN starting state, so the default
 * initial condition is used unless the caller supplies `init_sol`.
 *
 * THE RATE SCHEDULE IS DETECTED HERE, as `getTranAvg.m:76` detects it: a
 * transient of a model with a non-homogeneous source follows the intensity
 * exactly rather than its time average. A caller that has already filled
 * `opt.nhpp_sched` keeps its own list, so the nominal can still be asked for.
 *
 * `out_grid` REPLACES the uniform grid when a caller needs the trajectory at
 * points of its own choosing. It exists because interpolating a trajectory
 * cannot recover resolution it never had: SolverENV sums an exit average
 * against a sojourn CDF, and over a horizon of 1e3 read through an Exp(1) clock
 * a uniform 1001-point grid carries six samples where the whole weight lives.
 * LSODA takes an arbitrary increasing output vector, so asking for the points
 * that matter costs nothing and removes the interpolation entirely.
 */
template <class T>
std::vector<FluidTranPoint> solver_fluid_transient(const qn::NetworkStruct<T>& sn,
                                                   const FluidOptions& opt, double t_end,
                                                   std::size_t points = 101,
                                                   const std::vector<double>& out_grid =
                                                       std::vector<double>()) {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "solver_fluid_transient: the fluid drift is integrated by LSODA, which is double "
            "precision by construction; rerun with --arith double");
    if (!(t_end > 0.0)) throw InputError("solver_fluid_transient: t_end must be positive");
    if (points < 2) throw InputError("solver_fluid_transient: need at least two output points");

    FluidOptions o = opt;
    if (o.nhpp_sched.empty()) o.nhpp_sched = fluid_detect_nhpp(sn);
    FluidOdeSystem sys = fluid_ode_system(sn);
    sys.ratemult = detail::fluid_ratemult(sn, sys, o);
    const FluidLayout& L = sys.layout;
    if (L.nstates == 0)
        throw InputError("solver_fluid_transient: no station serves any class");

    std::vector<double> y0 =
        o.init_sol.empty() ? detail::fluid_default_initsol(sn, L) : o.init_sol;
    if (y0.size() != L.nstates)
        throw InputError("solver_fluid_transient: init_sol has the wrong length");

    std::vector<double> grid = out_grid;
    if (grid.empty()) {
        grid.resize(points);
        for (std::size_t j = 0; j < points; ++j)
            grid[j] = t_end * static_cast<double>(j) / static_cast<double>(points - 1);
    }

    LsodaOptions lopt;
    lopt.rtol = o.tol;
    lopt.atol = o.tol;
    const std::function<void(double, const double*, double*)> tdrift = fluid_drift(sys);
    const LsodaSolution sol = fluid_integrate_grid(tdrift, y0, grid, lopt);

    std::vector<FluidTranPoint> out;
    out.reserve(sol.y.size());
    for (std::size_t j = 0; j < sol.y.size(); ++j) {
        std::vector<double> xs = sol.y[j];
        for (double& v : xs)
            if (v < 0.0) v = 0.0;
        FluidTranPoint pt;
        pt.t = sol.t[j];
        Matrix<double> R;
        fluid_closing_metrics(sn, sys, std::string("closing"), xs, pt.QN, pt.UN, R, pt.TN);
        detail::fluid_snap_all(pt.QN, pt.UN, R, pt.TN);
        out.push_back(pt);
    }
    return out;
}

/**
 * The horizon a transient runs to when the caller gives none.
 *
 * FACTORED OUT OF `solver_fluid_tran_avg` because `dae` has a transient of its
 * own (fluid_dae.h) and must reach it by the SAME rule: a horizon rule that two
 * methods each computed for themselves is a horizon rule that can differ
 * between them, and then two trajectories of one model are read at different
 * times for no stated reason.
 *
 * `options.timespan` defaults to [0, Inf] and the reference does NOT integrate
 * to infinity for it. The rule lives in `@NetworkSolver/getTranAvg.m`, not in
 * the fluid analyzer: an unspecified end time becomes `30/minrate` with
 * `minrate = min(sn.rates(isfinite(sn.rates)))`, i.e. thirty mean events of the
 * SLOWEST RATE IN sn.rates -- arrival rates included, since the source's rate
 * sits in that same table. That is not the analyzer's own horizon-extension
 * rule (ten mean events of the slowest transition per pass), which governs how
 * far `solver_fluid` integrates while hunting the fixed point and is invisible
 * to the getter. A caller that sets `timespan_end` is integrated over exactly
 * that instead.
 *
 * MATLAB drops NaN (its disabled marker) through `isfinite`; the port carries a
 * separate `disabled` flag and stores zero, so the zero is skipped by the flag.
 */
template <class T>
double fluid_default_horizon(const qn::NetworkStruct<T>& sn, const FluidOptions& opt) {
    if (std::isfinite(opt.timespan_end) && opt.timespan_end > 0.0) return opt.timespan_end;
    double min_rate = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (sn.disabled[i][r]) continue;
            const double rate = num_traits<T>::to_double(sn.rates(i, r));
            if (std::isfinite(rate) && rate > opt.tol) min_rate = std::min(min_rate, rate);
        }
    if (!std::isfinite(min_rate)) min_rate = 1.0;
    return 30.0 / min_rate;
}

/**
 * `getTranAvg` on the first-order closing drift, over that horizon.
 *
 * THE METHOD IS NOT CONSULTED, as it is not in the reference: matrix, pnorm and
 * the smoothed variants are steady-state devices with no trajectory of their
 * own, so `getTranAvg.m` substitutes `closing` for them and warns. `dae` is the
 * exception the reference itself makes, and `solver_fluid_run_transient`
 * (fluid_runner.h) is where that routing lives -- it cannot live here, because
 * this header is below fluid_dae.h in the include order.
 */
template <class T>
std::vector<FluidTranPoint> solver_fluid_tran_avg(const qn::NetworkStruct<T>& sn,
                                                  const FluidOptions& opt,
                                                  std::size_t points = 101) {
    return solver_fluid_transient(sn, opt, fluid_default_horizon(sn, opt), points);
}

/**
 * The Jacobian of the fluid drift at a state, by central differences.
 *
 * The reference's `getJacobian` builds this SYMBOLICALLY and hands the
 * expression to a SAGE backend over HTTP; there is no symbolic engine here, so
 * this is the numerical counterpart. It is what the symbolic form is used for
 * in practice -- local stability of the fixed point, through the eigenvalues of
 * J -- and it needs no external service.
 */
template <class T>
Matrix<double> fluid_jacobian(const qn::NetworkStruct<T>& sn, const std::vector<double>& x) {
    const FluidOdeSystem sys = fluid_ode_system(sn);
    const std::size_t n = sys.layout.nstates;
    if (x.size() != n) throw InputError("fluid_jacobian: the state has the wrong length");
    const std::function<void(double, const double*, double*)> f = fluid_drift(sys);
    Matrix<double> J(n, n, 0.0);
    std::vector<double> xp(x), xm(x), fp(n, 0.0), fm(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        // A step scaled to the component, floored so a zero entry still moves.
        const double h = 1e-6 * std::max(1.0, std::fabs(x[j]));
        xp = x;
        xm = x;
        xp[j] += h;
        xm[j] -= h;
        f(0.0, xp.data(), fp.data());
        f(0.0, xm.data(), fm.data());
        for (std::size_t i = 0; i < n; ++i) J(i, j) = (fp[i] - fm[i]) / (2.0 * h);
    }
    return J;
}


/**
 * The joint probability of the per-class populations at station `i` (0-based)
 * under the linear noise approximation solved by the moment closure, the
 * `local_gaussian_cell` of the reference `@@SolverFLD/getProbAggr`.
 *
 * The state coordinates of class r at the station are `moments.class_block[i][r]`
 * (one per service phase), so the class population is their sum: its mean is the
 * reported `QN(i,r)` and the class-to-class covariance is the sum of the
 * corresponding block of `moments.Sigma`. The integer count n is then read off
 * the continuous law as the unit cell [n-1/2, n+1/2], with the two ends extended
 * to infinity at the boundaries of the state space, so that the mass the normal
 * puts on negative populations lands on the empty station and the mass above a
 * closed population lands on the full one. Those cells tile the state space, so
 * the probabilities sum to one over the reachable states.
 */
template <class T>
double fluid_prob_aggr_gaussian(const qn::NetworkStruct<T>& sn, const FluidSolution& sol,
                                std::size_t i, const std::vector<double>& nir,
                                double* logp_out = nullptr) {
    const std::size_t K = sn.nclasses;
    const std::vector<std::vector<std::size_t>>& cb = sol.moments.class_block[i];

    std::vector<std::size_t> idx;
    std::vector<double> m, a, b;
    for (std::size_t r = 0; r < K; ++r) {
        if (cb[r].empty()) {
            // the class has no service process here, so it has no coordinate:
            // any positive count is impossible rather than improbable
            if (nir[r] > 0.0) {
                if (logp_out) *logp_out = -std::numeric_limits<double>::infinity();
                return 0.0;
            }
            continue;
        }
        idx.push_back(r);
        m.push_back(sol.QN(i, r));
        a.push_back(nir[r] <= 0.0 ? -std::numeric_limits<double>::infinity() : nir[r] - 0.5);
        const double pop = sn.classes[r].population;
        b.push_back((std::isfinite(pop) && nir[r] >= pop) ? std::numeric_limits<double>::infinity()
                                                          : nir[r] + 0.5);
    }

    if (idx.empty()) {
        if (logp_out) *logp_out = 0.0;
        return 1.0;
    }

    const std::size_t nr = idx.size();
    Matrix<double> C(nr, nr, 0.0);
    for (std::size_t u = 0; u < nr; ++u)
        for (std::size_t v = u; v < nr; ++v) {
            double acc = 0.0;
            for (std::size_t p = 0; p < cb[idx[u]].size(); ++p)
                for (std::size_t q = 0; q < cb[idx[v]].size(); ++q)
                    acc += sol.moments.Sigma(cb[idx[u]][p], cb[idx[v]][q]);
            C(u, v) = acc;
            C(v, u) = acc;
        }

    const double p = fluid_mvn_rectangle(m, C, a, b);
    if (logp_out) *logp_out = p > 0.0 ? std::log(p) : -std::numeric_limits<double>::infinity();
    return p;
}

/**
 * Port of `@@SolverFLD/getProbAggr`: the probability that station `ist` holds
 * the marginal population of the model's default state.
 *
 * The fluid solver has no state space, so the probability is FITTED to the
 * mean queue lengths it does produce: a binomial for each closed class
 * (Schmidt) and, for the open ones, the BCMP marginal -- Poisson at an
 * infinite server, multinomial-geometric at a queue. The two contributions are
 * ADDED in log space, which is what lets a mixed model be evaluated at all;
 * the MVA port's version picks one branch because its callers are never mixed.
 *
 * `ist` is 1-based, as in the reference.
 */
template <class T>
double fluid_prob_aggr(const qn::NetworkStruct<T>& sn, const FluidSolution& sol, std::size_t ist,
                       double* logp_out = nullptr) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    if (ist == 0 || ist > M)
        throw InputError("fluid_prob_aggr: station number exceeds the number of stations");
    const std::size_t i = ist - 1;

    // The marginal of the DEFAULT state: a closed class sits at its reference
    // station, an open one holds nothing. This is what `State.toMarginal`
    // returns for the state the model starts in.
    std::vector<double> nir(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        const double pop = sn.classes[r].population;
        if (std::isfinite(pop) && sn.classes[r].refstat == ist) nir[r] = pop;
    }

    double logp = 0.0;
    bool minus_inf = false;
    const lang::SchedStrategy sc = sn.stations[i].sched;

    // The moment closure supplies the JOINT law of the per-class populations, so
    // the answer is the probability its multivariate normal assigns to the unit
    // cell around the state, correlation between the classes included. Three
    // exclusions, each structural rather than defensive:
    //   - a Source coordinate is a normalisation constant, not a population, and
    //     `fluid_moment_terms` projects it out of the covariance;
    //   - an OPEN class already has an EXACT first-order answer here (the BCMP
    //     marginal below), and a normal approximation of it would only lose: on
    //     M/M/1 at rho = 0.5 the cell returns 0.391 for the empty queue against
    //     an exact 0.500;
    //   - without moments there is no second moment anywhere in FLD.
    if (sol.has_moments && !sol.moments.class_block.empty() && sc != lang::SchedStrategy::EXT) {
        bool open_here = false;
        for (std::size_t r = 0; r < K; ++r)
            if (!sol.moments.class_block[i][r].empty() && !std::isfinite(sn.classes[r].population))
                open_here = true;
        if (!open_here)
            return fluid_prob_aggr_gaussian(sn, sol, i, nir, logp_out);
    }

    // ---- open classes ------------------------------------------------------
    bool any_open = false;
    for (std::size_t r = 0; r < K; ++r)
        if (!std::isfinite(sn.classes[r].population)) any_open = true;
    if (any_open && sc == lang::SchedStrategy::INF) {
        for (std::size_t r = 0; r < K; ++r) {
            if (std::isfinite(sn.classes[r].population)) continue;
            const double q = sol.QN(i, r);
            if (q > 0.0)
                logp += nir[r] * std::log(q) - q - std::lgamma(nir[r] + 1.0);
            else if (nir[r] > 0.0)
                minus_inf = true;
        }
    } else if (any_open && sc != lang::SchedStrategy::EXT) {
        double rho_total = 0.0, n_total = 0.0;
        for (std::size_t r = 0; r < K; ++r) {
            if (std::isfinite(sn.classes[r].population)) continue;
            rho_total += sol.UN(i, r);
            n_total += nir[r];
        }
        if (rho_total < 1.0) {
            logp += std::log(1.0 - rho_total) + std::lgamma(n_total + 1.0);
            for (std::size_t r = 0; r < K; ++r) {
                if (std::isfinite(sn.classes[r].population) || !(nir[r] > 0.0)) continue;
                const double rho_r = sol.UN(i, r);
                if (rho_r > 0.0)
                    logp += nir[r] * std::log(rho_r) - std::lgamma(nir[r] + 1.0);
                else
                    minus_inf = true;
            }
        } else {
            minus_inf = true;  // a saturated station has no stationary marginal
        }
    }

    // ---- closed classes: the Schmidt binomial ------------------------------
    for (std::size_t r = 0; r < K; ++r) {
        const double N = sn.classes[r].population;
        if (!std::isfinite(N)) continue;
        const double q = sol.QN(i, r);
        const double p = (N > 0.0) ? q / N : 0.0;
        // nchoosekln(N, nir)
        logp += std::lgamma(N + 1.0) - std::lgamma(nir[r] + 1.0) - std::lgamma(N - nir[r] + 1.0);
        if (p > 0.0) {
            logp += nir[r] * std::log(p);
        } else if (nir[r] > 0.0) {
            minus_inf = true;
        }
        if (p < 1.0) {
            logp += (N - nir[r]) * std::log(1.0 - p);
        } else if (N - nir[r] > 0.0) {
            minus_inf = true;
        }
    }

    if (minus_inf) {
        if (logp_out) *logp_out = -std::numeric_limits<double>::infinity();
        return 0.0;
    }
    if (logp_out) *logp_out = logp;
    return std::exp(logp);
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_SOLVER_FLUID_H

