/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * line-cli: multiprecision C++ front end, flag-compatible with
 * jar/src/main/java/jline/cli/LineCLI.java.
 *
 * The flag surface is the Java one plus --arith and --list-api. Anything not
 * yet ported is refused explicitly, naming what is missing; nothing is
 * silently approximated or answered from a partial implementation.
 *
 * ONE BINARY FOR BOTH MODEL KINDS, as the Java reference is: a Network model
 * (-i json) reaches solve_model_dispatch and a layered one (-i lqnx|xml)
 * reaches solve_lqn_dispatch. They were two executables until the LQN path was
 * folded in here; the split had no user-visible justification, since the Java
 * CLI has always taken both from one entry point.
 */
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "line/util/line_console.h"
#include "line/api/sn/sn_node_metrics.h"
#include "line/solvers/solver_node_tables.h"
#include "line/api/sn/sn_state.h"
#include "line/api/sym/sym_engines.h"
#include "line/io/docker_image.h"
#include "line/io/environment_reader.h"
#include "line/io/jsim_reader.h"
#include "line/io/network_reader.h"
#include "line/io/pnml.h"
#include "line/io/lqn_json_reader.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/num/number.h"
#include "line/reg/api_dispatch.h"
#include "line/reg/registry.h"
#include "line/util/method_type.h"
#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/auto/solver_auto.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/solver_default_cdf.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cdf.h"
#include "line/solvers/ctmc/solver_ctmc_cftp.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ctmc/solver_ctmc_mdd_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_reward.h"
#include "line/solvers/ctmc/solver_ctmc_sample.h"
#include "line/solvers/ctmc/solver_ctmc_sens.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/env/env_dispatch.h"
#include "line/solvers/fluid/fluid_jacobian.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/ldes/ldes_ln_engine.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/solvers/ag/ag_dispatch.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_prob.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_busyp.h"
#include "line/solvers/nc/solver_nc_cdf.h"
#include "line/solvers/nc/solver_nc_prob.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/sens/solver_sens_table.h"
#include "line/solvers/solver_chain_tables.h"
#include "line/solvers/ssa/solver_ssa_getters.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/uq/uq_dispatch.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"
#include "line/solvers/wrappers/jmt/jmt_logs.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/solvers/wrappers/lqns/solver_lqns.h"
#include "line/solvers/wrappers/qns/solver_qns.h"
#include "line/util/error.h"
#include "line/util/websocket.h"

namespace {

const char* kVersion = "0.1.0";

/**
 * Whether `-a avg` prints its table as JSON, i.e. `-o json`.
 *
 * A FILE-SCOPE FLAG, and the one place in this port that has one. The library
 * under `include/line/` keeps no global mutable state, deliberately, so a
 * pybind11 or MEX host can call it reentrantly; `line_cli.cpp` is a `main()`
 * translation unit and the output format is a property of ONE process
 * invocation, decided before any solve begins and never changed after. The
 * alternative is threading a bool through `solve_model_{mva,nc,mam,ba}`,
 * `solve_ctmc_avg`, `solve_fluid_avg` and `solve_ssa_avg`, three of which do
 * not take the `Knobs` struct at all, to carry a value none of them varies.
 */
bool g_json_output = false;

/**
 * Solve a Network model.json with SolverMVA and print its AvgTable.
 *
 * The columns and their order -- Station, JobClass, then QLen, Util, RespT,
 * ResidT, ArvR, Tput -- are the ones `parity/compare_parity.py` parses, so the
 * printed table is directly diffable against the MATLAB and Python rows. A
 * (station, class) pair with no presence at all (every metric zero) is dropped,
 * as `getAvgTable` drops its unvisited rows.
 */
/**
 * The solver knobs the CLI exposes, gathered so the dispatcher can check them
 * against the solver actually chosen.
 *
 * A sentinel means "not given" rather than a default, because the distinction
 * is what lets an option be REFUSED for a solver that has no such setting
 * instead of being accepted and dropped. That silent acceptance is the defect
 * this struct exists to prevent: `--samples 1e6` was previously taken without
 * complaint and the run still used 10000.
 */
struct Knobs {
    std::string method;
    double tol = -1.0;         // < 0 = not given
    double iter_tol = -1.0;
    int iter_max = -1;
    // --max-states, `options.config.maxStates`: the truncation level of an OPEN
    // agent's queue-length dimension in SolverAG. < 0 = not given, so AgOptions
    // keeps its own 100. It is a TRUNCATION and therefore part of the ANSWER,
    // not a budget: a run truncated at 100 that the caller asked to truncate at
    // 500 is a different number reported as theirs, which is why it travels
    // rather than being dropped the way the execution backends are.
    long long max_states = -1;
    // --multiserver, the AMVA rule that decides WHICH algorithm serves a
    // multiserver model. Empty = not given, so the solver keeps its default:
    // a rule this CLI silently dropped made every delegated solve answer under
    // 'default' while reporting the caller's choice.
    std::string multiserver;
    // --fork-join, `options.config.fork_join`: WHICH fork-join arm the shared
    // fixed point takes on a model with a Fork. Empty = not given, so the solver
    // keeps 'default' (the MMT transform); 'ht' is Heidelberger-Trivedi, which
    // is a different answer to the same model rather than a faster one.
    std::string fork_join;
    // The solver console has NO knob of its own: it IS VerboseLevel::DEBUG,
    // so `-v debug` is what asks for the running progress log.
    std::size_t samples = 0;   // 0 = not given
    unsigned long seed = 0;
    double cutoff = -1.0;      // < 0 = not given
    // --cutoff AS A MATRIX, `r1c1,r1c2;r2c1,r2c2`: the reference's
    // `options.cutoff` is a (station x class) table wherever a model needs a
    // different truncation per station, and reading only the scalar form left
    // `atof` silently taking the first number, i.e. truncating every station at
    // the first station's first class. Empty = not given.
    std::vector<std::vector<std::size_t>> cutoff_mat;
    /** True when `--cutoff` was given in EITHER spelling. */
    bool has_cutoff() const { return cutoff >= 0.0 || !cutoff_mat.empty(); }
    // --force, `options.force`: downgrade SolverCTMC's memory pre-gate from a
    // refusal to a warning. The gate exists because the alternative to refusing
    // is the OOM killer, so this is opt-in and never a default.
    bool force = false;
    // --fj-accuracy / --fj-tmode, `options.config.fj_accuracy` and
    // `options.config.fj_tmode` of solver_mam_fj.m. The first is the FJ_codes
    // truncation C of the queue-length difference between the two branches and
    // is the accuracy knob of that approximation; the second picks the route to
    // the T matrix. Kept apart from --tol/--iter_max for the same reason
    // --mdd-tol is: C is a STATE-SPACE bound, not a convergence threshold.
    int fj_accuracy = 0;       // 0 = not given
    std::string fj_tmode;      // empty = not given
    // --timescale, `options.config.timescale` of sn_is_discrete_time.m. "auto"
    // lets the distributions decide whether the model is slotted; "discrete"
    // and "continuous" force the reading, the first raising when the model
    // mixes lattice and non-lattice laws rather than solving the wrong time
    // scale. The slot itself comes from the SHARED `--slotlength` below, the
    // same one SolverNC's discrete product form reads.
    std::string timescale;     // empty = not given
    // --mdd-tol / --mdd-maxiter, the level iteration of `-s ctmc --method mdd`.
    // They are DELIBERATELY not --tol / --iter_tol: that iteration is an INNER
    // numerical solve whose fixed point is verified against the model's
    // population invariant at 1e-6, so an AMVA-sized tolerance converges short
    // of it and trips the guard. The reference keeps them apart for the same
    // reason (options.config.mdd_tol, not options.iter_tol).
    double mdd_tol = -1.0;     // < 0 = not given
    int mdd_maxiter = 0;       // 0 = not given
    // --level, the hierarchy level of the pbh/cbh/sib families and the
    // iteration count of pbk/bjbk. 0 = not given, so BaOptions keeps its 2.
    int level = 0;
    // --busyperiod / --busyperiod-subnet, the orders and the subnetwork of
    // `-a busyperiod`. The flag names are `ldes_cli`'s, so ONE spelling drives
    // the transform (solver_nc_busyp) and the simulation. The subnetwork has no
    // default -- a busy period is defined for a named set of stations and
    // choosing one here would answer about a subnetwork the caller never
    // named -- while the orders default to the ordinary busy period, 1.
    std::vector<std::size_t> busy_orders;
    std::vector<std::size_t> busy_subnet;
    // --qrf-params / --qrf-alpha, the blocking parameterisation and the
    // load-dependent scaling of the QRF arms of SolverBA. Both are JSON, given
    // inline or as a path; empty = not given.
    std::string qrf_params;
    std::string qrf_alpha;
    // --tspan, the horizon the transient CTMC analyses integrate over. There is
    // no default: pi(t) on an unstated horizon is not a quantity, and picking
    // one here would answer a question the caller did not ask.
    double t0 = 0.0, t1 = -1.0;  // t1 < 0 = not given
    std::size_t node = 0;        // --node, 1-based; 0 = not given
    // --class and --marg-states, the second and third arguments of
    // `@@SolverMVA/getProbMarg`. The class is 1-based and 0 = not given, i.e.
    // every class; the state list is the reference's `state_m` and empty = not
    // given, i.e. the default range each of the three laws picks for itself.
    // They are NOT folded into --node: a marginal is indexed by a PAIR, and one
    // flag carrying both would make "station 2" and "class 2" the same token.
    std::size_t jobclass = 0;       // --class, 1-based; 0 = not given
    std::vector<long> marg_states;  // --marg-states; empty = not given
    // --warmupfrac, `options.config.warmupfrac`: the leading fraction of an SSA
    // path discarded before the means are taken. The JAR CLI has carried it
    // since the SSA branch existed and this port did not, so a delegated solve
    // asking for a warmup discard silently kept the whole transient.
    // < 0 = not given, so the engine keeps its own 0.
    double warmupfrac = -1.0;
    // --pstar, `options.config.pstar`: the exponent of the fluid p-norm
    // smoothing (Ruuskanen et al., PEVA 151 (2021), eq. (26)). It selects the
    // DRIFT the matrix method integrates, so a solve that never receives it
    // returns the unsmoothed mean-field fixed point under the caller's choice.
    // < 0 = not given.
    double pstar = -1.0;
    // --notation, which form of the exported ODE document is wanted; empty = not
    // given, and only `-a odes` has one.
    std::string notation;
    // --cdf-algorithm, `options.config.algorithm` of `@@SolverNC/getCdfRespT`:
    // 'exact' is the pfqn_stdf sojourn-time inversion, 'rd' the pfqn_stdf_heur
    // reduction. Empty = not given, so the solver keeps the reference's 'exact'.
    // It is NOT --method: the ladder that computes the constants and the
    // algorithm that inverts the sojourn law are separate choices, and folding
    // them would make one name silently select the other.
    std::string cdf_algorithm;
    // -s ctmc -a firstpasst: the two state sets of `getCdfFirstPassT(A, B)`.
    // Each is either a 1-based index list into the state space ("3,5") or
    // semicolon-separated state rows ("0,2;1,1"), resolved against the space
    // the engine enumerates -- rows travel across the boundary because the two
    // enumerations need not order (or even purge) states identically.
    std::string passage_from;
    std::string passage_into;
    // "expm" (default) or "lt", `options.config.passage_method` of the reference
    std::string passage_method;
    // -s ctmc -a firstpasstmom: the highest moment order, `nmax` of
    // `getFirstPassTMoments(A, B, nmax)`. 0 means the reference's default of 3.
    std::size_t passage_orders = 0;
    // --perm-engine, the permanent estimator of `@@SolverNC/getProbSysMarg.m`:
    // 'exact' is Ryser's expansion with column multiplicities, and 'spm',
    // 'bethe', 'heur', 'huberlaw', 'adapart' are the five approximations. It is
    // NOT --method: the ladder that computes the normalizing constant and the
    // estimator that evaluates the permanent are separate choices. The five
    // approximations REFUSE a demand matrix with a structural zero rather than
    // flooring it, since they need full support. 'spm' is the only one that does
    // not expand the matrix to order sum(N), so it is the one whose cost does
    // not grow with the population and whose error falls as it grows.
    std::string method_perm = "exact";
    // --symbolic, `options.config.symbolic` of `@@SolverFLD/getJacobian`: auto
    // to search for a line-sage-rest backend, a URL, an image name, or none to
    // stay with the locally differentiated Jacobian. Empty = not given, i.e.
    // auto. --equilibria is the reference's fourth output, which is REQUESTED
    // and not implied: solving f(x) = 0 needs the backend, so implying it would
    // turn a Jacobian this port answers on its own into one that fails without
    // a container.
    std::string symbolic;
    bool equilibria = false;
    // ---- the layered path's own knobs, on the same not-given discipline ----
    // They are refused on the Network path rather than dropped, exactly as the
    // ones above are refused for a solver that has no such setting.
    bool no_interlocking = false;  // --no-interlocking was passed
    int repeat = 0;                // --repeat; 0 = not given, i.e. one run
    std::string layer_solver;      // --layer-solver; empty = not given
    /**
     * `--stage-solver`, which solver runs each stage of an ENVIRONMENT.
     *
     * IT IS NOT `--layer-solver`, and folding the two would be wrong: a LAYER of
     * an LQN is solved in STEADY STATE and a STAGE of a random environment is
     * solved TRANSIENTLY, so their admissible solver sets are different sets for
     * different reasons. Empty = not given, i.e. the coupling's own default
     * (fluid for the mean-field one, ctmc for the state-vector one).
     */
    std::string stage_solver;
    // --ln-transient / --ln-transient-channels, the coupling of the layered
    // transient and which inter-layer channels it injects. Empty = not given.
    std::string ln_transient;
    std::string ln_transient_channels;
    // --sens-method / --sens-scheme / --sens-step, the name-value contract of
    // getSensitivityTable. They are NOT --method: --method names the LN update
    // (default / moment3 / mwba.*), and the branch that differentiates it is a
    // separate choice. Empty / <= 0 = not given.
    std::string sens_method;
    std::string sens_scheme;
    double sens_step = -1.0;
    // --uq-solver, the engine SolverUQ runs at each design point. Empty is not
    // a default here but a missing argument: UQ computes nothing itself, so
    // there is no engine to fall back to, and `-s uq` without it is refused.
    std::string uq_solver;
    // --tran-points, the resolution of the uniform transient grid the ENV
    // mean-field coupling sums its exit metrics over. It is that coupling's
    // accuracy knob, not a cosmetic one: the sum is a Riemann-Stieltjes
    // quadrature against the holding-time CDF, so the answer moves with the
    // grid. 0 = not given, i.e. EnvOptions' own default.
    std::size_t tran_points = 0;
    // ---- the LQNS wrapper's own knobs -------------------------------------
    // `options.keep` and `options.verbose` of the reference wrapper, plus the
    // two that pick a REMOTE lqns. They apply to `-s lqns` alone and are
    // refused elsewhere: nothing else in this CLI runs a child process whose
    // working directory a caller might want to inspect.
    bool keep = false;
    bool verbose = false;
    bool remote = false;
    std::string remote_url;   // empty = not given, i.e. the wrapper's default
    int timeout_seconds = 0;  // 0 = not given, i.e. no deadline
    // ---- the simulator's own knobs, all `--ldes-*` -------------------------
    // They are PREFIXED rather than folded into the shared names because they
    // are an external engine's settings and not this port's: --ldes-tranfilter
    // is the warmup filter of a simulation run and has nothing to do with
    // --tol, and an unprefixed --warmupfrac would read as a knob every solver
    // has. Same not-given discipline as the rest -- empty, <= 0 or false means
    // not given -- so the engine keeps its own default and the command line
    // stays minimal, which is what an older AOT native binary can still parse.
    std::string ldes_tranfilter;  // mser5 | fixed | none
    double ldes_warmupfrac = -1.0;
    std::string ldes_cimethod;  // obm | bm | spectral | none
    bool ldes_cnvgon = false;
    double ldes_cnvgtol = -1.0;
    bool ldes_slotted = false;
    double ldes_slotlength = -1.0;
    /**
     * `--slotted` / `--slotlength`: the discrete time scale for the
     * ANALYTICAL solvers, distinct from the `--ldes-*` pair above, which
     * configures the simulator. SolverNC reads it and routes to the
     * discrete-time product form.
     */
    bool slotted = false;
    double slotlength = -1.0;
    int ldes_replications = 0;
    int ldes_numthreads = 0;
    double ldes_maxtime = -1.0;
    std::vector<double> ldes_initsol;  // station-major warm-start placement
    std::string ldes_rest_url;
    // ---- the JAR CLI's own five, ported so `line-cli` answers every question
    // `jline.cli.LineCLI` answers ------------------------------------------
    /**
     * `--state`, the state vector `-a prob` asks about.
     *
     * WITHOUT IT THE QUERY IS ABOUT THE MODEL'S DEFAULT INITIAL STATE, which is
     * what this CLI reported before and remains the default. The JAR takes an
     * explicit one because `getProb(node, state)` is a different question from
     * `getProb(node)`: the second names a state the model already holds, the
     * first names any state of the node's own space, and a caller sweeping a
     * marginal law needs the first. Empty = not given.
     */
    std::vector<long> state;
    /**
     * `--events`, the length of a sampled trajectory, `options.samples` of
     * `@@SolverSSA/sample`. It is NOT `--samples`: the JAR keeps them apart
     * because `--samples` is a simulation run length or a Monte Carlo draw
     * count and reaches a solver's options, while this is the number of EVENTS
     * one `sample` call walks. Folding them would make `-s ssa -a avg --events`
     * silently lengthen the run. 0 = not given, i.e. the reference's 1000.
     */
    std::size_t events = 0;
    /**
     * `--timestep`, the fixed output step of a transient analysis. Without it
     * the grid is whatever the integrator chose, which is the reference's
     * adaptive default; with it the trajectory is resampled onto a uniform
     * lattice of that step, which is what a caller diffing two transients needs.
     * <= 0 = not given.
     */
    double timestep = -1.0;
    /**
     * `--transient-method`, `options.config.transient_method` of
     * `solver_ctmc_transient_analyzer.m`: "ode" integrates the forward
     * equation, "fau" marches fast adaptive uniformization over the output
     * grid. Empty = not given, so the solver keeps the reference's "ode".
     *
     * It is NOT `--method`: the state-space path and the way the forward
     * equation is advanced on it are separate choices, and the reference keeps
     * this one out of its valid-method list because it changes no stationary
     * answer.
     */
    std::string transient_method;
    /**
     * `--fau-epsilon` and `--fau-delta`, the two tolerances of the "fau"
     * transient: the total probability mass the grid may discard, and the
     * occupancy below which a state is dropped from the support. <= 0 = not
     * given, i.e. the reference's 1e-6 and 1e-12.
     */
    double fau_epsilon = -1.0;
    double fau_delta = -1.0;
    /**
     * `--percentiles`, the levels `getPerctRespT` is read at. Empty = not
     * given, i.e. the reference's `pers_stored` {0.50, 0.90, 0.95, 0.99}.
     * Accepted as fractions (0.9) or as percents (90), told apart by magnitude:
     * a level above 1 cannot be a probability.
     */
    std::vector<double> percentiles;
    /**
     * `-v/--verbosity`: silent | standard | debug.
     *
     * IT WAS ACCEPTED AND DISCARDED, which is the silent-acceptance defect this
     * struct exists to prevent one layer up: a caller who asked for `silent`
     * still received every warning the arms print on stderr. It gates them now,
     * and nothing else -- the tables on stdout are the answer and are printed
     * whatever the level, exactly as the JAR prints them.
     */
    std::string verbosity;
    /**
     * `--reward-name`, which declared reward `-a reward-value` returns the
     * value function of. Empty = not given, and `-a reward-value` without it is
     * refused rather than defaulted to the first reward: the value functions of
     * two rewards are different objects, and picking one silently would label
     * the wrong matrix with the caller's question.
     */
    std::string reward_name;
};

/**
 * `--cutoff` written as a per-(station,class) matrix, `';'`-separated rows of
 * `','`-separated non-negative counts. Empty on anything that is not one.
 *
 * The spelling is the JAR CLI's, so one example pins one string for both
 * engines. A zero entry is legal and means the station may not hold that class
 * at all, which is how the reference bounds a queue in the classes it does not
 * serve.
 */
std::vector<std::vector<std::size_t>> parse_cutoff_matrix(const std::string& s) {
    std::vector<std::vector<std::size_t>> out;
    std::string::size_type pos = 0;
    while (pos <= s.size()) {
        const std::string::size_type semi = s.find(';', pos);
        const std::string row = s.substr(pos, semi == std::string::npos ? std::string::npos
                                                                       : semi - pos);
        std::vector<std::size_t> cells;
        std::string::size_type cp = 0;
        while (cp <= row.size()) {
            const std::string::size_type comma = row.find(',', cp);
            const std::string cell = row.substr(cp, comma == std::string::npos ? std::string::npos
                                                                              : comma - cp);
            if (cell.empty()) return std::vector<std::vector<std::size_t>>();
            for (std::string::size_type i = 0; i < cell.size(); ++i)
                if (!std::isdigit(static_cast<unsigned char>(cell[i])))
                    return std::vector<std::vector<std::size_t>>();
            cells.push_back(static_cast<std::size_t>(std::atol(cell.c_str())));
            if (comma == std::string::npos) break;
            cp = comma + 1;
        }
        if (cells.empty()) return std::vector<std::vector<std::size_t>>();
        if (!out.empty() && cells.size() != out[0].size())
            return std::vector<std::vector<std::size_t>>();
        out.push_back(cells);
        if (semi == std::string::npos) break;
        pos = semi + 1;
    }
    return out;
}

/**
 * The model text piped on stdin, drained ONCE and kept.
 *
 * `-s auto` parses the model twice: once to choose an engine and once for the
 * engine to solve. A stream can only be drained once, so without this buffer the
 * second parse would see nothing and the chooser would be unusable on a piped
 * model. It is deliberately not a function template -- a static local inside one
 * would be per-instantiation, and the two parses need not share an arithmetic.
 */
const std::string& stdin_model_text() {
    static std::string buf;
    static bool loaded = false;
    if (!loaded) {
        buf.assign(std::istreambuf_iterator<char>(std::cin), std::istreambuf_iterator<char>());
        loaded = true;
    }
    return buf;
}

/**
 * Whether `-i` selected a JMT document rather than a model.json.
 *
 * A FILE-SCOPE FLAG AND NOT A PARAMETER, because `read_model` is called from
 * every solver arm and threading the format through all of them would touch
 * fifty signatures to carry one bit that main already knows. It is set once,
 * before any arm runs, and never again.
 */
bool g_jsim_input = false;

/**
 * Whether `-i` selected a PNML place/transition net.
 *
 * A second flag rather than a format string for the same reason `g_jsim_input`
 * is one: `read_model` is called from every solver arm, and the three input
 * kinds it can be handed are mutually exclusive, so two booleans say what a
 * threaded enum would and touch no signature. Set once in main.
 */
bool g_pnml_input = false;

/**
 * `-v/--verbosity`, hoisted to file scope for the same reason as `g_jsim_input`.
 *
 * The knobs carry it too, but the priority warning below is raised where the
 * model is READ rather than inside a solver arm, and that function receives no
 * knobs. Set once in main, before any arm runs.
 */
std::string g_verbosity = "standard";

/**
 * Say so when class priorities were declared and no station will read them.
 *
 * The twin of the block at the tail of MATLAB's `@MNetwork/refreshStruct.m` and
 * of `Network.refreshStruct` in the JAR. It is a WARNING and not a refusal: a
 * priority at a station whose discipline ignores one is a legitimate model --
 * `prio_hol_open` is built on it -- and only becomes a defect when NO station
 * reads it, at which point the metrics are not the priority ones the caller is
 * about to read them as. Priority-awareness is a property of the declared
 * policy and is never inferred from the data; see network_struct.h.
 */
template <class T>
void warn_priorities_ignored(line::qn::Network<T>& net) {
    if (g_verbosity == "silent") return;
    // `raw_struct` and not `get_struct`: the priorities and the disciplines are
    // written when the classes and stations are added, and nothing in the
    // refresh chain touches either, so reading them here costs no refresh.
    if (!net.raw_struct().priorities_ignored()) return;
    std::fprintf(stderr,
                 "Warning: Priority classes are specified but no priority-aware scheduling "
                 "policy (PSPRIO, DPSPRIO, GPSPRIO, HOL, FCFSPRIO, FCFSPRPRIO, FCFSPIPRIO, "
                 "LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO, SRPTPRIO) is used in the model. "
                 "Priorities will be ignored.\n");
}

/** Read the model named by `file`, or stdin when it is empty. */
template <class T>
line::qn::Network<T> read_model(const std::string& file) {
    if (g_pnml_input) {
        // The PNML reader walks a DOM and has no stdin form, so a piped document
        // is staged to a temporary file, as the JSIM arm below does.
        if (!file.empty()) {
            line::qn::Network<T> net = line::io::pnml_load<T>(file);
            warn_priorities_ignored(net);
            return net;
        }
        const std::string tmp = line::io::jsim_stage_stdin(stdin_model_text());
        line::qn::Network<T> net = line::io::pnml_load<T>(tmp);
        std::remove(tmp.c_str());
        warn_priorities_ignored(net);
        return net;
    }
    if (g_jsim_input) {
        // The JSIM reader walks a DOM and has no stdin form, so a piped JMT
        // document is staged to a temporary file rather than refused: `cat
        // model.jsimg | line-cli -i jsimg` is how the JAR CLI is used and the
        // Docker image documents it.
        if (!file.empty()) {
            line::qn::Network<T> net = line::io::read_jsim<T>(file);
            warn_priorities_ignored(net);
            return net;
        }
        const std::string tmp = line::io::jsim_stage_stdin(stdin_model_text());
        line::qn::Network<T> net = line::io::read_jsim<T>(tmp);
        std::remove(tmp.c_str());
        warn_priorities_ignored(net);
        return net;
    }
    if (!file.empty()) {
        line::qn::Network<T> net = line::io::read_network_json<T>(file);
        warn_priorities_ignored(net);
        return net;
    }
    std::istringstream in(stdin_model_text());
    line::io::detail::json root;
    in >> root;
    line::qn::Network<T> net = line::io::build_network_from_json<T>(root);
    warn_priorities_ignored(net);
    return net;
}

/** Read the Environment model.json named by `file`, or stdin when it is empty. */
template <class T>
line::env::Environment<T> read_env_model(const std::string& file) {
    if (!file.empty()) return line::io::read_environment_json<T>(file);
    std::istringstream in(stdin_model_text());
    line::io::detail::json root;
    in >> root;
    return line::io::build_environment_from_json<T>(root);
}

/** One row of an average table, already reduced to double. */
struct AvgRow {
    double q, u, r, w, a, t;
};

/**
 * Render an average table, as text or as the JSON the hosts parse.
 *
 * THE ONE PLACE `-o json` IS HONOURED, reached by every solver arm through a
 * per-(station,class) accessor. The SSA and fluid arms used to hand-roll their
 * own printer -- their solution types are `SsaSolution`/`FluidSolution`, not
 * `mva::AvgResult<T>` -- and neither consulted `g_json_output`, so `-o json` was
 * ACCEPTED AND IGNORED for `-s ssa` and `-s fluid`: the caller got the readable
 * table with no brace in it, which is what the Python `lang='cpp'` bridge hit as
 * "no JSON object found in solver output". Accepting a flag and not applying it
 * is the silent-acceptance defect this CLI refuses everywhere else, so the
 * rendering is shared rather than reimplemented per solution type.
 *
 * THE JSON FORM IS THE JAR's, key for key. `jline.cli.LineCLI -a avg -o json`
 * emits {"avg": {"type": "AvgTable", "Station": [...], "JobClass": [...],
 * "QLen": [...], ...}}, column-oriented, and the Python wrapper's
 * `station_matrices_via_jar` parses exactly that shape. Emitting anything else
 * would force a second parser into the wrapper for a table that is the same
 * table, so the host can treat `lang='cpp'` and `lang='java'` as one transport
 * with two binaries. The JAR's `data` key (the rendered text) is NOT reproduced:
 * no caller reads it, and a second rendering of the same numbers is one more
 * thing that can disagree with the first.
 *
 * The rows are also byte-identical across the solvers that share it, which
 * matters more than it looks: `parity/compare_parity.py` diffs these tables
 * column by column, and a solver whose printer drifted by a space would read as
 * a parity failure in a harness that is supposed to be measuring the numbers.
 */
template <class Get>
void emit_avg_table_named(const std::vector<std::string>& stations,
                          const std::vector<std::string>& classes, const char* arith,
                          const std::string& method, Get get, const line::reg::Json& extra,
                          const line::reg::Json& envelope) {
    line::reg::Json rows = line::reg::Json::object();
    for (const char* key : {"Station", "JobClass", "QLen", "Util", "RespT", "ResidT", "ArvR",
                            "Tput"})
        rows[key] = line::reg::Json::array();
    rows["type"] = "AvgTable";

    if (!g_json_output)
        std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen",
                    "Util", "RespT", "ResidT", "ArvR", "Tput");
    for (std::size_t i = 0; i < stations.size(); ++i) {
        for (std::size_t c = 0; c < classes.size(); ++c) {
            const AvgRow v = get(i, c);
            if (v.q == 0.0 && v.u == 0.0 && v.r == 0.0 && v.w == 0.0 && v.a == 0.0 && v.t == 0.0)
                continue;
            if (!g_json_output) {
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                            stations[i].c_str(), classes[c].c_str(), v.q, v.u, v.r, v.w, v.a, v.t);
                continue;
            }
            // FULL PRECISION on the JSON path, against the readable path's
            // %12.6g. The table is for a human to read; the JSON is for a host
            // to compare against another codebase's answer, and six digits
            // would cap any parity check at six digits.
            rows["Station"].push_back(stations[i]);
            rows["JobClass"].push_back(classes[c]);
            rows["QLen"].push_back(v.q);
            rows["Util"].push_back(v.u);
            rows["RespT"].push_back(v.r);
            rows["ResidT"].push_back(v.w);
            rows["ArvR"].push_back(v.a);
            rows["Tput"].push_back(v.t);
        }
    }
    if (g_json_output) {
        // `extra` NESTS INSIDE the "avg" payload, for emit_analysis's reason: a
        // solver-specific key at envelope level can collide with an analysis's
        // own name, and a per-list cost qualifies THIS answer.
        for (line::reg::Json::const_iterator it = extra.begin(); it != extra.end(); ++it)
            rows[it.key()] = it.value();
        line::reg::Json out = line::reg::Json::object();
        out["avg"] = rows;
        out["arith"] = arith;
        out["method"] = method;
        // `envelope`, unlike `extra`, sits BESIDE "avg" rather than inside it.
        // Provenance of the solve as a whole belongs there: the iteration count
        // and the convergence flag qualify the answer, not any one row, and a
        // host reads them where the JAR CLI puts them.
        for (line::reg::Json::const_iterator it = envelope.begin(); it != envelope.end(); ++it)
            out[it.key()] = it.value();
        std::printf("%s\n", out.dump().c_str());
    }
}

/**
 * The same table, labelled from a struct.
 *
 * THE NAMES, NOT THE STRUCT, are what the printer needs, and one solver has no
 * struct to give it: `-s ldes` forwards the model document to an external engine
 * and is labelled from the names that engine reports, precisely so a model this
 * port cannot itself parse still prints the same table. Everything else calls
 * this overload and nothing about those rows changes.
 */
template <class T, class Get>
void emit_avg_table(const line::qn::NetworkStruct<T>& sn, const std::string& method, Get get,
                    const line::reg::Json& extra = line::reg::Json::object(),
                    const line::reg::Json& envelope = line::reg::Json::object()) {
    std::vector<std::string> stations, classes;
    stations.reserve(sn.nstations);
    classes.reserve(sn.nclasses);
    for (std::size_t i = 0; i < sn.nstations; ++i) stations.push_back(sn.stations[i].name);
    for (std::size_t c = 0; c < sn.nclasses; ++c) classes.push_back(sn.classes[c].name);
    emit_avg_table_named(stations, classes, line::num_traits<T>::name(), method, get, extra,
                         envelope);
}

/**
 * Print one analysis that is NOT the average table, as the JSON a host parses.
 *
 * ONE ENVELOPE FOR EVERY ANALYSIS: the payload sits under a key named after the
 * `-a` it answers, and the arithmetic and the resolved method sit beside it,
 * exactly as `emit_avg_table` places them beside "avg". A host therefore reads
 * the provenance the same way whatever it asked for, and a caller that asked for
 * `-a cdf` and received a "states" key knows the answer is not its own instead of
 * misreading the numbers as its own.
 *
 * `method` MAY BE EMPTY, in which case the key is omitted rather than filled
 * with the requested name: `-a reward` returns rewards and no solved chain, so
 * there is no resolved method to report, and echoing back "default" would claim
 * a resolution that never happened.
 *
 * EVERY INDEX IN A PAYLOAD IS 0-BASED, against the readable tables' 1-based
 * columns, and each payload carries `indexBase` so a host cannot get it wrong
 * silently. The two conventions are deliberate: the table is read by a human
 * diffing it against MATLAB, whose indices start at 1, while the JSON is
 * consumed by code that will index a numpy array or a std::vector with it.
 */
template <class T>
void emit_analysis(const char* key, const line::reg::Json& payload, const std::string& method,
                   const line::reg::Json& extra = line::reg::Json::object()) {
    line::reg::Json body = payload;
    // `extra` GOES INSIDE THE PAYLOAD, not beside it. At envelope level a
    // solver-specific key can collide with the analysis's own name -- the CTMC
    // state count is "states" and so is the `-a states` payload, and the merge
    // silently replaced the whole answer with the integer 3. Nesting it makes
    // that class of collision unrepresentable, and it is where the value belongs
    // anyway: a cutoff qualifies THIS answer.
    for (line::reg::Json::const_iterator it = extra.begin(); it != extra.end(); ++it)
        body[it.key()] = it.value();
    line::reg::Json out = line::reg::Json::object();
    out[key] = body;
    out["arith"] = line::num_traits<T>::name();
    if (!method.empty()) out["method"] = method;
    std::printf("%s\n", out.dump().c_str());
}

/** A Matrix as a row-major array of arrays, each entry reduced to double. */
template <class T>
line::reg::Json matrix_json(const line::Matrix<T>& M) {
    line::reg::Json rows = line::reg::Json::array();
    for (std::size_t i = 0; i < M.rows(); ++i) {
        line::reg::Json row = line::reg::Json::array();
        for (std::size_t j = 0; j < M.cols(); ++j)
            row.push_back(line::num_traits<T>::to_double(M(i, j)));
        rows.push_back(row);
    }
    return rows;
}

/** A vector of field elements as a JSON array of doubles. */
template <class T>
line::reg::Json vector_json(const std::vector<T>& v) {
    line::reg::Json a = line::reg::Json::array();
    for (std::size_t i = 0; i < v.size(); ++i) a.push_back(line::num_traits<T>::to_double(v[i]));
    return a;
}

/** A vector of sizes as a JSON array, unchanged: a width is not a measurement. */
line::reg::Json index_json(const std::vector<std::size_t>& v) {
    line::reg::Json a = line::reg::Json::array();
    for (std::size_t i = 0; i < v.size(); ++i) a.push_back(v[i]);
    return a;
}

/**
 * The per-Cache result block, as `-a avg` carries it inside the "avg" payload.
 *
 * PER-CACHE RESULTS RIDE WITH THE AVG TABLE, because a cache's hit, miss and
 * delayed-hit fractions are a SOLVER RESULT and a host that solves through this
 * CLI has no other way to get them back onto its own Cache node. Emitted per
 * node, and each vector is OMITTED when the solver computed none: absent must
 * CLEAR the host's copy, and a zero-filled vector would instead assert that
 * nothing hits.
 *
 * FACTORED OUT of `print_avg_table` because the SSA and Fluid arms build the
 * table themselves rather than through it, and so carried no block at all:
 * `CPPLINE.restoreCacheResults` then cleared MATLAB's Cache node, refreshed the
 * visits, and reported link()'s offered 1/2-1/2 for a split both engines had
 * measured. On cache_replc_rr that is hit 1.0 against 1.1246 (fluid) and
 * 1.1460 (ssa). A second copy of this serialization here would be free to drop
 * a field again, so there is one.
 */
template <class T>
line::reg::Json cache_extra_json(const line::solvers::CacheMetrics<T>& cache) {
    line::reg::Json caches = line::reg::Json::array();
    for (std::size_t c = 0; c < cache.caches.size(); ++c) {
        const line::solvers::CacheNodeMetrics<T>& m = cache.caches[c];
        line::reg::Json e = line::reg::Json::object();
        // NAME FIRST, because the index is not portable. `node` is an index
        // into THIS process's node order, which is not the model.json
        // declaration order: on retrieval_simple the JSON declares Source,
        // Cache, Queue, Sink and this struct holds Source, Queue, Sink,
        // Cache, so the Cache is 2 to the host and 4 here. A host matching
        // on the index wrote onto its Sink, found no Cache and silently kept
        // the PREVIOUS solver's numbers. `node` stays for provenance.
        e["name"] = m.name;
        e["node"] = m.node;
        if (!m.hitprob.empty()) e["HitProb"] = vector_json(m.hitprob);
        if (!m.missprob.empty()) e["MissProb"] = vector_json(m.missprob);
        if (!m.delayedprob.empty()) e["DelayedHitProb"] = vector_json(m.delayedprob);
        if (!m.latency.empty()) e["ResidT"] = vector_json(m.latency);
        if (!m.listcost.empty()) e["ListCost"] = vector_json(m.listcost);
        if (!m.delayedhitqlen.empty()) {
            e["DelayedHitQLen"] = vector_json(m.delayedhitqlen);
            e["DelayedHitQLenFull"] = vector_json(m.delayedhitqlenfull);
        }
        caches.push_back(e);
    }
    return caches;
}

/**
 * Print an AvgResult as the parity table.
 *
 * Shared by every solver that returns `mva::AvgResult` -- MVA, NC, MAM and BA
 * -- so the rows stay byte-identical across them.
 */
template <class T>
void print_avg_table(const line::qn::NetworkStruct<T>& sn, const line::mva::AvgResult<T>& r,
                     const line::reg::Json& extra_in = line::reg::Json::object()) {
    // THE JSON FORM IS THE JAR's, key for key. `jline.cli.LineCLI -a avg -o json`
    // emits {"avg": {"type": "AvgTable", "Station": [...], "JobClass": [...],
    // "QLen": [...], ...}}, column-oriented, and the Python wrapper's
    // `station_matrices_via_jar` parses exactly that shape. Emitting anything
    // else here would force a second parser into the wrapper for a table that is
    // the same table, so the host can treat `lang='cpp'` and `lang='java'` as one
    // transport with two binaries. The JAR's `data` key (the rendered text) is
    // NOT reproduced: no caller reads it, and a second rendering of the same
    // numbers is one more thing that can disagree with the first.
    // A CARRIED WARNING IS ONLY A WARNING IF SOMEONE SEES IT. `AvgResult`
    // carries the reference's text verbatim for the cases where the answer is
    // usable but not one the reference stands behind (the SJN starvation cap,
    // immediate feedback approximated as re-queueing), and until now no CLI
    // path printed it -- so the table looked authoritative exactly where the
    // reference declines. On STDERR, not stdout: stdout is the table a parity
    // harness diffs and the JSON a wrapper parses, and a warning line in either
    // would be read as data. Python's `lang='cpp'` bridge re-raises anything on
    // stderr as a Python warning, so it reaches that caller too.
    if (!r.warning.empty())
        std::fprintf(stderr, "Warning: %s\n", r.warning.c_str());

    // Mean per-list cache storage cost, the ListCost column of getAvgCacheTable.
    // Present only on a cache model carrying item sizes, and omitted entirely
    // otherwise rather than emitted empty, so a host can test for the key.
    line::reg::Json extra = line::reg::Json::object();
    if (!r.listcost.empty()) {
        line::reg::Json lc = line::reg::Json::array();
        for (std::size_t j = 0; j < r.listcost.size(); ++j)
            lc.push_back(line::num_traits<T>::to_double(r.listcost[j]));
        extra["ListCost"] = lc;
    }

    // WHAT THE CALLER ADDS, merged rather than replaced. The JMT arm carries the
    // finite-capacity-region rows this way: they are metric rows past the last
    // station, and the station table drops them, so without a channel of their
    // own a host reading `-a avg` cannot see a region at all (MATLAB's
    // `getAvgNodeTable` then filtered the FCR node out for having no numbers).
    if (extra_in.is_object())
        for (line::reg::Json::const_iterator it = extra_in.begin(); it != extra_in.end(); ++it)
            extra[it.key()] = it.value();

    // See `cache_extra_json`: the split is a solver result and the host has no
    // other channel back to its own Cache node.
    if (!r.cache.empty()) extra["Cache"] = cache_extra_json<T>(r.cache);

    // Provenance of the solve, beside "avg" and keyed as the JAR CLI keys it, so
    // one host-side reduction reads both backends. `converged` is OMITTED when
    // the handler reports none: absent is not false, and there the count is the
    // signal a caller may fall back on.
    line::reg::Json envelope = line::reg::Json::object();
    envelope["iter"] = r.iter;
    if (r.converged.has_value()) envelope["converged"] = r.converged.value();
    // `@SolverNC/getProbNormConstAggr`, keyed as the reference names the field it
    // stores it in. OMITTED for every solver that computes no constant: log G = 0
    // is the constant of an empty network, so a zero here would be a claim.
    if (r.lognormconst.has_value()) envelope["logNormConstAggr"] = r.lognormconst.value();

    emit_avg_table<T>(sn, r.actualmethod, [&](std::size_t i, std::size_t c) {
        AvgRow row;
        row.q = line::num_traits<T>::to_double(r.QN(i, c));
        row.u = line::num_traits<T>::to_double(r.UN(i, c));
        row.r = line::num_traits<T>::to_double(r.RN(i, c));
        row.w = line::num_traits<T>::to_double(r.WN(i, c));
        row.a = line::num_traits<T>::to_double(r.AN(i, c));
        row.t = line::num_traits<T>::to_double(r.TN(i, c));
        return row;
    }, extra, envelope);
}

template <class T>
int solve_model_mva(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::mva::MvaOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    if (!k.fork_join.empty()) opt.fork_join = k.fork_join;
    line::Matrix<T> init;
    const line::mva::AvgResult<T> r =
        line::mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    std::printf("SolverMVA arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                r.actualmethod.c_str(),
                line::util::method_type("MVA", r.actualmethod).c_str());
    print_avg_table<T>(sn, r);
    return 0;
}


/**
 * Solve a Network model.json with SolverJMT: write the JMT document, run the
 * external engine, print the same table every other solver prints.
 *
 * DOUBLE ONLY, and refused rather than relabelled under another arithmetic:
 * JMT simulates in double and reports in double, so an answer tagged
 * `real:128` would name a precision that never touched the computation. The
 * same rule the LDES arm applies, for the same reason.
 */
int solve_model_jmt(const std::string& file, const Knobs& k, const std::string& analysis) {
    line::qn::Network<double> net = read_model<double>(file);
    const line::qn::NetworkStruct<double>& sn = net.get_struct();

    line::jmt::JmtOptions o;
    if (!k.method.empty() && k.method != "default") o.method = k.method;
    if (k.samples > 0) o.samples = static_cast<double>(k.samples);
    if (k.seed != 0) o.seed = static_cast<long>(k.seed);
    o.keep = k.keep;
    if (k.t1 >= 0.0) o.max_simulated_time = k.t1;
    o.verbose = k.verbose;

    // `-a prob`: `getProbAggr` per station and `getProbSysAggr`, weighed off ONE
    // instrumented run. The DETAILED pair -- `getProb` and `getProbSys` -- has
    // no counterpart here and is OMITTED rather than aliased to the aggregate:
    // a JMT log records per-class job counts at a node and nothing about the
    // buffer order or the service phase, so the encoding those two are
    // probabilities of is not observed at all.
    if (analysis == "prob") {
        std::size_t target = 0;
        if (k.node) {
            if (k.node <= sn.nodes.size()) target = sn.nodes[k.node - 1].station;
            if (target == 0)
                throw line::InputError(
                    "--node " + std::to_string(k.node) +
                    " is not a station, so it holds no per-class job count to ask about");
        }
        const line::jmt::JmtProbAggr r = line::jmt::jmt_prob_aggr(
            sn, o, target, std::vector<double>(k.state.begin(), k.state.end()));
        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "ProbAggr";
            p["indexBase"] = 0;
            p["ProbSysAggr"] = r.sys;
            // WHETHER THE STATE OCCURRED AT ALL, beside the number. On an exact
            // solver a zero probability is a property of the model; on a
            // simulation it is far more often a property of the run length, and
            // a caller cannot tell the two apart from the zero alone.
            p["SysStateSeen"] = r.sys_seen;
            line::reg::Json st = line::reg::Json::array(), pa = line::reg::Json::array(),
                            sv = line::reg::Json::array();
            for (std::size_t i = 0; i < sn.nstations; ++i) {
                st.push_back(sn.stations[i].name);
                pa.push_back(r.station[i]);
                sv.push_back(static_cast<bool>(r.station_seen[i]));
            }
            p["Station"] = st;
            p["ProbAggr"] = pa;
            p["StateSeen"] = sv;
            // ALWAYS jsim, whatever `--method` asked for: this answer is read
            // off a simulated trajectory, and JMVA computes means from a
            // product form and logs nothing, so labelling it with the
            // requested method would name an engine that never ran.
            emit_analysis<double>("prob", p, std::string("jsim"));
            return 0;
        }
        std::printf("SolverJMT arith=double method=jsim\n");
        std::printf("ProbSysAggr %.10g%s\n", r.sys, r.sys_seen ? "" : "  (state never observed)");
        std::printf("%-16s %14s\n", "Station", "ProbAggr");
        for (std::size_t i = 0; i < sn.nstations; ++i)
            std::printf("%-16s %14.10g%s\n", sn.stations[i].name.c_str(), r.station[i],
                        r.station_seen[i] ? "" : "  (state never observed)");
        return 0;
    }

    if (analysis == "cdf" || analysis == "trancdf" || analysis == "trancdfpasst") {
        // -a cdf preloads the rounded steady-state queue lengths, the seeded
        // getCdfRespT pipeline; the transient names start from the default
        // initial state, getTranCdfRespT's contract
        const std::map<std::pair<std::size_t, std::size_t>,
                       std::vector<std::pair<double, double> > >
            rd = line::jmt::jmt_get_cdf_resp_t(sn, o, analysis == "cdf");
        line::reg::Json cdf = line::reg::Json::object();
        for (std::map<std::pair<std::size_t, std::size_t>,
                      std::vector<std::pair<double, double> > >::const_iterator it = rd.begin();
             it != rd.end(); ++it) {
            line::reg::Json rows = line::reg::Json::array();
            for (std::size_t i = 0; i < it->second.size(); ++i) {
                line::reg::Json row = line::reg::Json::array();
                row.push_back(it->second[i].first);
                row.push_back(it->second[i].second);
                rows.push_back(row);
            }
            const std::string key =
                sn.nodes[sn.station_to_node[it->first.first - 1] - 1].name + "/" +
                sn.classes[it->first.second - 1].name;
            cdf[key] = rows;
            if (!g_json_output)
                std::printf("%-24s %8zu points  respT(max)=%.6g\n", key.c_str(),
                            it->second.size(), it->second.back().second);
        }
        if (g_json_output) {
            line::reg::Json out = line::reg::Json::object();
            out["cdf"] = cdf;
            std::printf("%s\n", out.dump(2).c_str());
        }
        return 0;
    }

    const line::jmt::JmtResult<double> r = line::jmt::solver_jmt_run_analyzer(sn, o);
    std::printf("SolverJMT arith=double method=%s\n", r.avg.actualmethod.c_str());

    // The metric matrices carry `nregions` EXTRA rows past the stations; the
    // station table takes the first `nstations` of them and the regions are
    // reported separately, as the LDES arm reports its own FCR block.
    line::mva::AvgResult<double> table = r.avg;
    if (table.QN.rows() > sn.nstations) {
        const std::size_t M = sn.nstations, K = sn.nclasses;
        line::Matrix<double>* dst[6] = {&table.QN, &table.UN, &table.RN,
                                        &table.TN, &table.AN, &table.WN};
        const line::Matrix<double>* src[6] = {&r.avg.QN, &r.avg.UN, &r.avg.RN,
                                              &r.avg.TN, &r.avg.AN, &r.avg.WN};
        for (int m = 0; m < 6; ++m) {
            line::Matrix<double> t(M, K, 0.0);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < K; ++c) t(i, c) = (*src[m])(i, c);
            *dst[m] = t;
        }
    }
    // THE REGION ROWS TRAVEL BESIDE THE STATION TABLE, in the same shape the
    // node view prints them: a host solving through `-a avg` has no other way to
    // reach them, and MATLAB's `getAvgNodeTable` needs `result.Avg` to carry
    // `nstations + nregions` rows or it filters the FCR node out for having no
    // numbers at all (fcr_mm1waitq[M2C], 'row FCR1 missing').
    line::reg::Json fcr_extra = line::reg::Json::object();
    if (r.avg.QN.rows() > sn.nstations && sn.regions.size() > 0) {
        const std::size_t M = sn.nstations, K = sn.nclasses;
        line::reg::Json names = line::reg::Json::array();
        line::reg::Json q = line::reg::Json::array(), u = line::reg::Json::array();
        line::reg::Json rt = line::reg::Json::array(), w = line::reg::Json::array();
        line::reg::Json a = line::reg::Json::array(), t = line::reg::Json::array();
        for (std::size_t f = 0; f < sn.regions.size() && M + f < r.avg.QN.rows(); ++f) {
            names.push_back(f < sn.regions.size() && !sn.regions[f].name.empty()
                                ? sn.regions[f].name
                                : "FCR" + std::to_string(f + 1));
            for (std::size_t c = 0; c < K; ++c) {
                q.push_back(r.avg.QN(M + f, c));
                u.push_back(r.avg.UN(M + f, c));
                rt.push_back(r.avg.RN(M + f, c));
                w.push_back(r.avg.WN(M + f, c));
                a.push_back(r.avg.AN(M + f, c));
                t.push_back(r.avg.TN(M + f, c));
            }
        }
        fcr_extra["Region"] = names;
        fcr_extra["QLen"] = q;
        fcr_extra["Util"] = u;
        fcr_extra["RespT"] = rt;
        fcr_extra["ResidT"] = w;
        fcr_extra["ArvR"] = a;
        fcr_extra["Tput"] = t;
    }
    line::reg::Json avg_extra = line::reg::Json::object();
    if (!fcr_extra.empty()) avg_extra["FCR"] = fcr_extra;
    print_avg_table<double>(sn, table, avg_extra);

    if (!g_json_output && r.TNfcr.rows() > 0) {
        std::printf("%-16s %-14s %12s %12s\n", "Region", "JobClass", "Tput", "DropRate");
        for (std::size_t f = 0; f < r.TNfcr.rows(); ++f)
            for (std::size_t c = 0; c < sn.nclasses; ++c)
                // THE REGION'S DECLARED NAME, not the JSIM document's internal
                // `FCRegion<n>` label: that spelling is what the writer puts in
                // the exported model, and reporting it back renamed the user's
                // own region (`FCR1` in every other codebase's node table).
                std::printf("%-16s %-14s %12.6g %12.6g\n",
                            (f < sn.regions.size() && !sn.regions[f].name.empty()
                                 ? sn.regions[f].name
                                 : "FCR" + std::to_string(f + 1))
                                .c_str(),
                            sn.classes[c].name.c_str(), r.TNfcr(f, c), r.DropRateNfcr(f, c));
    }
    if (!g_json_output && !r.cache_hit_prob.empty())
        for (std::map<std::size_t, std::vector<double> >::const_iterator it =
                 r.cache_hit_prob.begin();
             it != r.cache_hit_prob.end(); ++it)
            for (std::size_t c = 0; c < it->second.size(); ++c)
                std::printf("%-16s %-14s hitProb=%12.6g\n", sn.nodes[it->first - 1].name.c_str(),
                            sn.classes[c].name.c_str(), it->second[c]);
    return 0;
}

/**
 * Solve a Network model.json with SolverNC and print the same table.
 *
 * `--samples` and `--seed` are wired because NC has stochastic methods that read
 * them ('mci', 'imci', 'ls', 'is', 'sampling' and 'mcmc'), and a run length nobody
 * can set is a method nobody can drive. `highvar` still has no flag and keeps its
 * SolverOptions('NC') default rather than being invented here.
 */
template <class T>
int solve_model_nc(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::nc::NcSolverOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (k.samples) opt.samples = k.samples;
    if (k.seed) opt.seed = k.seed;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    if (!k.fork_join.empty()) opt.fork_join = k.fork_join;
    if (k.slotted) opt.slotted = true;
    if (k.slotlength > 0.0) {
        opt.slotted = true;
        opt.slotlength = k.slotlength;
    }
    const line::mva::AvgResult<T> r = line::nc::solver_nc_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    // `lognormconst=` rides ON THE BANNER LINE rather than on one of its own,
    // as the fluid arm's `iters=` does: the banner is provenance and a new line
    // between it and the table is one more thing a table parser has to skip.
    // It is `getProbNormConstAggr`, which no other `-a` reports.
    std::printf("SolverNC arith=%s method=%s type=%s lognormconst=%.10g\n",
                line::num_traits<T>::name(), r.actualmethod.c_str(),
                line::util::method_type("NC", r.actualmethod).c_str(),
                r.lognormconst.has_value() ? r.lognormconst.value() : 0.0);
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * Solve a Network model.json with SolverMAM and print the same table.
 *
 * DOUBLE ONLY, and refused by name otherwise in the dispatcher: the analyzer
 * fits phase-type representations (aph_fit), which static_asserts on
 * transcendental arithmetic, so an exact instantiation would fail to COMPILE
 * rather than refuse at run time. `MamOptions`' tol, iter_max, space_max and
 * preserveDet keep their SolverOptions('MAM') defaults.
 */
template <class T>
int solve_model_mam(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::mam::MamOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    const line::mva::AvgResult<T> r = line::mam::solver_mam_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    std::printf("SolverMAM arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                r.actualmethod.c_str(),
                line::util::method_type("MAM", r.actualmethod).c_str());
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * Solve a Network model.json with SolverAG and print the same table.
 *
 * The RCAT/INAP arm, reachable from the CLI rather than API-only. It was the
 * one engine `run_avg_engine` already served for `-a node` and its three
 * sibling views while `-a avg` -- the view every parity row and every wrapper
 * asks for first -- refused the token outright, so `-s ag` reported an argument
 * error where the solver was present and working.
 *
 * `--method` picks between inap, inapplus, inapinf and exact; `default`
 * resolves to inap inside the analyzer, as it does in every codebase. `--tol`
 * and `--iter_max` are the fixed point's, matching -s mam.
 */
/**
 * The AgOptions a command line asks for.
 *
 * ONE PLACE, because there are three call sites (-a avg, -a cdf and the
 * run_avg_engine views) and a knob added to only some of them is a knob that
 * works or not depending on which view was asked for.
 */
void apply_ag_knobs(const Knobs& k, line::ag::AgOptions& opt) {
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (k.max_states > 0) opt.max_states = static_cast<std::size_t>(k.max_states);
}

template <class T>
int solve_model_ag(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ag::AgOptions opt;
    apply_ag_knobs(k, opt);
    const line::mva::AvgResult<T> r = line::ag::solver_ag_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    std::printf("SolverAG arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                r.actualmethod.c_str(),
                line::util::method_type("AG", r.actualmethod).c_str());
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * The station a per-node MAM query is about: `--node` when given, and otherwise
 * the model's only Queue.
 *
 * Defaulting is legitimate here and nowhere else: `getProb`, `getProbMarg` and
 * `getMAMResult` all run `require_single_queue`, so a model that reaches them
 * has exactly ONE queue and there is nothing to choose. `--node` is still
 * accepted, because naming the node one means is clearer than relying on that.
 */
template <class T>
std::size_t mam_query_node(const line::qn::NetworkStruct<T>& sn, const Knobs& k) {
    if (k.node) return k.node;
    for (std::size_t a = 0; a < sn.nof_nodes(); ++a)
        if (sn.nodes[a].nodetype == line::qn::NodeType::Queue) return a + 1;
    throw line::InputError(
        "the MAM per-node analyses report a queue's internals and this model has no Queue node");
}

/** The MamOptions every `-s mam` entry point builds from the CLI knobs. */
line::mam::MamOptions mam_options(const Knobs& k) {
    line::mam::MamOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (k.cutoff >= 0.0) opt.cutoff = static_cast<std::size_t>(k.cutoff);
    if (k.fj_accuracy > 0) opt.fj_accuracy = static_cast<std::size_t>(k.fj_accuracy);
    if (!k.fj_tmode.empty()) opt.fj_tmode = k.fj_tmode;
    if (!k.timescale.empty()) opt.timescale = k.timescale;
    if (k.slotlength > 0.0) opt.slotlength = k.slotlength;
    if (k.t1 >= 0.0) {
        opt.timespan_start = k.t0;
        opt.timespan_end = k.t1;
    }
    return opt;
}

/**
 * `-s mam -a prob`: `@@SolverMAM/getProb` and `@@SolverMAM/getProbMarg`.
 *
 * BOTH, in one answer, because they are two views of the same queue-length law:
 * the joint (level, phase) table the first returns, and the per-class marginal
 * P(n jobs of class r) the second does. Reporting only one would leave the other
 * unreachable again, which is the state this wiring closes.
 *
 * `--cutoff` is the level truncation an OPEN model needs -- its queue length is
 * unbounded, so the table has to stop somewhere -- and is passed straight
 * through as `options.cutoff`; a closed model bounds itself by its population
 * and ignores it.
 */
template <class T>
int solve_model_mam_prob(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::mam::MamOptions opt = mam_options(k);
    // The reference's getters read `self.getAvg()` first, so the solve comes
    // before the query and a caller cannot reach them on an unsolved model.
    const line::mva::AvgResult<T> avg = line::mam::solver_mam_run_analyzer(sn, opt);
    const std::size_t node = mam_query_node<T>(sn, k);
    const std::size_t ist = sn.nodes[node - 1].station;
    const line::mam::ProbTable<T> P = line::mam::solver_mam_get_prob(sn, opt, node, avg);
    std::vector<std::vector<T> > marg(sn.nclasses);
    for (std::size_t r = 0; r < sn.nclasses; ++r)
        marg[r] = line::mam::solver_mam_get_prob_marg(sn, opt, ist, r + 1, avg);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbTable";
        p["indexBase"] = 0;
        p["node"] = node - 1;
        p["Node"] = sn.nodes[node - 1].name;
        p["levels"] = P.P.rows();
        p["phases"] = P.P.cols();
        line::reg::Json joint = line::reg::Json::array();
        for (std::size_t n = 0; n < P.P.rows(); ++n) {
            line::reg::Json row = line::reg::Json::array();
            for (std::size_t j = 0; j < P.P.cols(); ++j)
                row.push_back(line::num_traits<T>::to_double(P.P(n, j)));
            joint.push_back(row);
        }
        p["joint"] = joint;
        line::reg::Json mj = line::reg::Json::array();
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            line::reg::Json e = line::reg::Json::object();
            e["JobClass"] = sn.classes[r].name;
            e["jobclass"] = r;
            e["P"] = vector_json(marg[r]);
            mj.push_back(e);
        }
        p["marginal"] = mj;
        emit_analysis<T>("prob", p, avg.actualmethod);
        return 0;
    }
    std::printf("SolverMAM arith=%s method=%s node=%s levels=%zu phases=%zu\n",
                line::num_traits<T>::name(), avg.actualmethod.c_str(),
                sn.nodes[node - 1].name.c_str(), P.P.rows(), P.P.cols());
    std::printf("%-8s %-8s %16s\n", "Level", "Phase", "Prob");
    for (std::size_t n = 0; n < P.P.rows(); ++n)
        for (std::size_t j = 0; j < P.P.cols(); ++j)
            std::printf("%-8zu %-8zu %16.10g\n", n, j + 1,
                        line::num_traits<T>::to_double(P.P(n, j)));
    std::printf("%-14s %-8s %16s\n", "JobClass", "Jobs", "Prob");
    for (std::size_t r = 0; r < sn.nclasses; ++r)
        for (std::size_t n = 0; n < marg[r].size(); ++n)
            std::printf("%-14s %-8zu %16.10g\n", sn.classes[r].name.c_str(), n,
                        line::num_traits<T>::to_double(marg[r][n]));
    return 0;
}

/**
 * The levels `getPerctRespT` is read at: `--percentiles`, or the reference's
 * `pers_stored` when the caller named none.
 */
std::vector<double> percentile_levels(const Knobs& k) {
    if (!k.percentiles.empty()) return k.percentiles;
    std::vector<double> pcts;
    pcts.push_back(0.50);
    pcts.push_back(0.90);
    pcts.push_back(0.95);
    pcts.push_back(0.99);
    return pcts;
}

/**
 * `-s mam -a cdf`: `@@SolverMAM/getCdfRespT` (and its aliases getSjrnT / sjrnT),
 * with `@@SolverMAM/getPerctRespT` beside it.
 *
 * THE PERCENTILE LEVELS COME FROM `--percentiles`, defaulting to the
 * {0.50, 0.90, 0.95, 0.99} that `solver_mam_fj.m` stores as `pers_stored`: the
 * reference takes them as an argument to `getPerctRespT`, and the flag is that
 * argument. On the CDF path they are a READING of the same curve --
 * `mam_percentiles_from_cdf` inverts the CDF that is printed above them -- so a
 * level the caller names costs nothing beyond the inversion.
 *
 * A FORK-JOIN MODEL HAS NO CDF HERE, and that is the reference's structure
 * rather than a gap in this port: `getPerctRespT.m` reads the table
 * `solver_mam_fj.m` stored, and its own fallback comment records that
 * `getCdfRespT` is "not available for this model". The two are separate methods
 * in MATLAB and only this CLI bundles them, so the fork-join case prints the
 * percentiles alone and says why. The condition is tested explicitly, not
 * discovered by catching the CDF's refusal: a caught exception cannot tell "no
 * fork-join route exists" from "the passage-time engine failed on this model".
 */
template <class T>
int solve_model_mam_cdf(const std::string& file, const Knobs& k, const char* key,
                        const char* type) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::mam::MamOptions opt = mam_options(k);
    const std::vector<double> pcts = percentile_levels(k);

    if (line::mam::mam_has_fj_percentiles(sn, opt)) {
        const std::vector<std::vector<T> > perc =
            line::mam::solver_mam_get_perct_respt(sn, opt, pcts);
        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "PerctRespT";
            p["indexBase"] = 0;
            p["source"] = "fjcodes";
            line::reg::Json arr = line::reg::Json::array();
            for (std::size_t r = 0; r < perc.size(); ++r) {
                line::reg::Json e = line::reg::Json::object();
                e["JobClass"] = sn.classes[r].name;
                e["jobclass"] = r;
                e["percentileLevels"] = pcts;
                e["percentiles"] = vector_json(perc[r]);
                arr.push_back(e);
            }
            p["respt"] = arr;
            emit_analysis<T>(key, p, opt.method);
            return 0;
        }
        std::printf("SolverMAM arith=%s method=%s classes=%zu\n", line::num_traits<T>::name(),
                    opt.method.c_str(), sn.nclasses);
        std::printf("# the response-time CDF has no fork-join route in the reference; these are "
                    "the FJ_codes percentiles getPerctRespT reads\n");
        std::printf("%-14s %14s %14s\n", "JobClass", "Percentile", "RespT");
        for (std::size_t r = 0; r < perc.size(); ++r)
            for (std::size_t j = 0; j < perc[r].size(); ++j)
                std::printf("%-14s %14.4g %14.10g\n", sn.classes[r].name.c_str(), pcts[j],
                            line::num_traits<T>::to_double(perc[r][j]));
        return 0;
    }

    const std::vector<line::mam::RespTCdf<T> > rd = line::mam::solver_mam_get_cdf_respt(sn, opt);
    std::vector<std::vector<T> > perc;
    for (std::size_t r = 0; r < rd.size(); ++r)
        perc.push_back(line::mam::mam_percentiles_from_cdf(rd[r], pcts));

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = type;
        p["indexBase"] = 0;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t r = 0; r < rd.size(); ++r) {
            // An empty curve is not a degenerate one: the class has no passage
            // through a queue here, and it is OMITTED rather than sent as a flat
            // zero law.
            if (rd[r].X.empty()) continue;
            line::reg::Json e = line::reg::Json::object();
            e["JobClass"] = sn.classes[r].name;
            e["jobclass"] = r;
            e["t"] = vector_json(rd[r].X);
            e["F"] = vector_json(rd[r].F);
            e["percentileLevels"] = pcts;
            e["percentiles"] = vector_json(perc[r]);
            arr.push_back(e);
        }
        p["respt"] = arr;
        emit_analysis<T>(key, p, opt.method);
        return 0;
    }
    std::printf("SolverMAM arith=%s method=%s classes=%zu\n", line::num_traits<T>::name(),
                opt.method.c_str(), sn.nclasses);
    std::printf("%-14s %14s %14s\n", "JobClass", "Time", "F(t)");
    for (std::size_t r = 0; r < rd.size(); ++r)
        for (std::size_t j = 0; j < rd[r].X.size(); ++j)
            std::printf("%-14s %14.8g %14.10g\n", sn.classes[r].name.c_str(),
                        line::num_traits<T>::to_double(rd[r].X[j]),
                        line::num_traits<T>::to_double(rd[r].F[j]));
    std::printf("%-14s %14s %14s\n", "JobClass", "Percentile", "RespT");
    for (std::size_t r = 0; r < perc.size(); ++r)
        for (std::size_t j = 0; j < perc[r].size(); ++j)
            std::printf("%-14s %14.4g %14.10g\n", sn.classes[r].name.c_str(), pcts[j],
                        line::num_traits<T>::to_double(perc[r][j]));
    return 0;
}

/**
 * `-s mam -a perct-respt`: `@@SolverMAM/getPerctRespT` ALONE.
 *
 * The percentiles ride beside the curve under `-a cdf` as well, and this arm
 * exists because they are a separate METHOD in the reference and a separate
 * answer to a caller: a service-level question ("what is the 99th percentile")
 * wants four numbers, not the whole law printed above them. On a fork-join
 * model it is the only route -- `getCdfRespT` has none there -- and on every
 * other model the levels are inverted from the CDF the other arm prints, so
 * the two never disagree.
 */
template <class T>
int solve_model_mam_perct(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::mam::MamOptions opt = mam_options(k);
    const std::vector<double> pcts = percentile_levels(k);

    std::vector<std::vector<T> > perc;
    const bool fj = line::mam::mam_has_fj_percentiles(sn, opt);
    if (fj) {
        perc = line::mam::solver_mam_get_perct_respt(sn, opt, pcts);
    } else {
        const std::vector<line::mam::RespTCdf<T> > rd =
            line::mam::solver_mam_get_cdf_respt(sn, opt);
        for (std::size_t r = 0; r < rd.size(); ++r)
            perc.push_back(line::mam::mam_percentiles_from_cdf(rd[r], pcts));
    }

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "PerctRespT";
        p["indexBase"] = 0;
        p["source"] = fj ? "fjcodes" : "cdf";
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t r = 0; r < perc.size(); ++r) {
            if (perc[r].empty()) continue;
            line::reg::Json e = line::reg::Json::object();
            e["JobClass"] = sn.classes[r].name;
            e["jobclass"] = r;
            e["percentileLevels"] = pcts;
            e["percentiles"] = vector_json(perc[r]);
            arr.push_back(e);
        }
        p["respt"] = arr;
        emit_analysis<T>("perct", p, opt.method);
        return 0;
    }
    std::printf("SolverMAM arith=%s method=%s classes=%zu source=%s\n",
                line::num_traits<T>::name(), opt.method.c_str(), sn.nclasses,
                fj ? "fjcodes" : "cdf");
    std::printf("%-14s %14s %14s\n", "JobClass", "Percentile", "RespT");
    for (std::size_t r = 0; r < perc.size(); ++r)
        for (std::size_t j = 0; j < perc[r].size(); ++j)
            std::printf("%-14s %14.4g %14.10g\n", sn.classes[r].name.c_str(), pcts[j],
                        line::num_traits<T>::to_double(perc[r][j]));
    return 0;
}

/**
 * `-s mam -a tran`: `@@SolverMAM/getTranAvg`.
 *
 * The reference FORCES `options.method = 'ldqbd'` before delegating, so the
 * transient path is not the caller's method and `--method` does not select it;
 * which engine runs -- the Laplace-domain transient QBD or the QBD fast path --
 * is decided from the model by `mam_transient_qbd_applicable`. `--tspan` is
 * required: a transient curve on an unstated horizon is not a quantity.
 */
template <class T>
int solve_model_mam_tran(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::mam::MamOptions opt = mam_options(k);
    const line::mam::TranResult<T> tr = line::mam::solver_mam_get_tran_avg(sn, opt);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranAvgTable";
        p["indexBase"] = 0;
        p["t0"] = opt.timespan_start;
        p["t1"] = opt.timespan_end;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < tr.Qt.size(); ++i)
            for (std::size_t r = 0; r < tr.Qt[i].size(); ++r) {
                if (tr.Qt[i][r].times.empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[r].name;
                e["station"] = i;
                e["jobclass"] = r;
                e["t"] = tr.Qt[i][r].times;
                e["QLen"] = vector_json(tr.Qt[i][r].values);
                e["Util"] = vector_json(tr.Ut[i][r].values);
                e["Tput"] = vector_json(tr.Tt[i][r].values);
                arr.push_back(e);
            }
        p["curves"] = arr;
        emit_analysis<T>("tran", p, "ldqbd");
        return 0;
    }
    std::printf("SolverMAM arith=%s method=ldqbd tspan=[%g,%g]\n", line::num_traits<T>::name(),
                opt.timespan_start, opt.timespan_end);
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "Time", "QLen", "Util",
                "Tput");
    for (std::size_t i = 0; i < tr.Qt.size(); ++i)
        for (std::size_t r = 0; r < tr.Qt[i].size(); ++r)
            for (std::size_t j = 0; j < tr.Qt[i][r].times.size(); ++j)
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n",
                            sn.stations[i].name.c_str(), sn.classes[r].name.c_str(),
                            tr.Qt[i][r].times[j],
                            line::num_traits<T>::to_double(tr.Qt[i][r].values[j]),
                            line::num_traits<T>::to_double(tr.Ut[i][r].values[j]),
                            line::num_traits<T>::to_double(tr.Tt[i][r].values[j]));
    return 0;
}

/**
 * `-s mam -a internals`: `@@SolverMAM/getMAMResult`, the M/G/1-type internals of
 * a single queue.
 *
 * It is NOT a metric table and is deliberately not folded into `-a avg`: the
 * answer is the matrix-analytic machinery itself -- the randomized blocks, the
 * G matrix, the drift, the decay rate, the level probabilities -- which is what
 * a caller checking a queue's stability or tail decay asks for, and what a
 * cross-codebase comparison of the QBD assembly needs to see.
 */
template <class T>
int solve_model_mam_internals(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    (void)k;
    const line::qsys::BmapM1Result<T> r = line::mam::solver_mam_get_mam_result(sn);
    const double q = line::num_traits<T>::to_double(r.q);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "MAMResult";
        p["indexBase"] = 0;
        p["lambda"] = line::num_traits<T>::to_double(r.lambda);
        p["rho"] = line::num_traits<T>::to_double(r.rho);
        p["uniformization"] = q;
        p["drift"] = line::num_traits<T>::to_double(r.drift);
        p["decayRate"] = r.decayRate;
        p["pi0"] = line::num_traits<T>::to_double(r.pi0);
        p["QLen"] = line::num_traits<T>::to_double(r.meanQueueLength);
        p["Util"] = line::num_traits<T>::to_double(r.utilization);
        p["Tput"] = line::num_traits<T>::to_double(r.throughput);
        p["truncLevel"] = r.truncLevel;
        p["truncError"] = r.truncError;
        p["gConverged"] = r.gConverged;
        p["theta"] = vector_json(r.theta);
        p["alpha"] = vector_json(r.alpha);
        line::reg::Json lv = line::reg::Json::array();
        for (std::size_t n = 0; n < r.levelProb.rows(); ++n) {
            line::reg::Json row = line::reg::Json::array();
            for (std::size_t j = 0; j < r.levelProb.cols(); ++j)
                row.push_back(line::num_traits<T>::to_double(r.levelProb(n, j)));
            lv.push_back(row);
        }
        p["levelProb"] = lv;
        emit_analysis<T>("internals", p, std::string());
        return 0;
    }
    std::printf("SolverMAM arith=%s stations=%zu\n", line::num_traits<T>::name(), sn.nstations);
    std::printf("lambda=%.10g rho=%.10g q=%.10g drift=%.10g decayRate=%.10g\n",
                line::num_traits<T>::to_double(r.lambda), line::num_traits<T>::to_double(r.rho), q,
                line::num_traits<T>::to_double(r.drift), r.decayRate);
    std::printf("pi0=%.10g QLen=%.10g Util=%.10g Tput=%.10g truncLevel=%zu truncError=%.3g "
                "gConverged=%s\n",
                line::num_traits<T>::to_double(r.pi0),
                line::num_traits<T>::to_double(r.meanQueueLength),
                line::num_traits<T>::to_double(r.utilization),
                line::num_traits<T>::to_double(r.throughput), r.truncLevel, r.truncError,
                r.gConverged ? "yes" : "no");
    std::printf("%-8s %16s\n", "Level", "Prob");
    for (std::size_t n = 0; n < r.levelProb.rows(); ++n) {
        double s = 0.0;
        for (std::size_t j = 0; j < r.levelProb.cols(); ++j)
            s += line::num_traits<T>::to_double(r.levelProb(n, j));
        std::printf("%-8zu %16.10g\n", n, s);
    }
    return 0;
}

/**
 * Solve a Network model.json with SolverBA and print the same table.
 *
 * A BOUND, not an estimate: `--method` names which side of which hierarchy is
 * wanted (`aba.upper`, `gb.lower`, ...) and `default` resolves to `gb.upper`,
 * as in the reference. `options.level` keeps its default of 2. The table is
 * therefore not comparable with an exact solver's row except as a bracket,
 * which is why the parity row is separate.
 */
/** A JSON flag value, given inline or as the path of a file holding it. */
line::reg::Json read_json_arg(const std::string& spec, const char* flag) {
    std::string text = spec;
    const std::size_t at = spec.find_first_not_of(" \t\r\n");
    if (at == std::string::npos || (spec[at] != '{' && spec[at] != '[')) {
        std::ifstream in(spec.c_str());
        if (!in)
            throw line::InputError(std::string(flag) + " is neither inline JSON nor a readable "
                                   "file (got '" + spec + "')");
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    try {
        return line::reg::Json::parse(text);
    } catch (const line::reg::Json::parse_error& e) {
        throw line::InputError(std::string("malformed ") + flag + " JSON: " + e.what());
    }
}

/** A JSON array of arrays as an integer table. */
std::vector<std::vector<int> > json_int_table(const line::reg::Json& j, const char* name) {
    if (!j.is_array())
        throw line::InputError(std::string("--qrf-params ") + name + " must be an array of rows");
    std::vector<std::vector<int> > out;
    for (std::size_t m = 0; m < j.size(); ++m) {
        if (!j[m].is_array())
            throw line::InputError(std::string("--qrf-params ") + name + " must be an array of "
                                   "rows");
        std::vector<int> row;
        for (std::size_t c = 0; c < j[m].size(); ++c) row.push_back(j[m][c].get<int>());
        out.push_back(row);
    }
    return out;
}

/**
 * Decode `--qrf-params` into `BaOptions::qrf_params`.
 *
 * Required: f, MR, BB, MM, MM1, ZZ, exactly what `sn_to_qrf_params` demands.
 * F is optional and falls back to sn.cap. ZM is DERIVED from ZZ and a supplied
 * one is ignored: it is max(ZZ) by definition, and a larger one empties the
 * polytope through THM3I instead of failing cleanly.
 */
void decode_qrf_params(const std::string& spec, line::ba::BaOptions& opt) {
    const line::reg::Json j = read_json_arg(spec, "--qrf-params");
    if (!j.is_object()) throw line::InputError("--qrf-params must be a JSON object");
    const char* required[] = {"f", "MR", "BB", "MM", "MM1", "ZZ"};
    std::string missing;
    for (std::size_t i = 0; i < 6; ++i)
        if (!j.contains(required[i]))
            missing += (missing.empty() ? "" : ", ") + std::string(required[i]);
    if (!missing.empty())
        throw line::InputError("--qrf-params is missing the field(s) " + missing +
                               "; required are f, MR, BB, MM, MM1, ZZ (F is optional, ZM is "
                               "derived from ZZ)");
    line::ba::BaOptions::QrfParams p;
    p.supplied = true;
    p.f = j["f"].get<int>();
    p.MR = j["MR"].get<int>();
    p.BB = json_int_table(j["BB"], "BB");
    p.MM = json_int_table(j["MM"], "MM");
    p.MM1 = json_int_table(j["MM1"], "MM1");
    if (!j["ZZ"].is_array()) throw line::InputError("--qrf-params ZZ must be an array");
    for (std::size_t i = 0; i < j["ZZ"].size(); ++i) p.ZZ.push_back(j["ZZ"][i].get<int>());
    if (j.contains("F"))
        for (std::size_t i = 0; i < j["F"].size(); ++i) p.F.push_back(j["F"][i].get<int>());
    opt.qrf_params = p;
}

/** Decode `--qrf-alpha`, the (nstations x N) load-dependent scaling. */
void decode_qrf_alpha(const std::string& spec, line::ba::BaOptions& opt) {
    const line::reg::Json j = read_json_arg(spec, "--qrf-alpha");
    if (!j.is_array() || j.empty() || !j[0].is_array())
        throw line::InputError("--qrf-alpha must be a JSON array of rows, one per station");
    line::Matrix<double> a(j.size(), j[0].size(), 0.0);
    for (std::size_t i = 0; i < j.size(); ++i) {
        if (!j[i].is_array() || j[i].size() != j[0].size())
            throw line::InputError("--qrf-alpha rows must all have the same length");
        for (std::size_t n = 0; n < j[i].size(); ++n) a(i, n) = j[i][n].get<double>();
    }
    opt.qrf_alpha = a;
}

/** The knobs the two SolverBA arms share. */
void apply_ba_knobs(const Knobs& k, line::ba::BaOptions& opt) {
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.level > 0) opt.level = k.level;
    if (!k.qrf_params.empty()) decode_qrf_params(k.qrf_params, opt);
    if (!k.qrf_alpha.empty()) decode_qrf_alpha(k.qrf_alpha, opt);
}

template <class T>
int solve_model_ba(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ba::BaOptions opt;
    apply_ba_knobs(k, opt);
    const line::mva::AvgResult<T> r = line::ba::solver_ba_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    std::printf("SolverBA arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                r.actualmethod.c_str(),
                line::util::method_type("BA", r.actualmethod).c_str());
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * `-s ba -a bounds`: `SolverBA.getBounds` and its `getBoundsTable` row filter.
 *
 * NOT THE SAME ANSWER AS `-a avg`, which is why it is its own analysis. `-a avg`
 * reports ONE side of ONE family -- whichever `--method` named -- so a caller who
 * wants the bracket has to run the solver twice and know which two method names
 * pair up. `getBounds` takes the family (the method's prefix before the first
 * dot) and re-runs both sides under the caller's FULL option set, so a
 * hierarchical family tightens with `--level` as it should.
 *
 * A ONE-SIDED FAMILY REPORTS NaN ON THE SIDE IT LACKS, never zero: `cub` is
 * upper-only and `mbjb`/`ldbcmp` are lower-only, and a zero there would read as
 * a lower bound of zero rather than as the absence of one.
 */
template <class T>
int solve_model_ba_bounds(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ba::BaOptions opt;
    apply_ba_knobs(k, opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::ba::BaBounds<T> b = line::ba::ba_bounds(sn, opt);
    const std::string am = line::ba::resolve_method(opt.method);

    auto d = [](const T& v) { return line::num_traits<T>::to_double(v); };
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "BoundsTable";
        p["indexBase"] = 0;
        p["family"] = am.substr(0, am.find('.'));
        p["hasLower"] = b.has_lower;
        p["hasUpper"] = b.has_upper;
        for (const char* key : {"Station", "JobClass", "QLower", "QUpper", "TLower", "TUpper"})
            p[key] = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                if (!b.keep[i][c]) continue;
                p["Station"].push_back(sn.stations[i].name);
                p["JobClass"].push_back(sn.classes[c].name);
                p["QLower"].push_back(d(b.Qlower(i, c)));
                p["QUpper"].push_back(d(b.Qupper(i, c)));
                p["TLower"].push_back(d(b.Tlower(i, c)));
                p["TUpper"].push_back(d(b.Tupper(i, c)));
            }
        emit_analysis<T>("bounds", p, am);
        return 0;
    }
    std::printf("SolverBA arith=%s method=%s type=%s family=%s sides=%s\n",
                line::num_traits<T>::name(), am.c_str(),
                line::util::method_type("BA", am).c_str(), am.substr(0, am.find('.')).c_str(),
                b.has_lower ? (b.has_upper ? "lower,upper" : "lower") : "upper");
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "QLower", "QUpper",
                "TLower", "TUpper");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (!b.keep[i][c]) continue;
            std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n", sn.stations[i].name.c_str(),
                        sn.classes[c].name.c_str(), d(b.Qlower(i, c)), d(b.Qupper(i, c)),
                        d(b.Tlower(i, c)), d(b.Tupper(i, c)));
        }
    return 0;
}

/**
 * `-s qns`: the model handed to the external `qnsolver` binary.
 *
 * The only wrapper on the model-solving path. The banner names it separately
 * from the table so a reader can tell an independent tool's numbers from the
 * port's own -- which is the whole reason the wrapper exists.
 */
template <class T>
int solve_model_qns(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::qns::QnsOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    // --samples is deliberately NOT forwarded: `options.samples` reaches the
    // JMVA document as `maxSamples`, which is JMT's Monte Carlo cap and which
    // qnsolver reads past. The shared knob ladder below refuses it for every
    // solver that draws nothing, and QNS is one of them.
    opt.timeout = k.timeout_seconds;
    opt.keep = k.keep;
    const line::mva::AvgResult<T> r = line::qns::solver_qns_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    if (!g_json_output)
        std::printf("SolverQNS arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                    r.actualmethod.c_str(),
                    line::util::method_type("QNS", r.actualmethod).c_str());
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * Emit the base-class exponential response-time CDF fallback, the same
 * `CdfRespT` document every real distributional arm emits.
 *
 * This is `@@NetworkSolver/getCdfRespT.m`: the solvers without a
 * distributional result of their own (MVA, QNS, BA, AG) inherit an exponential
 * law with the right mean in the reference, and refusing `-a cdf` for them
 * diverged from it. The curve says nothing about the tail; the banner names
 * the solver so the caller knows which mean it wraps.
 */
template <class T>
int emit_default_cdf(const char* solver_name, const line::qn::NetworkStruct<T>& sn,
                     const line::mva::AvgResult<T>& r) {
    const std::vector<std::vector<line::solvers::DefaultCdfCurve>> RD =
        line::solvers::solver_default_cdf_respt<T>(sn, r.RN);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "CdfRespT";
        p["indexBase"] = 0;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < RD.size(); ++i)
            for (std::size_t c = 0; c < RD[i].size(); ++c) {
                if (RD[i][c].t.empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = line::reg::Json(RD[i][c].t);
                e["F"] = line::reg::Json(RD[i][c].F);
                arr.push_back(e);
            }
        p["respt"] = arr;
        emit_analysis<T>("cdf", p, std::string());
        return 0;
    }
    std::printf("%s arith=%s method=%s (exponential fallback with the solver's mean)\n",
                solver_name, line::num_traits<T>::name(), r.actualmethod.c_str());
    std::printf("%-16s %-14s %14s %14s\n", "Station", "JobClass", "Time", "F(t)");
    for (std::size_t i = 0; i < RD.size(); ++i)
        for (std::size_t c = 0; c < RD[i].size(); ++c)
            for (std::size_t j = 0; j < RD[i][c].t.size(); ++j)
                std::printf("%-16s %-14s %14.8g %14.10g\n", sn.stations[i].name.c_str(),
                            sn.classes[c].name.c_str(), RD[i][c].t[j], RD[i][c].F[j]);
    return 0;
}

/** `-s mva -a cdf`: the inherited exponential fallback over the MVA means. */
template <class T>
int solve_model_mva_cdf(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::mva::MvaOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    if (!k.fork_join.empty()) opt.fork_join = k.fork_join;
    line::Matrix<T> init;
    const line::mva::AvgResult<T> r = line::mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
    return emit_default_cdf<T>("SolverMVA", net.get_struct(), r);
}

/** `-s ag -a cdf`: the inherited exponential fallback over the RCAT means. */
template <class T>
int solve_model_ag_cdf(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ag::AgOptions opt;
    apply_ag_knobs(k, opt);
    const line::mva::AvgResult<T> r = line::ag::solver_ag_run_analyzer(net.get_struct(), opt);
    return emit_default_cdf<T>("SolverAG", net.get_struct(), r);
}

/** `-s ba -a cdf`: the inherited exponential fallback over the bound means. */
template <class T>
int solve_model_ba_cdf(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ba::BaOptions opt;
    apply_ba_knobs(k, opt);
    const line::mva::AvgResult<T> r = line::ba::solver_ba_run_analyzer(net.get_struct(), opt);
    return emit_default_cdf<T>("SolverBA", net.get_struct(), r);
}

/** `-s qns -a cdf`: the inherited exponential fallback over qnsolver's means. */
int solve_model_qns_cdf(const std::string& file, const Knobs& k) {
    line::qn::Network<double> net = read_model<double>(file);
    line::qns::QnsOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    opt.timeout = k.timeout_seconds;
    opt.keep = k.keep;
    const line::mva::AvgResult<double> r = line::qns::solver_qns_run_analyzer(net.get_struct(), opt);
    return emit_default_cdf<double>("SolverQNS", net.get_struct(), r);
}


/**
 * `-a interval`: the support-only range, the reference's `getIntervalTable`.
 *
 * A SEPARATE FUNCTION because it is printed before the ensemble exists: on the
 * exact path there is no design to solve, so the interval arrives without a
 * `UqSolution` beside it.
 */
template <class T>
int print_uq_interval(const line::qn::NetworkStruct<T>& sn, const line::uq::UqInterval<T>& iv,
                      const std::string& stage, std::size_t npriors) {
    auto v = [](const line::Matrix<T>& M, std::size_t i, std::size_t c) {
        return M.empty() ? 0.0 : line::num_traits<T>::to_double(M(i, c));
    };
    line::reg::Json p = line::reg::Json::object();
    p["type"] = "IntervalTable";
    p["indexBase"] = 0;
    p["exact"] = iv.exact;
    p["intervalMethod"] = iv.method;
    if (!iv.why.empty()) p["why"] = iv.why;
    if (iv.has_totals) {
        p["X"] = {line::num_traits<T>::to_double(iv.Xlo),
                  line::num_traits<T>::to_double(iv.Xup)};
        p["Rtot"] = {line::num_traits<T>::to_double(iv.Rtot_lo),
                     line::num_traits<T>::to_double(iv.Rtot_up)};
    }
    for (const char* key : {"Station", "JobClass", "QLen_lo", "QLen_up", "Util_lo", "Util_up",
                            "RespT_lo", "RespT_up", "Tput_lo", "Tput_up"})
        p[key] = line::reg::Json::array();
    if (!g_json_output) {
        std::printf("SolverUQ arith=%s interval=%s exact=%s priors=%zu stage=%s\n",
                    line::num_traits<T>::name(), iv.method.c_str(), iv.exact ? "yes" : "no",
                    npriors, stage.c_str());
        // A RANGE THAT IS NOT AN ENCLOSURE MUST SAY SO. The sampled path
        // spans the design points only, so on a continuous Prior it lies
        // strictly inside the true range; printing it beside an exact hull
        // without the reason would make the two indistinguishable.
        if (!iv.exact)
            std::fprintf(stderr,
                         "Warning: exact interval MVA does not apply (%s); the range below is "
                         "over the solved design points and is not an enclosure.\n",
                         iv.why.c_str());
        std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s %12s %12s\n", "Station",
                    "JobClass", "QLen_lo", "QLen_up", "Util_lo", "Util_up", "RespT_lo",
                    "RespT_up", "Tput_lo", "Tput_up");
    }
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            // The reference's row filter: an upper endpoint of zero on every
            // presence metric means the class never visits the station.
            if (v(iv.Qup, i, c) <= 0.0 && v(iv.Uup, i, c) <= 0.0 && v(iv.Tup, i, c) <= 0.0)
                continue;
            if (!g_json_output) {
                std::printf(
                    "%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                    sn.stations[i].name.c_str(), sn.classes[c].name.c_str(), v(iv.Qlo, i, c),
                    v(iv.Qup, i, c), v(iv.Ulo, i, c), v(iv.Uup, i, c), v(iv.Rlo, i, c),
                    v(iv.Rup, i, c), v(iv.Tlo, i, c), v(iv.Tup, i, c));
                continue;
            }
            p["Station"].push_back(sn.stations[i].name);
            p["JobClass"].push_back(sn.classes[c].name);
            p["QLen_lo"].push_back(v(iv.Qlo, i, c));
            p["QLen_up"].push_back(v(iv.Qup, i, c));
            p["Util_lo"].push_back(v(iv.Ulo, i, c));
            p["Util_up"].push_back(v(iv.Uup, i, c));
            p["RespT_lo"].push_back(v(iv.Rlo, i, c));
            p["RespT_up"].push_back(v(iv.Rup, i, c));
            p["Tput_lo"].push_back(v(iv.Tlo, i, c));
            p["Tput_up"].push_back(v(iv.Tup, i, c));
        }
    if (g_json_output) emit_analysis<T>("interval", p, iv.method);
    return 0;
}

/**
 * Solve a Network model.json carrying a Prior with SolverUQ.
 *
 * TWO FLAGS MEAN SOMETHING ELSE HERE, and both are UQ's own rather than the
 * stage solver's: `--method` names the DESIGN (quadrature or montecarlo, the
 * reference's `options.method`), and `--samples` the number of nodes per
 * continuous Prior (`options.samples`, 11 by default and not the simulation
 * default, since each node is a full solver run). The engine that runs at each
 * point is `--uq-solver`, and it keeps its own defaults for everything except
 * the convergence knobs, which UQ does not have and therefore passes through.
 * A caller who wants an SSA run length AND a UQ design cannot state both, so
 * the SSA stage keeps its default sample count; that is stated in --help rather
 * than resolved by giving one flag two meanings.
 *
 * `-a posterior` prints the design itself -- every point, its weight and its
 * metrics -- because the expectation alone hides whether it averaged two nearby
 * models or two wildly different ones, and that spread IS the answer to an
 * uncertainty question.
 */
template <class T>
int solve_model_uq(const std::string& file, const Knobs& k, const std::string& analysis) {
    line::qn::Network<T> net = read_model<T>(file);
    line::uq::UqOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (k.samples) opt.samples = k.samples;
    if (k.seed) opt.seed = k.seed;
    line::uq::UqStageOptions so;
    so.solver = k.uq_solver;
    so.tol = k.tol;
    so.iter_tol = k.iter_tol;
    so.iter_max = k.iter_max;
    so.cutoff = k.cutoff;

    if (analysis == "interval") {
        // BEFORE THE ENSEMBLE, because the exact path does not need one: it is
        // 2*(m+2) MVA calls over the demand box, and the reference's
        // `getInterval` likewise reaches `intervalByMVA` without touching
        // `self.results`. Only the sampling fallback solves the design.
        const line::uq::UqInterval<T> iv =
            line::uq::uq_interval_run<T>(net, line::uq::uq_stage_solver<T>(so), opt);
        const line::qn::NetworkStruct<T>& isn = net.get_struct();
        return print_uq_interval<T>(isn, iv, so.solver, line::uq::uq_detect_priors(isn).size());
    }

    const line::uq::UqSolution<T> r =
        line::uq::solver_uq_run_analyzer<T>(net, line::uq::uq_stage_solver<T>(so), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    if (analysis == "avg") {
        if (!g_json_output)
            std::printf("SolverUQ arith=%s design=%s points=%zu priors=%zu stage=%s method=%s\n",
                        line::num_traits<T>::name(), r.method.c_str(), r.points.size(),
                        r.sites.size(), so.solver.c_str(), r.avg.actualmethod.c_str());
        print_avg_table<T>(sn, r.avg);
        return 0;
    }

    // -a posterior: the per-design-point table, the reference's getPosteriorTable.
    line::reg::Json p = line::reg::Json::object();
    p["type"] = "PosteriorTable";
    p["indexBase"] = 0;
    p["design"] = r.method;
    p["stage"] = so.solver;
    line::reg::Json sites = line::reg::Json::array();
    for (std::size_t l = 0; l < r.sites.size(); ++l) {
        line::reg::Json s = line::reg::Json::object();
        s["node"] = sn.nodes[r.sites[l].node - 1].name;
        s["class"] = sn.classes[r.sites[l].cls - 1].name;
        s["kind"] = r.sites[l].arrival ? "arrival" : "service";
        sites.push_back(s);
    }
    p["priors"] = sites;
    for (const char* key : {"Point", "Weight", "Station", "JobClass", "QLen", "Util", "RespT",
                            "Tput"})
        p[key] = line::reg::Json::array();
    line::reg::Json means = line::reg::Json::array();

    if (!g_json_output) {
        std::printf("SolverUQ arith=%s design=%s points=%zu priors=%zu stage=%s\n",
                    line::num_traits<T>::name(), r.method.c_str(), r.points.size(), r.sites.size(),
                    so.solver.c_str());
        // THE SUBSTITUTED MEANS ARE PROVENANCE, not decoration: a design point
        // is a model, and this line is what says which one.
        std::printf("%-6s %12s  substituted means\n", "Point", "Weight");
        for (std::size_t e = 0; e < r.points.size(); ++e) {
            std::printf("%-6zu %12.6g ", e + 1, line::num_traits<T>::to_double(r.weights[e]));
            for (std::size_t l = 0; l < r.sites.size(); ++l)
                std::printf(" %s@%s=%.6g", sn.nodes[r.sites[l].node - 1].name.c_str(),
                            sn.classes[r.sites[l].cls - 1].name.c_str(),
                            line::num_traits<T>::to_double(r.design[e].dists[l].mean));
            std::printf("\n");
        }
        std::printf("%-6s %12s %-16s %-14s %12s %12s %12s %12s\n", "Point", "Weight", "Station",
                    "JobClass", "QLen", "Util", "RespT", "Tput");
    }
    for (std::size_t e = 0; e < r.points.size(); ++e) {
        line::reg::Json row = line::reg::Json::array();
        for (std::size_t l = 0; l < r.sites.size(); ++l)
            row.push_back(line::num_traits<T>::to_double(r.design[e].dists[l].mean));
        means.push_back(row);
        const line::mva::AvgResult<T>& a = r.points[e];
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                const double q = line::num_traits<T>::to_double(a.QN(i, c));
                const double u = line::num_traits<T>::to_double(a.UN(i, c));
                const double rr = line::num_traits<T>::to_double(a.RN(i, c));
                const double t = line::num_traits<T>::to_double(a.TN(i, c));
                // getPosteriorTable's row filter, the PRESENCE metrics only: a
                // response time alone does not put a class at a station, and
                // reading it as presence kept rows the reference drops.
                if (q <= 0.0 && u <= 0.0 && t <= 0.0) continue;
                if (!g_json_output) {
                    std::printf("%-6zu %12.6g %-16s %-14s %12.6g %12.6g %12.6g %12.6g\n", e + 1,
                                line::num_traits<T>::to_double(r.weights[e]),
                                sn.stations[i].name.c_str(), sn.classes[c].name.c_str(), q, u, rr,
                                t);
                    continue;
                }
                p["Point"].push_back(e);
                p["Weight"].push_back(line::num_traits<T>::to_double(r.weights[e]));
                p["Station"].push_back(sn.stations[i].name);
                p["JobClass"].push_back(sn.classes[c].name);
                p["QLen"].push_back(q);
                p["Util"].push_back(u);
                p["RespT"].push_back(rr);
                p["Tput"].push_back(t);
            }
    }
    if (g_json_output) {
        p["substitutedMean"] = means;
        emit_analysis<T>("posterior", p, r.avg.actualmethod);
    }
    return 0;
}

/**
 * The banner every `-s ctmc` analysis opens with.
 *
 * It carries the state count and the cutoff because a CTMC number is not the
 * model's number without them: the space is what was enumerated, and on an open
 * model the cutoff is what truncated it.
 */
template <class T>
void print_ctmc_banner(const line::ctmc::CtmcSolution<T>& d) {
    std::size_t cut = 0;
    for (std::size_t i = 0; i < d.cutoff.size(); ++i) cut = std::max(cut, d.cutoff[i]);
    if (cut)
        std::printf("SolverCTMC arith=%s method=%s type=%s states=%zu cutoff=%zu\n",
                    line::num_traits<T>::name(), d.actualmethod.c_str(),
                    line::util::method_type("CTMC", d.actualmethod).c_str(), d.chain.space.size(),
                    cut);
    else
        std::printf("SolverCTMC arith=%s method=%s type=%s states=%zu\n",
                    line::num_traits<T>::name(), d.actualmethod.c_str(),
                    line::util::method_type("CTMC", d.actualmethod).c_str(), d.chain.space.size());
    // The library never writes to stderr, so an unseeded reducible mixture is
    // reported here or not at all -- and it is the one case where the printed
    // distribution is not the model's.
    if (!d.warning.empty()) std::fprintf(stderr, "warning: %s\n", d.warning.c_str());
}

/**
 * The banner's own content as JSON, added to every `-s ctmc` payload.
 *
 * A CTMC number is not the model's number without them, which is why the
 * readable banner carries them: the space is what was actually enumerated, and
 * on an open model the cutoff is what truncated it. A host that reported the
 * occupancy of a truncated chain as the model's would be reporting a different
 * model's answer, so the two travel with the payload rather than only above it.
 */
template <class T>
line::reg::Json ctmc_meta(const line::ctmc::CtmcSolution<T>& d) {
    line::reg::Json m = line::reg::Json::object();
    m["states"] = d.chain.space.size();
    std::size_t cut = 0;
    for (std::size_t i = 0; i < d.cutoff.size(); ++i) cut = std::max(cut, d.cutoff[i]);
    if (cut) m["cutoff"] = cut;
    return m;
}

/**
 * `-a avg`: the AvgTable, on whichever path the model's regions require.
 *
 * `solver_ctmc_analyzer_any` and not `solver_ctmc_analyzer`, so a WAITQ region
 * reaches the augmented walk that carries its token FIFO instead of the
 * lattice analyzer, which refuses it. Every other region rule keeps refusing.
 */
template <class T>
int solve_ctmc_avg(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    const line::ctmc::CtmcAnySolution<T> a = line::ctmc::solver_ctmc_analyzer_any(sn, opt);
    const line::mva::AvgResult<T> r = line::ctmc::solver_ctmc_avg_table(sn, a.sol, opt.method);
    print_ctmc_banner<T>(a.sol);
    print_avg_table<T>(sn, r);
    // A PARKED JOB IS IN NO QLen COLUMN: it left its station and sits in the
    // region's FIFO, so the model's population balances only once this is read
    // alongside the table rather than instead of it.
    if (a.waitq) {
        std::printf("%-16s %-14s %12s\n", "Region", "JobClass", "Parked");
        for (std::size_t c = 0; c < sn.nclasses && c < a.parked.size(); ++c)
            std::printf("%-16s %-14s %12.6g\n", "(all)", sn.classes[c].name.c_str(),
                        line::num_traits<T>::to_double(a.parked[c]));
    }
    return 0;
}

/**
 * `-a avg` under `--method mdd`: the AvgTable from the level aggregation.
 *
 * A SEPARATE ARM, and not a branch inside `solve_ctmc_avg`, because the method
 * never forms the |S|-state generator: there is no state space to print a state
 * count from and no cutoff to report, so the banner is the diagram's own -- the
 * cardinality it counts symbolically, what the levels actually hold, and whether
 * the answer is certified exact. Reporting the generator banner here would name
 * a chain that was never built.
 */
template <class T>
int solve_ctmc_mdd_avg(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                       const line::mdd::MddMcdOptions& mcdopt) {
    // EVERY BACKEND, `exact` included: the level solve picks Householder QR in
    // floating point and `line::lstsq` under exact arithmetic, so nothing here
    // takes a square root that a rational field has no answer for.
    const line::ctmc::CtmcMddSolution<T> s = line::ctmc::solver_ctmc_mdd_analyzer(sn, opt, mcdopt);
    line::ctmc::CtmcSolution<T> d;
    d.avg = s.avg;
    d.actualmethod = s.actualmethod;
    const line::mva::AvgResult<T> r = line::ctmc::solver_ctmc_avg_table(sn, d, opt.method);
    std::size_t held = 0;
    for (std::size_t k = 0; k < s.level_sizes.size(); ++k) held += s.level_sizes[k];
    std::printf(
        "SolverCTMC arith=%s method=%s type=%s states=%lld held=%zu levels=%zu iters=%d "
        "encoding=%s exact=%s\n",
        line::num_traits<T>::name(), s.actualmethod.c_str(),
        line::util::method_type("CTMC", s.actualmethod).c_str(), s.num_states, held,
        s.level_sizes.size(), s.iters, s.encoding.c_str(),
        // "certified" is not "exact": a product-form model is exact however
        // much its diagram shares, so the false case says only that the
        // STRUCTURAL test did not fire.
        s.no_aggregation ? "certified" : "product-form-only");
    print_avg_table<T>(sn, r);
    return 0;
}

/**
 * `-a avg` under `--method cftp` / `cftp.approx`: the AvgTable from perfect
 * sampling.
 *
 * THE BANNER NAMES THE RUN LENGTH because these numbers carry Monte Carlo error
 * and are not comparable to an exact solver's at solver tolerance. The mean
 * coalescence horizon travels with it: it is the cost the exact sampler paid,
 * has no a-priori bound, and is the one number that says whether the draw was
 * cheap or the model is nearly saturated.
 */
template <class T>
int solve_ctmc_cftp_avg(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                        const line::ctmc::CtmcCftpOptions& cftpopt) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the cftp methods draw random states and form the station balance functions in the "
            "log domain, neither of which exists in exact rational arithmetic; rerun with --arith "
            "double or --arith real");
    } else {
        const line::ctmc::CtmcCftpSolution<T> s =
            line::ctmc::solver_ctmc_cftp(sn, opt, cftpopt);
        line::ctmc::CtmcSolution<T> d;
        d.avg = s.avg;
        d.actualmethod = s.actualmethod;
        const line::mva::AvgResult<T> r = line::ctmc::solver_ctmc_avg_table(sn, d, opt.method);
        double horizon = 0.0;
        for (std::size_t i = 0; i < s.horizon.size(); ++i)
            horizon += static_cast<double>(s.horizon[i]);
        if (!s.horizon.empty()) horizon /= static_cast<double>(s.horizon.size());
        std::printf(
            "SolverCTMC arith=%s method=%s type=%s samples=%zu seed=%lu distinct=%zu "
            "meanhorizon=%.6g\n",
            line::num_traits<T>::name(), s.actualmethod.c_str(),
            line::util::method_type("CTMC", s.actualmethod).c_str(), cftpopt.samples,
            cftpopt.seed, s.distinct_states.size(), horizon);
        print_avg_table<T>(sn, r);
        return 0;
    }
}

/**
 * `-a prob`: the four SolverCTMC probability queries over the model's default
 * initial state -- `getProbSys`, `getProbSysAggr`, and `getProb`/`getProbAggr`
 * per station.
 *
 * ALL FOUR AND NOT ONE, because the joint and the aggregate answer different
 * questions and the pair is what makes either readable: `getProbSys` is the
 * probability of exactly that state, phases and buffer arrangement included,
 * while `getProbSysAggr` sums over every arrangement realizing the same
 * per-class counts. On a single-phase model with no buffer ordering the two
 * coincide, and where they do not the ratio is what the encoding added.
 *
 * NOT THE SAME NUMBERS AS `-s mva -a prob` OR `-s nc -a prob`, and that is the
 * point of having all three: MVA fits a binomial to its own means and NC takes a
 * ratio of normalizing constants under the product form, whereas these are the
 * stationary law of the chain itself and are exact for any model the chain
 * represents, product-form or not.
 */
template <class T>
int solve_ctmc_prob(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                    const Knobs& k) {
    const line::ctmc::CtmcSolution<T> d = line::ctmc::solver_ctmc_analyzer(sn, opt);
    line::qn::NetState<T> init;
    if (!line::ctmc::analyzer_detail::default_init_state(sn, init))
        throw line::UnsupportedError(
            "-a prob reports the probability of the model's DEFAULT INITIAL STATE, and this "
            "model's initial marking admits no state; check the class populations against their "
            "reference stations");
    // `--state` is `getProb(node, state)`'s second argument: the ENCODED ROW of
    // one stateful node's own state space, not a per-class job count. The
    // reference substitutes it into `sn.state{node}` and leaves every other
    // node at its default, which is what happens here -- so the station
    // marginals below are the requested state's, and the two system
    // probabilities are that state's joined with the rest of the default
    // marking, exactly as `setState` followed by `getProbSys` reports it.
    if (!k.state.empty()) {
        if (!k.node)
            throw line::InputError(
                "--state is the state of ONE node and needs --node to say which; a bare state "
                "vector cannot be matched against a network whose nodes have different widths");
        const std::size_t isf = sn.stateful_index(k.node);
        if (isf == 0)
            throw line::InputError("--node " + std::to_string(k.node) +
                                   " is not a stateful node, so it holds no state to ask about");
        // THE WIDTH MUST MATCH EXACTLY, and a short vector is refused rather
        // than padded. Every row of a node's block is stored at the node's
        // widest encoding, so a padded row IS a state -- just not the one the
        // caller named: on an FCFS queue holding two jobs the encoding is
        // (buffer, in-service phase count) = (1,1), and padding `--state 2` to
        // (0,2) names a state the chain never visits, which would answer 0
        // where the caller expected the marginal. The width is reported so the
        // next attempt can be right.
        const std::size_t w = d.chain.space.empty() ? k.state.size()
                                                    : d.chain.space[0].local[isf - 1].size();
        if (k.state.size() != w)
            throw line::InputError(
                "--state has " + std::to_string(k.state.size()) + " entries but node " +
                std::to_string(k.node) + " encodes its state in " + std::to_string(w) +
                "; getProb(node, state) takes the node's whole encoded row, and a shorter one "
                "padded with zeros is a different state rather than a partial one (use -a marg "
                "for a per-class job-count marginal)");
        std::vector<T> row(w, line::num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < k.state.size(); ++i)
            row[i] = line::num_traits<T>::from_int(k.state[i]);
        init.local[isf - 1] = row;
    }
    const T psys = line::ctmc::solver_ctmc_joint(sn, d, init);
    const T psysaggr = line::ctmc::solver_ctmc_jointaggr(sn, d, init);
    const std::vector<T> pmarg = line::ctmc::solver_ctmc_marg(sn, d, init);
    const std::vector<T> pmargaggr = line::ctmc::solver_ctmc_margaggr(sn, d, init);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbAggr";
        p["indexBase"] = 0;
        p["ProbSys"] = line::num_traits<T>::to_double(psys);
        p["ProbSysAggr"] = line::num_traits<T>::to_double(psysaggr);
        line::reg::Json st = line::reg::Json::array(), pm = line::reg::Json::array(),
                        pa = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            st.push_back(sn.stations[i].name);
            pm.push_back(line::num_traits<T>::to_double(pmarg[i]));
            pa.push_back(line::num_traits<T>::to_double(pmargaggr[i]));
        }
        p["Station"] = st;
        p["Prob"] = pm;
        p["ProbAggr"] = pa;
        emit_analysis<T>("prob", p, d.actualmethod, ctmc_meta<T>(d));
        return 0;
    }
    print_ctmc_banner<T>(d);
    std::printf("ProbSys %.10g\n", line::num_traits<T>::to_double(psys));
    std::printf("ProbSysAggr %.10g\n", line::num_traits<T>::to_double(psysaggr));
    std::printf("%-16s %14s %14s\n", "Station", "Prob", "ProbAggr");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        std::printf("%-16s %14.10g %14.10g\n", sn.stations[i].name.c_str(),
                    line::num_traits<T>::to_double(pmarg[i]),
                    line::num_traits<T>::to_double(pmargaggr[i]));
    return 0;
}

/**
 * `-a gen`: `getInfGen`, as Q in sparse triplets plus the synchronization list
 * its event filtration is indexed by.
 *
 * Q IS PRINTED SPARSELY. A generator's rows hold one entry per enabled
 * transition and the space can run to thousands of states, so the dense form
 * would be quadratic in a quantity that is linear in the model.
 */
/**
 * The derived START/PREEMPT filtration as one JSON block per (station, class),
 * each carrying its own 0-based Station and Class beside the usual From/To/Rate
 * triplets. An all-zero block is omitted: on a model with no preemption that is
 * every block of the PREEMPT filtration.
 */
template <class T>
line::reg::Json aux_filt_json(const std::vector<std::vector<line::Matrix<T> > >& filt,
                              std::size_t n) {
    line::reg::Json blocks = line::reg::Json::array();
    for (std::size_t i = 0; i < filt.size(); ++i)
        for (std::size_t r = 0; r < filt[i].size(); ++r) {
            line::reg::Json bfrom = line::reg::Json::array(), bto = line::reg::Json::array(),
                            brate = line::reg::Json::array();
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b) {
                    const double q = line::num_traits<T>::to_double(filt[i][r](a, b));
                    if (q == 0.0) continue;
                    bfrom.push_back(a);
                    bto.push_back(b);
                    brate.push_back(q);
                }
            if (bfrom.empty()) continue;
            line::reg::Json e = line::reg::Json::object();
            e["Station"] = i;
            e["Class"] = r;
            e["From"] = bfrom;
            e["To"] = bto;
            e["Rate"] = brate;
            blocks.push_back(e);
        }
    return blocks;
}

template <class T>
int solve_ctmc_gen(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    line::ctmc::CtmcOptions o = opt;
    o.keep_filtration = true;  // the filtration is half of what getInfGen returns
    const line::ctmc::CtmcSolution<T> d = line::ctmc::solver_ctmc_analyzer(sn, o);
    const line::ctmc::CtmcGenerator<T> g = line::ctmc::ctmc_get_infgen(sn, d);
    const std::size_t n = g.Q.rows();
    std::size_t nnz = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (line::num_traits<T>::to_double(g.Q(i, j)) != 0.0) ++nnz;
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "InfGen";
        p["indexBase"] = 0;
        p["size"] = n;
        p["nnz"] = nnz;
        // SPARSE HERE TOO, for the reason the table is: a generator has one entry
        // per enabled transition, so the dense form is quadratic in a quantity
        // that is linear in the model. A host rebuilds Q with one scatter.
        line::reg::Json from = line::reg::Json::array(), to = line::reg::Json::array(),
                        rate = line::reg::Json::array();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const double q = line::num_traits<T>::to_double(g.Q(i, j));
                if (q == 0.0) continue;
                from.push_back(i);
                to.push_back(j);
                rate.push_back(q);
            }
        p["From"] = from;
        p["To"] = to;
        p["Rate"] = rate;
        // The synchronization list is HALF of what getInfGen returns: without it
        // the filtration's index has no meaning, so it travels with Q. Its node
        // and class are the port's own 1-BASED indices, with 0 the LOCAL dummy
        // -- `indexBase` above governs the state indices, which is where a host
        // would otherwise be indexing a chain by an off-by-one.
        line::reg::Json sync = line::reg::Json::array();
        for (std::size_t a = 0; a < g.sync.size() && a < g.filt.size(); ++a) {
            // THE FILTER ITSELF, not only its nnz as the readable table reports:
            // `eventFilt` is half of what getInfGen returns, and a host handed
            // only Q cannot recover which synchronization contributed a rate --
            // Q's entries have already summed every one of them.
            std::size_t fnz = 0;
            line::reg::Json ffrom = line::reg::Json::array(), fto = line::reg::Json::array(),
                            frate = line::reg::Json::array();
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) {
                    const double q = line::num_traits<T>::to_double(g.filt[a](i, j));
                    if (q == 0.0) continue;
                    ++fnz;
                    ffrom.push_back(i);
                    fto.push_back(j);
                    frate.push_back(q);
                }
            line::reg::Json e = line::reg::Json::object();
            e["From"] = ffrom;
            e["To"] = fto;
            e["Rate"] = frate;
            e["activeEvent"] = line::lang::event_to_text(g.sync[a].active.event);
            e["activeNode"] = g.sync[a].active.node;
            e["activeClass"] = g.sync[a].active.cls;
            e["passiveEvent"] = line::lang::event_to_text(g.sync[a].passive.event);
            e["passiveNode"] = g.sync[a].passive.node;
            e["passiveClass"] = g.sync[a].passive.cls;
            e["nnz"] = fnz;
            sync.push_back(e);
        }
        p["sync"] = sync;
        // The DERIVED filtrations travel under their own keys, one block per
        // (station, class), because they are NOT synchronizations: a START
        // rides on an arc `sync` already carries, so folding them in would make
        // a host summing the filtration double-count the generator.
        p["startFilt"] = aux_filt_json<T>(g.start_filt, n);
        p["preemptFilt"] = aux_filt_json<T>(g.preempt_filt, n);
        emit_analysis<T>("gen", p, d.actualmethod, ctmc_meta<T>(d));
        return 0;
    }
    print_ctmc_banner<T>(d);
    std::printf("InfGen events=%zu nnz=%zu\n", g.sync.size(), nnz);
    std::printf("%8s %8s %16s\n", "From", "To", "Rate");
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double q = line::num_traits<T>::to_double(g.Q(i, j));
            if (q != 0.0) std::printf("%8zu %8zu %16.10g\n", i + 1, j + 1, q);
        }
    std::printf("%6s %-10s %6s %6s %-10s %6s %6s %8s\n", "Event", "ActEvent", "ActNode", "ActCls",
                "PasEvent", "PasNode", "PasCls", "Nnz");
    for (std::size_t a = 0; a < g.sync.size() && a < g.filt.size(); ++a) {
        std::size_t fnz = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                if (line::num_traits<T>::to_double(g.filt[a](i, j)) != 0.0) ++fnz;
        std::printf("%6zu %-10s %6zu %6zu %-10s %6zu %6zu %8zu\n", a + 1,
                    line::lang::event_to_text(g.sync[a].active.event), g.sync[a].active.node,
                    g.sync[a].active.cls, line::lang::event_to_text(g.sync[a].passive.event),
                    g.sync[a].passive.node, g.sync[a].passive.cls, fnz);
    }
    return 0;
}

/** `-a states`: `getStateSpace` and `getStateSpaceAggr`, side by side. */
template <class T>
int solve_ctmc_states(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    const line::ctmc::CtmcSolution<T> d = line::ctmc::solver_ctmc_analyzer(sn, opt);
    const line::ctmc::CtmcStateSpace<T> s = line::ctmc::ctmc_get_state_space(sn, d);
    const line::Matrix<T> A = line::ctmc::ctmc_get_state_space_aggr(sn, d);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "StateSpace";
        p["indexBase"] = 0;
        // The per-node widths are the ONLY thing that makes a flat row decodable:
        // they are where one stateful node's block ends and the next begins.
        p["NodeWidths"] = index_json(s.node_width);
        p["space"] = matrix_json(s.flat);
        p["spaceAggr"] = matrix_json(A);
        // `localStateSpace` of `[stateSpace, localStateSpace] = getStateSpace`.
        // One entry per stateful node, in stateful-node order, holding that
        // node's distinct local rows. A host with only the flat space can
        // return the second output as a SINGLE cell holding the whole space,
        // which is what `getStateSpace.m` does under `lang='cpp'`, and a caller
        // indexing it per node then reads the global space for every node.
        line::reg::Json loc = line::reg::Json::array();
        for (std::size_t f = 0; f < s.local.size(); ++f) loc.push_back(matrix_json(s.local[f]));
        p["localSpace"] = loc;
        // pi travels with the space because a row of the space is not an answer:
        // the pair (state, probability) is, and reading them from two invocations
        // would risk pairing one solve's states with another solve's law.
        p["pi"] = vector_json(d.pi);
        emit_analysis<T>("states", p, d.actualmethod, ctmc_meta<T>(d));
        return 0;
    }
    print_ctmc_banner<T>(d);
    // The per-node widths are printed because the flat row is only decodable
    // with them: they are where one node's block ends and the next begins.
    std::printf("NodeWidths");
    for (std::size_t f = 0; f < s.node_width.size(); ++f) std::printf(" %zu", s.node_width[f]);
    std::printf("\n");
    std::printf("%8s %12s   %s\n", "State", "Prob", "Detailed | Aggregate");
    for (std::size_t i = 0; i < s.flat.rows(); ++i) {
        std::printf("%8zu %12.6g  ", i + 1, line::num_traits<T>::to_double(d.pi[i]));
        for (std::size_t c = 0; c < s.flat.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(s.flat(i, c)));
        std::printf(" |");
        for (std::size_t c = 0; c < A.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(A(i, c)));
        std::printf("\n");
    }
    return 0;
}

/**
 * `-a sens`: `getSensitivityRanking` over the exponential service rates.
 *
 * THE REWARD IS STATED, NOT INFERRED. The reference makes the caller supply one,
 * and a model.json carries no reward function, so the CLI has to name the reward
 * it ranks against or the numbers mean nothing: it is the mean number of jobs at
 * the QUEUEING stations, which excludes the Source (whose column is the infinite
 * reservoir) and the Delay stations (whose population is think time, not work).
 *
 * Only an exponential service pair becomes a parameter. Perturbing any other
 * distribution would mean replacing it with an exponential of the new rate,
 * which changes the model's shape rather than one of its parameters; those pairs
 * are listed as skipped rather than silently differenced.
 */
template <class T>
int solve_ctmc_sens(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const line::ctmc::CtmcSolution<T> d = line::ctmc::solver_ctmc_analyzer(sn, opt);
    const line::Matrix<T> A = line::ctmc::ctmc_state_space_aggr(sn, d.chain.space);

    std::vector<bool> queueing(M, false);
    for (std::size_t i = 0; i < M; ++i)
        queueing[i] = sn.stations[i].nodetype != line::lang::NodeType::Source &&
                      sn.stations[i].sched != line::lang::SchedStrategy::INF &&
                      sn.stations[i].sched != line::lang::SchedStrategy::EXT;

    std::vector<T> reward(d.chain.space.size(), line::num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < reward.size(); ++s)
        for (std::size_t i = 0; i < M; ++i)
            if (queueing[i])
                for (std::size_t c = 0; c < K; ++c) reward[s] += A(s, i * K + c);

    std::vector<line::ctmc::CtmcSensParam<T> > params;
    std::vector<std::string> skipped;
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].nodetype == line::lang::NodeType::Source) continue;
        for (std::size_t c = 0; c < K; ++c) {
            if (sn.disabled[i][c]) continue;
            const double mu = line::num_traits<T>::to_double(sn.rates(i, c));
            if (!(mu > 0.0)) continue;
            const std::string nm =
                "mu(" + sn.stations[i].name + "," + sn.classes[c].name + ")";
            if (sn.service[i][c].type != line::lang::ProcessType::EXP) {
                skipped.push_back(nm);
                continue;
            }
            line::ctmc::CtmcSensParam<T> p;
            p.name = nm;
            p.value = mu;
            p.set = [i, c](line::qn::NetworkStruct<T>& s, double v) {
                s.service[i][c] = line::lang::Distrib<T>::exp_rate(line::num_traits<T>::from_double(v));
                s.refresh_rates();
            };
            params.push_back(p);
        }
    }
    if (params.empty())
        throw line::UnsupportedError(
            "-a sens found no exponential service rate to differentiate: every enabled "
            "(station, class) pair carries a non-exponential distribution, and replacing one with "
            "an exponential of the perturbed rate would change the model rather than a parameter "
            "of it");

    const std::vector<line::ctmc::CtmcSensRank<T> > rank =
        line::ctmc::solver_ctmc_sensitivity_ranking(sn, opt, params, reward);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "SensRanking";
        // THE REWARD IS STATED, NOT INFERRED, on this path too: a ranking is a
        // ranking against something, and a host that read these numbers without
        // knowing what was differentiated would be reading a sensitivity of an
        // unnamed functional.
        p["reward"] = "mean number of jobs at the queueing stations";
        line::reg::Json par = line::reg::Json::array(), val = line::reg::Json::array(),
                        S = line::reg::Json::array(), SS = line::reg::Json::array();
        for (std::size_t l = 0; l < rank.size(); ++l) {
            par.push_back(rank[l].parameter);
            val.push_back(rank[l].value);
            S.push_back(line::num_traits<T>::to_double(rank[l].S));
            // null, not NaN and not 0: a zero mean reward makes the scaled form
            // undefined, which is a property of the model, and JSON has no NaN.
            if (rank[l].scaled_valid)
                SS.push_back(line::num_traits<T>::to_double(rank[l].SS));
            else
                SS.push_back(line::reg::Json());
        }
        p["Parameter"] = par;
        p["Value"] = val;
        p["Sens"] = S;
        p["ScaledSens"] = SS;
        line::reg::Json sk = line::reg::Json::array();
        for (std::size_t l = 0; l < skipped.size(); ++l) sk.push_back(skipped[l]);
        p["Skipped"] = sk;
        emit_analysis<T>("sens", p, d.actualmethod, ctmc_meta<T>(d));
        return 0;
    }
    print_ctmc_banner<T>(d);
    std::printf("Reward mean number of jobs at the queueing stations\n");
    std::printf("%-28s %14s %16s %16s\n", "Parameter", "Value", "Sens", "ScaledSens");
    for (std::size_t l = 0; l < rank.size(); ++l) {
        if (rank[l].scaled_valid)
            std::printf("%-28s %14.6g %16.8g %16.8g\n", rank[l].parameter.c_str(), rank[l].value,
                        line::num_traits<T>::to_double(rank[l].S),
                        line::num_traits<T>::to_double(rank[l].SS));
        else
            // MATLAB reports NaN here; the reason is printed instead, since a
            // zero mean reward is a property of the model and not a failure.
            std::printf("%-28s %14.6g %16.8g %16s\n", rank[l].parameter.c_str(), rank[l].value,
                        line::num_traits<T>::to_double(rank[l].S), "undefined(E[r]=0)");
    }
    for (std::size_t l = 0; l < skipped.size(); ++l)
        std::printf("Skipped %s: service is not exponential\n", skipped[l].c_str());
    return 0;
}

/** `-a reward`: `getAvgReward`, the steady-state expectation of each reward. */
template <class T>
int solve_ctmc_reward(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    std::vector<std::string> names;
    const std::vector<T> r = line::ctmc::solver_ctmc_avg_reward(sn, opt, &names);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgReward";
        line::reg::Json nm = line::reg::Json::array();
        for (std::size_t l = 0; l < names.size(); ++l) nm.push_back(names[l]);
        p["Reward"] = nm;
        p["E"] = vector_json(r);
        // NO "method": this arm returns the expectations and no solved chain, so
        // there is no resolved method to report and none is invented.
        emit_analysis<T>("reward", p, std::string());
        return 0;
    }
    std::printf("SolverCTMC arith=%s rewards=%zu\n", line::num_traits<T>::name(), r.size());
    std::printf("%-28s %16s\n", "Reward", "E[r]");
    for (std::size_t l = 0; l < r.size(); ++l)
        std::printf("%-28s %16.10g\n", names[l].c_str(), line::num_traits<T>::to_double(r[l]));
    return 0;
}

/**
 * `-a reward-value`: `@@SolverCTMC/getRewardValueFunction`, V^k(s).
 *
 * A DIFFERENT OBJECT FROM BOTH `-a reward` AND `-a tranreward`. `-a reward`
 * returns one number per reward, the steady-state E[r]; `-a tranreward` the
 * expected rate along a horizon; this returns the VALUE FUNCTION of ONE named
 * reward -- the reward accumulated over k uniformized steps, from every state
 * of the chain, as a (Tmax+1 x nstates) matrix. It is the object a policy
 * evaluation reads, and it is indexed by state, not by station.
 *
 * `--reward-name` IS REQUIRED and not defaulted to the first declared reward:
 * the value functions of two rewards are different matrices, and labelling one
 * with the caller's question would be a wrong answer rather than a missing one.
 */
template <class T>
int solve_ctmc_reward_value(const line::qn::NetworkStruct<T>& sn,
                            const line::ctmc::CtmcOptions& opt, const Knobs& k) {
    if (k.reward_name.empty())
        throw line::InputError(
            "-a reward-value returns the value function of ONE reward and needs --reward-name to "
            "say which; -a reward returns the steady-state expectation of every declared reward");
    const line::ctmc::CtmcReward<T> rr = line::ctmc::solver_ctmc_reward(sn, opt);
    std::size_t which = rr.names.size();
    for (std::size_t l = 0; l < rr.names.size(); ++l)
        if (rr.names[l] == k.reward_name) which = l;
    if (which == rr.names.size()) {
        std::string avail;
        for (std::size_t l = 0; l < rr.names.size(); ++l)
            avail += (l ? ", " : "") + rr.names[l];
        throw line::InputError("--reward-name '" + k.reward_name +
                               "' is not declared by this model; it declares: " +
                               (avail.empty() ? std::string("(none)") : avail));
    }
    const line::Matrix<T>& V = rr.V[which];

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "RewardValueFunction";
        p["indexBase"] = 0;
        p["Reward"] = rr.names[which];
        p["steps"] = V.rows();
        p["states"] = V.cols();
        p["t"] = vector_json(rr.t);
        p["V"] = matrix_json(V);
        p["stateSpaceAggr"] = matrix_json(rr.state_space_aggr);
        emit_analysis<T>("rewardvalue", p, std::string());
        return 0;
    }
    std::printf("SolverCTMC arith=%s reward=%s steps=%zu states=%zu\n",
                line::num_traits<T>::name(), rr.names[which].c_str(), V.rows(), V.cols());
    std::printf("%-10s %-10s %20s\n", "Step", "State", "V");
    for (std::size_t i = 0; i < V.rows(); ++i)
        for (std::size_t j = 0; j < V.cols(); ++j)
            std::printf("%-10zu %-10zu %20.10g\n", i, j, line::num_traits<T>::to_double(V(i, j)));
    return 0;
}

/**
 * `-a tranreward`: `getTranReward`, E[r(X(t))] over the --tspan horizon.
 *
 * A DIFFERENT QUANTITY FROM `-a reward`, not a formatting of it. `-a reward`
 * returns the steady-state expectation, one number per reward; this returns the
 * expected reward RATE along the trajectory, which converges to that number but
 * is not it at any finite t. It is also not the accumulated reward `V`, which
 * the same header computes and which diverges -- the header calls that the
 * easiest mistake to make here, so the two are kept on separate flags.
 *
 * THE HORIZON IS REQUIRED, as it is for `-a tranprob`: E[r(X(t))] on an
 * unstated horizon is not a quantity, and the reference refuses an infinite one
 * rather than picking a bound.
 */
template <class T>
int solve_ctmc_tran_reward(const line::qn::NetworkStruct<T>& sn,
                           const line::ctmc::CtmcOptions& opt, const Knobs& k) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        (void)sn; (void)opt; (void)k;
        throw line::UnsupportedError(
            "-a tranreward integrates the forward equation, which needs transcendental "
            "arithmetic; rerun with --arith double or --arith real");
    } else {
        if (k.t1 < 0.0)
            throw line::InputError(
                "-a tranreward integrates E[r(X(t))] and needs a horizon: pass --tspan <t0>:<t1>");
        std::vector<std::string> names;
        std::vector<T> t;
        const std::vector<std::vector<T> > r = line::ctmc::solver_ctmc_tran_reward(
            sn, opt, line::num_traits<T>::from_double(k.t0),
            line::num_traits<T>::from_double(k.t1), &t, &names);
        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "TranReward";
            p["indexBase"] = 0;
            line::reg::Json nm = line::reg::Json::array();
            for (std::size_t l = 0; l < names.size(); ++l) nm.push_back(names[l]);
            p["Reward"] = nm;
            p["t"] = vector_json<T>(t);
            line::reg::Json e = line::reg::Json::array();
            for (std::size_t l = 0; l < r.size(); ++l) e.push_back(vector_json<T>(r[l]));
            p["E"] = e;
            emit_analysis<T>("tranreward", p, std::string());
            return 0;
        }
        std::printf("SolverCTMC arith=%s rewards=%zu points=%zu tspan=[%g,%g]\n",
                    line::num_traits<T>::name(), r.size(), t.size(), k.t0, k.t1);
        std::printf("%-16s", "t");
        for (std::size_t l = 0; l < names.size(); ++l) std::printf(" %16s", names[l].c_str());
        std::printf("\n");
        for (std::size_t i = 0; i < t.size(); ++i) {
            std::printf("%-16.10g", line::num_traits<T>::to_double(t[i]));
            for (std::size_t l = 0; l < r.size(); ++l)
                std::printf(" %16.10g", line::num_traits<T>::to_double(r[l][i]));
            std::printf("\n");
        }
        return 0;
    }
}

/**
 * `-a tranprob`: pi(t) over the --tspan horizon, labelled by the whole network
 * (`getTranProbSys` / `getTranProbSysAggr`) or by one node when `--node` names
 * it (`getTranProb` / `getTranProbAggr`).
 *
 * The forward equation is integrated ONCE and both label sets are taken off the
 * same `CtmcTransient`. Calling the (sn, opt, ...) overloads twice would solve
 * and integrate the chain twice for two views of one answer.
 */
template <class T>
int solve_ctmc_tranprob(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                        const Knobs& k) {
    if (k.t1 < 0.0)
        throw line::InputError(
            "-a tranprob integrates pi(t) and needs a horizon: pass --tspan <t0>:<t1>");
    const line::ctmc::CtmcTransient<T> tr = line::ctmc::solver_ctmc_transient_analyzer(
        sn, opt, line::num_traits<T>::from_double(k.t0), line::num_traits<T>::from_double(k.t1));
    const line::ctmc::CtmcTranProb<T> det =
        k.node ? line::ctmc::ctmc_get_tran_prob(sn, tr, k.node)
               : line::ctmc::ctmc_get_tran_prob_sys(sn, tr);
    const line::ctmc::CtmcTranProb<T> agg =
        k.node ? line::ctmc::ctmc_get_tran_prob_aggr(sn, tr, k.node)
               : line::ctmc::ctmc_get_tran_prob_sys_aggr(sn, tr);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranProb";
        p["indexBase"] = 0;
        p["scope"] = k.node ? sn.nodes[k.node - 1].name : std::string("(system)");
        if (k.node) p["node"] = k.node - 1;
        line::reg::Json span = line::reg::Json::array();
        span.push_back(k.t0);
        span.push_back(k.t1);
        // THE HORIZON IS PART OF THE ANSWER: pi(t) on an unstated span is not a
        // quantity, and a host that took the last row for "the" occupancy without
        // it would be quoting a time it does not know.
        p["tspan"] = span;
        p["t"] = vector_json(det.t);
        // BOTH VIEWS COME OFF ONE INTEGRATION, as on the readable path: the
        // detailed labels answer getTranProb(Sys) and the aggregate ones
        // getTranProb(Sys)Aggr, and a host asking for both would otherwise
        // integrate the same forward equation twice for one answer.
        p["labels"] = matrix_json(det.labels);
        p["labelsAggr"] = matrix_json(agg.labels);
        p["pit"] = matrix_json(det.pit);
        p["pitAggr"] = matrix_json(agg.pit);
        emit_analysis<T>("tranprob", p, tr.chain.actualmethod, ctmc_meta<T>(tr.chain));
        return 0;
    }
    print_ctmc_banner<T>(tr.chain);
    std::printf("TranProb times=%zu tspan=%g:%g scope=%s\n", det.t.size(), k.t0, k.t1,
                k.node ? sn.nodes[k.node - 1].name.c_str() : "(system)");
    // The labels come FIRST and the occupancy after, because pi(t) is a row per
    // time over columns that mean nothing until the state they index is named.
    std::printf("%8s   %s\n", "State", "Detailed | Aggregate");
    for (std::size_t s = 0; s < det.labels.rows(); ++s) {
        std::printf("%8zu  ", s + 1);
        for (std::size_t c = 0; c < det.labels.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(det.labels(s, c)));
        std::printf(" |");
        for (std::size_t c = 0; c < agg.labels.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(agg.labels(s, c)));
        std::printf("\n");
    }
    std::printf("%14s", "Time");
    for (std::size_t s = 0; s < det.pit.cols(); ++s) std::printf(" %12zu", s + 1);
    std::printf("\n");
    for (std::size_t i = 0; i < det.t.size(); ++i) {
        std::printf("%14.8g", line::num_traits<T>::to_double(det.t[i]));
        for (std::size_t s = 0; s < det.pit.cols(); ++s)
            std::printf(" %12.6g", line::num_traits<T>::to_double(det.pit(i, s)));
        std::printf("\n");
    }
    return 0;
}

/**
 * `-s ctmc -a tran`: `getTranAvg`, the transient MEANS Q(t), U(t) and X(t) over
 * the --tspan horizon.
 *
 * NOT `-a tranprob`, AND THE TWO ARE NOT REDUCIBLE TO ONE ANOTHER FOR A CALLER.
 * `tranprob` sends pi(t) with the labels that index it, from which a host COULD
 * form these means -- and that is exactly the computation that must not happen
 * in a host: the utilization is not a linear functional of the labels (the PS
 * and DPS shares divide by the state's own total, and every other discipline
 * takes min(n_k, c)/c with the reference's warning attached), so a host folding
 * pi(t) itself would be reimplementing `solver_ctmc_transient_analyzer`'s
 * discipline switch and would diverge from it silently. The analyzer already
 * computes all three trajectories on the way to pi(t); this arm reports them.
 *
 * The payload is the SAME `TranAvgTable` the fluid and MAM transients emit, key
 * for key, so one host reader serves every solver that answers `-a tran`.
 */
template <class T>
int solve_ctmc_tran(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                    const Knobs& k) {
    if (k.t1 < 0.0)
        throw line::InputError(
            "-a tran integrates the forward equation and needs a horizon: pass --tspan <t0>:<t1>");
    const line::ctmc::CtmcTransient<T> tr = line::ctmc::solver_ctmc_transient_analyzer(
        sn, opt, line::num_traits<T>::from_double(k.t0), line::num_traits<T>::from_double(k.t1));
    const std::size_t M = sn.nstations, K = sn.nclasses, nt = tr.t.size();

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranAvgTable";
        p["indexBase"] = 0;
        p["t0"] = k.t0;
        p["t1"] = k.t1;
        line::reg::Json ts = line::reg::Json::array();
        for (std::size_t j = 0; j < nt; ++j) ts.push_back(line::num_traits<T>::to_double(tr.t[j]));
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < K; ++c) {
                // A DISABLED PAIR IS OMITTED, not sent as zeros: the reference
                // leaves its cell empty and the host turns an absent curve into
                // the disabled handle's NaN. Zeros would read as a station that
                // is genuinely idle for that class.
                if (sn.disabled[i][c]) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = ts;
                line::reg::Json q = line::reg::Json::array(), u = line::reg::Json::array(),
                                 x = line::reg::Json::array();
                for (std::size_t j = 0; j < nt; ++j) {
                    q.push_back(line::num_traits<T>::to_double(tr.QNt[i][c][j]));
                    u.push_back(line::num_traits<T>::to_double(tr.UNt[i][c][j]));
                    x.push_back(line::num_traits<T>::to_double(tr.TNt[i][c][j]));
                }
                e["QLen"] = q;
                e["Util"] = u;
                e["Tput"] = x;
                arr.push_back(e);
            }
        p["curves"] = arr;
        emit_analysis<T>("tran", p, tr.chain.actualmethod, ctmc_meta<T>(tr.chain));
        return 0;
    }
    print_ctmc_banner<T>(tr.chain);
    std::printf("TranAvg times=%zu tspan=%g:%g\n", nt, k.t0, k.t1);
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "Time", "QLen", "Util",
                "Tput");
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            if (sn.disabled[i][c]) continue;
            for (std::size_t j = 0; j < nt; ++j)
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n",
                            sn.stations[i].name.c_str(), sn.classes[c].name.c_str(),
                            line::num_traits<T>::to_double(tr.t[j]),
                            line::num_traits<T>::to_double(tr.QNt[i][c][j]),
                            line::num_traits<T>::to_double(tr.UNt[i][c][j]),
                            line::num_traits<T>::to_double(tr.TNt[i][c][j]));
        }
    return 0;
}

/**
 * `-a sample`: `sampleSys` and `sampleSysAggr`, one marked trajectory; with
 * `--node`, also that node's own block (`sample`) and per-class counts
 * (`sampleAggr`), which is the view MATLAB's per-node sampler returns.
 */
template <class T>
int solve_ctmc_sample(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt,
                      const Knobs& k) {
    // `--events` NAMES THIS NUMBER, `--samples` only stands in for it. The
    // JAR keeps the two apart because a sampled trajectory is walked for a
    // number of EVENTS while `--samples` is a solver-wide run length, so a
    // caller who set the run length and then asked for a trajectory would
    // silently get one of that length. Both are honoured, --events first,
    // and the reference default of 1000 stands when neither is given.
    const std::size_t nevents = k.events ? k.events : (k.samples ? k.samples : 1000);
    const unsigned long seed = k.seed ? k.seed : 23000;
    const line::ctmc::CtmcSamplePath<T> path =
        line::ctmc::solver_ctmc_sample_sys<T>(sn, opt, nevents, seed);
    const line::Matrix<T> A = line::ctmc::solver_ctmc_sample_sys_aggr(sn, path);
    line::Matrix<T> L, LA;
    if (k.node) {
        L = line::ctmc::solver_ctmc_sample(sn, path, k.node);
        LA = line::ctmc::solver_ctmc_sample_aggr(sn, path, k.node);
    }
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "SamplePath";
        p["indexBase"] = 0;
        // The seed and the requested length are part of the ANSWER, as on the SSA
        // path: two runs are the same trace only if both are stated.
        p["events"] = nevents;
        p["seed"] = seed;
        p["drawn"] = path.state.size();
        p["scope"] = k.node ? sn.nodes[k.node - 1].name : std::string("(system)");
        if (k.node) p["node"] = k.node - 1;
        p["t"] = vector_json(path.t);
        p["state"] = index_json(path.state);
        line::reg::Json ev = line::reg::Json::array();
        for (std::size_t i = 0; i < path.state.size(); ++i) {
            // null where the readable table prints "absorb": the walk left no
            // state, so there is no synchronization index, and any integer here
            // would name an event that did not fire.
            if (i < path.event.size() && path.event[i] != static_cast<std::size_t>(-1))
                ev.push_back(path.event[i]);
            else
                ev.push_back(line::reg::Json());
        }
        p["event"] = ev;
        p["sysAggr"] = matrix_json(A);
        // THE STATE SPACE TRAVELS WITH THE TRAJECTORY, so a host can turn the
        // visited indices into the states themselves without enumerating the
        // chain a second time in another process -- which would also be a second
        // chance for the two enumerations to disagree while looking paired.
        {
            const line::ctmc::CtmcStateSpace<T> s =
                line::ctmc::ctmc_get_state_space(sn, path.chain);
            p["space"] = matrix_json(s.flat);
            p["NodeWidths"] = index_json(s.node_width);
        }
        if (k.node) {
            p["nodeState"] = matrix_json(L);
            p["nodeAggr"] = matrix_json(LA);
        }
        emit_analysis<T>("sample", p, path.chain.actualmethod, ctmc_meta<T>(path.chain));
        return 0;
    }
    // The seed and the requested length are part of the ANSWER, as on the SSA
    // path: two runs are the same trace only if both are stated.
    std::printf("SolverCTMC arith=%s states=%zu events=%zu seed=%lu drawn=%zu scope=%s\n",
                line::num_traits<T>::name(), path.chain.chain.space.size(), nevents, seed,
                path.state.size(), k.node ? sn.nodes[k.node - 1].name.c_str() : "(system)");
    std::printf("%14s %8s %8s   %s\n", "Time", "State", "Event",
                k.node ? "SysAggregate | NodeState | NodeAggregate" : "SysAggregate");
    for (std::size_t i = 0; i < path.state.size(); ++i) {
        const std::size_t ev = i < path.event.size() ? path.event[i] : static_cast<std::size_t>(-1);
        std::printf("%14.8g %8zu ", line::num_traits<T>::to_double(path.t[i]), path.state[i] + 1);
        if (ev == static_cast<std::size_t>(-1))
            std::printf("%8s  ", "absorb");
        else
            std::printf("%8zu  ", ev + 1);
        for (std::size_t c = 0; c < A.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(A(i, c)));
        if (k.node) {
            std::printf(" |");
            for (std::size_t c = 0; c < L.cols(); ++c)
                std::printf(" %g", line::num_traits<T>::to_double(L(i, c)));
            std::printf(" |");
            for (std::size_t c = 0; c < LA.cols(); ++c)
                std::printf(" %g", line::num_traits<T>::to_double(LA(i, c)));
        }
        std::printf("\n");
    }
    return 0;
}

/** `-a cdf`: `getCdfRespT` per (station, class), and `getCdfSysRespT` per chain. */
template <class T>
int solve_ctmc_cdf(const line::qn::NetworkStruct<T>& sn, const line::ctmc::CtmcOptions& opt) {
    const std::vector<std::vector<line::ctmc::CdfCurve<T> > > RD =
        line::ctmc::solver_ctmc_cdf_respt(sn, opt);
    const std::vector<line::ctmc::CdfCurve<T> > RS = line::ctmc::solver_ctmc_cdf_sys_respt(sn, opt);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "CdfRespT";
        p["chains"] = sn.nchains;
        // ONE OBJECT PER CURVE rather than four parallel columns: the grids are
        // per-pair and of different lengths, so a column-oriented form would have
        // to repeat the station and the class on every sample and leave the host
        // to re-group them.
        line::reg::Json rd = line::reg::Json::array();
        for (std::size_t i = 0; i < RD.size(); ++i)
            for (std::size_t c = 0; c < RD[i].size(); ++c) {
                // An empty curve is not a degenerate one: the chain never visits
                // that pair, so there is no arrival event to condition on, and it
                // is OMITTED rather than sent as a flat zero law.
                if (RD[i][c].empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = vector_json(RD[i][c].t);
                e["F"] = vector_json(RD[i][c].F);
                rd.push_back(e);
            }
        p["respt"] = rd;
        line::reg::Json rs = line::reg::Json::array();
        for (std::size_t c = 0; c < RS.size(); ++c) {
            if (RS[c].empty()) continue;
            line::reg::Json e = line::reg::Json::object();
            e["chain"] = c;
            e["t"] = vector_json(RS[c].t);
            e["F"] = vector_json(RS[c].F);
            rs.push_back(e);
        }
        p["sysrespt"] = rs;
        p["indexBase"] = 0;
        // NO "method", as on the reward arm: the response-time laws are computed
        // from tagged chains this function does not return, so there is no
        // resolved method to report and the requested one is not it.
        emit_analysis<T>("cdf", p, std::string());
        return 0;
    }
    std::printf("SolverCTMC arith=%s chains=%zu\n", line::num_traits<T>::name(), sn.nchains);
    std::printf("%-16s %-14s %14s %14s\n", "Station", "JobClass", "Time", "F(t)");
    for (std::size_t i = 0; i < RD.size(); ++i)
        for (std::size_t c = 0; c < RD[i].size(); ++c) {
            // An empty curve is not a degenerate one: the chain never visits
            // that pair, so there is no arrival event to condition on.
            if (RD[i][c].empty()) continue;
            for (std::size_t j = 0; j < RD[i][c].t.size(); ++j)
                std::printf("%-16s %-14s %14.8g %14.10g\n", sn.stations[i].name.c_str(),
                            sn.classes[c].name.c_str(),
                            line::num_traits<T>::to_double(RD[i][c].t[j]),
                            line::num_traits<T>::to_double(RD[i][c].F[j]));
        }
    std::printf("%-16s %-14s %14s %14s\n", "System", "Chain", "Time", "F(t)");
    for (std::size_t c = 0; c < RS.size(); ++c) {
        if (RS[c].empty()) continue;
        for (std::size_t j = 0; j < RS[c].t.size(); ++j)
            std::printf("%-16s %-14zu %14.8g %14.10g\n", "(system)", c + 1,
                        line::num_traits<T>::to_double(RS[c].t[j]),
                        line::num_traits<T>::to_double(RS[c].F[j]));
    }
    return 0;
}

/**
 * A `--passage-from`/`--passage-into` spec as a Matrix<double>: either a flat
 * comma list ("3,5", one row the resolver reads as 1-based indices when every
 * entry is one) or semicolon-separated state rows ("0,2;1,1"), which resolve
 * against the enumerated space by content.
 */
inline line::Matrix<double> parse_passage_set(const std::string& spec, const char* flag) {
    if (spec.empty()) return line::Matrix<double>(0, 0);
    std::vector<std::vector<double>> rows;
    std::size_t pos = 0;
    while (pos <= spec.size()) {
        std::size_t semi = spec.find(';', pos);
        if (semi == std::string::npos) semi = spec.size();
        std::string rowtxt = spec.substr(pos, semi - pos);
        std::vector<double> row;
        std::size_t p2 = 0;
        while (p2 <= rowtxt.size()) {
            std::size_t comma = rowtxt.find(',', p2);
            if (comma == std::string::npos) comma = rowtxt.size();
            std::string cell = rowtxt.substr(p2, comma - p2);
            if (!cell.empty()) {
                char* endp = 0;
                const double v = std::strtod(cell.c_str(), &endp);
                if (endp == cell.c_str() || *endp != '\0')
                    throw line::InputError(std::string(flag) + ": '" + cell +
                                           "' is not a number");
                row.push_back(v);
            }
            p2 = comma + 1;
        }
        if (!row.empty()) rows.push_back(row);
        pos = semi + 1;
    }
    if (rows.empty()) return line::Matrix<double>(0, 0);
    for (std::size_t i = 1; i < rows.size(); ++i)
        if (rows[i].size() != rows[0].size())
            throw line::InputError(std::string(flag) +
                                   ": every state row must have the same width");
    line::Matrix<double> out(rows.size(), rows[0].size());
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (std::size_t j = 0; j < rows[i].size(); ++j) out(i, j) = rows[i][j];
    return out;
}

/**
 * `-s ctmc -a firstpasst`: `@@SolverCTMC/getCdfFirstPassT(A, B)`, the first
 * passage time between two state sets the caller names. `--passage-into` is
 * required; an empty `--passage-from` starts from the conditional stationary
 * law on the complement of the target, as the reference does.
 */
template <class T>
int solve_ctmc_firstpasst(const line::qn::NetworkStruct<T>& sn,
                          const line::ctmc::CtmcOptions& opt, const Knobs& k) {
    if (k.passage_into.empty())
        throw line::InputError(
            "-a firstpasst times the passage INTO a state set and needs --passage-into; name it "
            "as 1-based rows of the state space ('3,5') or as state rows ('0,2;1,1')");
    const line::Matrix<double> A = parse_passage_set(k.passage_from, "--passage-from");
    const line::Matrix<double> B = parse_passage_set(k.passage_into, "--passage-into");
    const std::string method = k.passage_method.empty() ? "expm" : k.passage_method;

    const line::ctmc::CtmcFirstPassage fp =
        line::ctmc::ctmc_cdf_firstpasst<T>(sn, opt, A, B, method);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "CdfFirstPassT";
        p["indexBase"] = 0;
        p["t"] = line::reg::Json(fp.t);
        p["F"] = line::reg::Json(fp.F);
        p["f"] = line::reg::Json(fp.f);
        line::reg::Json src = line::reg::Json::array(), tgt = line::reg::Json::array();
        for (std::size_t i : fp.source) src.push_back(static_cast<double>(i));
        for (std::size_t i : fp.target) tgt.push_back(static_cast<double>(i));
        p["source"] = src;
        p["target"] = tgt;
        emit_analysis<T>("firstpasst", p, method);
        return 0;
    }
    std::printf("SolverCTMC arith=%s method=%s getCdfFirstPassT\n",
                line::num_traits<T>::name(), method.c_str());
    std::printf("source states: %zu%s, target states: %zu\n", fp.source.size(),
                fp.source.empty() ? " (conditional stationary law)" : "", fp.target.size());
    std::printf("%14s %14s %14s\n", "Time", "F(t)", "f(t)");
    for (std::size_t j = 0; j < fp.t.size(); j += 111)
        std::printf("%14.8g %14.10g %14.10g\n", fp.t[j], fp.F[j], fp.f[j]);
    std::printf("%14.8g %14.10g %14.10g\n", fp.t.back(), fp.F.back(), fp.f.back());
    return 0;
}

/**
 * `-s ctmc -a firstpasstmom`: `@@SolverCTMC/getFirstPassTMoments(A, B, nmax)`,
 * the moments of the same passage `firstpasst` gives the curve of.
 *
 * These are EXACT and cost one linear solve per order, so a caller who wants a
 * variance or a skewness should ask for them here rather than integrate the
 * truncated curve the other arm returns. `--passage-orders` is the nmax; it
 * defaults to the reference's 3.
 */
template <class T>
int solve_ctmc_firstpasst_moments(const line::qn::NetworkStruct<T>& sn,
                                  const line::ctmc::CtmcOptions& opt, const Knobs& k) {
    if (k.passage_into.empty())
        throw line::InputError(
            "-a firstpasstmom times the passage INTO a state set and needs --passage-into; name "
            "it as 1-based rows of the state space ('3,5') or as state rows ('0,2;1,1')");
    const line::Matrix<double> A = parse_passage_set(k.passage_from, "--passage-from");
    const line::Matrix<double> B = parse_passage_set(k.passage_into, "--passage-into");
    const std::size_t nmax = (k.passage_orders > 0) ? k.passage_orders : 3;

    const line::ctmc::CtmcFirstPassageMoments<T> fm =
        line::ctmc::ctmc_firstpasst_moments<T>(sn, opt, A, B, nmax);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "FirstPassTMoments";
        p["indexBase"] = 0;
        line::reg::Json m = line::reg::Json::array();
        for (std::size_t i = 0; i < fm.m.size(); ++i)
            m.push_back(line::num_traits<T>::to_double(fm.m[i]));
        p["m"] = m;
        line::reg::Json mall = line::reg::Json::array();
        for (std::size_t i = 0; i < fm.mall.rows(); ++i) {
            line::reg::Json row = line::reg::Json::array();
            for (std::size_t j = 0; j < fm.mall.cols(); ++j)
                row.push_back(line::num_traits<T>::to_double(fm.mall(i, j)));
            mall.push_back(row);
        }
        p["mall"] = mall;
        line::reg::Json src = line::reg::Json::array(), tgt = line::reg::Json::array();
        for (std::size_t i : fm.source) src.push_back(static_cast<double>(i));
        for (std::size_t i : fm.target) tgt.push_back(static_cast<double>(i));
        p["source"] = src;
        p["target"] = tgt;
        emit_analysis<T>("firstpasstmom", p, "moments");
        return 0;
    }
    std::printf("SolverCTMC arith=%s getFirstPassTMoments\n", line::num_traits<T>::name());
    std::printf("source states: %zu%s, target states: %zu\n", fm.source.size(),
                fm.source.empty() ? " (conditional stationary law)" : "", fm.target.size());
    std::printf("%8s %20s\n", "Order", "Moment");
    for (std::size_t i = 0; i < fm.m.size(); ++i)
        std::printf("%8zu %20.10g\n", i + 1, line::num_traits<T>::to_double(fm.m[i]));
    return 0;
}

/**
 * Solve a Network model.json with SolverCTMC and print what `-a` asked for.
 *
 * EXACT, and the only ported solver that is exact on a non-product-form model:
 * it enumerates the state space and solves pi Q = 0, so the numbers are the
 * chain's own and not an approximation of them. The price is the state space,
 * which is why an OPEN model needs `--cutoff`: without a bound on the open
 * population the chain is infinite. The banner reports the cutoff that was used,
 * because a truncated chain's answer is not the model's answer without it.
 *
 * Every arithmetic backend runs the stationary analyses: the generator assembly
 * and the stationary solve are field operations throughout, so `--arith exact`
 * returns the exact rational stationary law of a chain with rational rates. The
 * transient ones refuse by name under exact, since a forward integration, an
 * exponential clock and a matrix exponential are all transcendental.
 */
template <class T>
int solve_model_ctmc(const std::string& file, const Knobs& k, const std::string& analysis) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ctmc::CtmcOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (k.cutoff >= 0.0) opt.cutoff = k.cutoff;
    opt.cutoff_mat = k.cutoff_mat;
    opt.force = k.force;
    if (k.timestep > 0.0) opt.timestep = k.timestep;  // `--timestep`, the fixed output grid
    // `--transient-method` and its two tolerances, `options.config` of the
    // reference's transient analyzer. The analyzer validates the name.
    if (!k.transient_method.empty()) opt.transient_method = k.transient_method;
    if (k.fau_epsilon > 0.0) opt.fau_epsilon = k.fau_epsilon;
    if (k.fau_delta >= 0.0) opt.fau_delta = k.fau_delta;
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    // The two methods that never build the generator are served by their own
    // analyzers, ahead of every state-space path below. The dispatcher has
    // already refused every analysis but `avg` for them, since neither produces
    // a state space, a filtration or a trajectory to answer one from.
    if (opt.method == "mdd") {
        line::mdd::MddMcdOptions mcdopt;
        if (k.mdd_tol > 0.0) mcdopt.tol = k.mdd_tol;
        if (k.mdd_maxiter > 0) mcdopt.maxiter = k.mdd_maxiter;
        return solve_ctmc_mdd_avg<T>(sn, opt, mcdopt);
    }
    if (opt.method == "cftp" || opt.method == "cftp.approx") {
        line::ctmc::CtmcCftpOptions cftpopt;
        if (k.samples) cftpopt.samples = k.samples;
        if (k.seed) cftpopt.seed = k.seed;
        return solve_ctmc_cftp_avg<T>(sn, opt, cftpopt);
    }

    // The field-arithmetic analyses: every step from the generator to the answer
    // is a ring operation, so exact returns the exact rational quantity.
    if (analysis == "avg") return solve_ctmc_avg<T>(sn, opt);
    if (analysis == "prob") return solve_ctmc_prob<T>(sn, opt, k);
    if (analysis == "gen") return solve_ctmc_gen<T>(sn, opt);
    if (analysis == "states") return solve_ctmc_states<T>(sn, opt);
    if (analysis == "sens") return solve_ctmc_sens<T>(sn, opt);
    if (analysis == "reward") return solve_ctmc_reward<T>(sn, opt);
    if (analysis == "rewardvalue") return solve_ctmc_reward_value<T>(sn, opt, k);

    // The rest integrate a forward equation, draw an exponential clock or take a
    // matrix exponential, none of which exists in a rational field. Their
    // static_asserts are behind if-constexpr so the refusal is a message rather
    // than a compile error in the exact instantiation.
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the -a " + analysis +
            " analysis integrates the forward equation, draws exponential clocks or takes a matrix "
            "exponential, none of which exists in exact rational arithmetic; rerun with --arith "
            "double or --arith real");
    } else {
        if (analysis == "tran") return solve_ctmc_tran<T>(sn, opt, k);
        if (analysis == "tranprob") return solve_ctmc_tranprob<T>(sn, opt, k);
        if (analysis == "tranreward") return solve_ctmc_tran_reward<T>(sn, opt, k);
        if (analysis == "sample") return solve_ctmc_sample<T>(sn, opt, k);
        if (analysis == "firstpasst") return solve_ctmc_firstpasst<T>(sn, opt, k);
        if (analysis == "firstpasstmom") return solve_ctmc_firstpasst_moments<T>(sn, opt, k);
        return solve_ctmc_cdf<T>(sn, opt);  // the dispatcher admitted no other name
    }
}

/**
 * `-s ssa -a prob`: the four SolverSSA probability queries over the model's
 * DEFAULT INITIAL STATE, the same state `-s ctmc -a prob` reports on.
 *
 * THE PAIR IS THE POINT. The CTMC answer is the stationary law of the chain and
 * this one is a time average of a finite sample path, so running both on a model
 * small enough for the chain measures the simulation error directly instead of
 * inferring it. The banner therefore carries the run length and the seed, as
 * every simulated number on this CLI does.
 *
 * `seen` TRAVELS WITH EACH PROBABILITY because a zero here has two meanings: the
 * path visited the state and left immediately, or it never got there at all. The
 * reference warns on the second; a machine-readable answer has to carry the
 * distinction rather than print it.
 */
template <class T>
int solve_model_ssa_prob(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    line::ssa::SsaSerialOptions opt;
    // `serial`, whatever `-m` said, and for the reference's own reason: the NRM
    // simulates per-(node, class, phase) counts rather than the state ENCODING,
    // so it has no row to compare a requested state against.
    // `@@SolverSSA/getProb.m` rewrites the method the same way.
    opt.method = "serial";
    if (k.samples) opt.samples = k.samples;
    if (k.seed) opt.seed = k.seed;
    if (k.warmupfrac >= 0.0) opt.warmupfrac = k.warmupfrac;
    if (k.cutoff >= 0.0) opt.cutoff = k.cutoff;
    const line::ssa::SsaProbReport r = line::ssa::solver_ssa_prob<T>(sn, opt);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbAggr";
        p["indexBase"] = 0;
        p["samples"] = r.samples;
        p["seed"] = r.seed;
        p["ProbSys"] = r.sys.prob;
        p["ProbSysAggr"] = r.sys_aggr.prob;
        p["ProbSysSeen"] = r.sys.seen;
        p["ProbSysAggrSeen"] = r.sys_aggr.seen;
        line::reg::Json st = line::reg::Json::array(), pm = line::reg::Json::array(),
                        pa = line::reg::Json::array(), sm = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            st.push_back(sn.stations[i].name);
            pm.push_back(r.marg[i].prob);
            pa.push_back(r.aggr[i].prob);
            sm.push_back(r.marg[i].seen);
        }
        p["Station"] = st;
        p["Prob"] = pm;
        p["ProbAggr"] = pa;
        p["Seen"] = sm;
        emit_analysis<T>("prob", p, "serial");
        return 0;
    }
    std::printf("SolverSSA arith=%s method=serial samples=%zu seed=%lu time=%.6g\n",
                line::num_traits<T>::name(), r.samples, r.seed, r.simulated_time);
    std::printf("ProbSys      = %.8g%s\n", r.sys.prob, r.sys.seen ? "" : "  (state never visited)");
    std::printf("ProbSysAggr  = %.8g%s\n", r.sys_aggr.prob,
                r.sys_aggr.seen ? "" : "  (state never visited)");
    std::printf("%-20s %14s %14s\n", "Station", "Prob", "ProbAggr");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        std::printf("%-20s %14.8g %14.8g\n", sn.stations[i].name.c_str(), r.marg[i].prob,
                    r.aggr[i].prob);
    return 0;
}

/**
 * `-s ssa -a sample`: `sampleSys` and `sampleSysAggr`, one simulated trajectory;
 * with `--node`, also that node's own block (`sample`) and per-class counts
 * (`sampleAggr`).
 *
 * The SAME shape `-s ctmc -a sample` emits, deliberately: the CTMC sampler walks
 * the jump chain of an enumerated generator and this one walks the network's own
 * encoding, and a host that can read one trajectory should be able to read the
 * other. The event column indexes the synchronization list, so the two are
 * comparable only within a solver -- which is why it is printed and not
 * interpreted here.
 */
template <class T>
int solve_model_ssa_sample(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    line::ssa::SsaSerialOptions opt;
    opt.method = "serial";
    opt.samples = k.events ? k.events : (k.samples ? k.samples : 1000);  // see -a sample above
    if (k.seed) opt.seed = k.seed;
    if (k.warmupfrac >= 0.0) opt.warmupfrac = k.warmupfrac;
    if (k.cutoff >= 0.0) opt.cutoff = k.cutoff;
    const line::ssa::SsaSerialSolution<T> sim = line::ssa::solver_ssa_serial_analyzer(sn, opt);
    const line::ssa::SsaSamplePath<T> sys = line::ssa::ssa_sample_sys(sn, sim.run);
    line::ssa::SsaSamplePath<T> nodep;
    if (k.node) nodep = line::ssa::ssa_sample_node(sn, sim.run, k.node);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "SamplePath";
        p["indexBase"] = 0;
        p["events"] = opt.samples;
        p["seed"] = sys.seed;
        p["drawn"] = sys.t.size();
        p["scope"] = k.node ? sn.nodes[k.node - 1].name : std::string("(system)");
        if (k.node) p["node"] = k.node - 1;
        p["t"] = vector_json(sys.t);
        p["event"] = index_json(sys.event);
        p["state"] = matrix_json(sys.state);
        // THE BLOCK BOUNDARIES TRAVEL, as under `-s ctmc -a sample`. A row of
        // `state` is the stateful nodes' local encodings laid end to end, and
        // `sampleSys` reports them one node at a time; a host that had to guess
        // the widths would cut the row at the wrong columns on any model with a
        // phase-type service, and read a phase index as a job count.
        {
            line::reg::Json w = line::reg::Json::array();
            for (std::size_t f = 0; f < sn.stateful_nodes.size(); ++f)
                w.push_back(sim.run.space.empty() ? 0 : sim.run.space[0].local[f].size());
            p["NodeWidths"] = w;
        }
        p["sysAggr"] = matrix_json(sys.aggr);
        if (k.node) {
            p["nodeState"] = matrix_json(nodep.state);
            p["nodeAggr"] = matrix_json(nodep.aggr);
        }
        emit_analysis<T>("sample", p, "serial");
        return 0;
    }
    std::printf("SolverSSA arith=%s method=serial events=%zu seed=%lu drawn=%zu scope=%s\n",
                line::num_traits<T>::name(), opt.samples, sys.seed, sys.t.size(),
                k.node ? sn.nodes[k.node - 1].name.c_str() : "(system)");
    std::printf("%14s %8s   %s\n", "Time", "Event",
                k.node ? "SysAggregate | NodeState | NodeAggregate" : "SysAggregate");
    for (std::size_t i = 0; i < sys.t.size(); ++i) {
        std::printf("%14.8g %8zu  ", sys.t[i], sys.event[i]);
        for (std::size_t c = 0; c < sys.aggr.cols(); ++c)
            std::printf(" %g", line::num_traits<T>::to_double(sys.aggr(i, c)));
        if (k.node) {
            std::printf(" |");
            for (std::size_t c = 0; c < nodep.state.cols(); ++c)
                std::printf(" %g", line::num_traits<T>::to_double(nodep.state(i, c)));
            std::printf(" |");
            for (std::size_t c = 0; c < nodep.aggr.cols(); ++c)
                std::printf(" %g", line::num_traits<T>::to_double(nodep.aggr(i, c)));
        }
        std::printf("\n");
    }
    return 0;
}

/**
 * A `Matrix<double>` read as this arithmetic's matrix.
 *
 * SSA and Fluid return plain-double solutions whatever `T` the CLI was asked
 * for, so every shared routine that takes a `Matrix<T>` -- the residence-time
 * conversion, the chain aggregation -- needs this one lift. Both are refused
 * outside `--arith double` anyway, so it is a type bridge and not a precision
 * claim.
 */
template <class T>
line::Matrix<T> to_matrix(const line::Matrix<double>& m) {
    line::Matrix<T> out(m.rows(), m.cols(), line::num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            out(i, j) = line::num_traits<T>::from_double(m(i, j));
    return out;
}

/**
 * Solve a Network model.json with SolverSSA and print the same table.
 *
 * A SIMULATION: its numbers carry Monte Carlo error, so the row is compared
 * against the other codebases' SSA rows and not against an exact solver's.
 *
 * `--samples` and `--seed` set the run length and the stream; without them the
 * defaults are 10000 firings and seed 23000, which on a two-station model is
 * roughly 3300 job cycles and lands a few percent from the analytical answer.
 * Diffing THAT against an exact solver reads as a defect and is not one: a
 * measured -3.12% at 1e4 on an M/M/1 falls to +0.06% at 2.56e6. The banner
 * therefore carries both numbers, so a row quoting an SSA figure carries the
 * conditions that produced it. Double only, refused by name in the dispatcher.
 *
 * `-m` REACHES THE SOLVER'S OWN DISPATCHER, `ssa::solver_ssa`, and not one
 * engine's entry. Calling `solver_ssa_nrm_analyzer` here would run the NRM
 * whatever `-m` said, so `-m serial` would silently answer with a different
 * estimator than the one asked for -- and the three names the NRM entry cannot
 * serve (`serial`, `para`, `parallel`) are all methods the library honours.
 */
template <class T>
int solve_model_ssa(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::ssa::SsaOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.samples) opt.samples = k.samples;
    if (k.seed) opt.seed = k.seed;
    if (k.warmupfrac >= 0.0) opt.warmupfrac = k.warmupfrac;
    // The cache write-back is COLLECTED HERE, not only on the `-a node` path:
    // the realized hit and miss shares are what the sample path measured, and
    // `-a avg` is the arm MATLAB's lang='cpp' bridge calls. Without it
    // `CPPLINE.restoreCacheResults` found no block, cleared the Cache node and
    // refreshed the visits back to link()'s offered 1/2-1/2.
    std::vector<line::ssa::SsaCacheRatio> cacheratio;
    const line::ssa::SsaSolution r = line::ssa::solver_ssa(net.get_struct(), opt, &cacheratio);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    // The seed and the sample count are part of the ANSWER, not of the
    // invocation: two runs of a simulation are the same measurement only if
    // both are stated, so the banner carries them and a parity row that quotes
    // an SSA number carries them with it.
    //
    // `r.method` is the engine that ACTUALLY ran, never `k.method`: `-m
    // parallel` on an NRM-eligible model reports `nrm`, because that is what
    // produced the numbers below and the banner may not claim otherwise.
    std::printf("SolverSSA arith=%s method=%s type=%s samples=%zu seed=%lu time=%.6g\n",
                line::num_traits<T>::name(), r.method.c_str(),
                line::util::method_type("SSA", r.method).c_str(), r.samples, opt.seed,
                r.simulated_time);
    // ResidT IS NOT RespT UNLESS EVERY STATION IS VISITED ONCE PER CYCLE.
    // `sn_get_residt_from_respt` is the reference's own per-visit -> per-job
    // conversion and a pure function of `sn` and RN, so a solver that reports no
    // residence time of its own still owes the caller this one: reporting RespT
    // in its place was a factor of 3 out on sdroute_closed and 17 on Queue1 of
    // init_state_ps, both multi-visit closed models.
    //
    // TAKEN ON THE MEASURED CACHE SPLIT, not on the offered one: the visits this
    // conversion divides by are a function of the routing, and a cache's routing
    // is a RESULT. On tut06_cache_lru_zipf the base struct still carried
    // `link()`'s even hit/miss share, so both classes came back at exactly half
    // their response time (0.1 and 0.5 against 0.16475 and 0.17625) -- a number
    // that is not the residence time of any model.
    const line::qn::NetworkStruct<T> snw = line::ssa::sn_with_ssa_cache_split<T>(sn, cacheratio);
    const line::Matrix<T> WN = line::mva::sn_get_residt_from_respt<T>(snw, to_matrix<T>(r.RN));
    line::reg::Json extra = line::reg::Json::object();
    const line::solvers::CacheMetrics<T> cache = line::ssa::cache_metrics_of_ssa<T>(sn, cacheratio);
    if (!cache.empty()) extra["Cache"] = cache_extra_json<T>(cache);
    // ArvR IS NOT Tput, and reading it off the throughput column was wrong
    // wherever the two differ -- most visibly at a JOIN, which takes in one
    // sibling per branch and fires once per parent, so its arrival rate is the
    // fork degree times its throughput. `ssa_fj_foldback` already divides the
    // Join's QLen by the DERIVED rate to get its response time, so reporting
    // Tput in the ArvR column left the printed row self-inconsistent
    // (0.62372/1.02229 is 0.610, not the 0.3038 beside it). Taken on the
    // measured-cache-split struct for the same reason ResidT is.
    const line::Matrix<T> AN = line::mva::sn_get_arvr_from_tput<T>(snw, to_matrix<T>(r.TN));
    emit_avg_table<T>(sn, r.method, [&](std::size_t i, std::size_t c) {
        AvgRow row;
        row.q = r.QN(i, c);
        row.u = r.UN(i, c);
        row.r = r.RN(i, c);
        row.w = line::num_traits<T>::to_double(WN(i, c));
        row.t = r.TN(i, c);
        // A Source has no arrivals TO ITSELF, so its ArvR is 0 while its Tput is
        // the arrival rate.
        row.a = sn.stations[i].sched == line::lang::SchedStrategy::EXT
                    ? 0.0
                    : line::num_traits<T>::to_double(AN(i, c));
        return row;
    }, extra);
    return 0;
}

/** The knobs the fluid solver reads, in one place so every fluid arm reads the
 *  same set: an arm that quietly dropped one would answer a different model. */
inline line::fluid::FluidOptions fluid_options(const Knobs& k) {
    line::fluid::FluidOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = static_cast<std::size_t>(k.iter_max);
    if (k.t1 >= 0.0) opt.timespan_end = k.t1;
    // `pstar_set` is what makes the smoothing active under a method that did
    // not ask for it by name, exactly as `options.config.pstar` does in the
    // reference: setting the exponent alone would leave `-a avg` integrating
    // the hard-min drift while reporting the caller's choice.
    if (k.pstar > 0.0) {
        opt.pstar = k.pstar;
        opt.pstar_set = true;
    }
    return opt;
}

/**
 * Solve a Network model.json with the fluid solver and print the same table as
 * the MVA path, so the parity harness can diff the two rows unchanged.
 *
 * ResidT and ArvR are reported as the per-visit response time and the
 * throughput: the fluid analyzer works at station level and, unlike the MVA
 * runner, has no chain-visit conversion behind it. That is what MATLAB's
 * fluid `getAvgTable` shows for these columns on a single-visit model.
 *
 * Only `double` reaches here -- the drift is integrated by LSODA -- so the
 * other backends are refused by name in `solve_model_dispatch` rather than
 * being narrowed silently.
 */
template <class T>
int solve_model_fluid(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    // `solver_fluid_run_analyzer` is runAnalyzer's resolution over the analyzer, so a
    // Cache model reaches the rmf branch, a DPS model the closing drift, and
    // anything the moment closure accepts reaches `minnormal`.
    // The converged hit/miss split is COLLECTED, for the reason the SSA arm
    // above collects its own: `-a avg` is what MATLAB's lang='cpp' bridge calls,
    // and a missing block there CLEARS the host's Cache node.
    line::solvers::CacheMetrics<T> cache;
    // THE REFRESHED STRUCT IS TAKEN, not dropped: on a cache model the routing
    // the analyzer converged to carries the ACTUAL hit/miss split, where
    // `net.get_struct()` still carries link()'s offered one. The arrival rates
    // below are read off that routing, so the offered split reported 0.5/0.5
    // where the model converged to 0.4/0.6 (cache_replc_routing).
    line::qn::NetworkStruct<T> refreshed;
    const line::fluid::FluidSolution r = line::fluid::solver_fluid_run_analyzer(
        net.get_struct(), opt, static_cast<line::qn::NetworkStruct<T>*>(nullptr),
        &refreshed, &cache);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::qn::NetworkStruct<T>& snflow =
        refreshed.nstations == sn.nstations && refreshed.nclasses == sn.nclasses ? refreshed : sn;

    std::printf("SolverFluid arith=%s method=%s type=%s iters=%zu\n", line::num_traits<T>::name(),
                r.method.c_str(),
                line::util::method_type("FLD", r.method).c_str(), r.iters);
    // The residence-time conversion the SSA arm above applies, for the same
    // reason: neither analyzer produces a per-job residence time, and
    // `sn_get_residt_from_respt` derives one from RN and the visit ratios.
    const line::Matrix<T> WN = line::mva::sn_get_residt_from_respt<T>(sn, to_matrix<T>(r.RN));
    // THE ARRIVAL RATE IS A FLOW, NOT A COPY OF THE THROUGHPUT. `runAnalyzer`
    // takes it from `sn_get_arvr_from_tput`, i.e. from the class-expanded
    // routing, and the two agree only where every job a station serves it also
    // completes. They part on a station a job LEAVES by another route: on
    // cache_replc_routing the fluid solution puts zero throughput on the two
    // Delay stations while 0.4 and 0.6 arrive at them, so copying the
    // throughput made both rows all-zero and the table dropped them.
    const line::Matrix<T> AN = line::mva::sn_get_arvr_from_tput<T>(snflow, to_matrix<T>(r.TN));
    line::reg::Json extra = line::reg::Json::object();
    if (!cache.empty()) extra["Cache"] = cache_extra_json<T>(cache);
    emit_avg_table<T>(sn, r.method, [&](std::size_t i, std::size_t c) {
        AvgRow row;
        row.q = r.QN(i, c);
        row.u = r.UN(i, c);
        row.r = r.RN(i, c);
        row.w = line::num_traits<T>::to_double(WN(i, c));
        row.t = r.TN(i, c);
        // A Source has no arrivals TO ITSELF, so its ArvR is 0 while its Tput is
        // the arrival rate -- the same rule the SSA arm above applies. Without it
        // the two CLI paths disagreed on one column of the same model:
        // gallery_mm1 reported Source ArvR 0 under -s mva and 1 under -s fluid.
        row.a = sn.stations[i].sched == line::lang::SchedStrategy::EXT
                    ? 0.0
                    : line::num_traits<T>::to_double(AN(i, c));
        return row;
    }, extra);
    return 0;
}

/**
 * `-s fluid -a statevec`: the converged FLUID STATE VECTOR, `result.odeStateVec`.
 *
 * NOT A METRIC AND NOT INDEXED LIKE ONE. The ODE state carries one coordinate
 * per (station, class, PHASE), so a two-phase Erlang service contributes two
 * entries where the AvgTable contributes one number, and the sum over a
 * station's phases is its mean queue length. It is what a caller needs to
 * restart an integration, to seed another solver, or to read the phase
 * occupancy the means average away -- which is why the JAR exposes it as its
 * own `-a statevec` rather than as a column.
 *
 * The index layout is `fluid_state_layout`'s and is emitted BESIDE the vector,
 * because a bare list of numbers cannot be related back to a station without it.
 */
template <class T>
int solve_model_fluid_statevec(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::fluid::FluidSolution r = line::fluid::solver_fluid_run_analyzer(net.get_struct(), opt);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    if (r.xvec.empty())
        throw line::UnsupportedError(
            "-a statevec reports the converged ODE state and this solve produced none; the rmf "
            "and closing branches integrate a drift and fill it, so a branch that returns means "
            "directly has no state vector to report");

    // The (station, class, phase) each coordinate belongs to, in the order the
    // ODE state is laid out: station-major, then class, then phase.
    std::vector<std::size_t> ist, cls, phs;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const std::size_t np = sn.phases_of(i + 1, c + 1);
            for (std::size_t j = 0; j < np; ++j) {
                ist.push_back(i);
                cls.push_back(c);
                phs.push_back(j);
            }
        }
    // The layout is a CLAIM about the solver's state ordering, so it is checked
    // rather than asserted in a comment: a mismatch means the labels below would
    // name the wrong station, which is worse than no labels at all.
    const bool labelled = ist.size() == r.xvec.size();

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "FluidStateVec";
        p["indexBase"] = 0;
        p["xvec"] = vector_json(r.xvec);
        p["labelled"] = labelled;
        if (labelled) {
            line::reg::Json st = line::reg::Json::array(), cl = line::reg::Json::array(),
                            ph = line::reg::Json::array();
            for (std::size_t j = 0; j < ist.size(); ++j) {
                st.push_back(sn.stations[ist[j]].name);
                cl.push_back(sn.classes[cls[j]].name);
                ph.push_back(phs[j]);
            }
            p["Station"] = st;
            p["JobClass"] = cl;
            p["Phase"] = ph;
        }
        emit_analysis<T>("statevec", p, r.method);
        return 0;
    }
    std::printf("SolverFluid arith=%s method=%s coords=%zu\n", line::num_traits<T>::name(),
                r.method.c_str(), r.xvec.size());
    if (!labelled) {
        std::printf("# the ODE state is %zu wide and the (station, class, phase) layout accounts "
                    "for %zu; the coordinates are printed unlabelled\n",
                    r.xvec.size(), ist.size());
        std::printf("%-10s %20s\n", "Index", "x");
        for (std::size_t j = 0; j < r.xvec.size(); ++j)
            std::printf("%-10zu %20.10g\n", j, r.xvec[j]);
        return 0;
    }
    std::printf("%-16s %-14s %-8s %20s\n", "Station", "JobClass", "Phase", "x");
    for (std::size_t j = 0; j < r.xvec.size(); ++j)
        std::printf("%-16s %-14s %-8zu %20.10g\n", sn.stations[ist[j]].name.c_str(),
                    sn.classes[cls[j]].name.c_str(), phs[j], r.xvec[j]);
    return 0;
}

/**
 * Solve an ENVIRONMENT model.json with SolverENV and print the same average
 * table every other arm prints.
 *
 * THE COLUMNS ARE THE REFERENCE'S, INCLUDING THE TWO IT LEAVES EMPTY.
 * `@@SolverENV/getEnsembleAvg` returns Q, U and T from the coupling, sets
 * `WNclass = QNclass ./ TNclass` and returns `RNclass` and `ANclass` as NaN --
 * ENV blends per-stage metrics over the environment process and computes no
 * response time or arrival rate at all. Printing Q/T under RespT here would
 * invent a number the reference declines to give, so RespT and ArvR are NaN and
 * ResidT carries the Little's-law ratio, exactly as MATLAB's table does.
 *
 * The station and class names come from stage 1. The mean-field coupling
 * already refuses an environment whose stages disagree on the station or class
 * count, so any stage names the same rows.
 */
/**
 * Run the coupling `o.method` names, at the arithmetic the caller asked for.
 *
 * WHY THIS IS NOT JUST `env::solver_env`. That entry instantiates BOTH
 * couplings, and the mean-field one solves each stage with the fluid transient
 * -- LSODA, hence double. Calling it at `Rational` does not merely give a worse
 * answer, it does not compile (`sqrt` on a rational), so the template below
 * carries only the state-vector coupling and the double overload beside it
 * carries the full dispatch. The dispatcher has already refused `-s env
 * --arith exact` without `--method statevec`, so a non-double run reaching here HAS
 * asked for the state-vector coupling and gets it, banner included.
 */
template <class T>
line::env::EnvAnalyzerSolution<T> env_run(line::env::Environment<T>& e,
                                          const line::env::EnvOptions& o) {
    line::env::EnvAnalyzerSolution<T> out;
    out.statevec =
        line::env::solver_env_statevec(e, line::env::dispatch_detail::env_statevec_options<T>(o));
    line::env::dispatch_detail::env_take_statevec(out);
    return out;
}

/** The double case, where both couplings are available. */
inline line::env::EnvAnalyzerSolution<double> env_run(line::env::Environment<double>& e,
                                                      const line::env::EnvOptions& o) {
    return line::env::solver_env(e, o);
}

template <class T>
int solve_model_env(const std::string& file, const Knobs& k) {
    line::env::Environment<T> e = read_env_model<T>(file);
    line::env::EnvOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (k.t1 >= 0.0) opt.timespan_end = k.t1;
    if (k.tran_points) opt.tran_points = k.tran_points;
    if (k.tol >= 0.0) opt.stage.tol = k.tol;
    // THE STAGE SOLVER FOLLOWS THE COUPLING, because each coupling has exactly
    // one. `EnvOptions::stage_solver` defaults to `fluid`, which is what the
    // mean-field coupling needs (a transient mean) and what the state-vector one
    // refuses (it needs an enumerated generator and state space, which only
    // SolverCTMC exposes). Leaving the default in place made every
    // `--method statevec` run die on "stage solver 'fluid' is not available",
    // with no flag to fix it. This is not a silent fallback: there is one
    // admissible stage solver per coupling, and picking the other would be the
    // error.
    if (opt.method == "statevec" || opt.method == "blend") opt.stage_solver = "ctmc";
    // `--stage-solver` OVERRIDES that default, and only the mean-field coupling
    // has a choice to make: it needs a transient mean, which both the fluid
    // analyzer and the enumerated CTMC produce, and the two are different
    // models rather than two routes to one answer (a chain holds whole jobs).
    // An ensemble built on SolverCTMC stages therefore has to say so, or the
    // engine answers the fluid ensemble under its name -- which is what the
    // hosts refused lang='cpp' for.
    if (!k.stage_solver.empty()) opt.stage_solver = k.stage_solver;
    if (k.cutoff >= 0.0) opt.stage_cutoff = k.cutoff;

    // UNQUALIFIED, so the double overload above wins for T = double: naming the
    // template explicitly would send every arithmetic to the state-vector
    // coupling and quietly ignore `--method meanfield`.
    const line::env::EnvAnalyzerSolution<T> r = env_run(e, opt);
    const line::qn::NetworkStruct<T>& sn = e.stage(0).model;

    // The horizon and the grid are part of the ANSWER on this path, for the
    // same reason the seed and the sample count are on the SSA one: the
    // mean-field exit metrics are a quadrature, and two runs are the same
    // measurement only if both knobs are stated.
    // `points` is the mean-field quadrature's and is printed only there: the
    // state-vector coupling carries the whole joint law across a switch and
    // never sums over that grid, so reporting a grid it did not use would
    // describe a computation that did not happen.
    // A closed-form limit reads neither knob: it solves each stage, or one
    // rate-averaged model, in STEADY STATE, so a horizon and an iteration count
    // would describe a transient and a fixed point that never ran.
    if (r.method == "avg" || r.method == "dec")
        std::printf("SolverENV arith=%s method=%s stages=%zu (closed-form limit)\n",
                    line::num_traits<T>::name(), r.method.c_str(), e.nstages());
    else if (r.method == "statevec")
        std::printf("SolverENV arith=%s method=%s stages=%zu horizon=%.6g iters=%d%s\n",
                    line::num_traits<T>::name(), r.method.c_str(), e.nstages(), opt.timespan_end,
                    r.iterations, r.converged ? "" : " (NOT CONVERGED)");
    else
        // Every remaining method -- meanfield, and the smp and statedep runs of
        // the same analyzer -- sums the same quadrature, so all of them report
        // the grid it was summed over.
        std::printf("SolverENV arith=%s method=%s stages=%zu horizon=%.6g points=%zu iters=%d%s\n",
                    line::num_traits<T>::name(), r.method.c_str(), e.nstages(), opt.timespan_end,
                    opt.tran_points, r.iterations, r.converged ? "" : " (NOT CONVERGED)");
    emit_avg_table<T>(sn, r.method, [&](std::size_t i, std::size_t c) {
        AvgRow row;
        row.q = line::num_traits<T>::to_double(r.QN(i, c));
        row.u = line::num_traits<T>::to_double(r.UN(i, c));
        row.t = line::num_traits<T>::to_double(r.TN(i, c));
        row.r = std::numeric_limits<double>::quiet_NaN();
        row.a = std::numeric_limits<double>::quiet_NaN();
        row.w = row.q / row.t;
        return row;
    });
    return 0;
}

/**
 * `-a var`: the SECOND moment of the queue length, which only the fluid solver
 * has and only through two of its methods.
 *
 * `minnormal` and `refined` report the STATIONARY covariance of the linear noise
 * approximation (`@@SolverFLD/getMoments`), `kp` the covariance integrated along
 * the trajectory (`@@SolverFLD/getTranAvgVar`); this prints the per-station,
 * per-class variance and its standard deviation, plus, on the JSON path, the full
 * state covariance so that cross-station terms survive rather than only the
 * per-block totals. Every other method carries a first moment only and is refused
 * by name -- a variance of zero would be a claim, not an absence.
 */
template <class T>
int solve_model_fluid_var(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::fluid::FluidSolution r = line::fluid::solver_fluid_run_analyzer(sn, opt);
    if (!r.has_moments)
        throw line::UnsupportedError(
            "-a var needs a fluid method that computes a second moment: 'minnormal', 'refined' or "
            "'dae' for the stationary covariance, 'kp' for the covariance along the trajectory. "
            "The '" +
            r.method + "' method integrates the mean only");

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "QueueLengthVariance";
        p["indexBase"] = 0;
        p["Station"] = line::reg::Json::array();
        p["JobClass"] = line::reg::Json::array();
        p["QVar"] = line::reg::Json::array();
        p["QStd"] = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                if (r.moments.QVar(i, c) == 0.0) continue;
                p["Station"].push_back(sn.stations[i].name);
                p["JobClass"].push_back(sn.classes[c].name);
                p["QVar"].push_back(r.moments.QVar(i, c));
                p["QStd"].push_back(r.moments.QStd(i, c));
            }
        p["Sigma"] = matrix_json<double>(r.moments.Sigma);
        emit_analysis<T>("var", p, r.method);
        return 0;
    }
    std::printf("SolverFluid arith=%s method=%s second moment\n", line::num_traits<T>::name(),
                r.method.c_str());
    std::printf("%-16s %-14s %12s %12s\n", "Station", "JobClass", "QVar", "QStd");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (r.moments.QVar(i, c) == 0.0) continue;
            std::printf("%-16s %-14s %12.6g %12.6g\n", sn.stations[i].name.c_str(),
                        sn.classes[c].name.c_str(), r.moments.QVar(i, c), r.moments.QStd(i, c));
        }
    return 0;
}

/**
 * `-a odes`: `@@SolverFLD/exportODEs`, the drift itself rather than its fixed
 * point.
 *
 * The output is the LaTeX document the reference writes, printed to stdout so
 * that it can be redirected. It carries a machine-readable comment header
 * naming every state variable and every event, which is what makes the document
 * diffable against MATLAB's rather than only readable.
 */
template <class T>
int solve_model_fluid_odes(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::fluid::FluidOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    // `--notation` REACHES THE EXPORTER, and is not fixed at "scalar" here: the
    // matrix form is a different document of the same drift, and hardcoding one
    // while accepting a flag naming the other is the silent-acceptance defect
    // this CLI refuses everywhere else. An unrecognised name is refused by
    // `export_odes_latex` itself rather than defaulted.
    const std::string notation = k.notation.empty() ? "scalar" : k.notation;
    const std::string tex = line::fluid::solver_fluid_export_odes(sn, opt, notation, sn.name);
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ODEs";
        p["notation"] = notation;
        // The document goes in a STRING VALUE, escaped by the JSON writer: it is
        // LaTeX, so it carries backslashes and newlines that a host reading raw
        // stdout would have to re-parse out of the surrounding table.
        p["latex"] = tex;
        emit_analysis<T>("odes", p,
                         line::fluid::detail::fluid_resolve_method(sn, opt.method, opt));
        return 0;
    }
    std::printf("%s\n", tex.c_str());
    return 0;
}

/**
 * `-a jacobian`: `@@SolverFLD/getJacobian`, d f_i / d x_j of the mean-field
 * drift, with the equilibria beside it when they are asked for.
 *
 * This is the fixed point's LOCAL BEHAVIOUR, which no integration reports: the
 * eigenvalues of J tell a stable fixed point from a limit cycle and give the
 * rate at which the fluid approximation converges to it.
 *
 * The method must be a smooth one. `fluid_symbolic_drift` refuses the min-scaled
 * drifts by the factor that carries the kink, before any backend is contacted.
 */
template <class T>
int solve_model_fluid_jacobian(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::fluid::FluidOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    // pstar REACHES THE SYSTEM ONLY UNDER `pnorm`, which is the rule the
    // integrator and the exporter both follow: `matrix` and `default` leave it
    // at zero, selecting the hard min. That is why they have no Jacobian here
    // and `pnorm` does, and it is the reference's rule too -- MATLAB reads
    // options.config.pstar, which is unset unless asked for.
    std::string m = opt.method;
    if (m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    const line::fluid::FluidSymSystem sys = line::fluid::fluid_symodes(
        sn, opt.method, (m == "pnorm" || opt.pstar_set) ? opt.pstar : 0.0, std::vector<double>());
    line::fluid::FluidSymbolicOptions symopt;
    if (!k.symbolic.empty()) symopt.backend = k.symbolic;
    symopt.equilibria = k.equilibria;
    const line::fluid::FluidJacobian jac = line::fluid::fluid_jacobian(sys, symopt);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "Jacobian";
        p["engine"] = jac.engine;
        p["vars"] = line::reg::Json(jac.vars);
        p["rhs"] = line::reg::Json(jac.rhs);
        line::reg::Json rows = line::reg::Json::array();
        for (std::size_t i = 0; i < jac.J.size(); ++i) rows.push_back(line::reg::Json(jac.J[i]));
        p["jacobian"] = rows;
        // `hasEquilibria` separates "asked and answered with none" from "never
        // asked"; an empty list alone would read as "this system has none".
        p["hasEquilibria"] = jac.has_equilibria;
        line::reg::Json eqs = line::reg::Json::array();
        for (std::size_t e = 0; e < jac.equilibria.size(); ++e) {
            line::reg::Json one = line::reg::Json::object();
            for (std::map<std::string, std::string>::const_iterator it = jac.equilibria[e].begin();
                 it != jac.equilibria[e].end(); ++it)
                one[it->first] = it->second;
            eqs.push_back(one);
        }
        p["equilibria"] = eqs;
        emit_analysis<T>("jacobian", p,
                         line::fluid::detail::fluid_resolve_method(sn, opt.method, opt));
        return 0;
    }

    std::printf("engine=%s states=%zu\n", jac.engine.c_str(), jac.vars.size());
    for (std::size_t i = 0; i < jac.rhs.size(); ++i)
        std::printf("d%s/dt = %s\n", jac.vars[i].c_str(), jac.rhs[i].c_str());
    for (std::size_t i = 0; i < jac.J.size(); ++i)
        for (std::size_t j = 0; j < jac.J[i].size(); ++j) {
            // A structurally zero entry is printed, not skipped: a reader must be
            // able to tell a zero derivative from a row this port never emitted.
            std::printf("J[%s,%s] = %s\n", jac.vars[i].c_str(), jac.vars[j].c_str(),
                        jac.J[i][j].c_str());
        }
    if (jac.has_equilibria) {
        if (jac.equilibria.empty())
            std::printf("equilibria: none in closed form (the solve found none, which is not a "
                        "proof that none exist)\n");
        for (std::size_t e = 0; e < jac.equilibria.size(); ++e)
            for (std::map<std::string, std::string>::const_iterator it = jac.equilibria[e].begin();
                 it != jac.equilibria[e].end(); ++it)
                std::printf("equilibrium %zu: %s = %s\n", e + 1, it->first.c_str(),
                            it->second.c_str());
    }
    return 0;
}

/**
 * `-s fluid -a tranvar`: `@@SolverFLD/getTranAvgVar`, the queue-length VARIANCE
 * along the trajectory, plus the full state covariance at each time point.
 *
 * NOT `-a var`, which reports the STATIONARY covariance of `minnormal` /
 * `refined` -- one number per (station, class) at the fixed point. This is the
 * diffusion limit of Ko and Pender integrated alongside the fluid limit, so it
 * has a value at every t, and only `--method kp` produces it. Asking any other
 * method for it is an error rather than a misleading zero, which is the header's
 * own rule and is left to the header to enforce so the two flags cannot drift.
 */
template <class T>
int solve_model_fluid_tranvar(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::fluid::FluidKpTransient tr = line::fluid::solver_fluid_tran_avg_var(sn, opt);
    if (tr.t.empty())
        throw line::UnsupportedError("-a tranvar produced no trajectory points");

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranAvgVarTable";
        p["indexBase"] = 0;
        p["t0"] = tr.t.front();
        p["t1"] = tr.t.back();
        line::reg::Json ts = line::reg::Json::array();
        for (std::size_t j = 0; j < tr.t.size(); ++j) ts.push_back(tr.t[j]);
        p["t"] = ts;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                if (sn.disabled[i][c]) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                line::reg::Json v = line::reg::Json::array();
                for (std::size_t j = 0; j < tr.QVar.size(); ++j) v.push_back(tr.QVar[j](i, c));
                e["QVar"] = v;
                arr.push_back(e);
            }
        p["curves"] = arr;
        // The full covariance is the answer's other half: the per-pair variances
        // are its diagonal, and a caller asking for the diffusion limit wants the
        // off-diagonal correlations the limit is about.
        line::reg::Json sig = line::reg::Json::array();
        for (std::size_t j = 0; j < tr.Sigma.size(); ++j)
            sig.push_back(matrix_json<double>(tr.Sigma[j]));
        p["Sigma"] = sig;
        emit_analysis<T>("tranvar", p, "kp");
        return 0;
    }
    std::printf("SolverFluid arith=%s method=kp tspan=[%g,%g] points=%zu dim=%zu\n",
                line::num_traits<T>::name(), tr.t.front(), tr.t.back(), tr.t.size(),
                tr.Sigma.empty() ? std::size_t(0) : tr.Sigma.front().rows());
    std::printf("%-16s %-14s %12s %12s %12s\n", "Station", "JobClass", "Time", "QVar", "QStd");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (sn.disabled[i][c]) continue;
            for (std::size_t j = 0; j < tr.QVar.size(); ++j) {
                const double var = tr.QVar[j](i, c);
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g\n", sn.stations[i].name.c_str(),
                            sn.classes[c].name.c_str(), tr.t[j], var,
                            var >= 0.0 ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN());
            }
        }
    return 0;
}

/**
 * `-s fluid -a tran`: `@@SolverFLD/getTranAvg`, the metrics ALONG the
 * trajectory rather than at its fixed point.
 *
 * `--tspan` is optional here and required on the MAM arm, and the difference is
 * the reference's: `options.timespan` defaults to [0, Inf] for the fluid solver,
 * which does not mean "integrate forever" but "integrate until the state stops
 * moving" -- `solver_fluid_tran_avg` reproduces the horizon that adaptive loop
 * converges at. A caller that names a horizon gets exactly that one.
 *
 * The reference forces `closing` for a transient (the matrix and smoothed
 * variants are steady-state devices) and warns when it does; the port forces it
 * too, and the banner names the method that actually integrated. `dae` is the
 * one exception the reference itself makes -- it has a trajectory of its own,
 * with conservation carried as an algebraic equation -- so
 * `solver_fluid_run_transient` routes it rather than substituting the
 * first-order drift under its name.
 */
template <class T>
int solve_model_fluid_tran(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const std::vector<line::fluid::FluidTranPoint> tr =
        line::fluid::solver_fluid_run_transient(sn, opt);
    if (tr.empty()) throw line::UnsupportedError("-a tran produced no trajectory points");

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranAvgTable";
        p["indexBase"] = 0;
        p["t0"] = 0.0;
        p["t1"] = tr.back().t;
        line::reg::Json ts = line::reg::Json::array();
        for (std::size_t j = 0; j < tr.size(); ++j) ts.push_back(tr[j].t);
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                if (sn.disabled[i][c]) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = ts;
                line::reg::Json q = line::reg::Json::array(), u = line::reg::Json::array(),
                                 x = line::reg::Json::array();
                for (std::size_t j = 0; j < tr.size(); ++j) {
                    q.push_back(tr[j].QN(i, c));
                    u.push_back(tr[j].UN(i, c));
                    x.push_back(tr[j].TN(i, c));
                }
                e["QLen"] = q;
                e["Util"] = u;
                e["Tput"] = x;
                arr.push_back(e);
            }
        p["curves"] = arr;
        emit_analysis<T>("tran", p, "closing");
        return 0;
    }
    std::printf("SolverFluid arith=%s method=closing tspan=[0,%g] points=%zu\n",
                line::num_traits<T>::name(), tr.back().t, tr.size());
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "Time", "QLen", "Util",
                "Tput");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (sn.disabled[i][c]) continue;
            for (std::size_t j = 0; j < tr.size(); ++j)
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n",
                            sn.stations[i].name.c_str(), sn.classes[c].name.c_str(), tr[j].t,
                            tr[j].QN(i, c), tr[j].UN(i, c), tr[j].TN(i, c));
        }
    return 0;
}

/**
 * `-s fluid -a prob`: `@@SolverFLD/getProbAggr`, the probability that a station
 * holds the marginal population of the model's default state.
 *
 * IT IS NOT THE MVA ARM'S ANSWER AND IS NOT MEANT TO BE. The fluid solver has
 * no state space, so the law is fitted to the means it does produce -- a
 * binomial per closed class, the BCMP marginal per open one -- and under a
 * moment closure it is instead the multivariate normal the closure supplies,
 * correlation between the classes included. Two solvers disagreeing here is the
 * approximation showing, not a defect.
 *
 * The log-probability is reported beside it because the fitted law underflows
 * on a large population, where the linear value is 0 and the log one is not.
 */
template <class T>
int solve_model_fluid_prob(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::fluid::FluidSolution r = line::fluid::solver_fluid_run_analyzer(sn, opt);

    std::vector<double> pr(sn.nstations, 0.0), lg(sn.nstations, 0.0);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        pr[i] = line::fluid::fluid_prob_aggr(sn, r, i + 1, &lg[i]);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbAggr";
        p["indexBase"] = 0;
        line::reg::Json st = line::reg::Json::array(), pa = line::reg::Json::array(),
                         lp = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            st.push_back(sn.stations[i].name);
            pa.push_back(pr[i]);
            lp.push_back(lg[i]);
        }
        p["Station"] = st;
        p["ProbAggr"] = pa;
        p["logProbAggr"] = lp;
        emit_analysis<T>("prob", p, r.method);
        return 0;
    }
    std::printf("SolverFluid arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                r.method.c_str(), line::util::method_type("FLD", r.method).c_str());
    std::printf("%-16s %14s %14s\n", "Station", "ProbAggr", "logProbAggr");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        std::printf("%-16s %14.10g %14.10g\n", sn.stations[i].name.c_str(), pr[i], lg[i]);
    return 0;
}

/**
 * `-s fluid -a cdf`: `@@SolverFLD/getCdfRespT`, the WHOLE response-time law per
 * (station, class) and not only its mean.
 *
 * The law is read off a second integration in which the jobs present at the
 * steady state are MARKED and followed to their departure, so it is the
 * stationary response-time distribution of the fluid model. The solve that
 * produces the state to mark in is run here, as the reference runs it: its
 * `getCdfRespT` clears the cached result and re-runs `getAvg` first.
 */
template <class T>
int solve_model_fluid_cdf(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::fluid::FluidOptions opt = fluid_options(k);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const std::vector<std::vector<line::fluid::FluidPassage> > RD =
        line::fluid::solver_fluid_cdf_respt(sn, opt);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "CdfRespT";
        p["indexBase"] = 0;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t i = 0; i < RD.size(); ++i)
            for (std::size_t c = 0; c < RD[i].size(); ++c) {
                // An empty curve is an ABSENT law -- a Source, or a class the
                // station does not serve -- and is omitted rather than sent as a
                // degenerate one.
                if (RD[i][c].t.empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = sn.stations[i].name;
                e["JobClass"] = sn.classes[c].name;
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = line::reg::Json(RD[i][c].t);
                e["F"] = line::reg::Json(RD[i][c].cdf);
                arr.push_back(e);
            }
        p["respt"] = arr;
        emit_analysis<T>("cdf", p, std::string());
        return 0;
    }
    std::printf("SolverFluid arith=%s\n", line::num_traits<T>::name());
    std::printf("%-16s %-14s %14s %14s\n", "Station", "JobClass", "Time", "F(t)");
    for (std::size_t i = 0; i < RD.size(); ++i)
        for (std::size_t c = 0; c < RD[i].size(); ++c)
            for (std::size_t j = 0; j < RD[i][c].t.size(); ++j)
                std::printf("%-16s %-14s %14.8g %14.10g\n", sn.stations[i].name.c_str(),
                            sn.classes[c].name.c_str(), RD[i][c].t[j], RD[i][c].cdf[j]);
    return 0;
}

/**
 * `-s fluid -a aoi`: `@@SolverFLD/getAvgAoI` and `getCdfAoI` in one answer, the
 * Age of Information and Peak AoI laws of a status-update system.
 *
 * ONLY THE `mfq` METHOD HAS THEM, and only on the topology the age laws are
 * defined for: one open class through Source -> Queue -> Sink, a single server,
 * capacity 1 (bufferless) or 2 (single buffer), FCFS/LCFS/LCFSPR. The topology
 * is tested first so a model that is not one is told WHICH condition it fails
 * rather than being handed a number computed for a different system.
 *
 * `--method` may only say `mfq` here: silently overriding a caller who asked for
 * another method would report the age laws under a method that does not produce
 * them.
 */
template <class T>
int solve_model_fluid_aoi(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    const line::fluid::AoiTopology top = line::fluid::aoi_is_aoi(sn);
    if (!top.ok)
        throw line::UnsupportedError(
            "-a aoi reports the age of a status-update system and needs the topology the age laws "
            "are defined for: " +
            (top.error.empty() ? std::string("this model is not one") : top.error));
    line::fluid::FluidOptions opt = fluid_options(k);
    const std::string requested = line::fluid::detail::fluid_unqualify(opt.method);
    if (requested != "default" && requested != "mfq")
        throw line::UnsupportedError(
            "-a aoi is the AoI branch of the 'mfq' method; '" + requested +
            "' integrates the mean-field drift and carries no age process");
    opt.method = "mfq";
    const line::fluid::FluidSolution r = line::fluid::solver_fluid_run_analyzer(sn, opt);
    if (!r.has_aoi)
        throw line::UnsupportedError(
            "-a aoi: the 'mfq' method did not take its AoI branch on this model");

    // The grid `getCdfAoI` builds when the caller names no time points: five
    // mean ages, which covers the bulk of both laws.
    const double base = (std::isfinite(r.aoi.aoi.mean) && r.aoi.aoi.mean > 0.0) ? r.aoi.aoi.mean : 1.0;
    const std::size_t np = 200;
    std::vector<double> tv(np), fa(np), fp(np);
    for (std::size_t j = 0; j < np; ++j) {
        tv[j] = 5.0 * base * static_cast<double>(j) / static_cast<double>(np - 1);
        fa[j] = line::fluid::aoi_cdf(r.aoi.aoi, tv[j]);
        fp[j] = line::fluid::aoi_cdf(r.aoi.paoi, tv[j]);
    }
    const double asd = std::sqrt(std::max(0.0, r.aoi.aoi.var));
    const double psd = std::sqrt(std::max(0.0, r.aoi.paoi.var));

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AoI";
        p["systemType"] = r.aoi.system_type;
        p["preemption"] = r.aoi.preemption;
        p["AoIMean"] = r.aoi.aoi.mean;
        p["AoIVar"] = r.aoi.aoi.var;
        p["AoIStd"] = asd;
        p["PAoIMean"] = r.aoi.paoi.mean;
        p["PAoIVar"] = r.aoi.paoi.var;
        p["PAoIStd"] = psd;
        p["t"] = line::reg::Json(tv);
        p["AoICdf"] = line::reg::Json(fa);
        p["PAoICdf"] = line::reg::Json(fp);
        // The (g, A, h) DENSITY triples, not only the curve evaluated above.
        // `getCdfAoI` takes an optional t_values, and a caller who names their
        // own grid cannot be served from a fixed 200-point one; with the triple
        // they evaluate the same law at their own abscissae. `solve_mfq_aoi`
        // normalizes g so that g*expm(A t)*h is the density, so the survival
        // function carries an extra inv(A) -- the MATLAB getter's own note.
        p["AoI_g"] = line::reg::Json(r.aoi.aoi.g);
        p["AoI_A"] = matrix_json(r.aoi.aoi.A);
        p["AoI_h"] = line::reg::Json(r.aoi.aoi.h);
        p["PAoI_g"] = line::reg::Json(r.aoi.paoi.g);
        p["PAoI_A"] = matrix_json(r.aoi.paoi.A);
        p["PAoI_h"] = line::reg::Json(r.aoi.paoi.h);
        emit_analysis<T>("aoi", p, r.method);
        return 0;
    }
    std::printf("SolverFluid arith=%s method=mfq system=%s preemption=%.6g\n",
                line::num_traits<T>::name(), r.aoi.system_type.c_str(), r.aoi.preemption);
    std::printf("%-8s %14s %14s %14s\n", "Metric", "Mean", "Var", "Std");
    std::printf("%-8s %14.10g %14.10g %14.10g\n", "AoI", r.aoi.aoi.mean, r.aoi.aoi.var, asd);
    std::printf("%-8s %14.10g %14.10g %14.10g\n", "PAoI", r.aoi.paoi.mean, r.aoi.paoi.var, psd);
    std::printf("%14s %14s %14s\n", "Time", "F_AoI(t)", "F_PAoI(t)");
    for (std::size_t j = 0; j < np; ++j)
        std::printf("%14.8g %14.10g %14.10g\n", tv[j], fa[j], fp[j]);
    return 0;
}

/**
 * Solve a Network model.json and print its aggregate state probabilities:
 * getProbSysAggr (the whole-system joint) and getProbAggr per station, over the
 * model's default initial state. These fit the MVA means and so need logarithms;
 * under exact/Rational they refuse by name.
 */
template <class T>
int solve_model_prob(const std::string& file, const Knobs& k) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the -a prob analysis fits a binomial/product-form law and needs transcendental "
            "arithmetic; rerun with --arith double or --arith real");
    } else {
        line::qn::Network<T> net = read_model<T>(file);
        line::mva::MvaOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        line::Matrix<T> init;
        const line::mva::AvgResult<T> r = line::mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
        const line::qn::NetworkStruct<T>& sn = net.get_struct();
        const line::mva::AggrResult<T> ps =
            line::mva::solver_mva_get_prob_sys_aggr(sn, r, opt.method);
        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "ProbAggr";
            p["indexBase"] = 0;
            p["ProbSysAggr"] = line::num_traits<T>::to_double(ps.P);
            line::reg::Json st = line::reg::Json::array(), pa = line::reg::Json::array();
            for (std::size_t i = 0; i < sn.nstations; ++i) {
                st.push_back(sn.stations[i].name);
                pa.push_back(line::num_traits<T>::to_double(
                    line::mva::solver_mva_get_prob_aggr(sn, r, i + 1, opt.method).P));
            }
            p["Station"] = st;
            p["ProbAggr"] = pa;
            emit_analysis<T>("prob", p, r.actualmethod);
            return 0;
        }
        std::printf("SolverMVA arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                    r.actualmethod.c_str(),
                    line::util::method_type("MVA", r.actualmethod).c_str());
        std::printf("ProbSysAggr %.10g\n", line::num_traits<T>::to_double(ps.P));
        std::printf("%-16s %14s\n", "Station", "ProbAggr");
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const line::mva::AggrResult<T> pa =
                line::mva::solver_mva_get_prob_aggr(sn, r, i + 1, opt.method);
            std::printf("%-16s %14.10g\n", sn.stations[i].name.c_str(),
                        line::num_traits<T>::to_double(pa.P));
        }
        return 0;
    }
}

/**
 * `-s mva -a marg`: `@@SolverMVA/getProbMarg`, P(n jobs of class r at station i).
 *
 * THE WHOLE GRID BY DEFAULT, one curve per (station, class): the reference takes
 * the station and the class as arguments and this CLI has no notion of a
 * "current" pair, so reporting every pair is the only reading that answers the
 * method rather than a choice this file would be making on the caller's behalf.
 * `--node` and `--class` narrow it to one node's station and one class, and
 * `--marg-states` is the reference's third argument `state_m`: the n values to
 * report, in place of the default range each case picks for itself (0..N_r for a
 * closed class, mean + 5 sigma for a Poisson, the 1e-10 tail for a geometric).
 *
 * A NODE THAT IS NOT A STATION IS AN ERROR, not an empty answer: a queue-length
 * law at a ClassSwitch is not a quantity, and defaulting to the whole network
 * after the caller narrowed it would report more than was asked for.
 */
template <class T>
int solve_model_marg(const std::string& file, const Knobs& k) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the -a marg analysis fits a binomial / Poisson / geometric law and needs "
            "transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
        line::qn::Network<T> net = read_model<T>(file);
        line::mva::MvaOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        line::Matrix<T> init;
        const line::mva::AvgResult<T> r = line::mva::solver_mva_run_analyzer(net.get_struct(), opt, init);
        const line::qn::NetworkStruct<T>& sn = net.get_struct();

        std::vector<std::size_t> ists;  // 1-based station indices to report
        if (k.node) {
            if (k.node > sn.nof_nodes())
                throw line::InputError("--node " + std::to_string(k.node) +
                                       " exceeds the number of nodes in the model (" +
                                       std::to_string(sn.nof_nodes()) + ")");
            const std::size_t ist = sn.nodes[k.node - 1].station;
            if (!ist)
                throw line::InputError("--node " + std::to_string(k.node) + " ('" +
                                       sn.nodes[k.node - 1].name +
                                       "') is not a station, and a queue-length distribution is "
                                       "reported per station");
            ists.push_back(ist);
        } else {
            for (std::size_t i = 0; i < sn.nstations; ++i) ists.push_back(i + 1);
        }
        std::vector<std::size_t> rs;  // 1-based class indices to report
        if (k.jobclass) {
            if (k.jobclass > sn.nclasses)
                throw line::InputError("--class " + std::to_string(k.jobclass) +
                                       " exceeds the number of classes in the model");
            rs.push_back(k.jobclass);
        } else {
            for (std::size_t c = 0; c < sn.nclasses; ++c) rs.push_back(c + 1);
        }

        // Every curve first: a pair the reference refuses must not leave a
        // banner and a column header standing above an answer that never came.
        std::vector<line::mva::MargResult<T> > curves;
        for (std::size_t a = 0; a < ists.size(); ++a)
            for (std::size_t b = 0; b < rs.size(); ++b)
                curves.push_back(line::mva::solver_mva_get_prob_marg(sn, r, ists[a], rs[b],
                                                                     k.marg_states, opt.method));

        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "ProbMarg";
            p["indexBase"] = 0;
            line::reg::Json arr = line::reg::Json::array();
            for (std::size_t a = 0, q = 0; a < ists.size(); ++a)
                for (std::size_t b = 0; b < rs.size(); ++b, ++q) {
                    const line::mva::MargResult<T>& m = curves[q];
                    line::reg::Json e = line::reg::Json::object();
                    e["station"] = ists[a] - 1;
                    e["Station"] = sn.stations[ists[a] - 1].name;
                    e["jobclass"] = rs[b] - 1;
                    e["JobClass"] = sn.classes[rs[b] - 1].name;
                    line::reg::Json jobs = line::reg::Json::array();
                    for (std::size_t n = 0; n < m.P.size(); ++n)
                        jobs.push_back(k.marg_states.empty() ? static_cast<long>(n)
                                                             : k.marg_states[n]);
                    e["Jobs"] = jobs;
                    e["P"] = vector_json(m.P);
                    e["logP"] = vector_json(m.logP);
                    arr.push_back(e);
                }
            p["marginal"] = arr;
            emit_analysis<T>("marg", p, r.actualmethod);
            return 0;
        }
        std::printf("SolverMVA arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                    r.actualmethod.c_str(),
                    line::util::method_type("MVA", r.actualmethod).c_str());
        std::printf("%-16s %-14s %-8s %16s\n", "Station", "JobClass", "Jobs", "ProbMarg");
        for (std::size_t a = 0, q = 0; a < ists.size(); ++a)
            for (std::size_t b = 0; b < rs.size(); ++b, ++q) {
                const line::mva::MargResult<T>& m = curves[q];
                for (std::size_t n = 0; n < m.P.size(); ++n)
                    std::printf("%-16s %-14s %-8ld %16.10g\n",
                                sn.stations[ists[a] - 1].name.c_str(),
                                sn.classes[rs[b] - 1].name.c_str(),
                                k.marg_states.empty() ? static_cast<long>(n) : k.marg_states[n],
                                line::num_traits<T>::to_double(m.P[n]));
            }
        return 0;
    }
}

/**
 * `-a normconst`: `@@SolverMVA/getProbNormConstAggr` and `@@SolverNC`'s.
 *
 * ONE ANALYSIS, TWO SOLVERS, AND THEY DO NOT COMPUTE IT THE SAME WAY. The NC
 * arm reads the constant its own solve already formed. The MVA arm RE-ENTERS the
 * analyzer at method='exact', as the reference does, because only the exact MVA
 * recursion carries a G: an AMVA solve has none, and reporting the requested
 * method's number would attribute the constant to an algorithm that never
 * produced one. That re-entry is why this is a separate `-a` and not a field on
 * the `-s mva` banner, where it would charge every average solve for a second
 * exact one.
 *
 * WHAT A MODEL WITH NO CONSTANT REPORTS IS THE ANALYZER'S OWN ANSWER, not a
 * substitution made here, and the two cases differ: the branches that form no G
 * at all -- MVAC, the LCFS chain -- set lG to NaN at the source and print nan,
 * while the open-queue closed forms report lG = 0 exactly as
 * solver_mva_qsys_analyzer.m:54,96,235 does. Neither is edited on the way out.
 */
/**
 * Mean busy period of a named subnetwork, Daduna (J. ACM 35(3), 1988).
 *
 * The transform `solver_nc_busyp` has been in the port since it landed, and
 * `ldes_cli` has answered `--busyperiod` all along, so the ONLY thing between
 * a caller and the analytical form was a `-a` token: asking this CLI for a busy
 * period meant simulating a quantity there is a closed form for.
 *
 * `--busyperiod-subnet` is 1-BASED, as every station index this CLI takes is,
 * and is required: a busy period is defined for a NAMED set of stations and
 * defaulting it would answer about a subnetwork the caller never chose. The
 * orders default to 1, the ordinary busy period.
 */
template <class T>
int solve_model_nc_busyp(const std::string& file, const Knobs& k) {
    if (k.busy_subnet.empty())
        throw line::InputError(
            "-a busyperiod needs --busyperiod-subnet: the busy period is defined for a named "
            "subnetwork of stations, and no default can choose one");
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::vector<std::size_t> subnet;
    for (std::size_t t = 0; t < k.busy_subnet.size(); ++t) {
        if (k.busy_subnet[t] > sn.nstations)
            throw line::InputError("--busyperiod-subnet names station " +
                                   std::to_string(k.busy_subnet[t]) + ", beyond the model's " +
                                   std::to_string(sn.nstations));
        subnet.push_back(k.busy_subnet[t] - 1);
    }
    std::vector<std::size_t> orders = k.busy_orders;
    if (orders.empty()) orders.push_back(1);
    const std::vector<double> b = line::nc::solver_nc_busyp(sn, subnet, orders);

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "BusyPeriod";
        p["indexBase"] = 0;
        line::reg::Json sj = line::reg::Json::array();
        for (std::size_t t = 0; t < subnet.size(); ++t) sj.push_back(subnet[t]);
        line::reg::Json oj = line::reg::Json::array();
        for (std::size_t t = 0; t < orders.size(); ++t) oj.push_back(orders[t]);
        line::reg::Json bj = line::reg::Json::array();
        for (std::size_t t = 0; t < b.size(); ++t) bj.push_back(b[t]);
        p["subnet"] = sj;
        p["orders"] = oj;
        p["b"] = bj;
        emit_analysis<T>("busyperiod", p, "daduna");
        return 0;
    }
    std::printf("SolverNC arith=%s busy period, subnetwork {", line::num_traits<T>::name());
    for (std::size_t t = 0; t < k.busy_subnet.size(); ++t)
        std::printf("%s%zu", t ? "," : "", k.busy_subnet[t]);
    std::printf("}\n");
    for (std::size_t t = 0; t < orders.size(); ++t)
        std::printf("  order %zu   %.10g\n", orders[t], b[t]);
    return 0;
}

template <class T>
int solve_model_normconst(const std::string& file, const Knobs& k, const std::string& solver) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    double lG = 0.0;
    std::string method;
    if (solver == "nc") {
        line::nc::NcSolverOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        const line::mva::AvgResult<T> r = line::nc::solver_nc_run_analyzer(sn, opt);
        lG = r.lognormconst.has_value() ? r.lognormconst.value()
                                        : std::numeric_limits<double>::quiet_NaN();
        method = r.actualmethod;
    } else {
        line::mva::MvaOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        lG = line::num_traits<T>::to_double(line::mva::solver_mva_get_prob_norm_const_aggr(sn, opt));
        method = "exact";
    }
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "NormConst";
        p["indexBase"] = 0;
        // A NaN rides as JSON null, which is the table's nan on the wire.
        p["logNormConstAggr"] = lG;
        emit_analysis<T>("normconst", p, method);
        return 0;
    }
    std::printf("Solver%s arith=%s method=%s lognormconst=%.10g\n",
                solver == "nc" ? "NC" : "MVA", line::num_traits<T>::name(), method.c_str(), lG);
    return 0;
}

/** The model's declared per-class placement, in `solver_nc_*`'s own container. */
template <class T>
line::nc::MarginalState nc_declared_marginal(const line::qn::NetworkStruct<T>& sn) {
    const line::Matrix<T> nir = line::api::sn_declared_marginal<T>(sn);
    line::nc::MarginalState out(sn.nstations, std::vector<int>(sn.nclasses, 0));
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            out[i][r] = static_cast<int>(std::llround(line::num_traits<T>::to_double(nir(i, r))));
    return out;
}

/**
 * The same two probabilities under SolverNC: `@@SolverNC/getProbSysAggr.m` and
 * `@@SolverNC/getProbAggr.m`, over the model's declared state.
 *
 * NOT THE SAME NUMBERS AS `-s mva -a prob`, and that is the point of having
 * both. SolverMVA fits a binomial to its own means (Schmidt 1997); these are a
 * ratio of normalizing constants and are the product-form model's own
 * probabilities exactly. A closed model therefore reports different figures
 * under the two solvers, and the NC ones are the reference.
 */
template <class T>
int solve_model_nc_prob(const std::string& file, const Knobs& k) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the -s nc -a prob analysis exponentiates a difference of log normalizing constants "
            "and needs transcendental arithmetic; rerun with --arith double or --arith real");
    } else {
        line::qn::Network<T> net = read_model<T>(file);
        line::nc::NcSolverOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        const line::qn::NetworkStruct<T>& sn = net.get_struct();

        // The state is the MODEL'S OWN, which `model.json` now carries: a
        // stateful node's declared row is decoded to the per-class counts the
        // reference reads with `State.toMarginal(sn, ist, state{isf})`. Where a
        // station declares none, the default marking is rebuilt for it -- every
        // closed class at its reference station -- which is what `initDefault`
        // would have put there.
        line::nc::MarginalState nir = nc_declared_marginal<T>(sn);
        // `--state` is `getProb(node, state)`'s second argument, decoded the way
        // the reference decodes it: it substitutes the row into `sn.state{isf}`
        // and takes `State.toMarginal` of the result, so what reaches the
        // probability is that node's PER-CLASS COUNTS and every other node's
        // declared ones. Passing the row through untouched would treat an
        // encoding as a job vector, and on a station with phase-type service the
        // two differ in both width and meaning.
        if (!k.state.empty()) {
            if (!k.node)
                throw line::InputError(
                    "--state is the state of ONE node and needs --node to say which");
            const std::size_t ist =
                k.node <= sn.nodes.size() ? sn.nodes[k.node - 1].station : 0;
            if (ist == 0)
                throw line::InputError("--node " + std::to_string(k.node) +
                                       " is not a station, so it has no queue-length state");
            std::vector<std::size_t> ph(sn.nclasses, 1), shift(sn.nclasses, 0);
            std::size_t w = 0;
            for (std::size_t c = 0; c < sn.nclasses; ++c) {
                ph[c] = sn.phases_of(ist, c + 1);
                shift[c] = w;
                w += ph[c];
            }
            std::vector<T> row(k.state.size());
            for (std::size_t i = 0; i < k.state.size(); ++i)
                row[i] = line::num_traits<T>::from_int(k.state[i]);
            const line::qn::Marginal<T> m =
                line::qn::to_marginal(sn, ist, row, ph, shift, sn.nvars_of(k.node));
            for (std::size_t c = 0; c < sn.nclasses; ++c)
                nir[ist - 1][c] =
                    static_cast<int>(std::llround(line::num_traits<T>::to_double(m.nir[c])));
        }

        // THE SOLVE'S lG IS NOT THIS lG, and handing it over here was wrong.
        // `solver_nc_solve` normalizes the SEIDMANN-REDUCED model -- a
        // multiserver station enters as demand/c with the residual folded into
        // the delay -- while the probability identity F_i G_{-i} / G needs the
        // constant of the load-dependent lattice mu(n) = min(n, c) that F_i and
        // G_{-i} are themselves computed on. Mixing the two scaled every
        // probability of a model with a multiserver station by one common
        // factor: on the 2-job Delay -> PS -> PS(c=2) chain the three stations
        // came back 0.17225 / 0.68900 / 0.32536 against the exact 0.18 / 0.72 /
        // 0.34, and the error is invisible on a single-server model because
        // there the two constants coincide. `solver_nc_margaggr` computes its
        // own, once, for every station -- which is the reference's
        // `logNormConstAggr` caching, not a per-station resolve.
        const line::nc::NcMargResult<T> mr = line::nc::solver_nc_margaggr(
            sn, opt, nir, std::numeric_limits<double>::quiet_NaN());
        const T ps = line::nc::solver_nc_getprob_sys_aggr(sn, opt, nir);
        // THE DETAILED PAIR TOO, because `getProb` and `getProbSys` are not the
        // aggregate ones with rounding: `solver_nc_marg` and `solver_nc_joint`
        // carry the class-within-chain split `lg0_i - lG0_i` that the aggregate
        // pair sums out, so on a multichain model they are different numbers
        // rather than the same one to more places.
        const line::nc::NcMargResult<T> mdet = line::nc::solver_nc_marg(
            sn, opt, nir, std::numeric_limits<double>::quiet_NaN());
        const T pjoint = line::nc::solver_nc_joint<T>(sn, opt, nir, nullptr);

        const line::nc::NcSolution<T> d = line::nc::solver_nc_solve(sn, opt);
        const std::string am = (opt.method == "default" && !d.actualmethod.empty() &&
                                d.actualmethod != "default")
                                   ? "default/" + d.actualmethod
                                   : d.actualmethod;
        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "ProbAggr";
            p["indexBase"] = 0;
            p["ProbSysAggr"] = line::num_traits<T>::to_double(ps);
            p["ProbSys"] = line::num_traits<T>::to_double(pjoint);
            line::reg::Json st = line::reg::Json::array(), pa = line::reg::Json::array(),
                            pm = line::reg::Json::array();
            for (std::size_t i = 0; i < sn.nstations; ++i) {
                st.push_back(sn.stations[i].name);
                pa.push_back(line::num_traits<T>::to_double(mr.P[i]));
                pm.push_back(line::num_traits<T>::to_double(mdet.P[i]));
            }
            p["Station"] = st;
            p["ProbAggr"] = pa;
            p["Prob"] = pm;
            emit_analysis<T>("prob", p, am);
            return 0;
        }
        std::printf("SolverNC arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                    am.c_str(), line::util::method_type("NC", am).c_str());
        std::printf("ProbSysAggr %.10g\n", line::num_traits<T>::to_double(ps));
        std::printf("ProbSys     %.10g\n", line::num_traits<T>::to_double(pjoint));
        std::printf("%-16s %14s %14s\n", "Station", "Prob", "ProbAggr");
        for (std::size_t i = 0; i < sn.nstations; ++i)
            std::printf("%-16s %14.10g %14.10g\n", sn.stations[i].name.c_str(),
                        line::num_traits<T>::to_double(mdet.P[i]),
                        line::num_traits<T>::to_double(mr.P[i]));
        return 0;
    }
}

/**
 * `-s nc -a sysmarg`: `@@SolverNC/getProbSysMarg.m`, the JOINT law of the
 * per-station total queue lengths.
 *
 * NEITHER `-a prob` NOR `-a marg`, and the three are worth telling apart.
 * `-a prob` fixes the PER-CLASS population of every station and is a product
 * form; `-a marg` is this law marginalized down to ONE station; this arm is the
 * joint over all of them, with the classes summed out. Each value is the sum of
 * `-a prob` over the whole fibre of per-class tables with these row sums, and
 * that fibre grows combinatorially, so it is evaluated as a permanent of the
 * demand matrix replicated once per job (Ryser 1963) rather than enumerated.
 *
 * The whole lattice of total states is swept, so the printed column sums to one
 * and the sweep pays for the normalizing constant once.
 */
template <class T>
int solve_model_nc_sysmarg(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::nc::NcSolverOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    double Ntot = 0.0;
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        const double pop = sn.classes[r].population;
        if (!std::isfinite(pop))
            throw line::UnsupportedError(
                "getProbSysMarg requires a closed model: the joint law of the total queue lengths "
                "is not defined when a class has an infinite population");
        Ntot += pop;
    }
    const std::vector<std::vector<int> > states = line::pfqn::multichoose_rows(
        static_cast<int>(sn.nstations), static_cast<int>(std::llround(Ntot)));

    std::vector<double> P(states.size(), 0.0);
    for (std::size_t j = 0; j < states.size(); ++j)
        P[j] = line::num_traits<T>::to_double(
            line::nc::solver_nc_getprob_sys_marg(sn, opt, states[j], k.method_perm));

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbSysMarg";
        p["indexBase"] = 0;
        p["engine"] = k.method_perm;
        line::reg::Json st = line::reg::Json::array(), arr = line::reg::Json::array();
        for (std::size_t i = 0; i < sn.nstations; ++i) st.push_back(sn.stations[i].name);
        for (std::size_t j = 0; j < states.size(); ++j) {
            line::reg::Json e = line::reg::Json::object();
            line::reg::Json n = line::reg::Json::array();
            for (std::size_t i = 0; i < sn.nstations; ++i) n.push_back(states[j][i]);
            e["n"] = n;
            e["P"] = P[j];
            arr.push_back(e);
        }
        p["Station"] = st;
        p["states"] = arr;
        emit_analysis<T>("sysmarg", p, opt.method);
        return 0;
    }
    std::printf("SolverNC arith=%s method=%s engine=%s\n", line::num_traits<T>::name(),
                opt.method.c_str(), k.method_perm.c_str());
    for (std::size_t i = 0; i < sn.nstations; ++i)
        std::printf("%10s", sn.stations[i].name.c_str());
    std::printf(" %14s\n", "ProbSysMarg");
    double total = 0.0;
    for (std::size_t j = 0; j < states.size(); ++j) {
        for (std::size_t i = 0; i < sn.nstations; ++i) std::printf("%10d", states[j][i]);
        std::printf(" %14.10g\n", P[j]);
        total += P[j];
    }
    std::printf("%*s %14.10g\n", static_cast<int>(10 * sn.nstations), "sum", total);
    return 0;
}

/**
 * `-s nc -a marg`: `@@SolverNC/getProbMarg.m`, the TOTAL queue-length law.
 *
 * NOT THE SAME QUANTITY AS `-s mva -a marg`, although the reference gives both
 * methods the same name. SolverMVA's getProbMarg is per (station, CLASS) and is
 * a binomial / Poisson / geometric fitted to the solver's own means; SolverNC's
 * is the TOTAL number of jobs at a station, summed over classes, and is exact --
 * a ratio of normalizing constants, obtained either from one `pfqn_procomom`
 * solve (`--method comom`) or by summing the aggregate marginal over the
 * per-class partitions of n. `--class` and `--marg-states` are therefore refused
 * for it rather than ignored: this law has no class argument and its support is
 * 0..sum(N), which the model fixes.
 *
 * Both P and log P are reported. The log is not a formatting of the other: the
 * enumeration forms it first and a probability that underflows to 0 in double
 * still has a finite log, so dropping it would lose the only number left.
 */
template <class T>
int solve_model_nc_marg(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::nc::NcSolverOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    const line::qn::NetworkStruct<T>& sn = net.get_struct();

    std::vector<std::size_t> ists;  // 1-based station indices to report
    if (k.node) {
        if (k.node > sn.nof_nodes())
            throw line::InputError("--node " + std::to_string(k.node) +
                                   " exceeds the number of nodes in the model (" +
                                   std::to_string(sn.nof_nodes()) + ")");
        const std::size_t ist = sn.nodes[k.node - 1].station;
        if (!ist)
            throw line::InputError("--node " + std::to_string(k.node) + " ('" +
                                   sn.nodes[k.node - 1].name +
                                   "') is not a station, and a queue-length distribution is "
                                   "reported per station");
        ists.push_back(ist);
    } else {
        for (std::size_t i = 0; i < sn.nstations; ++i) ists.push_back(i + 1);
    }

    // Every curve first, for solve_model_marg's reason: a station the reference
    // refuses must not leave a header standing above an answer that never came.
    std::vector<line::nc::NcQueueLengthDist<T> > curves;
    for (std::size_t a = 0; a < ists.size(); ++a)
        curves.push_back(line::nc::solver_nc_getprob_marg(sn, opt, ists[a]));

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "ProbMargAggr";
        p["indexBase"] = 0;
        line::reg::Json arr = line::reg::Json::array();
        for (std::size_t a = 0; a < ists.size(); ++a) {
            line::reg::Json e = line::reg::Json::object();
            e["station"] = ists[a] - 1;
            e["Station"] = sn.stations[ists[a] - 1].name;
            e["P"] = vector_json<T>(curves[a].P);
            e["logP"] = vector_json<T>(curves[a].logP);
            arr.push_back(e);
        }
        p["curves"] = arr;
        emit_analysis<T>("marg", p, opt.method);
        return 0;
    }
    std::printf("SolverNC arith=%s method=%s type=%s\n", line::num_traits<T>::name(),
                opt.method.c_str(), line::util::method_type("NC", opt.method).c_str());
    for (std::size_t a = 0; a < ists.size(); ++a) {
        std::printf("%-16s %-8s %14s %14s\n", "Station", "n", "P", "logP");
        for (std::size_t n = 0; n < curves[a].P.size(); ++n)
            std::printf("%-16s %-8zu %14.10g %14.10g\n", sn.stations[ists[a] - 1].name.c_str(), n,
                        line::num_traits<T>::to_double(curves[a].P[n]),
                        line::num_traits<T>::to_double(curves[a].logP[n]));
    }
    return 0;
}

/**
 * `-s nc -a cdf`: `@@SolverNC/getCdfRespT.m` and its aliases `getSjrnT`/`sjrnT`.
 *
 * THE WHOLE LAW, NOT ITS MEAN. `-a avg` reports E[R]; this reports F(t) per
 * (station, class) on one shared logarithmic grid, so a percentile or a tail
 * probability can be read off it. The algorithm is `pfqn_stdf` (`--method-cdf
 * exact`, the default) or the `pfqn_stdf_heur` reduction (`rd`), selected
 * through `options.config.algorithm` exactly as in the reference.
 *
 * FCFS ONLY, and the reference says so by WARNING and returning an empty
 * result rather than raising: the sojourn law of a processor-sharing or
 * infinite-server station is not the one this inversion computes. That warning
 * is carried through to stderr here and the analysis reports no curve, which is
 * distinguishable from a curve that is flat.
 */
template <class T>
int solve_model_nc_cdf(const std::string& file, const Knobs& k) {
    if constexpr (!line::num_traits<T>::has_transcendental) {
        throw line::UnsupportedError(
            "the -s nc -a cdf analysis evaluates the sojourn law on a logarithmic time grid and "
            "inverts a generating function; it needs transcendental arithmetic, so rerun with "
            "--arith double or --arith real");
    } else {
        line::qn::Network<T> net = read_model<T>(file);
        line::nc::NcSolverOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        if (!k.cdf_algorithm.empty()) opt.cdf_algorithm = k.cdf_algorithm;
        const line::qn::NetworkStruct<T>& sn = net.get_struct();
        const line::nc::CdfRespTResult<T> r = line::nc::solver_nc_cdf_respt(sn, opt);
        if (!r.warning.empty()) std::fprintf(stderr, "warning: %s\n", r.warning.c_str());

        if (g_json_output) {
            line::reg::Json p = line::reg::Json::object();
            p["type"] = "CdfRespT";
            p["indexBase"] = 0;
            p["algorithm"] = opt.cdf_algorithm;
            // ONE OBJECT PER CURVE, as on the CTMC arm: the pairs that carry a
            // law are a subset of the grid, so a column form would leave the
            // host to re-group them.
            line::reg::Json rd = line::reg::Json::array();
            for (std::size_t i = 0; i < r.RD.size(); ++i)
                for (std::size_t c = 0; c < r.RD[i].size(); ++c) {
                    // An empty entry is an ABSENT law -- the station is not FCFS
                    // or does not serve the class -- and is omitted rather than
                    // sent as a degenerate one.
                    if (r.RD[i][c].empty()) continue;
                    line::reg::Json e = line::reg::Json::object();
                    e["Station"] = sn.stations[i].name;
                    e["JobClass"] = sn.classes[c].name;
                    e["station"] = i;
                    e["jobclass"] = c;
                    line::reg::Json tt = line::reg::Json::array(), ff = line::reg::Json::array();
                    for (std::size_t j = 0; j < r.RD[i][c].rows(); ++j) {
                        ff.push_back(line::num_traits<T>::to_double(r.RD[i][c](j, 0)));
                        tt.push_back(line::num_traits<T>::to_double(r.RD[i][c](j, 1)));
                    }
                    e["t"] = tt;
                    e["F"] = ff;
                    rd.push_back(e);
                }
            p["respt"] = rd;
            p["tset"] = vector_json(r.tset);
            if (!r.warning.empty()) p["warning"] = r.warning;
            // NO "method": the law comes from the sojourn-time inversion and not
            // from the normalizing-constant ladder, so the requested method name
            // would not be the algorithm that produced these numbers.
            emit_analysis<T>("cdf", p, std::string());
            return 0;
        }
        std::printf("SolverNC arith=%s algorithm=%s grid=%zu\n", line::num_traits<T>::name(),
                    opt.cdf_algorithm.c_str(), r.tset.size());
        std::printf("%-16s %-14s %14s %14s\n", "Station", "JobClass", "Time", "F(t)");
        for (std::size_t i = 0; i < r.RD.size(); ++i)
            for (std::size_t c = 0; c < r.RD[i].size(); ++c) {
                if (r.RD[i][c].empty()) continue;
                for (std::size_t j = 0; j < r.RD[i][c].rows(); ++j)
                    std::printf("%-16s %-14s %14.8g %14.10g\n", sn.stations[i].name.c_str(),
                                sn.classes[c].name.c_str(),
                                line::num_traits<T>::to_double(r.RD[i][c](j, 1)),
                                line::num_traits<T>::to_double(r.RD[i][c](j, 0)));
            }
        return 0;
    }
}

/**
 * The clean-up a sensitivity table applies before printing, for a quantity that
 * may legitimately be NEGATIVE.
 *
 * `ln_sanitize` below tests `x <= FineTol`, which is right for a queue length or
 * a utilization -- every metric it was written for is nonnegative, so that test
 * reads as "negligible". A DERIVATIVE is not: raising a service rate lowers the
 * response time, the queue length and the utilization, so the whole sensitivity
 * table is negative by construction and the unsigned test would print an exact
 * zero for every one of those columns. Only the MAGNITUDE decides negligibility
 * here. The NC and the layered tables share it, which is why it sits above both.
 */
double sens_sanitize_signed(double x) {
    if (std::fabs(x) <= line::lang::GlobalConstants::FineTol) return 0.0;
    return x;
}

/**
 * `-s nc -a sens`: `@@NetworkSolver/getSensitivityTable.m` under SolverNC.
 *
 * ONE ROW PER (station, class) carrying dTput/dRate, dRespT/dRate, dQLen/dRate
 * and dUtil/dRate, i.e. the derivative of that row's means with respect to that
 * row's service RATE. Two branches produce them and the banner names the one
 * that ran: `exact` differentiates the product-form recursion analytically
 * (`pfqn_sens` at chain level for a closed model, the closed-form BCMP
 * derivatives for an open one), `fd` re-solves rate-perturbed copies of the
 * model with THIS solver and forms the quotient.
 *
 * NC IS ONE OF THE TWO ENGINES THAT CAN TAKE THE EXACT BRANCH, which is what
 * `@@SolverNC/supportsExactSensitivity.m` returns true for, so `auto` resolves
 * to `exact` whenever the model is in its scope (single-server queues plus
 * delays, not mixed) and only falls back to differences outside it. Asking for
 * `--sens-method exact` outside that scope is refused by name rather than
 * silently downgraded: the two branches answer to different precision.
 *
 * The struct is COPIED rather than referenced because the fd branch writes a
 * scaled service process into it between solves; `net.get_struct()` hands out a
 * const reference to the model's own, which must not move under the caller.
 */
template <class T>
int solve_model_nc_sens(const std::string& file, const Knobs& k) {
    line::qn::Network<T> net = read_model<T>(file);
    line::nc::NcSolverOptions opt;
    if (!k.method.empty() && k.method != "default") opt.method = k.method;
    if (k.tol >= 0.0) opt.tol = k.tol;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    line::qn::NetworkStruct<T> sn = net.get_struct();

    line::sens::SensOptions so;
    if (!k.sens_method.empty()) so.method = k.sens_method;
    if (!k.sens_scheme.empty()) so.scheme = k.sens_scheme;
    if (k.sens_step > 0.0) so.step = k.sens_step;
    so.simulation = false;  // the normalizing-constant path is deterministic

    // `getAvg`, not the raw analyzer: the reference's difference quotient is
    // taken on the metrics the solver reports, which are the filtered ones.
    const line::sens::SensTable<T> tbl = line::sens::solver_sensitivity_table<T>(
        sn, so, /*exact_available=*/true, [&sn, &opt]() {
            const line::mva::AvgResult<T> a = line::nc::solver_nc_run_analyzer(sn, opt);
            line::mva::MvaSolution<T> s;
            s.Q = a.QN;
            s.U = a.UN;
            s.R = a.RN;
            s.Tp = a.TN;
            s.C = a.CN;
            s.X = a.XN;
            s.method = a.actualmethod;
            s.iter = a.iter;
            return s;
        });

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "SensitivityTable";
        p["indexBase"] = 0;
        p["branch"] = tbl.method;
        line::reg::Json rows = line::reg::Json::array();
        for (const line::sens::SensRow<T>& r : tbl.rows) {
            line::reg::Json o = line::reg::Json::object();
            o["Station"] = r.station;
            o["JobClass"] = r.jobclass;
            o["dTput_dRate"] = sens_sanitize_signed(line::num_traits<T>::to_double(r.dTput));
            o["dRespT_dRate"] = sens_sanitize_signed(line::num_traits<T>::to_double(r.dRespT));
            o["dQLen_dRate"] = sens_sanitize_signed(line::num_traits<T>::to_double(r.dQLen));
            o["dUtil_dRate"] = sens_sanitize_signed(line::num_traits<T>::to_double(r.dUtil));
            rows.push_back(o);
        }
        p["rows"] = rows;
        // The branch rides in "branch" because it is not a normalizing-constant
        // method name. The method slot is EMPTY on the exact branch, which
        // differentiates the recursion in closed form and never runs a solve:
        // naming the NC method there would attribute the numbers to an
        // algorithm that did not produce them.
        emit_analysis<T>("sens", p, tbl.method == "fd" ? opt.method : std::string());
        return 0;
    }

    std::printf("SolverNC arith=%s branch=%s method=%s rows=%zu\n", line::num_traits<T>::name(),
                tbl.method.c_str(), tbl.method == "fd" ? opt.method.c_str() : "-",
                tbl.rows.size());
    std::printf("%-16s %-14s %14s %14s %14s %14s\n", "Station", "JobClass", "dTput_dRate",
                "dRespT_dRate", "dQLen_dRate", "dUtil_dRate");
    for (const line::sens::SensRow<T>& r : tbl.rows)
        std::printf("%-16s %-14s %14.6g %14.6g %14.6g %14.6g\n", r.station.c_str(),
                    r.jobclass.c_str(),
                    sens_sanitize_signed(line::num_traits<T>::to_double(r.dTput)),
                    sens_sanitize_signed(line::num_traits<T>::to_double(r.dRespT)),
                    sens_sanitize_signed(line::num_traits<T>::to_double(r.dQLen)),
                    sens_sanitize_signed(line::num_traits<T>::to_double(r.dUtil)));
    return 0;
}

/** The @@Solver method each `-a` stands for, which is what the chooser keys on. */
// ===================== SolverLDES, the simulator ==========================

/**
 * The knobs one `-s ldes` invocation resolves to.
 *
 * `--samples` is the SERVICE-COMPLETION budget the engine stops at, `--seed` the
 * stream, and both are part of the answer rather than of the invocation, which is
 * why the banner carries them. Everything else is `--ldes-*`.
 *
 * A `--ldes-initsol` placement is NOT accompanied by a forced `fixed` warmup
 * filter here, deliberately: `initFromSolver` sets `tranfilter='fixed'` with
 * `warmupfrac=0` because the placement it computes IS a steady state, while a
 * placement handed in on the command line may equally be the start of a
 * transient. Pass `--ldes-tranfilter fixed --ldes-warmupfrac 0` alongside it to
 * reproduce `initFromSolver` exactly.
 */
inline line::ldes::LdesOptions ldes_options(const Knobs& k) {
    line::ldes::LdesOptions o;
    if (k.samples) o.samples = k.samples;
    if (k.seed) o.seed = static_cast<long>(k.seed);
    if (!k.method.empty()) o.method = k.method;
    if (!k.ldes_tranfilter.empty()) o.tranfilter = k.ldes_tranfilter;
    if (k.ldes_warmupfrac >= 0.0) o.warmupfrac = k.ldes_warmupfrac;
    if (!k.ldes_cimethod.empty()) o.cimethod = k.ldes_cimethod;
    if (k.ldes_cnvgon) o.cnvgon = true;
    if (k.ldes_cnvgtol > 0.0) o.cnvgtol = k.ldes_cnvgtol;
    if (k.ldes_slotted) o.slotted = true;
    if (k.ldes_slotlength > 0.0) {
        o.slotted = true;
        o.slot_length = k.ldes_slotlength;
    }
    if (k.ldes_replications > 0) o.replications = k.ldes_replications;
    if (k.ldes_numthreads > 0) o.numthreads = k.ldes_numthreads;
    if (k.ldes_maxtime > 0.0) o.timeout = k.ldes_maxtime;
    if (!k.ldes_initsol.empty()) o.init_sol = k.ldes_initsol;
    if (!k.ldes_rest_url.empty()) o.rest_url = k.ldes_rest_url;
    o.verbose = k.verbose;
    return o;
}

/**
 * The model.json text the engine is handed.
 *
 * FORWARDED BYTE FOR BYTE, and never through `read_model`: the reader is scoped
 * to the subset the analytical solvers need, and a round trip through it would
 * degrade exactly the models LDES exists for. `-a reward` is the one arm that
 * also parses the document, because a reward DECLARATION is what it needs.
 */
inline std::string ldes_document(const std::string& file) {
    return file.empty() ? stdin_model_text() : line::ldes::detail::read_file(file);
}

/** A reported entry, or 0 where the engine reported no such row. */
inline double ldes_at(const line::Matrix<double>& M, std::size_t i, std::size_t j) {
    return i < M.rows() && j < M.cols() ? M(i, j) : 0.0;
}

/**
 * One LDES run.
 *
 * A HARD TIMEOUT IS A REFUSAL HERE, not an empty result. The two other clients
 * return an empty result flagged `timedOut` and warn, which suits a caller that
 * can inspect the flag; a CLI's caller reads a table, and a table of zeros that
 * means "the run was killed" is the silent-wrong-number outcome this CLI refuses
 * everywhere else.
 */
inline line::ldes::LdesResult ldes_run(const std::string& file,
                                      const line::ldes::LdesOptions& o,
                                      const std::vector<std::string>& extra) {
    const line::ldes::LdesResult r = line::ldes::solver_ldes_text(ldes_document(file), o, extra);
    if (r.timed_out)
        throw line::NumericError(
            "SolverLDES exceeded its wall-clock budget (--ldes-maxtime) and was terminated before "
            "it wrote a result; raise the budget or lower --samples");
    if (r.station_names.empty())
        throw line::NumericError(
            "SolverLDES: the engine reported no station names, so its metrics cannot be labelled; "
            "the run produced no result document");
    return r;
}

/**
 * The provenance line every LDES arm prints first.
 *
 * `engine=` is not decoration: the AOT native image and the jar are two builds of
 * one engine and the first can lag the sources, so a number quoted from an LDES
 * run has to say which produced it. `stopping=` is the reason the run ended --
 * a `max_events` stop at a low `--samples` is a wide confidence interval and a
 * `max_time` one is a truncated run, and neither is visible in the means.
 */
inline void ldes_banner(const line::ldes::LdesResult& r, const line::ldes::LdesOptions& o) {
    std::printf("SolverLDES arith=double method=%s type=%s engine=%s samples=%zu seed=%ld "
                "time=%.6g events=%lld stopping=%s\n",
                r.method.c_str(), line::util::method_type("LDES", r.method).c_str(),
                r.engine.c_str(), o.events ? o.events : o.samples, o.seed, r.runtime,
                r.total_simulated_events, r.stopping_reason.c_str());
}

/** The envelope keys that qualify an LDES solve as a whole. */
inline line::reg::Json ldes_envelope(const line::ldes::LdesResult& r,
                                     const line::ldes::LdesOptions& o) {
    line::reg::Json e = line::reg::Json::object();
    e["engine"] = r.engine;
    e["samples"] = o.events ? o.events : o.samples;
    e["seed"] = o.seed;
    e["converged"] = r.converged;
    e["stoppingReason"] = r.stopping_reason;
    e["totalSimulatedEvents"] = r.total_simulated_events;
    e["runtime"] = r.runtime;
    return e;
}

/**
 * `-s ldes -a avg`: the steady-state table, the engine's `getAvg`.
 *
 * THE FINITE-CAPACITY-REGION ROWS DO NOT JOIN THE STATION TABLE, unlike MATLAB's
 * `getAvgTable`, which appends them after the stations. A region is not a station
 * and the "Station" column of this CLI's table is read by a parity harness that
 * pairs rows with another codebase's stations; a region row there would pair with
 * nothing. They are printed as their own table and carried under `avg.fcr`, which
 * is the same information without the collision.
 */
int solve_model_ldes_avg(const std::string& file, const Knobs& k) {
    const line::ldes::LdesOptions o = ldes_options(k);
    const line::ldes::LdesResult r = ldes_run(file, o, std::vector<std::string>());
    if (!g_json_output) ldes_banner(r, o);

    line::reg::Json extra = line::reg::Json::object();
    // The confidence intervals are the half-widths the engine reports, one per
    // metric; a simulation that quoted a mean without them would be quoting a
    // point estimate as if it were exact.
    line::reg::Json ci = line::reg::Json::object();
    if (!r.QNCI.empty()) ci["QNCI"] = matrix_json<double>(r.QNCI);
    if (!r.UNCI.empty()) ci["UNCI"] = matrix_json<double>(r.UNCI);
    if (!r.RNCI.empty()) ci["RNCI"] = matrix_json<double>(r.RNCI);
    if (!r.TNCI.empty()) ci["TNCI"] = matrix_json<double>(r.TNCI);
    if (!r.ANCI.empty()) ci["ANCI"] = matrix_json<double>(r.ANCI);
    if (!r.WNCI.empty()) ci["WNCI"] = matrix_json<double>(r.WNCI);
    if (!ci.empty()) extra["CI"] = ci;
    if (!r.QNfcr.empty()) {
        line::reg::Json f = line::reg::Json::object();
        f["nregions"] = r.nregions;
        f["QNfcr"] = matrix_json<double>(r.QNfcr);
        f["RNfcr"] = matrix_json<double>(r.RNfcr);
        f["TNfcr"] = matrix_json<double>(r.TNfcr);
        f["WNfcr"] = matrix_json<double>(r.WNfcr);
        if (!r.WeightNfcr.empty()) f["WeightNfcr"] = matrix_json<double>(r.WeightNfcr);
        if (!r.MemOccNfcr.empty()) f["MemOccNfcr"] = matrix_json<double>(r.MemOccNfcr);
        if (!r.DropRateNfcr.empty()) f["DropRateNfcr"] = matrix_json<double>(r.DropRateNfcr);
        extra["fcr"] = f;
    }
    if (!r.DropRateJoin.empty()) extra["DropRateJoin"] = matrix_json<double>(r.DropRateJoin);
    if (!r.cache_metrics.empty()) {
        line::reg::Json cm = line::reg::Json::object();
        for (std::map<std::string, line::ldes::LdesCacheMetrics>::const_iterator it =
                 r.cache_metrics.begin();
             it != r.cache_metrics.end(); ++it) {
            line::reg::Json c = line::reg::Json::object();
            if (!it->second.hit.empty()) c["hit"] = matrix_json<double>(it->second.hit);
            if (!it->second.delayed.empty()) c["delayed"] = matrix_json<double>(it->second.delayed);
            if (!it->second.miss.empty()) c["miss"] = matrix_json<double>(it->second.miss);
            if (!it->second.latency.empty()) c["latency"] = matrix_json<double>(it->second.latency);
            if (!it->second.hitList.empty()) c["hitList"] = matrix_json<double>(it->second.hitList);
            if (!it->second.itemProb.empty())
                c["itemProb"] = matrix_json<double>(it->second.itemProb);
            if (!it->second.listCost.empty())
                c["listCost"] = matrix_json<double>(it->second.listCost);
            cm[it->first] = c;
        }
        extra["cacheMetrics"] = cm;
    }

    // THE RESIDENCE TIME IS DERIVED HERE, not taken from the engine. The engine
    // reports WN = RN because it counts one visit per station, which is only
    // true when every visit ratio is 1; `getAvg.m:204` therefore discards the
    // WN a solver returned and recomputes `sn_get_residt_from_respt(sn, RN)`,
    // and every other C++ solver already routes through the same helper. On
    // cqn_repairmen, whose Queue1 is visited 0.3 times per cycle, the engine's
    // WN came out 11.8136 against the reference's 3.5205 -- the response time
    // reported as if the station were visited once.
    //
    // MATCHED BY NAME. The engine's station order is its own; a Cache or a
    // Source can sit at a different index in the struct, and pairing the two
    // off positionally would scale one station's time by another's visits.
    line::Matrix<double> WNd(r.station_names.size(), r.class_names.size(), 0.0);
    {
        line::qn::Network<double> net = read_model<double>(file);
        const line::qn::NetworkStruct<double>& sn = net.get_struct();
        std::vector<std::size_t> st_of(r.station_names.size(), 0);  // 1-based, 0 = unmatched
        for (std::size_t i = 0; i < r.station_names.size(); ++i)
            for (std::size_t j = 0; j < sn.nstations; ++j)
                if (sn.stations[j].name == r.station_names[i]) { st_of[i] = j + 1; break; }
        std::vector<std::size_t> cl_of(r.class_names.size(), 0);
        for (std::size_t c = 0; c < r.class_names.size(); ++c)
            for (std::size_t k = 0; k < sn.nclasses; ++k)
                if (sn.classes[k].name == r.class_names[c]) { cl_of[c] = k + 1; break; }
        line::Matrix<double> RNs(sn.nstations, sn.nclasses, 0.0);
        for (std::size_t i = 0; i < r.station_names.size(); ++i)
            for (std::size_t c = 0; c < r.class_names.size(); ++c)
                if (st_of[i] && cl_of[c]) RNs(st_of[i] - 1, cl_of[c] - 1) = ldes_at(r.RN, i, c);
        const line::Matrix<double> WNs = line::mva::sn_get_residt_from_respt(sn, RNs);
        for (std::size_t i = 0; i < r.station_names.size(); ++i)
            for (std::size_t c = 0; c < r.class_names.size(); ++c)
                WNd(i, c) = (st_of[i] && cl_of[c]) ? WNs(st_of[i] - 1, cl_of[c] - 1)
                                                   : ldes_at(r.WN, i, c);
    }

    emit_avg_table_named(r.station_names, r.class_names, "double", r.method,
                         [&](std::size_t i, std::size_t c) {
                             AvgRow v;
                             v.q = ldes_at(r.QN, i, c);
                             v.u = ldes_at(r.UN, i, c);
                             v.r = ldes_at(r.RN, i, c);
                             v.w = WNd(i, c);
                             v.a = ldes_at(r.AN, i, c);
                             v.t = ldes_at(r.TN, i, c);
                             return v;
                         },
                         extra, ldes_envelope(r, o));

    if (!g_json_output && !r.QNfcr.empty()) {
        std::printf("%-16s %-14s %12s %12s %12s %12s %12s\n", "Region", "JobClass", "QLen", "RespT",
                    "Tput", "Weight", "MemOcc");
        for (std::size_t i = 0; i < r.QNfcr.rows(); ++i)
            for (std::size_t c = 0; c < r.class_names.size(); ++c)
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                            ("Region" + std::to_string(i + 1)).c_str(), r.class_names[c].c_str(),
                            ldes_at(r.QNfcr, i, c), ldes_at(r.RNfcr, i, c),
                            ldes_at(r.TNfcr, i, c), ldes_at(r.WeightNfcr, i, c),
                            ldes_at(r.MemOccNfcr, i, c));
    }
    if (!g_json_output && !r.cache_metrics.empty()) {
        std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Cache", "JobClass", "Hit", "Delayed",
                    "Miss", "Latency");
        for (std::map<std::string, line::ldes::LdesCacheMetrics>::const_iterator it =
                 r.cache_metrics.begin();
             it != r.cache_metrics.end(); ++it)
            for (std::size_t c = 0; c < r.class_names.size(); ++c)
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n", it->first.c_str(),
                            r.class_names[c].c_str(), ldes_at(it->second.hit, 0, c),
                            ldes_at(it->second.delayed, 0, c), ldes_at(it->second.miss, 0, c),
                            ldes_at(it->second.latency, 0, c));
    }
    return 0;
}

/**
 * `-s ldes -a tran`: `getTranAvg`, the per-bucket QNt / UNt / TNt series over
 * `--tspan`.
 *
 * THE HORIZON IS REQUIRED. `options.timespan` is what turns the engine's run into
 * a transient one, and there is no default: a trajectory over an unstated horizon
 * is not a quantity. A SINGLE PATH IS NOT E[N](t) either -- there is no time
 * ergodicity at fixed t -- so `--ldes-replications` is how an ensemble mean is
 * asked for, exactly as `runAnalyzer.m` passes `--replications` for the same
 * reason.
 *
 * The series are indexed by STATION, as `LDESResultIO` writes them
 * (`result.QNt = new Matrix[numStations][numClasses]`).
 */
int solve_model_ldes_tran(const std::string& file, const Knobs& k) {
    line::ldes::LdesOptions o = ldes_options(k);
    o.has_timespan = true;
    o.t0 = k.t0;
    o.t1 = k.t1;
    std::vector<std::string> extra;
    extra.push_back("--trajectory");
    const line::ldes::LdesResult r = ldes_run(file, o, extra);
    if (r.t.empty() || r.QNt.empty())
        throw line::NumericError(
            "SolverLDES -a tran produced no trajectory: the engine ran but recorded no bucket over "
            "[" + line::ldes::detail::shortest(k.t0) + "," +
            line::ldes::detail::shortest(k.t1) + "]");

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "TranAvgTable";
        p["indexBase"] = 0;
        p["t0"] = k.t0;
        p["t1"] = k.t1;
        line::reg::Json curves = line::reg::Json::array();
        for (std::size_t i = 0; i < r.QNt.size(); ++i)
            for (std::size_t c = 0; c < r.QNt[i].size(); ++c) {
                // An empty series is an ABSENT one -- the class does not visit
                // the station -- and is omitted rather than sent as a flat zero.
                if (r.QNt[i][c].empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = i < r.station_names.size() ? r.station_names[i]
                                                          : "Station" + std::to_string(i);
                e["JobClass"] =
                    c < r.class_names.size() ? r.class_names[c] : "Class" + std::to_string(c);
                e["station"] = i;
                e["jobclass"] = c;
                line::reg::Json tt = line::reg::Json::array(), q = line::reg::Json::array(),
                                u = line::reg::Json::array(), x = line::reg::Json::array();
                for (std::size_t j = 0; j < r.QNt[i][c].rows(); ++j) {
                    tt.push_back(r.QNt[i][c](j, 1));
                    q.push_back(r.QNt[i][c](j, 0));
                }
                if (i < r.UNt.size() && c < r.UNt[i].size())
                    for (std::size_t j = 0; j < r.UNt[i][c].rows(); ++j)
                        u.push_back(r.UNt[i][c](j, 0));
                if (i < r.TNt.size() && c < r.TNt[i].size())
                    for (std::size_t j = 0; j < r.TNt[i][c].rows(); ++j)
                        x.push_back(r.TNt[i][c](j, 0));
                e["t"] = tt;
                e["QLen"] = q;
                e["Util"] = u;
                e["Tput"] = x;
                curves.push_back(e);
            }
        p["curves"] = curves;
        p["tset"] = vector_json(r.t);
        // The envelope is built ONCE into a local: `begin()` and `end()` taken
        // from two different temporaries are iterators into two different
        // objects, which is undefined behaviour and not a style point.
        const line::reg::Json env = ldes_envelope(r, o);
        for (line::reg::Json::const_iterator it = env.begin(); it != env.end(); ++it)
            p[it.key()] = it.value();
        emit_analysis<double>("tran", p, r.method);
        return 0;
    }
    ldes_banner(r, o);
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "Time", "QLen", "Util",
                "Tput");
    for (std::size_t i = 0; i < r.QNt.size(); ++i)
        for (std::size_t c = 0; c < r.QNt[i].size(); ++c) {
            if (r.QNt[i][c].empty()) continue;
            for (std::size_t j = 0; j < r.QNt[i][c].rows(); ++j) {
                const bool hu = i < r.UNt.size() && c < r.UNt[i].size() &&
                                j < r.UNt[i][c].rows();
                const bool hx = i < r.TNt.size() && c < r.TNt[i].size() &&
                                j < r.TNt[i][c].rows();
                std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n",
                            (i < r.station_names.size() ? r.station_names[i].c_str() : "?"),
                            (c < r.class_names.size() ? r.class_names[c].c_str() : "?"),
                            r.QNt[i][c](j, 1), r.QNt[i][c](j, 0),
                            hu ? r.UNt[i][c](j, 0) : 0.0, hx ? r.TNt[i][c](j, 0) : 0.0);
            }
        }
    return 0;
}

/**
 * `-s ldes -a cdf`: `getCdfRespT`, the EMPIRICAL response-time law.
 *
 * A SIMULATOR MUST REPORT WHAT IT OBSERVED. The base solver's fallback fabricates
 * an exponential law with the right mean, which says nothing about the tail; the
 * engine records every per-job response time under `--respt-samples`, and the
 * curve here is the ecdf of those samples with repeated observations collapsed to
 * their largest F, exactly as `@@SolverLDES/getCdfRespT.m` builds it.
 *
 * `--respt-samples` POSTDATES the prebuilt AOT image, which is why the runner
 * order flips for it (see `ldes_runners`): on that image the flag is accepted and
 * ignored, and the arm would refuse for want of samples that were never asked for.
 */
int solve_model_ldes_cdf(const std::string& file, const Knobs& k, const char* key,
                         const char* type) {
    const line::ldes::LdesOptions o = ldes_options(k);
    std::vector<std::string> extra;
    extra.push_back("--respt-samples");
    const line::ldes::LdesResult r = ldes_run(file, o, extra);
    if (r.respTimeSamples.empty())
        throw line::NumericError(
            "SolverLDES -a cdf needs the per-job response times the engine records under "
            "--respt-samples and the run returned none; raise --samples so completions are "
            "observed at all");

    // The ecdf of each (station, class): sorted observations, F = i/n, and one
    // pair per DISTINCT value carrying the largest F at it.
    std::vector<std::vector<std::vector<double>>> tt(r.respTimeSamples.size()), ff(
        r.respTimeSamples.size());
    for (std::size_t i = 0; i < r.respTimeSamples.size(); ++i) {
        tt[i].resize(r.respTimeSamples[i].size());
        ff[i].resize(r.respTimeSamples[i].size());
        for (std::size_t c = 0; c < r.respTimeSamples[i].size(); ++c) {
            std::vector<double> x = r.respTimeSamples[i][c];
            if (x.empty()) continue;
            std::sort(x.begin(), x.end());
            const double n = static_cast<double>(x.size());
            for (std::size_t j = 0; j < x.size(); ++j) {
                if (j + 1 < x.size() && x[j + 1] == x[j]) continue;
                tt[i][c].push_back(x[j]);
                ff[i][c].push_back(static_cast<double>(j + 1) / n);
            }
        }
    }

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = type;
        p["indexBase"] = 0;
        p["algorithm"] = "empirical";
        line::reg::Json rd = line::reg::Json::array();
        for (std::size_t i = 0; i < tt.size(); ++i)
            for (std::size_t c = 0; c < tt[i].size(); ++c) {
                if (tt[i][c].empty()) continue;
                line::reg::Json e = line::reg::Json::object();
                e["Station"] = i < r.station_names.size() ? r.station_names[i]
                                                          : "Station" + std::to_string(i);
                e["JobClass"] =
                    c < r.class_names.size() ? r.class_names[c] : "Class" + std::to_string(c);
                e["station"] = i;
                e["jobclass"] = c;
                e["t"] = vector_json(tt[i][c]);
                e["F"] = vector_json(ff[i][c]);
                e["samples"] = r.respTimeSamples[i][c].size();
                rd.push_back(e);
            }
        p["respt"] = rd;
        emit_analysis<double>(key, p, std::string());
        return 0;
    }
    ldes_banner(r, o);
    std::printf("%-16s %-14s %14s %14s\n", "Station", "JobClass", "Time", "F(t)");
    for (std::size_t i = 0; i < tt.size(); ++i)
        for (std::size_t c = 0; c < tt[i].size(); ++c)
            for (std::size_t j = 0; j < tt[i][c].size(); ++j)
                std::printf("%-16s %-14s %14.8g %14.10g\n",
                            (i < r.station_names.size() ? r.station_names[i].c_str() : "?"),
                            (c < r.class_names.size() ? r.class_names[c].c_str() : "?"),
                            tt[i][c][j], ff[i][c][j]);
    return 0;
}

/**
 * `-s ldes -a sample`: `sampleSys` / `sampleSysAggr`, one simulated trajectory.
 *
 * The horizon is `[0, --samples]`, which is `runTransientJson`'s: the event budget
 * doubles as the transient horizon there because the engine ignores the budget in
 * transient mode, so the number the caller gave has to name the horizon or name
 * nothing. `--tspan` names a horizon in its own right and belongs to `-a tran`.
 *
 * The state is the per-class queue length at each station, which is why there is
 * no separate `sampleSysAggr` column: an LDES trajectory is ALREADY per class, so
 * the aggregate view is the row sum and the reference's two getters return the
 * same data with one flag flipped.
 */
int solve_model_ldes_sample(const std::string& file, const Knobs& k) {
    line::ldes::LdesOptions o = ldes_options(k);
    o.has_timespan = true;
    o.t0 = 0.0;
    o.t1 = static_cast<double>(o.events ? o.events : o.samples);
    std::vector<std::string> extra;
    extra.push_back("--trajectory");
    const line::ldes::LdesResult r = ldes_run(file, o, extra);
    if (r.t.empty() || r.QNt.empty())
        throw line::NumericError(
            "SolverLDES -a sample produced no trajectory over [0," +
            line::ldes::detail::shortest(o.t1) + "]");

    const std::size_t M = r.QNt.size(), K = r.class_names.size(), n = r.t.size();
    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "SamplePath";
        p["indexBase"] = 0;
        p["scope"] = "(system)";
        p["drawn"] = n;
        p["t"] = vector_json(r.t);
        line::reg::Json st = line::reg::Json::array();
        for (std::size_t j = 0; j < n; ++j) {
            line::reg::Json row = line::reg::Json::array();
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < K; ++c)
                    row.push_back(c < r.QNt[i].size() && j < r.QNt[i][c].rows()
                                      ? r.QNt[i][c](j, 0)
                                      : 0.0);
            st.push_back(row);
        }
        p["state"] = st;
        line::reg::Json cols = line::reg::Json::array();
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < K; ++c)
                cols.push_back((i < r.station_names.size() ? r.station_names[i] : "?") + "," +
                               (c < K ? r.class_names[c] : "?"));
        p["columns"] = cols;
        // The envelope is built ONCE into a local: `begin()` and `end()` taken
        // from two different temporaries are iterators into two different
        // objects, which is undefined behaviour and not a style point.
        const line::reg::Json env = ldes_envelope(r, o);
        for (line::reg::Json::const_iterator it = env.begin(); it != env.end(); ++it)
            p[it.key()] = it.value();
        emit_analysis<double>("sample", p, r.method);
        return 0;
    }
    ldes_banner(r, o);
    std::printf("%14s   %s\n", "Time", "SysState (station-major, per class)");
    for (std::size_t j = 0; j < n; ++j) {
        std::printf("%14.8g  ", r.t[j]);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < K; ++c)
                std::printf(" %g", c < r.QNt[i].size() && j < r.QNt[i][c].rows()
                                       ? r.QNt[i][c](j, 0)
                                       : 0.0);
        std::printf("\n");
    }
    return 0;
}

/**
 * `-s ldes -a reward`: `getAvgReward`, E[r] on the EXACT joint-state histogram.
 *
 * The engine exports, under `--export-histogram`, the residence time of every
 * joint state it visited, in the aggregate layout `ctmc_state_space_aggr` builds.
 * The rewards are then evaluated here, state by state, so E[r] = sum_s (t_s /
 * sum t) r(state_s) is correct for a NONLINEAR reward too -- which is the whole
 * point of the histogram over the means: E[n^2] cannot be recovered from E[n].
 *
 * THIS IS THE ONE ARM THAT ALSO PARSES THE DOCUMENT, because a reward is a
 * DECLARATION and the declarations live in the model, not in the result. A model
 * outside this port's reader therefore reaches every other LDES arm and not this
 * one, and says so by the reader's own refusal.
 */
int solve_model_ldes_reward(const std::string& file, const Knobs& k) {
    line::qn::Network<double> net = read_model<double>(file);
    const line::qn::NetworkStruct<double>& sn = net.get_struct();
    if (sn.reward.empty())
        throw line::InputError(
            "-s ldes -a reward needs a reward declared on the model (set_reward(name, fn), the "
            "`rewards` block of model.json); there is nothing to average");

    const line::ldes::LdesOptions o = ldes_options(k);
    std::vector<std::string> extra;
    extra.push_back("--export-histogram");
    const line::ldes::LdesResult r = ldes_run(file, o, extra);
    if (r.histogram_space.empty() || r.histogram_time.empty())
        throw line::NumericError(
            "SolverLDES -a reward needs the joint-state residence-time histogram the engine "
            "exports under --export-histogram and the run returned none");

    double total = 0.0;
    for (std::size_t s = 0; s < r.histogram_time.rows(); ++s)
        for (std::size_t c = 0; c < r.histogram_time.cols(); ++c) total += r.histogram_time(s, c);
    if (!(total > 0.0))
        throw line::NumericError(
            "SolverLDES -a reward: the state histogram carries no residence time, so no state "
            "distribution can be formed from it");

    const std::size_t ns = r.histogram_space.rows(), w = r.histogram_space.cols();
    std::vector<double> E(sn.reward.size(), 0.0);
    std::vector<std::string> names(sn.reward.size());
    for (std::size_t l = 0; l < sn.reward.size(); ++l) {
        names[l] = sn.reward[l].name;
        for (std::size_t s = 0; s < ns; ++s) {
            std::vector<double> row(w);
            for (std::size_t c = 0; c < w; ++c) row[c] = r.histogram_space(s, c);
            const double t = s < r.histogram_time.rows() && r.histogram_time.cols() > 0
                                 ? r.histogram_time(s, 0)
                                 : 0.0;
            E[l] += (t / total) * sn.reward[l].fn(row);
        }
    }

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgReward";
        line::reg::Json nm = line::reg::Json::array();
        for (std::size_t l = 0; l < names.size(); ++l) nm.push_back(names[l]);
        p["Reward"] = nm;
        p["E"] = vector_json(E);
        p["states"] = ns;
        // The envelope is built ONCE into a local: `begin()` and `end()` taken
        // from two different temporaries are iterators into two different
        // objects, which is undefined behaviour and not a style point.
        const line::reg::Json env = ldes_envelope(r, o);
        for (line::reg::Json::const_iterator it = env.begin(); it != env.end(); ++it)
            p[it.key()] = it.value();
        // NO "method": these are expectations over an empirical distribution and
        // no chain was solved, so there is no resolved method to report.
        emit_analysis<double>("reward", p, std::string());
        return 0;
    }
    ldes_banner(r, o);
    std::printf("%-28s %16s\n", "Reward", "E[r]");
    for (std::size_t l = 0; l < E.size(); ++l)
        std::printf("%-28s %16.10g\n", names[l].c_str(), E[l]);
    return 0;
}

/**
 * The `-s ldes` entry: one analysis per getter of the reference's surface.
 *
 * `-a prob` IS REFUSED BY NAME, and the refusal is about this port's model layer
 * rather than about LDES. `getProb` / `getProbSys` weigh the simulated trajectory
 * by how long it spent in the model's CURRENT state, `sn.state{isf}` per stateful
 * node; `qn::NetworkStruct` carries no such row -- a Cache's `initstate`, a
 * Place's `initmarking` and the `stateprior`/`statespace` pair are declarations of
 * a different thing -- so there is no state to weigh against and any answer would
 * be about a state the caller never named.
 */
int solve_model_ldes(const std::string& file, const Knobs& k, const std::string& analysis) {
    if (analysis == "avg") return solve_model_ldes_avg(file, k);
    if (analysis == "tran") return solve_model_ldes_tran(file, k);
    // THE SAME ECDF UNDER FOUR NAMES, which is the reference's own structure and
    // not a shortcut here: `@@SolverLDES/getCdfRespT`, `getTranCdfRespT` and
    // `getTranCdfPassT` all read `respTimeSamples`, and the JAR's
    // getTranCdfPassT is literally `return getTranCdfRespT(R)`. A simulator
    // observes one per-job passage and there is no second measurement to make;
    // emitting the curve under the key the caller asked for is what tells them
    // apart, and the payload's `type` says which question it answers.
    if (analysis == "cdf") return solve_model_ldes_cdf(file, k, "cdf", "CdfRespT");
    if (analysis == "cdfpasst") return solve_model_ldes_cdf(file, k, "cdfpasst", "CdfPassT");
    if (analysis == "trancdf") return solve_model_ldes_cdf(file, k, "trancdf", "TranCdfRespT");
    if (analysis == "trancdfpasst")
        return solve_model_ldes_cdf(file, k, "trancdfpasst", "TranCdfPassT");
    if (analysis == "sample") return solve_model_ldes_sample(file, k);
    if (analysis == "reward") return solve_model_ldes_reward(file, k);
    if (analysis == "prob")
        throw line::UnsupportedError(
            "-s ldes -a prob is not ported: getProb/getProbSys weigh the trajectory by the time "
            "spent in the model's CURRENT state, and the C++ NetworkStruct carries no such state "
            "row to compare against (only a Cache initstate, a Place initmarking and the "
            "statePrior/space pair). Use -s ctmc -a prob for an exact marginal, or -s ssa -a prob "
            "for a simulated one");
    throw line::UnsupportedError(
        "SolverLDES ports -a avg (getAvg) and its four views -a node, -a sys, -a chain and "
        "-a nodechain, -a tran (getTranAvg), -a cdf / cdf-passt / "
        "tran-cdf-respt / tran-cdf-passt (the empirical passage law, one measurement under the "
        "four names the reference gives it), -a sample (sampleSys) and -a reward "
        "(getAvgReward); got '" + analysis + "'");
}

std::string auto_getter_of_analysis(const std::string& analysis) {
    if (analysis == "prob") return "getProbSysAggr";
    if (analysis == "marg") return "getProbMarg";
    if (analysis == "sysmarg") return "getProbSysMarg";
    if (analysis == "normconst") return "getProbNormConstAggr";
    if (analysis == "tranprob") return "getTranProbSysAggr";
    if (analysis == "sample") return "sampleSys";
    if (analysis == "cdf") return "getCdfRespT";
    if (analysis == "gen") return "getInfGen";
    if (analysis == "states") return "getStateSpace";
    if (analysis == "reward") return "getAvgReward";
    if (analysis == "sens") return "getSensitivityTable";
    if (analysis == "tran") return "getTranAvg";
    if (analysis == "internals") return "getMAMResult";
    // `-a bounds` and `-a tranreward` deliberately keep the AvgTable getter:
    // `getBoundsTable` and `getTranReward` are absent from chooseSolverHeur's
    // own method lists, so naming them here would ask the chooser about a getter
    // it does not rank. Each is served by exactly one engine anyway (BA and
    // CTMC), which refuses by name when `-s auto` sends the run elsewhere.
    if (analysis == "node") return "getAvgNodeTable";
    return "getAvgTable";
}

/**
 * What `-s auto` resolved to: the engines to try, in order, and the method the
 * first of them runs.
 */
struct AutoPlan {
    std::vector<std::string> order;  ///< CLI solver method names, chosen first
    std::string method;              ///< the method the chosen engine runs, "" for its default
    std::string note;                ///< what the ranking preferred and this port cannot build
};

/** The CLI method name of a method family, or a refusal naming what the family needs. */
std::string auto_cli_token_of_family(const std::string& fam) {
    if (fam == "mva" || fam == "nc" || fam == "ctmc" || fam == "mam" || fam == "ag" ||
        fam == "ssa" || fam == "ba" || fam == "uq" || fam == "env")
        return fam;
    if (fam == "fluid") return "fluid";
    if (fam == "ldes") {
        // The engine is not built here, it is RUN here, so the family resolves
        // whenever the machine has one and refuses -- naming what is missing --
        // when it does not, rather than diverting to a different simulator.
        if (!line::ldes::ldes_is_available())
            throw line::UnsupportedError(
                "--method ldes names the discrete-event engine, and no engine was found beside "
                "this binary (common/ldes or common/ldes.jar, or $LINE_LDES_DIR); -s ssa is the "
                "simulator this port builds in process");
        return "ldes";
    }
    if (fam == "jmt") {
        // IT IS WRAPPED NOW. This refused the token outright until the JMT
        // client landed, and stayed behind: `-s jmt` drives jsim and jmva
        // through `solver_jmt_run_analyzer`, so refusing `--method jmt` denied
        // under `-s auto` what the very same binary answers under `-s jmt`.
        // Availability is left to the wrapper, which names what is missing (a
        // JVM, common/JMT.jar or a REST endpoint) rather than guessing here.
        return "jmt";
    }
    if (fam == "qns") return "qns";
    if (fam == "lqns")
        throw line::UnsupportedError(
            "--method lqns names the external layered binary, which solves a LayeredNetwork and "
            "not a Network; reach it as -i lqnx -s lqns");
    if (fam == "ln")
        throw line::UnsupportedError(
            "--method ln names the layered solver, which takes a LayeredNetwork: pass the model "
            "as -i lqnx -s ln rather than as a Network");
    throw line::InputError("SolverAUTO: no engine stands behind method family '" + fam + "'");
}

/**
 * Is this model.json an Environment envelope rather than a Network?
 *
 * `-s auto` has to know before it parses: the two readers take different
 * documents, and the reference's chooser has an Environment arm that is
 * unreachable if every auto run is assumed to hold a Network. A malformed
 * document answers false, so the real reader reports the parse error.
 */
bool model_is_environment(const std::string& file) {
    try {
        line::io::detail::json root;
        if (file.empty()) {
            std::istringstream in(stdin_model_text());
            in >> root;
        } else {
            std::ifstream in(file.c_str());
            if (!in) return false;
            in >> root;
        }
        // The envelope may wrap the model, exactly as build_environment_from_json
        // unwraps it; `type` is what that reader keys on, so this reads the same
        // field rather than a second convention of its own.
        if (!root.is_object()) return false;
        const line::io::detail::json& model = root.contains("model") ? root.at("model") : root;
        return model.is_object() && model.value("type", std::string()) == "Environment";
    } catch (...) {
        return false;
    }
}

/**
 * `-s auto`: which engine answers, by `chooseSolver.m` through solver_auto.h.
 *
 * THREE THINGS DECIDE, in the reference's own order. `--method` is resolved
 * first, because a method FAMILY ('nc', 'nc.comom') names the engine outright
 * and bypasses every ranking, while a selection INTENT ('exact', 'sim', 'fast',
 * 'accurate') picks which ranking runs. Then the model class: an Environment
 * envelope takes the Environment arm, a Network the Network one. Then the
 * GETTER the caller asked for, since the reference keys the ranking on the
 * metric family and not on the model alone.
 *
 * THE ORDER IS A RETRY LIST, not a single name. `delegate.m` tries the chosen
 * solver and then every feasible candidate, so a refusal moves to the next
 * engine instead of ending the run; the caller sees which one answered.
 *
 * The choice reads structure only -- traits, feature sets, product form,
 * populations -- so it is made in the arithmetic the run will use and costs one
 * extra parse of the model, nothing more.
 */
template <class T>
AutoPlan choose_auto_plan(const std::string& file, const std::string& analysis,
                          const std::string& method_token) {
    const line::autosolver::AutoToken tok = line::autosolver::auto_resolve_token(method_token);
    AutoPlan plan;
    if (!tok.is_intent) {
        plan.order.push_back(auto_cli_token_of_family(tok.family));
        if (tok.submethod != "default") plan.method = tok.submethod;
        return plan;
    }

    const std::string getter = auto_getter_of_analysis(analysis);
    if (model_is_environment(file)) {
        const line::autosolver::AutoEnvChoice ec =
            line::autosolver::auto_choose_env_solver(getter, tok.mode);
        plan.order.push_back("env");
        for (std::size_t i = 0; i < ec.skipped.size(); ++i)
            plan.note += std::string(i ? ", " : " (the ranking preferred ") +
                         line::autosolver::auto_env_name(ec.skipped[i]);
        if (!ec.skipped.empty()) plan.note += ", which this port does not build)";
        return plan;
    }

    line::qn::Network<T> net = read_model<T>(file);
    const line::autosolver::AutoChoice c =
        line::autosolver::auto_choose_solver_mode(net.get_struct(), getter, tok.mode);
    plan.method = c.method;
    const std::vector<line::autosolver::AutoSolver> proposed =
        line::autosolver::auto_proposed_solvers(net.get_struct(), getter, tok.mode);
    for (std::size_t i = 0; i < proposed.size(); ++i)
        plan.order.push_back(line::autosolver::auto_solver_name(proposed[i]));
    for (std::size_t i = 0; i < c.skipped.size(); ++i)
        plan.note += std::string(i ? ", " : " (the ranking preferred ") +
                     line::autosolver::auto_solver_name(c.skipped[i]);
    if (!c.skipped.empty()) plan.note += ", which this port does not build)";
    return plan;
}

AutoPlan choose_auto_plan_dispatch(const std::string& arith, const std::string& file,
                                   const std::string& analysis, const std::string& method_token) {
    if (arith == "exact") return choose_auto_plan<line::Rational>(file, analysis, method_token);
    if (arith == "real:16") return choose_auto_plan<line::Real<16> >(file, analysis, method_token);
    if (arith == "real" || arith == "real:32")
        return choose_auto_plan<line::Real<32> >(file, analysis, method_token);
    if (arith == "real:64") return choose_auto_plan<line::Real<64> >(file, analysis, method_token);
    if (arith == "real:128")
        return choose_auto_plan<line::Real<128> >(file, analysis, method_token);
    if (arith == "real:256")
        return choose_auto_plan<line::Real<256> >(file, analysis, method_token);
    return choose_auto_plan<double>(file, analysis, method_token);
}

/**
 * The numeric clean-up SolverLN.getAvgTable applies before printing.
 *
 * It is reproduced here because the reference's reported table IS the
 * sanitized one, so a row-by-row comparison against it has to compare like
 * with like. Two rules: snap a value to one decimal place when it is already
 * within CoarseTol of it relatively, and snap anything at or below FineTol to
 * zero. The second rule is what turns the residual queue length of a chain of
 * Immediate classes -- of order 1e-8 by construction, since Immediate has rate
 * 1e8 -- into the exact zero the reference prints.
 */
double ln_sanitize(double x) {
    const double r = std::round(x * 10.0);
    if (std::fabs(x * 10.0 - r) < line::lang::GlobalConstants::CoarseTol * x * 10.0) x = r / 10.0;
    if (x <= line::lang::GlobalConstants::FineTol) x = 0.0;
    return x;
}

const char* ln_element_kind(const line::lqn::LqnStruct<double>& l, std::size_t i) {
    switch (l.type[i]) {
        case line::lang::LqnElement::HOST: return "Processor";
        case line::lang::LqnElement::TASK: return l.isref[i] ? "RefTask" : "Task";
        case line::lang::LqnElement::ENTRY: return "Entry";
        default: return "Activity";
    }
}

/**
 * Print every layer's stations, classes and routing, in a form a MATLAB dump of
 * `solver.ensemble{k}` can be diffed against line by line.
 *
 * A layered result that is close but not equal across codebases is almost never
 * a difference in the MVA call; it is a layer that was BUILT differently -- a
 * class that is present in one and not the other, a population, a routing
 * probability. Comparing the final AvgTable cannot tell those apart, so the
 * structure has to be observable directly.
 */
template <class T>
void ln_dump_layers(const line::ln::SolverLN<T>& solver) {
    using namespace line;
    const std::vector<qn::Layer<T> >& ens = solver.layers();
    for (std::size_t k = 0; k < ens.size(); ++k) {
        const qn::Layer<T>& L = ens[k];
        std::printf("LAYER %zu %s nstations=%zu nclasses=%zu nchains=%zu\n", k + 1, L.name.c_str(),
                    L.stations.size(), L.classes.size(), L.nchains);
        for (std::size_t i = 0; i < L.stations.size(); ++i)
            std::printf("  STATION %zu %s sched=%s nservers=%g\n", i + 1, L.stations[i].name.c_str(),
                        lang::sched_to_text(L.stations[i].sched), L.stations[i].nservers);
        for (std::size_t r = 0; r < L.classes.size(); ++r)
            std::printf("  CLASS %zu %s pop=%.17g refstat=%zu completes=%d\n", r + 1,
                        L.classes[r].name.c_str(), L.classes[r].population, L.classes[r].refstat,
                        int(L.classes[r].completes));
        for (std::size_t i = 0; i < L.stations.size(); ++i)
            for (std::size_t r = 0; r < L.classes.size(); ++r) {
                if (L.disabled.empty() || L.disabled[i][r]) continue;
                std::printf("  RATE %s %s %.17g scv=%.17g\n", L.stations[i].name.c_str(),
                            L.classes[r].name.c_str(), num_traits<T>::to_double(L.rates(i, r)),
                            num_traits<T>::to_double(L.scv(i, r)));
            }
        for (const auto& kv : L.P) {
            const Matrix<T>& B = kv.second;
            for (std::size_t i = 0; i < B.rows(); ++i)
                for (std::size_t j = 0; j < B.cols(); ++j) {
                    const double p = num_traits<T>::to_double(B(i, j));
                    if (p == 0.0) continue;
                    std::printf("  ROUTE %s->%s %s->%s %.17g\n",
                                L.classes[kv.first.first - 1].name.c_str(),
                                L.classes[kv.first.second - 1].name.c_str(),
                                L.nodes[i].name.c_str(), L.nodes[j].name.c_str(), p);
                }
        }
    }
}

/**
 * Solve a .lqnx layered queueing network with SolverLN and print its AvgTable.
 *
 * The output is one row per LQN element, with its queue length, utilization,
 * response time, residence time and throughput, in the same element order as
 * MATLAB's getAvgTable, so a row-by-row numeric comparison against the
 * reference is a plain diff.
 *
 * `--repeat` re-runs the whole solve K times and reports the best wall-clock
 * time, which is what the arithmetic-backend benchmark reads.
 */
/** The LnOptions a set of CLI knobs describes; arithmetic-independent. */
inline line::ln::LnOptions ln_options_from(const Knobs& k) {
    using namespace line;
    // The reference's own LnOptions defaults stand unless the caller overrode
    // them; the CLI does not restate them, so an untouched knob keeps whatever
    // SolverLN itself considers default.
    ln::LnOptions opt;
    if (k.iter_max >= 0) opt.iter_max = k.iter_max;
    if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
    if (k.no_interlocking) opt.interlocking = false;
    if (!k.layer_solver.empty()) opt.layer_solver = k.layer_solver;
    if (!k.method.empty()) opt.method = k.method;
    if (k.samples) opt.layer_ssa.samples = k.samples;
    if (k.seed) opt.layer_ssa.seed = k.seed;
    if (k.t1 >= 0.0) opt.timespan_end = k.t1;
    if (k.tran_points) opt.tran_points = k.tran_points;
    if (!k.ln_transient.empty()) opt.ln_transient = k.ln_transient;
    if (!k.ln_transient_channels.empty()) opt.ln_transient_channels = k.ln_transient_channels;
    return opt;
}

/** The engine name the banner reports, which is never a hardcoded MVA. */
inline const char* ln_layer_engine_name(const std::string& layer_solver) {
    if (layer_solver == "fluid") return "Fluid";
    if (layer_solver == "nc") return "NC";
    if (layer_solver == "ssa") return "SSA";
    return "MVA";
}

/**
 * `-a tran`: the layered transient, one block of time series per layer.
 *
 * The blocks are printed rather than the block-diagonal cell array the
 * reference assembles, because the off-diagonal blocks of that array are empty
 * by construction and each layer keeps its own grid.
 */
template <class T>
int run_ln_tran(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<T> model = io::read_layered_model<T>(file);
    ln::LnOptions opt = ln_options_from(k);
    ln::SolverLN<T> solver(model, opt);
    const ln::LnTranSolution tr = solver.get_tran_avg();
    const std::vector<qn::Layer<T>>& layers = solver.layers();

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = num_traits<T>::name();
        j["mode"] = tr.mode;
        j["iterations"] = tr.iterations;
        j["gap"] = ln_sanitize(tr.gap);
        reg::Json ls = reg::Json::array();
        for (std::size_t e = 0; e < tr.layers.size(); ++e) {
            reg::Json le = reg::Json::object();
            le["layer"] = layers[e].name;
            reg::Json tt = reg::Json::array();
            for (double t : tr.layers[e].t) tt.push_back(ln_sanitize(t));
            le["t"] = tt;
            reg::Json series = reg::Json::array();
            for (std::size_t i = 0; i < tr.layers[e].QN.size(); ++i)
                for (std::size_t r = 0; r < tr.layers[e].QN[i].size(); ++r) {
                    reg::Json s = reg::Json::object();
                    s["station"] = layers[e].stations[i].name;
                    s["class"] = layers[e].classes[r].name;
                    auto arr = [&](const std::vector<double>& v) {
                        reg::Json a = reg::Json::array();
                        for (double x : v) a.push_back(ln_sanitize(x));
                        return a;
                    };
                    s["QLen"] = arr(tr.layers[e].QN[i][r]);
                    s["Util"] = arr(tr.layers[e].UN[i][r]);
                    s["Tput"] = arr(tr.layers[e].TN[i][r]);
                    series.push_back(s);
                }
            le["series"] = series;
            ls.push_back(le);
        }
        j["layers"] = ls;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLN(Solver%s) getTranAvg arith=%s mode=%s layers=%zu iterations=%ld gap=%.3e\n",
                ln_layer_engine_name(opt.layer_solver), num_traits<T>::name(), tr.mode.c_str(),
                tr.layers.size(), tr.iterations, tr.gap);
    for (std::size_t e = 0; e < tr.layers.size(); ++e) {
        const std::vector<double>& t = tr.layers[e].t;
        if (t.empty()) continue;
        std::printf("\nLayer %s  (%zu points on [%.6g, %.6g])\n", layers[e].name.c_str(), t.size(),
                    t.front(), t.back());
        std::printf("%-30s %-24s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen(0)",
                    "QLen(end)", "Util(end)", "Tput(end)");
        for (std::size_t i = 0; i < tr.layers[e].QN.size(); ++i)
            for (std::size_t r = 0; r < tr.layers[e].QN[i].size(); ++r) {
                const std::vector<double>& q = tr.layers[e].QN[i][r];
                if (q.empty()) continue;
                std::printf("%-30s %-24s %12.6g %12.6g %12.6g %12.6g\n",
                            layers[e].stations[i].name.c_str(), layers[e].classes[r].name.c_str(),
                            ln_sanitize(q.front()), ln_sanitize(q.back()),
                            ln_sanitize(tr.layers[e].UN[i][r].back()),
                            ln_sanitize(tr.layers[e].TN[i][r].back()));
            }
    }
    return 0;
}

/** `-a sens`: the layer sensitivity tables under a leading Layer column. */
template <class T>
int run_ln_sens(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<T> model = io::read_layered_model<T>(file);
    ln::LnOptions opt = ln_options_from(k);
    ln::SolverLN<T> solver(model, opt);
    sens::SensOptions so;
    if (!k.sens_method.empty()) so.method = k.sens_method;
    if (!k.sens_scheme.empty()) so.scheme = k.sens_scheme;
    if (k.sens_step > 0.0) so.step = k.sens_step;
    const ln::LnSensTable<T> tbl = solver.get_sensitivity_table(so);

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = num_traits<T>::name();
        j["method"] = tbl.method;
        reg::Json rows = reg::Json::array();
        for (const auto& r : tbl.rows) {
            reg::Json o = reg::Json::object();
            o["Layer"] = r.layer;
            o["Station"] = r.station;
            o["JobClass"] = r.jobclass;
            o["dTput_dRate"] = sens_sanitize_signed(num_traits<T>::to_double(r.dTput));
            o["dRespT_dRate"] = sens_sanitize_signed(num_traits<T>::to_double(r.dRespT));
            o["dQLen_dRate"] = sens_sanitize_signed(num_traits<T>::to_double(r.dQLen));
            o["dUtil_dRate"] = sens_sanitize_signed(num_traits<T>::to_double(r.dUtil));
            rows.push_back(o);
        }
        j["rows"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLN(Solver%s) getSensitivityTable arith=%s method=%s rows=%zu\n",
                ln_layer_engine_name(opt.layer_solver), num_traits<T>::name(), tbl.method.c_str(),
                tbl.rows.size());
    std::printf("%-28s %-28s %-20s %14s %14s %14s %14s\n", "Layer", "Station", "JobClass",
                "dTput_dRate", "dRespT_dRate", "dQLen_dRate", "dUtil_dRate");
    for (const auto& r : tbl.rows)
        std::printf("%-28s %-28s %-20s %14.6g %14.6g %14.6g %14.6g\n", r.layer.c_str(),
                    r.station.c_str(), r.jobclass.c_str(),
                    sens_sanitize_signed(num_traits<T>::to_double(r.dTput)),
                    sens_sanitize_signed(num_traits<T>::to_double(r.dRespT)),
                    sens_sanitize_signed(num_traits<T>::to_double(r.dQLen)),
                    sens_sanitize_signed(num_traits<T>::to_double(r.dUtil)));
    return 0;
}

/** `-a cdf`: the per-entry response-time distribution of the moment3 method. */
template <class T>
int run_ln_cdf(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<T> model = io::read_layered_model<T>(file);
    ln::LnOptions opt = ln_options_from(k);
    // getCdfRespT.m runs the ensemble under moment3 whatever the caller asked
    // for, and restores the method afterwards: the mean-based update forms no
    // distribution at all, so there is nothing else to report.
    opt.method = "moment3";
    ln::SolverLN<T> solver(model, opt);
    const std::vector<ln::LnCdf> cdf = solver.get_cdf_respt();
    const lqn::LqnStruct<double> names = io::read_layered_model<double>(file);

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = num_traits<T>::name();
        reg::Json rows = reg::Json::array();
        for (std::size_t e = 1; e <= model.nentries && e < cdf.size(); ++e) {
            reg::Json o = reg::Json::object();
            o["entry"] = names.names[model.eshift + e];
            reg::Json tt = reg::Json::array(), ff = reg::Json::array();
            for (std::size_t p = 0; p < cdf[e].t.size(); ++p) {
                tt.push_back(ln_sanitize(cdf[e].t[p]));
                ff.push_back(ln_sanitize(cdf[e].cdf[p]));
            }
            o["t"] = tt;
            o["F"] = ff;
            rows.push_back(o);
        }
        j["entries"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLN(Solver%s) getCdfRespT arith=%s method=moment3 entries=%zu\n",
                ln_layer_engine_name(opt.layer_solver), num_traits<T>::name(), model.nentries);
    for (std::size_t e = 1; e <= model.nentries && e < cdf.size(); ++e) {
        if (cdf[e].t.empty()) {
            std::printf("%-40s (no distribution: the entry has no fitted term)\n",
                        names.names[model.eshift + e].c_str());
            continue;
        }
        // Quartiles, which is what a CDF is read for; the full grid goes to JSON.
        auto quantile = [&](double p) {
            for (std::size_t i = 0; i < cdf[e].cdf.size(); ++i)
                if (cdf[e].cdf[i] >= p) return cdf[e].t[i];
            return cdf[e].t.back();
        };
        std::printf("%-40s p25=%12.6g p50=%12.6g p75=%12.6g p95=%12.6g points=%zu\n",
                    names.names[model.eshift + e].c_str(), ln_sanitize(quantile(0.25)),
                    ln_sanitize(quantile(0.50)), ln_sanitize(quantile(0.75)),
                    ln_sanitize(quantile(0.95)), cdf[e].t.size());
    }
    return 0;
}

template <class T>
int run_ln(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<T> model = io::read_layered_model<T>(file);

    ln::LnOptions opt = ln_options_from(k);
    const int repeat = k.repeat > 0 ? k.repeat : 1;

    double best = 1e300;
    ln::LnSolution<T> sol;
    std::size_t nlayers = 0;
    for (int rep = 0; rep < repeat; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        ln::SolverLN<T> solver(model, opt);
        sol = solver.get_ensemble_avg();
        nlayers = solver.nlayers();
        if (output == "layers") {
            ln_dump_layers(solver);
            return 0;
        }
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
    }

    // element names and kinds are arithmetic-independent, so read them once
    const lqn::LqnStruct<double> names = io::read_layered_model<double>(file);

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = num_traits<T>::name();
        j["layers"] = nlayers;
        j["iterations"] = sol.iterations;
        j["converged"] = sol.converged;
        j["seconds"] = best;
        reg::Json rows = reg::Json::array();
        for (std::size_t i = 1; i <= model.nidx; ++i) {
            reg::Json r = reg::Json::object();
            r["node"] = names.names[i];
            r["type"] = ln_element_kind(names, i);
            auto put = [&](const char* key, const std::vector<T>& v, const std::vector<bool>& d) {
                if (d[i]) r[key] = ln_sanitize(num_traits<T>::to_double(v[i]));
                else if (sol.is_bound) r[key] = 0.0;  // see the note in the table below
                else r[key] = nullptr;
            };
            put("QLen", sol.QN, sol.defined_Q);
            put("Util", sol.UN, sol.defined_U);
            put("RespT", sol.RN, sol.defined_R);
            put("ResidT", sol.WN, sol.defined_W);
            put("Tput", sol.TN, sol.defined_T);
            rows.push_back(r);
        }
        j["rows"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    // The banner names the LAYER solver actually used, not a hardcoded MVA: the
    // two converge to different fixed points, so a reader (or a parity row) that
    // cannot tell them apart is reading numbers it cannot attribute.
    std::printf(
        "SolverLN(Solver%s) arith=%s type=%s layers=%zu iterations=%d converged=%d time=%.4fs\n",
        ln_layer_engine_name(opt.layer_solver), num_traits<T>::name(),
        line::util::method_type("LN", opt.method).c_str(), nlayers, sol.iterations,
        int(sol.converged), best);
    std::printf("%-62s %-10s %12s %12s %12s %12s %12s\n", "Node", "NodeType", "QLen", "Util",
                "RespT", "ResidT", "Tput");
    for (std::size_t i = 1; i <= model.nidx; ++i) {
        // NaN IS THE UNDEFINED MARKER EVERYWHERE EXCEPT UNDER A BOUND. MATLAB's
        // getAvgTable prints NaN for a measure the element does not have (a
        // processor has no queue length), and this reproduces that. A BOUND is
        // the one case where the reference prints 0 instead: `mwba.*` defines
        // throughput and processor utilization only, and both MATLAB and the JAR
        // report the rest as zero rather than as absent (the JAR maps the NaN
        // explicitly, SolverLN.java:3294-3299). The `defined_*` flags still say
        // undefined to any caller of the API; only the printed table follows the
        // reference, so that a numeric parity row compares like with like.
        //
        // This table is for a human; a cross-codebase comparison must read the
        // `-o json` above, which carries the raw double, and quantize it itself.
        auto fmt = [&](const std::vector<T>& v, const std::vector<bool>& d, char* buf) {
            if (!d[i] && sol.is_bound) std::snprintf(buf, 24, "%12.6g", 0.0);
            else if (!d[i]) std::snprintf(buf, 24, "%12s", "NaN");
            else std::snprintf(buf, 24, "%12.6g", ln_sanitize(num_traits<T>::to_double(v[i])));
        };
        char q[24], u[24], rr[24], w[24], t[24];
        fmt(sol.QN, sol.defined_Q, q);
        fmt(sol.UN, sol.defined_U, u);
        fmt(sol.RN, sol.defined_R, rr);
        fmt(sol.WN, sol.defined_W, w);
        fmt(sol.TN, sol.defined_T, t);
        std::printf("%-62s %-10s %s %s %s %s %s\n", names.names[i].c_str(),
                    ln_element_kind(names, i), q, u, rr, w, t);
    }
    return 0;
}

/**
 * `-s ldes`: the layered model simulated directly, by the IN-PROCESS engine.
 *
 * IT IS NOT THE SUBPROCESS `-s ldes` OF THE FLAT PATH. That arm hands a
 * `model.json` to `common/ldes`, and the engine behind that wire refuses a
 * layered document outright ("LDES currently supports Network models only"), so
 * there is nothing to forward. This arm calls `ldes_ln_engine_solve` in process
 * -- the C++ twin of `Solver_ssj_ln.java` -- which simulates entries,
 * activities, task threads and synchronous calls directly instead of
 * decomposing the model into layers. So it is not a noisier route to `-s ln`:
 * SolverLN's decomposition is an APPROXIMATION and this is a sample path of the
 * model itself, which is what makes it the reference the layered solvers are
 * checked against.
 *
 * THE NaN MASK OF THE LAYERED TABLE IS PART OF THE ANSWER, and it belongs to
 * the table rather than to whichever solver filled it: a processor has no queue
 * length, a task no response time, an entry no residence, and `-s ln` and
 * `-s lqns` both print NaN there. This engine MEASURES more than that -- a
 * processor's completion rate is sitting in `LnResult::TLN` -- and printing it
 * under a column the other two arms leave empty would make one column mean
 * different things depending on who filled it. That is the divergence the JAR
 * removed from `getLNAvgTable` on 2026-08-21, and masking here rather than in
 * the engine keeps it removed on this side too. Nothing is discarded: the
 * unmasked measurements stay on `LnResult` for a programmatic caller, and only
 * the shared table is masked, so that it can be diffed row by row.
 *
 * THE NUMBERS ARE NOT `ln_sanitize`d, for the reason `run_lqns` states below:
 * that helper snaps to a tenth and floors at FineTol, which is right for a
 * fixed point this port computed and wrong for a measurement it made. A
 * simulated utilization of 1e-9 is a rare event that was observed, not a
 * negligible residue of an iteration, and the JAR's own layered LDES table
 * reports the raw estimate too -- so snapping here would put the two out of
 * step on exactly the models a parity row is read on.
 */
/**
 * `-i lqnx -s ldes -a cdf`: the SIMULATED response time distribution per entry.
 *
 * The measured counterpart of `run_ln_cdf`, which fits an APH to three moments of
 * a fluid passage time and convolves. Here the engine timed every invocation from
 * the instant the request reached the entry to its reply -- the interval `RLN`
 * averages -- so this law's mean reproduces that row and its tail is observed
 * rather than extrapolated.
 */
int run_ln_ldes_cdf(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<double> model = io::read_layered_model<double>(file);
    ldes::LdesOptions o;
    if (k.samples) o.samples = k.samples;
    if (k.seed) o.seed = static_cast<long>(k.seed);

    const ldes::engine::LnResult r = ldes::ldes_ln_engine_solve(model, o);

    // The per-entry ecdf, through the same API function the other codebases
    // call getCdfRespTLN.
    const std::vector<ldes::engine::LnEntryCdf> cdf =
        ldes::engine::ldes_ln_cdf_respt(r, model.nentries);
    std::vector<std::vector<double> > tt(model.nentries), ff(model.nentries);
    for (std::size_t e = 0; e < model.nentries; ++e) {
        tt[e] = cdf[e].t;
        ff[e] = cdf[e].F;
    }

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["engine"] = "native-ln";
        reg::Json rows = reg::Json::array();
        for (std::size_t e = 0; e < model.nentries; ++e) {
            reg::Json ob = reg::Json::object();
            ob["entry"] = model.names[model.eshift + e + 1];
            reg::Json at = reg::Json::array(), af = reg::Json::array();
            for (std::size_t i = 0; i < tt[e].size(); ++i) {
                at.push_back(tt[e][i]);
                af.push_back(ff[e][i]);
            }
            ob["t"] = at;
            ob["F"] = af;
            ob["observations"] = static_cast<double>(
                e < r.entry_resp_samples.size() ? r.entry_resp_samples[e].size() : 0);
            rows.push_back(ob);
        }
        j["entries"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLDES(native LN engine) getCdfRespT entries=%zu\n", model.nentries);
    for (std::size_t e = 0; e < model.nentries; ++e) {
        if (tt[e].empty()) {
            std::printf("%-40s (no observation)\n", model.names[model.eshift + e + 1].c_str());
            continue;
        }
        // Quartiles plus p95, which is what a measured law is read for.
        struct Q {
            const std::vector<double>& t;
            const std::vector<double>& f;
            double operator()(double p) const {
                for (std::size_t i = 0; i < f.size(); ++i)
                    if (f[i] >= p) return t[i];
                return t.back();
            }
        } q = {tt[e], ff[e]};
        std::printf("%-40s p25=%12.6g p50=%12.6g p75=%12.6g p95=%12.6g n=%zu\n",
                    model.names[model.eshift + e + 1].c_str(), q(0.25), q(0.50), q(0.75),
                    q(0.95), r.entry_resp_samples[e].size());
    }
    return 0;
}

int run_ln_ldes(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnStruct<double> model = io::read_layered_model<double>(file);

    // The engine reads THREE settings and no more (`o.samples`, `o.events`,
    // `o.seed`); the dispatcher refuses the rest of the --ldes-* family rather
    // than letting this function drop them silently.
    ldes::LdesOptions o;
    if (k.samples) o.samples = k.samples;
    if (k.seed) o.seed = static_cast<long>(k.seed);

    // --repeat is a TIMING loop and stays honest here only because the stream is
    // seeded: every replication of one command walks the same sample path, so
    // the table is that run's table and `time=` is the only thing that varies.
    // This is why -s lqns refuses the flag and this arm does not -- lqsim seeds
    // itself, so re-running it would report different numbers under one banner.
    const int repeat = k.repeat > 0 ? k.repeat : 1;
    double best = 1e300;
    ldes::engine::LnResult r;
    for (int rep = 0; rep < repeat; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        r = ldes::ldes_ln_engine_solve(model, o);
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
    }

    // The mask lives beside the data it masks (`ldes::engine::ln_defined`), so
    // that the doctest suite can pin it against `LnSolution::defined_*` without
    // reaching into this file. Spelling it here instead would put the rule that
    // decides what the table MEANS inside a printer.
    typedef ldes::engine::LnColumn Col;
    const auto def = [&](std::size_t i, Col c) {
        return ldes::engine::ln_defined(model, r, i, c);
    };

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = "double";
        j["solver"] = "ldes";
        // WHICH ENGINE ANSWERED, on the same grounds the flat arm's `engine=` is
        // not decoration: `-s ldes` names two different simulators depending on
        // whether the model is layered, and a number quoted from one must not be
        // read as the other's.
        j["engine"] = "native-ln";
        j["samples"] = o.samples;
        j["seed"] = o.seed;
        j["simulatedTime"] = r.simulated_time;
        j["completions"] = r.completions;
        j["seconds"] = best;
        reg::Json rows = reg::Json::array();
        for (std::size_t i = 1; i <= model.nidx; ++i) {
            reg::Json row = reg::Json::object();
            row["node"] = model.names[i];
            row["type"] = ln_element_kind(model, i);
            auto put = [&](const char* key, double v, bool defined) {
                if (defined) row[key] = v;
                else row[key] = nullptr;
            };
            put("QLen", r.QLN(i, 0), def(i, Col::QLen));
            put("Util", r.ULN(i, 0), def(i, Col::Util));
            put("RespT", r.RLN(i, 0), def(i, Col::RespT));
            put("ResidT", r.WLN(i, 0), def(i, Col::ResidT));
            put("Tput", r.TLN(i, 0), def(i, Col::Tput));
            rows.push_back(row);
        }
        j["rows"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLDES(native LN engine) arith=double type=%s samples=%zu seed=%ld "
                "simtime=%.6g completions=%lld time=%.4fs\n",
                line::util::method_type("LDES", o.method).c_str(), o.samples, o.seed,
                r.simulated_time, r.completions, best);
    std::printf("%-62s %-10s %12s %12s %12s %12s %12s\n", "Node", "NodeType", "QLen", "Util",
                "RespT", "ResidT", "Tput");
    for (std::size_t i = 1; i <= model.nidx; ++i) {
        // Six digits, as on the two tables around it, so the three arms diff.
        auto fmt = [&](double v, bool defined, char* buf) {
            if (!defined) std::snprintf(buf, 24, "%12s", "NaN");
            else std::snprintf(buf, 24, "%12.6g", v);
        };
        char q[24], u[24], rr[24], w[24], t[24];
        fmt(r.QLN(i, 0), def(i, Col::QLen), q);
        fmt(r.ULN(i, 0), def(i, Col::Util), u);
        fmt(r.RLN(i, 0), def(i, Col::RespT), rr);
        fmt(r.WLN(i, 0), def(i, Col::ResidT), w);
        fmt(r.TLN(i, 0), def(i, Col::Tput), t);
        std::printf("%-62s %-10s %s %s %s %s %s\n", model.names[i].c_str(),
                    ln_element_kind(model, i), q, u, rr, w, t);
    }
    return 0;
}

/**
 * `-s lqns`: the same layered model, solved by the external binary.
 *
 * The table has the SAME columns as run_ln's so that the two can be diffed row
 * by row, but the numbers are not sanitized the same way: `ln_sanitize` also
 * snaps anything at or below FineTol to zero, which is right for a fixed point
 * this port computed and wrong for a measurement it did not -- lqns reports no
 * residence time and no arrival rate at all, and a zero there would read as a
 * computed zero. Those two columns print NaN, and only the snap-to-tenth of the
 * reference's getAvgTable is applied.
 */
template <class T>
int run_lqns(const std::string& file, const std::string& output, const Knobs& k) {
    using namespace line;
    const lqn::LqnModel<T> model = lqn::read_lqnx_model<T>(file);

    lqns::LqnsOptions opt;
    if (!k.method.empty()) opt.method = k.method;
    if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
    if (k.samples) opt.samples = static_cast<double>(k.samples);
    opt.verbose = k.verbose;
    opt.keep = k.keep;
    opt.remote = k.remote;
    if (!k.remote_url.empty()) opt.remote_url = k.remote_url;
    opt.timeout_seconds = k.timeout_seconds;

    lqns::SolverLQNS<T> solver(model, opt);
    const lqns::LqnsSolution<T> sol = solver.get_ensemble_avg();
    const lqn::LqnStruct<T>& sn = solver.get_struct();
    const lqn::LqnStruct<double> names = io::read_layered_model<double>(file);

    if (output == "json") {
        reg::Json j = reg::Json::object();
        j["model"] = file;
        j["arith"] = num_traits<T>::name();
        j["solver"] = lqns::SolverLQNS<T>::is_stochastic_method(opt.method) ? "lqsim" : "lqns";
        j["method"] = opt.method;
        j["iterations"] = sol.iterations;
        j["seconds"] = solver.runtime();
        reg::Json rows = reg::Json::array();
        for (std::size_t i = 1; i <= sn.nidx; ++i) {
            reg::Json r = reg::Json::object();
            r["node"] = names.names[i];
            r["type"] = ln_element_kind(names, i);
            auto put = [&](const char* key, const std::vector<T>& v, const std::vector<bool>& d) {
                if (d[i]) r[key] = lqns::detail::snap_to_tenth(num_traits<T>::to_double(v[i]));
                else r[key] = nullptr;
            };
            put("QLen", sol.QN, sol.defined_Q);
            put("Util", sol.UN, sol.defined_U);
            put("RespT", sol.RN, sol.defined_R);
            put("ResidT", sol.WN, sol.defined_W);
            put("Tput", sol.TN, sol.defined_T);
            rows.push_back(r);
        }
        j["rows"] = rows;
        std::printf("%s\n", j.dump(1).c_str());
        return 0;
    }

    std::printf("SolverLQNS(%s) arith=%s type=%s iterations=%d time=%.4fs\n",
                lqns::SolverLQNS<T>::version().c_str(), num_traits<T>::name(),
                line::util::method_type("LQNS", opt.method).c_str(), sol.iterations,
                solver.runtime());
    std::printf("%-62s %-10s %12s %12s %12s %12s %12s\n", "Node", "NodeType", "QLen", "Util",
                "RespT", "ResidT", "Tput");
    for (std::size_t i = 1; i <= sn.nidx; ++i) {
        // Six digits, as on the SolverLN table above and for the same reason.
        auto fmt = [&](const std::vector<T>& v, const std::vector<bool>& d, char* buf) {
            if (!d[i]) std::snprintf(buf, 24, "%12s", "NaN");
            else
                std::snprintf(buf, 24, "%12.6g",
                              lqns::detail::snap_to_tenth(num_traits<T>::to_double(v[i])));
        };
        char q[24], u[24], rr[24], w[24], t[24];
        fmt(sol.QN, sol.defined_Q, q);
        fmt(sol.UN, sol.defined_U, u);
        fmt(sol.RN, sol.defined_R, rr);
        fmt(sol.WN, sol.defined_W, w);
        fmt(sol.TN, sol.defined_T, t);
        std::printf("%-62s %-10s %s %s %s %s %s\n", names.names[i].c_str(),
                    ln_element_kind(names, i), q, u, rr, w, t);
    }
    return 0;
}

/**
 * `-i lqnx|xml`: the layered path, sibling to solve_model_dispatch.
 *
 * `-s auto` CONSULTS THE LAYERED ARM of the chooser, which is keyed on the
 * metric and on one trait of the model: a cache task promotes the NC layer
 * engine, because the cache layer is where NC beats MVA. LQNS leads two of
 * those rankings and IS wrapped, so `-s auto` selects it wherever the binary is
 * installed -- which makes the choice machine-dependent, exactly as
 * `chooseAvgSolverHeur.m` makes it, since LINE ships no LQNS binary. The banner
 * always names what answered.
 *
 * THE TOKENS NAME THE LAYER ENGINE, as the Java CLI's do: `ln.mva` runs the
 * layers under SolverMVA and `ln.comom` under SolverNC. The bare `ln` keeps the
 * MVA layers this port has always given it -- the Java CLI reads it as NC, but
 * changing it here would re-baseline every existing layered result under a token
 * whose meaning nothing in this tree states -- so `ln.comom` is the way to ask
 * for NC layers. `lqns` is not one of them: it does not solve LAYERS at all, it
 * hands the whole model to another program, so it takes the wrapper's own knobs
 * and none of SolverLN's.
 */
/** The solver method name as the console names it, e.g. "mva" -> "MVA". */
std::string upper_tag(const std::string& s) {
    std::string out = s;
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
    return out;
}

int solve_lqn_dispatch(const std::string& arith, const std::string& solver,
                       const std::string& analysis, const std::string& output,
                       const std::string& file, const Knobs& k) {
    if (file.empty())
        throw line::InputError(
            "a layered model is read from a file: pass -f <model.lqnx> or -f <model.json> (neither "
            "layered reader has a stdin form)");
    if (analysis != "avg" && analysis != "tran" && analysis != "sens" && analysis != "cdf")
        throw line::UnsupportedError(
            "the layered path ports -a avg (getAvgTable), -a tran (getTranAvg), -a sens "
            "(getSensitivityTable) and -a cdf (getCdfRespT); got '" + analysis + "'");
    const std::string s = solver.empty() ? "auto" : solver;
    if (s != "auto" && s != "ln" && s != "ln.mva" && s != "ln.comom" && s != "lqns" &&
        s != "ldes")
        throw line::UnsupportedError(
            "the layered path takes -s ln, ln.mva, ln.comom, ldes, lqns or auto (got '" + s +
            "'); a Network solver cannot be applied to a LayeredNetwork directly");

    // Solver console: the layered path has its own dispatcher and never
    // reaches solve_model_dispatch, so it opens the narrated run here. The
    // guard's destructor closes it, on an exception too.
    line::util::LineConsole::Run consoleRun(upper_tag(s), "", true);

    Knobs kk = k;
    // What will actually answer: the token, or what `-s auto` resolves to.
    std::string engine = s;
    // `-s auto` with no explicit layer engine consults chooseLayeredSolver. An
    // explicit --layer-solver is a choice the caller already made, so the
    // chooser does not overrule it.
    if (s == "auto" && kk.layer_solver.empty()) {
        std::string getter = "getAvgTable";
        if (analysis == "tran") getter = "getTranAvg";
        else if (analysis == "cdf") getter = "getCdfRespT";
        else if (analysis == "sens") getter = "getSensitivityTable";
        // `iscache` is the one trait the ranking reads, and it costs one parse
        // of the .lqnx: the same document the solve parses again below.
        const line::lqn::LqnStruct<double> probe = line::io::read_layered_model<double>(file);
        bool has_cache_task = false;
        for (std::size_t i = 0; i < probe.iscache.size(); ++i)
            if (probe.iscache[i]) has_cache_task = true;
        const line::autosolver::AutoLayeredChoice lc =
            line::autosolver::auto_choose_layered_solver(getter, has_cache_task);
        const std::string token = line::autosolver::auto_layered_name(lc.solver);
        if (token == "lqns") engine = "lqns";
        else if (token == "ln.comom") kk.layer_solver = "nc";
        else if (token == "ln.fluid") kk.layer_solver = "fluid";
        else kk.layer_solver = "mva";
        std::string note;
        for (std::size_t i = 0; i < lc.skipped.size(); ++i)
            note += std::string(i ? ", " : " (the ranking preferred ") +
                    line::autosolver::auto_layered_name(lc.skipped[i]);
        if (!lc.skipped.empty()) note += ", which is not available here)";
        std::printf("SolverAUTO selected %s%s\n", token.c_str(), note.c_str());
    }

    // ---- the wrapper, BEFORE the SolverLN knob ladder ---------------------
    // It is a wrapper, not a layer engine: --samples is the lqsim run length
    // rather than a simulated LAYER's, and --keep, --remote and --timeout
    // describe a child process no SolverLN run has. Asking the ladder below
    // about them would answer for the wrong solver.
    if (engine == "lqns") {
        if (analysis != "avg")
            throw line::UnsupportedError(
                "SolverLQNS reports the mean table its binary computes; it has no transient, no "
                "sensitivity and no response-time distribution here, so it takes -a avg (got '" +
                analysis + "')");
        if (!kk.layer_solver.empty())
            throw line::InputError(
                "--layer-solver names the engine SolverLN runs on each layer; -s lqns solves no "
                "layers, it hands the whole model to the lqns binary");
        if (output == "layers")
            throw line::UnsupportedError(
                "-o layers dumps the stations and routing SolverLN BUILT from the model; lqns "
                "builds its own submodels inside another process and this port never sees them");
        if (kk.iter_tol >= 0.0 || kk.iter_max > 0)
            throw line::UnsupportedError(
                "--iter_tol and --iter_max are SolverLN's layer-iteration knobs; lqns runs its own "
                "iteration and takes neither (its --iteration-limit is unreliable as of 6.2.27, "
                "which is why the reference stopped passing it)");
        if (kk.seed)
            throw line::UnsupportedError(
                "--seed sets the stream of a simulator this port drives; lqsim seeds itself and "
                "the wrapper passes no seed, exactly as the reference does not");
        if (kk.repeat > 1)
            throw line::UnsupportedError(
                "--repeat times a solve by re-running it; re-running lqsim would report a "
                "different answer under the same banner");
        if (kk.no_interlocking || !kk.ln_transient.empty() || !kk.ln_transient_channels.empty() ||
            !kk.sens_method.empty() || !kk.sens_scheme.empty() || kk.sens_step > 0.0)
            throw line::UnsupportedError(
                "--no-interlocking, --ln-transient*, and --sens-* are SolverLN options; -s lqns "
                "has none of them");
        if (arith != "double")
            throw line::UnsupportedError(
                "SolverLQNS reads a result file another program wrote in decimal double "
                "precision; there is no higher precision to carry, so rerun with --arith double "
                "(got '" + arith + "')");
        // The Network-only knobs the SolverLN ladder below refuses are refused
        // HERE TOO. This branch returns before that ladder runs, so a knob left
        // out of it is silently DROPPED rather than refused -- and which branch
        // a bare `.lqnx` path takes depends on whether an lqns binary is
        // installed, so the same command line would be refused on one machine
        // and quietly ignored on another.
        if (kk.has_cutoff())
            throw line::UnsupportedError(
                "--cutoff bounds the open population of a CTMC state space and applies to -s "
                "ctmc; the layered path enumerates no states");
        if (kk.node)
            throw line::UnsupportedError(
                "--node selects the stateful node a CTMC query is labelled by; the layered path "
                "reports every LQN element");
        return run_lqns<double>(file, output, kk);
    }
    if (k.keep || k.verbose || k.remote || !k.remote_url.empty() || k.timeout_seconds)
        throw line::UnsupportedError(
            "--keep, --verbose, --remote, --remote-url and --timeout describe the child process "
            "of an external solver and apply to -s lqns only");
    // ---- the native LN simulator, BEFORE the SolverLN knob ladder ----------
    // It solves no LAYERS, so the ladder below asks its questions of a
    // decomposition this arm never builds: --iter_tol and --iter_max bound a
    // fixed point it does not iterate, --layer-solver names an engine it does
    // not run, and --samples and --seed -- which the ladder refuses outright
    // unless a layer is simulated -- are precisely this arm's two settings.
    if (engine == "ldes") {
        if (analysis != "avg" && analysis != "cdf")
            throw line::UnsupportedError(
                "the native LN engine measures a sample path: it takes -a avg for the mean table "
                "and -a cdf for the per-entry response time distribution, and has no transient "
                "and no sensitivity here (got '" + analysis + "')");
        if (arith != "double")
            throw line::UnsupportedError(
                "the native LN engine accumulates its estimators in double, so there is no higher "
                "precision to carry; rerun with --arith double (got '" + arith + "')");
        if (output == "layers")
            throw line::UnsupportedError(
                "-o layers dumps the stations and routing SolverLN BUILT from the model; -s ldes "
                "simulates the layered semantics directly and builds no submodels");
        if (!kk.layer_solver.empty())
            throw line::InputError(
                "--layer-solver names the engine SolverLN runs on each layer; -s ldes solves no "
                "layers, it simulates entries, activities and calls directly");
        if (kk.iter_tol >= 0.0 || kk.iter_max > 0 || kk.no_interlocking)
            throw line::UnsupportedError(
                "--iter_tol, --iter_max and --no-interlocking are SolverLN's layer-iteration "
                "knobs; a simulated sample path converges by run length, which is --samples");
        if (!kk.ln_transient.empty() || !kk.ln_transient_channels.empty() ||
            !kk.sens_method.empty() || !kk.sens_scheme.empty() || kk.sens_step > 0.0)
            throw line::UnsupportedError(
                "--ln-transient* and --sens-* are SolverLN options; -s ldes has none of them");
        if (!kk.method.empty() && kk.method != "default")
            throw line::UnsupportedError(
                "--method on the layered path names the LN UPDATE (default, moment3, mwba.*); "
                "-s ldes performs no update, it simulates the model (got '" + kk.method + "')");
        // The --ldes-* family is the SUBPROCESS engine's settings. The native LN
        // engine reads samples, events and seed and nothing else, so a warmup
        // filter or a CI estimator passed here would be DROPPED rather than
        // honoured -- and a dropped `--ldes-tranfilter none` reads as a run with
        // no warmup removal that in fact removed one.
        if (!kk.ldes_tranfilter.empty() || kk.ldes_warmupfrac >= 0.0 ||
            !kk.ldes_cimethod.empty() || kk.ldes_cnvgon || kk.ldes_cnvgtol > 0.0 ||
            kk.ldes_slotted || kk.ldes_slotlength > 0.0 || kk.ldes_replications > 0 ||
            kk.ldes_numthreads > 0 || kk.ldes_maxtime > 0.0 || !kk.ldes_initsol.empty() ||
            !kk.ldes_rest_url.empty())
            throw line::UnsupportedError(
                "the --ldes-* flags configure the SUBPROCESS engine that answers -s ldes on a "
                "Network (warmup filter, CI estimator, slot lattice, replications, warm-start "
                "placement); the native LN engine behind -i lqnx -s ldes reads --samples and "
                "--seed only");
        // The Network-only knobs the SolverLN ladder refuses below are refused
        // HERE TOO: this branch returns before that ladder runs, so a knob left
        // out of it would be silently dropped rather than refused.
        if (kk.has_cutoff())
            throw line::UnsupportedError(
                "--cutoff bounds the open population of a CTMC state space and applies to -s "
                "ctmc; the layered path enumerates no states");
        if (kk.node)
            throw line::UnsupportedError(
                "--node selects the stateful node a CTMC query is labelled by; the layered path "
                "reports every LQN element");
        if (kk.t1 >= 0.0)
            throw line::UnsupportedError(
                "--tspan sets the horizon of a transient analysis; the native LN engine runs to "
                "a completion budget, which is --samples");
        if (kk.tol >= 0.0)
            throw line::UnsupportedError(
                "--tol is not an LDES option; the run length is set with --samples");
        if (analysis == "cdf") return run_ln_ldes_cdf(file, output, kk);
        return run_ln_ldes(file, output, kk);
    }
    // The solver method name and --layer-solver name the same choice, so they may not
    // disagree: silently letting one win would report the other in the banner.
    if (s == "ln.comom") {
        if (!kk.layer_solver.empty() && kk.layer_solver != "nc")
            throw line::InputError("-s ln.comom already selects NC layers, but --layer-solver says '" +
                                   kk.layer_solver + "'");
        kk.layer_solver = "nc";
    } else if (s == "ln.mva") {
        if (!kk.layer_solver.empty() && kk.layer_solver != "mva")
            throw line::InputError("-s ln.mva already selects MVA layers, but --layer-solver says '" +
                                   kk.layer_solver + "'");
        kk.layer_solver = "mva";
    }
    // --layer-solver is this port's own flag, the C++ spelling of the
    // reference's solver FACTORY: `LN(model, @(m) MVA(m))` against
    // `LN(model, @(m) Fluid(m))`. They converge to DIFFERENT fixed points.
    if (!kk.layer_solver.empty() && kk.layer_solver != "mva" && kk.layer_solver != "nc" &&
        kk.layer_solver != "fluid" && kk.layer_solver != "ssa")
        throw line::InputError("--layer-solver takes mva, nc, fluid or ssa (got '" +
                               kk.layer_solver + "')");
    // ---- knobs the layered path does not have are refused, not dropped ----
    if ((k.samples || k.seed) && kk.layer_solver != "ssa")
        throw line::UnsupportedError(
            "--samples and --seed set the run length and the stream of a SIMULATED layer; the "
            "layered path draws no random numbers unless --layer-solver ssa is in force");
    if (k.has_cutoff())
        throw line::UnsupportedError(
            "--cutoff bounds the open population of a CTMC state space and applies to -s ctmc; the "
            "layered path enumerates no states");
    if (k.t1 >= 0.0 && analysis != "tran")
        throw line::UnsupportedError(
            "--tspan sets the horizon of a transient analysis and applies to the layered path "
            "only with -a tran");
    if (analysis == "tran" && !(k.t1 >= 0.0))
        throw line::InputError(
            "-a tran integrates each layer's drift and needs a horizon: pass --tspan <t0>:<t1>");
    if (k.node)
        throw line::UnsupportedError(
            "--node selects the stateful node a CTMC query is labelled by; the layered path "
            "reports every LQN element");
    if (k.tol >= 0.0)
        throw line::UnsupportedError(
            "--tol is not a SolverLN option (LnOptions carries iter_tol and iter_max); "
            "use --iter_tol");
    // --method now names the LN UPDATE, which is a different question from the
    // layer engine: `moment3` reports a distribution the default never forms,
    // and the two bound requests report a bound instead of a fixed point.
    if (!k.method.empty() && k.method != "default" && k.method != "moment3" &&
        k.method != "mwba.upper" && k.method != "mwba.lower")
        throw line::UnsupportedError(
            "--method on the layered path takes default, moment3, mwba.upper or mwba.lower "
            "(got '" + k.method + "'); the per-layer engine is chosen with --layer-solver");
    if ((k.method == "mwba.upper" || k.method == "mwba.lower") && analysis != "avg")
        throw line::UnsupportedError(
            "--method mwba.* reports a throughput and utilization BOUND and solves no layer, so "
            "it has no transient, no sensitivity and no response-time law; use -a avg");
    if (!k.sens_method.empty() && analysis != "sens")
        throw line::UnsupportedError("--sens-method applies to -a sens");
    if (!k.ln_transient.empty() && analysis != "tran")
        throw line::UnsupportedError("--ln-transient applies to -a tran");

    if (analysis == "tran") {
        if (arith != "double")
            throw line::UnsupportedError(
                "the layered transient integrates each layer's drift with LSODA, which is double "
                "precision by construction; rerun with --arith double (got '" + arith + "')");
        return run_ln_tran<double>(file, output, kk);
    }
    if (analysis == "cdf") {
        if (arith != "double")
            throw line::UnsupportedError(
                "-a cdf fits an APH to a fluid passage time, integrated by LSODA in double "
                "precision; rerun with --arith double (got '" + arith + "')");
        return run_ln_cdf<double>(file, output, kk);
    }
    if (analysis == "sens") {
        if (arith == "double") return run_ln_sens<double>(file, output, kk);
        if (arith == "exact") return run_ln_sens<line::Rational>(file, output, kk);
        if (arith == "real:16") return run_ln_sens<line::Real<16> >(file, output, kk);
        if (arith == "real" || arith == "real:32")
            return run_ln_sens<line::Real<32> >(file, output, kk);
        if (arith == "real:64") return run_ln_sens<line::Real<64> >(file, output, kk);
        if (arith == "real:128") return run_ln_sens<line::Real<128> >(file, output, kk);
        if (arith == "real:256") return run_ln_sens<line::Real<256> >(file, output, kk);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }

    if (arith == "double") return run_ln<double>(file, output, kk);
    if (arith == "exact") return run_ln<line::Rational>(file, output, kk);
    // precision-ladder rationale (real:16 rung): see _kb/14-cpp-multiprecision.md
    if (arith == "real:16") return run_ln<line::Real<16> >(file, output, kk);
    if (arith == "real" || arith == "real:32") return run_ln<line::Real<32> >(file, output, kk);
    if (arith == "real:64") return run_ln<line::Real<64> >(file, output, kk);
    if (arith == "real:128") return run_ln<line::Real<128> >(file, output, kk);
    if (arith == "real:256") return run_ln<line::Real<256> >(file, output, kk);
    throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
}

/**
 * `-a node`: `@@NetworkSolver/getAvgNodeTable`, the means per NODE.
 *
 * A DIFFERENT INDEX SPACE FROM `-a avg`, not a relabelling of it. The AvgTable
 * is indexed by STATION, so a ClassSwitch, a Router, a Fork, a Join and a Sink
 * are absent from it entirely -- they hold no jobs, so they have no row -- yet
 * jobs flow through them and the flow is what a caller sizing a link or an
 * interconnect needs. This table has one row per node and reports the arrival
 * rate and the throughput at every one of them, which is the only place those
 * two numbers exist for a non-station node.
 *
 * QLen, Util, RespT and ResidT ARE the station numbers, scattered to the node
 * indices and left at zero elsewhere; that is the reference's own construction
 * and not a gap, because a node that is not a station holds no jobs and serves
 * nothing. ArvR and Tput are the two the reference recomputes, through
 * `sn_get_node_arvr_from_tput` and `sn_get_node_tput_from_tput`.
 *
 * THE F REGION PSEUDO-NODE ROWS OF THE REFERENCE ARE NOT EMITTED. MATLAB
 * appends one row per finite-capacity region, filled from `result.Avg` rows
 * M+1..M+F; the C++ `AvgResult` carries no per-region queue length or
 * utilization, so those rows have no data source here and are omitted rather
 * than fabricated as zeros, which would read as an empty region.
 */
/* `avg_result_from_sim` MOVED to line/solvers/solver_node_tables.h, where the
 * example corpus can reach it too: four `cache_replc_*` twins print the NODE
 * table their references print, and it is the same view of the same numbers.
 * Named unqualified below, as it was when it was defined here. */
using line::solvers::avg_result_from_sim;

/**
 * The LDES result document mapped onto the station AvgResult, BY NAME.
 *
 * The engine reports its own station and class order, and pairing the two off
 * positionally would put one station's numbers on another's row; `-a avg`
 * already matches by name for exactly that reason and this is the same match.
 * `WN` is recomputed rather than taken, for the reason the `-a avg` arm states:
 * the engine counts one visit per station, so its residence time is the response
 * time whenever a visit ratio is not 1.
 *
 * REGION ROWS RIDE PAST THE STATIONS, in the (M+F) layout `jmt_map_measures`
 * already uses, so a caller reads both engines' finite-capacity rows the same
 * way.
 */
inline line::mva::AvgResult<double> avg_result_from_ldes(
    const line::qn::NetworkStruct<double>& sn, const line::ldes::LdesResult& a) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const std::size_t F = a.QNfcr.empty() ? 0 : a.QNfcr.rows();
    line::mva::AvgResult<double> r;
    line::Matrix<double>* dst[6] = {&r.QN, &r.UN, &r.RN, &r.TN, &r.AN, &r.WN};
    const line::Matrix<double>* src[6] = {&a.QN, &a.UN, &a.RN, &a.TN, &a.AN, &a.WN};
    const line::Matrix<double>* fcr[6] = {&a.QNfcr, &a.UNfcr, &a.RNfcr,
                                          &a.TNfcr, &a.ANfcr, &a.WNfcr};
    std::vector<std::size_t> st_of(a.station_names.size(), 0);  // 1-based, 0 = unmatched
    for (std::size_t i = 0; i < a.station_names.size(); ++i)
        for (std::size_t j = 0; j < M; ++j)
            if (sn.stations[j].name == a.station_names[i]) { st_of[i] = j + 1; break; }
    std::vector<std::size_t> cl_of(a.class_names.size(), 0);
    for (std::size_t c = 0; c < a.class_names.size(); ++c)
        for (std::size_t j = 0; j < K; ++j)
            if (sn.classes[j].name == a.class_names[c]) { cl_of[c] = j + 1; break; }
    for (int m = 0; m < 6; ++m) {
        *dst[m] = line::Matrix<double>(M + F, K, 0.0);
        for (std::size_t i = 0; i < a.station_names.size(); ++i)
            for (std::size_t c = 0; c < a.class_names.size(); ++c)
                if (st_of[i] && cl_of[c])
                    (*dst[m])(st_of[i] - 1, cl_of[c] - 1) = ldes_at(*src[m], i, c);
        for (std::size_t f = 0; f < F; ++f)
            for (std::size_t c = 0; c < a.class_names.size(); ++c)
                if (cl_of[c]) (*dst[m])(M + f, cl_of[c] - 1) = ldes_at(*fcr[m], f, c);
    }
    line::Matrix<double> RNs(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) RNs(i, c) = r.RN(i, c);
    const line::Matrix<double> WNs = line::mva::sn_get_residt_from_respt(sn, RNs);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) r.WN(i, c) = WNs(i, c);
    for (std::size_t c = 0; c < K; ++c) {
        r.CN.push_back(ldes_at(a.CN, 0, c));
        r.XN.push_back(ldes_at(a.XN, 0, c));
    }
    r.method = a.method;
    r.actualmethod = a.method;
    return r;
}

/**
 * Solve for the station AvgResult with the named engine.
 *
 * SHARED BY EVERY @@NetworkSolver TABLE THAT IS NOT THE AvgTable -- `-a node`,
 * `-a sys`, `-a chain`, `-a nodechain` -- because each of them is a VIEW of the
 * same solved result and differs only in how it is indexed and aggregated.
 * Giving each arm its own engine ladder would let the five drift apart in which
 * knobs they honour, which is exactly the silent divergence the Knobs struct
 * exists to prevent one layer up.
 *
 * SSA and Fluid reach it through `avg_result_from_sim`, and only under
 * `--arith double`: both integrate transcendental quantities (exponential
 * clocks, an LSODA drift), so the dispatcher refuses the other backends by name
 * before the ladder is entered, exactly as their own `-a avg` arms do.
 *
 * THE TWO WRAPPERS ARE HERE FOR THE SAME REASON THE TWO SIMULATORS ARE.
 * `getAvgNodeTable` is @@NetworkSolver's and is a VIEW of whatever AvgResult a
 * solver returned; JMT and LDES both return one (`JmtResult::avg`, the LDES
 * result document), so refusing them the four views said "no C++ simulation
 * engine" about engines this port drives. `file` is what LDES needs and the
 * others ignore: it hands the model DOCUMENT to an external process rather
 * than reading the struct.
 */
template <class T>
line::mva::AvgResult<T> run_avg_engine(const line::qn::NetworkStruct<T>& sn, const Knobs& k,
                                       const std::string& s, std::string& banner,
                                       std::string* suffix = nullptr,
                                       const std::string* file = nullptr) {
    line::mva::AvgResult<T> r;
    if (s == "jmt" || s == "ldes") {
        // COMPILE-TIME, as for the two simulators: both wrappers report in
        // double and the dispatcher has already refused every other arithmetic
        // by name, so the discarded branch is unreachable rather than narrowed.
        if constexpr (std::is_same_v<T, double>) {
            if (s == "jmt") {
                line::jmt::JmtOptions o;
                if (!k.method.empty() && k.method != "default") o.method = k.method;
                if (k.samples > 0) o.samples = static_cast<double>(k.samples);
                if (k.seed != 0) o.seed = static_cast<long>(k.seed);
                o.keep = k.keep;
                if (k.t1 >= 0.0) o.max_simulated_time = k.t1;
                o.verbose = k.verbose;
                const line::jmt::JmtResult<double> a = line::jmt::solver_jmt_run_analyzer(sn, o);
                r = a.avg;
                banner = "SolverJMT";
                // THE SEED IS PART OF THE ANSWER, as it is for SSA: two runs of
                // a simulation are the same measurement only if both state it.
                if (suffix) {
                    char buf[128];
                    std::snprintf(buf, sizeof(buf), " samples=%g seed=%ld", o.samples, o.seed);
                    *suffix = buf;
                }
            } else {
                const line::ldes::LdesOptions o = ldes_options(k);
                const line::ldes::LdesResult a =
                    ldes_run(file ? *file : std::string(), o, std::vector<std::string>());
                r = avg_result_from_ldes(sn, a);
                banner = "SolverLDES";
                if (suffix) {
                    char buf[160];
                    std::snprintf(buf, sizeof(buf), " engine=%s samples=%zu seed=%ld",
                                  a.engine.empty() ? "?" : a.engine.c_str(), o.samples, o.seed);
                    *suffix = buf;
                }
            }
        } else {
            throw line::UnsupportedError("-s " + s + " runs under --arith double only");
        }
    } else if (s == "ssa" || s == "fluid") {
        // COMPILE-TIME, not just run-time: both runners static_assert on
        // transcendental arithmetic inside (an exponential clock, a square root
        // in the Cox refit), so instantiating them at Rational is a hard error
        // and not a refusal. The dispatcher has already rejected every arith
        // but double by name, so the discarded branch is unreachable rather
        // than silently narrowed.
        if constexpr (std::is_same_v<T, double>) {
            if (s == "ssa") {
                line::ssa::SsaOptions opt;
                if (!k.method.empty() && k.method != "default") opt.method = k.method;
                if (k.samples) opt.samples = k.samples;
                if (k.seed) opt.seed = k.seed;
                // The cache write-back rides beside the metric table, for the
                // reason `node_metrics` states: the realized hit and miss shares
                // are what the simulation MEASURED, and without them the node
                // table falls back to the split `link()` offered.
                std::vector<line::ssa::SsaCacheRatio> cache;
                const line::ssa::SsaSolution a = line::ssa::solver_ssa(sn, opt, &cache);
                // ResidT and ArvR are derived from the VISITS, and a cache's
                // split is routing, so both are taken on the struct carrying the
                // measured hit/miss shares rather than on `link()`'s even offer.
                r = avg_result_from_sim<T>(line::ssa::sn_with_ssa_cache_split<T>(sn, cache),
                                           a.QN, a.UN, a.RN, a.TN, a.CN, a.XN, a.method);
                r.cache = line::ssa::cache_metrics_of_ssa<T>(sn, cache);
                // THE SEED AND THE SAMPLE COUNT ARE PART OF THE ANSWER, not of
                // the invocation, so they ride in the banner here as they do in
                // `-a avg`: two runs of a simulation are the same measurement
                // only if both are stated, and a parity row that quotes one of
                // these tables has to carry them with it.
                banner = "SolverSSA";
                // APPENDED, NOT PREFIXED. Every consumer recognises a banner by
                // `Solver<name> arith=`, so a fact wedged between the two makes
                // the table belong to no solver at all -- which is how the
                // parity harness lost the whole SSA section.
                if (suffix) {
                    char buf[128];
                    std::snprintf(buf, sizeof(buf), " samples=%zu seed=%lu time=%.6g", a.samples,
                                  static_cast<unsigned long>(opt.seed), a.simulated_time);
                    *suffix = buf;
                }
            } else {
                const line::fluid::FluidOptions opt = fluid_options(k);
                // The cache decomposition renormalizes the self-switch at the
                // converged split; `-a node` needs that struct or it reports the
                // 1/2-1/2 `link()` left behind (see `node_metrics`).
                line::qn::NetworkStruct<T> refreshed;
                // The null must be typed: a bare `nullptr` is `std::nullptr_t`
                // and blocks deduction of T from the fourth argument.
                // THE SPLIT ITSELF IS TAKEN TOO, not only the struct it
                // renormalized: `node_metrics` prefers the stated split over the
                // visit ratios, and `-a cache` is built from nothing else.
                line::solvers::CacheMetrics<T> cache;
                const line::fluid::FluidSolution a = line::fluid::solver_fluid_run_analyzer(
                    sn, opt, static_cast<line::qn::NetworkStruct<T>*>(nullptr), &refreshed, &cache);
                // Non-empty only where the cache branch ran, which is the same
                // test `node_metrics` makes on the pointer. IT IS ALSO THE
                // STRUCT THE ARRIVAL RATE IS READ FROM: that column is derived
                // from the class-expanded routing, and the base struct still
                // holds the 1/2-1/2 self-switch `link()` offered, so deriving it
                // there reports 0.5/0.5 where cache_replc_routing's Delay1 sees
                // 0.4/0.6. Every other column is indexed by station and is the
                // same in both structs.
                const bool has_ref = !refreshed.nodes.empty();
                r = avg_result_from_sim<T>(has_ref ? refreshed : sn, a.QN, a.UN, a.RN, a.TN, a.CN,
                                           a.XN, a.method);
                if (has_ref) r.refreshed_struct.reset(new line::qn::NetworkStruct<T>(refreshed));
                r.cache = cache;
                r.iter = static_cast<int>(a.iters);
                banner = "SolverFluid";
                if (suffix) {
                    char buf[64];
                    std::snprintf(buf, sizeof(buf), " iters=%zu", a.iters);
                    *suffix = buf;
                }
            }
        } else {
            throw line::UnsupportedError("-s " + s + " runs under --arith double only");
        }
    } else if (s == "nc") {
        line::nc::NcSolverOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        if (!k.fork_join.empty()) opt.fork_join = k.fork_join;
        r = line::nc::solver_nc_run_analyzer(sn, opt);
        banner = "SolverNC";
    } else if (s == "mam") {
        line::mam::MamOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        r = line::mam::solver_mam_run_analyzer(sn, opt);
        banner = "SolverMAM";
    } else if (s == "ag") {
        line::ag::AgOptions opt;
        apply_ag_knobs(k, opt);
        r = line::ag::solver_ag_run_analyzer(sn, opt);
        banner = "SolverAG";
    } else if (s == "ba") {
        line::ba::BaOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        r = line::ba::solver_ba_run_analyzer(sn, opt);
        banner = "SolverBA";
    } else if (s == "ctmc") {
        line::ctmc::CtmcOptions opt;
        if (!k.method.empty()) opt.method = k.method;
        if (k.cutoff >= 0.0) opt.cutoff = k.cutoff;
        opt.cutoff_mat = k.cutoff_mat;
        opt.force = k.force;
        const line::ctmc::CtmcAnySolution<T> a = line::ctmc::solver_ctmc_analyzer_any(sn, opt);
        r = line::ctmc::solver_ctmc_avg_table(sn, a.sol, opt.method);
        banner = "SolverCTMC";
    } else {
        line::mva::MvaOptions opt;
        if (!k.method.empty() && k.method != "default") opt.method = k.method;
        if (k.tol >= 0.0) opt.tol = k.tol;
        if (k.iter_tol >= 0.0) opt.iter_tol = k.iter_tol;
        if (k.iter_max >= 0) opt.iter_max = k.iter_max;
        // `-a avg` has honoured --multiserver since the flag existed; these
        // views are the SAME solve indexed differently, so ignoring it here made
        // `-a chain` answer a different model than `-a avg` for the same command
        // line. Measured on cqn_repairmen_multi, where softmin and the default
        // rule differ by 39% on the Delay queue length.
        if (!k.multiserver.empty()) opt.multiserver = k.multiserver;
        if (!k.fork_join.empty()) opt.fork_join = k.fork_join;
        line::Matrix<T> init;
        r = line::mva::solver_mva_run_analyzer(sn, opt, init);
        banner = "SolverMVA";
    }
    return r;
}

/* `NodeMetrics` / `node_metrics` MOVED to line/solvers/solver_node_tables.h;
 * see the note above `run_avg_engine`. Used unchanged by both arms below. */
using line::solvers::NodeMetrics;
using line::solvers::node_metrics;

template <class T>
int solve_model_node(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);

    const std::size_t I = sn.nodes.size(), R = sn.nclasses;
    const NodeMetrics<T> nm = node_metrics<T>(sn, r);
    const line::Matrix<T>&QNn = nm.QN, &UNn = nm.UN, &RNn = nm.RN, &WNn = nm.WN, &ANn = nm.AN,
                         &TNn = nm.TN;

    auto d = [](const T& v) { return line::num_traits<T>::to_double(v); };
    // A FINITE CAPACITY REGION IS NOT A NODE, and the reference still prints it
    // in this table: `getAvgNodeTable` appends one row per region past the
    // nodes, because a WAITQ region holds jobs that are in no station's QLen and
    // the model's population only balances once they are read. The rows ride
    // past the stations in the returned AvgResult -- the (M+F) layout both
    // wrappers report -- and no analytical solver fills them, so this block is
    // empty for every engine that does not measure a region.
    //
    // Util AND ArvR ARE NaN, NOT ZERO. A region has no server to be busy and no
    // arrival process of its own; the reference reports both as missing, and a
    // zero there would be a number the run never measured.
    const std::size_t F =
        r.QN.rows() > sn.nstations ? r.QN.rows() - sn.nstations : static_cast<std::size_t>(0);
    const double region_nan = std::numeric_limits<double>::quiet_NaN();
    auto region_name = [&](std::size_t f) {
        return (f < sn.regions.size() && !sn.regions[f].name.empty())
                   ? sn.regions[f].name
                   : "FCR" + std::to_string(f + 1);
    };
    // The reference's own filter, the region twin of the all-zero row test: a
    // region no job ever entered is absent rather than a row of zeros.
    auto region_empty = [&](std::size_t f, std::size_t c) {
        return !(d(r.QN(sn.nstations + f, c)) > 0.0 || d(r.RN(sn.nstations + f, c)) > 0.0 ||
                 d(r.TN(sn.nstations + f, c)) > 0.0);
    };
    if (g_json_output) {
        // The banner under `-o json` too, for print_chain_table's reason: the
        // envelope names the arithmetic and the method but never the SOLVER.
        std::printf("%s arith=%s method=%s nodes=%zu%s\n", banner.c_str(), line::num_traits<T>::name(),
                    r.actualmethod.c_str(), I, suffix.c_str());
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgNodeTable";
        p["indexBase"] = 0;
        for (const char* key : {"Node", "JobClass", "QLen", "Util", "RespT", "ResidT", "ArvR",
                                "Tput"})
            p[key] = line::reg::Json::array();
        for (std::size_t i = 0; i < I; ++i)
            for (std::size_t c = 0; c < R; ++c) {
                // The reference's own row filter: a node a class never reaches
                // is absent, not a row of zeros, exactly as in the AvgTable.
                if (d(QNn(i, c)) == 0.0 && d(UNn(i, c)) == 0.0 && d(RNn(i, c)) == 0.0 &&
                    d(WNn(i, c)) == 0.0 && d(ANn(i, c)) == 0.0 && d(TNn(i, c)) == 0.0)
                    continue;
                p["Node"].push_back(sn.nodes[i].name);
                p["JobClass"].push_back(sn.classes[c].name);
                p["QLen"].push_back(d(QNn(i, c)));
                p["Util"].push_back(d(UNn(i, c)));
                p["RespT"].push_back(d(RNn(i, c)));
                p["ResidT"].push_back(d(WNn(i, c)));
                p["ArvR"].push_back(d(ANn(i, c)));
                p["Tput"].push_back(d(TNn(i, c)));
            }
        for (std::size_t f = 0; f < F; ++f)
            for (std::size_t c = 0; c < R; ++c) {
                if (region_empty(f, c)) continue;
                p["Node"].push_back(region_name(f));
                p["JobClass"].push_back(sn.classes[c].name);
                p["QLen"].push_back(d(r.QN(sn.nstations + f, c)));
                p["Util"].push_back(region_nan);
                p["RespT"].push_back(d(r.RN(sn.nstations + f, c)));
                p["ResidT"].push_back(d(r.WN(sn.nstations + f, c)));
                p["ArvR"].push_back(region_nan);
                p["Tput"].push_back(d(r.TN(sn.nstations + f, c)));
            }
        emit_analysis<T>("node", p, r.actualmethod);
        return 0;
    }
    std::printf("%s arith=%s method=%s nodes=%zu%s\n", banner.c_str(), line::num_traits<T>::name(),
                r.actualmethod.c_str(), I, suffix.c_str());
    std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Node", "JobClass", "QLen", "Util",
                "RespT", "ResidT", "ArvR", "Tput");
    for (std::size_t i = 0; i < I; ++i)
        for (std::size_t c = 0; c < R; ++c) {
            if (d(QNn(i, c)) == 0.0 && d(UNn(i, c)) == 0.0 && d(RNn(i, c)) == 0.0 &&
                d(WNn(i, c)) == 0.0 && d(ANn(i, c)) == 0.0 && d(TNn(i, c)) == 0.0)
                continue;
            std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                        sn.nodes[i].name.c_str(), sn.classes[c].name.c_str(), d(QNn(i, c)),
                        d(UNn(i, c)), d(RNn(i, c)), d(WNn(i, c)), d(ANn(i, c)), d(TNn(i, c)));
        }
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t c = 0; c < R; ++c) {
            if (region_empty(f, c)) continue;
            std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                        region_name(f).c_str(), sn.classes[c].name.c_str(),
                        d(r.QN(sn.nstations + f, c)), region_nan, d(r.RN(sn.nstations + f, c)),
                        d(r.WN(sn.nstations + f, c)), region_nan, d(r.TN(sn.nstations + f, c)));
        }
    return 0;
}

/** NaN as this arithmetic spells it, the reference's "not computed" marker. */
template <class T>
T cache_nan() {
    return line::num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
}

/** `v[r]` when the vector reaches r, NaN otherwise -- the reference's nanGetAt. */
template <class T>
double cache_at(const std::vector<T>& v, std::size_t r) {
    if (r >= v.size()) return std::numeric_limits<double>::quiet_NaN();
    return line::num_traits<T>::to_double(v[r]);
}

/**
 * `-a cache`: `@@NetworkSolver/getAvgCacheTable`, per Cache node and READ class.
 *
 * ONE TOTAL ROW PER (node, read class), plus one row per cache list where the
 * solver reported a per-list breakdown and the cache has more than one list.
 * The `List` column tells them apart: 0 is the total over every list, l is
 * list l. Only the total row carries the delayed-hit and miss columns, because
 * a miss is a property of the cache and not of any one list.
 *
 * THE ArvR COLUMN IS THE RETRIEVAL FLOW, `arvr * (missprob + delayedprob)`, on
 * the total row and the raw read rate on a list row. That is the reference's
 * choice and it is Little-consistent with ResidT: the residence time reported
 * beside it is the expected retrieval latency, which only the requests that
 * actually retrieve wait for.
 *
 * A READ CLASS IS ONE WITH A HIT CLASS DEFINED. A class that never reads the
 * cache has no row at all rather than a row of zeros, which is the same rule
 * the AvgTable applies to a class that never visits a station.
 */
template <class T>
int solve_model_cache(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);
    if (r.cache.empty())
        throw line::UnsupportedError(
            "-a cache reports the per-Cache hit and miss table and this model has no Cache node, "
            "or the solver that ran analyzes none; SolverNC's cache branches are what fill it");

    // The read-class arrival rate is that class's SOURCE throughput: every read
    // request enters the cache, so this holds across solvers, including the
    // simulators where a delayed hit is not folded into the hit throughput.
    const NodeMetrics<T> nm = node_metrics<T>(sn, r);
    std::size_t srcnode = 0;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i)
        if (sn.nodes[i].nodetype == line::lang::NodeType::Source) { srcnode = i + 1; break; }

    struct Row {
        std::string node, cls;
        double list, listcap, items, hitp, dhitp, missp, hitr, dhitr, missr, arvr, residt, cost;
    };
    std::vector<Row> rows;
    const double dnan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t c = 0; c < r.cache.caches.size(); ++c) {
        const line::solvers::CacheNodeMetrics<T>& m = r.cache.caches[c];
        const typename std::map<std::size_t, line::qn::CacheParam<T> >::const_iterator it =
            sn.nodeparam.find(m.node);
        if (it == sn.nodeparam.end()) continue;
        const std::vector<std::size_t>& hitclass = it->second.hitclass;
        const std::size_t h = m.itemcap.size();
        double totcap = 0.0;
        for (std::size_t l = 0; l < h; ++l) totcap += m.itemcap[l];
        double totcost = dnan;
        if (!m.listcost.empty()) {
            totcost = 0.0;
            for (std::size_t l = 0; l < m.listcost.size(); ++l)
                totcost += line::num_traits<T>::to_double(m.listcost[l]);
        }

        for (std::size_t cl = 0; cl < sn.nclasses; ++cl) {
            if (cl >= hitclass.size() || hitclass[cl] == 0) continue;  // not a read class
            double ph = cache_at(m.hitprob, cl), pm = cache_at(m.missprob, cl),
                   pd = cache_at(m.delayedprob, cl);
            if (std::isnan(ph) && std::isnan(pm) && std::isnan(pd)) continue;
            if (std::isnan(ph)) ph = 0.0;
            if (std::isnan(pm)) pm = 0.0;
            if (std::isnan(pd)) pd = 0.0;
            const double arvr =
                srcnode ? line::num_traits<T>::to_double(nm.TN(srcnode - 1, cl)) : 0.0;
            const double lat = cache_at(m.latency, cl);

            Row t;
            t.node = sn.nodes[m.node - 1].name;
            t.cls = sn.classes[cl].name;
            t.list = 0;
            t.listcap = totcap;
            t.items = static_cast<double>(m.nitems);
            t.hitp = ph;
            t.dhitp = pd;
            t.missp = pm;
            t.hitr = arvr * ph;
            t.dhitr = arvr * pd;
            t.missr = arvr * pm;
            t.arvr = arvr * (pm + pd);
            t.residt = lat;
            t.cost = totcost;
            rows.push_back(t);

            // Per-list rows, only where a genuine multi-list breakdown exists.
            bool any = false;
            if (h > 1 && cl < m.hitproblist.rows())
                for (std::size_t l = 0; l < m.hitproblist.cols(); ++l)
                    if (!std::isnan(line::num_traits<T>::to_double(m.hitproblist(cl, l))))
                        any = true;
            if (!any) continue;
            for (std::size_t l = 0; l < h; ++l) {
                double phl = l < m.hitproblist.cols()
                                 ? line::num_traits<T>::to_double(m.hitproblist(cl, l))
                                 : dnan;
                if (std::isnan(phl)) phl = 0.0;
                Row u;
                u.node = t.node;
                u.cls = t.cls;
                u.list = static_cast<double>(l + 1);
                u.listcap = m.itemcap[l];
                u.items = t.items;
                u.hitp = phl;
                u.dhitp = dnan;
                u.missp = dnan;
                u.hitr = arvr * phl;
                u.dhitr = dnan;
                u.missr = dnan;
                u.arvr = arvr;
                u.residt = dnan;
                u.cost = l < m.listcost.size()
                             ? line::num_traits<T>::to_double(m.listcost[l])
                             : dnan;
                rows.push_back(u);
            }
        }
    }

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgCacheTable";
        p["indexBase"] = 0;
        for (const char* key : {"Node", "JobClass", "List", "ListCap", "Items", "HitProb",
                                "DelayedHitProb", "MissProb", "HitRate", "DelayedHitRate",
                                "MissRate", "ArvR", "ResidT", "ListCost"})
            p[key] = line::reg::Json::array();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            p["Node"].push_back(rows[i].node);
            p["JobClass"].push_back(rows[i].cls);
            p["List"].push_back(rows[i].list);
            p["ListCap"].push_back(rows[i].listcap);
            p["Items"].push_back(rows[i].items);
            p["HitProb"].push_back(rows[i].hitp);
            p["DelayedHitProb"].push_back(rows[i].dhitp);
            p["MissProb"].push_back(rows[i].missp);
            p["HitRate"].push_back(rows[i].hitr);
            p["DelayedHitRate"].push_back(rows[i].dhitr);
            p["MissRate"].push_back(rows[i].missr);
            p["ArvR"].push_back(rows[i].arvr);
            p["ResidT"].push_back(rows[i].residt);
            p["ListCost"].push_back(rows[i].cost);
        }
        emit_analysis<T>("cache", p, r.actualmethod);
        return 0;
    }
    std::printf("%s arith=%s method=%s caches=%zu\n", banner.c_str(), line::num_traits<T>::name(),
                r.actualmethod.c_str(), r.cache.caches.size());
    std::printf("%-14s %-12s %5s %8s %6s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "Node",
                "JobClass", "List", "ListCap", "Items", "HitProb", "DHitProb", "MissProb",
                "HitRate", "DHitRate", "MissRate", "ArvR", "ResidT", "ListCost");
    for (std::size_t i = 0; i < rows.size(); ++i)
        std::printf(
            "%-14s %-12s %5g %8g %6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g "
            "%10.6g\n",
            rows[i].node.c_str(), rows[i].cls.c_str(), rows[i].list, rows[i].listcap,
            rows[i].items, rows[i].hitp, rows[i].dhitp, rows[i].missp, rows[i].hitr,
            rows[i].dhitr, rows[i].missr, rows[i].arvr, rows[i].residt, rows[i].cost);
    return 0;
}

/**
 * `-a item`: `@@NetworkSolver/getAvgItemTable`, one row per (Cache, item, list).
 *
 * THE PER-ITEM OCCUPANCY, which only a solver that computes a genuine per-item
 * distribution has: the exact NC cache recursions and the delayed-hit retrieval
 * algorithms. Every other branch measures the aggregate hit probability and
 * never forms the item law, and the arm refuses rather than filling the column
 * with the uniform guess that would reproduce the same aggregate.
 *
 * `Cost` is `Size * Prob`, so summing it over the items of a list reproduces
 * that list's ListCost in the cache table -- which is what makes the two tables
 * checkable against each other.
 */
template <class T>
int solve_model_item(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);
    if (r.cache.empty())
        throw line::UnsupportedError(
            "-a item reports the per-item cache occupancy and this model has no Cache node, or "
            "the solver that ran analyzes none");

    const double dnan = std::numeric_limits<double>::quiet_NaN();
    struct Row {
        std::string node;
        double item, list, listcap, size, prob, cost, dhq, dhqf;
    };
    std::vector<Row> rows;
    for (std::size_t c = 0; c < r.cache.caches.size(); ++c) {
        const line::solvers::CacheNodeMetrics<T>& m = r.cache.caches[c];
        const std::size_t h = m.itemcap.size();
        // EITHER measurement earns the item its rows. A solver may form the
        // per-item occupancy (the NC and MVA cache recursions) or the per-item
        // delayed-hit queue length (the exact chain) and not the other, and
        // requiring both would drop the CTMC's whole table.
        if (h == 0 || (m.itemprob.rows() == 0 && m.delayedhitqlen.empty())) continue;
        const std::size_t nit =
            m.itemprob.rows() > 0 ? m.itemprob.rows() : m.delayedhitqlen.size();
        for (std::size_t i = 0; i < nit; ++i)
            for (std::size_t l = 0; l < h; ++l) {
                Row t;
                t.node = sn.nodes[m.node - 1].name;
                t.item = static_cast<double>(i + 1);
                t.list = static_cast<double>(l + 1);
                t.listcap = m.itemcap[l];
                t.size = i < m.itemsize.size() ? m.itemsize[i] : dnan;
                // Column 0 of `itemprob` is the MISS column, so list l is column
                // l+1; reading it as l would report every item one list too low.
                t.prob = (l + 1) < m.itemprob.cols()
                             ? line::num_traits<T>::to_double(m.itemprob(i, l + 1))
                             : dnan;
                t.cost = t.size * t.prob;
                t.dhq = i < m.delayedhitqlen.size()
                            ? line::num_traits<T>::to_double(m.delayedhitqlen[i])
                            : dnan;
                t.dhqf = i < m.delayedhitqlenfull.size()
                             ? line::num_traits<T>::to_double(m.delayedhitqlenfull[i])
                             : dnan;
                rows.push_back(t);
            }
    }
    if (rows.empty())
        throw line::UnsupportedError(
            "-a item needs a per-item occupancy law and this solve produced none; the NC/MVA "
            "cache recursions (isolated and integrated alike) and the delayed-hit retrieval "
            "algorithms compute the embedded one, SolverCTMC the time-weighted one, and the "
            "simulators none");

    if (g_json_output) {
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgItemTable";
        p["indexBase"] = 0;
        for (const char* key : {"Node", "Item", "List", "ListCap", "Size", "Prob", "Cost",
                                "DelayedHitQLen", "DelayedHitQLenFull"})
            p[key] = line::reg::Json::array();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            p["Node"].push_back(rows[i].node);
            p["Item"].push_back(rows[i].item);
            p["List"].push_back(rows[i].list);
            p["ListCap"].push_back(rows[i].listcap);
            p["Size"].push_back(rows[i].size);
            p["Prob"].push_back(rows[i].prob);
            p["Cost"].push_back(rows[i].cost);
            p["DelayedHitQLen"].push_back(rows[i].dhq);
            p["DelayedHitQLenFull"].push_back(rows[i].dhqf);
        }
        emit_analysis<T>("item", p, r.actualmethod);
        return 0;
    }
    std::printf("%s arith=%s method=%s rows=%zu\n", banner.c_str(), line::num_traits<T>::name(),
                r.actualmethod.c_str(), rows.size());
    std::printf("%-14s %6s %6s %8s %10s %12s %12s %14s %18s\n", "Node", "Item", "List", "ListCap",
                "Size", "Prob", "Cost", "DelayedHitQLen", "DelayedHitQLenFull");
    for (std::size_t i = 0; i < rows.size(); ++i)
        std::printf("%-14s %6g %6g %8g %10g %12.8g %12.8g %14.8g %18.8g\n", rows[i].node.c_str(),
                    rows[i].item, rows[i].list, rows[i].listcap, rows[i].size, rows[i].prob,
                    rows[i].cost, rows[i].dhq, rows[i].dhqf);
    return 0;
}

/**
 * `-a sys`: `@@NetworkSolver/getAvgSysTable`, one row per CHAIN.
 *
 * SysRespT is the chain's CYCLE TIME and SysTput the flow that completes it,
 * both measured at the chain's reference station -- not a column of the
 * AvgTable summed up. On a closed chain the two are tied by Little's law and
 * the table is the standard capacity-planning view: N = X * R.
 */
template <class T>
int solve_model_sys(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);
    const line::solvers::SysResult<T> sys = line::solvers::solver_get_avg_sys<T>(sn, r);
    const std::vector<std::string> cn = line::solvers::chain_names(sn.nchains);
    const std::vector<std::string> cc = line::solvers::chain_class_labels<T>(sn);
    auto d = [](const T& v) { return line::num_traits<T>::to_double(v); };

    if (g_json_output) {
        // The banner under `-o json` too, for print_chain_table's reason: the
        // envelope names the arithmetic and the method but never the SOLVER.
        std::printf("%s arith=%s method=%s chains=%zu%s\n", banner.c_str(), line::num_traits<T>::name(),
                    r.actualmethod.c_str(), cn.size(), suffix.c_str());
        line::reg::Json p = line::reg::Json::object();
        p["type"] = "AvgSysTable";
        p["indexBase"] = 0;
        for (const char* key : {"Chain", "JobClasses", "SysRespT", "SysTput"})
            p[key] = line::reg::Json::array();
        for (std::size_t c = 0; c < sn.nchains; ++c) {
            p["Chain"].push_back(cn[c]);
            p["JobClasses"].push_back(cc[c]);
            p["SysRespT"].push_back(d(sys.CN[c]));
            p["SysTput"].push_back(d(sys.XN[c]));
        }
        emit_analysis<T>("sys", p, r.actualmethod);
        return 0;
    }
    std::printf("%s arith=%s method=%s chains=%zu%s\n", banner.c_str(),
                line::num_traits<T>::name(), r.actualmethod.c_str(), sn.nchains, suffix.c_str());
    std::printf("%-10s %-24s %14s %14s\n", "Chain", "JobClasses", "SysRespT", "SysTput");
    for (std::size_t c = 0; c < sn.nchains; ++c)
        std::printf("%-10s %-24s %14.6g %14.6g\n", cn[c].c_str(), cc[c].c_str(), d(sys.CN[c]),
                    d(sys.XN[c]));
    return 0;
}

/** Render a station- or node-level chain table, as text or as the host's JSON. */
template <class T>
void print_chain_table(const char* key, const char* type, const char* rowlabel,
                       const std::vector<std::string>& rows,
                       const std::vector<std::string>& chains,
                       const std::vector<std::string>& classes,
                       const line::solvers::ChainResult<T>& t, const std::string& method,
                       const char* banner, const char* arith, const char* suffix = "") {
    auto d = [](const T& v) { return line::num_traits<T>::to_double(v); };
    if (g_json_output) {
        // THE BANNER IS PRINTED UNDER `-o json` TOO, as the AvgTable path does:
        // the envelope names the arithmetic and the method but never the
        // SOLVER, so a host that pairs a table with another codebase's by the
        // solver in its banner cannot attribute a bannerless one at all. Its
        // absence here made `-o json` unusable for these four analyses.
        std::printf("%s arith=%s method=%s chains=%zu%s\n", banner, arith, method.c_str(),
                    chains.size(), suffix);
        line::reg::Json p = line::reg::Json::object();
        p["type"] = type;
        p["indexBase"] = 0;
        for (const char* c : {rowlabel, "Chain", "JobClasses", "QLen", "Util", "RespT", "ResidT",
                              "ArvR", "Tput"})
            p[c] = line::reg::Json::array();
        // ROW-MAJOR OVER (row, chain), the reference's `(ist-1)*C+c` ordering,
        // so a host reading the two codebases' tables side by side indexes them
        // the same way.
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (std::size_t c = 0; c < chains.size(); ++c) {
                p[rowlabel].push_back(rows[i]);
                p["Chain"].push_back(chains[c]);
                p["JobClasses"].push_back(classes[c]);
                p["QLen"].push_back(d(t.QN(i, c)));
                p["Util"].push_back(d(t.UN(i, c)));
                p["RespT"].push_back(d(t.RN(i, c)));
                p["ResidT"].push_back(d(t.WN(i, c)));
                p["ArvR"].push_back(d(t.AN(i, c)));
                p["Tput"].push_back(d(t.TN(i, c)));
            }
        emit_analysis<T>(key, p, method);
        return;
    }
    std::printf("%s arith=%s method=%s chains=%zu%s\n", banner, arith, method.c_str(),
                chains.size(), suffix);
    std::printf("%-16s %-10s %-20s %12s %12s %12s %12s %12s %12s\n", rowlabel, "Chain",
                "JobClasses", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (std::size_t c = 0; c < chains.size(); ++c)
            std::printf("%-16s %-10s %-20s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                        rows[i].c_str(), chains[c].c_str(), classes[c].c_str(), d(t.QN(i, c)),
                        d(t.UN(i, c)), d(t.RN(i, c)), d(t.WN(i, c)), d(t.AN(i, c)),
                        d(t.TN(i, c)));
}

/**
 * `-a chain`: `@@NetworkSolver/getAvgChainTable`, the station table by CHAIN.
 *
 * EVERY ROW IS EMITTED, including the all-zero ones, unlike the AvgTable and the
 * AvgNodeTable. The reference builds this table with a full (M x C) grid and no
 * row filter, and a chain that is absent from a station is information -- it is
 * the shape of the routing -- where a class absent from a station in the
 * AvgTable is only the class's own scope.
 */
template <class T>
int solve_model_chain(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);
    const line::solvers::ChainResult<T> t = line::solvers::solver_get_avg_chain<T>(sn, r);
    std::vector<std::string> rows;
    for (std::size_t i = 0; i < sn.nstations; ++i) rows.push_back(sn.stations[i].name);
    print_chain_table<T>("chain", "AvgChainTable", "Station", rows,
                         line::solvers::chain_names(sn.nchains),
                         line::solvers::chain_class_labels<T>(sn), t, r.actualmethod,
                         banner.c_str(), line::num_traits<T>::name(), suffix.c_str());
    return 0;
}

/** `-a nodechain`: `@@NetworkSolver/getAvgNodeChainTable`, the node table by CHAIN. */
template <class T>
int solve_model_nodechain(const std::string& file, const Knobs& k, const std::string& s) {
    line::qn::Network<T> net = read_model<T>(file);
    const line::qn::NetworkStruct<T>& sn = net.get_struct();
    std::string banner, suffix;
    const line::mva::AvgResult<T> r = run_avg_engine<T>(sn, k, s, banner, &suffix, &file);
    const NodeMetrics<T> nm = node_metrics<T>(sn, r);
    const line::solvers::ChainResult<T> t =
        line::solvers::solver_get_avg_node_chain<T>(sn, nm.QN, nm.UN, nm.RN, nm.WN, nm.AN, nm.TN);
    std::vector<std::string> rows;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) rows.push_back(sn.nodes[i].name);
    print_chain_table<T>("nodechain", "AvgNodeChainTable", "Node", rows,
                         line::solvers::chain_names(sn.nchains),
                         line::solvers::chain_class_labels<T>(sn), t, r.actualmethod,
                         banner.c_str(), line::num_traits<T>::name(), suffix.c_str());
    return 0;
}

/** The four @@NetworkSolver tables that are VIEWS of one solved AvgResult. */
inline bool is_avg_view(const std::string& analysis) {
    return analysis == "node" || analysis == "sys" || analysis == "chain" ||
           analysis == "nodechain";
}

/**
 * Dispatch one of those four views at a fixed arithmetic.
 *
 * The wrapper arms reach the views through this rather than through the shared
 * ladder further down: they return before it, having validated their own knobs
 * against an engine the ladder knows nothing about. The arithmetic is fixed
 * because both wrappers report in double and have already refused every other
 * backend by name.
 */
template <class T>
int solve_avg_view(const std::string& file, const Knobs& k, const std::string& s,
                   const std::string& analysis) {
    if (analysis == "node") return solve_model_node<T>(file, k, s);
    if (analysis == "sys") return solve_model_sys<T>(file, k, s);
    if (analysis == "chain") return solve_model_chain<T>(file, k, s);
    return solve_model_nodechain<T>(file, k, s);
}

int solve_model_dispatch(const std::string& arith, const std::string& solver,
                         const std::string& analysis, const std::string& file, const Knobs& k) {
    std::string s = solver.empty() ? "auto" : solver;
    if (s != "mva" && s != "auto" && s != "fluid" && s != "fld" && s != "nc" && s != "mam" &&
        s != "ag" && s != "ba" && s != "ssa" && s != "ctmc" && s != "uq" && s != "env" &&
        s != "qns" && s != "ldes" && s != "jmt")
        throw line::UnsupportedError(
            "the model-solving path ports -s mva, nc, ctmc, mam, ag, ba, ssa, fluid, ldes, jmt, "
            "uq, env and qns (got '" + s + "'); other solvers remain API-only");
    // Every solver on this path but JMT and QNS runs in-process, so the flags
    // that describe an external solver's child process have nothing to act on.
    // Those two do run one, and both write a scratch directory --keep names:
    // `solve_model_jmt` forwards it to JmtOptions.keep, which is what leaves
    // model.jsim behind, and refusing it here made the one document a parity
    // difference has to be read from unobtainable.
    if ((k.verbose && s != "ldes") || k.remote || !k.remote_url.empty() ||
        (k.keep && s != "qns" && s != "jmt"))
        throw line::UnsupportedError(
            "--keep, --verbose, --remote and --remote-url describe the child process of an "
            "external solver; on this path -s jmt and -s qns run one and take --keep, -s ldes "
            "runs one and takes --verbose (which echoes the resolved engine command line), and "
            "no path takes --remote or --remote-url");
    if (k.timeout_seconds && s != "qns")
        throw line::UnsupportedError(
            "--timeout is the deadline of an external solver's child process; on this path only "
            "-s qns runs one");
    // --fork-join names an arm of the fork-join FIXED POINT, which only the
    // mean-value arms drive: a simulator walks the fork on its sample path and
    // a CTMC enumerates it, so neither has a transform to choose. Refused by
    // name rather than ignored, which would report the default arm's numbers
    // under the caller's choice.
    if (!k.fork_join.empty() && s != "mva" && s != "nc")
        throw line::UnsupportedError(
            "--fork-join selects the fork-join transform of the shared mean-value fixed point "
            "and is read by -s mva and -s nc; -s " + s +
            " either simulates or enumerates the fork and applies no transform");
    // REFUSED BY NAME RATHER THAN IGNORED, on the same grounds as --fork-join
    // above: a knob silently dropped reports the DEFAULT arm's numbers under
    // the caller's choice, which is the one outcome stating the flag exists to
    // rule out.
    if (k.warmupfrac >= 0.0 && s != "ssa")
        throw line::UnsupportedError(
            "--warmupfrac discards a leading fraction of a SIMULATED path before the means are "
            "taken and is read by -s ssa; -s " + s +
            " has no path to discard (the LDES engine takes --ldes-warmupfrac)");
    if (k.pstar > 0.0 && s != "fluid" && s != "fld")
        throw line::UnsupportedError(
            "--pstar is the exponent of the fluid p-norm smoothing of the drift and is read by "
            "-s fluid; -s " + s + " integrates no drift");
    if ((!k.busy_orders.empty() || !k.busy_subnet.empty()) && analysis != "busyperiod")
        throw line::UnsupportedError(
            "--busyperiod and --busyperiod-subnet name the orders and the subnetwork of "
            "-a busyperiod; got -a " + analysis);
    // Solver console: this dispatcher is the single point every model-solving
    // arm passes through, so the narrated run is opened here and closed by the
    // guard's destructor -- on an exception too, so a failed analysis still
    // reports what it had reached. The model name is not known before the file
    // is read, so the header names the file's model once the struct compiles.
    line::util::LineConsole::Run consoleRun(upper_tag(s), "", true);

    // ---- SolverUQ, BEFORE the per-solver knob ladder ----------------------
    // It is a wrapper, not an engine: `--method`, `--samples` and `--seed`
    // describe its DESIGN and the convergence knobs belong to whatever
    // `--uq-solver` names, so the ladder below -- which asks "does THIS solver
    // have a sample count" -- answers about the wrong solver here.
    if (s == "uq") {
        if (analysis != "avg" && analysis != "posterior" && analysis != "interval")
            throw line::UnsupportedError(
                "SolverUQ ports -a avg (the prior-weighted expectation), -a posterior (the "
                "per-design-point table) and -a interval (the support-only range); got '" +
                analysis + "'");
        if (k.uq_solver.empty())
            throw line::InputError(
                "-s uq needs --uq-solver: UQ computes nothing itself, it expands the Prior and "
                "runs another solver at each design point (the C++ spelling of UQ(model, "
                "@SolverMVA)). Naming one here by default would attribute the numbers to an "
                "engine the caller never chose");
        if (k.t1 >= 0.0 || k.node || !k.notation.empty())
            throw line::UnsupportedError(
                "--tspan, --node and --notation name a transient horizon, a stateful node and an "
                "ODE document; SolverUQ reports steady-state means over a design of models and "
                "has none of the three");
        if (k.no_interlocking || k.repeat > 0 || !k.layer_solver.empty())
            throw line::UnsupportedError(
                "--no-interlocking, --repeat and --layer-solver are options of the layered solver "
                "and apply to -i lqnx; a Network model has no layers to interlock");
        if (arith == "double") return solve_model_uq<double>(file, k, analysis);
        if (arith == "exact") return solve_model_uq<line::Rational>(file, k, analysis);
        if (arith == "real:16") return solve_model_uq<line::Real<16> >(file, k, analysis);
        if (arith == "real" || arith == "real:32")
            return solve_model_uq<line::Real<32> >(file, k, analysis);
        if (arith == "real:64") return solve_model_uq<line::Real<64> >(file, k, analysis);
        if (arith == "real:128") return solve_model_uq<line::Real<128> >(file, k, analysis);
        if (arith == "real:256") return solve_model_uq<line::Real<256> >(file, k, analysis);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    // ---- SolverENV, also BEFORE the per-solver knob ladder ----------------
    // It reads a DIFFERENT MODEL TYPE: an Environment envelope, whose stages
    // each hold a Network. The ladder below asks its questions of one Network
    // model -- `--cutoff` per open class, `--node` on a stateful node -- and
    // the model-level knobs it would validate belong to the STAGE solver here,
    // so ENV states its own refusals and routes around it.
    if (s == "env") {
        if (analysis != "avg")
            throw line::UnsupportedError(
                "SolverENV reports -a avg, the environment-blended means; getEnsembleAvg is its "
                "only metric entry in the reference too (got '" + analysis + "')");
        if (k.samples || k.seed)
            throw line::UnsupportedError(
                "--samples and --seed describe a simulation; SolverENV iterates a fixed point over "
                "transient stage solves and draws nothing");
        // `--cutoff` IS ADMITTED WITH `--stage-solver ctmc`, and only then: every
        // stage is then enumerated, and an open stage's chain has to be
        // truncated somewhere. It stays refused for a fluid ensemble, which
        // enumerates nothing.
        if (k.has_cutoff() && k.stage_solver != "ctmc")
            throw line::UnsupportedError(
                "--cutoff bounds the open population of an enumerated state space and applies to "
                "-s env only beside --stage-solver ctmc; the fluid stages of this ensemble "
                "enumerate no states");
        if (k.node || !k.notation.empty())
            throw line::UnsupportedError(
                "--node and --notation name a stateful node and an ODE document of ONE network; "
                "an Environment holds a network per stage and "
                "SolverENV reports the blend over them");
        if (k.no_interlocking || k.repeat > 0 || !k.layer_solver.empty())
            throw line::UnsupportedError(
                "--no-interlocking, --repeat and --layer-solver are options of the layered solver "
                "and apply to -i lqnx; an Environment has stages, not layers");
        if (k.t0 != 0.0)
            throw line::UnsupportedError(
                "--tspan on the ENV path states the transient HORIZON each stage solve integrates "
                "to, and every stage starts from its entry state at 0; a nonzero t0 would name a "
                "start the coupling has no state for");
        // The mean-field coupling integrates the stage drift with LSODA, which
        // is double; the state-vector one uniformizes a CTMC and carries the
        // whole ladder. Narrowing silently would report an `exact` banner over
        // a double solve, so the refusal names the coupling that decided it.
        const std::string coupling =
            (k.method.empty() || k.method == "default") ? "meanfield" : k.method;
        if (k.tran_points && (coupling == "statevec" || coupling == "blend"))
            throw line::UnsupportedError(
                "--tran-points is the mean-field coupling's quadrature grid; the state-vector "
                "coupling carries the whole joint law across a switch and sums over no such grid, "
                "so the value would be accepted and never used");
        // `statedep` is a C++-API method and not a file one: it needs a rate
        // function PER ARC (`Environment::set_env_rate_reset`), which is a
        // function of the stage exit metrics and has no representation in
        // model.json -- the reference cannot serialize `resetEnvRatesFun`
        // either. Reaching it from a file would find no hook and refuse deeper
        // in, with a message about an environment the caller never wrote.
        if (coupling == "statedep")
            throw line::UnsupportedError(
                "--method statedep makes each environment transition depend on the state its "
                "stage is left in, through a rate function per arc that no model.json can carry "
                "(the reference cannot serialize resetEnvRatesFun either); it is reachable from "
                "the C++ API, through Environment::set_env_rate_reset");
        if ((k.tran_points || k.t1 >= 0.0) && (coupling == "avg" || coupling == "dec"))
            throw line::UnsupportedError(
                "--tran-points and --tspan state the grid and the horizon of a TRANSIENT stage "
                "solve; the closed-form limits --method avg and --method dec solve in steady "
                "state and carry nothing across a switch, so both would be accepted and never "
                "used");
        if (arith != "double" && coupling != "statevec" && coupling != "blend")
            throw line::UnsupportedError(
                "SolverENV solves a stage with the fluid analyzer on every method but the "
                "state-vector one -- the mean-field coupling transiently, the avg and dec limits "
                "in steady state -- and that analyzer is LSODA's and therefore double; --arith " +
                arith +
                " reaches ENV only through the state-vector coupling (--method statevec or --method blend)");
        if (arith == "double") return solve_model_env<double>(file, k);
        if (arith == "exact") return solve_model_env<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_env<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32")
            return solve_model_env<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_env<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_env<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_env<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (!k.uq_solver.empty())
        throw line::UnsupportedError(
            "--uq-solver names the engine SolverUQ runs at each design point and applies to -s uq; "
            "'" + s + "' solves one model, not a design of them");
    if (k.tran_points)
        throw line::UnsupportedError(
            "--tran-points is the resolution of the transient grid SolverENV sums its stage exit "
            "metrics over and applies to -s env; '" + s + "' has no such quadrature");
    // ---- SolverLDES, BEFORE the per-solver knob ladder --------------------
    // It hands the model DOCUMENT to an external engine instead of reading it,
    // so the ladder below -- which validates knobs against a parsed struct --
    // asks its questions of a model this path never parses. LDES states its own
    // refusals here for the same reason UQ and ENV state theirs.
    const bool ldes_knob = !k.ldes_tranfilter.empty() || k.ldes_warmupfrac >= 0.0 ||
                           !k.ldes_cimethod.empty() || k.ldes_cnvgon || k.ldes_slotted ||
                           k.ldes_replications > 0 || k.ldes_numthreads > 0 ||
                           k.ldes_maxtime > 0.0 || !k.ldes_initsol.empty() ||
                           !k.ldes_rest_url.empty();
    if (ldes_knob && s != "ldes" && s != "auto")
        throw line::UnsupportedError(
            "the --ldes-* flags are the discrete-event engine's own settings (warmup filter, "
            "confidence-interval estimator, slot lattice, replications, warm-start placement) and "
            "apply to -s ldes; '" + s + "' has none of them");
    if (s == "jmt") {
        if (arith != "double")
            throw line::UnsupportedError(
                "SolverJMT is a client of the Java Modelling Tools engine, which simulates in "
                "double and reports in double; --arith " + arith +
                " would label a double answer with an arithmetic that never touched it");
        if (analysis != "avg" && analysis != "cdf" && analysis != "trancdf" &&
            analysis != "trancdfpasst" && analysis != "prob" && !is_avg_view(analysis))
            throw line::UnsupportedError(
                "-s jmt reports -a avg (the JSIM or JMVA mean table), its four views -a node, "
                "-a sys, -a chain and -a nodechain, -a cdf (the empirical response-time law read "
                "back from the JMT logs, preloaded at the rounded steady-state queue lengths), "
                "-a tran-cdf-respt and -a tran-cdf-passt (the same logged run from the default "
                "initial state, so the samples cover the transient) and -a prob (the time each "
                "declared state is held for along the logged trajectory); -a tranprob needs one "
                "run per replication and is not exposed here");
        const std::vector<std::string> valid = line::jmt::jmt_list_valid_methods();
        if (!k.method.empty() &&
            std::find(valid.begin(), valid.end(), k.method) == valid.end())
            throw line::UnsupportedError(
                "SolverJMT methods are default, jsim and the jmva family (jmva, jmva.amva, "
                "jmva.mva, jmva.recal, jmva.comom, jmva.chow, jmva.bs, jmva.aql, jmva.lin, "
                "jmva.dmlin); got '" + k.method + "'");
        if (k.tol >= 0.0 || k.iter_tol >= 0.0 || k.iter_max >= 0 || !k.multiserver.empty())
            throw line::UnsupportedError(
                "--tol, --iter_tol, --iter_max and --multiserver configure a fixed-point "
                "iteration; JSIM simulates a sample path and JMVA takes its tolerance from the "
                "exported document. The simulation's stopping rule is --samples");
        if (k.has_cutoff())
            throw line::UnsupportedError(
                "--cutoff truncates an enumerated state space; JMT enumerates none");
        // `--node` NAMES THE STATION `--state` OVERRIDES, and nothing else: the
        // JMT arms report tables over every station and class, and `-a prob`
        // does too. The pair is `getProbAggr(node, state_a)`'s two arguments,
        // where the reference substitutes the given counts into that station's
        // row and leaves every other station at its declared state.
        if ((k.node && analysis != "prob") || k.jobclass || !k.marg_states.empty())
            throw line::UnsupportedError(
                "--node, --class and --marg-states select the marginal law of one (node, class); "
                "the JMT arms report tables over every station and class, and -a prob takes "
                "--node only to say which station --state overrides");
        if (!k.state.empty() && analysis != "prob")
            throw line::UnsupportedError(
                "--state names the state a probability is asked about and applies to -a prob");
        if (!k.state.empty() && !k.node)
            throw line::UnsupportedError(
                "--state is the per-class job count of ONE station and needs --node to say "
                "which; a bare count vector cannot be matched against a whole network");
        if (!k.notation.empty() || !k.symbolic.empty() || k.equilibria)
            throw line::UnsupportedError(
                "--notation, --symbolic and --equilibria describe an exported ODE document; JMT "
                "integrates no ODE");
        if (!k.cdf_algorithm.empty())
            throw line::UnsupportedError(
                "--cdf-algorithm selects between the two sojourn-time INVERSIONS of -s nc; the "
                "JMT response-time law is the ecdf of the passages its loggers recorded and is "
                "not computed from a transform");
        if (is_avg_view(analysis)) return solve_avg_view<double>(file, k, "jmt", analysis);
        return solve_model_jmt(file, k, analysis);
    }
    if (s == "ldes") {
        if (arith != "double")
            throw line::UnsupportedError(
                "SolverLDES is a client of the SSJ engine, which simulates in double and reports "
                "in double; --arith " + arith +
                " would label a double answer with an arithmetic that never touched it");
        // 'parallel' asks the engine for INDEPENDENT REPLICATIONS and the mean
        // over them, which is what its parallel analyzer is; it is not a second
        // engine. It resolves to a replication count here, taking --ldes-
        // replications when given and 8 otherwise -- the same default the SSA
        // parallel analyzer uses. Mirrors SolverLDES.listValidMethods in every
        // codebase, which advertises exactly {default, parallel}.
        if (!k.method.empty() && k.method != "default" && k.method != "parallel")
            throw line::UnsupportedError(
                "SolverLDES has two methods, 'default' and 'parallel' (listValidMethods returns "
                "exactly those in every codebase); got '" + k.method + "'");

        if (k.tol >= 0.0 || k.iter_tol >= 0.0 || k.iter_max >= 0 || !k.multiserver.empty())
            throw line::UnsupportedError(
                "--tol, --iter_tol, --iter_max and --multiserver are the knobs of a fixed-point "
                "iteration; LDES simulates a sample path and iterates nothing. Its stopping rule "
                "is --samples, or --ldes-cnvgon with --ldes-cnvgtol");
        if (k.has_cutoff())
            throw line::UnsupportedError(
                "--cutoff truncates an enumerated state space; a simulator visits the states the "
                "sample path reaches and enumerates none");
        if (k.node || k.jobclass || !k.marg_states.empty())
            throw line::UnsupportedError(
                "--node, --class and --marg-states select a marginal law of one (node, class); the "
                "LDES arms report tables over every station and class, and -a prob is refused for "
                "want of a target state");
        if (!k.notation.empty() || !k.symbolic.empty() || k.equilibria)
            throw line::UnsupportedError(
                "--notation, --symbolic and --equilibria describe an exported ODE document; LDES "
                "integrates no ODE");
        if (!k.cdf_algorithm.empty())
            throw line::UnsupportedError(
                "--cdf-algorithm selects between the two sojourn-time INVERSIONS of -s nc; the "
                "LDES response-time law is the ecdf of the samples the engine recorded and is not "
                "computed from a transform");
        if (k.no_interlocking || k.repeat > 0 || !k.layer_solver.empty())
            throw line::UnsupportedError(
                "--no-interlocking, --repeat and --layer-solver are options of the layered solver "
                "and apply to -i lqnx; the LDES layered path runs in the JAR's own ensemble "
                "backend and has no JSON interface to reach from here");
        if (!k.sens_method.empty() || !k.sens_scheme.empty() || k.sens_step >= 0.0)
            throw line::UnsupportedError(
                "--sens-method, --sens-scheme and --sens-step configure the layered sensitivity "
                "table; LDES reports no sensitivity");
        if (analysis == "tran" && !(k.t1 >= 0.0))
            throw line::InputError(
                "-s ldes -a tran needs --tspan <t1> or --tspan <t0>:<t1>: a trajectory over an unstated "
                "horizon is not a quantity, and the engine only records buckets once a timespan "
                "makes the run transient");
        if (k.t1 >= 0.0 && analysis != "tran")
            throw line::UnsupportedError(
                "--tspan names the horizon of -a tran; -a sample runs over [0, --samples], which "
                "is the horizon runTransientJson uses, and the other arms are steady state");
        // `-a tran-cdf-*` NAMES THE TRANSIENT LAW AND IS THE STEADY-STATE ONE
        // HERE, because a simulator has only the samples it observed: the
        // reference's `getTranCdfRespT` reads the same `respTimeSamples` its
        // `getCdfRespT` reads, and its `getTranCdfPassT` is a one-line delegation
        // to `getTranCdfRespT`. Warned rather than refused, so a script written
        // against the JAR runs and its author is told what the curve is.
        if ((analysis == "trancdf" || analysis == "trancdfpasst") &&
            k.verbosity != "silent")
            std::fprintf(stderr,
                         "Warning: -a %s is the ecdf of the per-job response times the run "
                         "observed, the same curve -a cdf reports; the reference's LDES "
                         "getTranCdfRespT reads the same samples\n",
                         analysis.c_str());
        // method='parallel' resolved to its replication count, after every
        // other knob has been validated against the caller's own Knobs.
        Knobs kldes = k;
        if (kldes.method == "parallel" && kldes.ldes_replications <= 1) kldes.ldes_replications = 8;
        if (is_avg_view(analysis)) return solve_avg_view<double>(file, kldes, "ldes", analysis);
        return solve_model_ldes(file, kldes, analysis);
    }
    // SolverAUTO resolves to a real engine BEFORE the knob checks below, so a
    // model the chooser sends to SSA accepts --samples and one it sends to CTMC
    // accepts --cutoff: the checks must see the solver that will actually run.
    if (s == "auto") {
        // `chooseSolverHeur` picks an engine from a GETTER, and the age laws are
        // not one of its getters: a model whose table it would send to MVA does
        // not thereby have an AoI answer. Refusing by name here beats letting
        // the chosen engine refuse an analysis it was never asked about.
        if (analysis == "aoi")
            throw line::UnsupportedError(
                "-a aoi is the AoI branch of the fluid 'mfq' method and no other engine reports "
                "it, so SolverAUTO does not choose for it: ask for it by name with -s fluid");
        const AutoPlan plan = choose_auto_plan_dispatch(arith, file, analysis, k.method);
        // `delegate.m` runs the chosen solver and, when it fails, every other
        // feasible candidate in slot order. Reproduced here by re-entering this
        // dispatch with a CONCRETE method name, so each attempt is validated against
        // the knobs of the solver that will actually run it. Only the first
        // attempt carries the method the ranking gated on ('exact'): a method
        // name is a solver's own vocabulary and does not travel to the next.
        if (plan.order.size() == 1) {
            // A forced method name (a method family, or the Environment envelope)
            // leaves ONE solver, and its own diagnostic is then the whole
            // story: the reference rethrows it rather than wrapping it.
            Knobs kk = k;
            kk.method = plan.method;
            std::printf("SolverAUTO selected %s%s\n", plan.order[0].c_str(), plan.note.c_str());
            return solve_model_dispatch(arith, plan.order[0], analysis, file, kk);
        }
        std::string first_error;
        for (std::size_t i = 0; i < plan.order.size(); ++i) {
            Knobs kk = k;
            kk.method = (i == 0) ? plan.method : std::string();
            if (i == 0)
                std::printf("SolverAUTO selected %s%s\n", plan.order[i].c_str(),
                            plan.note.c_str());
            else
                std::printf("SolverAUTO retrying with %s\n", plan.order[i].c_str());
            try {
                return solve_model_dispatch(arith, plan.order[i], analysis, file, kk);
            } catch (const line::UnsupportedError& e) {
                if (first_error.empty()) first_error = plan.order[i] + ": " + e.what();
                std::printf("SolverAUTO: %s cannot serve this run (%s)\n", plan.order[i].c_str(),
                            e.what());
            }
        }
        throw line::UnsupportedError(
            "SolverAUTO: every candidate refused this run. The chosen engine reported -- " +
            first_error);
    }
    // ---- knobs the CHOSEN solver does not have are refused, not dropped ----
    // Accepting an option and discarding it is the one place the port would
    // answer a question it was not asked: the caller believes a value the
    // solver never saw. Every refusal below names the option and the solver.
    const bool is_sim = (s == "ssa");
    // `-s ctmc -a sample` walks the chain with an exponential clock, so it has a
    // run length and a stream in the same sense a simulation does; every other
    // CTMC analysis is a solve and still refuses both.
    // The cftp methods draw iid stationary states, so they too have a run length
    // and a stream; unlike -a sample the draw is the ANSWER there, which is why
    // --samples is required rather than defaulted for them.
    const bool is_cftp =
        (s == "ctmc" && (k.method == "cftp" || k.method == "cftp.approx"));
    const bool draws_samples = is_sim || (s == "ctmc" && analysis == "sample") || is_cftp;
    if (k.samples && !draws_samples)
        throw line::UnsupportedError("--samples applies to the simulation solver (-s ssa), to "
                                     "-s ctmc -a sample and to -s ctmc --method cftp; '" + s +
                                     "' has no sample count");
    if (k.seed && !draws_samples)
        throw line::UnsupportedError("--seed applies to the simulation solver (-s ssa), to "
                                     "-s ctmc -a sample and to -s ctmc --method cftp; '" + s +
                                     "' draws no random numbers");
    // --mdd-tol / --mdd-maxiter drive the level iteration, which only the mdd
    // method runs. Accepting them elsewhere would let a caller believe a
    // tolerance was applied to a solve that has no iteration in it.
    if ((k.mdd_tol > 0.0 || k.mdd_maxiter > 0) && !(s == "ctmc" && k.method == "mdd"))
        throw line::UnsupportedError(
            "--mdd-tol and --mdd-maxiter set the coupled level iteration of -s ctmc --method mdd; "
            "'" + s + " / " + (k.method.empty() ? std::string("default") : k.method) +
            "' iterates no levels");
    // --tspan names a transient horizon, which only the CTMC transient analyses
    // have; accepting it elsewhere would let a caller believe a horizon was used.
    // The fluid solver integrates a forward equation too, and its horizon is what
    // `kp` reports its covariance AT, so --tspan reaches it as well; every other
    // fluid method restarts from its own end state until the moved mass stops
    // changing, and a horizon there caps that iteration rather than naming a time.
    if (k.t1 >= 0.0 &&
        !(s == "ctmc" &&
          (analysis == "tran" || analysis == "tranprob" || analysis == "tranreward")) &&
        !(s == "fluid" || s == "fld") && !(s == "mam" && analysis == "tran"))
        throw line::UnsupportedError(
            "--tspan sets the horizon of a transient analysis and applies to -s ctmc -a tran, "
            "-a tranprob and -a tranreward, to -s mam -a tran and to -s fluid; '" + s + " / " +
            analysis + "' integrates no forward equation");
    // --node narrows an answer to one node's block of the state, which only the
    // per-node CTMC queries and the per-node MAM queue-length law have; every
    // other analysis reports the whole network.
    // `-s ctmc -a prob` takes it TOGETHER WITH --state and only then: the arm
    // reports every station's marginal, so a bare --node would narrow nothing,
    // while --state names a row of ONE node's own space and needs --node to say
    // whose.
    if (k.node && !(s == "ctmc" && (analysis == "tranprob" || analysis == "sample")) &&
        !(s == "ctmc" && analysis == "prob" && !k.state.empty()) &&
        !(s == "nc" && analysis == "prob" && !k.state.empty()) &&
        !(s == "ssa" && analysis == "sample") && !(s == "mam" && analysis == "prob") &&
        !((s == "mva" || s == "nc") && analysis == "marg"))
        throw line::UnsupportedError(
            "--node selects the stateful node a state query is labelled by and applies to -s ctmc "
            "-a tranprob and -a sample, to -s ctmc|nc -a prob beside --state, to -s ssa -a sample, "
            "to -s mam -a prob and to -s mva|nc -a marg; '" + s + " / " + analysis +
            "' reports the whole network");
    // --class and --marg-states are the remaining two arguments of getProbMarg,
    // and nothing else in the surface takes either: every other analysis reports
    // all classes, and no other law is evaluated at a caller-chosen job count.
    if ((k.jobclass || !k.marg_states.empty()) && !(s == "mva" && analysis == "marg"))
        throw line::UnsupportedError(
            "--class and --marg-states are the job class and the state list of getProbMarg and "
            "apply to -s mva -a marg; '" + s + " / " + analysis +
            "' reports every class over its own range");
    // --notation selects which document the ODE export writes, and only the
    // export writes one; every other analysis reports numbers, which have no
    // notation to choose.
    if (!k.notation.empty() && !((s == "fluid" || s == "fld") && analysis == "odes"))
        throw line::UnsupportedError(
            "--notation selects the form of the exported ODE document and applies to -s fluid -a "
            "odes; '" + s + " / " + analysis + "' exports no equations");
    // --symbolic and --equilibria select the computer-algebra backend and ask it
    // to solve f(x) = 0; only the Jacobian consults one.
    if ((!k.symbolic.empty() || k.equilibria) &&
        !((s == "fluid" || s == "fld") && analysis == "jacobian"))
        throw line::UnsupportedError(
            "--symbolic selects the computer-algebra backend and --equilibria asks it for the "
            "solutions of f(x) = 0; both apply to -s fluid -a jacobian, and '" + s + " / " +
            analysis + "' consults no backend");
    // ---- the five JAR knobs, each refused where it would be accepted and never
    // read. Same discipline as every knob above: a caller who passed one to an
    // arm that does not consult it would believe a setting had been applied.
    // `-s jmt -a prob` reads it too, and reads it DIFFERENTLY: an exact chain is
    // indexed by the ENCODED row of a node's own state space, while a JMT log
    // records per-class job counts and nothing else, so there the vector is
    // `getProbAggr(node, state_a)`'s per-class count. Both are "the state this
    // probability is about"; which encoding it is in follows the solver.
    if (!k.state.empty() &&
        !(analysis == "prob" && (s == "ctmc" || s == "auto" || s == "jmt" || s == "nc")))
        throw line::UnsupportedError(
            "--state names the state `getProb(node, state)` asks about and applies to -a prob "
            "under -s ctmc, -s nc and -s jmt, the arms whose answer is indexed by a state; '" + s +
            " / " + analysis + "' reports a mean or a law over all of them");
    if (k.events && analysis != "sample")
        throw line::UnsupportedError(
            "--events is the length of ONE sampled trajectory and applies to -a sample; use "
            "--samples for a solver's run length ('" + s + " / " + analysis + "')");
    if (!k.percentiles.empty() && !(s == "mam" && (analysis == "cdf" || analysis == "cdfpasst" ||
                                                   analysis == "perct")))
        throw line::UnsupportedError(
            "--percentiles names the levels getPerctRespT is read at and applies to -s mam -a "
            "perct-respt (and to the percentiles printed beside -a cdf); '" + s + " / " +
            analysis + "' inverts no response-time law");
    if (!k.reward_name.empty() && analysis != "rewardvalue")
        throw line::UnsupportedError(
            "--reward-name selects which declared reward -a reward-value returns the value "
            "function of; -a reward returns every reward's steady-state expectation and needs no "
            "name ('" + s + " / " + analysis + "')");
    if (k.timestep > 0.0 &&
        !(s == "ctmc" &&
          (analysis == "tran" || analysis == "tranprob" || analysis == "tranreward")))
        throw line::UnsupportedError(
            "--timestep is the fixed output grid of a transient CTMC solve, `options.timestep` of "
            "ctmc_transient.m, and applies to -s ctmc -a tran, -a tranprob and -a tranreward; the "
            "fluid "
            "and simulated transients report the points their own integrator or engine produced "
            "('" + s + " / " + analysis + "')");
    if ((!k.transient_method.empty() || k.fau_epsilon > 0.0 || k.fau_delta >= 0.0) &&
        !(s == "ctmc" &&
          (analysis == "tran" || analysis == "tranprob" || analysis == "tranreward")))
        throw line::UnsupportedError(
            "--transient-method (and --fau-epsilon / --fau-delta) selects how the CTMC forward "
            "equation is advanced, `options.config.transient_method` of "
            "solver_ctmc_transient_analyzer.m, and applies to -s ctmc -a tran, -a tranprob and "
            "-a tranreward; every other analysis solves no forward equation ('" +
            s + " / " + analysis + "')");
    // --cdf-algorithm selects how the sojourn law is inverted, which only the NC
    // response-time distribution does; the CTMC one is read off tagged chains
    // and has no such choice.
    if (!k.cdf_algorithm.empty() && !(s == "nc" && analysis == "cdf"))
        throw line::UnsupportedError(
            "--cdf-algorithm selects the sojourn-time inversion of the NC response-time "
            "distribution and applies to -s nc -a cdf; '" + s + " / " + analysis +
            "' inverts no generating function");
    // The passage flags name the two state sets of getCdfFirstPassT and of
    // getFirstPassTMoments, which are the only two arms that time a state-set
    // passage. --passage-method selects the inversion and so belongs to the
    // curve alone; the moments involve no inversion at all.
    const bool passage_arm =
        (s == "ctmc" && (analysis == "firstpasst" || analysis == "firstpasstmom"));
    if ((!k.passage_from.empty() || !k.passage_into.empty() || k.passage_orders > 0) &&
        !passage_arm)
        throw line::UnsupportedError(
            "--passage-from, --passage-into and --passage-orders name the state sets and the "
            "moment order of -s ctmc -a firstpasst / firstpasstmom; '" + s + " / " + analysis +
            "' times no state-set passage");
    if (!k.passage_method.empty() && !(s == "ctmc" && analysis == "firstpasst"))
        throw line::UnsupportedError(
            "--passage-method selects the transform inversion of -s ctmc -a firstpasst; '" + s +
            " / " + analysis + "' inverts none (the moments arm solves for them directly)");
    // --perm-engine selects the permanent estimator, which only the NC joint law
    // of the per-station totals uses; nothing else in the tree evaluates one.
    if (k.method_perm != "exact" && !(s == "nc" && analysis == "sysmarg"))
        throw line::UnsupportedError(
            "--perm-engine selects the permanent estimator of the NC joint total-queue-length "
            "law and applies to -s nc -a sysmarg; '" + s + " / " + analysis +
            "' evaluates no permanent");
    // The layered path's own knobs, refused here for the same reason every other
    // knob above is: a caller who passed one to a Network solve would believe a
    // setting was applied that no Network solver has. --sens-* is the exception:
    // getSensitivityTable is a @@NetworkSolver method, so it selects the branch
    // of `-s nc -a sens` as much as of the layered one.
    const bool nc_sens = s == "nc" && analysis == "sens";
    if (!nc_sens && (!k.sens_method.empty() || !k.sens_scheme.empty() || k.sens_step > 0.0))
        throw line::UnsupportedError(
            "--sens-method, --sens-scheme and --sens-step select the branch of a sensitivity "
            "table and apply to -i lqnx -a sens or to -s nc -a sens; '" + s + " / " + analysis +
            "' differentiates nothing");
    if (k.no_interlocking || k.repeat > 0 || !k.layer_solver.empty() ||
        !k.ln_transient.empty() || !k.ln_transient_channels.empty())
        throw line::UnsupportedError(
            "--no-interlocking, --repeat, --layer-solver and --ln-transient* are options "
            "of the layered solver and apply to -i lqnx; a Network model has no layers to "
            "interlock");
    if (s == "ba" && (k.tol >= 0.0 || k.iter_tol >= 0.0 || k.iter_max >= 0))
        throw line::UnsupportedError(
            "--tol, --iter_tol and --iter_max do not apply to -s ba: a bound is a closed form, "
            "with nothing to converge");
    // Silent acceptance is the defect these guard against: a caller who passed
    // a QRF table to another solver would believe a parameterisation was
    // applied that nothing read.
    if (s != "ba" && (!k.qrf_params.empty() || !k.qrf_alpha.empty()))
        throw line::UnsupportedError(
            "--qrf-params and --qrf-alpha parameterise the QRF reduction bounds and apply to "
            "-s ba; '" + s + "' solves no reduction program");
    if (s != "ba" && k.level > 0)
        throw line::UnsupportedError(
            "--level is the hierarchy level of the SolverBA bound families and applies to -s ba; "
            "'" + s + "' has no bound hierarchy");
    if (is_sim && (k.tol >= 0.0 || k.iter_tol >= 0.0 || k.iter_max >= 0))
        throw line::UnsupportedError(
            "--tol, --iter_tol and --iter_max do not apply to -s ssa: a sample path is not an "
            "iteration; use --samples to set its length");
    if (s == "mam" && k.iter_tol >= 0.0)
        throw line::UnsupportedError(
            "--iter_tol is not a SolverMAM option (MamOptions carries tol and iter_max); "
            "use --tol");
    if (s == "ag" && k.iter_tol >= 0.0)
        throw line::UnsupportedError(
            "--iter_tol is not a SolverAG option (AgOptions carries tol and iter_max, the "
            "tolerance and the sweep budget of the reversed-rate fixed point); use --tol");
    // --max-states BOUNDS AN OPEN AGENT'S QUEUE-LENGTH DIMENSION, and only the
    // RCAT agents have one: a closed class is bounded by its own population
    // instead, and 'inapinf' ignores the level entirely and solves the open
    // agents on the infinite state space. Accepting it elsewhere would be the
    // silent-acceptance defect these guards exist for -- a caller who passed it
    // to -s ctmc would believe a truncation applied that nothing truncated.
    if (k.max_states >= 0 && s != "ag")
        throw line::UnsupportedError(
            "--max-states truncates the queue-length dimension of a SolverAG agent and applies "
            "to -s ag; '" + s + "' truncates no agent (use --cutoff for a CTMC state space)");
    // The cutoff BOUNDS A STATE SPACE, and only SolverCTMC and the MAM
    // queue-length law have one -- the latter because an OPEN queue's level
    // process is unbounded and `getProb` has to stop somewhere. Accepting it
    // elsewhere would be the silent-acceptance defect: a caller who passed it to
    // -s mva would believe the answer was truncated when nothing truncated it.
    if (k.has_cutoff() && s != "ctmc" && !(s == "mam" && analysis == "prob") && s != "env")
        throw line::UnsupportedError(
            "--cutoff bounds the open population of a CTMC state space, and the level truncation "
            "of -s mam -a prob; '" + s + " / " + analysis + "' enumerates no states");
    // THE MATRIX SPELLING IS NARROWER THAN THE SCALAR ONE. Only the CTMC state
    // space is enumerated per station, so only it can honour a per-station
    // bound; the MAM level truncation and an environment's stage cutoff are one
    // number each. Refused rather than reduced to a maximum, because that
    // silently answers a LARGER chain than the caller asked for.
    if (!k.cutoff_mat.empty() && s != "ctmc")
        throw line::UnsupportedError(
            "--cutoff as a per-(station,class) matrix bounds an enumerated state space per "
            "station and applies to -s ctmc; '" + s + "' takes one number");
    // `--stage-solver` names the solver each STAGE of an environment is run
    // with, and nothing else has stages: a layer of an LQN takes
    // `--layer-solver`, which is a different set for a different reason (a
    // layer is solved in steady state, a stage transiently).
    if (!k.stage_solver.empty() && s != "env")
        throw line::UnsupportedError(
            "--stage-solver names the solver each stage of a random environment is run with and "
            "applies to -s env; '" + s + "' has no stages");
    if (!k.stage_solver.empty() && k.stage_solver != "fluid" && k.stage_solver != "ctmc" &&
        k.stage_solver != "mam")
        throw line::UnsupportedError(
            "--stage-solver '" + k.stage_solver +
            "' is not available: the environment coupling needs a TRANSIENT stage solve, and only "
            "the fluid analyzer, the enumerated CTMC and the flattened LD-QBD provide one in this "
            "port");
    // `mam` IS THE STATE-VECTOR COUPLING'S BACKEND ONLY. The mean-field one
    // carries marginal means and reads them off a transient mean the LD-QBD
    // reduction does not produce; accepting it there would run the CTMC
    // ensemble under the MAM name.
    if (k.stage_solver == "mam" && k.method != "statevec")
        throw line::UnsupportedError(
            "--stage-solver mam applies to -s env --method statevec: the LD-QBD backend flattens "
            "its blocks into a generator the state-vector coupling propagates a distribution "
            "across, and the mean-field coupling carries marginal MEANS instead");
    // The two FJ_codes knobs configure ONE analyzer, solver_mam_fj, which the
    // MAM dispatch reaches on a homogeneous fork-join model. Accepting them
    // anywhere else would let a caller believe an accuracy setting had been
    // honoured by a solver that never read it.
    if ((k.fj_accuracy > 0 || !k.fj_tmode.empty()) && s != "mam")
        throw line::UnsupportedError(
            "--fj-accuracy and --fj-tmode configure the FJ_codes fork-join approximation of "
            "solver_mam_fj.m and apply to -s mam; '" + s + "' does not run it");
    // --timescale gates the slotted branch of the MAM dispatch alone. SolverNC
    // has its own discrete product form and takes --slotted for it, so a
    // caller that names a time scale for any other solver is told rather than
    // silently answered on the continuous one.
    if (!k.timescale.empty() && s != "mam")
        throw line::UnsupportedError(
            "--timescale selects the time scale of the MAM discrete-time path and applies to "
            "-s mam; '" + s + "' does not read it (SolverNC takes --slotted)");
    // ---- `-a node`, ahead of the per-solver whitelists ---------------------
    // getAvgNodeTable is @@NetworkSolver's, not any one solver's: it is the
    // station table scattered to the node index space plus the two flow columns
    // recomputed from it, so every solver that produces an AvgResult can answer
    // it and none of them needs its own arm.
    // `-a node`, `-a sys`, `-a chain` and `-a nodechain` are the four
    // @@NetworkSolver tables that are VIEWS of one solved AvgResult -- scattered
    // to nodes, aggregated to chains, or reduced to the reference station -- so
    // they share the engine whitelist and the arithmetic ladder. Adding a
    // per-arm copy of either would let the four drift on which solvers and
    // which arithmetics they accept, for tables built from the same numbers.
    // `-s ssa` and `-s fluid` return their own solution types rather than an
    // AvgResult; `run_avg_engine` bridges them (`avg_result_from_sim`), so the
    // reference's rule holds here too -- a solver that reports an AvgTable
    // reports its four views. The two external wrappers obey the same rule and
    // are dispatched in their own arms above, which run before this one: they
    // validate knobs the ladder here knows nothing about, and LDES is handed the
    // model DOCUMENT rather than a parsed struct.
    if (analysis == "node" || analysis == "sys" || analysis == "chain" ||
        analysis == "nodechain" || analysis == "cache" || analysis == "item") {
        const bool is_sim_engine = (s == "ssa" || s == "fluid");
        // The two cache tables are read off the SAME solved result, so they
        // belong to the same group; SolverNC is the only engine here whose
        // branches fill `AvgResult::cache`, and the arms say so when it is
        // empty rather than being whitelisted to nc alone -- `-s auto` on a
        // cache model resolves to nc, and refusing the token would refuse the
        // model.
        //
        // BOTH SIMULATORS ARE ADMITTED TO `-a cache`, and only there.
        // `run_avg_engine` fills `r.cache` for `-s ssa` from
        // `cache_metrics_of_ssa` -- the realized hit, delayed-hit and miss
        // SHARES the sample path measured -- and for `-s fluid` from the
        // cacheqn decomposition's converged split, which is the same quantity
        // its refreshed struct is renormalized at. Both are how
        // `SSA(model).getAvgCacheTable()` and `Fluid(model).getAvgCacheTable()`
        // answer in the reference, which reads them off the node the analyzer
        // wrote. `-a item` stays refused for both: the per-item occupancy is a
        // recursion of the
        // NC/MVA cache branches and no simulator forms it, so admitting it
        // would report an empty table for a quantity that was never measured.
        const bool cache_table = (analysis == "cache" || analysis == "item");
        const bool sim_cache_ok = (analysis == "cache");
        // `-a cache` and `-a item` are the CACHE tables, which RCAT does not
        // form -- so `ag` joins the AvgResult views and not those two.
        const bool ag_view = (s == "ag" && !cache_table);
        if ((s != "mva" && s != "auto" && s != "nc" && s != "mam" && s != "ba" && s != "ctmc" &&
             !ag_view && !is_sim_engine) ||
            (cache_table && is_sim_engine && !sim_cache_ok))
            throw line::UnsupportedError(
                // `-a node`'s own wording is kept verbatim: it is the message a
                // caller has been reading since the arm existed, and the group
                // it now shares does not change what it says.
                (analysis == "node"
                     ? std::string("-a node reports the per-node table of -s mva, nc, mam, ag, ba, "
                                   "ctmc, ssa, fluid, jmt, ldes and auto; '")
                     : "-a " + analysis +
                           " is a view of the station AvgResult and is reported by -s mva, nc, "
                           "mam, " + (analysis == "item" ? "" : "ag, ") + "ba, ctmc" +
                           (analysis == "item" ? "" : ", ssa, fluid, jmt, ldes") + " and auto; '") +
                s + "' does not return the station AvgResult it is built from");
        // The same refusal their own `-a avg` arms raise, for the same reason:
        // an SSA sample path is generated from exponential clocks and a fluid
        // trajectory is integrated by LSODA, so neither is carried by an exact
        // or an extended-precision backend. Refused BY NAME here rather than
        // narrowed silently in the ladder below.
        if (is_sim_engine && arith != "double")
            throw line::UnsupportedError(
                "-a " + analysis + " under -s " + s +
                " is read off a " + (s == "ssa" ? "sample path" : "fluid trajectory") +
                ", which is transcendental; rerun with --arith double (got '" + arith + "')");
        const std::string eng = (s == "auto") ? std::string("mva") : s;
#define LINE_CLI_TABLE_LADDER(FN)                                                     \
    do {                                                                              \
        if (arith == "double") return FN<double>(file, k, eng);                        \
        if (arith == "exact") return FN<line::Rational>(file, k, eng);                 \
        if (arith == "real:16") return FN<line::Real<16> >(file, k, eng);              \
        if (arith == "real" || arith == "real:32") return FN<line::Real<32> >(file, k, eng); \
        if (arith == "real:64") return FN<line::Real<64> >(file, k, eng);              \
        if (arith == "real:128") return FN<line::Real<128> >(file, k, eng);            \
        if (arith == "real:256") return FN<line::Real<256> >(file, k, eng);            \
    } while (0)
        if (analysis == "node") LINE_CLI_TABLE_LADDER(solve_model_node);
        if (analysis == "sys") LINE_CLI_TABLE_LADDER(solve_model_sys);
        if (analysis == "chain") LINE_CLI_TABLE_LADDER(solve_model_chain);
        if (analysis == "nodechain") LINE_CLI_TABLE_LADDER(solve_model_nodechain);
        if (analysis == "cache") LINE_CLI_TABLE_LADDER(solve_model_cache);
        if (analysis == "item") LINE_CLI_TABLE_LADDER(solve_model_item);
#undef LINE_CLI_TABLE_LADDER
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (s == "ctmc") {
        // The CTMC surface the port reaches: the AvgTable, plus the @@SolverCTMC
        // methods that are not means -- the generator and the state space
        // (getInfGen, getStateSpace), the transient occupancy (getTranProbSysAggr),
        // a marked trajectory (sampleSys), the declared rewards (getAvgReward),
        // the response-time laws (getCdfRespT, getCdfSysRespT) and the parametric
        // sensitivity (getSensitivityRanking).
        if (analysis != "avg" && analysis != "prob" && analysis != "gen" && analysis != "states" &&
            analysis != "tran" && analysis != "tranprob" && analysis != "tranreward" &&
            analysis != "sample" && analysis != "reward" && analysis != "rewardvalue" &&
            analysis != "cdf" && analysis != "sens" && analysis != "firstpasst" &&
            analysis != "firstpasstmom")
            throw line::UnsupportedError(
                "the CTMC solver ports -a avg, node, sys, chain, nodechain, prob, gen, states, "
                "tran, tranprob, tranreward, sample, reward, reward-value, cdf, first-passt, "
                "first-passt-moments and sens (got '" + analysis + "')");
        // The generator-free methods answer MEANS and nothing else: mdd holds
        // the reachable set in a diagram and cftp never enumerates it at all, so
        // there is no state space to list, no filtration to split and no
        // trajectory to walk. Refused by name rather than served from the
        // enumerated chain, which would report an answer under a method that did
        // not produce it.
        if ((k.method == "mdd" || k.method == "cftp" || k.method == "cftp.approx") &&
            analysis != "avg")
            throw line::UnsupportedError(
                "the '" + k.method +
                "' method never builds the explicit generator, so it serves -a avg only (got '" +
                analysis + "'); use --method default for the state-space analyses");
        if (k.tol >= 0.0 || k.iter_tol >= 0.0 || k.iter_max >= 0)
            throw line::UnsupportedError(
                "--tol, --iter_tol and --iter_max do not apply to -s ctmc: the stationary vector "
                "is obtained by a direct solve of pi Q = 0, with nothing to converge. The mdd "
                "method's level iteration has --mdd-tol and --mdd-maxiter of its own");
        // Every step from the generator to the means is a field operation, so
        // there is no arithmetic to refuse: exact returns the exact rational
        // stationary law. The transient analyses refuse inside, by name.
        if (arith == "double") return solve_model_ctmc<double>(file, k, analysis);
        if (arith == "exact") return solve_model_ctmc<line::Rational>(file, k, analysis);
        if (arith == "real:16") return solve_model_ctmc<line::Real<16> >(file, k, analysis);
        if (arith == "real" || arith == "real:32")
            return solve_model_ctmc<line::Real<32> >(file, k, analysis);
        if (arith == "real:64") return solve_model_ctmc<line::Real<64> >(file, k, analysis);
        if (arith == "real:128") return solve_model_ctmc<line::Real<128> >(file, k, analysis);
        if (arith == "real:256") return solve_model_ctmc<line::Real<256> >(file, k, analysis);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }

    // ---- the solvers that carry their own arithmetic restriction -----------
    // Each refuses BY NAME rather than being narrowed silently, and the reason
    // is the solver's own, not a limitation of the CLI.
    if (s == "mam") {
        // The @@SolverMAM surface the port reaches: the AvgTable, the
        // queue-length law (getProb / getProbMarg), the response-time law
        // (getCdfRespT and its getSjrnT / sjrnT aliases, plus getPerctRespT),
        // the transient means (getTranAvg) and the M/G/1-type internals of the
        // queue (getMAMResult).
        if (analysis != "avg" && analysis != "prob" && analysis != "cdf" &&
            analysis != "cdfpasst" && analysis != "perct" && analysis != "tran" &&
            analysis != "internals")
            throw line::UnsupportedError(
                "the MAM solver ports -a avg, node, prob, cdf, cdf-passt, perct-respt, tran and "
                "internals (got '" + analysis + "')");
        if (arith != "double")
            throw line::UnsupportedError(
                "the MAM solver fits phase-type representations, whose fitter requires "
                "transcendental arithmetic; rerun with --arith double (got '" + arith + "')");
        if (analysis == "tran" && k.t1 < 0.0)
            throw line::UnsupportedError(
                "-s mam -a tran integrates the transient queue length over a horizon and there is "
                "no default for it; pass --tspan t0 t1");
        if (analysis == "prob") return solve_model_mam_prob<double>(file, k);
        if (analysis == "cdf") return solve_model_mam_cdf<double>(file, k, "cdf", "CdfRespT");
        // `getCdfPassT` IS `getCdfRespT` here, and that is the reference's
        // construction rather than an alias invented in the CLI: SolverMAM.java
        // computes both from the SAME `solver_mam_passage_time(sn, sn.proc,
        // options)` call. It is emitted under its own key so a caller that
        // asked the passage-time question is answered it, and the payload's
        // `type` records which of the two names produced the curve.
        if (analysis == "cdfpasst")
            return solve_model_mam_cdf<double>(file, k, "cdfpasst", "CdfPassT");
        if (analysis == "perct") return solve_model_mam_perct<double>(file, k);
        if (analysis == "tran") return solve_model_mam_tran<double>(file, k);
        if (analysis == "internals") return solve_model_mam_internals<double>(file, k);
        return solve_model_mam<double>(file, k);
    }
    if (s == "ssa") {
        // `-a cdf` is refused BY NAME rather than falling into the generic
        // message, because the refusal is the reference's own answer and not a
        // port gap: `@@SolverSSA/getCdfRespT.m` raises the same error, since SSA
        // samples state trajectories and not per-job sojourn times.
        if (analysis == "cdf") line::ssa::ssa_cdf_respt_refuse();
        if (analysis != "avg" && analysis != "prob" && analysis != "sample")
            throw line::UnsupportedError("the SSA solver ports -a avg, -a prob and -a sample (got '" +
                                         analysis + "')");
        if (arith != "double")
            throw line::UnsupportedError(
                "an SSA sample path is generated from exponential clocks, which are "
                "transcendental; rerun with --arith double (got '" + arith + "')");
        if (analysis == "prob") return solve_model_ssa_prob<double>(file, k);
        if (analysis == "sample") return solve_model_ssa_sample<double>(file, k);
        return solve_model_ssa<double>(file, k);
    }
    if (s == "nc") {
        if (analysis != "avg" && analysis != "prob" && analysis != "marg" &&
            analysis != "sysmarg" && analysis != "cdf" && analysis != "sens" &&
            analysis != "normconst" && analysis != "busyperiod")
            throw line::UnsupportedError(
                "the NC solver ports -a avg, -a node, -a prob, -a marg, -a sysmarg, -a cdf, "
                "-a sens, -a normconst and -a busyperiod (got '" + analysis + "')");
        if (analysis == "busyperiod") {
            if (arith == "double") return solve_model_nc_busyp<double>(file, k);
            if (arith == "exact") return solve_model_nc_busyp<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_busyp<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_busyp<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_busyp<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_busyp<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_busyp<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "sysmarg") {
            if (arith == "double") return solve_model_nc_sysmarg<double>(file, k);
            if (arith == "exact") return solve_model_nc_sysmarg<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_sysmarg<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_sysmarg<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_sysmarg<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_sysmarg<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_sysmarg<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "marg") {
            if (arith == "double") return solve_model_nc_marg<double>(file, k);
            if (arith == "exact") return solve_model_nc_marg<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_marg<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_marg<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_marg<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_marg<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_marg<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        // Served under NC too so `-s auto -a normconst`, which chooseSolverHeur
        // sends to NC on a product-form model, lands on an arm that answers.
        if (analysis == "normconst") {
            if (arith == "double") return solve_model_normconst<double>(file, k, s);
            if (arith == "exact") return solve_model_normconst<line::Rational>(file, k, s);
            if (arith == "real:16") return solve_model_normconst<line::Real<16> >(file, k, s);
            if (arith == "real" || arith == "real:32")
                return solve_model_normconst<line::Real<32> >(file, k, s);
            if (arith == "real:64") return solve_model_normconst<line::Real<64> >(file, k, s);
            if (arith == "real:128") return solve_model_normconst<line::Real<128> >(file, k, s);
            if (arith == "real:256") return solve_model_normconst<line::Real<256> >(file, k, s);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "sens") {
            if (arith == "double") return solve_model_nc_sens<double>(file, k);
            if (arith == "exact") return solve_model_nc_sens<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_sens<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_sens<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_sens<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_sens<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_sens<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "cdf") {
            if (arith == "double") return solve_model_nc_cdf<double>(file, k);
            if (arith == "exact") return solve_model_nc_cdf<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_cdf<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_cdf<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_cdf<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_cdf<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_cdf<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "prob") {
            if (arith == "double") return solve_model_nc_prob<double>(file, k);
            if (arith == "exact") return solve_model_nc_prob<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_nc_prob<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_nc_prob<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_nc_prob<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_nc_prob<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_nc_prob<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (arith == "double") return solve_model_nc<double>(file, k);
        if (arith == "exact") return solve_model_nc<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_nc<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32") return solve_model_nc<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_nc<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_nc<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_nc<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (s == "ag") {
        // The @@SolverAG surface the port reaches: the AvgTable and its four
        // views, plus -a cdf as the inherited base-class exponential fallback.
        // RCAT converges a fixed point over the synchronization rates and
        // reports mean measures; it forms no state probability and no
        // transient, so there is nothing else to expose.
        if (analysis != "avg" && analysis != "cdf")
            throw line::UnsupportedError(
                "the AG solver ports -a avg (with its views -a node, -a sys, -a chain and "
                "-a nodechain) and -a cdf, the inherited exponential fallback: RCAT converges a "
                "fixed point over the synchronization rates and reports mean measures, forming no "
                "state probability or transient (got '" + analysis + "')");
        if (analysis == "cdf") {
            if (arith == "double") return solve_model_ag_cdf<double>(file, k);
            if (arith == "exact") return solve_model_ag_cdf<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_ag_cdf<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_ag_cdf<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_ag_cdf<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_ag_cdf<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_ag_cdf<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (arith == "double") return solve_model_ag<double>(file, k);
        if (arith == "exact") return solve_model_ag<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_ag<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32") return solve_model_ag<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_ag<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_ag<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_ag<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (s == "ba") {
        if (analysis != "avg" && analysis != "bounds" && analysis != "cdf")
            throw line::UnsupportedError("the BA solver ports -a avg, -a node, -a bounds and "
                                         "-a cdf, the inherited exponential fallback (got '" +
                                         analysis + "')");
        if (analysis == "cdf") {
            if (arith == "double") return solve_model_ba_cdf<double>(file, k);
            if (arith == "exact") return solve_model_ba_cdf<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_ba_cdf<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_ba_cdf<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_ba_cdf<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_ba_cdf<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_ba_cdf<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (analysis == "bounds") {
            if (arith == "double") return solve_model_ba_bounds<double>(file, k);
            if (arith == "exact") return solve_model_ba_bounds<line::Rational>(file, k);
            if (arith == "real:16") return solve_model_ba_bounds<line::Real<16> >(file, k);
            if (arith == "real" || arith == "real:32")
                return solve_model_ba_bounds<line::Real<32> >(file, k);
            if (arith == "real:64") return solve_model_ba_bounds<line::Real<64> >(file, k);
            if (arith == "real:128") return solve_model_ba_bounds<line::Real<128> >(file, k);
            if (arith == "real:256") return solve_model_ba_bounds<line::Real<256> >(file, k);
            throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
        }
        if (arith == "double") return solve_model_ba<double>(file, k);
        if (arith == "exact") return solve_model_ba<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_ba<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32") return solve_model_ba<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_ba<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_ba<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_ba<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (s == "qns") {
        if (analysis != "avg" && analysis != "cdf")
            throw line::UnsupportedError(
                "SolverQNS reports -a avg and -a cdf, the inherited exponential fallback: "
                "qnsolver returns one chain-level table of means and computes no state "
                "probability and no transient (got '" + analysis + "')");
        // The numbers arrive as the decimal text qnsolver printed, so every
        // digit past double is one this port invented; the ladder is refused
        // rather than run at a width the answer does not have.
        if (arith != "double")
            throw line::UnsupportedError(
                "SolverQNS reads its results back as the decimal text an external binary printed, "
                "which is double at best; --arith " + arith +
                " would report a precision the tool never produced");
        if (analysis == "cdf") return solve_model_qns_cdf(file, k);
        return solve_model_qns<double>(file, k);
    }
    if (s == "fluid" || s == "fld") {
        if (analysis != "avg" && analysis != "odes" && analysis != "var" &&
            analysis != "tranvar" && analysis != "jacobian" && analysis != "tran" &&
            analysis != "prob" && analysis != "cdf" && analysis != "aoi" &&
            analysis != "statevec")
            throw line::UnsupportedError(
                "the fluid solver ports -a avg, -a tran, -a tranvar, -a prob, -a cdf, -a aoi, "
                "-a odes, -a statevec, -a var and -a jacobian (got '" + analysis + "')");
        // The drift is integrated by LSODA, whose coefficients assume double;
        // a higher-precision request is refused rather than quietly narrowed.
        if (arith != "double")
            throw line::UnsupportedError(
                "the fluid solver integrates its drift with LSODA, which is double precision by "
                "construction; rerun with --arith double (got '" + arith + "')");
        if (analysis == "odes") return solve_model_fluid_odes<double>(file, k);
        if (analysis == "statevec") return solve_model_fluid_statevec<double>(file, k);
        if (analysis == "jacobian") return solve_model_fluid_jacobian<double>(file, k);
        if (analysis == "var") return solve_model_fluid_var<double>(file, k);
        if (analysis == "tranvar") return solve_model_fluid_tranvar<double>(file, k);
        if (analysis == "tran") return solve_model_fluid_tran<double>(file, k);
        if (analysis == "prob") return solve_model_fluid_prob<double>(file, k);
        if (analysis == "cdf") return solve_model_fluid_cdf<double>(file, k);
        if (analysis == "aoi") return solve_model_fluid_aoi<double>(file, k);
        return solve_model_fluid<double>(file, k);
    }
    if (analysis == "prob") {
        if (arith == "double") return solve_model_prob<double>(file, k);
        if (arith == "exact") return solve_model_prob<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_prob<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32") return solve_model_prob<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_prob<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_prob<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_prob<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (analysis == "marg") {
        if (arith == "double") return solve_model_marg<double>(file, k);
        if (arith == "exact") return solve_model_marg<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_marg<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32") return solve_model_marg<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_marg<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_marg<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_marg<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (analysis == "normconst") {
        if (arith == "double") return solve_model_normconst<double>(file, k, s);
        if (arith == "exact") return solve_model_normconst<line::Rational>(file, k, s);
        if (arith == "real:16") return solve_model_normconst<line::Real<16> >(file, k, s);
        if (arith == "real" || arith == "real:32")
            return solve_model_normconst<line::Real<32> >(file, k, s);
        if (arith == "real:64") return solve_model_normconst<line::Real<64> >(file, k, s);
        if (arith == "real:128") return solve_model_normconst<line::Real<128> >(file, k, s);
        if (arith == "real:256") return solve_model_normconst<line::Real<256> >(file, k, s);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (analysis == "cdf") {
        // The inherited base-class exponential fallback over the MVA means,
        // as @@NetworkSolver/getCdfRespT.m serves it for SolverMVA
        if (arith == "double") return solve_model_mva_cdf<double>(file, k);
        if (arith == "exact") return solve_model_mva_cdf<line::Rational>(file, k);
        if (arith == "real:16") return solve_model_mva_cdf<line::Real<16> >(file, k);
        if (arith == "real" || arith == "real:32")
            return solve_model_mva_cdf<line::Real<32> >(file, k);
        if (arith == "real:64") return solve_model_mva_cdf<line::Real<64> >(file, k);
        if (arith == "real:128") return solve_model_mva_cdf<line::Real<128> >(file, k);
        if (arith == "real:256") return solve_model_mva_cdf<line::Real<256> >(file, k);
        throw line::InputError("--arith '" + arith + "' is not a model-solve backend");
    }
    if (analysis != "avg")
        throw line::UnsupportedError(
            "the model-solving path ports -a avg, -a node, -a prob, -a marg, -a normconst and "
            "-a cdf (got '" +
            analysis + "')");
    // The mvaDispatch ladder's transcendental-only analyzers (open-queue closed
    // forms, DPS-exact, Marie, size-based) refuse by name under exact/real via
    // their if-constexpr guards, so the field-arithmetic branches (product-form
    // MVA, LD scaling) stay exact while a transcendental model is refused rather
    // than silently degraded.
    if (arith == "double") return solve_model_mva<double>(file, k);
    if (arith == "exact") return solve_model_mva<line::Rational>(file, k);
    if (arith == "real:16") return solve_model_mva<line::Real<16> >(file, k);
    if (arith == "real" || arith == "real:32") return solve_model_mva<line::Real<32> >(file, k);
    if (arith == "real:64") return solve_model_mva<line::Real<64> >(file, k);
    if (arith == "real:128") return solve_model_mva<line::Real<128> >(file, k);
    if (arith == "real:256") return solve_model_mva<line::Real<256> >(file, k);
    throw line::InputError(
        "--arith '" + arith +
        "' is not a model-solve backend; use double, exact or real:<16|32|64|128|256>");
}

/**
 * The no-argument and `-h` message: the flags a first run actually needs.
 *
 * The full reference below is ~460 lines, which is not a thing a human reads at
 * a prompt; it stays one flag away under `--help-all` rather than being the
 * first thing the binary says.
 */
void print_brief_help() {
    std::printf(
        "LINE solver (C++), version %s\n"
        "\n"
        "Usage: line-cli -f <model> [-s <solver>] [-a <analysis>] [-o <format>]\n"
        "       cat model.json | line-cli\n"
        "\n"
        "Common options:\n"
        "  -f, --file <path>      model file: .json (network), .lqnx (layered),\n"
        "                         .jsimg (JMT), .pnml (Petri net); stdin if omitted\n"
        "  -s, --solver <name>    auto (default), mva, nc, ctmc, mam, fluid, ssa,\n"
        "                         ldes, jmt, ag, ba, uq; ln for layered models, env\n"
        "                         for random environments\n"
        "  -a, --analysis <type>  avg (default), node, sys, chain, tran, prob, cdf,\n"
        "                         states, sample, normconst, bounds, ... (comma list)\n"
        "  -o, --output <fmt>     readable (default) | json\n"
        "      --method <name>    algorithm within the chosen solver\n"
        "      --samples <n>      simulation run length (ssa, ldes); default 10000\n"
        "      --seed <n>         random seed; default 23000\n"
        "  -v, --verbosity <lvl>  silent | standard | debug; debug turns on\n"
        "                         the solver console, a running progress log\n"
        "      --find-solver [m]  which solvers and methods can analyze this model,\n"
        "                         optionally only those answering measure <m>\n"
        "                         (avg, tran, cdf, prob, sample, ...); reports and exits\n"
        "      --find-solver-all [m]  the same, keeping the refused pairs and the\n"
        "                         reason each was refused\n"
        "  -h, --help             this message\n"
        "      --help-all         every option, solver by solver\n"
        "  -V, --version          version string\n"
        "      --install          environment check: which optional backends\n"
        "                         (Java/JMT, LQNS, qnsolver, SageMath) are reachable\n"
        "\n"
        "Examples:\n"
        "  line-cli -f model.json                       solve, letting auto pick the solver\n"
        "  line-cli -f model.json -s mva -a avg         mean queue lengths, MVA\n"
        "  line-cli -f model.lqnx -s ln                 solve a layered model\n"
        "  line-cli -f model.json -s ssa --samples 1e6 -o json\n"
        "  line-cli -f model.json --find-solver         what can solve this model\n"
        "  line-cli -f model.json --find-solver cdf     ... and return a passage-time law\n"
        "\n"
        "Each solver has flags of its own (tolerances, horizons, engine choices):\n"
        "run `line-cli --help-all` for the full reference.\n",
        kVersion);
}

void print_help() {
    std::printf(
        "LINE multiprecision solver (C++), version %s\n"
        "\n"
        "Usage: line-cli [OPTIONS]\n"
        "       cat model.json | line-cli -i json [OPTIONS]\n"
        "       line-cli model.lqnx [OPTIONS]\n"
        "\n"
        "Options (flag-compatible with jline.cli.LineCLI):\n"
        "  -f, --file <path>       model file; stdin when omitted (json only)\n"
        "  -i, --input <fmt>       input format: json | jsim | jsimg | jsimw |\n"
        "                          lqnx | xml | pnml. Without it a .lqnx or .xml\n"
        "                          path is read as a layered model, a\n"
        "                          .jsim/.jsimg/.jsimw path as a JMT simulation\n"
        "                          document, a .pnml path as a place/transition\n"
        "                          net (ISO/IEC 15909-2), and anything else as a\n"
        "                          Network model.json. The three jsim spellings\n"
        "                          name ONE format, as they do in JMT\n"
        "  -o, --output <fmt>      output format: readable | json (| layers, lqnx).\n"
        "                          json is honoured by EVERY analysis, not only\n"
        "                          -a avg: each answers under a key named after\n"
        "                          its -a, with the arithmetic and the resolved\n"
        "                          method beside it, and every index inside a\n"
        "                          payload is 0-based against the tables' 1-based\n"
        "                          columns (each payload states its indexBase)\n"
        "  -s, --solver <name>     Network: auto, mva, nc, ctmc, mam, ba, ssa, fluid, uq,\n"
        "                          ldes (the SSJ discrete-event engine, run as a\n"
        "                          subprocess on common/ldes or common/ldes.jar),\n"
        "                          jmt (the Java Modelling Tools engine),\n"
        "                          qns (the external qnsolver binary)\n"
        "                          layered: auto, ln, ln.mva, ln.comom, lqns,\n"
        "                          ldes (the native in-process LN simulator, a\n"
        "                          sample path of the layered model itself and\n"
        "                          not a decomposition into layers; it takes\n"
        "                          --samples and --seed, and NOT the --ldes-*\n"
        "                          family, which configures the subprocess\n"
        "                          engine that answers -s ldes on a Network)\n"
        "                          environment: env (an Environment model.json)\n"
        "  -a, --analysis <type>   analysis: avg, node, sys, chain, nodechain,\n"
        "                          stage, cache, item, prob, marg, cdf, cdf-passt,\n"
        "                          perct-respt, tran, tran-cdf-respt,\n"
        "                          tran-cdf-passt, tranprob, tranreward,\n"
        "                          reward-value, normconst, gen, states, sample,\n"
        "                          reward, sens, first-passt, odes, statevec,\n"
        "                          var, tranvar, busyperiod,\n"
        "                          jacobian, aoi, internals, bounds, posterior,\n"
        "                          interval. A COMMA LIST runs several in order\n"
        "                          (`-a avg,sys`), emitting one -o json envelope\n"
        "                          per analysis rather than one merged object.\n"
        "                          The JAR CLI's own spellings are accepted as\n"
        "                          aliases -- cdf-respt, prob-sys-aggr, tran-avg,\n"
        "                          generator, reward-steady, all -- and collapse\n"
        "                          onto the arm that already answers them whole:\n"
        "                          -a prob reports getProbSys, getProbSysAggr and\n"
        "                          the per-station pair together, -a sample walks\n"
        "                          all four samplers at once, and -a stage IS\n"
        "                          -a avg on a Network (one implicit stage).\n"
        "                          Which names a\n"
        "                          solver serves is stated by its own refusal;\n"
        "                          ssa serves avg, prob and sample (prob and\n"
        "                          sample run the SERIAL engine whatever -m\n"
        "                          said, since the NRM simulates counts rather\n"
        "                          than the state encoding)\n"
        "  -v, --verbosity <lvl>   silent | standard | debug; debug turns on\n"
        "                          the solver console, a running progress log\n"
        "                          of every solver run\n"
        "  -d, --seed <n>          random seed (SSA, ctmc --method cftp); default\n"
        "                          23000. -d is the JAR CLI's spelling\n"
        "      --warmupfrac <f>    SSA: leading fraction of the path discarded\n"
        "                          before the means are taken, in [0,1)\n"
        "      --pstar <p>         fluid: exponent of the p-norm smoothing of the\n"
        "                          drift; without it the hard min() is integrated\n"
        "      --busyperiod <n,..> -a busyperiod: the orders wanted; default 1\n"
        "      --busyperiod-subnet <i,..>\n"
        "                          -a busyperiod: the 1-based stations forming the\n"
        "                          subnetwork. Required -- a busy period is defined\n"
        "                          for a NAMED set and no default can choose one\n"
        "      --method <name>     algorithm within the solver\n"
        "      --samples <n>       simulation run length (SSA); default 10000.\n"
        "                          Accepts 1e6 as well as 1000000. A simulation\n"
        "                          figure is only a measurement WITH this number,\n"
        "                          which is why the SSA banner reports it back.\n"
        "                          REQUIRED by ctmc --method cftp, where the draw\n"
        "                          is the answer rather than a run length.\n"
        "      --mdd-tol <x>       level-iteration tolerance of ctmc --method mdd;\n"
        "                          default 1e-12. NOT --tol: that iteration is an\n"
        "                          inner solve whose fixed point is checked\n"
        "                          against the population invariant at 1e-6, so a\n"
        "                          solver-sized tolerance stops short of it\n"
        "      --mdd-maxiter <n>   coupled sweeps before ctmc --method mdd is\n"
        "                          declared non-convergent; default 500\n"
        "      --level <n>         hierarchy level of the ba pbh/cbh/sib families\n"
        "                          and the iteration count of pbk/bjbk; default 2\n"
        "      --qrf-params <j>    JSON (inline or a path) with the QRF blocking\n"
        "                          tables of -s ba --method qrf.bas|qrf.rsrd:\n"
        "                          f, MR, BB, MM, MM1, ZZ and optionally F, the\n"
        "                          fields sn_to_qrf_params assembles. ZM is\n"
        "                          derived from ZZ. There is no default: assuming\n"
        "                          no blocking puts the bound ~31x farther from\n"
        "                          exact, so its absence is refused\n"
        "      --qrf-alpha <j>     JSON (nstations x N) load-dependent scaling of\n"
        "                          the ba qrf.mmi.ld, qrf.mmi.linear and qrf.rsrd\n"
        "                          arms; default all ones\n"
        "      --tol <x>           convergence tolerance (mva, nc, mam, ag, fluid)\n"
        "      --iter_tol <x>      outer-loop tolerance (mva, nc, fluid)\n"
        "      --iter_max <n>      iteration cap (mva, nc, mam, ag, fluid)\n"
        "      --max-states <n>    truncation level of an OPEN agent's queue-length\n"
        "                          dimension (ag), options.config.maxStates;\n"
        "                          default 100. A closed class is bounded by its\n"
        "                          own population instead, and --method inapinf\n"
        "                          ignores the level and solves the open agents on\n"
        "                          the infinite state space\n"
        "      --fork-join <arm>   which fork-join transform the mean-value fixed\n"
        "                          point takes (mva, nc): default|mmt|fjt is the\n"
        "                          MMT transform, ht|heidelberger-trivedi the\n"
        "                          Heidelberger-Trivedi one, which is CLOSED\n"
        "                          models only and a different answer to the same\n"
        "                          model rather than a faster route to one\n"
        "      --cutoff <n|matrix> open jobs per class in the CTMC state space;\n"
        "                          a matrix is per (station,class), '1,1,0;3,3,0;0,0,3'\n"
        "                          (ctmc); without it the reference's\n"
        "                          ceil(6000^(1/(M*K))) is used and reported.\n"
        "                          Also the level truncation of mam -a prob,\n"
        "                          whose open queue has no bound of its own\n"
        "      --fj-accuracy <n>   FJ_codes truncation C of the queue-length\n"
        "                          difference between the two fork-join\n"
        "                          branches (mam, homogeneous fork-join);\n"
        "                          default 100, larger is more accurate\n"
        "      --fj-tmode <mode>   how that approximation solves for its T\n"
        "                          matrix: NARE (default) or Sylves\n"
        "      --timescale <mode>  auto (default), discrete or continuous: how\n"
        "                          -s mam reads the time scale. auto lets the\n"
        "                          distributions decide; discrete raises rather\n"
        "                          than solve a model that mixes lattice and\n"
        "                          non-lattice laws. The slot is --slotlength\n"
        "      --tspan <t0>:<t1>   transient horizon (ctmc -a tranprob,\n"
        "                          mam -a tran, fluid); a bare <t1> starts at 0.\n"
        "                          For the CTMC and MAM there is no default:\n"
        "                          pi(t) on an unstated horizon is not a\n"
        "                          quantity. For the fluid solver it bounds the\n"
        "                          integration, and is what -s fluid --method kp\n"
        "                          reports its covariance AT\n"
        "  -n, --node <n>          1-based stateful node a state query is\n"
        "                          labelled by (ctmc -a tranprob, -a sample and\n"
        "                          ssa -a sample), the\n"
        "                          queue mam -a prob reports, or the station\n"
        "                          mva -a marg reports; without it the\n"
        "                          whole network is reported (the MAM queries\n"
        "                          take the model's only Queue)\n"
        "  -c, --class <r>         1-based job class of mva -a marg; without it\n"
        "                          every class is reported\n"
        "                          NOTE ON THE INDEX BASE: -n and -c are 1-BASED\n"
        "                          here, as every station index this CLI takes\n"
        "                          is, and 0-BASED in jline.cli.LineCLI, which\n"
        "                          indexes as Java does. The short spellings are\n"
        "                          accepted so one command line parses in both,\n"
        "                          but the SAME number names a different node --\n"
        "                          the two bridges (cpp_dispatch, jar_dispatch)\n"
        "                          each convert for their own CLI\n"
        "      --marg-states <ns>  comma-separated job counts the mva -a marg\n"
        "                          curve is evaluated at, the reference's\n"
        "                          state_m; without it each law takes its own\n"
        "                          default range (0..N_r closed, mean + 5 sigma\n"
        "                          Poisson, the 1e-10 tail geometric)\n"
        "      --notation <form>   scalar (default) | matrix, the form the ODE\n"
        "                          export writes (fluid -a odes only)\n"
        "      --symbolic <b>      computer-algebra backend of fluid -a\n"
        "                          jacobian: auto (default, searches for a\n"
        "                          line-sage-rest service), a URL, an image\n"
        "                          name, or none to differentiate locally\n"
        "      --equilibria        also ask the backend for the solutions of\n"
        "                          f(x) = 0 (fluid -a jacobian). Needs a\n"
        "                          backend: solving is not differentiating\n"
        "      --cdf-algorithm <a> exact (default, pfqn_stdf) | rd\n"
        "                          (pfqn_stdf_heur), how the sojourn law is\n"
        "                          inverted (nc -a cdf only)\n"
        "      --perm-engine <e>   exact (default, Ryser) | spm | bethe | heur |\n"
        "                          huberlaw | adapart, the permanent estimator\n"
        "                          of nc -a sysmarg. The five approximations\n"
        "                          refuse a demand matrix with a zero entry;\n"
        "                          spm is the saddle point, whose cost does not\n"
        "                          grow with the population\n"
        "      --tran-points <n>   points on the uniform transient grid the ENV\n"
        "                          mean-field coupling sums its stage exit\n"
        "                          metrics over (env only); default 1001\n"
        "      --state <n,...>     the ENCODED state row of the node named by\n"
        "                          --node that -a prob asks about, i.e.\n"
        "                          getProb(node, state)'s second argument; without\n"
        "                          it the query is about the model's default\n"
        "                          initial state\n"
        "      --events <n>        length of ONE sampled trajectory (-a sample);\n"
        "                          default 1000. NOT --samples, which is a\n"
        "                          solver's run length\n"
        "      --timestep <dt>     fixed output step of a transient CTMC solve\n"
        "      --transient-method <m>  ode (default) or fau, how the CTMC forward\n"
        "                          equation is advanced over the output grid\n"
        "      --fau-epsilon <e>   fau: probability mass the grid may discard\n"
        "      --fau-delta <d>     fau: occupancy below which a state is dropped\n"
        "                          (-a tranprob, -a tranreward), options.timestep\n"
        "                          of ctmc_transient.m; without it the grid is the\n"
        "                          integrator's own adaptive one. It changes WHERE\n"
        "                          the solution is reported, not how it is\n"
        "                          computed: the grid points are read off the same\n"
        "                          interpolant\n"
        "      --percentiles <p,..>  levels getPerctRespT is read at (-s mam),\n"
        "                          as fractions (0.9) or percents (90); default\n"
        "                          0.50,0.90,0.95,0.99, the reference's pers_stored\n"
        "      --reward-name <nm>  which declared reward -a reward-value returns\n"
        "                          the value function of. REQUIRED there and never\n"
        "                          defaulted: two rewards have different value\n"
        "                          functions, and picking one would mislabel it\n"
        "  -p, --port <n>          run as a solve SERVER on this port, speaking\n"
        "                          LineWebSocketServer's protocol: one WebSocket\n"
        "                          text message per connection, its first line the\n"
        "                          comma-separated argument list and its remainder\n"
        "                          the model document; the CLI's output comes back\n"
        "                          as one text message. Plaintext, one connection\n"
        "                          at a time, and -f is refused beside it\n"
        "  -m, --maxreq <n>        quit after serving n requests; without it the\n"
        "                          server runs until interrupted\n"
        "  -h, --help              the short message: the flags a first run needs\n"
        "      --help-all          this message, every option solver by solver\n"
        "  -V, --version           version string\n"
        "      --install           environment check: report which optional backends\n"
        "                          (Java/JMT, LQNS, qnsolver, SageMath) are reachable\n"
        "\n"
        "Layered models (-i lqnx|xml) additionally take:\n"
        "      --layer-solver <s>  solver run in each layer: mva (default)|nc|\n"
        "                          fluid|ssa, the C++ spelling of the reference's\n"
        "                          factory argument: LN(model, @(m) MVA(m))\n"
        "                          against LN(model, @(m) Fluid(m)). They converge\n"
        "                          to DIFFERENT fixed points, not to the same one\n"
        "                          by different routes, because each layer's\n"
        "                          results feed the next outer iteration's demands.\n"
        "                          fluid is double only and refuses a layer with a\n"
        "                          fork; ssa is NOISY, so the outer loop switches\n"
        "                          to the Robbins-Monro / Polyak-Ruppert controller\n"
        "      --method <name>     the LN update: default | moment3 | mwba.upper |\n"
        "                          mwba.lower. moment3 fits an APH to each layer's\n"
        "                          response-time CDF and convolves along the entry,\n"
        "                          which is what makes -a cdf possible; mwba.*\n"
        "                          reports Majumdar-Woodside box BOUNDS on\n"
        "                          throughput and processor utilization and solves\n"
        "                          no layer at all (every other metric is NaN)\n"
        "      --ln-transient <m>  coupled (default) | decoupled, how -a tran\n"
        "                          couples the layers. decoupled freezes the\n"
        "                          inter-layer demands at the fixed point; coupled\n"
        "                          relaxes time-varying demands through the fluid\n"
        "                          rate schedule until the trajectories settle\n"
        "      --ln-transient-channels <c>  both (default) | thinkt | callservt,\n"
        "                          which inter-layer coupling the relaxation\n"
        "                          injects, for isolating one channel's share\n"
        "      --sens-method <m>   auto (default) | exact | fd, the branch each\n"
        "                          LAYER's sensitivity table takes (-a sens);\n"
        "                          the same three under -s nc -a sens\n"
        "      --sens-scheme <s>   forward (default) | central, the difference\n"
        "                          quotient of the fd branch\n"
        "      --sens-step <h>     relative rate perturbation of the fd branch,\n"
        "                          in (0,1); default 1e-4, or 1e-2 for ssa layers\n"
        "      --no-interlocking   disable the interlocking correction\n"
        "      --repeat <k>        re-solve k times, report the best wall clock\n"
        "      -o layers           dump every layer's stations, classes, rates\n"
        "                          and routing instead of the AvgTable\n"
        "  -a takes avg (getAvgTable), tran (getTranAvg, needs --tspan and fluid\n"
        "  layers), sens (getSensitivityTable) and cdf (getCdfRespT, moment3).\n"
        "  --iter_max and --iter_tol set the outer LN loop; --tol does not apply,\n"
        "  and neither does any Network solver token.\n"
        "\n"
        "The external LQNS binary (-s lqns) additionally takes:\n"
        "      --method <name>     default | lqns | srvn | exactmva |\n"
        "                          srvn.exactmva | sim | lqsim | lqnsdefault.\n"
        "                          sim and lqsim run lqsim, the SIMULATOR;\n"
        "                          lqnsdefault is lqns with no pragma at all,\n"
        "                          which is a different fixed point and not a\n"
        "                          synonym for default\n"
        "      --multiserver <p>   conway|rolia|zhou|suri|reiser|schmidt|default\n"
        "                          (= rolia), the -Pmultiserver= pragma. Not\n"
        "                          passed to lqsim, which has no MVA to configure\n"
        "      --samples <n>       lqsim run length (-A); default 10000\n"
        "      --timeout <s>       kill the child after s seconds; without it the\n"
        "                          wrapper waits\n"
        "      --keep              keep the working directory with model.lqnx and\n"
        "                          model.lqxo instead of removing it\n"
        "      --verbose           echo the command line and the binary's output\n"
        "      --remote[-url <u>]  solve on a host running lqns-rest instead of\n"
        "                          locally; -url implies --remote. LINE ships no\n"
        "                          LQNS binary, so this is the other way to reach\n"
        "                          one\n"
        "  -a takes avg only: lqns computes no transient, no sensitivity and no\n"
        "  response-time distribution. QLen is the element utilization, Util its\n"
        "  processor utilization per server, RespT its phase-1 service time;\n"
        "  ResidT and ArvR print NaN because lqns reports neither.\n"
        "\n"
        "The external qnsolver binary (-s qns) additionally takes:\n"
        "      --method <name>     default (= rolia) | conway | rolia | zhou |\n"
        "                          reiser. The multiserver approximation, passed\n"
        "                          as qnsolver -m and only when the model HAS a\n"
        "                          multiserver station. suri and schmidt are\n"
        "                          listed by the reference but reach the tool\n"
        "                          through its SolverLQNS branch, which this port\n"
        "                          does not carry, so they are refused by name\n"
        "      --multiserver <p>   the same choice under its config spelling;\n"
        "                          --method wins when it names one\n"
        "      --timeout <s>       kill the child after s seconds; without it the\n"
        "                          wrapper waits\n"
        "      --keep              keep the working directory with model.jmva and\n"
        "                          result.jmva instead of removing it\n"
        "  -a takes avg only. The model is marshalled to the JMVA interchange\n"
        "  format at CHAIN level and the chain results are de-aggregated back to\n"
        "  classes, so only Queue, Delay and Source stations are expressible; any\n"
        "  other station is refused rather than dropped. A closed model that is\n"
        "  NOT product-form is refused too: the reference converts it with QN2LQN\n"
        "  and delegates to SolverLQNS, and QN2LQN is not ported. --arith is\n"
        "  double only, since the results arrive as the text an external binary\n"
        "  printed.\n"
        "\n"
        "Discrete-event simulation (-s ldes) additionally takes:\n"
        "      --ldes-tranfilter <f>  warmup filter: mser5 (default), fixed, none\n"
        "      --ldes-warmupfrac <x>  fraction the fixed filter discards (0.2)\n"
        "      --ldes-cimethod <m>    CI estimator: obm (default), bm, spectral,\n"
        "                             none\n"
        "      --ldes-cnvgon          stop on relative precision instead of on\n"
        "                             the --samples budget\n"
        "      --ldes-cnvgtol <x>     that precision target (0.05); implies\n"
        "                             --ldes-cnvgon\n"
        "      --slotted              run the analytical solver on a discrete\n"
        "                             (slotted) time scale; SolverNC routes to the\n"
        "                             discrete-time product form and refuses a model\n"
        "                             outside it\n"
        "      --slotlength <x>       the slot in model time units; implies --slotted\n"
        "      --ldes-slotted         run on a discrete (slotted) time scale; a\n"
        "                             sample off the lattice is an error, never\n"
        "                             rounded\n"
        "      --ldes-slotlength <x>  the slot (1.0); implies --ldes-slotted\n"
        "      --ldes-replications <n>  independent runs, averaged. A single path\n"
        "                             is NOT E[N](t): -a tran over an ensemble\n"
        "                             mean needs this\n"
        "      --ldes-numthreads <n>  workers for those replications\n"
        "      --ldes-maxtime <s>     wall-clock budget; the engine stops early\n"
        "                             and reports stopping=max_time\n"
        "      --ldes-initsol <v,..>  warm-start placement, station-major\n"
        "                             [st0_cl0, st0_cl1, ...]. Add --ldes-tranfilter\n"
        "                             fixed --ldes-warmupfrac 0 to reproduce\n"
        "                             initFromSolver, which assumes the placement\n"
        "                             is already a steady state\n"
        "      --ldes-rest-url <u>    solve on an LDES REST server instead of a\n"
        "                             local binary; same wire format, same numbers\n"
        "  The model.json is forwarded to the engine BYTE FOR BYTE, so a model\n"
        "  this port cannot itself parse (a cache with retrieval, an SPN, a\n"
        "  polling server) is simulated exactly as the MATLAB and Python clients\n"
        "  simulate it. -a avg, tran, cdf (the empirical response-time law),\n"
        "  sample and reward are served; -a prob is refused, because getProbSys\n"
        "  weighs the trajectory against the model's current state and the C++\n"
        "  NetworkStruct carries no such row. --arith is double only.\n"
        "\n"
        "Uncertainty quantification (-s uq) additionally takes:\n"
        "      --uq-solver <s>     the engine run at each design point: mva, nc,\n"
        "                          mam, ba, ctmc, fluid or ssa. REQUIRED: UQ\n"
        "                          computes nothing itself, and defaulting it\n"
        "                          would attribute the numbers to an engine the\n"
        "                          caller never chose\n"
        "  A model whose service or arrival process is a Prior is a FAMILY of\n"
        "  models. UQ discretizes each Prior, solves the tensor product of the\n"
        "  alternatives, and reports the prior-weighted expectation. Under -s uq\n"
        "  three flags describe the DESIGN and not the engine: --method is\n"
        "  quadrature (default, and the alias of discrete) or montecarlo,\n"
        "  --samples the nodes per continuous Prior (11), --seed the Monte Carlo\n"
        "  stream. The stage solver therefore keeps its own sample count and its\n"
        "  own seed; --tol, --iter_tol, --iter_max and --cutoff pass through to\n"
        "  it, since UQ has no convergence of its own. -a posterior prints every\n"
        "  design point, its weight and the means it substituted, which is what\n"
        "  says whether the expectation averaged two nearby models or two very\n"
        "  different ones. Every other solver REFUSES a model carrying a Prior\n"
        "  rather than lowering it to its mixture moments.\n"
        "  -a interval answers the OTHER epistemic question, in which a\n"
        "  parameter is bounded but not distributed: it drops the weights and\n"
        "  keeps the endpoints. On a single-class closed model of LI\n"
        "  single-server queues and delays it is the EXACT hull of MVA over the\n"
        "  demand box (2*(m+2) MVA calls, no design solved at all); otherwise it\n"
        "  falls back to the range over the solved design points, which for a\n"
        "  continuous Prior lies strictly inside the true range. The table says\n"
        "  which, and the fallback warns on stderr: a range that is not an\n"
        "  enclosure must not read like one.\n"
        "\n"
        "Random environments (-s env) read an Environment model.json, whose\n"
        "  stages each hold a Network and whose transitions carry the stage\n"
        "  holding times. The stages are solved TRANSIENTLY and coupled: each\n"
        "  stage starts from the queue lengths the previous one left, and the\n"
        "  reported means are the per-stage sojourn averages blended by the\n"
        "  environment probabilities. --method selects the coupling: meanfield\n"
        "  (default, the reference's, carries the marginal means across a\n"
        "  switch) or statevec|blend (carries the whole joint distribution).\n"
        "  meanfield solves each stage with the fluid transient and is double\n"
        "  only; statevec uniformizes a CTMC and takes the whole --arith ladder.\n"
        "  --method avg|dec asks instead for a closed-form limit, which carries\n"
        "  nothing across a switch and iterates nothing: avg solves ONE model\n"
        "  whose modulated rates are their probEnv-weighted averages (exact as\n"
        "  the environment gets fast), dec solves each stage in steady state and\n"
        "  blends by probEnv (exact as it gets slow). A model with an\n"
        "  environment-declared node breakdown is read from the nodeFailures\n"
        "  block, in the expanded or the one-stage macro form.\n"
        "  --tspan bounds the stage horizon and --tran-points its grid; --iter_tol\n"
        "  and --iter_max drive the fixed point. RespT and ArvR print as nan\n"
        "  because ENV computes neither -- the reference returns them as NaN too,\n"
        "  and ResidT carries QLen/Tput.\n"
        "\n"
        "Additions specific to this port:\n"
        "      --arith <mode>      double (default) | exact | real:<digits>\n"
        "      --list-api          list the API functions ported so far\n"
        "      --api <name>        invoke one API function directly\n"
        "      --args <path>       JSON arguments for --api; stdin when omitted\n"
        "\n"
        "Solvers and what they honour. Every model-solving solver reads -a avg,\n"
        "-f/-i and --method; nothing else is wired, so tolerances, iteration\n"
        "caps, seeds and sample counts keep their SolverOptions defaults rather\n"
        "than being invented here. mva and nc additionally read -a prob -- mva\n"
        "fits a binomial to its own means, nc returns the exact product-form\n"
        "probability, so the two disagree by construction. mva also reads -a\n"
        "marg, @SolverMVA's getProbMarg: P(n jobs of class r at station i) for\n"
        "every (station, class) pair, narrowed by --node / --class and evaluated\n"
        "at --marg-states. A closed class takes the Schmidt binomial fitted to\n"
        "Q(i,r); an open one takes the station's exact BCMP marginal (Poisson at\n"
        "an infinite server, multinomial-geometric at a queue). Both solvers read\n"
        "-a normconst, getProbNormConstAggr: nc reports the constant its solve\n"
        "already formed, mva RE-ENTERS its analyzer at method='exact', since only\n"
        "the exact recursion carries a G -- a model whose branch forms none, an\n"
        "open or mixed one above all, reports nan, as the reference does. nc also\n"
        "reads -a cdf,\n"
        "@SolverNC's getCdfRespT: the whole response-time law per (station,\n"
        "class) on one logarithmic grid, FCFS stations only, with\n"
        "--cdf-algorithm exact (pfqn_stdf, the default) or rd (pfqn_stdf_heur).\n"
        "nc reads -a sens as well, @NetworkSolver's getSensitivityTable: the\n"
        "derivative of each row's means with respect to its own service rate,\n"
        "selected with --sens-method / --sens-scheme / --sens-step. NC is one of\n"
        "the two engines whose exact branch differentiates the product-form\n"
        "recursion analytically, so auto takes it wherever the model is in its\n"
        "scope (single-server queues plus delays, not mixed) and falls back to\n"
        "finite differences elsewhere. NOT the -s ctmc -a sens analysis, which\n"
        "is getSensitivityRanking, a ranking of rate perturbations and not a\n"
        "table of derivatives.\n"
        "-a node is @NetworkSolver's getAvgNodeTable and is served by mva, nc,\n"
        "mam, ba and ctmc, the model solvers whose runner returns the station\n"
        "AvgResult it is built from. It is a DIFFERENT INDEX SPACE from -a avg,\n"
        "not a relabelling: the AvgTable has one row per STATION, so a\n"
        "ClassSwitch, Router, Fork, Join or Sink never appears in it, yet jobs\n"
        "flow through all of them. QLen, Util, RespT and ResidT are the station\n"
        "numbers scattered to their node indices and zero elsewhere -- a node\n"
        "that is not a station holds no jobs -- while ArvR and Tput are\n"
        "recomputed for every node by sn_get_node_arvr_from_tput and\n"
        "sn_get_node_tput_from_tput. The reference's finite-capacity-region\n"
        "pseudo-node rows are NOT emitted: this port's AvgResult carries no\n"
        "per-region queue length or utilization to fill them with.\n"
        "\n"
        "auto picks the engine\n"
        "with chooseSolverHeur and prints the name it picked; a branch selecting\n"
        "JMT or LDES refuses by name rather than substituting another. ctmc reads\n"
        "--cutoff and ports -a avg, prob, gen, states, tranprob, sample, reward,\n"
        "cdf, first-passt and sens, which are @SolverCTMC's\n"
        "getProbSys/getProbSysAggr and the per-station getProb/getProbAggr,\n"
        "getInfGen, getStateSpace, getTranProbSysAggr, sampleSys, getAvgReward,\n"
        "getCdfRespT, getCdfFirstPassT (state sets via --passage-from and\n"
        "--passage-into, as 1-based space rows '3,5' or state rows '0,2;1,1'),\n"
        "first-passt-moments (-a firstpasstmom, the same two sets plus\n"
        "--passage-orders: exact moments by one linear solve per order, so a\n"
        "variance or a skewness costs no truncated curve)\n"
        "and getSensitivityRanking. Its --method also takes mdd, which holds the\n"
        "reachable set in a decision diagram and solves K coupled level-CTMCs\n"
        "instead of the |S|-state generator, and cftp / cftp.approx, which draw\n"
        "iid states from the exact stationary law by coupling from the past; all\n"
        "three serve -a avg only, having no explicit chain to answer the rest\n"
        "from, and the cftp rows carry Monte Carlo error. mam ports -a avg,\n"
        "prob, cdf, tran and internals, which are @SolverMAM's getProb and\n"
        "getProbMarg (the joint (level, phase) law of the queue and its\n"
        "per-class marginals, truncated at --cutoff when the model is open),\n"
        "getCdfRespT with getPerctRespT beside it, getTranAvg over --tspan (the\n"
        "reference forces method ldqbd there, so --method does not select it),\n"
        "and getMAMResult, the M/G/1-type internals of a single queue. fluid\n"
        "additionally reads -a tran, -a prob, -a cdf, -a var and -a aoi, which\n"
        "are @SolverFLD's getTranAvg (the metrics along the trajectory, over\n"
        "--tspan or, without one, the horizon the reference's own adaptive loop\n"
        "converges at), getProbAggr (a law FITTED to the fluid means, so it does\n"
        "not agree with the mva or nc answer by construction), getCdfRespT (the\n"
        "response-time law per station and class, read off a marked-fluid\n"
        "integration started from the steady state), getMoments/getTranAvgVar\n"
        "and getAvgAoI with getCdfAoI beside it -- the last needing method mfq\n"
        "and the Source/Queue/Sink topology the age laws are defined for. It\n"
        "further reads -a odes, which is\n"
        "@SolverFLD/exportODEs, and --notation for the form it writes, and -a\n"
        "jacobian, which is @SolverFLD/getJacobian: d f_i / d x_j of the drift,\n"
        "differentiated locally over the structure of the system, with the\n"
        "equilibria beside it under --equilibria, which needs the\n"
        "line-sage-rest backend --symbolic names. Only the smooth methods have\n"
        "a Jacobian: min(n,S) has none where the regime switches, so the\n"
        "min-scaled drifts are refused by the factor that carries the kink.\n"
        "nc, ba and ctmc run\n"
        "under every --arith backend, ctmc's cftp method excepted: its sampler\n"
        "works in the log domain, so it refuses 'exact' by name rather than\n"
        "answering in a field it does not live in. mdd does run under 'exact'\n"
        "(its level solve drops Householder for a rational least squares there),\n"
        "but the LEVEL AGGREGATION is still an approximation away from product\n"
        "form: exact arithmetic pins the fixed point, not the model. mam is\n"
        "double only (its phase-type fitter\n"
        "needs transcendental arithmetic) and so is ssa (its sample path is\n"
        "generated from exponential clocks) and fluid (LSODA). ba reports a\n"
        "BOUND, not an estimate, and ssa a simulation carrying Monte Carlo\n"
        "error; neither is comparable with an exact solver except as such.\n"
        "\n"
        "Arithmetic: 'exact' computes in arbitrary-precision rationals and\n"
        "reports numerator and denominator alongside the double value; 'real'\n"
        "computes in fixed high-precision binary floating point, at 50, 100 or\n"
        "200 digits (a request in between is rounded up to the next tier).\n"
        "\n"
        "--api arguments are a JSON object keyed by the MATLAB parameter names,\n"
        "e.g. {\"L\": [[0.6,0.4]], \"N\": [2,1], \"Z\": [1,0.5]}: a 2-D array is a\n"
        "matrix (row-major), a 1-D array a row vector, a bare number a scalar.\n"
        "A JSON number is read as its shortest round-tripping decimal, so 0.6 is\n"
        "3/5 in exact arithmetic; pass a string such as \"1/3\" for anything else.\n",
        kVersion);
}

/**
 * `--install`: the environment check, the C++ twin of MATLAB's `lineInstall`,
 * the JAR's `jline.cli.LineInstall` and Python's `line-install`.
 *
 * NOTHING IT LOOKS FOR IS REQUIRED. The C++ edition is header-only and its
 * native solvers stand alone, so every dependency probed here backs one
 * optional wrapper or backend. A miss is therefore a warning on stderr that
 * names the solvers it disables and how to install it, never an error: the
 * point of the command is to tell a fresh checkout which solvers it can
 * actually reach, not to refuse to run.
 *
 * @return true when nothing warned
 */
bool install_check() {
    bool has_warnings = false;
    // stdout is block-buffered when the check is piped or redirected while
    // stderr is not, so an unflushed progress line would surface AFTER the
    // warning it belongs to. Flush before every warning so the transcript reads
    // in the order the checks ran.
    const auto warn = [&](const std::string& text) {
        std::fflush(stdout);
        std::fprintf(stderr, "%s\n", text.c_str());
        std::fflush(stderr);
        has_warnings = true;
    };

    std::printf("Checking LINE (C++)...\n");
    std::printf("  line-cli %s\n", kVersion);

    std::printf("Checking Java runtime (JMT wrapper)...\n");
    const std::string java = line::jmt::detail::find_java();
    if (java.empty()) {
        warn("WARNING: no Java runtime was found in LINE_JAVA, JAVA_HOME or on PATH, so the "
             "JMT wrapper (-s jmt) cannot run. Install a JRE 8 or later.");
    } else {
        std::printf("  %s\n", java.c_str());
    }

    std::printf("Checking JMT...\n");
    const std::string jmt_dir = line::jmt::detail::jmt_jar_path();
    if (jmt_dir.empty()) {
        warn("WARNING: JMT.jar was not found, this is required by the JMT wrapper. Download it "
             "from https://line-solver.sourceforge.net/latest/JMT.jar into common/, or point "
             "LINE_JMT_DIR at the folder holding it.");
    } else {
        std::printf("  %s/JMT.jar\n", jmt_dir.c_str());
    }

    std::printf("Checking LQNS...\n");
    if (!line::lqns::lqns_is_available()) {
        const std::string banner = line::lqns::lqns_version();
        if (banner.empty())
            warn("WARNING: lqns is not installed, this is required by the LQNS wrapper (-s lqns) "
                 "for layered models. Download it at: https://github.com/layeredqueuing/dist");
        else
            warn("WARNING: the installed lqns is too old for LINE, which needs release 6 or "
                 "later; it reports '" + banner + "'. Upgrade it from: "
                 "https://github.com/layeredqueuing/dist");
    } else {
        std::printf("  %s\n", line::lqns::lqns_version().c_str());
    }

    std::printf("Checking QNS (qnsolver)...\n");
    if (!line::qns::is_available())
        warn("WARNING: qnsolver is not installed, this is required by the QNS wrapper (-s qns). "
             "It ships with LQNS: https://github.com/layeredqueuing/dist");

    std::printf("Checking symbolic backend (line-sage-rest)...\n");
    if (!line::io::docker_daemon_available())
        warn(std::string("WARNING: Docker is not available, so the SageMath symbolic backend "
                         "cannot start. It is the only computer algebra system this edition "
                         "reaches and is required by the symbolic methods of "
                         "SolverCTMC/SolverFluid. Install Docker, then run: "
                         "docker run -d -p 8080:8080 ") + line::sym::SYM_DOCKER_IMAGE);
    else if (line::sym::sym_find_image().empty())
        warn(std::string("WARNING: the line-sage-rest image is not present locally, this may be "
                         "required by some LINE methods. Pull it with: "
                         "docker run -d -p 8080:8080 ") + line::sym::SYM_DOCKER_IMAGE);

    if (has_warnings)
        std::printf("Completed. LINE has warnings.\n");
    else
        std::printf("Success. LINE is ready to use.\n");
    return !has_warnings;
}

void list_api() {
    const auto& reg = line::api_registry();
    std::printf("%-24s %-8s %-24s %s\n", "function", "domain", "arithmetic", "ported from");
    for (const auto& e : reg) {
        std::string modes;
        for (std::size_t k = 0; k < e.arith.size(); ++k) {
            if (k) modes += ",";
            modes += line::arith_name(e.arith[k]);
        }
        std::printf("%-24s %-8s %-24s %s\n", e.name.c_str(), e.domain.c_str(), modes.c_str(),
                    e.reference.c_str());
    }
    std::printf("\n%zu of ~480 API functions ported.\n", reg.size());
}

/**
 * Read the --api argument object: from the file named by --args, or from stdin
 * when --args is absent. A parse failure names the source and the position, so
 * a malformed file is a legible error rather than an empty argument set.
 */
line::reg::Json read_api_args(const std::string& path) {
    std::string text;
    std::string source;
    if (path.empty()) {
        source = "standard input";
        text.assign(std::istreambuf_iterator<char>(std::cin), std::istreambuf_iterator<char>());
    } else {
        source = "'" + path + "'";
        std::ifstream in(path.c_str());
        if (!in) throw line::InputError("cannot open the --args file " + source);
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    // No --args and nothing on stdin means the caller forgot the arguments.
    if (text.find_first_not_of(" \t\r\n") == std::string::npos)
        throw line::InputError("no --api arguments given: pass --args <path> or a JSON object on "
                               "standard input (read from " +
                               source + ")");
    try {
        return line::reg::Json::parse(text);
    } catch (const line::reg::Json::parse_error& e) {
        throw line::InputError("malformed --api arguments in " + source + ": " + e.what());
    }
}

/**
 * The JAR CLI's `-a` spelling, mapped onto this port's token.
 *
 * `jline.cli.LineCLI` names its analyses with hyphenated, getter-shaped tokens
 * (`cdf-respt`, `prob-sys-aggr`, `tran-avg`) and this port names them by the
 * question (`cdf`, `prob`, `tran`), because one arm here answers what the JAR
 * splits across several: `-a prob` reports `getProbSys`, `getProbSysAggr` and
 * the per-station `getProb`/`getProbAggr` in ONE table, and `-a sample` walks
 * `sampleSys`, `sampleSysAggr` and the per-node pair in one trajectory. Both
 * spellings are therefore accepted and the JAR's collapse onto the arm that
 * already contains the answer -- a script written against the JAR CLI keeps
 * working, and nothing here is renamed to make that true.
 *
 * AN UNKNOWN METHOD NAME IS RETURNED UNCHANGED, not rejected here: the per-solver
 * whitelists downstream refuse by name and say which analyses that solver
 * serves, which is a better message than a table of every token in the CLI.
 */
std::string normalize_analysis(const std::string& a) {
    // Same question, different spelling.
    if (a == "cdf-respt" || a == "cdfrespt") return "cdf";
    if (a == "cdf-passt" || a == "cdfpasst") return "cdfpasst";
    if (a == "first-passt" || a == "cdf-firstpasst" || a == "cdffirstpasst") return "firstpasst";
    if (a == "first-passt-moments" || a == "firstpasst-moments" || a == "firstpasstmoments")
        return "firstpasstmom";
    if (a == "perct-respt" || a == "perctrespt") return "perct";
    if (a == "tran-avg" || a == "tranavg") return "tran";
    if (a == "tran-cdf-respt" || a == "trancdfrespt") return "trancdf";
    if (a == "tran-cdf-passt" || a == "trancdfpasst") return "trancdfpasst";
    if (a == "tran-prob" || a == "tranprob-sys-aggr") return "tranprob";
    if (a == "generator") return "gen";
    if (a == "state-space" || a == "statespace") return "states";
    if (a == "reward-steady" || a == "rewardsteady") return "reward";
    if (a == "reward-value" || a == "rewardvalue") return "rewardvalue";
    if (a == "node-chain" || a == "node-chain-table") return "nodechain";
    // The JAR's four probability getters and its four samplers, each answered
    // whole by one arm here. Collapsing them is not a loss: the arm emits every
    // one of the four, so a caller asking for the aggregate receives it beside
    // the joint rather than instead of it.
    if (a == "prob-aggr" || a == "prob-sys" || a == "prob-sys-aggr") return "prob";
    if (a == "prob-marg" || a == "probmarg") return "marg";
    if (a == "prob-sys-marg" || a == "probsysmarg" || a == "sys-marg") return "sysmarg";
    if (a == "sample-aggr" || a == "sample-sys" || a == "sample-sys-aggr") return "sample";
    // `getStageTable` IS `getAvgTable` on a Network model, in the reference and
    // in the JAR both: a network has one implicit stage, and only an
    // Environment has several. Mapped here rather than given an arm of its own,
    // because an arm would be a second name for one table and free to drift
    // from it.
    if (a == "stage") return "avg";
    return a;
}

/** `-a` split on commas, each token normalized; never empty. */
std::vector<std::string> analysis_list(const std::string& spec) {
    std::vector<std::string> out;
    std::size_t at = 0;
    while (at <= spec.size()) {
        const std::size_t comma = spec.find(',', at);
        std::string tok =
            spec.substr(at, comma == std::string::npos ? std::string::npos : comma - at);
        // A stray space around a comma is a typo, not a different analysis.
        while (!tok.empty() && std::isspace(static_cast<unsigned char>(tok.front())))
            tok.erase(tok.begin());
        while (!tok.empty() && std::isspace(static_cast<unsigned char>(tok.back())))
            tok.pop_back();
        if (tok.empty())
            throw line::InputError("-a takes a comma-separated list of analyses and one entry of '" +
                                   spec + "' is empty");
        // `all` is the JAR's composite of the station table and the system one,
        // expanded HERE rather than inside a solver arm so every downstream
        // whitelist sees the two analyses it already knows.
        if (normalize_analysis(tok) == "all") {
            out.push_back("avg");
            out.push_back("sys");
        } else {
            out.push_back(normalize_analysis(tok));
        }
        if (comma == std::string::npos) break;
        at = comma + 1;
    }
    if (out.empty()) throw line::InputError("-a takes at least one analysis");
    return out;
}

struct Options {
    std::string file, input = "json", output = "readable", solver = "auto", analysis = "avg";
    std::string arith = "double", api, args;
    bool help = false, help_all = false, version = false, list = false;
    /**
     * `--find-solver [metric]`: report which solvers and methods can analyze the
     * model named by -f, and exit without solving it.
     *
     * `find_solver_all` is `--find-solver-all`, which keeps the refused pairs
     * and the reason each was refused. The report is arithmetic-independent --
     * it asks feature sets and shapes, not numbers -- so it always reads the
     * model at double and ignores --arith.
     */
    bool find_solver = false, find_solver_all = false;
    std::string find_solver_metric;
    /** `--install`: run the environment check and exit, solving nothing. */
    bool install = false;
    /**
     * Whether -i was actually passed.
     *
     * Without it a `.lqnx` path could not be recognised: `input` defaults to
     * json, and a defaulted json is indistinguishable from an explicit one, so
     * the extension sniff below would either never fire or would override a
     * caller who said `-i json` deliberately.
     */
    bool input_given = false;
    /**
     * `-p/--port` and `-m/--maxreq`: server mode.
     *
     * `port == 0` is "not given" and not "port 0": binding port 0 asks the
     * kernel for an ephemeral one, which a caller who typed no port did not
     * ask for. `maxreq == 0` is unbounded, matching the JAR's documented
     * "quit after this many requests" with no cap by default.
     */
    int port = 0;
    int maxreq = 0;
    Knobs knobs;
};

/** Whether `file` ends in one of JMT's three simulation-document extensions. */
bool has_jsim_extension(const std::string& file) {
    const std::string::size_type dot = file.find_last_of('.');
    if (dot == std::string::npos) return false;
    std::string ext = file.substr(dot + 1);
    for (std::size_t i = 0; i < ext.size(); ++i)
        ext[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(ext[i])));
    return ext == "jsim" || ext == "jsimg" || ext == "jsimw";
}

/** Whether `file` ends in the PNML extension. */
bool has_pnml_extension(const std::string& file) {
    const std::string::size_type dot = file.find_last_of('.');
    if (dot == std::string::npos) return false;
    std::string ext = file.substr(dot + 1);
    for (std::size_t i = 0; i < ext.size(); ++i)
        ext[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(ext[i])));
    return ext == "pnml";
}

/** Whether `file` ends in one of the layered model's two extensions. */
bool has_lqn_extension(const std::string& file) {
    const std::string::size_type dot = file.find_last_of('.');
    if (dot == std::string::npos) return false;
    std::string ext = file.substr(dot + 1);
    for (std::size_t i = 0; i < ext.size(); ++i)
        ext[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(ext[i])));
    return ext == "lqnx" || ext == "xml";
}

Options parse_args(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* what) -> std::string {
            if (i + 1 >= argc) throw line::InputError(std::string("missing value after ") + what);
            return argv[++i];
        };
        if (a == "-h" || a == "--help") o.help = true;
        else if (a == "--help-all" || a == "--help-full") o.help_all = true;
        else if (a == "-V" || a == "--version") o.version = true;
        else if (a == "--install") o.install = true;
        else if (a == "--list-api") o.list = true;
        else if (a == "--find-solver" || a == "--find-method" || a == "--help-model") {
            o.find_solver = true;
            // The metric is OPTIONAL, so it is taken only when the next token is
            // not itself a flag: `--find-solver -f m.json` must not swallow -f.
            if (i + 1 < argc && argv[i + 1][0] != '-') o.find_solver_metric = argv[++i];
        } else if (a == "--find-solver-all" || a == "--find-method-all") {
            o.find_solver = true;
            o.find_solver_all = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') o.find_solver_metric = argv[++i];
        }
        else if (a == "-f" || a == "--file") o.file = next("-f");
        else if (a == "-i" || a == "--input") { o.input = next("-i"); o.input_given = true; }
        else if (a == "-o" || a == "--output") o.output = next("-o");
        else if (a == "-s" || a == "--solver") o.solver = next("-s");
        else if (a == "-a" || a == "--analysis") o.analysis = next("-a");
        else if (a == "--arith") o.arith = next("--arith");
        else if (a == "-p" || a == "--port") {
            const std::string v = next("-p");
            const long n = std::atol(v.c_str());
            if (n < 1 || n > 65535)
                throw line::InputError("-p takes a TCP port in 1..65535 (got '" + v + "')");
            o.port = static_cast<int>(n);
        } else if (a == "-m" || a == "--maxreq") {
            const std::string v = next("-m");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError(
                    "-m is the number of requests the server serves before quitting and must be "
                    "positive; omit it to serve indefinitely (got '" + v + "')");
            o.maxreq = static_cast<int>(n);
        }
        else if (a == "--api") o.api = next("--api");
        else if (a == "--args") o.args = next("--args");
        else if (a == "--method") o.knobs.method = next("--method");
        else if (a == "--qrf-params") o.knobs.qrf_params = next("--qrf-params");
        else if (a == "--qrf-alpha") o.knobs.qrf_alpha = next("--qrf-alpha");
        else if (a == "--level") {
            const std::string v = next("--level");
            const int lv = std::atoi(v.c_str());
            if (lv < 1) throw line::InputError("--level must be a positive integer (got '" + v + "')");
            o.knobs.level = lv;
        }
        else if (a == "--samples") {
            const std::string v = next("--samples");
            const double d = std::atof(v.c_str());  // accepts 1e6 as well as 1000000
            if (!(d >= 1.0))
                throw line::InputError("--samples must be a positive count (got '" + v + "')");
            o.knobs.samples = static_cast<std::size_t>(d);
        } else if (a == "-d" || a == "--seed") {
            const std::string v = next("--seed");
            o.knobs.seed = std::strtoul(v.c_str(), nullptr, 10);
            if (o.knobs.seed == 0)
                throw line::InputError("--seed must be a positive integer (got '" + v + "')");
        } else if (a == "--warmupfrac") {
            const std::string v = next("--warmupfrac");
            const double f = std::atof(v.c_str());
            if (!(f >= 0.0 && f < 1.0))
                throw line::InputError(
                    "--warmupfrac is the fraction of the path discarded before the means are "
                    "taken and must lie in [0,1) (got '" + v + "')");
            o.knobs.warmupfrac = f;
        } else if (a == "--pstar") {
            const std::string v = next("--pstar");
            const double ps = std::atof(v.c_str());
            if (!(ps > 0.0))
                throw line::InputError(
                    "--pstar is the exponent of the fluid p-norm smoothing and must be positive "
                    "(got '" + v + "')");
            o.knobs.pstar = ps;
        } else if (a == "--busyperiod" || a == "--busyperiod-subnet") {
            // Same all-or-nothing parse as --marg-states: an entry dropped from
            // the list is a DIFFERENT report, not a shorter one.
            const bool orders = (a == "--busyperiod");
            const std::string v = next(orders ? "--busyperiod" : "--busyperiod-subnet");
            std::vector<std::size_t>& into = orders ? o.knobs.busy_orders : o.knobs.busy_subnet;
            std::size_t at = 0;
            while (at <= v.size()) {
                const std::size_t comma = v.find(',', at);
                const std::string tok =
                    v.substr(at, comma == std::string::npos ? std::string::npos : comma - at);
                if (tok.empty() || tok.find_first_not_of("0123456789") != std::string::npos ||
                    std::atol(tok.c_str()) < 1)
                    throw line::InputError(
                        std::string(orders ? "--busyperiod takes a comma-separated list of "
                                             "positive orders"
                                           : "--busyperiod-subnet takes a comma-separated list of "
                                             "1-based station indexes") +
                        " (got '" + v + "')");
                into.push_back(static_cast<std::size_t>(std::atol(tok.c_str())));
                if (comma == std::string::npos) break;
                at = comma + 1;
            }
        } else if (a == "--tol") o.knobs.tol = std::atof(next("--tol").c_str());
        else if (a == "--iter_tol") o.knobs.iter_tol = std::atof(next("--iter_tol").c_str());
        else if (a == "--iter_max") o.knobs.iter_max = std::atoi(next("--iter_max").c_str());
        else if (a == "--max-states") {
            const std::string v = next("--max-states");
            const long long n = std::atoll(v.c_str());
            if (n <= 0)
                throw line::InputError(
                    "--max-states truncates an open agent's queue-length dimension and takes a "
                    "positive state count (got '" + v + "')");
            o.knobs.max_states = n;
        }
        else if (a == "--multiserver") o.knobs.multiserver = next("--multiserver");
        else if (a == "--fork-join" || a == "--fork_join")
            o.knobs.fork_join = next("--fork-join");
        else if (a == "--tran-points" || a == "--tran_points") {
            const std::string v = next("--tran-points");
            const long n = std::atol(v.c_str());
            if (n < 2)
                throw line::InputError(
                    "--tran-points is the number of points on the transient grid and needs at "
                    "least two, a start and an end (got '" + v + "')");
            o.knobs.tran_points = static_cast<std::size_t>(n);
        }
        else if (a == "--mdd-tol" || a == "--mdd_tol") {
            const std::string v = next("--mdd-tol");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0))
                throw line::InputError("--mdd-tol must be a positive tolerance (got '" + v + "')");
            o.knobs.mdd_tol = d;
        } else if (a == "--mdd-maxiter" || a == "--mdd_maxiter") {
            const std::string v = next("--mdd-maxiter");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError(
                    "--mdd-maxiter must be a positive sweep count (got '" + v + "')");
            o.knobs.mdd_maxiter = static_cast<int>(n);
        }
        else if (a == "--fj-accuracy") {
            const std::string v = next("--fj-accuracy");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError(
                    "--fj-accuracy is the FJ_codes truncation C of the queue-length difference "
                    "between the two fork-join branches and must be at least 1 (got '" + v + "')");
            o.knobs.fj_accuracy = static_cast<int>(n);
        } else if (a == "--fj-tmode") {
            const std::string v = next("--fj-tmode");
            if (v != "NARE" && v != "Sylves")
                throw line::InputError(
                    "--fj-tmode selects how computeT.m solves for the T matrix and is 'NARE' (the "
                    "Riccati route, the default) or 'Sylves' (the fixed-point iteration); got '" +
                    v + "'");
            o.knobs.fj_tmode = v;
        } else if (a == "--timescale") {
            const std::string v = next("--timescale");
            if (v != "auto" && v != "discrete" && v != "continuous")
                throw line::InputError(
                    "--timescale decides whether the model is read on a slot lattice and is "
                    "'auto' (the default), 'discrete' or 'continuous'; got '" + v + "'");
            o.knobs.timescale = v;
        }
        else if (a == "--force") {
            o.knobs.force = true;
        }
        else if (a == "--cutoff") {
            const std::string v = next("--cutoff");
            if (v.find(',') != std::string::npos || v.find(';') != std::string::npos) {
                o.knobs.cutoff_mat = parse_cutoff_matrix(v);
                if (o.knobs.cutoff_mat.empty())
                    throw line::InputError(
                        "--cutoff takes a number or a per-(station,class) matrix written "
                        "'r1c1,r1c2;r2c1,r2c2' (got '" + v + "')");
            } else {
                const double d = std::atof(v.c_str());
                if (!(d >= 1.0))
                    throw line::InputError(
                        "--cutoff must be a positive job count per open class (got '" + v + "')");
                o.knobs.cutoff = d;
            }
        } else if (a == "--tspan" || a == "--timespan") {
            const std::string v = next("--tspan");
            // BOTH SEPARATORS. This port has always written the horizon
            // `t0:t1`, and `jline.cli.LineCLI --timespan` has always written it
            // `t0,t1`; accepting only one made a command line that names a
            // horizon unportable between the two CLIs even after the flag names
            // were reconciled.
            std::string::size_type sep = v.find(':');
            if (sep == std::string::npos) sep = v.find(',');
            // A bare value is the END of the horizon and starts at 0, which is
            // what a transient from the initial state means; `t0:t1` states both.
            const double lo = sep == std::string::npos ? 0.0 : std::atof(v.substr(0, sep).c_str());
            const double hi = std::atof(
                (sep == std::string::npos ? v : v.substr(sep + 1)).c_str());
            // An INFINITE horizon is refused here rather than integrated to:
            // pi(t) on [0, Inf) is the stationary vector, which -a avg reports.
            if (!(hi > lo) || !(lo >= 0.0) || !std::isfinite(hi))
                throw line::InputError(
                    "--tspan must be a finite horizon 0 <= t0 < t1, given as <t1>, <t0>:<t1> or "
                    "<t0>,<t1> (got '" + v + "')");
            o.knobs.t0 = lo;
            o.knobs.t1 = hi;
        } else if (a == "-n" || a == "--node") {
            const std::string v = next("--node");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError("--node must be a positive 1-based node index (got '" + v +
                                       "')");
            o.knobs.node = static_cast<std::size_t>(n);
        } else if (a == "-c" || a == "--class") {
            const std::string v = next("--class");
            const long c = std::atol(v.c_str());
            if (c < 1)
                throw line::InputError("--class must be a positive 1-based class index (got '" + v +
                                       "')");
            o.knobs.jobclass = static_cast<std::size_t>(c);
        } else if (a == "--marg-states" || a == "--marg_states") {
            // Every entry must parse: a dropped one shortens the curve silently.
            const std::string v = next("--marg-states");
            std::size_t at = 0;
            while (at <= v.size()) {
                const std::size_t comma = v.find(',', at);
                const std::string tok =
                    v.substr(at, comma == std::string::npos ? std::string::npos : comma - at);
                if (tok.empty() || tok.find_first_not_of("0123456789") != std::string::npos)
                    throw line::InputError(
                        "--marg-states takes a comma-separated list of non-negative job counts "
                        "(got '" + v + "')");
                o.knobs.marg_states.push_back(std::atol(tok.c_str()));
                if (comma == std::string::npos) break;
                at = comma + 1;
            }
        } else if (a == "--state") {
            // Same all-or-nothing discipline as --marg-states: a state vector
            // with one entry dropped is a DIFFERENT state, not a shorter one,
            // and the length is checked against the node's own space later.
            const std::string v = next("--state");
            std::size_t at = 0;
            while (at <= v.size()) {
                const std::size_t comma = v.find(',', at);
                const std::string tok =
                    v.substr(at, comma == std::string::npos ? std::string::npos : comma - at);
                if (tok.empty() || tok.find_first_not_of("0123456789") != std::string::npos)
                    throw line::InputError(
                        "--state is the state vector -a prob asks about and takes a "
                        "comma-separated list of non-negative counts (got '" + v + "')");
                o.knobs.state.push_back(std::atol(tok.c_str()));
                if (comma == std::string::npos) break;
                at = comma + 1;
            }
        } else if (a == "--events") {
            const std::string v = next("--events");
            const double d = std::atof(v.c_str());  // accepts 5e3 as well as 5000
            if (!(d >= 1.0))
                throw line::InputError(
                    "--events is the length of a sampled trajectory and must be a positive event "
                    "count (got '" + v + "')");
            o.knobs.events = static_cast<std::size_t>(d);
        } else if (a == "--timestep") {
            const std::string v = next("--timestep");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0) || !std::isfinite(d))
                throw line::InputError(
                    "--timestep is the fixed output step of a transient analysis and must be a "
                    "positive finite time (got '" + v + "')");
            o.knobs.timestep = d;
        } else if (a == "--percentiles") {
            const std::string v = next("--percentiles");
            std::size_t at = 0;
            while (at <= v.size()) {
                const std::size_t comma = v.find(',', at);
                const std::string tok =
                    v.substr(at, comma == std::string::npos ? std::string::npos : comma - at);
                if (tok.empty())
                    throw line::InputError(
                        "--percentiles takes a comma-separated list of levels (got '" + v + "')");
                double p = std::atof(tok.c_str());
                // A LEVEL ABOVE 1 IS A PERCENT, below it a probability. The JAR
                // documents `--percentiles 50,90,95,99` and MATLAB stores
                // `pers_stored` as fractions, so both spellings reach this port
                // and the magnitude is what tells them apart. 1 itself is read
                // as the fraction: P(T <= t) = 1 is a level, 1% is not one
                // anybody asks for beside 50, 90 and 99.
                if (p > 1.0) p /= 100.0;
                if (!(p > 0.0) || !(p < 1.0))
                    throw line::InputError(
                        "--percentiles levels lie strictly inside (0,1) as fractions or (0,100) "
                        "as percents; the 100th percentile of an unbounded law is not finite "
                        "(got '" + tok + "')");
                o.knobs.percentiles.push_back(p);
                if (comma == std::string::npos) break;
                at = comma + 1;
            }
        } else if (a == "--reward-name" || a == "--reward_name") {
            o.knobs.reward_name = next("--reward-name");
        } else if (a == "--notation") {
            const std::string v = next("--notation");
            // The exporter validates the name; it is not defaulted here, so an
            // unrecognised notation is refused rather than answered with scalar.
            o.knobs.notation = v;
        } else if (a == "--symbolic") {
            // `sym_resolve` validates the name: auto, none, a URL or an image.
            // It is not defaulted here, so `--symbolic none` stays local rather
            // than being read as "not given" and searching anyway.
            o.knobs.symbolic = next("--symbolic");
        } else if (a == "--equilibria") o.knobs.equilibria = true;
        else if (a == "--perm-engine") {
            // The analyzer validates the name; an unrecognised engine is
            // refused rather than answered with exact.
            o.knobs.method_perm = next("--perm-engine");
        } else if (a == "--transient-method") {
            // The analyzer validates the name, so an unrecognised one is
            // refused rather than answered with the default.
            o.knobs.transient_method = next("--transient-method");
        } else if (a == "--fau-epsilon") {
            const std::string v = next("--fau-epsilon");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0) || !std::isfinite(d))
                throw line::InputError(
                    "--fau-epsilon is the probability mass the transient grid may discard and "
                    "must be a positive finite number (got '" + v + "')");
            o.knobs.fau_epsilon = d;
        } else if (a == "--fau-delta") {
            const std::string v = next("--fau-delta");
            const double d = std::atof(v.c_str());
            if (!(d >= 0.0) || !std::isfinite(d))
                throw line::InputError(
                    "--fau-delta is the occupancy below which a state is dropped and must be a "
                    "nonnegative finite number (got '" + v + "')");
            o.knobs.fau_delta = d;
        } else if (a == "--cdf-algorithm") {
            // The analyzer validates the name; it is not defaulted here, so an
            // unrecognised algorithm is refused rather than answered with exact.
            o.knobs.cdf_algorithm = next("--cdf-algorithm");
        } else if (a == "--passage-from") o.knobs.passage_from = next("--passage-from");
        else if (a == "--passage-into") o.knobs.passage_into = next("--passage-into");
        else if (a == "--passage-method") o.knobs.passage_method = next("--passage-method");
        else if (a == "--passage-orders")
            o.knobs.passage_orders = static_cast<std::size_t>(std::stoul(next("--passage-orders")));
        else if (a == "--no-interlocking") o.knobs.no_interlocking = true;
        else if (a == "--layer-solver") o.knobs.layer_solver = next("--layer-solver");
        else if (a == "--stage-solver") o.knobs.stage_solver = next("--stage-solver");
        else if (a == "--ln-transient") o.knobs.ln_transient = next("--ln-transient");
        else if (a == "--ln-transient-channels")
            o.knobs.ln_transient_channels = next("--ln-transient-channels");
        else if (a == "--sens-method") o.knobs.sens_method = next("--sens-method");
        else if (a == "--sens-scheme") o.knobs.sens_scheme = next("--sens-scheme");
        else if (a == "--sens-step") {
            const std::string v = next("--sens-step");
            const double h = std::atof(v.c_str());
            if (!(h > 0.0) || !(h < 1.0))
                throw line::InputError(
                    "--sens-step is the RELATIVE rate perturbation and must lie in (0,1) (got '" +
                    v + "')");
            o.knobs.sens_step = h;
        }
        else if (a == "--uq-solver") o.knobs.uq_solver = next("--uq-solver");
        else if (a == "--keep") o.knobs.keep = true;
        else if (a == "--verbose") o.knobs.verbose = true;
        else if (a == "--remote") o.knobs.remote = true;
        else if (a == "--remote-url") {
            // Implies --remote: a URL given and then ignored because the flag
            // was forgotten would solve LOCALLY and report nothing about it.
            o.knobs.remote_url = next("--remote-url");
            o.knobs.remote = true;
        }
        else if (a == "--timeout") {
            const std::string v = next("--timeout");
            const long s = std::atol(v.c_str());
            if (s < 1)
                throw line::InputError("--timeout is a deadline in seconds and must be positive "
                                       "(got '" + v + "')");
            o.knobs.timeout_seconds = static_cast<int>(s);
        }
        else if (a == "--repeat") {
            const std::string v = next("--repeat");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError("--repeat must be a positive run count (got '" + v + "')");
            o.knobs.repeat = static_cast<int>(n);
        }
        // ---- the simulator's own knobs ------------------------------------
        else if (a == "--ldes-tranfilter") {
            const std::string v = next("--ldes-tranfilter");
            if (v != "mser5" && v != "fixed" && v != "none")
                throw line::InputError(
                    "--ldes-tranfilter selects the warmup filter and is mser5, fixed or none (got '" +
                    v + "')");
            o.knobs.ldes_tranfilter = v;
        }
        else if (a == "--ldes-warmupfrac") {
            const std::string v = next("--ldes-warmupfrac");
            const double d = std::atof(v.c_str());
            if (!(d >= 0.0 && d < 1.0))
                throw line::InputError(
                    "--ldes-warmupfrac is the fraction of the run the fixed filter discards and "
                    "lies in [0,1) (got '" + v + "')");
            o.knobs.ldes_warmupfrac = d;
        }
        else if (a == "--ldes-cimethod") {
            const std::string v = next("--ldes-cimethod");
            if (v != "obm" && v != "bm" && v != "spectral" && v != "none")
                throw line::InputError(
                    "--ldes-cimethod selects the confidence-interval estimator and is obm, bm, "
                    "spectral or none (got '" + v + "')");
            o.knobs.ldes_cimethod = v;
        }
        else if (a == "--ldes-cnvgon") o.knobs.ldes_cnvgon = true;
        else if (a == "--ldes-cnvgtol") {
            // Implies --ldes-cnvgon: a tolerance given and then ignored because
            // the switch was forgotten would run the full budget and say nothing.
            const std::string v = next("--ldes-cnvgtol");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0 && d < 1.0))
                throw line::InputError(
                    "--ldes-cnvgtol is a RELATIVE precision target and lies in (0,1) (got '" + v +
                    "')");
            o.knobs.ldes_cnvgtol = d;
            o.knobs.ldes_cnvgon = true;
        }
        else if (a == "--ldes-slotted") o.knobs.ldes_slotted = true;
        else if (a == "--slotted") o.knobs.slotted = true;
        else if (a == "--slotlength") {
            // Implies --slotted, as --ldes-slotlength does for the simulator.
            const std::string v = next("--slotlength");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0))
                throw line::InputError(
                    "--slotlength is the slot of the discrete time scale and must be positive "
                    "(got '" + v + "')");
            o.knobs.slotlength = d;
            o.knobs.slotted = true;
        }
        else if (a == "--ldes-slotlength") {
            // Implies --ldes-slotted, as --slotlength does on the engine's CLI.
            const std::string v = next("--ldes-slotlength");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0))
                throw line::InputError(
                    "--ldes-slotlength is the slot of the discrete time scale and must be positive "
                    "(got '" + v + "')");
            o.knobs.ldes_slotlength = d;
            o.knobs.ldes_slotted = true;
        }
        else if (a == "--ldes-replications") {
            const std::string v = next("--ldes-replications");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError(
                    "--ldes-replications is a positive count of independent runs (got '" + v + "')");
            o.knobs.ldes_replications = static_cast<int>(n);
        }
        else if (a == "--ldes-numthreads") {
            const std::string v = next("--ldes-numthreads");
            const long n = std::atol(v.c_str());
            if (n < 1)
                throw line::InputError(
                    "--ldes-numthreads is a positive worker count (got '" + v + "')");
            o.knobs.ldes_numthreads = static_cast<int>(n);
        }
        else if (a == "--ldes-maxtime") {
            const std::string v = next("--ldes-maxtime");
            const double d = std::atof(v.c_str());
            if (!(d > 0.0))
                throw line::InputError(
                    "--ldes-maxtime is a wall-clock budget in seconds and must be positive (got '" +
                    v + "')");
            o.knobs.ldes_maxtime = d;
        }
        else if (a == "--ldes-initsol") {
            // A STATION-MAJOR placement, [st0_cl0, st0_cl1, ...]; the engine
            // reads it as the initial state and skips its default placement.
            const std::string v = next("--ldes-initsol");
            std::size_t b = 0;
            while (b <= v.size()) {
                const std::size_t e = v.find(',', b);
                const std::string tok =
                    v.substr(b, e == std::string::npos ? std::string::npos : e - b);
                if (tok.empty())
                    throw line::InputError(
                        "--ldes-initsol is a comma-separated placement with no empty entry (got '" +
                        v + "')");
                o.knobs.ldes_initsol.push_back(std::atof(tok.c_str()));
                if (e == std::string::npos) break;
                b = e + 1;
            }
        }
        else if (a == "--ldes-rest-url") o.knobs.ldes_rest_url = next("--ldes-rest-url");
        else if (a == "-v" || a == "--verbosity") {
            const std::string v = next(a.c_str());
            // The JAR names two levels and this port's help documents three;
            // all five spellings are accepted, and an unknown one is refused
            // rather than read as `standard`, which would silence nothing while
            // reporting that it had.
            if (v != "silent" && v != "standard" && v != "normal" && v != "debug" &&
                v != "verbose")
                throw line::InputError(
                    "-v takes silent, standard (the JAR spells it normal) or debug; got '" + v +
                    "'");
            o.knobs.verbosity = (v == "normal") ? "standard" : v;
        }
        else if (!a.empty() && a[0] == '-')
            throw line::InputError("unknown option: " + a);
        else
            o.file = a;
    }
    return o;
}

}  // namespace

/**
 * Everything one invocation does once the arguments are in hand.
 *
 * SPLIT OUT OF `main` FOR SERVER MODE, which runs it once per request with a
 * different `-f` and a captured stdout. Keeping one body means a request served
 * over the socket takes exactly the path the same command line takes at the
 * shell -- the failure this avoids is a server that answers slightly differently
 * from the CLI it is supposed to BE.
 */
/**
 * `--find-solver`: which solvers and methods can analyze the model named by -f.
 *
 * It reports rather than solves, so it stops before the solver ladder in
 * `solve_model_dispatch` and before every knob that describes a run. The answer
 * is arithmetic-independent -- `auto_find_solver` asks feature sets, shapes and
 * gates, never numbers -- so the model is read at double whatever --arith says,
 * and a caller who passed one is told rather than silently obeyed.
 *
 * The layered path is not covered: `auto_find_solver` narrows the flat Network
 * families, and a LayeredNetwork's are `ln` and `lqns`, which this port reaches
 * through `solve_lqn_dispatch` and not through an AUTO of its own.
 */
int find_solver_report(const Options& o) {
    if (o.file.empty())
        throw line::InputError("--find-solver reports on a model; name one with -f");
    if (!o.input_given && (has_lqn_extension(o.file) || line::io::is_layered_json(o.file)))
        throw line::UnsupportedError(
            "--find-solver reports on a flat Network model; a layered one is solved by -s ln "
            "and -s lqns, which this port reaches through the -i lqnx path");
    line::qn::Network<double> net = read_model<double>(o.file);
    const std::vector<line::autosolver::SolverCandidate> rows = line::autosolver::auto_find_solver(
        net.get_struct(), o.find_solver_metric, o.find_solver_all);
    std::printf("%s", line::autosolver::auto_find_solver_table(rows).c_str());
    return 0;
}

int run_invocation(Options o) {
        // Validated on every path, not only --api: an unrecognised --arith is a
        // caller error whatever else the invocation asks for.
        line::reg::parse_arith(o.arith);
        if (!o.api.empty()) {
            if (o.output != "readable" && o.output != "json")
                throw line::InputError("unknown -o '" + o.output +
                                       "'; accepted forms are: readable, json");
            const line::reg::Json result =
                line::reg::api_invoke(o.api, o.arith, read_api_args(o.args));
            if (o.output == "json")
                std::printf("%s\n", result.dump(2).c_str());
            else
                std::printf("%s", line::reg::api_render_readable(result).c_str());
            return 0;
        }
        // A .lqnx/.xml path with no -i is taken as a layered model, so naming
        // the file is enough to solve it. An explicit -i always wins, so a
        // caller who says -i json about an oddly-named file still gets json.
        if (!o.input_given && has_lqn_extension(o.file)) o.input = "lqnx";
        // The same courtesy for a JMT document: naming the file is enough.
        if (!o.input_given && has_jsim_extension(o.file)) o.input = "jsimg";
        // And for a PNML document.
        if (!o.input_given && has_pnml_extension(o.file)) o.input = "pnml";
        // A model.json can carry EITHER model kind, and `linemodel_save` writes
        // `.json` for both, so the extension cannot decide it. The content can:
        // a `LayeredNetwork` type takes the layered path whatever `-i` says,
        // because handing it to the Network reader only produces "model type
        // 'LayeredNetwork' is not a Network" one frame further down.
        if (!o.file.empty() && o.input != "lqnx" && o.input != "xml" && o.input != "pnml" &&
            !has_pnml_extension(o.file) &&
            !has_jsim_extension(o.file) && o.input.compare(0, 4, "jsim") != 0 &&
            line::io::is_layered_json(o.file))
            o.input = "lqnx";
        if (o.input == "lqnx" || o.input == "xml") {
            if (o.output != "readable" && o.output != "json" && o.output != "layers")
                throw line::InputError("unknown -o '" + o.output +
                                       "'; accepted forms on the layered path are: readable, "
                                       "json, layers");
            // ONE ENVELOPE PER ANALYSIS, in the order asked for. `-a avg,sens`
            // is the JAR's comma list, and the JAR merges the results into one
            // JSON object; here each arm prints as it computes, so the multi
            // form emits a sequence of the SAME envelopes a single `-a` emits
            // (JSON Lines). That keeps one parser for both spellings, where a
            // merged object would need a second one for the multi case alone.
            const std::vector<std::string> as = analysis_list(o.analysis);
            for (std::size_t i = 0; i + 1 < as.size(); ++i) {
                const int rc = solve_lqn_dispatch(o.arith, o.solver, as[i], o.output, o.file,
                                                  o.knobs);
                if (rc != 0) return rc;
            }
            return solve_lqn_dispatch(o.arith, o.solver, as.back(), o.output, o.file, o.knobs);
        }
        // `-i jsim | jsimg | jsimw`: a JMT simulation document, read by
        // `read_jsim`. The three extensions name ONE format -- JMT writes the
        // same `<sim>` document under all three -- and are accepted separately
        // because the JAR CLI accepts them separately and a script naming the
        // wrong one is not describing a different model.
        if (o.input == "jsim" || o.input == "jsimg" || o.input == "jsimw")
            g_jsim_input = true;
        else if (o.input == "pnml")
            g_pnml_input = true;
        else if (o.input != "json")
            throw line::UnsupportedError(
                "the model-solving path reads -i json for a Network model, -i jsim|jsimg|jsimw "
                "for a JMT simulation document, -i pnml for a place/transition net and "
                "-i lqnx|xml for a layered one (got '" + o.input + "')");
        // `-o json` IS WIRED FOR EVERY ANALYSIS, not only `-a avg`. It was
        // avg-only, and the restriction was real rather than nominal: each other
        // arm printed its own readable shape and nothing else, so accepting
        // `json` for one would have handed a caller a table it asked to receive
        // as JSON -- the silent-acceptance failure the Knobs struct exists to
        // prevent, one layer up. The refusal is therefore removed only now that
        // every arm emits a payload (`emit_analysis`), and each answers under a
        // key named after its own `-a` so a caller can tell which question was
        // answered. On `-a avg` the solver BANNER still precedes the object on
        // stdout, as the JAR's does; it contains no brace, so the wrapper's
        // first-`{` scan is unaffected. The other arms print no banner: what it
        // carried -- the arithmetic, the resolved method, the state count and the
        // cutoff -- is inside the payload, where a host reads it as data instead
        // of scraping it from a line above.
        if (o.output != "readable" && o.output != "json")
            throw line::InputError("unknown -o '" + o.output +
                                   "'; accepted forms are: readable, json");
        g_json_output = (o.output == "json");
        // See the layered path above for why a comma list emits one envelope
        // per analysis rather than one merged object.
        const std::vector<std::string> as = analysis_list(o.analysis);
        for (std::size_t i = 0; i + 1 < as.size(); ++i) {
            const int rc = solve_model_dispatch(o.arith, o.solver, as[i], o.file, o.knobs);
            if (rc != 0) return rc;
        }
        return solve_model_dispatch(o.arith, o.solver, as.back(), o.file, o.knobs);
}

/**
 * Run one invocation with stdout captured, and return what it printed.
 *
 * SERVER MODE'S ONE PIECE OF MACHINERY. Every arm of this CLI writes its answer
 * with `printf`, which is the right thing for a command-line tool and leaves a
 * server nothing to send. Redirecting fd 1 around the call means the arms need
 * no server-aware variant and cannot drift from the command-line behaviour;
 * stderr is deliberately NOT captured, so a warning still reaches the operator's
 * console rather than being folded into the client's answer.
 */
std::string run_invocation_captured(Options o, int& rc) {
    const char* tmpdir = std::getenv("TMPDIR");
    std::string path = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/line-cli-out-XXXXXX";
    std::vector<char> buf(path.begin(), path.end());
    buf.push_back('\0');
    const int tfd = ::mkstemp(&buf[0]);
    if (tfd < 0) {
        rc = 3;
        return "line-cli: cannot create a capture file for the response\n";
    }
    path.assign(&buf[0]);
    std::fflush(stdout);
    const int saved = ::dup(1);
    ::dup2(tfd, 1);
    std::string err;
    try {
        rc = run_invocation(o);
    } catch (const line::Error& e) {
        rc = 2;
        err = std::string("line-cli: ") + e.what() + "\n";
    } catch (const std::exception& e) {
        rc = 3;
        err = std::string("line-cli: unexpected failure: ") + e.what() + "\n";
    }
    std::fflush(stdout);
    ::dup2(saved, 1);
    ::close(saved);
    ::close(tfd);
    std::ifstream in(path.c_str());
    std::string out((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    std::remove(path.c_str());
    // THE ERROR IS THE ANSWER when the solve failed: a client that receives an
    // empty message cannot tell a refusal from a model with no rows.
    return err.empty() ? out : out + err;
}

/**
 * `-p/--port`: serve solve requests over a WebSocket, as `LineWebSocketServer`
 * does.
 *
 * THE PROTOCOL IS THE JAR's, unchanged: one text message per connection, whose
 * FIRST LINE is the comma-separated argument list and whose remainder is the
 * model document. The JAR overwrites the first two arguments with `--file` and
 * the path it staged the document at, so the client's own first two tokens are
 * placeholders; the same substitution happens here, which is what lets an
 * existing client talk to this server without knowing which binary answered.
 */
int run_server(const Options& base) {
    line::ws::Server server(base.port);
    std::printf("--------------------------------------------------------------------\n");
    std::printf("LINE Solver - Command Line Interface (C++)\n");
    std::printf("Copyright (c) 2012-2026, QORE Lab, Imperial College London\n");
    std::printf("Version %s. All rights reserved.\n", kVersion);
    std::printf("--------------------------------------------------------------------\n");
    std::printf("Running in server mode on port %d.\n", base.port);
    if (base.maxreq)
        std::printf("Quitting after %d request(s).\n", base.maxreq);
    std::fflush(stdout);

    int served = 0;
    while (base.maxreq == 0 || served < base.maxreq) {
        const bool ok = server.serve_one([&](const std::string& msg) -> std::string {
            const std::string::size_type nl = msg.find('\n');
            if (nl == std::string::npos)
                return "line-cli: the request's first line is the argument list and its "
                       "remainder is the model document; this message has no newline\n";
            const std::string argline = msg.substr(0, nl);
            const std::string model = msg.substr(nl + 1);

            const char* tmpdir = std::getenv("TMPDIR");
            std::string path =
                std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/line-cli-req-XXXXXX";
            std::vector<char> nb(path.begin(), path.end());
            nb.push_back('\0');
            const int mfd = ::mkstemp(&nb[0]);
            if (mfd < 0) return "line-cli: cannot stage the client model\n";
            ::close(mfd);
            path.assign(&nb[0]);
            {
                std::ofstream mf(path.c_str());
                mf << model;
            }

            // The argument list, with the first two tokens replaced by the
            // staged path exactly as `LineWebSocketServer.onMessage` replaces
            // them. A list SHORTER than two is the client's error and is
            // reported rather than padded, since padding would solve the
            // default model instead of the one it sent.
            std::vector<std::string> toks;
            std::string::size_type at = 0;
            while (at <= argline.size()) {
                const std::string::size_type comma = argline.find(',', at);
                toks.push_back(argline.substr(
                    at, comma == std::string::npos ? std::string::npos : comma - at));
                if (comma == std::string::npos) break;
                at = comma + 1;
            }
            std::string result;
            if (toks.size() < 2) {
                result = "line-cli: the argument list needs at least two tokens; the first two "
                         "are replaced by --file and the staged model path\n";
            } else {
                toks[0] = "--file";
                toks[1] = path;
                std::vector<char*> argv;
                std::vector<std::string> store;
                store.push_back("line-cli");
                for (std::size_t i = 0; i < toks.size(); ++i) store.push_back(toks[i]);
                for (std::size_t i = 0; i < store.size(); ++i)
                    argv.push_back(const_cast<char*>(store[i].c_str()));
                int rc = 0;
                try {
                    Options ro = parse_args(static_cast<int>(argv.size()), &argv[0]);
                    // The server's own flags never travel into a request: a
                    // client that sent `-p` would otherwise make the server
                    // recurse into a second listener on the same process.
                    ro.port = 0;
                    ro.maxreq = 0;
                    result = run_invocation_captured(ro, rc);
                } catch (const line::Error& e) {
                    result = std::string("line-cli: ") + e.what() + "\n";
                }
            }
            std::remove(path.c_str());
            return result;
        });
        // A dropped client is not a reason to stop serving, and it does not
        // count against --maxreq either: the JAR counts REQUESTS, and a peer
        // that vanished before sending one made none.
        if (ok) ++served;
    }
    return 0;
}

int main(int argc, char** argv) {
    try {
        Options o = parse_args(argc, argv);
        // Before any arm runs: `read_model` raises the priority warning and has
        // no knobs of its own to read the level from.
        if (!o.knobs.verbosity.empty()) g_verbosity = o.knobs.verbosity;
        // Solver console: set for the whole process, so that the model compile
        // narrates too -- LineConsole::writes() falls back to the session level
        // when no run is open, and reading the model happens before any solver
        // exists. `-v debug` (and its `verbose` spelling) is the ONLY way in.
        line::util::LineConsole::set_verbose(
            g_verbosity == "silent"   ? line::util::VerboseLevel::SILENT
            : (g_verbosity == "debug" || g_verbosity == "verbose")
                                      ? line::util::VerboseLevel::DEBUG
                                      : line::util::VerboseLevel::STD);
        if (o.help_all) {
            print_help();
            return 0;
        }
        if (o.help || argc == 1) {
            print_brief_help();
            return 0;
        }
        if (o.version) {
            std::printf("line-cli %s\n", kVersion);
            return 0;
        }
        // A warning is not a failure: every dependency the check probes is
        // optional, so the command exits 0 whenever the check itself ran.
        if (o.install) {
            install_check();
            return 0;
        }
        if (o.list) {
            list_api();
            return 0;
        }
        if (o.find_solver) return find_solver_report(o);
        if (o.port) {
            // `-f` and `-p` together would be two model sources for one run;
            // the request carries the model in server mode, so a file named on
            // the command line could only be ignored.
            if (!o.file.empty())
                throw line::InputError(
                    "-p runs the solver as a server, where each request carries its own model; "
                    "-f names a model on the command line and the two cannot both be the source");
            return run_server(o);
        }
        return run_invocation(o);
    } catch (const line::Error& e) {
        std::fprintf(stderr, "line-cli: %s\n", e.what());
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "line-cli: unexpected failure: %s\n", e.what());
        return 3;
    }
}
