/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The --api dispatch table.
 *
 * One entry per exposed function, hand written: it pulls its named arguments
 * out of the parsed JSON and calls the templated function at the arithmetic
 * the caller selected. Registry-generated dispatch is the eventual goal, but a
 * generated table cannot know that pfqn_conv takes callables it must refuse or
 * that pfqn_comom's Z is a vector where pfqn_comomrm's is a matrix, so the
 * table is explicit while the argument shapes are being pinned down.
 *
 * Compiled once into a static library rather than instantiated per translation
 * unit: three arithmetics times a dozen algorithms is a lot of template
 * instantiation to repeat in the CLI, the tests and a future binding.
 */

#include "line/reg/api_dispatch.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/aoi/aoi_fcfs_dm1.h"
#include "line/api/aoi/aoi_fcfs_gim1.h"
#include "line/api/aoi/aoi_fcfs_md1.h"
#include "line/api/aoi/aoi_fcfs_mgi1.h"
#include "line/api/aoi/aoi_fcfs_mm1.h"
#include "line/api/aoi/aoi_lcfsd_gim1.h"
#include "line/api/aoi/aoi_lcfsd_mgi1.h"
#include "line/api/aoi/aoi_lcfspr_dm1.h"
#include "line/api/aoi/aoi_lcfspr_gim1.h"
#include "line/api/aoi/aoi_lcfspr_md1.h"
#include "line/api/aoi/aoi_lcfspr_mgi1.h"
#include "line/api/aoi/aoi_lcfspr_mm1.h"
#include "line/api/aoi/aoi_lcfss_gim1.h"
#include "line/api/aoi/aoi_lcfss_mgi1.h"
#include "line/api/aoi/aoi_lst_det.h"
#include "line/api/aoi/aoi_lst_erlang.h"
#include "line/api/aoi/aoi_lst_exp.h"
#include "line/api/aoi/aoi_lst_ph.h"
#include "line/api/fj/fj_bounds.h"
#include "line/api/fj/fj_char_max.h"
#include "line/api/fj/fj_gk_bound.h"
#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_quantile.h"
#include "line/api/fj/fj_ordstat_exp.h"
#include "line/api/fj/fj_quorum_moments.h"
#include "line/api/fj/fj_respt_2way.h"
#include "line/api/fj/fj_respt_nt.h"
#include "line/api/fj/fj_respt_varki.h"
#include "line/api/fj/fj_respt_vm.h"
#include "line/api/fj/fj_rmax.h"
#include "line/api/fj/fj_rmax_erlang.h"
#include "line/api/fj/fj_rmax_evd.h"
#include "line/api/fj/fj_sm_tput.h"
#include "line/api/fj/fj_synch_delay.h"
#include "line/api/fj/fj_xmax_2.h"
#include "line/api/fj/fj_xmax_approx.h"
#include "line/api/fj/fj_xmax_emma.h"
#include "line/api/fj/fj_xmax_erlang.h"
#include "line/api/fj/fj_xmax_exp.h"
#include "line/api/fj/fj_xmax_hyperexp.h"
#include "line/api/fj/fj_xmax_normal.h"
#include "line/api/fj/fj_xmax_pareto.h"
#include "line/api/moment/moment_binomial_from_factorial.h"
#include "line/api/moment/moment_binomial_from_negbinomial.h"
#include "line/api/moment/moment_binotrans.h"
#include "line/api/moment/moment_binotransinv.h"
#include "line/api/moment/moment_central_from_raw.h"
#include "line/api/moment/moment_factorial_from_binomial.h"
#include "line/api/moment/moment_factorial_from_raw.h"
#include "line/api/moment/moment_factorial_from_upfactorial.h"
#include "line/api/moment/moment_lah.h"
#include "line/api/moment/moment_negbinomial_from_binomial.h"
#include "line/api/moment/moment_negbinomial_from_upfactorial.h"
#include "line/api/moment/moment_raw_from_central.h"
#include "line/api/moment/moment_raw_from_factorial.h"
#include "line/api/moment/moment_raw_from_upfactorial.h"
#include "line/api/moment/moment_stirling1.h"
#include "line/api/moment/moment_stirling2.h"
#include "line/api/moment/moment_stirlingcycle.h"
#include "line/api/moment/moment_upfactorial_from_factorial.h"
#include "line/api/moment/moment_upfactorial_from_negbinomial.h"
#include "line/api/moment/moment_upfactorial_from_raw.h"
#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_fau.h"
#include "line/api/mc/ctmc_foxglynn.h"
#include "line/api/mc/ctmc_saddlepoint.h"
#include "line/api/mc/ctmc_bicgstab.h"
#include "line/api/mc/ctmc_gmres.h"
#include "line/api/mc/ctmc_gmres_multi.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_multi.h"
#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_sens.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/ctmc_solve_reducible_blkdecomp.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/api/mc/ctmc_transient_sens.h"
#include "line/api/mc/dtmc_makestochastic.h"
#include "line/api/mc/dtmc_solve_reducible.h"
#include "line/api/mc/stronglyconncomp.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comom.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_comomrm_ms.h"
#include "line/api/pfqn/pfqn_conv.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_procomom.h"
#include "line/api/pfqn/pfqn_recal.h"
#include "line/api/sim/sim_firquest.h"
#include "line/api/qsys/qsys_bmapphnn_retrial.h"
#include "line/api/qsys/qsys_dmc.h"
#include "line/api/dqsys/dqsys_geogeo1.h"
#include "line/api/dqsys/dqsys_geoxgeo1.h"
#include "line/api/qsys/qsys_gg1.h"
#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gig1_approx_gelenbe.h"
#include "line/api/qsys/qsys_gig1_approx_heyman.h"
#include "line/api/qsys/qsys_gig1_approx_kimura.h"
#include "line/api/qsys/qsys_gig1_approx_klb.h"
#include "line/api/qsys/qsys_gig1_approx_kobayashi.h"
#include "line/api/qsys/qsys_gig1_approx_marchal.h"
#include "line/api/qsys/qsys_gig1_approx_myskja.h"
#include "line/api/qsys/qsys_gig1_approx_myskja2.h"
#include "line/api/qsys/qsys_gig1_lbnd.h"
#include "line/api/qsys/qsys_gig1_rqt.h"
#include "line/api/qsys/qsys_gig1_ubnd_kingman.h"
#include "line/api/qsys/qsys_gigk_approx.h"
#include "line/api/qsys/qsys_gigk_approx_cosmetatos.h"
#include "line/api/qsys/qsys_gigk_approx_kingman.h"
#include "line/api/qsys/qsys_gigk_approx_whitt.h"
#include "line/api/qsys/qsys_ggnm_diffusion.h"
#include "line/api/qsys/qsys_gig1_bnds_extremal.h"
#include "line/api/qsys/qsys_mgisrgi_whitt.h"
#include "line/api/qsys/qsys_mmk_qed.h"
#include "line/api/qsys/qsys_gigk_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt_gamma.h"
#include "line/api/qsys/qsys_gm1.h"
#include "line/api/qsys/qsys_mapd1.h"
#include "line/api/qsys/qsys_mapdc.h"
#include "line/api/qsys/qsys_mapg1.h"
#include "line/api/qsys/qsys_mapm1.h"
#include "line/api/qsys/qsys_mapmap1.h"
#include "line/api/qsys/qsys_mapmc.h"
#include "line/api/qsys/qsys_mapph1.h"
#include "line/api/qsys/qsys_mapphc.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mg1_fb.h"
#include "line/api/qsys/qsys_mg1_lrpt.h"
#include "line/api/qsys/qsys_mg1_prio.h"
#include "line/api/qsys/qsys_mg1_psjf.h"
#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_mg1_srpt.h"
#include "line/api/qsys/qsys_mg1k_loss_mgs.h"
#include "line/api/qsys/qsys_mginf.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mm1_dps.h"
#include "line/api/qsys/qsys_mm1k_loss.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/api/qsys/qsys_mmcc_retrial_fp.h"
#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_mxm1.h"
#include "line/api/qsys/qsys_phm1.h"
#include "line/api/qsys/qsys_phmc.h"
#include "line/api/qsys/qsys_phph1.h"
#include "line/api/sim/sim_fquest.h"
#include "line/api/sim/sim_quest_options.h"
#include "line/api/sim/sim_shapirowilk.h"
#include "line/api/sim/sim_sts_quantile_areas.h"
#include "line/api/sim/sim_vonneumann.h"

namespace line {
namespace reg {

// Arithmetic selection

std::string ArithSpec::str() const {
    switch (mode) {
        case Arith::Double: return "double";
        case Arith::Exact: return "exact";
        case Arith::Real: return "real:" + std::to_string(digits);
    }
    return "?";
}

ArithSpec parse_arith(const std::string& text) {
    ArithSpec s;
    if (text == "double") {
        s.mode = Arith::Double;
        return s;
    }
    if (text == "exact") {
        s.mode = Arith::Exact;
        return s;
    }
    if (text.compare(0, 5, "real:") == 0 || text == "real") {
        unsigned want = 50;
        if (text != "real") {
            const std::string d = text.substr(5);
            if (d.empty() || d.find_first_not_of("0123456789") != std::string::npos)
                throw InputError("--arith real:<digits> needs a decimal digit count, got '" + d +
                                 "'");
            want = static_cast<unsigned>(std::strtoul(d.c_str(), nullptr, 10));
            if (want == 0) throw InputError("--arith real:<digits> needs at least one digit");
        }
        s.mode = Arith::Real;
        // Round UP to an instantiated tier: never give less precision than asked.
        if (want <= 50)
            s.digits = 50;
        else if (want <= 100)
            s.digits = 100;
        else if (want <= 200)
            s.digits = 200;
        else
            throw UnsupportedError(
                "--arith real:" + std::to_string(want) +
                " exceeds the precision tiers this build instantiates (50, 100, 200 digits); "
                "a wider tier has to be compiled in, it cannot be selected at run time.");
        return s;
    }
    throw InputError("unknown --arith '" + text +
                     "'; accepted forms are: double, exact, real:<digits> "
                     "(digits rounded up to 50, 100 or 200)");
}

// ---------------------------------------------------------------------------
// The table
// ---------------------------------------------------------------------------

namespace {

/**
 * Run an operation at the selected arithmetic. Op is a struct with a static
 * member template run<T>(Args&) returning the results object; the arithmetic
 * gate against the registry has already passed when this is reached, so every
 * branch here is a mode the function genuinely supports.
 */
template <class Op>
Json run_at(const ArithSpec& s, Args& a) {
    switch (s.mode) {
        case Arith::Double: return Op::template run<double>(a);
        case Arith::Exact: return Op::template run<Rational>(a);
        case Arith::Real:
            if (s.digits <= 50) return Op::template run<Real50>(a);
            if (s.digits <= 100) return Op::template run<Real100>(a);
            return Op::template run<Real200>(a);
    }
    throw InputError("unreachable arithmetic selection");
}

/**
 * Wrap an Op whose registry entry lists Double and Real but NOT Exact.
 *
 * `run_at` instantiates all three backends because the mode is a runtime value,
 * so an Op that calls exp, log or sqrt would not COMPILE at Rational even
 * though it can never be REACHED there -- `api_invoke` refuses the mode against
 * the registry before dispatching. The if-constexpr makes that same fact known
 * to the compiler. The throw is the arm the gate has already excluded, kept as
 * a real refusal rather than a placeholder so that a registry entry mistakenly
 * widened to Exact reports the reason instead of returning a wrong number.
 */
template <class Op>
struct NeedsTranscendental {
    template <class T>
    static Json run(Args& a) {
        if constexpr (num_traits<T>::has_transcendental) {
            return Op::template run<T>(a);
        } else {
            (void)a;
            throw UnsupportedError(
                "this function evaluates transcendental terms (exp, log, sqrt or a "
                "matrix exponential), which the exact rational field does not contain; "
                "rerun with --arith double or --arith real:<digits>.");
        }
    }
};

/** G and lG, the shape most of the normalizing-constant family returns. */
template <class T>
Json nc_results(const pfqn::NcResult<T>& r) {
    Json out;
    out["G"] = encode_scalar(r.G);
    out["lG"] = r.lG;
    return out;
}

struct OpPfqnCa {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        a.done();
        return nc_results(pfqn::pfqn_ca(L, N, Z));
    }
};

struct OpPfqnConv {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        a.unsupported("cdscaling",
                      "the class-dependence scaling argument is a vector of callables, which has "
                      "no JSON representation; call pfqn_conv from the library for that case");
        a.unsupported("options", "solver options are not carried over the --api boundary");
        a.done();
        return nc_results(pfqn::pfqn_conv(L, N, Z));
    }
};

struct OpPfqnRecal {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        const std::vector<int> m0 = a.ints_or_empty("m0");
        a.done();
        return nc_results(pfqn::pfqn_recal(L, N, Z, m0));
    }
};

struct OpPfqnGld {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> mu = a.matrix<T>("mu");
        a.unsupported("options", "solver options are not carried over the --api boundary");
        a.done();
        return nc_results(pfqn::pfqn_gld(L, N, mu));
    }
};

struct OpPfqnMva {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        const std::vector<int> mi = a.ints_or_empty("mi");
        a.done();
        const pfqn::MvaResult<T> r = pfqn::pfqn_mva(L, N, Z, mi);
        Json out;
        out["XN"] = encode_vector(r.XN);
        out["QN"] = encode_matrix(r.QN);
        out["UN"] = encode_matrix(r.UN);
        out["CN"] = encode_matrix(r.CN);
        out["G"] = encode_scalar(r.G);
        out["lG"] = r.lG;
        return out;
    }
};

/** G, lG and the CoMoM basis. */
template <class T>
Json comom_results(const pfqn::ComomResult<T>& r) {
    Json out;
    out["G"] = encode_scalar(r.G);
    out["lG"] = r.lG;
    out["basis"] = encode_vector(r.basis);
    return out;
}

struct OpPfqnComom {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const std::vector<T> Z = a.vector_or_empty<T>("Z");
        const T atol = a.scalar<T>("atol", num_traits<T>::from_int(0));
        a.done();
        return comom_results(pfqn::pfqn_comom(L, N, Z, atol));
    }
};

struct OpPfqnComomrm {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        const int m = a.integer("m", 1);
        a.done();
        return comom_results(pfqn::pfqn_comomrm(L, N, Z, m));
    }
};

struct OpPfqnComomrmOrig {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        const T atol = a.scalar<T>("atol", num_traits<T>::from_int(0));
        a.done();
        return comom_results(pfqn::pfqn_comomrm_orig(L, N, Z, atol));
    }
};

struct OpPfqnComomrmMs {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const Matrix<T> Z = a.matrix_or_empty<T>("Z");
        const int m = a.integer("m", 1);
        const int S = a.integer("S", 1);
        a.done();
        const pfqn::ComomRmResult<T> r = pfqn::pfqn_comomrm_ms(L, N, Z, m, S);
        Json out;
        out["G"] = encode_scalar(r.G);
        out["lG"] = r.lG;
        out["prob"] = encode_vector(r.prob);
        return out;
    }
};

struct OpPfqnProcomom {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> L = a.matrix<T>("L");
        const std::vector<int> N = a.ints("N");
        const std::vector<T> Z = a.vector_or_empty<T>("Z");
        const T atol = a.scalar<T>("atol", num_traits<T>::from_int(0));
        a.done();
        const pfqn::ProcomomResult<T> r = pfqn::pfqn_procomom(L, N, Z, atol);
        Json out;
        out["Pr"] = encode_matrix(r.Pr);
        out["Q"] = encode_vector(r.Q);
        out["rankdef"] = r.rankdef;
        return out;
    }
};

struct OpCtmcSolve {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        a.unsupported("options", "solver options are not carried over the --api boundary");
        a.done();
        Json out;
        out["p"] = encode_vector(mc::ctmc_solve(Q));
        return out;
    }
};

struct OpDtmcSolve {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> P = a.matrix<T>("P");
        a.unsupported("options", "solver options are not carried over the --api boundary");
        a.done();
        Json out;
        out["PROB"] = encode_vector(mc::dtmc_solve(P));
        return out;
    }
};

struct OpCtmcMakeinfgen {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        a.done();
        Json out;
        out["Q"] = encode_matrix(mc::ctmc_makeinfgen(Q));
        return out;
    }
};

// ---------------------------------------------------------------------------
// api/sim, the simulation output analysis family
// ---------------------------------------------------------------------------

/*
 * EXACT ARITHMETIC IS REFUSED BY EVERY ENTRY BELOW, and the refusal has to be a
 * runtime branch even though the registry already gates it: `run_at` names all
 * three arithmetics in one switch, so the Rational instantiation is COMPILED
 * whether or not it can be reached, and the family's `static_assert` would fire
 * at build time. The `if constexpr` guard is what keeps the translation unit
 * compilable, exactly as the mdd path in the CLI does.
 *
 * There is no `--api sym_*` counterpart. `api/sym` is a backend, not a family of
 * numeric functions: it has no MATLAB `sym_*` entry points to mirror, and what
 * it serves reaches the user through the solvers that call it (`--symbolic` on
 * fluid, the CTMC symbolic getters). Exposing a REST client as an `--api`
 * function would invent a signature that exists in no other codebase.
 */

/** Option overrides and whether the caller gave any, which FIRQUEST reads. */
struct QuestOptionsRead {
    sim::QuestOptions opt;
    bool supplied = false;
};

/** The QUEST option overrides, read from the same keys MATLAB's struct uses. */
QuestOptionsRead read_quest_options(Args& a) {
    QuestOptionsRead out;
    const char* keys[] = {"b0", "m0", "s", "beta", "eta", "theta", "weight", "force"};
    for (std::size_t k = 0; k < sizeof(keys) / sizeof(keys[0]); ++k)
        if (a.has(keys[k])) out.supplied = true;
    out.opt.b0 = static_cast<long>(a.integer("b0", static_cast<int>(out.opt.b0)));
    out.opt.m0 = static_cast<long>(a.integer("m0", static_cast<int>(out.opt.m0)));
    const std::vector<int> s = a.ints_or_empty("s");
    if (!s.empty()) {
        out.opt.s.clear();
        for (std::size_t i = 0; i < s.size(); ++i) out.opt.s.push_back(static_cast<long>(s[i]));
    }
    out.opt.beta = a.scalar<double>("beta", out.opt.beta);
    out.opt.eta = a.scalar<double>("eta", out.opt.eta);
    out.opt.theta = a.scalar<double>("theta", out.opt.theta);
    out.opt.weight = a.scalar<double>("weight", out.opt.weight);
    out.opt.force = a.integer("force", out.opt.force ? 1 : 0) != 0;
    return out;
}

/** The shared shape of `sim_fquest` and `sim_firquest`. */
template <class T>
Json quest_results(const sim::QuestResult<T>& r) {
    Json out;
    out["estimate"] = encode_scalar(r.estimate);
    out["lower"] = encode_scalar(r.lower);
    out["upper"] = encode_scalar(r.upper);
    out["halfwidth"] = encode_scalar(r.halfwidth);
    out["b"] = static_cast<long long>(r.b);
    out["m"] = static_cast<long long>(r.m);
    out["n"] = static_cast<long long>(r.n);
    out["R"] = static_cast<long long>(r.R);
    out["truncated"] = static_cast<long long>(r.truncated);
    out["Ap"] = encode_scalar(r.Ap);
    out["Np"] = encode_scalar(r.Np);
    out["Vp"] = encode_scalar(r.Vp);
    out["heuristic"] = r.heuristic;
    out["warnings"] = r.warnings;
    return out;
}

/** The message every exact-arithmetic refusal in this family shares. */
std::string sim_exact_refusal(const char* fn, const char* why) {
    return std::string(fn) + ": " + why +
           ", so exact rational arithmetic has nothing to preserve and is refused; rerun with "
           "--arith double or --arith real";
}

struct OpSimVonneumann {
    template <class T>
    static Json run(Args& a) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                sim_exact_refusal("sim_vonneumann", "the p-value is a normal tail"));
        } else {
            const std::vector<T> x = vector_from_json<T>(a.get("x"), "sim_vonneumann: x");
            const double alpha = a.scalar<double>("alpha", 0.05);
            a.done();
            const sim::VonNeumannResult<T> r = sim::sim_vonneumann(x, alpha);
            Json out;
            out["ratio"] = encode_scalar(r.ratio);
            out["zscore"] = r.zscore;
            out["pvalue"] = r.pvalue;
            out["reject"] = r.reject;
            out["nobs"] = static_cast<long long>(r.nobs);
            return out;
        }
    }
};

struct OpSimShapirowilk {
    template <class T>
    static Json run(Args& a) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(sim_exact_refusal(
                "sim_shapirowilk", "the weights and the p-value are normal-order statistics"));
        } else {
            const std::vector<T> x = vector_from_json<T>(a.get("x"), "sim_shapirowilk: x");
            const double alpha = a.scalar<double>("alpha", 0.05);
            a.done();
            const sim::ShapiroWilkResult<T> r = sim::sim_shapirowilk(x, alpha);
            Json out;
            out["W"] = encode_scalar(r.W);
            out["pvalue"] = r.pvalue;
            out["zscore"] = r.zscore;
            out["reject"] = r.reject;
            out["nobs"] = static_cast<long long>(r.nobs);
            return out;
        }
    }
};

struct OpSimStsQuantileAreas {
    template <class T>
    static Json run(Args& a) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(sim_exact_refusal(
                "sim_sts_quantile_areas", "the areas carry the irrational weight sqrt(12)"));
        } else {
            const std::vector<T> Y =
                vector_from_json<T>(a.get("Y"), "sim_sts_quantile_areas: Y");
            const int b = a.integer("b", 0);
            const int m = a.integer("m", 0);
            const double p = a.scalar<double>("p", 0.5);
            const double weight = a.scalar<double>("weight", std::sqrt(12.0));
            a.done();
            if (b <= 0 || m <= 0)
                throw InputError(
                    "sim_sts_quantile_areas: the batch count b and batch size m are required and "
                    "must be positive integers");
            const sim::StsQuantileStats<T> r = sim::sim_sts_quantile_areas(
                Y, static_cast<std::size_t>(b), static_cast<std::size_t>(m), p, weight);
            Json out;
            out["areas"] = encode_vector(r.areas);
            out["bqe"] = encode_vector(r.bqe);
            out["quantile"] = encode_scalar(r.quantile);
            out["Ap"] = encode_scalar(r.Ap);
            out["Np"] = encode_scalar(r.Np);
            out["Vp"] = encode_scalar(r.Vp);
            out["b"] = static_cast<long long>(r.b);
            out["m"] = static_cast<long long>(r.m);
            out["n"] = static_cast<long long>(r.n);
            return out;
        }
    }
};

struct OpSimFquest {
    template <class T>
    static Json run(Args& a) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(sim_exact_refusal(
                "sim_fquest", "the interval is a t quantile times a square root"));
        } else {
            const std::vector<T> Y = vector_from_json<T>(a.get("Y"), "sim_fquest: Y");
            const double p = a.scalar<double>("p", 0.5);
            const double alpha = a.scalar<double>("alpha", 0.05);
            const QuestOptionsRead o = read_quest_options(a);
            a.done();
            return quest_results(sim::sim_fquest(Y, p, alpha, o.opt));
        }
    }
};

struct OpSimFirquest {
    template <class T>
    static Json run(Args& a) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(sim_exact_refusal(
                "sim_firquest", "the interval is a t quantile times a square root"));
        } else {
            // MATLAB passes an (n x R) matrix, one column per replication; the
            // C++ signature is replication-major, so the columns are transposed
            // here rather than reinterpreted -- pooling order is what the
            // procedure's stage tests act on.
            const Matrix<T> Ym = a.matrix<T>("Y");
            const double p = a.scalar<double>("p", 0.5);
            const double alpha = a.scalar<double>("alpha", 0.05);
            const QuestOptionsRead o = read_quest_options(a);
            a.done();
            std::vector<std::vector<T> > Y(Ym.cols(), std::vector<T>(Ym.rows()));
            for (std::size_t i = 0; i < Ym.rows(); ++i)
                for (std::size_t r = 0; r < Ym.cols(); ++r) Y[r][i] = Ym(i, r);
            // A missing option set is not the same as the FQUEST defaults here:
            // FIRQUEST's own defaults depend on R, and nullptr is how they are
            // asked for.
            return quest_results(
                sim::sim_firquest(Y, p, alpha, o.supplied ? &o.opt : nullptr));
        }
    }
};

// ---------------------------------------------------------------------------
// qsys: the closed-form queueing systems
// ---------------------------------------------------------------------------
//
// The most mechanical part of the surface -- scalar rates and coefficients of
// variation in, a named result struct out -- so it is exposed wholesale rather
// than one function at a time. What is NOT exposed is exactly the set whose
// MATLAB signature takes a FUNCTION: qsys_gig1_rq takes the interarrival law as
// a handle, qsys_mg1k_loss a density, qsys_mapg1k / qsys_mapg1k_perflow /
// qsys_mmapg1k a service law and qsys_ldps_workload a workload law. A callable
// has no JSON representation, so those six report that from api_invoke's
// not-exposed arm rather than being given an argument shape that silently
// stands in for the caller's function.
//
// AN ABSENT OPTIONAL ARGUMENT CALLS THE PORT'S OWN SHORT OVERLOAD, never a
// default written out here. Several of these functions carry a tuned default
// (qsys_mapd1's 4096 arrival truncation, qsys_phmc's 50000 iterations,
// qsys_phm1's 1e-16) and a plausible-looking value restated at this boundary
// would answer a different question from the library and the reference.
//
// The KEYS ARE THE MATLAB PARAMETER NAMES, which is why qsys_mapmap1 reads its
// ARRIVAL process from (C0, C1) and its SERVICE process from (D0, D1): that is
// the order matlab/src/api/qsys/qsys_mapmap1.m declares, and the opposite
// reading silently swaps the two processes of a stable queue for an unstable
// one. The returned keys are the port's struct field names.

/** [W, rhohat], what most of the family returns. */
template <class T>
Json qsys_results(const qsys::QsysResult<T>& r) {
    Json out;
    out["W"] = encode_scalar(r.W);
    // MATLAB spells this `rho` in qsys_mm1/qsys_mmk and `rhohat` everywhere
    // else; it is one quantity, the modified utilization, and it is reported
    // under the port's single name so a host does not key on the file it came
    // from.
    out["rhohat"] = encode_scalar(r.rhohat);
    return out;
}

#define LINE_QSYS_GIG1(OpName, fn)                                    \
    struct OpName {                                                   \
        template <class T>                                            \
        static Json run(Args& a) {                                    \
            const T lambda = a.number<T>("lambda");                   \
            const T mu = a.number<T>("mu");                           \
            const T ca = a.number<T>("ca");                           \
            const T cs = a.number<T>("cs");                           \
            a.done();                                                 \
            return qsys_results(qsys::fn(lambda, mu, ca, cs));        \
        }                                                             \
    }

LINE_QSYS_GIG1(OpQsysGig1AllenCunneen, qsys_gig1_approx_allencunneen);
LINE_QSYS_GIG1(OpQsysGig1Gelenbe, qsys_gig1_approx_gelenbe);
LINE_QSYS_GIG1(OpQsysGig1Heyman, qsys_gig1_approx_heyman);
LINE_QSYS_GIG1(OpQsysGig1Kimura, qsys_gig1_approx_kimura);
LINE_QSYS_GIG1(OpQsysGig1Klb, qsys_gig1_approx_klb);
LINE_QSYS_GIG1(OpQsysGig1Kobayashi, qsys_gig1_approx_kobayashi);
LINE_QSYS_GIG1(OpQsysGig1Marchal, qsys_gig1_approx_marchal);
LINE_QSYS_GIG1(OpQsysGig1Lbnd, qsys_gig1_lbnd);
LINE_QSYS_GIG1(OpQsysGig1UbndKingman, qsys_gig1_ubnd_kingman);
#undef LINE_QSYS_GIG1

#define LINE_QSYS_GIGK(OpName, fn)                                    \
    struct OpName {                                                   \
        template <class T>                                            \
        static Json run(Args& a) {                                    \
            const T lambda = a.number<T>("lambda");                   \
            const T mu = a.number<T>("mu");                           \
            const T ca = a.number<T>("ca");                           \
            const T cs = a.number<T>("cs");                           \
            const unsigned k = a.uinteger("k");                       \
            a.done();                                                 \
            return qsys_results(qsys::fn(lambda, mu, ca, cs, k));     \
        }                                                             \
    }

LINE_QSYS_GIGK(OpQsysGigkApprox, qsys_gigk_approx);
LINE_QSYS_GIGK(OpQsysGigkCosmetatos, qsys_gigk_approx_cosmetatos);
LINE_QSYS_GIGK(OpQsysGigkKingman, qsys_gigk_approx_kingman);
LINE_QSYS_GIGK(OpQsysGigkWhitt, qsys_gigk_approx_whitt);
#undef LINE_QSYS_GIGK

// The abandonment and QED families. Only the entry points whose arguments are
// JSON-expressible are exposed here: qsys_mgisrgi_whitt with a general patience
// law, qsys_ggisgi_fluid and qsys_mtginf all take CALLABLES (a hazard, a ccdf,
// an arrival rate), which no argument object can carry, so they stay
// library-only. qsys_erlanga is the exponential-patience case of the first, and
// its whole model is five numbers.
struct OpQsysErlangA {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T theta = a.number<T>("theta");
        const unsigned s = a.uinteger("s");
        // JSON has no infinity, so an absent r is the unbounded waiting room.
        const double r = a.has("r") ? num_traits<double>::to_double(a.number<double>("r"))
                                    : std::numeric_limits<double>::infinity();
        qsys::MgisrgiOptions opts;
        opts.wPoints = a.vector_or_empty<double>("wPoints");
        opts.maxQueue = a.count("maxQueue", qsys::MgisrgiOptions().maxQueue);
        a.done();
        const qsys::QsysAbandonResult<T> res = qsys::qsys_erlanga(lambda, mu, theta, s, r, opts);
        Json out;
        out["queueLengthDist"] = encode_vector(res.queueLengthDist);
        out["probLoss"] = encode_scalar(res.probLoss);
        out["probNoWait"] = encode_scalar(res.probNoWait);
        out["probServed"] = encode_scalar(res.probServed);
        out["probAbandon"] = encode_scalar(res.probAbandon);
        out["meanNumber"] = encode_scalar(res.meanNumber);
        out["varNumber"] = encode_scalar(res.varNumber);
        out["meanQueueLength"] = encode_scalar(res.meanQueueLength);
        out["varQueueLength"] = encode_scalar(res.varQueueLength);
        out["utilization"] = encode_scalar(res.utilization);
        out["throughput"] = encode_scalar(res.throughput);
        out["abandonRate"] = encode_scalar(res.abandonRate);
        out["meanWaitServed"] = encode_scalar(res.meanWaitServed);
        out["varWaitServed"] = encode_scalar(res.varWaitServed);
        out["meanWaitAbandon"] = encode_scalar(res.meanWaitAbandon);
        out["varWaitAbandon"] = encode_scalar(res.varWaitAbandon);
        out["meanWait"] = encode_scalar(res.meanWait);
        out["secondMomentWait"] = encode_scalar(res.secondMomentWait);
        out["numWaitingSpaces"] = static_cast<double>(res.numWaitingSpaces);
        if (!res.waitPoints.empty()) {
            out["waitPoints"] = encode_vector(res.waitPoints);
            out["cdfWaitServed"] = encode_vector(res.cdfWaitServed);
            out["cdfWaitAbandon"] = encode_vector(res.cdfWaitAbandon);
            out["cdfWait"] = encode_vector(res.cdfWait);
        }
        return out;
    }
};

struct OpQsysGgnmDiffusion {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned n = a.uinteger("n");
        // JSON has no infinity, so an absent m is the unbounded waiting room.
        const double m = a.has("m") ? a.number<double>("m")
                                    : std::numeric_limits<double>::infinity();
        const T ca = a.number<T>("ca");
        const T cs = a.number<T>("cs");
        a.done();
        // The service ccdf is a callable, so the CLI can only offer the
        // exponential default; the library entry point takes the general law.
        const qsys::QsysGgnmResult<T> r = qsys::qsys_ggnm_diffusion<T>(lambda, mu, n, m, ca, cs);
        Json out;
        out["beta"] = encode_scalar(r.beta);
        out["peakedness"] = encode_scalar(r.peakedness);
        out["variability"] = encode_scalar(r.variability);
        out["probDelay"] = encode_scalar(r.probDelay);
        out["probBlock"] = encode_scalar(r.probBlock);
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanNumber"] = encode_scalar(r.meanNumber);
        out["meanWait"] = encode_scalar(r.meanWait);
        out["utilization"] = encode_scalar(r.utilization);
        out["throughput"] = encode_scalar(r.throughput);
        return out;
    }
};

struct OpQsysGig1BndsExtremal {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T ca = a.number<T>("ca");
        const T cs = a.number<T>("cs");
        const std::size_t K = a.count("K", 4000);
        const std::size_t N = a.count("N", 2000);
        const bool skipTight = a.boolean("skipTight", false);
        a.done();
        const qsys::Gig1ExtremalResult<T> r =
            qsys::qsys_gig1_bnds_extremal(lambda, mu, ca, cs, K, N, skipTight);
        Json out;
        out["trafficIntensity"] = encode_scalar(r.trafficIntensity);
        out["lowerBound"] = encode_scalar(r.lowerBound);
        out["upperBound"] = encode_scalar(r.upperBound);
        out["upperBoundClosed"] = encode_scalar(r.upperBoundClosed);
        out["upperBoundDaley"] = encode_scalar(r.upperBoundDaley);
        out["upperBoundKingman"] = encode_scalar(r.upperBoundKingman);
        out["heavyTraffic"] = encode_scalar(r.heavyTraffic);
        out["delta"] = encode_scalar(r.delta);
        out["relativeWidth"] = encode_scalar(r.relativeWidth);
        out["tightComputed"] = r.tightComputed;
        return out;
    }
};

struct OpQsysMmkQed {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned s = a.uinteger("s");
        a.done();
        const qsys::QsysQedResult<T> r = qsys::qsys_mmk_qed(lambda, mu, s);
        Json out;
        out["offeredLoad"] = encode_scalar(r.offeredLoad);
        out["trafficIntensity"] = encode_scalar(r.trafficIntensity);
        out["beta"] = encode_scalar(r.beta);
        out["probDelay"] = encode_scalar(r.probDelay);
        out["meanWaitDelayed"] = encode_scalar(r.meanWaitDelayed);
        out["meanWait"] = encode_scalar(r.meanWait);
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanNumber"] = encode_scalar(r.meanNumber);
        out["utilization"] = encode_scalar(r.utilization);
        return out;
    }
};

struct OpQsysMmkQedAlpha {
    template <class T>
    static Json run(Args& a) {
        const T beta = a.number<T>("beta");
        a.done();
        Json out;
        out["alpha"] = encode_scalar(qsys::qsys_mmk_qed_alpha(beta));
        return out;
    }
};

struct OpQsysMmkQedStaffing {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const std::string crit = a.text("criterion", "delay");
        const T target = a.has("target") ? a.number<T>("target") : num_traits<T>::from_int(0);
        const T deadline = a.has("deadline") ? a.number<T>("deadline") : num_traits<T>::from_int(0);
        const T level = a.has("level") ? a.number<T>("level") : num_traits<T>::from_int(0);
        const bool exact = a.boolean("exact", false);
        a.done();
        qsys::QedCriterion c = qsys::QedCriterion::Delay;
        if (crit == "meanwait") {
            c = qsys::QedCriterion::MeanWait;
        } else if (crit == "servicelevel") {
            c = qsys::QedCriterion::ServiceLevel;
        } else if (crit != "delay") {
            throw InputError("qsys_mmk_qed_staffing: unknown criterion '" + crit + "'");
        }
        const qsys::QsysQedStaffingResult<T> r =
            qsys::qsys_mmk_qed_staffing(lambda, mu, target, c, deadline, level, exact);
        Json out;
        out["numServers"] = static_cast<double>(r.numServers);
        out["beta"] = encode_scalar(r.beta);
        out["betaTarget"] = encode_scalar(r.betaTarget);
        out["offeredLoad"] = encode_scalar(r.offeredLoad);
        out["probDelay"] = encode_scalar(r.probDelay);
        out["meanWait"] = encode_scalar(r.meanWait);
        if (c == qsys::QedCriterion::ServiceLevel)
            out["serviceLevel"] = encode_scalar(r.serviceLevel);
        return out;
    }
};

// The RQT family takes uncertainty-set variability parameters rather than
// coefficients of variation, and returns the worst case alongside the bound.
struct OpQsysGigkRqt {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T Gamma_a = a.number<T>("Gamma_a");
        const T Gamma_s = a.number<T>("Gamma_s");
        const unsigned k = a.uinteger("k");
        const T alpha_a = a.number<T>("alpha_a");
        const T alpha_s = a.number<T>("alpha_s");
        a.done();
        const qsys::GigkRqtResult<T> r =
            qsys::qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, k, alpha_a, alpha_s);
        Json out;
        out["W"] = encode_scalar(r.W);
        out["rhohat"] = encode_scalar(r.rhohat);
        out["Sworst"] = encode_scalar(r.Sworst);
        return out;
    }
};

struct OpQsysGig1Rqt {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T Gamma_a = a.number<T>("Gamma_a");
        const T Gamma_s = a.number<T>("Gamma_s");
        const T alpha_a = a.number<T>("alpha_a");
        const T alpha_s = a.number<T>("alpha_s");
        a.done();
        const qsys::GigkRqtResult<T> r =
            qsys::qsys_gig1_rqt(lambda, mu, Gamma_a, Gamma_s, alpha_a, alpha_s);
        Json out;
        out["W"] = encode_scalar(r.W);
        out["rhohat"] = encode_scalar(r.rhohat);
        out["Sworst"] = encode_scalar(r.Sworst);
        return out;
    }
};

struct OpQsysGigkRqtGamma {
    template <class T>
    static Json run(Args& a) {
        const T rho = a.number<T>("rho");
        const T mu = a.number<T>("mu");
        const T Gamma_a = a.number<T>("Gamma_a");
        const T sigma_s = a.number<T>("sigma_s");
        const unsigned k = a.uinteger("k");
        const T alpha_a = a.number<T>("alpha_a");
        const std::string regime = a.text("regime", "independent");
        a.done();
        Json out;
        out["Gamma_s"] =
            encode_scalar(qsys::qsys_gigk_rqt_gamma(rho, mu, Gamma_a, sigma_s, k, alpha_a, regime));
        return out;
    }
};

#define LINE_QSYS_MYSKJA(OpName, fn)                                            \
    struct OpName {                                                             \
        template <class T>                                                      \
        static Json run(Args& a) {                                              \
            const T lambda = a.number<T>("lambda");                             \
            const T mu = a.number<T>("mu");                                     \
            const T ca = a.number<T>("ca");                                     \
            const T cs = a.number<T>("cs");                                     \
            const T q0 = a.number<T>("q0");                                     \
            const T qa = a.number<T>("qa");                                     \
            a.done();                                                           \
            return qsys_results(qsys::fn(lambda, mu, ca, cs, q0, qa));          \
        }                                                                       \
    }

LINE_QSYS_MYSKJA(OpQsysGig1Myskja, qsys_gig1_approx_myskja);
LINE_QSYS_MYSKJA(OpQsysGig1Myskja2, qsys_gig1_approx_myskja2);
#undef LINE_QSYS_MYSKJA

struct OpQsysMm1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        a.done();
        return qsys_results(qsys::qsys_mm1(lambda, mu));
    }
};

struct OpQsysMmk {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned k = a.uinteger("k");
        a.done();
        return qsys_results(qsys::qsys_mmk(lambda, mu, k));
    }
};

struct OpQsysMg1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T cs = a.number<T>("cs");
        a.done();
        return qsys_results(qsys::qsys_mg1(lambda, mu, cs));
    }
};

struct OpQsysGg1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T ca2 = a.number<T>("ca2");
        const T cs2 = a.number<T>("cs2");
        a.done();
        return qsys_results(qsys::qsys_gg1(lambda, mu, ca2, cs2));
    }
};

/** The one member of the family returning a bare scalar, MATLAB's `W`. */
struct OpQsysGm1 {
    template <class T>
    static Json run(Args& a) {
        const T sigma = a.number<T>("sigma");
        const T mu = a.number<T>("mu");
        a.done();
        Json out;
        out["W"] = encode_scalar(qsys::qsys_gm1(sigma, mu));
        return out;
    }
};

struct OpQsysMginf {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const bool has_k = a.has("k");
        const unsigned k = has_k ? a.uinteger("k") : 0u;
        a.done();
        const qsys::MginfResult<T> r =
            has_k ? qsys::qsys_mginf(lambda, mu, k) : qsys::qsys_mginf(lambda, mu);
        Json out;
        out["L"] = encode_scalar(r.L);
        out["Lq"] = encode_scalar(r.Lq);
        out["W"] = encode_scalar(r.W);
        out["Wq"] = encode_scalar(r.Wq);
        out["p0"] = encode_scalar(r.p0);
        // pk IS REPORTED ONLY WHEN IT EXISTS: the reference returns it for the
        // k the caller named and leaves it undefined otherwise, and a zero
        // there would read as a probability rather than as the absence of one.
        if (r.has_pk) out["pk"] = encode_scalar(r.pk);
        return out;
    }
};

struct OpQsysMmck {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned c = a.uinteger("c");
        const unsigned K = a.uinteger("K");
        a.done();
        const qsys::MmckResult<T> r = qsys::qsys_mmck(lambda, mu, c, K);
        Json out;
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanQueueLengthQ"] = encode_scalar(r.meanQueueLengthQ);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["utilization"] = encode_scalar(r.utilization);
        out["throughput"] = encode_scalar(r.throughput);
        out["lossProbability"] = encode_scalar(r.lossProbability);
        out["queueLengthDist"] = encode_vector(r.queueLengthDist);
        return out;
    }
};

struct OpQsysMm1kLoss {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned K = a.uinteger("K");
        a.done();
        const qsys::Mm1kLossResult<T> r = qsys::qsys_mm1k_loss(lambda, mu, K);
        Json out;
        out["lossProbability"] = encode_scalar(r.lossProbability);
        out["utilization"] = encode_scalar(r.utilization);
        return out;
    }
};

struct OpQsysMg1kLossMgs {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const T mu_scv = a.number<T>("mu_scv");
        const unsigned K = a.uinteger("K");
        a.done();
        const qsys::Mg1kLossMgsResult<T> r = qsys::qsys_mg1k_loss_mgs(lambda, mu, mu_scv, K);
        Json out;
        out["lossProbability"] = encode_scalar(r.lossProbability);
        out["utilization"] = encode_scalar(r.utilization);
        return out;
    }
};

/**
 * The batch-arrival M[X]/M/1. MATLAB's fourth argument is
 * `E_X2_or_Var_X`, whose reading depends on a trailing 'variance' flag; the two
 * readings are separate C++ overloads, so the boundary takes the two names
 * apart and refuses both at once rather than guessing which was meant.
 */
struct OpQsysMxm1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda_batch = a.number<T>("lambda_batch");
        const T mu = a.number<T>("mu");
        const T E_X = a.number<T>("E_X");
        const bool has_m2 = a.has("E_X2"), has_var = a.has("Var_X");
        if (has_m2 == has_var)
            throw InputError(
                "qsys_mxm1: give exactly one of 'E_X2' (the second moment of the batch size) "
                "and 'Var_X' (its variance); MATLAB selects between them with a trailing "
                "'variance' flag, which has no place in a named-argument object");
        const T m2 = has_m2 ? a.number<T>("E_X2")
                            : T(a.number<T>("Var_X") + E_X * E_X);
        a.done();
        const qsys::MxM1Result<T> r = qsys::qsys_mxm1(lambda_batch, mu, E_X, m2);
        Json out;
        out["W"] = encode_scalar(r.W);
        out["Wq"] = encode_scalar(r.Wq);
        out["U"] = encode_scalar(r.U);
        out["Q"] = encode_scalar(r.Q);
        return out;
    }
};

struct OpQsysDmc {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda_arr");
        const T mu = a.number<T>("mu");
        const unsigned c = a.uinteger("c");
        const bool tuned = a.has("truncation") || a.has("quadSteps");
        const unsigned truncation = a.uinteger("truncation", 0u);
        const unsigned quadSteps = a.uinteger("quadSteps", 200u);
        a.done();
        const qsys::DmcResult<T> r = tuned ? qsys::qsys_dmc(lambda, mu, c, truncation, quadSteps)
                                           : qsys::qsys_dmc(lambda, mu, c);
        Json out;
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanWaitingQueue"] = encode_scalar(r.meanWaitingQueue);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["utilization"] = encode_scalar(r.utilization);
        return out;
    }
};

struct OpQsysMmccRetrialFp {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        const unsigned c = a.uinteger("c");
        const bool tuned = a.has("tol") || a.has("maxiter");
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        const std::size_t maxiter = a.count("maxiter", 0);
        a.done();
        const qsys::MmccRetrialFpResult<T> r =
            tuned ? qsys::qsys_mmcc_retrial_fp(lambda, mu, c, tol, maxiter)
                  : qsys::qsys_mmcc_retrial_fp(lambda, mu, c);
        Json out;
        out["blockingProbability"] = encode_scalar(r.blockingProbability);
        out["retrialRate"] = encode_scalar(r.retrialRate);
        out["iterations"] = encode_count(r.iterations);
        out["converged"] = r.converged;
        return out;
    }
};

struct OpQsysMm1Dps {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> lambda = a.vector<T>("lambda");
        const std::vector<T> mu = a.vector<T>("mu");
        const std::vector<T> w = a.vector<T>("w");
        const bool tuned = a.has("tol") || a.has("maxCutoff");
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        const unsigned maxCutoff = a.uinteger("maxCutoff", 0u);
        a.done();
        const qsys::Mm1DpsResult<T> r = tuned
                                            ? qsys::qsys_mm1_dps(lambda, mu, w, tol, maxCutoff)
                                            : qsys::qsys_mm1_dps(lambda, mu, w);
        Json out;
        // The port names the per-class response times `T_` because `T` is the
        // arithmetic parameter; MATLAB's output is `T`, which is the name a
        // caller reading the reference will look for, so that is the wire key.
        out["T"] = encode_vector(r.T_);
        out["rho"] = encode_scalar(r.rho);
        return out;
    }
};

/** The six M/G/1 disciplines, all {W per class, rhohat}. */
#define LINE_QSYS_MG1_DISC(OpName, fn, ResultT)                                 \
    struct OpName {                                                             \
        template <class T>                                                      \
        static Json run(Args& a) {                                              \
            const std::vector<T> lambda = a.vector<T>("lambda");                \
            const std::vector<T> mu = a.vector<T>("mu");                        \
            const std::vector<T> cs = a.vector<T>("cs");                        \
            a.done();                                                           \
            const qsys::ResultT<T> r = qsys::fn(lambda, mu, cs);                \
            Json out;                                                           \
            out["W"] = encode_vector(r.W);                                      \
            out["rhohat"] = encode_scalar(r.rhohat);                            \
            return out;                                                         \
        }                                                                       \
    }

LINE_QSYS_MG1_DISC(OpQsysMg1Fb, qsys_mg1_fb, Mg1DisciplineResult);
LINE_QSYS_MG1_DISC(OpQsysMg1Lrpt, qsys_mg1_lrpt, Mg1DisciplineResult);
LINE_QSYS_MG1_DISC(OpQsysMg1Psjf, qsys_mg1_psjf, Mg1DisciplineResult);
LINE_QSYS_MG1_DISC(OpQsysMg1Setf, qsys_mg1_setf, Mg1DisciplineResult);
LINE_QSYS_MG1_DISC(OpQsysMg1Srpt, qsys_mg1_srpt, Mg1DisciplineResult);
LINE_QSYS_MG1_DISC(OpQsysMg1Prio, qsys_mg1_prio, Mg1PrioResult);
#undef LINE_QSYS_MG1_DISC

/**
 * The two slotted systems. `convention` selects LAS-DA (late arrival, delayed
 * access; the reference's default) or EAS (early arrival). They are DIFFERENT
 * queues, not two readings of one, so an unknown name is refused rather than
 * falling back to the default.
 */
inline dqsys::GeoConvention read_geo_convention(Args& a) {
    const std::string c = a.text("convention", "LAS_DA");
    if (c == "LAS_DA" || c == "las_da" || c == "LASDA") return dqsys::GeoConvention::LAS_DA;
    if (c == "EAS" || c == "eas") return dqsys::GeoConvention::EAS;
    throw InputError("'convention' must be LAS_DA or EAS, got '" + c + "'");
}

inline const char* geo_convention_name(dqsys::GeoConvention c) {
    return c == dqsys::GeoConvention::EAS ? "EAS" : "LAS_DA";
}

struct OpQsysGeoGeo1 {
    template <class T>
    static Json run(Args& a) {
        const T arrival = a.number<T>("a");
        const T service = a.number<T>("s");
        const dqsys::GeoConvention conv = read_geo_convention(a);
        a.done();
        const dqsys::GeoGeo1Result<T> r = dqsys::dqsys_geogeo1(arrival, service, conv);
        Json out;
        out["convention"] = geo_convention_name(r.convention);
        out["arrivalProb"] = encode_scalar(r.arrivalProb);
        out["serviceProb"] = encode_scalar(r.serviceProb);
        out["utilization"] = encode_scalar(r.utilization);
        out["throughput"] = encode_scalar(r.throughput);
        out["emptyProb"] = encode_scalar(r.emptyProb);
        out["ratio"] = encode_scalar(r.ratio);
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanWaitingQueue"] = encode_scalar(r.meanWaitingQueue);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanServiceTime"] = encode_scalar(r.meanServiceTime);
        return out;
    }
};

struct OpQsysGeoxGeo1 {
    template <class T>
    static Json run(Args& a) {
        const T arrival = a.number<T>("a");
        const T beta = a.number<T>("beta");
        const T service = a.number<T>("s");
        const dqsys::GeoConvention conv = read_geo_convention(a);
        a.done();
        const dqsys::GeoXGeo1Result<T> r = dqsys::dqsys_geoxgeo1(arrival, beta, service, conv);
        Json out;
        out["convention"] = geo_convention_name(r.convention);
        out["batchArrivalProb"] = encode_scalar(r.batchArrivalProb);
        out["batchMean"] = encode_scalar(r.batchMean);
        out["batchSecondFactorialMoment"] = encode_scalar(r.batchSecondFactorialMoment);
        out["serviceProb"] = encode_scalar(r.serviceProb);
        out["arrivalRate"] = encode_scalar(r.arrivalRate);
        out["throughput"] = encode_scalar(r.throughput);
        out["utilization"] = encode_scalar(r.utilization);
        out["boundaryEmptyProb"] = encode_scalar(r.boundaryEmptyProb);
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanWaitingQueue"] = encode_scalar(r.meanWaitingQueue);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanServiceTime"] = encode_scalar(r.meanServiceTime);
        return out;
    }
};

/** The four means plus the level distribution, shared by the MAP-driven queues. */
template <class R>
Json qsys_map_results(const R& r) {
    Json out;
    out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
    out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
    out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
    out["utilization"] = encode_scalar(r.utilization);
    out["queueLengthDist"] = encode_vector(r.queueLengthDist);
    return out;
}

/** A MAP as the (D0, D1) pair the reference passes around, under any two keys. */
template <class T>
mam::Map<T> read_map(Args& a, const char* d0, const char* d1) {
    mam::Map<T> m;
    m.D0 = a.matrix<T>(d0);
    m.D1 = a.matrix<T>(d1);
    return m;
}

struct OpQsysMapm1 {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const T mu = a.number<T>("mu");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        return qsys_map_results(sized ? qsys::qsys_mapm1(arrival, mu, dist_size)
                                      : qsys::qsys_mapm1(arrival, mu));
    }
};

struct OpQsysMapmc {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const T mu = a.number<T>("mu");
        const unsigned c = a.uinteger("c");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        return qsys_map_results(sized ? qsys::qsys_mapmc(arrival, mu, c, dist_size)
                                      : qsys::qsys_mapmc(arrival, mu, c));
    }
};

struct OpQsysMapmap1 {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "C0", "C1");
        const mam::Map<T> service = read_map<T>(a, "D0", "D1");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        return qsys_map_results(sized ? qsys::qsys_mapmap1(arrival, service, dist_size)
                                      : qsys::qsys_mapmap1(arrival, service));
    }
};

struct OpQsysMapph1 {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const std::vector<T> sigma = a.vector<T>("sigma");
        const Matrix<T> S = a.matrix<T>("S");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        return qsys_map_results(sized ? qsys::qsys_mapph1(arrival, sigma, S, dist_size)
                                      : qsys::qsys_mapph1(arrival, sigma, S));
    }
};

/**
 * MAP/PH/c. The port offers exactly two overloads, the defaulted one and the
 * one taking all three tuning arguments, so naming ANY of them requires naming
 * ALL of them rather than mixing a caller value with defaults restated here.
 */
struct OpQsysMapphc {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const std::vector<T> alpha = a.vector<T>("alpha");
        const Matrix<T> S = a.matrix<T>("S");
        const unsigned c = a.uinteger("c");
        const bool any =
            a.has("dist_size") || a.has("num_w_moms") || a.has("w_points");
        const bool all =
            a.has("dist_size") && a.has("num_w_moms") && a.has("w_points");
        if (any && !all)
            throw InputError(
                "qsys_mapphc: 'dist_size', 'num_w_moms' and 'w_points' select one tuned "
                "overload together; give all three or none");
        const std::size_t dist_size = a.count("dist_size", 0);
        const std::size_t num_w_moms = a.count("num_w_moms", 0);
        const std::vector<T> w_points = a.vector_or_empty<T>("w_points");
        a.done();
        const qsys::MapPhcResult<T> r =
            all ? qsys::qsys_mapphc(arrival, alpha, S, c, dist_size, num_w_moms, w_points)
                : qsys::qsys_mapphc(arrival, alpha, S, c);
        Json out = qsys_map_results(r);
        out["waitingTimeMoments"] = encode_vector(r.waitingTimeMoments);
        out["waitingTimeCCDF"] = encode_vector(r.waitingTimeCCDF);
        out["waitingTimePoints"] = encode_vector(r.waitingTimePoints);
        out["probWait"] = encode_scalar(r.probWait);
        // The repeating configuration count, binomial(ms+c-1,c): it is what
        // sizes the solve, so a caller comparing runtimes needs it reported
        // rather than recomputed from ms and c.
        out["phaseCount"] = encode_count(r.phaseCount);
        return out;
    }
};

struct OpQsysPhph1 {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> alpha = a.vector<T>("alpha");
        const Matrix<T> Tm = a.matrix<T>("T");
        const std::vector<T> beta = a.vector<T>("beta");
        const Matrix<T> S = a.matrix<T>("S");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        return qsys_map_results(sized ? qsys::qsys_phph1(alpha, Tm, beta, S, dist_size)
                                      : qsys::qsys_phph1(alpha, Tm, beta, S));
    }
};

/**
 * The two deterministic-service MAP queues. Their four tuning arguments move
 * together -- the arrival truncation bounds a sum the level iteration then
 * inverts -- so naming ANY of them requires naming ALL of them, rather than
 * mixing one caller value with three library defaults from a different regime.
 */
inline bool qsys_map_d_tuned(Args& a, const char* who) {
    const bool any = a.has("dist_size") || a.has("max_arrivals") || a.has("max_levels") ||
                     a.has("tol");
    const bool all = a.has("dist_size") && a.has("max_arrivals") && a.has("max_levels") &&
                     a.has("tol");
    if (any && !all)
        throw InputError(std::string(who) +
                         ": 'dist_size', 'max_arrivals', 'max_levels' and 'tol' tune one "
                         "truncation together; give all four or none");
    return all;
}

struct OpQsysMapd1 {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const T s = a.number<T>("s");
        const bool tuned = qsys_map_d_tuned(a, "qsys_mapd1");
        const std::size_t dist_size = a.count("dist_size", 0);
        const unsigned max_arrivals = a.uinteger("max_arrivals", 0u);
        const std::size_t max_levels = a.count("max_levels", 0);
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        a.done();
        return qsys_map_results(
            tuned ? qsys::qsys_mapd1(arrival, s, dist_size, max_arrivals, max_levels, tol)
                  : qsys::qsys_mapd1(arrival, s));
    }
};

struct OpQsysMapdc {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const T s = a.number<T>("s");
        const unsigned c = a.uinteger("c");
        const bool tuned = qsys_map_d_tuned(a, "qsys_mapdc");
        const std::size_t dist_size = a.count("dist_size", 0);
        const unsigned max_arrivals = a.uinteger("max_arrivals", 0u);
        const std::size_t max_levels = a.count("max_levels", 0);
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        a.done();
        return qsys_map_results(
            tuned ? qsys::qsys_mapdc(arrival, s, c, dist_size, max_arrivals, max_levels, tol)
                  : qsys::qsys_mapdc(arrival, s, c));
    }
};

struct OpQsysMapg1 {
    template <class T>
    static Json run(Args& a) {
        const mam::Map<T> arrival = read_map<T>(a, "D0", "D1");
        const std::vector<T> moments = a.vector<T>("serviceMoments");
        const bool sized = a.has("dist_size");
        const std::size_t dist_size = a.count("dist_size", 0);
        a.done();
        const qsys::MapG1Result<T> r = sized ? qsys::qsys_mapg1(arrival, moments, dist_size)
                                             : qsys::qsys_mapg1(arrival, moments);
        Json out = qsys_map_results(r);
        // WHICH LAW WAS FITTED IS PART OF THE ANSWER: the moments select an
        // exponential, an Erlang, a hyperexponential or an acyclic PH, and a
        // caller comparing against MATLAB has to know which branch produced the
        // number before it can call a difference a discrepancy.
        const char* kind = "Acyclic";
        switch (r.fitKind) {
            case qsys::MapG1ServiceFit::Exponential: kind = "Exponential"; break;
            case qsys::MapG1ServiceFit::Erlang: kind = "Erlang"; break;
            case qsys::MapG1ServiceFit::Hyperexponential: kind = "Hyperexponential"; break;
            case qsys::MapG1ServiceFit::Acyclic: kind = "Acyclic"; break;
        }
        out["fitKind"] = kind;
        out["servicePhases"] = encode_count(r.servicePhases);
        out["serviceFitD0"] = encode_matrix(r.serviceFit.D0);
        out["serviceFitD1"] = encode_matrix(r.serviceFit.D1);
        return out;
    }
};

struct OpQsysPhm1 {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> alpha = a.vector<T>("alpha");
        const Matrix<T> Tm = a.matrix<T>("T");
        const T mu = a.number<T>("mu");
        const bool tuned = a.has("tol");
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        a.done();
        const qsys::PhM1Result<T> r = tuned ? qsys::qsys_phm1(alpha, Tm, mu, tol)
                                            : qsys::qsys_phm1(alpha, Tm, mu);
        Json out;
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanWaitingQueue"] = encode_scalar(r.meanWaitingQueue);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["utilization"] = encode_scalar(r.utilization);
        out["sigma"] = encode_scalar(r.sigma);
        return out;
    }
};

struct OpQsysPhmc {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> alpha = a.vector<T>("alpha");
        const Matrix<T> Tm = a.matrix<T>("T");
        const T mu = a.number<T>("mu");
        const unsigned c = a.uinteger("c");
        const bool any = a.has("maxIter") || a.has("tol");
        const bool all = a.has("maxIter") && a.has("tol");
        if (any && !all)
            throw InputError(
                "qsys_phmc: 'maxIter' and 'tol' bound one fixed-point iteration together; give "
                "both or neither");
        const unsigned maxIter = a.uinteger("maxIter", 0u);
        const T tol = a.scalar<T>("tol", num_traits<T>::from_int(0));
        a.done();
        const qsys::PhMcResult<T> r = all ? qsys::qsys_phmc(alpha, Tm, mu, c, maxIter, tol)
                                          : qsys::qsys_phmc(alpha, Tm, mu, c);
        Json out;
        out["meanQueueLength"] = encode_scalar(r.meanQueueLength);
        out["meanWaitingQueue"] = encode_scalar(r.meanWaitingQueue);
        out["meanWaitingTime"] = encode_scalar(r.meanWaitingTime);
        out["meanSojournTime"] = encode_scalar(r.meanSojournTime);
        out["utilization"] = encode_scalar(r.utilization);
        return out;
    }
};

struct OpQsysBmapphnnRetrial {
    template <class T>
    static Json run(Args& a) {
        const std::vector<Matrix<T> > D = a.matrices<T>("D");
        const std::vector<T> beta = a.vector<T>("beta");
        const Matrix<T> S = a.matrix<T>("S");
        const int N = a.required_integer("N");
        const T alpha = a.number<T>("alpha");
        const T gamma = a.number<T>("gamma");
        const T p = a.number<T>("p");
        const std::vector<int> Ri = a.ints("R");
        qsys::BmapPhNnRetrialOptions opt;
        opt.maxLevel = a.count("maxLevel", 0);
        a.done();
        std::vector<long> R(Ri.size());
        for (std::size_t i = 0; i < Ri.size(); ++i) R[i] = static_cast<long>(Ri[i]);
        const qsys::BmapPhNnRetrialResult<T> r =
            qsys::qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, opt);
        Json out;
        out["L_orbit"] = encode_scalar(r.L_orbit);
        out["N_server"] = encode_scalar(r.N_server);
        out["L_system"] = encode_scalar(r.L_system);
        out["utilization"] = encode_scalar(r.utilization);
        out["throughput"] = encode_scalar(r.throughput);
        out["P_idle"] = encode_scalar(r.P_idle);
        out["P_empty_orbit"] = encode_scalar(r.P_empty_orbit);
        out["P_empty_system"] = encode_scalar(r.P_empty_system);
        out["pi"] = encode_matrix(r.pi);
        out["truncLevel"] = encode_count(r.truncLevel);
        out["topLevelMass"] = encode_scalar(r.topLevelMass);
        // WHETHER THE LEVEL WAS CLIPPED is part of the answer: a clipped
        // truncation makes every reported mass a lower bound on the true one.
        out["clipped"] = r.clipped;
        return out;
    }
};

// ---------------------------------------------------------------------------
// moment: the moment-basis changes
// ---------------------------------------------------------------------------
//
// Every member is a pure change of basis on a moment sequence -- raw, central,
// factorial, up-factorial, binomial, negative-binomial -- so each takes one
// vector and returns one vector under the reference's own argument and output
// names. They are field operations throughout, which is why the registry lists
// all three arithmetics and none of them needs the transcendental guard.

#define LINE_MOMENT_VEC(OpName, fn, in_key, out_key)                     \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const std::vector<T> v = a.vector<T>(in_key);                \
            a.done();                                                    \
            Json out;                                                    \
            out[out_key] = encode_vector(moment::fn(v));                 \
            return out;                                                  \
        }                                                                \
    }

LINE_MOMENT_VEC(OpMomentBinomialFromFactorial, moment_binomial_from_factorial, "f", "b");
LINE_MOMENT_VEC(OpMomentBinomialFromNegbinomial, moment_binomial_from_negbinomial, "bm", "b");
LINE_MOMENT_VEC(OpMomentBinotrans, moment_binotrans, "x", "y");
LINE_MOMENT_VEC(OpMomentBinotransinv, moment_binotransinv, "y", "x");
LINE_MOMENT_VEC(OpMomentCentralFromRaw, moment_central_from_raw, "m", "mc");
LINE_MOMENT_VEC(OpMomentFactorialFromBinomial, moment_factorial_from_binomial, "b", "f");
LINE_MOMENT_VEC(OpMomentFactorialFromRaw, moment_factorial_from_raw, "m", "f");
LINE_MOMENT_VEC(OpMomentFactorialFromUpfactorial, moment_factorial_from_upfactorial, "fp", "f");
LINE_MOMENT_VEC(OpMomentNegbinomialFromBinomial, moment_negbinomial_from_binomial, "b", "bm");
LINE_MOMENT_VEC(OpMomentNegbinomialFromUpfactorial, moment_negbinomial_from_upfactorial, "fp",
                "bm");
LINE_MOMENT_VEC(OpMomentRawFromFactorial, moment_raw_from_factorial, "f", "m");
LINE_MOMENT_VEC(OpMomentRawFromUpfactorial, moment_raw_from_upfactorial, "fp", "m");
LINE_MOMENT_VEC(OpMomentUpfactorialFromFactorial, moment_upfactorial_from_factorial, "f", "fp");
LINE_MOMENT_VEC(OpMomentUpfactorialFromNegbinomial, moment_upfactorial_from_negbinomial, "bm",
                "fp");
LINE_MOMENT_VEC(OpMomentUpfactorialFromRaw, moment_upfactorial_from_raw, "m", "fp");
#undef LINE_MOMENT_VEC

/**
 * The inverse of moment_central_from_raw needs the MEAN as well: a central
 * sequence has lost it (mc(1) is 0 by construction), so it cannot be recovered
 * from the sequence alone and is a second argument rather than an option.
 */
struct OpMomentRawFromCentral {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> mc = a.vector<T>("mc");
        const T m1 = a.number<T>("m1");
        a.done();
        Json out;
        out["m"] = encode_vector(moment::moment_raw_from_central(mc, m1));
        return out;
    }
};

/** The four triangular coefficient matrices, each a function of the order n. */
#define LINE_MOMENT_TRI(OpName, fn, out_key)                             \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const int n = a.required_integer("n");                       \
            a.done();                                                    \
            Json out;                                                    \
            out[out_key] = encode_matrix(moment::fn<T>(n));              \
            return out;                                                  \
        }                                                                \
    }

LINE_MOMENT_TRI(OpMomentLah, moment_lah, "L");
LINE_MOMENT_TRI(OpMomentStirling1, moment_stirling1, "s");
LINE_MOMENT_TRI(OpMomentStirling2, moment_stirling2, "S");
LINE_MOMENT_TRI(OpMomentStirlingcycle, moment_stirlingcycle, "sigma");
#undef LINE_MOMENT_TRI

// ---------------------------------------------------------------------------
// mc: the Markov-chain solvers beyond ctmc_solve / dtmc_solve
// ---------------------------------------------------------------------------
//
// `ctmc_rand` is NOT exposed: it takes a random generator by reference, so its
// answer depends on a stream this boundary cannot carry, and a generator seeded
// here would return a matrix the caller cannot reproduce.
//
// AN ITERATIVE SOLVER REPORTS WHETHER IT CONVERGED, and that flag rides beside
// the vector rather than being folded into it. A caller that reads only the
// numbers from a run that hit its iteration cap would be reading a partial
// sweep as a stationary law.

/** The block partition MS/MSS: a list of state-index lists, 0-based. */
inline std::vector<std::vector<std::size_t> > read_partition(Args& a, const char* key) {
    const Json& j = a.get(key);
    const std::string where = std::string("partition '") + key + "'";
    if (!j.is_array()) throw InputError(where + " must be an array of index arrays");
    std::vector<std::vector<std::size_t> > out;
    for (std::size_t b = 0; b < j.size(); ++b) {
        if (!j[b].is_array())
            throw InputError(where + ": block " + std::to_string(b) + " must be an array");
        std::vector<std::size_t> blk;
        for (std::size_t i = 0; i < j[b].size(); ++i) {
            const std::string text = decimal_text(j[b][i], where);
            const long long v = std::strtoll(text.c_str(), nullptr, 10);
            if (v < 0)
                throw InputError(where + ": a state index cannot be negative, got " + text);
            blk.push_back(static_cast<std::size_t>(v));
        }
        out.push_back(blk);
    }
    return out;
}

/** A list of state indices, reported unchanged: an index is not a measurement. */
inline Json encode_indices(const std::vector<std::size_t>& v) {
    Json a = Json::array();
    for (std::size_t i = 0; i < v.size(); ++i) a.push_back(static_cast<std::uint64_t>(v[i]));
    return a;
}

/**
 * The reducible solvers' five outputs. `pis` is the per-component stationary
 * law and `pi0` the absorption probabilities into each component, so the two
 * together say WHERE the mass went as well as how it is spread once there;
 * reporting only `pi` would lose the decomposition the function exists for.
 */
template <class T>
Json reducible_results(const mc::ReducibleResult<T>& r) {
    Json out;
    out["pi"] = encode_vector(r.pi);
    out["pis"] = encode_matrix(r.pis);
    out["pi0"] = encode_matrix(r.pi0);
    out["scc"] = encode_indices(r.scc);
    Json isrec = Json::array();
    for (std::size_t i = 0; i < r.isrec.size(); ++i) isrec.push_back(bool(r.isrec[i]));
    out["isrec"] = isrec;
    out["Pl"] = encode_matrix(r.Pl);
    out["pil"] = encode_matrix(r.pil);
    return out;
}

/** The block-decomposition variant, which has no lumped chain to report. */
template <class T>
Json blkdecomp_results(const mc::BlkDecompResult<T>& r) {
    Json out;
    out["pi"] = encode_vector(r.pi);
    out["pis"] = encode_matrix(r.pis);
    out["pi0"] = encode_matrix(r.pi0);
    out["scc"] = encode_indices(r.scc);
    Json isrec = Json::array();
    for (std::size_t i = 0; i < r.isrec.size(); ++i) isrec.push_back(bool(r.isrec[i]));
    out["isrec"] = isrec;
    return out;
}

struct OpCtmcSolveReducible {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<T> pi0 = a.vector<T>("pi0");
        const double zeroColTol = a.scalar<double>("zeroColTol", 1e-12);
        a.done();
        const mc::ReducibleResult<T> r = mc::ctmc_solve_reducible(Q, pi0, zeroColTol);
        return reducible_results(r);
    }
};

struct OpCtmcSolveReducibleBlkdecomp {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<T> pin = a.vector<T>("pin");
        const double reachTol = a.scalar<double>("reachTol", 1e-15);
        const double zeroColTol = a.scalar<double>("zeroColTol", 1e-12);
        a.done();
        const mc::BlkDecompResult<T> r =
            mc::ctmc_solve_reducible_blkdecomp(Q, pin, reachTol, zeroColTol);
        return blkdecomp_results(r);
    }
};

struct OpDtmcSolveReducible {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> P = a.matrix<T>("P");
        const std::vector<T> pin = a.vector<T>("pin");
        const double zeroColTol = a.scalar<double>("zeroColTol", 1e-12);
        a.done();
        const mc::ReducibleResult<T> r = mc::dtmc_solve_reducible(P, pin, zeroColTol);
        return reducible_results(r);
    }
};

struct OpDtmcMakestochastic {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> P = a.matrix<T>("Pin");
        a.done();
        Json out;
        out["P"] = encode_matrix(mc::dtmc_makestochastic(P));
        return out;
    }
};

struct OpCtmcSens {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const Matrix<T> dQ = a.matrix<T>("dQ");
        const std::vector<T> pi = a.vector<T>("pi");
        a.done();
        Json out;
        out["dpi"] = encode_vector(mc::ctmc_sens(Q, dQ, pi));
        return out;
    }
};

struct OpStronglyConnComp {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> A = a.matrix<T>("A");
        a.done();
        const mc::SccResult r = mc::stronglyconncomp(A);
        Json out;
        out["scc"] = encode_indices(r.scc);
        Json rec = Json::array();
        for (std::size_t i = 0; i < r.recurrent.size(); ++i) rec.push_back(bool(r.recurrent[i]));
        out["recurrent"] = rec;
        Json mem = Json::array();
        for (std::size_t i = 0; i < r.members.size(); ++i)
            mem.push_back(encode_indices(r.members[i]));
        out["members"] = mem;
        return out;
    }
};

struct OpCtmcRandomization {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const bool given = a.has("q");
        const T q = a.scalar<T>("q", num_traits<T>::from_int(0));
        a.done();
        if (!given)
            throw InputError(
                "ctmc_randomization: 'q' is the uniformization rate and has no default here; "
                "MATLAB requires it too");
        const mc::RandomizationResult<T> r = mc::ctmc_randomization(Q, q);
        Json out;
        out["P"] = encode_matrix(r.P);
        out["q"] = encode_scalar(r.q);
        return out;
    }
};

/**
 * The three nearly-completely-decomposable solvers. `q` is the uniformization
 * rate and is OPTIONAL for courtois and multi, because the port carries the
 * reference's own derived default (courtois_default_rate) and restating it here
 * would fix a value the library computes from Q and MS.
 */
struct OpCtmcCourtois {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<std::vector<std::size_t> > MS = read_partition(a, "MS");
        const bool given = a.has("q");
        const T q = a.scalar<T>("q", num_traits<T>::from_int(0));
        a.done();
        const mc::CourtoisResult<T> r =
            given ? mc::ctmc_courtois(Q, MS, q) : mc::ctmc_courtois(Q, MS);
        Json out;
        out["p"] = encode_vector(r.p);
        out["v"] = encode_indices(r.v);
        out["Qperm"] = encode_matrix(r.Qperm);
        out["Qdec"] = encode_matrix(r.Qdec);
        out["P"] = encode_matrix(r.P);
        out["B"] = encode_matrix(r.B);
        out["C"] = encode_scalar(r.C);
        out["eps"] = encode_scalar(r.eps);
        out["epsRowMax"] = encode_scalar(r.epsRowMax);
        out["epsMAX"] = encode_scalar(r.epsMAX);
        out["q"] = encode_scalar(r.q);
        return out;
    }
};

/** KMS and Takahashi return the same six outputs, under the same names. */
template <class R>
Json kms_results(const R& r) {
    Json out;
    out["p"] = encode_vector(r.p);
    out["p_1"] = encode_vector(r.p_1);
    out["pcourt"] = encode_vector(r.pcourt);
    out["Qperm"] = encode_matrix(r.Qperm);
    out["eps"] = encode_scalar(r.eps);
    out["epsMAX"] = encode_scalar(r.epsMAX);
    return out;
}

struct OpCtmcKms {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<std::vector<std::size_t> > MS = read_partition(a, "MS");
        const std::size_t numSteps = a.count("numSteps");
        a.done();
        return kms_results(mc::ctmc_kms(Q, MS, numSteps));
    }
};

struct OpCtmcTakahashi {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<std::vector<std::size_t> > MS = read_partition(a, "MS");
        const std::size_t numSteps = a.count("numSteps");
        const double massTol = a.scalar<double>("massTol", 1e-14);
        a.done();
        return kms_results(mc::ctmc_takahashi(Q, MS, numSteps, massTol));
    }
};

struct OpCtmcMulti {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<std::vector<std::size_t> > MS = read_partition(a, "MS");
        const std::vector<std::vector<std::size_t> > MSS = read_partition(a, "MSS");
        const bool given = a.has("q");
        const T q = a.scalar<T>("q", num_traits<T>::from_int(0));
        a.done();
        const mc::MultiResult<T> r =
            given ? mc::ctmc_multi(Q, MS, MSS, q) : mc::ctmc_multi(Q, MS, MSS);
        Json out;
        out["p"] = encode_vector(r.p);
        out["pcourt"] = encode_vector(r.pcourt);
        out["Qperm"] = encode_matrix(r.Qperm);
        out["eps"] = encode_scalar(r.eps);
        out["epsMAX"] = encode_scalar(r.epsMAX);
        return out;
    }
};

/**
 * GMRES reports MATLAB's `flag` rather than a boolean: 0 is convergence and the
 * nonzero values distinguish an iteration cap from a stagnation, which are
 * different failures and are what a caller checks before using `x`.
 */
struct OpCtmcGmres {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> A = a.matrix<T>("A");
        const std::vector<T> b = a.vector<T>("b");
        const double tol = a.scalar<double>("tol", 1e-12);
        const int restart = a.integer("restart", 0);
        const int maxit = a.integer("maxit", 0);
        const std::vector<T> x0 = a.vector_or_empty<T>("x0");
        a.done();
        const mc::GmresResult<T> r = mc::ctmc_gmres(A, b, tol, restart, maxit, x0);
        Json out;
        out["x"] = encode_vector(r.x);
        out["flag"] = r.flag;
        out["relres"] = encode_scalar(r.relres);
        out["iter"] = static_cast<std::int64_t>(r.iter);
        return out;
    }
};

struct OpCtmcGmresMulti {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> A = a.matrix<T>("A");
        const Matrix<T> B = a.matrix<T>("B");
        const double tol = a.scalar<double>("tol", 1e-12);
        const int restart = a.integer("restart", 0);
        const int maxit = a.integer("maxit", 0);
        a.done();
        const mc::GmresMultiResult<T> r = mc::ctmc_gmres_multi(A, B, tol, restart, maxit);
        Json out;
        out["X"] = encode_matrix(r.X);
        out["flag"] = r.flag;
        return out;
    }
};

/**
 * BiCGSTAB reports MATLAB's `flag` on the same convention as GMRES, with the
 * addition of 4 for a scalar breakdown of the underlying Lanczos process, which
 * is a different failure from an iteration cap and is what a caller checks
 * before using `x`. `iter` counts matrix-vector products with A, so it is
 * directly comparable with the `iter` GMRES reports.
 */
struct OpCtmcBicgstab {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> A = a.matrix<T>("A");
        const std::vector<T> b = a.vector<T>("b");
        const double tol = a.scalar<double>("tol", 1e-12);
        const int maxit = a.integer("maxit", 0);
        const std::vector<T> x0 = a.vector_or_empty<T>("x0");
        a.done();
        const mc::BicgstabResult<T> r = mc::ctmc_bicgstab(A, b, tol, maxit, x0);
        Json out;
        out["x"] = encode_vector(r.x);
        out["flag"] = r.flag;
        out["relres"] = encode_scalar(r.relres);
        out["iter"] = static_cast<std::int64_t>(r.iter);
        return out;
    }
};

struct OpCtmcBicgstabMulti {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> A = a.matrix<T>("A");
        const Matrix<T> B = a.matrix<T>("B");
        const double tol = a.scalar<double>("tol", 1e-12);
        const int maxit = a.integer("maxit", 0);
        a.done();
        const mc::BicgstabMultiResult<T> r = mc::ctmc_bicgstab_multi(A, B, tol, maxit);
        Json out;
        out["X"] = encode_matrix(r.X);
        out["flag"] = r.flag;
        return out;
    }
};

struct OpCtmcFoxglynn {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> pi0 = a.vector<T>("pi0");
        const Matrix<T> Q = a.matrix<T>("Q");
        const T t = a.number<T>("t");
        const double tol = a.scalar<double>("tol", 1e-12);
        const int maxiter = a.integer("maxiter", -1);
        a.done();
        const mc::FoxGlynnResult<T> r = mc::ctmc_foxglynn(pi0, Q, t, tol, maxiter);
        Json out;
        out["pi"] = encode_vector(r.pi);
        // The Poisson weight window: `left` and `right` are the truncation
        // points the algorithm chose, and `w` the weights between them, so a
        // caller can see how much mass the truncation kept.
        out["left"] = static_cast<std::int64_t>(r.left);
        out["right"] = static_cast<std::int64_t>(r.right);
        out["w"] = encode_vector(r.w);
        return out;
    }
};

struct OpCtmcSaddlepoint {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> D0 = a.matrix<T>("D0");
        const Matrix<T> D1 = a.matrix<T>("D1");
        const std::vector<T> t = a.vector<T>("t");
        const std::vector<int> ki = a.ints("k");
        const std::string method = a.text("method", "daniels2");
        const std::vector<T> pi0 = a.vector_or_empty<T>("pi0");
        a.done();
        mc::SaddlepointMethod m = mc::SADDLEPOINT_DANIELS2;
        if (method == "daniels" || method == "sp1")
            m = mc::SADDLEPOINT_DANIELS;
        else if (method == "plain" || method == "bare")
            m = mc::SADDLEPOINT_PLAIN;
        else if (!(method == "daniels2" || method == "sp2"))
            throw InputError("ctmc_saddlepoint: unknown method '" + method +
                             "', expected daniels2, daniels or plain");
        std::vector<long> k(ki.size());
        for (std::size_t i = 0; i < ki.size(); ++i) k[i] = static_cast<long>(ki[i]);
        const mc::SaddlepointResult<T> r = mc::ctmc_saddlepoint(D0, D1, t, k, m, pi0);
        Json out;
        out["p"] = encode_vector(r.p);
        out["logp"] = encode_vector(r.logp);
        out["theta"] = encode_vector(r.theta);
        out["k2"] = encode_vector(r.k2);
        out["lambda"] = encode_scalar(r.lambda);
        // K2 = t*eta''(theta*) is the expansion parameter, so a caller that
        // ignores out_of_regime is reading a number outside its own validity
        out["outOfRegime"] = r.out_of_regime;
        return out;
    }
};

struct OpCtmcFau {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> pi0 = a.vector<T>("pi0");
        const Matrix<T> Q = a.matrix<T>("Q");
        const T t = a.number<T>("t");
        const double epsilon = a.scalar<double>("epsilon", 1e-6);
        const double delta = a.scalar<double>("delta", 1e-12);
        const int maxsteps = a.integer("maxsteps", -1);
        a.done();
        const mc::FauResult<T> r = mc::ctmc_fau(pi0, Q, t, epsilon, delta, maxsteps);
        Json out;
        out["pi"] = encode_vector(r.pit);
        // The whole error budget, since the method removes mass and never puts
        // any back: errorBound is the L1 distance to the exact distribution,
        // and the three components say where it went.
        out["steps"] = static_cast<std::int64_t>(r.steps);
        out["lambdaMin"] = r.lambdaMin;
        out["lambdaMax"] = r.lambdaMax;
        out["uniformRate"] = r.uniformRate;
        out["weightTail"] = num_traits<T>::to_double(r.weightTail);
        out["weightWindow"] = num_traits<T>::to_double(r.weightWindow);
        out["droppedMass"] = num_traits<T>::to_double(r.droppedMass);
        out["errorBound"] = num_traits<T>::to_double(r.errorBound);
        out["supportMax"] = static_cast<std::int64_t>(r.supportMax);
        out["supportFinal"] = static_cast<std::int64_t>(r.supportFinal);
        out["truncated"] = r.truncated;
        out["absorbed"] = r.absorbed;
        return out;
    }
};

struct OpCtmcTransient {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const std::vector<T> pi0 = a.vector<T>("pi0");
        const T t0 = a.number<T>("t0");
        const T t1 = a.number<T>("t1");
        const double rtol = a.scalar<double>("rtol", 1e-3);
        const double atol = a.scalar<double>("atol", 1e-6);
        a.done();
        const mc::TransientResult<T> r = mc::ctmc_transient(Q, pi0, t0, t1, rtol, atol);
        Json out;
        out["t"] = encode_vector(r.t);
        out["pi"] = encode_matrix(r.pi);
        return out;
    }
};

struct OpCtmcTransientSens {
    template <class T>
    static Json run(Args& a) {
        const Matrix<T> Q = a.matrix<T>("Q");
        const Matrix<T> dQ = a.matrix<T>("dQ");
        const std::vector<T> pi0 = a.vector<T>("pi0");
        const T t0 = a.number<T>("t0");
        const T t1 = a.number<T>("t1");
        const double rtol = a.scalar<double>("rtol", 1e-3);
        const double atol = a.scalar<double>("atol", 1e-6);
        a.done();
        const mc::TransientSensResult<T> r =
            mc::ctmc_transient_sens(Q, dQ, pi0, t0, t1, rtol, atol);
        Json out;
        out["t"] = encode_vector(r.t);
        out["pi"] = encode_matrix(r.pi);
        out["dpi"] = encode_matrix(r.dpi);
        return out;
    }
};

// ---------------------------------------------------------------------------
// aoi: the age-of-information laws
// ---------------------------------------------------------------------------
//
// The MATLAB signatures take the service law as a LAPLACE-STIELTJES TRANSFORM,
// which in C++ is a std::function and has no JSON representation. Rather than
// leave nine of the fifteen laws unreachable, the boundary takes a NAMED law
// and builds the transform with the port's own aoi_lst_* constructors -- the
// same four the reference offers -- so the caller states which law it means and
// nothing about the transform is invented here.
//
// `lstAoI` IS NOT REPORTED, because it is a transform and not a number; the
// header states that nothing in the family inverts one. `has_lst` says whether
// the law returned one at all, which is the part a caller can act on.

/** The service law of the G/M/1 and M/G/1 age formulas, by name. */
template <class T>
aoi::Lst<T> read_aoi_lst(Args& a, T& E1, T& E2, bool& has_E2) {
    const std::string kind = a.text("lst", "");
    if (kind.empty())
        throw InputError(
            "'lst' names the service law whose Laplace-Stieltjes transform this age formula "
            "needs: exp, erlang, det or ph. A transform itself cannot cross a JSON boundary");
    has_E2 = true;
    if (kind == "exp") {
        const T mu = a.number<T>("lst_mu");
        E1 = T(num_traits<T>::from_int(1) / mu);
        E2 = T(num_traits<T>::from_int(2) / (mu * mu));
        return aoi::aoi_lst_exp(mu);
    }
    if (kind == "erlang") {
        const unsigned k = a.uinteger("lst_k");
        const T mu = a.number<T>("lst_mu");
        const T kk = num_traits<T>::from_int(static_cast<int>(k));
        E1 = T(kk / mu);
        E2 = T(kk * (kk + num_traits<T>::from_int(1)) / (mu * mu));
        return aoi::aoi_lst_erlang(k, mu);
    }
    if (kind == "det") {
        const T d = a.number<T>("lst_d");
        E1 = d;
        E2 = T(d * d);
        return aoi::aoi_lst_det(d);
    }
    if (kind == "ph") {
        const std::vector<T> alpha = a.vector<T>("lst_alpha");
        const Matrix<T> Tmat = a.matrix<T>("lst_T");
        // The two moments of a PH law are -alpha inv(T) 1 and 2 alpha inv(T)^2 1;
        // the caller states them rather than having this boundary invert T, so
        // the number the formula uses is the caller's own and not a second
        // inversion that could disagree with the transform beside it.
        E1 = a.number<T>("E_1");
        E2 = a.number<T>("E_2");
        return aoi::aoi_lst_ph(alpha, Tmat);
    }
    throw InputError("'lst' must be one of exp, erlang, det, ph; got '" + kind + "'");
}

/** [meanAoI, varAoI, peakAoI], the closed-form arm of the family. */
template <class T>
Json aoi_results(const aoi::AoiResult<T>& r) {
    Json out;
    out["meanAoI"] = encode_scalar(r.meanAoI);
    out["varAoI"] = encode_scalar(r.varAoI);
    out["peakAoI"] = encode_scalar(r.peakAoI);
    return out;
}

/** The transform arm: the two means, and whether an LST came back with them. */
template <class T>
Json aoi_lst_results(const aoi::AoiLstResult<T>& r) {
    Json out;
    out["meanAoI"] = encode_scalar(r.meanAoI);
    out["peakAoI"] = encode_scalar(r.peakAoI);
    out["has_lst"] = r.has_lst;
    return out;
}

#define LINE_AOI_TWO(OpName, fn, k1, k2)                                 \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const T x = a.number<T>(k1);                                 \
            const T y = a.number<T>(k2);                                 \
            a.done();                                                    \
            return aoi_results(aoi::fn(x, y));                           \
        }                                                                \
    }

LINE_AOI_TWO(OpAoiFcfsDm1, aoi_fcfs_dm1, "tau", "mu");
LINE_AOI_TWO(OpAoiFcfsMd1, aoi_fcfs_md1, "lambda", "d");
LINE_AOI_TWO(OpAoiFcfsMm1, aoi_fcfs_mm1, "lambda", "mu");
LINE_AOI_TWO(OpAoiLcfsprDm1, aoi_lcfspr_dm1, "tau", "mu");
LINE_AOI_TWO(OpAoiLcfsprMd1, aoi_lcfspr_md1, "lambda", "d");
LINE_AOI_TWO(OpAoiLcfsprMm1, aoi_lcfspr_mm1, "lambda", "mu");
#undef LINE_AOI_TWO

/** The two M/GI/1 laws that need only the first two service moments. */
#define LINE_AOI_MGI1_MOMENTS(OpName, fn)                                \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const T lambda = a.number<T>("lambda");                      \
            const T E_H = a.number<T>("E_H");                            \
            const T E_H2 = a.number<T>("E_H2");                          \
            a.done();                                                    \
            return aoi_lst_results(aoi::fn(lambda, E_H, E_H2));          \
        }                                                                \
    }

LINE_AOI_MGI1_MOMENTS(OpAoiLcfsdMgi1, aoi_lcfsd_mgi1);
LINE_AOI_MGI1_MOMENTS(OpAoiLcfssMgi1, aoi_lcfss_mgi1);
#undef LINE_AOI_MGI1_MOMENTS

struct OpAoiFcfsMgi1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        T E1 = num_traits<T>::from_int(0), E2 = num_traits<T>::from_int(0);
        bool has_E2 = false;
        const aoi::Lst<T> H = read_aoi_lst<T>(a, E1, E2, has_E2);
        a.done();
        return aoi_lst_results(aoi::aoi_fcfs_mgi1(lambda, H, E1, E2));
    }
};

struct OpAoiLcfsprMgi1 {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        T E1 = num_traits<T>::from_int(0), E2 = num_traits<T>::from_int(0);
        bool has_E2 = false;
        const aoi::Lst<T> H = read_aoi_lst<T>(a, E1, E2, has_E2);
        a.done();
        return aoi_lst_results(aoi::aoi_lcfspr_mgi1(lambda, H, E1));
    }
};

struct OpAoiFcfsGim1 {
    template <class T>
    static Json run(Args& a) {
        const T mu = a.number<T>("mu");
        T E1 = num_traits<T>::from_int(0), E2 = num_traits<T>::from_int(0);
        bool has_E2 = false;
        const aoi::Lst<T> Y = read_aoi_lst<T>(a, E1, E2, has_E2);
        a.done();
        return aoi_lst_results(aoi::aoi_fcfs_gim1(Y, mu, E1, E2));
    }
};

#define LINE_AOI_GIM1_ONE(OpName, fn)                                    \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const T mu = a.number<T>("mu");                              \
            T E1 = num_traits<T>::from_int(0), E2 = num_traits<T>::from_int(0);  \
            bool has_E2 = false;                                         \
            const aoi::Lst<T> Y = read_aoi_lst<T>(a, E1, E2, has_E2);    \
            a.done();                                                    \
            return aoi_lst_results(aoi::fn(Y, mu, E1));                  \
        }                                                                \
    }

LINE_AOI_GIM1_ONE(OpAoiLcfsdGim1, aoi_lcfsd_gim1);
LINE_AOI_GIM1_ONE(OpAoiLcfsprGim1, aoi_lcfspr_gim1);
LINE_AOI_GIM1_ONE(OpAoiLcfssGim1, aoi_lcfss_gim1);
#undef LINE_AOI_GIM1_ONE

// ---------------------------------------------------------------------------
// fj: the fork-join order-statistic bounds and approximations
// ---------------------------------------------------------------------------
//
// `fj_order_stat` is NOT exposed: its last argument is the CDF of the branch
// service law as a callable, and unlike the age family there is no set of named
// laws the reference offers in its place, so a substitute would be this file's
// choice rather than the caller's.

/** The many members returning one bare number, under the reference's name. */
#define LINE_FJ_K_LAMBDA_MU(OpName, fn, out_key)                         \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const unsigned K = a.uinteger("K");                          \
            const T lambda = a.number<T>("lambda");                      \
            const T mu = a.number<T>("mu");                              \
            a.done();                                                    \
            Json out;                                                    \
            out[out_key] = encode_scalar(fj::fn(K, lambda, mu));         \
            return out;                                                  \
        }                                                                \
    }

LINE_FJ_K_LAMBDA_MU(OpFjResptNt, fj_respt_nt, "R");
LINE_FJ_K_LAMBDA_MU(OpFjResptVarki, fj_respt_varki, "R");
LINE_FJ_K_LAMBDA_MU(OpFjResptVm, fj_respt_vm, "R");
LINE_FJ_K_LAMBDA_MU(OpFjRmax, fj_rmax, "Rmax");
#undef LINE_FJ_K_LAMBDA_MU

#define LINE_FJ_K_MU(OpName, fn, out_key)                                \
    struct OpName {                                                      \
        template <class T>                                               \
        static Json run(Args& a) {                                       \
            const unsigned K = a.uinteger("K");                          \
            const T mu = a.number<T>("mu");                              \
            a.done();                                                    \
            Json out;                                                    \
            out[out_key] = encode_scalar(fj::fn(K, mu));                 \
            return out;                                                  \
        }                                                                \
    }

LINE_FJ_K_MU(OpFjSmTput, fj_sm_tput, "X");
LINE_FJ_K_MU(OpFjXmaxEmma, fj_xmax_emma, "Xmax");
LINE_FJ_K_MU(OpFjXmaxExp, fj_xmax_exp, "Xmax");
#undef LINE_FJ_K_MU

struct OpFjHarmonic {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        a.done();
        Json out;
        out["HK"] = encode_scalar(fj::fj_harmonic<T>(K));
        return out;
    }
};

struct OpFjQuantile {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T q = a.number<T>("q");
        a.done();
        Json out;
        out["x"] = encode_scalar(fj::fj_quantile(K, q));
        return out;
    }
};

struct OpFjResptTwoway {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        a.done();
        Json out;
        out["R"] = encode_scalar(fj::fj_respt_2way(lambda, mu));
        return out;
    }
};

struct OpFjSynchDelay {
    template <class T>
    static Json run(Args& a) {
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        a.done();
        Json out;
        out["D"] = encode_scalar(fj::fj_synch_delay(lambda, mu));
        return out;
    }
};

struct OpFjXmax2 {
    template <class T>
    static Json run(Args& a) {
        const T lambda1 = a.number<T>("lambda1");
        const T lambda2 = a.number<T>("lambda2");
        a.done();
        Json out;
        out["Xmax"] = encode_scalar(fj::fj_xmax_2(lambda1, lambda2));
        return out;
    }
};

struct OpFjXmaxErlang {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const unsigned k = a.uinteger("k");
        const T mu = a.number<T>("mu");
        a.done();
        Json out;
        out["Xmax"] = encode_scalar(fj::fj_xmax_erlang(K, k, mu));
        return out;
    }
};

struct OpFjRmaxErlang {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const unsigned k = a.uinteger("k");
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        a.done();
        Json out;
        out["Rmax"] = encode_scalar(fj::fj_rmax_erlang(K, k, lambda, mu));
        return out;
    }
};

struct OpFjXmaxHyperexp {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T p1 = a.number<T>("p1");
        const T mu1 = a.number<T>("mu1");
        const T mu2 = a.number<T>("mu2");
        a.done();
        Json out;
        out["Xmax"] = encode_scalar(fj::fj_xmax_hyperexp(K, p1, mu1, mu2));
        return out;
    }
};

/**
 * The EVD response-time fit. `calibrated` selects the version whose constants
 * were refitted, which is a DIFFERENT approximation and not a refinement of the
 * same one, so it defaults to the reference's own false rather than to the
 * one that happens to be more accurate on any given model.
 */
struct OpFjRmaxEvd {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T R = a.number<T>("R");
        const T sigma_R = a.number<T>("sigma_R");
        const bool calibrated = a.boolean("calibrated", false);
        a.done();
        Json out;
        out["Rmax"] = encode_scalar(fj::fj_rmax_evd(K, R, sigma_R, calibrated));
        return out;
    }
};

struct OpFjBounds {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T lambda = a.number<T>("lambda");
        const T mu = a.number<T>("mu");
        a.done();
        const fj::FJBoundsResult<T> r = fj::fj_bounds(K, lambda, mu);
        Json out;
        out["Rmax"] = encode_scalar(r.Rmax);
        out["Rmin"] = encode_scalar(r.Rmin);
        return out;
    }
};

struct OpFjCharMax {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T mu = a.number<T>("mu");
        a.done();
        const fj::FJCharMaxResult<T> r = fj::fj_char_max(K, mu);
        Json out;
        out["MK"] = encode_scalar(r.MK);
        out["mK"] = encode_scalar(r.mK);
        return out;
    }
};

struct OpFjGkBound {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        a.done();
        const fj::FJGKBoundResult<T> r = fj::fj_gk_bound<T>(K);
        Json out;
        out["K"] = encode_count(r.K);
        out["exponential"] = encode_scalar(r.exponential);
        out["uniform"] = encode_scalar(r.uniform);
        out["evd"] = encode_scalar(r.evd);
        out["upper_bound"] = encode_scalar(r.upper_bound);
        return out;
    }
};

struct OpFjOrdstatExp {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> ri = a.vector<T>("ri");
        const std::size_t k = a.count("k");
        a.done();
        Json out;
        out["m"] = encode_scalar(fj::fj_ordstat_exp<T>(ri, k));
        return out;
    }
};

struct OpFjQuorumMoments {
    template <class T>
    static Json run(Args& a) {
        const std::vector<T> branchMeans = a.vector<T>("branchMeans");
        const std::vector<T> branchVars = a.vector<T>("branchVars");
        const std::size_t k = a.count("k");
        a.done();
        const fj::FJQuorumMomentsResult<T> r =
            fj::fj_quorum_moments(branchMeans, branchVars, k);
        Json out;
        out["m"] = encode_scalar(r.m);
        out["v"] = encode_scalar(r.v);
        return out;
    }
};

struct OpFjXmaxApprox {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T mu_X = a.number<T>("mu_X");
        const T sigma_X = a.number<T>("sigma_X");
        const std::string t = a.text("type", "Exp");
        fj::FJDistType type = fj::FJDistType::Exp;
        if (t == "Exp") type = fj::FJDistType::Exp;
        else if (t == "Uniform") type = fj::FJDistType::Uniform;
        else if (t == "Evd") type = fj::FJDistType::Evd;
        else if (t == "Bound") type = fj::FJDistType::Bound;
        else throw InputError("'type' must be Exp, Uniform, Evd or Bound; got '" + t + "'");
        a.done();
        const fj::FJXmaxApproxResult<T> r = fj::fj_xmax_approx(K, mu_X, sigma_X, type);
        Json out;
        out["Xmax"] = encode_scalar(r.Xmax);
        out["GK"] = encode_scalar(r.GK);
        return out;
    }
};

struct OpFjXmaxNormal {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T mu = a.number<T>("mu");
        const T sigma = a.number<T>("sigma");
        const std::string m = a.text("method", "Johnson");
        fj::FJNormalMethod method = fj::FJNormalMethod::Johnson;
        if (m == "Johnson") method = fj::FJNormalMethod::Johnson;
        else if (m == "Arnold") method = fj::FJNormalMethod::Arnold;
        else if (m == "Corrected") method = fj::FJNormalMethod::Corrected;
        else throw InputError("'method' must be Johnson, Arnold or Corrected; got '" + m + "'");
        a.done();
        const fj::FJXmaxNormalResult<T> r = fj::fj_xmax_normal(K, mu, sigma, method);
        Json out;
        out["Xmax"] = encode_scalar(r.Xmax);
        out["Vmax"] = encode_scalar(r.Vmax);
        return out;
    }
};

struct OpFjXmaxPareto {
    template <class T>
    static Json run(Args& a) {
        const unsigned K = a.uinteger("K");
        const T beta = a.number<T>("beta");
        const T k = a.number<T>("k");
        a.done();
        const fj::FJXmaxParetoResult<T> r = fj::fj_xmax_pareto(K, beta, k);
        Json out;
        out["Xmax"] = encode_scalar(r.Xmax);
        out["MK"] = encode_scalar(r.MK);
        return out;
    }
};



using Invoker = std::function<Json(const ArithSpec&, Args&)>;

template <class Op>
Invoker make() {
    return [](const ArithSpec& s, Args& a) { return run_at<Op>(s, a); };
}

const std::map<std::string, Invoker>& dispatch_table() {
    static const std::map<std::string, Invoker> table = {
        {"pfqn_ca", make<OpPfqnCa>()},
        {"pfqn_conv", make<OpPfqnConv>()},
        {"pfqn_recal", make<OpPfqnRecal>()},
        {"pfqn_gld", make<OpPfqnGld>()},
        {"pfqn_mva", make<OpPfqnMva>()},
        {"pfqn_comom", make<OpPfqnComom>()},
        {"pfqn_comomrm", make<OpPfqnComomrm>()},
        {"pfqn_comomrm_orig", make<OpPfqnComomrmOrig>()},
        {"pfqn_comomrm_ms", make<OpPfqnComomrmMs>()},
        {"pfqn_procomom", make<OpPfqnProcomom>()},
        {"ctmc_solve", make<OpCtmcSolve>()},
        {"dtmc_solve", make<OpDtmcSolve>()},
        {"ctmc_makeinfgen", make<OpCtmcMakeinfgen>()},
        {"sim_vonneumann", make<OpSimVonneumann>()},
        {"sim_shapirowilk", make<OpSimShapirowilk>()},
        {"sim_sts_quantile_areas", make<OpSimStsQuantileAreas>()},
        {"sim_fquest", make<OpSimFquest>()},
        {"sim_firquest", make<OpSimFirquest>()},
        {"qsys_bmapphnn_retrial", make<OpQsysBmapphnnRetrial>()},
        {"aoi_fcfs_dm1", make<NeedsTranscendental<OpAoiFcfsDm1>>()},
        {"aoi_fcfs_gim1", make<NeedsTranscendental<OpAoiFcfsGim1>>()},
        {"aoi_fcfs_md1", make<NeedsTranscendental<OpAoiFcfsMd1>>()},
        {"aoi_fcfs_mgi1", make<NeedsTranscendental<OpAoiFcfsMgi1>>()},
        {"aoi_fcfs_mm1", make<OpAoiFcfsMm1>()},
        {"aoi_lcfsd_gim1", make<NeedsTranscendental<OpAoiLcfsdGim1>>()},
        {"aoi_lcfsd_mgi1", make<OpAoiLcfsdMgi1>()},
        {"aoi_lcfspr_dm1", make<NeedsTranscendental<OpAoiLcfsprDm1>>()},
        {"aoi_lcfspr_gim1", make<NeedsTranscendental<OpAoiLcfsprGim1>>()},
        {"aoi_lcfspr_md1", make<NeedsTranscendental<OpAoiLcfsprMd1>>()},
        {"aoi_lcfspr_mgi1", make<NeedsTranscendental<OpAoiLcfsprMgi1>>()},
        {"aoi_lcfspr_mm1", make<OpAoiLcfsprMm1>()},
        {"aoi_lcfss_gim1", make<NeedsTranscendental<OpAoiLcfssGim1>>()},
        {"aoi_lcfss_mgi1", make<OpAoiLcfssMgi1>()},
        {"fj_bounds", make<OpFjBounds>()},
        {"fj_char_max", make<NeedsTranscendental<OpFjCharMax>>()},
        {"fj_gk_bound", make<NeedsTranscendental<OpFjGkBound>>()},
        {"fj_harmonic", make<OpFjHarmonic>()},
        {"fj_quantile", make<NeedsTranscendental<OpFjQuantile>>()},
        {"fj_ordstat_exp", make<OpFjOrdstatExp>()},
        {"fj_quorum_moments", make<NeedsTranscendental<OpFjQuorumMoments>>()},
        {"fj_respt_2way", make<OpFjResptTwoway>()},
        {"fj_respt_nt", make<OpFjResptNt>()},
        {"fj_respt_varki", make<OpFjResptVarki>()},
        {"fj_respt_vm", make<OpFjResptVm>()},
        {"fj_rmax", make<OpFjRmax>()},
        {"fj_rmax_erlang", make<NeedsTranscendental<OpFjRmaxErlang>>()},
        {"fj_rmax_evd", make<NeedsTranscendental<OpFjRmaxEvd>>()},
        {"fj_sm_tput", make<OpFjSmTput>()},
        {"fj_synch_delay", make<OpFjSynchDelay>()},
        {"fj_xmax_2", make<OpFjXmax2>()},
        {"fj_xmax_approx", make<NeedsTranscendental<OpFjXmaxApprox>>()},
        {"fj_xmax_emma", make<NeedsTranscendental<OpFjXmaxEmma>>()},
        {"fj_xmax_erlang", make<NeedsTranscendental<OpFjXmaxErlang>>()},
        {"fj_xmax_exp", make<OpFjXmaxExp>()},
        {"fj_xmax_hyperexp", make<OpFjXmaxHyperexp>()},
        {"fj_xmax_normal", make<NeedsTranscendental<OpFjXmaxNormal>>()},
        {"fj_xmax_pareto", make<NeedsTranscendental<OpFjXmaxPareto>>()},
        {"ctmc_courtois", make<NeedsTranscendental<OpCtmcCourtois>>()},
        {"ctmc_fau", make<NeedsTranscendental<OpCtmcFau>>()},
        {"ctmc_foxglynn", make<NeedsTranscendental<OpCtmcFoxglynn>>()},
        {"ctmc_saddlepoint", make<NeedsTranscendental<OpCtmcSaddlepoint>>()},
        {"ctmc_bicgstab", make<NeedsTranscendental<OpCtmcBicgstab>>()},
        {"ctmc_bicgstab_multi", make<NeedsTranscendental<OpCtmcBicgstabMulti>>()},
        {"ctmc_gmres", make<NeedsTranscendental<OpCtmcGmres>>()},
        {"ctmc_gmres_multi", make<NeedsTranscendental<OpCtmcGmresMulti>>()},
        {"ctmc_kms", make<NeedsTranscendental<OpCtmcKms>>()},
        {"ctmc_multi", make<NeedsTranscendental<OpCtmcMulti>>()},
        {"ctmc_randomization", make<OpCtmcRandomization>()},
        {"ctmc_sens", make<OpCtmcSens>()},
        {"ctmc_solve_reducible", make<OpCtmcSolveReducible>()},
        {"ctmc_solve_reducible_blkdecomp", make<OpCtmcSolveReducibleBlkdecomp>()},
        {"ctmc_takahashi", make<NeedsTranscendental<OpCtmcTakahashi>>()},
        {"ctmc_transient", make<NeedsTranscendental<OpCtmcTransient>>()},
        {"ctmc_transient_sens", make<NeedsTranscendental<OpCtmcTransientSens>>()},
        {"dtmc_makestochastic", make<OpDtmcMakestochastic>()},
        {"dtmc_solve_reducible", make<OpDtmcSolveReducible>()},
        {"moment_binomial_from_factorial", make<OpMomentBinomialFromFactorial>()},
        {"moment_binomial_from_negbinomial", make<OpMomentBinomialFromNegbinomial>()},
        {"moment_binotrans", make<OpMomentBinotrans>()},
        {"moment_binotransinv", make<OpMomentBinotransinv>()},
        {"moment_central_from_raw", make<OpMomentCentralFromRaw>()},
        {"moment_factorial_from_binomial", make<OpMomentFactorialFromBinomial>()},
        {"moment_factorial_from_raw", make<OpMomentFactorialFromRaw>()},
        {"moment_factorial_from_upfactorial", make<OpMomentFactorialFromUpfactorial>()},
        {"moment_lah", make<OpMomentLah>()},
        {"moment_negbinomial_from_binomial", make<OpMomentNegbinomialFromBinomial>()},
        {"moment_negbinomial_from_upfactorial", make<OpMomentNegbinomialFromUpfactorial>()},
        {"moment_raw_from_central", make<OpMomentRawFromCentral>()},
        {"moment_raw_from_factorial", make<OpMomentRawFromFactorial>()},
        {"moment_raw_from_upfactorial", make<OpMomentRawFromUpfactorial>()},
        {"moment_stirling1", make<OpMomentStirling1>()},
        {"moment_stirling2", make<OpMomentStirling2>()},
        {"moment_stirlingcycle", make<OpMomentStirlingcycle>()},
        {"moment_upfactorial_from_factorial", make<OpMomentUpfactorialFromFactorial>()},
        {"moment_upfactorial_from_negbinomial", make<OpMomentUpfactorialFromNegbinomial>()},
        {"moment_upfactorial_from_raw", make<OpMomentUpfactorialFromRaw>()},
        {"stronglyconncomp", make<OpStronglyConnComp>()},
        {"qsys_dmc", make<NeedsTranscendental<OpQsysDmc>>()},
        {"dqsys_geogeo1", make<OpQsysGeoGeo1>()},
        {"dqsys_geoxgeo1", make<OpQsysGeoxGeo1>()},
        {"qsys_gg1", make<NeedsTranscendental<OpQsysGg1>>()},
        {"qsys_gig1_approx_allencunneen", make<OpQsysGig1AllenCunneen>()},
        {"qsys_gig1_approx_gelenbe", make<NeedsTranscendental<OpQsysGig1Gelenbe>>()},
        {"qsys_gig1_approx_heyman", make<OpQsysGig1Heyman>()},
        {"qsys_gig1_approx_kimura", make<OpQsysGig1Kimura>()},
        {"qsys_gig1_approx_klb", make<NeedsTranscendental<OpQsysGig1Klb>>()},
        {"qsys_gig1_approx_kobayashi", make<NeedsTranscendental<OpQsysGig1Kobayashi>>()},
        {"qsys_gig1_approx_marchal", make<OpQsysGig1Marchal>()},
        {"qsys_gig1_approx_myskja", make<NeedsTranscendental<OpQsysGig1Myskja>>()},
        {"qsys_gig1_approx_myskja2", make<NeedsTranscendental<OpQsysGig1Myskja2>>()},
        {"qsys_gig1_lbnd", make<OpQsysGig1Lbnd>()},
        {"qsys_gig1_ubnd_kingman", make<OpQsysGig1UbndKingman>()},
        {"qsys_gigk_approx", make<NeedsTranscendental<OpQsysGigkApprox>>()},
        {"qsys_gigk_approx_cosmetatos", make<NeedsTranscendental<OpQsysGigkCosmetatos>>()},
        {"qsys_gigk_approx_kingman", make<OpQsysGigkKingman>()},
        {"qsys_gigk_approx_whitt", make<NeedsTranscendental<OpQsysGigkWhitt>>()},
        {"qsys_erlanga", make<NeedsTranscendental<OpQsysErlangA>>()},
        {"qsys_ggnm_diffusion", make<NeedsTranscendental<OpQsysGgnmDiffusion>>()},
        {"qsys_gig1_bnds_extremal", make<NeedsTranscendental<OpQsysGig1BndsExtremal>>()},
        {"qsys_mmk_qed", make<NeedsTranscendental<OpQsysMmkQed>>()},
        {"qsys_mmk_qed_alpha", make<NeedsTranscendental<OpQsysMmkQedAlpha>>()},
        {"qsys_mmk_qed_staffing", make<NeedsTranscendental<OpQsysMmkQedStaffing>>()},
        {"qsys_gigk_rqt", make<NeedsTranscendental<OpQsysGigkRqt>>()},
        {"qsys_gig1_rqt", make<NeedsTranscendental<OpQsysGig1Rqt>>()},
        {"qsys_gigk_rqt_gamma", make<NeedsTranscendental<OpQsysGigkRqtGamma>>()},
        {"qsys_gm1", make<OpQsysGm1>()},
        {"qsys_mapd1", make<NeedsTranscendental<OpQsysMapd1>>()},
        {"qsys_mapdc", make<NeedsTranscendental<OpQsysMapdc>>()},
        {"qsys_mapg1", make<NeedsTranscendental<OpQsysMapg1>>()},
        {"qsys_mapm1", make<NeedsTranscendental<OpQsysMapm1>>()},
        {"qsys_mapmap1", make<NeedsTranscendental<OpQsysMapmap1>>()},
        {"qsys_mapmc", make<NeedsTranscendental<OpQsysMapmc>>()},
        {"qsys_mapph1", make<NeedsTranscendental<OpQsysMapph1>>()},
        {"qsys_mapphc", make<NeedsTranscendental<OpQsysMapphc>>()},
        {"qsys_mg1", make<OpQsysMg1>()},
        {"qsys_mg1_fb", make<NeedsTranscendental<OpQsysMg1Fb>>()},
        {"qsys_mg1_lrpt", make<NeedsTranscendental<OpQsysMg1Lrpt>>()},
        {"qsys_mg1_prio", make<OpQsysMg1Prio>()},
        {"qsys_mg1_psjf", make<NeedsTranscendental<OpQsysMg1Psjf>>()},
        {"qsys_mg1_setf", make<NeedsTranscendental<OpQsysMg1Setf>>()},
        {"qsys_mg1_srpt", make<NeedsTranscendental<OpQsysMg1Srpt>>()},
        {"qsys_mg1k_loss_mgs", make<NeedsTranscendental<OpQsysMg1kLossMgs>>()},
        {"qsys_mginf", make<NeedsTranscendental<OpQsysMginf>>()},
        {"qsys_mm1", make<OpQsysMm1>()},
        {"qsys_mm1_dps", make<NeedsTranscendental<OpQsysMm1Dps>>()},
        {"qsys_mm1k_loss", make<OpQsysMm1kLoss>()},
        {"qsys_mmcc_retrial_fp", make<NeedsTranscendental<OpQsysMmccRetrialFp>>()},
        {"qsys_mmck", make<OpQsysMmck>()},
        {"qsys_mmk", make<OpQsysMmk>()},
        {"qsys_mxm1", make<OpQsysMxm1>()},
        {"qsys_phm1", make<NeedsTranscendental<OpQsysPhm1>>()},
        {"qsys_phmc", make<NeedsTranscendental<OpQsysPhmc>>()},
        {"qsys_phph1", make<NeedsTranscendental<OpQsysPhph1>>()},
    };
    return table;
}

std::string join(const std::vector<std::string>& v, const char* sep) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) s += sep;
        s += v[i];
    }
    return s;
}

}  // namespace

std::vector<std::string> api_exposed_functions() {
    std::vector<std::string> names;
    for (const auto& kv : dispatch_table()) names.push_back(kv.first);
    return names;  // std::map already orders them
}

bool api_is_exposed(const std::string& name) {
    return dispatch_table().find(name) != dispatch_table().end();
}

Json api_invoke(const std::string& name, const std::string& arith, const Json& args) {
    const ApiEntry* entry = find_api(name);
    if (entry == nullptr)
        throw UnsupportedError("API function '" + name +
                               "' is not ported to C++ yet (--list-api shows what is).");

    const ArithSpec spec = parse_arith(arith);

    // arithmetic-gate-before-exposure-gate rationale: see _kb/14-cpp-multiprecision.md
    if (!api_supports(*entry, spec.mode)) {
        std::vector<std::string> modes;
        for (Arith m : entry->arith) modes.push_back(arith_name(m));
        throw UnsupportedError("API function '" + name + "' does not support --arith " +
                               spec.str() + "; it supports: " + join(modes, ", ") + ".");
    }

    auto it = dispatch_table().find(name);
    if (it == dispatch_table().end()) {
        // The two api domains that take a MODEL, not matrices, can never be
        // reached from here whatever the port does next, and telling a caller
        // to wait for an exposure that will never come is worse than telling it
        // nothing: it names a route that exists instead of the one that does.
        if (entry->domain == "sn" || entry->domain == "lqn")
            throw UnsupportedError(
                "API function '" + name +
                "' takes a NetworkStruct or a LayeredNetworkStruct, not matrices, so it has no "
                "named-JSON-argument form and --api cannot carry it in any future version. Reach "
                "the model-level answer with -s <solver> -a <analysis>, or call it from the "
                "library.");
        throw UnsupportedError(
            "API function '" + name +
            "' is ported to C++ but not yet exposed over --api; it is reachable from the library "
            "and its tests only. Exposed over --api so far: " +
            join(api_exposed_functions(), ", ") + ".");
    }

    Args reader(args, name);
    Json out;
    out["function"] = name;
    out["arith"] = spec.str();
    out["results"] = it->second(spec, reader);
    return out;
}

// ---------------------------------------------------------------------------
// Readable rendering
// ---------------------------------------------------------------------------

namespace {

void render_value(std::ostringstream& os, const Json& v, const std::string& indent) {
    // readable-vs-JSON rendering consistency: see _kb/14-cpp-multiprecision.md
    if (v.is_object() && v.contains("double")) {
        os << v["double"].dump();
        if (v.contains("num"))
            os << "  = " << v["num"].get<std::string>() << " / " << v["den"].get<std::string>();
        else if (v.contains("dec"))
            os << "  = " << v["dec"].get<std::string>();
        os << "\n";
        return;
    }
    if (v.is_array()) {
        if (!v.empty() && v[0].is_array()) {
            os << "\n";
            for (const Json& row : v) {
                os << indent << "  ";
                for (std::size_t k = 0; k < row.size(); ++k) {
                    if (k) os << "  ";
                    if (row[k].is_object() && row[k].contains("double"))
                        os << row[k]["double"].dump();
                    else
                        os << row[k].dump();
                }
                os << "\n";
            }
            return;
        }
        os << "[";
        for (std::size_t k = 0; k < v.size(); ++k) {
            if (k) os << ", ";
            if (v[k].is_object() && v[k].contains("double"))
                os << v[k]["double"].dump();
            else
                os << v[k].dump();
        }
        os << "]\n";
        return;
    }
    os << v.dump() << "\n";
}

}  // namespace

std::string api_render_readable(const Json& result) {
    std::ostringstream os;
    os << result["function"].get<std::string>() << " (arith: "
       << result["arith"].get<std::string>() << ")\n";
    const Json& res = result["results"];
    for (auto it = res.begin(); it != res.end(); ++it) {
        os << "  " << it.key() << " = ";
        render_value(os, it.value(), "  ");
    }
    return os.str();
}

}  // namespace reg
}  // namespace line
