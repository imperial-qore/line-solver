/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_NONMARKOV_TOPH_H
#define LINE_API_SN_SN_NONMARKOV_TOPH_H

/**
 * Replace every non-Markovian service and firing law by a Markovian surrogate.
 *
 * Templated port of matlab/src/api/sn/sn_nonmarkov_toph.m, mirrored by
 * jline.api.sn.SnNonmarkovToPh and the native Python sn_nonmarkov_toph. The
 * MAM, CTMC, Fluid and SSA analyzers all run it on their own copy of the struct
 * before any state space is built, so a Gamma, Weibull, Lognormal, Pareto,
 * Uniform or Det service law reaches an algorithm that only understands
 * generators.
 *
 * IT IS NOT `convertToMAP`. The struct refresh already lowers those families to
 * an Erlang of ceil(1/scv) phases (`dist_to_map` in lang/distribution.h); that
 * fit matches the mean and, above SCV 1, nothing else. This is the SOLVER-side
 * conversion, with an explicit phase budget (20 by default) and a choice of fit:
 *
 *  - phfit = Cme, the default: a concentrated matrix exponential convolved with
 *    an exponential (`dist_fit_me`), which matches the first TWO moments exactly
 *    whenever the SCV is in (0,1). It is what makes an M/G/1 mean come out
 *    right: on M/Gamma/1 at rho 0.5 the two-moment fit lands on the exact mean
 *    queue length where the density fit lands 2.7e-2 away.
 *  - phfit = Ph: the Bernstein density fit (`map_bernstein`), kept because it
 *    carries SHAPE information a two-moment fit cannot, and because SSA, Fluid
 *    and JMT cannot consume a matrix exponential at all.
 *
 * DET IS ERLANG, WHATEVER phfit ASKS FOR. A concentrated ME matches a Det's
 * moments far better (SCV 5.7e-3 against Erlang-20's 0.05) but is not a
 * generator: its off-diagonal entries are not rates, so a CTMC built from it
 * does not describe the model. Measured on the reference's test_OQN_DM1 (Det(1)
 * arrivals, Exp(2) service, rho = 0.5), the ME surrogate reported U = 0.993 and
 * a departure rate of 2.0 against a mean interarrival of 1.0, where JMT gives
 * 0.497, LDES 0.502 and the golden 0.500.
 *
 * WHAT IS NOT TOUCHED. The Markovian families, DISABLED, IMMEDIATE, and the
 * three time-inhomogeneous families NHPP / MAPt / PHt, whose whole content is
 * the schedule -- collapsing one to a single homogeneous surrogate would erase
 * the time dependence that is the reason it was declared.
 *
 * NO STATE SURGERY IS NEEDED HERE. The MATLAB reference rewrites sn.phases,
 * phasessz, phaseshift, mu, phi, pie and then splices extra columns into any
 * pre-initialized sn.state and sn.space, because those are stored arrays that
 * would otherwise disagree with the new phase count. In this port every one of
 * them is DERIVED from `sn.service[i][r]` on demand (`phases_of`, `dist_pie`,
 * `dist_to_map`), so replacing the distribution updates all of them at once and
 * there is nothing left to splice. Callers run this before generating states,
 * exactly as the reference's analyzers do.
 *
 * ARITHMETIC: transcendental. The fits evaluate densities and take square
 * roots, so this does not instantiate under Rational.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <vector>

#include "line/api/mam/cme.h"
#include "line/api/mam/hyperexp_fit_longtail.h"
#include "line/api/mam/map_bernstein.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/sn/sn_is_phasetype.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace api {

/** Which Markovian surrogate to fit; `options.config.phfit`. */
enum class PhFit {
    Cme,  ///< concentrated ME plus exponential tail: two moments, exactly
    Ph,   ///< Bernstein density fit: a genuine phase-type, shape-carrying
    /**
     * Mixture of exponentials fitted to the ccdf ITSELF across decades of time
     * scale (`hyperexp_fit_longtail`, Feldmann and Whitt 1998).
     *
     * The only one of the three that says anything about a LONG TAIL: a Pareto
     * with tail index below 2 has no finite variance, so a two-moment fit does
     * not exist at all, and even where the moments are finite they say nothing
     * about the several orders of magnitude over which such a law acts. It
     * applies only to the families with a closed-form tail below; a light-tailed
     * law falls through to the two-moment surrogate.
     */
    Hyperexp
};

/** `options.config.nonmkv` and friends. */
struct NonmarkovOptions {
    bool enabled = true;       ///< false is the reference's nonmkv = 'none'
    std::size_t order = 20;    ///< `nonmkvorder`, the phase budget
    PhFit phfit = PhFit::Cme;  ///< which surrogate family
    /**
     * Leave Det alone for the exact MAP/D/c branch, `options.config.preserveDet`.
     * The MAM analyzer sets it; every other caller converts.
     */
    bool preserve_det = false;
};

namespace detail {

/** True for the families that already carry a Markovian representation. */
inline bool nonmkv_is_markovian(lang::ProcessType p) {
    using lang::ProcessType;
    switch (p) {
        case ProcessType::EXP:
        case ProcessType::ERLANG:
        case ProcessType::HYPEREXP:
        case ProcessType::PH:
        case ProcessType::APH:
        case ProcessType::MAP:
        case ProcessType::DMAP:
        case ProcessType::MMAP:
        case ProcessType::BMAP:
        case ProcessType::ME:
        case ProcessType::RAP:
        case ProcessType::COXIAN:
        case ProcessType::COX2:
        case ProcessType::MMPP2:
        case ProcessType::IMMEDIATE:
        case ProcessType::DISABLED:
            return true;
        default:
            return false;
    }
}

/**
 * True for the families whose content IS a schedule.
 *
 * A single homogeneous surrogate has no way to carry a piecewise-constant
 * D0(t)/D1(t), so converting one would silently answer for a stationary model.
 */
inline bool nonmkv_is_scheduled(lang::ProcessType p) {
    using lang::ProcessType;
    return p == ProcessType::NHPP || p == ProcessType::MAPT || p == ProcessType::PHT;
}

/**
 * The density of the five families the reference fits by shape, or an empty
 * function for anything else (which then takes the moment-only surrogate).
 */
template <class T>
std::function<double(double)> nonmkv_density(const lang::Distrib<T>& d) {
    using lang::ProcessType;
    std::vector<double> p;
    for (std::size_t i = 0; i < d.params.size(); ++i)
        p.push_back(num_traits<T>::to_double(d.params[i]));
    switch (d.type) {
        case ProcessType::GAMMA: {  // params (shape, scale)
            if (p.size() < 2) break;
            const double k = p[0], th = p[1];
            return [k, th](double x) {
                if (!(x > 0.0)) return 0.0;
                return std::exp((k - 1.0) * std::log(x) - x / th - k * std::log(th) -
                                std::lgamma(k));
            };
        }
        case ProcessType::WEIBULL: {  // params (scale, shape), the C++ builder order
            if (p.size() < 2) break;
            const double a = p[0], r = p[1];
            return [a, r](double x) {
                if (!(x > 0.0)) return 0.0;
                return (r / a) * std::pow(x / a, r - 1.0) * std::exp(-std::pow(x / a, r));
            };
        }
        case ProcessType::LOGNORMAL: {  // params (mu, sigma) of the log
            if (p.size() < 2) break;
            const double mu = p[0], sg = p[1];
            return [mu, sg](double x) {
                if (!(x > 0.0)) return 0.0;
                const double z = (std::log(x) - mu) / sg;
                return std::exp(-0.5 * z * z) / (x * sg * std::sqrt(2.0 * 3.14159265358979323846));
            };
        }
        case ProcessType::PARETO: {  // params (shape alpha, scale k)
            if (p.size() < 2) break;
            const double al = p[0], k = p[1];
            return [al, k](double x) {
                if (!(x >= k) || !(x > 0.0)) return 0.0;
                return al * std::pow(k, al) / std::pow(x, al + 1.0);
            };
        }
        case ProcessType::UNIFORM: {  // params (a, b)
            if (p.size() < 2) break;
            const double a = p[0], b = p[1];
            return [a, b](double x) { return (x >= a && x <= b) ? 1.0 / (b - a) : 0.0; };
        }
        default:
            break;
    }
    return std::function<double(double)>();
}

/**
 * The declared law's complementary cdf, for the families whose tail the
 * long-tail fit is stated for. Empty for every other law, which is what makes
 * `phfit = Hyperexp` fall through to the two-moment surrogate there.
 *
 * Built on `dist_cdf`, which already carries the closed form of each of these
 * four families -- the Gamma from the regularized incomplete gamma, the Pareto
 * from 1 - (k/x)^alpha -- so there is no second copy of any of them here.
 */
template <class T>
std::function<double(double)> nonmkv_ccdf(const lang::Distrib<T>& d) {
    using lang::ProcessType;
    if (d.type != ProcessType::GAMMA && d.type != ProcessType::WEIBULL &&
        d.type != ProcessType::LOGNORMAL && d.type != ProcessType::PARETO)
        return std::function<double(double)>();
    const lang::Distrib<T> law = d;
    return [law](double x) {
        if (!(x > 0.0)) return 1.0;
        return 1.0 - num_traits<T>::to_double(lang::dist_cdf(law, num_traits<T>::from_double(x)));
    };
}

/**
 * A hyperexponential fitted to the ccdf across decades, as its MAP pair and
 * rescaled to the mean the struct carries.
 *
 * `ok` comes back false when the recursion declines the law: the components have
 * to dominate one another at their own time scales, which a light-tailed law
 * does not provide, and answering with a fit that does not hold is worse than
 * falling through to the two-moment surrogate.
 */
template <class T>
mam::Map<T> nonmkv_longtail(const std::function<double(double)>& ccdf, double mean, bool& ok) {
    ok = false;
    mam::Map<T> out;
    // The fit itself runs in double: it bisects for quantiles and takes
    // logarithms of the ccdf, which is a `double(double)` here whatever T is.
    mam::HyperexpLongtailResult<double> fit;
    try {
        fit = mam::hyperexp_fit_longtail<double>(ccdf);
    } catch (const std::exception&) {
        return out;
    }
    const std::size_t n = fit.p.size();
    if (n == 0) return out;
    Matrix<T> D0(n, n, num_traits<T>::from_int(0));
    Matrix<T> D1(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(fit.lambda[i]) || fit.lambda[i] <= 0.0) return out;
        D0(i, i) = num_traits<T>::from_double(-fit.lambda[i]);
        for (std::size_t j = 0; j < n; ++j)
            D1(i, j) = num_traits<T>::from_double(fit.lambda[i] * fit.p[j]);
    }
    mam::Map<T> m;
    m.D0 = D0;
    m.D1 = D1;
    out = mam::map_scale(m, num_traits<T>::from_double(mean));
    ok = true;
    return out;
}

/**
 * The moment-only surrogate: a concentrated ME under phfit = Cme, an Erlang of
 * the full budget otherwise. Port of the reference's fitConcentratedSurrogate.
 */
template <class T>
mam::Map<T> nonmkv_concentrated(double mean, double scv, std::size_t order, PhFit phfit) {
    if (phfit == PhFit::Cme) {
        // A Det has SCV 0, which no ME attains; the budget-limited branch of the
        // fitter then returns the most concentrated member that fits.
        const double s = scv > 1e-12 ? scv : 1e-12;
        if (s < 1.0) return mam::dist_fit_me<T>(mean, s < 1.0 - 1e-12 ? s : 1.0 - 1e-12, order);
    }
    return mam::map_erlang(num_traits<T>::from_double(mean), static_cast<unsigned>(order));
}

/** Install a fitted (D0,D1) as the law of a station-class pair or a firing mode. */
template <class T>
void nonmkv_install(lang::Distrib<T>& d, const mam::Map<T>& m) {
    std::vector<Matrix<T>> blocks;
    blocks.push_back(m.D0);
    blocks.push_back(m.D1);
    const bool isph = sn_is_phasetype(blocks, std::vector<T>());
    // The declared law is kept alongside the fit. Everything downstream reads
    // the surrogate, as it did; MMAP[K]/G[K]/1 evaluates the ORIGINAL transform
    // and would otherwise find only the fit it exists to avoid.
    std::shared_ptr<lang::Distrib<T>> declared(new lang::Distrib<T>(d));
    declared->declared.reset();
    // A concentrated-ME surrogate is NOT a phase-type, so it is tagged ME: the
    // CTMC then assembles a rational generator and the PH-only consumers refuse
    // it. sn_is_phasetype is the single test, exactly as in the reference.
    d = lang::Distrib<T>::map_dist(m.D0, m.D1,
                                   isph ? lang::ProcessType::PH : lang::ProcessType::ME);
    lang::dist_refresh_moments(d);
    d.declared = declared;
}

/**
 * True when mode `m`'s firing law is ALREADY a usable Markovian one and needs
 * no fit.
 *
 * A POSITIVE PHASE COUNT IS NOT ENOUGH, and reading it as such is what left
 * `spn_pareto_service` absorbing. `network_reader` records
 * `firingphases[m] = dist_to_map(proc).order()` at read time, so a Pareto mode
 * arrives with a phase count of the fit that WOULD be made while `firingproc[m]`
 * is still the Pareto, whose (D0,D1) are empty. The old guard read that count,
 * skipped the conversion, and `after_global_event` then found no firing matrix:
 * the chain reached the enabled state and stopped there, filling the place to
 * the cutoff and reporting Tput 0 with QLen equal to the cutoff at every cutoff.
 * The matrices themselves are the test.
 */
template <class T>
bool nonmkv_firing_ready(const qn::TransitionParam<T>& np, std::size_t m) {
    if (m >= np.firingphases.size() || np.firingphases[m] <= 0) return false;
    if (m >= np.firingproc.size()) return false;
    return np.firingproc[m].D0.rows() == static_cast<std::size_t>(np.firingphases[m]);
}

}  // namespace detail

/**
 * Whether any law in the struct would be replaced, so a caller can skip copying
 * the struct when there is nothing to convert.
 *
 * @param sn           the struct to inspect
 * @param preserve_det leave Det alone, as the MAM analyzer does
 */
template <class T>
bool sn_has_nonmarkov(const qn::NetworkStruct<T>& sn, bool preserve_det = false) {
    using lang::ProcessType;
    for (std::size_t ist = 0; ist < sn.nstations; ++ist)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (ist < sn.disabled.size() && r < sn.disabled[ist].size() && sn.disabled[ist][r])
                continue;
            const lang::Distrib<T>& d = sn.service[ist][r];
            if (d.disabled || d.is_prior()) continue;
            if (detail::nonmkv_is_markovian(d.type) || detail::nonmkv_is_scheduled(d.type))
                continue;
            if (d.type == ProcessType::DET && preserve_det) continue;
            return true;
        }
    for (typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
             sn.transparam.begin();
         it != sn.transparam.end(); ++it) {
        const qn::TransitionParam<T>& np = it->second;
        for (std::size_t m = 0; m < np.firingproc.size(); ++m) {
            if (detail::nonmkv_firing_ready(np, m)) continue;
            const lang::Distrib<T>& d = np.firingproc[m];
            if (d.disabled || d.is_prior()) continue;
            if (detail::nonmkv_is_markovian(d.type) || detail::nonmkv_is_scheduled(d.type))
                continue;
            if (d.type == ProcessType::DET && preserve_det) continue;
            return true;
        }
    }
    return false;
}

/**
 * @param sn   the struct to convert IN PLACE; callers pass their own copy
 * @param opts the conversion method, order and fit family
 */
template <class T>
void sn_nonmarkov_toph(qn::NetworkStruct<T>& sn, const NonmarkovOptions& opts = NonmarkovOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "sn_nonmarkov_toph evaluates densities and fits moments");
    using lang::ProcessType;
    if (!opts.enabled) return;
    if (opts.order == 0) throw InputError("sn_nonmarkov_toph: the phase budget must be positive");

    for (std::size_t ist = 0; ist < sn.nstations; ++ist) {
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (ist < sn.disabled.size() && r < sn.disabled[ist].size() && sn.disabled[ist][r])
                continue;
            lang::Distrib<T>& d = sn.service[ist][r];
            if (d.disabled || d.is_prior()) continue;
            const ProcessType p = d.type;
            if (detail::nonmkv_is_markovian(p) || detail::nonmkv_is_scheduled(p)) continue;
            if (p == ProcessType::DET && opts.preserve_det) continue;

            const double mean = num_traits<T>::to_double(d.mean);
            const double scv = num_traits<T>::to_double(d.scv);
            if (!std::isfinite(mean) || mean <= 0.0) continue;

            if (p == ProcessType::DET) {
                detail::nonmkv_install(
                    d, mam::map_erlang(num_traits<T>::from_double(mean),
                                       static_cast<unsigned>(opts.order)));
                continue;
            }

            // The long-tail fit, when it was asked for and the law has a tail
            // to fit. It matches the ccdf at points spread over decades rather
            // than matching two moments, so it is the route for a Pareto, a
            // Weibull with shape below one or a Lognormal with a large sigma.
            if (opts.phfit == PhFit::Hyperexp) {
                const std::function<double(double)> ccdf = detail::nonmkv_ccdf(d);
                if (ccdf) {
                    bool ok = false;
                    const mam::Map<T> lt = detail::nonmkv_longtail<T>(ccdf, mean, ok);
                    if (ok) {
                        detail::nonmkv_install(d, lt);
                        continue;
                    }
                }
            }

            const std::function<double(double)> pdf = detail::nonmkv_density(d);
            if (!pdf) {  // no density to fit by shape: moments only
                detail::nonmkv_install(d,
                                       detail::nonmkv_concentrated<T>(mean, scv, opts.order,
                                                                      opts.phfit));
                continue;
            }
            if (opts.phfit == PhFit::Cme && scv >= 0.0 && scv < 1.0) {
                detail::nonmkv_install(d,
                                       detail::nonmkv_concentrated<T>(mean, scv, opts.order,
                                                                      opts.phfit));
                continue;
            }
            const mam::Map<T> fit =
                mam::map_scale(mam::map_bernstein<T>(pdf, static_cast<unsigned>(opts.order)),
                               num_traits<T>::from_double(mean));
            detail::nonmkv_install(d, fit);
        }
    }

    // ---- SPN transition firing laws -------------------------------------
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind) {
        if (sn.nodes[ind].nodetype != qn::NodeType::Transition) continue;
        typename std::map<std::size_t, qn::TransitionParam<T>>::iterator it =
            sn.transparam.find(ind + 1);
        if (it == sn.transparam.end()) continue;
        qn::TransitionParam<T>* np = &it->second;
        for (std::size_t m = 0; m < np->firingproc.size(); ++m) {
            // refreshPetriNetNodes already gives a Markovian mode a valid
            // (D0,D1) and a positive phase count; those are left alone. The
            // MATRICES decide, not the count -- see `nonmkv_firing_ready`.
            if (detail::nonmkv_firing_ready(*np, m)) continue;
            lang::Distrib<T>& d = np->firingproc[m];
            if (d.disabled || d.is_prior()) continue;
            const ProcessType p = d.type;
            if (detail::nonmkv_is_markovian(p) || detail::nonmkv_is_scheduled(p)) continue;
            if (p == ProcessType::DET && opts.preserve_det) continue;

            const double mean = num_traits<T>::to_double(d.mean);
            const double scv = num_traits<T>::to_double(d.scv);
            mam::Map<T> fit;
            if (p == ProcessType::DET || !std::isfinite(mean) || mean <= 0.0) {
                // A firing law with no usable mean falls back to the unit-mean
                // Erlang, as the reference does rather than refusing the model.
                const double mu = (std::isfinite(mean) && mean > 0.0) ? mean : 1.0;
                fit = mam::map_erlang(num_traits<T>::from_double(mu),
                                      static_cast<unsigned>(opts.order));
            } else {
                const std::function<double(double)> pdf = detail::nonmkv_density(d);
                if (!pdf || (opts.phfit == PhFit::Cme && scv >= 0.0 && scv < 1.0))
                    fit = detail::nonmkv_concentrated<T>(mean, scv, opts.order, opts.phfit);
                else
                    fit = mam::map_scale(
                        mam::map_bernstein<T>(pdf, static_cast<unsigned>(opts.order)),
                        num_traits<T>::from_double(mean));
            }
            detail::nonmkv_install(d, fit);
            if (m < np->firingphases.size()) np->firingphases[m] = fit.order();
        }
    }
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_NONMARKOV_TOPH_H
