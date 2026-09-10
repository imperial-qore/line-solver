/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_PREDICATES_H
#define LINE_API_SN_SN_PREDICATES_H

/**
 * Ports of the sn_has_* / sn_is_* predicate family of matlab/src/api/sn.
 *
 * The reference keeps each predicate in its own one-line .m file; collecting
 * them in one header keeps the port readable and, more importantly, keeps the
 * DEFINITIONS in one place. Several of these predicates gate solver dispatch,
 * so a predicate that quietly differs from the reference does not raise an
 * error -- it routes the model to a different algorithm and returns a
 * different number.
 *
 * NetworkStruct carries member functions with some of these names. Where both
 * exist the free function here is the reference-faithful one and the member is
 * a convenience used inside the struct's own refresh; sn_has_product_form is
 * the case where the two DISAGREE, see the note on that function.
 * sn_has_blocking and sn_is_mm1k_loss are the other way round: the SOLVERS read
 * them through NetworkStruct::has_product_form, so the body lives on the struct
 * and these forward to it rather than carrying a second copy.
 *
 * ARITHMETIC: field. Every test is a comparison; nothing here is
 * transcendental, so the whole family instantiates under Rational.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/sn/sn_get_buffer_size.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

namespace detail {

/** True when some station runs `s`, the shape shared by the sn_has_<sched> family. */
template <class T>
bool sn_any_sched(const qn::NetworkStruct<T>& sn, qn::SchedStrategy s) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        if (sn.stations[i].sched == s) return true;
    return false;
}

}  // namespace detail

// ---- the sn_has_<discipline> family ------------------------------------

template <class T>
bool sn_has_dps(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::DPS);
}
template <class T>
bool sn_has_dps_prio(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::DPSPRIO);
}
template <class T>
bool sn_has_fcfs(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::FCFS);
}
template <class T>
bool sn_has_gps(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::GPS);
}
template <class T>
bool sn_has_gps_prio(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::GPSPRIO);
}
template <class T>
bool sn_has_hol(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::HOL);
}
template <class T>
bool sn_has_inf(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::INF);
}
template <class T>
bool sn_has_lcfs(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LCFS);
}
template <class T>
bool sn_has_lcfs_pi(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LCFSPI);
}
template <class T>
bool sn_has_lcfs_pr(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LCFSPR);
}
template <class T>
bool sn_has_lept(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LEPT);
}
template <class T>
bool sn_has_ljf(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LJF);
}
template <class T>
bool sn_has_lps(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::LPS);
}
template <class T>
bool sn_has_polling(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::POLLING);
}
template <class T>
bool sn_has_ps(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::PS);
}
template <class T>
bool sn_has_ps_prio(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::PSPRIO);
}
template <class T>
bool sn_has_sept(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::SEPT);
}
template <class T>
bool sn_has_setf(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::SETF);
}
template <class T>
bool sn_has_siro(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::SIRO);
}
template <class T>
bool sn_has_sjf(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::SJF);
}
template <class T>
bool sn_has_srpt(const qn::NetworkStruct<T>& sn) {
    return detail::sn_any_sched(sn, qn::SchedStrategy::SRPT);
}

/**
 * Port of sn_has_homogeneous_scheduling.
 *
 * The reference is `length(findstring(sn.sched, strategy)) == sn.nstations`,
 * and findstring compares STRINGS: against the numeric sched vector its strcmp
 * is false everywhere, so it returns the sentinel -1, whose length is 1. The
 * predicate therefore reduces to nstations == 1 whatever the disciplines are.
 * Reproduced and not corrected: this is the branch every codebase takes.
 */
template <class T>
bool sn_has_homogeneous_scheduling(const qn::NetworkStruct<T>& sn, qn::SchedStrategy) {
    return sn.nstations == 1;
}

// ---- class and chain counts -------------------------------------------

template <class T>
bool sn_has_open_classes(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (std::isinf(sn.classes[k].population)) return true;
    return false;
}

template <class T>
bool sn_has_closed_classes(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (std::isfinite(sn.classes[k].population)) return true;
    return false;
}

template <class T>
bool sn_has_multiple_closed_classes(const qn::NetworkStruct<T>& sn) {
    std::size_t n = 0;
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (std::isfinite(sn.classes[k].population)) ++n;
    return n > 1;
}

template <class T>
bool sn_has_mixed_classes(const qn::NetworkStruct<T>& sn) {
    return sn_has_closed_classes(sn) && sn_has_open_classes(sn);
}

template <class T>
bool sn_has_multi_chain(const qn::NetworkStruct<T>& sn) {
    return sn.nchains > 1;
}
template <class T>
bool sn_has_single_chain(const qn::NetworkStruct<T>& sn) {
    return sn.nchains == 1;
}
template <class T>
bool sn_has_multi_class(const qn::NetworkStruct<T>& sn) {
    return sn.nclasses > 1;
}
template <class T>
bool sn_has_single_class(const qn::NetworkStruct<T>& sn) {
    return sn.nclasses == 1;
}
template <class T>
bool sn_has_class_switching(const qn::NetworkStruct<T>& sn) {
    return sn.nclasses != sn.nchains;
}

/** `all(isinf(sn.njobs))`: EVERY class is open, which a mixed model fails. */
template <class T>
bool sn_is_open_model(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (!std::isinf(sn.classes[k].population)) return false;
    return true;
}

/** `all(isfinite(sn.njobs))`: EVERY class is closed. */
template <class T>
bool sn_is_closed_model(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (!std::isfinite(sn.classes[k].population)) return false;
    return true;
}

template <class T>
bool sn_is_mixed_model(const qn::NetworkStruct<T>& sn) {
    return sn_has_mixed_classes(sn);
}

// ---- structural predicates --------------------------------------------

template <class T>
bool sn_has_priorities(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k)
        if (sn.classes[k].prio > 0) return true;
    return false;
}

template <class T>
bool sn_has_multi_server(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        if (std::isfinite(sn.stations[i].nservers) && sn.stations[i].nservers > 1.0) return true;
    return false;
}

/** `any(sn.njobs ~= round(sn.njobs))`; an infinite population rounds to itself. */
template <class T>
bool sn_has_fractional_populations(const qn::NetworkStruct<T>& sn) {
    for (std::size_t k = 0; k < sn.classes.size(); ++k) {
        const double n = sn.classes[k].population;
        if (std::isfinite(n) && n != std::floor(n + 0.5)) return true;
    }
    return false;
}

/** `size(sn.lldscaling,2)>0`: some station carries a limited-load scaling row. */
template <class T>
bool sn_has_load_dependence(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i)
        if (!sn.stations[i].lldscaling.empty()) return true;
    return false;
}

/**
 * `any(sn.fj(:) > 0)`: the model has a Fork closed by a Join.
 *
 * `sn.fj` is the (fork, join) pair list here, so non-empty is the predicate.
 * It is populated on the UN-augmented struct too, exactly as MATLAB's fj
 * matrix is, so this agrees with the reference before and after fj_tag.
 */
template <class T>
bool sn_has_fork_join(const qn::NetworkStruct<T>& sn) {
    return !sn.fj.empty();
}

/**
 * `sn_has_sd_routing`: some node dispatches on the STATE of the network.
 *
 * RROBIN, WRROBIN, JSQ and SQ. RROBIN matters even though the refresh spreads
 * its probabilities uniformly: the uniform matrix is the right MEAN and the
 * wrong higher moments, so a round-robin model is not product form and must
 * not be handed to exact MVA or to the convolution algorithm.
 */
template <class T>
bool sn_has_sd_routing(const qn::NetworkStruct<T>& sn) {
    for (std::size_t a = 0; a < sn.nodes.size(); ++a)
        for (std::size_t r = 0; r < sn.nodes[a].routing.size(); ++r) {
            const qn::RoutingStrategy s = sn.nodes[a].routing[r];
            if (s == qn::RoutingStrategy::RROBIN || s == qn::RoutingStrategy::WRROBIN ||
                s == qn::RoutingStrategy::JSQ || s == qn::RoutingStrategy::SQ ||
                s == qn::RoutingStrategy::SDR)
                return true;
        }
    return false;
}

/** `sum(sn.rates(i,:)>0)>1` at some FCFS station: more than one class is served there. */
template <class T>
bool sn_has_multi_class_fcfs(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched != qn::SchedStrategy::FCFS) continue;
        std::size_t n = 0;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            if (sn.rates(i, r) > zero) ++n;
        }
        if (n > 1) return true;
    }
    return false;
}

/**
 * `sn_has_multi_class_heter_fcfs`: an FCFS station whose per-class rates differ.
 *
 * MATLAB filters the row with `isfinite(rates) & ~isnan(rates)`, which drops
 * BOTH the NaN of a disabled pair and the Inf of an Immediate service. The
 * disabled flag covers the first; the isfinite test below covers the second,
 * without which an Immediate-served class would make every FCFS station look
 * heterogeneous and take product form away from the model.
 */
template <class T>
bool sn_has_multi_class_heter_fcfs(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched != qn::SchedStrategy::FCFS) continue;
        bool any = false;
        double lo = 0.0, hi = 0.0;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            const double v = num_traits<T>::to_double(sn.rates(i, r));
            if (!std::isfinite(v)) continue;
            if (!any) {
                lo = hi = v;
                any = true;
            } else {
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
        }
        if (any && hi - lo > 0.0) return true;
    }
    return false;
}

/** An FCFS station with heterogeneous rates whose SCVs are all one (exponential). */
template <class T>
bool sn_has_multi_class_heter_exp_fcfs(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched != qn::SchedStrategy::FCFS) continue;
        // MATLAB tests `range(sn.rates(i,:)) > 0`, and range's max and min skip
        // NaN but keep Inf -- unlike sn_has_multi_class_heter_fcfs, which
        // filters on isfinite. The two predicates really do differ there.
        bool any = false;
        double lo = 0.0, hi = 0.0;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            const double v = num_traits<T>::to_double(sn.rates(i, r));
            if (std::isnan(v)) continue;
            if (!any) {
                lo = hi = v;
                any = true;
            } else {
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
        }
        if (!(any && hi - lo > 0.0)) continue;
        // max/min over the scv row also skip NaN, so a disabled pair neither
        // qualifies nor disqualifies the station; an all-disabled row leaves
        // max([]) empty and the reference's `if` false, hence `seen`.
        bool allone = true, seen = false;
        for (std::size_t r = 0; r < sn.nclasses && allone; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            const double c = num_traits<T>::to_double(sn.scv(i, r));
            if (std::isnan(c)) continue;
            seen = true;
            if (!(c < 1.0 + qn::GlobalConstants::FineTol && c > 1.0 - qn::GlobalConstants::FineTol))
                allone = false;
        }
        if (seen && allone) return true;
    }
    return false;
}

/** The disciplines the reference admits in sn_has_product_form. */
template <class T>
bool sn_sched_is_product_form(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        const qn::SchedStrategy s = sn.stations[i].sched;
        if (!(s == qn::SchedStrategy::INF || s == qn::SchedStrategy::PS ||
              s == qn::SchedStrategy::FCFS || s == qn::SchedStrategy::LCFSPR ||
              s == qn::SchedStrategy::LCFS || s == qn::SchedStrategy::EXT))
            return false;
    }
    return true;
}

/** Defined below, after sn_is_mm1k_loss, which it exempts. */
template <class T>
bool sn_has_blocking(const qn::NetworkStruct<T>& sn);

/**
 * Port of sn_has_product_form.
 *
 * NOTE. `NetworkStruct::has_product_form` omits the fork-join and the
 * state-dependent-routing conjuncts the reference carries, so it answers TRUE
 * on a model with a surviving Fork or a round-robin dispatcher. That is not a
 * cosmetic difference: the predicate gates exact MVA, the NC convolution and
 * the AUTO solver tree, so the omission silently hands a non-product-form
 * model to an exact algorithm. This free function is the reference-faithful
 * one and callers in the api layer must use it.
 */
template <class T>
bool sn_has_product_form(const qn::NetworkStruct<T>& sn) {
    // BCMP asks for infinite buffers. Nothing here read sn.cap/sn.classcap/
    // sn.droprule, so a BAS-blocked station or any binding finite buffer passed the
    // gate and the network read as product form while its truncation couples the
    // station occupancies.
    return sn_sched_is_product_form(sn) && !sn_has_multi_class_heter_fcfs(sn) &&
           !sn_has_priorities(sn) && !sn_has_fork_join(sn) && !sn_has_sd_routing(sn) &&
           !sn_has_blocking(sn) && sn.has_exponential_fcfs();
}

/**
 * Port of sn_has_product_form_not_het_fcfs: LCFS is out, and every enabled
 * FCFS pair with a finite positive SCV must have that SCV equal to one, with
 * the service means agreeing across the classes served there.
 *
 * CHECK_MEANS drops that second half; pass false only for an algorithm that
 * models class-dependent FCFS itself (ab, schmidt, schmidt-ext).
 */
template <class T>
bool sn_has_product_form_not_het_fcfs(const qn::NetworkStruct<T>& sn, bool check_means = true) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        const qn::SchedStrategy s = sn.stations[i].sched;
        if (!(s == qn::SchedStrategy::INF || s == qn::SchedStrategy::PS ||
              s == qn::SchedStrategy::FCFS || s == qn::SchedStrategy::LCFSPR ||
              s == qn::SchedStrategy::EXT))
            return false;
    }
    if (sn_has_priorities(sn) || sn_has_fork_join(sn) || sn_has_sd_routing(sn)) return false;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched != qn::SchedStrategy::FCFS) continue;
        // BCMP type 1 asks the FCFS service to be exponential AND
        // class-independent, so the means are checked alongside the SCVs: with
        // unequal means the product-form solve returns a wait proportional to
        // each class's own demand where FCFS makes every class wait behind the
        // same queue. The mean comparison is between CHAIN service times
        // (visit-weighted over the classes that actually visit the station): a
        // class that never visits cannot break product form, and within-chain
        // heterogeneity is invisible to both the product-form and the qd
        // branch, which deaggregate a chain result proportionally to each
        // class's own demand, so only between-chain heterogeneity warrants the
        // divert. LN layers carry seeded rates for classes with zero visits,
        // which a raw per-class comparison mistakes for heterogeneity.
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            const double v = num_traits<T>::to_double(sn.scv(i, r));
            if (std::isfinite(v) && v > 0.0 &&
                !(v > 1.0 - qn::GlobalConstants::FineTol && v < 1.0 + qn::GlobalConstants::FineTol))
                return false;
        }
        if (!check_means || sn.visits.empty()) continue;
        const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
        double stmin = 0.0, stmax = 0.0;
        bool anyserved = false;
        for (std::size_t c = 0; c < sn.nchains && c < sn.visits.size(); ++c) {
            double num = 0.0, den = 0.0;
            for (std::size_t r = 0; r < sn.nclasses; ++r) {
                if (!sn.chains.empty() && !sn.chains[c][r]) continue;
                if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
                const double w = num_traits<T>::to_double(sn.visits[c](isf, r));
                const double rate = num_traits<T>::to_double(sn.rates(i, r));
                if (w > qn::GlobalConstants::Zero && std::isfinite(rate) && rate > 0.0) {
                    num += w / rate;
                    den += w;
                }
            }
            if (den > 0.0) {
                const double st = num / den;
                if (!anyserved) {
                    stmin = stmax = st;
                    anyserved = true;
                } else {
                    stmin = std::min(stmin, st);
                    stmax = std::max(stmax, st);
                }
            }
        }
        if (anyserved && stmax - stmin > qn::GlobalConstants::CoarseTol * stmax) return false;
    }
    return true;
}

/**
 * `sn_is_population_model`: every station is population dependent only, which
 * is what lets the fluid and mean-field limits close on the queue lengths.
 */
template <class T>
bool sn_is_population_model(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        const qn::SchedStrategy s = sn.stations[i].sched;
        if (!(s == qn::SchedStrategy::INF || s == qn::SchedStrategy::PS ||
              s == qn::SchedStrategy::PSPRIO || s == qn::SchedStrategy::DPS ||
              s == qn::SchedStrategy::GPS || s == qn::SchedStrategy::GPSPRIO ||
              s == qn::SchedStrategy::DPSPRIO || s == qn::SchedStrategy::EXT))
            return false;
    }
    return !sn_has_priorities(sn) && !sn_has_fork_join(sn);
}

/**
 * `sn_is_bas_model`: a single-class closed model with blocking-after-service.
 */
template <class T>
bool sn_is_bas_model(const qn::NetworkStruct<T>& sn) {
    if (sn.nclasses != 1) return false;
    double closed = 0.0;
    for (std::size_t k = 0; k < sn.classes.size(); ++k) {
        if (std::isinf(sn.classes[k].population)) return false;  // open class present
        closed += sn.classes[k].population;
    }
    if (!(closed > 0.0)) return false;
    // `sn.droprule`, the REFRESHED (nstations x nclasses) table, not the
    // per-station user declaration: refresh_bas_blocking propagates a
    // destination-declared BAS rule into the table, and reading the
    // declaration would miss exactly those models
    for (std::size_t i = 0; i < sn.droprule.size(); ++i)
        for (std::size_t r = 0; r < sn.droprule[i].size(); ++r)
            if (sn.droprule[i][r] == qn::DropStrategy::BAS) return true;
    return false;
}

/**
 * `sn_is_mm1k_loss`: the three-node Source/Queue/Sink model of an M/M/1/K with
 * loss, the shape MVA answers in closed form instead of iterating.
 */
template <class T>
bool sn_is_mm1k_loss(const qn::NetworkStruct<T>& sn) {
    return sn.is_mm1k_loss();
}

/**
 * Port of sn_has_blocking: some station can REFUSE a job.
 *
 * Either its own buffer binds -- Kendall's K below the population that can
 * reach it, whatever the drop rule (WAITQ, DROP, BAS, BBS, RSRD) -- or a finite
 * capacity region caps a set of stations jointly. Such a network is not product
 * form: the truncation couples the station occupancies, so no BCMP
 * factorization of the equilibrium distribution exists.
 *
 * Only a buffer that can actually BIND counts, which is what
 * sn_get_buffer_size decides; refresh_capacity derives a finite classcap (the
 * chain population) at every station of every closed model, so a plain
 * finiteness test would call every closed model blocking.
 *
 * A Cache is exempt: it builds its own capped retrieval queues (classcap 1),
 * which the cache analyzers solve rather than treat as a buffer constraint. So
 * is the single-station M/M/1/K loss system, whose truncated geometric
 * distribution is a product form over its one station.
 */
template <class T>
bool sn_has_blocking(const qn::NetworkStruct<T>& sn) {
    return sn.has_blocking();
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_PREDICATES_H
