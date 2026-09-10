/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SPN_SPN_MDD_H
#define LINE_API_SPN_SPN_MDD_H

/**
 * Decision-diagram reachable set and Kronecker rate descriptor of a stochastic
 * Petri net, so that `mdd::mdd_mcd` can analyse it.
 *
 * Port of matlab/src/api/spn/spn_mdd.m, jline.api.spn.Spn_mdd and
 * python/line_solver/api/spn/mdd.py.
 *
 * Levels are of two kinds. PLACE LEVELS hold a token count, one per Place.
 * PHASE LEVELS hold the phase the running server of a multi-phase mode
 * occupies, one per such mode.
 *
 * ONE LEVEL PER PLACE, NOT PER (PLACE, CLASS). The other three codebases carry
 * (nnodes x nclasses) arc matrices and so give a multiclass net P*R place
 * levels. This port reads `NetworkStruct::transparam`, whose enabling, firing
 * and inhibiting arcs are per (mode, NODE) with no class dimension -- which is
 * also how `state_events.h` evaluates them, summing the marking over classes
 * before every test. A multiclass net is therefore REJECTED here rather than
 * silently collapsed onto class-aggregated levels: the aggregation would be a
 * different model, not an approximation of this one. Single-class nets, which
 * is what an SPN normally is, agree level for level with the other codebases.
 *
 * The rate structure factorises exactly under single-server firing semantics: a
 * mode fires at a constant rate whenever every input level holds its enabling
 * multiplicity and no inhibitor level has reached its threshold, so
 * W_l^e[i, i + fire(l) - enab(l)] = 1 for enab(l) <= i < inhib(l) at every place
 * level. A phase-type mode contributes two event families on its phase level,
 * the internal phase changes D0 (marking unchanged) and the firings D1 (marking
 * moved), each gated by the same per-level enabling indicators. Both are
 * products of per-level terms, which is what Eq. 1 of the paper requires.
 *
 * PHASE-TYPE FIRING AND THE MEMORY POLICY. LINE discards a running server's
 * phase when its mode becomes disabled, i.e. preemptive repeat. Resetting a
 * mode's phase is then triggered by a JOINT condition on the place levels, which
 * is not a product of per-level terms and has no Kronecker form. What this
 * descriptor encodes is preemptive resume: a disabled mode's phase freezes and
 * continues when the mode is re-enabled. The two policies coincide exactly when
 * a mode is never disabled while running, so reachability records, for free,
 * whether any phase-type mode was ever found disabled. phmemory "exact" (the
 * default) errors when one was; "resume" proceeds deliberately with the resume
 * semantics.
 *
 * OTHER RESTRICTIONS, each an error and never a silent approximation: no
 * immediate transitions (they make vanishing states, which must be eliminated
 * before a Kronecker rate descriptor exists) and no marking-dependent firing
 * rates. A multi-server mode is accepted only when its enabling touches ONE
 * level, because the enabling degree min_l floor(m(l)/enab(l)) is otherwise not
 * a product of per-level terms.
 *
 * The initial marking comes from the reference station of each closed class
 * (`SpnOptions::init` overrides it): unlike the object-graph codebases, a
 * NetworkStruct carries no per-place state to read instead.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_types.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace spn {

/** One (transition, mode) pair of the net, in level coordinates. */
template <class T>
struct SpnMode {
    /** 1-based node index of the transition. */
    std::size_t trans = 0;
    /** Mode index within the transition, 0-based. */
    std::size_t mode = 0;
    /** Enabling multiplicity per place level. */
    std::vector<double> enab;
    /** Inhibition threshold per place level; infinite when absent. */
    std::vector<double> inhib;
    /** Firing outcome per place level. */
    std::vector<double> fire;
    Matrix<T> D0;
    Matrix<T> D1;
    std::vector<T> pie;
    std::size_t nph = 1;
    double srv = 1.0;
    /** Marking-dependent firing-rate multiplier; empty for the unit one. */
    std::function<T(const std::vector<T>&)> dep;
};

/** Everything the caller needs alongside the descriptor. */
template <class T>
struct SpnInfo {
    /** 1-based node indices of the places. */
    std::vector<std::size_t> places;
    std::vector<std::string> placenames;
    /** 1 for a place level, 2 for a phase level. */
    std::vector<int> levelkind;
    std::vector<std::string> levelname;
    std::vector<SpnMode<T>> modes;
    std::vector<int> init;
    mdd::MDD diagram;
    /** 1-based phase level of each mode, 0 when the mode has one phase. */
    std::vector<std::size_t> phaseof;
    /** Whether each mode was ever found disabled in a reachable marking. */
    std::vector<bool> ever_disabled;
    std::size_t nplacelevels = 0;
    /** Node count of the model, so a firingdep argument can be rebuilt. */
    std::size_t nnodes = 0;
    /** Whether the Kronecker descriptor was built. */
    bool descriptor = true;

    SpnInfo() : diagram(std::vector<int>(1, 1)) {}
};

/** Descriptor, diagram and metadata returned together. */
template <class T>
struct SpnResult {
    mdd::MddStruct mdds;
    mdd::MddDescriptor<T> desc;
    SpnInfo<T> info;
};

/** Options of the translation. */
struct SpnOptions {
    /** Per-place-level token bound; empty infers it from a place invariant. */
    std::vector<double> bound;
    /** Initial marking per place level; empty takes it from the reference stations. */
    std::vector<double> init;
    /** "exact" (default) or "resume". */
    std::string phmemory = "exact";
    /**
     * Build the Kronecker rate descriptor (default true). Pass false for the
     * MDD-rec route, which reads only the reachable set: the restrictions that
     * exist purely because a Kronecker form must factorise per level
     * (marking-dependent firing rates, multi-server modes drawing from several
     * places) are then lifted, in exchange for the firing times having to be
     * exponential.
     */
    bool descriptor = true;
};

namespace detail {

/**
 * Local matrix of one mode at place level l.
 *
 * Move the count by net, but only from local states that satisfy this level's
 * enabling and inhibition. apply_degree scales each row by the enabling degree
 * min(floor(m/enab), srv), i.e. the number of concurrently firing servers; it is
 * carried by the single enabling level of a multi-server mode. scale carries the
 * scalar firing rate on whichever level the caller chose to put it.
 */
template <class T>
mdd::MddLocalMatrix<T> spn_placemat(std::size_t l, const SpnMode<T>& mde, double net, int d,
                                    bool apply_degree, const T& scale) {
    typename mdd::MddLocalMatrix<T>::Builder bld(static_cast<std::size_t>(d));
    for (int i = 0; i < d; ++i) {
        if (!(static_cast<double>(i) >= mde.enab[l] && static_cast<double>(i) < mde.inhib[l]))
            continue;
        const int j = i + static_cast<int>(net);
        if (j < 0 || j > d - 1) continue;
        double val = 1.0;
        if (apply_degree && mde.enab[l] > 0) {
            const double deg = std::floor(static_cast<double>(i) / mde.enab[l]);
            val = deg < mde.srv ? deg : mde.srv;
        }
        bld.add(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                T(scale * num_traits<T>::from_double(val)));
    }
    return bld.build();
}

/** gcd of the integral entries of both vectors, 0 when any is not integral. */
inline long spn_gcd_vec(const std::vector<double>& a, const std::vector<double>& b) {
    long g = 0;
    const std::vector<double>* all[2] = {&a, &b};
    for (int i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < all[i]->size(); ++k) {
            const double x = (*all[i])[k];
            if (std::fabs(x - std::rint(x)) > 1e-9) return 0;
            long y = std::labs(static_cast<long>(std::rint(x)));
            while (y != 0) {
                const long t = g % y;
                g = y;
                y = t;
            }
        }
    return g;
}

/**
 * Drop every row whose support strictly contains another's, which is what leaves
 * the minimal supports and stops the pair expansion from blowing up.
 */
inline void spn_minimal_support(std::vector<std::vector<double>>& M,
                                std::vector<std::vector<double>>& B, std::size_t L) {
    const std::size_t n = B.size();
    std::vector<bool> drop(n, false);
    for (std::size_t i = 0; i < n; ++i) {
        if (drop[i]) continue;
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j || drop[j]) continue;
            bool contained = true, strict = false;
            for (std::size_t l = 0; l < L && contained; ++l) {
                const bool si = std::fabs(B[i][l]) > 1e-12;
                const bool sj = std::fabs(B[j][l]) > 1e-12;
                if (sj && !si) contained = false;
                if (si && !sj) strict = true;
            }
            if (contained && strict) {
                drop[i] = true;
                break;
            }
        }
    }
    std::vector<std::vector<double>> Mk, Bk;
    for (std::size_t i = 0; i < n; ++i)
        if (!drop[i]) {
            Mk.push_back(M[i]);
            Bk.push_back(B[i]);
        }
    M.swap(Mk);
    B.swap(Bk);
}

/**
 * A strictly positive place invariant w, i.e. w >= 0 with netm w = 0.
 *
 * Scanning a null-space BASIS for a positive vector is not enough, and the
 * fork-join net is the counterexample: its two minimal-support invariants are
 * (1,1,0,1) and (1,0,1,1), neither positive, while their sum (2,1,1,2) is. Which
 * basis a codebase's null space routine returns then decides whether the net is
 * accepted, which is how MATLAB and the JAR came to disagree on it. So the
 * non-negative generators are computed directly, by Farkas' algorithm on
 * [netm' | I] -- the same construction `spn_sinvariants` uses on the incidence
 * matrix -- and summed.
 */
inline bool spn_place_invariant(const std::vector<std::vector<double>>& netm, std::size_t L,
                                std::vector<double>& w) {
    const std::size_t E = netm.size();
    bool conservative = true;
    for (std::size_t e = 0; e < E; ++e) {
        double s = 0;
        for (std::size_t l = 0; l < L; ++l) s += netm[e][l];
        if (std::fabs(s) > 1e-12) {
            conservative = false;
            break;
        }
    }
    if (conservative) {
        w.assign(L, 1.0);
        return true;
    }
    if (E == 0) return false;

    std::vector<std::vector<double>> M(L, std::vector<double>(E, 0.0));  // level l's column
    std::vector<std::vector<double>> B(L, std::vector<double>(L, 0.0));  // its combination
    for (std::size_t l = 0; l < L; ++l) {
        for (std::size_t e = 0; e < E; ++e) M[l][e] = netm[e][l];
        B[l][l] = 1.0;
    }

    for (std::size_t e = 0; e < E; ++e) {
        std::vector<std::vector<double>> Mn, Bn;
        for (std::size_t i = 0; i < M.size(); ++i)
            if (std::fabs(M[i][e]) < 1e-12) {
                Mn.push_back(M[i]);
                Bn.push_back(B[i]);
            }
        for (std::size_t a = 0; a < M.size(); ++a) {
            if (!(M[a][e] > 1e-12)) continue;
            for (std::size_t b = 0; b < M.size(); ++b) {
                if (!(M[b][e] < -1e-12)) continue;
                const double ca = -M[b][e], cb = M[a][e];
                std::vector<double> c(E, 0.0), d(L, 0.0);
                for (std::size_t k = 0; k < E; ++k) c[k] = ca * M[a][k] + cb * M[b][k];
                for (std::size_t k = 0; k < L; ++k) d[k] = ca * B[a][k] + cb * B[b][k];
                const long g = spn_gcd_vec(c, d);
                if (g > 0) {
                    for (std::size_t k = 0; k < E; ++k) c[k] /= static_cast<double>(g);
                    for (std::size_t k = 0; k < L; ++k) d[k] /= static_cast<double>(g);
                }
                Mn.push_back(c);
                Bn.push_back(d);
            }
        }
        spn_minimal_support(Mn, Bn, L);
        M.swap(Mn);
        B.swap(Bn);
    }

    if (B.empty()) return false;
    std::vector<double> s(L, 0.0);
    for (std::size_t i = 0; i < B.size(); ++i)
        for (std::size_t l = 0; l < L; ++l) s[l] += B[i][l];
    double mn = std::numeric_limits<double>::max();
    for (std::size_t l = 0; l < L; ++l) {
        if (!(s[l] > 1e-9)) return false;
        mn = std::min(mn, s[l]);
    }
    w.assign(L, 0.0);
    for (std::size_t l = 0; l < L; ++l) w[l] = s[l] / mn;
    return true;
}

/** Successors of s, recording which modes were found disabled here. */
template <class T>
void spn_successors(const std::vector<int>& s, const std::vector<SpnMode<T>>& md,
                    const std::vector<std::vector<double>>& netm,
                    const std::vector<std::size_t>& phaseof, const std::vector<int>& domain,
                    std::size_t L, std::vector<std::vector<int>>& out,
                    std::vector<bool>& ever_disabled) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t e = 0; e < md.size(); ++e) {
        const SpnMode<T>& mde = md[e];
        bool enabled = true;
        for (std::size_t l = 0; l < L; ++l)
            if (static_cast<double>(s[l]) < mde.enab[l] ||
                static_cast<double>(s[l]) >= mde.inhib[l]) {
                enabled = false;
                break;
            }
        if (!enabled) {
            ever_disabled[e] = true;
            continue;
        }
        if (phaseof[e] == 0) {
            std::vector<int> t = s;
            bool ok = true;
            for (std::size_t l = 0; l < L; ++l) {
                const int nv = s[l] + static_cast<int>(netm[e][l]);
                if (nv < 0 || nv > domain[l] - 1) {
                    ok = false;
                    break;
                }
                t[l] = nv;
            }
            if (ok) out.push_back(t);
        } else {
            const std::size_t q = phaseof[e] - 1;
            const int ph = s[q];
            for (std::size_t j = 0; j < mde.nph; ++j)
                if (static_cast<int>(j) != ph && mde.D0(ph, j) != zero) {
                    std::vector<int> t = s;
                    t[q] = static_cast<int>(j);
                    out.push_back(t);
                }
            std::vector<int> moved = s;
            bool ok = true;
            for (std::size_t l = 0; l < L; ++l) {
                const int nv = s[l] + static_cast<int>(netm[e][l]);
                if (nv < 0 || nv > domain[l] - 1) {
                    ok = false;
                    break;
                }
                moved[l] = nv;
            }
            if (ok)
                for (std::size_t j = 0; j < mde.nph; ++j)
                    if (mde.D1(ph, j) != zero) {
                        std::vector<int> t = moved;
                        t[q] = static_cast<int>(j);
                        out.push_back(t);
                    }
        }
    }
}

}  // namespace detail

/**
 * Build the reachable set and Kronecker descriptor of a stochastic Petri net.
 *
 * @param sn a NetworkStruct holding Places and Transitions
 * @param options translation options
 */
template <class T>
SpnResult<T> spn_mdd(const qn::NetworkStruct<T>& sn, const SpnOptions& options = SpnOptions()) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const double inf = std::numeric_limits<double>::infinity();

    // ONE LEVEL PER PLACE is the decomposition this translation is built on
    // (see the header comment), so a COLOURED net -- whose level would have to
    // be a (place, class) pair -- has no descriptor here. The struct does carry
    // the class dimension since 2026-08-12; what is missing is the per-(place,
    // class) level, and the refusal names that rather than the storage.
    if (sn.nclasses > 1)
        throw UnsupportedError(
            "spn_mdd: the net has " + std::to_string(sn.nclasses) +
            " classes, and this translation puts ONE LEVEL PER PLACE; a coloured net needs a "
            "level per (place, class) pair, whose Kronecker descriptor is not the one built "
            "here. Solve the single-class net, or use SolverCTMC / SolverLDES, which evaluate "
            "the per-class arcs directly.");

    std::vector<std::size_t> places, transitions;
    for (std::size_t i = 1; i <= sn.nodes.size(); ++i) {
        if (sn.nodes[i - 1].nodetype == lang::NodeType::Place) places.push_back(i);
        else if (sn.nodes[i - 1].nodetype == lang::NodeType::Transition) transitions.push_back(i);
    }
    if (places.empty() || transitions.empty())
        throw InputError("spn_mdd: the model holds no Place or no Transition node");
    const std::size_t P = places.size();
    const std::size_t L = P;  // one level per place; see the header comment

    // ---- collect the (transition, mode) pairs
    std::vector<SpnMode<T>> md;
    for (std::size_t t = 0; t < transitions.size(); ++t) {
        const std::size_t ind = transitions[t];
        const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
            sn.transparam.find(ind);
        if (it == sn.transparam.end()) continue;
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            if (m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE)
                throw UnsupportedError("spn_mdd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " is IMMEDIATE; vanishing states must be eliminated before "
                                       "the net has a Kronecker rate descriptor");
            if (m < tp.firingdep.size() && tp.firingdep[m] && options.descriptor)
                throw UnsupportedError("spn_mdd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " has a marking-dependent firing rate; g(marking) is not a "
                                       "product of per-level terms");
            if (m >= tp.firingproc.size() || tp.firingproc[m].disabled)
                throw InputError("spn_mdd: mode " + std::to_string(m + 1) + " of node " +
                                 std::to_string(ind) + " has no firing process");
            SpnMode<T> e;
            e.trans = ind;
            e.mode = m;
            e.enab.assign(L, 0.0);
            e.inhib.assign(L, inf);
            e.fire.assign(L, 0.0);
            // The net is single class here (refused above), so the class-summed
            // arc IS the arc; `arc_total` is what says that out loud.
            const std::vector<T> en_t = qn::TransitionParam<T>::arc_total(tp.enabling, m);
            const std::vector<T> ih_t = qn::TransitionParam<T>::inhibit_total(tp.inhibiting, m);
            const std::vector<T> fi_t = qn::TransitionParam<T>::arc_total(tp.firing, m);
            for (std::size_t pp = 0; pp < P; ++pp) {
                const std::size_t q = places[pp] - 1;  // arcs are indexed by node
                if (q < en_t.size()) e.enab[pp] = std::max(0.0, num_traits<T>::to_double(en_t[q]));
                if (q < ih_t.size()) e.inhib[pp] = num_traits<T>::to_double(ih_t[q]);
                if (q < fi_t.size()) e.fire[pp] = std::max(0.0, num_traits<T>::to_double(fi_t[q]));
            }
            const mam::Map<T> proc = lang::dist_to_map(tp.firingproc[m]);
            e.D0 = proc.D0;
            e.D1 = proc.D1;
            e.nph = proc.order();
            e.pie = mdd::mdd_entry_law(std::vector<T>(), e.D1, e.nph, m, "spn_mdd");
            e.srv = m < tp.nmodeservers.size() ? tp.nmodeservers[m] : 1.0;
            if (m < tp.firingdep.size()) e.dep = tp.firingdep[m];
            if (!options.descriptor && e.nph > 1)
                throw UnsupportedError("spn_mdd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " has a phase-type firing time; the reachable-set-only "
                                       "mode carries no phase level, and a product-form marking "
                                       "process must be memoryless in the marking alone");
            md.push_back(e);
        }
    }
    const std::size_t E = md.size();
    for (std::size_t e = 0; options.descriptor && e < E; ++e) {
        std::size_t nz = 0;
        for (std::size_t l = 0; l < L; ++l)
            if (md[e].enab[l] != 0) ++nz;
        if (md[e].srv != 1 && nz > 1)
            throw UnsupportedError("spn_mdd: mode " + std::to_string(md[e].mode + 1) +
                                   " of node " + std::to_string(md[e].trans) + " has " +
                                   std::to_string(md[e].srv) + " servers and draws from " +
                                   std::to_string(nz) +
                                   " levels; the enabling degree min_l floor(m(l)/enab(l)) is then "
                                   "not a product of per-level terms and admits no Kronecker form");
        if (md[e].srv != 1 && nz == 0)
            throw InputError("spn_mdd: mode " + std::to_string(md[e].mode + 1) + " of node " +
                             std::to_string(md[e].trans) + " has " + std::to_string(md[e].srv) +
                             " servers but consumes from no place, so its enabling degree is "
                             "unbounded and its firing rate undefined");
    }

    // ---- phase levels for the multi-phase modes
    std::vector<std::size_t> phaseof(E, 0);  // 0 = none; else 1-based level
    std::size_t Q = 0;
    for (std::size_t e = 0; options.descriptor && e < E; ++e)
        if (md[e].nph > 1) {
            ++Q;
            phaseof[e] = L + Q;
        }
    const std::size_t K = L + Q;

    std::vector<std::vector<double>> netm(E, std::vector<double>(L, 0.0));
    for (std::size_t e = 0; e < E; ++e)
        for (std::size_t l = 0; l < L; ++l) netm[e][l] = md[e].fire[l] - md[e].enab[l];

    // ---- initial marking and per-level bounds
    std::vector<double> init0(L, 0.0);
    if (!options.init.empty()) {
        if (options.init.size() != L)
            throw InputError("spn_mdd: options.init must hold one token count per place");
        init0 = options.init;
    } else {
        for (std::size_t r = 0; r < sn.classes.size(); ++r) {
            const double njobs = sn.classes[r].population;
            if (!std::isfinite(njobs))
                throw UnsupportedError("spn_mdd: class " + std::to_string(r + 1) +
                                       " is open; an unbounded token population has no finite "
                                       "place level");
            const std::size_t ref_node = sn.station_to_node[sn.classes[r].refstat - 1];
            for (std::size_t pp = 0; pp < P; ++pp)
                if (places[pp] == ref_node) init0[pp] += njobs;
        }
    }
    std::vector<double> winv;
    const bool has_inv = detail::spn_place_invariant(netm, L, winv);
    double vinv = 0;
    if (has_inv)
        for (std::size_t l = 0; l < L; ++l) vinv += winv[l] * init0[l];

    std::vector<double> bound(L, 0.0);
    if (!options.bound.empty()) {
        for (std::size_t l = 0; l < L; ++l)
            bound[l] = options.bound.size() == 1 ? options.bound[0] : options.bound[l];
    } else if (has_inv) {
        for (std::size_t l = 0; l < L; ++l)
            bound[l] = winv[l] > 0 ? std::floor(vinv / winv[l]) : vinv;
    } else {
        throw InputError("spn_mdd: the net has no place invariant with positive weights, so the "
                         "marking is not bounded a priori; pass options.bound");
    }

    std::vector<int> domain(K, 1);
    for (std::size_t l = 0; l < L; ++l) domain[l] = static_cast<int>(bound[l]) + 1;
    for (std::size_t e = 0; e < E; ++e)
        if (phaseof[e] > 0) domain[phaseof[e] - 1] = static_cast<int>(md[e].nph);

    std::vector<int> init(K, 0);
    for (std::size_t l = 0; l < L; ++l) init[l] = static_cast<int>(init0[l]);
    for (std::size_t e = 0; e < E; ++e)
        if (phaseof[e] > 0)
            for (std::size_t a = 0; a < md[e].nph; ++a)
                if (md[e].pie[a] > zero) {
                    init[phaseof[e] - 1] = static_cast<int>(a);
                    break;
                }

    // ---- reachable set; the closure records which modes were ever disabled, so
    // the phase-memory question is answered without a second pass over |S|
    std::vector<bool> ever_disabled(E, false);
    mdd::MDD diagram(domain);
    diagram.insert(init);
    std::vector<std::vector<int>> frontier;
    frontier.push_back(init);
    std::size_t head = 0;
    while (head < frontier.size()) {
        const std::vector<int> s = frontier[head];
        ++head;
        std::vector<std::vector<int>> succ;
        detail::spn_successors(s, md, netm, phaseof, domain, L, succ, ever_disabled);
        for (std::size_t r = 0; r < succ.size(); ++r)
            if (!diagram.member(succ[r])) {
                diagram.insert(succ[r]);
                frontier.push_back(succ[r]);
            }
        if (head > 1024 && 2 * head > frontier.size()) {
            frontier.erase(frontier.begin(), frontier.begin() + static_cast<long>(head));
            head = 0;
        }
    }
    diagram.compact();

    for (std::size_t e = 0; e < E; ++e)
        if (phaseof[e] > 0 && options.descriptor && ever_disabled[e] &&
            options.phmemory != "resume" &&
            options.phmemory != "RESUME")
            throw UnsupportedError(
                "spn_mdd: mode " + std::to_string(md[e].mode + 1) + " of node " +
                std::to_string(md[e].trans) +
                " has a phase-type firing time AND is disabled in some reachable marking. LINE "
                "discards the phase on disabling (preemptive repeat) but that reset is a joint "
                "condition on the place levels and has no Kronecker form, so this descriptor would "
                "encode preemptive resume instead and disagree with SolverCTMC. Pass "
                "phmemory=\"resume\" to accept the resume semantics.");

    // ---- Kronecker event matrices
    mdd::MddDescriptor<T> desc;
    for (std::size_t e = 0; options.descriptor && e < E; ++e) {
        const SpnMode<T>& mde = md[e];
        std::set<std::size_t> gate, touched;
        for (std::size_t l = 0; l < L; ++l)
            if (mde.enab[l] > 0 || std::isfinite(mde.inhib[l])) {
                gate.insert(l);
                touched.insert(l);
            }
        for (std::size_t l = 0; l < L; ++l)
            if (netm[e][l] != 0) touched.insert(l);
        const std::vector<std::size_t> touched_sorted(touched.begin(), touched.end());
        long degl = -1;
        if (mde.srv != 1) {
            for (std::size_t l = 0; l < L; ++l)
                if (mde.enab[l] > 0) {
                    degl = static_cast<long>(l);  // level carrying the degree
                    break;
                }
        }
        if (phaseof[e] == 0) {
            if (touched_sorted.empty()) continue;
            mdd::MddEvent<T> event;
            event.a = mde.trans;
            event.b = mde.mode;
            for (std::size_t t = 0; t < touched_sorted.size(); ++t) {
                const std::size_t l = touched_sorted[t];
                event.lev.push_back(l);
                event.W.push_back(detail::spn_placemat(l, mde, netm[e][l], domain[l],
                                                       static_cast<long>(l) == degl, one));
            }
            // the scalar firing rate rides on the first touched level
            const std::size_t l0 = touched_sorted[0];
            event.W[0] = detail::spn_placemat(l0, mde, netm[e][l0], domain[l0],
                                              static_cast<long>(l0) == degl, mde.D1(0, 0));
            desc.events.push_back(event);
        } else {
            const std::size_t q = phaseof[e] - 1;
            // (1) internal phase changes: marking unchanged, gated by enabling
            std::size_t nnz_off = 0;
            for (std::size_t a = 0; a < mde.nph; ++a)
                for (std::size_t b = 0; b < mde.nph; ++b)
                    if (a != b && mde.D0(a, b) != zero) ++nnz_off;
            if (nnz_off > 0) {
                mdd::MddEvent<T> event;
                event.a = mde.trans;
                event.b = mde.mode;
                for (std::set<std::size_t>::const_iterator g = gate.begin(); g != gate.end();
                     ++g) {
                    event.lev.push_back(*g);
                    event.W.push_back(
                        detail::spn_placemat(*g, mde, 0.0, domain[*g], false, one));
                }
                typename mdd::MddLocalMatrix<T>::Builder bld(mde.nph);
                for (std::size_t a = 0; a < mde.nph; ++a)
                    for (std::size_t b = 0; b < mde.nph; ++b)
                        if (a != b) bld.add(a, b, mde.D0(a, b));
                event.lev.push_back(q);
                event.W.push_back(bld.build());
                desc.events.push_back(event);
            }
            // (2) firings: marking moved, phase redrawn through D1
            mdd::MddEvent<T> event;
            event.a = mde.trans;
            event.b = mde.mode;
            for (std::size_t t = 0; t < touched_sorted.size(); ++t) {
                const std::size_t l = touched_sorted[t];
                event.lev.push_back(l);
                event.W.push_back(detail::spn_placemat(l, mde, netm[e][l], domain[l],
                                                       static_cast<long>(l) == degl, one));
            }
            typename mdd::MddLocalMatrix<T>::Builder bld(mde.nph);
            for (std::size_t a = 0; a < mde.nph; ++a)
                for (std::size_t b = 0; b < mde.nph; ++b) bld.add(a, b, mde.D1(a, b));
            event.lev.push_back(q);
            event.W.push_back(bld.build());
            desc.events.push_back(event);
        }
    }

    desc.K = K;
    desc.N = 0;  // the invariant below replaces it
    desc.domain = domain;
    if (has_inv) {
        desc.invariant_weights.assign(K, 0.0);
        for (std::size_t l = 0; l < L; ++l) desc.invariant_weights[l] = winv[l];
        desc.invariant_value = vinv;
    }

    // ---- descriptive information
    SpnInfo<T> info;
    info.places = places;
    for (std::size_t pp = 0; pp < P; ++pp) info.placenames.push_back(sn.nodes[places[pp] - 1].name);
    info.levelkind.assign(K, 2);
    for (std::size_t l = 0; l < L; ++l) info.levelkind[l] = 1;
    info.levelname.assign(K, std::string());
    for (std::size_t pp = 0; pp < P; ++pp) info.levelname[pp] = info.placenames[pp];
    for (std::size_t e = 0; e < E; ++e)
        if (phaseof[e] > 0)
            info.levelname[phaseof[e] - 1] = "phase(" + sn.nodes[md[e].trans - 1].name + ".m" +
                                             std::to_string(md[e].mode + 1) + ")";
    info.modes = md;
    info.init = init;
    info.diagram = diagram;
    info.phaseof = phaseof;
    info.ever_disabled = ever_disabled;
    info.nplacelevels = L;
    info.nnodes = sn.nodes.size();
    info.descriptor = options.descriptor;

    SpnResult<T> res;
    res.mdds = diagram.to_struct();
    res.desc = desc;
    res.info = info;
    return res;
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_MDD_H
