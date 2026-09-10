#pragma once
/**
 * @file spn_lpbnd.h
 * @brief Linear-programming bounds on the mean marking and the throughputs of a
 *        stochastic timed Petri net.
 *
 * The stationary chain is relaxed to a MOMENT POLYTOPE: the uniformized
 * evolution equation is written for E[X_p], E[X_p^2] and E[X_p1 X_p2], which
 * gives linear equalities among the mean marking x, the enabling probabilities
 * q and the products y(p,t) = E[X_p e_t]; behavioural and probabilistic
 * inequalities are added on top; and every reported measure is then obtained by
 * minimising and maximising its linear form over that polytope. Any stationary
 * point of the true chain satisfies every row, so the two optima BRACKET the
 * exact value whatever the polytope leaves out.
 *
 * This is the Petri-net sibling of the QRF bounds in SolverBA: same technique,
 * a different index space, and a LINEAR objective, so there is no stationary
 * point to escape from and the answer is a property of the model alone.
 *
 * VARIABLES, over place levels l and modes e: x(l) the mean tokens, q(e) the
 * probability that mode e is enabled, th(e) its throughput, u(e) the
 * state-equation firing counts, and y(l,e) = E[X_l e_e] on the Markovian side
 * only. u is EXISTENTIAL and is not reported: E[X] is a convex combination of
 * reachable markings, each of which is m0 + C h for some nonnegative integer h,
 * so the mean satisfies m0 + C u for some nonnegative real u.
 *
 * ONE LEVEL PER PLACE, NOT PER (PLACE, CLASS), and a multiclass net is REFUSED.
 * The other three codebases carry (nnodes x nclasses) arc matrices and so give
 * a multiclass net P*R place levels. This port reads
 * `NetworkStruct::transparam`, whose arcs are per (mode, NODE) with no class
 * dimension, the same restriction `spn_mdd` states and enforces. A coloured net
 * is a different model, not an approximation of this one, so it is refused
 * rather than collapsed. `spn_sinvariants` is class-summed here for the same
 * reason, which makes the invariant family (8) below the coarser one -- still a
 * genuine invariant, just not the finest.
 *
 * THE INITIAL MARKING COMES FROM THE REFERENCE STATION of each closed class, or
 * from `SpnLpOptions::init`: unlike the object-graph codebases, a NetworkStruct
 * carries no per-place state to read instead. A net whose tokens do not all
 * start at the reference station must pass `init`.
 *
 * THE TOKEN COUNTS CREATED BY A FIRING ARE DETERMINISTIC IN LINE, which removes
 * a whole branch of the reference: it allows sigma_(t,p)(n) to be random and
 * splits the covariance family into an independent case (its eq. 7) and a
 * selective one (its eq. 8). A firing outcome is an integer weight, so
 * E[sigma^2] = sigma^2 and E[sigma_p1 sigma_p2] = sigma_p1 sigma_p2 hold
 * exactly and eq. (7) is the correct form. Eq. (8) has no LINE model behind it
 * and is deliberately absent.
 *
 * LIVENESS IS OFF BY DEFAULT, AND THAT IS DELIBERATE. The reference's two
 * liveness rows (sum_t q_t >= 1 and x_p <= sum_t y_(p,t)) hold only on a live
 * net, and liveness is not something this function can cheaply certify -- an
 * inhibitor arc alone is enough to deadlock a net that looks well formed. A
 * bound that silently assumed it would be wrong rather than loose on exactly
 * the models where a bound is most wanted, so the rows are opt-in.
 *
 * WHAT THE ROWS ARE WORTH, MEASURED. They are the whole of the lower side. On
 * the reference's own Table 2 (its Fig. 2b production line, five rate vectors)
 * `assumelive` reproduces its published l.b. column to four decimals -- 1.1653
 * against 1.165, 1.8288 against 1.829, 1.5814 against 1.581, 1.3592 against
 * 1.359, 1.3497 against 1.350 -- while without them the Markovian lower bound
 * collapses onto the OPERATIONAL one on four of the five. The upper side needs
 * neither row and matches the published u.b.2 either way.
 *
 * Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
 * Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
 * 1014-1030. The constraint families are its Table 1, p. 1022; the bracket
 * statement is its Theorem 3, p. 1021.
 *
 * MATLAB twin: `spn_lpbnd.m`. JAR twin: `Spn_lpbnd.java`. Python twin:
 * `api/spn/lpbnd.py`.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/spn/spn_sinvariants.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/lp_highs.h"
#include "line/util/simplex.h"

namespace line {
namespace spn {

/** Options of the relaxation. */
struct SpnLpOptions {
    /**
     * true (the default) uses the second-moment, covariance and Little's law
     * families, which need exponential firing times; false drops them and the
     * whole y block, leaving the operational bound, which needs only a mean
     * firing time and so admits any phase-type law.
     */
    bool markovian = true;
    /** true adds the two liveness rows, valid only on a live net. */
    bool assumelive = false;
    /** Initial tokens per place level; empty takes the reference stations. */
    std::vector<double> init;
    /** Slack added to the inequality sides. */
    double tol = 0.0;
};

/** One (transition, mode) pair over place levels. */
struct SpnLpMode {
    std::size_t trans = 0;   ///< 1-based node index of the transition
    std::size_t mode = 0;    ///< mode index within the transition, 0-based
    std::vector<double> enab;
    std::vector<double> inhib;
    std::vector<double> fire;
    double rate = 0.0;
};

/** The brackets; each vector pair holds the minimum then the maximum. */
struct SpnLpBounds {
    std::vector<std::size_t> places;     ///< 1-based node indices, in level order
    std::vector<std::string> levelname;
    std::vector<SpnLpMode> modes;
    std::vector<double> tokens_lo, tokens_hi;
    std::vector<double> place_tput_lo, place_tput_hi;
    std::vector<double> mode_tput_lo, mode_tput_hi;
    std::vector<double> mode_util_lo, mode_util_hi;
    std::vector<double> bound;           ///< a priori per-level cap, inf where none
    std::size_t nplacelevels = 0;
    bool markovian = true;
    std::size_t nvars = 0;
    std::size_t nrows = 0;
};

namespace detail {

/** One LP. An infeasible or unbounded program answers NaN, as MATLAB does. */
inline double spn_lp_opt(lp::LpModel<double>& m, const std::vector<double>& c, bool minimize) {
    for (std::size_t j = 0; j < m.num_vars(); ++j) m.set_cost(j, c[j]);
    m.set_maximize(!minimize);
    const lp::LpSolution<double> sol = lp::lp_solve(m);
    if (!sol.ok() || !std::isfinite(sol.objective)) return std::numeric_limits<double>::quiet_NaN();
    return sol.objective;
}

}  // namespace detail

/**
 * Bracket the mean tokens and the throughputs of a stochastic Petri net.
 *
 * @param sn a NetworkStruct holding Places and Transitions, single class
 * @param options the relaxation options
 * @return the brackets, per place level and per mode
 */
template <class T>
SpnLpBounds spn_lpbnd(const qn::NetworkStruct<T>& sn,
                      const SpnLpOptions& options = SpnLpOptions()) {
    const double inf = std::numeric_limits<double>::infinity();

    // ONE LEVEL PER PLACE, so a coloured net is a different model. Same refusal
    // spn_mdd makes, and for the same reason.
    if (sn.nclasses > 1)
        throw UnsupportedError(
            "spn_lpbnd: the net has " + std::to_string(sn.nclasses) +
            " classes, and this translation puts ONE LEVEL PER PLACE; a coloured net needs a "
            "level per (place, class) pair, which the moment relaxation would then be written "
            "over. Solve the single-class net, or use the MATLAB, JAR or python twin, which "
            "carry the class dimension.");

    std::vector<std::size_t> places, transitions;
    for (std::size_t i = 1; i <= sn.nodes.size(); ++i) {
        if (sn.nodes[i - 1].nodetype == lang::NodeType::Place) places.push_back(i);
        else if (sn.nodes[i - 1].nodetype == lang::NodeType::Transition) transitions.push_back(i);
    }
    if (places.empty() || transitions.empty())
        throw InputError("spn_lpbnd: the model holds no Place or no Transition node");
    const std::size_t P = places.size();
    const std::size_t L = P;

    // ---- the (transition, mode) table. spn_mdd builds the same one, but only
    // as a step of reachable-set construction, which is the cost this bound
    // exists to avoid.
    std::vector<SpnLpMode> md;
    for (std::size_t t = 0; t < transitions.size(); ++t) {
        const std::size_t ind = transitions[t];
        const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
            sn.transparam.find(ind);
        if (it == sn.transparam.end()) continue;
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            if (m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE)
                throw UnsupportedError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " is IMMEDIATE; the moment relaxation is written for a net "
                                       "whose transitions all have finite rates, so vanishing "
                                       "states must be eliminated first");
            if (m < tp.firingdep.size() && tp.firingdep[m])
                throw UnsupportedError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " has a marking-dependent firing rate; the uniformization "
                                       "step needs one rate per mode");
            const double srv = m < tp.nmodeservers.size() ? tp.nmodeservers[m] : 1.0;
            if (srv != 1.0)
                throw UnsupportedError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) + " has " + std::to_string(srv) +
                                       " servers; the relaxation is derived under single-server "
                                       "semantics, where the firing rate is mu*q. Its "
                                       "infinite-server form needs the K-fold transition "
                                       "expansion of the reference's Section 7, which is not "
                                       "implemented");
            if (m >= tp.firingproc.size() || tp.firingproc[m].disabled)
                throw InputError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                 std::to_string(ind) + " has no firing process");
            SpnLpMode e;
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
            // A PHASE-TYPE FIRING LAW IS WHERE THE TWO VARIANTS PART. The mean
            // of a (D0,D1) pair is all the operational bound needs; the
            // Markovian one needs the marking alone to be the state, which a
            // multi-phase mode breaks.
            const mam::Map<T> proc = lang::dist_to_map(tp.firingproc[m]);
            const std::size_t nph = proc.order();
            if (options.markovian && nph > 1)
                throw UnsupportedError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                       std::to_string(ind) +
                                       " has a phase-type firing time; the relaxation is written "
                                       "over the marking alone, and a phase-type mode needs the "
                                       "state-machine expansion of the reference's Section 7, "
                                       "which is not implemented. Use the operational bound, "
                                       "which needs only the mean");
            const double mean = num_traits<T>::to_double(mam::map_mean(proc));
            e.rate = nph == 1 ? num_traits<T>::to_double(proc.D1(0, 0)) : 1.0 / mean;
            if (!(e.rate > 0) || !std::isfinite(e.rate))
                throw InputError("spn_lpbnd: mode " + std::to_string(m + 1) + " of node " +
                                 std::to_string(ind) + " has mean firing rate " +
                                 std::to_string(e.rate) + "; a bound needs a finite positive one");
            md.push_back(e);
        }
    }
    if (md.empty()) throw InputError("spn_lpbnd: the net has no firing mode");
    const std::size_t E = md.size();

    std::vector<double> mu(E, 0.0);
    std::vector<std::vector<double>> net(E, std::vector<double>(L, 0.0));
    for (std::size_t e = 0; e < E; ++e) {
        mu[e] = md[e].rate;
        for (std::size_t l = 0; l < L; ++l) net[e][l] = md[e].fire[l] - md[e].enab[l];
    }

    // Per-level a priori bounds and the conserved sums, both off the same
    // minimal-support P-invariant basis. The reference writes its "cycle
    // population" family for UNWEIGHTED cycles; spn_sinvariants returns the
    // weighted invariants S m = V, which are equally linear and strictly
    // tighter, so those are what is emitted.
    const SpnInvariants inv = spn_sinvariants(sn, options.init);
    const std::vector<std::vector<long long>>& S = inv.S;
    const std::vector<long long>& V = inv.V;
    const std::vector<long long>& m0 = inv.m0;
    std::vector<double> B(L, inf);
    for (std::size_t i = 0; i < S.size(); ++i)
        for (std::size_t l = 0; l < L; ++l)
            if (S[i][l] > 0)
                B[l] = std::min(B[l], std::floor(static_cast<double>(V[i]) /
                                                 static_cast<double>(S[i][l])));

    // ---- variable layout
    const std::size_t ix = 0;
    const std::size_t iq = L;
    const std::size_t ith = L + E;
    const std::size_t iu = L + 2 * E;
    const std::size_t iy = L + 3 * E;
    const std::size_t nv = options.markovian ? L + 3 * E + L * E : L + 3 * E;

    lp::LpModel<double> lpm(nv);
    for (std::size_t e = 0; e < E; ++e) lpm.set_bounds(iq + e, 0.0, 1.0);
    for (std::size_t l = 0; l < L; ++l) {
        if (std::isfinite(B[l])) {
            lpm.set_bounds(ix + l, 0.0, B[l]);
            if (options.markovian)
                for (std::size_t e = 0; e < E; ++e) lpm.set_bounds(iy + l * E + e, 0.0, B[l]);
        }
    }

    const double tol = options.tol;

    // ---- (1) throughput: th_e = mu_e q_e
    for (std::size_t e = 0; e < E; ++e) {
        lpm.row_clear();
        lpm.row_add(ith + e, 1.0);
        lpm.row_add(iq + e, -mu[e]);
        lpm.emit_eq(0.0);
    }

    // ---- (2) flow balance: tokens are created at a level at the rate they are
    // consumed there. Holds for any stable net, Markovian or not.
    for (std::size_t l = 0; l < L; ++l) {
        lpm.row_clear();
        for (std::size_t e = 0; e < E; ++e) lpm.row_add(iq + e, mu[e] * net[e][l]);
        lpm.emit_eq(0.0);
    }

    // ---- (3)+(4) second moment and population covariance, from the
    // stationarity of E[X_l1 X_l2] under the uniformized chain. Table 1 writes
    // the q side as four sums over set intersections; since the memberships are
    // exactly "sigma > 0" and "pi > 0", those collapse to
    //     -(sigma_1 - pi_1)(sigma_2 - pi_2) = -net_1 net_2
    // per mode, which also makes the l1 == l2 case reduce to the second-moment
    // family with no separate derivation.
    if (options.markovian) {
        for (std::size_t l1 = 0; l1 < L; ++l1)
            for (std::size_t l2 = l1; l2 < L; ++l2) {
                lpm.row_clear();
                for (std::size_t e = 0; e < E; ++e) {
                    // at l1 == l2 the two y terms address the same column and
                    // row_add SUMS them, which is the factor of two the
                    // reference's (6) carries
                    lpm.row_add(iy + l1 * E + e, mu[e] * net[e][l2]);
                    lpm.row_add(iy + l2 * E + e, mu[e] * net[e][l1]);
                    lpm.row_add(iq + e, mu[e] * net[e][l1] * net[e][l2]);
                }
                lpm.emit_eq(0.0);
            }
    }

    // ---- (5) liveness, only when the caller vouches for it
    if (options.assumelive) {
        lpm.row_clear();
        for (std::size_t e = 0; e < E; ++e) lpm.row_add(iq + e, 1.0);
        lpm.emit_ge(1.0 - tol);
        if (options.markovian)
            for (std::size_t l = 0; l < L; ++l) {
                lpm.row_clear();
                lpm.row_add(ix + l, 1.0);
                for (std::size_t e = 0; e < E; ++e) lpm.row_add(iy + l * E + e, -1.0);
                lpm.emit_le(tol);
            }
    }

    // ---- (6) conflicting transitions: a mode that consumes no more and is
    // inhibited no sooner is enabled whenever the other is
    for (std::size_t e1 = 0; e1 < E; ++e1)
        for (std::size_t e2 = 0; e2 < E; ++e2) {
            if (e1 == e2) continue;
            bool dominated = true;
            for (std::size_t l = 0; l < L && dominated; ++l)
                if (!(md[e1].enab[l] <= md[e2].enab[l] && md[e1].inhib[l] >= md[e2].inhib[l]))
                    dominated = false;
            if (dominated) {
                lpm.row_clear();
                lpm.row_add(iq + e1, 1.0);
                lpm.row_add(iq + e2, -1.0);
                lpm.emit_ge(-tol);
            }
        }

    // ---- (7) boundedness, per level; and (8) cycle population, as the
    // weighted invariant equalities and their y companions
    if (options.markovian)
        for (std::size_t l = 0; l < L; ++l) {
            if (!std::isfinite(B[l])) continue;
            for (std::size_t e = 0; e < E; ++e) {
                lpm.row_clear();
                lpm.row_add(iy + l * E + e, 1.0);
                lpm.row_add(iq + e, -B[l]);
                lpm.emit_le(tol);
                lpm.row_clear();
                lpm.row_add(ix + l, 1.0);
                lpm.row_add(iy + l * E + e, -1.0);
                lpm.row_add(iq + e, B[l]);
                lpm.emit_le(B[l] + tol);
                if (B[l] > 0) {
                    lpm.row_clear();
                    lpm.row_add(ix + l, 1.0 - 1.0 / B[l]);
                    lpm.row_add(iy + l * E + e, -1.0);
                    lpm.row_add(iq + e, 1.0);
                    lpm.emit_ge(-tol);
                }
            }
        }
    for (std::size_t i = 0; i < S.size(); ++i) {
        lpm.row_clear();
        for (std::size_t l = 0; l < L; ++l) lpm.row_add(ix + l, static_cast<double>(S[i][l]));
        lpm.emit_eq(static_cast<double>(V[i]));
        if (options.markovian)
            for (std::size_t e = 0; e < E; ++e) {
                lpm.row_clear();
                for (std::size_t l = 0; l < L; ++l)
                    lpm.row_add(iy + l * E + e, static_cast<double>(S[i][l]));
                lpm.row_add(iq + e, -static_cast<double>(V[i]));
                lpm.emit_eq(0.0);
            }
    }

    // ---- (9) reachable marking: the mean lies in the state-equation cone
    for (std::size_t l = 0; l < L; ++l) {
        lpm.row_clear();
        lpm.row_add(ix + l, 1.0);
        for (std::size_t e = 0; e < E; ++e) lpm.row_add(iu + e, -net[e][l]);
        lpm.emit_eq(static_cast<double>(m0[l]));
    }

    // ---- (10) sample-path comparisons
    if (options.markovian) {
        double mutot = 0;
        for (std::size_t e = 0; e < E; ++e) mutot += mu[e];
        for (std::size_t l = 0; l < L; ++l) {
            for (std::size_t e = 0; e < E; ++e) {
                lpm.row_clear();
                lpm.row_add(iy + l * E + e, 1.0);
                lpm.row_add(ix + l, -1.0);
                lpm.emit_le(tol);
                if (md[e].enab[l] > 0) {
                    lpm.row_clear();
                    lpm.row_add(iy + l * E + e, 1.0);
                    lpm.row_add(iq + e, -md[e].enab[l]);
                    lpm.emit_ge(-tol);
                }
                if (std::isfinite(md[e].inhib[l])) {
                    lpm.row_clear();
                    lpm.row_add(iy + l * E + e, 1.0);
                    lpm.row_add(iq + e, -(md[e].inhib[l] - 1.0));
                    lpm.emit_le(tol);
                }
            }
            lpm.row_clear();
            lpm.row_add(ix + l, mutot);
            for (std::size_t e = 0; e < E; ++e) lpm.row_add(iy + l * E + e, -mu[e]);
            lpm.emit_ge(-tol);
        }
        for (std::size_t e = 0; e < E; ++e) {
            std::size_t ent = 0;
            std::size_t nent = 0;
            bool inhibited = false;
            for (std::size_t l = 0; l < L; ++l) {
                if (md[e].enab[l] > 0) {
                    ent = l;
                    ++nent;
                }
                if (std::isfinite(md[e].inhib[l])) inhibited = true;
            }
            if (nent == 1 && !inhibited) {
                lpm.row_clear();
                lpm.row_add(ix + ent, 1.0);
                lpm.row_add(iy + ent * E + e, -1.0);
                lpm.emit_le(md[e].enab[ent] - 1.0 + tol);
            }
        }
    }

    // ---- (11) enabling bounds, from Chernoff's inequality on the marking
    for (std::size_t e = 0; e < E; ++e) {
        std::vector<std::size_t> ent, inh;
        for (std::size_t l = 0; l < L; ++l) {
            if (md[e].enab[l] > 0) ent.push_back(l);
            if (std::isfinite(md[e].inhib[l])) inh.push_back(l);
        }
        const std::size_t d = ent.size() + inh.size();
        if (d == 0) continue;
        bool entBounded = !ent.empty();
        for (std::size_t j = 0; j < ent.size(); ++j)
            if (!std::isfinite(B[ent[j]])) entBounded = false;
        if (entBounded) {
            lpm.row_clear();
            lpm.row_add(iq + e, 1.0);
            double rhs = 1.0;
            bool ok = true;
            for (std::size_t j = 0; j < ent.size() && ok; ++j) {
                const std::size_t l = ent[j];
                const double den = B[l] - md[e].enab[l] + 1.0;
                if (den <= 0) {
                    ok = false;
                    break;
                }
                lpm.row_add(ix + l, -1.0 / den);
                rhs -= B[l] / den;
            }
            if (ok) {
                for (std::size_t j = 0; j < inh.size(); ++j)
                    lpm.row_add(ix + inh[j], 1.0 / md[e].inhib[inh[j]]);
                lpm.emit_ge(rhs - tol);
            } else {
                lpm.row_clear();
            }
        }
        // Upper side. Every term of the sum over input levels carries a "min"
        // operator, and Table 1's convention is that either operand may be
        // taken; each choice is a valid row and the whole set is the tightest
        // linear relaxation, so all of them are emitted while the count stays
        // small.
        bool ok = true;
        for (std::size_t j = 0; j < inh.size(); ++j) {
            const std::size_t l = inh[j];
            if (!std::isfinite(B[l]) || B[l] - md[e].inhib[l] + 1.0 <= 0) ok = false;
        }
        if (!ok) continue;
        const std::size_t nc = ent.size();
        std::vector<unsigned long> combos;
        if (nc <= 4) {
            for (unsigned long c = 0; c < (1UL << nc); ++c) combos.push_back(c);
        } else {
            combos.push_back(0);
            combos.push_back((1UL << nc) - 1);
        }
        for (std::size_t ci = 0; ci < combos.size(); ++ci) {
            const unsigned long c = combos[ci];
            lpm.row_clear();
            lpm.row_add(iq + e, static_cast<double>(d));
            double rhs = 0.0;
            for (std::size_t j = 0; j < inh.size(); ++j) {
                const std::size_t l = inh[j];
                const double den = B[l] - md[e].inhib[l] + 1.0;
                rhs += B[l] / den;
                lpm.row_add(ix + l, 1.0 / den);
            }
            for (std::size_t j = 0; j < nc; ++j) {
                const std::size_t l = ent[j];
                if (((c >> j) & 1UL) == 0) {
                    lpm.row_add(ix + l, -1.0 / md[e].enab[l]);
                } else {
                    rhs += 1.0;
                }
            }
            lpm.emit_le(rhs + tol);
        }
    }

    // ---- (12) Little's law at each level: the mean sojourn time of a token is
    // at least the mean minimum firing time of the modes that can remove it
    if (options.markovian)
        for (std::size_t l = 0; l < L; ++l) {
            double out = 0;
            for (std::size_t e = 0; e < E; ++e)
                if (md[e].enab[l] > 0) out += mu[e];
            if (out <= 0) continue;
            lpm.row_clear();
            lpm.row_add(ix + l, out);
            for (std::size_t e = 0; e < E; ++e) lpm.row_add(iq + e, -mu[e] * md[e].fire[l]);
            lpm.emit_ge(-tol);
        }

    SpnLpBounds out;
    out.tokens_lo.assign(L, 0.0);
    out.tokens_hi.assign(L, 0.0);
    out.place_tput_lo.assign(L, 0.0);
    out.place_tput_hi.assign(L, 0.0);
    std::vector<double> c(nv, 0.0);
    for (std::size_t l = 0; l < L; ++l) {
        std::fill(c.begin(), c.end(), 0.0);
        c[ix + l] = 1.0;
        out.tokens_lo[l] = detail::spn_lp_opt(lpm, c, true);
        out.tokens_hi[l] = detail::spn_lp_opt(lpm, c, false);
        std::fill(c.begin(), c.end(), 0.0);
        bool any = false;
        for (std::size_t e = 0; e < E; ++e)
            if (md[e].enab[l] > 0) {
                c[iq + e] = mu[e] * md[e].enab[l];
                any = true;
            }
        out.place_tput_lo[l] = any ? detail::spn_lp_opt(lpm, c, true) : 0.0;
        out.place_tput_hi[l] = any ? detail::spn_lp_opt(lpm, c, false) : 0.0;
    }
    out.mode_tput_lo.assign(E, 0.0);
    out.mode_tput_hi.assign(E, 0.0);
    out.mode_util_lo.assign(E, 0.0);
    out.mode_util_hi.assign(E, 0.0);
    for (std::size_t e = 0; e < E; ++e) {
        std::fill(c.begin(), c.end(), 0.0);
        c[ith + e] = 1.0;
        out.mode_tput_lo[e] = detail::spn_lp_opt(lpm, c, true);
        out.mode_tput_hi[e] = detail::spn_lp_opt(lpm, c, false);
        std::fill(c.begin(), c.end(), 0.0);
        c[iq + e] = 1.0;
        out.mode_util_lo[e] = detail::spn_lp_opt(lpm, c, true);
        out.mode_util_hi[e] = detail::spn_lp_opt(lpm, c, false);
    }

    out.places = places;
    out.levelname.resize(L);
    for (std::size_t pp = 0; pp < P; ++pp) out.levelname[pp] = sn.nodes[places[pp] - 1].name;
    out.modes = md;
    out.bound = B;
    out.nplacelevels = L;
    out.markovian = options.markovian;
    out.nvars = nv;
    out.nrows = lpm.num_rows();
    return out;
}

}  // namespace spn
}  // namespace line
