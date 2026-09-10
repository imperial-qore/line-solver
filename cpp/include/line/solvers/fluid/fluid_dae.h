/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_DAE_H
#define LINE_SOLVERS_FLUID_FLUID_DAE_H

/**
 * The min-normal closure as a DIFFERENTIAL-ALGEBRAIC system: `solver_fluid_dae.m`.
 *
 * SOLVER_FLUID_MOMENTS already solves a differential system (the mean) coupled
 * to an algebraic one (the covariance). It solves them by SUCCESSIVE
 * SUBSTITUTION: integrate the mean to its fixed point at a held variance, solve
 * the Lyapunov equation there, extract sigma2, repeat. This states the same
 * closure as one system and solves it as one system. The closure itself is
 * unchanged -- the drift, the rate factors and the Lyapunov equation are taken
 * from FLUID_MOMENT_TERMS and FLUID_LYAPUNOV untouched -- only the way the
 * coupled equations are discharged.
 *
 * STEADY STATE is an algebraic system, not an integration:
 *
 *     0 = D r(x, sigma2)               drift residual, nstate rows
 *     0 = C x - Nchain                 population conservation, one row per
 *                                      closed chain
 *     0 = sigma2 - sigmaOf(x, sigma2)  closure consistency, one row per
 *                                      closable station
 *
 * solved simultaneously by a damped projected Newton. Integrating a stable ODE
 * until it stops moving is a poor way to solve f(x)=0: the cost is set by the
 * slowest mode of the model rather than by the accuracy wanted, which is why
 * the stiff models are expensive on the substitution route.
 *
 * TRANSIENT is an index-1 DAE with a SINGULAR MASS MATRIX, integrated by the
 * vendored RODAS (third_party/rodas.hpp):
 *
 *     d/dt x = D r(x, sigma2)    differential rows
 *     0      = C x - Nchain      algebraic rows
 *
 * RODAS IS THE DEFAULT AND THE ONLY INTEGRATOR HERE. The rest of the fluid
 * solver runs LSODA (line/util/lsoda.h) so that MATLAB, the JAR, Python and C++
 * agree in the last digits, but LSODA integrates y' = f and cannot carry a
 * singular mass matrix at all, so it is not a candidate for this path. RODAS is
 * a Rosenbrock method of order (3)4 for M y' = f with singular M, and being a
 * fixed sequence of six linear solves rather than an iteration it has no
 * convergence history for a future port to diverge on -- which is what makes it
 * the right choice for a route that the other three codebases do not have yet.
 *
 * WHY THE CONSTRAINT IS WRITTEN DIFFERENTLY IN THE TWO MODES. At a fixed point
 * every flow already balances, so the differentiated form d/dt(Cx) = 0 is
 * satisfied by anything and pins nothing; the steady state therefore uses
 * `C x = N` directly. The transient needs the opposite: the mass matrix zeroes
 * one row per chain and the constraint residual is written there, which is the
 * index-1 form RODAS integrates.
 *
 * WHAT THE ALGEBRAIC CONSTRAINT BUYS. Population conservation otherwise holds
 * only to integrator tolerance: it is a consequence of the drift (the rows of D
 * sum to zero on a closed chain), never an equation. Writing it as a constraint
 * also makes the Newton system solvable, because the drift Jacobian is singular
 * along exactly the conserved directions -- the same singularity FLUID_LYAPUNOV
 * works around by projecting onto range(D) -- so the constraint rows supply the
 * missing rank instead of a pseudo-inverse hiding it.
 *
 * WHY THE COVARIANCE IS NOT A NEWTON UNKNOWN. Sigma is nstate^2 entries, so a
 * Jacobian over it is quartic work -- strictly worse than the cubic Lyapunov
 * solves it would replace. Sigma is LINEAR in itself for a held x, so it is
 * eliminated by one Lyapunov solve per residual evaluation and only sigma2, M
 * numbers, joins x in the unknown vector.
 *
 * @see fluid_moments.h - the closure this solves, and the substitution route
 * @see third_party/rodas.hpp - the DAE integrator
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/solvers/fluid/fluid_moments.h"
#include "line/solvers/fluid/fluid_nonhyperbolic.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"
#include "line/util/lstsq.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "rodas.hpp"

namespace line {
namespace fluid {

/**
 * Finite capacity regions as linear admission constraints on the fluid state.
 *
 * EVERY FORM THE REGION CAN TAKE IS ONE FAMILY. A region carries a global job
 * cap, per-class caps, a memory budget with per-class sizes, and an optional
 * linear pair (A,b); all four are the same object once written against the
 * state, `Arow x <= b` with Arow the per-class weight on the region's
 * coordinates. Collapsing them here means the solver carries one mechanism
 * rather than four.
 *
 * ONLY WAITQ IS A CONSTRAINT ON THIS DRIFT. Under a waiting queue a blocked job
 * waits outside the region and is admitted later, so the population is
 * conserved and only the admission FLOW is throttled. DROP destroys the job,
 * BAS/BBS/RSRD hold a server upstream, and the retrial rules move it to an
 * orbit: each changes the event set itself, so each needs a different drift
 * rather than a constraint on this one.
 */
struct FluidDaeConstraints {
    Matrix<double> A;                        ///< (ncon x nstate)
    Matrix<double> As;                       ///< (ncon x nstaging), see fluid_dae_extend
    std::vector<double> b;
    std::vector<std::size_t> region;         ///< which region, npos for a station cap
    std::vector<std::size_t> station;        ///< which station, npos for a region cap
    std::vector<std::size_t> klass_row;      ///< which class, npos when several
    std::vector<bool> staged;                ///< true where the job waits in a room
    std::vector<std::string> label;
    std::vector<std::vector<bool> > member;  ///< (nregions x nstate)
    std::vector<std::size_t> coord_class;
    std::size_t nregions = 0;
    bool empty() const { return b.empty(); }
    static std::size_t none() { return static_cast<std::size_t>(-1); }
};

/**
 * The largest `row x` the population can produce, ignoring the coupling.
 *
 * An upper bound is what is wanted: too loose only costs a constraint that stays
 * inactive, while too tight would discard a cap that does bind. THE BOUND IS PER
 * CLASS, and it has to be -- the heaviest weight times the whole population never
 * prunes a PER-CLASS cap set to its own class population, exactly the row the
 * struct refresh derives at every station of every closed model.
 */
inline double fluid_dae_reach(const std::vector<double>& row,
                              const std::vector<std::size_t>& coord_class,
                              const std::vector<double>& njobs, std::size_t K) {
    double total = 0.0;
    for (std::size_t k = 0; k < K; ++k) {
        double w = 0.0;
        bool any = false;
        for (std::size_t s = 0; s < row.size(); ++s)
            if (coord_class[s] == k) { w = std::max(w, row[s]); any = true; }
        if (!any || w <= 0.0) continue;
        const double nk = (k < njobs.size()) ? njobs[k] : std::numeric_limits<double>::infinity();
        if (!std::isfinite(nk)) return std::numeric_limits<double>::infinity();
        total += w * nk;
    }
    return total;
}

template <typename T>
inline FluidDaeConstraints fluid_dae_constraints(const qn::NetworkStruct<T>& sn,
                                                 const FluidMomentTerms& terms) {
    FluidDaeConstraints con;
    const std::size_t nstate = terms.nstate;
    const std::size_t M = terms.class_block.size();
    const std::size_t K = M ? terms.class_block[0].size() : 0;
    const std::size_t NONE = FluidDaeConstraints::none();
    con.nregions = sn.regions.size();
    con.coord_class.assign(nstate, 0);
    std::vector<std::size_t> coord_station(nstate, NONE);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t s : terms.class_block[i][k]) {
                con.coord_class[s] = k;
                coord_station[s] = i;
            }
    con.member.assign(std::max<std::size_t>(con.nregions, 1), std::vector<bool>(nstate, false));

    const double UNB = -1.0;
    std::vector<std::vector<double> > rows;
    std::vector<std::size_t> rgs, sts, kls;
    std::vector<bool> sdg;
    for (std::size_t f = 0; f < con.nregions; ++f) {
        const typename qn::NetworkStruct<T>::Region& R = sn.regions[f];
        std::vector<bool> inR(nstate, false);
        bool any_member = false;
        for (std::size_t s = 0; s < nstate; ++s) {
            const std::size_t i = coord_station[s];
            if (i != NONE && i < R.members.size() && R.members[i]) { inR[s] = true; any_member = true; }
        }
        if (!any_member) continue;
        con.member[f] = inR;

        for (std::size_t k = 0; k < K && k < R.rule.size(); ++k)
            if (R.rule[k] != lang::DropStrategy::WAITQ)
                throw UnsupportedError(
                    "fluid_dae_constraints: region " + std::to_string(f + 1) +
                    " applies a drop rule other than a waiting queue to class " +
                    std::to_string(k + 1) +
                    ". Only a waiting queue conserves the population and throttles the admission "
                    "flow, which is what an algebraic equation on this drift can express.");
        for (std::size_t k = 0; k < K && k < R.weight.size(); ++k)
            if (std::fabs(static_cast<double>(R.weight[k]) - 1.0) > 1e-12)
                throw UnsupportedError(
                    "fluid_dae_constraints: region " + std::to_string(f + 1) +
                    " sets per-class admission weights, which decide WHICH blocked class enters "
                    "when capacity frees up. The throttle carries no such priority and would "
                    "ignore them silently.");

        std::size_t member_row = 0;
        for (std::size_t i = 0; i < R.members.size(); ++i)
            if (R.members[i]) { member_row = i; break; }

        // 1. region-global job cap: column K of `cap`
        if (member_row < R.cap.size() && R.cap[member_row].size() > K) {
            const double g = static_cast<double>(R.cap[member_row][K]);
            if (std::isfinite(g) && g != UNB && g >= 0.0) {
                std::vector<double> row(nstate, 0.0);
                for (std::size_t s = 0; s < nstate; ++s) if (inR[s]) row[s] = 1.0;
                rows.push_back(row); con.b.push_back(g);
                rgs.push_back(f); sts.push_back(NONE); kls.push_back(NONE); sdg.push_back(true);
                con.label.push_back("region " + std::to_string(f + 1) + " global job cap");
            }
        }
        // 2. per-class job caps
        if (member_row < R.cap.size())
            for (std::size_t k = 0; k < K && k < R.cap[member_row].size(); ++k) {
                const double c = static_cast<double>(R.cap[member_row][k]);
                if (!std::isfinite(c) || c == UNB || c < 0.0) continue;
                std::vector<double> row(nstate, 0.0);
                bool any = false;
                for (std::size_t s = 0; s < nstate; ++s)
                    if (inR[s] && con.coord_class[s] == k) { row[s] = 1.0; any = true; }
                if (!any) continue;
                rows.push_back(row); con.b.push_back(c);
                rgs.push_back(f); sts.push_back(NONE); kls.push_back(k); sdg.push_back(true);
                con.label.push_back("region " + std::to_string(f + 1) + " class " +
                                    std::to_string(k + 1) + " job cap");
            }
        // 3. region-global memory budget, weighted by each class's size
        if (member_row < R.maxmem.size()) {
            const double mem = static_cast<double>(R.maxmem[member_row]);
            if (std::isfinite(mem) && mem != UNB && mem >= 0.0) {
                std::vector<double> row(nstate, 0.0);
                bool any = false;
                for (std::size_t s = 0; s < nstate; ++s) {
                    if (!inR[s]) continue;
                    const std::size_t k = con.coord_class[s];
                    const double w = (k < R.size.size()) ? static_cast<double>(R.size[k]) : 1.0;
                    if (w != 0.0) { row[s] = w; any = true; }
                }
                if (any) {
                    rows.push_back(row); con.b.push_back(mem);
                    rgs.push_back(f); sts.push_back(NONE); kls.push_back(NONE); sdg.push_back(true);
                    con.label.push_back("region " + std::to_string(f + 1) + " memory budget");
                }
            }
        }
        // 4. explicit linear constraints
        for (std::size_t c = 0; c < R.lincon_A.rows(); ++c) {
            std::vector<double> row(nstate, 0.0);
            bool any = false;
            for (std::size_t k = 0; k < K && k < R.lincon_A.cols(); ++k) {
                const double a = static_cast<double>(R.lincon_A(c, k));
                if (a == 0.0) continue;
                for (std::size_t s = 0; s < nstate; ++s)
                    if (inR[s] && con.coord_class[s] == k) { row[s] = a; any = true; }
            }
            if (!any) continue;
            rows.push_back(row);
            con.b.push_back(static_cast<double>(R.lincon_b[c]));
            rgs.push_back(f); sts.push_back(NONE); kls.push_back(NONE); sdg.push_back(true);
            con.label.push_back("region " + std::to_string(f + 1) + " linear constraint " +
                                std::to_string(c + 1));
        }
    }

    // 5. per-station buffers. A STATION CAP IS THE ONE-STATION CASE OF THE SAME
    // ROW, and every fluid method other than this one ignores it outright: nothing
    // in the fluid tree reads sn.cap or sn.classcap, so a capped station was
    // integrated as an unbounded one and the table reported more jobs in the buffer
    // than the buffer holds.
    for (std::size_t i = 0; i < M; ++i) {
        if (terms.is_ext[i]) continue;   // a source holds no jobs
        bool any_at = false;
        for (std::size_t s = 0; s < nstate; ++s) any_at |= coord_station[s] == i;
        if (!any_at) continue;
        if (i < sn.cap.size()) {
            const double g = static_cast<double>(sn.cap[i]);
            if (std::isfinite(g) && g >= 0.0) {
                std::vector<double> row(nstate, 0.0);
                for (std::size_t s = 0; s < nstate; ++s) if (coord_station[s] == i) row[s] = 1.0;
                rows.push_back(row); con.b.push_back(g);
                rgs.push_back(NONE); sts.push_back(i); kls.push_back(NONE); sdg.push_back(false);
                con.label.push_back("station " + std::to_string(i + 1) + " buffer");
            }
        }
        if (i < sn.classcap.size())
            for (std::size_t k = 0; k < K && k < sn.classcap[i].size(); ++k) {
                const double c = static_cast<double>(sn.classcap[i][k]);
                if (!std::isfinite(c) || c < 0.0) continue;
                std::vector<double> row(nstate, 0.0);
                bool any = false;
                for (std::size_t s = 0; s < nstate; ++s)
                    if (coord_station[s] == i && con.coord_class[s] == k) { row[s] = 1.0; any = true; }
                if (!any) continue;
                rows.push_back(row); con.b.push_back(c);
                rgs.push_back(NONE); sts.push_back(i); kls.push_back(k); sdg.push_back(false);
                con.label.push_back("station " + std::to_string(i + 1) + " class " +
                                    std::to_string(k + 1) + " buffer");
            }
    }

    const std::size_t ncand = rows.size();
    std::vector<bool> keep(ncand, true);
    const std::vector<double> njobs = sn.njobs();

    // A CAP THE POPULATION CANNOT REACH IS NOT A CAP, and dropping it here keeps it
    // out of the active-set loop, where it would be tested on every pass and never
    // bind while its multiplier stayed an unknown with no equation to pin it.
    for (std::size_t c = 0; c < ncand; ++c)
        if (fluid_dae_reach(rows[c], con.coord_class, njobs, K) <= con.b[c] + 1e-14)
            keep[c] = false;

    // A CAP ON COORDINATES THE IMMEDIATE REDUCTION FOLDED AWAY IS VACUOUS, NOT
    // MALFORMED. `hide_immediate` stochastic-complements an Immediate-rate
    // coordinate out of the event set -- the MMT transform's zero-service Join is
    // one, until the fork-join fixed point gives it a synchronisation delay -- and
    // the reduced drift then holds no mass there and has no event landing on it.
    // Left in, such a row reaches `fluid_dae_gates` with no gating event and is
    // reported as a limit the model cannot approach, which is the right message for
    // a station that really does hold jobs and the wrong one here: this station
    // holds none, so its buffer is satisfied identically.
    for (std::size_t c = 0; c < ncand; ++c) {
        if (!keep[c]) continue;
        bool live = false;
        for (std::size_t s = 0; s < nstate && !live; ++s)
            if (rows[c][s] != 0.0 && !fluid_coord_eliminated(terms, s)) live = true;
        if (!live) keep[c] = false;
    }

    // THE SAME ROW TWICE IS A SINGULAR NEWTON SYSTEM, not a redundancy the least
    // squares absorbs: two identical rows both bind, each takes a multiplier, and
    // nothing distinguishes them. The struct refresh derives classcap from cap, so a
    // single-class model declares the station total and the class buffer as the same
    // row. The TIGHTER bound survives; on a tie the region row does.
    for (std::size_t c = 0; c < ncand; ++c) {
        if (!keep[c]) continue;
        for (std::size_t d = c + 1; d < ncand; ++d) {
            if (!keep[d]) continue;
            bool same = true;
            for (std::size_t s = 0; s < nstate && same; ++s)
                same = std::fabs(rows[c][s] - rows[d][s]) <= 1e-14;
            if (!same) continue;
            const bool takeover = con.b[d] < con.b[c] - 1e-14 ||
                (std::fabs(con.b[d] - con.b[c]) <= 1e-14 && sdg[d] && !sdg[c]);
            if (takeover) {
                con.b[c] = con.b[d]; rgs[c] = rgs[d]; sts[c] = sts[d]; kls[c] = kls[d];
                sdg[c] = sdg[d]; con.label[c] = con.label[d];
            }
            keep[d] = false;
        }
    }

    // A ROW ITS OWN PER-CLASS ROWS ALREADY IMPLY IS RANK, NOT INFORMATION: a
    // two-class station capped 3 and 3 also declares a total of 6, exactly the sum
    // of the two class rows, and all three then bind together with rank 2. Only an
    // exact implication is dropped, so a total TIGHTER than the sum of its parts
    // survives.
    for (std::size_t c = 0; c < ncand; ++c) {
        if (!keep[c] || kls[c] != NONE) continue;
        std::vector<std::size_t> classes;
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t s = 0; s < nstate; ++s)
                if (con.coord_class[s] == k && coord_station[s] != NONE && rows[c][s] > 0.0) {
                    classes.push_back(k);
                    break;
                }
        if (classes.empty()) continue;
        bool implied = true;
        double budget = 0.0;
        for (std::size_t ci = 0; ci < classes.size() && implied; ++ci) {
            const std::size_t k = classes[ci];
            std::size_t part = NONE;
            for (std::size_t d = 0; d < ncand; ++d) {
                if (!keep[d] || d == c || kls[d] != k || rgs[d] != rgs[c] || sts[d] != sts[c])
                    continue;
                bool covers = true;
                for (std::size_t s = 0; s < nstate && covers; ++s)
                    if (con.coord_class[s] == k && coord_station[s] != NONE &&
                        rows[c][s] > 0.0 && rows[d][s] < rows[c][s] - 1e-14)
                        covers = false;
                if (covers) { part = d; break; }
            }
            if (part == NONE) implied = false;
            else budget += con.b[part];
        }
        if (implied && budget <= con.b[c] + 1e-14) keep[c] = false;
    }

    std::vector<std::vector<double> > kept_rows;
    std::vector<double> kept_b;
    std::vector<std::string> kept_lab;
    for (std::size_t c = 0; c < ncand; ++c) {
        if (!keep[c]) continue;
        kept_rows.push_back(rows[c]);
        kept_b.push_back(con.b[c]);
        kept_lab.push_back(con.label[c]);
        con.region.push_back(rgs[c]);
        con.station.push_back(sts[c]);
        con.klass_row.push_back(kls[c]);
        con.staged.push_back(sdg[c]);
    }
    con.b = kept_b;
    con.label = kept_lab;
    con.A = Matrix<double>(kept_rows.size(), nstate, 0.0);
    for (std::size_t r = 0; r < kept_rows.size(); ++r)
        for (std::size_t s = 0; s < nstate; ++s) con.A(r, s) = kept_rows[r][s];
    con.As = Matrix<double>(kept_rows.size(), 0, 0.0);

    // WHICH STATION RULES ARE A CONSTRAINT ON THIS DRIFT, and which are a different
    // event set. THE RULE IS NOT WHAT DECIDES THE SEMANTICS -- the class type is,
    // exactly as State.arrivalIsLost decides it for every other solver: a closed job
    // is never lost (the upstream departure is disabled instead) and an open one is
    // never held. So a waiting queue and a drop are BOTH constraints here; what is
    // refused is the rules that add STATE. A rule is only a contradiction where the
    // cap can BIND, so this runs on the surviving rows.
    for (std::size_t c = 0; c < con.b.size(); ++c) {
        const std::size_t i = con.station[c];
        if (i == NONE || i >= sn.droprule.size()) continue;
        for (std::size_t k = 0; k < K && k < sn.droprule[i].size(); ++k) {
            bool weighs = false;
            for (std::size_t s = 0; s < nstate; ++s)
                if (coord_station[s] == i && con.coord_class[s] == k && con.A(c, s) > 0.0)
                    weighs = true;
            if (!weighs) continue;
            const lang::DropStrategy rule = sn.droprule[i][k];
            if (rule != lang::DropStrategy::WAITQ && rule != lang::DropStrategy::DROP)
                throw UnsupportedError(
                    "fluid_dae_constraints: station " + std::to_string(i + 1) +
                    " applies a blocking or retrial rule to class " + std::to_string(k + 1) +
                    ", and its buffer binds. Only a waiting queue or a drop is a constraint on "
                    "this drift: BAS/BBS/RSRD add a blocked-server state to the upstream station "
                    "and the retrial rules add an orbit, so each needs a different drift rather "
                    "than an algebraic equation on this one.");
        }
    }
    return con;
}

/**
 * Which events each cap throttles, and what happens to the mass it stops.
 *
 * THE THREE ANSWERS ARE THE MODEL, AND THEY ARE NOT INTERCHANGEABLE:
 *   STAGED  a finite capacity region under a waiting queue. The job COMPLETES
 *           upstream service and waits outside the region, which is what JMT and
 *           LDES simulate.
 *   HELD    a station buffer reached by a CLOSED class. LINE refuses to lose a
 *           closed job and disables the upstream departure instead, so the job is
 *           still at the upstream station and the whole event is scaled.
 *   LOSS    a station buffer reached by an OPEN class. The arrival fires and only
 *           the carried flow is admitted, so the ARRIVAL leg alone is scaled.
 */
struct FluidDaeGates {
    std::vector<std::vector<bool> > gate;   ///< (ncon x nevents)
    std::vector<std::vector<bool> > held;
    std::vector<std::vector<bool> > loss;
    Matrix<double> Dn, DnExt, Dp;           ///< the jump matrix split in three
    bool has = false;
};

template <typename T>
inline FluidDaeGates fluid_dae_gates(const qn::NetworkStruct<T>& sn, const FluidMomentTerms& t,
                                     const FluidDaeConstraints& con) {
    FluidDaeGates g;
    const std::size_t ncon = con.b.size(), nev = t.D.cols(), nstate = t.nstate;
    g.gate.assign(ncon, std::vector<bool>(nev, false));
    g.held.assign(ncon, std::vector<bool>(nev, false));
    g.loss.assign(ncon, std::vector<bool>(nev, false));
    if (ncon == 0) return g;
    g.has = true;
    g.Dn = Matrix<double>(nstate, nev, 0.0);
    g.DnExt = Matrix<double>(nstate, nev, 0.0);
    g.Dp = Matrix<double>(nstate, nev, 0.0);
    std::vector<bool> ext_coord(nstate, false);
    for (std::size_t i = 0; i < t.station_block.size(); ++i)
        if (t.is_ext[i])
            for (std::size_t s : t.station_block[i]) ext_coord[s] = true;
    for (std::size_t s = 0; s < nstate; ++s)
        for (std::size_t e = 0; e < nev; ++e) {
            const double d = t.D(s, e);
            // A LOST ARRIVAL IS RETURNED TO THE SOURCE POOL: the EXT coordinate is a
            // normalisation and not a population, so scaling only the arrival leg of
            // a lost event would unbalance its row by exactly the loss.
            if (d < 0.0) { if (ext_coord[s]) g.DnExt(s, e) = d; else g.Dn(s, e) = d; }
            if (d > 0.0) g.Dp(s, e) = d;
        }
    const double tol = 1e-7;
    const std::vector<double> njobs = sn.njobs();
    for (std::size_t c = 0; c < ncon; ++c) {
        bool any = false;
        for (std::size_t e = 0; e < nev; ++e) {
            double delta = 0.0;
            for (std::size_t s = 0; s < nstate; ++s) delta += con.A(c, s) * t.D(s, e);
            if (delta <= tol) continue;
            g.gate[c][e] = true;
            any = true;
            if (con.staged[c]) continue;
            // the class is read off the coordinate the mass LANDS on, inside the
            // capped station: a class switch on entry would otherwise ask the class
            // the job is leaving behind whether it may be lost
            std::size_t k = static_cast<std::size_t>(-1);
            for (std::size_t s = 0; s < nstate; ++s)
                if (g.Dp(s, e) > tol && con.A(c, s) > 0.0) { k = con.coord_class[s]; break; }
            const bool isopen = k != static_cast<std::size_t>(-1) && k < njobs.size() &&
                                !std::isfinite(njobs[k]);
            g.loss[c][e] = isopen;
            g.held[c][e] = !isopen;
        }
        if (!any)
            throw UnsupportedError(
                "fluid_dae_gates: no event increases " + con.label[c] +
                ", so the cap can never be approached and there is no admission flow for the "
                "constraint to throttle. This is a malformed limit rather than a solvable one.");
    }
    return g;
}


/**
 * The waiting room outside a capped region, as fluid coordinates.
 *
 * WHY THE CONSTRAINT ALONE IS NOT ENOUGH. Throttling a region's admission events
 * does hold its population at the cap, but it holds it by slowing the UPSTREAM
 * STATION'S COMPLETIONS -- an admission event IS that station finishing a job --
 * so blocked mass piles up at a station it has already finished being served by.
 * Where that station is a delay the error is visible as a broken Little's law.
 * A waiting queue means the job COMPLETES upstream service and then waits; it is
 * somewhere else, and the model needs somewhere else to put it.
 *
 * Each region gains one coordinate per class, and every admission splits in two:
 *   upstream -> staging   at the nominal rate, so the upstream station empties
 *                         exactly as it would with no region;
 *   staging  -> region    at theta_f * s_{f,c}, the throttled leg.
 * The two jumps sum to the original, so only where the mass rests changes.
 *
 * ONE THROTTLE PER REGION, hence one binding cap per region: a region has a
 * single admission control -- how fast its queue drains -- so it can satisfy
 * exactly one equality.
 */
struct FluidDaeStaging {
    std::size_t n = 0;
    std::vector<std::size_t> region, klass;      ///< per staging coordinate
    std::vector<bool> adm;                       ///< per event: an admission?
    std::vector<std::size_t> adm_region, adm_stage;
    std::vector<std::vector<bool> > gated_by;    ///< (ncon x n) which rows gate which room
    Matrix<double> Dn, Dp;                       ///< negative and positive parts of D
};

inline FluidDaeStaging fluid_dae_staging(const FluidMomentTerms& t,
                                         const FluidDaeConstraints& con) {
    FluidDaeStaging stg;
    const std::size_t nev = t.D.cols();
    const std::size_t ncon = con.b.size();
    stg.adm.assign(nev, false);
    stg.adm_region.assign(nev, 0);
    stg.adm_stage.assign(nev, 0);
    stg.gated_by.assign(ncon, std::vector<bool>());
    bool any_staged = false;
    for (std::size_t c = 0; c < ncon; ++c) any_staged = any_staged || con.staged[c];
    // A STATION BUFFER GETS NO ROOM, and that is not an omission: LINE disables the
    // upstream departure rather than moving the job out, so the blocked mass is
    // still at the upstream station and still counted there.
    if (con.empty() || con.nregions == 0 || !any_staged) return stg;

    const double tol = 1e-9;
    stg.Dn = Matrix<double>(t.nstate, nev, 0.0);
    stg.Dp = Matrix<double>(t.nstate, nev, 0.0);
    for (std::size_t s = 0; s < t.nstate; ++s)
        for (std::size_t e = 0; e < nev; ++e) {
            const double d = t.D(s, e);
            if (d < 0.0) stg.Dn(s, e) = d;
            if (d > 0.0) stg.Dp(s, e) = d;
        }

    std::vector<bool> staged_region(con.nregions, false);
    for (std::size_t c = 0; c < ncon; ++c)
        if (con.staged[c] && con.region[c] < con.nregions) staged_region[con.region[c]] = true;

    const std::size_t K = t.class_block.empty() ? 0 : t.class_block[0].size();
    std::vector<std::vector<std::size_t> > idx(con.nregions, std::vector<std::size_t>(K, 0));
    for (std::size_t f = 0; f < con.nregions; ++f) {
        if (f >= con.member.size() || !staged_region[f]) continue;
        for (std::size_t e = 0; e < nev; ++e) {
            // net change of this region's population: positive means the event
            // brings mass in from outside, which is what the queue feeds
            double delta = 0.0;
            for (std::size_t s = 0; s < t.nstate; ++s)
                if (con.member[f][s]) delta += t.D(s, e);
            if (delta <= tol) continue;
            // the class is read off the coordinate the mass LANDS on, inside the
            // region: a class switch on entry would otherwise stage the job under
            // the class it is leaving behind
            std::size_t k = K;
            for (std::size_t s = 0; s < t.nstate; ++s)
                if (con.member[f][s] && stg.Dp(s, e) > tol) { k = con.coord_class[s]; break; }
            if (k >= K) continue;
            if (idx[f][k] == 0) {
                stg.region.push_back(f);
                stg.klass.push_back(k);
                idx[f][k] = stg.region.size();  // 1-based, 0 means none
            }
            stg.adm[e] = true;
            stg.adm_region[e] = f;
            stg.adm_stage[e] = idx[f][k] - 1;
        }
    }
    stg.n = stg.region.size();

    // WHICH ROWS GATE WHICH ROOM. A region-global cap gates every room of its
    // region, a per-class cap only the room of its class. A room gated by several
    // ACTIVE rows drains at the harmonic composition of their rates, which is what
    // lets two caps of one region bind at once.
    stg.gated_by.assign(ncon, std::vector<bool>(stg.n, false));
    for (std::size_t c = 0; c < ncon; ++c) {
        if (!con.staged[c] || con.region[c] >= con.nregions) continue;
        for (std::size_t j = 0; j < stg.n; ++j) {
            if (stg.region[j] != con.region[c]) continue;
            double w = 0.0;
            for (std::size_t s = 0; s < t.nstate; ++s)
                if (con.member[con.region[c]][s] && con.coord_class[s] == stg.klass[j])
                    w = std::max(w, con.A(c, s));
            stg.gated_by[c][j] = w > 0.0;
        }
    }
    return stg;
}

/**
 * Extend every cap to the staging coordinates that hold mass INSIDE it.
 *
 * A waiting room is outside the region it feeds, which is the whole point of it --
 * but it is not outside every OTHER limit. Where two regions overlap, an admission
 * into the inner one is an INTERNAL move of the outer one: the job leaves a station
 * of the outer region, waits, and re-enters a station of the same outer region,
 * never having left it. Counting only the state coordinates would take that mass
 * out of the outer cap for as long as it waits, and the Newton system that results
 * is inconsistent rather than merely inexact.
 */
inline void fluid_dae_extend(FluidDaeConstraints& con, const FluidDaeStaging& stg,
                             const FluidMomentTerms& t) {
    const std::size_t ncon = con.b.size();
    con.As = Matrix<double>(ncon, stg.n, 0.0);
    if (ncon == 0 || stg.n == 0) return;
    const std::size_t nev = t.D.cols();
    const double tol = 1e-9;
    std::vector<std::vector<std::size_t> > feeds(stg.n);
    std::vector<std::size_t> dest(stg.n, static_cast<std::size_t>(-1));
    for (std::size_t e = 0; e < nev; ++e) {
        if (!stg.adm[e]) continue;
        const std::size_t j = stg.adm_stage[e];
        for (std::size_t s = 0; s < t.nstate; ++s) {
            if (t.D(s, e) < -tol &&
                std::find(feeds[j].begin(), feeds[j].end(), s) == feeds[j].end())
                feeds[j].push_back(s);
            if (dest[j] == static_cast<std::size_t>(-1) && t.D(s, e) > tol &&
                con.member[stg.region[j]][s])
                dest[j] = s;
        }
    }
    for (std::size_t c = 0; c < ncon; ++c)
        for (std::size_t j = 0; j < stg.n; ++j) {
            if (dest[j] == static_cast<std::size_t>(-1) || feeds[j].empty()) continue;
            const double w = con.A(c, dest[j]);
            if (w <= 0.0) continue;
            bool all_inside = true;
            for (std::size_t a = 0; a < feeds[j].size() && all_inside; ++a)
                all_inside = con.A(c, feeds[j][a]) > 0.0;
            if (all_inside) con.As(c, j) = w;
        }
}

/** The two legs of every event under the active caps, and the waiting rooms. */
struct FluidDaeLegs {
    std::vector<double> rup;    ///< the rate each event FIRES at
    std::vector<double> rin;    ///< the rate mass LANDS at
    std::vector<double> ds;     ///< the derivative of each room
    std::vector<double> drain;  ///< each room's total outflow
};

/**
 * Shared by the steady-state residual and the transient right-hand side so that
 * the two solve the SAME model and not two spellings of it. What differs is only
 * what a STAGED cap's multiplier means: a drain RATE for the steady state, where
 * the room mass is pinned by the cap; the admitted FLOW for the transient, because
 * at the instant a region fills its room is EMPTY and a rate times zero mass cannot
 * hold the cap. The two rules agree at a fixed point, where the room mass is
 * proportional to its inflow.
 *
 * A held or lost cap composes as a PRODUCT of fractions either way, which is what
 * independent blocking gives and what keeps every active cap present in the
 * Jacobian. Several staged caps gating one room compose HARMONICALLY, because the
 * waits a job serves in turn add.
 */
inline FluidDaeLegs fluid_dae_legs(const FluidMomentTerms& t, const FluidDaeGates& gates,
                                   const FluidDaeStaging& stg, const FluidDaeConstraints& con,
                                   const std::vector<std::size_t>& active,
                                   const std::vector<double>& sg, const std::vector<double>& r,
                                   const std::vector<double>& mult, bool staged_flow) {
    const std::size_t nev = r.size(), ns = stg.n, nact = active.size();
    FluidDaeLegs out;
    std::vector<double> whole(nev, 1.0), entry(nev, 1.0);
    std::vector<double> inv_theta(ns, 0.0), flow_of(ns, 0.0);
    std::vector<bool> throttled(ns, false);
    for (std::size_t k = 0; k < nact; ++k) {
        const std::size_t c = active[k];
        const double m = mult[k];
        if (con.staged[c]) {
            for (std::size_t j = 0; j < ns; ++j)
                if (stg.gated_by[c][j]) {
                    throttled[j] = true;
                    if (!staged_flow)
                        inv_theta[j] += (m > 1e-14) ? 1.0 / m
                                                    : std::numeric_limits<double>::infinity();
                }
        } else {
            for (std::size_t e = 0; e < nev; ++e) {
                if (gates.held[c][e]) whole[e] *= m;
                if (gates.loss[c][e]) entry[e] *= m;
            }
        }
    }
    std::vector<bool> split(nev, false);
    for (std::size_t e = 0; e < nev; ++e)
        if (ns && stg.adm[e] && throttled[stg.adm_stage[e]]) split[e] = true;

    out.rup.assign(nev, 0.0);
    out.rin.assign(nev, 0.0);
    for (std::size_t e = 0; e < nev; ++e) {
        // A staged event is NOT suppressed upstream even when a held cap gates it:
        // the job completes upstream service into the waiting room, so the held
        // fraction moves to the room's exit leg below.
        out.rup[e] = r[e] * (split[e] ? 1.0 : whole[e]);
        out.rin[e] = out.rup[e] * entry[e];
    }
    std::vector<double> inflow(ns, 0.0), R(ns, 0.0);
    for (std::size_t e = 0; e < nev; ++e)
        if (split[e]) {
            inflow[stg.adm_stage[e]] += out.rup[e];
            R[stg.adm_stage[e]] += out.rup[e];
        }
    if (staged_flow) {
        for (std::size_t k = 0; k < nact; ++k) {
            const std::size_t c = active[k];
            if (!con.staged[c]) continue;
            double mass = 0.0, tot = 0.0;
            std::size_t nrooms = 0;
            for (std::size_t j = 0; j < ns; ++j)
                if (stg.gated_by[c][j]) { mass += sg[j]; tot += inflow[j]; ++nrooms; }
            if (nrooms == 0) continue;
            for (std::size_t j = 0; j < ns; ++j) {
                if (!stg.gated_by[c][j]) continue;
                double w;
                if (mass > 1e-8) w = sg[j] / mass;
                else if (tot > 1e-8) w = inflow[j] / tot;
                else w = 1.0 / static_cast<double>(nrooms);
                flow_of[j] += mult[k] * w;
            }
        }
    } else {
        for (std::size_t j = 0; j < ns; ++j) {
            const double theta = inv_theta[j] > 0.0 ? 1.0 / inv_theta[j] : 0.0;
            flow_of[j] = theta * sg[j];
        }
    }
    out.drain = flow_of;
    std::vector<double> left(ns, 0.0);
    for (std::size_t e = 0; e < nev; ++e) {
        if (!split[e]) continue;
        const std::size_t j = stg.adm_stage[e];
        const double q = R[j] > 1e-8 ? flow_of[j] * out.rup[e] / R[j] : 0.0;
        // the held fraction gates the room's EXIT: what it stops stays in the room.
        // The lost fraction gates the ARRIVAL: that mass leaves and is destroyed.
        out.rin[e] = q * whole[e] * entry[e];
        left[j] += q * whole[e];
    }
    for (std::size_t j = 0; j < ns; ++j)
        if (R[j] > 1e-8) out.drain[j] = left[j];
    out.ds.assign(ns, 0.0);
    for (std::size_t j = 0; j < ns; ++j)
        // a waiting room with no cap above it must be EMPTY, not merely balanced
        out.ds[j] = throttled[j] ? inflow[j] - out.drain[j] : sg[j];
    return out;
}

/** Population conservation, one row per CLOSED chain, in state space. */
struct FluidDaeConservation {
    Matrix<double> C;         ///< (nchain x nstate)
    std::vector<double> N;    ///< chain populations
};

/**
 * The conserved chains as equations.
 *
 * An open chain has no conserved population and contributes nothing. The EXT
 * coordinates are excluded because the closing representation holds unit mass
 * there as a normalisation constant, not as a job count.
 */
template <typename T>
inline FluidDaeConservation fluid_dae_conservation(const qn::NetworkStruct<T>& sn,
                                                   const FluidMomentTerms& terms,
                                                   const FluidDaeStaging& stg) {
    const std::size_t nstate = terms.nstate;
    const std::size_t M = terms.class_block.size();
    const std::size_t K = M ? terms.class_block[0].size() : 0;

    std::vector<std::size_t> coord_class(nstate, 0), coord_station(nstate, 0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            for (std::size_t s : terms.class_block[i][k]) {
                coord_class[s] = k;
                coord_station[s] = i;
            }

    std::vector<std::vector<double> > rows;
    std::vector<double> nvec;
    const std::size_t nchains = sn.chains.size();
    const std::vector<double> njobs = sn.njobs();
    for (std::size_t ch = 0; ch < nchains; ++ch) {
        std::vector<std::size_t> inch;
        double Nch = 0.0;
        bool finite = true;
        for (std::size_t k = 0; k < K; ++k) {
            if (k >= sn.chains[ch].size() || !sn.chains[ch][k]) continue;
            inch.push_back(k);
            const double nj = (k < njobs.size()) ? njobs[k] : 0.0;
            if (!std::isfinite(nj)) finite = false;
            Nch += nj;
        }
        if (inch.empty() || !finite || Nch <= 0.0) continue;  // open chain
        std::vector<double> row(nstate + stg.n, 0.0);
        bool any = false;
        for (std::size_t s = 0; s < nstate; ++s) {
            if (terms.is_ext[coord_station[s]]) continue;
            if (std::find(inch.begin(), inch.end(), coord_class[s]) == inch.end()) continue;
            row[s] = 1.0;
            any = true;
        }
        // A JOB IN THE WAITING QUEUE IS STILL IN THE CHAIN.
        for (std::size_t j = 0; j < stg.n; ++j)
            if (std::find(inch.begin(), inch.end(), stg.klass[j]) != inch.end())
                row[nstate + j] = 1.0;
        if (!any) continue;
        rows.push_back(row);
        nvec.push_back(Nch);
    }

    FluidDaeConservation out;
    out.C = Matrix<double>(rows.size(), nstate + stg.n, 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (std::size_t s = 0; s < nstate + stg.n; ++s) out.C(r, s) = rows[r][s];
    out.N = nvec;
    return out;
}



/** Options that only the DAE route reads. */
struct FluidDaeOptions {
    /**
     * The simultaneous solve carries one Lyapunov solve per residual evaluation
     * and takes a finite-difference Jacobian over nstate + nclosable unknowns,
     * so its cost is cubic per evaluation and quartic overall. That is
     * affordable at the scale the moment methods run at, but the crossover is
     * lower than the 200 `moment_maxstate` permits, so it gets its own limit.
     */
    std::size_t maxstate = 100;
    std::size_t newton_max = 50;
    /**
     * Largest covariance dimension integrated ALONGSIDE the mean on the
     * transient. The covariance adds nc^2 differential states and the Jacobian
     * is formed by finite differences over all of them, so the cost grows as
     * nc^4. Above this the mean is still integrated as a DAE -- population
     * conservation stays an algebraic equation -- but the variance is held at
     * its stationary value, which is what `minnormal` does for the whole of its
     * transient anyway. Twin of `options.config.dae_maxcov`.
     */
    std::size_t maxcov = 25;
};

/**
 * The controls the DAE route actually reads: the struct a caller pinned, with
 * whatever `options.config` set on top of it.
 *
 * WHY THIS EXISTS. `dae_maxstate` and `dae_maxcov` are `options.config` entries
 * in MATLAB, the JAR and native Python -- they are what a user reaches for to
 * raise a refusal -- and they arrive here on `FluidOptions`, while the route
 * reads `FluidDaeOptions`. Without this step the fields were unreachable from
 * every ordinary entry point (`solver_fluid_run_analyzer` hands the defaults), so the
 * two limits were fixed at 100 and 25 whatever the options said. Zero on
 * `FluidOptions` means NOT SET, so an explicit struct still wins where the
 * options are silent.
 *
 * THE NEWTON CAP IS NOT A KNOB BUT A DIVERGENCE. `max(50, iter_max)` is what
 * all three references solve with, and `iter_max` defaults to 200, so leaving
 * it at 50 here was not a missing option: it was a different solver, one that
 * refuses a fixed point the reference reaches on its 60th step.
 */
inline FluidDaeOptions fluid_dae_options(const FluidOptions& opt,
                                         const FluidDaeOptions& dopt) {
    FluidDaeOptions d = dopt;
    if (opt.dae_maxstate > 0) d.maxstate = opt.dae_maxstate;
    if (opt.dae_maxcov > 0) d.maxcov = opt.dae_maxcov;
    if (opt.iter_max > d.newton_max) d.newton_max = opt.iter_max;
    return d;
}

/**
 * Stations whose variance enters the drift.
 *
 * A delay or a source has no min() to close. A station that cannot fill its
 * servers has min(n,c) = n on its whole support, so closing it is not an
 * improvement but an error (see FLUID_MOMENT_TERMS). Those are held at zero and
 * are not unknowns, which keeps the Newton system as small as the closure is.
 */
inline std::vector<std::size_t> fluid_dae_closable(const FluidMomentTerms& t) {
    std::vector<std::size_t> idx;
    const std::size_t M = t.station_block.size();
    for (std::size_t i = 0; i < M; ++i) {
        if (t.is_ext[i] || t.min_exact[i] || t.station_block[i].empty()) continue;
        if (t.sys.sched[i] == lang::SchedStrategy::INF) continue;
        if (!std::isfinite(t.S[i])) continue;
        idx.push_back(i);
    }
    return idx;
}

/** Project a state-level covariance onto the per-station variances the drift reads. */
inline std::vector<double> fluid_dae_sigma_from(const Matrix<double>& Sigma,
                                                const FluidMomentTerms& t,
                                                const std::vector<std::size_t>& cidx) {
    std::vector<double> s2(t.station_block.size(), 0.0);
    for (std::size_t j = 0; j < cidx.size(); ++j) {
        const std::vector<std::size_t>& blk = t.station_block[cidx[j]];
        double acc = 0.0;
        for (std::size_t a = 0; a < blk.size(); ++a)
            for (std::size_t b = 0; b < blk.size(); ++b) acc += Sigma(blk[a], blk[b]);
        s2[cidx[j]] = std::max(0.0, acc);
    }
    return s2;
}

/** Everything the residual needs, gathered so the Newton can stay generic. */
struct FluidDaeSystem {
    const FluidMomentTerms* terms = nullptr;
    const FluidDaeConservation* cons = nullptr;
    const FluidDaeConstraints* con = nullptr;
    const FluidDaeStaging* stg = nullptr;
    const FluidDaeGates* gates = nullptr;
    std::vector<std::size_t> cidx;        ///< closable stations
    std::vector<std::size_t> active;      ///< binding capacity constraints
    std::size_t nstate = 0;
    /**
     * The tangent space of the caps that CLAMP, or an empty matrix. A cap that
     * holds the job upstream or loses it fixes its own combination of the state
     * for as long as it binds, so that combination does not fluctuate: WITHOUT the
     * projection the Lyapunov solve has no solution at all, because with the
     * multiplier held constant the drift along the constrained direction is
     * neutral. A STAGED cap is not clamped -- its room is a state and supplies the
     * restoring force -- and must not be projected.
     */
    Matrix<double> clampT = Matrix<double>(0, 0, 0.0);
    /** u = [x; staging; sigma2; mult] */
    std::size_t nstaging() const { return stg ? stg->n : 0; }
};

/**
 * Orthogonal projector onto the subspace the CLAMPING caps leave free:
 * `I - R'(RR')^-1 R` over the covariance coordinates. Empty when no active cap
 * clamps, which is every model without a station buffer and every region model.
 */
inline Matrix<double> fluid_dae_clamp_tangent(const FluidDaeConstraints& con,
                                              const std::vector<std::size_t>& active,
                                              const FluidMomentTerms& t) {
    std::vector<std::size_t> rows;
    for (std::size_t k = 0; k < active.size(); ++k)
        if (!con.staged[active[k]]) rows.push_back(active[k]);
    if (rows.empty()) return Matrix<double>(0, 0, 0.0);
    const std::size_t nc = t.cov_idx.size();
    Matrix<double> R(rows.size(), nc, 0.0);
    bool any = false;
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (std::size_t j = 0; j < nc; ++j) {
            R(i, j) = con.A(rows[i], t.cov_idx[j]);
            if (std::fabs(R(i, j)) > 1e-14) any = true;
        }
    if (!any) return Matrix<double>(0, 0, 0.0);
    // (R R')^-1 by a small Gauss-Jordan: the block is one row per clamping cap
    const std::size_t m = rows.size();
    Matrix<double> G(m, 2 * m, 0.0);
    for (std::size_t a = 0; a < m; ++a) {
        for (std::size_t b = 0; b < m; ++b) {
            double acc = 0.0;
            for (std::size_t j = 0; j < nc; ++j) acc += R(a, j) * R(b, j);
            G(a, b) = acc;
        }
        G(a, m + a) = 1.0;
    }
    for (std::size_t a = 0; a < m; ++a) {
        std::size_t piv = a;
        for (std::size_t b = a + 1; b < m; ++b)
            if (std::fabs(G(b, a)) > std::fabs(G(piv, a))) piv = b;
        if (std::fabs(G(piv, a)) < 1e-300) return Matrix<double>(0, 0, 0.0);
        if (piv != a)
            for (std::size_t j = 0; j < 2 * m; ++j) std::swap(G(a, j), G(piv, j));
        const double d = G(a, a);
        for (std::size_t j = 0; j < 2 * m; ++j) G(a, j) /= d;
        for (std::size_t b = 0; b < m; ++b) {
            if (b == a) continue;
            const double f = G(b, a);
            if (f == 0.0) continue;
            for (std::size_t j = 0; j < 2 * m; ++j) G(b, j) -= f * G(a, j);
        }
    }
    Matrix<double> T(nc, nc, 0.0);
    for (std::size_t i = 0; i < nc; ++i) T(i, i) = 1.0;
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < nc; ++j) {
            double acc = 0.0;
            for (std::size_t a = 0; a < m; ++a)
                for (std::size_t b = 0; b < m; ++b) acc += R(a, i) * G(a, m + b) * R(b, j);
            T(i, j) -= acc;
        }
    return T;
}

/**
 * The coupled algebraic system, stacked: drift, conservation, closure consistency.
 *
 * Returns false when the closure cannot be evaluated at this iterate -- the
 * Lyapunov solve throws on a non-hyperbolic fixed point -- so that the line
 * search can back off rather than the whole solve failing.
 *
 * THE UNKNOWNS ARE [x; staging; sigma2; mult], with MULT one multiplier per ACTIVE
 * cap: a FRACTION of the admissions allowed for a cap that holds the job upstream
 * or loses it, the RATE its waiting room drains at for one that stages it.
 */
inline bool fluid_dae_residual(const FluidDaeSystem& sysd, const std::vector<double>& u,
                               std::vector<double>& G, std::vector<double>* rates_out,
                               bool rethrow = false, std::vector<double>* fire_out = nullptr) {
    const FluidMomentTerms& t = *sysd.terms;
    const std::size_t n = sysd.nstate, ns = sysd.nstaging();
    const std::size_t ncl = sysd.cidx.size(), nact = sysd.active.size();
    const std::vector<double> x(u.begin(), u.begin() + n);
    const std::vector<double> sg(u.begin() + n, u.begin() + n + ns);

    FluidClosure cl;
    cl.sigma2.assign(t.station_block.size(), 0.0);
    cl.cov.assign(t.station_block.size(), Matrix<double>(0, 0, 0.0));
    // Clamped below only. FLUID_DAE_NEWTON projects the iterate into the
    // feasible box, so this is a guard rather than the mechanism: clamping alone
    // would flatten the Jacobian at the boundary and pin the unknown there.
    for (std::size_t j = 0; j < ncl; ++j)
        cl.sigma2[sysd.cidx[j]] = std::max(0.0, u[n + ns + j]);
    std::vector<double> mult(nact, 0.0);
    for (std::size_t k = 0; k < nact; ++k) mult[k] = std::max(0.0, u[n + ns + ncl + k]);

    std::vector<double> r;
    FluidDaeLegs lg;
    Matrix<double> Sigma(0, 0, 0.0);
    try {
        r = fluid_moment_rates(t, x, cl);
        lg = fluid_dae_legs(t, *sysd.gates, *sysd.stg, *sysd.con, sysd.active, sg, r, mult, false);
        const Matrix<double> A = fluid_drift_jacobian(t, x, cl);
        // THE FIRING RATE, not the nominal one: a cap that holds the job upstream
        // suppresses the event itself, so the diffusion counts what happens.
        const Matrix<double>* cT = sysd.clampT.rows() ? &sysd.clampT : nullptr;
        Sigma = fluid_moment_lyapunov(t, A, lg.rup, cT);
    } catch (const std::exception&) {
        // The line search reads `false` as "back off", but the FIRST evaluation
        // has nowhere to back off to: there the underlying reason -- a
        // non-hyperbolic fixed point, or LAPACK missing under the Lyapunov
        // solve -- is the answer and must not be replaced by a generic one.
        if (rethrow) throw;
        return false;
    }
    const std::vector<double> s2new = fluid_dae_sigma_from(Sigma, t, sysd.cidx);

    const std::size_t nev = r.size();
    std::vector<double> drift(n, 0.0);
    for (std::size_t s = 0; s < n; ++s) {
        double acc = 0.0;
        for (std::size_t e = 0; e < nev; ++e) {
            if (!sysd.gates->has) acc += t.D(s, e) * lg.rup[e];
            else
                // THE THREE LEGS: the removal leaves at the rate the event FIRES and
                // the arrival lands at the rate mass actually ARRIVES. A LOST arrival
                // is returned to the source pool rather than destroyed -- the EXT
                // coordinate is a normalisation and its row a real equation -- while
                // a job lost leaving a REAL station is destroyed.
                acc += sysd.gates->Dn(s, e) * lg.rup[e] + sysd.gates->DnExt(s, e) * lg.rin[e] +
                       sysd.gates->Dp(s, e) * lg.rin[e];
        }
        drift[s] = acc;
    }

    G.clear();
    G.insert(G.end(), drift.begin(), drift.end());
    G.insert(G.end(), lg.ds.begin(), lg.ds.end());
    const Matrix<double>& C = sysd.cons->C;
    for (std::size_t c = 0; c < C.rows(); ++c) {
        double acc = 0.0;
        for (std::size_t s = 0; s < n; ++s) acc += C(c, s) * x[s];
        for (std::size_t j = 0; j < ns; ++j) acc += C(c, n + j) * sg[j];
        G.push_back(acc - sysd.cons->N[c]);
    }
    for (std::size_t j = 0; j < ncl; ++j)
        G.push_back(cl.sigma2[sysd.cidx[j]] - s2new[sysd.cidx[j]]);
    // THE UNDIFFERENTIATED CONSTRAINT, on purpose: at a fixed point every flow
    // already balances, so d/dt(Ax)=0 is satisfied by anything and pins no
    // multiplier. A x = b does pin it, through the fixed point's dependence on it.
    for (std::size_t k = 0; k < nact; ++k) {
        double acc = 0.0;
        for (std::size_t s = 0; s < n; ++s) acc += sysd.con->A(sysd.active[k], s) * x[s];
        for (std::size_t j = 0; j < ns; ++j) acc += sysd.con->As(sysd.active[k], j) * sg[j];
        G.push_back(acc - sysd.con->b[sysd.active[k]]);
    }

    // RATES_OUT is the flow that actually CROSSES each event, which is what the
    // throughput table reports; FIRE_OUT the rate the event fires at.
    if (rates_out) *rates_out = lg.rin;
    if (fire_out) *fire_out = lg.rup;
    return true;
}

/** Onto the feasible box: the state is free, the variances are not. */
inline void fluid_dae_project(std::vector<double>& u, std::size_t nfree) {
    for (std::size_t i = nfree; i < u.size(); ++i) u[i] = std::max(0.0, u[i]);
}

/** What the Newton reports about where it stopped. */
struct FluidDaeNewtonInfo {
    std::size_t iters = 0;
    double residual = std::numeric_limits<double>::infinity();
    bool converged = false;
};

/**
 * Damped PROJECTED Newton with a finite-difference Jacobian and an Armijo
 * backtrack on the residual norm.
 *
 * The step is solved in least squares: the drift block is rank deficient by
 * exactly the number of conserved chains, and the constraint rows restore that
 * rank, so the stacked system is consistent and overdetermined rather than
 * square.
 *
 * WHY PROJECTED, AND NOT MERELY CLAMPED. The unknowns past NFREE are variances,
 * which may not go negative, and the residual reads them through max(0,.).
 * Clamping inside the residual alone is a trap: once an iterate goes negative
 * the residual stops depending on it, so the finite-difference column is exactly
 * zero, there is no derivative to climb back on, and the unknown is pinned at
 * the boundary for good. Projecting the ITERATE keeps every evaluation inside
 * the box, where the forward difference across max(0,.) is live even at zero.
 */
inline FluidDaeNewtonInfo fluid_dae_newton(const FluidDaeSystem& sysd, std::vector<double>& u,
                                           double tol, std::size_t maxit) {
    FluidDaeNewtonInfo info;
    const std::size_t nfree = sysd.nstate;
    fluid_dae_project(u, nfree);
    std::vector<double> G;
    fluid_dae_residual(sysd, u, G, nullptr, /*rethrow=*/true);
    double resnorm = 0.0;
    for (double g : G) resnorm = std::max(resnorm, std::fabs(g));

    const std::size_t nun = u.size(), m = G.size();
    while (resnorm >= tol && info.iters < maxit) {
        ++info.iters;
        Matrix<double> J(m, nun, 0.0);
        std::vector<double> Gp;
        for (std::size_t j = 0; j < nun; ++j) {
            const double h = std::max(1e-7 * std::fabs(u[j]), 1e-9);
            std::vector<double> up = u;
            up[j] += h;
            if (!fluid_dae_residual(sysd, up, Gp, nullptr)) continue;  // column left at zero
            for (std::size_t i = 0; i < m; ++i) J(i, j) = (Gp[i] - G[i]) / h;
        }
        std::vector<double> du;
        try {
            du = lstsq(J, G).x;
        } catch (const std::exception&) {
            break;
        }
        for (double& d : du) d = -d;

        double lam = 1.0;
        bool stepped = false;
        for (int ls = 0; ls < 25; ++ls) {
            std::vector<double> un(nun);
            for (std::size_t j = 0; j < nun; ++j) un[j] = u[j] + lam * du[j];
            fluid_dae_project(un, nfree);  // project BEFORE evaluating, so the
            std::vector<double> Gn;        // accepted point and its residual agree
            if (fluid_dae_residual(sysd, un, Gn, nullptr)) {
                double rn = 0.0;
                bool ok = true;
                for (double g : Gn) {
                    if (!std::isfinite(g)) { ok = false; break; }
                    rn = std::max(rn, std::fabs(g));
                }
                if (ok && rn < resnorm * (1.0 - 1e-4 * lam)) {
                    u = un; G = Gn; resnorm = rn;
                    stepped = true;
                    break;
                }
            }
            lam *= 0.5;
        }
        if (!stepped) break;  // no descent along this direction
    }
    info.residual = resnorm;
    info.converged = resnorm < tol;
    return info;
}


namespace detail {

/**
 * The transient DAE, as RODAS sees it.
 *
 * RODAS's callbacks are plain C function pointers with no user-data channel,
 * so the system under integration is reached through this pointer. That is
 * consistent with the vendored solver itself, whose upstream COMMON blocks make
 * it single-integration-at-a-time regardless.
 */
struct FluidDaeTransient {
    const FluidMomentTerms* terms = nullptr;
    const FluidDaeConservation* cons = nullptr;
    FluidClosure closure;                  ///< the stationary variance, when held
    std::vector<std::size_t> alg_row;      ///< one differential row per chain, replaced
    std::size_t nstate = 0;
    /**
     * THE COVARIANCE, INTEGRATED RATHER THAN HELD. `withcov` selects between the
     * two: below `maxcov` the nc x nc block on `cov_idx` joins the state as
     * ordinary DIFFERENTIAL rows (mass 1, unlike the one algebraic row per
     * chain) and advances by dS = A S + S A' + D diag(r) D', so the drift is
     * read at the variance the trajectory actually HAS at each instant. Above
     * it the block is dropped and `closure` supplies the stationary variance
     * instead -- which is what `minnormal` uses for the whole of its transient,
     * so the fallback is the reference's own and not a different closure.
     */
    bool withcov = false;
    std::vector<std::size_t> cov_idx;
    std::size_t nc = 0;
    std::vector<std::size_t> closable;     ///< stations whose variance enters the drift
    /// The output grid, and the states RODAS's dense output was read at. The
    /// cursor is the next grid point still owed a value; it only advances, so a
    /// step that spans several grid points fills them all in one call.
    std::vector<double> grid;
    std::size_t cursor = 0;
    std::vector<std::vector<double> > out;

    /**
     * THE CAPS, AND THE MODE THIS SEGMENT IS IN. Under a cap the transient is a
     * HYBRID DAE: the system switches every time a cap starts or stops binding, so
     * the horizon is covered by SEGMENTS, each one index-1 with a fixed binding set,
     * ending at a located crossing. The multiplier of each cap is an ALGEBRAIC
     * unknown carried in the state with a zero mass row, and its equation is the
     * DIFFERENTIATED constraint d/dt(A x + As s) = 0 -- the undifferentiated form is
     * index 2, which RODAS does not solve. A staged cap's unknown is the admitted
     * FLOW, not the drain rate the steady state solves for: at the instant a region
     * fills its room is EMPTY, and a rate times zero mass cannot hold the cap.
     */
    const FluidDaeConstraints* con = nullptr;
    const FluidDaeGates* gates = nullptr;
    const FluidDaeStaging* stg = nullptr;
    std::vector<std::size_t> active;
    std::vector<char> armed;          ///< a staged cap whose room has really filled
    std::vector<char> gated_rooms;
    std::vector<double> inert;        ///< one for a fraction, zero for a flow
    std::size_t ncon = 0, nstg = 0, o_sg = 0, o_m = 0, o_cov = 0;
    /// the event functions at the last accepted step, and the crossing located
    std::vector<double> gprev;
    int hit = -1;
    double hit_t = 0.0;
    std::vector<double> hit_z;
};

/**
 * The event functions, all of them: a cap that is NOT active is watched for
 * REACHING its bound, one that IS active for stopping to bind -- a fraction above
 * one, or a room that has emptied.
 *
 * ARMED IS A LATCH, and it has to be: a staged cap activates with an EMPTY room, so
 * "the room emptied" is true at the instant it starts binding and the release would
 * fire immediately. The latch only ever goes false to true.
 */
inline std::vector<double> fluid_dae_cap_events(const FluidDaeTransient& S,
                                                const double* z) {
    std::vector<double> g(S.ncon, 1.0);
    for (std::size_t c = 0; c < S.ncon; ++c) {
        const bool is_active =
            std::find(S.active.begin(), S.active.end(), c) != S.active.end();
        if (is_active) {
            if (S.con->staged[c]) {
                double mass = 0.0;
                for (std::size_t j = 0; j < S.nstg; ++j)
                    if (S.stg->gated_by[c][j]) mass += z[S.o_sg + j];
                g[c] = S.armed[c] ? mass : 1.0;
            } else {
                g[c] = 1.0 - z[S.o_m + c];
            }
        } else {
            double val = 0.0;
            for (std::size_t s = 0; s < S.nstate; ++s) val += S.con->A(c, s) * z[s];
            for (std::size_t j = 0; j < S.nstg; ++j) val += S.con->As(c, j) * z[S.o_sg + j];
            g[c] = S.con->b[c] - val;
        }
    }
    return g;
}
inline FluidDaeTransient*& fluid_dae_active() {
    static FluidDaeTransient* p = nullptr;
    return p;
}

/** The closure to read the drift at: the trajectory's own variance, or the held one. */
inline FluidClosure fluid_dae_closure_at(const FluidDaeTransient& S,
                                         const rodas_impl::doublereal* y) {
    if (!S.withcov) return S.closure;
    FluidClosure cl;
    cl.sigma2.assign(S.terms->station_block.size(), 0.0);
    cl.cov.assign(S.terms->station_block.size(), Matrix<double>(0, 0, 0.0));
    Matrix<double> Sigma(S.terms->nstate, S.terms->nstate, 0.0);
    for (std::size_t a = 0; a < S.nc; ++a)
        for (std::size_t b = 0; b < S.nc; ++b)
            // the equation preserves symmetry; rounding does not
            Sigma(S.cov_idx[a], S.cov_idx[b]) =
                0.5 * (y[S.o_cov + a * S.nc + b] + y[S.o_cov + b * S.nc + a]);
    const std::vector<double> s2 = fluid_dae_sigma_from(Sigma, *S.terms, S.closable);
    for (std::size_t i = 0; i < s2.size(); ++i) cl.sigma2[i] = s2[i];
    return cl;
}

/** M y' = f: the drift, with the algebraic rows carrying the constraint residual. */
inline int fluid_dae_fcn(rodas_impl::integer*, rodas_impl::doublereal*,
                         rodas_impl::doublereal* y, rodas_impl::doublereal* f,
                         rodas_impl::doublereal*, rodas_impl::integer*) {
    FluidDaeTransient& S = *fluid_dae_active();
    const std::vector<double> x(y, y + S.nstate);
    const FluidClosure cl = fluid_dae_closure_at(S, y);
    std::vector<double> dx;
    FluidDaeLegs lg;
    if (S.ncon == 0) {
        dx = fluid_moment_drift(*S.terms, x, cl);
    } else {
        const std::vector<double> sg(y + S.o_sg, y + S.o_sg + S.nstg);
        std::vector<double> mact(S.active.size(), 0.0);
        for (std::size_t k = 0; k < S.active.size(); ++k) mact[k] = y[S.o_m + S.active[k]];
        const std::vector<double> r = fluid_moment_rates(*S.terms, x, cl);
        lg = fluid_dae_legs(*S.terms, *S.gates, *S.stg, *S.con, S.active, sg, r, mact, true);
        dx.assign(S.nstate, 0.0);
        for (std::size_t s = 0; s < S.nstate; ++s) {
            double acc = 0.0;
            for (std::size_t e = 0; e < r.size(); ++e)
                acc += S.gates->Dn(s, e) * lg.rup[e] + S.gates->DnExt(s, e) * lg.rin[e] +
                       S.gates->Dp(s, e) * lg.rin[e];
            dx[s] = acc;
        }
    }
    for (std::size_t i = 0; i < S.nstate; ++i) f[i] = dx[i];
    const Matrix<double>& C = S.cons->C;
    for (std::size_t c = 0; c < C.rows(); ++c) {
        double acc = 0.0;
        for (std::size_t s = 0; s < S.nstate; ++s) acc += C(c, s) * x[s];
        for (std::size_t j = 0; j < S.nstg; ++j) acc += C(c, S.nstate + j) * y[S.o_sg + j];
        f[S.alg_row[c]] = acc - S.cons->N[c];  // the mass matrix zeroed this row
    }
    if (S.ncon) {
        for (std::size_t j = 0; j < S.nstg; ++j)
            f[S.o_sg + j] = S.gated_rooms[j] ? lg.ds[j] : y[S.o_sg + j];
        for (std::size_t c = 0; c < S.ncon; ++c) {
            const bool is_active =
                std::find(S.active.begin(), S.active.end(), c) != S.active.end();
            if (!is_active) {
                f[S.o_m + c] = y[S.o_m + c] - S.inert[c];
                continue;
            }
            double acc = 0.0;
            for (std::size_t s = 0; s < S.nstate; ++s) acc += S.con->A(c, s) * dx[s];
            for (std::size_t j = 0; j < S.nstg; ++j) acc += S.con->As(c, j) * lg.ds[j];
            f[S.o_m + c] = acc;
        }
    }
    if (S.withcov) {
        // dS = A S + S A' + D diag(r) D', on the coordinates that carry a real
        // population. ORDINARY DIFFERENTIAL ROWS: only the algebraic rows above are
        // zeroed by the mass matrix, so these stay at 1.
        const Matrix<double> A = fluid_drift_jacobian(*S.terms, x, cl);
        const std::vector<double> r =
            S.ncon ? lg.rup : fluid_moment_rates(*S.terms, x, cl);
        for (std::size_t a = 0; a < S.nc; ++a)
            for (std::size_t b = 0; b < S.nc; ++b) {
                double acc = 0.0;
                for (std::size_t l = 0; l < S.nc; ++l)
                    acc += A(S.cov_idx[a], S.cov_idx[l]) * y[S.o_cov + l * S.nc + b]
                         + y[S.o_cov + a * S.nc + l] * A(S.cov_idx[b], S.cov_idx[l]);
                for (std::size_t e = 0; e < r.size(); ++e)
                    acc += S.terms->D(S.cov_idx[a], e) * r[e] * S.terms->D(S.cov_idx[b], e);
                f[S.o_cov + a * S.nc + b] = acc;
            }
    }
    return 0;
}

/** Analytic Jacobian: the drift's, with the algebraic rows replaced by C. */
inline int fluid_dae_jac(rodas_impl::integer*, rodas_impl::doublereal*,
                         rodas_impl::doublereal* y, rodas_impl::doublereal* dfy,
                         rodas_impl::integer* ldfy, rodas_impl::doublereal*,
                         rodas_impl::integer*) {
    FluidDaeTransient& S = *fluid_dae_active();
    const std::vector<double> x(y, y + S.nstate);
    const Matrix<double> A = fluid_drift_jacobian(*S.terms, x, S.closure);
    const int ld = *ldfy;
    for (std::size_t i = 0; i < S.nstate; ++i)
        for (std::size_t j = 0; j < S.nstate; ++j) dfy[i + j * ld] = A(i, j);
    const Matrix<double>& C = S.cons->C;
    for (std::size_t c = 0; c < C.rows(); ++c)
        for (std::size_t j = 0; j < S.nstate; ++j) dfy[S.alg_row[c] + j * ld] = C(c, j);
    return 0;
}

/** The singular mass matrix, banded with MLMAS=MUMAS=0, i.e. the diagonal. */
inline int fluid_dae_mas(rodas_impl::integer* n, rodas_impl::doublereal* am,
                         rodas_impl::integer* lmas, rodas_impl::doublereal*,
                         rodas_impl::integer*) {
    FluidDaeTransient& S = *fluid_dae_active();
    const int ld = *lmas;
    for (int j = 0; j < *n; ++j) am[0 + j * ld] = 1.0;
    for (std::size_t c = 0; c < S.alg_row.size(); ++c) am[0 + S.alg_row[c] * ld] = 0.0;
    // a room with no cap above it is held EMPTY by an algebraic row, and every
    // multiplier is algebraic: its equation is the differentiated constraint
    for (std::size_t j = 0; j < S.nstg; ++j)
        if (!S.gated_rooms[j]) am[0 + (S.o_sg + j) * ld] = 0.0;
    for (std::size_t c = 0; c < S.ncon; ++c) am[0 + (S.o_m + c) * ld] = 0.0;
    return 0;
}

/**
 * The trajectory, read off RODAS's own dense output rather than off its steps.
 *
 * A Rosenbrock method chooses its step from the local error, so its accepted
 * points are wherever the stiffness put them and never the ones a caller asked
 * for. `contro_` is the third-order interpolant RODAS carries for exactly this,
 * valid over the step just accepted, so asking for an arbitrary grid costs no
 * extra step and no interpolation of the caller's own.
 *
 * THE FIRST CALL IS MADE BEFORE ANY STEP (`rodas.f` calls SOLOUT once with
 * XOLD = X = t0 and NACCPT = 0), so `cont` holds no coefficients yet and the
 * state at that point is `y` itself. Reading `contro_` there would interpolate
 * uninitialised work space.
 */
inline int fluid_dae_solout(rodas_impl::integer* nr, rodas_impl::doublereal* xold,
                            rodas_impl::doublereal* x, rodas_impl::doublereal* y,
                            rodas_impl::doublereal* cont, rodas_impl::integer* lrc,
                            rodas_impl::integer*, rodas_impl::doublereal*,
                            rodas_impl::integer*, rodas_impl::integer* irtrn) {
    FluidDaeTransient& S = *fluid_dae_active();
    const std::size_t nz = S.o_cov + (S.withcov ? S.nc * S.nc : 0);
    double stop_at = std::numeric_limits<double>::infinity();
    if (S.ncon && *nr > 1) {
        // arm a staged release only once its room has really filled
        for (std::size_t k = 0; k < S.active.size(); ++k) {
            const std::size_t c = S.active[k];
            if (!S.con->staged[c] || S.armed[c]) continue;
            double mass = 0.0;
            for (std::size_t j = 0; j < S.nstg; ++j)
                if (S.stg->gated_by[c][j]) mass += y[S.o_sg + j];
            if (mass > 1e-8) S.armed[c] = 1;
        }
        const std::vector<double> gcur = fluid_dae_cap_events(S, y);
        int cross = -1;
        for (std::size_t c = 0; c < S.ncon; ++c)
            if (S.gprev.size() == S.ncon && S.gprev[c] > 0.0 && gcur[c] <= 0.0) {
                cross = static_cast<int>(c);
                break;
            }
        S.gprev = gcur;
        if (cross >= 0) {
            // BISECT ON THE INTERPOLANT rather than on the integration: RODAS carries
            // a third-order one over the step it has just accepted, so locating the
            // crossing costs no extra step and the state at it is the integrator's.
            double lo = *xold;
            double hi = *x;
            std::vector<double> ymid(nz, 0.0);
            for (int bit = 0; bit < 60; ++bit) {
                const double mid = 0.5 * (lo + hi);
                for (std::size_t i = 0; i < nz; ++i) {
                    rodas_impl::integer ii = static_cast<rodas_impl::integer>(i + 1);
                    double tq = mid;
                    ymid[i] = rodas_impl::contro_(&ii, &tq, cont, lrc);
                }
                const std::vector<double> gmid = fluid_dae_cap_events(S, ymid.data());
                if (gmid[static_cast<std::size_t>(cross)] > 0.0) lo = mid;
                else hi = mid;
                if (hi - lo <= 1e-12 * std::max(1.0, std::fabs(hi))) break;
            }
            std::vector<double> yhit(nz, 0.0);
            for (std::size_t i = 0; i < nz; ++i) {
                rodas_impl::integer ii = static_cast<rodas_impl::integer>(i + 1);
                double tq = hi;
                yhit[i] = rodas_impl::contro_(&ii, &tq, cont, lrc);
            }
            S.hit = cross;
            S.hit_t = hi;
            S.hit_z = yhit;
            stop_at = hi;
        }
    } else if (S.ncon) {
        S.gprev = fluid_dae_cap_events(S, y);
    }
    // THE OVERSHOOTING STEP IS NOT REPORTED: grid points beyond the located
    // crossing belong to the NEXT segment, which starts from it.
    while (S.cursor < S.grid.size() && S.grid[S.cursor] <= *x + 1e-13 &&
           S.grid[S.cursor] <= stop_at + 1e-13) {
        double tq = S.grid[S.cursor];
        std::vector<double> xs(nz, 0.0);
        for (std::size_t i = 0; i < nz; ++i) {
            if (*nr <= 1) {
                xs[i] = y[i];
            } else {
                rodas_impl::integer ii = static_cast<rodas_impl::integer>(i + 1);
                xs[i] = rodas_impl::contro_(&ii, &tq, cont, lrc);
            }
        }
        S.out.push_back(xs);
        ++S.cursor;
    }
    // RODAS READS THE STOP THROUGH IRTRN, not through the return value: rodas.f's
    // SOLOUT is a subroutine and the f2c translation keeps that convention, so a
    // located crossing has to be written into the argument or the integration walks
    // straight past the cap it just found.
    if (S.hit >= 0 && irtrn) *irtrn = -1;
    return 0;
}
inline int fluid_dae_dfx(rodas_impl::integer*, rodas_impl::doublereal*,
                         rodas_impl::doublereal*, rodas_impl::doublereal*,
                         rodas_impl::doublereal*, rodas_impl::integer*) { return 0; }

}  // namespace detail

/**
 * Integrate the closure as an index-1 DAE, with RODAS, and report the state at
 * every point of `grid`.
 *
 * One differential equation per closed chain is redundant -- the rows of D sum
 * to zero there -- so one is REPLACED by the constraint rather than added to it.
 * The row dropped is the one carrying the most mass at t=0, which keeps the
 * algebraic equation away from a coordinate that is identically zero. Chains
 * partition the classes, so no index is claimed twice.
 *
 * `grid` must be increasing and end at the horizon; the trajectory comes from
 * RODAS's own dense output (`detail::fluid_dae_solout`), so the grid costs no
 * extra step. The LAST entry is also the return value of the integration
 * proper, which is what a steady-state caller asking for a single point wants.
 *
 * @return one state vector per grid point, in grid order
 */
inline std::vector<std::vector<double> > fluid_dae_integrate(
    const FluidMomentTerms& terms, const FluidDaeConservation& cons,
    const FluidClosure& closure, const std::vector<double>& x0,
    const std::vector<double>& grid, double tol, bool withcov = false,
    const std::vector<std::size_t>& closable = std::vector<std::size_t>()) {
    using namespace rodas_impl;
    const std::size_t n = terms.nstate;

    detail::FluidDaeTransient S;
    S.terms = &terms;
    S.cons = &cons;
    S.closure = closure;
    S.nstate = n;
    S.withcov = withcov;
    S.cov_idx = terms.cov_idx;
    S.nc = terms.cov_idx.size();
    S.closable = closable;
    if (S.nc == 0) S.withcov = false;
    S.o_sg = n;
    S.o_m = n;
    S.o_cov = n;
    const std::size_t nz = n + (S.withcov ? S.nc * S.nc : 0);
    for (std::size_t c = 0; c < cons.C.rows(); ++c) {
        std::size_t best = n;
        double bestv = -1.0;
        for (std::size_t s = 0; s < n; ++s) {
            if (cons.C(c, s) == 0.0) continue;
            bool taken = false;
            for (std::size_t d = 0; d < S.alg_row.size(); ++d)
                if (S.alg_row[d] == s) taken = true;
            if (taken) continue;
            if (x0[s] > bestv) { bestv = x0[s]; best = s; }
        }
        if (best == n) throw NumericError("fluid_dae_integrate: a chain has no free coordinate");
        S.alg_row.push_back(best);
    }
    if (grid.empty()) throw InputError("fluid_dae_integrate: the output grid is empty");
    for (std::size_t j = 1; j < grid.size(); ++j)
        if (!(grid[j] > grid[j - 1]))
            throw InputError("fluid_dae_integrate: the output grid must be increasing");
    S.grid = grid;
    detail::fluid_dae_active() = &S;

    // LWORK per the documented formula N*(LJAC+LMAS+LE1+14)+20 with a full
    // Jacobian (LJAC=LE1=N) and a diagonal mass matrix (LMAS=1). The 14 already
    // covers CONT's 4*N, so dense output needs no extra room.
    const integer N = static_cast<integer>(nz);
    const std::size_t lwork = nz * (nz + 1 + nz + 14) + 20 + 32;
    const std::size_t liwork = nz + 20 + 32;
    std::vector<doublereal> y(nz, 0.0), work(lwork, 0.0), rpar(1, 0.0);
    for (std::size_t i = 0; i < n && i < x0.size(); ++i) y[i] = x0[i];
    // Sigma(0) = 0 is the consistent initialisation, and the physically right
    // one: the population at t=0 is a known deterministic state, so it has no
    // variance. C x0 = N holds by construction, so the algebraic rows are
    // satisfied at t=0 and RODAS needs no separate consistency solve.
    std::vector<integer> iwork(liwork, 0), ipar(1, 0);
    doublereal rtol = tol, atol = tol * 1e-2, x = 0.0, xend = grid.back(), h = 1e-6;
    // THE ANALYTIC JACOBIAN IS THE DRIFT'S, so it is only the whole Jacobian
    // while the covariance is HELD. Once the nc^2 covariance rows join the
    // state they have derivatives of their own, and supplying a Jacobian that
    // is right on one block and zero on the other is worse than supplying none:
    // RODAS would take it as exact. So the covariance run differences it
    // numerically, which is also what ode15s does for the same rows in the
    // MATLAB reference.
    integer itol = 0, ifcn = 0, ijac = S.withcov ? 0 : 1, mljac = N, mujac = N, idfx = 0;
    integer imas = 1, mlmas = 0, mumas = 0, iout = 1, idid = 0;
    integer lw = static_cast<integer>(lwork), liw = static_cast<integer>(liwork);

    try {
        rodas_(const_cast<integer*>(&N), (U_fp)detail::fluid_dae_fcn, &ifcn, &x, y.data(),
               &xend, &h, &rtol, &atol, &itol,
               (U_fp)detail::fluid_dae_jac, &ijac, &mljac, &mujac,
               (U_fp)detail::fluid_dae_dfx, &idfx,
               (U_fp)detail::fluid_dae_mas, &imas, &mlmas, &mumas,
               (U_fp)detail::fluid_dae_solout, &iout,
               work.data(), &lw, iwork.data(), &liw, rpar.data(), ipar.data(), &idid);
    } catch (...) {
        // The active pointer is a global, so an exception thrown out of a
        // callback must not leave it dangling for the next integration.
        detail::fluid_dae_active() = nullptr;
        throw;
    }
    detail::fluid_dae_active() = nullptr;

    if (idid != 1)
        throw NumericError("fluid_dae_integrate: RODAS returned idid=" + std::to_string(idid));
    // The horizon itself is filled from `y` rather than from the interpolant:
    // RODAS lands exactly on XEND, and the last accepted step's dense output is
    // evaluated at its own right endpoint, where the two agree to rounding.
    if (S.out.size() + 1 == grid.size()) S.out.push_back(std::vector<double>(y.begin(), y.end()));
    if (S.out.size() != grid.size())
        throw NumericError("fluid_dae_integrate: RODAS reported " + std::to_string(S.out.size()) +
                           " of the " + std::to_string(grid.size()) + " requested output points");
    S.out.back().assign(y.begin(), y.end());
    return S.out;
}

/**
 * The multipliers that hold the active caps at this state, by small Newton.
 *
 * Each active cap contributes one equation, d/dt(A x + As s) = 0, and one unknown
 * -- a fraction for a cap that holds or loses, an admitted flow for one that stages
 * -- so the system is square and small. This is what makes an event RESTART
 * consistent: RODAS needs the algebraic unknowns to satisfy their equations at the
 * initial point of an index-1 DAE. It is also the FEASIBILITY test at an
 * activation: a fraction above one means the cap would have to admit more than
 * arrives, so it is not binding after all.
 */
inline std::vector<double> fluid_dae_hold_multipliers(
    const FluidMomentTerms& terms, const FluidDaeGates& gates, const FluidDaeStaging& stg,
    const FluidDaeConstraints& con, const std::vector<std::size_t>& active,
    const std::vector<double>& x, const std::vector<double>& sg, const FluidClosure& cl,
    const std::vector<double>& m0, bool& ok) {
    std::vector<double> m = m0;
    ok = true;
    if (active.empty()) return m;
    const std::size_t n = terms.nstate;
    auto resid = [&](const std::vector<double>& mm) {
        const std::vector<double> r = fluid_moment_rates(terms, x, cl);
        const FluidDaeLegs lg = fluid_dae_legs(terms, gates, stg, con, active, sg, r, mm, true);
        std::vector<double> dx(n, 0.0);
        for (std::size_t s = 0; s < n; ++s) {
            double acc = 0.0;
            for (std::size_t e = 0; e < r.size(); ++e)
                acc += gates.Dn(s, e) * lg.rup[e] + gates.DnExt(s, e) * lg.rin[e] +
                       gates.Dp(s, e) * lg.rin[e];
            dx[s] = acc;
        }
        std::vector<double> F(active.size(), 0.0);
        for (std::size_t k = 0; k < active.size(); ++k) {
            double acc = 0.0;
            for (std::size_t s = 0; s < n; ++s) acc += con.A(active[k], s) * dx[s];
            for (std::size_t j = 0; j < stg.n; ++j) acc += con.As(active[k], j) * lg.ds[j];
            F[k] = acc;
        }
        return F;
    };
    auto inf_norm = [](const std::vector<double>& v) {
        double a = 0.0;
        for (double e : v) a = std::max(a, std::fabs(e));
        return a;
    };
    std::vector<double> F = resid(m);
    for (int it = 0; it < 40; ++it) {
        if (inf_norm(F) < std::max(1e-12, 1e-10 * (inf_norm(m) + 1.0))) break;
        Matrix<double> J(F.size(), m.size(), 0.0);
        for (std::size_t j = 0; j < m.size(); ++j) {
            const double h = std::max(1e-7 * std::fabs(m[j]), 1e-9);
            std::vector<double> mp = m;
            mp[j] += h;
            const std::vector<double> Fp = resid(mp);
            for (std::size_t i = 0; i < F.size(); ++i) J(i, j) = (Fp[i] - F[i]) / h;
        }
        std::vector<double> step;
        try {
            step = lstsq(J, F).x;
        } catch (const std::exception&) {
            break;
        }
        double lam = 1.0;
        bool stepped = false;
        for (int ls = 0; ls < 20; ++ls) {
            std::vector<double> mn(m.size(), 0.0);
            for (std::size_t j = 0; j < m.size(); ++j) mn[j] = std::max(0.0, m[j] - lam * step[j]);
            const std::vector<double> Fn = resid(mn);
            if (inf_norm(Fn) < inf_norm(F)) {
                m = mn;
                F = Fn;
                stepped = true;
                break;
            }
            lam *= 0.5;
        }
        if (!stepped) break;
    }
    ok = inf_norm(F) < 1e-6;
    return m;
}

/** What a hybrid transient did, beside the trajectory. */
struct FluidDaeSwitch {
    double t = 0.0;
    std::size_t row = 0;
    int kind = 0;   ///< 0 release, 1 activate, 2 a crossing the cap could not hold
};

/**
 * The transient UNDER CAPS: one index-1 DAE per segment, restarted at every located
 * crossing. See FluidDaeTransient for why the multiplier is an algebraic unknown and
 * why a staged cap's is a flow rather than a rate.
 */
inline std::vector<std::vector<double> > fluid_dae_integrate_hybrid(
    const FluidMomentTerms& terms, const FluidDaeConservation& cons,
    const FluidDaeConstraints& con, const FluidDaeGates& gates, const FluidDaeStaging& stg,
    const FluidClosure& closure, const std::vector<double>& x0In,
    const std::vector<double>& grid, double tol, std::vector<FluidDaeSwitch>* switches = nullptr) {
    using namespace rodas_impl;
    const std::size_t n = terms.nstate;
    const std::size_t ncon = con.b.size(), nstg = stg.n;

    detail::FluidDaeTransient S;
    S.terms = &terms;
    S.cons = &cons;
    S.closure = closure;
    S.nstate = n;
    S.withcov = false;
    S.cov_idx = terms.cov_idx;
    S.nc = terms.cov_idx.size();
    S.con = &con;
    S.gates = &gates;
    S.stg = &stg;
    S.ncon = ncon;
    S.nstg = nstg;
    S.o_sg = n;
    S.o_m = n + nstg;
    S.o_cov = n + nstg + ncon;
    S.armed.assign(std::max<std::size_t>(ncon, 1), 0);
    S.inert.assign(ncon, 1.0);
    for (std::size_t c = 0; c < ncon; ++c)
        if (con.staged[c]) S.inert[c] = 0.0;   // a flow is inert at zero
    const std::size_t nz = S.o_cov;

    // A STATE ABOVE A CAP IS NOT A STATE THE MODEL CAN BE IN: holding the cap from
    // there would freeze the violation for the whole horizon. Start on the cap, and
    // MOVE the excess rather than dropping it -- conservation is an algebraic row,
    // so an initial state that does not satisfy it is an inconsistent initialisation
    // and no index-1 solver may be handed one.
    std::vector<double> x0 = x0In;
    x0.resize(n, 0.0);
    std::vector<double> sg0(nstg, 0.0), excess(std::max<std::size_t>(ncon, 1), 0.0);
    std::vector<std::size_t> over;
    for (std::size_t c = 0; c < ncon; ++c) {
        double val = 0.0, tot = 0.0;
        for (std::size_t s = 0; s < n; ++s) {
            val += con.A(c, s) * x0[s];
            if (con.A(c, s) > 0.0) tot += x0[s];
        }
        if (val > con.b[c] + std::max(1e-9, tol) && con.b[c] > 0.0) {
            excess[c] = tot * (1.0 - con.b[c] / val);
            for (std::size_t s = 0; s < n; ++s)
                if (con.A(c, s) > 0.0) x0[s] *= con.b[c] / val;
            over.push_back(c);
        }
    }
    // THE INITIAL MODE: a cap the initial state sits ON is already binding, so the
    // first segment must carry it. Feasibility decides.
    std::vector<double> m_init = S.inert;
    if (!over.empty()) {
        std::vector<double> m0(over.size(), 0.0);
        for (std::size_t a = 0; a < over.size(); ++a) m0[a] = con.staged[over[a]] ? 0.0 : 1.0;
        bool ok = false;
        const std::vector<double> mm = fluid_dae_hold_multipliers(terms, gates, stg, con, over,
                                                                  x0, sg0, closure, m0, ok);
        bool feasible = ok;
        for (std::size_t a = 0; a < over.size(); ++a)
            if (!con.staged[over[a]] && mm[a] > 1.0 + 1e-9) feasible = false;
        if (feasible) {
            S.active = over;
            for (std::size_t a = 0; a < over.size(); ++a) m_init[over[a]] = mm[a];
        }
    }
    for (std::size_t oi = 0; oi < over.size(); ++oi) {
        const std::size_t c = over[oi];
        if (excess[c] <= 0.0) continue;
        std::size_t rooms = 0;
        const bool staged_active =
            con.staged[c] && std::find(S.active.begin(), S.active.end(), c) != S.active.end();
        if (staged_active)
            for (std::size_t j = 0; j < nstg; ++j)
                if (stg.gated_by[c][j]) ++rooms;
        if (rooms) {
            for (std::size_t j = 0; j < nstg; ++j)
                if (stg.gated_by[c][j]) sg0[j] += excess[c] / static_cast<double>(rooms);
            continue;
        }
        // otherwise onto the coordinates that FEED the capped stations, which is
        // where a held job waits
        std::vector<bool> pool(n, false);
        bool any_feeder = false;
        std::vector<bool> feeders(n, false);
        for (std::size_t e = 0; e < terms.D.cols(); ++e) {
            if (!gates.gate[c][e]) continue;
            for (std::size_t s = 0; s < n; ++s)
                if (gates.Dn(s, e) < 0.0 || gates.DnExt(s, e) < 0.0) feeders[s] = true;
        }
        for (std::size_t s = 0; s < n; ++s) {
            pool[s] = con.A(c, s) <= 0.0;
            any_feeder = any_feeder || (pool[s] && feeders[s]);
        }
        if (any_feeder)
            for (std::size_t s = 0; s < n; ++s) pool[s] = pool[s] && feeders[s];
        double w = 0.0;
        std::size_t cnt = 0;
        for (std::size_t s = 0; s < n; ++s)
            if (pool[s]) { w += x0[s]; ++cnt; }
        if (w > 1e-8) {
            for (std::size_t s = 0; s < n; ++s)
                if (pool[s]) x0[s] += excess[c] * x0[s] / w;
        } else if (cnt) {
            for (std::size_t s = 0; s < n; ++s)
                if (pool[s]) x0[s] += excess[c] / static_cast<double>(cnt);
        }
    }
    for (std::size_t k = 0; k < S.active.size(); ++k) {
        const std::size_t c = S.active[k];
        if (!con.staged[c]) continue;
        double mass = 0.0;
        for (std::size_t j = 0; j < nstg; ++j)
            if (stg.gated_by[c][j]) mass += sg0[j];
        if (mass > 1e-8) S.armed[c] = 1;
    }

    for (std::size_t c = 0; c < cons.C.rows(); ++c) {
        std::size_t best = n;
        double bestv = -1.0;
        for (std::size_t s = 0; s < n; ++s) {
            if (cons.C(c, s) == 0.0) continue;
            bool taken = false;
            for (std::size_t d = 0; d < S.alg_row.size(); ++d)
                if (S.alg_row[d] == s) taken = true;
            if (taken) continue;
            if (x0[s] > bestv) { bestv = x0[s]; best = s; }
        }
        if (best == n)
            throw NumericError("fluid_dae_integrate_hybrid: a chain has no free coordinate");
        S.alg_row.push_back(best);
    }
    if (grid.empty()) throw InputError("fluid_dae_integrate_hybrid: the output grid is empty");
    S.grid = grid;

    std::vector<double> z(nz, 0.0);
    for (std::size_t s = 0; s < n; ++s) z[s] = x0[s];
    for (std::size_t j = 0; j < nstg; ++j) z[S.o_sg + j] = sg0[j];
    for (std::size_t c = 0; c < ncon; ++c) z[S.o_m + c] = m_init[c];

    double tcur = 0.0;
    const double tend = grid.back();
    const std::size_t seg_max = 4 * ncon + 8;
    for (std::size_t seg = 0; seg < seg_max; ++seg) {
        S.gated_rooms.assign(nstg, 0);
        for (std::size_t k = 0; k < S.active.size(); ++k)
            if (con.staged[S.active[k]])
                for (std::size_t j = 0; j < nstg; ++j)
                    if (stg.gated_by[S.active[k]][j]) S.gated_rooms[j] = 1;
        S.hit = -1;
        S.hit_z.clear();
        S.gprev.clear();
        detail::fluid_dae_active() = &S;

        const integer N = static_cast<integer>(nz);
        const std::size_t lwork = nz * (nz + 1 + nz + 14) + 20 + 32;
        const std::size_t liwork = nz + 20 + 32;
        std::vector<doublereal> y(z.begin(), z.end()), work(lwork, 0.0), rpar(1, 0.0);
        std::vector<integer> iwork(liwork, 0), ipar(1, 0);
        doublereal rtol = tol, atol = tol * 1e-2, x = tcur, xend = tend, h = 1e-6;
        // THE JACOBIAN IS DIFFERENCED, not the drift's: under a cap the rows carry
        // the room and multiplier equations too, and a Jacobian right on one block
        // and zero on the others is worse than none -- RODAS would take it as exact.
        integer itol = 0, ifcn = 0, ijac = 0, mljac = N, mujac = N, idfx = 0;
        integer imas = 1, mlmas = 0, mumas = 0, iout = 1, idid = 0;
        integer lw = static_cast<integer>(lwork), liw = static_cast<integer>(liwork);
        try {
            rodas_(const_cast<integer*>(&N), (U_fp)detail::fluid_dae_fcn, &ifcn, &x, y.data(),
                   &xend, &h, &rtol, &atol, &itol,
                   (U_fp)detail::fluid_dae_jac, &ijac, &mljac, &mujac,
                   (U_fp)detail::fluid_dae_dfx, &idfx,
                   (U_fp)detail::fluid_dae_mas, &imas, &mlmas, &mumas,
                   (U_fp)detail::fluid_dae_solout, &iout,
                   work.data(), &lw, iwork.data(), &liw, rpar.data(), ipar.data(), &idid);
        } catch (...) {
            detail::fluid_dae_active() = nullptr;
            throw;
        }
        detail::fluid_dae_active() = nullptr;
        if (idid != 1 && idid != 2 && S.hit < 0)
            throw NumericError("fluid_dae_integrate_hybrid: RODAS returned idid=" +
                               std::to_string(static_cast<long>(idid)));
        if (S.hit < 0) break;   // the horizon was reached with this set binding

        const std::size_t c = static_cast<std::size_t>(S.hit);
        tcur = S.hit_t;
        z = S.hit_z;
        const std::vector<double> xh(z.begin(), z.begin() + n);
        const std::vector<double> sgh(z.begin() + S.o_sg, z.begin() + S.o_sg + nstg);
        const bool was_active =
            std::find(S.active.begin(), S.active.end(), c) != S.active.end();
        if (was_active) {
            S.active.erase(std::remove(S.active.begin(), S.active.end(), c), S.active.end());
            S.armed[c] = 0;
            z[S.o_m + c] = S.inert[c];
            if (switches) switches->push_back(FluidDaeSwitch{tcur, c, 0});
        } else {
            std::vector<std::size_t> trial = S.active;
            trial.push_back(c);
            std::sort(trial.begin(), trial.end());
            std::vector<double> m0(trial.size(), 0.0);
            for (std::size_t a = 0; a < trial.size(); ++a)
                m0[a] = std::find(S.active.begin(), S.active.end(), trial[a]) != S.active.end()
                            ? z[S.o_m + trial[a]]
                            : (con.staged[trial[a]] ? 0.0 : 1.0);
            bool ok = false;
            const std::vector<double> mm = fluid_dae_hold_multipliers(
                terms, gates, stg, con, trial, xh, sgh, closure, m0, ok);
            bool feasible = ok;
            for (std::size_t a = 0; a < trial.size(); ++a)
                if (!con.staged[trial[a]] && mm[a] > 1.0 + 1e-9) feasible = false;
            if (!feasible) {
                if (switches) switches->push_back(FluidDaeSwitch{tcur, c, 2});
                break;   // a crossing whose cap cannot be held
            }
            S.active = trial;
            for (std::size_t a = 0; a < trial.size(); ++a) z[S.o_m + trial[a]] = mm[a];
            if (switches) switches->push_back(FluidDaeSwitch{tcur, c, 1});
        }
        if (tcur >= tend - 1e-12) break;
    }

    // Whatever the switching did not reach is reported at the last state, so the
    // grid the caller asked for always comes back full.
    while (S.out.size() < grid.size()) S.out.push_back(z);
    std::vector<std::vector<double> > out;
    for (std::size_t i = 0; i < S.out.size(); ++i)
        out.push_back(std::vector<double>(S.out[i].begin(), S.out[i].begin() + n));
    return out;
}

/**
 * The metrics of one state, read exactly as the steady-state table reads them.
 *
 * Shared by the fixed point and by every point of a trajectory so that the two
 * cannot drift: `solver_fluid_transient` reuses `fluid_closing_metrics` for the
 * same reason, and the last point of a long enough run has to reproduce the
 * table or one of them is reading the drift differently.
 */
inline void fluid_dae_metrics(const FluidMomentTerms& terms, const std::vector<double>& x,
                              const FluidClosure& cl, Matrix<double>& QN, Matrix<double>& UN,
                              Matrix<double>& RN, Matrix<double>& TN,
                              const std::vector<double>* rates = nullptr) {
    const std::size_t M = terms.station_block.size();
    const std::size_t K = M ? terms.class_block[0].size() : 0;
    // RATES, when given, is the flow that actually CROSSES each event: a cap that
    // holds the job upstream suppresses that station's departures, and reading the
    // nominal vector here would report a throughput the model does not carry.
    const std::vector<double> r = rates ? *rates : fluid_moment_rates(terms, x, cl);
    const std::vector<double> g = fluid_moment_factors(terms, x, cl);
    QN = Matrix<double>(M, K, 0.0);
    UN = Matrix<double>(M, K, 0.0);
    RN = Matrix<double>(M, K, 0.0);
    TN = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const std::vector<std::size_t>& blk = terms.class_block[i][k];
            if (blk.empty()) continue;
            double q = 0.0, gg = 0.0;
            for (std::size_t s : blk) { q += x[s]; gg += g[s]; }
            QN(i, k) = q;
            UN(i, k) = (terms.sys.sched[i] == lang::SchedStrategy::INF) ? q : gg / terms.S[i];
            double tn = 0.0;
            for (std::size_t e = 0; e < r.size(); ++e)
                if (terms.ev_is_departure[e] && terms.ev_station[e] == i && terms.ev_class[e] == k)
                    tn += r[e];
            TN(i, k) = tn;
            // TN is zero only to the integrator's accuracy: a class that never visits leaves
            // a ~1e-20 residue in TN too, and a strict > 0 test then divides residue by residue.
            if (tn > lang::GlobalConstants::Zero) RN(i, k) = q / tn;
        }
}

namespace detail {

/**
 * The t=0 state of the DAE, in the layout `terms` indexes.
 *
 * `solver_fluid_dae.m:821-828` takes the first row of the SEED TRAJECTORY, with
 * the comment that this is certainly in the closing state layout while
 * `options.init_sol` need not be. The port's seed carries no trajectory -- only
 * its fixed point -- but the two are the same vector by construction:
 * `fluid_dispatch` starts from `opt.init_sol` when the caller gave one and from
 * `fluid_default_initsol` otherwise, so that is what the seed's first row held.
 *
 * `C x0 = N` therefore holds at t=0, which is what makes the algebraic rows
 * consistent from the start and lets RODAS begin without a projection step.
 */
template <class T>
inline std::vector<double> fluid_dae_init_state(const qn::NetworkStruct<T>& sn,
                                                const FluidMomentTerms& terms,
                                                const FluidOptions& opt) {
    std::vector<double> x0 = opt.init_sol.empty()
                                 ? detail::fluid_default_initsol(sn, terms.sys.layout)
                                 : opt.init_sol;
    if (x0.size() != terms.nstate)
        throw InputError("solver_fluid_dae: the initial condition has " +
                         std::to_string(x0.size()) + " entries where the closure state has " +
                         std::to_string(terms.nstate));
    return x0;
}

}  // namespace detail

/**
 * `solver_fluid_dae.m`: the min-normal closure solved as one system.
 *
 * Backs options.method = "dae". The answer is the same closure `minnormal`
 * computes; what differs is that the mean and the variance are solved
 * simultaneously rather than alternated, and that population conservation is an
 * equation rather than a consequence of the drift.
 */
template <class T>
FluidSolution solver_fluid_dae(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                               const FluidDaeOptions& dopt_in = FluidDaeOptions()) {
    const FluidDaeOptions dopt = fluid_dae_options(opt, dopt_in);
    const FluidMomentTerms terms = fluid_moment_terms(sn, opt);
    const std::size_t n = terms.nstate;
    const std::size_t M = terms.station_block.size();
    const std::size_t K = M ? terms.class_block[0].size() : 0;

    if (n > dopt.maxstate)
        throw UnsupportedError(
            "solver_fluid_dae: the dae method solves a " + std::to_string(n) +
            "-unknown algebraic system with a finite-difference Jacobian, above the limit of " +
            std::to_string(dopt.maxstate) +
            " set by options.config.dae_maxstate. Raise it, or use options.method='minnormal' for "
            "the same closure by successive substitution.");
    // DPS and GPS close on the covariance BETWEEN a station's class coordinates,
    // not on the station total, so their closure state is a matrix block rather
    // than the scalar this solves for. `minnormal` carries those blocks through
    // its outer iteration; refusing is better than silently dropping them.
    for (std::size_t i = 0; i < M; ++i)
        if (terms.sys.sched[i] == lang::SchedStrategy::DPS ||
            terms.sys.sched[i] == lang::SchedStrategy::GPS)
            throw UnsupportedError(
                "solver_fluid_dae: the dae method closes on the per-station variance only, and "
                "DPS/GPS close on the covariance between their class coordinates. Use "
                "options.method='minnormal'.");

    // Caps first: staging depends on them, and conservation must count the blocked
    // mass staging holds. GATES decides, per cap and per event, where the mass a
    // cap stops actually goes -- held upstream, lost, or staged.
    FluidDaeConstraints con = fluid_dae_constraints(sn, terms);
    const FluidDaeGates gates = fluid_dae_gates(sn, terms, con);
    const FluidDaeStaging stg = fluid_dae_staging(terms, con);
    fluid_dae_extend(con, stg, terms);
    const FluidDaeConservation cons = fluid_dae_conservation(sn, terms, stg);
    const std::size_t ncon = con.b.size();

    // C*D vanishes on a conserving event set, which is the structural fact that
    // makes the drift Jacobian singular; a leak here would mean the constraint
    // contradicts the drift rather than completing it.
    for (std::size_t c = 0; c < cons.C.rows(); ++c)
        for (std::size_t e = 0; e < terms.D.cols(); ++e) {
            double acc = 0.0;
            for (std::size_t s = 0; s < n; ++s) acc += cons.C(c, s) * terms.D(s, e);
            if (std::fabs(acc) > 1e-7)
                throw NumericError(
                    "solver_fluid_dae: the event set does not conserve a closed chain");
        }

    FluidDaeSystem sysd;
    sysd.terms = &terms;
    sysd.cons = &cons;
    sysd.con = &con;
    sysd.stg = &stg;
    sysd.gates = &gates;
    sysd.cidx = fluid_dae_closable(terms);
    sysd.nstate = n;

    // The seed: Newton needs a point in the basin, not an answer. One
    // first-order solve supplies it, at the cost of the single integration this
    // method exists to avoid repeating twenty times.
    FluidOptions mo = opt;
    mo.method = "closing";
    mo.closure = FluidClosure();
    mo.timespan_end = std::numeric_limits<double>::infinity();
    const FluidSolution seed = detail::fluid_dispatch(sn, mo);
    std::vector<double> xcur(seed.xvec.begin(), seed.xvec.end());
    xcur.resize(n, 0.0);

    // THE VARIANCE IS SEEDED POSITIVE, which is why this route needs no kink
    // probe. sigma2 = 0 is where min(n,c) has no derivative and a saturated
    // model's first-order fixed point sits exactly there. The station mean is an
    // O(N) starting value in the right units and costs no Lyapunov solve.
    std::vector<double> s2seed(sysd.cidx.size(), 0.0);
    for (std::size_t j = 0; j < sysd.cidx.size(); ++j) {
        double acc = 0.0;
        for (std::size_t sIdx : terms.station_block[sysd.cidx[j]]) acc += xcur[sIdx];
        s2seed[j] = std::max(1e-8, acc);
    }

    const double tol = (opt.tol > 0.0 && std::isfinite(opt.tol)) ? opt.tol : 1e-8;
    FluidDaeNewtonInfo info;
    std::vector<double> u, sgv(stg.n, 0.0), mult;

    // ACTIVE SET. A capacity constraint is an inequality, and an inequality has
    // no residual to hand a Newton solver -- only the caps that actually bind
    // become equations. Each pass is a complete simultaneous solve, so the loop
    // iterates over WHICH caps bind, not over the closure.
    //
    // ONE MULTIPLIER PER ACTIVE ROW, not per region. A region's waiting room used
    // to carry a single drain rate, so two caps of one region were two equalities
    // against one control and the case was refused by name. A room gated by
    // several active rows now drains at the harmonic composition of their rates
    // and a held cap composes as a product of fractions.
    //
    // THE SET STARTS EMPTY, and the first pass therefore asks for the UNCONSTRAINED
    // fixed point -- which is what makes the answer of a model whose caps bind
    // independent of how far outside the seed happened to land. It is also not
    // always solvable: an overloaded M/M/1/K has NO equilibrium without its cap, so
    // that pass fails rather than converging, and the caps the iterate violates are
    // seeded into the set instead (below). Doing that up front for every model
    // looked cheaper and changed answers: a region of 8 over two identical queues
    // settled 7.43/0.57 from the clipped seed and 4/4 from the unconstrained one.
    bool seeded_from_failure = false;

    // the last iterate that WAS a fixed point, and the answer read off it: a pass
    // that fails to converge leaves an iterate that is not a fixed point of
    // anything, and reading the next active set off it is how one bad pass turns
    // into a walk through unrelated capacity combinations
    std::vector<double> x_ok = xcur, s2_ok = s2seed;
    bool have_best = false, best_feasible = false;
    std::vector<double> best_x, best_sg, best_s2, best_mult;
    std::vector<std::size_t> best_active;
    FluidDaeNewtonInfo best_info;
    bool clamped = false, best_clamped = false;
    const std::size_t aset_max = std::max<std::size_t>(4, 2 * ncon + 2);
    for (std::size_t aset = 0; aset < aset_max; ++aset) {
        bool has_clamp = false;
        for (std::size_t k = 0; k < sysd.active.size(); ++k)
            has_clamp = has_clamp || !con.staged[sysd.active[k]];

        // Seed each waiting room with the mass that does not fit, and every
        // multiplier at unity -- an unthrottled fraction for a held or lost cap, and
        // the drain rate the region route has always started from.
        std::vector<double> sg0(stg.n, 0.0);
        for (std::size_t k = 0; k < sysd.active.size(); ++k) {
            const std::size_t c = sysd.active[k];
            if (!con.staged[c]) continue;
            double cur = 0.0;
            for (std::size_t sIdx = 0; sIdx < n; ++sIdx) cur += con.A(c, sIdx) * xcur[sIdx];
            const double excess = std::max(0.0, cur - con.b[c]);
            std::size_t cnt = 0;
            for (std::size_t j = 0; j < stg.n; ++j) if (stg.gated_by[c][j]) ++cnt;
            if (cnt)
                for (std::size_t j = 0; j < stg.n; ++j)
                    if (stg.gated_by[c][j]) sg0[j] = excess / static_cast<double>(cnt);
        }
        std::vector<double> u0;
        u0.assign(xcur.begin(), xcur.end());
        u0.insert(u0.end(), sg0.begin(), sg0.end());
        u0.insert(u0.end(), s2seed.begin(), s2seed.end());
        u0.insert(u0.end(), sysd.active.size(), 1.0);

        // THE CLAMPED COVARIANCE IS A FALLBACK, NOT THE DEFAULT. A cap that holds or
        // loses fixes its own combination of the state, so the honest linear-noise
        // approximation puts no fluctuation there -- but the truth is neither that
        // nor the unprojected variance: the population under a cap follows a
        // TRUNCATED distribution. The unprojected solve is what every other fluid
        // method computes, so it runs first; the projection is tried only when the
        // unprojected system has no stationary covariance at all, which is the
        // neutral case an overloaded loss station produces. Deciding once per PASS
        // matters: a projection switching between iterates would give Newton a
        // discontinuous system.
        bool solved = false;
        bool nohyp = false;
        // WHETHER THE FAILURE WAS THE NON-HYPERBOLIC ONE, kept apart from the
        // message. Rethrowing every Newton failure as a plain NumericError below
        // erases exactly the distinction FluidNonHyperbolicError exists to carry:
        // the runner's fallback ladder catches THAT type and only that type, so a
        // model whose fixed point is merely neutral would reach the caller as a
        // hard error instead of walking to the next rung.
        bool nohyp_typed = false;
        std::string nohyp_what;
        for (int attempt = 0; attempt < (has_clamp ? 2 : 1); ++attempt) {
            sysd.clampT = attempt == 1 ? fluid_dae_clamp_tangent(con, sysd.active, terms)
                                       : Matrix<double>(0, 0, 0.0);
            u = u0;
            try {
                info = fluid_dae_newton(sysd, u, tol, dopt.newton_max);
            } catch (const FluidNonHyperbolicError& e) {
                if (attempt == 1 || !has_clamp) {
                    nohyp = true;
                    nohyp_typed = true;
                    nohyp_what = e.what();
                    break;
                }
                continue;
            } catch (const std::exception& e) {
                if (attempt == 1 || !has_clamp) {
                    nohyp = true;
                    nohyp_what = e.what();
                    break;
                }
                continue;
            }
            clamped = attempt == 1;
            solved = true;
            if (info.converged) break;
        }
        if (nohyp) {
            // THE UNCONSTRAINED FIXED POINT NEED NOT EXIST. An overloaded open
            // station has no equilibrium until its buffer bounds it, so the pass
            // that asks for one fails and the caps the iterate violates -- or left
            // infinite, which no comparison catches -- are seeded into the active
            // set instead. Once: a second failure with caps already bound is the
            // model's answer and not a starting point to improve.
            std::vector<std::size_t> cand;
            for (std::size_t c = 0; c < ncon; ++c) {
                if (std::find(sysd.active.begin(), sysd.active.end(), c) != sysd.active.end())
                    continue;
                double acc = 0.0;
                bool bad = false;
                for (std::size_t s = 0; s < n; ++s) {
                    if (con.A(c, s) > 0.0 && !std::isfinite(xcur[s])) bad = true;
                    acc += con.A(c, s) * xcur[s];
                }
                if (bad || acc > con.b[c] + std::max(1e-9, tol)) cand.push_back(c);
            }
            if (seeded_from_failure || cand.empty()) {
                if (nohyp_typed) throw FluidNonHyperbolicError(nohyp_what);
                throw NumericError(nohyp_what);
            }
            for (std::size_t ci = 0; ci < cand.size(); ++ci) {
                const std::size_t c = cand[ci];
                double val = 0.0;
                std::size_t cnt = 0;
                for (std::size_t s = 0; s < n; ++s) {
                    val += con.A(c, s) * xcur[s];
                    if (con.A(c, s) > 0.0) ++cnt;
                }
                if (!cnt || !(con.b[c] > 0.0)) continue;
                if (std::isfinite(val) && val > con.b[c]) {
                    const double f = con.b[c] / val;
                    for (std::size_t s = 0; s < n; ++s)
                        if (con.A(c, s) > 0.0) xcur[s] *= f;
                } else if (!std::isfinite(val)) {
                    for (std::size_t s = 0; s < n; ++s)
                        if (con.A(c, s) > 0.0) xcur[s] = con.b[c] / static_cast<double>(cnt);
                }
            }
            for (std::size_t s = 0; s < n; ++s)
                if (!std::isfinite(xcur[s])) xcur[s] = 0.0;
            for (std::size_t ci = 0; ci < cand.size(); ++ci) sysd.active.push_back(cand[ci]);
            std::sort(sysd.active.begin(), sysd.active.end());
            sysd.active.erase(std::unique(sysd.active.begin(), sysd.active.end()),
                              sysd.active.end());
            seeded_from_failure = true;
            continue;
        }
        if (!solved) throw NumericError("solver_fluid_dae: the closure could not be evaluated");
        xcur.assign(u.begin(), u.begin() + n);
        sgv.assign(u.begin() + n, u.begin() + n + stg.n);
        for (std::size_t j = 0; j < sysd.cidx.size(); ++j)
            s2seed[j] = std::max(0.0, u[n + stg.n + j]);
        mult.assign(u.begin() + n + stg.n + sysd.cidx.size(), u.end());
        if (ncon == 0) break;

        std::vector<double> slack(ncon, 0.0);
        bool feasible = true;
        for (std::size_t c = 0; c < ncon; ++c) {
            double acc = 0.0;
            for (std::size_t sIdx = 0; sIdx < n; ++sIdx) acc += con.A(c, sIdx) * xcur[sIdx];
            for (std::size_t j = 0; j < stg.n; ++j) acc += con.As(c, j) * sgv[j];
            slack[c] = con.b[c] - acc;
            if (slack[c] < -std::max(1e-9, tol)) feasible = false;
        }
        if (info.converged) {
            x_ok = xcur;
            s2_ok = s2seed;
            have_best = true;
            best_x = xcur; best_sg = sgv; best_s2 = s2seed; best_mult = mult;
            best_active = sysd.active; best_info = info; best_clamped = clamped;
            best_feasible = feasible;
        }
        std::vector<std::size_t> violated;
        for (std::size_t c = 0; c < ncon; ++c) {
            if (std::find(sysd.active.begin(), sysd.active.end(), c) != sysd.active.end()) continue;
            if (slack[c] < -std::max(1e-9, tol)) violated.push_back(c);
        }
        // THE RELEASE SIGNAL IS THE MULTIPLIER'S OWN UNITS, and the two kinds do not
        // share them. A held or lost cap throttles by a FRACTION, so one that came
        // back above one was holding the flow down for no reason. A staged cap
        // throttles by a RATE, which has no such scale -- there the signal is a
        // waiting room with no blocked mass at all.
        std::vector<std::size_t> released;
        for (std::size_t k = 0; k < sysd.active.size(); ++k) {
            const std::size_t c = sysd.active[k];
            if (con.staged[c]) {
                double held = 0.0;
                std::size_t rooms = 0;
                for (std::size_t j = 0; j < stg.n; ++j)
                    if (stg.gated_by[c][j]) { held += sgv[j]; ++rooms; }
                if (!rooms || held < 1e-9) released.push_back(c);
            } else if (mult[k] > 1.0 + std::max(1e-9, tol)) {
                released.push_back(c);
            }
        }
        if (!info.converged) {
            // A FAILED PASS SAYS NOTHING ABOUT WHICH CAPS BIND. Its release signal is
            // still information, but its state is not, so no row is ADDED from it and
            // the next pass restarts from the last point that was a fixed point.
            violated.clear();
            xcur = x_ok;
            s2seed = s2_ok;
            if (released.empty()) break;
        }
        if (violated.empty() && released.empty()) break;
        std::vector<std::size_t> next;
        for (std::size_t a : sysd.active)
            if (std::find(released.begin(), released.end(), a) == released.end()) next.push_back(a);
        for (std::size_t v : violated) next.push_back(v);
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        sysd.active = next;
    }
    // A CONVERGED POINT BEATS THE LAST ITERATE. The loop can end on a pass that did
    // not converge -- a cap the closure cannot hold at any multiplier is added,
    // fails and is released for as long as the loop runs.
    if (have_best && !info.converged) {
        xcur = best_x; sgv = best_sg; s2seed = best_s2; mult = best_mult;
        sysd.active = best_active; info = best_info; clamped = best_clamped;
        // THIS PORT THROWS WHERE THE REFERENCE WARNS, and only because it has no
        // warning channel: MATLAB, the JAR and native Python report the converged
        // point beside a warning naming the caps it exceeds. Returning a
        // cap-violating point silently is the one option none of them takes.
        if (!best_feasible)
            throw NumericError(
                "solver_fluid_dae: no fixed point of the closure satisfies every cap. The closure "
                "wants more jobs there than the cap allows and no admission multiplier holds it. "
                "Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES for this model.");
    }
    sysd.clampT = clamped ? fluid_dae_clamp_tangent(con, sysd.active, terms)
                          : Matrix<double>(0, 0, 0.0);

    std::vector<double> x = xcur;
    FluidClosure cl;
    cl.sigma2.assign(M, 0.0);
    cl.cov.assign(M, Matrix<double>(0, 0, 0.0));
    for (std::size_t j = 0; j < sysd.cidx.size(); ++j)
        cl.sigma2[sysd.cidx[j]] = std::max(0.0, u[n + stg.n + j]);

    // The transient is the same closure integrated as an index-1 DAE, held at
    // the variance the steady state converged to -- so the trajectory and the
    // table are read off one drift rather than two.
    //
    // THE INITIAL CONDITION IS THE MODEL'S, NOT THE SEED'S ANSWER. `seed` is the
    // first-order solve run to its fixed point, so `seed.xvec` is a STEADY
    // STATE; starting the integration there makes every horizon report the
    // answer it already had and hides the trajectory entirely. The reference
    // takes `xfall(1,:)`, the first row of that seed's own trajectory, which is
    // the t=0 state -- and that is `fluid_default_initsol`, the initial
    // condition the seed run itself started from, in the layout `terms` indexes.
    if (std::isfinite(opt.timespan_end) && opt.timespan_end > 0.0) {
        const std::vector<double> x0 = detail::fluid_dae_init_state(sn, terms, opt);
        const std::vector<double> grid(1, opt.timespan_end);
        // Held at the converged variance HERE, and only here: this call wants
        // the state at one horizon to read a table off, not a path, so the
        // covariance rows would be integrated and thrown away.
        //
        // UNDER A CAP THIS IS A HYBRID DAE, integrated segment by segment with the
        // binding set updated at each located crossing: what the steady state
        // settles once with an active-set loop, the trajectory settles again at
        // every fill and every drain.
        const std::vector<double> zend =
            con.empty()
                ? fluid_dae_integrate(terms, cons, cl, x0, grid, tol).back()
                : fluid_dae_integrate_hybrid(terms, cons, con, gates, stg, cl, x0, grid, tol)
                      .back();
        x.assign(zend.begin(), zend.begin() + n);
    }

    // THE FIRING RATES, and the same covariance treatment the converged pass used:
    // reading the nominal vector here would count admissions a held cap suppressed,
    // and dropping the clamp would ask for a covariance the constrained fixed point
    // does not have (an overloaded loss station is neutral along its cap).
    const std::vector<double> rnom = fluid_moment_rates(terms, x, cl);
    const FluidDaeLegs lgf =
        fluid_dae_legs(terms, gates, stg, con, sysd.active, sgv, rnom, mult, false);
    const std::vector<double> r = lgf.rin;
    const Matrix<double> A = fluid_drift_jacobian(terms, x, cl);
    const Matrix<double>* cTf = sysd.clampT.rows() ? &sysd.clampT : nullptr;
    const Matrix<double> Sigma = fluid_moment_lyapunov(terms, A, lgf.rup, cTf);

    FluidSolution out;
    fluid_dae_metrics(terms, x, cl, out.QN, out.UN, out.RN, out.TN, &r);

    // UTILIZATION MUST BE READ FROM THE FLOW THAT ACTUALLY CROSSES once a cap
    // binds. Away from a constraint the in-service fluid sum(g)/s and the carried
    // utilization T/(mu s) are the same number, because the drift balances; under
    // an ACTIVE cap they are not -- the multiplier throttles the departures (TN)
    // and leaves the in-service fluid alone, so the two columns disagreed: a closed
    // tandem capped at 1 reported Util 0.688 at a station whose own Tput/mu was
    // 0.550, and the exact answer is neither. Every other solver reports the
    // CARRIED utilization (Util = X*D, the utilization law), and SolverCTMC gives
    // 0.444 for the same station, so that is the convention the throttled point has
    // to keep. It is applied HERE and not inside fluid_dae_metrics because that
    // helper is shared with every point of the trajectory, whose active set is not
    // this one. Unconstrained runs are bit-identical: the loop is entered only when
    // the active set is non-empty. see _kb/06-solver-catalog.md
    if (!sysd.active.empty()) {
        for (std::size_t i = 0; i < M; ++i) {
            if (terms.sys.sched[i] == lang::SchedStrategy::INF || terms.is_ext[i]) continue;
            for (std::size_t k = 0; k < K; ++k) {
                if (terms.class_block[i][k].empty()) continue;
                const double mu = num_traits<T>::to_double(sn.rates(i, k));
                if (std::isfinite(mu) && mu > 0.0)
                    out.UN(i, k) = out.TN(i, k) / (mu * terms.S[i]);
            }
        }
    }
    out.CN.assign(K, 0.0);
    out.XN.assign(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t rs = (k < sn.classes.size()) ? sn.classes[k].refstat : 0;
        if (rs >= 1 && rs <= M) out.XN[k] = out.TN(rs - 1, k);
        double q = 0.0;
        for (std::size_t i = 0; i < M; ++i) q += out.QN(i, k);
        if (out.XN[k] > 0.0) out.CN[k] = q / out.XN[k];
    }

    out.xvec = x;
    out.iters = info.iters + seed.iters;
    out.method = "dae";
    out.closure = cl;
    out.has_moments = true;
    out.moments.Sigma = Sigma;
    out.moments.sigma2 = cl.sigma2;
    out.moments.outer_iters = info.iters;
    out.moments.class_block = terms.class_block;
    out.moments.QVar = Matrix<double>(M, K, 0.0);
    out.moments.QStd = Matrix<double>(M, K, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const std::vector<std::size_t>& blk = terms.class_block[i][k];
            double acc = 0.0;
            for (std::size_t a = 0; a < blk.size(); ++a)
                for (std::size_t b = 0; b < blk.size(); ++b) acc += Sigma(blk[a], blk[b]);
            out.moments.QVar(i, k) = std::max(0.0, acc);
            out.moments.QStd(i, k) = std::sqrt(out.moments.QVar(i, k));
        }
    return out;
}

/**
 * `@@SolverFLD/getTranAvg` for the DAE route: the metrics ALONG the trajectory.
 *
 * The counterpart of `solver_fluid_transient` (solver_fluid.h), and different
 * from it in exactly the way this method is different. That one forces the
 * method to `closing` -- the reference does the same, because matrix and the
 * smoothed variants are steady-state devices -- and integrates y' = f with
 * LSODA, so population conservation holds only to integrator tolerance. This
 * integrates M y' = f with a SINGULAR M, so conservation is an algebraic
 * equation satisfied at every reported point rather than a quantity that drifts.
 *
 * THE VARIANCE IS HELD AT ITS STATIONARY VALUE, which is the narrowing the port
 * makes against `solver_fluid_dae.m`. The reference integrates the covariance
 * alongside the mean when the closable state is small enough (nc <= dae_maxcov)
 * and holds it above that; this always holds it. What is held is the variance
 * the steady-state Newton converged to, so the trajectory and the table are read
 * off ONE drift -- and a held variance is what `minnormal` uses for the whole of
 * its own transient anyway, so this is the reference's own fallback, not a
 * different closure.
 *
 * The steady state is solved FIRST and is not optional: it is where that
 * variance comes from. A caller that wants only the trajectory still pays for
 * the fixed point.
 */
template <class T>
std::vector<FluidTranPoint> solver_fluid_dae_transient(
    const qn::NetworkStruct<T>& sn, const FluidOptions& opt, double t_end,
    std::size_t points = 101, const std::vector<double>& out_grid = std::vector<double>(),
    const FluidDaeOptions& dopt_in = FluidDaeOptions()) {
    if (!(t_end > 0.0)) throw InputError("solver_fluid_dae_transient: t_end must be positive");
    if (points < 2) throw InputError("solver_fluid_dae_transient: need at least two output points");
    const FluidDaeOptions dopt = fluid_dae_options(opt, dopt_in);

    FluidOptions os = opt;
    os.timespan_end = std::numeric_limits<double>::infinity();
    const FluidMomentTerms terms = fluid_moment_terms(sn, os);
    FluidDaeConstraints con = fluid_dae_constraints(sn, terms);
    const FluidDaeGates gates = fluid_dae_gates(sn, terms, con);
    const FluidDaeStaging stg = fluid_dae_staging(terms, con);
    fluid_dae_extend(con, stg, terms);
    const FluidDaeConservation cons = fluid_dae_conservation(sn, terms, stg);

    // The steady state, for the variance the trajectory STARTS from -- and, above
    // `maxcov`, is held at. Asked for over an infinite horizon so the solve is
    // the algebraic one and this does not recurse into an integration.
    const FluidSolution ss = solver_fluid_dae(sn, os, dopt);

    // THE COVARIANCE IS INTEGRATED ALONGSIDE THE MEAN where it is small enough
    // to afford, which is what makes this and `kp` the only fluid methods with a
    // time-varying second moment; `minnormal` evaluates its whole transient at
    // the single STATIONARY variance. The cost is the reason for the cap: nc^2
    // extra differential states differenced numerically is nc^4 work. Above it
    // the mean is still integrated as a DAE -- conservation stays an equation --
    // and the variance falls back to the stationary one, which is the
    // reference's own fallback rather than a different closure.
    const std::vector<std::size_t> closable = fluid_dae_closable(terms);
    const std::size_t nc = terms.cov_idx.size();
    const bool withcov = nc > 0 && nc <= dopt.maxcov;

    std::vector<double> grid = out_grid;
    if (grid.empty()) {
        grid.resize(points);
        for (std::size_t j = 0; j < points; ++j)
            grid[j] = t_end * static_cast<double>(j) / static_cast<double>(points - 1);
    }
    // t=0 is the initial condition, which no integrator has to be asked for.
    const bool has_zero = !grid.empty() && grid.front() <= 0.0;
    std::vector<double> inner(grid.begin() + (has_zero ? 1 : 0), grid.end());

    const double tol = (opt.tol > 0.0 && std::isfinite(opt.tol)) ? opt.tol : 1e-8;
    const std::vector<double> x0 = detail::fluid_dae_init_state(sn, terms, opt);
    const std::size_t nz = terms.nstate + (withcov ? nc * nc : 0);
    std::vector<std::vector<double> > xs;
    if (has_zero) {
        std::vector<double> z0(nz, 0.0);
        for (std::size_t i = 0; i < terms.nstate; ++i) z0[i] = x0[i];
        xs.push_back(z0);
    }
    if (!inner.empty()) {
        // UNDER A CAP THIS IS A HYBRID DAE, integrated segment by segment with the
        // binding set updated at each located crossing. The covariance is held there
        // rather than integrated: a segment restart would have to carry it across
        // the switch, and what the reference holds above `maxcov` it holds here for
        // the whole capped run.
        const std::vector<std::vector<double> > got =
            con.empty()
                ? fluid_dae_integrate(terms, cons, ss.closure, x0, inner, tol, withcov, closable)
                : fluid_dae_integrate_hybrid(terms, cons, con, gates, stg, ss.closure, x0, inner,
                                             tol);
        xs.insert(xs.end(), got.begin(), got.end());
    }

    std::vector<FluidTranPoint> out;
    out.reserve(xs.size());
    const std::size_t M = terms.station_block.size();
    const std::size_t K = M ? terms.class_block[0].size() : 0;
    for (std::size_t j = 0; j < xs.size(); ++j) {
        std::vector<double> x(xs[j].begin(), xs[j].begin() + terms.nstate);
        // The interpolant is a polynomial and does not know the state is a
        // population; a point that undershoots zero between two steps is
        // rounding, not a negative queue.
        for (double& v : x)
            if (v < 0.0) v = 0.0;
        // The rate factors are read at the variance the trajectory HAD at this
        // instant, not at the stationary one, whenever the covariance was
        // carried; that is the whole difference from `minnormal`'s transient.
        FluidClosure cl = ss.closure;
        Matrix<double> Sigma(terms.nstate, terms.nstate, 0.0);
        if (withcov) {
            for (std::size_t a = 0; a < nc; ++a)
                for (std::size_t b = 0; b < nc; ++b)
                    Sigma(terms.cov_idx[a], terms.cov_idx[b]) =
                        0.5 * (xs[j][terms.nstate + a * nc + b]
                               + xs[j][terms.nstate + b * nc + a]);
            cl.sigma2.assign(M, 0.0);
            cl.cov.assign(M, Matrix<double>(0, 0, 0.0));
            const std::vector<double> s2 = fluid_dae_sigma_from(Sigma, terms, closable);
            for (std::size_t i = 0; i < s2.size(); ++i) cl.sigma2[i] = s2[i];
        }
        FluidTranPoint pt;
        pt.t = grid[j];
        Matrix<double> R;
        fluid_dae_metrics(terms, x, cl, pt.QN, pt.UN, R, pt.TN);
        detail::fluid_snap_all(pt.QN, pt.UN, R, pt.TN);
        if (withcov) {
            pt.QVar = Matrix<double>(M, K, 0.0);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t k = 0; k < K; ++k) {
                    const std::vector<std::size_t>& blk = terms.class_block[i][k];
                    double acc = 0.0;
                    for (std::size_t a = 0; a < blk.size(); ++a)
                        for (std::size_t b = 0; b < blk.size(); ++b)
                            acc += Sigma(blk[a], blk[b]);
                    pt.QVar(i, k) = std::max(0.0, acc);
                }
        }
        out.push_back(pt);
    }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_DAE_H
