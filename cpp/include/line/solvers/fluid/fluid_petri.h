/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_PETRI_H
#define LINE_SOLVERS_FLUID_PETRI_H

/**
 * Fluid analysis of a stochastic Petri net: one simultaneous algebraic solve per
 * active set.
 *
 * Port of `matlab/src/solvers/FLD/solver_fluid_petri.m`, together with
 * `fluid_petri_conservation.m`, `fluid_petri_constraints.m`,
 * `fluid_petri_immediate.m` and `fluid_petri_applicable.m`. Cross-checked
 * against `jar/src/main/java/jline/solvers/fluid/petri/PetriSolver.java` and the
 * python `solver_fld/methods/petri.py`.
 *
 * A net's fluid limit is NOT the queueing drift with places in place of
 * stations. Its immediate transitions have no rate at all: their limit is a
 * FLOW, an algebraic unknown pinned by the constraint that the input place they
 * bind holds no mass. The marking, the diffusion covariance, those flows, the
 * multi-server phase latches and the capacity gates therefore solve
 * SIMULTANEOUSLY, as one algebraic system per active set, rather than by
 * integrating an ODE to its fixed point.
 *
 * THE ACTIVE SET IS WHAT ITERATES. A pass solves the system for a fixed choice
 * of which immediate modes fire, which coordinate each pins, and which
 * capacities bind; the answer says whether that choice was right (a negative
 * flow, a negative marking, a violated or a released capacity), and the next
 * pass makes ONE move. The moves are ordered so that a failure of each
 * invalidates the next.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "line/solvers/fluid/fluid_moments.h"
#include "line/solvers/fluid/fluid_petri_system.h"
#include "line/solvers/fluid/fluid_petri_terms.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"
#include "line/util/ode.h"
#include "line/util/svd.h"

namespace line {
namespace fluid {
namespace petri {

// ============================ conservation ===============================

/** The conserved quantities of a net, as equations. */
struct PetriConservation {
    Matrix<double> C;
    std::vector<double> N;
    double leak = 0.0;
    std::vector<std::string> label;
};

namespace cons_detail {

/**
 * A rational basis of the null space, as MATLAB's `null(A,'r')` returns.
 *
 * D is INTEGRAL -- arc multiplicities and unit phase moves -- so the reduced row
 * echelon form gives exact rational rows, which keeps each conservation row
 * readable as a statement about named places instead of an arbitrary orthogonal
 * mixture. An SVD basis would span the same space and say nothing.
 */
inline Matrix<double> null_rational(const Matrix<double>& A) {
    const std::size_t rows = A.rows(), cols = A.cols();
    if (rows == 0 || cols == 0) return Matrix<double>(cols, 0, 0.0);
    Matrix<double> R = A;
    std::vector<std::size_t> piv;
    std::size_t r = 0;
    const double tol = 1e-12;
    for (std::size_t c = 0; c < cols && r < rows; ++c) {
        std::size_t k = r;
        double best = std::fabs(R(r, c));
        for (std::size_t i = r + 1; i < rows; ++i)
            if (std::fabs(R(i, c)) > best) {
                best = std::fabs(R(i, c));
                k = i;
            }
        if (best <= tol) {
            for (std::size_t i = r; i < rows; ++i) R(i, c) = 0.0;
            continue;
        }
        if (k != r)
            for (std::size_t j = 0; j < cols; ++j) std::swap(R(r, j), R(k, j));
        const double p = R(r, c);
        for (std::size_t j = 0; j < cols; ++j) R(r, j) /= p;
        for (std::size_t i = 0; i < rows; ++i) {
            if (i == r) continue;
            const double f = R(i, c);
            if (f == 0.0) continue;
            for (std::size_t j = 0; j < cols; ++j) R(i, j) -= f * R(r, j);
        }
        piv.push_back(c);
        ++r;
    }
    std::vector<std::size_t> free;
    for (std::size_t c = 0; c < cols; ++c)
        if (std::find(piv.begin(), piv.end(), c) == piv.end()) free.push_back(c);
    Matrix<double> Z(cols, free.size(), 0.0);
    for (std::size_t i = 0; i < free.size(); ++i) {
        Z(free[i], i) = 1.0;
        for (std::size_t rr = 0; rr < piv.size(); ++rr) Z(piv[rr], i) = -R(rr, free[i]);
    }
    return Z;
}

/** %g formatting, so an integral weight prints without a decimal tail. */
inline std::string trim(double v) {
    if (v == std::rint(v) && std::fabs(v) < 1e15)
        return std::to_string(static_cast<long long>(v));
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%g", v);
    return std::string(buf);
}

/** One conservation row, written out over the coordinates it touches. */
inline std::string row_label(const PetriTerms& t, const Matrix<double>& C, std::size_t r) {
    std::string out;
    for (std::size_t s = 0; s < C.cols(); ++s) {
        const double w = C(r, s);
        if (w == 0.0) continue;
        std::string nm;
        if (s < t.nm) {
            nm = t.names_node[t.coord_node[s]] + "(class " +
                 std::to_string(t.coord_class[s] + 1) + ")";
        } else {
            nm = "phase";
            for (std::size_t j = 0; j < t.modes.size(); ++j)
                for (std::size_t q = 0; q < t.modes[j].zblk.size(); ++q)
                    if (t.modes[j].zblk[q] == s)
                        nm = t.modes[j].label + " phase " + std::to_string(q + 1);
        }
        if (!out.empty()) out += " + ";
        out += (w == 1.0) ? nm : (trim(w) + "*" + nm);
    }
    return out;
}

}  // namespace cons_detail

/**
 * The conserved quantities, as equations: `u'D = 0` implies `u'x` is constant.
 *
 * On the marking coordinates those u are the net's P-invariants; on a mode's
 * phase block the all-ones vector is one of them, which is the statement that
 * the phase coordinates are a distribution. Both come out of the SAME null
 * space, so the phase normalisation needs no separate row.
 *
 * AN OPEN NET LOSES THE ROWS ITS ARRIVALS BREAK, automatically: the arrival
 * columns are part of D, so a u an arrival moves is not in the null space.
 */
inline PetriConservation petri_conservation(const PetriTerms& t) {
    PetriConservation cons;
    if (t.D.rows() == 0 || t.D.cols() == 0 || t.nev == 0) {
        cons.C = Matrix<double>(0, t.nstate, 0.0);
        return cons;
    }
    Matrix<double> Dt(t.D.cols(), t.D.rows(), 0.0);
    for (std::size_t i = 0; i < t.D.rows(); ++i)
        for (std::size_t j = 0; j < t.D.cols(); ++j) Dt(j, i) = t.D(i, j);
    const Matrix<double> Z = cons_detail::null_rational(Dt);
    Matrix<double> C(Z.cols(), Z.rows(), 0.0);
    for (std::size_t i = 0; i < Z.rows(); ++i)
        for (std::size_t j = 0; j < Z.cols(); ++j) C(j, i) = Z(i, j);
    // A rational basis has exact zeros; anything below this is a rounding
    // artefact of the elimination, not a coefficient.
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < C.cols(); ++j)
            if (std::fabs(C(i, j)) < 1e-12) C(i, j) = 0.0;
    // Scale each row by its smallest nonzero magnitude, so a row reads as a
    // statement about named places with small integer weights.
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < C.rows(); ++i) {
        double smallest = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < C.cols(); ++j)
            if (C(i, j) != 0.0) smallest = std::min(smallest, std::fabs(C(i, j)));
        if (!std::isfinite(smallest)) continue;
        for (std::size_t j = 0; j < C.cols(); ++j) C(i, j) /= smallest;
        keep.push_back(i);
    }
    cons.C = Matrix<double>(keep.size(), t.nstate, 0.0);
    for (std::size_t r = 0; r < keep.size(); ++r)
        for (std::size_t j = 0; j < t.nstate; ++j) cons.C(r, j) = C(keep[r], j);
    cons.N.assign(cons.C.rows(), 0.0);
    for (std::size_t r = 0; r < cons.C.rows(); ++r) {
        double v = 0.0;
        for (std::size_t j = 0; j < t.nstate; ++j) v += cons.C(r, j) * t.x0[j];
        cons.N[r] = v;
    }
    for (std::size_t r = 0; r < cons.C.rows(); ++r)
        for (std::size_t c = 0; c < t.nev; ++c) {
            double v = 0.0;
            for (std::size_t j = 0; j < t.nstate; ++j) v += cons.C(r, j) * t.D(j, c);
            cons.leak = std::max(cons.leak, std::fabs(v));
        }
    for (std::size_t r = 0; r < cons.C.rows(); ++r)
        cons.label.push_back(cons_detail::row_label(t, cons.C, r));
    return cons;
}

// ============================ capacities =================================

/** Every finite place capacity as a linear row `A x <= b`. */
struct PetriConstraints {
    Matrix<double> A;
    std::vector<double> b;
    std::vector<std::string> label;
    std::vector<std::vector<bool>> cover;
};

/**
 * THE GATE IS A LOSS ON THE DEPOSIT: LINE loses the tokens a firing would push
 * past a place's capacity, so the fluid analogue scales the DEPOSIT leg of every
 * event adding mass to the capped place and leaves the removal leg alone. The
 * rows here name which coordinates are capped; the scaling lives in the tangent
 * clamp.
 */
template <class T>
inline PetriConstraints petri_constraints(const qn::NetworkStruct<T>& sn, const PetriTerms& t) {
    std::vector<std::vector<double>> rows;
    std::vector<double> bs;
    std::vector<std::string> labels;
    for (std::size_t pi = 0; pi < t.places.size(); ++pi) {
        const std::size_t ind = t.places[pi];
        const std::size_t ist = sn.nodes[ind - 1].station;
        std::vector<std::size_t> slots;
        for (std::size_t k = 0; k < t.K; ++k)
            if (t.pidx[ind - 1][k] >= 0)
                slots.push_back(static_cast<std::size_t>(t.pidx[ind - 1][k]));
        if (slots.empty()) continue;
        if (ist >= 1 && ist <= sn.stations.size()) {
            const double cap = sn.stations[ist - 1].cap;
            if (std::isfinite(cap)) {
                std::vector<double> row(t.nstate, 0.0);
                for (std::size_t s : slots) row[s] = 1.0;
                rows.push_back(row);
                bs.push_back(cap);
                labels.push_back("capacity " + cons_detail::trim(cap) + " of place " +
                                 t.names_node[ind - 1]);
            }
            for (std::size_t k = 0; k < t.K; ++k) {
                if (t.pidx[ind - 1][k] < 0) continue;
                if (k >= sn.stations[ist - 1].classcap.size()) continue;
                const double ccap = sn.stations[ist - 1].classcap[k];
                if (!std::isfinite(ccap)) continue;
                std::vector<double> row(t.nstate, 0.0);
                row[static_cast<std::size_t>(t.pidx[ind - 1][k])] = 1.0;
                rows.push_back(row);
                bs.push_back(ccap);
                labels.push_back("class-" + std::to_string(k + 1) + " capacity " +
                                 cons_detail::trim(ccap) + " of place " + t.names_node[ind - 1]);
            }
        }
    }
    // TWO ROWS THAT SAY THE SAME THING ARE A SINGULAR NEWTON SYSTEM, not a
    // redundancy the least squares absorbs, so an exact duplicate is pruned and
    // the tighter bound survives.
    std::vector<bool> keep(rows.size(), true);
    for (std::size_t c = 0; c < rows.size(); ++c) {
        if (!keep[c]) continue;
        for (std::size_t d = c + 1; d < rows.size(); ++d) {
            if (!keep[d] || rows[c] != rows[d]) continue;
            if (bs[d] < bs[c]) {
                keep[c] = false;
                break;
            }
            keep[d] = false;
        }
    }
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < keep.size(); ++i)
        if (keep[i]) idx.push_back(i);
    PetriConstraints con;
    con.A = Matrix<double>(idx.size(), t.nstate, 0.0);
    con.b.assign(idx.size(), 0.0);
    con.cover.assign(idx.size(), std::vector<bool>(t.nstate, false));
    for (std::size_t r = 0; r < idx.size(); ++r) {
        for (std::size_t s = 0; s < t.nstate; ++s) {
            con.A(r, s) = rows[idx[r]][s];
            con.cover[r][s] = rows[idx[r]][s] > 0.0;
        }
        con.b[r] = bs[idx[r]];
        con.label.push_back(labels[idx[r]]);
    }
    return con;
}

// ============================ the immediates =============================

/** The active set of the immediate modes and the equations pinning their flows. */
struct PetriImmediate {
    enum Kind { PIN = 0, RATIO = 1, ZERO = 2 };
    struct Row {
        Kind kind = PIN;
        std::size_t a = 0, b = 0;
        double wa = 0.0, wb = 0.0;
    };
    std::size_t n = 0;
    std::vector<bool> active;
    std::vector<std::ptrdiff_t> bind;
    std::vector<std::size_t> pins;
    std::vector<Row> rows;
    bool initialized = false;
};

namespace imm_detail {

/**
 * True when an inhibitor arc of this mode has reached its threshold.
 *
 * A HARD TEST ON THE MEAN, AND A KNOWN WRONG ANSWER WHEN THE MEAN SITS ON THE
 * THRESHOLD -- a timed mode closes the same indicator smoothly, as
 * `Phi((thr-m)/sd)` in `petri_theta`; this path does not. See
 * `_kb/06-solver-catalog.md` for the measurement and for what a real fix costs.
 */
inline bool inhibited(const PetriMode& md, const std::vector<double>& x) {
    for (std::size_t b = 0; b < md.inh_slot.size(); ++b)
        if (x[md.inh_slot[b]] >= md.inh_thr[b]) return true;
    return false;
}

}  // namespace imm_detail

/**
 * An immediate transition has no rate: its fluid limit is a FLOW, an algebraic
 * unknown pinned by the constraint that its binding input place holds no mass,
 *
 *     phi_j >= 0,   x_b = 0 for the coordinate b that binds mode j
 *
 * with the GSPN conflict rule supplying the extra equation when two modes drain
 * one place: `phi_j*weight_l = phi_l*weight_j` among the enabled modes of highest
 * firing priority, and `phi = 0` below it.
 *
 * THE COUNT IS SQUARE BY CONSTRUCTION: V pins plus (F-V) ratio rows is F
 * equations for F flows.
 */
inline void petri_immediate(const PetriTerms& t, const std::vector<double>& x,
                            PetriImmediate& imm) {
    const std::size_t n = t.imm_idx.size();
    if (!imm.initialized) {
        imm.n = n;
        imm.active.assign(n, true);
        imm.bind.assign(n, -1);
        imm.pins.clear();
        // An inhibited mode never fires, so it neither carries a flow nor empties
        // a place. AN EMPTY INPUT PLACE IS NOT A REASON TO DEACTIVATE -- that is
        // the normal state of an enabled immediate mode, and the whole content of
        // its pin.
        for (std::size_t k = 0; k < n; ++k) {
            const PetriMode& md = t.modes[t.imm_idx[k]];
            if (md.arc_slot.empty())
                throw InputError("fluid_petri_immediate: immediate mode " + md.label +
                                 " has no enabling arc, so nothing bounds its firing flow and "
                                 "the net has no fluid limit. Give it an input place, or make "
                                 "it timed");
            if (imm_detail::inhibited(md, x)) imm.active[k] = false;
        }
        imm.initialized = true;
    }

    // ---- the assignment: each active mode binds the input arc it is shortest of
    for (std::size_t k = 0; k < n; ++k) {
        if (!imm.active[k]) {
            imm.bind[k] = -1;
            continue;
        }
        const PetriMode& md = t.modes[t.imm_idx[k]];
        if (imm.bind[k] >= 0 &&
            std::find(md.arc_slot.begin(), md.arc_slot.end(),
                      static_cast<std::size_t>(imm.bind[k])) != md.arc_slot.end())
            continue;  // a binding the caller set explicitly is kept
        std::size_t best = md.arc_slot[0];
        double best_lev = x[md.arc_slot[0]] / md.arc_w[0];
        for (std::size_t a = 1; a < md.arc_slot.size(); ++a) {
            const double lev = x[md.arc_slot[a]] / md.arc_w[a];
            if (lev < best_lev) {
                best_lev = lev;
                best = md.arc_slot[a];
            }
        }
        imm.bind[k] = static_cast<std::ptrdiff_t>(best);
    }

    // ---- the equations
    std::set<std::size_t> pinset;
    for (std::size_t k = 0; k < n; ++k)
        if (imm.active[k] && imm.bind[k] >= 0)
            pinset.insert(static_cast<std::size_t>(imm.bind[k]));
    imm.pins.assign(pinset.begin(), pinset.end());
    imm.rows.clear();
    for (std::size_t p : imm.pins) {
        PetriImmediate::Row row;
        row.kind = PetriImmediate::PIN;
        row.a = p;
        imm.rows.push_back(row);
        std::vector<std::size_t> grp;
        for (std::size_t k = 0; k < n; ++k)
            if (imm.active[k] && imm.bind[k] == static_cast<std::ptrdiff_t>(p)) grp.push_back(k);
        if (grp.size() <= 1) continue;
        int topprio = std::numeric_limits<int>::min();
        for (std::size_t k : grp) topprio = std::max(topprio, t.modes[t.imm_idx[k]].prio);
        std::vector<std::size_t> top, low;
        for (std::size_t k : grp)
            (t.modes[t.imm_idx[k]].prio == topprio ? top : low).push_back(k);
        const double w0 = t.modes[t.imm_idx[top[0]]].weight;
        for (std::size_t i = 1; i < top.size(); ++i) {
            PetriImmediate::Row rr;
            rr.kind = PetriImmediate::RATIO;
            rr.a = top[0];
            rr.b = top[i];
            rr.wa = w0;
            rr.wb = t.modes[t.imm_idx[top[i]]].weight;
            imm.rows.push_back(rr);
        }
        for (std::size_t k : low) {
            PetriImmediate::Row rr;
            rr.kind = PetriImmediate::ZERO;
            rr.a = k;
            imm.rows.push_back(rr);
        }
    }
    for (std::size_t k = 0; k < n; ++k)
        if (!imm.active[k]) {
            PetriImmediate::Row rr;
            rr.kind = PetriImmediate::ZERO;
            rr.a = k;
            imm.rows.push_back(rr);
        }
}

// ============================ the refusals ===============================

/** Whether the fluid Petri route can answer this model, and why not. */
struct PetriVerdict {
    bool ok = true;
    std::string reason;
};

/**
 * A QUEUEING STATION IS THE ONE STRUCTURAL EXCLUSION: a net whose tokens also
 * visit a Queue or a Delay is two formalisms at once, and LINE has no reference
 * semantics for the hand-off.
 */
template <class T>
inline PetriVerdict petri_applicable(const qn::NetworkStruct<T>& sn) {
    PetriVerdict v;
    const std::size_t I = sn.nodes.size();
    bool has_transition = false;
    for (std::size_t i = 0; i < I; ++i)
        if (sn.nodes[i].nodetype == lang::NodeType::Transition) has_transition = true;
    if (!has_transition) {
        v.ok = false;
        v.reason = "the model has no Transition node, so it is not a Petri net";
        return v;
    }
    for (std::size_t i = 0; i < I; ++i) {
        const lang::NodeType nt = sn.nodes[i].nodetype;
        if (nt != lang::NodeType::Place && nt != lang::NodeType::Transition &&
            nt != lang::NodeType::Source && nt != lang::NodeType::Sink) {
            v.ok = false;
            v.reason = "node " + sn.nodes[i].name + " is a " +
                       std::string(lang::node_type_to_text(nt)) +
                       ". The fluid Petri route solves the marking of a Petri net, and a model "
                       "that also holds queueing stations is two formalisms at once with no "
                       "reference semantics for the hand-off; use SolverCTMC, SolverJMT, "
                       "SolverSSA or SolverLDES";
            return v;
        }
    }
    // A queueing place declares a service process, which is what turns it into a
    // station with an embedded queue and a depository.
    for (std::size_t i = 0; i < I; ++i) {
        if (sn.nodes[i].nodetype != lang::NodeType::Place) continue;
        const std::size_t ist = sn.nodes[i].station;
        if (ist < 1 || ist > sn.nstations) continue;
        for (std::size_t k = 0; k < sn.nclasses; ++k) {
            const double rate = num_traits<T>::to_double(sn.rates(ist - 1, k));
            if (!std::isnan(rate) && rate > 0) {
                v.ok = false;
                v.reason = "place " + sn.nodes[i].name +
                           " is a QUEUEING place (it declares a service process), whose embedded "
                           "queue this drift does not carry; use SolverLDES";
                return v;
            }
        }
    }
    return v;
}

// ============================ the driver =================================

/** What the Petri route computes that the station table has no column for. */
struct PetriReport {
    Matrix<double> marking, marking_var;
    std::vector<std::string> mode_label;
    std::vector<double> mode_flow, immediate_flow;
    std::vector<std::string> invariant_label;
    std::vector<double> invariant_value, invariant_error;
    std::vector<std::string> capacity_label;
    std::vector<std::size_t> capacity_active;
    std::vector<double> capacity_fraction;
    std::vector<std::size_t> pinned;
    Matrix<double> Sigma;
};

/** Everything `solver_fluid_petri` returns. */
struct PetriSolution {
    Matrix<double> QN, UN, RN, TN;
    std::vector<double> t;
    std::vector<std::vector<double>> xvec_t;
    std::vector<std::vector<std::vector<double>>> QNt, UNt, TNt;
    std::vector<double> x;
    std::size_t iters = 0;
    double resnorm = std::numeric_limits<double>::infinity();
    bool converged = false;
    double runtime = 0.0;
    Matrix<double> Sigma, QVar, QStd;
    PetriReport petri;
    std::vector<std::string> warnings;
};

/** Tuning of the outer solve. */
struct PetriOptions {
    std::size_t dae_maxstate = 100;
    double tol = 1e-8;
    std::size_t newton_max = 50;
};

namespace solver_detail {

/** Everything the residual needs that does not change within one active set. */
struct Ctx {
    const PetriTerms* terms = nullptr;
    const PetriImmediate* imm = nullptr;
    std::vector<std::size_t> active;
    const PetriConstraints* con = nullptr;
    Matrix<double> C;
    std::vector<double> N;
    Matrix<double> Dp, Dn;
};

/**
 * THE CONSERVATION ROWS A BINDING CAP BREAKS ARE DROPPED: a capped place loses
 * the tokens that do not fit, so a conserved quantity supported on it is not
 * conserved while the cap binds, and keeping its row would state an equation the
 * drift contradicts -- a singular Newton system rather than an inaccuracy.
 */
inline Ctx context(const PetriTerms& terms, const PetriConservation& cons,
                   const PetriConstraints& con, const PetriImmediate& imm,
                   const std::vector<std::size_t>& active) {
    Ctx ctx;
    ctx.terms = &terms;
    ctx.imm = &imm;
    ctx.active = active;
    ctx.con = &con;
    Matrix<double> C = cons.C;
    std::vector<double> N = cons.N;
    if (!active.empty() && C.rows() > 0) {
        std::vector<bool> hit(terms.nstate, false);
        for (std::size_t c : active)
            for (std::size_t s = 0; s < terms.nstate; ++s)
                if (con.cover[c][s]) hit[s] = true;
        std::vector<std::size_t> keep;
        for (std::size_t r = 0; r < C.rows(); ++r) {
            bool drop = false;
            for (std::size_t s = 0; s < terms.nstate && !drop; ++s)
                if (hit[s] && C(r, s) != 0.0) drop = true;
            if (!drop) keep.push_back(r);
        }
        Matrix<double> Ck(keep.size(), terms.nstate, 0.0);
        std::vector<double> Nk(keep.size(), 0.0);
        for (std::size_t i = 0; i < keep.size(); ++i) {
            for (std::size_t s = 0; s < terms.nstate; ++s) Ck(i, s) = C(keep[i], s);
            Nk[i] = N[keep[i]];
        }
        C = Ck;
        N = Nk;
    }
    ctx.C = C;
    ctx.N = N;
    ctx.Dp = Matrix<double>(terms.D.rows(), terms.D.cols(), 0.0);
    ctx.Dn = Matrix<double>(terms.D.rows(), terms.D.cols(), 0.0);
    for (std::size_t i = 0; i < terms.D.rows(); ++i)
        for (std::size_t j = 0; j < terms.D.cols(); ++j) {
            ctx.Dp(i, j) = std::max(terms.D(i, j), 0.0);
            ctx.Dn(i, j) = std::min(terms.D(i, j), 0.0);
        }
    return ctx;
}

/** The unknown vector, split into its blocks. */
struct Unpacked {
    std::vector<double> x, s2, phi, mu, zeta;
};

inline Unpacked unpack(const std::vector<double>& u, const PetriTerms& t,
                       const PetriImmediate& imm, std::size_t na) {
    const std::size_t n = t.nstate, npair = t.npair, ni = imm.n, nl = t.latch_mode.size();
    Unpacked up;
    up.x.assign(u.begin(), u.begin() + n);
    up.s2.assign(u.begin() + n, u.begin() + n + npair);
    up.phi.assign(u.begin() + n + npair, u.begin() + n + npair + ni);
    for (std::size_t i = 0; i < up.phi.size(); ++i) up.phi[i] = std::max(0.0, up.phi[i]);
    up.mu.assign(u.begin() + n + npair + ni, u.begin() + n + npair + ni + nl);
    up.zeta.assign(u.begin() + n + npair + ni + nl, u.begin() + n + npair + ni + nl + na);
    for (std::size_t i = 0; i < up.zeta.size(); ++i) up.zeta[i] = std::max(0.0, up.zeta[i]);
    return up;
}

inline std::ptrdiff_t index_of(const std::vector<std::size_t>& a, std::size_t v) {
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] == v) return static_cast<std::ptrdiff_t>(i);
    return -1;
}

/**
 * The reduction of the fluctuation onto the manifold the fast and clamped
 * directions leave free.
 *
 * AN IMMEDIATE PIN REDUCES OBLIQUELY, ALONG THE FAST REACTION ITSELF, and this
 * is the one place where the orthogonal projector the queueing twin uses is
 * WRONG rather than merely different: a slow event depositing into a pinned
 * place is answered instantly by the immediate transition, so its effective jump
 * is its own plus the immediate flow it triggers -- the token is forwarded, not
 * lost. An orthogonal projection deletes the deposit and destroys mass in the
 * diffusion.
 *
 *     P = I - Cf * G * E_B,   G = G0 * (E_B Cf G0)^-1
 *
 * A CAPACITY CAP REDUCES ORTHOGONALLY -- mass that does not fit is genuinely
 * lost, so there is nothing to forward it to. A SERVER LATCH REDUCES
 * ORTHOGONALLY TOO, on the LINEARISED row `[-dtheta_j/dm, 1 over the phase
 * block]`, which is why this projector depends on the iterate.
 *
 * THE PSEUDO-INVERSE IS THE RANK-REVEALING ONE, not a regularised normal-equation
 * solve: the oblique reduction must ANNIHILATE the pinned coordinate exactly, and
 * a Tikhonov term of 1e-12 leaves it at 1e-12, after which the reduced generator
 * keeps a marginal eigenvalue and `fluid_lyapunov` refuses the fixed point as
 * non-hyperbolic -- a failure to reduce reported as a property of the model.
 */
inline bool clamp_tangent(const PetriTerms& terms, const PetriImmediate& imm,
                          const PetriConstraints& con, const std::vector<std::size_t>& active,
                          const PetriTheta* th, Matrix<double>& T) {
    const std::vector<std::size_t>& idx = terms.cov_idx;
    const std::size_t nc = idx.size();
    T = Matrix<double>(nc, nc, 0.0);
    for (std::size_t i = 0; i < nc; ++i) T(i, i) = 1.0;

    std::vector<std::size_t> actk;
    for (std::size_t k = 0; k < imm.n; ++k)
        if (imm.active[k] && imm.bind[k] >= 0) actk.push_back(k);
    const std::vector<std::size_t>& B = imm.pins;
    if (!actk.empty() && !B.empty()) {
        Matrix<double> Cf(nc, actk.size(), 0.0);
        for (std::size_t a = 0; a < actk.size(); ++a) {
            const std::vector<double>& cv = terms.modes[terms.imm_idx[actk[a]]].cvec;
            for (std::size_t i = 0; i < nc; ++i) Cf(i, a) = cv[idx[i]];
        }
        Matrix<double> G0(actk.size(), B.size(), 0.0);
        for (std::size_t jb = 0; jb < B.size(); ++jb) {
            std::vector<std::size_t> grp;
            double wsum = 0.0;
            for (std::size_t q = 0; q < actk.size(); ++q)
                if (imm.bind[actk[q]] == static_cast<std::ptrdiff_t>(B[jb])) {
                    grp.push_back(q);
                    wsum += terms.modes[terms.imm_idx[actk[q]]].weight;
                }
            if (grp.empty()) continue;
            const bool uniform = !(wsum > 0);
            for (std::size_t q : grp) {
                const double w =
                    uniform ? 1.0 : terms.modes[terms.imm_idx[actk[q]]].weight;
                G0(q, jb) = w / (uniform ? static_cast<double>(grp.size()) : wsum);
            }
        }
        Matrix<double> EB(B.size(), nc, 0.0);
        for (std::size_t jb = 0; jb < B.size(); ++jb) {
            const std::ptrdiff_t at = index_of(idx, B[jb]);
            if (at >= 0) EB(jb, static_cast<std::size_t>(at)) = 1.0;
        }
        const Matrix<double> Mb = matmul(matmul(EB, Cf), G0);
        const Matrix<double> corr = matmul(matmul(Cf, matmul(G0, pinv(Mb))), EB);
        for (std::size_t i = 0; i < nc; ++i)
            for (std::size_t j = 0; j < nc; ++j) T(i, j) -= corr(i, j);
    }

    std::vector<std::vector<double>> R;
    for (std::size_t c : active) {
        std::vector<double> row(nc, 0.0);
        for (std::size_t i = 0; i < nc; ++i) row[i] = con.A(c, idx[i]);
        R.push_back(row);
    }
    if (th != nullptr) {
        for (std::size_t j : terms.latch_mode) {
            std::vector<double> row(nc, 0.0);
            for (std::size_t z : terms.modes[j].zblk) {
                const std::ptrdiff_t at = index_of(idx, z);
                if (at >= 0) row[static_cast<std::size_t>(at)] = 1.0;
            }
            for (std::size_t q = 0; q < th->dslot[j].size(); ++q) {
                const std::ptrdiff_t at = index_of(idx, th->dslot[j][q]);
                if (at >= 0) row[static_cast<std::size_t>(at)] -= th->dval[j][q];
            }
            R.push_back(row);
        }
    }
    if (!R.empty()) {
        bool any = false;
        for (const std::vector<double>& row : R)
            for (double v : row)
                if (std::fabs(v) > 1e-14) any = true;
        if (any) {
            Matrix<double> Rm(R.size(), nc, 0.0);
            for (std::size_t i = 0; i < R.size(); ++i)
                for (std::size_t j = 0; j < nc; ++j) Rm(i, j) = R[i][j];
            Matrix<double> Rt(nc, R.size(), 0.0);
            for (std::size_t i = 0; i < R.size(); ++i)
                for (std::size_t j = 0; j < nc; ++j) Rt(j, i) = Rm(i, j);
            const Matrix<double> RRt = matmul(Rm, Rt);
            const Matrix<double> corr = matmul(matmul(Rt, pinv(RRt)), Rm);
            Matrix<double> P(nc, nc, 0.0);
            for (std::size_t i = 0; i < nc; ++i)
                for (std::size_t j = 0; j < nc; ++j) P(i, j) = (i == j ? 1.0 : 0.0) - corr(i, j);
            T = matmul(P, T);
        }
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < nc; ++j)
            worst = std::max(worst, std::fabs(T(i, j) - (i == j ? 1.0 : 0.0)));
    return worst > 1e-14;
}

/**
 * One Lyapunov solve, over the marking coordinates.
 *
 * THE DIFFUSION COUNTS THE STOCHASTIC EVENTS ONLY: an immediate flow is not a
 * Poisson stream with an intensity but the limit of an infinitely fast one whose
 * fluctuation is slaved, and its pinned coordinate is projected out.
 */
inline Matrix<double> sigma_of(const PetriTerms& terms, const Matrix<double>& A,
                               const std::vector<double>& r, const Matrix<double>* clampT) {
    const std::vector<std::size_t>& idx = terms.cov_idx;
    const std::size_t nc = idx.size(), ns = terms.stoch_col.size();
    Matrix<double> Dc(nc, std::max<std::size_t>(ns, 1), 0.0);
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < ns; ++j) Dc(i, j) = terms.D(idx[i], terms.stoch_col[j]);
    Matrix<double> Am(nc, nc, 0.0);
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < nc; ++j) Am(i, j) = A(idx[i], idx[j]);
    if (clampT != nullptr) {
        // BOTH the jump directions and the generator are reduced: reducing Dc
        // alone would fix the subspace but leave the generator's orthogonal
        // component on it, which is not the reduced dynamics when the reduction
        // is oblique.
        Dc = matmul(*clampT, Dc);
        Am = matmul(*clampT, Am);
    }
    Matrix<double> Q(nc, nc, 0.0);
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < nc; ++j) {
            double v = 0.0;
            for (std::size_t e = 0; e < ns; ++e)
                v += Dc(i, e) * r[terms.stoch_col[e]] * Dc(j, e);
            Q(i, j) = v;
        }
    FluidLyapunovInfo info;
    const Matrix<double> Sc = fluid_lyapunov(Am, Q, Dc, info);
    Matrix<double> Sigma(terms.nstate, terms.nstate, 0.0);
    for (std::size_t i = 0; i < nc; ++i)
        for (std::size_t j = 0; j < nc; ++j) Sigma(idx[i], idx[j]) = Sc(i, j);
    return Sigma;
}

/** The stacked residual, the rates it was evaluated at, and the covariance. */
struct Residual {
    bool ok = false;
    std::vector<double> G, r;
    PetriTheta th;
    Matrix<double> Sigma;
};

/**
 * The coupled algebraic system, stacked.
 *
 * `ok` is false when the closure cannot be evaluated at this iterate, so the line
 * search can back off; the first evaluation of a pass runs with `quiet=false`,
 * where a genuine failure surfaces.
 */
inline Residual residual(const std::vector<double>& u, const Ctx& ctx, bool quiet) {
    const PetriTerms& terms = *ctx.terms;
    const PetriImmediate& imm = *ctx.imm;
    const std::size_t n = terms.nstate;
    const Unpacked up = unpack(u, terms, imm, ctx.active.size());

    Residual res;
    res.th = petri_theta(terms, up.x, up.s2);
    res.r = petri_rates(terms, up.x, up.phi, up.mu, res.th);

    // the deposit gate of every binding capacity, as a product of fractions
    std::vector<double> gain(n, 1.0);
    for (std::size_t k = 0; k < ctx.active.size(); ++k)
        for (std::size_t s = 0; s < n; ++s)
            if (ctx.con->cover[ctx.active[k]][s]) gain[s] *= up.zeta[k];
    std::vector<double> drift(n, 0.0);
    for (std::size_t s = 0; s < n; ++s) {
        double neg = 0.0, pos = 0.0;
        for (std::size_t e = 0; e < terms.nev; ++e) {
            neg += ctx.Dn(s, e) * res.r[e];
            pos += ctx.Dp(s, e) * res.r[e];
        }
        drift[s] = neg + gain[s] * pos;
    }

    try {
        const Matrix<double> A = petri_jacobian(terms, res.th);
        Matrix<double> T;
        const bool reduced =
            clamp_tangent(terms, imm, *ctx.con, ctx.active, &res.th, T);
        res.Sigma = sigma_of(terms, A, res.r, reduced ? &T : nullptr);
    } catch (const Error&) {
        if (!quiet) throw;
        res.ok = false;
        return res;
    }

    std::vector<double> G;
    G.reserve(n + ctx.C.rows() + terms.npair + imm.rows.size() + terms.latch_mode.size() +
              ctx.active.size());
    for (double d : drift) G.push_back(d);
    for (std::size_t r = 0; r < ctx.C.rows(); ++r) {
        double v = 0.0;
        for (std::size_t s = 0; s < n; ++s) v += ctx.C(r, s) * up.x[s];
        G.push_back(v - ctx.N[r]);
    }
    for (std::size_t p = 0; p < terms.npair; ++p)
        G.push_back(up.s2[p] - res.Sigma(terms.cov_pairs[p].first, terms.cov_pairs[p].second));
    for (const PetriImmediate::Row& row : imm.rows) {
        if (row.kind == PetriImmediate::PIN) G.push_back(up.x[row.a]);
        else if (row.kind == PetriImmediate::RATIO)
            G.push_back(up.phi[row.a] * row.wb - up.phi[row.b] * row.wa);
        else G.push_back(up.phi[row.a]);
    }
    for (std::size_t j : terms.latch_mode) {
        double s = 0.0;
        for (std::size_t z : terms.modes[j].zblk) s += up.x[z];
        G.push_back(s - res.th.theta[j]);
    }
    for (std::size_t k = 0; k < ctx.active.size(); ++k) {
        double v = 0.0;
        for (std::size_t s = 0; s < n; ++s) v += ctx.con->A(ctx.active[k], s) * up.x[s];
        G.push_back(v - ctx.con->b[ctx.active[k]]);
    }
    res.G = G;
    res.ok = true;
    return res;
}

inline double inf_norm(const std::vector<double>& v) {
    double m = 0.0;
    for (double d : v) m = std::max(m, std::fabs(d));
    return m;
}

/**
 * An iterate projected onto its feasible box, the lower bound only.
 *
 * The bound is a VECTOR, one entry per unknown, never a count: a scalar "nfree"
 * is indistinguishable from a one-unknown bound vector, which is how the MATLAB
 * twin crashed on the simplest net in the tree.
 */
inline std::vector<double> project(const std::vector<double>& u, const std::vector<double>& lb) {
    if (lb.empty()) return u;
    if (lb.size() != u.size())
        throw InputError("fluid_petri: the bound vector has " + std::to_string(lb.size()) +
                         " entries for " + std::to_string(u.size()) + " unknowns");
    std::vector<double> w = u;
    for (std::size_t i = 0; i < w.size(); ++i)
        if (std::isfinite(lb[i])) w[i] = std::max(lb[i], w[i]);
    return w;
}

/** Forward-difference Jacobian of the residual. */
inline Matrix<double> fdjac(const Ctx& ctx, const std::vector<double>& u,
                            const std::vector<double>& G) {
    const std::size_t n = u.size(), m = G.size();
    Matrix<double> J(std::max<std::size_t>(m, 1), std::max<std::size_t>(n, 1), 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const double h = 1e-7 * std::max(1.0, std::fabs(u[k]));
        std::vector<double> up = u;
        up[k] += h;
        Residual rp = residual(up, ctx, true);
        if (rp.ok) {
            for (std::size_t i = 0; i < m; ++i) J(i, k) = (rp.G[i] - G[i]) / h;
            continue;
        }
        up[k] = u[k] - h;
        Residual rm = residual(up, ctx, true);
        if (!rm.ok) continue;
        for (std::size_t i = 0; i < m; ++i) J(i, k) = (G[i] - rm.G[i]) / h;
    }
    return J;
}

struct NewtonResult {
    std::vector<double> u;
    std::size_t iterations = 0;
    bool converged = false;
    double resnorm = std::numeric_limits<double>::infinity();
};

/** Damped projected Newton with a finite-difference Jacobian. */
inline NewtonResult newton(const Ctx& ctx, const std::vector<double>& u0, double tol,
                           std::size_t maxit, const std::vector<double>& lb) {
    NewtonResult out;
    std::vector<double> u = project(u0, lb);
    Residual res = residual(u, ctx, false);
    if (!res.ok) {
        out.u = u;
        return out;
    }
    std::vector<double> G = res.G;
    double resnorm = inf_norm(G);
    std::size_t it = 0;
    for (it = 1; it <= maxit; ++it) {
        if (resnorm <= tol) {
            out.u = u;
            out.iterations = it - 1;
            out.converged = true;
            out.resnorm = resnorm;
            return out;
        }
        const Matrix<double> J = fdjac(ctx, u, G);
        Matrix<double> rhs(G.size(), 1, 0.0);
        for (std::size_t i = 0; i < G.size(); ++i) rhs(i, 0) = -G[i];
        const Matrix<double> du = matmul(pinv(J), rhs);
        double lam = 1.0;
        bool improved = false;
        for (int b = 0; b < 30; ++b) {
            std::vector<double> un(u.size(), 0.0);
            for (std::size_t i = 0; i < u.size(); ++i) un[i] = u[i] + lam * du(i, 0);
            un = project(un, lb);
            Residual rn = residual(un, ctx, true);
            if (rn.ok) {
                const double v = inf_norm(rn.G);
                if (v < resnorm) {
                    u = un;
                    G = rn.G;
                    resnorm = v;
                    improved = true;
                    break;
                }
            }
            lam *= 0.5;
        }
        if (!improved) break;
    }
    out.u = u;
    out.iterations = it;
    out.converged = resnorm <= tol;
    out.resnorm = resnorm;
    return out;
}

/**
 * The first-order drift: the same rates at zero variance, with an immediate mode
 * firing at LAM times its enabling degree and its firing weight, and the server
 * latch relaxed at LAM towards the enabling degree instead of solved. Both are
 * approximations of an algebraic constraint by a fast reaction, and both are
 * confined to the seed.
 */
inline std::vector<double> seed_drift(const PetriTerms& terms, const std::vector<double>& xin,
                                      double lam) {
    std::vector<double> x = xin;
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::max(0.0, x[i]);
    const PetriTheta th =
        petri_theta(terms, x, std::vector<double>(std::max<std::size_t>(terms.npair, 1), 0.0));
    std::vector<double> phi(terms.imm_idx.size(), 0.0);
    for (std::size_t k = 0; k < phi.size(); ++k) {
        const std::size_t j = terms.imm_idx[k];
        phi[k] = lam * terms.modes[j].weight * th.theta[j];
    }
    std::vector<double> mu(terms.latch_mode.size(), 0.0);
    for (std::size_t q = 0; q < mu.size(); ++q) {
        const std::size_t j = terms.latch_mode[q];
        double s = 0.0;
        for (std::size_t z : terms.modes[j].zblk) s += x[z];
        mu[q] = lam * (th.theta[j] - s);
    }
    const std::vector<double> r = petri_rates(terms, x, phi, mu, th);
    std::vector<double> d(terms.nstate, 0.0);
    for (std::size_t s = 0; s < terms.nstate; ++s) {
        double v = 0.0;
        for (std::size_t e = 0; e < terms.nev; ++e) v += terms.D(s, e) * r[e];
        d[s] = v;
    }
    return d;
}

}  // namespace solver_detail

/**
 * Fluid analysis of a stochastic Petri net.
 *
 * @param sn  a model whose nodes are Places, Transitions, Sources and Sinks only
 * @param opt the state-size cap, the Newton tolerance and its iteration cap
 */
template <class T>
inline PetriSolution solver_fluid_petri(const qn::NetworkStruct<T>& sn,
                                        const PetriOptions& opt = PetriOptions()) {
    const std::size_t M = sn.nstations, K = sn.nclasses;

    const PetriVerdict v = petri_applicable(sn);
    if (!v.ok)
        throw UnsupportedError("solver_fluid_petri: the fluid Petri route cannot solve this "
                               "model: " + v.reason + ".");

    const PetriTerms terms = petri_build_terms(sn);
    const std::size_t n = terms.nstate, npair = terms.npair;
    if (n > opt.dae_maxstate)
        throw UnsupportedError(
            "solver_fluid_petri: the fluid Petri route solves a " + std::to_string(n) +
            "-unknown algebraic system with a finite-difference Jacobian, above the limit of " +
            std::to_string(opt.dae_maxstate) +
            " set by options.config.dae_maxstate. Raise that limit, or use SolverSSA for a net "
            "of this size");

    const PetriConservation cons = petri_conservation(terms);
    if (cons.leak > 1e-7)
        throw NumericError("solver_fluid_petri: the conserved directions and the jump matrix "
                           "disagree: the largest leak per unit rate is " +
                           std::to_string(cons.leak) + ", where it must be zero");
    const PetriConstraints con = petri_constraints(sn, terms);
    const std::size_t ncon = con.b.size();

    // ---- seed ---------------------------------------------------------------
    // Newton needs a point in the basin, not an answer. The immediate modes get a
    // large FINITE rate here and only here, scaled to the model's own timescale
    // rather than taken from GlobalConstants.Immediate: 1e8 against a rate of
    // order one is a stiffness the seed does not need, and the answer does not
    // depend on the seed's accuracy.
    double rmax = 0.0;
    for (std::size_t j : terms.timed_idx) {
        double s = 0.0;
        for (double d : terms.modes[j].d1) s += d;
        rmax = std::max(rmax, s);
    }
    for (std::size_t e = 0; e < terms.nev; ++e)
        if (terms.ev_kind[e] == 3) rmax = std::max(rmax, terms.rate_base[e]);
    if (!(rmax > 0)) rmax = 1.0;
    const double lam = std::min(1e8, 1e4 * rmax);
    double mass = 0.0;
    for (std::size_t s = 0; s < terms.nm; ++s) mass += terms.x0[s];
    double Thor = 50.0 * (mass + 1.0) / rmax;

    const auto drift_fn = [&terms, lam](const double&, const std::vector<double>& y) {
        return solver_detail::seed_drift(terms, y, lam);
    };
    std::vector<double> tseed;
    std::vector<std::vector<double>> xseed;
    for (int attempt = 0; attempt < 6; ++attempt) {
        OdeOptions<double> oo;
        oo.rtol = 1e-6;
        oo.atol = 1e-10;
        oo.store_trajectory = true;
        oo.max_steps = 200000;
        const OdeSolution<double> sol =
            ode_rosenbrock4(drift_fn, 0.0, Thor, terms.x0, oo);
        tseed = sol.t;
        xseed = sol.y;
        const std::vector<double> d = solver_detail::seed_drift(terms, xseed.back(), lam);
        if (solver_detail::inf_norm(d) <= 1e-6 * std::max(1.0, rmax * (mass + 1.0))) break;
        Thor *= 4.0;
    }
    std::vector<double> x = xseed.back();

    PetriImmediate imm;
    petri_immediate(terms, x, imm);
    const std::size_t nimm = imm.n;

    // THE VARIANCE IS SEEDED POSITIVE: sigma2 = 0 is where min() has no
    // derivative, and a saturated net's first-order fixed point sits there.
    std::vector<double> s2(npair, 0.0);
    std::vector<bool> ondiag(npair, false);
    for (std::size_t p = 0; p < npair; ++p) {
        ondiag[p] = (terms.cov_pairs[p].first == terms.cov_pairs[p].second);
        if (ondiag[p]) s2[p] = std::max(petri_fine_tol(), x[terms.cov_pairs[p].first]);
    }

    // ---- steady state -------------------------------------------------------
    std::vector<std::size_t> active;
    std::vector<double> phi(nimm, 0.0);
    const std::size_t nlatch = terms.latch_mode.size();
    std::vector<double> mu(nlatch, 0.0), zeta;
    std::size_t iters = 0;
    const std::size_t aset_max = std::max<std::size_t>(6, 2 * (ncon + nimm) + 2);
    bool converged = false;
    double resnorm = std::numeric_limits<double>::infinity();
    bool have_best = false;
    std::vector<double> bx, bs2, bphi, bmu, bzeta;
    std::vector<std::size_t> bactive;
    PetriImmediate bimm;
    double bres = std::numeric_limits<double>::infinity();
    bool bconv = false;

    for (std::size_t sweep = 0; sweep < aset_max; ++sweep) {
        const solver_detail::Ctx ctx =
            solver_detail::context(terms, cons, con, imm, active);
        std::vector<double> u0;
        u0.insert(u0.end(), x.begin(), x.end());
        u0.insert(u0.end(), s2.begin(), s2.end());
        u0.insert(u0.end(), phi.begin(), phi.end());
        u0.insert(u0.end(), mu.begin(), mu.end());
        u0.insert(u0.end(), active.size(), 1.0);
        std::vector<double> lb(u0.size(), -std::numeric_limits<double>::infinity());
        for (std::size_t k = 0; k < nimm; ++k) lb[n + npair + k] = 0.0;
        for (std::size_t k = 0; k < active.size(); ++k)
            lb[n + npair + nimm + nlatch + k] = 0.0;
        for (std::size_t p = 0; p < npair; ++p)
            if (ondiag[p]) lb[n + p] = 0.0;

        const solver_detail::NewtonResult nr =
            solver_detail::newton(ctx, u0, opt.tol, opt.newton_max, lb);
        iters += nr.iterations;
        resnorm = nr.resnorm;
        converged = nr.converged;
        const solver_detail::Unpacked up =
            solver_detail::unpack(nr.u, terms, imm, active.size());
        x = up.x;
        s2 = up.s2;
        phi = up.phi;
        mu = up.mu;
        zeta = up.zeta;

        if (!have_best || resnorm < bres) {
            bx = x; bs2 = s2; bphi = phi; bmu = mu; bzeta = zeta;
            bactive = active; bimm = imm; bres = resnorm; bconv = converged;
            have_best = true;
        }

        // The active-set moves, in the order that a failure of each invalidates
        // the next: a negative flow means the mode does not fire at all, a
        // negative marking means the wrong coordinate was pinned, and only then
        // is it worth asking which capacity rows bind.
        bool moved = false;
        if (nimm > 0) {
            double scale = 1.0;
            for (double p : phi) scale = std::max(scale, std::fabs(p));
            const double thr = -std::max(opt.tol, 1e-10) * scale;
            for (std::size_t k = 0; k < nimm && !moved; ++k)
                if (phi[k] < thr) {
                    imm.active[k] = false;
                    petri_immediate(terms, x, imm);
                    moved = true;
                }
        }
        if (!moved && nimm > 0) {
            const double thr = -std::max(opt.tol, 1e-10);
            for (std::size_t s = 0; s < terms.nm && !moved; ++s) {
                if (x[s] >= thr) continue;
                for (std::size_t k = 0; k < nimm; ++k) {
                    const PetriMode& md = terms.modes[terms.imm_idx[k]];
                    if (imm.active[k] &&
                        std::find(md.arc_slot.begin(), md.arc_slot.end(), s) !=
                            md.arc_slot.end() &&
                        imm.bind[k] != static_cast<std::ptrdiff_t>(s)) {
                        imm.bind[k] = static_cast<std::ptrdiff_t>(s);
                        petri_immediate(terms, x, imm);
                        moved = true;
                        break;
                    }
                }
            }
        }
        if (!moved && ncon > 0) {
            std::vector<std::size_t> over;
            for (std::size_t c = 0; c < ncon; ++c) {
                double val = 0.0;
                for (std::size_t s = 0; s < n; ++s) val += con.A(c, s) * x[s];
                if (val > con.b[c] + std::max(1e-9, opt.tol) &&
                    std::find(active.begin(), active.end(), c) == active.end())
                    over.push_back(c);
            }
            if (!over.empty()) {
                active.insert(active.end(), over.begin(), over.end());
                moved = true;
            } else {
                std::vector<std::size_t> keep;
                for (std::size_t i = 0; i < active.size(); ++i)
                    if (!(i < zeta.size() && zeta[i] > 1 + std::max(1e-9, opt.tol)))
                        keep.push_back(active[i]);
                if (keep.size() != active.size()) {
                    active = keep;
                    moved = true;
                }
            }
        }
        if (!moved) break;
    }

    PetriSolution out;
    if (have_best && !converged && bconv) {
        x = bx; s2 = bs2; phi = bphi; mu = bmu; zeta = bzeta;
        active = bactive; imm = bimm; converged = true; resnorm = bres;
    }
    if (!converged)
        out.warnings.push_back(
            "The simultaneous closure solve stopped at residual " + std::to_string(resnorm) +
            " after " + std::to_string(iters) + " Newton steps without reaching " +
            std::to_string(opt.tol) + ". The reported point is the last iterate.");
    // A NEGATIVE MARKING IS NOT ROUNDING: the fixed point wanted mass a place
    // cannot supply and no immediate mode could be rebound to pin it, so the
    // answer is outside the model's own state space and is reported as such.
    for (std::size_t s = 0; s < terms.nm; ++s)
        if (x[s] < -std::max(opt.tol, 1e-10)) {
            out.warnings.push_back(
                "The fixed point holds " + std::to_string(x[s]) + " tokens at " +
                terms.names_node[terms.coord_node[s]] +
                ", which is negative: no immediate transition could be rebound to pin that "
                "place at zero. Use SolverCTMC, SolverSSA or SolverLDES for this net.");
            break;
        }

    const solver_detail::Ctx ctx = solver_detail::context(terms, cons, con, imm, active);
    std::vector<double> ufin;
    ufin.insert(ufin.end(), x.begin(), x.end());
    ufin.insert(ufin.end(), s2.begin(), s2.end());
    ufin.insert(ufin.end(), phi.begin(), phi.end());
    ufin.insert(ufin.end(), mu.begin(), mu.end());
    ufin.insert(ufin.end(), zeta.begin(), zeta.end());
    const solver_detail::Residual fin = solver_detail::residual(ufin, ctx, false);

    // ---- the metrics --------------------------------------------------------
    // A PLACE IS AN INF STATION: its queue length is its mean token count and its
    // utilization is the same number. Its throughput is the rate at which TOKENS
    // leave it, so a consuming mode contributes its firing rate times the arc
    // multiplicity -- SolverCTMC's convention, and the one Little's law needs.
    out.QN = Matrix<double>(M, K, 0.0);
    out.UN = Matrix<double>(M, K, 0.0);
    out.RN = Matrix<double>(M, K, 0.0);
    out.TN = Matrix<double>(M, K, 0.0);
    for (std::size_t s = 0; s < terms.nm; ++s) {
        if (terms.coord_station[s] < 0) continue;
        const std::size_t i = static_cast<std::size_t>(terms.coord_station[s]);
        const std::size_t k = terms.coord_class[s];
        out.QN(i, k) += x[s];
        out.UN(i, k) = out.QN(i, k);
    }
    for (std::map<std::size_t, std::vector<std::size_t>>::const_iterator it =
             terms.consumers.begin();
         it != terms.consumers.end(); ++it) {
        const std::size_t i = it->first / K, k = it->first % K;
        const std::vector<double>& w = terms.consumer_w.at(it->first);
        double val = 0.0;
        for (std::size_t q = 0; q < it->second.size(); ++q) val += w[q] * fin.r[it->second[q]];
        out.TN(i, k) += val;
    }
    for (std::map<std::size_t, std::vector<std::size_t>>::const_iterator it =
             terms.producers.begin();
         it != terms.producers.end(); ++it) {
        const std::size_t i = it->first / K, k = it->first % K;
        double val = 0.0;
        for (std::size_t e : it->second) val += fin.r[e];
        out.TN(i, k) += val;
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (out.TN(i, k) > 1e-14) out.RN(i, k) = out.QN(i, k) / out.TN(i, k);

    out.t = tseed;
    out.xvec_t = xseed;
    out.x = x;
    out.iters = iters;
    out.resnorm = resnorm;
    out.converged = converged;
    out.Sigma = fin.Sigma;
    out.QVar = Matrix<double>(M, K, 0.0);
    out.QStd = Matrix<double>(M, K, 0.0);
    for (std::size_t s = 0; s < terms.nm; ++s) {
        if (terms.coord_station[s] < 0) continue;
        const std::size_t i = static_cast<std::size_t>(terms.coord_station[s]);
        out.QVar(i, terms.coord_class[s]) = std::max(0.0, fin.Sigma(s, s));
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) out.QStd(i, k) = std::sqrt(out.QVar(i, k));

    // ---- the report ---------------------------------------------------------
    out.petri.marking = Matrix<double>(terms.I, terms.K, 0.0);
    out.petri.marking_var = Matrix<double>(terms.I, terms.K, 0.0);
    for (std::size_t s = 0; s < terms.nm; ++s) {
        out.petri.marking(terms.coord_node[s], terms.coord_class[s]) = x[s];
        out.petri.marking_var(terms.coord_node[s], terms.coord_class[s]) =
            std::max(0.0, fin.Sigma(s, s));
    }
    for (const PetriMode& md : terms.modes) out.petri.mode_label.push_back(md.label);
    out.petri.mode_flow.assign(terms.modes.size(), 0.0);
    for (std::size_t j = 0; j < terms.modes.size(); ++j) {
        double val = 0.0;
        for (std::size_t e = 0; e < terms.nev; ++e)
            if (terms.ev_mode[e] == static_cast<int>(j) &&
                (terms.ev_kind[e] == 1 || terms.ev_kind[e] == 4))
                val += fin.r[e];
        out.petri.mode_flow[j] = val;
    }
    out.petri.immediate_flow = phi;
    out.petri.invariant_label = cons.label;
    out.petri.invariant_value = cons.N;
    out.petri.invariant_error.assign(cons.C.rows(), 0.0);
    for (std::size_t c = 0; c < cons.C.rows(); ++c) {
        double val = 0.0;
        for (std::size_t s = 0; s < terms.nstate; ++s) val += cons.C(c, s) * x[s];
        out.petri.invariant_error[c] = val - cons.N[c];
    }
    out.petri.capacity_label = con.label;
    out.petri.capacity_active = active;
    out.petri.capacity_fraction = zeta;
    out.petri.pinned = imm.pins;
    out.petri.Sigma = fin.Sigma;
    return out;
}

}  // namespace petri
}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_PETRI_H
