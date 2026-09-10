/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_PN_FIRING_RATES_H
#define LINE_API_SN_SN_PN_FIRING_RATES_H

/**
 * Ports of matlab/src/api/sn/sn_pn_firing_rates.m and sn_pn_avg_rates.m.
 *
 * A solver reports a PLACE's throughput as the rate at which its tokens are
 * consumed; what a stochastic Petri net is actually driven by is the firing
 * rate of each transition MODE. These two functions invert that: they recover
 * the mode firing rates x from the reported place throughputs, and then rewrite
 * the place-level throughput, arrival rate and response time so that the
 * reported tables are consistent with the token balance.
 *
 * THE SYSTEM. Two families of equations, both linear in x.
 *  (a) measurement: for each (place, class) whose consumption is non-empty,
 *      the consumed rate equals the throughput the solver reported. With
 *      tputIsTokens the consumption is weighted by the arc multiplicity, and
 *      without it only by whether the mode touches the place at all -- the
 *      difference between "tokens per unit time" and "firings per unit time",
 *      which is what the two callers respectively hold.
 *  (b) balance: production equals consumption at every place, so the marking
 *      is stationary.
 * Immediate modes are excluded from the measurement rows (they fire in zero
 * time, so no measured rate is theirs) but not from the balance rows.
 * The least-squares solution is the pseudo-inverse, as the reference has it,
 * and a solution with a materially negative rate is REJECTED rather than
 * clamped: a negative firing rate means the reported throughputs are not
 * consistent with any marking-stationary firing vector, and the caller must
 * keep its own numbers rather than be handed a repaired impossibility.
 *
 * GUARDS. The reference refuses an open net (a Source or a Sink present) and a
 * net in which a place is routed to anything other than a transition, because
 * neither shape satisfies the balance it is about to impose.
 *
 * PER-CLASS ARCS. The reference's enabling and firing tables are (nnodes x
 * nclasses); this port's `TransitionParam` sums the class dimension away when
 * the model is read, so a multi-class net cannot be answered here and is
 * refused by name rather than answered with a class-blind arc.
 *
 * ARITHMETIC: double. The pseudo-inverse is an SVD.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/svd.h"

namespace line {
namespace api {

/** What sn_pn_firing_rates returns; `x` empty means "no answer", as in the reference. */
struct SnPnFiringRates {
    std::vector<double> x;                      ///< firing rate per mode
    /** consumed[mode][place * nclasses + class], the reference's (mm, pp, k). */
    std::vector<std::vector<double>> consumed;
    /** produced[mode][place * nclasses + class], the reference's (mm, pp, k). */
    std::vector<std::vector<double>> produced;
    std::vector<std::size_t> place_nodes;       ///< 1-based Place node indices
    std::size_t nclasses = 0;                   ///< the stride of the two tables above
};

/**
 * @param TN            station throughput table, (nstations x nclasses)
 * @param tput_is_tokens whether TN counts tokens (true) or firings (false)
 */
template <class T>
SnPnFiringRates sn_pn_firing_rates(const qn::NetworkStruct<T>& sn, const Matrix<T>& TN,
                                   bool tput_is_tokens) {
    SnPnFiringRates out;
    const std::size_t R = sn.nclasses, I = sn.nodes.size();
    std::vector<std::size_t> places, trans;
    for (std::size_t a = 1; a <= I; ++a) {
        if (sn.nodes[a - 1].nodetype == qn::NodeType::Place) places.push_back(a);
        if (sn.nodes[a - 1].nodetype == qn::NodeType::Transition) trans.push_back(a);
    }
    if (places.empty() || trans.empty() || TN.rows() == 0) return out;
    for (std::size_t a = 0; a < I; ++a)
        if (sn.nodes[a].nodetype == qn::NodeType::Source ||
            sn.nodes[a].nodetype == qn::NodeType::Sink)
            return out;
    // every place must route only to transitions, or the balance does not hold
    const std::size_t S = sn.nof_stateful();
    if (sn.rt.rows() != S * R) return out;
    for (std::size_t pp = 0; pp < places.size(); ++pp) {
        const std::size_t sfp = sn.stateful_index(places[pp]);
        if (sfp == 0) return out;
        for (std::size_t sfj = 1; sfj <= S; ++sfj) {
            if (sfj == sfp) continue;
            bool touches = false;
            for (std::size_t r = 0; r < R && !touches; ++r)
                for (std::size_t s = 0; s < R; ++s)
                    if (sn.rt((sfp - 1) * R + r, (sfj - 1) * R + s) >
                            num_traits<T>::from_int(0) ||
                        sn.rt((sfj - 1) * R + s, (sfp - 1) * R + r) >
                            num_traits<T>::from_int(0)) {
                        touches = true;
                        break;
                    }
            if (touches && sn.nodes[sn.stateful_nodes[sfj - 1] - 1].nodetype !=
                               qn::NodeType::Transition)
                return out;
        }
    }

    std::vector<std::size_t> mode_trans, mode_idx;
    std::vector<bool> mode_timed;
    for (std::size_t tt = 0; tt < trans.size(); ++tt) {
        const std::size_t ind = trans[tt];
        typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
            sn.transparam.find(ind);
        if (it == sn.transparam.end() || it->second.nmodes == 0) return out;
        for (std::size_t m = 0; m < it->second.nmodes; ++m) {
            mode_trans.push_back(ind);
            mode_idx.push_back(m);
            bool timed = true;
            if (m < it->second.timing.size())
                timed = it->second.timing[m] != lang::TimingStrategy::IMMEDIATE;
            mode_timed.push_back(timed);
        }
    }
    const std::size_t nModes = mode_trans.size();
    if (nModes == 0) return out;

    out.place_nodes = places;
    out.nclasses = R;
    out.consumed.assign(nModes, std::vector<double>(places.size() * R, 0.0));
    out.produced.assign(nModes, std::vector<double>(places.size() * R, 0.0));
    for (std::size_t mm = 0; mm < nModes; ++mm) {
        const qn::TransitionParam<T>& tp = sn.transparam.at(mode_trans[mm]);
        const Matrix<T>& enab = tp.enabling[mode_idx[mm]];
        const Matrix<T>& fire = tp.firing[mode_idx[mm]];
        for (std::size_t pp = 0; pp < places.size(); ++pp) {
            const std::size_t p = places[pp] - 1;
            for (std::size_t k = 0; k < R; ++k) {
                if (p < enab.rows() && k < enab.cols())
                    out.consumed[mm][pp * R + k] =
                        std::max(0.0, num_traits<T>::to_double(enab(p, k)));
                if (p < fire.rows() && k < fire.cols())
                    out.produced[mm][pp * R + k] =
                        std::max(0.0, num_traits<T>::to_double(fire(p, k)));
            }
        }
    }

    const std::size_t nEq = 2 * places.size() * R;
    Matrix<double> A(nEq, nModes, 0.0);
    std::vector<double> b(nEq, 0.0);
    std::size_t row = 0, nMeasured = 0;
    for (std::size_t pp = 0; pp < places.size(); ++pp) {
        const std::size_t ist = sn.nodes[places[pp] - 1].station;
        for (std::size_t k = 0; k < R; ++k) {
            std::vector<double> arow(nModes, 0.0);
            bool any = false;
            for (std::size_t mm = 0; mm < nModes; ++mm) {
                const double c = out.consumed[mm][pp * R + k];
                double v = tput_is_tokens ? c : (c > 0.0 ? 1.0 : 0.0);
                if (!mode_timed[mm]) v = 0.0;
                arow[mm] = v;
                if (v != 0.0) any = true;
            }
            if (any) {
                for (std::size_t mm = 0; mm < nModes; ++mm) A(row, mm) = arow[mm];
                b[row] = ist == 0 ? 0.0 : num_traits<T>::to_double(TN(ist - 1, k));
                ++row;
                ++nMeasured;
            }
            for (std::size_t mm = 0; mm < nModes; ++mm)
                A(row, mm) = out.produced[mm][pp * R + k] - out.consumed[mm][pp * R + k];
            b[row] = 0.0;
            ++row;
        }
    }
    if (nMeasured == 0) return out;

    Matrix<double> At(row, nModes, 0.0);
    for (std::size_t a = 0; a < row; ++a)
        for (std::size_t m = 0; m < nModes; ++m) At(a, m) = A(a, m);
    const Matrix<double> Ap = pinv(At);
    std::vector<double> xfit(nModes, 0.0);
    for (std::size_t m = 0; m < nModes; ++m) {
        double acc = 0.0;
        for (std::size_t a = 0; a < row; ++a) acc += Ap(m, a) * b[a];
        xfit[m] = acc;
    }
    double mx = 1.0;
    for (std::size_t m = 0; m < nModes; ++m) mx = std::max(mx, std::fabs(xfit[m]));
    for (std::size_t m = 0; m < nModes; ++m)
        if (xfit[m] < -1e-6 * mx) return out;
    out.x = xfit;
    return out;
}

/** What sn_pn_avg_rates rewrites in place; empty tables are left empty. */
template <class T>
struct SnPnAvgRates {
    Matrix<T> TN, AN, RN;
};

/**
 * Port of sn_pn_avg_rates: rewrite the place rows of TN, AN and RN so that
 * they agree with the recovered mode firing rates.
 *
 * A place's throughput becomes the rate its tokens are CONSUMED at, its
 * arrival rate the rate they are PRODUCED at, and its response time the
 * quotient the queue length and that throughput define. A model with no Place
 * node, or one whose firing rates could not be recovered, is returned
 * untouched.
 */
template <class T>
SnPnAvgRates<T> sn_pn_avg_rates(const qn::NetworkStruct<T>& sn, const Matrix<T>& QN,
                                const Matrix<T>& TN, const Matrix<T>& AN, const Matrix<T>& RN) {
    SnPnAvgRates<T> out;
    out.TN = TN;
    out.AN = AN;
    out.RN = RN;
    if (TN.rows() == 0) return out;
    bool has_place = false;
    for (std::size_t a = 0; a < sn.nodes.size(); ++a)
        if (sn.nodes[a].nodetype == qn::NodeType::Place) has_place = true;
    if (!has_place) return out;
    const SnPnFiringRates fr = sn_pn_firing_rates(sn, TN, false);
    if (fr.x.empty()) return out;
    const std::size_t R = sn.nclasses;
    for (std::size_t pp = 0; pp < fr.place_nodes.size(); ++pp) {
        const std::size_t ist = sn.nodes[fr.place_nodes[pp] - 1].station;
        if (ist == 0) continue;
        for (std::size_t k = 0; k < R; ++k) {
            double tk = 0.0, ak = 0.0;
            for (std::size_t m = 0; m < fr.x.size(); ++m) {
                tk += fr.consumed[m][pp * R + k] * fr.x[m];
                ak += fr.produced[m][pp * R + k] * fr.x[m];
            }
            out.TN(ist - 1, k) = num_traits<T>::from_double(tk);
            if (out.AN.rows() != 0) out.AN(ist - 1, k) = num_traits<T>::from_double(ak);
            if (out.RN.rows() != 0)
                out.RN(ist - 1, k) =
                    tk > 0.0 ? T(QN(ist - 1, k) / num_traits<T>::from_double(tk))
                             : num_traits<T>::from_int(0);
        }
    }
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_PN_FIRING_RATES_H
