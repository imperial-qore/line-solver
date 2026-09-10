/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_MAP_MODULATION_H
#define LINE_API_SN_SN_MAP_MODULATION_H

/**
 * Port of matlab/src/api/sn/sn_map_modulation.m.
 *
 * Every MAP-like process in the model, collected as the modulating environment
 * it really is: `map2renv` turns this list into a random-environment model, in
 * which the MAP's phase process is the environment stage and the network runs
 * at the rate that stage selects.
 *
 * MARKED PROCESSES. When a station's arrival is a marked MAP (an MMAP), the
 * classes it marks share ONE modulating phase process, so they are reported as
 * a SINGLE entry carrying (D0, D1^(1), ..., D1^(C)) over the marked class set,
 * not as one entry per class. Splitting them would give each class its own
 * independent environment and lose exactly the correlation the MMAP encodes.
 *
 * `is_mmpp` records whether every mark block is diagonal: a diagonal D1 means
 * an arrival never changes the phase, which is the MMPP special case whose
 * environment is the phase process alone.
 *
 * ARITHMETIC: field. A Frobenius norm comparison, no transcendental.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** One modulating process: which station and classes it drives, and its blocks. */
template <class T>
struct SnMapModulation {
    std::size_t ist = 0;   ///< 1-based station
    std::size_t node = 0;  ///< 1-based node
    bool arrival = false;  ///< true at a Source, false for a service process
    std::vector<std::size_t> classes;  ///< 1-based classes this process drives
    Matrix<T> D0;
    std::vector<Matrix<T>> D1;  ///< one block per entry of `classes`
    std::size_t order = 0;      ///< phases of the modulating process
    bool is_mmpp = false;       ///< every D1 block is diagonal
};

namespace detail {

/** `norm(D1 - diag(diag(D1)),'fro') <= Zero * max(1, norm(D1,'fro'))`. */
template <class T>
bool sn_map_is_diagonal(const Matrix<T>& D1) {
    double off = 0.0, all = 0.0;
    for (std::size_t a = 0; a < D1.rows(); ++a)
        for (std::size_t b = 0; b < D1.cols(); ++b) {
            const double v = num_traits<T>::to_double(D1(a, b));
            all += v * v;
            if (a != b) off += v * v;
        }
    return std::sqrt(off) <= lang::GlobalConstants::Zero * std::max(1.0, std::sqrt(all));
}

}  // namespace detail

template <class T>
std::vector<SnMapModulation<T>> sn_map_modulation(const qn::NetworkStruct<T>& sn) {
    std::vector<SnMapModulation<T>> mods;
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t nd = sn.station_to_node[ist - 1];
        const bool arrival =
            nd != 0 && sn.nodes[nd - 1].nodetype == qn::NodeType::Source;
        const std::vector<std::size_t>& marks = sn.stations[ist - 1].marked_classes;
        std::vector<bool> done(sn.nclasses, false);
        for (std::size_t r = 1; r <= sn.nclasses; ++r) {
            if (done[r - 1]) continue;
            const lang::ProcessType pt = sn.service[ist - 1][r - 1].type;
            if (!(pt == lang::ProcessType::MAP || pt == lang::ProcessType::MMPP2 ||
                  pt == lang::ProcessType::MMAP))
                continue;
            const lang::Distrib<T>& d = sn.service[ist - 1][r - 1];
            if (d.disabled || d.D0.rows() == 0) continue;
            bool is_marked = false;
            for (std::size_t c : marks)
                if (c == r) is_marked = true;
            SnMapModulation<T> m;
            m.ist = ist;
            m.node = nd;
            m.arrival = arrival;
            if (is_marked && !marks.empty()) {
                // one entry for the whole marked set, carried by the block list
                // of the FIRST marked class -- the reference's "carrier"
                const lang::Distrib<T>& carrier = sn.service[ist - 1][marks[0] - 1];
                if (carrier.Dmark.size() < marks.size())
                    throw InputError(
                        "sn_map_modulation: the marked arrival process at station " +
                        std::to_string(ist) + " carries " +
                        std::to_string(carrier.Dmark.size()) + " mark matrices for " +
                        std::to_string(marks.size()) +
                        " marked classes; the (D0,D1,D1^(1),...,D1^(C)) form is required");
                m.classes = marks;
                m.D0 = carrier.D0;
                m.is_mmpp = true;
                for (std::size_t k = 0; k < marks.size(); ++k) {
                    m.D1.push_back(carrier.Dmark[k]);
                    if (!detail::sn_map_is_diagonal(carrier.Dmark[k])) m.is_mmpp = false;
                    done[marks[k] - 1] = true;
                }
                m.order = carrier.D0.rows();
            } else {
                m.classes.push_back(r);
                m.D0 = d.D0;
                m.D1.push_back(d.D1);
                m.order = d.D0.rows();
                m.is_mmpp = detail::sn_map_is_diagonal(d.D1);
                done[r - 1] = true;
            }
            mods.push_back(m);
        }
    }
    return mods;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_MAP_MODULATION_H
