/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_EXAMPLES_EXAMPLE_NODE_TABLE_H
#define LINE_EXAMPLES_EXAMPLE_NODE_TABLE_H

/**
 * `getAvgNodeTable()` and `getAvgChainTable()` for a ported example.
 *
 * NOT IN examples_common.h, AND DELIBERATELY SO. This needs the station->node
 * scatter, which needs `mva::AvgResult`, which instantiates the solver template
 * stack; examples_common.h is included by every one of the ~35 example
 * translation units and only the handful that print a node table should pay for
 * that. Those already include the solver stack for their own solve.
 *
 * WHY AN EXAMPLE WOULD WANT THIS TABLE. `getAvgTable` is indexed by STATION,
 * and a cache model's answer is not at a station: the hit and miss rates live at
 * the Cache node and at the ClassSwitch below it, neither of which is one. The
 * four `cache_replc_*` references call `avg_node_table()` for that reason, and a
 * twin that prints the station table instead drops every row that carries the
 * result -- 71 golden cells across the four.
 *
 * The scatter itself is `line::solvers::node_metrics`, which is also what
 * `line-cli -a node` and `-a nodechain` call, so an example and the CLI cannot
 * report different node tables for one model.
 */

#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/solvers/solver_chain_tables.h"
#include "line/solvers/solver_node_tables.h"

namespace line {
namespace examples {

/** The header of the node table, which names its first column `Station`. */
inline void node_header() {
    // `Station`, NOT `Node`, and it is not a mislabel: every consumer of these
    // tables -- the recorder's own labels, the shared goldens' `Station` key --
    // uses one name for the first label column, and the reference's own
    // `getAvgNodeTable` prints its index under whatever the printer calls it.
    // Renaming the column here would fork the two label spaces for one table.
    std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen",
                "Util", "RespT", "ResidT", "ArvR", "Tput");
}

/**
 * `getAvgNodeTable()`: one row per (node, class) that carries a metric.
 *
 * The all-zero row is dropped exactly as `avg_rows` drops it and as the
 * reference drops it -- a class that never reaches a node has no row, rather
 * than a row of zeros.
 */
template <class T, class Res>
void print_avg_node(const qn::NetworkStruct<T>& sn, const Res& r) {
    const solvers::NodeMetrics<T> m = solvers::node_metrics<T>(sn, r);
    const std::size_t I = sn.nodes.size(), R = sn.nclasses;
    namespace parity = line::examples::parity;
    node_header();
    if (parity::enabled()) parity::begin_table("avg", "Station", "JobClass");
    auto emit = [&](const std::string& name, const std::string& cls, double q, double u, double rr,
                    double w, double a, double t) {
        std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n", name.c_str(),
                    cls.c_str(), q, u, rr, w, a, t);
        if (!parity::enabled()) return;
        std::vector<parity::Cell> cells;
        cells.push_back(parity::Cell{"QLen", q});
        cells.push_back(parity::Cell{"Util", u});
        cells.push_back(parity::Cell{"RespT", rr});
        cells.push_back(parity::Cell{"ResidT", w});
        cells.push_back(parity::Cell{"ArvR", a});
        cells.push_back(parity::Cell{"Tput", t});
        parity::add_row(name, cls, cells);
    };
    for (std::size_t i = 0; i < I; ++i)
        for (std::size_t c = 0; c < R; ++c) {
            const double q = num_traits<T>::to_double(m.QN(i, c));
            const double u = num_traits<T>::to_double(m.UN(i, c));
            const double rr = num_traits<T>::to_double(m.RN(i, c));
            const double w = num_traits<T>::to_double(m.WN(i, c));
            const double a = num_traits<T>::to_double(m.AN(i, c));
            const double t = num_traits<T>::to_double(m.TN(i, c));
            if (q == 0.0 && u == 0.0 && rr == 0.0 && w == 0.0 && a == 0.0 && t == 0.0) continue;
            emit(sn.nodes[i].name, sn.classes[c].name, q, u, rr, w, a, t);
        }
    // A FINITE CAPACITY REGION IS NOT A NODE, and the reference still prints it
    // in this table: a WAITQ region holds jobs that are in no station's QLen, and
    // the model's population only balances once they are read. The rows ride past
    // the stations in the returned result -- the (M+F) layout both wrappers
    // report -- so a solver that does not measure a region contributes none.
    //
    // Util AND ArvR ARE NaN, NOT ZERO. A region has no server to be busy and no
    // arrival process of its own; the reference reports both as missing, and a
    // zero there would be a number the run never measured.
    const std::size_t F =
        r.QN.rows() > sn.nstations ? r.QN.rows() - sn.nstations : static_cast<std::size_t>(0);
    const double region_nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t c = 0; c < R; ++c) {
            const double q = num_traits<T>::to_double(r.QN(sn.nstations + f, c));
            const double rr = num_traits<T>::to_double(r.RN(sn.nstations + f, c));
            const double w = num_traits<T>::to_double(r.WN(sn.nstations + f, c));
            const double t = num_traits<T>::to_double(r.TN(sn.nstations + f, c));
            // The reference's own filter, the region twin of the all-zero row
            // test: a region no job ever entered is absent rather than a row of
            // zeros.
            if (!(q > 0.0 || rr > 0.0 || t > 0.0)) continue;
            const std::string name = (f < sn.regions.size() && !sn.regions[f].name.empty())
                                         ? sn.regions[f].name
                                         : "FCR" + std::to_string(f + 1);
            emit(name, sn.classes[c].name, q, region_nan, rr, w, region_nan, t);
        }
}

/**
 * `getAvgNodeTable()` for a WRAPPER solve, which reports the six matrices flat.
 *
 * JMT and LDES return their tables through `WrapperAvg` rather than an
 * `mva::AvgResult`, so the scatter has no refreshed struct and no cache split to
 * read -- neither wrapper reports one, and neither is ever the engine behind a
 * cache model in this corpus.
 */
inline void print_avg_node(const qn::NetworkStruct<double>& sn, const WrapperAvg& w) {
    mva::AvgResult<double> r;
    r.QN = w.QN;
    r.UN = w.UN;
    r.RN = w.RN;
    r.WN = w.WN;
    r.AN = w.AN;
    r.TN = w.TN;
    print_avg_node<double>(sn, r);
}

/**
 * `getAvgChainTable()`: the station table aggregated by CHAIN.
 *
 * EVERY ROW IS EMITTED, including the all-zero ones, unlike the AvgTable and the
 * AvgNodeTable. The reference builds this table on a full (M x C) grid with no
 * row filter, and a chain that is absent from a station is information -- it is
 * the shape of the routing -- where a class absent from a station in the
 * AvgTable is only that class's own scope. `line-cli -a chain` prints it the
 * same way, from the same `solver_get_avg_chain`.
 *
 * Chains are named `Chain1..ChainC`, which is what `solvers::chain_names` gives
 * and what the goldens key on: a chain has no name of its own in the model.
 */
template <class T, class Res>
void print_avg_chain(const qn::NetworkStruct<T>& sn, const Res& r) {
    const solvers::ChainResult<T> t = solvers::solver_get_avg_chain<T>(sn, r);
    const std::vector<std::string> chains = solvers::chain_names(sn.nchains);
    const std::vector<std::string> classes = solvers::chain_class_labels<T>(sn);
    namespace parity = line::examples::parity;
    std::printf("%-16s %-10s %-20s %12s %12s %12s %12s %12s %12s\n", "Station", "Chain",
                "JobClasses", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
    if (parity::enabled()) parity::begin_table("chain", "Station", "JobClass");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < chains.size(); ++c) {
            const double q = num_traits<T>::to_double(t.QN(i, c));
            const double u = num_traits<T>::to_double(t.UN(i, c));
            const double rr = num_traits<T>::to_double(t.RN(i, c));
            const double w = num_traits<T>::to_double(t.WN(i, c));
            const double a = num_traits<T>::to_double(t.AN(i, c));
            const double tp = num_traits<T>::to_double(t.TN(i, c));
            std::printf("%-16s %-10s %-20s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
                        sn.stations[i].name.c_str(), chains[c].c_str(), classes[c].c_str(), q, u,
                        rr, w, a, tp);
            if (!parity::enabled()) continue;
            std::vector<parity::Cell> cells;
            cells.push_back(parity::Cell{"QLen", q});
            cells.push_back(parity::Cell{"Util", u});
            cells.push_back(parity::Cell{"RespT", rr});
            cells.push_back(parity::Cell{"ResidT", w});
            cells.push_back(parity::Cell{"ArvR", a});
            cells.push_back(parity::Cell{"Tput", tp});
            parity::add_row(sn.stations[i].name, chains[c], cells);
        }
}

}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_EXAMPLE_NODE_TABLE_H
