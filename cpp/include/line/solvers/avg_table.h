/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AVG_TABLE_H
#define LINE_SOLVERS_AVG_TABLE_H

/**
 * The result tables a solver returns.
 *
 * `AvgTable` is what Python's `getAvgTable()` returns and what this port spells
 * `avg_table()`: one row per (station, class) that carries a metric, plus the
 * per-class system columns `getAvgSysTable()` reports. The columns are the
 * reference's own -- QLen, Util, RespT, ResidT, ArvR, Tput -- and a
 * (station, class) pair with no row is a zero, not a gap.
 *
 * THESE TYPES CARRY NO SOLVER MACHINERY ON PURPOSE. The solver templates are a
 * heavy instantiation, so this header stays free of them and every example,
 * test and CLI translation unit can name a result without paying for the stack
 * that produced it; `line/solvers/solver.h` declares the solvers themselves and
 * their bodies are compiled once into `line_mp_api`.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace line {

/** `getAvgTable`, one row per (station, class) that carries a metric. */
struct AvgTable {
    std::vector<std::string> Station, JobClass;
    std::vector<double> QLen, Util, RespT, ResidT, ArvR, Tput;
    /** `getAvgSysTable`: system response time and throughput, per class. */
    std::vector<std::string> SysClass;
    std::vector<double> SysRespT, SysTput;
    std::string solver, method;
    int iter = 0;
    bool has_lognormconst = false;
    double lognormconst = 0.0;
    /** `getAvgCacheTable`'s ListCost column; empty on a model without item sizes. */
    std::vector<double> ListCost;
    /** The reference's own warning text, empty when it did not warn. */
    std::string warning;

    /** One cell of the table, by station and class NAME; NaN when absent. */
    double get(const std::string& column, const std::string& station,
               const std::string& jobclass) const;
    /** The column of a station over every class it has a row for. */
    std::vector<double> column(const std::string& column) const;
    /** The number of rows the table carries. */
    std::size_t size() const { return Station.size(); }
    bool empty() const { return Station.empty(); }
};

/** `SolverBA(model, method).getBoundsTable()`. */
struct BoundsTable {
    std::vector<std::string> Station, JobClass;
    std::vector<double> Qlower, Qupper, Tlower, Tupper;
    std::string method;
};

/** One response-time CDF curve: `F` the CDF value, `t` the time it is reached. */
struct CdfCurve {
    std::vector<double> F, t;
};

/** `getTranAvg`: the transient mean queue length per (station, class). */
struct TranAvg {
    std::vector<double> t;                    ///< the time axis
    std::vector<std::vector<double> > QNt;    ///< [step][station*class]
    std::vector<std::string> label;           ///< the (station, class) of each column
};

}  // namespace line

#endif  // LINE_SOLVERS_AVG_TABLE_H
