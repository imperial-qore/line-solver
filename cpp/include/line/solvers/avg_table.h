/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AVG_TABLE_H
#define LINE_SOLVERS_AVG_TABLE_H

/**
 * @file
 * @ingroup line_solvers
 * @ingroup line_public
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
#include <iosfwd>
#include <string>
#include <vector>

namespace line {

class JobClass;
class Node;

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
    /** A complete metric column. */
    std::vector<double> column(const std::string& column) const;

    /** Rows whose station or class has `name`, MATLAB's one-argument `filterBy`. */
    AvgTable filter_by(const std::string& name) const;
    /** Rows at one station for one class, MATLAB's two-argument `filterBy`. */
    AvgTable filter_by(const std::string& station, const std::string& jobclass) const;
    AvgTable filter_by(const Node& node) const;
    AvgTable filter_by(const ::line::JobClass& jobclass) const;
    AvgTable filter_by(const Node& node, const ::line::JobClass& jobclass) const;
    AvgTable filter_by(const ::line::JobClass& jobclass, const Node& node) const;

    /** MATLAB-compatible filtering aliases, preserving C++ snake-case naming. */
    AvgTable get(const Node& node) const;
    AvgTable get(const ::line::JobClass& jobclass) const;
    AvgTable get(const Node& node, const ::line::JobClass& jobclass) const;
    AvgTable get(const ::line::JobClass& jobclass, const Node& node) const;
    AvgTable tget(const Node& node) const;
    AvgTable tget(const ::line::JobClass& jobclass) const;
    AvgTable tget(const Node& node, const ::line::JobClass& jobclass) const;
    AvgTable tget(const ::line::JobClass& jobclass, const Node& node) const;

    /** Direct object indexing: `table(queue, jobs)`. */
    AvgTable operator()(const Node& node) const;
    AvgTable operator()(const ::line::JobClass& jobclass) const;
    AvgTable operator()(const Node& node, const ::line::JobClass& jobclass) const;
    AvgTable operator()(const ::line::JobClass& jobclass, const Node& node) const;

    /** Print the labelled station-class table. */
    void print(std::ostream& out) const;
    void print() const;
    /** The number of rows the table carries. */
    std::size_t size() const { return Station.size(); }
    bool empty() const { return Station.empty(); }
};

std::ostream& operator<<(std::ostream& out, const AvgTable& table);

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

/**
 * `getSymbolicSolution`: the stationary law as a function of the rate symbols.
 *
 * The entries are expression strings, not numbers, and are NOT comparable with
 * another codebase's by text: the symbol numbering x1..xE follows event
 * enumeration order and the printed normal form depends on the engine version.
 * Substitute rates and compare numbers instead.
 */
struct SymbolicSolution {
    std::vector<std::string> pi;       ///< stationary probability of each state
    std::vector<std::string> num;      ///< numerator of each entry over `den`
    std::string den = "1";             ///< common denominator of the vector
    std::vector<std::string> symbols;  ///< x1..xE, empty where an event has no positive rate
    std::vector<double> rate0;         ///< nominal value of each symbol
    std::string engine;                ///< backend that answered, e.g. "sage"
};

/** `getTranAvg`: the transient mean queue length per (station, class). */
struct TranAvg {
    std::vector<double> t;                    ///< the time axis
    std::vector<std::vector<double> > QNt;    ///< [step][station*class]
    std::vector<std::string> label;           ///< the (station, class) of each column
};

/**
 * `sampleSysAggr` / `sampleAggr`: ONE simulated trajectory, not a mean.
 *
 * It has the shape of a `TranAvg` and is deliberately a type of its own,
 * because the two must never be read for each other: `TranAvg` is E[N](t)
 * estimated over independent runs, this is a single sample path, and averaging
 * a column of `state` over `t` is a time average rather than an ensemble one.
 */
struct SamplePath {
    std::vector<double> t;                      ///< event times, ascending
    std::vector<std::vector<double> > state;    ///< [step][column]
    std::vector<std::string> label;             ///< what each column counts
};

}  // namespace line

#endif  // LINE_SOLVERS_AVG_TABLE_H
