/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#ifndef LINE_EXAMPLES_PARITY_RECORDER_H
#define LINE_EXAMPLES_PARITY_RECORDER_H

#include <string>
#include <vector>

/**
 * @file parity_recorder.h
 * @brief Capture what a twin computed, with the solver that computed it.
 *
 * WHY THIS EXISTS. Cross-codebase parity is asserted against one shared golden
 * per example (`goldens/baselines/*.json`), keyed by SOLVER NAME. Until
 * 2026-08-19 the only way to recover that key was to scrape the banner a twin
 * printed above each table -- one regex dialect per codebase, and a value
 * truncated to the six digits `%12.6g` shows. The recorder supplies the same
 * attribution BY CONSTRUCTION and at full precision. It is the C++ twin of
 * `python/line_solver/result_recorder.py` and
 * `matlab/src/io/LineResultRecorder.m`.
 *
 * WHAT MAKES IT SMALL. The twins already declare the solver: `section("MVA")`
 * exists precisely to say which solver the following table belongs to, and
 * `avg_rows` is the ONE place every average table is materialised. So four hook
 * points cover the whole corpus -- `section`, `na`, `avg_rows` and the two
 * printers that do not go through it -- and not one of the 388 registered twins
 * changes.
 *
 * It is OFF unless `--record <path>` asked for it, and then it costs one
 * `enabled()` test per printed row.
 */
namespace line {
namespace examples {
namespace parity {

/** One cell of a recorded table: a metric name and its value. */
struct Cell {
    std::string metric;
    double value;
};

/** One recorded table, with what produced it. */
struct Record {
    std::string solver;   ///< the golden's key: what `section()` declared
    std::string method;   ///< the method that solver resolved, or "default"
    std::string view;     ///< which table this is: "avg", "sys", "scalar"
    std::vector<std::string> labels;   ///< the two label column names
    /// One row per (station, class): the two labels then its metric cells.
    struct Row {
        std::string row_label;
        std::string col_label;
        std::vector<Cell> cells;
    };
    std::vector<Row> rows;
    long seq;             ///< call order, 0-based, as the other recorders number
    bool derived;         ///< compared at the golden's own written precision
};

/** True when a run asked to record. Nothing below costs anything when false. */
bool enabled();

/** Turn recording on and direct the dump at `path`. Called by `--record`. */
void enable(const std::string& path);

/**
 * Declare the solver the following tables belong to.
 *
 * `section()` calls this, so the declaration is the one the twin already makes
 * for its reader. A name that is not a golden key ("MVA (exact)") is kept
 * verbatim: the comparator reconciles spellings, and inventing a mapping here
 * would hide which twin produced what.
 */
void set_solver(const std::string& name);

/** Declare the method the current solver resolved, when the twin knows it. */
void set_method(const std::string& method);

/**
 * The solver `section()` last declared, or "" when none has (or when recording
 * is off). It is how a printer that derives a quantity files it under the solver
 * that produced the numbers behind it, instead of naming one by hand.
 */
const std::string& current_solver();

/**
 * Declare the solver from a RESULT TABLE that names its own, not from a banner.
 *
 * An `AvgTable` carries the solver that produced it, which is the same
 * attribution `section()` makes and is available even where a reference prints
 * no banner at all -- and a good part of the corpus prints none, because the
 * reference script is `print(MVA(model).getAvgTable())` and nothing more. Nine
 * examples recorded no table for exactly that reason.
 *
 * IT YIELDS TO `section()` AND NOT THE OTHER WAY ROUND. A twin that declared
 * `CTMC (cutoff 4)` or `LN(NC)` meant that name -- it is telling its reader
 * WHICH of several runs this is -- and a table's plain `CTMC` must not quietly
 * overwrite it. But it does replace a name IT set: `dispatch_open` prints an MVA
 * table and then an SSA one with no banner between them, and holding the first
 * name would file the simulation under the analytical solver's key.
 */
void set_solver_from_table(const std::string& name);

/** Open a new table under the current solver. Rows follow until the next open. */
void begin_table(const std::string& view, const std::string& row_label,
                 const std::string& col_label);

/** Add one row to the open table. */
void add_row(const std::string& row_label, const std::string& col_label,
             const std::vector<Cell>& cells);

/**
 * Record a labelled scalar the example DERIVED.
 *
 * Twenty-nine of the goldens hold a quantity no result table carries -- a state
 * probability, a phase-type moment, a cache hit rate. The twin prints them
 * through `print_scalar`, and they are recorded as one-row tables so everything
 * downstream handles a single shape.
 */
void add_scalar(const std::string& label, double value);

/**
 * Record a derived quantity under the GOLDEN'S OWN (row, column) key.
 *
 * `add_scalar` files a value under ('<label>', '0'), which is the shape the
 * Python recorder produces from a hooked getter. EIGHTEEN OF THE GOLDENS ARE
 * NOT SOLVER TABLES: they hold a workflow's phase-type moments, a state
 * probability, a cache hit rate, a reward component, a mean integrated off a
 * CDF -- and each keys those by a pair the EXAMPLE chooses, ('PH', 'mean'),
 * ('OptResult', 'Value'), ('QueueLength', 'Reward'). Python reassembles that
 * shape afterwards, in `python/tests/parity/derived.py`, by reading the lines
 * the example printed: no hook reaches an example's own arithmetic there. In
 * C++ the example IS the code, so it names the key directly and nothing has to
 * be scraped back out of stdout.
 *
 * `solver` is the golden's key for the quantity -- a solver's name where one
 * produced it ('CTMC'), or the shape's own key where none did ('WF', 'OPT',
 * 'CACHE', 'DEP'). Rows accumulate into ONE table per solver, so an example may
 * declare its keys one at a time as it prints them.
 *
 * `metric` is the column the golden holds the value in; every derived golden
 * but `cdf_respt_populations` uses `QLen`, which is why it is the default.
 */
void add_derived(const std::string& solver, const std::string& row_label,
                 const std::string& col_label, double value,
                 const std::string& metric = "QLen");

/**
 * Record a refusal the twin made BY NAME (`na("MAM", "...")`).
 *
 * A refusal is a fact about the port and is a named skip; a table that simply
 * never arrived is a failure. The two must not look alike downstream, which is
 * why this is recorded rather than only printed.
 */
void note_refusal(const std::string& solver, const std::string& why);

/** Everything recorded so far, in call order. */
const std::vector<Record>& records();

/** The refusals recorded so far, one sentence each. */
const std::vector<std::string>& notes();

/** Write the JSON dump to the path `enable()` was given. No-op when off. */
bool dump();

}  // namespace parity
}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_PARITY_RECORDER_H
