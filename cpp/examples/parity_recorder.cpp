/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#include "parity_recorder.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>

#include "json.hpp"

namespace line {
namespace examples {
namespace parity {
namespace {

using json = nlohmann::json;

/**
 * The whole recorder state, function-local so it is initialised on first use.
 *
 * A twin runs in one process and prints in one thread, so nothing here is
 * synchronised; the sequence number is the call order the other recorders also
 * report, which is what lets a consumer keep the FIRST of several models.
 */
struct State {
    bool on = false;
    std::string path;
    std::string solver;
    std::string method = "default";
    long seq = 0;
    std::vector<Record> records;
    std::vector<std::string> notes;
    /// Index of the table rows are currently appended to, or -1 when none.
    long open = -1;
    /// Index of the open DERIVED table per solver key; see `add_derived`.
    std::map<std::string, long> derived_tables;
    /// True when `solver` came from a result table rather than from `section()`.
    bool solver_from_table = false;
};

State& state() {
    static State s;
    return s;
}

}  // namespace

bool enabled() { return state().on; }

void enable(const std::string& path) {
    State& s = state();
    s.on = true;
    s.path = path;
    s.records.clear();
    s.notes.clear();
    s.seq = 0;
    s.open = -1;
    s.solver.clear();
    s.method = "default";
    s.derived_tables.clear();
    s.solver_from_table = false;
}

void set_solver(const std::string& name) {
    State& s = state();
    if (!s.on) return;
    s.solver = name;
    s.solver_from_table = false;
    // A NEW SOLVER CLOSES THE OPEN TABLE. `section()` is the twin's own
    // declaration that what follows belongs to a different solver, so carrying
    // rows across it would file one solver's answer under another's name.
    s.method = "default";
    s.open = -1;
}

void set_method(const std::string& method) {
    State& s = state();
    if (!s.on) return;
    s.method = method.empty() ? std::string("default") : method;
}

const std::string& current_solver() { return state().solver; }

void set_solver_from_table(const std::string& name) {
    State& s = state();
    if (!s.on || name.empty()) return;
    if (!s.solver.empty() && !s.solver_from_table) return;
    if (s.solver == name) return;
    s.solver = name;
    s.solver_from_table = true;
    s.method = "default";
    s.open = -1;
}

void begin_table(const std::string& view, const std::string& row_label,
                 const std::string& col_label) {
    State& s = state();
    if (!s.on || s.solver.empty()) return;
    Record rec;
    rec.solver = s.solver;
    rec.method = s.method;
    rec.view = view;
    rec.labels.push_back(row_label);
    rec.labels.push_back(col_label);
    rec.seq = s.seq++;
    rec.derived = (view == "scalar");
    s.records.push_back(rec);
    s.open = static_cast<long>(s.records.size()) - 1;
}

void add_row(const std::string& row_label, const std::string& col_label,
             const std::vector<Cell>& cells) {
    State& s = state();
    if (!s.on || s.open < 0) return;
    Record::Row row;
    row.row_label = row_label;
    row.col_label = col_label;
    row.cells = cells;
    s.records[static_cast<std::size_t>(s.open)].rows.push_back(row);
}

void add_scalar(const std::string& label, double value) {
    State& s = state();
    if (!s.on) return;
    // A scalar belongs to whichever solver last declared itself; where none has,
    // it is the example's own arithmetic and is filed under the reserved key
    // the derived goldens use.
    const std::string owner = s.solver.empty() ? std::string("DERIVED") : s.solver;
    Record rec;
    rec.solver = owner;
    rec.method = s.method;
    rec.view = "scalar";
    rec.labels.push_back("Quantity");
    rec.labels.push_back("Index");
    rec.seq = s.seq++;
    rec.derived = true;
    Record::Row row;
    row.row_label = label;
    row.col_label = "0";
    Cell cell;
    cell.metric = "QLen";
    cell.value = value;
    row.cells.push_back(cell);
    rec.rows.push_back(row);
    s.records.push_back(rec);
    // A scalar does not open a table: the next `add_row` belongs to whatever
    // table was open before it, not to this one-row record.
}

void add_derived(const std::string& solver, const std::string& row_label,
                 const std::string& col_label, double value, const std::string& metric) {
    State& s = state();
    if (!s.on) return;
    // ONE TABLE PER SOLVER KEY, appended to, so an example may declare its keys
    // as it prints them. It does NOT touch `open`: a derived quantity may be
    // printed between two rows of an ordinary table, and stealing the cursor
    // would file the rest of that table under this one.
    std::map<std::string, long>::iterator it = s.derived_tables.find(solver);
    if (it == s.derived_tables.end()) {
        Record rec;
        rec.solver = solver;
        rec.method = s.method;
        rec.view = "derived";
        rec.labels.push_back("Station");
        rec.labels.push_back("JobClass");
        rec.seq = s.seq++;
        // COMPARED AT THE GOLDEN'S OWN WRITTEN PRECISION. A derived cell was
        // written by the example's format string rather than by a result table,
        // so the half-unit the comparator allows is a DECIMAL place and not a
        // significant figure.
        rec.derived = true;
        s.records.push_back(rec);
        it = s.derived_tables.insert(std::make_pair(
            solver, static_cast<long>(s.records.size()) - 1)).first;
    }
    Record::Row row;
    row.row_label = row_label;
    row.col_label = col_label;
    Cell cell;
    cell.metric = metric;
    cell.value = value;
    row.cells.push_back(cell);
    s.records[static_cast<std::size_t>(it->second)].rows.push_back(row);
}

void note_refusal(const std::string& solver, const std::string& why) {
    State& s = state();
    if (!s.on) return;
    s.notes.push_back(solver + ": " + why);
}

const std::vector<Record>& records() { return state().records; }

const std::vector<std::string>& notes() { return state().notes; }

bool dump() {
    State& s = state();
    if (!s.on) return true;
    json out;
    out["records"] = json::array();
    for (std::size_t i = 0; i < s.records.size(); ++i) {
        const Record& rec = s.records[i];
        json jr;
        jr["solver"] = rec.solver;
        jr["method"] = rec.method;
        jr["view"] = rec.view;
        jr["labels"] = rec.labels;
        jr["seq"] = rec.seq;
        jr["derived"] = rec.derived;
        jr["rows"] = json::array();
        for (std::size_t k = 0; k < rec.rows.size(); ++k) {
            const Record::Row& row = rec.rows[k];
            json jrow;
            jrow[rec.labels[0]] = row.row_label;
            jrow[rec.labels[1]] = row.col_label;
            for (std::size_t c = 0; c < row.cells.size(); ++c) {
                const double v = row.cells[c].value;
                // A NON-FINITE CELL IS WRITTEN AS A STRING, not as bare `nan`,
                // which is not JSON and which every reader downstream would
                // reject as a malformed document rather than as a missing value.
                if (std::isnan(v))
                    jrow[row.cells[c].metric] = "NaN";
                else if (std::isinf(v))
                    jrow[row.cells[c].metric] = v > 0 ? "Infinity" : "-Infinity";
                else
                    jrow[row.cells[c].metric] = v;
            }
            jr["rows"].push_back(jrow);
        }
        out["records"].push_back(jr);
    }
    out["notes"] = s.notes;
    std::ofstream fh(s.path.c_str());
    if (!fh) {
        std::fprintf(stderr, "parity recorder: cannot write %s\n", s.path.c_str());
        return false;
    }
    fh << out.dump(1);
    return static_cast<bool>(fh);
}

}  // namespace parity
}  // namespace examples
}  // namespace line
