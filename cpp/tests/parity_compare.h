/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#ifndef LINE_TESTS_PARITY_COMPARE_H
#define LINE_TESTS_PARITY_COMPARE_H

#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "json.hpp"

/**
 * @file parity_compare.h
 * @brief Score one recorded row against its shared golden.
 *
 * The golden is keyed by SOLVER NAME, then by `(Station, JobClass)`, then by
 * metric. `line-examples --record` supplies the same shape, so the comparison is
 * a join on names -- no text is parsed anywhere. This is the C++ twin of
 * `python/tests/parity/compare.py` and
 * `line-test.git/test/testsParity/ParityCompare.m`, and the rules below were
 * learned from real false passes:
 *
 *   - ONLY AN ENSEMBLE RECORD MAY BE ADOPTED onto a golden key it does not
 *     match. A recorded table always names its solver, so the adoption path that
 *     existed for the scraper is gone.
 *   - A ROW THAT COMPARED NO CELL IS NOT A PASS. `measure_row` returns the
 *     count, and the caller must treat zero as a named skip or a failure.
 *   - A solver the engine REFUSED BY NAME is excused; a solver that merely
 *     failed to appear is a failure.
 */
namespace line {
namespace tests {
namespace parity {

using json = nlohmann::json;

/** The metrics a golden may carry, in the order the goldens list them. */
inline const std::vector<std::string>& metrics() {
    static const std::vector<std::string> m = {"QLen", "Util",  "RespT",
                                               "ResidT", "ArvR", "Tput"};
    return m;
}

/**
 * THE SHARED GOLDENS CARRY FIVE SIGNIFICANT DIGITS AND THE RECORDER RETURNS ALL
 * OF THEM. The goldens were generated from printed tables, where a cell shows
 * five digits (0.45505) and the recorder returns 0.4550451573262394. Compared
 * directly that is 1.1e-5 relative -- eleven times the deterministic tolerance
 * -- so every deterministic cell in the corpus would fail on ROUNDING.
 *
 * The fix is NOT a looser tolerance, which would hide real defects of the same
 * size. It is that A GOLDEN IS ONLY KNOWN TO ITS LAST WRITTEN DIGIT.
 */
static const int GOLDEN_SIGFIGS = 5;

/** How many significant digits this golden value actually carries. */
inline int sigfigs(double value) {
    if (value == 0.0 || std::isnan(value) || std::isinf(value)) return GOLDEN_SIGFIGS;
    char buf[64];
    for (int k = 1; k <= 17; ++k) {
        std::snprintf(buf, sizeof(buf), "%.*g", k, value);
        if (std::strtod(buf, nullptr) == value) return k;
    }
    return 17;
}

/** Decimal places in the value's shortest round-tripping form; -1 if exponential. */
inline int decimals(double value) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", value);
    std::string text(buf);
    for (int k = 1; k <= 17; ++k) {
        std::snprintf(buf, sizeof(buf), "%.*g", k, value);
        if (std::strtod(buf, nullptr) == value) {
            text = buf;
            break;
        }
    }
    if (text.find('e') != std::string::npos || text.find('E') != std::string::npos) return -1;
    const std::size_t dot = text.find('.');
    return dot == std::string::npos ? 0 : static_cast<int>(text.size() - dot - 1);
}

/**
 * Half a unit in the last digit the golden actually carries.
 *
 * For an ordinary table cell the printer's precision is five significant digits,
 * so the unit comes from the value's MAGNITUDE and not from how many digits it
 * happens to show: a golden written `1.0` came from a printer that would have
 * shown `1.0004`, so it is known to 5e-5 and not to 0.5.
 *
 * A DERIVED cell was written by the example's own format string rather than by a
 * result table, and there the unit is a DECIMAL place.
 *
 * An exact zero gets nothing: `0.0` is what a printer shows for a cell that
 * really is zero, and a half-unit allowance there would absorb an actual 0.04.
 */
inline double golden_slack(double expected, bool derived) {
    if (expected == 0.0 || std::isnan(expected) || std::isinf(expected)) return 0.0;
    if (derived) {
        const int places = decimals(expected);
        return places < 0 ? 0.0 : 0.5 * std::pow(10.0, -places);
    }
    const int digits = std::max(GOLDEN_SIGFIGS, sigfigs(expected));
    const double exponent = std::floor(std::log10(std::fabs(expected)));
    return 0.5 * std::pow(10.0, exponent - digits + 1);
}

/**
 * `MVA` for both `MVA` and `MVA:schmidt`: the family a golden key names.
 *
 * A golden's solver key is either a bare FAMILY, which means that family's
 * DEFAULT path, or `FAMILY:method`, which names one runnable method of it. The
 * two classify alike -- a swept `SSA:nrm` is no less a simulator for being
 * written with a colon, and reading the whole key would hold it to the 1e-6 a
 * closed-form solver gets.
 */
inline std::string solver_family(const std::string& solver) {
    const std::size_t colon = solver.find(':');
    return colon == std::string::npos ? solver : solver.substr(0, colon);
}

/** The tolerance policy, read from `goldens/tolerance_policy.json`. */
class Policy {
  public:
    explicit Policy(const std::string& golden_root) {
        std::ifstream fh((golden_root + "/tolerance_policy.json").c_str());
        if (!fh) throw std::runtime_error("cannot read tolerance_policy.json");
        fh >> doc_;
        std::ifstream th((golden_root + "/test_tolerances.json").c_str());
        if (th) th >> per_test_;
    }

    /**
     * (rel, abs) for an (example, solver) pair; a negative value means "unset".
     *
     * An override wins outright -- including one that sets only `abs`, which
     * then leaves `rel` unset and makes the absolute test the whole gate. That
     * is deliberate: `cdf_respt_populations` is goldened for an ODE truncation
     * residual, not for a relative agreement.
     *
     * AN OVERRIDE IS MATCHED ON THE EXACT KEY, a swept `FAMILY:method`
     * included, while the CLASS falls back to the FAMILY. A class says what
     * KIND of solver this is, which every method of a family shares; an
     * override is a measurement of ONE arm, and spreading it over the family's
     * other method names would widen gates nobody measured. To grant slack to one
     * method, name it in full: `"solver": "MVA:schmidt"`.
     */
    void tolerance_for(const std::string& example, const std::string& solver,
                       const std::string& solver_class, double* rel, double* abs_tol) const {
        *rel = -1.0;
        *abs_tol = -1.0;
        for (const json& entry : doc_.at("overrides")) {
            if (entry.value("example", std::string()) != example) continue;
            if (entry.value("solver", std::string()) != solver) continue;
            if (entry.contains("rel")) *rel = entry.at("rel").get<double>();
            if (entry.contains("abs")) *abs_tol = entry.at("abs").get<double>();
            return;
        }
        bool stochastic = solver_class == "stochastic";
        if (!stochastic)
            for (const json& n : doc_.at("stochasticSolvers"))
                if (n.get<std::string>() == solver_family(solver)) stochastic = true;
        const json& cls = doc_.at("classes").at(stochastic ? "stochastic" : "deterministic");
        if (cls.contains("rel")) *rel = cls.at("rel").get<double>();
        if (cls.contains("abs")) *abs_tol = cls.at("abs").get<double>();
    }

    /**
     * The tolerance with the generated per-example slack folded in.
     *
     * `test_tolerances.json` records the smallest absolute slack under which
     * every codebase matches that example's golden, so an example the codebases
     * agree on exactly carries 0.0 and grants nothing.
     */
    void effective_tolerance(const std::string& example, const std::string& solver,
                             const std::string& solver_class, double* rel,
                             double* abs_tol) const {
        tolerance_for(example, solver, solver_class, rel, abs_tol);
        double extra = 0.0;
        if (per_test_.is_object() && per_test_.contains(example) &&
            per_test_.at(example).is_number())
            extra = per_test_.at(example).get<double>();
        if (extra > 0.0) {
            extra = extra * (1.0 + 1e-9) + 1e-12;
            *abs_tol = (*abs_tol < 0.0) ? extra : std::max(*abs_tol, extra);
        }
    }

    /** (ok, observed difference) for one cell. NaN on either side is skipped. */
    bool values_match(double expected, double actual, double rel, double abs_tol,
                      double* diff) const {
        *diff = 0.0;
        if (std::isnan(expected) || std::isnan(actual)) return true;
        if (expected == actual) return true;
        double d = std::fabs(expected - actual);
        if (abs_tol >= 0.0 && d <= abs_tol) return true;
        const double denom = std::max(std::fabs(expected), std::fabs(actual));
        const double threshold = doc_.at("nearZero").at("threshold").get<double>();
        if (denom < threshold) {
            *diff = d;
            return d <= doc_.at("nearZero").at("abs").get<double>();
        }
        *diff = d / denom;
        return rel >= 0.0 && *diff <= rel;
    }

  private:
    json doc_;
    json per_test_;
};


/** One (station, class) key. */
typedef std::pair<std::string, std::string> Key;
/** One solver's cells, by key then by metric. */
typedef std::map<Key, std::map<std::string, double> > Table;

/** One recorded table, kept with the labels that say what its keys mean. */
struct View {
    std::vector<std::string> labels;
    Table table;
    std::string view;
    bool derived = false;
    /// The recorder's call order, which is what "the FIRST of several" means.
    long seq = 0;
};

/** A golden, in the shape the comparison wants. */
struct Golden {
    std::vector<std::string> solvers;          ///< keys AS WRITTEN
    std::vector<Table> tables;                 ///< aligned with `solvers`
    std::vector<std::string> classes;          ///< tolerance class per solver
};

/** A golden value that arrived as the wire's NaN spelling. */
inline double number_of(const json& v) {
    if (v.is_string()) {
        const std::string t = v.get<std::string>();
        if (t == "NaN") return std::numeric_limits<double>::quiet_NaN();
        if (t == "Infinity") return std::numeric_limits<double>::infinity();
        if (t == "-Infinity") return -std::numeric_limits<double>::infinity();
        return std::strtod(t.c_str(), nullptr);
    }
    return v.get<double>();
}

/** Read `goldens/baselines/<example>.json` into the shape above. */
inline Golden load_golden(const std::string& golden_root, const std::string& example) {
    std::ifstream fh((golden_root + "/baselines/" + example + ".json").c_str());
    if (!fh) throw std::runtime_error("no golden for " + example);
    json doc;
    fh >> doc;
    Golden g;
    for (json::const_iterator it = doc.at("solvers").begin(); it != doc.at("solvers").end();
         ++it) {
        Table table;
        for (const json& row : it.value()) {
            const Key key(row.value("Station", std::string()), row.value("JobClass", std::string()));
            std::map<std::string, double>& cells = table[key];
            for (std::size_t m = 0; m < metrics().size(); ++m)
                if (row.contains(metrics()[m])) cells[metrics()[m]] = number_of(row.at(metrics()[m]));
        }
        g.solvers.push_back(it.key());
        g.tables.push_back(table);
        std::string cls = "deterministic";
        if (doc.contains("class") && doc.at("class").contains(it.key()))
            cls = doc.at("class").at(it.key()).get<std::string>();
        g.classes.push_back(cls);
    }
    return g;
}

/**
 * The recorded document as views, by solver.
 *
 * VIEWS ARE KEPT APART rather than merged into one key space. An example that
 * asks a solver for both its cache table and its item table produces
 * ('Cache', 'InitClass') from the first and ('Cache', '1') ... from the second,
 * and merging them makes the station 'Cache' look as though it carried four job
 * classes -- which silently disables the default-class reconciliation the
 * golden's 'Jobs' row depends on.
 *
 * `first_wins` keeps the FIRST value where one solver supplies a key twice: an
 * example that builds SEVERAL models produces one table per model and its golden
 * holds the first (see `goldens/corpus.json`).
 */
inline std::map<std::string, std::vector<View> > records_to_views(const json& doc,
                                                                 bool first_wins) {
    std::map<std::string, std::vector<View> > out;
    if (!doc.contains("records")) return out;
    for (const json& rec : doc.at("records")) {
        View v;
        v.view = rec.value("view", std::string());
        v.derived = rec.value("derived", false);
        v.seq = rec.value("seq", static_cast<long>(0));
        for (const json& l : rec.at("labels")) v.labels.push_back(l.get<std::string>());
        if (v.labels.size() < 2) continue;
        for (const json& row : rec.at("rows")) {
            const Key key(row.value(v.labels[0], std::string()),
                          row.value(v.labels[1], std::string()));
            std::map<std::string, double>& cells = v.table[key];
            for (json::const_iterator it = row.begin(); it != row.end(); ++it) {
                if (it.key() == v.labels[0] || it.key() == v.labels[1]) continue;
                // A NON-FINITE CELL ARRIVES AS THE WIRE'S SPELLING ("NaN"), the
                // same one the goldens use, and is read back as the number. It
                // matters for the count rather than the verdict: a NaN on either
                // side is skipped by `values_match`, but a cell DROPPED here
                // would not be counted at all, and the count is what tells a row
                // that asserted nothing from one that agreed.
                if (it.value().is_number())
                    cells[it.key()] = it.value().get<double>();
                else if (it.value().is_string())
                    cells[it.key()] = number_of(it.value());
            }
        }
        if (v.table.empty()) continue;
        out[rec.value("solver", std::string())].push_back(v);
    }
    if (first_wins)
        for (std::map<std::string, std::vector<View> >::iterator it = out.begin();
             it != out.end(); ++it)
            std::reverse(it->second.begin(), it->second.end());
    return out;
}

/** The keys of a table whose station matches. */
inline std::vector<Key> keys_for_station(const Table& table, const std::string& station) {
    std::vector<Key> hits;
    for (Table::const_iterator it = table.begin(); it != table.end(); ++it)
        if (it->first.first == station) hits.push_back(it->first);
    return hits;
}

/** The two labels joined by a space, as a printed line would have shown them. */
inline std::string join_key(const Key& key) {
    std::string text = key.first + " " + key.second;
    while (!text.empty() && text[text.size() - 1] == ' ') text.erase(text.size() - 1);
    while (!text.empty() && text[0] == ' ') text.erase(0, 1);
    return text;
}

/**
 * Recover a golden key whose two labels were split at the wrong space.
 *
 * SOME GOLDENS RECORD WHERE THE SPACE FELL, NOT WHERE THE COLUMN DID: a model
 * whose station is named `Source 1` and whose class is `Class A` was written
 * down as `Station='Source 1 Class', JobClass='A'`. `oqn_cs_routing` is goldened
 * that way for all four of its solvers. The pair still identifies the same row,
 * because joining recovers the original line exactly, so a produced key is
 * renamed to the golden's spelling when the joined forms are equal AND unique on
 * both sides -- which makes it an identification and not a guess.
 */
inline void rejoin_split_labels(const std::vector<Key>& expected, Table* table) {
    std::map<std::string, std::vector<Key> > joined_exp, joined_act;
    for (std::size_t i = 0; i < expected.size(); ++i)
        joined_exp[join_key(expected[i])].push_back(expected[i]);
    for (Table::const_iterator it = table->begin(); it != table->end(); ++it)
        joined_act[join_key(it->first)].push_back(it->first);
    for (std::map<std::string, std::vector<Key> >::const_iterator it = joined_exp.begin();
         it != joined_exp.end(); ++it) {
        if (it->second.size() != 1) continue;
        std::map<std::string, std::vector<Key> >::const_iterator act = joined_act.find(it->first);
        if (act == joined_act.end() || act->second.size() != 1) continue;
        const Key& want = it->second[0];
        const Key& have = act->second[0];
        if (want == have || table->count(want)) continue;
        (*table)[want] = (*table)[have];
        table->erase(have);
    }
}

/**
 * This view's table with default-class labels reconciled against the golden.
 *
 * 'Jobs' IS A PLACEHOLDER FOR WHATEVER THE SECOND LABEL COLUMN HOLDS: a
 * single-class model is printed with the class named 'Jobs' by one codebase and
 * by its own name in another. THE UNIQUENESS REQUIREMENT IS WHAT MAKES THAT
 * SAFE -- renaming happens by STATION and only where that station has EXACTLY
 * ONE entry on the other side, so it can never merge two classes into one.
 */
inline Table reconcile(const std::vector<Key>& expected, const View& view) {
    Table out = view.table;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (expected[i].second != "Jobs") continue;
        const Key target(expected[i].first, "Jobs");
        if (out.count(target)) continue;
        const std::vector<Key> here = keys_for_station(out, expected[i].first);
        if (here.size() == 1) {
            out[target] = out[here[0]];
            out.erase(here[0]);
        }
    }
    std::vector<Key> present;
    for (Table::const_iterator it = out.begin(); it != out.end(); ++it)
        present.push_back(it->first);
    for (std::size_t i = 0; i < present.size(); ++i) {
        if (present[i].second != "Jobs") continue;
        std::vector<Key> there;
        for (std::size_t k = 0; k < expected.size(); ++k)
            if (expected[k].first == present[i].first) there.push_back(expected[k]);
        if (there.size() == 1 && !out.count(there[0])) {
            out[there[0]] = out[present[i]];
            out.erase(present[i]);
        }
    }
    rejoin_split_labels(expected, &out);
    return out;
}

/** One lookup table for a solver, from its views. Later views win. */
inline void flatten(const std::vector<Key>& expected, const std::vector<View>& views,
                    Table* out, std::set<Key>* derived) {
    out->clear();
    derived->clear();
    for (std::size_t v = 0; v < views.size(); ++v) {
        const Table table = reconcile(expected, views[v]);
        for (Table::const_iterator it = table.begin(); it != table.end(); ++it) {
            std::map<std::string, double>& cells = (*out)[it->first];
            for (std::map<std::string, double>::const_iterator c = it->second.begin();
                 c != it->second.end(); ++c)
                cells[c->first] = c->second;
            if (views[v].derived)
                derived->insert(it->first);
            else
                derived->erase(it->first);
        }
    }
}

/** ('LN', 'NC') for 'LN(NC)'; ('MVA', '') for 'MVA'. */
inline void split_ensemble(const std::string& name, std::string* base, std::string* member) {
    *base = name;
    member->clear();
    const std::size_t open = name.find('(');
    if (open != std::string::npos && !name.empty() && name[name.size() - 1] == ')') {
        *base = name.substr(0, open);
        *member = name.substr(open + 1, name.size() - open - 2);
    }
}

/**
 * The golden keys an ensemble record may be scored against.
 *
 * A layered solve is recorded as `LN(<layer>)`, while the goldens key it by
 * whatever the EXAMPLE printed. These are the same computation, so the spellings
 * may be reconciled -- but ONLY against the layer the row actually drove, which
 * the recorder knows and a scraper had to infer. THE BARE 'LN' IS A CANDIDATE
 * FOR EVERY MEMBER: a golden keyed plain 'LN' says only "the layered solve".
 */
inline std::vector<std::string> aliases(const std::string& base, const std::string& member) {
    std::vector<std::string> out;
    const char* none[] = {0};
    (void)none;
    if (base == "LN") {
        if (member == "NC") out = {"NC", "LN(NC)", "COMOM", "LN(COMOM)", "LN"};
        else if (member == "COMOM") out = {"COMOM", "LN(COMOM)", "NC", "LN(NC)", "LN"};
        else if (member == "MVA") out = {"MVA", "LN(MVA)", "LN"};
        else if (member == "FLD") out = {"FLD", "LN(FLD)", "LN(FLUID)", "LN"};
        else if (member == "SSA") out = {"SSA", "LN(SSA)", "LN"};
        else if (member == "CTMC") out = {"CTMC", "LN(CTMC)", "LN"};
        else if (member == "MAM") out = {"MAM", "LN(MAM)", "LN"};
        else if (member.empty()) out = {"LN"};
    } else if (base == "ENV") {
        if (member == "FLD") out = {"FLD", "ENV(FLD)", "ENV(FLUID)", "ENV"};
        else if (member == "CTMC") out = {"CTMC", "ENV(CTMC)", "ENV"};
        else if (member == "MVA") out = {"MVA", "ENV(MVA)", "ENV"};
        else if (member.empty()) out = {"ENV"};
    } else if (base == "UQ") {
        if (member == "MVA") out = {"MVA", "UQ(MVA)", "UQ"};
        else if (member == "NC") out = {"NC", "UQ(NC)", "UQ"};
        else if (member == "CTMC") out = {"CTMC", "UQ(CTMC)", "UQ"};
        else if (member == "FLD") out = {"FLD", "UQ(FLD)", "UQ"};
        else if (member == "SSA") out = {"SSA", "UQ(SSA)", "UQ"};
        else if (member.empty()) out = {"UQ"};
    }
    return out;
}

/**
 * Reconcile recorded solver labels against the golden's keys.
 *
 * Exact names match first. Then each ensemble record is offered to the ONE
 * golden key that names the member it drove, if exactly one such key is left
 * unmatched. Nothing else is adopted: a record that names its solver is that
 * solver's, and pairing it with a different solver's golden compares two
 * algorithms, which no tolerance can fix.
 */
inline std::map<std::string, std::vector<View> > align_solvers(
    const std::map<std::string, std::vector<View> >& produced,
    const std::vector<std::string>& ref_names) {
    std::map<std::string, std::vector<View> > aligned;
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        std::map<std::string, std::vector<View> >::const_iterator it = produced.find(ref_names[i]);
        if (it != produced.end()) aligned[ref_names[i]] = it->second;
    }
    std::vector<std::string> left_ref;
    for (std::size_t i = 0; i < ref_names.size(); ++i)
        if (!aligned.count(ref_names[i])) left_ref.push_back(ref_names[i]);
    // IN RECORD ORDER, NOT ALPHABETICAL. `produced` is a std::map, so iterating
    // it offers `LN(MVA)` before `LN(NC)` -- and where an example solves the
    // same model with BOTH, the golden holds whichever ran FIRST. On
    // `lcq_threehosts`, which runs NC layers and then MVA layers, the map order
    // handed the single `LN` key to the MVA solve and scored it against the NC
    // golden: MVA and NC layers are different fixed points there (cache hit 0.5
    // against 0.48331), so all 22 of that row's cells failed on an alignment
    // rather than on a number. The Python comparator has always taken the first
    // RECORDED one, because its `produced` is a dict in record order; the two
    // must agree, and record order is the one that carries the meaning.
    std::vector<std::pair<long, std::string> > order;
    for (std::map<std::string, std::vector<View> >::const_iterator it = produced.begin();
         it != produced.end(); ++it) {
        long first = 0;
        for (std::size_t v = 0; v < it->second.size(); ++v)
            if (v == 0 || it->second[v].seq < first) first = it->second[v].seq;
        order.push_back(std::make_pair(first, it->first));
    }
    std::sort(order.begin(), order.end());
    for (std::size_t o = 0; o < order.size(); ++o) {
        const std::string& name = order[o].second;
        if (std::find(ref_names.begin(), ref_names.end(), name) != ref_names.end()) continue;
        std::string base, member;
        split_ensemble(name, &base, &member);
        const std::vector<std::string> cand = aliases(base, member);
        if (cand.empty()) continue;
        std::vector<std::string> hits;
        for (std::size_t k = 0; k < left_ref.size(); ++k)
            if (std::find(cand.begin(), cand.end(), left_ref[k]) != cand.end())
                hits.push_back(left_ref[k]);
        if (hits.size() != 1) continue;
        aligned[hits[0]] = produced.find(name)->second;
        left_ref.erase(std::find(left_ref.begin(), left_ref.end(), hits[0]));
    }
    return aligned;
}

/** How one key reads in a failure message. */
inline std::string describe(const Key& key) {
    return "('" + key.first + "', '" + key.second + "')";
}

/**
 * (failures, cells compared) for one row.
 *
 * CELLS is what makes a PASS mean something: a row is MEASURED only where a
 * golden solver and a produced table agree BY NAME, and the count is the
 * evidence that anything happened at all. Without it, an example whose every
 * golden solver was excused returns an empty failure list and reports PASS
 * having asserted nothing.
 *
 * `excused` names solvers the engine REPORTED as refused. Those are absent by
 * declaration rather than missing.
 */
inline std::vector<std::string> measure_row(const std::string& example, const Policy& policy,
                                            const std::map<std::string, std::vector<View> >& produced,
                                            const Golden& golden,
                                            const std::set<std::string>& excused, long* cells) {
    std::vector<std::string> failures;
    *cells = 0;
    const std::map<std::string, std::vector<View> > aligned = align_solvers(produced, golden.solvers);
    for (std::size_t s = 0; s < golden.solvers.size(); ++s) {
        const std::string& solver = golden.solvers[s];
        std::map<std::string, std::vector<View> >::const_iterator it = aligned.find(solver);
        if (it == aligned.end()) {
            if (!excused.count(solver))
                failures.push_back("solver " + solver + " missing from the recorded results");
            continue;
        }
        std::vector<Key> expected;
        for (Table::const_iterator k = golden.tables[s].begin(); k != golden.tables[s].end(); ++k)
            expected.push_back(k->first);
        Table actual;
        std::set<Key> derived;
        flatten(expected, it->second, &actual, &derived);
        double rel = -1.0, abs_tol = -1.0;
        policy.effective_tolerance(example, solver, golden.classes[s], &rel, &abs_tol);
        for (Table::const_iterator k = golden.tables[s].begin(); k != golden.tables[s].end(); ++k) {
            Table::const_iterator got = actual.find(k->first);
            if (got == actual.end()) {
                failures.push_back(solver + ": row " + describe(k->first) + " missing");
                continue;
            }
            for (std::map<std::string, double>::const_iterator m = k->second.begin();
                 m != k->second.end(); ++m) {
                std::map<std::string, double>::const_iterator have = got->second.find(m->first);
                if (have == got->second.end()) continue;
                ++(*cells);
                // THE TWO ALLOWANCES ADD. The declared tolerance says how far two
                // codebases may legitimately differ; the slack says how precisely
                // the golden RECORDED the answer. Taking the larger loses the
                // knife edge a generated per-example tolerance sits on.
                const double slack = golden_slack(m->second, derived.count(k->first) > 0);
                const double cell_abs = slack + (abs_tol < 0.0 ? 0.0 : abs_tol);
                double diff = 0.0;
                if (!policy.values_match(m->second, have->second, rel, cell_abs, &diff)) {
                    char buf[512];
                    std::snprintf(buf, sizeof(buf),
                                  "%s %s.%s: expected %.10g got %.10g (rel %.3g, tol abs=%.3g)",
                                  solver.c_str(), describe(k->first).c_str(), m->first.c_str(),
                                  m->second, have->second, diff, cell_abs);
                    failures.push_back(buf);
                }
            }
        }
    }
    return failures;
}

/**
 * Why this row checked nothing, or "" when it checked something.
 *
 * Every golden solver EXCUSED is a named skip. A golden solver NOT excused that
 * simply did not appear has already been reported as missing, so this can only
 * ever turn a false PASS into a named SKIP.
 */
inline std::string unmeasured_reason(const Golden& golden, const std::set<std::string>& excused) {
    if (golden.solvers.empty()) return std::string();
    std::string names;
    for (std::size_t i = 0; i < golden.solvers.size(); ++i) {
        if (!excused.count(golden.solvers[i])) return std::string();
        names += (names.empty() ? "" : ", ") + golden.solvers[i];
    }
    return "the row checked nothing: every solver in the golden (" + names +
           ") was excused, so no number was compared";
}

}  // namespace parity
}  // namespace tests
}  // namespace line

#endif  // LINE_TESTS_PARITY_COMPARE_H
