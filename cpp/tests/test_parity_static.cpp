/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

#include <cstdio>
#include <cstdlib>
#include <sys/wait.h>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "doctest.h"
#include "parity_compare.h"

/**
 * @file test_parity_static.cpp
 * @brief Cross-codebase parity, asserted from inside the C++ suite.
 *
 * THE ROW THIS SUITE OWNS IS `CPP`: every goldened example run by its C++ twin,
 * with the result recorder on, checked cell by cell against the shared golden in
 * `goldens/baselines`. The PYTHON, P2J and P2C rows belong to the python suite,
 * MATLAB/M2J/M2P/M2C to the MATLAB suite and JAVA to the jar suite, because each
 * row is already asserted independently against the golden -- the standalone
 * harness batched them for SCHEDULING and never compared one row against
 * another, so splitting them across suites is a semantic identity.
 *
 * The twin runs as a SUBPROCESS (`line-examples --record`), which is what gives
 * isolation between examples: an example sets global state and several build
 * large models, so running 209 of them in one process would make each one's
 * result depend on its predecessors. That is also how the standalone harness has
 * always run them, which keeps the recorded values comparable to the goldens
 * those runs produced.
 *
 * SKIPPING IS BY NAME, NEVER BY DEFAULT. An unbuilt `line-examples` turns into
 * one named skip per example, and a twin that does not register the name into a
 * named skip for that example -- never into a silent pass. A refusal the twin
 * made BY NAME (`na("MAM", ...)`) excuses that solver and nothing else.
 *
 * Narrow a run with LINE_PARITY_EXAMPLES="cqn_oneline oqn_oneline".
 *
 * IT RUNS BY DEFAULT (since 2026-08-26). It was opt-in while it was new, on the
 * grounds that one subprocess per golden is a cost whoever runs phase 2 should
 * choose rather than discover. That reasoning inverted once this suite BECAME
 * the CPP row: the standalone `parity-static` harness is no longer the thing
 * that measures it, so a suite that skips by default leaves the row unmeasured
 * EVERYWHERE, which is the one outcome the named-skip discipline exists to
 * prevent. Set LINE_CPP_PARITY=0 to skip it -- for a quick unrelated run, not as
 * a way to make a red row go away.
 */
namespace {

namespace parity = line::tests::parity;
namespace fs = std::filesystem;

/** The shared golden directory: $LINE_GOLDENS, else $LINE_DEV, else in-tree. */
std::string golden_root() {
    if (const char* env = std::getenv("LINE_GOLDENS")) return env;
    if (const char* dev = std::getenv("LINE_DEV")) return std::string(dev) + "/goldens";
    return std::string(LINE_MP_REPO_ROOT) + "/goldens";
}

/** `common/line-examples`, or "" when this checkout has not built it. */
std::string twin_binary() {
    if (const char* env = std::getenv("LINE_EXAMPLES")) return env;
    const std::string path = std::string(LINE_MP_REPO_ROOT) + "/common/line-examples";
    std::error_code ec;
    return fs::exists(path, ec) ? path : std::string();
}

/** Every example that HAS a golden, sorted -- what parametrizes the suite. */
std::vector<std::string> goldened_examples() {
    std::vector<std::string> names;
    std::error_code ec;
    for (const fs::directory_entry& e : fs::directory_iterator(golden_root() + "/baselines", ec))
        if (e.path().extension() == ".json") names.push_back(e.path().stem().string());
    std::sort(names.begin(), names.end());
    return names;
}

/** The examples whose golden holds the FIRST of several models. */
std::set<std::string> multi_model_examples() {
    std::set<std::string> out;
    std::ifstream fh((golden_root() + "/corpus.json").c_str());
    if (!fh) return out;
    parity::json doc;
    fh >> doc;
    for (const parity::json& n : doc.at("multiModelExamples")) out.insert(n.get<std::string>());
    return out;
}

/** LINE_PARITY_EXAMPLES, split on whitespace; empty means every example. */
std::set<std::string> selected_examples() {
    std::set<std::string> out;
    const char* env = std::getenv("LINE_PARITY_EXAMPLES");
    if (!env) return out;
    std::string token;
    for (const char* p = env;; ++p) {
        if (*p == '\0' || *p == ' ' || *p == ',' || *p == '\t') {
            if (!token.empty()) out.insert(token);
            token.clear();
            if (*p == '\0') break;
        } else {
            token += *p;
        }
    }
    return out;
}

/** What one twin run produced. */
struct Run {
    bool ran = false;          ///< the binary executed and wrote a document
    bool registered = true;    ///< the twin knows this example's name
    parity::json doc;
    std::set<std::string> excused;
};

/**
 * The `FAMILY:method` keys this row cannot produce, excused BY NAME.
 *
 * A qualified key is a NON-DEFAULT method of a family, goldened by the METHOD
 * SWEEP (see `goldens/README.md`), and a row produces one by asking
 * `SolverAUTO(model, "mva.schmidt")` -- which the EXAMPLE never does. A twin in
 * `cpp/examples` runs its body and returns no model, so this suite has nothing
 * to replay the method against: excusing the key is the honest reading, where
 * "solver MVA:schmidt missing" would blame the row for a gap it reports
 * correctly.
 *
 * THIS IS A NAMED GAP, NOT A SILENT ONE, and it is printed. What closes it is
 * giving the twin a model handle, after which this row replays the keys exactly
 * as `python/tests/parity/rows.py` and `ParityRows.sweepModel` already do. Until
 * then the numbers ARE still checked, by the generator: its CPP row drives the
 * same pairs through `line-cli` and its gate refuses to write a key the C++
 * engine disagreed on.
 *
 * The bare family keys are unaffected and still measured here.
 */
std::set<std::string> unswept(const std::string& example,
                              const std::vector<std::string>& golden_keys) {
    std::set<std::string> keys;
    for (std::size_t i = 0; i < golden_keys.size(); ++i)
        if (golden_keys[i].find(':') != std::string::npos) keys.insert(golden_keys[i]);
    if (!keys.empty()) {
        std::string names;
        for (std::set<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it)
            names += (names.empty() ? "" : ", ") + *it;
        MESSAGE("UNSWEPT " << example << ": the C++ twin cannot replay a non-default "
                              "method, so " << names << " are excused here; the "
                              "generator's CPP row checks them.");
    }
    return keys;
}

/**
 * Excuse the solver a twin's note names, without ever excusing more than that.
 *
 * A note reads "<solver>: <why>". Taking the text before the FIRST colon was
 * right while every golden key was a bare family and is not any more: a swept
 * key is `MVA:schmidt`, so a note about it would excuse `MVA` -- the whole
 * family, and every bare-key cell the golden holds for it. That is the
 * false-pass shape the named-skip discipline exists to prevent, so the LONGEST
 * golden key the note names wins, and a note naming none excuses nothing.
 */
void add_excuse(const std::string& note, const std::vector<std::string>& keys,
                std::set<std::string>* excused) {
    std::string best;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        const std::string& key = keys[i];
        if (note.compare(0, key.size(), key) == 0 && note.size() > key.size() &&
            note[key.size()] == ':' && key.size() > best.size())
            best = key;
    }
    if (!best.empty()) {
        excused->insert(best);
        return;
    }
    const std::size_t colon = note.find(':');
    if (colon != std::string::npos) excused->insert(note.substr(0, colon));
}

/**
 * Run one twin with the recorder on.
 *
 * Exit status 2 is `line-examples`'s own "no such example", which is a NAMED
 * skip and not a failure -- the C++ tree does not carry a twin for every golden.
 * Any other non-zero status still yields a document, because the recorder writes
 * its dump even when a twin threw: a run that produced some tables and then died
 * is exactly the case that must be told apart from one that produced none.
 */
Run run_twin(const std::string& binary, const std::string& example,
             const std::vector<std::string>& golden_keys) {
    Run run;
    const fs::path out = fs::temp_directory_path() /
                         ("parity_cpp_" + example + "_" + std::to_string(::getpid()) + ".json");
    const std::string cmd = "'" + binary + "' --record '" + out.string() + "' '" + example +
                            "' > /dev/null 2>&1";
    const int status = std::system(cmd.c_str());
    const int code = (status == -1) ? -1 : WEXITSTATUS(status);
    if (code == 2) {
        run.registered = false;
        std::error_code ec;
        fs::remove(out, ec);
        return run;
    }
    std::ifstream fh(out.string().c_str());
    if (fh) {
        fh >> run.doc;
        run.ran = true;
    }
    std::error_code ec;
    fs::remove(out, ec);
    if (run.ran && run.doc.contains("notes"))
        for (const parity::json& n : run.doc.at("notes"))
            add_excuse(n.get<std::string>(), golden_keys, &run.excused);
    return run;
}

}  // namespace

TEST_CASE("parity: the CPP row against the shared goldens" * doctest::skip(false)) {
    // OPT-OUT, NOT OPT-IN: this suite IS the CPP row, so silence here is the row
    // going unmeasured rather than a cost avoided.
    const char* opt_out = std::getenv("LINE_CPP_PARITY");
    if (opt_out && std::string(opt_out) == "0") {
        MESSAGE("SKIP: LINE_CPP_PARITY=0 turned the CPP parity row off. Nothing "
                "about the CPP row is asserted by this run.");
        return;
    }
    const std::string root = golden_root();
    REQUIRE_MESSAGE(fs::exists(root), "shared goldens not found at " << root
                                      << "; set LINE_GOLDENS or LINE_DEV");
    const std::string binary = twin_binary();
    const std::vector<std::string> examples = goldened_examples();
    REQUIRE_MESSAGE(!examples.empty(), "no goldens under " << root << "/baselines");

    if (binary.empty()) {
        // ONE NAMED SKIP, NOT 209 SILENT PASSES. The count is printed so a reader
        // sees exactly how much went unmeasured.
        MESSAGE("common/line-examples is not built, so all "
                << examples.size()
                << " CPP parity rows are unmeasured; build it with `cd cpp && ./make.sh -O`");
        return;
    }

    const parity::Policy policy(root);
    const std::set<std::string> multi = multi_model_examples();
    const std::set<std::string> only = selected_examples();
    long checked = 0, skipped = 0;

    for (std::size_t i = 0; i < examples.size(); ++i) {
        const std::string& example = examples[i];
        if (!only.empty() && !only.count(example)) continue;
        CAPTURE(example);
        const parity::Golden golden = parity::load_golden(root, example);
        const Run run = run_twin(binary, example, golden.solvers);
        if (!run.registered) {
            MESSAGE("SKIP " << example << ": no C++ twin registers that name");
            ++skipped;
            continue;
        }
        if (!run.ran) {
            FAIL_CHECK(example << ": the twin wrote no result document");
            continue;
        }
        const std::map<std::string, std::vector<parity::View> > produced =
            parity::records_to_views(run.doc, multi.count(example) > 0);
        std::set<std::string> excused = run.excused;
        const std::set<std::string> skipped_keys = unswept(example, golden.solvers);
        excused.insert(skipped_keys.begin(), skipped_keys.end());
        long cells = 0;
        const std::vector<std::string> failures =
            parity::measure_row(example, policy, produced, golden, excused, &cells);

        // A ROW THAT COMPARED NO CELL IS NOT A PASS. Reaching here with no
        // failures and no cells means every golden solver was excused, so the row
        // asserted nothing -- and reporting a pass there is how an unmeasured
        // defect stays green.
        if (failures.empty() && cells == 0) {
            const std::string reason = parity::unmeasured_reason(golden, excused);
            if (!reason.empty()) {
                MESSAGE("SKIP " << example << ": " << reason);
                ++skipped;
                continue;
            }
            FAIL_CHECK(example << ": the row compared no cell against the golden, "
                                  "so it asserted nothing");
            continue;
        }
        for (std::size_t f = 0; f < failures.size(); ++f) FAIL_CHECK(example << ": " << failures[f]);
        if (failures.empty()) ++checked;
    }
    MESSAGE("CPP parity: " << checked << " examples agreed, " << skipped << " skipped by name");
}
