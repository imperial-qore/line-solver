/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_TESTS_EXAMPLE_SUITE_H
#define LINE_TESTS_EXAMPLE_SUITE_H

/**
 * @file example_suite.h
 * @brief Shared machinery for the shipped-example regression suites.
 *
 * WHY THIS LIVES IN cpp/ AND NOT IN line-test.git. `line-dev.git` keeps the
 * examples and the tests that exercise them (line-test.git/README.md). These
 * two translation units are the C++ half of that rule, and they are the reason
 * cpp/tests/ is no longer empty: a checkout with no sibling still gets a
 * regression suite over the examples it ships.
 *
 * WHAT IT ASSERTS: every registered example runs to completion without
 * throwing, and prints something. That is the whole invariant, and it is
 * deliberately narrow.
 *
 * WHAT IT USED TO ASSERT, AND WHY THAT WAS WRONG. A first version armed the
 * parity recorder and required every recorded value to be FINITE, and every
 * example to record a table or a refusal. Both are false invariants, and the
 * 2026-09-12 run proved it:
 *   - NaN IS AN EXPECTED VALUE. It is LINE's "not applicable for this station
 *     and metric" sentinel, and the goldens carry it as such --
 *     `goldens/baselines/fcr_mm1waitq.json` records `"Util": "NaN"` for station
 *     FCR1, and lcq_singlehost/lcq_threehosts carry 22 and 32 NaN cells. A
 *     finiteness rule fails on CORRECT output.
 *   - NOT EVERY EXAMPLE RECORDS. The recorder hooks `section`, `na`, `avg_rows`
 *     and two printers (parity_recorder.h); an example that prints its own
 *     table with std::printf, as busyp_subnetwork and the networkCalculus pair
 *     do, computes and prints real results while recording nothing.
 * So the recorder is not a general-purpose oracle: it exists to key results by
 * SOLVER for cross-codebase goldens, and it covers exactly the corpus that has
 * goldens. Comparing values against those goldens is parity-static's job, it
 * has the goldens to do it properly, and it stays there. Re-deriving a weaker
 * golden-less version here bought nothing and rejected correct examples.
 *
 * This suite answers the narrower question that has to be answerable WITHOUT
 * the sibling checkout: does every shipped example still run.
 */

#include <fcntl.h>
#include <unistd.h>

#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "doctest.h"

#include "examples_common.h"

namespace line_tests {

/**
 * Capture an example's own output for the duration of one run, and measure it.
 *
 * An example's whole body is a printed table; 212 of them would bury the
 * doctest report under tens of thousands of lines, and that report is the only
 * thing a CI reader actually looks at. The bytes are counted rather than
 * discarded because "it printed nothing at all" is the one cheap, UNIVERSAL
 * sign that an example did no work -- true for the examples that feed the
 * parity recorder and equally for the ones that print their own tables.
 *
 * THE REDIRECT IS AT THE FILE DESCRIPTOR, not at std::cout. The examples print
 * through std::printf, and C stdio and iostreams carry separate buffers, so
 * swapping cout's streambuf would leave every printf untouched. Both are
 * flushed on the way in and on the way out so no buffered line crosses the swap
 * and lands in the wrong place. Failure to redirect is not fatal: `captured()`
 * then reports false and the caller asserts nothing about the size, because a
 * suite that cannot run would be worse than a verbose one.
 */
class StdoutCapture {
 public:
    StdoutCapture() : m_saved(-1), m_sink(NULL), m_bytes(0), m_ok(false) {
        std::fflush(stdout);
        std::cout.flush();
        m_sink = std::tmpfile();
        if (m_sink == NULL) return;
        m_saved = ::dup(STDOUT_FILENO);
        if (m_saved < 0) {
            std::fclose(m_sink);
            m_sink = NULL;
            return;
        }
        m_ok = (::dup2(::fileno(m_sink), STDOUT_FILENO) >= 0);
    }

    ~StdoutCapture() { restore(); }

    /** Put stdout back, recording how much was written. Idempotent. */
    void restore() {
        if (m_saved < 0) return;
        std::fflush(stdout);
        std::cout.flush();
        // The offset AFTER the flush is exactly what the example wrote.
        const long pos = ::lseek(::fileno(m_sink), 0, SEEK_CUR);
        m_bytes = (pos > 0) ? static_cast<std::size_t>(pos) : 0;
        ::dup2(m_saved, STDOUT_FILENO);
        ::close(m_saved);
        m_saved = -1;
        // tmpfile() is already unlinked, so closing it is the whole cleanup.
        std::fclose(m_sink);
        m_sink = NULL;
    }

    /**
     * Whether the redirect actually took. This must NOT be inferred from the
     * byte count: "the capture failed" and "the example printed nothing" are
     * the two cases the size assertion has to tell apart, and both leave
     * bytes() == 0.
     */
    bool captured() const { return m_ok; }
    std::size_t bytes() const { return m_bytes; }

 private:
    StdoutCapture(const StdoutCapture&);
    StdoutCapture& operator=(const StdoutCapture&);

    int m_saved;
    std::FILE* m_sink;
    std::size_t m_bytes;
    bool m_ok;
};

/** The registered examples whose group starts with `prefix` ("basic/"). */
inline std::vector<line::examples::Example> examples_under(const std::string& prefix) {
    std::vector<line::examples::Example> out;
    const std::vector<line::examples::Example>& all = line::examples::registry();
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (all[i].group.size() >= prefix.size() &&
            all[i].group.compare(0, prefix.size(), prefix) == 0) {
            out.push_back(all[i]);
        }
    }
    return out;
}

/** Outcome of one example, collected with stdout captured and reported after. */
struct RunOutcome {
    bool threw;
    std::string error;
    bool measured;          ///< the redirect took, so out_bytes is meaningful
    std::size_t out_bytes;  ///< what the example printed
};

/** Run one example with its output captured, reporting nothing itself. */
inline RunOutcome run_one(const line::examples::Example& ex) {
    RunOutcome out;
    out.threw = false;
    out.measured = false;
    out.out_bytes = 0;

    {
        StdoutCapture capture;
        // Running an example IN PROCESS is safe here in a way it is not on the
        // Java side, where ExampleRunnerSupport must maintain EXCLUDED_MAINS
        // because a main() that calls System.exit takes the whole fork with it.
        // No C++ example calls exit, abort or terminate, and none writes a file
        // into the working directory, so there is nothing to exclude.
        try {
            ex.run();
        } catch (const std::exception& e) {
            out.threw = true;
            out.error = e.what();
        } catch (...) {
            out.threw = true;
            out.error = "unknown exception";
        }
        capture.restore();
        out.measured = capture.captured();
        out.out_bytes = capture.bytes();
    }  // stdout restored here, BEFORE any assertion below can report

    return out;
}

/**
 * Run every example under `prefix` and require each to complete.
 *
 * `expected_min` is a LINK GUARD, and it is why this asserts a count at all.
 * The examples reach the suite through an OBJECT library because each registers
 * itself from a static initializer that nothing references; were that ever to
 * become a static archive, the linker would drop the objects, the registry
 * would come back EMPTY, and a loop over nothing would report a green suite
 * having exercised no example. That silent pass is the one failure this suite
 * must not have. It is a floor and not an equality so that ADDING an example
 * never fails the suite that runs it.
 */
inline void run_examples_under(const std::string& prefix, std::size_t expected_min) {
    const std::vector<line::examples::Example> group = examples_under(prefix);
    REQUIRE_MESSAGE(group.size() >= expected_min,
                    "the registry holds " << group.size() << " examples under '" << prefix
                                          << "', expected at least " << expected_min
                                          << ": the static registrations were dropped at "
                                          << "link time, so this suite would pass over "
                                          << "nothing");

    for (std::size_t i = 0; i < group.size(); ++i) {
        const line::examples::Example& ex = group[i];
        const RunOutcome out = run_one(ex);
        const std::string id = ex.group + "/" + ex.name;

        CHECK_MESSAGE(!out.threw, id << " threw: " << out.error);
        if (out.threw) continue;

        // Only when the redirect took: an unmeasured run is silent about size
        // rather than failing on it. Every example ends in a printed table, so
        // zero bytes means it returned without doing its work -- the one
        // "produced nothing" signal that holds for the whole corpus, including
        // the examples that never touch the parity recorder.
        if (out.measured) {
            CHECK_MESSAGE(out.out_bytes > 0,
                          id << " printed nothing: it produced no result");
        }
    }
}

}  // namespace line_tests

#endif  // LINE_TESTS_EXAMPLE_SUITE_H
