/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <streambuf>
#include <string>
#include <thread>

#if defined(_WIN32)
#include <io.h>
#define LINE_MP_ISATTY() (_isatty(2) != 0)
#else
#include <unistd.h>
#define LINE_MP_ISATTY() (isatty(STDERR_FILENO) != 0)
#endif

// The suite is ~2800 test cases and runs for ~25 minutes at the -O0 that is the
// rule in this tree, so a run with no output until the summary is
// indistinguishable from a hung one. The progress meter below is a doctest
// LISTENER, not a reporter: it runs ALONGSIDE the console reporter rather than
// replacing it, so the failure blocks and the "[doctest] test cases:" tail that
// run-tests.sh quotes are byte-for-byte what they were.
//
// It writes to STDERR, never to the reporter's stream, and it flushes that
// stream first: doctest's console reporter does not flush per test case, so a
// piped run would otherwise land the progress lines out of order in the log.
// On a terminal it rewrites one line in place; when stderr is not a tty (the
// suite driver redirects the phase into a log) it emits one line per 1% step
// instead, which is 100 lines rather than a smeared carriage-return soup.
//
// A PERCENTAGE STEP ALONE CANNOT BOUND THE GAP between two lines, however
// small the step: the cases are wildly uneven, and a single one that runs for
// a quarter of an hour produces no listener callback at all while it runs, so
// the log goes silent exactly where the reader most needs to know what is
// holding it. A HEARTBEAT THREAD closes that: it wakes on a timer and re-emits
// the meter whenever nothing has been written for LINE_CPP_TEST_PROGRESS_INTERVAL
// seconds (default 120), naming the case in flight and how long it has been in
// it, so the slow case identifies itself instead of being bracketed. The
// heartbeat does NOT flush the reporter's stream -- that stream is written from
// the test thread and flushing it from here would be a data race -- so a
// heartbeat line can sort slightly early against reporter output; the step
// lines, which come from the test thread, keep the flush and the ordering.
//
// The denominator is the REGISTERED test count, and skipped cases advance the
// counter too, so a filtered run (-tc=..., --first/--last) still reaches 100%.
namespace {

std::string mp_format_duration(double seconds) {
    if (!(seconds >= 0)) return "--";
    const long total = static_cast<long>(seconds + 0.5);
    char buf[32];
    if (total >= 3600)
        std::snprintf(buf, sizeof(buf), "%ldh%02ldm", total / 3600, (total % 3600) / 60);
    else if (total >= 60)
        std::snprintf(buf, sizeof(buf), "%ldm%02lds", total / 60, total % 60);
    else
        std::snprintf(buf, sizeof(buf), "%lds", total);
    return std::string(buf);
}

bool mp_progress_wanted(const doctest::ContextOptions& opts) {
    if (opts.quiet) return false;
    const char* env = std::getenv("LINE_CPP_TEST_PROGRESS");
    if (env == nullptr) return true;
    return !(std::strcmp(env, "0") == 0 || std::strcmp(env, "off") == 0 ||
             std::strcmp(env, "no") == 0);
}

// A non-negative integer from the environment, or the default when the variable
// is unset, empty or not a number. 0 is a legal value and means "off" for both
// callers, so it is never silently promoted to the default.
long mp_env_long(const char* name, long fallback) {
    const char* env = std::getenv(name);
    if (env == nullptr || *env == '\0') return fallback;
    char* end = nullptr;
    const long value = std::strtol(env, &end, 10);
    if (end == env || *end != '\0' || value < 0) return fallback;
    return value;
}

bool mp_stderr_tag_wanted() {
    const char* env = std::getenv("LINE_CPP_TEST_STDERR_TAG");
    if (env == nullptr) return true;
    return !(std::strcmp(env, "0") == 0 || std::strcmp(env, "off") == 0 ||
             std::strcmp(env, "no") == 0);
}

// An engine that writes a diagnostic to std::cerr -- lsoda's "10000000 steps
// taken before reaching tout" is the recurring one -- lands it in the suite log
// with NO test case attached: doctest's console reporter names a case only when
// it fails, so a warning raised by a PASSING case is anonymous and the reader is
// left bracketing it between two progress markers. This streambuf closes that
// gap by prefixing every std::cerr LINE with the case that was running when the
// line completed, e.g. "[LN methods: the transient ...] [lsoda] ...".
//
// It buffers to the newline and NEVER emits on sync(). std::cerr is unitbuf, so
// each `<<` flushes; emitting on sync would split a line assembled from several
// insertions across several prefixes. A partial line left over at the end of a
// case is flushed then, so it cannot be attributed to the next one.
//
// The progress meter is unaffected: it writes through fprintf(stderr) on the C
// stream, below this buffer, so its lines stay untagged and unwrapped.
class TaggingStderrBuf : public std::streambuf {
public:
    explicit TaggingStderrBuf(std::streambuf* inner) : m_inner(inner) {}

    ~TaggingStderrBuf() override {
        std::lock_guard<std::mutex> guard(m_mutex);
        flush_partial();
    }

    void set_tag(const char* name) {
        std::lock_guard<std::mutex> guard(m_mutex);
        flush_partial();
        m_tag = name != nullptr ? name : "";
    }

protected:
    int_type overflow(int_type ch) override {
        if (traits_type::eq_int_type(ch, traits_type::eof())) return traits_type::not_eof(ch);
        std::lock_guard<std::mutex> guard(m_mutex);
        put(traits_type::to_char_type(ch));
        return ch;
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        std::lock_guard<std::mutex> guard(m_mutex);
        for (std::streamsize i = 0; i < n; ++i) put(s[i]);
        return n;
    }

    int sync() override { return m_inner->pubsync(); }

private:
    void put(char c) {
        m_line.push_back(c);
        if (c == '\n') emit();
    }

    // Caller holds m_mutex.
    void emit() {
        std::string out;
        if (!m_tag.empty()) out = "[" + m_tag + "] ";
        out += m_line;
        m_inner->sputn(out.data(), static_cast<std::streamsize>(out.size()));
        m_inner->pubsync();
        m_line.clear();
    }

    // Caller holds m_mutex.
    void flush_partial() {
        if (m_line.empty()) return;
        m_line.push_back('\n');
        emit();
    }

    std::streambuf* m_inner;
    std::mutex m_mutex;
    std::string m_line;
    std::string m_tag;
};

class ProgressListener : public doctest::IReporter {
public:
    explicit ProgressListener(const doctest::ContextOptions& opts)
        : m_opts(opts), m_enabled(mp_progress_wanted(opts)), m_tty(LINE_MP_ISATTY()) {}

    void report_query(const doctest::QueryData&) override {}

    void test_run_start() override {
        install_tagger();
        if (!m_enabled) return;
        {
            std::lock_guard<std::mutex> guard(m_mu);
            m_total = doctest::detail::getRegisteredTests().size();
            m_start = clock_type::now();
            m_case_start = m_start;
            m_last_emit = m_start;
            m_done = 0;
            m_next_step = 0;
            m_last_pct = -1;
            m_current.clear();
        }
        start_heartbeat();
    }

    void test_run_end(const doctest::TestRunStats&) override {
        stop_heartbeat();
        if (m_enabled) {
            std::lock_guard<std::mutex> guard(m_mu);
            m_current.clear();
            if (m_tty) {
                std::fputs("\r\033[K", stderr);
                std::fflush(stderr);
            } else if (m_last_pct != 100) {
                emit(true, false);
            }
        }
        uninstall_tagger();
    }

    void test_case_start(const doctest::TestCaseData& in) override {
        {
            std::lock_guard<std::mutex> guard(m_mu);
            m_current = in.m_name != nullptr ? in.m_name : "";
            m_case_start = clock_type::now();
        }
        set_tag(in.m_name);
        render();
    }

    void test_case_reenter(const doctest::TestCaseData&) override {}

    void test_case_end(const doctest::CurrentTestCaseStats&) override {
        finish_case();
    }

    void test_case_skipped(const doctest::TestCaseData&) override {
        finish_case();
    }

    void test_case_exception(const doctest::TestCaseException&) override {}
    void subcase_start(const doctest::SubcaseSignature&) override {}
    void subcase_end() override {}
    void log_assert(const doctest::AssertData&) override {}
    void log_message(const doctest::MessageData&) override {}

private:
    typedef std::chrono::steady_clock clock_type;

    // Not installed when the reporter itself writes to std::cerr (--out), since
    // tagging the reporter's own blocks would corrupt the output run-tests.sh
    // greps, nor when LINE_CPP_TEST_STDERR_TAG says no.
    void install_tagger() {
        if (m_tagger != nullptr) return;
        if (!mp_stderr_tag_wanted()) return;
        if (m_opts.cout == &std::cerr) return;
        m_saved_cerr = std::cerr.rdbuf();
        m_tagger.reset(new TaggingStderrBuf(m_saved_cerr));
        std::cerr.rdbuf(m_tagger.get());
    }

    void uninstall_tagger() {
        if (m_tagger == nullptr) return;
        m_tagger->set_tag(nullptr); // flushes any partial line under its own tag
        std::cerr.rdbuf(m_saved_cerr);
        m_tagger.reset();
        m_saved_cerr = nullptr;
    }

    void set_tag(const char* name) {
        if (m_tagger != nullptr) m_tagger->set_tag(name);
    }

    void finish_case() {
        {
            std::lock_guard<std::mutex> guard(m_mu);
            ++m_done;
            m_current.clear();
        }
        set_tag(nullptr);
        render();
    }

    // Called from the TEST thread, once per case boundary.
    void render() {
        if (!m_enabled) return;
        std::lock_guard<std::mutex> guard(m_mu);
        if (m_total == 0) return;
        const int pct = static_cast<int>((100.0 * m_done) / m_total);
        if (m_tty) {
            emit(false, false);
        } else if ((m_step > 0 && pct >= m_next_step) || due(clock_type::now())) {
            if (m_step > 0) m_next_step = ((pct / m_step) + 1) * m_step;
            emit(false, false);
        }
    }

    // Caller holds m_mu. True once the meter has been silent for the heartbeat
    // interval; always false when the heartbeat is switched off.
    bool due(clock_type::time_point now) const {
        if (m_interval <= 0) return false;
        return std::chrono::duration<double>(now - m_last_emit).count() >= m_interval;
    }

    void start_heartbeat() {
        if (!m_enabled || m_interval <= 0) return;
        m_stop = false;
        m_heart = std::thread(&ProgressListener::heartbeat_loop, this);
    }

    void stop_heartbeat() {
        if (!m_heart.joinable()) return;
        {
            std::lock_guard<std::mutex> guard(m_mu);
            m_stop = true;
        }
        m_cv.notify_all();
        m_heart.join();
    }

    // Wakes at the instant the current line goes stale rather than polling, so
    // an idle suite costs one wakeup per interval and nothing else.
    void heartbeat_loop() {
        std::unique_lock<std::mutex> lock(m_mu);
        while (!m_stop) {
            const clock_type::duration gap =
                std::chrono::duration_cast<clock_type::duration>(
                    std::chrono::duration<double>(m_interval));
            if (m_cv.wait_until(lock, m_last_emit + gap, [this] { return m_stop; })) return;
            if (m_total != 0 && due(clock_type::now())) emit(false, true);
        }
    }

    // Caller holds m_mu.
    void emit(bool final_line, bool from_heartbeat) {
        if (m_total == 0) return;
        // The reporter's stream belongs to the test thread; flushing it from the
        // heartbeat thread would race with the reporter writing to it.
        if (m_opts.cout != nullptr && !from_heartbeat) m_opts.cout->flush();

        const clock_type::time_point now = clock_type::now();
        const double elapsed = std::chrono::duration<double>(now - m_start).count();
        const int pct = static_cast<int>((100.0 * m_done) / m_total);
        m_last_pct = pct;
        m_last_emit = now;

        std::string tail;
        if (!m_opts.no_time_in_output) {
            tail = " " + mp_format_duration(elapsed);
            if (m_done > 0 && m_done < m_total)
                tail += " eta " + mp_format_duration(elapsed * (m_total - m_done) / m_done);
        }

        char head[96];
        std::snprintf(head, sizeof(head), "[cpp tests] %3d%% (%zu/%zu)", pct, m_done, m_total);

        // The case in flight, and how long it has been in flight: on a heartbeat
        // that is the whole point of the line, since the reader is waiting on
        // exactly this case and the log names it nowhere else until it fails.
        std::string running;
        if (!m_current.empty()) {
            running = m_current;
            if (running.size() > kNameWidth) running = running.substr(0, kNameWidth - 3) + "...";
        }

        if (m_tty && !final_line) {
            std::fprintf(stderr, "\r\033[K%s%s  %s", head, tail.c_str(), running.c_str());
        } else if (from_heartbeat && !running.empty()) {
            const double in_case = std::chrono::duration<double>(now - m_case_start).count();
            std::fprintf(stderr, "%s%s  in %s %s\n", head, tail.c_str(), running.c_str(),
                         mp_format_duration(in_case).c_str());
        } else {
            std::fprintf(stderr, "%s%s\n", head, tail.c_str());
        }
        std::fflush(stderr);
    }

    static const size_t kNameWidth = 48;

    const doctest::ContextOptions& m_opts;
    bool m_enabled;
    bool m_tty;
    // Percentage step for the redirected path, and the ceiling on the silence
    // between two lines, in seconds. 0 switches the respective rule off; with
    // the step off the heartbeat alone paces the log.
    int m_step = static_cast<int>(mp_env_long("LINE_CPP_TEST_PROGRESS_STEP", 1));
    double m_interval =
        static_cast<double>(mp_env_long("LINE_CPP_TEST_PROGRESS_INTERVAL", 120));
    mutable std::mutex m_mu; // guards every counter below and the emitting itself
    std::condition_variable m_cv;
    std::thread m_heart;
    bool m_stop = false;
    size_t m_total = 0;
    size_t m_done = 0;
    int m_next_step = 0;
    int m_last_pct = -1;
    std::string m_current;
    clock_type::time_point m_start;
    clock_type::time_point m_case_start;
    clock_type::time_point m_last_emit;
    std::unique_ptr<TaggingStderrBuf> m_tagger;
    std::streambuf* m_saved_cerr = nullptr;
};

} // namespace

DOCTEST_REGISTER_LISTENER("line_progress", 0, ProgressListener);
