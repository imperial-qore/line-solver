#pragma once

/**
 * @file line_console.h
 * @brief Running progress log of a LINE solver run (the "solver console").
 *
 * The console narrates what a solver is doing while it does it: reading the
 * model, compiling the network structure, computing chains and demands,
 * resolving the method, iterating, and closing with the figures of merit. Each
 * line carries the elapsed time since the run started:
 *
 * @code
 * [   0.014s] compiling the network structure of model 'cqn'
 * [   0.031s]   computing the routing table and the chains
 * [   0.052s] recognized a closed queueing network: 3 stations, 1 class, 1 chain
 * [   0.061s]   AMVA sweep 10: queue-length residual 1.11e-01, X = 1.2225
 * @endcode
 *
 * It prints no tables: the result table stays the caller's own printer. THE
 * CONSOLE IS VerboseLevel::DEBUG: it narrates exactly when the run is at DEBUG
 * and is silent at every lower level, and it never alters a numerical result.
 * There is no separate console switch -- the console was one, until it became
 * clear that a running progress log IS what a debug verbosity is for, and two
 * switches for one channel only let a session ask for DEBUG and get nothing.
 * line-cli asks for it with `-v debug`; a library caller with set_verbose().
 *
 * Nested runs (an inner solver driven by an ensemble) do not narrate: only the
 * outermost run writes. Use is_active() to suppress a legacy print, owns_log()
 * to emit.
 *
 * Mirrors MATLAB matlab/src/io/LineConsole.m, jline.io.LineConsole and python
 * line_solver/api/io/console.py.
 */

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

namespace line {
namespace util {

/**
 * Session verbosity, the twin of MATLAB's VerboseLevel, jline.VerboseLevel and
 * python's line_solver.constants.VerboseLevel. The ordering matters and matches
 * theirs: SILENT < STD < DEBUG.
 *
 * This port had none until the console became DEBUG. A per-run `bool verbose`
 * still travels with each solver's options and still means STD-or-SILENT;
 * DEBUG is a SESSION level, because the console must also narrate the model
 * compile, which happens before any solver options exist.
 */
enum class VerboseLevel { SILENT = 0, STD = 1, DEBUG = 2 };

class LineConsole {
public:
    // ---------------------------------------------------------------- state

    /** Set the session verbosity. DEBUG is what switches the console on. */
    static void set_verbose(VerboseLevel level) {
        state().verbose = level;
        reset();
    }

    /** The session verbosity. */
    static VerboseLevel get_verbose() { return state().verbose; }

    /** Forget any open run (used after an interrupted solve). */
    static void reset() {
        State& st = state();
        st.depth = 0;
        st.active = false;
        st.tag.clear();
        st.model_name.clear();
        st.t0 = clock_type::now();
        st.tsetup = -1.0;
        st.quiet = 0;
        st.muted = 0;
        st.force_detail = false;
        st.compiling_own = false;
        st.header_pending = false;
        st.last_loop.clear();
        st.shown = 0;
        st.trunc = false;
        st.detail_last.clear();
        st.detail_shapes.clear();
        st.detail_shape_count.clear();
        st.detail_total = 0;
    }

    /**
     * Resolve whether a run should narrate.
     *
     * THE CONSOLE IS DEBUG, and this is the whole rule: a run narrates when the
     * session is at DEBUG and at no lower level. The run's own `verbose` flag
     * still vetoes -- it spells SILENT when false, and a run asked to stay
     * silent stays silent whatever the session asks for -- but it cannot switch
     * the console ON, because `true` there means STD, not DEBUG.
     */
    static bool wanted(bool verbose = true) {
        return state().verbose == VerboseLevel::DEBUG && verbose;
    }

    /** True while a run is narrating; gates SUPPRESSION of legacy prints. */
    static bool is_active() { return state().active; }

    /** True only inside the OUTERMOST open run; gates EMISSION. */
    static bool owns_log() {
        const State& st = state();
        return st.active && st.depth <= 1 && st.muted == 0;
    }

    /**
     * True when a progress line should be printed: either the outermost run is
     * narrating, or no run is open and the session is at DEBUG -- the second
     * case is what lets model construction narrate before any solver exists.
     */
    static bool writes() {
        const State& st = state();
        if (st.muted > 0) return false;
        if (st.depth > 0) return owns_log();
        return st.verbose == VerboseLevel::DEBUG;
    }

    // ------------------------------------------------------------ lifecycle

    /**
     * Open a console run. Pair with end_run(), or use the Run guard below so
     * that a failed analysis still reports what it had reached.
     *
     * @param tag        short solver name, e.g. "MVA"
     * @param model_name the model under study
     * @param verbose    the run's own verbosity flag
     */
    static bool begin_run(const std::string& tag, const std::string& model_name,
                          bool verbose = true) {
        State& st = state();
        if (st.depth == 0 && !wanted(verbose)) {
            // A run that must stay silent MUTES the console for its whole
            // duration, so nothing it triggers leaks out at depth 0.
            st.muted++;
            return false;
        }
        st.depth++;
        if (st.depth > 1) return false;  // nested: the outer analyzer owns the log
        st.active = true;
        st.tag = tag;
        st.model_name = model_name;
        st.t0 = clock_type::now();
        st.tsetup = -1.0;
        st.last_loop.clear();
        st.shown = 0;
        st.trunc = false;
        st.detail_last.clear();
        st.detail_shapes.clear();
        st.detail_shape_count.clear();
        st.detail_total = 0;
        // The CLI knows the solver method name before it reads the model file, so the
        // header waits for the model name rather than printing an empty one; it
        // is emitted by the first compile, or at close if none happens.
        if (model_name.empty()) {
            st.header_pending = true;
        } else {
            emit_header();
        }
        return true;
    }

    /** Mark the end of the setup phase, so the closing line can split the time. */
    static void setup_done() {
        State& st = state();
        if (st.active && st.depth == 1) st.tsetup = elapsed(st.t0);
    }

    /** Close the innermost open run, writing the closing line. */
    static void end_run() {
        State& st = state();
        if (st.depth <= 0) {
            if (st.muted > 0) st.muted--;
            return;
        }
        st.depth--;
        if (st.depth > 0 || !st.active) return;
        if (st.header_pending) {
            st.header_pending = false;
            emit_header(false); // a header deferred to the close keeps the clock
        }
        const double total = elapsed(st.t0);
        if (st.tsetup < 0.0) {
            step("DONE in %.4f s", total);
        } else {
            step("DONE in %.4f s (setup %.4f s, analysis %.4f s)", total, st.tsetup,
                 std::max(0.0, total - st.tsetup));
        }
        // the mute of an enclosing silenced run outlives this run's own state
        const int outer_mute = st.muted;
        reset();
        state().muted = outer_mute;
    }

    /** RAII guard: opens a run on construction, closes it on scope exit. */
    class Run {
    public:
        Run(const std::string& tag, const std::string& model_name, bool verbose = true) {
            begin_run(tag, model_name, verbose);
        }
        ~Run() { end_run(); }
        Run(const Run&) = delete;
        Run& operator=(const Run&) = delete;
    };

    /** RAII guard that suppresses structure-compile detail while it lives. */
    class Quiet {
    public:
        Quiet() { state().quiet++; }
        ~Quiet() { state().quiet = std::max(0, state().quiet - 1); }
        Quiet(const Quiet&) = delete;
        Quiet& operator=(const Quiet&) = delete;
    };

    // ------------------------------------------------------------- emission

    /** Write one progress line. */
    static void step(const char* fmt, ...) {
        if (!writes()) return;
        va_list args;
        va_start(args, fmt);
        const std::string text = vformat(fmt, args);
        va_end(args);
        emit("", text);
    }

    /** Write one indented progress line. */
    static void substep(const char* fmt, ...) {
        if (!writes()) return;
        va_list args;
        va_start(args, fmt);
        const std::string text = vformat(fmt, args);
        va_end(args);
        emit("  ", text);
    }

    /**
     * One stage line of a structure compile: silenced inside a Quiet scope,
     * and inside an open run, where the structures compiled are those of
     * auxiliary models rather than of the model under study.
     */
    static void compile_detail(const char* fmt, ...) {
        const State& st = state();
        if (st.quiet > 0) return;
        if (st.depth > 0 && !st.force_detail && !st.compiling_own) return;
        if (!writes()) return;
        va_list args;
        va_start(args, fmt);
        const std::string text = vformat(fmt, args);
        va_end(args);
        emit("  ", text);
    }

    /**
     * Announce the compilation of a model structure.
     *
     * A run opened before the model file was read carries no model name yet
     * (the CLI knows the solver method name first), so the FIRST compile inside such
     * a run names the run: it is the model under study, not an auxiliary one.
     */
    static void compiling(const std::string& name) {
        State& st = state();
        if (st.depth > 0 && st.model_name.empty()) {
            st.model_name = name;
            if (st.header_pending) {
                st.header_pending = false;
                emit_header();
            }
        }
        if (st.depth > 0 && name != st.model_name) {
            // an ensemble rebuilds the same submodel once per stage, so these
            // go through detail() and collapse to one line
            st.compiling_own = false;
            detail("refreshing the auxiliary submodel '" + name + "'");
        } else {
            st.compiling_own = true;
            if (st.depth == 0) { // a compile outside any run opens its own timeline
                st.clock = clock_type::now();
                st.clock_started = true;
            }
            step("compiling the network structure of model '%s'", name.c_str());
        }
    }

    /**
     * Report a solver's own debug message as a substep. Consecutive repeats are
     * dropped, at most three messages of the same SHAPE (the text with its
     * numbers masked) are reported, and the channel is capped per run.
     */
    static void detail(const std::string& raw) {
        if (!owns_log()) return;
        const std::string text = trim(raw);
        State& st = state();
        if (text.empty() || text == st.detail_last) return;
        const std::string shape = mask_numbers(text);
        for (std::size_t i = 0; i < st.detail_shapes.size(); ++i) {
            if (st.detail_shapes[i] == shape) {
                if (++st.detail_shape_count[i] > kMaxPerShape) return;
                st.detail_last = text;
                bump_and_emit(text);
                return;
            }
        }
        st.detail_shapes.push_back(shape);
        st.detail_shape_count.push_back(1);
        st.detail_last = text;
        bump_and_emit(text);
    }

    /**
     * Announce an iteration loop and reset its reporting budget. Re-announcing
     * the SAME text (a solver that restarts its loop) neither reprints the
     * header nor refills the budget.
     */
    static void loop(const char* fmt, ...) {
        if (!owns_log()) return;
        va_list args;
        va_start(args, fmt);
        const std::string text = vformat(fmt, args);
        va_end(args);
        State& st = state();
        if (text == st.last_loop) return;
        st.last_loop = text;
        st.shown = 0;
        st.trunc = false;
        emit("", text);
    }

    /**
     * Report iteration @p k of the current loop. The first 20 iterations report
     * in full, then every 10th, and the loop stops after 30 lines so a long run
     * cannot bury the rest of the narration.
     */
    static void iter(long k, const char* fmt, ...) {
        if (!owns_log()) return;
        if (k > 20 && (k % 10) != 0) return;
        State& st = state();
        if (st.shown >= kMaxIterLines) {
            if (!st.trunc) {
                st.trunc = true;
                emit("  ", "further iterations of this loop not reported");
            }
            return;
        }
        va_list args;
        va_start(args, fmt);
        const std::string text = vformat(fmt, args);
        va_end(args);
        st.shown++;
        emit("  ", text);
    }

    /** Count with an agreeing noun, e.g. "1 chain" or "3 chains". */
    static std::string plural(long n, const std::string& singular,
                              const std::string& plural_form) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%ld ", n);
        return std::string(buf) + (n == 1 ? singular : plural_form);
    }

    /** Lift the compile detail rule while the run compiles its OWN model. */
    class OwnModelCompile {
    public:
        OwnModelCompile() { state().force_detail = true; }
        ~OwnModelCompile() { state().force_detail = false; }
        OwnModelCompile(const OwnModelCompile&) = delete;
        OwnModelCompile& operator=(const OwnModelCompile&) = delete;
    };

private:
    typedef std::chrono::steady_clock clock_type;

    static const int kMaxIterLines = 30;
    static const int kMaxDetailLines = 200;
    static const int kMaxPerShape = 3;

    struct State {
        VerboseLevel verbose = VerboseLevel::STD;
        int depth = 0;
        bool active = false;
        std::string tag;
        std::string model_name;
        clock_type::time_point t0 = clock_type::now();
        double tsetup = -1.0;
        int quiet = 0;
        int muted = 0;
        bool force_detail = false;
        bool compiling_own = false;
        bool header_pending = false;
        std::string last_loop;
        int shown = 0;
        bool trunc = false;
        std::string detail_last;
        std::vector<std::string> detail_shapes;
        std::vector<int> detail_shape_count;
        int detail_total = 0;
        clock_type::time_point clock = clock_type::now(); // restarted per run
        bool clock_started = false;
    };

    static State& state() {
        static State st;
        return st;
    }

    static void emit_header(bool restart_clock = true) {
        State& st = state();
        if (restart_clock) { // each run's timeline starts at zero
            st.clock = clock_type::now();
            st.clock_started = true;
        }
        if (writes()) { // the opening row is set off from whatever preceded it
            std::printf("\n");
        }
        step("LINE: Solver%s starting on model '%s' (lang cpp)", st.tag.c_str(),
             st.model_name.empty() ? "(unnamed)" : st.model_name.c_str());
    }

    static double elapsed(const clock_type::time_point& from) {
        return std::chrono::duration<double>(clock_type::now() - from).count();
    }

    static void emit(const char* indent, const std::string& text) {
        State& st = state();
        if (!st.clock_started) {
            st.clock = clock_type::now();
            st.clock_started = true;
        }
        // a top-level row opens with a capital, an indented substep stays lowercase
        std::string row = text;
        if (indent[0] == '\0' && !row.empty()) {
            row[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(row[0])));
        }
        std::printf("[%8.3fs] %s%s\n", elapsed(st.clock), indent, row.c_str());
        std::fflush(stdout);
    }

    static void bump_and_emit(const std::string& text) {
        State& st = state();
        if (++st.detail_total > kMaxDetailLines) {
            if (st.detail_total == kMaxDetailLines + 1)
                emit("  ", "further solver detail not reported");
            return;
        }
        emit("  ", lower_first(text));
    }

    static std::string vformat(const char* fmt, va_list args) {
        va_list copy;
        va_copy(copy, args);
        const int n = std::vsnprintf(nullptr, 0, fmt, copy);
        va_end(copy);
        if (n <= 0) return std::string(fmt);
        std::vector<char> buf(static_cast<std::size_t>(n) + 1);
        std::vsnprintf(buf.data(), buf.size(), fmt, args);
        return std::string(buf.data(), static_cast<std::size_t>(n));
    }

    static std::string trim(const std::string& s) {
        std::size_t b = s.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) return std::string();
        std::size_t e = s.find_last_not_of(" \t\r\n");
        return s.substr(b, e - b + 1);
    }

    /** The text with every number replaced by '#', so repeats collapse. */
    static std::string mask_numbers(const std::string& s) {
        std::string out;
        out.reserve(s.size());
        bool in_number = false;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const char c = s[i];
            const bool digit = (c >= '0' && c <= '9');
            const bool part = digit || (in_number && (c == '.' || c == 'e' || c == 'E' ||
                                                      ((c == '+' || c == '-') && i > 0 &&
                                                       (s[i - 1] == 'e' || s[i - 1] == 'E'))));
            if (part) {
                if (!in_number) {
                    out.push_back('#');
                    in_number = true;
                }
            } else {
                in_number = false;
                out.push_back(c);
            }
        }
        return out;
    }

    static std::string lower_first(const std::string& s) {
        if (s.size() >= 2) {
            const bool acronym = (s[0] >= 'A' && s[0] <= 'Z') && (s[1] >= 'A' && s[1] <= 'Z');
            if (!acronym && s[0] >= 'A' && s[0] <= 'Z') {
                std::string out = s;
                out[0] = static_cast<char>(s[0] - 'A' + 'a');
                return out;
            }
        }
        return s;
    }
};

}  // namespace util
}  // namespace line
