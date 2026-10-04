/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_CLI_CLI_RUN_H
#define LINE_CLI_CLI_RUN_H

/**
 * @file
 * @ingroup line_public
 * The CLI's argument vector, reachable in-process.
 *
 * `line-cli` is the only place in this port that serves the WHOLE vocabulary:
 * every solver, every `-a` analysis, every input format and the `--api`
 * registry. The in-process facade (`line/solvers/solver.h`) is deliberately a
 * facade -- `double` only, `avg_table()` and a handful of siblings -- so a host
 * that wants a response-time CDF, a layered solve or a state probability has
 * had no route but to spawn the binary. This header is that route.
 *
 * THE VOCABULARY IS argv, AND THAT IS THE POINT. `parse_args` accepts about a
 * hundred and twenty flags and gains more with every solver; a struct mirroring
 * them would be a second parser to keep in step, and every host binding would
 * have to learn it. Passing the argument vector verbatim means a host that can
 * build a command line can reach everything the command line reaches, on the
 * day the flag lands rather than on the day someone widens a struct.
 *
 * WHAT THIS ADDS OVER `system("line-cli ...")` is the three things a subprocess
 * cannot give: the JSON documents arrive as separate strings rather than as
 * text to be scanned for a brace, stdout and stderr are captured instead of
 * escaping into the host's console, and an error arrives with the stable
 * identifier of `line/io/marshal.h` rather than as prose on fd 2.
 *
 * REENTRANCY IS THE CALLER'S TO SERIALIZE. `line_cli.cpp` holds file-scope
 * state (the output format, the two input-format flags, the verbosity and the
 * buffered stdin text) and the capture redirects process-global file
 * descriptors, so two concurrent `run()` calls would interfere. `run()` resets
 * every one of those before it dispatches, so SEQUENTIAL calls are independent;
 * concurrent ones need a lock the caller holds. The C ABI in `r/capi/` takes
 * one, and that is the intended arrangement rather than a workaround.
 */

#include <string>
#include <vector>

namespace line {
namespace cli {

/** One invocation: the argument vector, plus what a pipe would have carried. */
struct Request {
    /**
     * The argument vector WITHOUT argv[0], e.g. {"-f", "m.json", "-s", "mva"}.
     * A leading program name is not expected and would be parsed as a stray
     * positional argument.
     */
    std::vector<std::string> argv;

    /**
     * The model document, as if piped to the binary's stdin.
     *
     * `has_stdin` rather than an empty check, because an empty document is a
     * legitimate thing to hand the reader and must produce the reader's own
     * refusal rather than a read from the host's real stdin.
     */
    std::string stdin_text;
    bool has_stdin = false;

    /** Collect the JSON documents as they are emitted. */
    bool collect_json = true;
    /** Redirect fd 1 for the duration of the call into Response::text. */
    bool capture_stdout = true;
    /** Redirect fd 2 for the duration of the call into Response::diagnostics. */
    bool capture_stderr = true;
    /**
     * Permit `-p`, which serves until its request budget is exhausted.
     *
     * OFF BY DEFAULT because a host that passes a user's argument string
     * through would otherwise hand that user a way to block the calling thread
     * indefinitely and open a listening socket. It is a capability, not a
     * safety check: a caller who means it says so.
     */
    bool allow_server = false;

    /**
     * Polled between analyses and at each document emission; non-zero aborts.
     *
     * A HOST CANNOT SIMPLY LONGJMP OUT OF HERE. R's `R_CheckUserInterrupt`
     * does exactly that, and doing it from a frame with C++ destructors above
     * it skips every one of them, leaking the capture's file descriptors and
     * its temporary files. So the host polls its own flag inside this callback
     * and says so by return value, and the abort unwinds normally.
     *
     * WHERE IT IS POLLED IS DELIBERATELY NARROW, and saying so is the point: a
     * half-working interrupt that is advertised as working is worse than none.
     * It fires between the analyses of a comma-separated `-a` list and at each
     * JSON emission. It does NOT reach inside one long solve, because the
     * solver loops take no callback and threading one through them is an
     * engine-wide change this seam has no business making.
     */
    int (*interrupt)(void* user) = nullptr;
    void* interrupt_user = nullptr;
};

/** What one invocation produced. */
struct Response {
    /** The CLI's own exit status: 0 ok, 2 `line::Error`, 3 other. */
    int exit_code = 0;

    /**
     * Each JSON document the run emitted, in emission order, AS EMITTED.
     *
     * STRINGS AND NOT PARSED VALUES, deliberately. The eleven emission sites
     * use three different indents and the numbers are printed at full
     * precision so a cross-codebase diff is not capped
     * (`line_cli.cpp`, `emit_avg_table_named`). Parsing here and re-dumping in
     * the host would reformat both, and the reformatting would be invisible
     * until someone diffed a golden. A host that wants a value parses this
     * string itself, once, at the boundary where it builds its own objects.
     *
     * PARTIAL ON FAILURE BY DESIGN: `-a avg,sens` whose `sens` arm throws has
     * already produced the `avg` envelope, and the binary printed it. Dropping
     * it here would lose an answer the equivalent command line gives. Read
     * `exit_code` for whether the run as a whole succeeded.
     */
    std::vector<std::string> documents;

    /** Everything written to fd 1, with the collected documents still in it. */
    std::string text;
    /** Everything written to fd 2: warnings, and the failure message. */
    std::string diagnostics;

    /** Empty on success; the `what()` of the exception otherwise. */
    std::string error;
    /**
     * `line::io::error_id`'s stable identifier for `error`.
     *
     * One of `line:input`, `line:numeric`, `line:unsupported`, `line:error`,
     * or `line:internal` for a `std::exception` that is not a `line::Error`
     * and which `error_id` therefore cannot classify.
     */
    std::string error_id;
};

/**
 * Run one invocation and return what it produced.
 *
 * Throws nothing: every failure the binary would report on fd 2 and an exit
 * status arrives here as `error`, `error_id` and `exit_code`.
 */
Response run(const Request& req);

/** The binary's `main`, so `line_cli_main.cpp` stays ten lines. */
int main_body(int argc, char** argv);

}  // namespace cli
}  // namespace line

#endif  // LINE_CLI_CLI_RUN_H
