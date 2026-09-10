/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Running progress log of a LINE solver run (the "solver console").
 *
 * <p>The console narrates what a solver is doing while it does it: reading the
 * model, compiling the network structure, computing chains, visits and demands,
 * resolving the method, iterating, and closing with the figures of merit. Each
 * line carries the elapsed time since the run started:</p>
 *
 * <pre>
 * [   0.014s] compiling the network structure of model 'cqn'
 * [   0.031s]   solving the routing DTMC for the visit ratios
 * [   0.052s] recognized a closed queueing network: 3 stations, 1 class, 1 chain
 * [   0.061s]   AMVA sweep 10: queue-length residual 1.11e-01, X = 1.2225
 * </pre>
 *
 * <p>It prints no tables: the result table stays the caller's own
 * {@code getAvgTable}. THE CONSOLE IS {@link VerboseLevel#DEBUG}: it narrates
 * exactly when the run is at DEBUG and is silent at every lower level, and it
 * never alters a numerical result. There is no separate console switch -- the
 * console was one, until it became clear that a running progress log IS what a
 * debug verbosity is for, and two switches for one channel only let a session
 * ask for DEBUG and get nothing.</p>
 *
 * <p>Typical use: {@code GlobalConstants.setVerbose(VerboseLevel.DEBUG)} for
 * the session, or {@code options.verbose = VerboseLevel.DEBUG} for one run.</p>
 *
 * <p>Nested runs (an inner solver driven by SolverLN, SolverENV or the UQ
 * ensemble) do not narrate: only the outermost run writes, so several inner
 * solves cannot interleave their lines. Use {@link #isActive} to suppress a
 * legacy print, {@link #ownsLog} to emit.</p>
 *
 * <p>Mirrors MATLAB {@code matlab/src/io/LineConsole.m} and python
 * {@code line_solver/api/io/console.py}.</p>
 */
public class LineConsole {

    private LineConsole() {
    }

    // ------------------------------------------------------------- state

    private static int depth = 0;             // nesting level of open runs
    private static boolean active = false;    // an outermost run is narrating
    private static String tag = "";           // solver name of the open run
    private static String modelName = "";     // model under study
    private static long t0 = 0L;              // nanoTime at the start of the run
    private static double tsetup = Double.NaN;// seconds before the analysis began
    private static int quiet = 0;             // >0 while compile detail is off
    private static int muted = 0;             // >0 while a silenced run executes
    private static boolean forceDetail = false;
    private static String lastLoop = "";
    private static int shown = 0;
    private static boolean trunc = false;
    private static long clock = 0L;           // console clock, restarted per run
    // legacy completion lines held until the run closes
    private static final List<String> pending = new ArrayList<String>();

    // routed line_debug bookkeeping
    private static String detailLast = "";
    private static final List<String> detailShapes = new ArrayList<String>();
    private static final List<Integer> detailShapeCount = new ArrayList<Integer>();
    private static int detailTotal = 0;

    private static final int MAX_ITER_LINES = 30;
    private static final int MAX_DETAIL_LINES = 200;
    private static final int MAX_PER_SHAPE = 3;

    /** Forget any open run (used after an interrupted solve). */
    public static synchronized void reset() {
        depth = 0;
        active = false;
        tag = "";
        modelName = "";
        t0 = 0L;
        tsetup = Double.NaN;
        quiet = 0;
        muted = 0;
        forceDetail = false;
        lastLoop = "";
        shown = 0;
        trunc = false;
        pending.clear();
        resetDetail();
    }

    private static void resetDetail() {
        detailLast = "";
        detailShapes.clear();
        detailShapeCount.clear();
        detailTotal = 0;
    }

    /**
     * Resolve whether a run with these options should narrate.
     *
     * <p>THE CONSOLE IS DEBUG, and this is the whole rule: a run narrates when
     * it is at {@link VerboseLevel#DEBUG} and at no lower level. The run's own
     * {@code options.verbose} decides when it carries one, otherwise the
     * session level does; there is no third switch that could put the two out
     * of step.</p>
     */
    public static synchronized boolean wanted(SolverOptions options) {
        if (options != null && options.verbose != null) {
            return options.verbose == VerboseLevel.DEBUG;
        }
        return GlobalConstants.getVerbose() == VerboseLevel.DEBUG;
    }

    /** @return true while a run is narrating; gates SUPPRESSION of legacy prints */
    public static synchronized boolean isActive() {
        return active;
    }

    /**
     * Queue a legacy line for after the run closes.
     *
     * <p>The console owns the log while a run narrates, so a solver's standard
     * completion message is held here and written once the closing DONE line
     * has gone out, reading as it would with the console off.</p>
     */
    public static synchronized void deferPrint(String fmt, Object... args) {
        final String text = String.format(fmt, args);
        if (!active) {
            System.out.print(text);
            return;
        }
        if (depth > 1) {
            return; // a nested run does not narrate, and does not report
        }
        pending.add(text);
    }

    /** @return true only inside the OUTERMOST open run; gates EMISSION */
    public static synchronized boolean ownsLog() {
        return active && depth <= 1 && muted == 0;
    }

    /**
     * @return true when a progress line should be printed: either the outermost
     * run is narrating, or no run is open and the session is at DEBUG -- the
     * second case is what lets model construction narrate before any solver
     * exists.
     */
    public static synchronized boolean writes() {
        if (muted > 0) {
            return false;
        }
        if (depth > 0) {
            return ownsLog();
        }
        return wanted(null);
    }

    // --------------------------------------------------------- lifecycle

    /**
     * Open a console run.
     *
     * <p>The caller MUST pair this with {@link #closeRun} in a finally block, so
     * that a failed analysis still reports what it had reached.</p>
     *
     * @param solver  the solver about to run (read for its model and result)
     * @param options its options
     * @return true when this call opened the outermost run
     */
    public static synchronized boolean beginRun(Object solver, SolverOptions options) {
        if (depth == 0 && !wanted(options)) {
            // A run that must stay silent MUTES the console for its whole
            // duration, so that nothing it triggers leaks out at depth 0.
            muted++;
            return false;
        }
        depth++;
        if (depth > 1) { // nested run: the outer analyzer owns the log
            return false;
        }
        active = true;
        tag = solverTag(solver);
        modelName = modelName(solver);
        t0 = System.nanoTime();
        tsetup = Double.NaN;
        lastLoop = "";
        shown = 0;
        trunc = false;
        resetDetail();
        openingLines(solver, options);
        tsetup = (System.nanoTime() - t0) / 1e9;
        return true;
    }

    /**
     * Close the innermost open run, writing the closing lines.
     *
     * @param solver the solver whose result is reported, or null
     */
    public static synchronized void closeRun(Object solver) {
        if (depth <= 0) {
            if (muted > 0) {
                muted--;
            }
            return;
        }
        depth--;
        if (depth > 0 || !active) {
            return;
        }
        if (solver != null) {
            closingLines(solver);
        }
        // the standard completion message follows the closing DONE line
        for (String line : pending) {
            System.out.print(line);
        }
        System.out.flush();
        // the mute of an enclosing silenced run outlives this run's own state
        final int outerMute = muted;
        reset();
        muted = outerMute;
    }

    // ---------------------------------------------------------- emission

    private static void emit(String indent, String text) {
        if (clock == 0L) {
            clock = System.nanoTime();
        }
        // a top-level row opens with a capital, an indented substep stays lowercase
        String row = text;
        if (indent.isEmpty() && !row.isEmpty()) {
            row = Character.toUpperCase(row.charAt(0)) + row.substring(1);
        }
        System.out.printf("[%8.3fs] %s%s%n", (System.nanoTime() - clock) / 1e9, indent, row);
    }

    /** Write one progress line. */
    public static synchronized void step(String fmt, Object... args) {
        if (!writes()) {
            return;
        }
        emit("", format(fmt, args));
    }

    /** Write one indented progress line. */
    public static synchronized void substep(String fmt, Object... args) {
        if (!writes()) {
            return;
        }
        emit("  ", format(fmt, args));
    }

    /**
     * One stage line of a structure compile.
     *
     * <p>Silenced inside a {@link #pushQuiet} scope, and inside an open run,
     * where the structures being compiled are those of auxiliary models
     * (SolverLN layers) rather than of the model under study.</p>
     */
    public static synchronized void compileDetail(String fmt, Object... args) {
        if (quiet > 0 || (depth > 0 && !forceDetail)) {
            return;
        }
        substep(fmt, args);
    }

    /** Announce the compilation of a model structure. */
    public static synchronized void compiling(String name) {
        if (depth > 0 && !name.equals(modelName)) {
            // an ensemble rebuilds the same submodel once per stage or per
            // iteration, so these go through detail() and collapse to one line
            detail("refreshing the auxiliary submodel '" + name + "'");
        } else {
            if (depth == 0) { // a compile outside any run opens its own timeline
                clock = System.nanoTime();
            }
            step("compiling the network structure of model '%s'", name);
        }
    }

    /**
     * Report a solver's own debug message as a substep.
     *
     * <p>{@code InputOutput.line_debug} routes here while a run narrates.
     * Consecutive repeats are dropped, at most three messages of the same SHAPE
     * (the text with its numbers masked) are reported, and the channel is capped
     * per run, since a message inside a loop would bury the narration.</p>
     */
    public static synchronized void detail(String text) {
        if (!ownsLog()) {
            return;
        }
        if (text == null) {
            return;
        }
        text = text.trim();
        if (text.isEmpty() || text.equals(detailLast)) {
            return;
        }
        String shape = text.replaceAll("[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?", "#");
        int hit = detailShapes.indexOf(shape);
        if (hit >= 0) {
            int n = detailShapeCount.get(hit) + 1;
            detailShapeCount.set(hit, n);
            if (n > MAX_PER_SHAPE) {
                return;
            }
        } else {
            detailShapes.add(shape);
            detailShapeCount.add(1);
        }
        detailLast = text;
        detailTotal++;
        if (detailTotal > MAX_DETAIL_LINES) {
            if (detailTotal == MAX_DETAIL_LINES + 1) {
                emit("  ", "further solver detail not reported");
            }
            return;
        }
        emit("  ", lowerFirst(text));
    }

    /**
     * Announce an iteration loop and reset its reporting budget.
     *
     * <p>Re-announcing the SAME text (a solver that restarts its loop) neither
     * reprints the header nor refills the budget.</p>
     */
    public static synchronized void loop(String fmt, Object... args) {
        if (!ownsLog()) {
            return;
        }
        String text = format(fmt, args);
        if (text.equals(lastLoop)) {
            return;
        }
        lastLoop = text;
        shown = 0;
        trunc = false;
        emit("", text);
    }

    /**
     * Report iteration {@code k} of the current loop.
     *
     * <p>Lines are decimated: the first 20 iterations report in full, then every
     * 10th, and the loop stops reporting after 30 lines, so that a long run
     * cannot bury the rest of the narration.</p>
     */
    public static synchronized void iter(long k, String fmt, Object... args) {
        if (!ownsLog()) {
            return;
        }
        if (k > 20 && k % 10 != 0) {
            return;
        }
        if (shown >= MAX_ITER_LINES) {
            if (!trunc) {
                trunc = true;
                emit("  ", "further iterations of this loop not reported");
            }
            return;
        }
        shown++;
        emit("  ", format(fmt, args));
    }

    /** Suppress structure-compile detail until {@link #popQuiet}. */
    public static synchronized void pushQuiet() {
        quiet++;
    }

    /** End one {@link #pushQuiet} scope. */
    public static synchronized void popQuiet() {
        quiet = Math.max(0, quiet - 1);
    }

    // -------------------------------------------------------- narration

    private static void openingLines(Object solver, SolverOptions options) {
        clock = System.nanoTime(); // each run's timeline starts at zero
        if (writes()) { // the opening row is set off from whatever preceded it
            System.out.println();
        }
        step("LINE %s: Solver%s starting on model '%s' (lang java)",
                GlobalConstants.Version, tag, modelName);
        readModel(solver);
        compileStruct(solver);
        recognizeModel(solver);
        reportMethod(solver, options);
    }

    private static void closingLines(Object solver) {
        Object res = field(solver, "result"); // jline.solvers.SolverResult
        double runtime = (System.nanoTime() - t0) / 1e9;
        if (res != null) {
            Object method = field(res, "method");
            Object iter = field(res, "iter");
            int iterations = iter instanceof Number ? ((Number) iter).intValue() : 0;
            String name = method == null ? "default" : method.toString();
            if (iterations > 1) {
                step("solved by %s in %d iterations", name, iterations);
            } else {
                step("solved by %s", name);
            }
            resultLines(structOf(solver), res);
        }
        if (Double.isNaN(tsetup)) {
            step("DONE in %.4f s", runtime);
        } else {
            step("DONE in %.4f s (setup %.4f s, analysis %.4f s)",
                    runtime, tsetup, Math.max(0.0, runtime - tsetup));
        }
    }

    private static void readModel(Object solver) {
        Object model = modelOf(solver);
        if (model == null) {
            return;
        }
        Object nnodes = call(model, "getNumberOfNodes");
        Object nclasses = call(model, "getNumberOfClasses");
        if (nnodes instanceof Number && nclasses instanceof Number) {
            step("reading the model: %s, %s",
                    plural(((Number) nnodes).intValue(), "node", "nodes"),
                    plural(((Number) nclasses).intValue(), "job class", "job classes"));
        }
    }

    private static void compileStruct(Object solver) {
        Object model = modelOf(solver);
        if (model == null) {
            return;
        }
        forceDetail = true; // this compile is the run's own model
        try {
            call(model, "getStruct");
        } finally {
            forceDetail = false;
        }
    }

    private static void recognizeModel(Object solver) {
        Object sn = structOf(solver);
        if (sn == null) {
            return;
        }
        Integer nstations = intField(sn, "nstations");
        Integer nclasses = intField(sn, "nclasses");
        Integer nchains = intField(sn, "nchains");
        if (nstations == null || nclasses == null || nchains == null) {
            return;
        }
        String kind = modelKind(sn, nclasses.intValue());
        step("recognized %s %s: %s, %s, %s", article(kind), kind,
                plural(nstations.intValue(), "station", "stations"),
                plural(nclasses.intValue(), "class", "classes"),
                plural(nchains.intValue(), "chain", "chains"));
        String sched = schedMix(sn, nstations.intValue());
        if (sched != null && !sched.isEmpty()) {
            substep("scheduling: %s", sched);
        }
    }

    private static void reportMethod(Object solver, SolverOptions options) {
        if (options == null) {
            return;
        }
        String requested = options.method == null ? "default" : options.method;
        // %g pads a tolerance to nine digits; the plain string form reads as
        // the user wrote it
        step("method '%s', tolerance %s, iteration cap %d",
                requested, Double.toString(options.tol), options.iter_max);
        if (tag.equals("SSA") || tag.equals("LDES") || tag.equals("JMT")) {
            substep("sample budget %d, seed %d", options.samples, options.seed);
        }
    }

    private static void resultLines(Object sn, Object res) {
        if (res == null) {
            return;
        }
        Matrix q = matrixField(res, "QN");
        Matrix u = matrixField(res, "UN");
        Matrix x = matrixField(res, "XN");
        if (x != null && x.length() > 0) {
            substep("system throughput %.4f", x.elementSum());
        }
        if (u != null && u.getNumRows() > 0 && sn != null) {
            // Only the first nstations rows are stations: a fork-join model's
            // result carries the transformed network's extra rows (and an FCR
            // model appends region rows) past that point, and indexing the
            // struct with one of those is out of range.
            Integer nst = intField(sn, "nstations");
            final int rows = nst == null ? u.getNumRows() : Math.min(u.getNumRows(), nst.intValue());
            int best = -1;
            double bestVal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < rows; i++) {
                if (isDelayOrSource(sn, i)) {
                    continue; // neither saturates, so neither is a candidate
                }
                double row = 0;
                for (int r = 0; r < u.getNumCols(); r++) {
                    row += u.get(i, r);
                }
                if (row > bestVal) {
                    bestVal = row;
                    best = i;
                }
            }
            if (best >= 0) {
                substep("busiest queueing station %s at utilization %.4f", stationName(sn, best), bestVal);
            }
        }
        if (q != null && q.length() > 0) {
            substep("mean jobs in the network %.4f", q.elementSum());
        }
    }

    // ---------------------------------------------------------- helpers

    private static String solverTag(Object solver) {
        String name = solver == null ? "" : solver.getClass().getSimpleName();
        return name.startsWith("Solver") ? name.substring(6) : name;
    }

    private static Object modelOf(Object solver) {
        Object model = field(solver, "model");
        return model;
    }

    private static String modelName(Object solver) {
        Object model = modelOf(solver);
        if (model == null) {
            return "(unnamed)";
        }
        Object name = call(model, "getName");
        return name == null ? "(unnamed)" : name.toString();
    }

    private static Object structOf(Object solver) {
        Object model = modelOf(solver);
        if (model == null) {
            return null;
        }
        return call(model, "getStruct");
    }

    private static String modelKind(Object sn, int nclasses) {
        Matrix njobs = matrixField(sn, "njobs");
        String base = "queueing network";
        int nopen = 0;
        if (njobs != null) {
            for (int i = 0; i < njobs.length(); i++) {
                if (Double.isInfinite(njobs.get(i))) {
                    nopen++;
                }
            }
        }
        if (nopen == 0) {
            return "closed " + base;
        }
        if (nopen == nclasses) {
            return "open " + base;
        }
        return "mixed " + base;
    }

    private static String schedMix(Object sn, int nstations) {
        Object sched = field(sn, "sched");
        if (!(sched instanceof java.util.Map)) {
            return "";
        }
        java.util.Map<?, ?> map = (java.util.Map<?, ?>) sched;
        List<String> names = new ArrayList<String>();
        List<Integer> counts = new ArrayList<Integer>();
        for (Object value : map.values()) {
            String name = String.valueOf(value).toUpperCase();
            int idx = names.indexOf(name);
            if (idx < 0) {
                names.add(name);
                counts.add(1);
            } else {
                counts.set(idx, counts.get(idx) + 1);
            }
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < names.size(); i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(names.get(i)).append(" x").append(counts.get(i));
        }
        return sb.toString();
    }

    private static boolean isDelayOrSource(Object sn, int stationIndex) {
        // sched is keyed by Station, so the station INDEX is resolved through
        // the struct's own station list rather than by map iteration order,
        // which carries no index meaning.
        Object stations = field(sn, "stations");
        Object sched = field(sn, "sched");
        if (!(sched instanceof java.util.Map) || !(stations instanceof List)) {
            return false;
        }
        List<?> list = (List<?>) stations;
        if (stationIndex < 0 || stationIndex >= list.size()) {
            return false;
        }
        Object value = ((java.util.Map<?, ?>) sched).get(list.get(stationIndex));
        return value == SchedStrategy.INF || value == SchedStrategy.EXT;
    }

    private static String stationName(Object sn, int stationIndex) {
        Matrix toNode = matrixField(sn, "stationToNode");
        Object names = field(sn, "nodenames");
        if (toNode != null && stationIndex >= 0 && stationIndex < toNode.length()
                && names instanceof List) {
            int node = (int) toNode.get(stationIndex);
            List<?> list = (List<?>) names;
            if (node >= 0 && node < list.size()) {
                return String.valueOf(list.get(node));
            }
        }
        return "station " + (stationIndex + 1);
    }

    private static String lowerFirst(String s) {
        if (s.length() >= 2 && !s.substring(0, 2).equals(s.substring(0, 2).toUpperCase())) {
            return Character.toLowerCase(s.charAt(0)) + s.substring(1);
        }
        return s;
    }

    private static String article(String word) {
        char c = Character.toLowerCase(word.charAt(0));
        return "aeiou".indexOf(c) >= 0 ? "an" : "a";
    }

    /** Count with an agreeing noun, e.g. "1 chain" or "3 chains". */
    public static String plural(int n, String singular, String plural) {
        return n + " " + (n == 1 ? singular : plural);
    }

    private static String format(String fmt, Object... args) {
        if (args == null || args.length == 0) {
            return fmt;
        }
        return String.format(fmt, args);
    }

    private static String str(Object o) {
        return o == null ? null : o.toString();
    }

    // Reflection keeps this class independent of the solver hierarchy: it is
    // called from jline.io, which the solver packages depend on and not the
    // other way round.
    private static Object call(Object target, String method) {
        if (target == null) {
            return null;
        }
        try {
            java.lang.reflect.Method m = target.getClass().getMethod(method);
            m.setAccessible(true);
            return m.invoke(target);
        } catch (Exception e) {
            return null;
        }
    }

    private static Object field(Object target, String name) {
        if (target == null) {
            return null;
        }
        try {
            java.lang.reflect.Field f = target.getClass().getField(name);
            return f.get(target);
        } catch (Exception e) {
            return null;
        }
    }

    private static Integer intField(Object target, String name) {
        Object v = field(target, name);
        if (v instanceof Number) {
            return Integer.valueOf(((Number) v).intValue());
        }
        return null;
    }

    private static Matrix matrixField(Object target, String name) {
        Object v = field(target, name);
        if (v instanceof Matrix) {
            return (Matrix) v;
        }
        Object viaGetter = call(target, "get" + Character.toUpperCase(name.charAt(0)) + name.substring(1));
        return viaGetter instanceof Matrix ? (Matrix) viaGetter : null;
    }
}
