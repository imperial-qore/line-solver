/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

/**
 * Capture the result tables a run produces, with the solver that produced them.
 *
 * <p>WHY THIS EXISTS. Cross-codebase parity is asserted against one shared
 * golden per example ({@code goldens/baselines/*.json}), keyed by SOLVER NAME.
 * The only way to recover that key used to be scraping the banner an example
 * printed above each table -- one regex dialect per codebase, and a value
 * truncated to the digits the printer showed. The recorder supplies the same
 * attribution BY CONSTRUCTION and at full precision. It is the Java twin of
 * {@code python/line_solver/result_recorder.py},
 * {@code matlab/src/io/LineResultRecorder.m} and
 * {@code cpp/examples/parity_recorder.h}.
 *
 * <p>WHERE IT HOOKS. Every result table in the JAR is materialised by one of a
 * dozen public getters on {@link jline.solvers.NetworkSolver} and its ensemble
 * counterparts, and each of those wraps its body in {@link #around}. So an
 * example that merely calls {@code solver.getAvgTable()} is recorded without
 * being touched -- which is what lets the 100-odd twins already in
 * {@code jline.examples.java} answer for the JAVA parity row unchanged.
 *
 * <p>ONE ENTRY PER OUTERMOST CALL. The ensemble solvers drive a layer solver's
 * getter of the same name, and recording both would file every layer's table
 * beside the ensemble's own. {@link #around} counts the nesting and records only
 * as the outermost call returns.
 *
 * <p>IT IS OFF unless {@link #enable} was called or {@code LINE_RECORD_RESULTS=1}
 * is set, and then it costs one boolean test per getter.
 *
 * <p>THE STATE IS PROCESS-WIDE AND UNSYNCHRONISED, which is correct for what
 * asks for it: a twin runs one example in one process on one thread, and the
 * nesting count only means anything within that thread. A caller that drove
 * several solves concurrently in one process would interleave their records; do
 * not enable it there.
 */
public final class LineResultRecorder {

    /** A getter body, so a caller can wrap one in {@link #around} in a single line. */
    public interface Body<T> {
        T get();
    }

    /** One recorded table, with what produced it. */
    private static final class Record {
        String solver;
        String method;
        String view;
        String rowLabel;
        String colLabel;
        boolean derived;
        long seq;
        final List<Map<String, Object>> rows = new ArrayList<Map<String, Object>>();
    }

    /**
     * The metrics a golden may carry.
     *
     * <p>Measured over all 209 baselines, these six and the two label columns are
     * the ONLY keys any golden holds, so a table is recorded through them and a
     * column no golden carries costs nothing to leave out.
     */
    private static final String[] METRICS = {"QLen", "Util", "RespT", "ResidT", "ArvR", "Tput"};

    /** Accessor on the table for each entry of {@link #METRICS}, in the same order. */
    private static final String[] METRIC_GETTERS = {"getQLen", "getUtil", "getRespT",
            "getResidT", "getArvR", "getTput"};

    private static boolean on = false;
    private static String path = null;
    private static int depth = 0;
    private static long seq = 0;
    private static final List<Record> RECORDS = new ArrayList<Record>();
    private static final List<String> NOTES = new ArrayList<String>();

    private LineResultRecorder() {
    }

    /** True when a run asked to record. Nothing below costs anything when false. */
    public static boolean isEnabled() {
        return on;
    }

    /** Turn recording on and direct the dump at {@code file}. */
    public static void enable(String file) {
        on = true;
        path = file;
        depth = 0;
        seq = 0;
        RECORDS.clear();
        NOTES.clear();
    }

    /** Turn recording on when {@code LINE_RECORD_RESULTS=1} asked for it. */
    public static void enableFromEnv() {
        String flag = System.getenv("LINE_RECORD_RESULTS");
        if (flag == null || !flag.equals("1")) {
            return;
        }
        String out = System.getenv("LINE_RECORD_PATH");
        enable(out == null ? "line-record.json" : out);
    }

    /**
     * Run one result getter, recording what it returned.
     *
     * <p>The nesting count is what keeps an ensemble solve one record rather than
     * one per layer: only the outermost getter files a table.
     *
     * @param solver the solver whose getter this is
     * @param view   which table this is: {@code avg}, {@code node}, {@code chain}, ...
     * @param body   the getter's own body
     * @param <T>    the table type the getter returns
     * @return exactly what {@code body} returned
     */
    public static <T> T around(Object solver, String view, Body<T> body) {
        if (!on) {
            return body.get();
        }
        depth++;
        T table;
        try {
            table = body.get();
        } finally {
            depth--;
        }
        if (depth == 0 && table != null) {
            capture(solver, view, table);
        }
        return table;
    }

    /**
     * Record a labelled scalar the example DERIVED.
     *
     * <p>Some goldens hold a quantity no result table carries -- a state
     * probability, a phase-type moment, a cache hit rate. They are recorded as
     * one-row tables so everything downstream handles a single shape, and marked
     * {@code derived} so the comparator scores them at the precision the example's
     * own format string wrote rather than at a table's five printed digits.
     *
     * @param solver the golden's key for this quantity
     * @param row    the golden's first label column for it
     * @param col    the golden's second label column for it
     * @param value  the quantity
     */
    public static void scalar(String solver, String row, String col, double value) {
        scalar(solver, row, col, "QLen", value);
    }

    /**
     * The same, under the metric column the golden files it in.
     *
     * <p>Almost every derived golden is filed under {@code QLen}, which is what
     * the scraper that generated them used for a labelled scalar. A few are not
     * -- {@code cdf_respt_populations} holds response times under {@code RespT} --
     * and the caller must name the column the golden actually uses, because the
     * comparator joins on it.
     */
    public static void scalar(String solver, String row, String col, String metric,
                              double value) {
        if (!on) {
            return;
        }
        Record rec = new Record();
        rec.solver = solver;
        rec.method = "default";
        rec.view = "scalar";
        rec.rowLabel = "Station";
        rec.colLabel = "JobClass";
        rec.derived = true;
        rec.seq = seq++;
        Map<String, Object> entry = new LinkedHashMap<String, Object>();
        entry.put("Station", row);
        entry.put("JobClass", col);
        entry.put(metric, Double.valueOf(value));
        rec.rows.add(entry);
        RECORDS.add(rec);
    }

    /**
     * Record a bare scalar an example prints with no label at all.
     *
     * <p>The goldens key these {@code OPT / OptResult / Value}, which is the shape
     * the scraper that generated them read a lone float in; recording that shape
     * is a transcription of the golden and not a new naming scheme.
     */
    public static void bareScalar(double value) {
        scalar("OPT", "OptResult", "Value", value);
    }

    /**
     * Record a refusal the run made BY NAME.
     *
     * <p>A refusal is a fact about the port and is a named skip downstream; a
     * table that simply never arrived is a failure. The two must not look alike,
     * which is why this is recorded rather than only printed.
     */
    public static void note(String solver, String why) {
        if (on) {
            NOTES.add(solver + ": " + why);
        }
    }

    /** Everything recorded so far, as the wire document. Public for tests. */
    public static JsonObject document() {
        JsonObject out = new JsonObject();
        JsonArray records = new JsonArray();
        for (int i = 0; i < RECORDS.size(); i++) {
            records.add(toJson(RECORDS.get(i)));
        }
        out.add("records", records);
        JsonArray notes = new JsonArray();
        for (int i = 0; i < NOTES.size(); i++) {
            notes.add(NOTES.get(i));
        }
        out.add("notes", notes);
        return out;
    }

    /** Write the JSON dump to the path {@link #enable} was given. No-op when off. */
    public static boolean dump() {
        if (!on) {
            return true;
        }
        Writer fh = null;
        try {
            fh = new FileWriter(path);
            new GsonBuilder().serializeNulls().create().toJson(document(), fh);
            return true;
        } catch (IOException e) {
            System.err.println("parity recorder: cannot write " + path + ": " + e.getMessage());
            return false;
        } finally {
            if (fh != null) {
                try {
                    fh.close();
                } catch (IOException ignored) {
                    // the dump is already written or already lost; nothing to add here
                }
            }
        }
    }

    private static JsonObject toJson(Record rec) {
        JsonObject jr = new JsonObject();
        jr.addProperty("solver", rec.solver);
        jr.addProperty("method", rec.method);
        jr.addProperty("view", rec.view);
        JsonArray labels = new JsonArray();
        labels.add(rec.rowLabel);
        labels.add(rec.colLabel);
        jr.add("labels", labels);
        jr.addProperty("seq", Long.valueOf(rec.seq));
        jr.addProperty("derived", Boolean.valueOf(rec.derived));
        JsonArray rows = new JsonArray();
        for (int i = 0; i < rec.rows.size(); i++) {
            JsonObject jrow = new JsonObject();
            Map<String, Object> row = rec.rows.get(i);
            for (Map.Entry<String, Object> cell : row.entrySet()) {
                Object v = cell.getValue();
                if (v instanceof Double) {
                    double d = ((Double) v).doubleValue();
                    // A NON-FINITE CELL IS WRITTEN AS A STRING, not as bare NaN, which is
                    // not JSON and which every reader downstream rejects as malformed
                    // rather than as a missing value.
                    if (Double.isNaN(d)) {
                        jrow.addProperty(cell.getKey(), "NaN");
                    } else if (Double.isInfinite(d)) {
                        jrow.addProperty(cell.getKey(), d > 0 ? "Infinity" : "-Infinity");
                    } else {
                        jrow.addProperty(cell.getKey(), Double.valueOf(d));
                    }
                } else {
                    jrow.addProperty(cell.getKey(), String.valueOf(v));
                }
            }
            rows.add(jrow);
        }
        jr.add("rows", rows);
        return jr;
    }

    // ---------------------------------------------------------------- capture

    /**
     * The two label columns of each view, in the order the goldens carry them.
     *
     * <p>The comparator joins a recorded row to a golden row BY POSITION, so what
     * matters is which accessor supplies the first label and which the second; the
     * names are carried for the failure message.
     */
    private static String[] labelSpec(String view) {
        if (view.equals("node") || view.equals("cache") || view.equals("item")) {
            return new String[]{"Node", "getNodeNames", "JobClass", "getClassNames"};
        }
        if (view.equals("chain")) {
            return new String[]{"Station", "getStationNames", "Chain", "getChainNames"};
        }
        if (view.equals("nodechain")) {
            return new String[]{"Node", "getNodeNames", "Chain", "getChainNames"};
        }
        if (view.equals("sys")) {
            return new String[]{"Chain", "getChainNames", "JobClass", "getInChainNames"};
        }
        if (view.equals("layered")) {
            // A layered table names its nodes and nothing else; every codebase prints
            // the class column of an LQN row as 'Jobs', which is the golden's spelling.
            return new String[]{"Station", "getNodeNames", "JobClass", null};
        }
        return new String[]{"Station", "getStationNames", "JobClass", "getClassNames"};
    }

    private static void capture(Object solver, String view, Object table) {
        String label = solverLabel(solver);
        if (label == null) {
            return;
        }
        String[] spec = labelSpec(view);
        List<String> rowNames = strings(table, spec[1]);
        List<String> colNames = strings(table, spec[3]);
        if (rowNames == null) {
            return;
        }
        List<List<Double>> columns = new ArrayList<List<Double>>();
        List<String> names = new ArrayList<String>();
        for (int m = 0; m < METRICS.length; m++) {
            List<Double> col = doubles(table, METRIC_GETTERS[m]);
            if (col != null) {
                columns.add(col);
                names.add(METRICS[m]);
            }
        }
        if (columns.isEmpty()) {
            return;
        }
        Record rec = new Record();
        rec.solver = label;
        rec.method = solverMethod(solver);
        rec.view = view;
        rec.rowLabel = spec[0];
        rec.colLabel = spec[2];
        rec.derived = false;
        rec.seq = seq++;
        for (int i = 0; i < rowNames.size(); i++) {
            Map<String, Object> row = new LinkedHashMap<String, Object>();
            row.put(spec[0], rowNames.get(i));
            row.put(spec[2], colNames == null || i >= colNames.size() ? "Jobs" : colNames.get(i));
            for (int c = 0; c < columns.size(); c++) {
                List<Double> col = columns.get(c);
                if (i < col.size() && col.get(i) != null) {
                    row.put(names.get(c), col.get(i));
                }
            }
            rec.rows.add(row);
        }
        if (rec.rows.isEmpty()) {
            return;
        }
        RECORDS.add(rec);
    }

    @SuppressWarnings("unchecked")
    private static List<String> strings(Object table, String getter) {
        Object v = call(table, getter);
        return v instanceof List ? (List<String>) v : null;
    }

    @SuppressWarnings("unchecked")
    private static List<Double> doubles(Object table, String getter) {
        Object v = call(table, getter);
        return v instanceof List ? (List<Double>) v : null;
    }

    /**
     * One no-argument accessor, or null when this table does not carry it.
     *
     * <p>Reflection rather than a cast per table class: the tables do not share an
     * interface that declares these, and a hand-kept list of casts is a list that
     * silently stops covering a table type added later.
     */
    private static Object call(Object table, String getter) {
        if (getter == null) {
            return null;
        }
        try {
            Method m = table.getClass().getMethod(getter);
            return m.invoke(table);
        } catch (Exception e) {
            return null;
        }
    }

    // ------------------------------------------------------------ solver keys

    /**
     * The golden's key for this solver, or null when it is not one we key by.
     *
     * <p>A LAYERED or ENVIRONMENT solve is qualified by the member solver it drove
     * -- {@code LN(NC)}, {@code ENV(FLD)} -- because that is how several goldens
     * spell it and because the two are genuinely different computations. The
     * comparator reconciles the qualified and bare spellings against the golden's
     * own key; recording the member is what gives it the evidence to do so safely.
     */
    public static String solverLabel(Object solver) {
        String base = baseLabel(solver.getClass());
        if (base == null) {
            return null;
        }
        if (!base.equals("LN") && !base.equals("ENV") && !base.equals("UQ")) {
            return base;
        }
        // AN ENSEMBLE'S OWN METHOD NAMES IT WHEN IT HAS ONE, because that is the
        // distinction the goldens draw: lqn_moment3 solves one model with the same NC
        // layers twice, default and moment3, and its golden holds the first.
        String method = solverMethod(solver);
        if (method != null && !method.equals("default") && !isLabel(method.toUpperCase())) {
            return base + "(" + method + ")";
        }
        String member = memberLabel(solver);
        return member == null ? base : base + "(" + member + ")";
    }

    /** Whether a name is already one of the labels a golden keys a solver by. */
    private static boolean isLabel(String name) {
        return baseLabelOf(name) != null;
    }

    /**
     * The label for a solver class, walking up to a parent when the class is an alias.
     *
     * <p>{@code jline.solvers.mva.MVA} is the short alias of {@code SolverMVA}, and
     * {@code LINE} is the alias of {@code SolverAUTO}; both key as their parent.
     */
    private static String baseLabel(Class<?> klass) {
        for (Class<?> k = klass; k != null; k = k.getSuperclass()) {
            String label = baseLabelOf(k.getSimpleName());
            if (label != null) {
                return label;
            }
        }
        return null;
    }

    private static String baseLabelOf(String simpleName) {
        String name = simpleName.startsWith("Solver") ? simpleName.substring(6) : simpleName;
        if (name.equals("MVA") || name.equals("NC") || name.equals("CTMC") || name.equals("SSA")
                || name.equals("MAM") || name.equals("JMT") || name.equals("LDES")
                || name.equals("LQNS") || name.equals("QNS") || name.equals("BA")
                || name.equals("AG") || name.equals("AUTO") || name.equals("LN")
                || name.equals("ENV") || name.equals("UQ")) {
            return name;
        }
        if (name.equals("FLD") || name.equals("Fluid")) {
            return "FLD";
        }
        if (name.equals("LINE")) {
            return "AUTO";
        }
        return null;
    }

    /** The method this solver resolved, or {@code default} when it pinned none. */
    public static String solverMethod(Object solver) {
        try {
            Object options = solver.getClass().getField("options").get(solver);
            if (options == null) {
                return "default";
            }
            Object method = options.getClass().getField("method").get(options);
            String text = method == null ? null : String.valueOf(method);
            return text == null || text.isEmpty() ? "default" : text;
        } catch (Exception e) {
            return "default";
        }
    }

    /**
     * The layer or stage solver an ensemble actually ran, as a golden label.
     *
     * <p>THE FACTORY IS THE DECLARATION and it hides the class, so the answer is
     * read off the solvers the ensemble BUILT rather than off the factory. Getting
     * it wrong is not a spelling difference: MVA layers and NC layers are different
     * fixed points.
     */
    private static String memberLabel(Object solver) {
        // EACH ENSEMBLE SPELLS THE FIELD ITS OWN WAY, AND ONE OF THEM LIES.
        // SolverLN fills the inherited `solvers`; SolverENV keeps its stage
        // solvers in `envSolvers` and leaves `solvers` allocated but FULL OF
        // NULLS until the macro-ensemble path assigns it. So the choice cannot be
        // made on the field's presence, or even on its length -- an ENV solve
        // picked the null-filled `solvers`, found no member in it, and went out
        // labelled bare "ENV", which aliases to no golden key and reported the
        // golden's FLD row missing. Every candidate is scanned, and the first
        // one holding an actual solver wins.
        String[] names = {"solvers", "envSolvers"};
        for (int n = 0; n < names.length; n++) {
            for (Class<?> k = solver.getClass(); k != null; k = k.getSuperclass()) {
                String label = labelFromField(solver, k, names[n]);
                if (label != null) {
                    return label;
                }
            }
        }
        return null;
    }

    /** The label of the first non-null solver in {@code klass.field}, or null. */
    private static String labelFromField(Object solver, Class<?> klass, String name) {
        try {
            java.lang.reflect.Field field = klass.getDeclaredField(name);
            field.setAccessible(true);
            Object value = field.get(solver);
            if (value == null || !value.getClass().isArray()) {
                return null;
            }
            int n = java.lang.reflect.Array.getLength(value);
            for (int i = 0; i < n; i++) {
                Object member = java.lang.reflect.Array.get(value, i);
                if (member != null) {
                    return baseLabel(member.getClass());
                }
            }
            return null;
        } catch (Exception e) {
            return null;
        }
    }
}
