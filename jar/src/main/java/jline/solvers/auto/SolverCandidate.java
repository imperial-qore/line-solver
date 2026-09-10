package jline.solvers.auto;

import java.util.ArrayList;
import java.util.List;

/**
 * One row of {@link SolverAUTO#findSolver}: a (family, method) pair this model
 * can be asked for, whether it runs, what kind of answer it returns and which
 * measures it can report.
 *
 * <p>The fields are the columns of the MATLAB table returned by
 * {@code @SolverAUTO/findSolver.m} and of the pandas frame returned by the
 * native python twin, under the same names, so a report can be compared across
 * the three codebases row for row.
 */
public final class SolverCandidate {

    /** Answer kinds a method can return; see {@link SolverAUTO#methodClass}. */
    public static final String CLASS_EXACT = "exact";
    /** An approximation of the model's true measures. */
    public static final String CLASS_APPROX = "approx";
    /** A bracket rather than an estimate. */
    public static final String CLASS_BOUND = "bound";
    /** A seed-dependent estimate from a sample path. */
    public static final String CLASS_SIMULATION = "simulation";

    /** The method family: "mva", "ctmc", "ldes", ... */
    public final String solver;
    /** The method name to pass as a method name, e.g. "mva.exact". */
    public final String method;
    /** True when the model passes this method's own support gate. */
    public final boolean runnable;
    /** One of the four CLASS_ constants. */
    public final String methodClass;
    /** The measure groups the family answers, comma-joined. */
    public final String metrics;
    /** Why a refused pair was refused; empty when {@link #runnable}. */
    public final String reason;

    /**
     * @param solver      the method family
     * @param method      the qualified method name
     * @param runnable    whether the model passes this method's gate
     * @param methodClass one of the CLASS_ constants
     * @param metrics     comma-joined measure groups
     * @param reason      the refusal reason, empty when runnable
     */
    public SolverCandidate(String solver, String method, boolean runnable, String methodClass,
                           String metrics, String reason) {
        this.solver = solver;
        this.method = method;
        this.runnable = runnable;
        this.methodClass = methodClass;
        this.metrics = metrics;
        this.reason = reason;
    }

    /** @return the method family */
    public String getSolver() {
        return solver;
    }

    /** @return the qualified method name */
    public String getMethod() {
        return method;
    }

    /** @return whether the model passes this method's gate */
    public boolean isRunnable() {
        return runnable;
    }

    /** @return one of the CLASS_ constants */
    public String getMethodClass() {
        return methodClass;
    }

    /** @return the comma-joined measure groups */
    public String getMetrics() {
        return metrics;
    }

    /** @return the refusal reason, empty when runnable */
    public String getReason() {
        return reason;
    }

    @Override
    public String toString() {
        return method + " [" + methodClass + (runnable ? "" : "; refused: " + reason) + "]";
    }

    /**
     * The rows as an aligned text table, the form the CLI and a console caller
     * want. The Reason column is last and unpadded, since it is the only one
     * whose width is unbounded.
     *
     * @param rows the rows to render, may be empty
     * @return the table as a printable string, with a trailing newline
     */
    public static String toTable(List<SolverCandidate> rows) {
        if (rows == null || rows.isEmpty()) {
            return "No solver method can analyze this model.\n";
        }
        List<String[]> cells = new ArrayList<String[]>();
        cells.add(new String[]{"Solver", "Method", "Runnable", "Class", "Metrics", "Reason"});
        for (SolverCandidate r : rows) {
            cells.add(new String[]{r.solver, r.method, r.runnable ? "true" : "false",
                    r.methodClass, r.metrics, r.reason});
        }
        int[] width = new int[6];
        for (String[] row : cells) {
            for (int c = 0; c < 5; c++) {
                if (row[c] != null && row[c].length() > width[c]) {
                    width[c] = row[c].length();
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        for (String[] row : cells) {
            for (int c = 0; c < 5; c++) {
                sb.append(pad(row[c], width[c])).append("  ");
            }
            sb.append(row[5] == null ? "" : row[5]);
            // trailing blanks on a reason-less row would be invisible padding
            while (sb.length() > 0 && sb.charAt(sb.length() - 1) == ' ') {
                sb.setLength(sb.length() - 1);
            }
            sb.append('\n');
        }
        return sb.toString();
    }

    private static String pad(String s, int width) {
        String v = s == null ? "" : s;
        StringBuilder sb = new StringBuilder(v);
        while (sb.length() < width) {
            sb.append(' ');
        }
        return sb.toString();
    }
}
