package jline.util;

/**
 * Session-level cooperative wall-clock budget (port of MATLAB
 * lineTimeoutExceeded's global-deadline form). A per-solver runAnalyzer with a
 * finite options.timeout arms the deadline and clears it in a finally block;
 * deep utility loops with no options argument (e.g. multichoose during CTMC
 * state-space generation) poll {@link #checkpoint(String)} at safe points and
 * abort with a RuntimeException, since a partially enumerated state space
 * would silently yield wrong results.
 */
public final class LineTimeout {

    private static volatile long deadlineNanos = 0L; // 0 = no budget armed
    private static int tochk = 0; // amortization counter (approximate under races, by design)

    private LineTimeout() {
    }

    /** Arm the deadline; a non-finite or non-positive budget disarms it. */
    public static void set(double budgetSeconds) {
        if (Double.isInfinite(budgetSeconds) || Double.isNaN(budgetSeconds) || budgetSeconds <= 0) {
            deadlineNanos = 0L;
        } else {
            deadlineNanos = System.nanoTime() + (long) (budgetSeconds * 1e9);
        }
    }

    /** Disarm the deadline. */
    public static void clear() {
        deadlineNanos = 0L;
    }

    /** True if a deadline is armed and has passed. */
    public static boolean exceeded() {
        long d = deadlineNanos;
        return d != 0L && System.nanoTime() > d;
    }

    /**
     * Amortized checkpoint for hot loops: polls the clock once every 256
     * calls and throws when the armed deadline has passed. The stride is low
     * because callers do heavy matrix work between calls (multichoose
     * concatenations), so a large stride overshoots the budget by minutes.
     */
    public static void checkpoint(String context) {
        if (++tochk >= 256) {
            tochk = 0;
            if (exceeded()) {
                throw new RuntimeException(context
                        + " exceeded the wall-clock time budget (options.timeout).");
            }
        }
    }
}
