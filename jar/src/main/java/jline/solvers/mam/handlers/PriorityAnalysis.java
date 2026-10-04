package jline.solvers.mam.handlers;

/** Analysis result indicating priority configuration support level. */
public final class PriorityAnalysis {
    public final boolean isIdentical;
    public final boolean isAllDistinct;
    public final boolean isSupported;
    public final String message;

    public PriorityAnalysis(boolean isIdentical, boolean isAllDistinct, boolean isSupported, String message) {
        this.isIdentical = isIdentical;
        this.isAllDistinct = isAllDistinct;
        this.isSupported = isSupported;
        this.message = message;
    }
}
