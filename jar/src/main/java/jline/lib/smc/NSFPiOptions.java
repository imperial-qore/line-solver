package jline.lib.smc;

/** Options for NSF_pi solver. */
public final class NSFPiOptions {
    private final int maxNumComp;
    private final boolean verbose;
    private final boolean firstBlockRow;

    public NSFPiOptions(int maxNumComp, boolean verbose, boolean firstBlockRow) {
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
        this.firstBlockRow = firstBlockRow;
    }

    public NSFPiOptions() {
        this(1000, false, false);
    }

    public NSFPiOptions(int maxNumComp, boolean verbose) {
        this(maxNumComp, verbose, false);
    }

    public int getMaxNumComp() { return maxNumComp; }
    public boolean isVerbose() { return verbose; }
    public boolean isFirstBlockRow() { return firstBlockRow; }
}
