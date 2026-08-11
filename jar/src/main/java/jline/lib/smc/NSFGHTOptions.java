package jline.lib.smc;

/** Options for NSF_GHT solver. */
public final class NSFGHTOptions {
    private final int maxNumIt;
    private final int verbose;
    private final boolean firstBlockRow;

    public NSFGHTOptions(int maxNumIt, int verbose, boolean firstBlockRow) {
        this.maxNumIt = maxNumIt;
        this.verbose = verbose;
        this.firstBlockRow = firstBlockRow;
    }

    public NSFGHTOptions() {
        this(10000, 0, false);
    }

    public NSFGHTOptions(int verbose) {
        this(10000, verbose, false);
    }

    public int getMaxNumIt() { return maxNumIt; }
    public int getVerbose() { return verbose; }
    public boolean isFirstBlockRow() { return firstBlockRow; }
}
