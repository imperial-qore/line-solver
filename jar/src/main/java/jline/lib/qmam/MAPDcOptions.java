package jline.lib.qmam;

/** Options for MAP/D/c queue analysis. */
public final class MAPDcOptions {
    private final int maxNumComp;
    private final int verbose;
    private final int numSteps;

    public MAPDcOptions(int maxNumComp, int verbose, int numSteps) {
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
        this.numSteps = numSteps;
    }

    public MAPDcOptions() {
        this(1000, 0, 1);
    }

    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }
    public int getNumSteps() { return numSteps; }
}
