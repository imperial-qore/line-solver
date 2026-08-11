package jline.lib.qmam;

/** Options for MAP/M/c queue analysis. */
public final class MAPMcOptions {
    public final String mode;
    public final int maxNumComp;
    public final int verbose;

    public MAPMcOptions() { this("SylvesCR", 1000, 0); }

    public MAPMcOptions(String mode, int maxNumComp, int verbose) {
        this.mode = mode;
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
    }

    public String getMode() { return mode; }
    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }
}
