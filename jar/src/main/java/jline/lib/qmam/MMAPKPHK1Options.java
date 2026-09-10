package jline.lib.qmam;

/** Options for MMAP[K]/PH[K]/1 queue analysis. */
public final class MMAPKPHK1Options {
    public final String mode;
    public final int maxNumComp;
    public final int verbose;

    public MMAPKPHK1Options() { this("Sylves", 1000, 0); }

    public MMAPKPHK1Options(String mode, int maxNumComp, int verbose) {
        this.mode = mode;
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
    }

    public String getMode() { return mode; }
    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }
}
