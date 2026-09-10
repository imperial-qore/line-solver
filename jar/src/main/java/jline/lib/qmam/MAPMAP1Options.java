/**
 * @file Options for MAP/MAP/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Objects;

/**
 * Options for MAP/MAP/1 queue analysis.
 */
public final class MAPMAP1Options {
    public final String mode;
    public final int maxNumComp;
    public final int verbose;

    public MAPMAP1Options(String mode, int maxNumComp, int verbose) {
        this.mode = mode;
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
    }

    public MAPMAP1Options() {
        this("SylvesCR", 1000, 0);
    }

    public String getMode() { return mode; }
    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }

    public String component1() { return mode; }
    public int component2() { return maxNumComp; }
    public int component3() { return verbose; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MAPMAP1Options)) return false;
        MAPMAP1Options that = (MAPMAP1Options) o;
        return maxNumComp == that.maxNumComp
                && verbose == that.verbose
                && Objects.equals(mode, that.mode);
    }

    @Override
    public int hashCode() {
        return Objects.hash(mode, maxNumComp, verbose);
    }

    @Override
    public String toString() {
        return "MAPMAP1Options(mode=" + mode + ", maxNumComp=" + maxNumComp + ", verbose=" + verbose + ")";
    }
}
