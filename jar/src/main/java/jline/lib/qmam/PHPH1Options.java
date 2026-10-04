/**
 * @file Options for PH/PH/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Objects;

/**
 * Options for PH/PH/1 queue analysis.
 */
public final class PHPH1Options {
    public final int maxNumComp;
    public final int verbose;

    public PHPH1Options(int maxNumComp, int verbose) {
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
    }

    public PHPH1Options() {
        this(1000, 0);
    }

    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }

    public int component1() { return maxNumComp; }
    public int component2() { return verbose; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PHPH1Options)) return false;
        PHPH1Options that = (PHPH1Options) o;
        return maxNumComp == that.maxNumComp && verbose == that.verbose;
    }

    @Override
    public int hashCode() {
        return Objects.hash(maxNumComp, verbose);
    }

    @Override
    public String toString() {
        return "PHPH1Options(maxNumComp=" + maxNumComp + ", verbose=" + verbose + ")";
    }
}
