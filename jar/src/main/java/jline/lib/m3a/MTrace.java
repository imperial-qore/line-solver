/**
 * @file MTrace data class for multiclass trace representation
 *
 * @since LINE 3.0
 */
package jline.lib.m3a;

import java.util.Arrays;

/**
 * Data structure for multiclass trace representation.
 */
public final class MTrace {
    private final double[] S;
    private final int[] C;
    private final int numClasses;

    public MTrace(double[] S, int[] C, int numClasses) {
        this.S = S;
        this.C = C;
        this.numClasses = numClasses;
    }

    public double[] getS() { return S; }
    public int[] getC() { return C; }
    public int getNumClasses() { return numClasses; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MTrace)) return false;
        MTrace that = (MTrace) o;
        return Arrays.equals(S, that.S) && Arrays.equals(C, that.C) && numClasses == that.numClasses;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(S);
        result = 31 * result + Arrays.hashCode(C);
        result = 31 * result + numClasses;
        return result;
    }

    @Override
    public String toString() {
        return "MTrace(S=" + Arrays.toString(S) + ", C=" + Arrays.toString(C)
                + ", numClasses=" + numClasses + ")";
    }
}
