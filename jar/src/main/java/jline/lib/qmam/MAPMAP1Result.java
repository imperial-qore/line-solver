/**
 * @file Result of MAP/MAP/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of MAP/MAP/1 queue analysis.
 */
public final class MAPMAP1Result {
    public final Matrix queueLength; // Queue length distribution
    public final Matrix sojAlpha;    // Sojourn time PH alpha vector
    public final Matrix waitAlpha;   // Waiting time PH alpha vector
    public final Matrix Smat;        // Service time matrix

    public MAPMAP1Result(Matrix queueLength, Matrix sojAlpha, Matrix waitAlpha, Matrix Smat) {
        this.queueLength = queueLength;
        this.sojAlpha = sojAlpha;
        this.waitAlpha = waitAlpha;
        this.Smat = Smat;
    }

    public Matrix getQueueLength() { return queueLength; }
    public Matrix getSojAlpha() { return sojAlpha; }
    public Matrix getWaitAlpha() { return waitAlpha; }
    public Matrix getSmat() { return Smat; }

    public Matrix component1() { return queueLength; }
    public Matrix component2() { return sojAlpha; }
    public Matrix component3() { return waitAlpha; }
    public Matrix component4() { return Smat; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MAPMAP1Result)) return false;
        MAPMAP1Result that = (MAPMAP1Result) o;
        return Objects.equals(queueLength, that.queueLength)
                && Objects.equals(sojAlpha, that.sojAlpha)
                && Objects.equals(waitAlpha, that.waitAlpha)
                && Objects.equals(Smat, that.Smat);
    }

    @Override
    public int hashCode() {
        return Objects.hash(queueLength, sojAlpha, waitAlpha, Smat);
    }

    @Override
    public String toString() {
        return "MAPMAP1Result(queueLength=" + queueLength + ", sojAlpha=" + sojAlpha
                + ", waitAlpha=" + waitAlpha + ", Smat=" + Smat + ")";
    }
}
