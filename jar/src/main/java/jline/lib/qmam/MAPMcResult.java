package jline.lib.qmam;

import jline.util.matrix.Matrix;

/**
 * Result of MAP/M/c queue analysis.
 */
public final class MAPMcResult {
    public final Matrix queueLength;
    public final Matrix waitAlpha; // nullable
    public final Matrix Smat;      // nullable

    public MAPMcResult(Matrix queueLength, Matrix waitAlpha, Matrix Smat) {
        this.queueLength = queueLength;
        this.waitAlpha = waitAlpha;
        this.Smat = Smat;
    }

    public Matrix getQueueLength() { return queueLength; }
    public Matrix getWaitAlpha() { return waitAlpha; }
    public Matrix getSmat() { return Smat; }
}
