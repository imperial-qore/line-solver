package jline.lib.butools.fitting;

import jline.util.matrix.Matrix;

/**
 * Result of MAP fitting.
 */
public class MAPFitResult {
    public final Matrix D0;
    public final Matrix D1;
    public final double logli;

    public MAPFitResult(Matrix D0, Matrix D1, double logli) {
        this.D0 = D0;
        this.D1 = D1;
        this.logli = logli;
    }
}
