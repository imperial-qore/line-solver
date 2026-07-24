package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * Result of returnWait
 */
public final class ReturnWaitResult {
    public final Matrix wait_alpha;
    public final Matrix wait_Smat;
    public final double prob_wait;
    public final Matrix alfa;

    public ReturnWaitResult(Matrix wait_alpha, Matrix wait_Smat, double prob_wait, Matrix alfa) {
        this.wait_alpha = wait_alpha;
        this.wait_Smat = wait_Smat;
        this.prob_wait = prob_wait;
        this.alfa = alfa;
    }

    public Matrix getWait_alpha() { return wait_alpha; }
    public Matrix getWait_Smat() { return wait_Smat; }
    public double getProb_wait() { return prob_wait; }
    public Matrix getAlfa() { return alfa; }
}
