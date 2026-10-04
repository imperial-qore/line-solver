package jline.solvers.ag.handlers;

import jline.util.matrix.Matrix;

/**
 * The matrix-geometric tail of one open RCAT component with more than one phase
 * per level: pi_(n+1) = pi_n R, closed by the boundary equations of levels 0
 * and 1. BUSY is sum_(n>=1) pi_n = pi_1 (I - R)^-1, the row vector every tail
 * moment is read off, and QLEN is E[N] = pi_1 (I - R)^-2 e.
 */
public final class QbdTail {
    public final Matrix R;
    public final Matrix pi0;
    public final Matrix pi1;
    public final Matrix busy;
    public final double qlen;

    public QbdTail(Matrix R, Matrix pi0, Matrix pi1, Matrix busy, double qlen) {
        this.R = R;
        this.pi0 = pi0;
        this.pi1 = pi1;
        this.busy = busy;
        this.qlen = qlen;
    }
}
