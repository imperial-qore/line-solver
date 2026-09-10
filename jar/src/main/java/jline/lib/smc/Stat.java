package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class Stat {
    private Stat() {}

    public static Matrix stat(Matrix A) {
        int S = A.getNumRows();
        Matrix e = Matrix.ones(S, 1);
        Matrix B = Matrix.concatColumns(A.add(-1.0, Matrix.eye(S)), e, null);
        Matrix y = Matrix.concatColumns(new Matrix(1, S), Matrix.singleton(1.0), null);
        // MATLAB: theta = y / B solves X * B = y for X (row vector)
        return y.rightMatrixDivide(B);
    }
}
