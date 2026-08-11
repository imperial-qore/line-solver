package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class QBD_ParsePara {
    private QBD_ParsePara() {}

    public static void QBD_ParsePara(Matrix A0, Matrix A1, Matrix A2) {
        if (A0.getNumCols() != A0.getNumRows()) {
            throw new RuntimeException("A0 is not a square matrix");
        }
        if (A1.getNumCols() != A1.getNumRows()) {
            throw new RuntimeException("A1 is not a square matrix");
        }
        if (A2.getNumCols() != A2.getNumRows()) {
            throw new RuntimeException("A2 is not a square matrix");
        }
        if (A0.getNumCols() != A1.getNumCols()) {
            throw new RuntimeException("The matrices A0 and A1 do not have the same dimension");
        }
        if (A0.getNumCols() != A2.getNumCols()) {
            throw new RuntimeException("The matrices A0 and A2 do not have the same dimension");
        }
        if (A0.elementMin() <= -Math.pow(-10.0, -14.0)) {
            throw new RuntimeException("The matrix A0 contains negative data");
        }
        if (A1.elementMin() <= -Math.pow(-10.0, -14.0)) {
            throw new RuntimeException("The matrix A1 contains negative data");
        }
        if (A2.elementMin() <= -Math.pow(-10.0, -14.0)) {
            throw new RuntimeException("he matrix A2 contains negative data");
        }
        if (A0.add(1.0, A1).add(1.0, A2).sumRows().elementMax() > 1 + Math.pow(10.0, -14.0)) {
            throw new RuntimeException("The matrix A0+A1+A2 has to be (sub)stochastic");
        }
    }
}
