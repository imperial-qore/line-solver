package jline.lib.smc;

import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

public final class QBD_EG {
    private QBD_EG() {}

    public static Map<String, Matrix> QBD_EG(Matrix A0, Matrix A1, Matrix A2, boolean optVerbose) {
        Matrix G = new Matrix(0, 0, 0);
        Matrix R = new Matrix(0, 0, 0);
        Matrix U = new Matrix(0, 0, 0);
        int m = A1.getNumRows();
        Matrix theta = Stat.stat(A0.add(1.0, A1).add(1.0, A2));
        double drift = theta.mult(A0.sumRows()).get(0) - theta.mult(A2.sumRows()).get(0);
        if (drift > 0) {
            if (A0.rank() == 1) {
                int non_zero_row = 0;
                Matrix row_sum = A0.sumRows();
                for (int i = 0; i < A0.getNumRows(); i++) {
                    if (row_sum.get(i) > 0) {
                        non_zero_row = i;
                        break;
                    }
                }
                Matrix beta = Matrix.scaleMult(
                        Matrix.extractRows(A0, non_zero_row, non_zero_row + 1, null),
                        1 / A0.sumRows(non_zero_row));
                G = Matrix.ones(m, 1).mult(beta);
                R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv());
            } else if (A2.rank() == 1) {
                double eta = QBD_CAUDAL.QBD_CAUDAL(A0, A1, A2);
                R = A2.mult(Matrix.eye(m).add(-1.0, A1).add(-eta, A0).inv());
                G = Matrix.eye(m).add(-1.0, A1.add(1.0, R.mult(A0))).inv().mult(A0);
            }
        } else if (drift < 0) {
            if (A2.rank() == 1) {
                Matrix alpha = A2.mult(Matrix.ones(m, 1));
                R = Matrix.scaleMult(alpha.mult(theta), theta.mult(alpha).get(0));
                G = Matrix.eye(m).add(-1.0, A1.add(1.0, R.mult(A0))).inv().mult(A0);
            } else if (A0.rank() == 0) {
                Matrix A0hat = Matrix.diag(theta.elementPow(-1.0).toArray1D()).mult(A2.transpose())
                        .mult(Matrix.diag(theta.toArray1D()));
                Matrix A1hat = Matrix.diag(theta.elementPow(-1.0).toArray1D()).mult(A1.transpose())
                        .mult(Matrix.diag(theta.toArray1D()));
                Matrix A2hat = Matrix.diag(theta.elementPow(-1.0).toArray1D()).mult(A0.transpose())
                        .mult(Matrix.diag(theta.toArray1D()));
                double etahat = QBD_CAUDAL.QBD_CAUDAL(A0hat, A1hat, A2hat);
                G = Matrix.diag(theta.elementPow(-1.0).toArray1D())
                        .mult(A2hat.mult(Matrix.eye(m).add(-1.0, A1hat).add(-etahat, A0hat).inv()).transpose())
                        .mult(Matrix.diag(theta.toArray1D()));
                R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv());
            }
        }

        if (!R.isEmpty()) {
            U = A1.add(1.0, R.mult(A0));
        }

        Map<String, Matrix> result = new HashMap<String, Matrix>();
        result.put("G", G);
        result.put("R", R);
        result.put("U", U);
        return result;
    }
}
