/**
 * @file QBD departure process joint moments computation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_pi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_depproc_jointmom {
    private Qbd_depproc_jointmom() {}

    /**
     * Compute joint moments of consecutive inter-departure times.
     */
    public static double[] qbd_depproc_jointmom(MatrixCell MAPa, MatrixCell MAPs, Matrix iset) {
        int na = MAPa.get(0).getNumRows();
        int ns = MAPs.get(0).getNumRows();
        int lvlsz = ns * na;

        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);
        Matrix F = MAPa.get(1).kron(IS);
        Matrix L = MAPa.get(0).kron(IS).add(IA.kron(MAPs.get(0)));
        Matrix B = IA.kron(MAPs.get(1));
        Matrix L0 = MAPa.get(0).kron(IS);

        Map<String, Matrix> qbdResult = QBD_CR.QBD_CR(B, L, F, null, null, null, null);
        Matrix R = qbdResult.get("R");
        if (R == null) {
            throw new RuntimeException("QBD_CR failed to compute R matrix");
        }

        Matrix pi = QBD_pi.QBD_pi(B, L0, R, 100, 0, null, 0);

        Matrix v0 = new Matrix(1, lvlsz);
        for (int j = 0; j < lvlsz; j++) {
            v0.set(0, j, pi.get(0, j));
        }

        double lambdaS = Map_lambda.map_lambda(MAPs);
        double invLambdaS = 1.0 / lambdaS;
        // departure epochs are the B transitions, so the embedded vector weighs the level
        // probabilities by B and not by the arrival matrix F
        Matrix v0R = v0.mult(R);
        Matrix v0D = v0R.mult(B).scale(invLambdaS);

        Matrix v0R2 = v0R.mult(R);
        Matrix v1D = v0R2.mult(B).scale(invLambdaS);

        Matrix v0R3 = v0R2.mult(R);
        Matrix ImR = Matrix.eye(R.getNumRows()).add(-1.0, R);
        Matrix ImRinv = ImR.inv();
        Matrix v2Dp = v0R3.mult(ImRinv).mult(B).scale(invLambdaS);

        Matrix z = new Matrix(1, 3 * lvlsz);
        for (int j = 0; j < lvlsz; j++) {
            z.set(0, j, v0D.get(0, j));
            z.set(0, lvlsz + j, v1D.get(0, j));
            z.set(0, 2 * lvlsz + j, v2Dp.get(0, j));
        }

        double zSum = z.elementSum();
        z.scaleEq(1.0 / zSum);

        int dim3 = 3 * lvlsz;
        Matrix M0 = Matrix.zeros(dim3, dim3);
        M0.insertSubMatrix(0, 0, lvlsz, lvlsz, L0);
        M0.insertSubMatrix(0, lvlsz, lvlsz, 2 * lvlsz, F);
        M0.insertSubMatrix(lvlsz, lvlsz, 2 * lvlsz, 2 * lvlsz, L);
        M0.insertSubMatrix(lvlsz, 2 * lvlsz, 2 * lvlsz, 3 * lvlsz, F);
        M0.insertSubMatrix(2 * lvlsz, 2 * lvlsz, 3 * lvlsz, 3 * lvlsz, L.add(F));

        Matrix M1 = Matrix.zeros(dim3, dim3);
        M1.insertSubMatrix(lvlsz, 0, 2 * lvlsz, lvlsz, B);
        M1.insertSubMatrix(2 * lvlsz, lvlsz, 3 * lvlsz, 2 * lvlsz, B);

        Matrix onesVec = Matrix.ones(dim3, 1);

        Matrix negM0 = M0.neg();
        Matrix negM0inv = negM0.inv();

        int numMoments = iset.getNumRows();
        double[] JM = new double[numMoments];

        for (int k = 0; k < numMoments; k++) {
            int iOrder = (int) iset.get(k, 0);
            int jOrder = (int) iset.get(k, 1);

            Matrix negM0invPowI1 = Matrix.eye(dim3);
            for (int p = 0; p < iOrder + 1; p++) {
                negM0invPowI1 = negM0invPowI1.mult(negM0inv);
            }

            Matrix negM0invPowJ = Matrix.eye(dim3);
            for (int p = 0; p < jOrder; p++) {
                negM0invPowJ = negM0invPowJ.mult(negM0inv);
            }

            double factI = factorial(iOrder);
            double factJ = factorial(jOrder);

            Matrix temp1 = z.mult(negM0invPowI1).scale(factI);
            Matrix temp2 = temp1.mult(M1);
            Matrix temp3 = temp2.mult(negM0invPowJ).scale(factJ);
            Matrix result = temp3.mult(onesVec);

            JM[k] = result.get(0, 0);
        }

        return JM;
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= (double) i;
        }
        return result;
    }
}
