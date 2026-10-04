/**
 * @file Quasi-Birth-Death process MAP/MAP/1 queue analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_pi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_mapmap1 {
    private Qbd_mapmap1() {}

    public static QbdMapMap1Result qbd_mapmap1(MatrixCell MAPa, MatrixCell MAPs) {
        return qbd_mapmap1(MAPa, MAPs, null);
    }

    /**
     * Analyze MAP/MAP/1 queue using QBD methods.
     */
    public static QbdMapMap1Result qbd_mapmap1(MatrixCell MAPa, MatrixCell MAPs, Double util) {
        int na = MAPa.get(0).getNumRows();
        int ns = MAPs.get(0).getNumRows();

        MatrixCell scaledMAPs = MAPs;
        if (util != null) {
            double lambdaA0 = Map_lambda.map_lambda(MAPa);
            scaledMAPs = Map_scale.map_scale(MAPs, util / lambdaA0);
        }

        double lambdaA = Map_lambda.map_lambda(MAPa);
        double lambdaS = Map_lambda.map_lambda(scaledMAPs);
        double actualUtil = lambdaA / lambdaS;

        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);
        Matrix A1 = MAPa.get(1).kron(IS);
        Matrix A0 = MAPa.get(0).kron(IS).add(IA.kron(scaledMAPs.get(0)));
        Matrix A_1 = IA.kron(scaledMAPs.get(1));
        Matrix A0bar = MAPa.get(0).kron(IS);

        Map<String, Matrix> qbdResult = QBD_CR.QBD_CR(A_1, A0, A1, null, null, null, null);
        Matrix G = qbdResult.get("G");
        Matrix R = qbdResult.get("R");
        Matrix U = qbdResult.get("U");
        if (G == null || R == null || U == null) {
            throw new RuntimeException("QBD_CR failed");
        }

        int n = na * ns;

        Matrix pi = QBD_pi.QBD_pi(A_1, A0bar, R, 100, 0, null, 0);
        int numLevels = pi.getNumCols() / n;

        Matrix pqueue = new Matrix(numLevels, n);
        for (int i = 0; i < numLevels; i++) {
            for (int j = 0; j < n; j++) {
                pqueue.set(i, j, pi.get(0, i * n + j));
            }
        }

        double busyProb = 0.0;
        for (int i = 1; i < numLevels; i++) {
            for (int j = 0; j < n; j++) {
                busyProb += pqueue.get(i, j);
            }
        }
        if (busyProb < actualUtil * 0.99) {
            Matrix pi2 = QBD_pi.QBD_pi(A_1, A0bar, R, 20000, 0, null, 0);
            int numLevels2 = pi2.getNumCols() / n;
            pqueue = new Matrix(numLevels2, n);
            for (int i = 0; i < numLevels2; i++) {
                for (int j = 0; j < n; j++) {
                    pqueue.set(i, j, pi2.get(0, i * n + j));
                }
            }
        }

        double UN;
        double QN;

        if (na == 1 && ns == 1) {
            UN = 1.0 - pqueue.get(0, 0);
            QN = 0.0;
            for (int i = 0; i < pqueue.getNumRows(); i++) {
                QN += (double) i * pqueue.get(i, 0);
            }
        } else {
            double sum0 = 0.0;
            for (int j = 0; j < n; j++) {
                sum0 += pqueue.get(0, j);
            }
            UN = 1.0 - sum0;
            QN = 0.0;
            for (int i = 0; i < pqueue.getNumRows(); i++) {
                double levelProb = 0.0;
                for (int j = 0; j < n; j++) {
                    levelProb += pqueue.get(i, j);
                }
                QN += (double) i * levelProb;
            }
        }

        double XN = lambdaA;

        Matrix eta = new Matrix(1, 1);
        eta.set(0, 0, spectralRadiusMapmap1(R));

        MatrixCell resultMAPs = new MatrixCell(2);
        resultMAPs.set(0, scaledMAPs.get(0));
        resultMAPs.set(1, scaledMAPs.get(1));

        return new QbdMapMap1Result(XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, resultMAPs);
    }

    private static double spectralRadiusMapmap1(Matrix A) {
        int n = A.getNumRows();
        org.ejml.data.DMatrixRMaj dm = new org.ejml.data.DMatrixRMaj(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dm.set(i, j, A.get(i, j));
            }
        }
        org.ejml.interfaces.decomposition.EigenDecomposition_F64<org.ejml.data.DMatrixRMaj> evd =
                org.ejml.dense.row.factory.DecompositionFactory_DDRM.eig(n, false);
        evd.decompose(dm);
        double maxAbs = 0.0;
        for (int i = 0; i < evd.getNumberOfEigenvalues(); i++) {
            org.ejml.data.Complex_F64 ev = evd.getEigenvalue(i);
            double absVal = Math.sqrt(ev.real * ev.real + ev.imaginary * ev.imaginary);
            if (absVal > maxAbs) {
                maxAbs = absVal;
            }
        }
        return maxAbs;
    }
}
