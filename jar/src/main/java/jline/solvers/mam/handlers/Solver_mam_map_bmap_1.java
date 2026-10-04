package jline.solvers.mam.handlers;

import java.util.List;

import jline.api.mam.Map_prob;
import jline.lib.smc.GIM1_ETAQA;
import jline.util.matrix.Matrix;

/**
 * MAP/BMAP/1 Queue Solver using GI/M/1 type analysis with ETAQA.
 */
public final class Solver_mam_map_bmap_1 {
    private Solver_mam_map_bmap_1() {}

    public static MAPBMAP1Result solver_mam_map_bmap_1(Matrix C0, Matrix C1, List<Matrix> D) {
        if (D.isEmpty()) throw new IllegalArgumentException("BMAP must have at least D0 matrix");
        if (D.size() < 2) throw new IllegalArgumentException("BMAP must have at least D0 and D1 matrices");

        int K = D.size() - 1;
        int ma = C0.getNumRows();
        int ms = D.get(0).getNumRows();
        int m = ma * ms;
        if (C0.getNumRows() != ma || C0.getNumCols() != ma) throw new IllegalArgumentException("C0 must be " + ma + "x" + ma);
        if (C1.getNumRows() != ma || C1.getNumCols() != ma) throw new IllegalArgumentException("C1 must be " + ma + "x" + ma);
        for (int i = 0; i < D.size(); i++) {
            if (D.get(i).getNumRows() != ms || D.get(i).getNumCols() != ms) {
                throw new IllegalArgumentException("All BMAP matrices must be " + ms + "x" + ms);
            }
        }

        Matrix piC = Map_prob.map_prob(C0, C1);
        Matrix eA = Matrix.ones(ma, 1);
        double lambdaArr = piC.mult(C1).mult(eA).get(0, 0);

        Matrix D1Total = new Matrix(ms, ms);
        for (int k = 1; k <= K; k++) D1Total = D1Total.add(D.get(k));
        Matrix piD = Map_prob.map_prob(D.get(0), D1Total);
        Matrix eS = Matrix.ones(ms, 1);

        double muTotal = 0.0;
        for (int k = 1; k <= K; k++) {
            double rateK = piD.mult(D.get(k)).mult(eS).get(0, 0);
            muTotal += k * rateK;
        }
        double batchRate = piD.mult(D1Total).mult(eS).get(0, 0);
        double meanBatchSize = (batchRate > 0) ? muTotal / batchRate : 0.0;
        double rho = lambdaArr / muTotal;
        if (rho >= 1.0) System.err.println("Warning: System is unstable (rho = " + rho + " >= 1). Results may be invalid.");

        Matrix A = new Matrix(m * (K + 2), m);
        Matrix Ims = Matrix.eye(ms);
        Matrix A0 = C1.kron(Ims);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) A.set(i, j, A0.get(i, j));
        }
        Matrix Ima = Matrix.eye(ma);
        Matrix A1 = C0.kron(Ims).add(Ima.kron(D.get(0)));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) A.set(m + i, j, A1.get(i, j));
        }
        for (int k = 1; k <= K; k++) {
            Matrix Ak = Ima.kron(D.get(k));
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) A.set((k + 1) * m + i, j, Ak.get(i, j));
            }
        }

        Matrix B = new Matrix(m * (K + 2), m);
        Matrix B1 = C0.kron(Ims).add(Ima.kron(D.get(0)));
        for (int k = 1; k <= K; k++) {
            Matrix Dk = Ima.kron(D.get(k));
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) B1.set(i, j, B1.get(i, j) + Dk.get(i, j));
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) B.set(i, j, B1.get(i, j));
        }
        for (int j = 1; j <= K; j++) {
            Matrix Bj = new Matrix(m, m);
            for (int k = j; k <= K; k++) {
                Matrix Dk = Ima.kron(D.get(k));
                for (int i = 0; i < m; i++) {
                    for (int jj = 0; jj < m; jj++) Bj.set(i, jj, Bj.get(i, jj) + Dk.get(i, jj));
                }
            }
            for (int i = 0; i < m; i++) {
                for (int jj = 0; jj < m; jj++) B.set(m + (j - 1) * m + i, jj, Bj.get(i, jj));
            }
        }

        // MAMSolver's GIM1_R_ETAQA / GIM1_pi_ETAQA / GIM1_qlen_ETAQA. The
        // reference asks for the FIRST moment only on this side, and its
        // GIM1_qlen_ETAQA carries the scalar-A(3) defect that can make the
        // reported mean negative for more than one phase; see GIM1_ETAQA.
        Matrix R = GIM1_ETAQA.gim1_r_etaqa(A);
        Matrix pi = GIM1_ETAQA.gim1_pi_etaqa(B, A, R, A0);
        double meanQueueLength = GIM1_ETAQA.gim1_qlen_etaqa(B, A, R, pi, 1, A0);
        double meanResponseTime = (lambdaArr > 0) ? meanQueueLength / lambdaArr : 0.0;
        return new MAPBMAP1Result(meanQueueLength, rho, meanResponseTime, lambdaArr, pi, R, meanBatchSize);
    }

}
