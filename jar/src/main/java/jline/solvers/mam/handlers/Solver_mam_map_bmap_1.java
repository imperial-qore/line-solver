package jline.solvers.mam.handlers;

import java.util.List;

import jline.api.mam.Map_prob;
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

        Matrix R = gim1_r_etaqa(A);
        Matrix pi = gim1_pi_etaqa(B, A, R, A0);
        double meanQueueLength = gim1_qlen_etaqa(B, A, R, pi, 1);
        double meanResponseTime = (lambdaArr > 0) ? meanQueueLength / lambdaArr : 0.0;
        return new MAPBMAP1Result(meanQueueLength, rho, meanResponseTime, lambdaArr, pi, R, meanBatchSize);
    }

    private static Matrix gim1_r_etaqa(Matrix A) {
        int m = A.getNumCols();
        int numBlocks = A.getNumRows() / m;
        Matrix[] blocks = new Matrix[numBlocks];
        for (int i = 0; i < numBlocks; i++) blocks[i] = A.extractRows(i * m, (i + 1) * m);
        Matrix A1 = blocks[1];
        double minDiag = Double.MAX_VALUE;
        for (int i = 0; i < m; i++) if (A1.get(i, i) < minDiag) minDiag = A1.get(i, i);
        Matrix[] uniformized;
        if (minDiag < 0) {
            double lamb = -minDiag;
            uniformized = new Matrix[numBlocks];
            for (int i = 0; i < numBlocks; i++) {
                Matrix scaled = blocks[i].scale(1.0 / lamb);
                if (i == 1) {
                    for (int j = 0; j < m; j++) scaled.set(j, j, scaled.get(j, j) + 1.0);
                }
                uniformized[i] = scaled;
            }
        } else {
            uniformized = blocks;
        }
        Matrix R = new Matrix(m, m);
        int maxIter = 200;
        double tol = 1e-14;
        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Rnew = uniformized[0].copy();
            Matrix Rpow = R.copy();
            for (int k = 1; k < numBlocks; k++) {
                Rnew = Rnew.add(uniformized[k].mult(Rpow));
                Rpow = Rpow.mult(R);
            }
            double diff = Rnew.sub(R).normFrobenius();
            if (diff < tol) return Rnew;
            R = Rnew;
        }
        return R;
    }

    private static Matrix gim1_pi_etaqa(Matrix B, Matrix A, Matrix R, Matrix B0) {
        int m = R.getNumRows();
        Matrix Iminr = Matrix.eye(m).sub(R);
        Matrix IminrInv = Iminr.inv();
        Matrix B1 = B.extractRows(0, m);
        Matrix temp = B1.add(B0.mult(R).mult(IminrInv));
        Matrix pi0 = stat(temp);
        double norm = pi0.mult(IminrInv).mult(Matrix.ones(m, 1)).get(0, 0);
        Matrix pi0Norm = (norm > 0) ? pi0.scale(1.0 / norm) : pi0;
        Matrix pi1 = pi0Norm.mult(R);
        Matrix piStar = pi1.mult(R);
        Matrix result = new Matrix(1, 3 * m);
        for (int i = 0; i < m; i++) {
            result.set(0, i, pi0Norm.get(0, i));
            result.set(0, m + i, pi1.get(0, i));
            result.set(0, 2 * m + i, piStar.get(0, i));
        }
        return result;
    }

    private static double gim1_qlen_etaqa(Matrix B, Matrix A, Matrix R, Matrix pi, int n) {
        int m = R.getNumRows();
        Matrix pi1 = pi.extractCols(m, 2 * m);
        Matrix piStar = pi.extractCols(2 * m, 3 * m);
        Matrix e = Matrix.ones(m, 1);
        if (n == 1) return pi1.mult(e).get(0, 0) + 2.0 * piStar.mult(e).get(0, 0);
        double moment = pi1.mult(e).get(0, 0);
        double piStarSum = piStar.mult(e).get(0, 0);
        for (int k = 2; k < 100; k++) {
            double levelProb = piStarSum * Math.pow(1.0 - piStarSum, (double) (k - 2));
            moment += Math.pow((double) k, (double) n) * levelProb;
            if (levelProb < 1e-15) break;
        }
        return moment;
    }

    private static Matrix stat(Matrix A) {
        int n = A.getNumRows();
        Matrix B = new Matrix(n, n + 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) B.set(i, j, A.get(j, i));
            B.set(i, n, 1.0);
        }
        Matrix y = new Matrix(1, n + 1);
        y.set(0, n, 1.0);
        Matrix pi = y.mult(B.pinv());
        for (int i = 0; i < n; i++) {
            if (pi.get(0, i) < 0) pi.set(0, i, 0.0);
        }
        double sum = pi.elementSum();
        if (sum > 0) {
            for (int i = 0; i < n; i++) pi.set(0, i, pi.get(0, i) / sum);
        }
        return pi;
    }
}
