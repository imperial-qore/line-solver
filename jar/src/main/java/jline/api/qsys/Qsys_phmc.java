/**
 * @file PH/M/c queueing system analysis (matrix-geometric)
 *
 * Solves the GI/M/c QBD via Neuts' matrix-geometric method:
 *   R^2 A2 + R A1 + A0 = 0,  A0 = (-T 1) alpha,  A1 = T - c mu I,  A2 = c mu I
 * Stationary level vectors satisfy pi_n = pi_c R^{n-c} for n >= c. Boundary
 * states pi_0..pi_c are determined by linear balance + normalization.
 */
package jline.api.qsys;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public final class Qsys_phmc {
    private Qsys_phmc() {}

    public static PhMcResult qsys_phmc(double[] alpha, double[][] T, double mu, int c, int maxIter, double tol) {
        if (!(mu > 0.0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(c >= 1)) throw new IllegalArgumentException("Number of servers c must be >= 1");
        int k = T.length;
        if (!(k > 0)) throw new IllegalArgumentException("T must be non-empty");
        for (int i = 0; i < k; i++) {
            if (T[i].length != k) throw new IllegalArgumentException("T must be square");
        }
        if (alpha.length != k) throw new IllegalArgumentException("alpha length must match T dimension");

        Array2DRowRealMatrix Tm = new Array2DRowRealMatrix(T, false);
        Array2DRowRealMatrix negTm = (Array2DRowRealMatrix) Tm.scalarMultiply(-1.0);
        double[] ones = new double[k];
        for (int i = 0; i < k; i++) ones[i] = 1.0;
        double[] tVec = negTm.operate(ones);
        double[][] D1 = new double[k][k];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                D1[i][j] = tVec[i] * alpha[j];
            }
        }
        Array2DRowRealMatrix D1m = new Array2DRowRealMatrix(D1, false);

        LUDecomposition negTmLU = new LUDecomposition(negTm);
        if (!negTmLU.getSolver().isNonSingular()) throw new IllegalArgumentException("PH sub-generator T is singular");
        double[] negTmInvOnes = negTmLU.getSolver().solve(new ArrayRealVector(ones)).toArray();
        double meanIa = 0.0;
        for (int i = 0; i < k; i++) meanIa += alpha[i] * negTmInvOnes[i];
        if (meanIa <= 0.0) throw new IllegalArgumentException("Non-positive mean inter-arrival: " + meanIa);
        double lambda = 1.0 / meanIa;
        double rho = lambda / (c * mu);
        if (rho >= 1.0 - 1e-12) throw new IllegalArgumentException("Load rho=" + rho + " must be strictly less than 1");

        RealMatrix Ik = MatrixUtils.createRealIdentityMatrix(k);
        Array2DRowRealMatrix A0 = D1m;
        RealMatrix A1 = Tm.subtract(Ik.scalarMultiply(c * mu));
        RealMatrix A2 = Ik.scalarMultiply(c * mu);

        Array2DRowRealMatrix R = new Array2DRowRealMatrix(k, k);  // zeros
        Array2DRowRealMatrix negA0 = (Array2DRowRealMatrix) A0.scalarMultiply(-1.0);
        for (int it = 0; it < maxIter; it++) {
            Array2DRowRealMatrix M = (Array2DRowRealMatrix) A1.add(R.multiply(A2));
            LUDecomposition lu = new LUDecomposition(M);
            if (!lu.getSolver().isNonSingular()) break;
            RealMatrix Rinv = lu.getSolver().getInverse();
            Array2DRowRealMatrix Rnew = (Array2DRowRealMatrix) negA0.multiply(Rinv);
            double maxDiff = 0.0;
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    maxDiff = Math.max(maxDiff, Math.abs(Rnew.getEntry(i, j) - R.getEntry(i, j)));
                }
            }
            R = Rnew;
            if (maxDiff < tol) break;
        }

        // Boundary linear system for pi_0, ..., pi_c (each row vec dim k).
        // Stack as a single vector of length (c+1)*k and build (c+1)*k x (c+1)*k system M_sys.
        int nVar = (c + 1) * k;
        double[][] Msys = new double[nVar][nVar];

        // Level 0: pi_0 T + pi_1 (mu I) = 0  → for each j ∈ [0, k), eq row j
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < k; i++) Msys[j][blockCol(0, k) + i] += T[i][j];
            Msys[j][blockCol(1, k) + j] += mu;
        }
        // Levels 1..c-1
        for (int n = 1; n < c; n++) {
            int eqBase = n * k;
            // A1n = T - n mu I
            for (int j = 0; j < k; j++) {
                int row = eqBase + j;
                for (int i = 0; i < k; i++) {
                    Msys[row][blockCol(n - 1, k) + i] += D1[i][j];
                    Msys[row][blockCol(n, k) + i] += T[i][j];
                    if (i == j) Msys[row][blockCol(n, k) + i] += -n * mu;
                }
                Msys[row][blockCol(n + 1, k) + j] += (n + 1) * mu;
            }
        }
        // Level c: pi_{c-1} D1 + pi_c (T - c mu I + R (c mu I)) = 0
        RealMatrix A1c_plus_RA2 = A1.add(R.scalarMultiply(c * mu));
        int eqBaseC = c * k;
        for (int j = 0; j < k; j++) {
            int row = eqBaseC + j;
            for (int i = 0; i < k; i++) {
                Msys[row][blockCol(c - 1, k) + i] += D1[i][j];
                Msys[row][blockCol(c, k) + i] += A1c_plus_RA2.getEntry(i, j);
            }
        }
        // Replace last row with normalization: sum_{n<c} pi_n*1 + pi_c*(I-R)^{-1}*1 = 1
        RealMatrix IR = Ik.subtract(R);
        LUDecomposition luIR = new LUDecomposition(IR);
        if (!luIR.getSolver().isNonSingular()) throw new IllegalArgumentException("I-R is singular (rho=" + rho + ")");
        double[] sumGeom = luIR.getSolver().solve(new ArrayRealVector(ones)).toArray();
        double[] normRow = new double[nVar];
        for (int n = 0; n < c; n++) {
            for (int i = 0; i < k; i++) normRow[blockCol(n, k) + i] = 1.0;
        }
        for (int i = 0; i < k; i++) normRow[blockCol(c, k) + i] = sumGeom[i];
        Msys[nVar - 1] = normRow;
        double[] rhs = new double[nVar];
        rhs[nVar - 1] = 1.0;

        double[] piVec = new LUDecomposition(new Array2DRowRealMatrix(Msys, false)).getSolver()
                .solve(new ArrayRealVector(rhs)).toArray();
        double[][] pis = new double[c + 1][k];
        for (int n = 0; n <= c; n++) {
            for (int i = 0; i < k; i++) pis[n][i] = piVec[blockCol(n, k) + i];
        }
        double[] piC = pis[c];

        RealMatrix IRinv = luIR.getSolver().getInverse();
        Array2DRowRealMatrix IRinv2 = (Array2DRowRealMatrix) IRinv.multiply(IRinv);

        // Lq = pi_c R IRinv2 ones
        Array2DRowRealMatrix tmp1 = (Array2DRowRealMatrix) R.multiply(IRinv2);
        double[] Lq_vec = tmp1.operate(ones);
        double Lq = 0.0;
        for (int i = 0; i < k; i++) Lq += piC[i] * Lq_vec[i];

        // L_bulk = pi_c (c IRinv + R IRinv2) ones
        Array2DRowRealMatrix cIRinv = (Array2DRowRealMatrix) IRinv.scalarMultiply((double) c);
        Array2DRowRealMatrix LbulkMat = (Array2DRowRealMatrix) cIRinv.add(tmp1);
        double[] Lbulk_vec = LbulkMat.operate(ones);
        double Lbulk = 0.0;
        for (int i = 0; i < k; i++) Lbulk += piC[i] * Lbulk_vec[i];

        double L = Lbulk;
        for (int n = 0; n < c; n++) {
            double s = 0.0;
            for (int i = 0; i < k; i++) s += pis[n][i];
            L += n * s;
        }

        double Wq = Lq / lambda;
        double W = Wq + 1.0 / mu;
        return new PhMcResult(L, Lq, Wq, W, rho);
    }

    public static PhMcResult qsys_phmc(double[] alpha, double[][] T, double mu, int c, int maxIter) {
        return qsys_phmc(alpha, T, mu, c, maxIter, 1e-14);
    }

    public static PhMcResult qsys_phmc(double[] alpha, double[][] T, double mu, int c) {
        return qsys_phmc(alpha, T, mu, c, 50000, 1e-14);
    }

    private static int blockCol(int n, int k) {
        return n * k;
    }
}
