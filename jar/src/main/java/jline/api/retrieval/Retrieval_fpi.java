/**
 * @file Retrieval_fpi.java
 * @brief Fixed-point (FPI) approximation of delayed-hit cache metrics.
 *
 * Port of matlab/src/api/retrieval/retrieval_fpi.m. Iteration t -> t+1:
 *   F_{s,i}    = 1 + sum_{k!=i} phi_{s,k}
 *   D_i        = 1 + lambda_i*eta_{0,i} + sum_s lambda_i*eta_{s,i}*F_{s,i}
 *   theta_{ij} = gamma_{ij}/D_i
 *   xi_j       = m_j / sum_k theta_{kj}(1 - sum_l pi_{kl})
 *   pi_{ij}    = theta_{ij} xi_j / (1 + sum_l theta_{il} xi_l)
 *   pi_{i0}    = (1 - sum_j pi_{ij})/D_i
 *   phi_{0i}   = lambda_i*eta_{0,i} pi_{i0};  phi_{si} = lambda_i*eta_{s,i} F_{s,i} pi_{i0}
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

public final class Retrieval_fpi {
    private Retrieval_fpi() {}

    /** Result: pmiss (1 x n), phit (h x n), pdh ((r+1) x n). */
    public static final class Result {
        public final double[] pmiss;
        public final double[][] phit;
        public final double[][] pdh;
        Result(double[] pmiss, double[][] phit, double[][] pdh) {
            this.pmiss = pmiss; this.phit = phit; this.pdh = pdh;
        }
    }

    public static Result retrieval_fpi(double[] m, double[] lambda, double[][] eta, double[][] gamma) {
        return retrieval_fpi(m, lambda, eta, gamma, 1000, 1e-6);
    }

    public static Result retrieval_fpi(double[] m, double[] lambda, double[][] eta, double[][] gamma,
                                       int maxIter, double tol) {
        int n = lambda.length;
        int h = m.length;
        int r = eta[0].length - 1;

        double[] eta0 = new double[n];
        double[][] etaPS = new double[n][r];
        for (int i = 0; i < n; i++) {
            eta0[i] = eta[i][0];
            for (int s = 0; s < r; s++) etaPS[i][s] = eta[i][s + 1];
        }

        double[][] phi = new double[r + 1][n];
        double[][] pij = new double[h][n];
        double[] pi0 = new double[n];
        double initPhi = 1.0 / ((h + 1.0) * (r + 2.0));
        double initPij = 1.0 / (h + 1.0);
        for (double[] row : phi) java.util.Arrays.fill(row, initPhi);
        for (double[] row : pij) java.util.Arrays.fill(row, initPij);
        java.util.Arrays.fill(pi0, initPhi);

        for (int t = 0; t < maxIter; t++) {
            // F(s,i) = 1 + sum_{k!=i} phi_{s,k}
            double[][] F = new double[r][n];
            for (int s = 0; s < r; s++) {
                double sumPhis = 0;
                for (int i = 0; i < n; i++) sumPhis += phi[s + 1][i];
                for (int i = 0; i < n; i++) F[s][i] = 1 + (sumPhis - phi[s + 1][i]);
            }
            // D_i
            double[] D = new double[n];
            for (int i = 0; i < n; i++) {
                D[i] = 1 + lambda[i] * eta0[i];
                for (int s = 0; s < r; s++) D[i] += lambda[i] * etaPS[i][s] * F[s][i];
            }
            // theta_{ij} = gamma_{ij}/D_i
            double[][] theta = new double[n][h];
            for (int i = 0; i < n; i++) for (int j = 0; j < h; j++) theta[i][j] = gamma[i][j] / D[i];
            // oneminus_k = 1 - sum_l pij_{l,k}
            double[] oneminus = new double[n];
            for (int k = 0; k < n; k++) {
                double sp = 0; for (int j = 0; j < h; j++) sp += pij[j][k];
                oneminus[k] = 1 - sp;
            }
            // xi_j = m_j / sum_k theta_{kj} oneminus_k
            double[] xi = new double[h];
            for (int j = 0; j < h; j++) {
                double den = 0; for (int k = 0; k < n; k++) den += theta[k][j] * oneminus[k];
                xi[j] = m[j] / den;
            }
            // pij_new_{j,i} = theta_ij xi_j / (1 + sum_l theta_il xi_l)
            double[][] pijNew = new double[h][n];
            double[] pi0New = new double[n];
            double[][] phiNew = new double[r + 1][n];
            for (int i = 0; i < n; i++) {
                double denom = 1;
                for (int l = 0; l < h; l++) denom += theta[i][l] * xi[l];
                double sumPij = 0;
                for (int j = 0; j < h; j++) {
                    pijNew[j][i] = theta[i][j] * xi[j] / denom;
                    sumPij += pijNew[j][i];
                }
                pi0New[i] = (1 - sumPij) / D[i];
                phiNew[0][i] = lambda[i] * eta0[i] * pi0New[i];
                for (int s = 0; s < r; s++) phiNew[s + 1][i] = lambda[i] * etaPS[i][s] * F[s][i] * pi0New[i];
            }
            double delta = Math.max(reldiff1(pi0New, pi0), Math.max(reldiff2(pijNew, pij), reldiff2(phiNew, phi)));
            pij = pijNew; pi0 = pi0New; phi = phiNew;
            if (!Double.isFinite(delta)) break;
            if (delta < tol) break;
        }
        return new Result(pi0, pij, phi);
    }

    private static double reldiff1(double[] a, double[] b) {
        double denom = 0, num = 0;
        for (int i = 0; i < a.length; i++) { denom = Math.max(denom, Math.abs(b[i])); num = Math.max(num, Math.abs(a[i] - b[i])); }
        if (denom == 0) denom = 1;
        return num / denom;
    }

    private static double reldiff2(double[][] a, double[][] b) {
        double denom = 0, num = 0;
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[i].length; j++) { denom = Math.max(denom, Math.abs(b[i][j])); num = Math.max(num, Math.abs(a[i][j] - b[i][j])); }
        if (denom == 0) denom = 1;
        return num / denom;
    }
}
