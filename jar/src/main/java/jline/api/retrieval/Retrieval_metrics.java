/**
 * @file Retrieval_metrics.java
 * @brief Exact miss/hit/delayed-hit metrics of a delayed-hit (list-based) cache.
 *
 * Port of matlab/src/api/retrieval/retrieval_metrics.m (paper
 * prop:performance_measures): with E(m)=E(0,m) and E_i the constant without item i,
 *   pi_{i,0} = E_i(m)/E,                                   (miss)
 *   pi_{i,j} = m_j*gamma_{i,j}*E_i(m-1_j)/E,               (hit in list j)
 *   phi_{0,i} = lambda_i*eta_{0,i}*E_i(m)/E,               (delayed hit, IS)
 *   phi_{s,i} = lambda_i*eta_{s,i}*E_i(1_s,m)/E.           (delayed hit, PS s)
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

public final class Retrieval_metrics {
    private Retrieval_metrics() {}

    /** Result: pmiss (1 x n), phit (h x n), pdh ((r+1) x n). */
    public static final class Result {
        public final double[] pmiss;
        public final double[][] phit;
        public final double[][] pdh;
        Result(double[] pmiss, double[][] phit, double[][] pdh) {
            this.pmiss = pmiss; this.phit = phit; this.pdh = pdh;
        }
    }

    public static Result retrieval_metrics(double[] m, double[] lambda, double[][] eta, double[][] gamma) {
        int n = lambda.length;
        int h = m.length;
        int r = eta[0].length - 1;
        double[] v0 = new double[r];

        double E = Retrieval_nc.retrieval_nc(v0, m, lambda, eta, gamma);

        double[] pmiss = new double[n];
        double[][] phit = new double[h][n];
        double[][] pdh = new double[r + 1][n];

        for (int i = 0; i < n; i++) {
            // keep = all items but i
            double[] lambda_i = new double[n - 1];
            double[][] eta_i = new double[n - 1][];
            double[][] gamma_i = new double[n - 1][];
            int idx = 0;
            for (int k = 0; k < n; k++) {
                if (k == i) continue;
                lambda_i[idx] = lambda[k];
                eta_i[idx] = eta[k];
                gamma_i[idx] = gamma[k];
                idx++;
            }

            double Ei = Retrieval_nc.retrieval_nc(v0, m, lambda_i, eta_i, gamma_i);
            pmiss[i] = Ei / E;
            pdh[0][i] = lambda[i] * eta[i][0] * Ei / E;

            for (int s = 1; s <= r; s++) {
                double[] vs = new double[r];
                vs[s - 1] = 1;
                double Eis = Retrieval_nc.retrieval_nc(vs, m, lambda_i, eta_i, gamma_i);
                pdh[s][i] = lambda[i] * eta[i][s] * Eis / E;
            }

            for (int j = 0; j < h; j++) {
                if (m[j] > 0) {
                    double[] mp = m.clone();
                    mp[j] -= 1;
                    double Eij = Retrieval_nc.retrieval_nc(v0, mp, lambda_i, eta_i, gamma_i);
                    phit[j][i] = m[j] * gamma[i][j] * Eij / E;
                }
            }
        }
        return new Result(pmiss, phit, pdh);
    }
}
