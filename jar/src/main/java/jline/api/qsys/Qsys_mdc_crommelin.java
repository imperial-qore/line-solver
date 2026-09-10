/**
 * @file M/D/c queueing system analysis (Poisson arrivals, deterministic service)
 *
 * Computes the exact queue-length distribution and Lq for an M/D/c queue
 * using Crommelin's embedded DTMC at multiples of the service time.
 *
 * The DTMC is X_{n+1} = max(0, X_n - c) + A_n where A_n ~ Poisson(lambda*s),
 * which is the "synchronized service start" embedding used for M/D/c with FCFS.
 * Truncates the state space at a level chosen so that the geometric tail
 * contribution falls below the requested tolerance.
 */
package jline.api.qsys;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;

public final class Qsys_mdc_crommelin {
    private Qsys_mdc_crommelin() {}

    /**
     * Solve M/D/c via the Crommelin embedded DTMC.
     *
     * @param lambdaArr Arrival rate (Poisson)
     * @param s Deterministic service time (>0)
     * @param c Number of servers (>=1)
     * @param truncation Truncation level. If <=0, chosen automatically.
     */
    public static MDcCrommelinResult qsys_mdc_crommelin(double lambdaArr, double s, int c, int truncation) {
        if (!(lambdaArr > 0.0)) throw new IllegalArgumentException("Arrival rate must be positive");
        if (!(s > 0.0)) throw new IllegalArgumentException("Service time must be positive");
        if (!(c >= 1)) throw new IllegalArgumentException("Number of servers must be >= 1");

        double a = lambdaArr * s;
        double rho = a / c;
        if (!(rho < 1.0 - 1e-12)) throw new IllegalArgumentException("Load rho=" + rho + " must be strictly less than 1");

        // see _kb/03-api-layer.md for rationale
        int autoN = Math.max(200, Math.min(2500, (int) (10.0 / (1.0 - rho)) + 200));
        int nMax = (truncation > 0) ? truncation : autoN;

        // Precompute Poisson(a) PMF up to nMax (log-space to avoid overflow)
        double logA = Math.log(a);
        double[] logPmf = new double[nMax + 1];
        double logFact = 0.0;
        for (int k = 0; k <= nMax; k++) {
            logPmf[k] = -a + k * logA - logFact;
            logFact += Math.log((double) (k + 1));
        }
        double[] pmf = new double[nMax + 1];
        for (int k = 0; k <= nMax; k++) pmf[k] = Math.exp(logPmf[k]);

        // Build (P^T - I), replace last row with all-ones for normalization.
        // P[i, j] = pmf[j - max(0, i - c)] when j >= max(0, i - c) else 0.
        int n = nMax + 1;
        double[][] data = new double[n][n];
        for (int i = 0; i < n; i++) {
            int base = (i <= c) ? 0 : i - c;
            for (int j = base; j < n; j++) {
                int k = j - base;
                if (k <= nMax) data[j][i] = pmf[k];  // P^T entry
            }
        }
        for (int i = 0; i < n; i++) data[i][i] -= 1.0;
        // Replace last row with sum(pi) = 1 normalization
        for (int j = 0; j < n; j++) data[n - 1][j] = 1.0;

        double[] rhs = new double[n];
        rhs[n - 1] = 1.0;

        // Use Apache Commons LU directly to bypass Matrix.solve's bogus
        // singularity check (det() underflows for large n even on well-conditioned A).
        LUDecomposition lu = new LUDecomposition(new Array2DRowRealMatrix(data, false));
        double[] piArr = lu.getSolver().solve(new ArrayRealVector(rhs)).toArray();

        double meanN = 0.0;
        double Lq = 0.0;
        for (int i = 0; i < n; i++) {
            double p = piArr[i];
            meanN += i * p;
            if (i > c) Lq += (i - c) * p;
        }

        double Wq = Lq / lambdaArr;
        double W = Wq + s;
        return new MDcCrommelinResult(meanN, Lq, Wq, W, rho);
    }

    public static MDcCrommelinResult qsys_mdc_crommelin(double lambdaArr, double s, int c) {
        return qsys_mdc_crommelin(lambdaArr, s, c, -1);
    }
}
