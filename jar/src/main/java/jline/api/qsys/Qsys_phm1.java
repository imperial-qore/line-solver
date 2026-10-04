/**
 * @file PH/M/1 queueing system analysis
 *
 * Solves the GI/M/1 fundamental equation
 *     sigma = psi_A(mu * (1 - sigma))
 * with the LST of a PH inter-arrival distribution
 *     psi_A(s) = alpha (s I - T)^{-1} (-T 1).
 * Returns time-average performance metrics.
 */
package jline.api.qsys;

import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Map;

public final class Qsys_phm1 {
    private Qsys_phm1() {}

    public static PhM1Result qsys_phm1(double[] alpha, double[][] T, double mu, double tol, int maxIter) {
        if (!(mu > 0.0)) throw new IllegalArgumentException("Service rate mu must be positive");
        int k = T.length;
        if (!(k > 0)) throw new IllegalArgumentException("T must be non-empty");
        for (int i = 0; i < k; i++) {
            if (T[i].length != k) throw new IllegalArgumentException("T must be square");
        }
        if (alpha.length != k) throw new IllegalArgumentException("alpha length must match T dimension");

        final RealMatrix Tm = new Array2DRowRealMatrix(T, false);
        final RealMatrix negTm = Tm.scalarMultiply(-1.0);

        double[] ones = new double[k];
        for (int i = 0; i < k; i++) ones[i] = 1.0;
        double[] tVec = negTm.operate(ones);  // -T 1

        LUDecomposition negTmLU = new LUDecomposition(negTm);
        if (!negTmLU.getSolver().isNonSingular()) {
            throw new IllegalArgumentException("PH sub-generator T is singular");
        }
        double[] negTmInvOnes = negTmLU.getSolver().solve(new ArrayRealVector(ones)).toArray();
        double meanIa = dotProduct(alpha, negTmInvOnes);
        if (meanIa <= 0.0) throw new IllegalArgumentException("Non-positive mean inter-arrival: " + meanIa);
        double lambda = 1.0 / meanIa;
        double rho = lambda / mu;
        if (rho >= 1.0 - 1e-12) throw new IllegalArgumentException("Load rho=" + rho + " must be strictly less than 1");

        final RealMatrix I = MatrixUtils.createRealIdentityMatrix(k);

        // Bracket and bisect/Brent the root in (0, 1) excluding the trivial sigma=1.
        double left = 1e-12;
        double right = 1.0 - 1e-12;
        double fL = fEval(left, alpha, Tm, I, tVec, mu);
        double fR = fEval(right, alpha, Tm, I, tVec, mu);
        double sigma;
        if (fL * fR > 0.0) {
            // Fall back to fixed-point iteration
            sigma = 0.5;
            for (int it = 0; it < maxIter; it++) {
                double sNew = lst(mu * (1.0 - sigma), alpha, Tm, I, tVec);
                if (Math.abs(sNew - sigma) < tol) {
                    sigma = sNew;
                    break;
                }
                sigma = sNew;
            }
        } else {
            // Bisection
            double lo = left;
            double hi = right;
            double fLo = fL;
            for (int it = 0; it < 200; it++) {
                double mid = 0.5 * (lo + hi);
                double fMid = fEval(mid, alpha, Tm, I, tVec, mu);
                if (Math.abs(fMid) < tol || (hi - lo) < tol) {
                    lo = mid; hi = mid;
                    break;
                }
                if (fLo * fMid < 0.0) {
                    hi = mid;
                } else {
                    lo = mid; fLo = fMid;
                }
            }
            sigma = 0.5 * (lo + hi);
        }

        double L = rho / (1.0 - sigma);
        double Lq = rho * sigma / (1.0 - sigma);
        double Wq = Lq / lambda;
        double W = Wq + 1.0 / mu;
        return new PhM1Result(L, Lq, Wq, W, rho, sigma);
    }

    public static PhM1Result qsys_phm1(double[] alpha, double[][] T, double mu, double tol) {
        return qsys_phm1(alpha, T, mu, tol, 200);
    }

    public static PhM1Result qsys_phm1(double[] alpha, double[][] T, double mu) {
        return qsys_phm1(alpha, T, mu, 1e-12, 200);
    }

    private static double lst(double s, double[] alpha, RealMatrix Tm, RealMatrix I, double[] tVec) {
        RealMatrix M = I.scalarMultiply(s).subtract(Tm);
        LUDecomposition lu = new LUDecomposition(M);
        double[] v = lu.getSolver().solve(new ArrayRealVector(tVec)).toArray();
        return dotProduct(alpha, v);
    }

    private static double fEval(double sigma, double[] alpha, RealMatrix Tm, RealMatrix I, double[] tVec, double mu) {
        return sigma - lst(mu * (1.0 - sigma), alpha, Tm, I, tVec);
    }

    private static double dotProduct(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    }

    /**
     * Extract a (alpha, T) PH representation for the qsys_phm1 sigma-root from the
     * sn.proc map at a given station/class. Recognizes the (D0, D1) MAP storage
     * form used by LINE's NetworkStruct and computes alpha via map_pie.
     * Returns null if the station does not have a usable PH representation.
     */
    public static Pair<double[], double[][]> extractPhPairForPhm1(NetworkStruct sn, int stationIdx, int classIdx) {
        try {
            if (stationIdx < 0 || stationIdx >= sn.stations.size()) return null;
            Station station = sn.stations.get(stationIdx);
            if (classIdx < 0 || classIdx >= sn.jobclasses.size()) return null;
            JobClass jobClass = sn.jobclasses.get(classIdx);
            Map<JobClass, MatrixCell> procStation = sn.proc.get(station);
            if (procStation == null) return null;
            MatrixCell cell = procStation.get(jobClass);
            if (cell == null) return null;
            if (cell.size() < 2) return null;
            Matrix D0 = cell.get(0);
            Matrix D1 = cell.get(1);
            if (D0 == null || D1 == null) return null;
            int k = D0.getNumRows();
            if (D0.getNumCols() != k || D1.getNumRows() != k || D1.getNumCols() != k) return null;
            // alpha = stationary entry distribution = map_pie(D0, D1)
            Matrix pieMatrix = Map_pie.map_pie(new MatrixCell(D0, D1));
            double[] alpha = new double[k];
            for (int i = 0; i < k; i++) alpha[i] = pieMatrix.get(0, i);
            double[][] T = new double[k][k];
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    T[i][j] = D0.get(i, j);
                }
            }
            return new Pair<double[], double[][]>(alpha, T);
        } catch (Exception e) {
            return null;
        }
    }

    public static Pair<double[], double[][]> extractPhPairForPhm1(NetworkStruct sn, int stationIdx) {
        return extractPhPairForPhm1(sn, stationIdx, 0);
    }
}
