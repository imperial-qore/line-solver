package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Andersen-Nielsen MAP fitting algorithm.
 */
public final class Map_anfit {
    private Map_anfit() {}

    public static MatrixCell map_anfit(double ls, double rho, double H, int n, int ds) {
        return map_anfit(ls, rho, H, n, ds, null, null, 100, 1e-9);
    }

    public static MatrixCell map_anfit(double ls, double rho, double H, int n, int ds,
                                       double[] SA, int[] SAlags) {
        return map_anfit(ls, rho, H, n, ds, SA, SAlags, 100, 1e-9);
    }

    public static MatrixCell map_anfit(final double ls, final double rho, double H, int n, int ds,
                                       double[] SA, int[] SAlags, int iter_max, double iter_tol) {
        final boolean LSQFIT = SA != null && SAlags != null;
        final double[] SAvals;
        if (LSQFIT) {
            SAvals = new double[SAlags.length];
            for (int i = 0; i < SAlags.length; i++) SAvals[i] = SA[SAlags[i]];
        } else {
            SAvals = null;
        }
        double beta = 2.0 - 2.0 * H;
        Map<Integer, Double> phi = new HashMap<Integer, Double>();
        int d = ds;
        int d0 = 0;
        double a = Math.pow(10.0, (double) n / (double) (d - 1));
        phi.put(d, 1.0);
        int i = 1;

        outer:
        while (true) {
            double S = 0.0;
            for (int j = 0; j < i; j++) {
                Double phiVal = phi.get(d - j);
                if (phiVal == null) phiVal = 0.0;
                S += phiVal * phiVal * Math.exp(1.0 - Math.pow(a, (double) (i - j)));
            }
            double D = Math.pow(a, (double) i * beta) - S;
            if (D < 0) {
                phi.put(d - i, 0.0);
                d0++;
                if (ds > d - d0) {
                    d++;
                    d0 = 0;
                    a = Math.pow(10.0, (double) n / (double) (d - 1));
                    phi.clear();
                    phi.put(d, 1.0);
                    i = 1;
                    continue outer;
                } else {
                    i++;
                    if (i == d) break outer;
                    else continue outer;
                }
            } else {
                phi.put(d - i, Math.sqrt(D));
                i++;
                if (i == d) break outer;
                else continue outer;
            }
        }

        final double[] k1 = new double[d];
        k1[0] = 0.8;
        for (int idx = 2; idx <= d; idx++) k1[idx - 1] = Math.pow(a, 1.0 - (double) idx) * k1[0];
        if (!(k1[0] < 1.0 && rho < 0.5)) {
            System.err.println("warning: if necessary adjust k(2,1) and/or rho");
        }

        double Seta = 0.0;
        for (int idx = 1; idx <= d; idx++) {
            double kappa = k1[idx - 1];
            double e = Math.exp(-kappa);
            Double phiVal = phi.get(idx);
            if (phiVal == null) phiVal = 0.0;
            Seta += phiVal * phiVal * Math.pow(kappa, -2.0)
                    * ((1.0 - e) * (1.0 - e) - 2.0 * rho * (kappa - (1.0 - e)));
        }
        double eta = Math.sqrt(4.0 * rho * ls) / Math.sqrt(Seta);

        double sumPhi = 0.0;
        for (int idx = 1; idx <= d; idx++) {
            Double v = phi.get(idx);
            if (v == null) v = 0.0;
            sumPhi += v;
        }
        double L = eta * sumPhi / 2.0;

        final double[] c1 = new double[d];
        final double[] c2 = new double[d];
        final double[] l = new double[d];
        final double lP;
        if (ls < L) {
            lP = 0.0;
            for (int idx = 1; idx <= d; idx++) {
                Double phiVal = phi.get(idx);
                if (phiVal == null) phiVal = 0.0;
                c1[idx - 1] = L * L / (ls * ls + L * L) * k1[idx - 1];
                c2[idx - 1] = k1[idx - 1] - c1[idx - 1];
                l[idx - 1] = phiVal * (ls * ls + L * L) / (ls * sumPhi);
            }
        } else {
            lP = ls - L;
            for (int idx = 1; idx <= d; idx++) {
                Double phiVal = phi.get(idx);
                if (phiVal == null) phiVal = 0.0;
                c2[idx - 1] = 0.5 * k1[idx - 1];
                c1[idx - 1] = c2[idx - 1];
                l[idx - 1] = eta * phiVal;
            }
        }

        final int dF = d;

        if (!LSQFIT) {
            MatrixCell MAP = buildPoissonMAP(lP);
            for (int idx = 0; idx < d; idx++) {
                MatrixCell IPP = buildIPP(c1[idx], c2[idx], l[idx]);
                MAP = Map_super.map_super(MAP, IPP);
            }
            return Map_normalize.map_normalize(MAP);
        }

        final double[] k0 = new double[d];
        final double[] lsLocal = new double[d];
        for (int idx = 0; idx < d; idx++) {
            double cSum = c1[idx] + c2[idx];
            k0[idx] = l[idx] * l[idx] * c1[idx] * c2[idx] / (cSum * cSum * cSum);
            lsLocal[idx] = c2[idx] * l[idx] / cSum;
        }

        double[] r0 = new double[d];
        for (int idx = 0; idx < d; idx++) {
            r0[idx] = 1.0 + k0[idx] * k1[idx] / (lsLocal[idx] * lsLocal[idx]);
        }

        final Matrix lagsMatrix = new Matrix(1, SAlags.length, SAlags.length);
        for (int idx = 0; idx < SAlags.length; idx++) lagsMatrix.set(idx, (double) SAlags[idx]);

        final double EPSTOL = 100.0 * iter_tol;
        final double PENALTY = 1e10;
        final double lP_ = lP;
        final int dM = dF;
        MultivariateFunction objective = new MultivariateFunction() {
            @Override
            public double value(double[] r) {
                for (int idx = 0; idx < dM; idx++) {
                    if (r[idx] <= EPSTOL) return PENALTY;
                }
                for (int idx = 0; idx < dM; idx++) {
                    double cv = -(lsLocal[idx] - Math.sqrt(k0[idx] * k1[idx] / r[idx]));
                    if (cv > 0) return PENALTY + cv * PENALTY;
                }
                MatrixCell MAP0 = buildPoissonMAP(lP_);
                for (int idx = 0; idx < dM; idx++) {
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, 0.0);
                    D0.set(0, 1, k1[idx] * r[idx] / (1.0 + r[idx]));
                    D0.set(1, 0, k1[idx] / (1.0 + r[idx]));
                    D0.set(1, 1, 0.0);
                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, lsLocal[idx] + Math.sqrt(k0[idx] * k1[idx] * r[idx]));
                    D1.set(0, 1, 0.0);
                    D1.set(1, 0, 0.0);
                    D1.set(1, 1, lsLocal[idx] - Math.sqrt(k0[idx] * k1[idx] / r[idx]));
                    MatrixCell IPP = new MatrixCell();
                    IPP.set(0, D0);
                    IPP.set(1, D1);
                    MAP0 = Map_super.map_super(MAP0, Map_normalize.map_normalize(IPP));
                }
                Matrix acfVals = Map_acf.map_acf(MAP0, lagsMatrix);
                double sumSq = 0.0;
                for (int idx = 0; idx < SAvals.length; idx++) {
                    double diff = acfVals.get(idx) - SAvals[idx];
                    sumSq += diff * diff;
                }
                return Math.sqrt(sumSq);
            }
        };

        MatrixCell bestMAP = null;
        try {
            SimplexOptimizer optimizer = new SimplexOptimizer(iter_tol, iter_tol);
            NelderMeadSimplex simplex = new NelderMeadSimplex(d, 0.1);
            PointValuePair result = optimizer.optimize(
                    new MaxEval(iter_max * 200),
                    new MaxIter(iter_max),
                    new ObjectiveFunction(objective),
                    GoalType.MINIMIZE,
                    new InitialGuess(r0),
                    simplex);
            double[] bestR = result.getPoint();
            MatrixCell MAP0 = buildPoissonMAP(lP);
            for (int idx = 0; idx < d; idx++) {
                Matrix D0 = new Matrix(2, 2, 4);
                D0.set(0, 0, 0.0);
                D0.set(0, 1, k1[idx] * bestR[idx] / (1.0 + bestR[idx]));
                D0.set(1, 0, k1[idx] / (1.0 + bestR[idx]));
                D0.set(1, 1, 0.0);
                Matrix D1 = new Matrix(2, 2, 4);
                D1.set(0, 0, lsLocal[idx] + Math.sqrt(k0[idx] * k1[idx] * bestR[idx]));
                D1.set(0, 1, 0.0);
                D1.set(1, 0, 0.0);
                D1.set(1, 1, lsLocal[idx] - Math.sqrt(k0[idx] * k1[idx] / bestR[idx]));
                MatrixCell IPP = new MatrixCell();
                IPP.set(0, D0);
                IPP.set(1, D1);
                MAP0 = Map_super.map_super(MAP0, Map_normalize.map_normalize(IPP));
            }
            bestMAP = MAP0;
        } catch (Exception e) {
            // Fallback below
        }

        if (bestMAP != null) return Map_normalize.map_normalize(bestMAP);

        MatrixCell MAP = buildPoissonMAP(lP);
        for (int idx = 0; idx < d; idx++) {
            MatrixCell IPP = buildIPP(c1[idx], c2[idx], l[idx]);
            MAP = Map_super.map_super(MAP, IPP);
        }
        return Map_normalize.map_normalize(MAP);
    }

    private static MatrixCell buildPoissonMAP(double rate) {
        MatrixCell MAP = new MatrixCell();
        Matrix D0 = new Matrix(1, 1, 1);
        D0.set(0, 0, -rate);
        Matrix D1 = new Matrix(1, 1, 1);
        D1.set(0, 0, rate);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return MAP;
    }

    private static MatrixCell buildIPP(double c1, double c2, double arrivalRate) {
        Matrix D0 = new Matrix(2, 2, 4);
        D0.set(0, 0, 0.0);
        D0.set(0, 1, c1);
        D0.set(1, 0, c2);
        D0.set(1, 1, 0.0);
        Matrix D1 = new Matrix(2, 2, 4);
        D1.set(0, 0, arrivalRate);
        D1.set(0, 1, 0.0);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, 0.0);
        MatrixCell IPP = new MatrixCell();
        IPP.set(0, D0);
        IPP.set(1, D1);
        return Map_normalize.map_normalize(IPP);
    }

    /** MAP Andersen-Nielsen fitting algorithms. */
    public static final class MapAnfitAlgo {}
}
