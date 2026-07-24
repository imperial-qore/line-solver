/**
 * Acyclic Markovian Arrival Process two-phase fitting with autocorrelation
 *
 * Fits AMAP(2) distributions to match moments and correlation characteristics.
 */
package jline.api.mam;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Top-level functions for AMAP2 gamma fitting (ported from Kotlin).
 */
public final class Amap2_fit_gamma {

    private Amap2_fit_gamma() {
    }

    /**
     * Triple holder for adjustment results.
     */
    public static class Triple {
        public final double first;
        public final double second;
        public final double third;
        public Triple(double a, double b, double c) {
            this.first = a;
            this.second = b;
            this.third = c;
        }
    }

    /**
     * Finds an AMAP(2) fitting the given characteristics.
     */
    public static Pair<MatrixCell, List<MatrixCell>> amap2_fit_gamma(double M1, double M2, double M3, double GAMMA) {
        if (Math.abs(M2 - 2 * M1 * M1) < 1e-6) {
            MatrixCell poissonMap = new MatrixCell(2);
            poissonMap.set(0, new Matrix(1, 1));
            poissonMap.get(0).set(0, 0, -1.0 / M1);
            poissonMap.set(1, new Matrix(1, 1));
            poissonMap.get(1).set(0, 0, 1.0 / M1);
            MatrixCell normalizedMap = Map_normalize.map_normalize(poissonMap);
            List<MatrixCell> single = new ArrayList<MatrixCell>();
            single.add(normalizedMap);
            return new Pair<MatrixCell, List<MatrixCell>>(normalizedMap, single);
        }

        List<MatrixCell> amaps = amap2_fitall_gamma(M1, M2, M3, GAMMA);
        List<MatrixCell> normalized = new ArrayList<MatrixCell>();
        for (MatrixCell mc : amaps) normalized.add(Map_normalize.map_normalize(mc));
        amaps = normalized;

        if (amaps.isEmpty()) {
            Triple adj = amap2_adjust_gamma(M1, M2, M3, GAMMA);
            amaps = amap2_fitall_gamma(M1, adj.first, adj.second, adj.third);
            normalized = new ArrayList<MatrixCell>();
            for (MatrixCell mc : amaps) normalized.add(Map_normalize.map_normalize(mc));
            amaps = normalized;
        }

        MatrixCell bestAmap = amaps.isEmpty() ? null : amaps.get(0);
        return new Pair<MatrixCell, List<MatrixCell>>(bestAmap, amaps);
    }

    /**
     * Finds all AMAP(2) solutions for given moments and correlation.
     */
    public static List<MatrixCell> amap2_fitall_gamma(double M1, double M2, double M3, double GAMMA) {
        List<MatrixCell> solutions = new ArrayList<MatrixCell>();
        double degentol = 1e-8;
        double r12tol = 1e-6;

        try {
            double SCV = (M2 - M1 * M1) / (M1 * M1);
            double M3lb = 3 * M1 * M1 * M1 * (3 * SCV - 1 + Math.sqrt(2.0) * Math.pow(1 - SCV, 1.5));

            double tmp0;
            if (SCV <= 1 && Math.abs(M3 - M3lb) < degentol) {
                tmp0 = 0.0;
            } else {
                double term1 = M3 * M3 / 9.0;
                double term2 = (8 * M1 * M1 * M1 / 3.0 - 2 * M2 * M1) * M3;
                double term3 = -3 * M1 * M1 * M2 * M2 + 2 * M2 * M2 * M2;
                tmp0 = term1 + term2 + term3;
            }

            if (tmp0 < 0) {
                return solutions;
            }

            double tmp1 = 3 * Math.sqrt(tmp0);
            double tmp2 = M3 - 3 * M1 * M2;
            double tmp3 = 6 * M2 - 12 * M1 * M1;

            int n = (tmp0 == 0.0) ? 1 : 2;

            double[] h1v = new double[n];
            double[] h2v = new double[n];

            if (n == 1) {
                h2v[0] = tmp2 / tmp3;
                h1v[0] = h2v[0];
            } else {
                h2v[0] = (tmp2 + tmp1) / tmp3;
                h2v[1] = (tmp2 - tmp1) / tmp3;
                h1v[1] = h2v[0];
                h1v[0] = h2v[1];
            }

            double minH2 = Double.POSITIVE_INFINITY;
            for (double v : h2v) if (v < minH2) minH2 = v;
            if (minH2 <= 0) {
                return solutions;
            }

            for (int j = 0; j < n; j++) {
                double h1 = h1v[j];
                double h2 = h2v[j];

                if (GAMMA >= 0) {
                    double z = M1 * M1 * GAMMA * GAMMA +
                            (2 * M1 * h1 + 2 * M1 * h2 - 4 * h1 * h2 - 2 * M1 * M1) * GAMMA +
                            M1 * M1 - 2 * M1 * h1 - 2 * M1 * h2 + h1 * h1 + 2 * h1 * h2 + h2 * h2;

                    if (Math.abs(z) < degentol) {
                        double r2 = (h1 - M1 + h2 + GAMMA * M1) / (2 * h1);
                        double r1 = (M1 - h1 - M1 * r2 + h1 * r2) / (h2 - M1 * r2);

                        if (isFeasible(r1, r2, r12tol)) {
                            double fixedR1 = Math.max(Math.min(r1, 1.0), 0.0);
                            double fixedR2 = Math.max(Math.min(r2, 1.0), 0.0);
                            solutions.add(amap2_assemble(h1, h2, fixedR1, fixedR2, 1));
                        }
                    } else if (z > 0) {
                        double[] r2v = new double[] {
                                (h1 - M1 + h2 - Math.sqrt(z) + GAMMA * M1) / (2 * h1),
                                (h1 - M1 + h2 + Math.sqrt(z) + GAMMA * M1) / (2 * h1)
                        };

                        for (int i = 0; i < r2v.length; i++) {
                            double r2 = r2v[i];
                            double r1 = (M1 - h1 - M1 * r2 + h1 * r2) / (h2 - M1 * r2);

                            if (isFeasible(r1, r2, r12tol)) {
                                double fixedR1 = Math.max(Math.min(r1, 1.0), 0.0);
                                double fixedR2 = Math.max(Math.min(r2, 1.0), 0.0);
                                solutions.add(amap2_assemble(h1, h2, fixedR1, fixedR2, 1));
                            }
                        }
                    }
                } else {
                    double r2 = (h1 - M1 + h2 + GAMMA * M1) / h1;
                    double r1 = (r2 + (h1 + h2 - h1 * r2) / M1 - 2) / (r2 - 1);

                    if (isFeasible(r1, r2, r12tol)) {
                        double fixedR1 = Math.max(Math.min(r1, 1.0), 0.0);
                        double fixedR2 = Math.max(Math.min(r2, 1.0), 0.0);
                        solutions.add(amap2_assemble(h1, h2, fixedR1, fixedR2, 2));
                    }
                }
            }

        } catch (Exception e) {
            // ignore
        }

        return solutions;
    }

    private static boolean isFeasible(double r1, double r2, double tol) {
        return !Double.isNaN(r1) && !Double.isNaN(r2) && Double.isFinite(r1) && Double.isFinite(r2) &&
                r1 >= -tol && r1 <= 1 + tol &&
                r2 >= -tol && r2 <= 1 + tol;
    }

    public static MatrixCell amap2_assemble(double l1, double l2, double p1, double p2, int form) {
        MatrixCell amap = new MatrixCell(2);
        Matrix D0 = new Matrix(2, 2);
        Matrix D1 = new Matrix(2, 2);

        if (form == 1) {
            D0.set(0, 0, -1.0 / l1);
            D0.set(0, 1, p1 / l1);
            D0.set(1, 0, 0.0);
            D0.set(1, 1, -1.0 / l2);

            D1.set(0, 0, (1 - p1) / l1);
            D1.set(0, 1, 0.0);
            D1.set(1, 0, (1 - p2) / l2);
            D1.set(1, 1, p2 / l2);
        } else if (form == 2) {
            D0.set(0, 0, -1.0 / l1);
            D0.set(0, 1, p1 / l1);
            D0.set(1, 0, 0.0);
            D0.set(1, 1, -1.0 / l2);

            D1.set(0, 0, 0.0);
            D1.set(0, 1, (1 - p1) / l1);
            D1.set(1, 0, (1 - p2) / l2);
            D1.set(1, 1, p2 / l2);
        } else {
            throw new IllegalArgumentException("Invalid form: should be either 1 (gamma > 0) or 2 (gamma < 0)");
        }

        amap.set(0, D0);
        amap.set(1, D1);
        return amap;
    }

    public static Triple amap2_adjust_gamma(double M1, double M2, double M3, double GAMMA) {
        return amap2_adjust_gamma(M1, M2, M3, GAMMA, new double[]{10.0, 1.0, 10.0}, 3, 2);
    }

    public static Triple amap2_adjust_gamma(double M1, double M2, double M3, double GAMMA,
                                             double[] weights, int method, int constraints) {
        double tol = 1e-2;
        switch (method) {
            case 1:
                return adjustMethod1(M1, M2, M3, GAMMA, weights, tol);
            case 2:
                return adjustMethod2(M1, M2, M3, GAMMA, weights, tol);
            case 3:
                return adjustMethod3(M1, M2, M3, GAMMA, tol);
            case 4:
                return adjustMethod4(M1, M2, M3, GAMMA, tol);
            default:
                throw new IllegalArgumentException("Invalid method for adjusting AMAP(2) characteristics");
        }
    }

    private static Triple adjustMethod3(double M1, double M2, double M3, double GAMMA, double tol) {
        Map<Integer, Double> adjusted = Aph2_adjust.aph2_adjust(M1, M2, M3, "simple");
        double M2a = adjusted.get(0);
        double M3a = adjusted.get(1);

        double[] bounds = computeGammaBounds(M1, M2a, M3a, tol);
        double lb = bounds[0];
        double ub = bounds[1];
        double GAMMAa = Math.max(lb, Math.min(GAMMA, ub));
        return new Triple(M2a, M3a, GAMMAa);
    }

    private static Triple adjustMethod4(double M1, double M2, double M3, double GAMMA, double tol) {
        double M1sq = M1 * M1;
        double scv = (M2 - M1sq) / M1sq;

        double M2a = (scv < 0.5) ? 1.5 * M1sq : M2;
        double scva = (M2a - M1sq) / M1sq;

        double M3_lb;
        double M3_ub;
        if (scva <= 1) {
            M3_lb = 3 * M1 * M1 * M1 * (3 * scva - 1 + Math.sqrt(2.0) * Math.pow(1 - scva, 1.5));
            M3_ub = 6 * M1 * M1 * M1 * scva;
        } else {
            M3_lb = 1.5 * M1 * M1 * M1 * (1 + scva) * (1 + scva);
            M3_ub = GlobalConstants.Inf;
        }

        double M3a;
        double GAMMAa;

        if (Math.abs(M3_lb - M3_ub) < tol) {
            M3a = (M3_lb + M3_ub) / 2;
            GAMMAa = 0.0;
        } else {
            M3a = optimizeM3ForGamma(M1, M2a, M3, GAMMA, M3_lb, M3_ub, tol);
            double[] bounds = computeGammaBounds(M1, M2a, M3a, tol);
            GAMMAa = Math.max(bounds[0], Math.min(GAMMA, bounds[1]));
        }

        return new Triple(M2a, M3a, GAMMAa);
    }

    private static Triple adjustMethod1(double M1, double M2, double M3, double GAMMA, double[] weights, double tol) {
        double[] target = new double[]{M2, M3, GAMMA};
        Triple bestSolution = new Triple(M2, M3, GAMMA);
        double bestObjective = Double.MAX_VALUE;

        double[] m2Range = generateRange(1.5 * M1 * M1, Math.max(M2 * 2, 5 * M1 * M1), 20);
        double[] m3Range = generateRange(M1 * M1 * M1, Math.max(M3 * 2, 10 * M1 * M1 * M1), 20);
        double[] gammaRange = generateRange(-0.99, 0.99, 20);

        for (double m2Test : m2Range) {
            for (double m3Test : m3Range) {
                for (double gammaTest : gammaRange) {
                    if (isAmap2Feasible(M1, m2Test, m3Test, gammaTest)) {
                        double[] candidate = new double[]{m2Test, m3Test, gammaTest};
                        double objective = computeWeightedDistance(candidate, target, weights);

                        if (objective < bestObjective) {
                            bestObjective = objective;
                            bestSolution = new Triple(m2Test, m3Test, gammaTest);
                        }
                    }
                }
            }
        }

        return bestSolution;
    }

    private static Triple adjustMethod2(double M1, double M2, double M3, double GAMMA, double[] weights, double tol) {
        double M2a = Math.max(1.5 * M1 * M1, M2);

        double n2 = M2a / (M1 * M1);
        double M3_LB;
        double M3_UB;
        if (n2 >= 1.5 && n2 < 2) {
            double p2 = 3 * (n2 - 2) / (3 * n2) * (-2 * Math.sqrt(3.0) / Math.sqrt(12 - 6 * n2) - 1);
            double a2 = (n2 - 2) / (p2 * (1 - n2) + Math.sqrt(p2 * p2 + 2 * p2 * (n2 - 2)));
            double l2 = 3 * (a2 + 1) / (a2 * p2 + 1) - 6 * a2 / (2 + a2 * p2 * (2 * a2 + 2));
            double u2 = 6 * (n2 - 1) / n2;
            M3_LB = l2 * M1 * M2a;
            M3_UB = u2 * M1 * M2a;
        } else {
            M3_LB = 1.5 * M2a * M2a / M1 + tol;
            M3_UB = GlobalConstants.Inf;
        }

        double M3a = Math.max(M3_LB, Math.min(M3_UB, M3));
        double[] bounds = computeGammaBounds(M1, M2a, M3a, tol);
        double GAMMAa = Math.max(bounds[0], Math.min(GAMMA, bounds[1]));

        return new Triple(M2a, M3a, GAMMAa);
    }

    private static double[] computeGammaBounds(double M1, double M2, double M3, double tol) {
        double n2 = M2 / (M1 * M1);
        double n3 = M3 / (M2 * M1);

        if (n2 < 2) {
            double lb = -(n2 * (n3 - 6) + 6) / (3 * n2 - 6);
            double ub = -(2 * Math.pow(0.5 * (n2 - 2) + 0.5 * Math.sqrt(n2 * n2 - 2 * n2 * n3 / 3), 2.0)) / (n2 - 2);
            return new double[]{lb, ub * (1 - tol)};
        } else if (n3 < 9 - 12 / n2) {
            double lb = -(n2 * (n3 - 6) + 6) / (3 * n2 - 6);
            return new double[]{lb, 1 - tol};
        } else {
            double x1 = Math.sqrt(n2 * (n2 * (18 * n2 + n3 * (n3 - 18) - 27) + 24 * n3));
            double x2 = n2 * (n3 - 9);
            double lb = (x2 - x1 + 12) / (x2 + x1 + 12);
            return new double[]{lb, 1 - tol};
        }
    }

    private static double[] generateRange(double min, double max, int count) {
        double[] arr = new double[count];
        for (int i = 0; i < count; i++) arr[i] = min + (max - min) * i / (count - 1);
        return arr;
    }

    private static double computeWeightedDistance(double[] candidate, double[] target, double[] weights) {
        double total = 0.0;
        for (int i = 0; i < candidate.length; i++) {
            double relativeError = (candidate[i] - target[i]) / target[i] * weights[i];
            total += relativeError * relativeError;
        }
        return total;
    }

    private static boolean isAmap2Feasible(double M1, double M2, double M3, double GAMMA) {
        try {
            return !amap2_fitall_gamma(M1, M2, M3, GAMMA).isEmpty();
        } catch (Exception e) {
            return false;
        }
    }

    private static double optimizeM3ForGamma(double M1, double M2, double M3, double GAMMA,
                                              double M3_lb, double M3_ub, double tol) {
        double phi = (1 + Math.sqrt(5.0)) / 2;
        double resphi = 2 - phi;

        double a = M3_lb + tol;
        double b = Double.isFinite(M3_ub) ? M3_ub : M3_lb + 10 * Math.abs(M3 - M3_lb);

        if (!Double.isFinite(b)) b = a + Math.abs(M3);

        int maxIter = 50;
        for (int i = 0; i < maxIter; i++) {
            double c = a + resphi * (b - a);
            double d = a + (1 - resphi) * (b - a);

            double objC = gammaObjective(M1, M2, c, GAMMA);
            double objD = gammaObjective(M1, M2, d, GAMMA);

            if (objC < objD) {
                b = d;
            } else {
                a = c;
            }

            if (Math.abs(b - a) < tol) break;
        }

        return (a + b) / 2;
    }

    private static double gammaObjective(double M1, double M2, double M3, double targetGamma) {
        double[] bounds = computeGammaBounds(M1, M2, M3, 1e-6);
        double adjustedGamma = Math.max(bounds[0], Math.min(targetGamma, bounds[1]));
        return Math.abs(adjustedGamma - targetGamma);
    }
}
