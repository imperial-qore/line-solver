/**
 * @file Absorbing Phase-type distribution moment adjustment
 *
 * Adjusts moments to ensure feasibility bounds for APH(2) fitting procedures.
 * Essential for stabilizing parameter estimation when input moments are infeasible.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.util.Pair;

public final class Aph2_adjust {
    private Aph2_adjust() {}

    /**
     * Adjusts the second and third moments (M2 and M3) of a distribution using the default "simple"
     * method for fitting an APH(2) distribution.
     *
     * @param M1 first moment (mean)
     * @param M2 second moment
     * @param M3 third moment
     * @return map with keys 0 and 1, representing the adjusted M2 and M3 respectively
     */
    public static Map<Integer, Double> aph2_adjust(double M1, double M2, double M3, String method) {
        double tol = 1e-4;
        double M2a = 0.0;
        double M3a = 0.0;
        if ("simple".equals(method)) {
            double scva;
            double M1sq = FastMath.pow(M1, 2);
            double scv = (M2 - M1sq) / M1sq;
            if (scv < 0.5) {
                M2a = 1.5 * FastMath.pow(M1, 2);
                scva = (M2a - M1sq) / M1sq;
            } else {
                M2a = M2;
                scva = scv;
            }
            if (scva < 1) {
                double lb = 3 * FastMath.pow(M1, 3) * (3 * scva - 1 + FastMath.sqrt(2.0) * FastMath.pow(1 - scva, 1.5));
                double ub = 6 * FastMath.pow(M1, 3) * scva;
                if (M3 < lb) {
                    M3a = lb;
                } else if (M3 > ub) {
                    M3a = ub;
                } else {
                    M3a = M3;
                }
            } else {
                double lb = 1.5 * FastMath.pow(M1, 3) * FastMath.pow(1 + scva, 2);
                if (M3 < lb) {
                    M3a = lb * (1 + tol);
                } else {
                    M3a = M3;
                }
            }
        } else if ("opt_char".equals(method) || "opt_char_gads".equals(method)) {
            Pair<Double, Double> results = optimizeInCharacteristicsSpace(M1, M2, M3, method.contains("gads"));
            M2a = results.getLeft();
            M3a = results.getRight();
        } else if ("opt_param".equals(method) || "opt_param_gads".equals(method)) {
            Pair<Double, Double> results = optimizeInParameterSpace(M1, M2, M3, method.contains("gads"));
            M2a = results.getLeft();
            M3a = results.getRight();
        } else {
            throw new IllegalArgumentException("Invalid method: " + method);
        }

        Map<Integer, Double> result = new HashMap<Integer, Double>();
        result.put(0, M2a);
        result.put(1, M3a);
        return result;
    }

    public static Map<Integer, Double> aph2_adjust(double M1, double M2, double M3) {
        return aph2_adjust(M1, M2, M3, "simple");
    }

    /**
     * Optimize in characteristics space (moment space).
     */
    private static Pair<Double, Double> optimizeInCharacteristicsSpace(double M1, double M2, double M3, boolean useGlobal) {
        int maxIter = useGlobal ? 100 : 50;
        double tolerance = 1e-6;

        Pair<Double, Double> result1 = optimizeForFitType(M1, M2, M3, 1, maxIter, tolerance);
        Pair<Double, Double> result2 = optimizeForFitType(M1, M2, M3, 2, maxIter, tolerance);

        double dist1 = Math.sqrt(Math.pow(result1.getLeft() - M2, 2) + Math.pow(result1.getRight() - M3, 2));
        double dist2 = Math.sqrt(Math.pow(result2.getLeft() - M2, 2) + Math.pow(result2.getRight() - M3, 2));

        return dist1 < dist2 ? result1 : result2;
    }

    /**
     * Optimize in parameter space.
     */
    private static Pair<Double, Double> optimizeInParameterSpace(double M1, double M2, double M3, boolean useGlobal) {
        int maxIter = useGlobal ? 100 : 50;
        double feastol = 1e-6;
        double degentol = 1e-8;

        double bestL2 = M1;
        double bestR1 = 0.5;
        double bestObj = Double.MAX_VALUE;

        double l2Min = feastol;
        double l2Max = M1 * 10;
        double r1Min = degentol;
        double r1Max = 1.0 - degentol;

        int gridSize = useGlobal ? 20 : 10;

        for (int i = 0; i <= gridSize; i++) {
            for (int j = 0; j <= gridSize; j++) {
                double l2 = l2Min + (l2Max - l2Min) * i / gridSize;
                double r1 = r1Min + (r1Max - r1Min) * j / gridSize;

                if (l2 * r1 > M1 - feastol) continue;

                double l1 = M1 - l2 * r1;
                if (l1 < feastol) continue;

                double xM2 = 2 * Math.pow(l1, 2) + 2 * r1 * l1 * l2 + 2 * r1 * Math.pow(l2, 2);
                double xM3 = 6 * Math.pow(l1, 3) + 6 * r1 * Math.pow(l1, 2) * l2 + 6 * r1 * l1 * Math.pow(l2, 2) + 6 * r1 * Math.pow(l2, 3);

                double obj = Math.sqrt(Math.pow(M2 - xM2, 2) + Math.pow(M3 - xM3, 2));

                if (obj < bestObj) {
                    bestObj = obj;
                    bestL2 = l2;
                    bestR1 = r1;
                }
            }
        }

        for (int iter = 0; iter < maxIter; iter++) {
            double stepSize = 0.01 * (1.0 - (double) iter / maxIter);
            boolean improved = false;

            double[] deltas = new double[]{-stepSize, stepSize};

            for (double dl2 : deltas) {
                for (double dr1 : deltas) {
                    double newL2 = Math.max(l2Min, Math.min(l2Max, bestL2 + dl2));
                    double newR1 = Math.max(r1Min, Math.min(r1Max, bestR1 + dr1));

                    if (newL2 * newR1 > M1 - feastol) continue;

                    double l1 = M1 - newL2 * newR1;
                    if (l1 < feastol) continue;

                    double xM2 = 2 * Math.pow(l1, 2) + 2 * newR1 * l1 * newL2 + 2 * newR1 * Math.pow(newL2, 2);
                    double xM3 = 6 * Math.pow(l1, 3) + 6 * newR1 * Math.pow(l1, 2) * newL2 + 6 * newR1 * l1 * Math.pow(newL2, 2) + 6 * newR1 * Math.pow(newL2, 3);

                    double obj = Math.sqrt(Math.pow(M2 - xM2, 2) + Math.pow(M3 - xM3, 2));

                    if (obj < bestObj) {
                        bestObj = obj;
                        bestL2 = newL2;
                        bestR1 = newR1;
                        improved = true;
                    }
                }
            }

            if (!improved && iter > 10) break;
        }

        double l1 = M1 - bestL2 * bestR1;
        double M2a = 2 * Math.pow(l1, 2) + 2 * bestR1 * l1 * bestL2 + 2 * bestR1 * Math.pow(bestL2, 2);
        double M3a = 6 * Math.pow(l1, 3) + 6 * bestR1 * Math.pow(l1, 2) * bestL2 + 6 * bestR1 * l1 * Math.pow(bestL2, 2) + 6 * bestR1 * Math.pow(bestL2, 3);

        return new Pair<Double, Double>(M2a, M3a);
    }

    /**
     * Optimize for a specific fit type (1 or 2).
     */
    private static Pair<Double, Double> optimizeForFitType(double M1, double M2, double M3, int fitType, int maxIter, double tolerance) {
        double bestM2 = M2;
        double bestM3 = M3;
        double bestObj = 0.0;

        double stepSize = Math.min(M2, M3) * 0.01;

        for (int iter = 0; iter < maxIter; iter++) {
            boolean improved = false;
            double currentStepSize = stepSize * (1.0 - (double) iter / maxIter);

            double[] deltas = new double[]{-currentStepSize, currentStepSize};

            for (double dM2 : deltas) {
                for (double dM3 : deltas) {
                    double newM2 = Math.max(tolerance, bestM2 + dM2);
                    double newM3 = Math.max(tolerance, bestM3 + dM3);

                    if (isFeasible(M1, newM2, newM3, fitType)) {
                        double obj = Math.sqrt(Math.pow(M2 - newM2, 2) + Math.pow(M3 - newM3, 2));

                        if (iter == 0 || obj < bestObj) {
                            bestObj = obj;
                            bestM2 = newM2;
                            bestM3 = newM3;
                            improved = true;
                        }
                    }
                }
            }

            if (!improved && iter > 10) break;
        }

        return new Pair<Double, Double>(bestM2, bestM3);
    }

    /**
     * Check feasibility of M2, M3 for APH(2) fitting.
     */
    private static boolean isFeasible(double M1, double M2, double M3, int fitType) {
        try {
            double tmp0 = (8 * Math.pow(M1, 3) * M3) / 3 - 3 * Math.pow(M1, 2) * Math.pow(M2, 2) - 2 * M1 * M2 * M3 + 2 * Math.pow(M2, 3) + Math.pow(M3, 2) / 9;

            if (tmp0 < 0) return false;

            double tmp1 = 3 * Math.sqrt(tmp0);
            double tmp2 = M3 - 3 * M1 * M2;
            double tmp3 = 6 * M2 - 12 * Math.pow(M1, 2);

            if (Math.abs(tmp3) < 1e-10) return false;

            double l1;
            double l2;
            if (fitType == 1) {
                l1 = (tmp2 + tmp1) / tmp3;
                l2 = (tmp2 - tmp1) / tmp3;
            } else {
                l1 = (tmp2 - tmp1) / tmp3;
                l2 = (tmp2 + tmp1) / tmp3;
            }

            if (l1 < 1e-6 || l2 < 1e-6) return false;

            double p1 = (M1 - l1) / l2;
            if (p1 < 0 || p1 > 1) return false;

            return true;
        } catch (Exception e) {
            return false;
        }
    }
}
