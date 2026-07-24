/**
 * @file Asymptotic Cache Miss Analysis
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

public final class Cache_miss_asy {
    private Cache_miss_asy() {}

    /**
     * Compute cache miss rates using asymptotic approximation (Fixed Point Iteration method).
     */
    public static double cache_miss_asy(Matrix gamma, Matrix m, int maxIter, double tolerance) {
        int n = gamma.getNumCols();
        int h = gamma.getNumRows();

        if (m.elementSum() == 0.0 || m.elementMin() < 0.0) {
            return 1.0;
        }

        Matrix pi = Matrix.ones(1, n).scale(1.0 / n);
        Matrix prevPi;

        for (int iter = 0; iter < maxIter; iter++) {
            prevPi = pi.copy();

            Matrix newPi = Matrix.zeros(1, n);

            for (int k = 0; k < n; k++) {
                double numerator = 0.0;
                double denominator = 0.0;

                for (int i = 0; i < h; i++) {
                    int cacheSizeAtLevel = (int) m.get(i);

                    if (cacheSizeAtLevel > 0) {
                        double probNotInCache = 1.0;

                        List<double[]> otherItems = new ArrayList<double[]>();
                        for (int j = 0; j < n; j++) {
                            if (j != k) {
                                otherItems.add(new double[]{(double) j, gamma.get(i, j) * (1.0 - prevPi.get(0, j))});
                            }
                        }

                        Collections.sort(otherItems, new Comparator<double[]>() {
                            @Override
                            public int compare(double[] a, double[] b) {
                                return Double.compare(b[1], a[1]);
                            }
                        });
                        int takeCount = (int) FastMath.min(cacheSizeAtLevel, otherItems.size());
                        List<double[]> topItems = otherItems.subList(0, takeCount);

                        if (topItems.size() < cacheSizeAtLevel) {
                            probNotInCache = 0.0;
                        } else {
                            double weakestInCache = topItems.get(topItems.size() - 1)[1];
                            double itemKPopularity = gamma.get(i, k) * (1.0 - prevPi.get(0, k));
                            if (itemKPopularity > weakestInCache) {
                                probNotInCache = 0.0;
                            }
                        }

                        numerator += gamma.get(i, k) * probNotInCache;
                        denominator += gamma.get(i, k);
                    }
                }

                if (denominator > GlobalConstants.Zero) {
                    newPi.set(0, k, numerator / denominator);
                } else {
                    newPi.set(0, k, 1.0);
                }
            }

            pi = newPi;

            double maxDiff = 0.0;
            for (int k = 0; k < n; k++) {
                maxDiff = FastMath.max(maxDiff, FastMath.abs(pi.get(0, k) - prevPi.get(0, k)));
            }
            if (maxDiff < tolerance) {
                break;
            }
        }

        double globalMissRate = 0.0;
        double totalRate = 0.0;
        for (int i = 0; i < h; i++) {
            for (int k = 0; k < n; k++) {
                globalMissRate += gamma.get(i, k) * pi.get(0, k);
                totalRate += gamma.get(i, k);
            }
        }

        return totalRate > GlobalConstants.Zero ? globalMissRate / totalRate : 1.0;
    }

    public static double cache_miss_asy(Matrix gamma, Matrix m) {
        return cache_miss_asy(gamma, m, 1000, 1e-8);
    }
}
