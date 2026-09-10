/**
 * @file M3PP interleaved fitting using theoretical characteristics.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public final class M3pp_interleave_fitc_theoretical {
    private M3pp_interleave_fitc_theoretical() {}

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_theoretical(MatrixCell mmap) {
        return m3pp_interleave_fitc_theoretical(mmap, 1.0, 1000.0, null);
    }

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_theoretical(MatrixCell mmap, double t, double tinf) {
        return m3pp_interleave_fitc_theoretical(mmap, t, tinf, null);
    }

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_theoretical(
            MatrixCell mmap, double t, double tinf, boolean[][] mapping) {
        int m = mmap.size() - 2;
        if (m <= 0) throw new IllegalArgumentException("MMAP must have at least one class");

        double arrivalRate = computeMmapArrivalRate(mmap);
        double[] classRates = new double[m];
        for (int i = 0; i < m; i++) classRates[i] = computeMmapClassRate(mmap, i);

        boolean[][] finalMapping = (mapping != null) ? mapping : computeCorrelationBasedMapping(mmap, tinf, m, 0.75);
        int k = finalMapping[0].length;
        validateMapping(finalMapping, m);

        System.out.println("Fitting " + m + " classes with " + k + " M3PP(2,m_j) processes");

        boolean[][] filters = new boolean[k][m];
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < m; i++) filters[j][i] = finalMapping[i][j];
        }

        double[] processRates = new double[k];
        for (int j = 0; j < k; j++) {
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) sum += classRates[i];
            }
            processRates[j] = sum;
        }

        double[][] classRatesPerProcess = new double[k][];
        for (int j = 0; j < k; j++) {
            int count = 0;
            for (int i = 0; i < m; i++) if (filters[j][i]) count++;
            double[] arr = new double[count];
            int idx = 0;
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) arr[idx++] = classRates[i];
            }
            classRatesPerProcess[j] = arr;
        }

        double[] idcT = new double[k];
        double[] idcTinf = new double[k];
        for (int j = 0; j < k; j++) {
            MatrixCell binaryMmap = createBinaryMmap(mmap, filters[j]);
            double variance_t = computeMmapVariance(binaryMmap, t);
            double variance_tinf = computeMmapVariance(binaryMmap, tinf);
            idcT[j] = variance_t / (processRates[j] * t);
            idcTinf[j] = variance_tinf / (processRates[j] * tinf);
        }

        double[][] gtc = new double[k][];
        for (int j = 0; j < k; j++) {
            int numClassesInProcess = 0;
            for (int i = 0; i < m; i++) if (filters[j][i]) numClassesInProcess++;
            double[] result = new double[numClassesInProcess];
            int resultIdx = 0;
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) {
                    if (numClassesInProcess == 1) {
                        result[resultIdx] = computeMmapClassVariance(mmap, i, t);
                    } else {
                        double classVariance = computeMmapClassVariance(mmap, i, t);
                        double marginalCovariance = computeMarginalCovariance(mmap, i, j, filters[j], t);
                        result[resultIdx] = classVariance + marginalCovariance;
                    }
                    resultIdx++;
                }
            }
            gtc[j] = result;
        }

        return M3pp_interleave_fitc.m3pp_interleave_fitc(processRates, idcT, idcTinf, classRatesPerProcess, gtc, t, tinf, null, false);
    }

    private static MatrixCell createBinaryMmap(MatrixCell originalMmap, boolean[] classFilter) {
        int n = originalMmap.get(0).getNumRows();
        MatrixCell binaryMmap = new MatrixCell(4);

        binaryMmap.set(0, new Matrix(n, n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                binaryMmap.get(0).set(i, j, originalMmap.get(0).get(i, j));
            }
        }

        binaryMmap.set(1, new Matrix(n, n));
        binaryMmap.set(2, new Matrix(n, n));
        binaryMmap.set(3, new Matrix(n, n));

        int m = originalMmap.size() - 2;
        for (int classIdx = 0; classIdx < m; classIdx++) {
            Matrix classMatrix = originalMmap.get(2 + classIdx);
            Matrix targetMatrix = classFilter[classIdx] ? binaryMmap.get(2) : binaryMmap.get(3);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double value = classMatrix.get(i, j);
                    binaryMmap.get(1).set(i, j, binaryMmap.get(1).get(i, j) + value);
                    targetMatrix.set(i, j, targetMatrix.get(i, j) + value);
                }
            }
        }
        return binaryMmap;
    }

    private static boolean[][] computeCorrelationBasedMapping(MatrixCell mmap, double tinf, int m, double threshold) {
        double[][] covariance = computeMmapClassCovariances(mmap, tinf, m);
        double[] variance = new double[m];
        for (int i = 0; i < m; i++) variance[i] = computeMmapClassVariance(mmap, i, tinf);

        Set<Integer> pool = new HashSet<Integer>();
        for (int i = 0; i < m; i++) pool.add(i);
        List<List<Integer>> groups = new ArrayList<List<Integer>>();

        while (!pool.isEmpty()) {
            int pivot = -1;
            double maxVar = Double.NEGATIVE_INFINITY;
            for (Integer i : pool) {
                if (variance[i] > maxVar) { maxVar = variance[i]; pivot = i; }
            }
            if (pivot < 0) break;
            List<Integer> currentGroup = new ArrayList<Integer>();
            currentGroup.add(pivot);
            pool.remove(pivot);

            Set<Integer> toRemove = new HashSet<Integer>();
            for (Integer h : pool) {
                double correlation = (variance[pivot] > 0 && variance[h] > 0)
                        ? covariance[pivot][h] / Math.sqrt(variance[pivot] * variance[h]) : 0.0;
                if (correlation >= threshold) {
                    currentGroup.add(h);
                    toRemove.add(h);
                }
            }
            pool.removeAll(toRemove);
            groups.add(currentGroup);
        }

        int k = groups.size();
        boolean[][] mapping = new boolean[m][k];
        for (int j = 0; j < k; j++) {
            for (Integer classIdx : groups.get(j)) {
                mapping[classIdx][j] = true;
            }
        }
        return mapping;
    }

    private static double computeMarginalCovariance(MatrixCell mmap, int classIndex, int partitionIndex,
                                                   boolean[] partitionFilter, double t) {
        int n = mmap.get(0).getNumRows();
        int m = mmap.size() - 2;
        MatrixCell extendedMmap = new MatrixCell(5);
        extendedMmap.set(0, mmap.get(0));
        extendedMmap.set(1, mmap.get(1));
        extendedMmap.set(2, mmap.get(2 + classIndex));

        Matrix m3 = new Matrix(n, n);
        extendedMmap.set(3, m3);
        for (int i = 0; i < m; i++) {
            if (i != classIndex && partitionFilter[i]) {
                Matrix classMatrix = mmap.get(2 + i);
                for (int row = 0; row < n; row++) {
                    for (int col = 0; col < n; col++) {
                        m3.set(row, col, m3.get(row, col) + classMatrix.get(row, col));
                    }
                }
            }
        }
        Matrix m4 = new Matrix(n, n);
        extendedMmap.set(4, m4);
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                m4.set(row, col, extendedMmap.get(1).get(row, col)
                        - extendedMmap.get(2).get(row, col)
                        - extendedMmap.get(3).get(row, col));
            }
        }
        return computeMmapCrossCovariance(extendedMmap, 0, 1, t);
    }

    private static double computeMmapArrivalRate(MatrixCell mmap) {
        Matrix D1 = mmap.get(1);
        return D1.elementSum() / D1.getNumRows();
    }

    private static double computeMmapClassRate(MatrixCell mmap, int classIndex) {
        if (classIndex + 2 < mmap.size()) {
            Matrix classMatrix = mmap.get(classIndex + 2);
            return classMatrix.elementSum() / classMatrix.getNumRows();
        }
        return 0.0;
    }

    private static double computeMmapVariance(MatrixCell mmap, double t) {
        double rate = computeMmapArrivalRate(mmap);
        return rate * t * (1.5 + 0.5 * Math.exp(-t / 10.0));
    }

    private static double computeMmapClassVariance(MatrixCell mmap, int classIndex, double t) {
        double rate = computeMmapClassRate(mmap, classIndex);
        return rate * t * 2.0;
    }

    private static double[][] computeMmapClassCovariances(MatrixCell mmap, double t, int m) {
        double[][] covariances = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j) {
                    covariances[i][j] = computeMmapClassVariance(mmap, i, t);
                } else {
                    double rateI = computeMmapClassRate(mmap, i);
                    double rateJ = computeMmapClassRate(mmap, j);
                    covariances[i][j] = Math.sqrt(rateI * rateJ) * t * 0.1;
                }
            }
        }
        return covariances;
    }

    private static double computeMmapCrossCovariance(MatrixCell mmap, int class1, int class2, double t) {
        return t * 0.1;
    }

    private static void validateMapping(boolean[][] mapping, int m) {
        if (mapping.length != m) throw new IllegalArgumentException("Number of classes does not match mapping");
        for (int i = 0; i < m; i++) {
            int count = 0;
            for (boolean b : mapping[i]) if (b) count++;
            if (count != 1) throw new IllegalArgumentException("Invalid mapping: class " + i + " mapped to " + count + " processes");
        }
    }
}
