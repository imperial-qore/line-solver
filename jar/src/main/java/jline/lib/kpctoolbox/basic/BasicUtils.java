package jline.lib.kpctoolbox.basic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

/**
 * Basic utility functions for the KPC-Toolbox.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/basic/
 */
public final class BasicUtils {
    private BasicUtils() {}

    /**
     * Finds the position of the minimum value in a vector.
     */
    public static int minpos(double[] v) {
        if (v.length == 0) {
            throw new IllegalArgumentException("Vector cannot be empty");
        }
        int minIdx = 0;
        double minVal = v[0];
        for (int i = 1; i < v.length; i++) {
            if (v[i] < minVal) {
                minVal = v[i];
                minIdx = i;
            }
        }
        return minIdx;
    }

    /**
     * Finds the positions of the n smallest values in a vector.
     */
    public static int[] minpos(double[] v, int n) {
        if (v.length == 0) {
            throw new IllegalArgumentException("Vector cannot be empty");
        }
        int count = Math.min(n, v.length);

        Integer[] idx = new Integer[v.length];
        for (int i = 0; i < v.length; i++) {
            idx[i] = i;
        }
        final double[] vRef = v;
        Arrays.sort(idx, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(vRef[a], vRef[b]);
            }
        });
        int[] result = new int[count];
        for (int i = 0; i < count; i++) {
            result[i] = idx[i];
        }
        return result;
    }

    /**
     * Finds the position of the maximum value in a vector.
     */
    public static int maxpos(double[] v) {
        if (v.length == 0) {
            throw new IllegalArgumentException("Vector cannot be empty");
        }
        int maxIdx = 0;
        double maxVal = v[0];
        for (int i = 1; i < v.length; i++) {
            if (v[i] > maxVal) {
                maxVal = v[i];
                maxIdx = i;
            }
        }
        return maxIdx;
    }

    /**
     * Finds the positions of the n largest values in a vector.
     */
    public static int[] maxpos(double[] v, int n) {
        if (v.length == 0) {
            throw new IllegalArgumentException("Vector cannot be empty");
        }
        int count = Math.min(n, v.length);

        Integer[] idx = new Integer[v.length];
        for (int i = 0; i < v.length; i++) {
            idx[i] = i;
        }
        final double[] vRef = v;
        Arrays.sort(idx, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(vRef[b], vRef[a]);
            }
        });
        int[] result = new int[count];
        for (int i = 0; i < count; i++) {
            result[i] = idx[i];
        }
        return result;
    }

    /**
     * Generates logarithmically spaced integers in [a, b].
     */
    public static int[] logspacei(double a, double b, int points) {
        if (points <= 0) {
            return new int[0];
        }
        if (points == 1) {
            return new int[]{(int) FastMath.round(a)};
        }

        double logA = FastMath.log10(a);
        double logB = FastMath.log10(b);
        double step = (logB - logA) / (points - 1);

        int[] result = new int[points];
        for (int i = 0; i < points; i++) {
            int value = (int) FastMath.round(FastMath.pow(10.0, logA + i * step));
            if (value < a) value = (int) FastMath.round(a);
            if (value > b) value = (int) FastMath.round(b);
            result[i] = value;
        }
        return result;
    }

    /**
     * Computes the spectral decomposition of a matrix.
     * Returns eigenvalues, eigenvectors, and projectors.
     */
    public static SpectralDecomposition spectd(Matrix A) {
        int n = A.getNumRows();

        RealMatrix realMatrix = MatrixUtils.createRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                realMatrix.setEntry(i, j, A.get(i, j));
            }
        }

        EigenDecomposition eigen = new EigenDecomposition(realMatrix);
        double[] realEigenvalues = eigen.getRealEigenvalues();

        RealMatrix vMatrix = eigen.getV();
        Matrix V = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                V.set(i, j, vMatrix.getEntry(i, j));
            }
        }

        RealMatrix vInverse = MatrixUtils.inverse(vMatrix);
        Matrix iV = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                iV.set(i, j, vInverse.getEntry(i, j));
            }
        }

        Matrix D = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            D.set(i, i, realEigenvalues[i]);
        }

        List<Matrix> projectors = new ArrayList<Matrix>();
        for (int k = 0; k < n; k++) {
            Matrix projector = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    projector.set(i, j, V.get(i, k) * iV.get(k, j));
                }
            }
            projectors.add(projector);
        }

        return new SpectralDecomposition(realEigenvalues, projectors, V, D);
    }

    /**
     * Creates a column vector of ones.
     */
    public static Matrix ones(int n) {
        Matrix result = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            result.set(i, 0, 1.0);
        }
        return result;
    }

    /**
     * Creates a vector e(n) = [1, 1, ..., 1]^T used in matrix operations.
     */
    public static Matrix e(int n) {
        return ones(n);
    }

    /**
     * Creates an identity matrix.
     */
    public static Matrix eye(int n) {
        return Matrix.eye(n);
    }

    /**
     * Creates a zero matrix.
     */
    public static Matrix zeros(int m, int n) {
        return new Matrix(m, n);
    }
}
