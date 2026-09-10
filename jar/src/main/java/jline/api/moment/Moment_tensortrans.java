package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Mode product of a joint moment array with a conversion matrix.
 *
 * <p>Applies a conversion matrix along one dimension of a joint moment array:
 * every fibre of the array along that dimension is replaced by T times that
 * fibre. Applying it once per dimension realises the Kronecker product of the
 * univariate conversions, which is the structure of every separable edge of the
 * house of moments.
 *
 * <p>Joint moment arrays are held flattened in ROW-MAJOR order, the last
 * dimension varying fastest, with dims[j] = n_j+1 the extent of dimension j.
 * The MATLAB counterpart uses native N-D arrays and therefore a column-major
 * layout; only the linearisation differs, never the value at a multi-index.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_tensortrans {
    private Moment_tensortrans() {}

    /**
     * Applies a conversion matrix along one dimension of a joint moment array.
     *
     * @param a flattened joint moment array in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @param T (dims[mode])x(dims[mode]) conversion matrix
     * @param mode zero-based dimension to transform
     * @return flattened array of the same extent, same layout
     */
    public static double[] moment_tensortrans(double[] a, int[] dims, Matrix T, int mode) {
        int d = dims.length;
        if (mode < 0 || mode >= d) {
            throw new IllegalArgumentException("moment_tensortrans: The mode must be in 0,...,"
                    + (d - 1) + ".");
        }
        int nel = 1;
        for (int i = 0; i < d; i++) {
            nel *= dims[i];
        }
        if (a.length != nel) {
            throw new IllegalArgumentException("moment_tensortrans: The array does not match dims.");
        }
        int n = dims[mode];
        if (T.getNumRows() != n || T.getNumCols() != n) {
            throw new IllegalArgumentException("moment_tensortrans: The matrix T must be " + n
                    + "x" + n + ".");
        }
        int stride = 1;
        for (int i = mode + 1; i < d; i++) {
            stride *= dims[i];
        }
        int outer = nel / (n * stride);
        double[] out = new double[nel];
        for (int o = 0; o < outer; o++) {
            for (int s = 0; s < stride; s++) {
                int base = o * n * stride + s;
                for (int i = 0; i < n; i++) {
                    double acc = 0.0;
                    for (int k = 0; k < n; k++) {
                        acc += T.get(i, k) * a[base + k * stride];
                    }
                    out[base + i * stride] = acc;
                }
            }
        }
        return out;
    }
}
