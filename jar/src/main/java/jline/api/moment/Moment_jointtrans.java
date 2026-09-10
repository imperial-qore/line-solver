package jline.api.moment;

/**
 * Separable joint moment conversion, applied dimension by dimension.
 *
 * <p>Applies the conversion matrix of one edge of the house of moments along
 * every dimension of a joint moment array. This is the Kronecker-product form
 * shared by all the joint conversions except the cumulant and the central ones.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_jointtrans {
    private Moment_jointtrans() {}

    /**
     * Applies one edge of the house of moments along every dimension.
     *
     * @param a flattened joint moment array in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @param edge edge label accepted by Moment_housematrix
     * @return flattened converted array, same layout
     */
    public static double[] moment_jointtrans(double[] a, int[] dims, String edge) {
        double[] out = a;
        for (int mode = 0; mode < dims.length; mode++) {
            out = Moment_tensortrans.moment_tensortrans(out, dims,
                    Moment_housematrix.moment_housematrix(edge, dims[mode] - 1), mode);
        }
        return out;
    }
}
