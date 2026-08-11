/**
 * @file Joint-dependent (non-product-form) scaling function evaluator
 *
 * Provides functionality to evaluate joint-dependent scaling functions eta_i(n)
 * in queueing networks. Unlike {@link Pfqn_cdfun} (product-form beta_{i,r}
 * depending on the own-class marginal n_{i,r}), eta may read the joint per-class
 * population vector arbitrarily and is therefore non-product-form: the analysis
 * is an approximation with no exactness/uniqueness guarantee. The numerical
 * evaluation matches {@link Pfqn_cdfun}; the distinction is semantic and is
 * carried by the separate sn.jdscaling field.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.List;

import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

public final class Pfqn_jdfun {
    private Pfqn_jdfun() {}

    /**
     * Evaluate joint-dependent (JD) scaling function for the first class.
     *
     * @param nvec      Per-class queue-length values. The values can be continuous.
     * @param jdscaling JD functions indexed by station index, or null if none
     * @param M         Number of stations
     * @return matrix of scaling factors
     */
    public static Matrix pfqn_jdfun(Matrix nvec,
                                    List<SerializableFunction<Matrix, Matrix>> jdscaling,
                                    int M) {
        return pfqn_jdfun(nvec, jdscaling, M, 0);
    }

    /**
     * Evaluate joint-dependent (JD) scaling function.
     *
     * <p>Returns, for every station i, the reciprocal of the joint-dependent
     * scaling eta_i(n_i1, ..., n_iR) evaluated at the per-class population vector
     * nvec(i,:), for class r = classIdx.</p>
     *
     * <p>jdscaling.get(i) is a function of the joint per-class population vector
     * at station i. It may return either a 1x1 matrix (a scaling shared by every
     * class, as in the flagship min(ni[0],c)), or a length-R row vector of
     * per-class scalings, of which element classIdx is taken (Sauer
     * chain-dependent rate mu_{r,i}(n)). Unlike Pfqn_cdfun, eta may read the
     * joint vector arbitrarily and is non-product-form.</p>
     *
     * <p>A null entry denotes a station with no joint dependence and is skipped,
     * leaving the neutral factor 1.</p>
     *
     * @param nvec      Per-class queue-length values. The values can be continuous.
     * @param jdscaling JD functions indexed by station index, or null if none
     * @param M         Number of stations
     * @param classIdx  Class index selecting eta_{i,r}
     * @return matrix of scaling factors
     */
    public static Matrix pfqn_jdfun(Matrix nvec,
                                    List<SerializableFunction<Matrix, Matrix>> jdscaling,
                                    int M,
                                    int classIdx) {
        Matrix r = new Matrix(M, 1);
        r.fill(1.0);
        if (!(jdscaling == null || jdscaling.isEmpty())) {
            for (int i = 0; i < M; i++) {
                SerializableFunction<Matrix, Matrix> func = (i < jdscaling.size()) ? jdscaling.get(i) : null;
                if (func != null) {
                    Matrix v = func.apply(Matrix.extractRows(nvec, i, i + 1, null));
                    double eta;
                    if (v.length() > 1) {
                        // per-class eta_{i,r}: select the requested class
                        eta = v.get(classIdx);
                    } else {
                        eta = v.get(0);
                    }
                    r.set(i, 0, 1.0 / eta);
                }
            }
        }
        return r;
    }
}
