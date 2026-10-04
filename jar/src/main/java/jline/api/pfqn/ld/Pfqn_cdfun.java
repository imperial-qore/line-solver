/**
 * @file Class-dependent scaling function evaluator for load-dependent queueing systems
 *
 * Provides functionality to evaluate class-dependent scaling functions in load-dependent
 * queueing networks. Calculates scaling factors based on queue-length dependent service
 * rates, supporting state-dependent service mechanisms where service capacity varies
 * with the number of jobs present.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.List;

import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

public final class Pfqn_cdfun {
    private Pfqn_cdfun() {}

    /**
     * Evaluate class-dependent (CD) scaling function for the first class.
     *
     * @param nvec      Per-class queue-length values. The values can be continuous.
     * @param cdscaling CD functions indexed by station index, or null if none
     * @param M         Number of stations
     * @return matrix of scaling factors
     */
    public static Matrix pfqn_cdfun(Matrix nvec,
                                    List<SerializableFunction<Matrix, Matrix>> cdscaling,
                                    int M) {
        return pfqn_cdfun(nvec, cdscaling, M, 0);
    }

    /**
     * Evaluate class-dependent (CD) scaling function.
     *
     * <p>Returns, for every station i, the reciprocal of the class-dependent
     * scaling beta_{i,r}(n_i1, ..., n_iR) evaluated at the per-class population
     * vector nvec(i,:), for class r = classIdx.</p>
     *
     * <p>cdscaling.get(i) is a function of the per-class population vector at
     * station i. It may return either a 1x1 matrix, i.e. a chain-independent
     * scaling beta_i(n) shared by every class (the common case), or a length-R
     * row vector, i.e. the per-class scalings [beta_{i,1}(n), ...,
     * beta_{i,R}(n)], of which element classIdx is taken. The per-class form
     * expresses Sauer's chain-dependent service rates mu_{r,i}(n) (Sauer 1983,
     * "Computational Algorithms for State-Dependent Queueing Networks",
     * eq. (40)), so a single class-dependence mechanism covers both the
     * chain-independent and the chain-specific cases.</p>
     *
     * <p>A null entry denotes a station with no class dependence and is skipped,
     * leaving the neutral factor 1. Such stations are deliberately not filled
     * with a constant function: under Sauer a constant 1 asserts that every
     * class completes at rate 1, which is not load independence.</p>
     *
     * @param nvec      Per-class queue-length values. The values can be continuous.
     * @param cdscaling CD functions indexed by station index, or null if none
     * @param M         Number of stations
     * @param classIdx  Class index selecting beta_{i,r}
     * @return matrix of scaling factors
     */
    public static Matrix pfqn_cdfun(Matrix nvec,
                                    List<SerializableFunction<Matrix, Matrix>> cdscaling,
                                    int M,
                                    int classIdx) {
        Matrix r = new Matrix(M, 1);
        r.fill(1.0);
        if (!(cdscaling == null || cdscaling.isEmpty())) {
            for (int i = 0; i < M; i++) {
                SerializableFunction<Matrix, Matrix> func = (i < cdscaling.size()) ? cdscaling.get(i) : null;
                if (func != null) {
                    Matrix v = func.apply(Matrix.extractRows(nvec, i, i + 1, null));
                    double beta;
                    if (v.length() > 1) {
                        // per-class beta_{i,r}: select the requested class
                        beta = v.get(classIdx);
                    } else {
                        beta = v.get(0);
                    }
                    r.set(i, 0, 1.0 / beta);
                }
            }
        }
        return r;
    }
}
