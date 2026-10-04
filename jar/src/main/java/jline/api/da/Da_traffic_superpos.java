/**
 * @file Asymptotic-method superposition of independent flows
 *
 * @since LINE 3.0
 */
package jline.api.da;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Asymptotic-method superposition of independent flows with given rates and
 * squared coefficients of variation: returns the rate-weighted SCV mixture of
 * the merged flow (Whitt's QNA stationary-interval formula). Entries with
 * non-finite rates are ignored.
 *
 * Mirrors matlab/src/api/da/da_traffic_superpos.m.
 */
public final class Da_traffic_superpos {
    private Da_traffic_superpos() {}

    public static double da_traffic_superpos(Matrix lambda, Matrix a2) {
        List<Integer> lambda_finite_idx = new ArrayList<Integer>();
        for (int i = 0; i < lambda.length(); i++) {
            if (Double.isFinite(lambda.get(i))) {
                lambda_finite_idx.add(Integer.valueOf(i));
            }
        }
        Matrix a2_new = new Matrix(1, lambda_finite_idx.size(), lambda_finite_idx.size());
        for (int i = 0; i < lambda_finite_idx.size(); i++) {
            a2_new.set(i, a2.get(lambda_finite_idx.get(i).intValue()));
        }
        Matrix lambda_new = new Matrix(1, lambda_finite_idx.size(), lambda_finite_idx.size());
        for (int i = 0; i < lambda_finite_idx.size(); i++) {
            lambda_new.set(i, lambda.get(lambda_finite_idx.get(i).intValue()));
        }
        return a2_new.mult(lambda_new.transpose()).get(0) / lambda_new.elementSum();
    }
}
