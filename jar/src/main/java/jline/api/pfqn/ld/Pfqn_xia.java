/**
 * @file Xia's asymptotic approximation for load-dependent normalizing constants
 *
 * Implements Xia's asymptotic approximation method for computing normalizing constants
 * in load-dependent closed queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_xia {
    private Pfqn_xia() {}

    private static double pfqn_xia_F(double u, double k) {
        double ret = 0.0;
        int j = 0;
        while (j < k) {
            ret += FastMath.pow(u, j) / CombinatoricsUtils.factorial(j);
            j++;
        }
        ret += FastMath.pow(u, k) / CombinatoricsUtils.factorial((int) k) / (1 - u / k);
        return ret;
    }

    public static double pfqn_xia(Matrix L, int N, Matrix s, SolverOptions options) {
        Matrix L_local = L.copy();
        Matrix rho = new Matrix(L_local.getNumRows(), L_local.getNumCols());
        for (int i = 0; i < s.getNumRows(); i++) {
            for (int j = 0; j < s.getNumCols(); j++) {
                rho.set(i, j, L_local.get(i, j) / s.get(i, j));
            }
        }
        double scalefactor = 1 / rho.elementMax();
        L_local.scaleEq(scalefactor);
        int M = L_local.getNumRows();
        rho.scaleEq(scalefactor);
        List<Integer> bnkset = new ArrayList<Integer>();
        List<Integer> nbnkset = new ArrayList<Integer>();
        for (int i = 0; i < rho.getNumElements(); i++) {
            if (rho.get(i) == rho.elementMax()) {
                bnkset.add(i);
            }
        }
        for (int i = 0; i < M; i++) {
            if (!bnkset.contains(i)) {
                nbnkset.add(i);
            }
        }
        int B = bnkset.size();
        double logGasy = -Maths.factln(B - 1) - N * FastMath.log(scalefactor);
        for (int b : bnkset) {
            logGasy = logGasy + s.get(b) * FastMath.log(L_local.get(b)) - Maths.factln((int) s.get(b));
        }
        for (int k : nbnkset) {
            double f = pfqn_xia_F(L_local.get(k), s.get(k));
            if (Double.isFinite(f)) {
                logGasy += FastMath.log(f);
            }
        }
        return logGasy;
    }
}
