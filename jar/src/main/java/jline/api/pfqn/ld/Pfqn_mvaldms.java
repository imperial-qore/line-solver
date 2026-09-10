/**
 * @file Load-dependent MVA wrapper for multi-server utilization adjustment
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_mvaldms {
    private Pfqn_mvaldms() {}

    /**
     * Wrapper for pfqn_mvaldmx that adjusts utilizations to account for multiservers.
     */
    public static Ret.pfqnMVA pfqn_mvaldms(Matrix lambda, Matrix D, Matrix N, Matrix Z, Matrix S) {
        int M = D.getNumRows();
        int R = D.getNumCols();
        double Nct = 0.0;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                double num = N.get(i, j);
                if (Double.isFinite(num)) {
                    Nct += num;
                }
            }
        }
        Matrix mu = Matrix.ones(M, (int) Nct);
        for (int ist = 0; ist < M; ist++) {
            for (int j = 0; j < mu.getNumCols(); j++) {
                mu.set(ist, j, Maths.min((double) (j + 1), S.get(ist)));
            }
        }
        if (Z.isEmpty()) {
            Z = new Matrix(1, R);
        }
        Ret.pfqnMVALDMX ret1 = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, D, N, Z, mu, S);
        List<Integer> openClasses = new ArrayList<Integer>();
        List<Integer> closedClasses = new ArrayList<Integer>();
        for (int i = 0; i < N.length(); i++) {
            if (Utils.isInf(N.get(i))) {
                openClasses.add(i);
            } else {
                closedClasses.add(i);
            }
        }
        Matrix XN = ret1.X;
        Matrix QN = ret1.Q;
        Matrix CN = ret1.R;
        double lGN = ret1.lG;
        Matrix UN = new Matrix(M, R);
        for (Integer r : closedClasses) {
            for (int ist = 0; ist < M; ist++) {
                UN.set(ist, r, XN.get(r) * D.get(ist, r) / S.get(ist));
            }
        }
        for (Integer r : openClasses) {
            for (int ist = 0; ist < M; ist++) {
                UN.set(ist, r, lambda.get(r) * D.get(ist, r) / S.get(ist));
            }
        }
        return new Ret.pfqnMVA(XN, QN, UN, CN, lGN);
    }
}
