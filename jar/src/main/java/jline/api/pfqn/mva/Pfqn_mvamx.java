/**
 * @file Mean Value Analysis for mixed open-closed queueing networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_mvamx {
    private Pfqn_mvamx() {}

    /**
     * Mean Value Analysis (MVA) method for open and mixed queueing networks with no multi-server nodes.
     */
    public static Ret.pfqnMVA pfqn_mvamx(Matrix lambda, Matrix D, Matrix N, Matrix Z, Matrix mi) {
        Matrix Zlocal = (Z == null || Z.isEmpty()) ? new Matrix(1, D.getNumCols()) : Z.copy();
        Matrix miLocal = (mi == null) ? null : mi.copy();
        for (int i = 0; i < lambda.getNumCols(); i++) {
            if (lambda.get(i) > 0 && N.get(i) > 0 && Double.isFinite(N.get(i))) {
                throw new RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.");
            }
        }
        int M = D.getNumRows();
        int R = D.getNumCols();
        if (miLocal == null) {
            miLocal = Matrix.ones(M, 1);
        }
        List<Integer> openClasses = new ArrayList<Integer>();
        List<Integer> closedClasses = new ArrayList<Integer>();
        for (int i = 0; i < N.length(); i++) {
            if (Utils.isInf(N.get(i))) {
                openClasses.add(i);
            } else {
                closedClasses.add(i);
            }
        }
        Matrix XN = new Matrix(1, R);
        Matrix UN = new Matrix(M, R);
        Matrix CN = new Matrix(M, R);
        Matrix QN = new Matrix(M, R);

        for (Integer r : openClasses) {
            for (int ist = 0; ist < M; ist++) {
                UN.set(ist, r, lambda.get(r) * D.get(ist, r));
            }
            XN.set(0, r, lambda.get(r));
        }

        Matrix UNt = UN.sumRows();
        Matrix Dc = new Matrix(D.getNumRows(), closedClasses.size());
        Matrix rep = UNt.repmat(1, closedClasses.size());
        for (int i = 0; i < Dc.getNumRows(); i++) {
            int j = 0;
            for (Integer closedClass : closedClasses) {
                Dc.set(i, j, D.get(i, closedClass) / (1 - rep.get(i, j)));
                j++;
            }
        }
        Matrix Nclosed = new Matrix(1, closedClasses.size());
        Matrix Zclosed = new Matrix(1, closedClasses.size());
        int idx = 0;
        for (Integer closedClass : closedClasses) {
            Nclosed.set(0, idx, N.get(closedClass));
            Zclosed.set(0, idx, Zlocal.get(closedClass));
            idx++;
        }
        Ret.pfqnMVA ret1 = Pfqn_mva.pfqn_mva(Dc, Nclosed, Zclosed, miLocal);
        for (int i = 0; i < closedClasses.size(); i++) {
            int cc = closedClasses.get(i);
            XN.set(cc, ret1.X.get(i));
            for (int j = 0; j < QN.getNumRows(); j++) {
                QN.set(j, cc, ret1.Q.get(j, i));
            }
            for (int j = 0; j < CN.getNumRows(); j++) {
                CN.set(j, cc, ret1.R.get(j, i));
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (Integer r : closedClasses) {
                UN.set(ist, r, XN.get(r) * D.get(ist, r));
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (Integer r : openClasses) {
                if (ret1.Q.isEmpty()) {
                    CN.set(ist, r, D.get(ist, r) / (1 - UNt.get(ist)));
                } else {
                    CN.set(ist, r, D.get(ist, r) * (1 + ret1.Q.sumRows(ist)) / (1 - UNt.get(ist)));
                }
                QN.set(ist, r, CN.get(ist, r) * XN.get(r));
            }
        }
        return new Ret.pfqnMVA(XN, QN, UN, CN, ret1.lGN);
    }

    public static Ret.pfqnMVA pfqn_mvamx(Matrix lambda, Matrix D, Matrix N, Matrix mi) {
        return pfqn_mvamx(lambda, D, N, new Matrix(1, D.getNumCols()), mi);
    }
}
