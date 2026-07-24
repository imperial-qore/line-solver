/**
 * @file Linearizer approximate MVA for mixed networks with multi-server stations
 *
 * Implements linearizer-based approximation methods for mixed queueing networks with
 * multi-server stations. Provides multiple solution techniques including standard linearizer,
 * general-form linearizer, and extended general-form linearizer with automatic method selection.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_linearizermx {
    private Pfqn_linearizermx() {}

    /**
     * Linearizer method for mixed models with multi-server stations.
     */
    public static Ret.pfqnAMVA pfqn_linearizermx(Matrix lambda,
                                                  Matrix L,
                                                  Matrix N,
                                                  Matrix Z,
                                                  Matrix nservers,
                                                  SchedStrategy[] type,
                                                  double tol,
                                                  int maxiter,
                                                  String method) {
        return pfqn_linearizermx(lambda, L, N, Z, nservers, type, tol, maxiter, method, null);
    }

    public static Ret.pfqnAMVA pfqn_linearizermx(Matrix lambda,
                                                  Matrix L,
                                                  Matrix N,
                                                  Matrix Z,
                                                  Matrix nservers,
                                                  SchedStrategy[] type,
                                                  double tol,
                                                  int maxiter,
                                                  String method,
                                                  Matrix QN0) {
        Matrix lambdaLocal = lambda.copy(); // Create local copy to avoid modifying the original
        Matrix LLocal = L.copy();           // Create local copy to avoid modifying the original

        for (int i = 0; i < lambdaLocal.getNumCols(); i++) {
            if (lambdaLocal.get(i) > 0 && N.get(i) > 0 && Double.isFinite(N.get(i))) {
                throw new RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.");
            }
        }
        int M = LLocal.getNumRows();
        int R = LLocal.getNumCols();

        for (int i = 0; i < lambdaLocal.getNumRows(); i++) {
            for (int j = 0; j < lambdaLocal.getNumCols(); j++) {
                if (Double.isNaN(lambdaLocal.get(i, j))) {
                    lambdaLocal.set(i, j, 0);
                }
            }
        }
        for (int i = 0; i < LLocal.getNumRows(); i++) {
            for (int j = 0; j < LLocal.getNumCols(); j++) {
                if (Double.isNaN(LLocal.get(i, j))) {
                    LLocal.set(i, j, 0);
                }
            }
        }
        for (int i = 0; i < Z.getNumRows(); i++) {
            for (int j = 0; j < Z.getNumCols(); j++) {
                if (Double.isNaN(Z.get(i, j))) {
                    Z.set(i, j, 0);
                }
            }
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
        Matrix WN = new Matrix(M, R);
        Matrix QN = new Matrix(M, R);
        Matrix CN = new Matrix(1, R);

        for (int r : openClasses) {
            for (int i = 0; i < M; i++) {
                UN.set(i, r, lambdaLocal.get(r) * LLocal.get(i, r));
            }
            XN.set(0, r, lambdaLocal.get(r));
        }

        Matrix UNt = UN.sumRows();

        if (Z.isEmpty()) {
            Z = new Matrix(1, R);
        } else {
            Z = Z.sumCols();
        }
        Matrix Dc = new Matrix(LLocal.getNumRows(), closedClasses.size());
        Matrix rep = UNt.repmat(1, closedClasses.size());
        for (int i = 0; i < Dc.getNumRows(); i++) {
            int j = 0;
            for (int closedClass : closedClasses) {
                Dc.set(i, j, LLocal.get(i, closedClass) / (1 - rep.get(i, j)));
                j++;
            }
        }

        Matrix Nclosed = new Matrix(1, closedClasses.size());
        Matrix Zclosed = new Matrix(1, closedClasses.size());
        int idx = 0;
        for (int closedClass : closedClasses) {
            Nclosed.set(0, idx, N.get(closedClass));
            Zclosed.set(0, idx, Z.get(closedClass));
            idx++;
        }

        // initial closed-class queue lengths (warm start) aligned to the algorithm layout
        Matrix QN0c = null;
        if (QN0 != null && !QN0.isEmpty()) {
            if (QN0.getNumRows() == M && QN0.getNumCols() == N.getNumCols()) {
                QN0c = new Matrix(M, closedClasses.size());
                int cc = 0;
                for (int closedClass : closedClasses) {
                    for (int i = 0; i < M; i++) {
                        QN0c.set(i, cc, QN0.get(i, closedClass));
                    }
                    cc++;
                }
            } else if (QN0.getNumRows() == M && QN0.getNumCols() == closedClasses.size()) {
                QN0c = QN0;
            }
        }

        Matrix QNc;
        Matrix UNc;
        Matrix WNc;
        Matrix CNc;
        Matrix XNc;
        int totiter;

        if (nservers.elementMax() == 1.0) {
            Ret.pfqnAMVA res;
            if ("lin".equals(method)) {
                res = Pfqn_linearizer.pfqn_linearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, QN0c);
            } else if ("gflin".equals(method)) {
                double linAlpha = 2.0;
                res = Pfqn_gflinearizer.pfqn_gflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, linAlpha, QN0c);
            } else {
                // Extended General Form Linearizer
                Matrix alphaM = new Matrix(1, Nclosed.getNumCols());
                // see _kb/03-api-layer.md for rationale
                for (int i = 0; i < Nclosed.getNumCols(); i++) {
                    alphaM.set(i, 0.6 + 1.4 * FastMath.exp(-8 * FastMath.exp(-0.8 * Nclosed.get(i))));
                }
                res = Pfqn_egflinearizer.pfqn_egflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, alphaM, QN0c);
            }
            QNc = res.Q;
            UNc = res.U;
            WNc = res.R;
            CNc = res.C;
            XNc = res.X;
            totiter = res.totiter;
        } else {
            List<SchedStrategy> typeMatrix = new ArrayList<SchedStrategy>();
            for (int i = 0; i < type.length; i++) {
                typeMatrix.add(type[i]);
            }
            Ret.pfqnAMVAMS res = Pfqn_linearizerms.pfqn_linearizerms(Dc, Nclosed, Zclosed, nservers, typeMatrix, tol, maxiter, QN0c);
            QNc = res.Q;
            UNc = res.U;
            WNc = res.R;
            CNc = res.C;
            XNc = res.X;
            totiter = res.totiter;
        }

        for (int i = 0; i < closedClasses.size(); i++) {
            XN.set(closedClasses.get(i), XNc.get(i));
            for (int j = 0; j < QN.getNumRows(); j++) {
                QN.set(j, closedClasses.get(i), QNc.get(j, i));
            }
            for (int j = 0; j < WN.getNumRows(); j++) {
                WN.set(j, closedClasses.get(i), WNc.get(j, i));
            }
            for (int j = 0; j < UN.getNumRows(); j++) {
                UN.set(j, closedClasses.get(i), UNc.get(j, i));
            }
            CN.set(closedClasses.get(i), CNc.get(i));
        }

        for (int i = 0; i < M; i++) {
            for (int r : closedClasses) {
                UN.set(i, r, XN.get(r) * LLocal.get(i, r));
            }
        }

        for (int i = 0; i < M; i++) {
            for (int r : openClasses) {
                if (QNc.isEmpty()) {
                    WN.set(i, r, LLocal.get(i, r) / (1 - UNt.get(i)));
                } else {
                    WN.set(i, r, LLocal.get(i, r) * (1 + QNc.sumRows(i)) / (1 - UNt.get(i)));
                }
                QN.set(i, r, WN.get(i, r) * XN.get(r));
            }
        }
        for (int r : openClasses) {
            CN.set(r, WN.sumCols(r));
        }
        return new Ret.pfqnAMVA(QN, UN, WN, null, CN, XN, totiter);
    }
}
