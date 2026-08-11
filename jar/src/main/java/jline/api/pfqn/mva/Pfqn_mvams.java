/**
 * @file General-purpose MVA for mixed networks with multi-server stations
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.api.pfqn.ld.Pfqn_mvald;
import jline.api.pfqn.ld.Pfqn_mvaldms;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_mvams {
    private Pfqn_mvams() {}

    /**
     * General purpose script to handle mixed Query Networks with multiserver nodes.
     */
    public static Ret.pfqnMVA pfqn_mvams(Matrix lambda, Matrix L, Matrix N, Matrix Z, Matrix mi, Matrix S) {
        Matrix Zlocal = (Z == null || Z.isEmpty()) ? new Matrix(1, L.getNumCols()) : Z.copy();
        if (Zlocal.getNumRows() > 1) {
            Zlocal = Zlocal.sumCols();
        }
        Matrix miLocal = (mi == null) ? null : mi.copy();
        Matrix Slocal = (S == null) ? null : S.copy();
        int M = L.getNumRows();
        double Ntot = 0.0;
        boolean NhasInf = false;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                double num = N.get(i, j);
                if (Double.isFinite(num)) {
                    Ntot += num;
                } else if (Utils.isInf(num)) {
                    NhasInf = true;
                }
            }
        }
        Matrix mu = Matrix.ones(M, (int) Ntot);
        if (Slocal == null) {
            Slocal = Matrix.ones(M, 1);
        }
        if (miLocal == null) {
            miLocal = Matrix.ones(M, 1);
        }
        for (int ist = 0; ist < M; ist++) {
            int j = 0;
            while (j < Ntot) {
                mu.set(ist, j, Maths.min((double) (j + 1), Slocal.get(ist)));
                j++;
            }
        }
        // see _kb/03-api-layer.md for rationale
        boolean hasMultiServer = false;
        for (int i = 0; i < Slocal.getNumRows() && !hasMultiServer; i++) {
            for (int j = 0; j < Slocal.getNumCols(); j++) {
                double num = Slocal.get(i, j);
                if (Double.isFinite(num) && num > 1.0) {
                    hasMultiServer = true;
                    break;
                }
            }
        }
        Ret.pfqnMVA returnObject;
        if (!hasMultiServer) {
            if (NhasInf) {
                Ret.pfqnMVA retMVAMX = Pfqn_mvamx.pfqn_mvamx(lambda, L, N, Zlocal, miLocal);
                returnObject = new Ret.pfqnMVA(retMVAMX.X, retMVAMX.Q, retMVAMX.U, retMVAMX.R, retMVAMX.lGN);
            } else {
                Ret.pfqnMVA retMVA = Pfqn_mva.pfqn_mva(L, N, Zlocal, miLocal);
                returnObject = new Ret.pfqnMVA(retMVA.X, retMVA.Q, retMVA.U, retMVA.R, retMVA.lGN);
            }
        } else {
            if (NhasInf) {
                if (miLocal.elementMax() == 1.0) {
                    double lG = Double.NaN;
                    Ret.pfqnMVA retMVALDMS = Pfqn_mvaldms.pfqn_mvaldms(lambda, L, N, Zlocal, Slocal);
                    returnObject = new Ret.pfqnMVA(retMVALDMS.X, retMVALDMS.Q, retMVALDMS.U, retMVALDMS.R, lG);
                } else {
                    throw new RuntimeException("pfqn_mvams: Queue replicas not available in exact MVA for mixed models.");
                }
            } else {
                Ret.pfqnMVALD retMVALD = Pfqn_mvald.pfqn_mvald(L, N, Zlocal, mu);
                double lG = retMVALD.lG.get(retMVALD.lG.size() - 1);
                // see _kb/03-api-layer.md for rationale
                int Mq = L.getNumRows();
                int Rq = L.getNumCols();
                Matrix CN = new Matrix(Mq, Rq);
                for (int r = 0; r < Rq; r++) {
                    double Nr = (N.getNumRows() == 1) ? N.get(0, r) : N.get(r, 0);
                    for (int i = 0; i < Mq; i++) {
                        if (Nr > 0) {
                            CN.set(i, r, retMVALD.Q.get(i, r) / retMVALD.X.get(0, r));
                        } else {
                            // an absent class, as in Pfqn_mva
                            CN.set(i, r, L.get(i, r) * miLocal.get(i));
                        }
                    }
                }
                returnObject = new Ret.pfqnMVA(retMVALD.X, retMVALD.Q, retMVALD.U, CN, lG);
            }
        }
        return returnObject;
    }

    public static Ret.pfqnMVA pfqn_mvams(Matrix lambda, Matrix L, Matrix N, Matrix mi, Matrix S) {
        return pfqn_mvams(lambda, L, N, new Matrix(1, L.getNumCols()), mi, S);
    }
}
