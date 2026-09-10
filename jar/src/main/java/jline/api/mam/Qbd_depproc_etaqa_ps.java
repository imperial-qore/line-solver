/**
 * @file QBD departure process ETAQA truncation for PS discipline
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;

import jline.lib.smc.QBD_CR;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_depproc_etaqa_ps {
    private Qbd_depproc_etaqa_ps() {}

    /**
     * Compute MAP departure process for MAP/MAP/1-PS via ETAQA truncation.
     */
    public static MatrixCell qbd_depproc_etaqa_ps(MatrixCell MAPa, MatrixCell MAPs, int n) {
        int na = MAPa.get(0).getNumRows();
        int ns = MAPs.get(0).getNumRows();
        int lvlsz = ns * na;

        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);
        Matrix F = MAPa.get(1).kron(IS);
        Matrix L = MAPa.get(0).kron(IS).add(IA.kron(MAPs.get(0)));
        Matrix B = IA.kron(MAPs.get(1));
        Matrix L0 = MAPa.get(0).kron(IS);

        Map<String, Matrix> qbdResult = QBD_CR.QBD_CR(B, L, F, null, null, null, null);
        Matrix R = qbdResult.get("R");
        if (R == null) {
            throw new RuntimeException("QBD_CR failed to compute R matrix");
        }

        Matrix G = L.add(R.mult(B)).neg().inv().mult(B);
        Matrix Lhat = F.add(L);
        Matrix Bbar = B.add(F.mult(G));
        Matrix Bhat = F.mult(G);

        int totalDim = (n + 1) * lvlsz;

        Matrix D0 = Matrix.zeros(totalDim, totalDim);

        D0.insertSubMatrix(0, 0, lvlsz, lvlsz, L0);
        if (n >= 1) {
            D0.insertSubMatrix(0, lvlsz, lvlsz, 2 * lvlsz, F);
        }

        for (int j = 1; j <= n; j++) {
            int rs = j * lvlsz;
            D0.insertSubMatrix(rs, rs, rs + lvlsz, rs + lvlsz, L);
            if (j < n) {
                D0.insertSubMatrix(rs, rs + lvlsz, rs + lvlsz, rs + 2 * lvlsz, F);
            }
        }

        if (n >= 1) {
            int rs = (n - 1) * lvlsz;
            D0.insertSubMatrix(rs, rs, rs + lvlsz, rs + lvlsz, Lhat);
        }

        if (n >= 1) {
            int rsN = n * lvlsz;
            int csN1 = (n - 1) * lvlsz;
            for (int i = 0; i < lvlsz; i++) {
                for (int j = 0; j < lvlsz; j++) {
                    D0.set(rsN + i, csN1 + j, D0.get(rsN + i, csN1 + j) + Bbar.get(i, j));
                }
            }
            double tailScale = (double) (n - 1) / (double) n;
            for (int i = 0; i < lvlsz; i++) {
                for (int j = 0; j < lvlsz; j++) {
                    D0.set(rsN + i, rsN + j, D0.get(rsN + i, rsN + j) + Bhat.get(i, j) * tailScale);
                }
            }
        }

        Matrix D1 = Matrix.zeros(totalDim, totalDim);

        if (n >= 1) {
            int rsN = n * lvlsz;
            int csN1 = (n - 1) * lvlsz;
            D1.insertSubMatrix(rsN, csN1, rsN + lvlsz, csN1 + lvlsz, Bbar);
            D1.insertSubMatrix(rsN, rsN, rsN + lvlsz, rsN + lvlsz, Bhat.scale(1.0 / (double) n));
        }

        for (int j = 1; j < n; j++) {
            int rs = j * lvlsz;
            int cs = (j - 1) * lvlsz;
            double invJ = 1.0 / (double) j;

            double nonDepartScale = 1.0 - invJ;
            for (int i = 0; i < lvlsz; i++) {
                for (int k = 0; k < lvlsz; k++) {
                    D0.set(rs + i, cs + k, D0.get(rs + i, cs + k) + B.get(i, k) * nonDepartScale);
                }
            }

            D1.insertSubMatrix(rs, cs, rs + lvlsz, cs + lvlsz, B.scale(invJ));
        }

        MatrixCell result = new MatrixCell(2);
        result.set(0, D0);
        result.set(1, D1);
        return result;
    }
}
