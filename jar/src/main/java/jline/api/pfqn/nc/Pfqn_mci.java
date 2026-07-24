/**
 * @file Monte Carlo Integration (MCI) methods for normalizing constant computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.Locale;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrixRMaj;

import jline.api.pfqn.mva.Pfqn_bs;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mci {
    private Pfqn_mci() {}

    public static Ret.pfqnNc pfqn_mci(Matrix D, Matrix N, Matrix Z, int I, String variant) {
        int M = D.getNumRows();
        int R = D.getNumCols();

        double lGn;
        if (D.isEmpty() || D.elementSum() < 1e-4) {
            lGn = -N.factln().elementSum() + N.scale(Math.log(Z.elementSum())).elementSum();
            double G = FastMath.exp(lGn);
            return new Ret.pfqnNc(G, lGn);
        }

        Matrix tput = new Matrix(1, R);
        Matrix util;
        Matrix gamma = new Matrix(1, M);

        String variantLower = variant.toLowerCase(Locale.getDefault());
        if (variantLower.equals("imci")) {
            tput = Pfqn_bs.pfqn_bs(D, N, Z).X;
            util = D.mult(tput.transpose());
            for (int i = 0; i < M; i++) {
                gamma.set(i, FastMath.max(0.01, 1 - util.get(i)));
            }
        } else if (variantLower.equals("mci")) {
            tput = Pfqn_bs.pfqn_bs(D, N, Z).X;
            util = D.mult(tput.transpose());
            for (int i = 0; i < M; i++) {
                if (util.get(i) > 0.9) {
                    gamma.set(i, 1 / FastMath.sqrt(N.elementMax()));
                } else {
                    gamma.set(i, 1 - util.get(i));
                }
            }
        } else if (variantLower.equals("rm")) {
            for (int r = 0; r < R; r++) {
                Matrix Dr = D.getColumn(r);
                tput.set(r, N.get(r) / (Dr.elementSum() + Z.get(r) + Dr.elementMax() * (N.elementSum() - 1)));
            }
            util = D.mult(tput.transpose());
            for (int i = 0; i < M; i++) {
                if (util.get(i) > 0.9) {
                    gamma.set(i, 1 / FastMath.sqrt(N.elementMax()));
                } else {
                    gamma.set(i, 1 - util.get(i));
                }
            }
        } else {
            util = D.mult(tput.transpose());
        }

        // see _kb/03-api-layer.md for rationale
        Matrix logfact = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double acc = 0.0;
            for (int n = 1; n <= (int) N.get(r); n++) {
                acc += FastMath.log((double) n);
            }
            logfact.set(r, acc);
        }

        // see _kb/03-api-layer.md for rationale
        Random rand = new Random();
        DMatrixRMaj Vdense = new DMatrixRMaj(I, M);
        for (int j = 0; j < M; j++) {
            double scale = -1.0 / gamma.get(j);
            for (int i = 0; i < I; i++) {
                Vdense.set(i, j, scale * FastMath.log(rand.nextDouble()));
            }
        }
        Matrix V = new Matrix(Vdense);

        // see _kb/03-api-layer.md for rationale
        Matrix VD = V.mult(D);
        int Rc = VD.getNumCols();
        DMatrixRMaj logVDZ = new DMatrixRMaj(I, Rc);
        for (int r = 0; r < Rc; r++) {
            double zr = Z.get(0, r);
            for (int i = 0; i < I; i++) {
                logVDZ.set(i, r, FastMath.log(VD.get(i, r) + zr));
            }
        }

        Matrix ones = Matrix.ones(1, M);
        Matrix lZ = ones.sub(gamma).scale(-1.0).mult(V.transpose());
        lZ.subEq(gamma.log().elementSum());
        lZ.subEq(logfact.elementSum());
        lZ.addEq(N.mult(new Matrix(logVDZ).transpose()));
        double lG = Maths.logmeanexp(lZ);

        return new Ret.pfqnNc(FastMath.exp(lG), lG);
    }
}
