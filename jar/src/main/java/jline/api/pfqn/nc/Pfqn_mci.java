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

    /**
     * Monte Carlo integration for the normalizing constant.
     *
     * <p>Variants: {@code mci}, {@code imci} (improved), {@code amci}, {@code lhsmci} and
     * {@code rm} (repairman). {@code amci} and {@code lhsmci} use the {@code imci} tilt
     * and differ only in how the uniforms are drawn.</p>
     *
     * <p>{@code amci} draws ANTITHETIC pairs (u, 1-u). This does NOT reliably reduce
     * variance here: the tilted integrand is not monotone in the exponential draws (the
     * tilt term -(1-gamma)V decreases while the N log(VD+Z) term increases), so the pair
     * correlation is not systematically negative; measured variance ratios against
     * {@code imci} range from 0.54 to 1.6 across models. It is kept because it is the
     * Ross-Wang construction, not because it is the better default. {@code lhsmci}
     * stratifies each coordinate by Latin hypercube sampling, which IS reliably
     * variance-reducing on the same models (ratios 0.0 to 0.48, exact quadrature in the
     * limit of one station) at O(I log I) extra cost.</p>
     *
     * @param D       service demand matrix (M x R)
     * @param N       population vector (1 x R)
     * @param Z       think time vector (1 x R)
     * @param I       number of samples
     * @param variant sampling variant
     * @return the normalizing constant and its logarithm
     */
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
        if (variantLower.equals("imci") || variantLower.equals("amci") || variantLower.equals("lhsmci")) {
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
        if (variantLower.equals("amci")) {
            // Antithetic pairs (u, 1-u), the Ross-Wang construction
            int Ih = (I + 1) / 2;
            for (int j = 0; j < M; j++) {
                double scale = -1.0 / gamma.get(j);
                for (int i = 0; i < Ih; i++) {
                    double u = rand.nextDouble();
                    Vdense.set(i, j, scale * FastMath.log(u));
                    if (Ih + i < I) {
                        Vdense.set(Ih + i, j, scale * FastMath.log(1.0 - u));
                    }
                }
            }
        } else if (variantLower.equals("lhsmci")) {
            // Latin hypercube: one sample per stratum in every coordinate, so no region
            // of the tilted density is over- or under-sampled by chance.
            int[] perm = new int[I];
            for (int j = 0; j < M; j++) {
                double scale = -1.0 / gamma.get(j);
                for (int i = 0; i < I; i++) {
                    perm[i] = i;
                }
                for (int i = I - 1; i > 0; i--) {
                    int k = rand.nextInt(i + 1);
                    int t = perm[i];
                    perm[i] = perm[k];
                    perm[k] = t;
                }
                for (int i = 0; i < I; i++) {
                    double u = (perm[i] + rand.nextDouble()) / I;
                    Vdense.set(i, j, scale * FastMath.log(u));
                }
            }
        } else {
            for (int j = 0; j < M; j++) {
                double scale = -1.0 / gamma.get(j);
                for (int i = 0; i < I; i++) {
                    Vdense.set(i, j, scale * FastMath.log(rand.nextDouble()));
                }
            }
        }
        Matrix V = new Matrix(Vdense);

        // see _kb/03-api-layer.md for rationale
        Matrix VD = V.mult(D);
        int Rc = VD.getNumCols();
        DMatrixRMaj logVDZ = new DMatrixRMaj(I, Rc);
        for (int r = 0; r < Rc; r++) {
            // Class r's think time is the COLUMN SUM of Z, which may carry one
            // row per delay node. The reference is explicit about it
            // (pfqn_mci.m:59 writes sum(Z,1)), so this reproduces it rather
            // than relying on every caller having summed first -- both current
            // callers do pass a single row, so nothing observable changes here.
            double zr = 0.0;
            for (int zi = 0; zi < Z.getNumRows(); zi++) zr += Z.get(zi, r);
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
