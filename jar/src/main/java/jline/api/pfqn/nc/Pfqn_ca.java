/**
 * Convolution Algorithm for Product-Form Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

public final class Pfqn_ca {
    private Pfqn_ca() {}

    public static Ret.pfqnNc pfqn_ca(Matrix L, Matrix N) {
        Matrix Z = N.copy();
        Z.zero();
        return pfqn_ca(L, N, Z);
    }

    public static Ret.pfqnNc pfqn_ca(Matrix L, Matrix N, Matrix Z) {
        Matrix Zlocal = Z;
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (M == 0) {
            Matrix tmp = new Matrix(1, N.length());
            for (int i = 0; i < N.length(); i++) {
                tmp.set(0, i, -Maths.factln(N.get(i)));
            }
            double lGn = tmp.sumRows().sumCols().value();

            Matrix tmp2 = Zlocal.sumCols();
            for (int i = 0; i < tmp2.length(); i++) {
                tmp2.set(i, FastMath.log(tmp2.get(i)));
            }
            if (N.length() == 1) {
                lGn += (N.get(0) * tmp2.sumRows().get(0));
            } else if (tmp2.length() == 1) {
                lGn += (tmp2.get(0) * N.sumRows().sumCols().value());
            } else {
                Matrix tmp3 = new Matrix(1, N.length());
                for (int i = 0; i < N.length(); i++) {
                    tmp3.set(i, N.get(i));
                }
                lGn += tmp3.elementMult(tmp2, null).sumRows().get(0);
            }
            double Gn = FastMath.exp(lGn);
            return new Ret.pfqnNc(Gn, lGn);
        }

        if (N.elementMin() < 0) {
            return new Ret.pfqnNc(0.0, GlobalConstants.NegInf);
        }

        if (N.sumRows().sumCols().get(0) == 0.0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        if (Zlocal.isEmpty()) {
            Matrix temp = new Matrix(1, R);
            temp.fill(0.0);
            Zlocal = temp;
        }

        // see _kb/03-api-layer.md for rationale
        double Nt = N.elementSum();
        Matrix Zsum = Zlocal.sumCols();
        // Each class independently takes whichever station -- or the delay -- gives it
        // its largest factor. The mixed state so named has term at least the product of
        // those factors, because a station holding several classes carries a multinomial
        // coefficient of at least one, so this is still a LOWER bound on log G. It
        // dominates the per-configuration maximum it replaces, which asked ONE station
        // (or the delay) to hold every class at once and so dropped the delay entirely
        // as soon as a single class had no think time. That collapse is what made the
        // scaling scale UP: on L=[1e-9,1], N=[99,1], Z=[1,0] the old estimate was the
        // all-at-the-queue -2051.6 against a true log G of -359.1, giving kscale=-30,
        // and Z/2^-30 = 1.07e9 overflowed the delay column Z^n/n! at n=[40,0].
        double lGest = 0.0;
        for (int r = 0; r < R; r++) {
            if (N.get(r) <= 0) {
                continue;
            }
            double best = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                if (L.get(i, r) > 0) {
                    best = Maths.max(best, N.get(r) * FastMath.log(L.get(i, r)));
                }
            }
            if (Zsum.get(r) > 0) {
                best = Maths.max(best, N.get(r) * FastMath.log(Zsum.get(r)) - Maths.factln(N.get(r)));
            }
            if (Double.isInfinite(best) || Double.isNaN(best)) {
                // no station and no delay can hold class r, so G(N) is exactly zero
                lGest = Double.NEGATIVE_INFINITY;
                break;
            }
            lGest += best;
        }
        int kscale;
        if (Double.isInfinite(lGest) || Double.isNaN(lGest)) {
            kscale = 0;
        } else {
            kscale = (int) Math.round(lGest / (Nt * FastMath.log(2)));
        }
        double cscale = Math.scalb(1.0, kscale);
        Matrix Lscaled = new Matrix(L);
        Lscaled.scaleEq(1.0 / cscale);
        Matrix Zscaled = new Matrix(Zlocal);
        Zscaled.scaleEq(1.0 / cscale);

        int product_N_plus_one = 1;
        for (int i = 0; i < N.length(); i++) {
            product_N_plus_one = (int) (product_N_plus_one * (N.get(i) + 1));
        }
        Matrix G = new Matrix(M + 1, product_N_plus_one);
        G.fill(1.0);
        Matrix n = PopulationLattice.pprod(N);

        while (FastMath.abs(n.sumRows().sumCols().get(0) + 1) > GlobalConstants.FineTol) {
            int idxn = PopulationLattice.hashpop(n, N);
            G.set(0, idxn, Pfqn_pff_delay.pfqn_pff_delay(Zscaled, n));
            for (int m = 1; m < M + 1; m++) {
                G.set(m, idxn, G.get(m - 1, idxn));
                for (int r = 0; r < R; r++) {
                    if (n.get(r) >= 1) {
                        n.set(r, n.get(r) - 1);
                        int idxn_1r = PopulationLattice.hashpop(n, N);
                        n.set(r, n.get(r) + 1);
                        double tmp_res = G.get(m, idxn) + Lscaled.get(m - 1, r) * G.get(m, idxn_1r);
                        G.set(m, idxn, tmp_res);
                    }
                }
            }
            n = PopulationLattice.pprod(n, N);
        }

        // see _kb/03-api-layer.md for rationale
        double rawG = G.get(M, G.getNumCols() - 1);
        double lGn = FastMath.log(rawG) + Nt * kscale * FastMath.log(2);
        double Gn = Math.scalb(rawG, (int) Math.round(Nt * kscale));
        return new Ret.pfqnNc(Gn, lGn);
    }
}
