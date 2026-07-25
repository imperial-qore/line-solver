/**
 * @file Main load-dependent normalizing constant computation with method selection.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.GlobalConstants;
import jline.api.pfqn.nc.*;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;

public final class Pfqn_ncld {
    private Pfqn_ncld() {}

    public static Ret.pfqnNc pfqn_ncld(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
        double lG = Double.NaN;
        double G = Double.NaN;
        String method = options.method;

        Matrix mu_new;
        if ((int) N.elementSum() >= mu.getNumCols()) {
            mu_new = mu.copy();
        } else {
            mu_new = new Matrix(mu.getNumRows(), 0);
            int i = 0;
            while (i < N.elementSum()) {
                Matrix mu_col_i = new Matrix(mu.getNumRows(), 1);
                Matrix.extract(mu, 0, mu.getNumRows(), i, i + 1, mu_col_i, 0, 0);
                mu_new = Matrix.concatColumns(mu_new, mu_col_i, null);
                i++;
            }
        }

        Matrix L_new = new Matrix(L.getNumRows(), 0);
        Matrix N_new = new Matrix(N.getNumRows(), 0);
        Matrix Z_new = new Matrix(Z.getNumRows(), 0);
        for (int i = 0; i < N.length(); i++) {
            if (FastMath.abs(N.get(i)) >= GlobalConstants.FineTol) {
                Matrix L_col_i = Matrix.extractColumn(L, i, null);
                Matrix N_col_i = Matrix.extractColumn(N, i, null);
                Matrix Z_col_i = Matrix.extractColumn(Z, i, null);
                L_new = Matrix.concatColumns(L_new, L_col_i, null);
                N_new = Matrix.concatColumns(N_new, N_col_i, null);
                Z_new = Matrix.concatColumns(Z_new, Z_col_i, null);
            }
        }

        int R = N_new.getNumCols();
        Matrix scalevec = new Matrix(1, R);
        scalevec.fill(1.0);
        for (int r = 0; r < R; r++) {
            Matrix L_col_r = Matrix.extractColumn(L_new, r, null);
            Matrix Z_col_r = Matrix.extractColumn(Z_new, r, null);
            scalevec.set(r, FastMath.max(L_col_r.elementMax(), Z_col_r.elementMax()));
        }

        for (int i = 0; i < L_new.getNumRows(); i++) {
            for (int j = 0; j < L_new.getNumCols(); j++) {
                L_new.set(i, j, L_new.get(i, j) / scalevec.get(j));
            }
        }

        for (int j = 0; j < Z_new.getNumCols(); j++) {
            Z_new.set(j, Z_new.get(j) / scalevec.get(j));
        }

        Matrix Lsum = new Matrix(L_new.getNumRows(), 1);
        Matrix Lmax = new Matrix(L_new.getNumRows(), 1);
        for (int i = 0; i < L_new.getNumRows(); i++) {
            Matrix L_row_i = new Matrix(1, L_new.getNumCols());
            Matrix.extract(L_new, i, i + 1, 0, L_new.getNumCols(), L_row_i, 0, 0);
            Lsum.set(i, L_new.sumRows(i));
            Lmax.set(i, L_row_i.elementMax());
        }

        List<Integer> demStations = new ArrayList<Integer>();
        Matrix L_tmp = new Matrix(0, L_new.getNumCols());
        Matrix mu_tmp = new Matrix(0, mu_new.getNumCols());

        for (int i = 0; i < L_new.getNumRows(); i++) {
            double ratio = Lmax.get(i) / Lsum.get(i);
            if (!Double.isNaN(ratio) && ratio > GlobalConstants.FineTol) {
                demStations.add(i);
                Matrix L_row_i = new Matrix(1, L_new.getNumCols());
                Matrix.extract(L_new, i, i + 1, 0, L_new.getNumCols(), L_row_i, 0, 0);
                L_tmp = Matrix.concatRows(L_tmp, L_row_i, null);
                Matrix mu_row_i = Matrix.extractRows(mu_new, i, i + 1, null);
                mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null);
            }
        }
        L_new = L_tmp.copy();
        mu_new = mu_tmp.copy();

        boolean flag = false;
        for (int i = 0; i < N_new.getNumCols(); i++) {
            if (FastMath.abs(L_new.sumCols(i) + Z_new.sumCols(i)) < GlobalConstants.FineTol
                    && N_new.get(i) > GlobalConstants.FineTol) {
                flag = true;
                break;
            }
        }

        if (flag) {
            System.out.println("pfqn_ncld warning: The model has no positive demands in any class.");
            if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
                lG = 0.0;
            } else {
                Matrix tmp1 = Z_new.sumCols();
                Matrix tmp2 = scalevec.copy();
                for (int i = 0; i < tmp1.length(); i++) {
                    tmp1.set(i, FastMath.log(tmp1.get(i)));
                    tmp2.set(i, FastMath.log(tmp2.get(i)));
                }
                lG = -Matrix.factln(N_new).elementSum()
                        + N_new.elementMult(tmp1, null).elementSum()
                        + N_new.mult(tmp2.transpose()).get(0);
            }
            G = Double.NaN;
            return new Ret.pfqnNc(G, lG, method);
        }

        int M = L_new.getNumRows();
        R = L_new.getNumCols();

        if (L_new.isEmpty() || L_new.elementSum() < options.tol) {
            if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
                lG = 0.0;
            } else {
                Matrix tmp1 = Z_new.sumCols();
                Matrix tmp2 = scalevec.copy();
                for (int i = 0; i < tmp1.length(); i++) {
                    tmp1.set(i, FastMath.log(tmp1.get(i)));
                    tmp2.set(i, FastMath.log(tmp2.get(i)));
                }
                lG = -Matrix.factln(N_new).elementSum()
                        + N_new.elementMult(tmp1, null).elementSum()
                        + N_new.mult(tmp2.transpose()).get(0);
            }
            return new Ret.pfqnNc(G, lG, method);
        } else if (M == 1 && (Z_new.isEmpty() || Z_new.elementSum() < options.tol)) {
            Matrix tmp1 = L_new.sumCols();
            Matrix tmp2 = scalevec.copy();
            for (int i = 0; i < tmp1.length(); i++) {
                tmp1.set(i, FastMath.log(tmp1.get(i)));
                tmp2.set(i, FastMath.log(tmp2.get(i)));
            }
            Matrix tmp3 = mu_new.copy();
            for (int i = 0; i < tmp3.length(); i++) {
                tmp3.set(i, FastMath.log(tmp3.get(i)));
            }
            lG = (Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum()
                    + N_new.elementMult(tmp1, null).elementSum()
                    + N_new.mult(tmp2.transpose()).get(0)) - tmp3.elementSum();
            return new Ret.pfqnNc(G, lG, method);
        }

        List<Integer> zeroDemandClasses = new ArrayList<Integer>();
        List<Integer> nonzeroDemandClasses = new ArrayList<Integer>();
        for (int i = 0; i < R; i++) {
            if (L_new.sumCols(i) < options.tol) {
                zeroDemandClasses.add(i);
            } else {
                nonzeroDemandClasses.add(i);
            }
        }

        double lGzdem;
        Matrix Nz;
        Matrix Zz = new Matrix(Z_new.getNumRows(), 0);
        for (Integer i : zeroDemandClasses) {
            Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
            Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i + 1, Z_col_i, 0, 0);
            Zz = Matrix.concatColumns(Zz, Z_col_i, null);
        }

        flag = true;
        for (int i = 0; i < Zz.getNumCols(); i++) {
            if (Zz.sumCols(i) >= options.tol) { flag = false; break; }
        }
        if (Z_new.isEmpty() || flag) {
            lGzdem = 0.0;
            Nz = new Matrix(1, 1);
            Nz.fill(0.0);
        } else if (zeroDemandClasses.isEmpty()) {
            lGzdem = 0.0;
            Nz = new Matrix(1, 1);
            Nz.fill(0.0);
        } else {
            Nz = new Matrix(1, 0);
            for (Integer i : zeroDemandClasses) {
                Matrix N_col_i = new Matrix(1, 1);
                Matrix.extract(N_new, 0, 1, i, i + 1, N_col_i, 0, 0);
                Nz = Matrix.concatColumns(Nz, N_col_i, null);
            }
            Matrix tmp1 = Zz.sumCols();
            Matrix tmp2 = new Matrix(1, 0);
            for (Integer i : zeroDemandClasses) {
                Matrix scalevec_col_i = new Matrix(1, 1);
                Matrix.extract(scalevec, 0, 1, i, i + 1, scalevec_col_i, 0, 0);
                tmp2 = Matrix.concatColumns(tmp2, scalevec_col_i, null);
            }
            for (int i = 0; i < tmp1.length(); i++) {
                tmp1.set(i, FastMath.log(tmp1.get(i)));
                tmp2.set(i, FastMath.log(tmp2.get(i)));
            }
            lGzdem = -Matrix.factln(Nz).elementSum()
                    + Nz.elementMult(tmp1, null).elementSum()
                    + Nz.mult(tmp2.transpose()).get(0);
        }

        L_tmp = new Matrix(L_new.getNumRows(), 0);
        Matrix N_tmp = new Matrix(1, 0);
        Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
        Matrix scalevecz = new Matrix(1, 0);
        for (Integer i : nonzeroDemandClasses) {
            Matrix L_col_i = new Matrix(L_new.getNumRows(), 1);
            Matrix N_col_i = new Matrix(1, 1);
            Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
            Matrix scalevec_col_i = new Matrix(1, 1);
            Matrix.extract(L_new, 0, L_new.getNumRows(), i, i + 1, L_col_i, 0, 0);
            Matrix.extract(N_new, 0, 1, i, i + 1, N_col_i, 0, 0);
            Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i + 1, Z_col_i, 0, 0);
            Matrix.extract(scalevec, 0, 1, i, i + 1, scalevec_col_i, 0, 0);
            L_tmp = Matrix.concatColumns(L_tmp, L_col_i, null);
            N_tmp = Matrix.concatColumns(N_tmp, N_col_i, null);
            Z_tmp = Matrix.concatColumns(Z_tmp, Z_col_i, null);
            scalevecz = Matrix.concatColumns(scalevecz, scalevec_col_i, null);
        }
        L_new = L_tmp;
        N_new = N_tmp;
        Z_new = Z_tmp;

        double lGnnzdem;
        if (N_new.elementMin() < 0.0) {
            lGnnzdem = 0.0;
        } else {
            Ret.pfqnNc ret = compute_norm_const_ld(L_new, N_new, Z_new, mu_new, options);
            lGnnzdem = ret.lG;
            method = ret.method;
        }

        Matrix tmp = scalevecz.copy();
        for (int i = 0; i < tmp.length(); i++) tmp.set(i, FastMath.log(tmp.get(i)));
        lG = lGnnzdem + lGzdem + N_new.mult(tmp.transpose()).get(0);
        G = FastMath.exp(lG);
        return new Ret.pfqnNc(G, lG, method);
    }

    public static Ret.pfqnNc compute_norm_const_ld(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        String method = options.method;
        Double lG = null;

        if ("default".equals(options.method) || "exact".equals(options.method)) {
            // In "default" mode prefer the Choudhury-Leung-Whitt generating-
            // function inversion (Pfqn_clw_lld) for genuinely multi-station,
            // low-class-count models. Two independent limits, both calibrated by
            // profiling the canonical JAR, gate its use (else fall back to
            // exact/gld):
            //  1) Speed. Cost is exactly prod_j 2*l_j*N_j contour points (l
            //     defaults [1,2,2,3,3,...]) at a steady ~1e7 points/s, so
            //     runtime ~= clwPredCost/1e7 s. It grows as the product of the
            //     per-class populations (degree R in the total population for
            //     balanced classes: ~2*Ntot^2 at R=2, ~Ntot^3 at R=3) and
            //     exponentially in R. Cap at a ~2s budget (CLW_MAX_COST).
            //  2) Numerical validity. The inversion sums 2*N_j alternating-sign
            //     contour points per class, so it loses accuracy to catastrophic
            //     cancellation as the population grows: clw_lld returns NaN
            //     intermittently above per-class ~450 (JAR; native ~700).
            //     Restrict to total population sum(N) <= CLW_MAX_POP = 200, well
            //     inside the reliable regime; larger models revert to exact.
            // The explicit "exact" method always uses the convolution path
            // below; "clw" forces the inversion irrespective of these caps.
            final int CLW_MAX_CLASSES = 5;   // class-count gate (cost ~exp in R)
            final double CLW_MAX_POP = 200;  // total-population cap (numerical; NaN onset ~450)
            final double CLW_MAX_COST = 2e7; // contour-point budget (~2s at ~1e7 pts/s)
            double clwPredCost = 1.0;
            double clwTotPop = N.elementSum();
            for (int j = 0; j < R; j++) {
                clwPredCost *= 2.0 * clwLattice(R, j) * N.get(j);
            }
            if ("default".equals(options.method) && M > 1 && R >= 2
                    && R <= CLW_MAX_CLASSES && clwTotPop <= CLW_MAX_POP
                    && clwPredCost <= CLW_MAX_COST) {
                lG = Pfqn_clw_lld.pfqn_clw_lld(L, N, Z.sumCols(), mu).lG;
                method = "clw";
            } else {
                Matrix Lz, muz;
                if (Z.elementSum() < GlobalConstants.FineTol) {
                    Lz = L; muz = mu;
                } else {
                    int D = Z.getNumRows();
                    Lz = Matrix.concatRows(L, Z, null);
                    Matrix tmp = new Matrix(1, mu.getNumCols());
                    int i = 0;
                    while (i < tmp.length()) {
                        tmp.set(i, i + 1.0);
                        i++;
                    }
                    muz = Matrix.concatRows(mu, tmp.repmat(D, 1), null);
                }
                if (R == 1) {
                    lG = Pfqn_gldsingle.pfqn_gldsingle(Lz, N, muz, options).lG;
                    method = "exact/gld";
                } else if (M == 1 && Z.elementMax() > 0) {
                    // see _kb/03-api-layer.md for rationale
                    Ret.pfqnComomrmLd ret = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, N, Z, mu, options);
                    lG = ret.lG;
                    method = "exact/comomld";
                } else if (M == 1 && Z.elementMax() < GlobalConstants.FineTol) {
                    // see _kb/03-api-layer.md for rationale
                    Matrix zeroZ = new Matrix(N.getNumRows(), N.getNumCols());
                    zeroZ.fill(0.0);
                    Ret.pfqnComomrmLd ret = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, N, zeroZ, mu, options);
                    lG = ret.lG;
                    method = "exact/comomld";
                } else {
                    Ret.pfqnNc ret = Pfqn_gld.pfqn_gld(Lz, N, muz, options);
                    lG = ret.lG;
                    method = "exact/gld";
                }
            }
        } else if ("is".equals(options.method)) {
            // see _kb/03-api-layer.md for rationale
            lG = Pfqn_ld_is.pfqn_ld_is(L, N, Z.sumCols(), mu, options).lG;
            method = "is";
        } else if ("clw".equals(options.method)) {
            // see _kb/03-api-layer.md for rationale
            lG = Pfqn_clw_lld.pfqn_clw_lld(L, N, Z.sumCols(), mu).lG;
            method = "clw";
        } else if ("panacea".equals(options.method) || "panaceald".equals(options.method)) {
            // Mitra-McKenna load-dependent PANACEA asymptotic expansion. Delay
            // terms may arrive either in Z or as mu(i,n)=n rows of L, both are
            // recognized by Pfqn_panaceald.
            lG = Pfqn_panaceald.pfqn_panaceald(L, N, Z.sumCols(), mu).lG;
            method = "panaceald";
            if (Double.isNaN(lG)) {
                // normal usage (1 - lambda_i/mu_i(Ntot) > 0 at every queueing
                // center) is the domain of the expansion, not a numerical failure
                throw new RuntimeException("The model is not in normal usage, so the \"panaceald\" asymptotic expansion does not apply. Use \"exact\", \"clw\" or an approximate load-dependent method instead.");
            }
        } else if ("rd".equals(options.method)) {
            lG = Pfqn_rd.pfqn_rd(L, N, Z, mu, options).lG;
        } else if ("nrp".equals(options.method)) {
            lG = Pfqn_nrp.pfqn_nrp(L, N, Z, mu, options);
        } else if ("nrl".equals(options.method)) {
            lG = Pfqn_nrl.pfqn_nrl(L, N, Z, mu, options);
        } else if ("comomld".equals(options.method)) {
            if (M <= 1 || Z.elementSum() <= GlobalConstants.Zero) {
                lG = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, N, Z, mu, options).lG;
            } else {
                InputOutput.line_warning(InputOutput.mfilename(new Object()),
                        "Load-dependent CoMoM is available only in models with a delay and m identical stations, running the \"rd\" algorithm instead.\n");
                lG = Pfqn_rd.pfqn_rd(L, N, Z, mu, options).lG;
                method = "rd";
            }
        } else {
            throw new RuntimeException("Unrecognized method: " + options.method);
        }

        double G = FastMath.exp(lG);
        return new Ret.pfqnNc(G, lG, method);
    }

    // Default CLW inner lattice parameter l_j (0-based chain index j out of p):
    // l_1=1, l_2=l_3=2, l_j>=4 = 3. Used to predict the clw_lld contour-point cost.
    private static int clwLattice(int p, int j) {
        if (j == 0) return 1;
        if (j == 1 || j == 2) return 2;
        return 3;
    }
}
