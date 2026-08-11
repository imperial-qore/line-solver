/**
 * Normalizing Constant Methods for Product-Form Networks
 *
 * Provides a comprehensive suite of normalizing constant algorithms.
 */
package jline.api.pfqn.nc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.InputOutput;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;

public final class Pfqn_nc {

    private Pfqn_nc() {
    }

    public static Ret.pfqnNcXQ pfqn_nc(Matrix lambda, Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        String method = "exact";
        N.length();

        Double lG = null;
        Matrix X = new Matrix(0, 0);
        Matrix Q = new Matrix(0, 0);
        Matrix lambda_new = lambda.copy();
        Matrix L_new = L.copy();
        Matrix N_new = N.copy();
        Matrix Z_new = Z.copy();

        boolean hasNegative = false;
        for (int i = 0; i < N_new.length(); i++) {
            if (N_new.get(i) < 0) {
                hasNegative = true;
                break;
            }
        }

        if (hasNegative || N_new.isEmpty()) {
            lG = GlobalConstants.NegInf;
            return new Ret.pfqnNcXQ(lG, X, Q, method);
        }

        if (N_new.elementSum() < GlobalConstants.FineTol) {
            lG = 0.0;
            return new Ret.pfqnNcXQ(lG, X, Q, method);
        }

        if (lambda_new.isEmpty()) {
            lambda_new = N_new.copy();
            lambda_new.fill(0.0);
        }

        Matrix Qopen = new Matrix(0, lambda_new.length());
        Matrix Ut = new Matrix(1, L_new.getNumRows());
        double lGopen = 0.0;
        for (int i = 0; i < L_new.getNumRows(); i++) {
            Matrix L_row_i = new Matrix(1, L_new.getNumCols());
            Matrix.extract(L_new, i, i + 1, 0, L_new.getNumCols(), L_row_i, 0, 0);
            Ut.set(i, 1 - lambda_new.mult(L_row_i.transpose()).get(0));
            if (Double.isNaN(Ut.get(i))) {
                Ut.set(i, 0.0);
            }
            for (int j = 0; j < L_row_i.length(); j++) {
                L_new.set(i, j, L_new.get(i, j) / Ut.get(i));
            }
            Matrix.extract(L_new, i, i + 1, 0, L_new.getNumCols(), L_row_i, 0, 0);
            for (int j = 0; j < L_row_i.length(); j++) {
                L_row_i.set(j, L_row_i.get(j) / Ut.get(i));
            }
            Matrix tmp = lambda_new.elementMult(L_row_i, null);
            Qopen = Matrix.concatRows(Qopen, tmp, null);
        }
        Qopen.removeNaN();

        List<Integer> ocl = new ArrayList<Integer>();
        for (int i = 0; i < N_new.length(); i++) {
            if (Utils.isInf(N_new.get(i))) {
                ocl.add(i);
            }
        }

        for (int i = 0; i < N_new.length(); i++) {
            if (Utils.isInf(N_new.get(i))) {
                N_new.set(i, 0.0);
            }
        }

        Matrix L_tmp = new Matrix(L_new.getNumRows(), 0);
        Matrix N_tmp = new Matrix(N_new.getNumRows(), 0);
        Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
        Matrix lambda_tmp = new Matrix(lambda_new.getNumRows(), 0);
        for (int i = 0; i < N_new.length(); i++) {
            if (FastMath.abs(N_new.get(i)) >= GlobalConstants.FineTol) {
                Matrix L_col = Matrix.extractColumn(L_new, i, null);
                Matrix N_col = Matrix.extractColumn(N_new, i, null);
                Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
                Matrix lambda_col = Matrix.extractColumn(lambda_new, i, null);
                L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
                N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
                Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
                lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null);
            }
        }
        L_new = L_tmp;
        N_new = N_tmp;
        Z_new = Z_tmp;

        int R = N_new.length();
        Matrix scalevec = new Matrix(1, R);
        scalevec.fill(1.0);
        for (int r = 0; r < R; r++) {
            Matrix L_col_r = new Matrix(L_new.getNumRows(), 1);
            if (L_new.getNumCols() > 0) {
                Matrix.extractColumn(L_new, r, L_col_r);
            } else {
                L_col_r.fill(0.0);
            }
            Matrix Z_col_r = new Matrix(Z_new.getNumRows(), 1);
            if (Z_new.getNumCols() > 0) {
                Matrix.extractColumn(Z_new, r, Z_col_r);
            } else {
                Z_col_r.fill(0.0);
            }
            scalevec.set(r, FastMath.max(L_col_r.elementMax(), Z_col_r.elementMax()));
        }

        for (int i = 0; i < L_new.getNumRows(); i++) {
            for (int j = 0; j < L_new.getNumCols(); j++) {
                L_new.set(i, j, L_new.get(i, j) / scalevec.get(j));
            }
        }

        for (int i = 0; i < Z_new.getNumCols(); i++) {
            Z_new.set(i, Z_new.get(i) / scalevec.get(i));
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
        List<Integer> noDemStations = new ArrayList<Integer>();
        L_tmp = new Matrix(0, L_new.getNumCols());

        for (int i = 0; i < L_new.getNumRows(); i++) {
            if (!Double.isNaN(Lmax.get(i) / Lsum.get(i)) && Lmax.get(i) / Lsum.get(i) > GlobalConstants.FineTol) {
                demStations.add(i);
                Matrix L_row_i = new Matrix(1, L_new.getNumCols());
                Matrix.extract(L_new, i, i + 1, 0, L_new.getNumCols(), L_row_i, 0, 0);
                L_tmp = Matrix.concatRows(L_tmp, L_row_i, null);
            } else {
                noDemStations.add(i);
            }
        }
        L_new = L_tmp;

        boolean flag = false;
        for (int i = 0; i < N_new.getNumCols(); i++) {
            if (FastMath.abs(L_new.sumCols(i) + Z_new.sumCols(i)) < GlobalConstants.FineTol && N_new.get(i) > GlobalConstants.FineTol) {
                flag = true;
                break;
            }
        }

        if (flag) {
            if (options.verbose != VerboseLevel.SILENT) {
                InputOutput.line_warning(InputOutput.mfilename(new Object() {}),
                        "pfqn_nc warning: The model has no positive demands in any class.");
            }
            if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
                lG = 0.0;
            } else {
                Matrix tmp1 = Z_new.sumCols();
                Matrix tmp2 = scalevec.copy();
                for (int i = 0; i < tmp1.length(); i++) {
                    tmp1.set(i, FastMath.log(tmp1.get(i)));
                    tmp2.set(i, FastMath.log(tmp2.get(i)));
                }
                lG = -Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum()
                        + N_new.mult(tmp2.transpose()).get(0);
            }
            return new Ret.pfqnNcXQ(lG, X, Q, method);
        }

        int M = L_new.getNumRows();
        R = L_new.getNumCols();

        if (L_new.isEmpty() || L_new.elementSum() < options.tol) {
            if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
                lG = lGopen;
            } else {
                Matrix tmp1 = Z_new.sumCols();
                Matrix tmp2 = scalevec.copy();
                for (int i = 0; i < tmp1.length(); i++) {
                    tmp1.set(i, FastMath.log(tmp1.get(i)));
                    tmp2.set(i, FastMath.log(tmp2.get(i)));
                }
                lG = (lGopen - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum() + N_new.mult(
                        tmp2.transpose()).get(0));
            }
            return new Ret.pfqnNcXQ(lG, X, Q, method);
        } else if (M == 1 && (Z_new.isEmpty() || Z_new.elementSum() < options.tol)) {
            Matrix tmp1 = L_new.sumCols();
            Matrix tmp2 = scalevec.copy();
            for (int i = 0; i < tmp1.length(); i++) {
                tmp1.set(i, FastMath.log(tmp1.get(i)));
                tmp2.set(i, FastMath.log(tmp2.get(i)));
            }
            lG = (Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null)
                    .elementSum() + N_new.mult(tmp2.transpose()).get(0));
            return new Ret.pfqnNcXQ(lG, X, Q, method);
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
        for (int i : zeroDemandClasses) {
            Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
            Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i + 1, Z_col_i, 0, 0);
            Zz = Matrix.concatColumns(Zz, Z_col_i, null);
        }

        flag = true;
        for (int i = 0; i < Zz.getNumCols(); i++) {
            if (Zz.sumCols(i) >= options.tol) {
                flag = false;
                break;
            }
        }
        if (Z_new.isEmpty() || flag) {
            lGzdem = 0.0;
            Nz = new Matrix(1, 1);
            Nz.fill(0.0);
        } else {
            if (zeroDemandClasses.isEmpty()) {
                lGzdem = 0.0;
                Nz = new Matrix(1, 1);
                Nz.fill(0.0);
            } else {
                Nz = new Matrix(1, 0);
                for (int i : zeroDemandClasses) {
                    Matrix N_col_i = new Matrix(1, 1);
                    Matrix.extract(N_new, 0, 1, i, i + 1, N_col_i, 0, 0);
                    Nz = Matrix.concatColumns(Nz, N_col_i, null);
                }

                Matrix tmp1 = Zz.sumCols();
                Matrix tmp2 = new Matrix(1, 0);
                for (int i : zeroDemandClasses) {
                    Matrix scalevec_col_i = new Matrix(1, 1);
                    Matrix.extract(scalevec, 0, 1, i, i + 1, scalevec_col_i, 0, 0);
                    tmp2 = Matrix.concatColumns(tmp2, scalevec_col_i, null);
                }
                for (int i = 0; i < tmp1.length(); i++) {
                    tmp1.set(i, FastMath.log(tmp1.get(i)));
                    tmp2.set(i, FastMath.log(tmp2.get(i)));
                }

                lGzdem = (-Matrix.factln(Nz).elementSum() + Nz.elementMult(tmp1, null)
                        .elementSum() + Nz.mult(tmp2.transpose()).get(0));
            }
        }

        L_tmp = new Matrix(L_new.getNumRows(), 0);
        N_tmp = new Matrix(1, 0);
        Z_tmp = new Matrix(Z_new.getNumRows(), 0);
        Matrix scalevecz = new Matrix(1, 0);

        for (int i : nonzeroDemandClasses) {
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

        Ret.pfqnNcXQ ret = compute_norm_const(L_new, N_new, Z_new, options);
        double lGnzdem = ret.lG;
        Matrix Xnnzdem = ret.X;
        Matrix _q = ret.Q;
        method = ret.method;

        if (Xnnzdem.isEmpty()) {
            X = new Matrix(0, 0);
            Q = new Matrix(0, 0);
        }

        Matrix tmp = scalevecz.copy();
        for (int i = 0; i < tmp.length(); i++) {
            tmp.set(i, FastMath.log(tmp.get(i)));
        }
        lG = lGopen + lGnzdem + lGzdem + N_new.mult(tmp.transpose()).get(0);

        return new Ret.pfqnNcXQ(lG, X, Q, method);
    }

    public static Ret.pfqnNcXQ compute_norm_const(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix X = new Matrix(0, 0);
        Matrix Q = new Matrix(0, 0);
        String method = options.method;
        Double lG = null;

        if ("ca".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Z.sumCols());
            lG = ret.lG;
        } else if ("clw".equals(options.method)) {
            // Choudhury-Leung-Whitt generating function inversion: each
            // single-server station is a multiplicity-1 queue, delay is the IS term
            Ret.pfqnNc ret = Pfqn_clw.pfqn_clw(L, N, Z.sumCols());
            lG = ret.lG;
        } else if ("default".equals(options.method) || "adaptive".equals(options.method)) {
            Matrix Z_colSum = Z.sumCols();

            if (M > 1) {
                if (N.elementSum() < 1e3) {
                    double Cmax = M * R * Math.pow(50.0, 3);
                    int maxOrder = Math.min((int) Math.ceil((N.elementSum() - 1) / 2.0), 16);

                    double totCost = 0.0;
                    int order = 0;

                    while (order < maxOrder) {
                        double nextCost = R * Maths.binomialCoeff(M + 2 * (order + 1), M - 1);
                        if (totCost + nextCost <= Cmax) {
                            order += 1;
                            totCost += nextCost;
                        } else {
                            break;
                        }
                    }

                    Ret.pfqnNc ret = Pfqn_cub.pfqn_cub(L, N, Z_colSum, order, GlobalConstants.FineTol);
                    lG = ret.lG;
                    method = "cub";
                } else {
                    Ret.pfqnNc ret = Pfqn_le.pfqn_le(L, N, Z_colSum);
                    lG = ret.lG;
                    method = "le";
                }
            } else {
                if (Z_colSum.getNumCols() == 1 && FastMath.abs(Z_colSum.get(0)) < GlobalConstants.FineTol) {
                    Matrix tmp = L.copy();
                    int i = 0;
                    while (i < tmp.length()) {
                        tmp.set(i, FastMath.log(tmp.get(i)));
                        i++;
                    }
                    lG = -N.mult(tmp.transpose()).get(0);
                    method = "exact";
                } else {
                    if (N.elementSum() < 10000) {
                        Ret.pfqnComomrm ret = Pfqn_comomrm.pfqn_comomrm(L, N, Z_colSum, 1, GlobalConstants.Zero);
                        lG = ret.lG;
                        method = "comom";
                    } else {
                        Ret.pfqnNc ret = Pfqn_le.pfqn_le(L, N, Z_colSum);
                        lG = ret.lG;
                        method = "le";
                    }
                }
            }
        } else if ("is".equals(options.method)) {
            // see _kb/03-api-layer.md for rationale
            Ret.pfqnNc ret = Pfqn_is.pfqn_is(L, N, Z.sumCols(), options);
            lG = ret.lG;
            method = "is";
        } else if ("sampling".equals(options.method)) {
            if (M == 1) {
                Ret.pfqnNc ret = Pfqn_mmsample2.pfqn_mmsample2(L, N, Z.sumCols(), options.samples);
                lG = ret.lG;
                method = "sampling";
            } else if (M > R) {
                Ret.pfqnNc ret = Pfqn_mci.pfqn_mci(L, N, Z.sumCols(), options.samples, "imci");
                lG = ret.lG;
                method = "imci";
            } else {
                Ret.pfqnNc ret = Pfqn_ls.pfqn_ls(L, N, Z.sumCols(), (long) options.samples, (long) options.seed);
                lG = ret.lG;
                method = "ls";
            }
        } else if ("cub".equals(options.method) || "gm".equals(options.method)) {
            int order = (int) FastMath.ceil((N.elementSum() - 1) / 2);
            Ret.pfqnNc ret = Pfqn_cub.pfqn_cub(L, N, Z.sumCols(), order, GlobalConstants.FineTol);
            lG = ret.lG;
        } else if ("kt".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_kt.pfqn_kt(L, N, Z.sumCols());
            lG = ret.lG;
            method = "kt";
        } else if ("mmint2".equals(options.method) || "gleint".equals(options.method)) {
            if (L.getNumRows() > 1) {
                throw new RuntimeException("The " + options.method + " method requires a model with a delay and a single queueing station.");
            } else {
                Ret.pfqnNc ret = Pfqn_mmint2_gausslegendre.pfqn_mmint2_gausslegendre(L, N, Z.sumCols(), null);
                lG = ret.lG;
            }
        } else if ("le".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_le.pfqn_le(L, N, Z.sumCols());
            lG = ret.lG;
        } else if ("ls".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_ls.pfqn_ls(L, N, Z.sumCols(), (long) options.samples, (long) options.seed);
            lG = ret.lG;
        } else if ("mci".equals(options.method) || "imci".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_mci.pfqn_mci(L, N, Z.sumCols(), options.samples, options.method);
            lG = ret.lG;
        } else if ("exact".equals(options.method)) {
            if (M >= R || N.elementSum() > 10 || Z.elementSum() > 0) {
                Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Z.sumCols());
                lG = ret.lG;
                method = "exact/ca";
            } else {
                Ret.pfqnNc ret = Pfqn_recal.pfqn_recal(L, N, Z.sumCols());
                lG = ret.lG;
                method = "exact/recal";
            }
        } else if ("comom".equals(options.method)) {
            if (R > 1) {
                // see _kb/03-api-layer.md for rationale
                if (M > 1) {
                    InputOutput.line_error("pfqn_nc", String.format(
                            "The 'comom' method supports a single queueing station, but this model has %d. "
                                    + "Use 'default' or 'ca' for an exact normalizing constant, or SolverJMT "
                                    + "with method 'jmva.comom'.", M));
                }
                try {
                    Ret.pfqnComomrm ret = Pfqn_comomrm.pfqn_comomrm(L, N, Z, 1, options.tol);
                    lG = ret.lG;
                } catch (Exception e) {
                    e.printStackTrace();
                    lG = Double.NaN;
                }
            } else {
                Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Z.sumCols());
                lG = ret.lG;
                method = "ca";
            }
        } else if ("comomld".equals(options.method)) {
            if (R > 1 && M <= 1) {
                try {
                    Ret.pfqnComomrm ret = Pfqn_comomrm.pfqn_comomrm(L, N, Z, 1, options.tol);
                    lG = ret.lG;
                    method = "comom";
                } catch (Exception e) {
                    e.printStackTrace();
                    Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Z.sumCols());
                    lG = ret.lG;
                    method = "ca";
                }
            } else {
                Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Z.sumCols());
                lG = ret.lG;
                method = "ca";
            }
        } else {
            InputOutput.line_warning("pfqn_nc", "unrecognized method \"%s\"", options.method);
            lG = Double.NaN;
            method = options.method;
        }
        return new Ret.pfqnNcXQ(lG, X, Q, method);
    }
}
