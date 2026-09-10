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

    /**
     * Normalizing constant of a closed or mixed product-form network, with the mean
     * values a few methods produce as a by-product.
     *
     * <p>Equivalent to {@code pfqn_nc(lambda, L, N, Z, options, true)}: the caller is
     * assumed to want X and Q when a method can supply them.</p>
     */
    public static Ret.pfqnNcXQ pfqn_nc(Matrix lambda, Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        return pfqn_nc(lambda, L, N, Z, options, true);
    }

    /**
     * @param wantXQ false when the caller wants lG ALONE. A method that supplies mean
     *        values by simulation ('mcmc') must not pay for a run nobody reads; the
     *        reference passes {@code nargout>=2} here.
     */
    public static Ret.pfqnNcXQ pfqn_nc(Matrix lambda, Matrix L, Matrix N, Matrix Z,
                                       SolverOptions options, boolean wantXQ) {
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

        Ret.pfqnNcXQ ret = compute_norm_const(L_new, N_new, Z_new, options, wantXQ);
        double lGnzdem = ret.lG;
        Matrix Xnnzdem = ret.X;
        Matrix _q = ret.Q;
        method = ret.method;

        if (Xnnzdem.isEmpty()) {
            X = new Matrix(0, 0);
            Q = new Matrix(0, 0);
        } else {
            // A method that produces mean values as a by-product (load concealment, and
            // the default path on many stations) hands them back in the REDUCED,
            // SCALED and REORDERED problem: only the demand-bearing stations, only
            // the nonzero-demand classes, in that order, on demands divided by
            // scalevec. Undo all three, or the caller silently reads a permuted
            // throughput of the wrong magnitude. Mirrors the MATLAB pfqn_nc.
            int Rfull = lambda_new.length();
            int Mfull = L.getNumRows();
            X = new Matrix(1, Rfull);
            X.fill(0.0);
            for (int j = 0; j < nonzeroDemandClasses.size(); j++) {
                int c = nonzeroDemandClasses.get(j);
                X.set(0, c, Xnnzdem.get(j) / scalevec.get(c));
            }
            for (int c : zeroDemandClasses) {
                // the whole class sits in the delay, so X = N / Z
                double zc = Z.getNumCols() > c ? Z.sumCols(c) : 0.0;
                X.set(0, c, zc > options.tol ? N.get(c) / zc : 0.0);
            }
            Q = new Matrix(Mfull, Rfull);
            Q.fill(0.0);
            for (int i = 0; i < demStations.size(); i++) {
                for (int j = 0; j < nonzeroDemandClasses.size(); j++) {
                    // Q is invariant under the per-class demand scaling, since
                    // it only ever sees the product L(i,r)*X(r)
                    Q.set(demStations.get(i), nonzeroDemandClasses.get(j), _q.get(i, j));
                }
            }
            for (int c : ocl) {
                X.set(0, c, lambda_new.get(c));
                for (int i = 0; i < Mfull && i < Qopen.getNumRows(); i++) {
                    Q.set(i, c, Qopen.get(i, c));
                }
            }
        }

        Matrix tmp = scalevecz.copy();
        for (int i = 0; i < tmp.length(); i++) {
            tmp.set(i, FastMath.log(tmp.get(i)));
        }
        lG = lGopen + lGnzdem + lGzdem + N_new.mult(tmp.transpose()).get(0);

        return new Ret.pfqnNcXQ(lG, X, Q, method);
    }

    public static Ret.pfqnNcXQ compute_norm_const(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        return compute_norm_const(L, N, Z, options, true);
    }

    /**
     * Auxiliary routine that computes lG after the initial filtering of L, N and Z.
     *
     * @param wantXQ false when the caller requested lG alone.
     */
    public static Ret.pfqnNcXQ compute_norm_const(Matrix L, Matrix N, Matrix Z,
                                                  SolverOptions options, boolean wantXQ) {
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
        } else if ("ger".equals(options.method)) {
            // Residue closed form of the same generating function "clw" inverts
            // numerically. A class eliminated by residues enters only as a pole ORDER,
            // so its population is free: this is the cheap route when one population
            // dwarfs the others, and the expensive one when the classes are many, since
            // the term count grows as C(S+M-1,M-1) per further elimination. The maxterms
            // cap REFUSES rather than truncating, so an oversized model errors here
            // instead of returning a wrong lG. The solver options are deliberately not
            // forwarded: pfqn_gerasimov's tol is a pole-merging threshold, not the
            // iterative tolerance options.tol carries.
            Ret.pfqnNc ret = Pfqn_gerasimov.pfqn_gerasimov(L, N, Z.sumCols());
            lG = ret.lG;
        } else if ("divdiff".equals(options.method)) {
            // Divided-difference closed form, Casale (SIGMETRICS 2017), Eqs. (15) and
            // (16). Load-independent single-server queues only: a think time needs the
            // integral form of Corollary 3.4, which is not implemented. Unlike the
            // default route below this one keeps Pfqn_explicit's warnings, since a
            // caller that named the method has no fallback.
            if (Z.sumCols().elementSum() > 0) {
                throw new IllegalArgumentException(
                        "pfqn_nc: the 'divdiff' method requires a model without think time, "
                                + "which needs the integral form of Corollary 3.4. Use 'ca' or "
                                + "'default'.");
            }
            Pfqn_explicit.Result ex = Pfqn_explicit.pfqn_explicit(L, N);
            lG = ex.lG;
            method = "divdiff/" + ex.method;
        } else if ("default".equals(options.method) || "adaptive".equals(options.method)) {
            Matrix Z_colSum = Z.sumCols();

            // ONE ESTIMATOR ANSWERS THE WHOLE FAMILY. The divided-difference closed
            // form of Casale (SIGMETRICS 2017) is exact here and was briefly tried first
            // on M>1 && R==1 && sum(Z)==0, but the default route does not serve a single
            // constant: the analyzer differences it at N-e_r for X and at the AUGMENTED
            // shape for Q, one extra class holding one job at station i. That shape has
            // R+1 classes, which the closed form refuses at any sizeable population
            // (the outer sum's cancellation), so it kept the cubature while G(N) turned
            // exact. Mixing the two costs more than either: on mqn_singleserver_ps the
            // closed-form G(N) under cubature numerators left sum_i Q_i at 99.500 of
            // N=100, and the conservation rescale then moved the entire cubature error
            // into X, 0.5% against the 0.06% the cubature ratio carries on its own.
            // 'divdiff' stays a NAMED method, where the caller owns the whole family.
            if (M > 1) {
                int order = -1;
                if (N.elementSum() < 1e3) {
                    double Cmax = M * R * Math.pow(50.0, 3);
                    int maxOrder = Math.min((int) Math.ceil((N.elementSum() - 1) / 2.0), 16);

                    double totCost = 0.0;
                    order = 0;

                    while (order < maxOrder) {
                        double nextCost = R * Maths.binomialCoeff(M + 2 * (order + 1), M - 1);
                        if (totCost + nextCost <= Cmax) {
                            order += 1;
                            totCost += nextCost;
                        } else {
                            break;
                        }
                    }
                }

                // Cmax prices neither the Grundmann-Moeller node count nor the
                // think-time v-integration, so the order is re-priced here and
                // lowered until it fits. Lowering the order keeps the cubature;
                // switching to le instead would hand these models to a Laplace
                // expansion whose mode sits on the simplex boundary when L has
                // near-zero rows, which is the flat-layer case that trips the budget
                while (order > 0 && Pfqn_cub.pfqn_cub_evals(M, order, Z_colSum) > Pfqn_cub.CUB_MAX_EVALS) {
                    order -= 1;
                }
                if (order >= 0) {
                    Ret.pfqnNc ret = Pfqn_cub.pfqn_cub(L, N, Z_colSum, order, GlobalConstants.FineTol);
                    lG = ret.lG;
                    method = "cub";
                } else {
                    // BLE on the default path: strictly better on lG and it
                    // cancels in G(N-e_r)/G(N). "le" stays the published form.
                    Ret.pfqnNc ret = Pfqn_ble.pfqn_ble(L, N, Z_colSum);
                    lG = ret.lG;
                    method = "ble";
                    // Birman-Kogan Algorithm 2 supplies the MEAN VALUES here.
                    // The caller's fallback differences lG at R+M*R reduced
                    // populations, which on many stations is both dearer and
                    // ~300x less accurate than the load concealment fixed point. Gated
                    // on the station count, since load concealment is mean
                    // field in M: see _kb/06-solver-catalog.md.
                    if (M >= 10 && R > 1 && N.elementMin() >= 0 && N.elementSum() > 0) {
                        Ret.pfqnBkLc thin =
                                Pfqn_bk.pfqn_bklc(L, N, Z_colSum, "mva", 1e-10, options.iter_max);
                        X = thin.X;
                        Q = thin.Q;
                        method = "ble/lc";
                    }
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
                        Ret.pfqnNc ret = Pfqn_ble.pfqn_ble(L, N, Z_colSum);
                        lG = ret.lG;
                        method = "ble";
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
        } else if ("bkt".equals(options.method)) {
            // KT minus the exact Stirling remainder of each Laplaced class; see _kb/03-api-layer.md
            Ret.pfqnNc ret = Pfqn_bkt.pfqn_bkt(L, N, Z.sumCols());
            lG = ret.lG;
            method = "bkt";
        } else if ("lekt".equals(options.method)) {
            // the estimator ble and bkt both compute, on the cheaper side; see _kb/03-api-layer.md
            Ret.pfqnNc ret = Pfqn_lekt.pfqn_lekt(L, N, Z.sumCols());
            lG = ret.lG;
            method = "lekt";
        } else if ("bk".equals(options.method)) {
            Ret.pfqnNc ret = Pfqn_bk.pfqn_bk(L, N, Z.sumCols());
            lG = ret.lG;
            method = "bk";
        } else if ("bkue".equals(options.method)) {
            // The uniform expansion is single chain by construction; the
            // multichain fallback is the saddle point of the same paper, which is
            // also how the analyzer reaches this branch, since it conditions on a
            // station population by augmenting the model with an auxiliary class.
            if (R > 1) {
                Ret.pfqnNc ret = Pfqn_bk.pfqn_bk(L, N, Z.sumCols());
                lG = ret.lG;
                method = "bkue/bk";
            } else {
                Matrix Zs = Z.sumCols();
                Ret.pfqnNc ret = Pfqn_bk.pfqn_bkue(L, N.get(0, 0), Zs.isEmpty() ? 0.0 : Zs.get(0, 0));
                lG = ret.lG;
                method = "bkue";
            }
        } else if ("lc".equals(options.method) || "lc.ue".equals(options.method)) {
            // the fixed point converges linearly and slowly, so a solver-level
            // reporting tolerance would stop it far from its own limit and at a
            // different sweep in each codebase: iterate to the method's accuracy
            Ret.pfqnBkLc thin = Pfqn_bk.pfqn_bklc(L, N, Z.sumCols(),
                    "lc.ue".equals(options.method) ? "ue" : "mva", 1e-10, options.iter_max);
            X = thin.X;
            Q = thin.Q;
            // Algorithm 2 returns mean values, not a multichain constant; the
            // saddle point that seeds it supplies lG on the same asymptotics
            lG = Pfqn_bk.pfqn_bk(L, N, Z.sumCols()).lG;
            method = options.method;
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
        } else if ("ble".equals(options.method)) {
            // LE plus the empirical eps->0 correction; see _kb/03-api-layer.md
            Ret.pfqnNc ret = Pfqn_ble.pfqn_ble(L, N, Z.sumCols());
            lG = ret.lG;
        } else if ("aghq".equals(options.method)) {
            // adaptive Gauss-Hermite over the simplex; q=1 would be "le".
            // options.config.aghq_nodes overrides the node count.
            int aghqNodes = 3;
            if (options.config != null && options.config.aghq_nodes != null) {
                aghqNodes = Math.max(1, options.config.aghq_nodes.intValue());
            }
            Ret.pfqnNc ret = Pfqn_aghq.pfqn_aghq(L, N, Z.sumCols(), aghqNodes);
            lG = ret.lG;
        } else if ("mcmc".equals(options.method)) {
            // Chen-O'Cinneide REGULARIZATION (TOMACS 8(3), 1998). The chain is simulated
            // on the regularized network, which shares the steady-state distribution of
            // the original one, so it returns X and Q directly through the pfqn_nc X/Q
            // channel, like 'lc'. What it does NOT return is the constant itself:
            // the algorithm estimates the RATIOS G(N-e_r)/G(N), never G, so lG here is
            // the BLE expansion and is not part of the paper. It cancels out of every
            // mean value reported by the analyzer; only getProbNormConstAggr reads it.
            if (wantXQ) {
                Ret.pfqnMcmc mc = Pfqn_mcmc.pfqn_mcmc(L, N, Z.sumCols(), options);
                X = mc.X;
                Q = mc.Q;
            }
            if (M > 1) {
                Ret.pfqnNc ret = Pfqn_ble.pfqn_ble(L, N, Z.sumCols());
                lG = ret.lG;
            } else {
                Ret.pfqnComomrm ret = Pfqn_comomrm.pfqn_comomrm(L, N, Z, 1, options.tol);
                lG = ret.lG;
            }
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
        } else if ("rgf".equals(options.method)) {
            // Recursion by generating functions. Single class: one sequence per
            // group of identically loaded stations (Coury-Harrison 1997, Property
            // 1). Multiclass: the residue recursion of Harrison-Coury 2002, Thm 1,
            // with think times carried by the Bertozzi-McKenna truncation that
            // neither RGF paper has. That sum is ALTERNATING, so Pfqn_rgfmc refuses
            // when the cancellation leaves no significant digits rather than
            // returning a wrong lG; the exact convolution answers those and the
            // reported method says so.
            if (R > 1) {
                Matrix Zs = Z.sumCols();
                double[] Nv = new double[R];
                double[] Zv = new double[R];
                for (int r = 0; r < R; r++) {
                    Nv[r] = N.get(0, r);
                    Zv[r] = Zs.isEmpty() ? 0.0 : Zs.get(0, r);
                }
                try {
                    lG = Pfqn_rgfmc.pfqn_rgfmc(L, Nv, Zv).lG;
                } catch (RuntimeException e) {
                    Ret.pfqnNc ret = Pfqn_ca.pfqn_ca(L, N, Zs);
                    lG = ret.lG;
                    method = "rgf/ca";
                }
            } else {
                Matrix Zsum = Z.sumCols();
                double Ztot = Zsum.isEmpty() ? 0.0 : Zsum.get(0, 0);
                lG = Pfqn_rgf.pfqn_rgf(L, N.get(0, 0), Ztot).lG;
            }
        } else {
            InputOutput.line_warning("pfqn_nc", "unrecognized method \"%s\"", options.method);
            lG = Double.NaN;
            method = options.method;
        }
        return new Ret.pfqnNcXQ(lG, X, Q, method);
    }
}
