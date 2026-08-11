package jline.solvers.nc.handlers;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.npfqn.Npfqn_nonexp_approx;
import jline.api.pfqn.nc.Pfqn_nc;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Solver_nc {
    private Solver_nc() {}

    public static SolverNC.SolverNCReturn solver_nc(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        Matrix nservers = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        java.util.Map<jline.lang.nodes.Station, SchedStrategy> sched = sn.sched;
        int C = sn.nchains;
        Matrix SCV = sn.scv;
        int K = sn.nclasses;
        long startTime = System.nanoTime();

        // Order-independent (OI) closed network: auto-detected when the closed
        // model consists solely of OI and delay stations. Solved exactly by the
        // balanced-fairness normalizing constant (pfqn_ncoi) with the OI
        // functional-server (pfqn_oi_fnc) identity for the mean queue lengths.
        // Intercept BEFORE the generic normalizing-constant path (e.g. the
        // comomld 2-station branch). There is no explicit 'oi' selector: the
        // OI path is reached only from 'default'/'exact', matching SolverMVA.
        boolean oiAuto = Solver_nc_oi.nc_is_oi_model(sn)
                && options.method != null
                && (options.method.equalsIgnoreCase("default") || options.method.equalsIgnoreCase("exact"));
        if (oiAuto) {
            return Solver_nc_oi.solver_nc_oi(sn, options);
        }

        // Importance sampling ('is') specialized to OI / pass-and-swap (P&S)
        // stations: the auto-normalized IS normalizing constant Pfqn_pas_is
        // (which reduces to the OI case when the swap graph is empty). Triggers
        // on the explicit IS selector 'is'; on 'sampling', which maps to 'is' in
        // the presence of OI/PAS stations rather than the generic Pfqn_mci /
        // Pfqn_ls estimators; and on 'default' for a P&S tandem with a non-empty
        // swap graph, which is reducible and has no exact path (a pure-OI tandem
        // on 'default'/'exact' is caught by the exact OI analyzer above).
        boolean pasAuto = Solver_nc_pas_is.nc_is_pas_model(sn)
                && options.method != null
                && (options.method.equalsIgnoreCase("default")
                    || options.method.equalsIgnoreCase("is")
                    || options.method.equalsIgnoreCase("sampling"));
        if (pasAuto) {
            return Solver_nc_pas_is.solver_nc_pas_is(sn, options);
        }

        // 'is' is the sample-an-ordering importance-sampling family. OI /
        // pass-and-swap stations are handled above (Pfqn_pas_is / Pfqn_oi_is);
        // every other CLOSED product-form network falls through to the standard
        // normalizing-constant path, where Pfqn_nc routes 'is' to Pfqn_is
        // (load-independent) and Pfqn_ncld routes it to Pfqn_ld_is
        // (load-dependent). The estimator has no open-class form.
        if (options.method != null && options.method.equalsIgnoreCase("is")) {
            boolean hasOpen = false;
            for (int r = 0; r < sn.njobs.getNumElements(); r++) {
                if (Double.isInfinite(sn.njobs.get(r))) {
                    hasOpen = true;
                    break;
                }
            }
            if (hasOpen) {
                throw new RuntimeException("The 'is' importance-sampling method requires a closed queueing network. Use 'sampling' (pfqn_mci/pfqn_ls) for open or mixed models.");
            }
        }

        List<Integer> lcfsStats = new ArrayList<Integer>();
        List<Integer> lcfsprStats = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sched.get(sn.stations.get(i));
            if (s == SchedStrategy.LCFS) lcfsStats.add(i);
            else if (s == SchedStrategy.LCFSPR) lcfsprStats.add(i);
        }
        if (!lcfsStats.isEmpty() && !lcfsprStats.isEmpty()) {
            if (lcfsStats.size() != 1 || lcfsprStats.size() != 1) {
                throw new RuntimeException("LCFS NC requires exactly one LCFS and one LCFS-PR station.");
            }
            for (int i = 0; i < NK.getNumRows(); i++) {
                if (Utils.isInf(NK.get(i))) throw new RuntimeException("LCFS NC requires a closed queueing network.");
            }
            return Solver_nc_lcfsqn.solver_nc_lcfsqn(sn, options, lcfsStats.get(0), lcfsprStats.get(0));
        } else if (!lcfsStats.isEmpty()) {
            throw new RuntimeException("LCFS scheduling requires a paired LCFS-PR station.");
        }

        Matrix V = new Matrix(sn.nstateful, K);
        for (int i = 0; i < sn.visits.size(); i++) V = V.add(1.0, sn.visits.get(i));
        Matrix rates = sn.rates;
        Matrix ST = rates.copy();
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                double rate = ST.get(i, j);
                if (Double.isNaN(rate)) ST.set(i, j, 0.0);
                else ST.set(i, j, 1.0 / rate);
            }
        }
        Matrix ST0 = ST.copy();
        Matrix Nchain = new Matrix(1, C);
        Nchain.fill(0.0);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            double sum = 0.0;
            for (int col = 0; col < inchain.getNumCols(); col++) sum += NK.get((int) inchain.get(0, col));
            Nchain.set(c, sum);
        }

        List<Integer> openChains = new ArrayList<Integer>();
        List<Integer> closedChains = new ArrayList<Integer>();
        for (int c = 0; c < C; c++) {
            if (Utils.isInf(Nchain.get(c))) openChains.add(c);
            else closedChains.add(c);
        }

        Matrix gamma = new Matrix(1, M);
        gamma.fill(0.0);
        Matrix eta_1 = new Matrix(1, M);
        eta_1.fill(0.0);
        Matrix eta = new Matrix(1, M);
        eta.fill(1.0);

        if (!sched.containsValue(SchedStrategy.FCFS)) options.iter_max = 1;

        int it = 0;
        Matrix tmp_eta = new Matrix(1, M);
        for (int i = 0; i < M; i++) tmp_eta.set(i, FastMath.abs(1 - eta.get(i) / eta_1.get(i)));

        Matrix lambda = null, Lchain = null, STchain = null, Vchain = null, alpha = null;
        Matrix Q = null, U = null, R = null, T = null, X = null, STeff = null;
        Double lG = null;
        String method = null;

        while (tmp_eta.elementMax() > options.iter_tol && it < options.iter_max) {
            it++;
            eta_1 = eta;

            if (it == 1) {
                lambda = new Matrix(1, C);
                lambda.fill(0.0);
                Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
                Lchain = ret.Dchain;
                STchain = ret.STchain;
                Vchain = ret.Vchain;
                alpha = ret.alpha;
                Nchain = ret.Nchain;
                for (int c = 0; c < C; c++) {
                    Matrix inchain = sn.inchain.get(c);
                    boolean isOpenChain = false;
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        if (Utils.isInf(NK.get((int) inchain.get(0, col)))) { isOpenChain = true; break; }
                    }
                    for (int i = 0; i < M; i++) {
                        if (isOpenChain && FastMath.abs(i - sn.refstat.get((int) inchain.get(0))) < 1e-6) {
                            lambda.set(c, 1.0 / STchain.get(i, c));
                        }
                    }
                }
            } else {
                for (int c = 0; c < C; c++) {
                    Matrix inchain = sn.inchain.get(c);
                    for (int i = 0; i < M; i++) {
                        Matrix ST_tmp = new Matrix(1, inchain.getNumCols());
                        Matrix alpha_tmp = new Matrix(1, inchain.getNumCols());
                        for (int col = 0; col < inchain.getNumCols(); col++) {
                            ST_tmp.set(col, ST.get(i, (int) inchain.get(0, col)));
                            alpha_tmp.set(col, alpha.get(i, (int) inchain.get(0, col)));
                        }
                        STchain.set(i, c, ST_tmp.mult(alpha_tmp.transpose()).get(0));
                        Lchain.set(i, c, Vchain.get(i, c) * STchain.get(i, c));
                    }
                }
            }
            STchain.removeInfinity();
            Lchain.removeInfinity();

            Matrix Lms = new Matrix(M, C);
            Matrix Z = new Matrix(M, C);
            Matrix Zms = new Matrix(M, C);
            Lms.fill(0.0); Z.fill(0.0); Zms.fill(0.0);
            List<Integer> infServers = new ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                if (Utils.isInf(nservers.get(i))) {
                    infServers.add(i);
                    for (int j = 0; j < C; j++) Z.set(i, j, Lchain.get(i, j));
                } else {
                    if (options.method != null && options.method.equalsIgnoreCase("exact") && nservers.get(i) > 1) {
                        if (options.verbose != VerboseLevel.SILENT) System.out.println("SolverNC does not support exact multiserver yet. Switching to approximate method.");
                    }
                    for (int j = 0; j < C; j++) {
                        Lms.set(i, j, Lchain.get(i, j) / nservers.get(i));
                        Zms.set(i, j, Lchain.get(i, j) * (nservers.get(i) - 1) / nservers.get(i));
                    }
                }
            }
            Matrix Z_new = new Matrix(1, C);
            for (int i = 0; i < C; i++) Z_new.set(i, Z.sumCols(i) + Zms.sumCols(i));
            Matrix Z_tmp_append_0 = new Matrix(1, Z_new.getNumCols() + 1);
            Z_tmp_append_0.set(Z_new.getNumCols(), 0.0);
            Matrix.extract(Z_new, 0, 1, 0, Z_new.length(), Z_tmp_append_0, 0, 0);

            Ret.pfqnNcXQ ret = Pfqn_nc.pfqn_nc(lambda, Lms, Nchain, Z_new, options);
            lG = ret.lG;
            Matrix Xchain = ret.X;
            Matrix Qchain = ret.Q;
            method = ret.method;

            if (Zms.sumCols().elementMin() > GlobalConstants.FineTol) {
                Xchain = new Matrix(0, 0);
                Qchain = new Matrix(0, 0);
            }
            if (lG == null) {
                double runtime = (System.nanoTime() - startTime) / 1.0e9;
                return new SolverNC.SolverNCReturn(Q, U, R, T, C, X, Double.NaN, STeff, it, runtime, method);
            }

            if (Xchain.isEmpty()) {
                Xchain = lambda.copy();
                Qchain = new Matrix(M, C);
                Qchain.fill(0.0);
                for (Integer r : closedChains) {
                    Matrix Nchain_tmp = Nchain.copy();
                    Nchain_tmp.set(r, Nchain_tmp.get(r) - 1);
                    Matrix Nchain_tmp_append_1 = new Matrix(1, Nchain.getNumCols() + 1);
                    Nchain_tmp_append_1.set(Nchain.getNumCols(), 1.0);
                    Matrix.extract(Nchain_tmp, 0, 1, 0, Nchain_tmp.length(), Nchain_tmp_append_1, 0, 0);
                    Xchain.set(r, FastMath.exp(Pfqn_nc.pfqn_nc(lambda, Lms, Nchain_tmp, Z_new, options).lG - lG));
                    for (int i = 0; i < M; i++) {
                        if (Lchain.get(i, r) > 1e-6) {
                            if (Utils.isInf(nservers.get(i))) {
                                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
                            } else {
                                Matrix lambda_tmp = new Matrix(1, lambda.length() + 1);
                                lambda_tmp.fill(0.0);
                                Matrix.extract(lambda, 0, 1, 0, lambda.length(), lambda_tmp, 0, 0);

                                Matrix L_tmp = new Matrix(0, Lms.getNumCols() + 1);
                                for (int row = 0; row < Lms.getNumRows(); row++) {
                                    if (row != i) {
                                        Matrix L_tmp_row = new Matrix(1, Lms.getNumCols() + 1);
                                        L_tmp_row.set(Lms.getNumCols(), 0.0);
                                        Matrix.extract(Lms, row, row + 1, 0, Lms.getNumCols(), L_tmp_row, 0, 0);
                                        L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null);
                                    }
                                }
                                Matrix L_tmp_row = new Matrix(1, Lms.getNumCols() + 1);
                                L_tmp_row.set(Lms.getNumCols(), 1.0);
                                Matrix.extract(Lms, i, i + 1, 0, Lms.getNumCols(), L_tmp_row, 0, 0);
                                L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null);

                                Ret.pfqnNcXQ ret_tmp = Pfqn_nc.pfqn_nc(lambda_tmp, L_tmp, Nchain_tmp_append_1, Z_tmp_append_0, options);
                                double res = ret_tmp.lG;
                                method = ret_tmp.method;
                                Qchain.set(i, r, Zms.get(i, r) * Xchain.get(r) + Lms.get(i, r) * Math.exp(res - lG));
                            }
                        }
                    }
                    Qchain.removeNaN();
                }

                for (Integer r : openChains) {
                    for (int i = 0; i < M; i++) {
                        Matrix lambda_open = new Matrix(1, openChains.size());
                        Matrix Lchain_i_open = new Matrix(1, openChains.size());
                        Matrix Qchain_i_closed = new Matrix(1, closedChains.size());
                        for (int j = 0; j < openChains.size(); j++) {
                            lambda_open.set(j, lambda.get(openChains.get(j)));
                            Lchain_i_open.set(j, Lchain.get(i, openChains.get(j)));
                        }
                        for (int j = 0; j < closedChains.size(); j++) {
                            Qchain_i_closed.set(j, Qchain.get(i, closedChains.get(j)));
                        }
                        Qchain.set(i, r,
                                lambda.get(r) * Lchain.get(i, r)
                                        / (1 - lambda_open.mult(Lchain_i_open.transpose()).get(0) / nservers.get(i))
                                        * (1 + Qchain_i_closed.elementSum()));
                    }
                }
            } else {
                for (int r = 0; r < C; r++) {
                    for (int i = 0; i < M; i++) {
                        if (Lchain.get(i, r) > 1e-6) {
                            if (Utils.isInf(nservers.get(i))) {
                                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
                            }
                        }
                    }
                }
            }

            boolean allNaN = true;
            for (int i = 0; i < Xchain.length(); i++) {
                if (!Double.isNaN(Xchain.get(i))) { allNaN = false; break; }
            }
            if (allNaN) {
                if (options.verbose != VerboseLevel.SILENT) System.out.println("Normalizing constant computations produced a floating-point range exception. Model is likely too large.");
            }

            Matrix Rchain = new Matrix(Qchain.getNumRows(), Qchain.getNumCols());
            for (int i = 0; i < Rchain.getNumRows(); i++) {
                for (int j = 0; j < Rchain.getNumCols(); j++) {
                    if (Qchain.get(i, j) < GlobalConstants.Zero) Rchain.set(i, j, 0.0);
                    else Rchain.set(i, j, Qchain.get(i, j) / Xchain.get(j) / Vchain.get(i, j));
                }
            }
            for (int i = 0; i < infServers.size(); i++) {
                int row = infServers.get(i);
                for (int j = 0; j < Rchain.getNumCols(); j++) Rchain.set(row, j, Lchain.get(row, j) / Vchain.get(row, j));
            }
            Matrix Tchain = new Matrix(Vchain.getNumRows(), Vchain.getNumCols());
            for (int i = 0; i < Tchain.getNumRows(); i++) {
                for (int j = 0; j < Tchain.getNumCols(); j++) Tchain.set(i, j, Xchain.get(j) * Vchain.get(i, j));
            }
            Ret.snDeaggregateChainResults ret1 = SnDeaggregateChainResults.snDeaggregateChainResults(
                    sn, Lchain, ST, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
            Q = ret1.Q;
            U = ret1.U;
            R = ret1.R;
            T = ret1.T;
            X = ret1.X;
            STeff = ST.copy();

            Ret.npfqnNonexpApprox npfqnRet = Npfqn_nonexp_approx.npfqn_nonexp_approx(
                    options.config.highvar == null ? "interp" : options.config.highvar,
                    sn, ST0, V, SCV, T, U, gamma, nservers);
            ST = npfqnRet.ST;
            gamma = npfqnRet.gamma;
            eta = npfqnRet.eta;
            for (int i = 0; i < M; i++) tmp_eta.set(i, FastMath.abs(1 - eta.get(i) / eta_1.get(i)));
        }
        if (Q != null) Q.absEq();
        if (R != null) R.absEq();
        if (X != null) X.absEq();
        if (U != null) U.absEq();
        if (X != null) X.removeInfinity();
        if (U != null) U.removeInfinity();
        if (Q != null) Q.removeInfinity();
        if (R != null) R.removeInfinity();

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            double Nchain_c = Nchain.get(c);
            if (Double.isFinite(Nchain_c) && Q != null) {
                double sumQ = 0.0;
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int classIdx = (int) inchain.get(0, col);
                    sumQ += Q.sumCols(classIdx);
                }
                if (sumQ > 0) {
                    double ratio = Nchain_c / sumQ;
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int classIdx = (int) inchain.get(0, col);
                        for (int i = 0; i < M; i++) Q.set(i, classIdx, ratio * Q.get(i, classIdx));
                        if (X != null) X.set(classIdx, ratio * X.get(classIdx));
                        if (T != null) {
                            for (int i = 0; i < M; i++) T.set(i, classIdx, ratio * T.get(i, classIdx));
                        }
                        // MATLAB solver_nc.m rescales U with the same population
                        // ratio; omitting it leaves Util uniformly off by 1/ratio
                        if (U != null) {
                            for (int i = 0; i < M; i++) U.set(i, classIdx, ratio * U.get(i, classIdx));
                        }
                    }
                }
            }
        }
        double runtime = (System.nanoTime() - startTime) / 1.0e9;
        return new SolverNC.SolverNCReturn(Q, U, R, T, C, X, lG != null ? lG : Double.NaN, STeff, it, runtime, method);
    }
}
