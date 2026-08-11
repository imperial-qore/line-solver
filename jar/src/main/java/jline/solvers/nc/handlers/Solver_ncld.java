package jline.solvers.nc.handlers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.api.npfqn.Npfqn_nonexp_approx;
import jline.api.pfqn.ld.Pfqn_fnc;
import jline.api.pfqn.ld.Pfqn_mushift;
import jline.api.pfqn.ld.Pfqn_mvaldmx;
import jline.api.pfqn.ld.Pfqn_ncld;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnGetProductFormParams;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Solver_ncld {
    private Solver_ncld() {}

    public static SolverNC.SolverNCLDReturn solver_ncld(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix nservers = sn.nservers;
        Matrix nserversFinite = nservers.copy();
        nserversFinite.removeInfinity();
        double minFiniteServer = Double.MAX_VALUE;
        for (int i = 0; i < nservers.getNumElements(); i++) {
            if (Double.isFinite(nservers.get(i)) && nservers.get(i) < minFiniteServer) {
                minFiniteServer = nservers.get(i);
            }
        }

        if (minFiniteServer < Double.MAX_VALUE && minFiniteServer > 1) {
            if (sn.lldscaling.isEmpty() && M == 2 && Double.isFinite(sn.njobs.elementMaxAbs())) {
                double Nt = sn.njobs.elementSum();
                sn.lldscaling = sn.lldscaling.concatCols(new Matrix(M, (int) Nt));
                for (int i = 0; i < M; i++) {
                    int j = 0;
                    while (j < Nt) {
                        sn.lldscaling.set(i, j, FastMath.min((double) (j + 1), sn.nservers.get(i)));
                        j++;
                    }
                }
            } else if (lldEncodesMultiserver(sn, nservers, M)) {
                // The caller (SolverNC) already expressed every multiserver as
                // mu(n)=min(n,c), which is what this guard asks for, so the model is
                // supported. nservers is deliberately left at c: utilization is the
                // fraction of the c servers busy and c cannot be read back from
                // lldscaling once the population is below it (min(1:Nt,c) is then 1:Nt).
            } else {
                throw new RuntimeException("The load-dependent solver does not support multi-server stations yet. Specify multi-server stations via limited load-dependence.");
            }
        }
        if ((sn.cdscaling != null && !sn.cdscaling.isEmpty())
                || (sn.jdscaling != null && !sn.jdscaling.isEmpty())) {
            // Class-dependent beta_{i,r}(n) or joint-dependent eta_i(n) service
            // rates: route to the convolution solver, which builds the station
            // factor X_m(n) by Sauer's
            // chain-dependent recurrence (Sauer 1983, Sec. 5.2, eq. (40)). This
            // covers both the chain-independent case (a scalar-valued function) and
            // the chain-specific one (a length-R vector), the latter being what the
            // flow-equivalent-server aggregation installs.
            //
            // The routing is unconditional on options.method: the algorithms below
            // read mu(n) from lldscaling only and have no way to apply beta, so
            // gating this on method='exact' made every other method (including the
            // default) return the UNSCALED network instead of erroring. Convolution
            // is the only exact algorithm for beta, and it is what the native
            // Python NC already does.
            return Solver_nc_conv.solver_nc_conv(sn, options);
        }
        Matrix NK = sn.njobs.transpose();
        // Mixed open/closed load-dependent networks are handled below through the
        // exact chain-level MVALDMX algorithm (Bruell-Balbo-Afshari effective
        // capacity); the closed-only convolution path is retained otherwise.
        int C = sn.nchains;
        Matrix SCV = sn.scv;
        Matrix gamma = new Matrix(M, 1);
        gamma.zero();
        Matrix V = new Matrix(sn.nstateful, K);
        for (int i = 0; i < sn.visits.size(); i++) {
            V = V.add(1.0, sn.visits.get(i));
        }
        Matrix ST = Matrix.ones(sn.rates.getNumRows(), sn.rates.getNumCols()).elementDiv(sn.rates);
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                if (Double.isNaN(ST.get(i, j))) {
                    ST.set(i, j, 0);
                }
            }
        }
        Matrix ST0 = ST.copy();
        Matrix lldscaling = sn.lldscaling;
        Matrix NKfinite = NK.copy();
        NKfinite.removeInfinity();
        double Nt = NKfinite.elementSum();
        if (lldscaling.isEmpty()) {
            lldscaling = Matrix.ones(M, (int) FastMath.ceil(Nt));
        }
        Ret.snGetDemands demandsChainReturn = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Vchain = demandsChainReturn.Vchain;
        Matrix alpha = demandsChainReturn.alpha;
        Matrix eta_1 = new Matrix(1, M);
        eta_1.zero();
        Matrix eta = Matrix.ones(1, M);
        if (!sn.sched.containsValue(SchedStrategy.FCFS)) {
            options.iter_max = 1;
        }
        int iter = 0;
        long Tstart = System.nanoTime();
        Ret.snDeaggregateChainResults snDeaggragatedChains = null;
        Matrix lambda = null;
        Matrix Lchain = null;
        Matrix STchain = null;
        Matrix Q = null;
        Matrix Cmat = null;
        Matrix U = null;
        Matrix R = null;
        Matrix T = null;
        Matrix X = null;
        Double lG = null;
        String method = null;

        while (Matrix.ones(eta.getNumRows(), eta.getNumCols()).sub(eta.elementDiv(eta_1))
                .elementMaxAbs() > options.iter_tol && iter < options.iter_max) {
            iter += 1;
            eta_1 = eta.copy();
            M = sn.nstations;
            K = sn.nclasses;
            C = sn.nchains;
            Lchain = new Matrix(M, C);
            Lchain.zero();
            STchain = new Matrix(M, C);
            STchain.zero();
            Matrix SCVchain = new Matrix(M, C);
            SCVchain.zero();
            Matrix Nchain = new Matrix(1, C);
            Nchain.zero();
            Matrix refstatchain = new Matrix(C, 1);
            refstatchain.zero();
            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                boolean isOpenChain = false;
                for (int i = 0; i < inchain.getNumElements(); i++) {
                    if (Utils.isInf(sn.njobs.get((int) inchain.get(i)))) {
                        isOpenChain = true;
                    }
                }
                for (int i = 0; i < M; i++) {
                    Matrix STinchain = new Matrix(1, inchain.getNumElements());
                    Matrix alphainchain = new Matrix(1, inchain.getNumElements());
                    Matrix SCVinchain = new Matrix(1, inchain.getNumElements());
                    for (int j = 0; j < inchain.getNumElements(); j++) {
                        STinchain.set(j, ST.get(i, (int) inchain.get(j)));
                        alphainchain.set(j, alpha.get(i, (int) inchain.get(j)));
                        SCVinchain.set(j, SCV.get(i, (int) inchain.get(j)));
                    }
                    double lchainVal = Vchain.get(i, c) * STinchain.mult(alphainchain.transpose()).toDouble();
                    Lchain.set(i, c, lchainVal);
                    STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble());
                    if (isOpenChain && (double) i == sn.refstat.get((int) inchain.get(0))) {
                        Matrix STinchainFinite = STinchain.copy();
                        STinchainFinite.removeInfinity();
                        STchain.set(i, c, STinchainFinite.elementSum());
                    } else {
                        STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble());
                    }
                    SCVchain.set(i, c, SCVinchain.mult(alphainchain.transpose()).toDouble());
                }
                Matrix NKinchain = new Matrix(1, inchain.getNumElements());
                for (int i = 0; i < inchain.getNumElements(); i++) {
                    NKinchain.set(i, NK.get((int) inchain.get(i)));
                }
                Nchain.set(c, NKinchain.elementSum());
                refstatchain.set(c, sn.refstat.get((int) inchain.get(0)));
                if ((sn.refstat.get((int) inchain.get(0)) - refstatchain.get(c)) != 0.0) {
                    throw new RuntimeException(String.format("Classes in chain %d have different reference station.", c));
                }
            }
            for (int i = 0; i < STchain.getNumRows(); i++) {
                for (int j = 0; j < STchain.getNumCols(); j++) {
                    if (!Double.isFinite(STchain.get(i, j))) {
                        STchain.set(i, j, 0);
                    }
                    if (!Double.isFinite(Lchain.get(i, j))) {
                        Lchain.set(i, j, 0);
                    }
                }
            }
            // chain-level arrival rates for open chains (source ST = 1/arrival rate)
            Matrix lambdaChain = new Matrix(1, C);
            lambdaChain.zero();
            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                boolean isOpen = false;
                for (int i = 0; i < inchain.getNumElements(); i++) {
                    if (Utils.isInf(sn.njobs.get((int) inchain.get(i)))) {
                        isOpen = true;
                    }
                }
                if (isOpen) {
                    int rst = (int) sn.refstat.get((int) inchain.get(0));
                    if (STchain.get(rst, c) > 0) {
                        lambdaChain.set(c, 1.0 / STchain.get(rst, c));
                    }
                }
            }
            Tstart = System.nanoTime();
            Matrix Nchainfinite = Nchain.copy();
            Nchainfinite.removeInfinity();
            Nt = Nchainfinite.elementSum();
            Matrix L = new Matrix(M, C);
            L.zero();
            Matrix mu = new Matrix(M, (int) FastMath.ceil(Nt));
            List<Integer> infServers = new ArrayList<Integer>();
            Matrix Z = L.copy();
            for (int i = 0; i < M; i++) {
                if (Utils.isInf(nservers.get(i))) {
                    infServers.add(i);
                    for (int j = 0; j < C; j++) {
                        L.set(i, j, Lchain.get(i, j));
                        Z.set(i, j, Lchain.get(i, j));
                    }
                    int j = 0;
                    while (j < Math.ceil(Nt)) {
                        mu.set(i, j, j + 1);
                        j++;
                    }
                } else {
                    if (options.method.equalsIgnoreCase("exact") && nservers.get(i) > 1) {
                        InputOutput.line_warning("solver_ncld", "SolverNC does not support exact multiserver yet. Switching to approximate method.");
                    }
                    for (int j = 0; j < C; j++) {
                        L.set(i, j, Lchain.get(i, j));
                    }
                    int j = 0;
                    while (j < Math.ceil(Nt)) {
                        mu.set(i, j, lldscaling.get(i, j));
                        j++;
                    }
                }
            }
            Matrix Qchain = new Matrix(M, C);
            Qchain.zero();
            Matrix Nchain0 = Nchain.copy();
            Nchain0.zero();
            List<Integer> openChains = new ArrayList<Integer>();
            for (int c = 0; c < C; c++) {
                if (Utils.isInf(Nchain.get(c))) {
                    openChains.add(c);
                }
            }
            Matrix Xchain;
            if (!openChains.isEmpty()) {
                // Mixed limited load-dependent network: exact chain-level MVALDMX.
                // Open chains enter via arrival rates lambda; their reference
                // (source) stations carry only the 1/lambda bookkeeping demand and
                // are excluded. Delay stations fold into the think-time vector; the
                // remaining queueing stations carry the lldscaling rates.
                java.util.Set<Integer> sourceSet = new java.util.HashSet<Integer>();
                for (int c : openChains) {
                    sourceSet.add((int) refstatchain.get(c));
                }
                List<Integer> delayStations = new ArrayList<Integer>();
                for (Integer i : infServers) {
                    if (!sourceSet.contains(i)) {
                        delayStations.add(i);
                    }
                }
                List<Integer> queueStations = new ArrayList<Integer>();
                for (int i = 0; i < M; i++) {
                    if (!sourceSet.contains(i) && !delayStations.contains(i)) {
                        queueStations.add(i);
                    }
                }
                int nq = queueStations.size();
                int ncol = Math.max(1, (int) Nchainfinite.elementSum());
                Matrix Zvec = new Matrix(1, C);
                Zvec.zero();
                for (Integer di : delayStations) {
                    for (int c = 0; c < C; c++) {
                        Zvec.set(c, Zvec.get(c) + Lchain.get(di, c));
                    }
                }
                Matrix Dq = new Matrix(nq, C);
                Matrix muq = Matrix.ones(nq, ncol);
                int availCols = Math.min(ncol, lldscaling.getNumCols());
                for (int qi = 0; qi < nq; qi++) {
                    int ist = queueStations.get(qi);
                    for (int c = 0; c < C; c++) {
                        Dq.set(qi, c, Lchain.get(ist, c));
                    }
                    for (int k = 0; k < availCols; k++) {
                        muq.set(qi, k, lldscaling.get(ist, k));
                    }
                }
                Ret.pfqnMVALDMX mret = Pfqn_mvaldmx.pfqn_mvaldmx(lambdaChain, Dq, Nchain, Zvec, muq, Matrix.ones(nq, 1));
                Xchain = mret.X;
                lG = mret.lG;
                method = "ncldmx";
                for (int qi = 0; qi < nq; qi++) {
                    int ist = queueStations.get(qi);
                    for (int c = 0; c < C; c++) {
                        Qchain.set(ist, c, mret.Q.get(qi, c));
                    }
                }
                for (Integer di : delayStations) {
                    for (int c = 0; c < C; c++) {
                        Qchain.set(di, c, Lchain.get(di, c) * Xchain.get(c));
                    }
                }
            } else {
            Ret.pfqnNc ret = Pfqn_ncld.pfqn_ncld(L, Nchain, Nchain0, mu, options);
            lG = ret.lG;
            method = ret.method;
            Xchain = new Matrix(1, 0);
            if (Xchain.isEmpty()) {
                Matrix lGr = new Matrix(1, C);
                Matrix lGhat_fnci = new Matrix(1, C);
                Matrix lGhatir = new Matrix(1, C);
                Matrix lGr_i = new Matrix(1, C);
                Matrix ldDemand = new Matrix(M, C);
                for (int r = 0; r < C; r++) {
                    Matrix Nchain_r = Matrix.oner(Nchain, new ArrayList<Integer>(Arrays.asList(r)));
                    lGr.set(r, Pfqn_ncld.pfqn_ncld(L, Nchain_r, Nchain0, mu, options).lG);
                    double xchainVal = Math.exp(lGr.get(r) - lG);
                    Xchain = Xchain.concatCols(Matrix.singleton(xchainVal));
                    for (int i = 0; i < M; i++) {
                        Qchain.set(i, r, 0);
                    }
                    Matrix CQchain_r = new Matrix(M, 1);
                    CQchain_r.zero();
                    if (M == 2 && Utils.isInf(sn.nservers.elementMaxAbs())) {
                        int firstDelay = -1;
                        for (int i = 0; i < sn.nservers.getNumElements(); i++) {
                            if (Utils.isInf(sn.nservers.get(i))) {
                                firstDelay = i;
                                break;
                            }
                        }
                        Qchain.set(firstDelay, r, Lchain.get(firstDelay, r) * Xchain.get(r));
                        for (int i = 0; i < M; i++) {
                            if (i != firstDelay) {
                                Qchain.set(i, r, Nchain.get(r) - Lchain.get(firstDelay, r) * Xchain.get(r));
                            }
                        }
                    } else {
                        for (int i = 0; i < M; i++) {
                            Matrix Lms_i = Matrix.extractRows(L, i, i + 1, null);
                            Matrix.extractRows(mu, i, i + 1, null);
                            Matrix muhati = Pfqn_mushift.pfqn_mushift(mu, i);
                            Ret.pfqnFnc fncret = Pfqn_fnc.pfqn_fnc(Matrix.extractRows(muhati, i, i + 1, null));
                            Matrix muhati_f = fncret.mu;
                            double cval = fncret.c.get(0);
                            if (Lchain.get(i, r) > 0) {
                                if (Utils.isInf(nservers.get(i))) {
                                    Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
                                } else {
                                    if (i == (M - 1) && nserversFinite.elementSum() == 1.0) {
                                        double Lchainsum = 0.0;
                                        double Qchainsum = 0.0;
                                        for (int j = 0; j < nservers.getNumElements(); j++) {
                                            if (Utils.isInf(nservers.get(j))) {
                                                Lchainsum += Lchain.get(j, r);
                                            } else {
                                                Qchainsum += Qchain.get(j, r);
                                            }
                                        }
                                        Qchain.set(i, r,
                                                FastMath.max(0.0, Nchain.get(r) - Lchainsum * Xchain.get(r) - Qchainsum));
                                    } else {
                                        lGhat_fnci.set(r,
                                                Pfqn_ncld.pfqn_ncld(Matrix.concatRows(L, Matrix.extractRows(L, i, i + 1, null), null),
                                                        Nchain_r, Nchain0,
                                                        Matrix.concatRows(muhati, muhati_f, null), options).lG);
                                        lGhatir.set(r, Pfqn_ncld.pfqn_ncld(L, Nchain_r, Nchain0, muhati, options).lG);
                                        lGr_i.set(r, Pfqn_ncld.pfqn_ncld(Lms_i, Nchain_r, Nchain0, muhati, options).lG);
                                        double dlGa = lGhat_fnci.get(r) - lGhatir.get(r);
                                        double dlG_i = lGr_i.get(r) - lGhatir.get(r);
                                        CQchain_r.set(i, (Math.exp(dlGa) - 1) + cval * (Math.exp(dlG_i) - 1));
                                        ldDemand.set(i, r,
                                                FastMath.log(L.get(i, r)) + lGhatir.get(r) - FastMath.log(mu.get(i, 0)) - lGr.get(r));
                                        Qchain.set(i, r,
                                                FastMath.exp(ldDemand.get(i, r)) * Xchain.get(r) * (1 + CQchain_r.get(i)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            } // close open/closed dispatch

            boolean allNaN = true;
            for (double v : Xchain.toArray1D()) {
                if (!Double.isNaN(v)) { allNaN = false; break; }
            }
            if (allNaN) {
                InputOutput.line_warning("solver_ncld", "Normalizing constant computations produced a floating-point range exception. Model is likely too large.");
            }
            Z = Z.sumCols();
            Matrix Rchain = Qchain.elementDivide(Xchain.repmat(M, 1)).elementDivide(Vchain);
            for (Integer i : infServers) {
                for (int j = 0; j < Rchain.getNumCols(); j++) {
                    Rchain.set(i, j, Lchain.get(i, j) / Vchain.get(i, j));
                }
            }
            Matrix Tchain = Xchain.repmat(M, 1).elementMult(Vchain, null);
            Tchain.elementMult(Lchain, null);
            Nchain.elementDiv(Xchain).sub(Z);
            // An EMPTY chain (Nchain(c) == 0) must report zero everywhere, mirroring
            // solver_ncld.m:270-272. This has to happen HERE, upstream of the
            // deaggregation, because the helper derives per-class Q from Rchain and
            // Xchain and then R from Q/T: zeroing after the fact leaves R alone.
            // Without it, the infinite-server overwrite just above
            // (Rchain = Lchain/Vchain, which does not consult the population) leaves
            // the delay's service demand in Rchain for an empty chain, and that value
            // survives into the reported RespT at EVERY station. Measured on a closed
            // model with a zero-population second chain: the JAR reported RespT = the
            // Delay's mean (1.0, and 3.0 when the mean was changed to 3.0) at both
            // stations where MATLAB reports 0. With no infinite server present Rchain
            // came from 0/0 -> NaN -> removeNaN -> 0 and the two agreed, which is why
            // the divergence only shows up when the model has a delay station.
            for (int c = 0; c < Nchain.getNumElements(); c++) {
                if (Nchain.get(c) == 0) {
                    Xchain.set(0, c, 0.0);
                    for (int i = 0; i < M; i++) {
                        Rchain.set(i, c, 0.0);
                        Tchain.set(i, c, 0.0);
                    }
                }
            }
            snDeaggragatedChains = SnDeaggregateChainResults.snDeaggregateChainResults(
                    sn, Lchain, ST, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
            Q = snDeaggragatedChains.Q;
            U = snDeaggragatedChains.U;
            R = snDeaggragatedChains.R;
            T = snDeaggragatedChains.T;
            Cmat = snDeaggragatedChains.C;
            X = snDeaggragatedChains.X;
            Ret.npfqnNonexpApprox NPFQNret = Npfqn_nonexp_approx.npfqn_nonexp_approx(
                    options.config.highvar == null ? "default" : options.config.highvar,
                    sn, ST0.copy(), V, SCV, T, U, gamma, nservers);
            ST = NPFQNret.ST;
            gamma = NPFQNret.gamma;
            eta = NPFQNret.eta.transpose();
        }
        SnGetProductFormParams.snGetProductFormParams(sn);
        double runtime = (System.nanoTime() - Tstart) / 1_000_000_000.0;
        Q = snDeaggragatedChains.Q;
        Q.absEq();
        R = snDeaggragatedChains.R;
        R.absEq();
        U = snDeaggragatedChains.U;
        U.absEq();
        X = snDeaggragatedChains.X;
        X.absEq();
        for (int i = 0; i < M; i++) {
            List<Integer> openClasses = new ArrayList<Integer>();
            for (int j = 0; j < NK.getNumElements(); j++) {
                if (Utils.isInf(NK.get(j))) {
                    openClasses.add(j);
                }
            }
            if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
                for (int r = 0; r < K; r++) {
                    int c = (int) Matrix.extractColumn(sn.chains, r, null).find().value();
                    if (!openClasses.contains(r) && snDeaggragatedChains.X.get(r) > 0) {
                        U.set(i, r,
                                X.get(r) * sn.visits.get(c).get(i, r) / sn.visits.get(c).get((int) sn.refstat.get(r), r)
                                        * ST.get(i, r) / sn.nservers.get(i));
                    } else if (openClasses.contains(r) && lambda != null && lambda.get(r) > 0) {
                        U.set(i, r,
                                lambda.get(r) * sn.visits.get(c).get(i, r) / sn.visits.get(c).get((int) sn.refstat.get(r), r)
                                        * ST.get(i, r) / sn.nservers.get(i));
                    }
                }
            } else if (Utils.isInf(sn.nservers.get(i))) {
                for (int r = 0; r < K; r++) {
                    int c = (int) Matrix.extractColumn(sn.chains, r, null).find().get(0);
                    if (!openClasses.contains(r) && snDeaggragatedChains.X.get(r) > 0) {
                        U.set(i, r,
                                X.get(r) * sn.visits.get(c).get(i, r) / sn.visits.get(c).get((int) sn.refstat.get(r), r) * ST.get(i, r));
                    } else if (openClasses.contains(r) && lambda != null && lambda.get(r) > 0) {
                        U.set(i, r,
                                lambda.get(r) * sn.visits.get(c).get(i, r) / sn.visits.get(c).get((int) sn.refstat.get(r), r) * ST.get(i, r));
                    }
                }
            } else {
                for (int j = 0; j < U.getNumCols(); j++) {
                    U.set(i, j, U.get(i, j) / Matrix.extractRows(lldscaling, i, i + 1, null).elementMax());
                }
                if (Matrix.extractRows(U, i, i + 1, null).elementSum() > 1) {
                    Matrix Uinotnan = Matrix.extractRows(U, i, i + 1, null);
                    Uinotnan.removeNaN();
                    for (int j = 0; j < U.getNumCols(); j++) {
                        U.set(i, j, U.get(i, j) / Uinotnan.elementSum());
                    }
                }
            }
        }
        X = snDeaggragatedChains.X;
        X.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        X.apply(Double.NaN, 0.0, "equal");
        U.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        U.apply(Double.NaN, 0.0, "equal");
        Q.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        Q.apply(Double.NaN, 0.0, "equal");
        R.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        R.apply(Double.NaN, 0.0, "equal");

        Matrix Nchain_final = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            double sum = 0.0;
            for (int col = 0; col < inchain.getNumCols(); col++) {
                sum += NK.get((int) inchain.get(0, col));
            }
            Nchain_final.set(c, sum);
        }

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            double Nchain_c = Nchain_final.get(c);
            if (Double.isFinite(Nchain_c)) {
                double sumQ = 0.0;
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int classIdx = (int) inchain.get(0, col);
                    sumQ += Q.sumCols(classIdx);
                }
                double ratio;
                if (sumQ > 0) {
                    ratio = Nchain_c / sumQ;
                } else {
                    ratio = 0.0;
                }
                if (sumQ > 0) {
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int classIdx = (int) inchain.get(0, col);
                        for (int i = 0; i < M; i++) {
                            Q.set(i, classIdx, ratio * Q.get(i, classIdx));
                        }
                        X.set(classIdx, ratio * X.get(classIdx));
                        for (int i = 0; i < M; i++) {
                            snDeaggragatedChains.T.set(i, classIdx, ratio * snDeaggragatedChains.T.get(i, classIdx));
                        }
                        for (int i = 0; i < M; i++) {
                            U.set(i, classIdx, ratio * U.get(i, classIdx));
                        }
                        for (int i = 0; i < M; i++) {
                            if (snDeaggragatedChains.T.get(i, classIdx) > 0) {
                                R.set(i, classIdx, Q.get(i, classIdx) / snDeaggragatedChains.T.get(i, classIdx));
                            }
                        }
                    }
                }
            }
        }

        return new SolverNC.SolverNCLDReturn(Q, U, R, snDeaggragatedChains.T, Cmat, X, lG, runtime, iter, method);
    }

    /**
     * True if every finite multi-server station already carries mu(n)=min(n,c) in
     * lldscaling, i.e. the multiserver is fully described by the load-dependent
     * rates and needs no further handling.
     */
    private static boolean lldEncodesMultiserver(NetworkStruct sn, Matrix nservers, int M) {
        if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
            return false;
        }
        double Ntot = 0;
        for (int r = 0; r < sn.njobs.length(); r++) {
            double n = sn.njobs.get(r);
            if (Double.isFinite(n)) {
                Ntot += n;
            }
        }
        if (!Double.isFinite(Ntot) || Ntot < 1 || sn.lldscaling.getNumCols() < (int) Ntot) {
            return false;
        }
        for (int i = 0; i < M; i++) {
            double c = nservers.get(i);
            if (Double.isFinite(c) && c > 1) {
                for (int j = 0; j < (int) Ntot; j++) {
                    if (sn.lldscaling.get(i, j) != FastMath.min((double) (j + 1), c)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}
