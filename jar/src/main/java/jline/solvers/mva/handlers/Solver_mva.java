package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_mvac;
import jline.api.pfqn.mva.Pfqn_mvams;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnHasProductForm;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

/**
 * Handler for the solver_mva function.
 */
public final class Solver_mva {
    private Solver_mva() {}

    public static MVAResult solver_mva(NetworkStruct sn, SolverOptions options) {
        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix Vchain = ret.Vchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;
        Matrix refstatchain = ret.refstatchain;
        Matrix nservers = sn.nservers;
        int M = sn.nstations;
        int K = sn.nchains;

        List<Integer> lcfsStats = new ArrayList<Integer>();
        List<Integer> lcfsprStats = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.LCFS) lcfsStats.add(i);
            else if (s == SchedStrategy.LCFSPR) lcfsprStats.add(i);
        }

        if (!lcfsStats.isEmpty() && !lcfsprStats.isEmpty()) {
            if (lcfsStats.size() != 1 || lcfsprStats.size() != 1) {
                throw new RuntimeException("LCFS MVA requires exactly one LCFS and one LCFS-PR station.");
            }
            for (int i = 0; i < Nchain.getNumCols(); i++) {
                if (Utils.isInf(Nchain.get(i))) {
                    throw new RuntimeException("LCFS MVA requires a closed queueing network.");
                }
            }
            Matrix rt = sn.rt;
            int nclasses = sn.nclasses;
            int[] istArr = new int[]{lcfsStats.get(0), lcfsprStats.get(0)};
            for (int ist : istArr) {
                for (int r = 0; r < nclasses; r++) {
                    int idx = ist * nclasses + r;
                    if (rt.get(idx, idx) > 0) throw new RuntimeException("LCFS MVA does not support self-loops at stations.");
                }
            }
            return Solver_mva_lcfsqn.solver_mva_lcfsqn(sn, options, lcfsStats.get(0), lcfsprStats.get(0));
        } else if (!lcfsStats.isEmpty()) {
            throw new RuntimeException("LCFS scheduling requires a paired LCFS-PR station.");
        }

        List<Integer> infSET = new ArrayList<Integer>();
        List<Integer> qSET = new ArrayList<Integer>();
        if (!SnHasProductForm.snHasProductForm(sn)) {
            throw new RuntimeException("Unsupported exact MVA analysis, the model does not have a product form");
        }
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.EXT) {
                // skip
            } else if (s == SchedStrategy.INF) {
                infSET.add(i);
            } else if (s == SchedStrategy.PS || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) {
                qSET.add(i);
            } else {
                throw new RuntimeException("Unsupported exact MVA analysis for " + SchedStrategy.toText(s) + " scheduling");
            }
        }
        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        Matrix C = new Matrix(1, K);
        Matrix Wchain = new Matrix(M, K);
        Matrix Qchain = new Matrix(M, K);
        Matrix lambda = new Matrix(1, K);

        List<Integer> ocl = new ArrayList<Integer>();
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Utils.isInf(Nchain.get(i))) ocl.add(i);
        }
        for (int r : ocl) {
            lambda.set(0, r, 1.0 / STchain.get((int) refstatchain.get(r), r));
            Qchain.set((int) refstatchain.get(r), r, GlobalConstants.Inf);
        }
        List<Integer> rset = new ArrayList<Integer>();
        for (int i = 0; i < K; i++) {
            if (Nchain.get(i) != 0.0) rset.add(i);
        }
        Matrix Lp = new Matrix(qSET.size(), STchain.getNumCols());
        for (int i = 0; i < qSET.size(); i++) {
            for (int j = 0; j < Lp.getNumCols(); j++) {
                Lp.set(i, j, STchain.get(qSET.get(i), j) * Vchain.get(qSET.get(i), j));
            }
        }
        Matrix Zp = new Matrix(infSET.size(), STchain.getNumCols());
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Zp.getNumCols(); j++) {
                Zp.set(i, j, STchain.get(infSET.get(i), j) * Vchain.get(infSET.get(i), j));
            }
        }
        Matrix nserversp = new Matrix(qSET.size(), 1);
        for (int i = 0; i < qSET.size(); i++) nserversp.set(i, nservers.get(qSET.get(i)));
        Ret.pfqnMVA ret1 = Pfqn_mvams.pfqn_mvams(lambda, Lp, Nchain, Zp, Matrix.ones(qSET.size(), 1), nserversp);
        Matrix Xchain = ret1.X;
        Matrix Qpf = ret1.Q;
        double lG = ret1.lGN;
        for (int i = 0; i < qSET.size(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) {
                Qchain.set(qSET.get(i), j, Qpf.get(i, j));
            }
        }

        Matrix Q2 = new Matrix(infSET.size(), Qchain.getNumCols());
        Matrix Xchainrep = Xchain.repmat(infSET.size(), 1);
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Q2.getNumCols(); j++) {
                Q2.set(i, j, Xchainrep.get(i, j) * STchain.get(infSET.get(i), j) * Vchain.get(infSET.get(i), j));
            }
        }
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) Qchain.set(infSET.get(i), j, Q2.get(i, j));
        }

        for (int r : rset) {
            for (int k : infSET) Wchain.set(k, r, STchain.get(k, r));
            for (int k : qSET) {
                if (Utils.isInf(nservers.get(k))) {
                    Wchain.set(k, r, STchain.get(k, r));
                } else {
                    if (Vchain.get(k, r) == 0.0 || Xchain.get(r) == 0.0) Wchain.set(k, r, 0);
                    else Wchain.set(k, r, Qchain.get(k, r) / (Xchain.get(r) * Vchain.get(k, r)));
                }
            }
        }
        for (int r : rset) {
            if (Matrix.extractColumn(Wchain, r, null).elementSum() == 0.0) {
                Xchain.set(r, 0.0);
            } else {
                if (Utils.isInf(Nchain.get(r))) {
                    Matrix vt = Matrix.extractColumn(Vchain, r, null);
                    Matrix wt = Matrix.extractColumn(Wchain, r, null);
                    C.set(r, vt.transpose().mult(wt).get(0));
                } else if (Nchain.get(r) == 0.0) {
                    Xchain.set(r, 0.0);
                    C.set(r, 0.0);
                } else {
                    Matrix vt = Matrix.extractColumn(Vchain, r, null);
                    Matrix wt = Matrix.extractColumn(Wchain, r, null);
                    C.set(r, vt.transpose().mult(wt).get(0));
                    Xchain.set(r, Nchain.get(r) / C.get(r));
                }
            }
            for (int k = 0; k < M; k++) {
                Qchain.set(k, r, Xchain.get(r) * Vchain.get(k, r) * Wchain.get(k, r));
                Tchain.set(k, r, Xchain.get(r) * Vchain.get(k, r));
            }
        }
        for (int k = 0; k < M; k++) {
            for (int r : rset) {
                if (Utils.isInf(nservers.get(k))) {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r));
                } else {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r) / nservers.get(k));
                }
            }
        }
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                if (Vchain.get(k, r) * STchain.get(k, r) > options.tol) {
                    SchedStrategy s = sn.sched.get(sn.stations.get(k));
                    if (s == SchedStrategy.FCFS || s == SchedStrategy.PS) {
                        Matrix Urow = Matrix.extractRows(Uchain, k, k + 1, null);
                        double UrowSum = Urow.elementSum();
                        if (UrowSum > 1 + options.tol) {
                            Matrix Vrow = Matrix.extractRows(Vchain, k, k + 1, null);
                            Matrix STrow = Matrix.extractRows(STchain, k, k + 1, null);
                            Uchain.set(k, r, Maths.min(1.0, UrowSum) * Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r)
                                    / Vrow.elementMult(STrow, null).mult(Xchain.columnMajorOrder()).get(0));
                        }
                    }
                }
            }
        }

        Matrix Vs = null;
        for (Object key : sn.nodevisits.keySet()) {
            if (Vs == null) Vs = sn.nodevisits.get(key);
            else Vs = Vs.add(1.0, sn.nodevisits.get(key));
        }
        List<Integer> sinks = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Sink) sinks.add(i);
        }

        Matrix Rchain = new Matrix(Qchain.getNumRows(), Qchain.getNumCols());
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                if (Qchain.get(i, j) < GlobalConstants.Zero) Rchain.set(i, j, 0.0);
                else Rchain.set(i, j, Qchain.get(i, j) / Tchain.get(i, j));
            }
        }
        for (int i = 0; i < Xchain.getNumRows(); i++) {
            for (int j = 0; j < Xchain.getNumCols(); j++) {
                if (!Double.isFinite(Xchain.get(i, j))) Xchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Uchain.getNumRows(); i++) {
            for (int j = 0; j < Uchain.getNumCols(); j++) {
                if (!Double.isFinite(Uchain.get(i, j))) Uchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Qchain.getNumRows(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) {
                if (!Double.isFinite(Qchain.get(i, j))) Qchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                if (!Double.isFinite(Rchain.get(i, j))) Rchain.set(i, j, 0);
            }
        }

        List<Integer> Nzero = new ArrayList<Integer>();
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Nchain.get(i) == 0.0) Nzero.add(i);
        }
        for (int j : Nzero) {
            Xchain.set(j, 0.0);
            for (int i = 0; i < Uchain.getNumRows(); i++) Uchain.set(i, j, 0.0);
            for (int i = 0; i < Qchain.getNumRows(); i++) Qchain.set(i, j, 0.0);
            for (int i = 0; i < Rchain.getNumRows(); i++) Rchain.set(i, j, 0.0);
            for (int i = 0; i < Tchain.getNumRows(); i++) Tchain.set(i, j, 0.0);
            for (int i = 0; i < Wchain.getNumRows(); i++) Wchain.set(i, j, 0.0);
        }

        Ret.snDeaggregateChainResults ret2 = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);

        MVAResult res = new MVAResult();
        res.QN = ret2.Q;
        res.UN = ret2.U;
        res.RN = ret2.R;
        res.TN = ret2.T;
        res.CN = ret2.C;
        res.XN = ret2.X;
        res.logNormConstAggr = lG;
        return res;
    }

    /**
     * Exact mean value analysis by chain (MVAC, Conway-de Souza e Silva-Lavenberg
     * 1989) via {@link Pfqn_mvac}. Closed product-form networks of single-server
     * fixed-rate (SSFR) queues and infinite-server centers only.
     */
    public static MVAResult solver_mvac(NetworkStruct sn, SolverOptions options) {
        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix Vchain = ret.Vchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;
        Matrix nservers = sn.nservers;
        int M = sn.nstations;
        int K = sn.nchains;

        if (!SnHasProductForm.snHasProductForm(sn)) {
            throw new RuntimeException("MVAC requires a product-form model.");
        }
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Utils.isInf(Nchain.get(i))) {
                throw new RuntimeException("MVAC supports closed models only; use method 'exact' for open/mixed networks.");
            }
        }

        List<Integer> infSET = new ArrayList<Integer>();
        List<Integer> qSET = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.EXT) {
                // skip
            } else if (s == SchedStrategy.INF) {
                infSET.add(i);
            } else if (s == SchedStrategy.PS || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) {
                qSET.add(i);
            } else {
                throw new RuntimeException("MVAC does not support " + SchedStrategy.toText(s) + " scheduling");
            }
        }
        for (int k : qSET) {
            if (nservers.get(k) != 1.0) {
                throw new RuntimeException("MVAC supports single-server (SSFR) queues only; use method 'exact' for multiserver stations.");
            }
        }

        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        Matrix C = new Matrix(1, K);
        Matrix Wchain = new Matrix(M, K);
        Matrix Qchain = new Matrix(M, K);

        List<Integer> rset = new ArrayList<Integer>();
        for (int i = 0; i < K; i++) {
            if (Nchain.get(i) != 0.0) rset.add(i);
        }
        Matrix Lp = new Matrix(qSET.size(), STchain.getNumCols());
        for (int i = 0; i < qSET.size(); i++) {
            for (int j = 0; j < Lp.getNumCols(); j++) {
                Lp.set(i, j, STchain.get(qSET.get(i), j) * Vchain.get(qSET.get(i), j));
            }
        }
        Matrix Zp = new Matrix(infSET.size(), STchain.getNumCols());
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Zp.getNumCols(); j++) {
                Zp.set(i, j, STchain.get(infSET.get(i), j) * Vchain.get(infSET.get(i), j));
            }
        }
        Ret.pfqnMVAC ret1 = Pfqn_mvac.pfqn_mvac(Lp, Nchain, Zp);
        Matrix Xchain = ret1.X;
        Matrix Qpf = ret1.Q;
        double lG = Double.NaN; // MVAC forms no normalizing constant
        for (int i = 0; i < qSET.size(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) {
                Qchain.set(qSET.get(i), j, Qpf.get(i, j));
            }
        }

        Matrix Q2 = new Matrix(infSET.size(), Qchain.getNumCols());
        Matrix Xchainrep = Xchain.repmat(infSET.size(), 1);
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Q2.getNumCols(); j++) {
                Q2.set(i, j, Xchainrep.get(i, j) * STchain.get(infSET.get(i), j) * Vchain.get(infSET.get(i), j));
            }
        }
        for (int i = 0; i < infSET.size(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) Qchain.set(infSET.get(i), j, Q2.get(i, j));
        }

        for (int r : rset) {
            for (int k : infSET) Wchain.set(k, r, STchain.get(k, r));
            for (int k : qSET) {
                if (Utils.isInf(nservers.get(k))) {
                    Wchain.set(k, r, STchain.get(k, r));
                } else {
                    if (Vchain.get(k, r) == 0.0 || Xchain.get(r) == 0.0) Wchain.set(k, r, 0);
                    else Wchain.set(k, r, Qchain.get(k, r) / (Xchain.get(r) * Vchain.get(k, r)));
                }
            }
        }
        for (int r : rset) {
            if (Matrix.extractColumn(Wchain, r, null).elementSum() == 0.0) {
                Xchain.set(r, 0.0);
            } else if (Nchain.get(r) == 0.0) {
                Xchain.set(r, 0.0);
                C.set(r, 0.0);
            } else {
                Matrix vt = Matrix.extractColumn(Vchain, r, null);
                Matrix wt = Matrix.extractColumn(Wchain, r, null);
                C.set(r, vt.transpose().mult(wt).get(0));
                Xchain.set(r, Nchain.get(r) / C.get(r));
            }
            for (int k = 0; k < M; k++) {
                Qchain.set(k, r, Xchain.get(r) * Vchain.get(k, r) * Wchain.get(k, r));
                Tchain.set(k, r, Xchain.get(r) * Vchain.get(k, r));
            }
        }
        for (int k = 0; k < M; k++) {
            for (int r : rset) {
                if (Utils.isInf(nservers.get(k))) {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r));
                } else {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r) / nservers.get(k));
                }
            }
        }
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                if (Vchain.get(k, r) * STchain.get(k, r) > options.tol) {
                    SchedStrategy s = sn.sched.get(sn.stations.get(k));
                    if (s == SchedStrategy.FCFS || s == SchedStrategy.PS) {
                        Matrix Urow = Matrix.extractRows(Uchain, k, k + 1, null);
                        double UrowSum = Urow.elementSum();
                        if (UrowSum > 1 + options.tol) {
                            Matrix Vrow = Matrix.extractRows(Vchain, k, k + 1, null);
                            Matrix STrow = Matrix.extractRows(STchain, k, k + 1, null);
                            Uchain.set(k, r, Maths.min(1.0, UrowSum) * Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r)
                                    / Vrow.elementMult(STrow, null).mult(Xchain.columnMajorOrder()).get(0));
                        }
                    }
                }
            }
        }

        Matrix Rchain = new Matrix(Qchain.getNumRows(), Qchain.getNumCols());
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                if (Qchain.get(i, j) < GlobalConstants.Zero) Rchain.set(i, j, 0.0);
                else Rchain.set(i, j, Qchain.get(i, j) / Tchain.get(i, j));
            }
        }
        for (int i = 0; i < Xchain.getNumRows(); i++) {
            for (int j = 0; j < Xchain.getNumCols(); j++) {
                if (!Double.isFinite(Xchain.get(i, j))) Xchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Uchain.getNumRows(); i++) {
            for (int j = 0; j < Uchain.getNumCols(); j++) {
                if (!Double.isFinite(Uchain.get(i, j))) Uchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Qchain.getNumRows(); i++) {
            for (int j = 0; j < Qchain.getNumCols(); j++) {
                if (!Double.isFinite(Qchain.get(i, j))) Qchain.set(i, j, 0);
            }
        }
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                if (!Double.isFinite(Rchain.get(i, j))) Rchain.set(i, j, 0);
            }
        }

        List<Integer> Nzero = new ArrayList<Integer>();
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Nchain.get(i) == 0.0) Nzero.add(i);
        }
        for (int j : Nzero) {
            Xchain.set(j, 0.0);
            for (int i = 0; i < Uchain.getNumRows(); i++) Uchain.set(i, j, 0.0);
            for (int i = 0; i < Qchain.getNumRows(); i++) Qchain.set(i, j, 0.0);
            for (int i = 0; i < Rchain.getNumRows(); i++) Rchain.set(i, j, 0.0);
            for (int i = 0; i < Tchain.getNumRows(); i++) Tchain.set(i, j, 0.0);
            for (int i = 0; i < Wchain.getNumRows(); i++) Wchain.set(i, j, 0.0);
        }

        Ret.snDeaggregateChainResults ret2 = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);

        MVAResult res = new MVAResult();
        res.QN = ret2.Q;
        res.UN = ret2.U;
        res.RN = ret2.R;
        res.TN = ret2.T;
        res.CN = ret2.C;
        res.XN = ret2.X;
        res.logNormConstAggr = lG;
        return res;
    }
}
