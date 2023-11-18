package jline.solvers.mva.handlers;

import jline.api.PFQN;
import jline.api.SN;
import jline.util.Maths;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;

import java.util.ArrayList;

/**
 * Handler for the solver_mva function.
 */
public class MVAHandler implements MVASolverHandler{
    @Override
    public SolverMVAResult solve(NetworkStruct sn, SolverOptions options) {
        SN.snGetDemandsChainReturn ret = SN.snGetDemandsChain(sn);
        Matrix Lchain = ret.Lchain;
        Matrix STchain = ret.STchain;
        Matrix Vchain = ret.Vchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;
        Matrix refstatchain = ret.refstatchain;
        Matrix nservers = sn.nservers;
        // We're using schedStrategies instead of schedid, as schedid is not available in NetworkStruct
        int M = sn.nstations;
        int K = sn.nchains;
        ArrayList<Integer> infSET = new ArrayList<>(); // set of infinite server stations
        ArrayList<Integer> qSET = new ArrayList<>(); // set of other product-form stations
        for(int i = 0; i < M; i++){
            switch (sn.sched.get(sn.stations.get(i))){
                case EXT:
                    // No-op
                    break;
                case INF:
                    infSET.add(i);
                    break;
                case PS: case LCFSPR:
                    qSET.add(i);
                    break;
                case FCFS: case SIRO:
                    double min = 0, max = 0;
                    for(int j = 0; j < sn.rates.getNumCols(); j++){
                        if(j == 0){
                            min = sn.rates.get(i, j);
                            max = sn.rates.get(i, j);
                        } else {
                            double current = sn.rates.get(i, j);
                            if(current > max)
                                max = current;
                            if(current < min)
                                min = current;
                        }
                    }
                    if(max - min == 0){
                        qSET.add(i);
                    }
                    break;
                default:
                    throw new RuntimeException("Unsupported exact MVA analysis for " +
                            SchedStrategy.toText(sn.sched.get(sn.stations.get(i))) + " scheduling");
            }
        }
        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        Matrix C = new Matrix(1, K);
        Matrix Wchain = new Matrix(M, K);
        Matrix Qchain = new Matrix(M, K);
        Matrix lambda = new Matrix(1, K);

        ArrayList<Integer> ocl = new ArrayList<>();
        for(int i = 0; i < Nchain.getNumCols(); i++){
            if(Double.isInfinite(Nchain.get(i))){
                ocl.add(i);
            }
        }
        for(int r : ocl){ // open classes
            lambda.set(0, r, 1 / STchain.get((int) refstatchain.get(r), r));
            Qchain.set((int) refstatchain.get(r), r, Double.POSITIVE_INFINITY);
        }
        ArrayList<Integer> rset = new ArrayList<>();
        for(int i = 0; i < K; i++){
            if(Nchain.get(i) != 0){
                rset.add(i);
            }
        }
        Matrix Lp = new Matrix(qSET.size(), STchain.getNumCols());
        for(int i = 0; i < qSET.size(); i++){
            for(int j = 0; j < Lp.getNumCols(); j++){
                Lp.set(i, j, STchain.get(qSET.get(i), j) * Vchain.get(qSET.get(i), j));
            }
        }
        Matrix Zp = new Matrix(infSET.size(), STchain.getNumCols());
        for(int i = 0; i < infSET.size(); i++){
            for(int j = 0; j < Zp.getNumCols(); j++){
                Zp.set(i, j, STchain.get(infSET.get(i), j) * Vchain.get(infSET.get(i), j));
            }
        }
        Matrix nserversp = new Matrix(qSET.size(), 1);
        for(int i = 0; i < qSET.size(); i++){
            nserversp.set(i, nservers.get(qSET.get(i)));
        }
        PFQN.pfqnMVAReturn ret1 = PFQN.pfqn_mvams(lambda, Lp, Nchain, Zp, Matrix.ones(qSET.size(), 1), nserversp);
        Matrix Xchain = ret1.XN;
        Matrix Qpf = ret1.QN;
        /* Uchain = ret1.UN; -> not needed */
        double lG = ret1.lGN;
        for(int i = 0; i < qSET.size(); i++){
            for(int j = 0; j < Qchain.getNumCols(); j++){
                Qchain.set(qSET.get(i), j, Qpf.get(i, j));
            }
        }

        Matrix Q2 = new Matrix(infSET.size(), Qchain.getNumCols());
        Matrix Xchainrep = Xchain.repmat(infSET.size(), 1);
        for(int i = 0; i < infSET.size(); i++){
            for(int j = 0; j < Q2.getNumCols(); j++){
                Q2.set(i, j, Xchainrep.get(i, j) * STchain.get(infSET.get(i), j) *
                        Vchain.get(infSET.get(i), j));
            }
        }

        for(int i = 0; i < infSET.size(); i++){
            for(int j = 0; j < Qchain.getNumCols(); j++){
                Qchain.set(infSET.get(i), j, Q2.get(i, j));
            }
        }

        ArrayList<Integer> ccl = new ArrayList<>();
        for(int i = 0; i < Nchain.getNumCols(); i++){
            if(Double.isFinite(Nchain.get(i))){
                ccl.add(i);
            }
        }
        for(int r : rset){
            for(int k : infSET){
                Wchain.set(k, r, STchain.get(k, r));
            }
            for(int k : qSET){
                if(Double.isInfinite(nservers.get(k))){ // Infinite server
                    Wchain.set(k, r, STchain.get(k, r));
                } else {
                    if(Vchain.get(k, r) == 0 || Xchain.get(r) == 0){
                        Wchain.set(k, r, 0);
                    } else {
                        Wchain.set(k, r, Qchain.get(k, r) / (Xchain.get(r) * Vchain.get(k, r)));
                    }
                }
            }
        }
        for(int r : rset){
            if(Matrix.extractColumn(Wchain, r, null).elementSum() == 0){
                Xchain.set(r, 0);
            } else {
                if(Double.isInfinite(Nchain.get(r))){
                    Matrix vt = Matrix.extractColumn(Vchain, r, null);
                    Matrix wt = Matrix.extractColumn(Wchain, r, null);
                    C.set(r, vt.transpose().mult(wt).get(0));
                } else if(Nchain.get(r) == 0){
                    Xchain.set(r, 0);
                    C.set(r, 0);
                } else {
                    Matrix vt = Matrix.extractColumn(Vchain, r, null);
                    Matrix wt = Matrix.extractColumn(Wchain, r, null);
                    C.set(r, vt.transpose().mult(wt).get(0));
                    Xchain.set(r, Nchain.get(r) / C.get(r));
                }
            }

            for(int k = 0; k < M; k++){
                Qchain.set(k, r, Xchain.get(r) * Vchain.get(k, r) * Wchain.get(k, r));
                Tchain.set(k, r, Xchain.get(r) * Vchain.get(k, r));
            }
        }
        for(int k = 0; k < M; k++){
            for(int r : rset){
                if(Double.isInfinite(nservers.get(k))){ // infinite server
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r));
                } else {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r) / nservers.get(k));
                }
            }
        }
        for(int k = 0; k < M; k++){
            for(int r = 0; r < K; r++){
                if(Vchain.get(k, r) * STchain.get(k, r) > options.tol){
                    switch(sn.sched.get(sn.stations.get(k))){
                        case FCFS: case PS:
                            Matrix Urow = Matrix.extractRows(Uchain, k, k+1, null);
                            double UrowSum = Urow.elementSum();
                            if(UrowSum > 1 + options.tol){
                                Matrix Vrow = Matrix.extractRows(Vchain, k, k+1, null);
                                Matrix STrow = Matrix.extractRows(STchain, k, k+1, null);
                                Uchain.set(k, r, Maths.min(1, UrowSum) * Vchain.get(k, r) *
                                        STchain.get(k, r) * Xchain.get(r) / Vrow.elementMult(STrow, null).mult(Xchain.columnMajorOrder()).get(0));
                            }
                    }
                }
            }
        }
        Matrix Vs = null;
        for(int key : sn.nodevisits.keySet()){
            if(Vs == null){
                Vs = sn.nodevisits.get(key);
            } else {
                Vs = Vs.add(1, sn.nodevisits.get(key));
            }
        }
        ArrayList<Integer> sinks = new ArrayList<>();
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Sink){
                sinks.add(i);
            }
        }
        Matrix Vsink = new Matrix(sinks.size(), Vs.getNumCols());
        for(int i = 0; i < sinks.size(); i++){
            for(int j = 0; j < Vsink.getNumCols(); j++){
                Vsink.set(i, j, Vs.get(sinks.get(i), j));
            }
        }
        for(int r = 0; r < Nchain.getNumCols(); r++){
            if(Double.isInfinite(Nchain.get(r))){ // open classes
                Xchain.set(r, Vsink.get(r) / STchain.get((int) refstatchain.get(r), r));
            }
        }

        Matrix Rchain = new Matrix(Qchain.getNumRows(), Qchain.getNumCols());
        for(int i = 0; i < Rchain.getNumRows(); i++){
            for(int j = 0; j < Rchain.getNumCols(); j++){
                Rchain.set(i, j, Qchain.get(i, j) / Tchain.get(i, j));
            }
        }

        for(int i = 0; i < Xchain.getNumRows(); i++){
            for(int j = 0; j < Xchain.getNumCols(); j++){
                if(!Double.isFinite(Xchain.get(i, j)))
                    Xchain.set(i, j, 0);
            }
        }

        for(int i = 0; i < Uchain.getNumRows(); i++){
            for(int j = 0; j < Uchain.getNumCols(); j++){
                if(!Double.isFinite(Uchain.get(i, j)))
                    Uchain.set(i, j, 0);
            }
        }

        for(int i = 0; i < Qchain.getNumRows(); i++){
            for(int j = 0; j < Qchain.getNumCols(); j++){
                if(!Double.isFinite(Qchain.get(i, j)))
                    Qchain.set(i, j, 0);
            }
        }

        for(int i = 0; i < Rchain.getNumRows(); i++){
            for(int j = 0; j < Rchain.getNumCols(); j++){
                if(!Double.isFinite(Rchain.get(i, j)))
                    Rchain.set(i, j, 0);
            }
        }

        ArrayList<Integer> Nzero = new ArrayList<>();
        for(int i = 0; i < Nchain.getNumCols(); i++){
            if(Nchain.get(i) == 0)
                Nzero.add(i);
        }
        for(int j : Nzero){
            Xchain.set(j, 0);
            for(int i = 0; i < Uchain.getNumRows(); i++){
                Uchain.set(i, j, 0);
            }
            for(int i = 0; i < Qchain.getNumRows(); i++){
                Qchain.set(i, j, 0);
            }
            for(int i = 0; i < Rchain.getNumRows(); i++){
                Rchain.set(i, j, 0);
            }
            for(int i = 0; i < Tchain.getNumRows(); i++){
                Tchain.set(i, j, 0);
            }
            for(int i = 0; i < Wchain.getNumRows(); i++){
                Wchain.set(i, j, 0);
            }
        }

        SN.snDeaggregateChainResultsReturn ret2 = SN.snDeaggregateChainResults(sn, Lchain, null,
                STchain, Vchain, alpha, null, null, Rchain, Tchain, null,
                Xchain);

        SolverMVAResult res = new SolverMVAResult();
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
