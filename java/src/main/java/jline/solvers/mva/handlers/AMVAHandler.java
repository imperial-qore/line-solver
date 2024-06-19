package jline.solvers.mva.handlers;

import jline.api.PFQN;
import jline.api.SN;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;
import org.apache.commons.lang3.NotImplementedException;

import java.util.ArrayList;

import static jline.api.PFQN.*;

/**
 * Handler for the solver_amva function.
 */
public class AMVAHandler implements MVASolverHandler {
    @Override
    public SolverMVAResult solve(NetworkStruct sn, SolverOptions options) {
        SN.snGetDemandsChainReturn chainReturn = SN.snGetDemandsChain(sn);
        Matrix Lchain = chainReturn.Lchain;
        Matrix STchain = chainReturn.STchain;
        Matrix Vchain = chainReturn.Vchain;
        Matrix alpha = chainReturn.alpha;
        Matrix Nchain = chainReturn.Nchain;
        Matrix SCVchain = chainReturn.SCVchain;
        Matrix refstatchain = chainReturn.refstatchain;
        if(options.config.np_priority == null){
            options.config.np_priority = "default";
        }
        if(options.config.multiserver == null){
            options.config.multiserver = "default";
        }
        if(options.config.highvar == null){
            options.config.highvar = "default";
        }
        switch(options.method){
            case "amva.qli":
                options.method = "qli";
                break;
            case "amva.qd": case "amva.qdamva": case "qdamva":
                options.method = "qd";
                break;
            case "amva.aql":
                options.method = "aql";
                break;
            case "amva.qdaql":
                options.method = "qdaql";
                break;
            case "amva.lin":
                options.method = "lin";
                break;
            case "amva.gflin":
                options.method = "gflin";
                break;
            case "amva.egflin":
                options.method = "egflin";
                break;
            case "amva.qdlin":
                options.method = "qdlin";
                break;
            case "amva.fli":
                options.method = "fli";
                break;
            case "amva.bs":
                options.method = "bs";
                break;
            case "default":
                double NchainSum = 0;
                boolean flag = false;
                for(int i = 0; i < Nchain.getNumRows(); i++){
                    for(int j = 0; j < Nchain.getNumCols(); j++){
                        NchainSum += Nchain.get(i, j);
                        if(Nchain.get(i, j) < 1)
                            flag = true;
                    }
                }
                if(NchainSum <= 2 || flag){
                    options.method = "qd"; // changing to bs degrades accuracy
                } else {
                    options.method = "egflin";
                    for (int i=0; i < sn.nstations; i++) {
                        if (sn.nservers.get(i,0)>1 && sn.nservers.get(i,0)<Integer.MAX_VALUE) {
                            // if multi-server
                            options.method = "lin"; // lin seems way worse than aql in test_LQN_8.xml
                        }
                    }
                }
                break;
        }

        // trivial models
        if(SN.snHasHomogeneousScheduling(sn, SchedStrategy.INF)){
            options.config.multiserver = "default";
            return new AMVALDHandler().solve(sn, options);
        }
        ArrayList<Integer> queueIdx = new ArrayList<>();
        ArrayList<Integer> delayIdx = new ArrayList<>();
        ArrayList<Integer> sourceIdx = new ArrayList<>();
        for(int i = 0; i < sn.nodetypes.size(); i++){
            switch (sn.nodetypes.get(i)){
                case Source:
                    sourceIdx.add(i);
                    break;
                case Queue:
                    queueIdx.add(i);
                    break;
                case Delay:
                    delayIdx.add(i);
                    break;
            }
        }
        // Run amva method
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = new Matrix(M, C);
        for(int s : sourceIdx){
            int i = (int) sn.nodeToStation.get(s);
            for(int c = 0; c < sn.nchains; c++){
                double rateSum = 0;
                for(int sp : sourceIdx){
                    for(int r = 0; r < sn.chains.getNumCols(); r++){
                        if(sn.chains.get(c, r) == 0){
                            continue;
                        }
                        if(!Double.isNaN(sn.rates.get((int) sn.nodeToStation.get(sp), r))){
                            rateSum += sn.rates.get((int) sn.nodeToStation.get(sp), r);
                        }
                    }
                }
                if(rateSum > 0){
                    V.set(i, c, 1);
                }
            }
        }
        Matrix Q = new Matrix(M, C);
        Matrix U = new Matrix(M, C);
        Matrix X = null;
        int totiter = 0;
        if(SN.snHasProductFormExceptMultiClassHeterExpFCFS(sn) && !SN.snHasLoadDependence(sn) &&
                (!SN.snHasOpenClasses(sn) || (SN.snHasProductForm(sn) && SN.snHasOpenClasses(sn) &&
                        options.method.equals("lin")))){
            SN.snGetProductFormChainParamsReturn ret = SN.snGetProductFormChainParams(sn);
            Matrix lambda = ret.lambda;
            Matrix L0 = new Matrix(ret.D);
            Matrix L = ret.D;
            Matrix N = ret.N;
            Matrix Z = ret.Z;
            Matrix Z0 = new Matrix(ret.Z);
            Matrix nservers = ret.S;
            int vIdx = 0;
            for(int i = 0; i < sn.nodetypes.size(); i++){
                if(sn.nodetypes.get(i) != NodeType.Queue && sn.nodetypes.get(i) != NodeType.Delay){
                    continue;
                }
                for(int j = 0; j < V.getNumCols(); j++){
                    V.set((int) sn.nodeToStation.get(i), j, ret.V.get(vIdx, j));
                }
                vIdx++;
            }
            switch(options.config.multiserver){
                case "default": case "seidmann":
                    Matrix nserversRep = nservers.columnMajorOrder().repmat(1, C);
                    // apply seidmann
                    for(int i = 0; i < L.getNumRows(); i++){
                        for(int j = 0; j < L.getNumCols(); j++){
                            L.set(i, j, L.get(i, j) / nserversRep.get(i, j));
                        }
                    }
                    for(int j = 0; j < L.getNumRows(); j++){
                        for(int k = 0; k < Z.getNumCols(); k++){
                            Z.set(0, k, Z.get(0, k) + L0.get(j, k) * (nservers.get(j) - 1) / nservers.get(j));
                        }
                    }
                    break;
                case "softmin":
                    return new AMVALDHandler().solve(sn, options);
            }
            switch(options.method){
                case "sqni":
                    throw new NotImplementedException("SQNI not implemented in AMVAHandler"); // TODO
                case "bs":
                    SchedStrategy[] schd = new SchedStrategy[queueIdx.size()];
                    for(int i = 0; i < queueIdx.size(); i++){
                        schd[i] =sn.sched.get(sn.stations.get(queueIdx.get(i)));
                    }
                    pfqnBSReturn bsret = PFQN.pfqn_bs(L, N, Z, options.tol, options.iter_max, null, schd);
                    X = bsret.XN;
                    int iResult = 0;
                    for(int i : queueIdx){
                        for(int j = 0; j < C; j++){
                            Q.set((int) sn.nodeToStation.get(i), j, bsret.QN.get(iResult, j));
                            U.set((int) sn.nodeToStation.get(i), j, bsret.UN.get(iResult, j));
                        }
                        iResult++;
                    }
                    totiter = bsret.it;
                    break;
                case "aql":
                    throw new NotImplementedException("AQL not implemented in AMVAHandler"); // TODO
                case "lin": case "gflin": case "egflin":
                    if(nservers.elementMax() == 1){
                        SchedStrategy[] schdi = new SchedStrategy[queueIdx.size() + delayIdx.size()];
                        int idx = 0;
                        for(int i = 0; i < sn.nodetypes.size(); i++){
                            if(sn.nodetypes.get(i) != NodeType.Queue && sn.nodetypes.get(i) != NodeType.Delay){
                                continue;
                            }
                            schdi[idx] = sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(i)));
                            idx++;
                        }

                        pfqnAMVAReturn res = pfqn_linearizermx(lambda, L, N, Z, nservers, schdi, options.tol, options.iter_max, options.method);

                        int iRes = 0;
                        for(int i : queueIdx){
                            for(int j = 0; j < C; j++){
                                Q.set((int) sn.nodeToStation.get(i), j, res.Q.get(iRes, j));
                                U.set((int) sn.nodeToStation.get(i), j, res.U.get(iRes, j));
                            }
                            iRes++;
                        }
                        X = res.X;
                        totiter = res.totiter;
                        break;
                    } else {
                        switch(options.config.multiserver){
                            case "conway":
                                throw new NotImplementedException("Conway not implemented in AMVAHandler"); // TODO
                            case "erlang":
                                throw new NotImplementedException("Erlang not implemented in AMVAHandler"); // TODO
                            case "krzesinski":
                                SchedStrategy[] schdi = new SchedStrategy[queueIdx.size() + delayIdx.size()];
                                int idx = 0;
                                for(int i = 0; i < sn.nodetypes.size(); i++){
                                    if(sn.nodetypes.get(i) != NodeType.Queue && sn.nodetypes.get(i) != NodeType.Delay){
                                        continue;
                                    }
                                    schdi[idx] = sn.sched.get(sn.stations.get((int)sn.nodeToStation.get(i)));
                                    idx++;
                                }
                                pfqnAMVAReturn res = pfqn_linearizermx(lambda, L, N, Z, nservers, schdi, options.tol, options.iter_max, null);
                                int iRes = 0;
                                for(int i : queueIdx){
                                    for(int j = 0; j < C; j++){
                                        Q.set((int) sn.nodeToStation.get(i), j, res.Q.get(iRes, j));
                                        U.set((int) sn.nodeToStation.get(i), j, res.U.get(iRes, j));
                                    }
                                    iRes++;
                                }
                                X = res.X;
                                totiter = res.totiter;
                                break;
                            case "default": case "softmin": case "seidmann":
                                return new AMVALDHandler().solve(sn, options);
                        }
                    }
                    break;
                default:
                    switch(options.config.multiserver){
                        case "conway": case "erlang": case "krzesinski":
                            options.config.multiserver = "default";
                    }
                    return new AMVALDHandler().solve(sn, options);
            }

            // Compute performance at delay, then unapply seidmann if needed
            for(int i = 0; i < Z0.getNumRows(); i++){
                Matrix mult = X.repmat(delayIdx.size(), 1).elementMult(Z, null);
                int zidx = 0;
                for(int d : delayIdx){
                    for(int j = 0; j < Q.getNumCols(); j++){
                        Q.set((int) sn.nodeToStation.get(d), j, mult.get(zidx, j));
                        U.set((int) sn.nodeToStation.get(d), j, mult.get(zidx, j));
                    }
                    zidx++;
                }
                switch(options.config.multiserver){
                    case "default": case "seidmann":
                        for(int j = 0; j < L.getNumRows(); j++){
                            if(i == 0 && nservers.get(j) > 1){
                                // un-apply seidmann from first delay and move it to the origin queue
                                for(int k = 0; k <= j; k++){
                                    if(k >= queueIdx.size())
                                        break;
                                    int jq = queueIdx.get(k);
                                    for(int l = 0; l < Q.getNumCols(); l++){
                                        Q.set(jq, l, Q.get(jq, l) + (L0.get(j, l) *
                                                (nservers.get(j) - 1) / nservers.get(j)) * X.get(l));
                                    }
                                }
                            }
                        }
                }
            }
            Matrix T = V.elementMult(X.repmat(M,1), null);
            Matrix R = new Matrix(Q.getNumRows(), Q.getNumCols());
            for(int i = 0; i < R.getNumRows(); i++){
                for(int j = 0; j < R.getNumCols(); j++){
                    R.set(i, j, Q.get(i, j) / T.get(i, j));
                }
            }
            Matrix Cm = new Matrix(N.getNumRows(), N.getNumCols());
            for(int i = 0; i < Cm.getNumRows(); i++){
                for(int j = 0; j < Cm.getNumCols(); j++){
                    Cm.set(i, j, N.get(i, j) / X.get(i, j) - Z.get(i, j));
                }
            }
            double lG = Double.NaN;
            SolverMVAResult result = new SolverMVAResult();
            result.QN = Q;
            result.UN = U;
            result.RN = R;
            result.TN = T;
            result.CN = Cm;
            result.XN = X;
            result.logNormConstAggr = lG;
            result.iter = totiter;
            if(SN.snHasClassSwitching(sn)){
                SN.snDeaggregateChainResultsReturn ret1 = SN.snDeaggregateChainResults(sn, Lchain, null,
                        STchain, Vchain, alpha, null, null, R, T, null, X);
                result.QN = ret1.Q;
                result.UN = ret1.U;
                result.RN = ret1.R;
                result.TN = ret1.T;
                result.CN = ret1.C;
                result.XN = ret1.X;
            }
            return result;
        } else {
            switch(options.config.multiserver){
                case "conway": case "erlang": case "krzesinski":
                    options.config.multiserver = "default";
            }
            return new AMVALDHandler().solve(sn, options);
        }

    }
}
