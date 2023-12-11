package jline.solvers.mva.handlers;

import jline.lang.NetworkStruct;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class QNAHandler implements MVASolverHandler{
    public SolverMVAResult solve(NetworkStruct sn, SolverOptions options){
        SolverOptions.Config config = options.config;
        config.space_max = 1;

        int K = sn.nclasses;
        Matrix rt = sn.rt.clone();
        Matrix S = sn.rates.element_power(-1);
        Matrix scv = sn.scv.clone();
        scv.removeNaN();

        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = Matrix.cellsum(sn.visits);
        Matrix Q = new Matrix(M,K,M*K);
        Matrix QN_1 = Q.clone();
        for(int i=0;i<QN_1.numRows;i++){
            for (int j=0;j<QN_1.numCols;j++){
                QN_1.set(i,j,Double.POSITIVE_INFINITY);
            }
        }

        Matrix U = new Matrix(M,K,M*K);
        Matrix R = new Matrix(M,K,M*K);
        Matrix T = new Matrix(M,K,M*K);
        Matrix X = new Matrix(1,K,K);

        Matrix lambda = new Matrix(1,C,C);

        int it = 0;

        for(int i=0;i<sn.njobs.length();i++){
            if(Double.isFinite(sn.njobs.get(i))){
                throw new RuntimeException("QNA does not support closed classes.");
            }
        }

        Matrix a1 = new Matrix(M,K,M*K);
        Matrix a2 = new Matrix(M,K,M*K);
        Matrix d2 = new Matrix(M,1,M);
        Matrix f2 = new Matrix(M*K,M*K,(int)Math.pow(M*K,2));
        for(int i=0;i<M;i++){
            for(int j=0;j<M;j++){
                if(sn.nodetypes.get((int)sn.stationToNode.get(j))!= NodeType.Source){
                    for(int r=0;r<K;r++){
                        for(int s=0;s<K;s++){
                            if(rt.get(i*K+r,j*K+s)>0){
                                f2.set(i*K+r,j*K+s,1);
                            }
                        }
                    }
                }
            }
        }
        Map<Integer,Matrix> lambdas_inchain = new HashMap<>();
        Map<Integer,Matrix> scvs_inchain = new HashMap<>();
        Matrix d2c = new Matrix(1,C,C);
        int last_source_idx = 0;
        for(int c=0;c<C;c++){
            Matrix inchain = sn.inchain.get(c);
            int sourceIdx = (int)sn.refstat.get((int)inchain.get(0));
            last_source_idx = sourceIdx;
            lambdas_inchain.put(c,new Matrix(1, inchain.length(), inchain.length()));
            for(int i=0;i< inchain.length();i++){
                lambdas_inchain.get(c).set(0,i,sn.rates.get(sourceIdx,(int)inchain.get(i)));
            }
            scvs_inchain.put(c,new Matrix(1, inchain.length(), inchain.length()));
            for (int i=0;i<inchain.length();i++){
                scvs_inchain.get(c).set(0,i,scv.get(sourceIdx,(int)inchain.get(i)));
            }
            Matrix lambdas_inchain_C = lambdas_inchain.get(c).clone();
            lambdas_inchain_C.removeINF();
            lambda.set(c,lambdas_inchain_C.elementSum());
            d2c.set(c,qna_superpos(lambdas_inchain.get(c),scvs_inchain.get(c)));
            boolean openChain = false;
            for(int i=0;i<inchain.length();i++){
                if(Double.isInfinite(sn.njobs.get((int)inchain.get(i)))){
                    openChain = true;
                }
            }
            if(openChain){
                for(int i=0;i< inchain.length();i++){
                    T.set(sourceIdx,(int)inchain.get(i),lambdas_inchain.get(c).get(i));
                }
            }
        }
        d2.set(last_source_idx,Matrix.extractRows(d2c,last_source_idx,last_source_idx+1,null).mult(lambda.transpose()).get(0)/lambda.elementSum());
        Matrix Q_diff = Q.add(-1,QN_1);
        Q_diff.abs();
        while (Q_diff.elementMax()> options.iter_tol&&it<= options.iter_max){
            it = it+1;
            QN_1 = Q.clone();

            if(it==1){
                for(int c=0;c<C;c++){
                    Matrix inchain = sn.inchain.get(c);
                    for(int m=0; m<M;m++){
                        for(int i=0;i<inchain.length();i++){
                            T.set(m,(int)inchain.get(i),V.get(m,(int) inchain.get(i))*lambda.get(c));
                        }
                    }
                }
            }
            for(int i=0;i<M;i++){
                for(int j=0;j<K;j++){
                    a1.set(i,j,0);
                    a2.set(i,j,0);
                }
                double lambda_i = T.sumRows(i);
                for (int j=0;j<M;j++){
                    for(int r=0;r<K;r++){
                        for(int s=0;s<K;s++){
                            a1.set(i,r,a1.get(i,r)+T.get(j,s)*rt.get(j*K+s,i*K+r));
                            a2.set(i,r,a2.get(i,r)+1/lambda_i*f2.get(j*K+s,i*K+r)*T.get(j,s)*rt.get(j*K+s,i*K+r));
                        }
                    }
                }
            }

            for(int ind = 0;ind<I;ind++){
                if(sn.isstation.get(ind)==1){
                    int ist = (int) sn.nodeToStation.get(ind);
                    if(sn.nodetypes.get(ind)!=NodeType.Join){
                        if(sn.sched.get(sn.stations.get(ist))==SchedStrategy.INF){
                            for(int i=0;i<M;i++){
                                for(int r=0;r<K;r++){
                                    for(int s=0;s<K;s++){
                                        d2.set(ist,s,a2.get(ist,s));
                                    }
                                }
                            }
                            for(int c=0;c<C;c++){
                                Matrix inchain = sn.inchain.get(c);
                                for(int k1=0;k1<inchain.length();k1++){
                                    int k = (int) inchain.get(k1);
                                    T.set(ist,k,a1.get(ist,k));
                                    U.set(ist,k,S.get(ist,k)*T.get(ist,k));
                                    Q.set(ist,k,T.get(ist,k)*S.get(ist,k)*V.get(ist,k));
                                    R.set(ist,k,Q.get(ist,k)/T.get(ist,k));
                                }
                            }
                        }else if(sn.sched.get(sn.stations.get(ist))==SchedStrategy.FCFS){
                            Matrix mu_ist = new Matrix(1,K,K);
                            for(int i=0;i<K;i++) {
                                mu_ist.set(i, sn.rates.get(ist, i));
                            }
                            mu_ist.removeNaN();
                            Matrix rho_ist_class = new Matrix(1,K,K);
                            for(int i=0;i<K;i++){
                                    rho_ist_class.set(i,a1.get(ist,i)/(GlobalConstants.FineTol+sn.rates.get(ist,i)));
                            }
                            rho_ist_class.removeNaN();
                            double lambda_ist = a1.sumRows(ist);
                            int mi = (int) sn.nservers.get(ist);
                            double rho_ist = rho_ist_class.elementSum()/mi;
                            double c2 = 0;
                            if(rho_ist<1-options.tol){
                                for (int k=0;k<K;k++){
                                    double alpha_mi;
                                    if(rho_ist>0.7){
                                        alpha_mi = (Math.pow(rho_ist,mi)+rho_ist)/2;
                                    }else {
                                        alpha_mi = Math.pow(rho_ist,(mi+1)/2.0);
                                    }
                                    double mubar = lambda_ist/rho_ist;
                                    c2 = -1;
                                    for (int r=0;r<K;r++){
                                        if(mu_ist.get(r)>0){
                                            c2 = c2+a1.get(ist,r)/lambda_ist*Math.pow(mubar/mi/mu_ist.get(r),2)*(scv.get(ist,r)+1);
                                        }
                                    }
                                    double Wiq = (alpha_mi/mubar)*1/(1-rho_ist)*(a2.sumRows(ist)+c2)/2;
                                    Q.set(ist,k,a1.get(ist,k)/mu_ist.get(k)+a1.get(ist,k)*Wiq);
                                }
                                d2.set(ist,1+Math.pow(rho_ist,2)*(c2-1)/Math.sqrt(mi)+(1-Math.pow(rho_ist,2))*(a2.sumRows(ist)-1));
                            }else {
                                for (int k=0;k<K;k++) {
                                    Q.set(ist,k,sn.njobs.get(k));
                                }
                                d2.set(ist,1);
                            }
                            for(int k=0;k<K;k++){
                                T.set(ist,k,a1.get(ist,k));
                                U.set(ist,k,T.get(ist,k)*S.get(ist,k)/sn.nservers.get(ist));
                                R.set(ist,k,Q.get(ist,k)/T.get(ist,k));
                            }
                        }
                    }
                }else {
                    if(sn.nodetypes.get(ind)==NodeType.Fork){
                        throw new RuntimeException("Fork nodes not supported yet by QNA solver.");
                    }
                }
            }

            for(int i=0;i<M;i++){
                for (int j=0;j<M;j++){
                    if(sn.nodetypes.get((int) sn.stationToNode.get(j))!=NodeType.Source){
                        for(int r=0;r<K;r++){
                            for(int s=0;s<K;s++){
                                if(rt.get(i*K+r,j*K+s)>0){
                                    f2.set(i*K+r,j*K+s,1+rt.get(i*K+r,j*K+s)*(d2.get(i)-1));
                                }
                            }
                        }
                    }
                }
            }
            Q_diff = Q.add(-1,QN_1);
            Q_diff.abs();
        }
        SolverMVAResult result = new SolverMVAResult();
        Matrix CN = R.sumCols();
        Q.abs();
        Q.removeNaN();
        U.removeNaN();
        R.removeNaN();
        CN.removeNaN();
        X.removeNaN();
        result.XN = X;
        result.QN = Q;
        result.CN = CN;
        result.UN = U;
        result.RN = R;
        result.TN = T;
        result.iter = it;
        result.logNormConstAggr = 0;
        return result;
    }
    public static double qna_superpos(Matrix lambda, Matrix a2){
        List<Integer> lambda_finite_idx = new ArrayList<>();
        for(int i=0;i< lambda.length();i++){
            if(Double.isFinite(lambda.get(i))){
                lambda_finite_idx.add(i);
            }
        }
        Matrix a2_new = new Matrix(1, lambda_finite_idx.size(), lambda_finite_idx.size());
        for(int i=0;i<lambda_finite_idx.size();i++){
            a2_new.set(i,a2.get(lambda_finite_idx.get(i)));
        }
        Matrix lambda_new = new Matrix(1, lambda_finite_idx.size(), lambda_finite_idx.size());
        for(int i=0;i<lambda_finite_idx.size();i++){
            lambda_new.set(i,lambda.get(lambda_finite_idx.get(i)));
        }
        return a2_new.mult(lambda_new.transpose()).get(0)/lambda_new.elementSum();
    }
}
