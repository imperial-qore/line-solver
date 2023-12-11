package jline.solvers.mam;


import jline.api.DTMC;
import jline.api.SN;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.*;
import jline.lang.distributions.APH;
import jline.lang.nodes.Station;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.Matrix;

import java.util.*;

import static jline.api.MAM.*;
import static jline.api.NPFQN.*;
import static jline.lib.KPCToolbox.*;
import static jline.lib.M3A.*;

import static jline.lib.thirdparty.BUTOOLS.*;



public class SolverMAM extends NetworkSolver {
    public SolverMAM(Network model){
        super(model,"SolverMAM",new SolverOptions(SolverType.MAM));
    }

    public SolverMAM(Network model, SolverOptions options){
        super(model,"SolverMAM",options);
    }

    public NetworkStruct getStruct(){
        return model.getStruct(true);
    }

    static public String listValidMethods(){
        return "default,dec.source,dec.mmap,dec.poisson";
    }

    @Override
    public void runAnalyzer(){
        double start = System.currentTimeMillis();
        NetworkStruct sn = getStruct();
        result = solver_mam_analyser(sn);
        int M = sn.nstations;
        int R = sn.nclasses;
        if(result.TN.length()>0){
            result.AN = new Matrix(M,R,M*R);
            for(int i=0;i<M;i++){
                for(int j=0;j<M;j++){
                    for(int k=0;k<R;k++){
                        for(int r=0;r<R;r++){
                            result.AN.set(i,k,result.AN.get(i,k)+result.TN.get(j,r)*sn.rt.get(j*R+r,i*R+k));
                        }
                    }
                }
            }
        }else {
            Matrix AN = new Matrix(0,0,0);
        }

        //setAvgResults();

        double finish = System.currentTimeMillis();
        result.runtime = (finish - start)/1000;
    }

    public SolverMAMResult solver_mam_analyser(NetworkStruct sn) {
        long start = System.nanoTime();
        options.config.merge = "super";
        options.config.compress = "mixture.order1";
        options.config.space_max = 128;
        SolverMAMResult result = new SolverMAMResult();
        if (options.method.equals("dec.mmap")) {
            result = (SolverMAMResult) solver_mam(sn);
        } else if (options.method.equals("default") || options.method.equals("dec.source")) {
            result = (SolverMAMResult) solver_mam_basic(sn);
        } else if (options.method.equals("poisson")) {
            options.config.space_max = 1;
            result = (SolverMAMResult) solver_mam_basic(sn);
        }else if(options.method.equals("qnamam")){
            boolean closed = false;
            boolean open = false;
            for(int i=0;i<sn.nclasses;i++){
                if(Double.isInfinite(sn.njobs.get(i))){
                    open = true;
                    break;
                }
            }
            for(int i=0;i<sn.nclasses;i++){
                if(Double.isFinite(sn.njobs.get(i))){
                    closed = true;
                    break;
                }
            }
            if(closed&&!open){
                result = (SolverMAMResult) solver_qna_mam_closed(sn);
            }else if (open &&!closed) {
                result = (SolverMAMResult) solver_qna_mam(sn);
            }else {
                throw new RuntimeException("Unsupport model for QNAMAM");
            }
        } else {
            throw new RuntimeException("Unknown method");
        }

        for(int i=0;i<sn.nstations;i++){
            if(sn.sched.get(sn.stations.get(i))==SchedStrategy.EXT){
                for(int j=0;j<result.TN.numCols;j++){
                    result.TN.set(i,j,sn.rates.get(i,j));
                }
            }
        }

        result.QN.removeNaN();
        result.CN.removeNaN();
        result.RN.removeNaN();
        result.UN.removeNaN();
        result.XN.removeNaN();
        result.TN.removeNaN();
        long finish = System.nanoTime();
        result.runtime = finish - start;
        return  result;
    }

    public SolverResult solver_qna_mam(NetworkStruct sn){
        SolverOptions.Config config = options.config;
        config.space_max = 1;

        int K = sn.nclasses;
        Matrix rt = sn.rt.clone();
        Matrix S = sn.rates.element_power(-1);
        Matrix scv = sn.scv.clone();
        scv.removeNaN();
        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = Matrix.cellsum(sn.visits);
        Matrix Q = new Matrix(M,K,M*K);
        Map<Integer,Map<Integer,Matrix>> pie = new HashMap<>();
        Map<Integer,Map<Integer,Matrix>> DO = new HashMap<>();


        Matrix U = new Matrix(M,K,M*K);
        Matrix R = new Matrix(M,K,M*K);
        Matrix T = new Matrix(M,K,M*K);
        Matrix X = new Matrix(1,K,K);
        for(int ist=0; ist<M;ist++){
            if(sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS||sn.sched.get(sn.stations.get(ist))==SchedStrategy.HOL||sn.sched.get(sn.stations.get(ist))==SchedStrategy.PS){
                pie.put(ist,new HashMap<Integer,Matrix>());
                DO.put(ist,new HashMap<Integer,Matrix>());
                for(int k = 0;k<K;k++){
                    PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k),map_scale(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1),S.get(ist,k)/sn.nservers.get(ist)));
                    pie.get(ist).put(k,map_pie(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));
                    DO.get(ist).put(k,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0));
                }
            }
        }

        Matrix lambda = new Matrix(1,C,C);

        int it = 0;

        Matrix a1 = new Matrix(M,K,M*K);
        Matrix a2 = new Matrix(M,K,M*K);
        Matrix a1_1 = a1.elementIncrease(Double.POSITIVE_INFINITY);

        Matrix a2_1 = a2.elementIncrease(Double.POSITIVE_INFINITY);

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

            for(int i=0;i< inchain.length();i++){
                T.set(sourceIdx,(int)inchain.get(i),lambdas_inchain.get(c).get(i));
            }
        }
        d2.set(last_source_idx,Matrix.extractRows(d2c,last_source_idx,last_source_idx+1,null).mult(lambda.transpose()).get(0)/lambda.elementSum());
        Matrix a1_diff = a1.add(-1,a1_1);
        a1_diff.abs();
        Matrix a2_diff = a2.add(-1,a2_1);
        a2_diff.abs();

        while ((a1_diff.elementMax()> options.iter_tol||a2_diff.elementMax()>options.iter_tol)&&it<= options.iter_max){
            it = it+1;
            a1_1 = a1.clone();
            a2_1 = a2.clone();

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
                            a2.set(i,r,a2.get(i,r)+(1/lambda_i)*f2.get(j*K+s,i*K+r)*T.get(j,s)*rt.get(j*K+s,i*K+r));
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


                                    double mubar = lambda_ist/rho_ist;
                                    c2 = -1;
                                    for (int r=0;r<K;r++){
                                        if(mu_ist.get(r)>0){
                                            c2 = c2+a1.get(ist,r)/lambda_ist*Math.pow(mubar/mi/mu_ist.get(r),2)*(scv.get(ist,r)+1);
                                        }
                                    }


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
                            }
                        }
                    }
                }else {
                    if(sn.nodetypes.get(ind)==NodeType.Fork){
                        // TODO: not implemented
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
            a1_diff = a1.add(-1,a1_1);
            a1_diff.abs();
            a2_diff = a2.add(-1,a2_1);
            a2_diff.abs();
        }

        for(int ind = 0;ind<I;ind++) {
            if (sn.isstation.get(ind) == 1) {
                int ist = (int) sn.nodeToStation.get(ind);
                if(sn.sched.get(sn.stations.get(ist))==SchedStrategy.FCFS){
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
                    if(rho_ist<1-options.tol) {
                        Map<Integer,Matrix> arri_class_anx = new HashMap<>();
                        Map<Integer,Matrix> arri_class = new HashMap<>();
                        Map<Integer,Matrix> arri_node = new HashMap<>();
                        for (int k = 0; k < K; k++) {
                            if (a1.get(ist, k) == 0) {
                                arri_class_anx = map_exponential(Double.POSITIVE_INFINITY);
                            } else {
                                arri_class_anx = APH.fitMeanAndSCV(1/a1.get(ist,k),a2.get(ist,k)).getRepres();
                            }
                            arri_class.put(0,arri_class_anx.get(0));
                            arri_class.put(1,arri_class_anx.get(1));
                            arri_class.put(2,arri_class_anx.get(1));
                            if (k == 0) {
                                arri_node.put(0,arri_class.get(0));
                                arri_node.put(1,arri_class.get(1));
                                arri_node.put(2,arri_class.get(2));
                            } else {
                                arri_node = mmap_super(arri_node,arri_class);
                            }


                        }
                        Map<Integer, Matrix> Qret = new HashMap<>();
                        Qret = MMAPPH1FCFS(mmap_shorten(arri_node), pie.get(ist), DO.get(ist), 1, null, null, null, false, false, null, null).get("ncMoms");
                        for (int i=0;i<Qret.size();i++){
                            Q.set(ist,i,Qret.get(i).get(0));
                        }
                    }else{
                        for(int k=0;k<K;k++){
                            Q.set(ist,k,sn.njobs.get(k));
                        }
                    }
                    for(int k=0;k<K;k++){
                        R.set(ist,k,Q.get(ist,k)/T.get(ist,k));
                    }

                }
            }
        }

        SolverMAMResult result = new SolverMAMResult();
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
    public SolverMAMResult solver_mam_basic(NetworkStruct sn){
        SolverOptions.Config conifg = options.config;
        double tol = options.tol;

        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        int I =  sn.nnodes;
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix N = sn.njobs.transpose();
        Matrix V = Matrix.cellsum(sn.visits);
        Matrix S = new Matrix(sn.rates);
        sn.rates.divide(1,S,false);
        Matrix Lchain = SN.snGetDemandsChain(sn).Lchain;

        Matrix QN = new Matrix(M,K,M*K);
        Matrix UN = new Matrix(M,K,M*K);
        Matrix RN = new Matrix(M,K,M*K);
        Matrix TN = new Matrix(M,K,M*K);
        Matrix WN = new Matrix(M,K,M*K);
        Matrix AN = new Matrix(M,K,M*K);
        Matrix CN = new Matrix(1,K,K);
        Matrix XN = new Matrix(1,K,K);

        Map<Integer,Map<Integer,Matrix>> pie = new HashMap<>();
        Map<Integer,Map<Integer,Matrix>> DO = new HashMap<>();

        Matrix lambda = new Matrix(1,C,C);
        Map<Integer, Map<Integer,Matrix>> chainSysArrivals = new HashMap<>();
        Matrix TN_1 = new Matrix(M,K,M*K);
        for(int i =0;i<M;i++){
            for(int j=0;j<K;j++){
                TN_1.set(i,j,Double.POSITIVE_INFINITY);
            }
        }

        int it =0;

        for(int ist=0; ist<M;ist++){
            if(sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS||sn.sched.get(sn.stations.get(ist))==SchedStrategy.HOL||sn.sched.get(sn.stations.get(ist))==SchedStrategy.PS){
                pie.put(ist,new HashMap<Integer,Matrix>());
                DO.put(ist,new HashMap<Integer,Matrix>());
                for(int k = 0;k<K;k++){
                    PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k),map_scale(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1),S.get(ist,k)/sn.nservers.get(ist)));
                    pie.get(ist).put(k,map_pie(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));
                    DO.get(ist).put(k,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0));
                }
            }
        }

        boolean isOpen = false;
        boolean isClosed = false;

        for(int i=0; i<sn.njobs.numRows;i++){
            for(int j=0; j< sn.njobs.numCols;j++){
                if(!Double.isFinite(sn.njobs.get(i,j))){
                    isOpen = true;
                    break;
                }
            }
        }

        for(int i=0; i<sn.njobs.numRows;i++){
            for(int j=0; j< sn.njobs.numCols;j++){
                if(Double.isFinite(sn.njobs.get(i,j))){
                    isClosed = true;
                    break;
                }
            }
        }

        boolean isMixed = isOpen && isClosed;

        if(isMixed){
            //todo line 65-68
        }

        HashMap<Integer, Matrix> lambdas_inchain = new HashMap<>(C);
        for(int c=0;c<C;c++){
            Matrix inchain = sn.inchain.get(c);
            lambdas_inchain.put(c, new Matrix(inchain.numCols,1,inchain.numCols));
            int ist = (int)sn.refstat.get((int)inchain.get(0,0),0);
            if(!PH.containsKey(sn.stations.get(ist))){
                PH.put(sn.stations.get(ist),new HashMap<>());
            }
            for(int j=0; j<inchain.numCols;j++) {
                lambdas_inchain.get(c).set(j,0,sn.rates.get(ist,(int)inchain.get(0,j)));
            }
            double sum = 0;
            for(int k=0;k<lambdas_inchain.get(c).numRows;k++){
                if(Double.isFinite(lambdas_inchain.get(c).get(k,0))){
                    sum += lambdas_inchain.get(c).get(k,0);
                }
            }
            lambda.set(0,c,sum);
            boolean openChain = false;
            for(int j=0;j<inchain.numCols;j++){
                if(Double.isInfinite(sn.njobs.get((int)inchain.get(0,j)))){
                    openChain = true;
                    break;
                }
            }
            if(openChain){
                for(int k=0;k<K;k++){
                    // if the matrix is NaN, then matrix has only one value and its is Double.NaN
                    if(Double.isNaN(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0).get(0))){
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k),map_exponential(Double.POSITIVE_INFINITY));
                    }
                }
                double k = inchain.get(0);
                chainSysArrivals.put(c,new HashMap<>());
                chainSysArrivals.get(c).put(0,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0));
                chainSysArrivals.get(c).put(1,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                chainSysArrivals.get(c).put(2,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                for(int ki=1;ki<inchain.length();ki++){
                    k = inchain.get(ki);
                    if(Double.isNaN(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0).get(0))){
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get((int) k),map_exponential(Double.POSITIVE_INFINITY));
                    }
                    Map<Integer,Map<Integer,Matrix>> MMAPS = new HashMap<>();
                    MMAPS.put(0,chainSysArrivals.get(c));
                    MMAPS.put(1,new HashMap<>());
                    MMAPS.get(1).put(0,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0));
                    MMAPS.get(1).put(1,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                    MMAPS.get(1).put(2,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                    chainSysArrivals.put(c,mmap_super_safe(MMAPS,conifg.space_max,"default"));
                }
                for(int i=0;i<inchain.numCols;i++){
                    TN.set(ist, (int) inchain.get(0,i),lambdas_inchain.get(c).get(i,0));
                }
            }
        }

        Matrix sd = new Matrix(sn.nservers.numRows,sn.nservers.numCols,sn.nservers.nz_length);
        for(int i=0; i<sn.nservers.numRows;i++){
            if(Double.isFinite(sn.nservers.get(i,0))){
                sd.set(i,0,1);
            }
        }

        Matrix dif_matrix  = TN.add(-1,TN_1);
        double dif = Math.abs(dif_matrix.elementMaxAbs());
        while(dif> tol && it<= options.iter_max){
            it++;
            TN_1 = TN.clone();
            Matrix Umax_matrix = new Matrix(M,1,M);
            for(int i=0;i<M;i++){
                if(sd.get(i,0)==0){
                    Umax_matrix.set(i,0,-Double.POSITIVE_INFINITY);
                }else {
                    double sum = 0;
                    for(int j=0;j<K;j++){
                        sum+= UN.get(i,j);
                    }
                    Umax_matrix.set(i,0,sum);
                }
            }

            double Umax = Umax_matrix.elementMax();
            if(Umax>1+ options.tol){
                lambda.divide(Umax,lambda,true);
            }else {
                for(int c=0;c<C;c++){
                    Matrix inchain = sn.inchain.get(c);
                    boolean openChain = false;
                    for(int j=0;j<inchain.numCols;j++){
                        if(Double.isInfinite(sn.njobs.get((int)inchain.get(0,j)))){
                            openChain = true;
                            break;
                        }
                    }
                    if(!openChain){
                        double Nc = 0;
                        for(int i=0; i< inchain.length();i++){
                            Nc = Nc + sn.njobs.get((int)inchain.get(i));
                        }
                        Matrix QN_omitnan = new Matrix(QN);
                        QN_omitnan.removeNaN();
                        Matrix QN_col_sum = new Matrix(inchain.length(),1,inchain.length());
                        double sum =0;
                        for(int i=0;i<inchain.length();i++){
                            QN_col_sum.set(i,0,QN_omitnan.sumCols((int)inchain.get(i)));
                        }
                        double QNc = QN_col_sum.elementSum();
                        QNc = Math.max(options.tol, QNc);
                        double TNlb = Nc/Lchain.sumCols(c);
                        if(it ==1){
                            lambda.set(0,c,TNlb);
                        }else {
                            lambda.set(0,c,lambda.get(0,c)*it/ options.iter_max+(Nc/QNc)*lambda.get(0,c)*(options.iter_max-it)/ options.iter_max);

                        }
                    }
                }
            }

            for(int c=0;c<C;c++){
                Matrix inchain = sn.inchain.get(c);
                Matrix lambda_c = new Matrix(1,1,1);
                lambda_c.set(0,0,lambda.get(c));
                chainSysArrivals.put(c,mmap_exponential(lambda_c));
                for(int m=0;m<M;m++){
                    for(int i=0; i<inchain.length();i++){
                        TN.set(m, (int)inchain.get(i),V.get(m,(int)inchain.get(i))*lambda.get(c));
                    }
                }
            }

            for(int ind=0;ind<I;ind++){
                if(sn.isstation.get(ind)==1){
                    int ist = (int) sn.nodeToStation.get(ind);
                    if(sn.nodetypes.get(ind)== NodeType.Join){
                        for(int c=0;c<C;c++){
                            Matrix inchain = sn.inchain.get(c);
                            for(int i=0;i<inchain.length();i++){
                                int fanin =0;
                                for(int j=0;j<sn.rtnodes.numRows;j++){
                                    if(sn.rtnodes.get(j,(int)((ind-1)*K+inchain.get(i)))!=0){
                                        fanin++;
                                    }
                                }
                                int k =(int)inchain.get(i);
                                TN.set(ist,k,lambda.get(c)*V.get(ist,k/fanin));
                                UN.set(ist,k,0);
                                QN.set(ist,k,0);
                                RN.set(ist,k,0);
                            }
                        }
                    }else {
                        if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    TN.set(ist, k, lambda.get(c) * V.get(ist, k));
                                    UN.set(ist, k, S.get(ist, k) * TN.get(ist, k));
                                    QN.set(ist, k, TN.get(ist, k) * S.get(ist, k) * V.get(ist, k));
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                                }
                            }

                        } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.PS) {
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    TN.set(ist, k, lambda.get(c) * V.get(ist, k));
                                    UN.set(ist, k, S.get(ist, k) * TN.get(ist, k));
                                }
                                double Uden = Math.min(1 - GlobalConstants.FineTol, UN.sumRows(ist));
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    QN.set(ist, k, UN.get(ist, k) / (1 - Uden));
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                                }
                            }
                        } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
                            Map<Integer, Map<Integer, Matrix>> chainArrivalAtNode = new HashMap<>();
                            Map<Integer, Map<Integer, Matrix>> rates = new HashMap<>();
                            Map<Integer, Matrix> aggrArrivalAtNode = new HashMap<>();
                            for (int c = 0; c < C; c++) {
                                Matrix a = Matrix.extractRows(V, ist, ist + 1, null);
                                a.scale(lambda.get(c));
                                if (c == 0) {
                                    rates.put(ist, new HashMap<>());
                                }
                                rates.get(ist).put(c, a);
                                Matrix inchain = sn.inchain.get(c);
                                Matrix markProb = new Matrix(1, inchain.length(), inchain.length());
                                for (int i = 0; i < inchain.length(); i++) {
                                    markProb.set(0, i, rates.get(ist).get(c).get((int) inchain.get(i)));
                                }
                                markProb.scale(1 / markProb.elementSum());
                                markProb.removeNaN();
                                chainArrivalAtNode.put(c, mmap_mark(chainSysArrivals.get(c), markProb));
                                chainArrivalAtNode.put(c, mmap_normalize(chainArrivalAtNode.get(c)));
                                Matrix b = new Matrix(rates.get(ist).get(c).numRows, rates.get(ist).get(c).numCols, rates.get(ist).get(c).numRows * rates.get(ist).get(c).numCols);
                                for(int i=0;i<b.numRows;i++){
                                    for (int j=0;j<b.numCols;j++){
                                        b.set(i,j,1/rates.get(ist).get(c).get(i,j));
                                    }
                                }
                                chainArrivalAtNode.put(c, mmap_scale(chainArrivalAtNode.get(c), b));
                                if (c == 0) {
                                    Map<Integer, Map<Integer, Matrix>> MMAPS = new HashMap<>();
                                    MMAPS.put(0, chainArrivalAtNode.get(c));
                                    MMAPS.put(1, mmap_exponential(new Matrix(0),1));
                                    Map<Integer, Matrix> aggrArrivalAtNode_temp = mmap_super_safe(MMAPS, conifg.space_max, "default");
                                    aggrArrivalAtNode.put(0, aggrArrivalAtNode_temp.get(0));
                                    aggrArrivalAtNode.put(1, aggrArrivalAtNode_temp.get(1));
                                    aggrArrivalAtNode.put(2, aggrArrivalAtNode_temp.get(1));
                                    double lc = map_lambda(chainArrivalAtNode.get(c).get(0), chainArrivalAtNode.get(c).get(1));
                                    if (lc > 0) {
                                        aggrArrivalAtNode = mmap_scale(aggrArrivalAtNode, new Matrix(1 / lc));
                                    }
                                } else {
                                    Map<Integer, Map<Integer, Matrix>> MMAPS = new HashMap<>();
                                    MMAPS.put(0, aggrArrivalAtNode);
                                    MMAPS.put(1, chainArrivalAtNode.get(c));
                                    aggrArrivalAtNode = mmap_super_safe(MMAPS, conifg.space_max, "default");
                                }
                            }
                            Map<Integer, Matrix> Qret = new HashMap<>();
                            boolean identical_prio = true;
                            for (int i = 0; i < sn.classprio.length(); i++) {
                                if (sn.classprio.get(i) != sn.classprio.get(0)) {
                                    identical_prio = false;
                                    break;
                                }
                            }
                            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL && !identical_prio) {
                                if (!sn.classprio.hasDuplicates()) {
                                    Qret = MMAPPH1NPPR(mmap_shorten(aggrArrivalAtNode), pie.get(ist), DO.get(ist), 1, null, null, null, null, null, null).get("ncMoms");
                                } else {
                                    throw new RuntimeException("Solver MAM requires either identical priorities or all distinct priorities");
                                }
                            } else {
                                Matrix sn_rates_ist_k = new Matrix(1, K, K);
                                for (int i = 0; i < K; i++) {
                                    sn_rates_ist_k.set(i, sn.rates.get(ist, i));
                                }
                                double aggrUtil = mmap_lambda(aggrArrivalAtNode).element_divide(sn_rates_ist_k.elementIncrease(GlobalConstants.FineTol)).elementSum();
                                if (aggrUtil < 1 - GlobalConstants.FineTol) {
                                    boolean closed = true;
                                    for (int i = 0; i < N.length(); i++) {
                                        if (Double.isInfinite(N.get(i))) {
                                            closed = false;
                                            break;
                                        }
                                    }
                                    if (!closed) {
                                        Qret = MMAPPH1FCFS(mmap_shorten(aggrArrivalAtNode), pie.get(ist), DO.get(ist), 1, null, null, null, false, false, null, null).get("ncMoms");
                                    } else {
                                        Matrix finite_N = N.clone();
                                        finite_N.removeINF();
                                        double maxLevel = finite_N.elementMax()+1;
                                        Map<Integer, Matrix> D = mmap_shorten(aggrArrivalAtNode);
                                        Map<Integer, Matrix> pdistr = new HashMap<>();
                                        if (map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                                            for (int k = 0; k < K; k++) {
                                                Matrix pdistrK = new Matrix(1, 2, 2);
                                                pdistrK.set(0, 1 - GlobalConstants.FineTol);
                                                pdistrK.set(1, GlobalConstants.FineTol);
                                                pdistr.put(k, pdistrK);
                                                Qret.put(k, new Matrix(GlobalConstants.FineTol / sn.rates.get(ist)));
                                            }
                                        } else {
                                            pdistr = MMAPPH1FCFS(D,pie.get(ist),DO.get(ist),null,(int)maxLevel,null,null,false,false,null,null).get("ncDistr");
                                            for (int k=0;k<K;k++){
                                                pdistr.put(k,Matrix.extractRows(pdistr.get(k).transpose(),0,(int) N.get(k)+1,null));
                                                pdistr.get(k).abs();
                                                double sum = 0;
                                                for(int i=0;i<pdistr.get(k).length()-1;i++){
                                                    sum = sum+pdistr.get(k).get(i);
                                                }
                                                pdistr.get(k).set(pdistr.get(k).length()-1,Math.abs(1-sum));
                                                pdistr.get(k).scale(1/pdistr.get(k).elementSum());
                                                Matrix  a = new Matrix(1,(int) N.get(k)+1,(int) N.get(k)+1);
                                                for(int i=0;i<a.length();i++){
                                                    a.set(i,i);
                                                }
                                                Matrix b = new Matrix(1,(int) N.get(k)+1,(int) N.get(k)+1);
                                                for(int i=0;i<a.length();i++){
                                                    b.set(i,pdistr.get(k).get(i));
                                                }
                                                Qret.put(k,new Matrix(Math.max(0,Math.min(N.get(k),a.mult(b.transpose()).get(0)))));
                                            }
                                        }
                                    }
                                }else {
                                    for (int k=0;k<K;k++){
                                        Qret.put(k,new Matrix(sn.njobs.get(k)));
                                    }
                                }
                            }
                            for (int i=0;i<Qret.size();i++){
                                QN.set(ist,i,Qret.get(i).get(0));
                            }
                            for (int k = 0; k < K; k++) {
                                int c = 0;
                                for(int i=0;i<sn.chains.numRows;i++){
                                    if(sn.chains.get(i,k)!=0){
                                        c = i;
                                        break;
                                    }
                                }
                                TN.set(ist,k,rates.get(ist).get(c).get(k));
                                UN.set(ist,k,TN.get(ist,k)*S.get(ist,k));
                                QN.set(ist,k,Qret.get(k).get(0));
                                QN.set(ist,k,QN.get(ist,k)+TN.get(ist,k)*(S.get(ist,k)*sn.nservers.get(ist))*(sn.nservers.get(ist)-1)/sn.nservers.get(ist));
                                RN.set(ist,k,QN.get(ist,k)/TN.get(ist,k));
                            }
                        }
                    }
                }else {
                    if(sn.nodetypes.get(ind)==NodeType.Fork){
                        // TODO: not implemented
                        throw new RuntimeException("Fork nodes not supported yet by MAM solver.");
                    }
                }
            }
            dif_matrix  = TN.add(-1,TN_1);
            dif = Math.abs(dif_matrix.elementMaxAbs());
        }

        CN = RN.sumCols();
        QN.abs();
        for (int i =1;i<=2;i ++){
            for(int c=0;c<C;c++){
                Matrix inchain = sn.inchain.get(c);
                double Nc = 0;
                for(int j=0;j<inchain.length();j++){
                    Nc = Nc + sn.njobs.get((int) inchain.get(j));
                }
                if(Double.isFinite(Nc)){
                    double QNc = 0 ;
                    for(int j=0;j<inchain.length();j++){
                        QNc = QNc + QN.sumCols((int) inchain.get(j));
                    }
                    for(int m=0;m<inchain.length();m++){
                        for(int n=0;n<QN.numRows;n++){
                            QN.set(n,(int)inchain.get(m),QN.get(n,(int)inchain.get(m))*(Nc/QNc));
                        }
                    }
                }

                for(int ind=0;ind<I;ind++){
                    for(int k=0;k<inchain.length();k++){
                        if(sn.isstation.get(ind)==1){
                            int ist = (int)sn.nodeToStation.get(ind);
                            if(V.get(ist,(int)inchain.get(k))>0){
                                if(Double.isInfinite(sn.nservers.get(ist))){
                                    RN.set(ist,(int)inchain.get(k),S.get(ist,(int)inchain.get(k)));
                                }else {
                                    RN.set(ist,(int) inchain.get(k),Math.max(S.get(ist,(int) inchain.get(k)),QN.get(ist,(int)inchain.get(k))/TN.get(ist,(int) inchain.get(k))));
                                }
                            }else {
                                RN.set(ist,(int)inchain.get(k),0);
                            }
                            QN.set(ist,(int) inchain.get(k), RN.get(ist,(int) inchain.get(k))*TN.get(ist,(int) inchain.get(k)));
                        }
                    }
                }
                if(Nc==0){
                    for(int k =0;k<QN.numRows;k++){
                        QN.set(k,c,0);
                    }
                    for(int k =0;k<UN.numRows;k++){
                        UN.set(k,c,0);
                    }
                    for(int k =0;k<RN.numRows;k++){
                        RN.set(k,c,0);
                    }
                    for(int k =0;k<TN.numRows;k++){
                        TN.set(k,c,0);
                    }
                    CN.set(0,c,0);
                    XN.set(0,c,0);
                }
            }
        }
        SolverMAMResult result = new SolverMAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.WN = WN;
        result.AN = AN;
        return result;
    }

    public SolverMAMResult solver_mam(NetworkStruct sn){
        SolverOptions.Config config = options.config;
        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix V = new Matrix(sn.visits.get(0));
        if (sn.nchains>1) {
            for (int i = 1; i < sn.nchains; i++) {
                V.add(1,sn.visits.get(i));
            }
        }

        Matrix QN = new Matrix(M,K,M*K);
        Matrix UN = new Matrix(M,K,M*K);
        Matrix RN = new Matrix(M,K,M*K);
        Matrix TN = new Matrix(M,K,M*K);
        Matrix AN = new Matrix(M,K,M*K);
        Matrix WN = new Matrix(M,K,M*K);
        Matrix CN = new Matrix(1,K,K);
        Matrix XN = new Matrix(1,K,K);

        Matrix lambda = new Matrix(1,K,K);
        for(int c=0;c<C;c++){
            Matrix inchain = sn.inchain.get(c);
            Matrix lambdas_inchain = new Matrix(1, inchain.length(), inchain.length());
            for(int i=0;i< inchain.length();i++){
                lambdas_inchain.set(0,i,sn.rates.get((int) sn.refstat.get((int) inchain.get(0)),(int) inchain.get(i)));
                double sum = 0;
                for(int j=0; j<lambdas_inchain.length();j++){
                    if(!Double.isInfinite(lambdas_inchain.get(j))){
                        sum = sum + lambdas_inchain.get(j);
                    }
               }
            }
        }

        Matrix chain = new Matrix(1,K,K);

        for(int k=0;k<K;k++){
            for(int i=0;i<chain.numCols;i++){
                if(sn.chains.get(i,k)!=0){
                    chain.set(0,k,i);
                    break;
                }
            }
        }

        boolean isopen = true;
        for(int i=0;i<sn.njobs.length();i++){
            if(Double.isFinite(sn.njobs.get(i))){
                isopen = false;
                break;
            }
        }

        if(isopen) {
            int last_it = 0;
            Map<Integer, Map<Integer, Matrix>> pie = new HashMap<>();
            Map<Integer, Map<Integer, Matrix>> D0 = new HashMap<>();
            for (int ist = 0; ist < M; ist++) {
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                    for (int i = 0; i < TN.numCols; i++) {
                        if (!Double.isNaN(sn.rates.get(ist, i))) {
                            TN.set(ist, i, sn.rates.get(ist, i));
                        }
                    }

                } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.PS) {
                    for (int k = 0; k < K; k++) {
                        Matrix D0_ = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0).clone();
                        D0_.scale(1 / sn.nservers.get(ist));
                        Matrix D1 = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1).clone();
                        D1.scale(1 / sn.nservers.get(ist));
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k), map_scale(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0), PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1), map_mean(D0_, D1)));
                        if (k == 0) {
                            pie.put(ist, new HashMap<>());
                            D0.put(ist, new HashMap<>());
                        }
                        pie.get(ist).put(k, map_pie(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0), PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));

                        D0.get(ist).put(k, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0));
                    }
                }
            }
            int it_max = options.iter_max;
            Map<Integer, Map<Integer, Map<Integer, Matrix>>> DEP = ph_reindex(PH,sn);
            for (int it = 0; it < it_max; it++) {
                last_it = it;
                if (it == 1) {
                    for (int ind = 1; ind < M; ind++) {
                        for (int r = 0; r < K; r++) {
                            int ist = (int) sn.nodeToStation.get(ind);
                            DEP.get(ind).put(r, map_scale(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0), PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(1), 1 / (lambda.get(r) * V.get(ind, r))));
                        }
                    }
                }

                Map<Integer, Map<Integer, Matrix>> ARV = solver_mam_traffic(sn, DEP, config);
                Matrix QN_1 = QN.clone();
                for (int ist = 0; ist < M; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    if (sn.nodetypes.get(ind) == NodeType.Queue) {
                        if (ARV.get(ind).get(0).length() > 50) { // config.space_max
                            System.out.println("Arrival process at node " + ind + " is now at " + ARV.get(ind).get(0).length() + " states. Compressing");
                            ARV.put(ind, mmap_compress(ARV.get(ind)));
                        }
                        Map<Integer, Matrix> D = new HashMap<>();
                        for (int i = 0; i < ARV.get(ind).size(); i++) {
                            if (i == 0) {
                                D.put(0, ARV.get(ind).get(0));
                            } else if (i > 1) {
                                D.put(i, ARV.get(ind).get(i - 1));
                            }
                        }
                        Map<Integer, Matrix> Qret = MMAPPH1FCFS(D, pie.get(ist), D0.get(ist), 1, 2, null, null, false, false, null, null).get("ncMoms");
                        for (int k = 0; k < K; k++) {
                            QN.set(ist, k, Qret.get(k).elementSum());
                        }
                        TN.insert_sub_matrix(ist, 0, ist + 1, TN.numCols, mmap_lambda(ARV.get(ind)));
                    }
                    for (int k=0;k<K;k++){
                        UN.set(ist,k,TN.get(ist,k)*map_mean(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));
                        QN.set(ist,k,QN.get(ist,k)+TN.get(ist,k)*map_mean(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));
                        RN.set(ist,k,QN.get(ist,k)/TN.get(ist,k));
                    }
                }

                if (it>=3){
                    Matrix Q_diff = QN.add(-1,QN_1);
                    Q_diff.abs();
                    if(Q_diff.element_divide(QN_1).elementMax()<options.iter_tol){
                        break;
                    }
                }

                for(int ist=0;ist<M;ist++){
                    int ind = (int) sn.stationToNode.get(ist);
                    if(sn.nodetypes.get(ind)==NodeType.Queue){
                        for(int r=0;r<K;r++){
                            Matrix types = new Matrix(1,K-1,K-1);
                            for(int i=0;i<K;i++){
                                if (i<r){
                                    types.set(i,i);
                                }else if (i>r){
                                    types.set(i,i+1);
                                }
                            }
                            Map<Integer,Matrix> A = mmap_hide(ARV.get(ind), types);
                            Matrix tA = A.get(1).sumRows();
                            Matrix pieA = map_pie(A.get(0),A.get(1));

                            Map<Integer,Matrix> S = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
                            Matrix pieS = map_pie(S.get(0),S.get(1));
                            Matrix tS = S.get(1).sumRows();

                            double rho = UN.sumRows(ist);

                            A = mmap_scale(A,new Matrix(map_mean(A.get(0),A.get(1))-map_mean(S.get(0),S.get(1))));
                            Matrix zAS = Matrix.createLike(tA.mult(pieS));
                            Matrix zSA = Matrix.createLike(tS.mult(pieA));
                            Matrix zA = Matrix.createLike(A.get(1));

                            Matrix DEP0ir = Matrix.concatRows(Matrix.concatColumns(S.get(0),Matrix.scale_mult(tS.mult(pieA),rho),null),Matrix.concatColumns(zAS,A.get(0),null),null);
                            Matrix DEP1ir = Matrix.concatRows(Matrix.concatColumns(Matrix.scale_mult(S.get(1),1-rho),zSA,null),Matrix.concatColumns(tA.mult(pieS),zA,null),null);
                            DEP.get(ind).put(r,map_normalize(DEP0ir,DEP1ir));
                            DEP.get(ind).put(r,map_scale(DEP.get(ind).get(r).get(0),DEP.get(ind).get(r).get(1),1/(lambda.get(r)*V.get(ind,r))));
                        }
                    }
                }
            }

            if(options.verbose== VerboseLevel.STD||options.verbose== VerboseLevel.DEBUG){
                System.out.println("MAM parametric decomposition completed in "+last_it+" iterations");
            }
        }else {
            System.out.println("This model is not supported by SolverMAM yet. Returning with no result.");
        }
        SolverMAMResult result = new SolverMAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.WN = WN;
        result.AN = AN;
        result.CN = CN;
        result.XN = XN;
        return result;
    }
    public SolverMAMResult solver_qna_mam_closed(NetworkStruct sn){
        SolverOptions.Config config = options.config;
        config.space_max = 1;

        int K = sn.nclasses;
        Matrix rt = sn.rt.clone();
        Matrix S = sn.rates.element_power(-1);
        Matrix scv = sn.scv.clone();
        scv.removeNaN();
        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix N = sn.njobs.clone();
        Matrix V = Matrix.cellsum(sn.visits);

        Map<Integer,Map<Integer,Matrix>> pie = new HashMap<>();
        Map<Integer,Map<Integer,Matrix>> DO = new HashMap<>();
        Matrix Q = new Matrix(M,K,M*K);
        Matrix U = new Matrix(M,K,M*K);
        Matrix R = new Matrix(M,K,M*K);
        Matrix T = new Matrix(M,K,M*K);
        Matrix X = new Matrix(1,K,K);


        for(int ist=0; ist<M;ist++){
            if(sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS){
                pie.put(ist,new HashMap<Integer,Matrix>());
                DO.put(ist,new HashMap<Integer,Matrix>());
                for(int k = 0;k<K;k++){
                    PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k),map_scale(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1),S.get(ist,k)/sn.nservers.get(ist)));
                    pie.get(ist).put(k,map_pie(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0),PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(1)));
                    DO.get(ist).put(k,PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0));
                }
            }
        }

        Matrix lambda = new Matrix(1,C,C);
        Matrix QNc = sn.njobs.clone();
        int it_out = 0;
        Matrix lambda_lb = new Matrix(1,K,K);
        Matrix lambda_ub = new Matrix(1,K,K);
        for (int k=0;k<K;k++){
            for(int i=0;i<I;i++){
                if(sn.nservers.get(i)!=Double.POSITIVE_INFINITY){
                    if(lambda_ub.get(k)==0){
                        lambda_ub.set(k,sn.rates.get(i,k));
                    }else {
                        lambda_ub.set(k,Math.min(lambda_ub.get(k),sn.rates.get(i,k)));
                    }
                }
            }
        }



        Map<Integer,Matrix> lambdas_inchain = new HashMap<>();
        Map<Integer,Matrix> scvs_inchain = new HashMap<>();
        Matrix d2c = new Matrix(1,C,C);



        Matrix QN = new Matrix(1,K,K);
        Matrix QN_diff = QN.add(-1,QNc);
        QN_diff.abs();

        while(QN_diff.elementMax()>options.iter_tol&&it_out<200) {
            it_out++;
            if(it_out==1){
                lambda = lambda_ub.clone();
            }else {
                for(int k=0;k<K;k++){
                    if(QN.get(k)<QNc.get(k)){
                        lambda_lb.set(k,lambda.get(k));
                    }else {
                        lambda_ub.set(k,lambda.get(k));
                    }
                    lambda.set(k,(lambda_ub.get(k)+lambda_lb.get(k))/2);
                }

            }
            int it=0;
            Q = new Matrix(M,K,M*K);
            U = new Matrix(M,K,M*K);
            R = new Matrix(M,K,M*K);
            T = new Matrix(M,K,M*K);
            X = new Matrix(1,K,K);
            Matrix a1 = new Matrix(M,K,M*K);
            Matrix a2 = new Matrix(M,K,M*K);
            Matrix a1_1 = a1.elementIncrease(Double.POSITIVE_INFINITY);

            Matrix a2_1 = a2.elementIncrease(Double.POSITIVE_INFINITY);
            Matrix a1_diff = a1.add(-1,a1_1);
            a1_diff.abs();
            Matrix a2_diff = a2.add(-1,a2_1);
            a2_diff.abs();
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
            while ((a1_diff.elementMax() > options.iter_tol || a2_diff.elementMax() > options.iter_tol) && it <= options.iter_max) {
                it = it + 1;
                a1_1 = a1.clone();
                a2_1 = a2.clone();

                if (it == 1) {
                    for (int c = 0; c < C; c++) {
                        Matrix inchain = sn.inchain.get(c);
                        for (int m = 0; m < M; m++) {
                            for (int i = 0; i < inchain.length(); i++) {
                                T.set(m, (int) inchain.get(i), V.get(m, (int) inchain.get(i)) * lambda.get(c));
                            }
                        }
                    }
                }
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < K; j++) {
                        a1.set(i, j, 0);
                        a2.set(i, j, 0);
                    }
                    double lambda_i = T.sumRows(i);
                    for (int j = 0; j < M; j++) {
                        for (int r = 0; r < K; r++) {
                            for (int s = 0; s < K; s++) {
                                a1.set(i, r, a1.get(i, r) + T.get(j, s) * rt.get(j * K + s, i * K + r));
                                a2.set(i, r, a2.get(i, r) + (1 / lambda_i) * f2.get(j * K + s, i * K + r) * T.get(j, s) * rt.get(j * K + s, i * K + r));
                            }
                        }
                    }
                }

                for (int ind = 0; ind < I; ind++) {
                    if (sn.isstation.get(ind) == 1) {
                        int ist = (int) sn.nodeToStation.get(ind);
                        if (sn.nodetypes.get(ind) != NodeType.Join) {
                            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                                for (int i = 0; i < M; i++) {
                                    for (int r = 0; r < K; r++) {
                                        for (int s = 0; s < K; s++) {
                                            d2.set(ist, s, a2.get(ist, s));
                                        }
                                    }
                                }
                                for (int c = 0; c < C; c++) {
                                    Matrix inchain = sn.inchain.get(c);
                                    for (int k1 = 0; k1 < inchain.length(); k1++) {
                                        int k = (int) inchain.get(k1);
                                        T.set(ist, k, a1.get(ist, k));
                                        U.set(ist, k, S.get(ist, k) * T.get(ist, k));
                                        Q.set(ist, k, T.get(ist, k) * S.get(ist, k) * V.get(ist, k));
                                        R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                                    }
                                }
                            } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
                                Matrix mu_ist = new Matrix(1, K, K);
                                for (int i = 0; i < K; i++) {
                                    mu_ist.set(i, sn.rates.get(ist, i));
                                }
                                mu_ist.removeNaN();
                                Matrix rho_ist_class = new Matrix(1, K, K);
                                for (int i = 0; i < K; i++) {
                                    rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)));
                                }
                                rho_ist_class.removeNaN();
                                double lambda_ist = a1.sumRows(ist);
                                int mi = (int) sn.nservers.get(ist);
                                double rho_ist = rho_ist_class.elementSum() / mi;
                                double c2 = 0;
                                if (rho_ist < 1 - options.tol) {
                                    for (int k = 0; k < K; k++) {


                                        double mubar = lambda_ist / rho_ist;
                                        c2 = -1;
                                        for (int r = 0; r < K; r++) {
                                            if (mu_ist.get(r) > 0) {
                                                c2 = c2 + a1.get(ist, r) / lambda_ist * Math.pow(mubar / mi / mu_ist.get(r), 2) * (scv.get(ist, r) + 1);
                                            }
                                        }


                                    }
                                    d2.set(ist, 1 + Math.pow(rho_ist, 2) * (c2 - 1) / Math.sqrt(mi) + (1 - Math.pow(rho_ist, 2)) * (a2.sumRows(ist) - 1));
                                } else {
                                    for (int k = 0; k < K; k++) {
                                        Q.set(ist, k, sn.njobs.get(k));
                                    }
                                    d2.set(ist, 1);
                                }
                                for (int k = 0; k < K; k++) {
                                    T.set(ist, k, a1.get(ist, k));
                                    U.set(ist, k, T.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist));
                                }
                            }
                        }
                    } else {
                        if (sn.nodetypes.get(ind) == NodeType.Fork) {
                            // TODO: not implemented
                            throw new RuntimeException("Fork nodes not supported yet by QNA solver.");
                        }
                    }
                }

                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        if (sn.nodetypes.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                            for (int r = 0; r < K; r++) {
                                for (int s = 0; s < K; s++) {
                                    if (rt.get(i * K + r, j * K + s) > 0) {
                                        f2.set(i * K + r, j * K + s, 1 + rt.get(i * K + r, j * K + s) * (d2.get(i) - 1));
                                    }
                                }
                            }
                        }
                    }
                }
                a1_diff = a1.add(-1, a1_1);
                a1_diff.abs();
                a2_diff = a2.add(-1, a2_1);
                a2_diff.abs();
            }

            for (int ind = 0; ind < I; ind++) {
                if (sn.isstation.get(ind) == 1) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {

                        Matrix rho_ist_class = new Matrix(1, K, K);
                        for (int i = 0; i < K; i++) {
                            rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)));
                        }
                        rho_ist_class.removeNaN();

                        int mi = (int) sn.nservers.get(ist);
                        double rho_ist = rho_ist_class.elementSum() / mi;
                        double c2 = 0;
                        if (rho_ist < 1 - options.tol) {
                            Map<Integer, Matrix> arri_class_anx = new HashMap<>();
                            Map<Integer, Matrix> arri_class = new HashMap<>();
                            Map<Integer, Matrix> arri_node = new HashMap<>();
                            for (int k = 0; k < K; k++) {
                                if (a1.get(ist, k) == 0) {
                                    arri_class_anx = map_exponential(Double.POSITIVE_INFINITY);
                                } else {
                                    arri_class_anx = APH.fitMeanAndSCV(1 / a1.get(ist, k), a2.get(ist, k)).getRepres();
                                }
                                arri_class.put(0, arri_class_anx.get(0));
                                arri_class.put(1, arri_class_anx.get(1));
                                arri_class.put(2, arri_class_anx.get(1));
                                if (k == 0) {
                                    arri_node.put(0, arri_class.get(0));
                                    arri_node.put(1, arri_class.get(1));
                                    arri_node.put(2, arri_class.get(2));
                                } else {
                                    arri_node = mmap_super(arri_node, arri_class);
                                }


                            }
                            Matrix finite_N = N.clone();
                            finite_N.removeINF();
                            double maxLevel = finite_N.elementMax()+1;
                            Map<Integer, Matrix> D = mmap_shorten(arri_node);
                            Map<Integer, Matrix> pdistr = new HashMap<>();
                            Map<Integer,Matrix> Qret = new HashMap<>();
                            if (map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                                for (int k = 0; k < K; k++) {
                                    Matrix pdistrK = new Matrix(1, 2, 2);
                                    pdistrK.set(0, 1 - GlobalConstants.FineTol);
                                    pdistrK.set(1, GlobalConstants.FineTol);
                                    pdistr.put(k, pdistrK);
                                    Qret.put(k, new Matrix(GlobalConstants.FineTol / sn.rates.get(ist)));
                                }
                            } else {
                                pdistr = MMAPPH1FCFS(D,pie.get(ist),DO.get(ist),null,(int)maxLevel,null,null,false,false,null,null).get("ncDistr");
                                for (int k=0;k<K;k++){
                                    pdistr.put(k,Matrix.extractRows(pdistr.get(k).transpose(),0,(int) N.get(k)+1,null));
                                    pdistr.get(k).abs();
                                    double sum = 0;
                                    for(int i=0;i<pdistr.get(k).length()-1;i++){
                                        sum = sum+pdistr.get(k).get(i);
                                    }
                                    pdistr.get(k).set(pdistr.get(k).length()-1,Math.abs(1-sum));
                                    pdistr.get(k).scale(1/pdistr.get(k).elementSum());
                                    Matrix  a = new Matrix(1,(int) N.get(k)+1,(int) N.get(k)+1);
                                    for(int i=0;i<a.length();i++){
                                        a.set(i,i);
                                    }
                                    Matrix b = new Matrix(1,(int) N.get(k)+1,(int) N.get(k)+1);
                                    for(int i=0;i<a.length();i++){
                                        b.set(i,pdistr.get(k).get(i));
                                    }
                                    Qret.put(k,new Matrix(Math.max(0,Math.min(N.get(k),a.mult(b.transpose()).get(0)))));
                                }
                            }
                            for (int i=0;i<Qret.size();i++){
                                Q.set(ist,i,Qret.get(i).get(0));
                            }
                        } else {
                            for (int k = 0; k < K; k++) {
                                Q.set(ist, k, sn.njobs.get(k));
                            }
                        }
                        for (int k = 0; k < K; k++) {
                            R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                        }

                    }
                }
            }
            QN = Q.sumCols();
            QN_diff = QN.add(-1,QNc);
            QN_diff.abs();
        }

        SolverMAMResult result = new SolverMAMResult();
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
        result.iter = it_out;

        return result;

    }

    private Map<Integer,Map<Integer,Matrix>> solver_mam_traffic(NetworkStruct sn, Map<Integer, Map<Integer, Map<Integer, Matrix>>> DEP, SolverOptions.Config config){
        int I = sn.nnodes;
        int C = sn.nchains;
        int R = sn.nclasses;
        Matrix non_cs_classes = new Matrix(1,I*R,I*R);
        Matrix isNCS = new Matrix(1,I,I);
        Matrix nodeToNCS = new Matrix(1,I,I);
        int end =0;
        for(int ind =0;ind<I;ind++){
            if(sn.nodetypes.get(ind)!=NodeType.ClassSwitch){
                for(int i=0;i<R;i++){
                    non_cs_classes.set(i+end,ind*R+i);
                }
                end = end+R;
                isNCS.set(ind,1);
                nodeToNCS.set(ind,isNCS.elementSum());
            }else {
                isNCS.set(ind,0);
            }
        }

        List<Integer> non_cs_classes_list = new ArrayList<>();
        for(int i=0;i<non_cs_classes.length();i++){
            non_cs_classes_list.add((int)non_cs_classes.get(i));
        }

        Matrix rtncs = DTMC.dtmc_stochcomp(sn.rtnodes,non_cs_classes_list);
        int Inc = I;
        for(int i=0;i<sn.nodetypes.size();i++){
            if(sn.nodetypes.get(i)==NodeType.ClassSwitch){
                Inc = Inc -1;
            }
        }

        for(int ist=0;ist<DEP.size();ist++){
            for(int r=0; r<DEP.get(ist).size();r++){
                if(DEP.get(ist).get(r).isEmpty()||DEP.get(ist).get(r).get(0).hasNaN()){
                    DEP.get(ist).get(r).put(0,new Matrix(1,1,0));
                    DEP.get(ist).get(r).put(1,new Matrix(1,1,0));
                    DEP.get(ist).get(r).put(2,new Matrix(1,1,0));
                }else {
                    DEP.get(ist).get(r).put(2,DEP.get(ist).get(ist).get(1));
                }
            }
        }


        Map<Integer,Map<Integer,Map<Integer,Matrix>>> DEP_new = new HashMap<>();
        Map<Integer,Map<Integer,Map<Integer,Matrix>>> LINKS = new HashMap<>();
        Map<Integer,Map<Integer,Matrix>> ARV = new HashMap<>();
        for(int ind=0;ind<I;ind++){
            if(isNCS.get(ind)==1){
                int inc = (int) nodeToNCS.get(ind);
                if(sn.nodetypes.get(ind)==NodeType.Source||sn.nodetypes.get(ind)==NodeType.Delay||sn.nodetypes.get(ind)==NodeType.Queue){
                    int ist = (int) sn.nodeToStation.get(ind);
                    if(R>1){
                        for(int r=0;r<R;r++) {
                            DEP_new.put(inc,new HashMap<>());
                            DEP_new.get(inc).put(r, mmap_super(DEP.get(ist).get(r)));
                        }
                    }else {
                        DEP_new.put(inc,new HashMap<>());
                        DEP_new.get(inc).put(0,DEP.get(ist).get(0));
                    }
                    Matrix Psplit = new Matrix(R,Inc*R,R*Inc*R);
                    for(int r=0;r<R;r++){
                        for(int jnd = 0; jnd<I;jnd++){
                            if(isNCS.get(jnd)==1){
                                int jnc = (int) nodeToNCS.get(jnd);
                                for (int s=0;s<R;s++){
                                    Psplit.set(r,(jnc-1)*R+s,rtncs.get((inc-1)*R+r,(jnc-1)*R+s));
                                }
                            }
                        }
                    }
                    Map<Integer,Map<Integer,Matrix>> Fsplit = npfqn_traffic_split_cs(DEP_new.get(inc).get(0),Psplit);
                    for(int jnc = 0; jnc<Inc;jnc++){
                        LINKS.put(inc,new HashMap<>());
                        LINKS.get(inc).put(jnc, Fsplit.get(jnc));
                        LINKS.get(inc).put(jnc,mmap_normalize(LINKS.get(inc).get(jnc)));
                    }
                }
            }
        }

        for(int ind =0; ind<I;ind++){
            Map<Integer,Map<Integer,Matrix>> FLOWS = new HashMap<>();
            if(isNCS.get(ind)==1 && sn.nodetypes.get(ind) != NodeType.Source){
                int inc = (int) nodeToNCS.get(ind);
                for(int jnd =0;jnd<Inc;jnd++){
                    if(!LINKS.get(jnd).get(inc).isEmpty()&&mmap_lambda(LINKS.get(jnd).get(inc)).elementSum()>GlobalConstants.FineTol){
                        FLOWS.put(FLOWS.size(),LINKS.get(jnd).get(inc));
                    }
                }
                if(FLOWS.size()>1){
                    ARV.put(ind,npfqn_traffic_merge(FLOWS,config.merge,config.compress));
                }else if(FLOWS.size()==1){
                    ARV.put(ind,FLOWS.get(0));
                }else {
                    ARV.put(ind,LINKS.get(Inc-1).get(0));
                }
            }else {
                ARV.put(ind,new HashMap<>());
            }
        }

        return ARV;
    }

    public Map<Integer,Map<Integer,Matrix>> solver_mam_passage_time(NetworkStruct sn, Map<Integer,Map<Integer,Map<Integer,Matrix>>> PH){
        Map<Integer,Map<Integer,Matrix>> RD = new HashMap<>();
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix N= sn.njobs.transpose();

        boolean open = true;
        for(int i=0;i< N.length();i++){
            if(Double.isFinite(N.get(i))){
                open = false;
                break;
            }
        }
        Map<Integer,Matrix> A = new HashMap<>();
        int idx_arv = 0;
        int idx_q = 0;
        if(M==2 && open){
            Map<Integer,Matrix> pie = new HashMap<>();
            Map<Integer,Matrix> S = new HashMap<>();
            for(int i=0;i<M;i++){
                if(sn.sched.get(sn.stations.get(i))==SchedStrategy.EXT){
                    A.put(0,PH.get(i).get(0).get(0));
                    A.put(1,PH.get(i).get(0).get(1));
                    A.put(2,PH.get(i).get(0).get(1));
                    for (int k=1;k<K;k++){
                        Map<Integer,Matrix> B = new HashMap<>();
                        A.put(0,PH.get(i).get(k).get(0));
                        A.put(1,PH.get(i).get(k).get(1));
                        A.put(2,PH.get(i).get(k).get(1));
                        A = mmap_super(A,B);
                    }
                    idx_arv = i;
                }else if(sn.sched.get(sn.stations.get(i))==SchedStrategy.FCFS||sn.sched.get(sn.stations.get(i))==SchedStrategy.HOL){
                    for (int k=0;k<K;k++) {
                        Matrix PH_i_K_D0 = PH.get(i).get(k).get(0);
                        PH_i_K_D0.scale(1/sn.nservers.get(i));
                        Matrix PH_i_K_D1 = PH.get(i).get(k).get(1);
                        PH_i_K_D1.scale(1/sn.nservers.get(i));
                        Map<Integer, Matrix> PH_new_i_k = map_scale(PH.get(i).get(k).get(0),PH.get(i).get(k).get(1),map_mean(PH_i_K_D0,PH_i_K_D1));
                        pie.put(k,map_pie(PH_new_i_k.get(0),PH_new_i_k.get(1)));
                        S.put(k,PH_new_i_k.get(0));
                    }
                    idx_q = i;
                }else {
                    throw new RuntimeException("Unsupported scheduling strategy");
                }
            }

            boolean identical_priority = true;
            for(int i=0;i<sn.classprio.length();i++){
                if(sn.classprio.get(i)!=sn.classprio.get(0)){
                    identical_priority = false;
                }
            }
            if(!identical_priority){
                throw new RuntimeException("Response time distribution in priority models not yet supported.");
            }else {
                for (int i=1;i<A.size()-1;i++){
                    A.put(i,A.get(i+1));
                }
                A.remove(A.size()-1);
                Map<Integer,Matrix> alpha = MMAPPH1FCFS(A,pie,S,null,null,null,null,false,true,null,null).get("stDistrPH_alpha");
                Map<Integer,Matrix> D0 = MMAPPH1FCFS(A,pie,S,null,null,null,null,false,true,null,null).get("stDistrPH_D0");
                Map<Integer,Map<Integer,Matrix>> RDph = new HashMap<>();
                Map<Integer,Double> sigma = new HashMap<>();
                Map<Integer,Double> mean = new HashMap<>();
                for (int k=0;k<K;k++) {
                    RDph.put(k, new HashMap<>());
                    Matrix neg_D0k = D0.get(k);
                    neg_D0k.scale(-1);
                    RDph.get(k).put(0,D0.get(k));
                    RDph.get(k).put(1,neg_D0k.mult(Matrix.ones(alpha.get(k).length(),1)).mult(alpha.get(k).transpose()));
                    sigma.put(k,Math.sqrt(map_var(RDph.get(k).get(0),RDph.get(k).get(1))));
                    mean.put(k,map_mean(RDph.get(k).get(0),RDph.get(k).get(1)));
                    int n=5;
                    Matrix point = new Matrix(mean.get(k)+n*sigma.get(k));

                    while (map_cdf(RDph.get(k).get(0),RDph.get(k).get(1),point).get(0)<1-GlobalConstants.FineTol){
                        n++;
                    }
                    Matrix X = new Matrix(10000,1,10000);
                    Matrix F = new Matrix(10000,1,10000);
                    double stepSize = mean.get(k)+n*sigma.get(k)/9999;
                    for(int i=0;i<10000;i++){
                        X.set(i,0,stepSize*i);
                        F.set(i,0,map_cdf(RDph.get(k).get(0),RDph.get(k).get(1),new Matrix(stepSize*i)).get(0));
                    }
                    RD.put(idx_arv,new HashMap<>());
                    RD.get(idx_arv).put(k,new Matrix(0,0,0));
                    RD.put(idx_q,new HashMap<>());
                    RD.get(idx_q).put(k,Matrix.concatColumns(F,X,null));
                }

            }

        }else {
            // TODO: not implemented
            throw new RuntimeException("This model is not supported by SolverMAM yet.");
        }

        return RD;
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.MAM);
    }

}



