package jline.lib.thirdparty;

import jline.api.CTMC;
import jline.util.Matrix;

import java.util.*;
import java.util.stream.IntStream;

import static jline.lib.thirdparty.SMCSolver.*;
import static org.apache.commons.math3.util.CombinatoricsUtils.factorial;

public class BUTOOLS {
    public static Map<String, Matrix> QBDFundamentalMatrices(Matrix B, Matrix L, Matrix F, Double precision_, Integer maxNumIt_, String method_, Integer Verbose_){
        double precision = 1e-14;
        if(precision_!=null){
            precision = precision_;
        }

        int maxNumIt = 50;
        if(maxNumIt_!=null){
            maxNumIt = maxNumIt_;
        }

        String method = "CR";
        if(method_!= null){
            method = method_;
        }

        int Verbose = 0;
        if(Verbose_!=null){
            Verbose = Verbose_;
        }

        if(method.equals("LR")){
            return QBD_LR(B,L,F,maxNumIt,Verbose,null,null);
        }

        if(method.equals("NI")){
            return QBD_NI(B,L,F,maxNumIt,Verbose,null,null);
        }

        if(method.equals("IS")) {
            return QBD_IS(B,L,F,maxNumIt,Verbose,null,null);
        }

        if(method.equals("FI")){
            return QBD_FI(B,L,F,maxNumIt,Verbose,null,null,null);
        }

        return QBD_CR(B,L,F,maxNumIt,Verbose,null,null);
    }

    public static Map<String,Matrix> FluidFundamentalMatrices(Matrix Fpp,Matrix Fpm, Matrix Fmp, Matrix Fmm, Double precision_, Integer maxNumIt_, String method_){
        Matrix Psi = new Matrix(0);
        double precision;
        if(precision_==null){
            precision = 1e-14;
        }else {
            precision = precision_;
        }

        int maxNumIt;
        if(maxNumIt_==null){
            maxNumIt = 150;
        }else {
            maxNumIt = maxNumIt_;
        }

        String method;
        if(method_ == null){
            method = "ADDA";
        }else {
            method = method_;
        }
        int numit = 0;
        if(Fpp.numRows==0){
            Psi = new Matrix(0,Fmm.numRows);
        }else if(method.equals("CR")){

        }else if(method.equals("ADDA")||method.equals("SDA")){
            Matrix A = Fpp.clone();
            A.scale(-1);
            Matrix B = Fpm.clone();
            Matrix C = Fmp.clone();
            Matrix D = Fmm.clone();
            D.scale(-1);
            Matrix diag_A = new Matrix(0);
            Matrix.extractDiag(A,diag_A);
            double gamma1 = diag_A.elementMax();
            Matrix diag_D = new Matrix(0);
            Matrix.extractDiag(D,diag_D);
            double gamma2 = diag_D.elementMax();
            if(method.equals("SDA")){
                gamma1 = Math.max(gamma1,gamma2);
                gamma2 = gamma1;
            }
            int sA = A.numRows;
            int sD = D.numRows;
            Matrix IA = Matrix.eye(sA);
            Matrix ID = Matrix.eye(sD);
            Matrix gamma2IA = IA.clone();
            gamma2IA.scale(gamma2);
            Matrix gamma1ID = ID.clone();
            gamma1ID.scale(gamma1);
            A = A.add(1,gamma2IA);
            D = D.add(1,gamma1ID);
            Matrix Dginv = D.inv();
            Matrix Vginv = D.add(-1,C.mult(A.inv()).mult(B)).inv();
            Matrix Wginv = A.add(-1,B.mult(Dginv).mult(C)).inv();
            Matrix gammaVginv = Vginv.clone();
            gammaVginv.scale(gamma1+gamma2);
            Matrix Eg = ID.add(-1,gammaVginv);
            Matrix gammaWginV = Wginv.clone();
            gammaWginV.scale(gamma1+gamma2);
            Matrix Fg = IA.add(-1,gammaWginV);
            Matrix Gg = Dginv.mult(C).mult(Wginv);
            Gg.scale(gamma1+gamma2);
            Matrix Hg = Wginv.mult(B).mult(Dginv);
            Hg.scale(gamma1+gamma2);

            double diff = 1;
            while (diff>precision && numit<maxNumIt){
                Vginv = Eg.mult(ID.add(-1,Gg.mult(Hg)).inv());
                Wginv = Fg.mult(IA.add(-1,Hg.mult(Gg)).inv());
                Gg = Gg.add(1,Vginv.mult(Gg).mult(Fg));
                Hg = Hg.add(1, Wginv.mult(Hg).mult(Eg));
                Eg = Vginv.mult(Eg);
                Fg = Wginv.mult(Fg);
                double neg = Matrix.first_norm(Eg);
                double nfg = Matrix.first_norm(Fg);
                if(method.equals("ADDA")){
                    double eta = Math.sqrt(nfg/neg);
                    Eg.scale(eta);
                    Fg.scale(1/eta);
                    diff = neg*nfg;
                }else {
                    diff = Math.min(neg,nfg);
                }
                numit++;
            }
            Psi = Hg;
        }

        if(numit==maxNumIt){
            System.out.print("Maximum Number of Iterations reached");
        }

        Map<String,Matrix> result = new HashMap<>();
        result.put("P",Psi);
        result.put("K",Fpp.add(1,Psi.mult(Fmp)));
        result.put("U",Fmm.add(1,Fmp.mult(Psi)));

        return result;
    }

    public static Matrix MomsFromFactorialMoms(Matrix fm){
        int n = fm.length();
        Matrix m = new Matrix(1,n,n);
        m.set(0,fm.get(0));

        for(int i=1;i<n;i++){
            Matrix a = new Matrix(1,i,i);
            for(int k=0;k<i;k++){
                a.set(k,k);
            }
            Matrix eh = poly(a);
            eh.scale(-1);
            Matrix eh_new = new Matrix(1,eh.numCols-2,eh.numCols-2);
            for(int k=0;k<eh.numCols-2;k++){
                eh_new.set(k,eh.get(eh.numCols-2-k));
            }
            m.set(i,fm.get(i)+eh_new.mult(Matrix.extractRows(m.transpose(),0,i,null)).get(0));
        }
        m.reshape(fm.numRows,fm.numCols);
        return m;
    }




    public static Matrix poly(Matrix x){
        int n = x.length();
        Matrix c= Matrix.concatColumns(new Matrix(1),new Matrix(1,n,n),null);
        for (int j=0;j<n;j++){
            for(int i=1;i<j+2;i++){
                c.set(i,c.get(i)-x.get(j)*c.get(i-1));
            }
        }
        return c;
    }

    public static Map<String, Map<Integer, Matrix>> MMAPPH1FCFS(Map<Integer,Matrix> D, Map<Integer,Matrix> sigma, Map<Integer,Matrix> S, Integer numOfQLMoms, Integer numOfQLProbs, Integer numOfSTMoms, Matrix stDistr, Boolean stDistrME, Boolean stDistrPH, Double prec, Matrix classes_){
        int K = D.size()-1;
        double precision = 1e-14;
        Matrix classes = new Matrix(1,K,K);
        for(int i=0;i<K;i++){
            classes.set(i,i);
        }
        if (prec!=null){
            precision = prec;
        }

        if(classes_!=null){
            classes = classes_;
        }

        Matrix D0 = D.get(0);
        int N = D0.numRows;
        Matrix Ia = Matrix.eye(N);
        Matrix Da = new Matrix(N,N,N*N);
        for(int q=0;q<K;q++){
            Da = Da.add(1,D.get(q+1));
        }
        Map<Integer,Matrix> beta = new HashMap<>();
        Matrix theta = CTMC.ctmc_solve(D0.add(1,Da));
        Matrix lambda = new Matrix(K,1,K);
        Matrix mu = new Matrix(K,1,K);
        Matrix Nsk = new Matrix(1,K,K);
        double ro = 0;
        for(int k=0;k<K;k++){
            lambda.set(k,theta.mult(D.get(k+1)).elementSum());
            beta.put(k,CTMC.ctmc_solve(S.get(k).add(-1,S.get(k).sumRows().mult(sigma.get(k)))));
            Matrix neg_sk = S.get(k).clone();
            neg_sk.scale(-1);
            mu.set(k,beta.get(k).mult(neg_sk).elementSum());
            Nsk.set(k,S.get(k).numRows);
            ro =ro+lambda.get(k)/mu.get(k);
        }

        Matrix alpha = theta.mult(Da);
        alpha.scale(1/lambda.elementSum());
        Matrix D0i = D0.inv().clone();
        D0i.scale(-1);

        Matrix Sa = S.get(0);
        Map<Integer,Matrix> sa = new HashMap<>();
        Map<Integer,Matrix> ba = new HashMap<>();
        Map<Integer,Matrix> sv = new HashMap<>();
        sa.put(0,sigma.get(0));
        ba.put(0,beta.get(0));
        Matrix sv0 = S.get(0).sumRows();
        sv0.scale(-1);
        sv.put(0,sv0);

        Map<Integer,Matrix> Pk = new HashMap<>();
        Pk.put(0,D0i.mult(D.get(1)));

        for(int q =1; q<K;q++){
            sa.put(q,new Matrix(1,sigma.get(0).length(),sigma.get(0).length()));
            ba.put(q,new Matrix(1,beta.get(0).length(),beta.get(0).length()));
            sv.put(q,new Matrix(sigma.get(0).length(),1,sigma.get(0).length()));
            Pk.put(q,D0i.mult(D.get(q+1)));
        }

        for (int k=1;k<K;k++){
            Sa = Sa.createBlockDiagonal(S.get(k));
            for(int q =0;q<K;q++){
                if(q==k){
                    sa.put(q,Matrix.concatColumns(sa.get(q),sigma.get(k),null));
                    ba.put(q,Matrix.concatColumns(ba.get(q),beta.get(k),null));
                    Matrix sk_neg_sum = S.get(k).sumRows();
                    sk_neg_sum.scale(-1);
                    sv.put(q,Matrix.concatRows(sv.get(q),sk_neg_sum,null));
                }else {
                    sa.put(q,Matrix.concatColumns(sa.get(q),new Matrix(sigma.get(k).numRows,sigma.get(k).numCols,0),null));
                    ba.put(q,Matrix.concatColumns(ba.get(q), new Matrix(beta.get(k).numRows,beta.get(k).numCols,0),null));
                    sv.put(q,Matrix.concatRows(sv.get(q),new Matrix(sigma.get(k).length(),1,0),null));
                }
            }
        }
        Matrix P = D0i.mult(Da);
        Matrix iVec = D.get(1).kron(sa.get(0));

        for(int k=1;k<K;k++){
            iVec = iVec.add(1,D.get(k+1).kron(sa.get(k)));
        }

        int Ns = Sa.numRows;
        Matrix Is = Matrix.eye(Ns);
        Matrix neg_Sa_row_sum = Sa.sumRows();
        neg_Sa_row_sum.scale(-1);
        Matrix Y0 = FluidFundamentalMatrices(Ia.kron(Sa),Ia.kron(neg_Sa_row_sum),iVec,D0,precision,null,null).get("P");
        Matrix T = Ia.kron(Sa).add(1,Y0.mult(iVec));
        Matrix pi0 = new Matrix(1,T.numRows,T.numRows);
        for(int k=0;k<K;k++){
            Matrix ba_mu = ba.get(k).clone();
            ba_mu.scale(1/mu.get(k));
            pi0 = pi0.add(1,theta.mult(D.get(k+1)).kron(ba_mu));
        }
        pi0 = pi0.mult(T);
        pi0.scale(-1);

        Matrix iT = T.inv();
        iT.scale(-1);
        Matrix oa = Matrix.ones(N,1);

        Map<String,Map<Integer,Matrix>> result = new HashMap<>();
        if(numOfSTMoms!=null){
            result.put("stNoms",new HashMap<>());
        }
        if(stDistr!=null){
            result.put("stDistr",new HashMap<>());
        }
        if(stDistrME){
            result.put("stDistrME_alpha",new HashMap<>());
            result.put("stDistrME_A",new HashMap<>());
        }
        if(stDistrPH){
            result.put("stDistrPH_alpha", new HashMap<>());
            result.put("stDistrPH_A", new HashMap<>());
        }

        if(numOfQLProbs!=null){
            result.put("ncDistr",new HashMap<>());
        }

        if(numOfQLMoms!=null){
            result.put("ncMoms",new HashMap<>());
        }



        for(int i=0;i<classes.length();i++){
            int k = (int)classes.get(i);
            Matrix clo = iT.mult(oa.kron(sv.get(k)));
            if(numOfSTMoms!=null){
                Matrix rtMoms = new Matrix(1,numOfSTMoms);
                for(int m=0;m<numOfSTMoms;m++){
                    double rtMoms_m = pi0.mult(Matrix.pow(iT,m)).mult(clo).get(0)/pi0.mult(clo).get(0);
                    rtMoms.set(m,factorial(m+1)*rtMoms_m);
                }
                result.get("stNoms").put(k,rtMoms);
            }

            if(stDistr!=null){
                Matrix cdf = new Matrix(0,0,0);
                for(int p=0;p<stDistr.length();p++){
                    Matrix Tt = T.clone();
                    Tt.scale(stDistr.get(p));
                    double pr = 1 - pi0.mult(Tt.expm()).mult(clo).get(0)/pi0.mult(clo).get(0);
                    Matrix cdf_ = new Matrix(1,1,1);
                    cdf_.set(0,0,pr);
                    if(p==0){
                        cdf = cdf_;
                    }else {
                        cdf = Matrix.concatColumns(cdf,cdf_,null);
                    }
                }
                result.get("stDistr").put(k,cdf);
            }

            if(stDistrME){
                Matrix Bm = SimilarityMatrixForVectors(clo.mult(pi0.mult(clo).inv()), Matrix.ones(N*Ns,1));
                Matrix Bmi = Bm.inv();
                Matrix A = Bm.mult(T).mult(Bmi);
                Matrix alpha_stDistrME = pi0.mult(Bmi);
                result.get("stDistrME_alpha").put(k,alpha);
                result.get("stDistrME_A").put(k,A);
            }


            if(stDistrPH){
                Matrix vv = pi0.mult(iT);

                List<Double> nz = new ArrayList<>();
                List<Integer> nz_index = new ArrayList<>();
                for(int n=0;n< vv.length();n++){
                    if(vv.get(n)>precision){
                        nz.add(vv.get(n));
                        nz_index.add(n);
                    }
                }
                Matrix delta = Matrix.diag(nz.stream().mapToDouble(Double::doubleValue).toArray());
                Matrix neg_T = T.clone();
                neg_T.scale(-1);
                Matrix cl = neg_T.mult(clo);
                cl.scale(pi0.mult(clo).get(0,0));
                Matrix alpha_stDistrPH = new Matrix(0,0,0);
                for(int n=0;n<nz_index.size();n++){
                    if(alpha_stDistrPH.length()==0){
                        alpha_stDistrPH = Matrix.extractRows(cl,nz_index.get(n),nz_index.get(n)+1,null);
                    }else {
                        Matrix.concatRows(alpha_stDistrPH,Matrix.extractRows(cl,nz_index.get(n),nz_index.get(n)+1,null),null);
                    }
                }
                Matrix A = new Matrix(nz_index.size(),nz_index.size(),(int)Math.pow(nz_index.size(),2));
                for(int n=0;n<nz_index.size();n++){
                    for(int m=0;m<nz_index.size();m++){
                        A.set(nz_index.get(n),nz_index.get(m),T.get(nz_index.get(n),nz_index.get(m)));
                    }
                }
                A = delta.inv().mult(A.transpose()).mult(delta);
                result.get("stDistrPH_alpha").put(k,alpha);
                result.get("stDistrPH_A").put(k,A);
            }

            if(numOfQLProbs!=null){
                Matrix value = new Matrix(1,numOfQLProbs,numOfQLProbs);
                Matrix jm = new Matrix(Ns,1,Ns);
                for(int n = (int)Nsk.sumRows(0,k).elementSum();n<Nsk.sumRows(0,k+1).elementSum();n++){
                    jm.set(n,0,1);
                }
                Matrix jmc = Matrix.ones(Ns,1);
                jmc = jmc.add(-1,jm);
                Matrix LmCurr = Matrix.lyap(T,D0.add(1,Da).add(-1,D.get(k+1)).kron(Is),Matrix.eye(N*Ns),null);
                value.set(0,1-ro+pi0.mult(LmCurr).mult(oa.kron(jmc)).get(0));
                for(int n=0;n<numOfQLProbs-1;n++){
                    Matrix LmPrev = LmCurr.clone();
                    LmCurr = Matrix.lyap(T,D0.add(1,Da).add(-1,D.get(k+1)).kron(Is),LmPrev.mult(D.get(k+1).kron(Is)),null);
                    value.set(n+1,pi0.mult(LmCurr).mult(oa.kron(jmc)).get(0)+pi0.mult(LmPrev).mult(oa.kron(jm)).get(0));
                }
                result.get("ncDistr").put(k,value);
            }

            if(numOfQLMoms!=null){
                Matrix jm = new Matrix(Ns,1,Ns);
                for(int n = (int)Nsk.sumRows(0,k).elementSum();n<Nsk.sumRows(0,k+1).elementSum();n++){
                    jm.set(n,0,1);
                }
                Map<Integer,Matrix> ELn = new HashMap<>();
                ELn.put(0,Matrix.lyap(T,D0.add(1,Da).kron(Is),Matrix.eye(N*Ns),null));
                Matrix qlMoms = new Matrix(1,numOfQLMoms,numOfQLMoms);
                for(int n=0;n<numOfQLMoms;n++){
                    int bino = 1;
                    Matrix Btag = new Matrix(N*Ns,N*Ns,(int) Math.pow(N*Ns,2));
                    for (int m=-1;m<=n-1;m++){
                        Btag = Btag.add(bino,ELn.get(m+1));
                        bino = bino*(n-m)/(m+2);
                    }
                    ELn.put(n+1,Matrix.lyap(T,D0.add(1,Da).kron(Is),Btag.mult(D.get(k+1).kron(Is)),null));
                    qlMoms.set(n,pi0.mult(ELn.get(n+1)).elementSum()+pi0.mult(Btag).mult(oa.kron(jm)).get(0));
                }
                result.get("ncMoms").put(k,qlMoms);
            }
        }


        return result;

    }

    public static Map<String,Map<Integer,Matrix>> MMAPPH1NPPR(Map<Integer,Matrix> D, Map<Integer,Matrix> sigma, Map<Integer,Matrix> S, Integer numOfQLMoms,Integer numOfQLProbs, Integer numOfSTMoms, Matrix stCdfPoints,Double prec, Integer erlMaxOrder_, Matrix classes_){
        int K = D.size()-1;
        int erlMaxOrder = 200;
        if(erlMaxOrder_!=null){
            erlMaxOrder = erlMaxOrder_;
        }
        double precision = 1e-14;
        if(prec!=null){
            precision = prec;
        }
        Matrix classes = new Matrix(1,K,K);
        for(int i=0;i<K;i++){
            classes.set(i,i);
        }
        if(classes_!=null){
            classes = classes_;
        }

        Matrix D0 = D.get(0);
        int N = D0.numRows;
        Matrix I = Matrix.eye(N);
        Matrix sD = new Matrix(N,N,N*N);

        for(int i=0;i<K+1;i++){
            sD = sD.add(1,D.get(i));
        }

        Matrix M = new Matrix(1,K,K);
        Map<Integer,Matrix> s = new HashMap<>();
        for(int i=0;i<K;i++){
            Matrix neg_si = S.get(i).clone();
            neg_si.scale(-1);
            s.put(i,neg_si.sumCols());
            M.set(i, M.get(i)+sigma.get(i).length());
        }
        Matrix QWMM = D0.clone();
        Matrix QWPP = new Matrix(N*(int)M.elementSum(),N*(int)M.elementSum(),N*(int)M.elementSum()*N*(int)M.elementSum());
        Matrix QWMP = new Matrix(N,N*(int)M.elementSum(),N*N*(int)M.elementSum());
        Matrix QWPM = new Matrix(N*(int)M.elementSum(),N,N*N*(int)M.elementSum());
        int kix =0;
        for(int i=0;i<K;i++){
            int bs = N*(int)M.get(i);
            QWPP.insert_sub_matrix(kix,kix,kix+bs,kix+bs,Matrix.eye(N).kron(S.get(i)));
            QWMP.insert_sub_matrix(0,kix,QWMP.numRows,kix+bs,D.get(i+1).kron(sigma.get(i)));
            QWPM.insert_sub_matrix(kix,0,kix+bs,QWPM.numCols,Matrix.eye(N).kron(s.get(i)));
            kix = kix+bs;
        }

        Map<String,Matrix> FluidFundamentalMatrices = FluidFundamentalMatrices(QWPP,QWPM,QWMP,QWMM,precision,null,null);

        Matrix Kw = FluidFundamentalMatrices.get("K");
        Matrix Uw = FluidFundamentalMatrices.get("U");
        Matrix neg_Kw = Kw.clone();
        neg_Kw.scale(-1);

        Matrix Ua = Matrix.ones(N,1).add(2,QWMP.mult(neg_Kw.inv()).sumRows());
        Matrix pm = new Matrix(0,0,0);
        Matrix.solve(Matrix.concatColumns(Uw,Ua,null).transpose(),Matrix.concatColumns(new Matrix(1,N,0),new Matrix(1),null).transpose(),pm);
        pm = pm.transpose();

        double ro = ((1-pm.elementSum()/2))/(pm.elementSum()+(1-pm.elementSum())/2);
        Matrix kappa = pm.clone();
        kappa.scale(1/pm.elementSum());

        Matrix pi = CTMC.ctmc_solve(sD);
        Matrix lambda = new Matrix(K,1,1);
        for(int i=0;i<K;i++){
            lambda.set(i,pi.mult(D.get(i+1)).elementSum());
        }

        Map<Integer,Matrix> Psiw = new HashMap<>();
        Map<Integer,Matrix> Qwmp = new HashMap<>();
        Map<Integer,Matrix> Qwzp = new HashMap<>();
        Map<Integer,Matrix> Qwpp = new HashMap<>();
        Map<Integer,Matrix> Qwmz = new HashMap<>();
        Map<Integer,Matrix> Qwpz = new HashMap<>();
        Map<Integer,Matrix> Qwzz = new HashMap<>();
        Map<Integer,Matrix> Qwmm = new HashMap<>();
        Map<Integer,Matrix> Qwpm = new HashMap<>();
        Map<Integer,Matrix> Qwzm = new HashMap<>();

        for(int k=0;k<K;k++){
            double Mlo;
            if(k==0){
                Mlo = 0;
            }else {
                Mlo = M.sumSubMatrix(0,1,0,k);
            }

            double Mhi = M.elementSum()-Mlo;
            int a = (int) (N*Mlo*Mhi+N*Mhi);
            Matrix Qkwpp = new Matrix(a,a,a*a);
            Matrix Qkwpz = new Matrix(a,N*(int)Mlo,a*N*(int)Mlo);
            Matrix Qkwpm = new Matrix(a,N,a*N);
            Matrix Qkwmz = new Matrix(N,N*(int)Mlo);
            Matrix Qkwmp = new Matrix(N,a,a*N);
            Matrix Dlo = D0.clone();
            for(int i=0;i<k-1;i++){
                Dlo = Dlo.add(1,D.get(i+1));
            }
            Matrix Qkwmm = Dlo;
            Matrix Qkwzp = new Matrix(N*(int)Mlo,a,N*(int)Mlo*a);
            Matrix Qkwzm = new Matrix(N*(int)Mlo,N,N*(int)Mlo*N);
            Matrix Qkwzz = new Matrix(N*(int)Mlo,N*(int)Mlo,N*(int)Mlo*N*(int)Mlo);
            kix = 0;
            for(int i=k;i<K;i++){
                int kix2=0;
                for(int j=0;j<k-1;j++){
                    int bs = (int) (N*M.get(j)*M.get(i));
                    int bs2= N*(int) M.get(j);
                    Qkwpp.insert_sub_matrix(kix,kix,kix+bs,kix+bs,Matrix.eye(N).kron(Matrix.eye((int)M.get(j)).kron(S.get(i))));
                    Qkwpz.insert_sub_matrix(kix,kix2,kix+bs,kix2+bs2,Matrix.eye(N).kron(Matrix.eye((int)M.get(j)).kron(s.get(i))));
                    Qkwzp.insert_sub_matrix(kix2,kix,kix2+bs2,kix+bs,D.get(i+1).kron(Matrix.eye((int)M.get(j)).kron(sigma.get(i))));
                    kix = kix+bs;
                    kix2 = kix2+bs2;
                }
            }
            for(int i=k;i<K;i++){
                int bs = N*(int) M.get(i);
                Qkwpp.insert_sub_matrix(kix,kix,kix+bs,kix+bs,Matrix.eye(N).kron(S.get(i)));
                Qkwpm.insert_sub_matrix(kix,0,kix+bs,Qkwpm.numCols,Matrix.eye(N).kron(s.get(i)));
                Qkwmp.insert_sub_matrix(0,kix,Qkwmp.numRows,kix+bs,D.get(i+1).kron(sigma.get(i)));
                kix = kix+bs;
            }
            kix = 0;
            for(int j=0;j<k-1;j++){
                int bs = N*(int)M.get(j);
                Qkwzz.insert_sub_matrix(kix,kix,kix+bs,kix+bs,Dlo.kron(Matrix.eye((int) M.get(j))).add(1,Matrix.eye(N).kron(S.get(j))));
                Qkwzm.insert_sub_matrix(kix,0,kix+bs,Qkwzm.numCols,Matrix.eye(N).kron(s.get(j)));
                kix = kix +bs;
            }
            Matrix neg_Qkwzz = Qkwzz.clone();
            neg_Qkwzz.scale(-1);
            Matrix Psikw = FluidFundamentalMatrices(Qkwpp.add(1,Qkwpz.mult(neg_Qkwzz.inv()).mult(Qkwzp)),Qkwpm.add(1,Qkwpz.mult(neg_Qkwzz.inv()).mult(Qkwzm)),Qkwmp,Qkwmm,precision,null,null).get("P");
            Psiw.put(k,Psikw);
            Qwzp.put(k,Qkwzp);
            Qwmp.put(k,Qkwmp);
            Qwpp.put(k,Qkwpp);
            Qwmz.put(k,Qkwmz);
            Qwpz.put(k,Qkwpz);
            Qwzz.put(k,Qkwzz);
            Qwmm.put(k,Qkwmm);
            Qwpm.put(k,Qkwpm);
            Qwzm.put(k,Qkwzm);
        }

        double lambdaS = lambda.elementSum();
        Map<Integer,Matrix> phi = new HashMap<>();
        Matrix neg_D0 = D0.clone();
        neg_D0.scale(-1);
        Matrix phi0 = kappa.mult(neg_D0);
        phi0.scale((1-ro)/lambdaS);
        phi.put(0,phi0);

        Map<Integer,Matrix> q0 = new HashMap<>();
        Map<Integer,Matrix> qL = new HashMap<>();
        q0.put(0, new Matrix(0,0,0));
        qL.put(0, new Matrix(0,0,0));

        for(int k=0;k<K-1;k++){
            Matrix sDk = D.get(0);
            for(int j=0;j<=k;j++){
                sDk = sDk.add(1,D.get(j+1));
            }
            double pk=0;
            for(int j=0;j<=k;j++){
                pk = pk+ lambda.get(j);
            }
            pk = pk/lambdaS - (1-ro)*kappa.mult(sDk.sumRows()).get(0)/lambdaS;
            Matrix Qwzpk = Qwzp.get(k+1);
            int vix= 0;
            Map<Integer,Matrix> Ak = new HashMap<>();
            for(int ii=0;ii<=k;ii++){
                int bs = (int)(N*M.get(ii));
                Matrix V1 = Matrix.extractRows(Qwzpk,vix,vix+bs,null);
                Matrix a = sDk.kron(Matrix.eye((int)M.get(ii)));
                a.scale(-1);
                a.add(-1,I.kron(S.get(ii)));
                Ak.put(ii,I.kron(sigma.get(ii))).mult(a.inv()).mult(I.kron(s.get(ii))).add(1,V1.mult(Psiw.get(k+1)));
                vix = vix +bs;
            }
            Matrix Qwmpk = Qwmp.get(k+1);
            Matrix Bk = Qwmpk.mult(Psiw.get(k+1));
            Matrix ztag = phi.get(0).mult(neg_D0.inv().mult(D.get(k+1)).mult(Ak.get(k)).add(-1,Ak.get(0)).add(1,D0.inv()).mult(Bk));
            for(int i=0;i<k-1;i++){
                ztag = ztag.add(1,phi.get(i+1).mult(Ak.get(i).add(-1,Ak.get(i+1)))).add(1,phi.get(0).mult(neg_D0.inv()).mult(D.get(i+1))).mult(Ak.get(i));
            }
            Matrix Mx = Matrix.eye(Ak.get(k).numCols).add(-1,Ak.get(k));
            Mx.insert_sub_matrix(0,0,Mx.numRows,1,Matrix.ones(N,1));
            phi.put(k+1,Matrix.concatColumns(new Matrix(pk),Matrix.extractRows(ztag.transpose(),1, ztag.length(), null).transpose(),null).mult(Mx.inv()));
            q0.put(k+1,phi.get(0).mult(neg_D0.inv()));
            qL.put(k+1,new Matrix(0,0,0));
            for(int ii=1;ii<k;ii++){
                Matrix a = sDk.kron(Matrix.eye((int)M.get(ii)));
                a.scale(-1);
                Matrix qLii = phi.get(ii+1).add(-1,phi.get(ii)).add(1,phi.get(0).mult(neg_D0.inv()).mult(D.get(ii+1))).mult(I.kron(sigma.get(ii))).mult(a.add(-1,I.kron(S.get(ii))));
                qL.put(k+1, Matrix.concatColumns(qL.get(k+1),qLii,null));
            }
        }

        Map<String,Map<Integer,Matrix>> result = new HashMap<>();
        if(numOfSTMoms!=null){
            result.put("stMoms",new HashMap<>());
        }
        if(stCdfPoints!=null){
            result.put("stDistr",new HashMap<>());
        }
        if(numOfQLMoms!=null){
            result.put("ncMoms",new HashMap<>());
        }
        if(numOfQLProbs!=null){
            result.put("ncDistr",new HashMap<>());
        }
        for(int g=0;g< classes.length();g++){
            int k = (int)classes.get(g);
            Matrix sD0k = D0;
            for(int i=0;i<=k-1;i++){
                sD0k = sD0k.add(1,D.get(i+1));
            }
            if(k<K){
                Matrix neg_Qwzz_k = Qwzz.get(k).clone();
                neg_Qwzz_k.scale(-1);
                Kw = Qwpp.get(k).add(1,Qwpz.get(k).mult(neg_Qwzz_k.inv()).mult(Qwzp.get(k))).add(1,Psiw.get(k).mult(Qwmp.get(k)));
                Matrix BM = new Matrix(0,0,0);
                Matrix CM = I.kron(s.get(0));
                Matrix DM = new Matrix(0,0,0);
                for(int i=0;i<=k-1;i++){
                    BM.createBlockDiagonal(I.kron(S.get(i)));
                    DM.createBlockDiagonal(D.get(k+1).kron(Matrix.eye((int)M.get(i))));
                    if(i!=0){
                        CM = Matrix.concatRows(CM,I.kron(s.get(i)),null);
                    }
                }
                Matrix Kwu = Matrix.concatRows(Matrix.concatColumns(Kw,Qwpz.get(k).add(1,Psiw.get(k).mult(Qwmz.get(k))).mult(Matrix.negative(Qwzz.get(k)).inv()).mult(DM),null),Matrix.concatColumns(new Matrix(BM.numRows,Kw.numRows,0),BM,null),null);
                Matrix Bwu = Matrix.concatRows(Psiw.get(k).mult(D.get(k+1)),CM,null);
                Matrix iniw;
                Matrix pwu;
                if(k>0){
                    iniw = Matrix.concatColumns(q0.get(k).mult(Qwmp.get(k)).add(1,qL.get(k).mult(Qwzp.get(k))),qL.get(k).mult(DM),null);
                    pwu = q0.get(k).mult(D.get(k+1));
                }else {
                    iniw = pm.mult(Qwmp.get(k));
                    pwu = pm.mult(D.get(k+1));
                }
                double norm = pwu.elementSum()+iniw.mult(Matrix.negative(Kwu).inv()).mult(Bwu).elementSum();
                pwu.scale(1/norm);
                iniw.scale(1/norm);
                int KN = Kwu.numRows;
                Matrix Qspp = new Matrix((int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()),(int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()),(int)Math.pow(KN+N*M.sumCols(k+1,M.numCols).elementSum(),2));
                Matrix Qspm = new Matrix((int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()),N,N*(int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()));
                Matrix Qsmp = new Matrix(N,(int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()),N*(int)(KN+N*M.sumCols(k+1,M.numCols).elementSum()));
                Matrix Qsmm = sD0k.add(1,D.get(k+1));
                kix=0;
                for(int i=k+1;i<K;i++){
                    int bs = N*(int)M.get(i);
                    Qspp.insert_sub_matrix(KN+kix,KN+kix,KN+kix+bs,KN+kix+bs,I.kron(S.get(i)));
                    Qspm.insert_sub_matrix(KN+kix,0,KN+kix+bs,Qspm.numCols,I.kron(s.get(i)));
                    Qsmp.insert_sub_matrix(0,KN+kix,Qsmp.numRows,KN+kix+bs,D.get(i+1).kron(sigma.get(i)));
                    kix = kix+bs;
                }
                Qspp.insert_sub_matrix(0,0,KN,KN,Kwu);
                Qspm.insert_sub_matrix(0,0,KN, Qspm.numCols, Bwu);
                Matrix inis = Matrix.concatColumns(iniw,new Matrix(1,N*(int)M.sumCols(K+1,M.numCols).elementSum()),null);

                Matrix Psis = FluidFundamentalMatrices(Qspp,Qspm,Qsmp,Qsmm,precision,null,null).get("P");

                if(numOfSTMoms!=null){
                    Map<Integer,Matrix> Pn = new HashMap<>();
                    Pn.put(0,Psis);
                    Matrix wtMons = new Matrix(1,numOfSTMoms,numOfSTMoms);
                    for (int n=0;n<numOfSTMoms;n++){
                        Matrix A = Qspp.add(1,Psis.mult(Qsmp));
                        Matrix B = Qsmm.add(1,Qsmp.mult(Psis));
                        Matrix C = Pn.get(n);
                        C.scale(-2*n);
                        int bino = 1;
                        for(int i=0;i<=n-1;i++){
                            bino = bino*(n-i+1)/(i+1);
                            C = C.add(bino,Pn.get(i+1).mult(Qsmp).mult(Pn.get(n-i+1)));
                        }
                        Matrix P = Matrix.lyap(A,B,C,null);
                        Pn.put(n+1,P);
                        wtMons.set(n,inis.mult(P).elementSum()*Math.pow(-1,n)/Math.pow(2,n));
                    }
                    Map<Integer,Matrix> Pnr = new HashMap<>();
                    Pnr.put(0,sigma.get(k));
                    Pnr.get(0).scale(inis.mult(Pn.get(1)).elementSum());
                    Matrix rtMoms = new Matrix(1,numOfSTMoms,numOfSTMoms);
                    for(int n=0;n<numOfSTMoms;n++){
                        Matrix P = Pnr.get(n).mult(Matrix.negative(S.get(k)).inv());
                        P.scale(n);
                        P.add(Math.pow(-1,n)*inis.mult(Pn.get(n+1)).elementSum()/Math.pow(2,n),sigma.get(k));
                    }
                    result.get("rtMoms").put(k,rtMoms);
                }

                if(stCdfPoints!=null){
                    Matrix res = new Matrix(1,0,0);
                    for(int o=0;o<stCdfPoints.length();o++){
                        int t = (int) stCdfPoints.get(o);
                        int L = erlMaxOrder;
                        double lambdae = L/(double)t/2.0;
                        Matrix Psie = FluidFundamentalMatrices(Qspp.add(-lambdae,Matrix.eye(Qspp.numRows)),Qspm,Qsmp,Qsmm.add(-lambdae,Matrix.eye(Qsmm.numRows)),precision,null,null).get("P");
                        Map<Integer,Matrix> Pn = new HashMap<>();
                        Pn.put(0,Psis);
                        double pr = (pwu.elementSum()+inis.mult(Psie).elementSum())*(1-sigma.get(k).mult(Matrix.pow(Matrix.eye(S.get(k).numRows).inv().add(-1/2.0/lambdae,S.get(k)),L)).elementSum());
                        for(int n=0;n<=L-1;n++){
                            Matrix A = Qspp.add(1,Psie.mult(Qsmp)).add(-lambdae,Matrix.eye(Qspp.length()));
                            Matrix B = Qsmm.add(1,Qsmp.mult(Psie)).add(-lambdae,Matrix.eye(Qsmm.length()));
                            Matrix C = Pn.get(n);
                            C.scale(2*lambdae);
                            for(int i=0;i<n-1;i++){
                                C = C.add(1,Pn.get(i+1).mult(Qsmp).mult(Pn.get(n-i+1)));
                            }
                            Matrix P = Matrix.lyap(A,B,C,null);
                            Pn.put(n+1,P);
                            pr = pr+ inis.mult(P).elementSum()*(1-sigma.get(k).mult(Matrix.pow(Matrix.eye(S.get(k).length()).inv().add(-1/2.0/lambdae,S.get(k)),L-n)).elementSum());
                        }
                        res = Matrix.concatColumns(res,new Matrix(pr),null);
                    }
                    result.get("stDistr").put(k,res);
                }
                if(numOfQLMoms!=null||numOfQLProbs!=null){
                    Matrix W = Matrix.negative(sD.add(-1,D.get(k+1)).kron(Matrix.eye((int)M.get(k)))).add(-1,I.kron(S.get(k))).inv().mult(D.get(k+1).kron(Matrix.eye((int)M.get(k))));
                    Matrix iW = Matrix.eye(W.numCols).add(-1,W).inv();
                    Matrix w = Matrix.eye(N).kron(sigma.get(k));
                    Matrix omega = Matrix.negative(sD.add(-1,D.get(k+1)).kron(Matrix.eye((int)M.get(k)))).add(-1,I.kron(S.get(k))).inv().mult(I.kron(s.get(k)));
                    if(numOfQLMoms!=null){
                        Map<Integer,Matrix> Psii = new HashMap<>();
                        Psii.put(0,Psis);
                        Map<Integer,Matrix> QLDPn = new HashMap<>();
                        QLDPn.put(0,inis.mult(Psii.get(0)).mult(w).mult(iW));
                        for (int n=0;n<numOfQLMoms;n++){
                            Matrix A = Qspp.add(1,Psis.mult(Qsmp));
                            Matrix B = Qsmm.add(1,Qsmp.mult(Psis));
                            Matrix C = Psii.get(n).mult(D.get(k+1));
                            C.scale(n);
                            int bino=1;
                            for(int i=0;i<=n-1;i++){
                                bino = bino*(n-i+1)/(i+1);
                                C = C.add(bino,Psii.get(i+1).mult(Qsmp).mult(Psii.get(n-i+1)));
                            }
                            Matrix P = Matrix.lyap(A,B,C,null);
                            Psii.put(n+1,P);
                            QLDPn.put(n+1,inis.mult(P).mult(w).mult(iW).add(n,QLDPn.get(n).mult(iW).mult(W)));
                        }
                        for(int n=-1;n<numOfQLMoms;n++){
                            QLDPn.put(n+1,QLDPn.get(n+1).add(1,pwu.mult(w).mult(Matrix.pow(iW,n+2)).mult(Matrix.pow(W,n+1)).mult(omega)));
                        }
                        Map<Integer,Matrix> QLPn = new HashMap<>();
                        QLPn.put(0,pi);
                        Matrix qlMOms = new Matrix(1,numOfQLMoms,numOfQLMoms);
                        Matrix iTerm = Matrix.ones(N,1).mult(pi).add(-1,sD).inv();
                        for(int n=0;n<numOfQLMoms;n++){
                            double sumP = QLDPn.get(n+1).elementSum()+n*(QLDPn.get(n).add(-1/lambda.get(k),QLPn.get(n).mult(D.get(k+1))).mult(iTerm).mult(D.get(k+1).sumRows())).get(0);
                            Matrix P = Matrix.scale_mult(pi,sumP).add(n,QLPn.get(n).mult(D.get(k+1)).add(-lambda.get(k),QLDPn.get(n)).mult(iTerm));
                            QLPn.put(n+1,P);
                            qlMOms.set(n,P.elementSum());
                        }
                        qlMOms = MomsFromFactorialMoms(qlMOms);
                        result.get("ncMoms").put(k,qlMOms);
                    }
                    if(numOfQLProbs!=null){
                        Matrix Psid = FluidFundamentalMatrices(Qspp,Qspm,Qsmp,sD0k,precision,null,null).get("P");
                        Map<Integer,Matrix> Pn = new HashMap<>();
                        Pn.put(0,Psid);
                        Matrix XDn = inis.mult(Psid).mult(w);
                        Matrix dqlProbs = XDn.add(1,pwu.mult(w)).mult(omega);
                        for(int n=0;n<numOfQLProbs-1;n++){
                            Matrix A = Qspp.add(1,Psid.mult(Qsmp));
                            Matrix B = sD0k.add(1,Qsmp.mult(Psid));
                            Matrix C = Pn.get(n).mult(D.get(k+1));
                            for(int i=0;i<=n-1;i++){
                                C = C.add(1,Pn.get(i+1).mult(Qsmp).mult(Pn.get(n-i+1)));
                            }
                            Matrix P = Matrix.lyap(A,B,C,null);
                            Pn.put(n+1,P);
                            XDn = XDn.mult(W).add(1,inis.mult(P).mult(w));
                            dqlProbs = Matrix.concatRows(dqlProbs,XDn.add(1,pwu.mult(w).mult(Matrix.pow(W,n))).mult(omega),null);
                        }
                        Matrix iTerm = Matrix.negative(sD.add(-1,D.get(k+1))).inv();
                        Matrix qlProbs = Matrix.scale_mult(Matrix.extractRows(dqlProbs,0,1,null).mult(iTerm),lambda.get(k));
                        for (int n=0;n<numOfQLProbs-1;n++){
                            Matrix P = Matrix.extractRows(qlProbs,n,n+1,null).mult(D.get(k+1)).add(1,Matrix.scale_mult(Matrix.extractRows(dqlProbs,n+1,n+2,null).add(-1,Matrix.extractRows(dqlProbs,n,n+1,null)),lambda.get(k))).mult(iTerm);
                            qlProbs = Matrix.concatRows(qlProbs,P,null);
                        }
                        result.get("ncDistr").put(k,qlProbs.sumRows().transpose());
                    }
                }
            }else if(k==K){
                if(numOfSTMoms!=null||stCdfPoints!=null){
                    Kw = Qwpp.get(k).add(1,Qwpz.get(k).mult(Matrix.negative(Qwzz.get(k)).inv()).mult(Qwzp.get(k))).add(1,Psiw.get(k).mult(Qwmp.get(k)));
                    Matrix AM = new Matrix(0,0,0);
                    Matrix BM = new Matrix(0,0,0);
                    Matrix CM = new Matrix(0,s.get(0).numCols,0);
                    Matrix DM = new Matrix(0,0,0);
                    for(int i=0;i<k-1;i++){
                        AM = AM.createBlockDiagonal(new Matrix(N,1).kron(Matrix.eye((int)M.get(i)).kron(s.get(k))));
                        BM = BM.createBlockDiagonal(S.get(i));
                        CM = Matrix.concatRows(CM,s.get(i),null);
                        DM = DM.createBlockDiagonal(D.get(k+1).kron(Matrix.eye((int)M.get(i))));
                    }
                    Matrix Z = Matrix.concatRows(Matrix.concatColumns(Kw,Matrix.concatRows(AM,new Matrix(N*(int)M.get(k), AM.numCols),null),null),Matrix.concatColumns(new Matrix(BM.numRows,Kw.numCols),BM,null),null);
                    Matrix z = Matrix.concatRows(Matrix.concatRows(new Matrix(AM.numRows,1,0),Matrix.ones(N,1).kron(s.get(k)),null),CM,null);
                    Matrix iniw = Matrix.concatColumns(q0.get(k).mult(Qwmp.get(k)).add(1,qL.get(k).mult(Qwzp.get(k))),new Matrix(1, BM.numRows),null);
                    Matrix zeta = Matrix.scale_mult(iniw,1/iniw.mult(Matrix.negative(Z).inv()).mult(z).elementSum());
                    if(numOfSTMoms!= null){
                        Matrix rtMomsH = new Matrix(1,numOfSTMoms,numOfSTMoms);
                        for(int i=0;i<numOfSTMoms;i++){
                            rtMomsH.set(i,factorial(i+1)*zeta.mult(Matrix.pow(Matrix.negative(Z).inv(),i+1)).mult(Z).get(0));
                        }
                        result.get("stMoms").put(k,rtMomsH);
                    }
                    if(stCdfPoints!=null){
                        Matrix rtDistr = zeta.mult(Matrix.negative(Z).inv()).mult(Matrix.eye(Z.numCols).add(-1,Z.mult(Matrix.scale_mult(Z,stCdfPoints.get(0)).expm()))).mult(z);
                        for(int i =1;i<stCdfPoints.length();i++){
                            rtDistr = Matrix.concatColumns(rtDistr,zeta.mult(Matrix.negative(Z).inv()).mult(Matrix.eye(Z.numCols).add(-1,Z.mult(Matrix.scale_mult(Z,stCdfPoints.get(i)).expm()))).mult(z),null);
                        }
                        result.get("stDistr").put(k,rtDistr);
                    }
                }

                if(numOfSTMoms!=null||stCdfPoints!=null){
                    Matrix L = new Matrix(N*(int) M.elementSum(),N*(int) M.elementSum(),(int) Math.pow(N*M.elementSum(),2));
                    Matrix B = new Matrix(N*(int) M.elementSum(),N*(int) M.elementSum(),(int) Math.pow(N*M.elementSum(),2));
                    Matrix F = new Matrix(N*(int) M.elementSum(),N*(int) M.elementSum(),(int) Math.pow(N*M.elementSum(),2));
                    kix = 0;
                    for(int i=0;i<K;i++){
                        int bs = N*(int)M.get(i);
                        F.insert_sub_matrix(kix,kix,kix+bs,kix+bs,D.get(k+1).kron(Matrix.eye((int)M.get(i))));
                        L.insert_sub_matrix(kix,kix,kix+bs,kix+bs,sD0k.kron(Matrix.eye((int)M.get(i)).add(1,I.kron(S.get(i)))));
                        if(i<K){
                            L.insert_sub_matrix(kix,1+N*(int) Matrix.extractRows(M,0,k,null).elementSum(),kix+bs,L.numCols,I.kron(s.get(i).mult(sigma.get(k))));
                        }else {
                            B.insert_sub_matrix(kix,1+N*(int) Matrix.extractRows(M,0,k,null).elementSum(),kix+bs,B.numCols,I.kron(s.get(i).mult(sigma.get(k))));
                        }
                        kix = kix+bs;
                    }
                    Matrix R = QBDFundamentalMatrices (B, L, F, precision,null ,null,null).get("R");
                    Matrix P0 = Matrix.concatColumns(qL.get(k),q0.get(k).mult(I.kron(sigma.get(k))),null);
                    P0.scale(1/P0.mult(Matrix.eye(R.numRows).add(-1,R).inv()).elementSum());

                    if(numOfQLMoms!=null){
                        Matrix qlMoms = new Matrix(1,numOfQLMoms,numOfQLMoms);
                        for(int i=0;i<numOfQLMoms;i++){
                            qlMoms.set(i,Matrix.scale_mult(P0.mult(Matrix.pow(R,i)).mult(Matrix.pow(Matrix.eye(R.numRows).add(-1,R),i+1).inv()),factorial(i)).elementSum());
                        }
                        result.get("ncMoms").put(k, MomsFromFactorialMoms(qlMoms));
                    }

                    if(numOfQLProbs!=null){
                        Matrix qlProbs = P0.clone();
                        for (int i=0;i<numOfQLProbs-1;i++) {
                            qlProbs = Matrix.concatRows(qlProbs, P0.mult(Matrix.pow(R, i)), null);
                        }
                        result.get("ncDistr").put(k,qlProbs.sumRows());
                    }
                }
            }
        }
        return result;
    }

    public static Matrix SimilarityMatrixForVectors (Matrix vecA, Matrix vecB){
        int m = vecA.length();
        Matrix neg_vecA = vecA.clone();
        neg_vecA.scale(-1);
        double[] neg_vecA_array = neg_vecA.toArray1D();
        Integer[] ix = IntStream.range(0, neg_vecA_array.length).boxed().sorted(Comparator.comparingDouble(i -> neg_vecA_array[i])).toArray(Integer[]::new);
        Matrix P = new Matrix(m,m,m);
        for(int i=0;i<m;i++){
            P.set(i,ix[i],1);
        }
        Matrix cp = P.mult(vecA);

        Matrix B = new Matrix(m,m,m*m);
        for(int i=0;i<m;i++){
            double cp_sum = 0 ;
            for(int j=0;j<=i;j++){
                cp_sum = cp_sum+cp.get(j,0);
            }
            for(int j=0;j<=i;j++){
                B.set(i,j,vecB.get(i)/cp_sum);
            }
        }
        return B.mult(P);
    }

}
