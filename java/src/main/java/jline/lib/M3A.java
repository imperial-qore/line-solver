package jline.lib;

import jline.util.Matrix;

import java.util.*;

import static jline.lib.KPCToolbox.*;
import static org.apache.commons.math3.util.CombinatoricsUtils.factorial;

public class M3A {
    public static Map<Integer, Matrix> aph_fit(double e1, double e2, double e3, int nmax){
        Map<Integer,Matrix> APH = new HashMap<>();
        if(Double.isInfinite(e2)||Double.isInfinite(e3)){
            return map_exponential(e1);
        }

        double n2 = e2/Math.pow(e1,2);
        double n3 = e3/e1/e2;

        boolean n2_feas = false;
        boolean n3_ubfeas = false;
        boolean n3_lbfeas = false;
        double n = 1;
        double un = 0;
        double un_1 = un;
        while ((!n2_feas|| !n3_lbfeas || !n3_ubfeas  ) && n < nmax){
            n++;
            double pn = ((n+1)*(n2-2)/(3*n2*(n-1)))*(-2*Math.sqrt(n+1)/Math.sqrt(4*(n+1)-3*n*n2) - 1);
            double an = (n2 - 2) / (pn*(1-n2) + Math.sqrt(Math.pow(pn,2)+pn*n*(n2-2)/(n-1)));
            double ln = ((3+an)*(n-1)+2*an)/((n-1)*(1+an*pn)) - (2*an*(n+1))/(2*(n-1)+an*pn*(n*an+2*n-2));

            un = (1/(Math.pow(n,2)*n2))*(2*(n-2)*(n*n2-n-1)*Math.sqrt(1+n*(n2-2)/(n-1))+(n+2)*(3*n*n2-2*n-2));

            if( n2 >= (n+1)/n && n2 <= (n+4)/(n+1)){
                n2_feas = true;
                if (n3 >= ln) {
                    n3_lbfeas = true;
                }
            }else if(n2 >= (n+4)/(n+1)){
                n2_feas = true;
                if (n3 >= n2*(n+1)/n) {
                    n3_lbfeas = true;
                }
            }

            if(n2 >= (n+1)/n && n2 <= n/(n-1)){
                n2_feas = true;
                if(n3 <= un){
                    n3_ubfeas = true;
                }
            }else if(n2 >= n/(n-1)){
                n2_feas = true;
                if(n3 < Double.POSITIVE_INFINITY){
                    n3_ubfeas = true;
                }
            }
        }

        if ((!n2_feas  || !n3_lbfeas  || ! n3_ubfeas  ) || (n == nmax)){
            System.out.print("'cannot match moment set exactly'");
            n2 = (n+1)/n;
            n3 =  2*n2-1;
        }

        if(n2 <= n/(n-1) || n3 <= 2*n2-1){
            double b = 2*(4-n*(3*n2-4))/(n2*(4+n-n*n3)+Math.sqrt(n*n2)*Math.sqrt(12*Math.pow(n2,2)*(n+1)+16*n3*(n+1)+n2*(n*(n3-15)*(n3+1)-8*(n3+3))));
            double a = (b*n2-2)*(n-1)*b/((b-1)*n);
            double p = (b-1)/a;

            double lambda = 1;
            double mu = lambda*(n-1)/a;
            Matrix alpha = new Matrix(1,(int)n,(int)n);
            alpha.set(0,0,p);
            alpha.set(0,(int)n-1,1-p);
            Matrix T = new Matrix((int)n,(int)n,(int)Math.pow(n,2));
            for(int i=0;i<n-1;i++){
                T.set(i,i,-mu);
                T.set(i,i+1,mu );

            }
            T.set((int)n-1,(int)n-1,-lambda);
            Matrix neg_T = T.clone();
            neg_T.scale(-1);
            Matrix one = new Matrix(1,(int)n-1,(int)n-1);
            for(int i=0;i<n-1;i++){
                one.set(0,i,1);
            }
            APH = map_scale(map_normalize(T,neg_T.mult(one).mult(alpha)).get(0),map_normalize(T,neg_T.mult(one).mult(alpha)).get(1),e1);
        }else if(n2 > n/(n-1) && n3 > un_1){
            double K1 = n-1;
            double K2 = n-2;
            double K3 = 3*n2-2*n3;
            double K4 = n3-3;
            double K5 = n-n2;
            double K6 = 1+n2-n3;
            double K7 = n+n2-n*n2;
            double K8 = 3+3*Math.pow(n2,2)+n3-3*n2*n3;
            double K9 = 108*Math.pow(K1,2)*(4*Math.pow(K2,2)*K3*Math.pow(n,2)*n2+Math.pow(K1,2)*K2*Math.pow(K4,2)*n*Math.pow(n2,2)+4*K1*K5*(Math.pow(K5,2)-3*K2*K6*n*n2)+Math.sqrt(-16*Math.pow(K1,2)*Math.pow(K7,6)+Math.pow(4*K1*Math.pow(K5,3)+Math.pow(K1,2)*K2*Math.pow(K4,2)*n*Math.pow(n2,2)+4*K2*n*n2*(K4*Math.pow(n,2)-3*K6*n2+K8*n),2)));
            double K10 = Math.pow(K4,2)/(4*Math.pow(K3,2)) - K5/(K1*K3*n2);
            double K11 = Math.pow(2,1.0/3)*(3*Math.pow(K5,2)+K2*(K3+2*K4)*n*n2)/(K3*Math.pow(K9,1.0/3)*n2);
            double K12 = Math.pow(K9,1.0/3) / (3*Math.pow(2,7.0/3)*Math.pow(K1,2)*K3*n2);
            double K13 = Math.sqrt(K10 + K11 + K12);
            double K14 = (6*K1*K3*K4*K5+4*K2*Math.pow(K3,2)*n-Math.pow(K1,2)*Math.pow(K4,2)*n2) / (4*Math.pow(K1,2)*Math.pow(K3,3)*K13*n2);
            double K15 = -K4/(2*K3);
            double K16 = Math.sqrt(2*K10 - K11 -K12 -K14);
            double K17 = Math.sqrt(2*K10 - K11 -K12 +K14);
            double K18 = 36*Math.pow(K5,3) + 36*K2*K4*K5*n*n2 + 9*K1*K2*Math.pow(K4,2)*n*Math.pow(n2,2) - Math.sqrt(81*Math.pow(4*Math.pow(K5,2)+4*K2*K4*K5*n*n2+K1*K2*Math.pow(K4,2)*n*Math.pow(n2,2),2)-48*Math.pow(3*Math.pow(K5,2)+2*K2*K4*n*n2,3));
            double K19 = -K5/(K1*K4*n2) -Math.pow(2,2.0/3)*(3*Math.pow(K5,2)+2*K2*K4*n*n2)/(Math.pow(3,1.0/3)*K1*K4*n2*Math.pow(K18,1.0/3)) - Math.pow(K18,1.0/3)/(Math.pow(6,2.0/3)*K1*K4*n2);
            double K20 = 6*K1*K3*K4*K5 + 4*K2*Math.pow(K3,2)*n - Math.pow(K1,2)*Math.pow(K4,2)*n2;
            double K21 = K11 + K12 + K5/(2*n*K1*K3);
            double K22 = Math.sqrt(3*Math.pow(K4,2)/(4*Math.pow(K3,2)) - 3*K5/(K1*K3*n2) + Math.sqrt(4*Math.pow(K21,2) - n*K2/(n2*Math.pow(K1,2)*K3)));
            double f = 0;
            if(n3 > un_1 && n3 < 3*n2/2){
                f = K13+K15-K17;
            }else if(n3 == 2*n2/2){
                f = K19;
            }else if(n3 > 3*n2/2 && K20 > 0){
                f = -K13+K15+K16;
            }else if(K20 == 0){
                f = K15+K22;
            }else if(K20<0){
                f = K13+K15+K17;
            }
            double a = 2*(f-1)*(n-1)/((n-1)*(n2*Math.pow(f,2)-2*f+2)-n);
            double p = (f-1)*a;
            double lambda = 1;
            double mu = lambda*(n-1)/a;
            Matrix alpha = new Matrix(1,(int)n,(int)n);
            alpha.set(0,0,p);
            alpha.set(0,1,1-p);
            Matrix T = new Matrix((int)n,(int)n,(int)Math.pow(n,2));
            for(int i=0;i<n-1;i++){
                T.set(i,i,-mu);
                T.set(i,i+1,mu );
            }
            T.set((int)n-1,(int)n-1,-mu);
            T.set(0,0,-lambda);
            T.set(0,1,lambda);
            Matrix neg_T = T.clone();
            neg_T.scale(-1);
            Matrix one = new Matrix(1,(int)n,(int)n);
            for(int i=0;i<n;i++){
                one.set(0,i,1);
            }
            APH = map_scale(map_normalize(T,neg_T.mult(one).mult(alpha)).get(0),map_normalize(T,neg_T.mult(one).mult(alpha)).get(1),e1);
        }else {
            System.out.print("moment set cannot be matched with an APH distribution");
        }

        return APH;
    }

    public static Map<Integer, Matrix> aph_fit(double e1, double e2, double e3){
        return aph_fit(e1,e2,e3,10);
    }

    public static Map<Integer, Map<Integer, Matrix>> aph2_fitall(double M1, double M2, double M3){
        Map<Integer,Map<Integer,Matrix>> APHS = new HashMap<>();
        double degentol = 1e-8;
        double SCV = (M2-Math.pow(M1,2))/Math.pow(M1,2);
        double M3lb = 3*Math.pow(M1,3)*(3*SCV-1+Math.sqrt(2)*Math.pow(1-SCV,1.5));
        double tmp0;
        if(SCV <= 1 && Math.abs(M3 - M3lb) < degentol){
            tmp0 =0;
        }else {
            tmp0 = Math.pow(M3,2.0/9) + ((8*Math.pow(M1,3))/3.0 - 2*M2*M1)*M3 - 3*Math.pow(M1,2)*Math.pow(M2,2) + 2*Math.pow(M2,3);
            if(tmp0 < 0){
                APHS.put(0, aph_fit(M1,M2,M3,2));
                return APHS;
            }
        }

        double tmp1 = 3*Math.sqrt(tmp0);
        double tmp2 = M3 - 3*M1*M2;
        double tmp3 = (6*M2 - 12*Math.pow(M1,2));
        int n;
        if(tmp0==0){
            n=1;
        }else {
            n=2;
        }

        Matrix h1v = new Matrix(n,1,n);
        Matrix h2v = new Matrix(n,1,n);
        Matrix r1v = new Matrix(n,1,n);

        if(n==1){
            h2v.set(0,0,tmp2/tmp3);
            h1v.set(0,0,tmp2/tmp3);
        }else {
            h2v.set(0,0,(tmp2 + tmp1)/tmp3);
            h2v.set(1,0,(tmp2 - tmp1)/tmp3);
            h1v.set(1,0,h2v.get(0));
            h1v.set(0,0,h2v.get(2));
        }

        for(int j=0;j<n;j++){
            double h1 = h1v.get(j);
            double h2 = h2v.get(j);
            r1v.set(j,(M1 - h1)/h2);
        }

        int idx = 0;
        for (int j=0;j<n;j++){
            double h1 = h1v.get(j);
            double h2 = h2v.get(j);
            double r1 = r1v.get(j);
            if(h1 > 0 && h2 > 0 && r1 >= -degentol && r1 <= (1+degentol)){
                r1 = Math.max(Math.min(r1,1),0);
                APHS.put(idx,aph2_assemble(h1,h2,r1));
                idx++;
            }
        }

        if(APHS.isEmpty()){
            APHS.put(0,aph_fit(M1,M2,M3,2));
        }

        return APHS;

    }

    public static Map<Integer,Matrix> aph2_assemble(double l1, double l2, double p1){
        Map<Integer,Matrix> APH = new HashMap<>();
        Matrix D0 = new Matrix(2,2,4);
        D0.set(0,0,-1/l1);
        D0.set(0,1,1/l1*p1);
        D0.set(1,0,0);
        D0.set(1,1,-1/l2);

        Matrix D1 = new Matrix(2,2,4);
        D1.set(0,0,1/l1*(1-p1));
        D1.set(0,1,0);
        D1.set(1,0,1/l2);
        D1.set(1,1,0);

        APH.put(0, D0);
        APH.put(1, D1);

        return APH;
    }


    public static class aph2_fit_return_type{
        public Map<Integer,Matrix> APH;
        public Map<Integer,Map<Integer,Matrix>> APHS;

        public aph2_fit_return_type(){
            APH = new HashMap<>();
            APHS = new HashMap<>();
        }
    }

    public static aph2_fit_return_type aph2_fit(double M1, double M2, double M3){
        aph2_fit_return_type result = new aph2_fit_return_type();
        result.APHS = aph2_fitall(M1,M2,M3);

        if(result.APHS.isEmpty()){
            result.APHS = aph2_fitall(M1,aph2_adjust(M1,M2,M3).get(0),aph2_adjust(M1,M2,M3).get(1));
            if(result.APHS.isEmpty()){
                throw new RuntimeException("Fitting APH(2): feasibility could not be restored");
            }
        }

        result.APH = result.APHS.get(0);

        return result;
    }

    public static Map<Integer,Double> aph2_adjust(double M1, double M2, double M3, String method){
        double tol = 1e-4;
        double M2a=0;
        double M3a=0;
        if(method.equals("simple") ){
            double scva;
            double M1sq = Math.pow(M1,2);
            double scv = ((M2-M1sq)/M1sq);
            if(scv<0.5){
                M2a = 1.5 * Math.pow(M1,2);
                scva = ((M2a-M1sq)/M1sq);
            }else {
                M2a = M2;
                scva = scv;
            }
            if(scva<1){
                double lb = 3*Math.pow(M1,3)*(3*scva-1+Math.sqrt(2)*Math.pow(1-scva,1.5));
                double ub = 6*Math.pow(M1,3)*scva;
                if(M3<lb){
                    M3a = lb;
                }else if(M3>ub){
                    M3a = ub;
                }else {
                    M3a = M3;
                }
            }else {
                double lb = 1.5*Math.pow(M1,3)*Math.pow(1+scva,2);
                if(M3<lb){
                    M3a = lb*(1+tol);
                }else {
                    M3a = M3;
                }
            }
        }else {
            //todo other adjust methods
        }

        Map<Integer,Double> result = new HashMap<>();
        result.put(0,M2a);
        result.put(1,M3a);
        return result;
    }

    public static Map<Integer,Double> aph2_adjust(double M1, double M2, double M3){
        return  aph2_adjust(M1, M2, M3, "simple");
    }

    public static Map<Integer, Matrix> mamap2m_fit_gamma_fb_mmap(Map<Integer,Matrix> mmap){
        double M1 = map_moment(mmap.get(0),mmap.get(1),1);
        double M2 = map_moment(mmap.get(0),mmap.get(1),2);
        double M3 = map_moment(mmap.get(0),mmap.get(1),3);

        Matrix P = mmap_pc(mmap);
        Matrix moments = new Matrix(1,1,1);
        moments.set(0,0,1);
        Matrix F = mmap_forward_moment(mmap,moments);
        Matrix B = mmap_backward_moment(mmap,moments);
        throw new RuntimeException("mamap2m_fit_gamma_fb_mmap has not been implemented");

    }

    /**
     * fits a order-n MMAP with given arrival rates lambda
     * @param lambda 1 by k Matrix describing the arrival rates of each jobclass
     * @param n numbers of MMAP states
     * @return MMAP
     */
    public static  Map<Integer,Matrix> mmap_exponential(Matrix lambda, int n){
        int K = lambda.length();
        Map<Integer,Matrix> MMAP = new HashMap<>(2+K);
        MMAP.put(0,new Matrix(n,n,n^2));
        MMAP.put(1, new Matrix(n,n,n^2));
        for(int k=0;k<K;k++){
            double a = lambda.get(0,k);
            Matrix m = new Matrix(n,n,n^2);
            for(int i=0;i<n;i++){
                m.set(i,n-1-i,a);
            }
            MMAP.put(2+k,m);
            MMAP.put(1,MMAP.get(1).add(1,m));
        }
        return mmap_normalize(MMAP);

    }

    public static Map<Integer,Matrix> mmap_exponential(Matrix lambda){
        return mmap_exponential(lambda,1);
    }

    /**
     * Fixes MMAP feasibility by setting negative values to zero and forcing the other conditions.
     * @param MMAP MMAP to be normalized
     * @return normalized MMAP
     */
    public static Map<Integer, Matrix> mmap_normalize(Map<Integer,Matrix> MMAP){
        if(MMAP.isEmpty()){
            return null;
        }
        int K = MMAP.get(0).numRows;
        int C = MMAP.size()-2;

        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++){
                if(i!=j){
                    MMAP.get(0).set(i,j,Math.max(MMAP.get(0).get(i,j),0));
                }
            }
        }

        MMAP.put(1, new Matrix(MMAP.get(0).numRows,MMAP.get(0).numCols,MMAP.get(0).numRows*MMAP.get(0).numCols));

        for(int c=0;c<C;c++){
            MMAP.get(2+c).removeNegative();
            if(Double.isNaN(MMAP.get(2+c).get(0))){
                MMAP.put(2+c,new Matrix(MMAP.get(2+c).numRows,MMAP.get(2+c).numCols,MMAP.get(2+c).numRows*MMAP.get(2+c).numCols));
            }
            MMAP.put(1,MMAP.get(1).add(1,MMAP.get(2+c)));
        }

        for(int k=0;k<K;k++){
            MMAP.get(0).set(k,k,0);
            MMAP.get(0).set(k,k,-MMAP.get(0).sumRows(k)-MMAP.get(1).sumRows(k));
        }

        return MMAP;
    }

    /**
     * takes a MMAP with K types and a probability
     * matrix with element PROB(k,s) giving the probability that a type-k
     * arrival should be marked as a class-s arrival and returns a new MMAP
     * that has R classes as output types.
     * @param MMAP the original MMAP with K tpyes
     * @param prob a K by S matrix describing the probability of a type-k arrival in original MMAP tp be marked as a
     *             type-S arrival in the new MMAP
     * @return a new MMAP with S tpyes
     */
    public static Map<Integer,Matrix> mmap_mark(Map<Integer,Matrix> MMAP, Matrix prob){
        int K = prob.numRows;
        int R = prob.numCols;
        Map<Integer,Matrix> mmap = new HashMap<>(2+R);
        mmap.put(0, MMAP.get(0).clone());
        mmap.put(1, MMAP.get(1).clone());
        for(int r=0;r<R;r++){
            mmap.put(2+r,new Matrix(MMAP.get(0).length(),MMAP.get(0).length(),MMAP.get(0).length()*MMAP.get(0).length()));
            for(int k=0;k<K;k++){
                Matrix a = MMAP.get(2+k).clone();
                a.scale(prob.get(k,r));
                mmap.put(2+r,mmap.get(2+r).add(1,a));
            }
        }

        return mmap;
    }

    /**
     * Changes the mean inter-arrival time of an MMAP.
     * @param MMAP MMAP to be scaled
     * @param M new mean
     * @return scaled MMAP
     */
    public static Map<Integer,Matrix> mmap_scale(Map<Integer,Matrix> MMAP, Matrix M){
        int C = MMAP.size()-2;
        Map<Integer,Matrix> SCALED = new HashMap<>(2+C);
        if(M.length()==1){
            double MOLD = map_mean(MMAP.get(0),MMAP.get(1));
            double ratio = MOLD/M.get(0);

            Matrix D0 = new Matrix(MMAP.get(0));
            D0.scale(ratio);
            SCALED.put(0,D0);
            Matrix D1 = new Matrix(MMAP.get(1));
            D1.scale(ratio);
            SCALED.put(1,D1);

            for(int c=0;c<C;c++){
                Matrix a = new Matrix(MMAP.get(2+c));
                a.scale(ratio);
                SCALED.put(2+c,a);
            }

        }else {
            SCALED.put(0,MMAP.get(0).clone());
            SCALED.put(1,new Matrix(MMAP.get(0).numRows,MMAP.get(0).numCols,MMAP.get(0).numRows*MMAP.get(0).numCols));
            Matrix l = mmap_count_lambda(MMAP);
            for(int c=0;c<C;c++){
                if(l.get(c)>0){
                    Matrix a = MMAP.get(2+c).clone();
                    a.scale((1/M.get(c))/l.get(c));
                    SCALED.put(2+c,a);
                    SCALED.put(1,SCALED.get(1).add(1,a));
                }else {
                    SCALED.put(2+c,new Matrix(MMAP.get(2+c).numRows,MMAP.get(2+c).numCols,MMAP.get(2+c).numRows*MMAP.get(2+c).numCols));
                }
            }

            SCALED = mmap_normalize(SCALED);
        }
        return SCALED;

        //todo iterative approximation
    }

    /**
     * Computes the arrival rate of the counting process, for the given Marked MAP.
     * @param mmap the Marked MAP
     * @return the vector with the rate for each job class
     */
    public static Matrix mmap_count_lambda(Map<Integer,Matrix> mmap){
        int n = mmap.get(0).numRows;
        int K = mmap.size()-2;

        // symbolic map haven't been implemented

        Matrix lk = new Matrix(K,1,K);

        Matrix theta = map_prob(mmap.get(0),mmap.get(1));

        for(int k=0;k<K;k++){
            lk.set(k,0,theta.mult(mmap.get(2+k)).elementSum());
        }

        return lk.transpose();
    }

    /**
     * same as mmap_count_lambda
     * @param MMAP
     * @return
     */
    public static Matrix mmap_lambda(Map<Integer,Matrix> MMAP){
        return mmap_count_lambda(MMAP);
    }

    /**
     * Combine two MMAPs in to one sup MMAP
     * @param MMAPa MMAP 1
     * @param MMAPb MMAP 2
     * @param opt "default" or "match", "default" treats two MMAPs are irrelevant and each class in MMAPa and MMAPb is a distinct class in SUP
     *            "match" requires MMAPa and MMAPb have the same number of classes and class c in both MMAPa and MMAPb is mapped both into class c of SUP
     * @return
     */
    public static Map<Integer,Matrix> mmap_super(Map<Integer,Matrix> MMAPa, Map<Integer,Matrix> MMAPb, String opt){
        Map<Integer,Matrix> sup = new HashMap<>();
        if(opt.equals("default")){
            int K1 = MMAPa.size()-2;
            int K2 = MMAPb.size()-2;
            int n1 = MMAPa.get(0).length();
            int n2 = MMAPb.get(0).length();

            sup.put(0,MMAPa.get(0).krons(MMAPb.get(0)));
            sup.put(1,MMAPa.get(1).krons(MMAPb.get(1)));

            for(int i=0;i<K1;i++){
                sup.put(2+i,MMAPa.get(2+i).krons(new Matrix(n2,n2,0)));
            }

            for(int j=0;j<K2;j++){
                Matrix a = new Matrix(n1,n1,0);
                sup.put(2+K1+j,a.krons(MMAPb.get(2+j)));
            }

        }else if(opt.equals("match")){
            int K1 = MMAPa.size();
            int K2 = MMAPb.size();

            if(K1!=K2){
                throw new RuntimeException("class matching failed: MMAPs have different number of classes");
            }

            for(int i=0;i<K1;i++){
                sup.put(i,MMAPa.get(i).krons(MMAPa.get(i)));
            }

        }else {
            throw new RuntimeException("unrecognized option");
        }

        return mmap_normalize(sup);
    }

    public static Map<Integer,Matrix> mmap_super(Map<Integer,Matrix> MMAPa, Map<Integer,Matrix> MMAPb){
        return mmap_super(MMAPa,MMAPb,"default");
    }

    public static Map<Integer,Matrix> mmap_super(Map<Integer,Matrix> MMAPa){
        for(Integer key:MMAPa.keySet()){
            if(MMAPa.get(key)==null){
                MMAPa.remove(key);
            }
        }
        Map<Integer,Matrix> SUP = new HashMap<>();
        SUP.put(0,MMAPa.get(0));
        for(int i=1;i< MMAPa.size();i++){
            Map<Integer,Matrix> MMAPb = new HashMap<>();
            MMAPb.put(0,MMAPa.get(i));
            SUP = mmap_super(SUP,MMAPb);
        }
        return SUP;
    }

    /**
     * Makes a subset of the MMAP type hidden. takes a MMAP with K arrival types and
     * % returns a new MMAP for the same stochastic process with those arrivals
     * % hidden from observation.
     * @param MMAP the original MMAP
     * @param types types to be hidden
     * @return
     */
    public static Map<Integer,Matrix> mmap_hide(Map<Integer,Matrix> MMAP, Matrix types){
        Map<Integer,Matrix> mmap = new HashMap<>();
        mmap.put(0,MMAP.get(0));
        mmap.put(1,MMAP.get(1));
        for(int i=0;i<types.length();i++){
            mmap.put((int) (2+types.get(2+i)),new Matrix(MMAP.get(0).numRows,MMAP.get(0).numRows));
        }

        return mmap_normalize(mmap);
    }

    /**
     * Returns K MAPs, one for each class of the MMAP[K] process.
     * @param MMAP MMAP
     * @return K MAPs
     */
    public static Map<Integer,Map<Integer,Matrix>> mmap_maps(Map<Integer,Matrix> MMAP){
        int K = MMAP.size()-2;
        Map<Integer,Map<Integer,Matrix>> Maps = new HashMap<>();
        for(int k=0;k<K;k++){
            Maps.put(k,new HashMap<>());
            Maps.get(k).put(0,MMAP.get(0).add(1,MMAP.get(1)).add(-1,MMAP.get(2+k)));
            Maps.get(k).put(1,MMAP.get(2+k));
        }
        return Maps;
    }
    public static mmap_mixture_fit_return_type mmap_mixture_fit_mmap(Map<Integer,Matrix> mmap){
        Map<Integer[],Matrix> P2 = mmap_sigma2(mmap);
        Matrix M1 = mmap_cross_moment(mmap,1);
        Matrix M2 = mmap_cross_moment(mmap,2);
        Matrix M3 = mmap_cross_moment(mmap,3);

        return mmap_mixture_fit(P2,M1,M2,M3);
    }

    public static class mmap_mixture_fit_return_type{
        Map<Integer,Matrix> MMAP;
        Map<Integer[],Map<Integer,Matrix>> PHs;

        public mmap_mixture_fit_return_type(){
            MMAP = new HashMap<>();
            PHs = new HashMap<>();
        }
    }

    public static mmap_mixture_fit_return_type mmap_mixture_fit(Map<Integer[],Matrix> P2, Matrix M1, Matrix M2, Matrix M3){
        mmap_mixture_fit_return_type result = new mmap_mixture_fit_return_type();
        int m = M1.numRows;
        for(int i=0;i<m;i++){
            for(int j=0;j<m;j++){
                result.PHs.put(new Integer[]{i,j},aph2_fit(M1.get(i,j),M2.get(i,j),M3.get(i,j)).APH);
            }
        }// todo need to correct
        result.MMAP.put(0, new Matrix(2*m*m,2*m*m, 4*(int)Math.pow(m,4)));
        result.MMAP.put(1, new Matrix(2*m*m,2*m*m, 4*(int)Math.pow(m,4)));
        for(int i =0;i<m;i++){
            result.MMAP.put(2+i,  new Matrix(2*m*m,2*m*m, 4*(int)Math.pow(m,4)));
        }
        int  k=1;
        for(int i=0;i<m;i++){
            for(int j=0;j<m;j++){
                //todo block
            }
        }
        return result;
    }

    public static Matrix mmap_cross_moment(Map<Integer,Matrix> mmap,int k){
        int C= mmap.size()-2;

        //todo sym matrix

        Map<Integer,Matrix> TG = new HashMap<>();
        Matrix MC = new Matrix(C,C,(int)Math.pow(C,2));

        for(int i=0;i<C;i++){
            Matrix a = mmap.get(0).inv().mult(mmap.get(2+i));
            a.scale(-1);
            TG.put(i,map_pie(mmap.get(0),mmap.get(1)).mult(a).sumCols());
        }

        for(int i=0;i<C;i++){
            Matrix a = mmap.get(0).inv().mult(mmap.get(2+i));
            a.scale(-1);
            Matrix start = map_pie(mmap.get(0),mmap.get(1)).mult(a).mult(TG.get(i).inv());
            for(int j=0;j<C;j++){
                Matrix neg_mmap0 = mmap.get(0).clone();
                neg_mmap0.scale(-1);
                MC.set(i,j,factorial(k)*start.mult((Matrix.pow(neg_mmap0.inv(),k+1)).mult(mmap.get(2+j))).elementSum());
                MC.set(i,j,MC.get(i,j)/start.mult(neg_mmap0.mult(mmap.get(2+j).inv())).elementSum());
            }
        }
        return MC;
    }

    public static Map<Integer[],Matrix> mmap_sigma2(Map<Integer,Matrix> mmap){
        int C = mmap.size()-2;
        //todo sym matrix
        Map<Integer[],Matrix> sigma = new HashMap<>();

        Matrix alpha = map_pie(mmap.get(0),mmap.get(1));

        for(int i=0;i<C;i++){
            Matrix starti_ = mmap.get(0).inv().mult(mmap.get(2+i));
            starti_.scale(-1);
            Matrix starti = alpha.mult(starti_);
            for(int j=0;j<C;j++){
                Matrix startj_ = mmap.get(0).inv().mult(mmap.get(2+j));
                starti_.scale(-1);
                Matrix startj = starti.mult(startj_);
                for(int h =0;h<C;h++){
                    Matrix starth_ = mmap.get(0).inv().mult(mmap.get(2+h));
                    starti_.scale(-1);
                    Matrix starth = startj.mult(starth_);
                    sigma.put(new Integer[]{i,j,h},starth.sumCols());//// todo need to correct
                }
            }
        }
        return sigma;
    }


    public static Map<Integer,Matrix> mmap_super_safe(Map<Integer,Map<Integer,Matrix>> MMAPS, int maxorder, String method){
        Map<Integer,Matrix> sup = new HashMap<>();
        //todo need to consider empty hashmap;
        List<Double> scv_unmarked = new ArrayList<>();
        for(int i=0;i< MMAPS.size();i++){
            scv_unmarked.add(map_scv(MMAPS.get(i).get(0),MMAPS.get(i).get(1)));
        }

        Integer[] indices = new Integer[scv_unmarked.size()];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }
        Arrays.sort(indices, Comparator.comparingDouble(scv_unmarked::get));
        int[] sortedIndices = Arrays.stream(indices).mapToInt(Integer::intValue).toArray();

        for(int i=0;i<sortedIndices.length;i++){
            int smallest_value = sortedIndices[i];
            if(sup.isEmpty()){
                sup = MMAPS.get(smallest_value);
                if(maxorder==1){
                    sup = mmap_exponential(mmap_lambda(MMAPS.get(smallest_value)));
                }
            }else {
                if(sup.get(0).length()*MMAPS.get(i).get(0).length()>maxorder){
                    if(sup.get(0).length()*2<maxorder){
                        sup = mmap_super(sup,mamap2m_fit_gamma_fb_mmap(MMAPS.get(smallest_value)),method);
                    }else {
                        sup = mmap_super(sup,mmap_exponential(mmap_lambda(MMAPS.get(smallest_value))),method);
                    }
                }else {
                    sup = mmap_super(sup,MMAPS.get(smallest_value),method);
                }
            }
        }


        return sup;
    }

    public static Map<Integer,Matrix> mmap_super_safe(Map<Integer,Map<Integer,Matrix>> MMAPS, int maxorder){
        return mmap_super_safe(MMAPS,maxorder,"default");
    }


    public static Matrix mmap_pc(Map<Integer,Matrix> MMAP){
        int m = MMAP.size()-2;

        //todo sym matrix
        Matrix neg_D0 = MMAP.get(0).clone();
        neg_D0.scale(-1);
        Matrix PC = new Matrix(m,1,m);
        for(int i=0;i<m;i++){
            PC.set(i,map_pie(MMAP.get(0),MMAP.get(1)).mult(neg_D0.inv().mult(MMAP.get(2+i))).elementSum());
        }
        return PC;
    }

    public static Matrix mmap_forward_moment(Map<Integer,Matrix> MMAP, Matrix ORDERS, int NORM){
        int C = MMAP.size()-2;
        int K = ORDERS.length();

        //todo sym matrix

        Matrix MOMENTS = new Matrix(C,K,C*K);
        Matrix pie = map_pie(MMAP.get(0),MMAP.get(1));

        Matrix neg_D0 = MMAP.get(0).clone();
        neg_D0.scale(-1);
        Matrix M = neg_D0.inv();

        for(int a =0;a<C;a++){
            double pa;
            if(NORM==1){
                pa = pie.mult(M.mult(MMAP.get(2+a))).elementSum();
            }else {
                pa = 1;
            }
            for (int h =0;h<ORDERS.length();h++){
                int k = (int) ORDERS.get(h);
                double fk = factorial(k);
                MOMENTS.set(a,h,fk/pa*pie.mult(M.mult(MMAP.get(2+a))).mult(Matrix.pow(M,k)).elementSum());
            }
        }
        return MOMENTS;
    }
    public static Matrix mmap_forward_moment(Map<Integer,Matrix> MMAP, Matrix ORDERS){
        return mmap_forward_moment(MMAP,ORDERS,1);
    }

    public static Matrix mmap_backward_moment(Map<Integer,Matrix> MMAP, Matrix ORDERS, int NORM){
        int C = MMAP.size()-2;
        int K = ORDERS.length();

        //todo sym matrix

        Matrix MOMENTS = new Matrix(C,K,C*K);
        Matrix pie = map_pie(MMAP.get(0),MMAP.get(1));

        Matrix neg_D0 = MMAP.get(0).clone();
        neg_D0.scale(-1);
        Matrix M = neg_D0.inv();

        for(int a =0;a<C;a++){
            double pa;
            if(NORM==1){
                pa = pie.mult(M.mult(MMAP.get(2+a))).elementSum();
            }else {
                pa = 1;
            }
            for (int h =0;h<ORDERS.length();h++){
                int k = (int) ORDERS.get(h);
                double fk = factorial(k);
                MOMENTS.set(a,h,fk/pa*pie.mult(Matrix.pow(M,k+1)).mult(MMAP.get(2+a)).elementSum());
            }
        }
        return MOMENTS;

    }
    public static Matrix mmap_backward_moment(Map<Integer,Matrix> MMAP, Matrix ORDERS){
        return mmap_backward_moment(MMAP,ORDERS,1);
    }

    public static Map<Integer,Matrix> mmap_mixture(Matrix alpha, Map<Integer,Map<Integer,Matrix>> MAPs){
        Map<Integer,Matrix> DK = new HashMap<>();
        int I = MAPs.size();

        for(int i=0;i<I;i++){
            if(MAPs.get(i).isEmpty()){
                MAPs.put(i,map_exponential(1e6));
            }
        }

        for(int i=0;i<I;i++){
            if(i==0){
                DK.put(0, MAPs.get(i).get(0));
            }
        }
        //todo block

        return DK;
    }

}
