package jline.lib;

import jline.api.CTMC;
import jline.util.Matrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class KPCToolbox {

    /**
     * Sanitizes the (D0,D1) matrices of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return MAP
     */
    public static Map<Integer, Matrix> map_normalize(Matrix D0, Matrix D1) {
        Map<Integer,Matrix> D = new HashMap<>();
        // This modifies the inputs
        D.put(0,D0);
        D.put(1,D1);
        long nPhases = D0.length();
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (D.get(0).get(i,j) < 0) {
                    D.get(0).set(i,j,0.0);
                }
                if (D.get(1).get(i,j) < 0) {
                    D.get(1).set(i,j,0.0);
                }
            }
        }

        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < nPhases; j++){
                if (j != i) {
                    rowSum = rowSum + D.get(0).get(i, j);
                }
                rowSum = rowSum + D.get(1).get(i,j);
                D.get(0).set(i,i,-rowSum);
            }
        }
        return D;
    }

    /**
     * CTMC underlying a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return CTMC infinitesimal generator Q
     */
    public static Matrix map_infgen(Matrix D0, Matrix D1) {
        return D0.add(1.0, D1);
    }

    /**
     * Steady-state vector of the CTMC underlying a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return CTMC steady state vector
     */
    public static Matrix map_prob(Matrix D0, Matrix D1) {
        return CTMC.ctmc_solve(map_infgen(D0,D1));
    }

    /**
     * Mean inter-arrival time of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return Mean inter-arrival time
     */
    public static double map_mean(Matrix D0, Matrix D1) {
        return 1/map_lambda(D0,D1);
    }

    /**
     * Arrival rate of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return Arrival rate
     */
    public static double map_lambda(Matrix D0, Matrix D1) {
        Matrix e = Matrix.ones(D0.getNumRows(),1);
        Matrix lambda = map_prob(D0,D1); // piq
        lambda=lambda.mult(D1);
        lambda=lambda.mult(e);
        return lambda.toDouble();
    }

    /**
     * Steady-state probability vector of the embedded DTMC of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return Embedded steady-state probability vector
     */
    public static Matrix map_pie(Matrix D0, Matrix D1) {
        Matrix e = Matrix.ones(D0.getNumRows(),1);
        Matrix A = map_prob(D0,D1).mult(D1); // piq*D1
        A.scale(1 / A.mult(e).toDouble()); // divide by A*ones(length(MAP{2}),1)
        return A;
    }

    /**
     * Raw moments of the inter-arrival times of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @param order Moment order, i.e., E[X^order]
     * @return Raw moment
     */
    public static double map_moment(Matrix D0, Matrix D1, int order) {
        Matrix pie = map_pie(D0,D1);
        Matrix iD0 = D0.clone(); iD0.scale(-1.0);
        iD0 = iD0.inv();
        Matrix iD0k = new Matrix(iD0);
        for (int i = 2; i<= order; i++) {
            iD0k = iD0k.mult(iD0); // scale(i) introduces the factorial
            iD0k.scale(i);
        }
        Matrix e = Matrix.ones(D0.getNumRows(),1);
        pie = pie.mult(iD0k);
        pie = pie.mult(e);
        return pie.toDouble();
    }

    /**
     * Cumulative distribution function (cdf) of the inter-arrival times of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @param points Points for computing the cdf
     * @return Cdf values
     */
    public static Matrix map_cdf(Matrix D0, Matrix D1, Matrix points) {
        Matrix CDFVals = new Matrix(1,points.length());
        Matrix pie = map_pie(D0, D1);
        Matrix e1 = Matrix.ones(D0.numRows,1);

        double nanVal = 0.0;
        for (int t = 0; t < points.length(); t++) {
            Matrix output = CTMC.ctmc_uniformization(pie,D0,points.get(t)).mult(e1);
            double val = 1 - output.get(0, 0);
            if (Double.isNaN(val)) {
                val = nanVal;
            } else { // after it finds the first non-zero, set nanVal to 1.0
                nanVal = 1.0;
            }
            CDFVals.set(0, t, val);
        }
        return CDFVals;
    }

    /**
     * Squared coefficient of variation (scv) of the inter-arrival times of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return Scv value
     */
    public static double map_scv(Matrix D0, Matrix D1) {
        double e1 = map_moment(D0,D1,1);
        double e2 = map_moment(D0,D1,2);

        double var = e2-e1*e1;
        double scv = var/e1/e1;
        return scv;
    }

    /**
     * Rescale mean inter-arrival time of a MAP
     *
     * @param D0 Hidden transition matrix
     * @param D1 Visible transition matrix
     * @return Scaled MAP
     */
    public static Map<Integer,Matrix> map_scale(Matrix D0, Matrix D1, double newMean) {
        Map<Integer,Matrix> D = new HashMap<>();
        D.put(0,D0.clone());
        D.put(1,D1.clone());
        map_normalize(D.get(0),D.get(1));
        double ratio = map_mean(D0,D1) / newMean;
        D.get(0).scale(ratio);
        D.get(1).scale(ratio);
        return D;
    }

    public static Map<Integer,Matrix> map_exponential(double mean){
        double mu = 1/mean;
        Map<Integer,Matrix> MAP = new HashMap<>();
        Matrix D0 = new Matrix(1,1,1);
        D0.set(0,0,-mu);
        Matrix D1 = new Matrix(1,1,1);
        D1.set(0,0,mu);
        MAP.put(0,D0);
        MAP.put(1,D1);
        return MAP;
    }
    public static double map_idc(Matrix D0, Matrix D1){
        Matrix e = new Matrix(D0.length(),1, D0.length());
        for (int i=0;i< D0.length();i++){
            e.set(i,0,1);
        }

        return  1+2*(map_lambda(D0,D1)-map_pie(D0,D1).mult(map_infgen(D0,D1).add(1,e.mult(map_prob(D0,D1))).inv()).mult(D1).mult(e).get(0));

    }

    public static class Map_fit_return_type{
        public Map<Integer,Matrix> MAP;
        public double error;
        public Map_fit_return_type(){
            MAP = new HashMap<>();
            error =0;
        }

    }
    public static Map_fit_return_type map2_fit(double e1, double e2, double e3, double g2){
        Map_fit_return_type result = new Map_fit_return_type();
        double r1 = e1;
        double r2 = e2/2;
        double h2 = (r2-Math.pow(r1,2))/Math.pow(r1,2);
        if(e3==-1){
            double scv = (e2-Math.pow(e1,2))/Math.pow(e1,2);
            if(1<=scv && scv<3){
                if(g2<0){
                    double h3 = h2-Math.pow(h2,2);
                    e3 = 12*Math.pow(e1,3)*h2+6*Math.pow(e1,3)*h3+6*Math.pow(e1,3)*(1+Math.pow(h2,2));
                }else {
                    e3 = 1.501*Math.pow(e2,2)/e1;
                }
            }else if(3<=scv){
                e3 = 1.501*Math.pow(e2,2)/e1;
            }else if(0<scv && scv<1){
                e3 = (1 + 1e-10) * (12 * Math.pow(e1, 3) * h2 + 6 * Math.pow(e1, 3) * (h2 * (1 - h2 - 2 * Math.sqrt(-h2))) + 6 * Math.pow(e1, 3) * (1 + Math.pow(h2, 2)));
            }
        }

        if (e3 == -2) {
            double scv = (e2 - Math.pow(e1, 2)) / Math.pow(e1, 2);

            if (scv >= 1) {
                e3 = (3.0 / 2 + 1e-6) * Math.pow(e2, 2) / e1;
            } else if (0 < scv && scv < 1) {
                double h3 = h2 * (1 - h2 - 2 * Math.sqrt(-h2));
                e3 = 6 * Math.pow(e1, 3) * (Math.pow(h2, 2) + h3);
            }
        }

        if (e3 == -3) {
            double scv = (e2 - Math.pow(e1, 2)) / Math.pow(e1, 2);

            if (scv >= 1) {
                e3 = Math.pow(10, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = Math.pow(-h2, 2);
                e3 = 6 * Math.pow(e1, 3) * (Math.pow(h2, 2) + h3);
            }
        }

        if (e3 == -4) {
            double scv = (e2 - Math.pow(e1, 2)) / Math.pow(e1, 2);
            double r = new Random().nextDouble();

            if (scv >= 1) {
                e3 = r * (3.0 / 2 + 1e-6) * Math.pow(e2, 2) / e1 + (1 - r) * Math.pow(10, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = r * Math.pow(-h2, 2) + (1 - r) * h2 * (1 - h2 - 2 * Math.sqrt(-h2));
                e3 = 6 * Math.pow(e1, 3) * (Math.pow(h2, 2) + h3);
            }
        }

        if (e3 > -1 && e3 < 0) {
            double scv = (e2 - Math.pow(e1, 2)) / Math.pow(e1, 2);
            double r = Math.abs(e3);

            if (scv >= 1) {
                e3 = r * (3.0 / 2 + 1e-6) * Math.pow(e2, 2) / e1 + (1 - r) * Math.pow(10, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = r * h2 * (1 - h2 - 2 * Math.sqrt(-h2)) + (1 - r) * Math.pow(-h2, 2);
                e3 = 6 * Math.pow(e1, 3) * (Math.pow(h2, 2) + h3);
            }
        }

        double r3 = e3 / 6;
        double h3 = (r3 * r1 - Math.pow(r2, 2)) / Math.pow(r1, 4);
        double b = h3 + Math.pow(r1, 2) - r1;
        double c = Math.sqrt(Math.pow(b, 2) + 4 * Math.pow(r1, 3));

        if(r1<=0){
            result.error = 10;
            return result;
        }

        if(h2==0){
            if(h3==0 && g2==0){
                result.MAP = map_exponential(e1);
            }else {
                result.error = 20;
                return result;
            }
        }

        if(h2>0 && h3>0){
            if(b>=0){
                if((b-c)/(b+c)<=g2){
                    Map<Integer,Matrix> MAP = new HashMap<>();
                    Matrix D0 = new Matrix(2,2,4);
                    D0.set(0,0,-(2*h2+b-c));
                    D0.set(0,1,0);
                    D0.set(1,0,0);
                    D0.set(1,1,-(2*h2+b+c));
                    D0.scale(1/(2*r1*h3));
                    MAP.put(0,D0);

                    Matrix D1 = new Matrix(2,2,4);
                    D1.set(0,0,(2*h2+b-c)*(1-b/c+g2*(1+b/c)));
                    D1.set(0,1,(2*h2+b-c)*(1+b/c)*(1-g2));
                    D1.set(1,0,(2*h2+b+c)*(1-b/c)*(1-g2));
                    D1.set(1,1,(2*h2+b+c)*(1+b/c+g2*(1-b/c)));
                    D1.scale(1/(4*r1*h3));
                    MAP.put(1,D1);

                    result.MAP = MAP;
                }else {
                    result.error = 51;
                }
            }else if(b<0){
                if(0<=g2 && g2<1){
                    Map<Integer,Matrix> MAP = new HashMap<>();
                    Matrix D0 = new Matrix(2,2,4);
                    D0.set(0,0,-(2*h2+b-c));
                    D0.set(0,1,0);
                    D0.set(1,0,0);
                    D0.set(1,1,-(2*h2+b+c));
                    D0.scale(1/(2*r1*h3));
                    MAP.put(0,D0);

                    Matrix D1 = new Matrix(2,2,4);
                    D1.set(0,0,(2*h2+b-c)*(1-b/c+g2*(1+b/c)));
                    D1.set(0,1,(2*h2+b-c)*(1+b/c)*(1-g2));
                    D1.set(1,0,(2*h2+b+c)*(1-b/c)*(1-g2));
                    D1.set(1,1,(2*h2+b+c)*(1+b/c+g2*(1-b/c)));
                    D1.scale(1/(4*r1*h3));
                    MAP.put(1,D1);

                    result.MAP = MAP;
                }else if(-(h3+Math.pow(h2,2))/h2 <= g2 && g2<0){
                    double a=(h3+Math.pow(h2,2))/h2;
                    double d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    double d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    Map<Integer,Matrix> MAP = new HashMap<>();
                    Matrix D0 = new Matrix(2,2,4);
                    D0.set(0,0,-(2*h2+b-c));
                    D0.set(0,1,(2*h2+b-c)*(1-a));
                    D0.set(1,0,0);
                    D0.set(1,1,-(2*h2+b+c));
                    D0.scale(1/(2*r1*h3));
                    MAP.put(0,D0);

                    Matrix D1 = new Matrix(2,2,4);
                    D1.set(0,0,(2*h2+b-c)*d1);
                    D1.set(0,1,(2*h2+b-c)*(a-d1));
                    D1.set(1,0,(2*h2+b+c)*d2);
                    D1.set(1,1,(2*h2+b+c)*(1-d2));
                    D1.scale(1/(2*r1*h3));
                    MAP.put(1,D1);

                    result.MAP = MAP;
                }else {
                    result.error = 52;
                }
            }
            if(result.MAP.size()>0&&!map_isfeasible(result.MAP)){
                result.error=-1;
            }
            return result;
        }else if(-0.25<=h2 && h2<0 && h2*(1-h2-2*Math.sqrt(-h2))<=h3 && h3<=-Math.pow(h2,2)){
            if(g2>=0){
                if(g2<=-Math.pow(h2+Math.sqrt(-h3),2)/h2){
                    double a=(2*h2+b-c)*(h2+Math.sqrt(-h3))/(2*h2*Math.sqrt(-h3));
                    c = -c;
                    double d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    double d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    Map<Integer,Matrix> MAP = new HashMap<>();
                    Matrix D0 = new Matrix(2,2,4);
                    D0.set(0,0,-(2*h2+b-c));
                    D0.set(0,1,(2*h2+b-c)*(1-a));
                    D0.set(1,0,0);
                    D0.set(1,1,-(2*h2+b+c));
                    D0.scale(1/(2*r1*h3));
                    MAP.put(0,D0);

                    Matrix D1 = new Matrix(2,2,4);
                    D1.set(0,0,(2*h2+b-c)*d1);
                    D1.set(0,1,(2*h2+b-c)*(a-d1));
                    D1.set(1,0,(2*h2+b+c)*d2);
                    D1.set(1,1,(2*h2+b+c)*(1-d2));
                    D1.scale(1/(2*r1*h3));
                    MAP.put(1,D1);

                    result.MAP = MAP;
                }else {
                    result.error=53;
                }
            }else if(g2<0){
                if(g2 >= -(h3+Math.pow(h2,2))/h2){
                    double a=(h3+Math.pow(h2,2))/h2;
                    c = -c;
                    double d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    double d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                    Map<Integer,Matrix> MAP = new HashMap<>();
                    Matrix D0 = new Matrix(2,2,4);
                    D0.set(0,0,-(2*h2+b-c));
                    D0.set(0,1,(2*h2+b-c)*(1-a));
                    D0.set(1,0,0);
                    D0.set(1,1,-(2*h2+b+c));
                    D0.scale(1/(2*r1*h3));
                    MAP.put(0,D0);

                    Matrix D1 = new Matrix(2,2,4);
                    D1.set(0,0,(2*h2+b-c)*d1);
                    D1.set(0,1,(2*h2+b-c)*(a-d1));
                    D1.set(1,0,(2*h2+b+c)*d2);
                    D1.set(1,1,(2*h2+b+c)*(1-d2));
                    D1.scale(1/(2*r1*h3));
                    MAP.put(1,D1);

                    result.MAP = MAP;
                }else {
                    result.error = 54;
                }
            }

            if(result.MAP.size()>0&&!map_isfeasible(result.MAP)){
                result.error=-1;
            }
            return result;
        }else {
            if(!(-0.25<=h2 && h2<0 && h2*(1-h2-2*Math.sqrt(-h2))<=h3 && h3<=-Math.pow(h2,2))){
                result.error = 30;
            }else if((h2>0 && h3<0)|| h2*(1-h2-2*Math.sqrt(-h2))>h3 || h3<=Math.pow(h2,2)){
                result.error = 40;
            }else {
                throw new RuntimeException("I lost an error");
            }
        }

        return result;
    }

    public static Map_fit_return_type map2_fit(double e1, double e2,double e3){
        return map2_fit(e1,e2,-1,e3);
    }

    public static boolean map_isfeasible(Map<Integer,Matrix> MAP, double TOL){
        return map_checkfeasible(MAP,TOL);
    }

    public static boolean map_isfeasible(Map<Integer,Matrix> MAP){

        int TOLMAGNITUDE = 15;
        for (int k = TOLMAGNITUDE; k >= 1; k--) {
            boolean check = map_checkfeasible(MAP, Math.pow(10, -k));
            if (check) {
                return k > map_feastol();
            }
        }
        return false;
    }

    public static boolean map_checkfeasible(Map<Integer,Matrix> MAP, double TOL){
        int n = MAP.get(0).length();
        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);
        if(D0.hasNaN()){
            return false;
        }
        boolean D0_has_inf = false;
        boolean D1_has_inf = false;

        for(int i=0;i<D0.numRows;i++){
            for (int j=0;j<D0.numCols;j++){
                if(Double.isInfinite(D0.get(i,j))){
                    D0_has_inf = true;
                    break;
                }
            }
        }

        for(int i=0;i<D1.numRows;i++){
            for (int j=0;j<D1.numCols;j++){
                if(Double.isInfinite(D1.get(i,j))){
                    D1_has_inf = true;
                    break;
                }
            }
        }

        if(D0.hasNaN()|| D1.hasNaN()||D0_has_inf||D1_has_inf){
            return false;
        }

        Matrix neg_D0 = D0.clone();
        neg_D0.scale(-1);

        Matrix P = neg_D0.inv().mult(D1);
        Matrix Q = D0.add(1,D1);

        for(int i=0;i<D0.numRows;i++){
            for (int j=0;j<D0.numCols;j++){
                if(Math.abs(D0.get(i,j))<TOL){
                    D0.set(i,j,0);
                }
            }
        }

        for(int i=0;i<D1.numRows;i++){
            for (int j=0;j<D1.numCols;j++){
                if(Math.abs(D1.get(i,j))<TOL){
                    D1.set(i,j,0);
                }
            }
        }

        for(int i=0;i<P.numRows;i++){
            for (int j=0;j<P.numCols;j++){
                if(Math.abs(P.get(i,j))<TOL){
                    P.set(i,j,0);
                }
            }
        }

        for(int i=0;i<Q.numRows;i++){
            for (int j=0;j<Q.numCols;j++){
                if(Math.abs(Q.get(i,j))<TOL){
                    Q.set(i,j,0);
                }
            }
        }


        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i!=j&&D0.get(i,j)<0){
                    return false;
                }
                if(i==j&&D0.get(i,j)>0){
                    return false;
                }
                if(D1.get(i,j)<0){
                    return false;
                }
                if(i!=j&&Q.get(i,j)<0){
                    return false;
                }
                if(i==j&&Q.get(i,j)>0){
                    return false;
                }
                if(P.get(i,j)<0){
                    return false;
                }
            }

            if(Math.abs(P.sumRows(i))<1-n*TOL){
                return false;
            }

            if(Math.abs(P.sumRows(i))>1+n*TOL){
                return false;
            }
            if(Math.abs(Q.sumRows(i))<0-n*TOL){
                return false;
            }
            if(Math.abs(Q.sumRows(i))<0+n*TOL){
                return false;
            }

        }

        if(n<map_largemap()){
            double [][] P_data = P.toArray2D();
            RealMatrix P_matrix = MatrixUtils.createRealMatrix(P_data);
            EigenDecomposition P_eigenDecomposition = new EigenDecomposition(P_matrix);
            double[] P_eigenvalues = P_eigenDecomposition.getRealEigenvalues();
            int P_sum = 0;
            for(double e:P_eigenvalues){
                if (e>1-TOL){
                    P_sum++;
                }
                if(P_sum>1){
                    return false;
                }
            }

            double [][] Q_data = Q.toArray2D();
            RealMatrix Q_matrix = MatrixUtils.createRealMatrix(Q_data);
            EigenDecomposition Q_eigenDecomposition = new EigenDecomposition(Q_matrix);
            double[] Q_eigenvalues = Q_eigenDecomposition.getRealEigenvalues();
            int Q_sum = 0;
            for(double e:Q_eigenvalues){
                if (e>0-TOL){
                    Q_sum++;
                }
                if(Q_sum>1){
                    return false;
                }
            }

        }

        return true;

    }

    public static int map_largemap(){
        return 100;
    }

    public static int map_feastol(){
        return 8;
    }
    public static Matrix map_embedded(Matrix D0, Matrix D1){
        Matrix neg_D0 = D0.clone();
        neg_D0.scale(-1);
        return neg_D0.inv().mult(D1);
    }

    public static Matrix map_acf(Matrix D0, Matrix D1,Matrix LAGS){
        Matrix P = map_embedded(D0,D1);
        Matrix x = map_prob(D0,D1);
        x.scale(map_lambda(D0,D1));

        //todo sym matrix

        Matrix neg_D0 = D0.clone();
        neg_D0.scale(-1);
        Matrix y = neg_D0.inv().sumRows();
        Matrix ACFCOEFFS = new Matrix(1,LAGS.length(), LAGS.length());
        for(int i=0;i< LAGS.length();i++){
            ACFCOEFFS.set(i,x.mult(Matrix.pow(P,(int)LAGS.get(i))).mult(y).get(0));
        }
        for(int i=0;i<ACFCOEFFS.length();i++){
            ACFCOEFFS.set(i,ACFCOEFFS.get(i)-1);
        }

        ACFCOEFFS.scale(1/map_scv(D0,D1));

        return ACFCOEFFS;
    }

    public static Matrix map_acf(Matrix D0, Matrix D1){
        Matrix LAGS = new Matrix(1,1,1);
        LAGS.set(0,0,1);
        return map_acf(D0,D1,LAGS);

    }

    public static double map_skew(Matrix D0, Matrix D1){
        Map<Integer,Double> m = new HashMap<>();
        for(int i=1;i<=3;i++){
            m.put(i,map_moment(D0,D1,i));
        }
        double M3 = m.get(3)-3*m.get(2)*m.get(1)+2*Math.pow(m.get(1),3);
        return M3/Math.pow(Math.sqrt(map_scv(D0,D1))*m.get(1),3);
    }

    public static double map_var(Matrix D0, Matrix D1){
        return map_moment(D0,D1,2)-Math.pow(map_mean(D0,D1),2);
    }

    public static Map<Integer, Matrix> mmpp2_fit1(double mean, double scv, double skew, double idc){
        double E1 = mean;
        double E2 = (1+scv)*Math.pow(E1,2);
        double g2=-(scv-idc)/(-1+idc);
        double E3;
        if(skew==-1){
            E3 = -1;
        }else {
            E3=-(2*Math.pow(E1,3)-3*E1*E2-skew*Math.pow(E2-Math.pow(E1,2),1.5));
        }

        return map2_fit(E1,E2,E3,g2).MAP;
    }

}
