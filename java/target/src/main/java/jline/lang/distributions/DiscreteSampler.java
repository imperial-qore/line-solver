package jline.lang.distributions;

import jline.util.Maths;
import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.lang3.NotImplementedException;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Class for discrete distributions specified from the probability mass function
 */
public class DiscreteSampler extends DiscreteDistribution{

    public DiscreteSampler(Matrix p){
        this(p, createX(p.length()));
    }

    /**
     * Constructs a discrete distribution from a finite probability vector p at the points specified in vector x
     * @param p - the probability of an item
     * @param x - the value of an item
     */
    public DiscreteSampler(Matrix p, Matrix x){
        super("DiscreteSampler", 3, new Pair<Double,Double>(x.elementMin(), x.elementMax()));
        setParam(1, "p", p.columnMajorOrder().transpose());
        setParam(2, "x", x.columnMajorOrder().transpose());
        Matrix f = p.columnMajorOrder().transpose().cumsumViaRow();
        double psum = p.elementSum();
        for(int i = 0; i < f.getNumCols(); i++){
            f.set(i, f.get(i) / psum);
        }
        setParam(3, "f", f);
    }

    @Override
    public List<Double> sample(long n) {
        return this.sample(n,new Random());
    }

    public List<Double> sample(long n, Random random) {
        Matrix x = (Matrix) this.getParam(2).getValue();
        Matrix f = (Matrix) this.getParam(3).getValue();
        Matrix r = new Matrix((int) n, 1);
        for(int i = 0; i < n; i++){
            r.set(i, random.nextDouble());
        }
        List<Double> retList = new ArrayList<>();
        Matrix rep = r.repmat(1, f.getNumCols());
        for(int i = 0; i < rep.getNumRows(); i++){
            for(int j = 0; j < rep.getNumCols(); j++){
                rep.set(i, j, Maths.min(rep.get(i, j) <= f.get(j) ? 1 : 0, 2));
            }
        }
        rep = rep.transpose();
        Matrix indexes = new Matrix(1, rep.getNumCols());
        for(int i = 0; i < rep.getNumCols(); i++){
            int max = Integer.MIN_VALUE;
            int max_index = -1;
            for(int j = 0; j < rep.getNumRows(); j++){
                if(indexes.get(j, i) > max){
                    max = (int) indexes.get(j, i);
                    max_index = j;
                }
            }
            indexes.set(i, max_index);
        }
        for(int i = 0; i < indexes.getNumCols(); i++){
            retList.add(x.get((int) indexes.get(i)));
        }
        return retList;
    }

    /**
     * Computes the distribution mean
     * @return - the mean of the distribution
     */
    @Override
    public double getMean() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        int n = p.length();
        return p.mult(createX(n).transpose()).elementSum();
    }

    @Override
    public double getRate() {
        throw new NotImplementedException("getRate() not implemented in DiscreteSampler!");
    }

    /**
     * Computes the distribution squared coefficient of variation (SCV = variance/mean^2)
     * @return
     */
    @Override
    public double getSCV() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        int n = p.length();
        double e2 = p.mult(createX(n).elementMult(createX(n), null).transpose()).elementSum();
        double ex = getMean();
        double var = e2 - ex * ex;
        return var / (ex * ex);
    }

    @Override
    public double getVar() {
        throw new NotImplementedException("getVar() not implemented in DiscreteSampler!");
    }

    @Override
    public double getSkew() {
        throw new NotImplementedException("getSkew() not implemented in DiscreteSampler!");
    }

    @Override
    public double evalCDF(double t) {
        Matrix f = (Matrix) this.getParam(3).getValue();
        if(t >= 0 && t < f.length()){
            return f.get((int) t);
        } else {
            return 0;
        }
    }

    public List<Double> evalPMF(){
        Matrix p = (Matrix) this.getParam(1).getValue();
        List<Double> retList = new ArrayList<>();
        for(int i = 0; i < p.getNumCols(); i++){
            retList.add(p.get(i));
        }
        return retList;
    }

    @Override
    public List<Double> evalPMF(List<Double> t){
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        List<Double> retList = new ArrayList<>();
        for(int i = 0; i < t.size(); i++){
            int j = 0;
            while(j < x.length() && x.get(j) != t.get(i)){
                j++;
            }
            retList.add(p.get(j));
        }
        return retList;
    }

    @Override
    public boolean isDisabled(){
        Matrix p = (Matrix) this.getParam(1).getValue();
        for(int i = 0; i < p.length(); i++){
            if(Double.isNaN(p.get(i)))
                return true;
        }
        return false;
    }

    @Override
    public double evalLST(double s) {
        throw new NotImplementedException("evalLST() not implemented in DiscreteSampler!");
    }

    private static Matrix createX(int length) {
        Matrix x = new Matrix(1, length);
        for(int i = 0; i < x.getNumCols(); i++){
            x.set(i, i+1);
        }
        return x;
    }
}
