package jline.lang.distributions;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;

public class Erlang extends MarkovianDistribution implements Serializable {
    public Erlang(double phaseRate, long nPhases) {
        super("Erlang", 2);
        this.setParam(1, "alpha", phaseRate);
        this.setParam(2, "r", nPhases);
        double alpha = (double)this.getParam(1).getValue();
        long r = (long) this.getParam(2).getValue();
    }

    @Override
    public List<Double> sample(long n) {
        throw new RuntimeException("Not implemented");
    }

    public List<Double> sample(long n, Random random) {
        return this.sample(n,random);
    }

    public long getNumberOfPhases() {
        return (long) this.getParam(2).getValue();
    }

    public double getMean() {
        double alpha = (double)this.getParam(1).getValue();
        long r = (long) this.getParam(2).getValue();
        return r/alpha;
    }

    public double getVar() {
        double alpha = (double)this.getParam(1).getValue();
        long r = (long) this.getParam(2).getValue();
        return r/Math.pow(alpha,2);
    }

    public double getSkew() {
        long r = (long) this.getParam(2).getValue();
        return 2.0/Math.sqrt(r);
    }

    public double getSCV() {
        long r = (long) this.getParam(2).getValue();
        return 1.0/r;
    }

    public double getRate() {
        return 1.0/getMean();
    }

    public double evalCDF(double t) {
        double alpha = (double)this.getParam(1).getValue();
        long r = (long) this.getParam(2).getValue();
        double ft = 1;

        for (int j = 0; j < r; j++) {
            int fac_j = 1;
            for (int k = 2; k <= j; k++) {
                fac_j *= k;
            }
            ft -= Math.exp(-alpha*t)*(alpha*t)*j/fac_j;
        }

        return ft;
    }

    public Map<Integer, Matrix> getPH()  {
        // This function has to be independent of getD0, getD1, getRepres
        double mu =  (double) this.getParam(1).getValue();
        long r =  (long) this.getParam(2).getValue();
        int size = (int) r;
		Matrix D0 = new Matrix(size, size);
        Matrix D1 = new Matrix(size, size);

		for(int i = 0; i < size - 1; i++) {
            D0.set(i, i, -mu);
			D0.set(i, i+1, mu);
		}
        D0.set(size-1, size-1, -mu);
		D1.set(size - 1, 0, mu);

		Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
		res.put(0, D0);
		res.put(1, D1);
        return res;
    }
    
    public double evalLST(double s) {
        double alpha = (double)this.getParam(1).getValue();
        long r = (long) this.getParam(2).getValue();
        return Math.pow(alpha/(alpha+s), r);
    }

    public static Erlang fitMeanAndSCV(double mean, double SCV) {
        long r = (long)  Math.ceil(1 / SCV);
        double alpha = r/mean;
        return new Erlang(alpha, r);
    }

    public static Erlang fitMeanAndStdDev(double mean, double stdDev) {
        return Erlang.fitMeanAndSCV(mean, (mean/Math.pow(stdDev,2)));
    }

    // Fit distribution with given mean and number of phases
    public static Erlang fitMeanAndOrder(double mean, long numPhases) {
        double SCV = (double) 1 / numPhases;
        long r = (long) Math.ceil(1 / SCV);
        double alpha = r / mean;
        Erlang er = new Erlang(alpha, r);
        er.immediate = mean < GlobalConstants.CoarseTol;
        return er;
    }
}
