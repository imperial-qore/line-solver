package jline.lang.distributions;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static jline.lib.KPCToolbox.map_hyperexp;


@SuppressWarnings("unchecked")
public class HyperExp extends MarkovianDistribution implements Serializable {

    private final long nPhases;

    public HyperExp(double p, double lambda1, double lambda2) {
        super("HyperExp", 1);

        this.setParam(1, "p", p);
        this.setParam(2, "lambda1", lambda1);
        this.setParam(3, "lambda2", lambda2);

        nPhases = 2;
    }

    public HyperExp(double p, double lambda) {
        this(p, lambda, lambda);
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public Matrix sample(long n) {
        return this.sample(n, new Random());
    }

    @Override
    public Matrix sample(long n, Random random) {
        throw new RuntimeException("Not implemented"); // TODO: not implemented
    }

    public long getNumberOfPhases() {
        return nPhases;
    }

    public double evalCDF(double t) {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return p * (1 - Math.exp(-mu1 * t)) + (1 - p) * (1 - Math.exp(-mu2 * t));
        } else {
            return super.evalCDF(t);
        }
    }

    public Map<Integer, Matrix> getPH() {
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();

        double p = (double) this.getParam(1).getValue();
        double mu1 = (double) this.getParam(2).getValue();
        double mu2 = (double) this.getParam(3).getValue();
        Matrix D0 = new Matrix(2, 2, 4);
        Matrix D1 = new Matrix(2, 2, 4);
        D0.set(0, 0, -mu1);
        D0.set(1, 1, -mu2);
        D1.set(0, 0, mu1 * p);
        D1.set(0, 1, mu1 * (1.0 - p));
        D1.set(1, 0, mu2 * p);
        D1.set(1, 1, mu2 * (1.0 - p));
        res.put(0, D0);
        res.put(1, D1);

        return res;
//	NOT USED (FOR N-PHASES)
//    		Matrix D0 = new Matrix(nPhases, nPhases, nPhases);
//    		for(int i = 0; i < nPhases; i++) 
//    			D0.set(i, i, -lambda.get(i));
//    		
//    		Matrix D1 = new Matrix(nPhases, nPhases);
//    		Matrix temP = new Matrix(nPhases, 1, nPhases);
//    		Matrix ones = new Matrix(1, nPhases, nPhases);
//    		CommonOps_DSCC.fill(ones, 1.0);
//    		for(int i = 0; i < nPhases; i++)
//    			temP.set(i, 0, p.get(i));
//    		CommonOps_DSCC.mult(CommonOps_DSCC.mult(D0, temP, null), ones, D1);
//    		CommonOps_DSCC.changeSign(D1, D1);
//    		res.put(0, D0);
//    		res.put(1, D1);
    }

    public double evalLST(double s) {
        return super.evalLST(s);
    }

    public double getSCV() {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return (2 * (p / Math.pow(mu1, 2) + (1 - p) / Math.pow(mu2, 2)) - Math.pow(p / mu1 + (1 - p) / mu2, 2)) / Math.pow(p / mu1 + (1 - p) / mu2, 2);
        } else {
            return super.getSCV();
        }
    }

    public double getRate() {
        return 1 / getMean();
    }

    public double getMean() {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return p / mu1 + (1 - p) / mu2;
        } else {
            return super.getMean();
        }
    }

    public double getVar() {
        return this.getSCV() * Math.pow(this.getMean(), 2);
    }

    public double getSkew() {
        return super.getSkew();
    }

    public String toString() {
        return String.format("jline.HyperExp(%f)", this.getRate());
    }

    public double getRateFromPhase(int phase) {
        if (phase > this.nPhases)
            throw new RuntimeException("Exceed the number of phases"); // TODO: not implemented

        return ((List<Double>) this.getParam(2)).get(phase - 1);
    }

    /**
     * Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
     */
    public static HyperExp fitMeanAndSCV(double mean, double scv) {
        double p, mu1, mu2;
        Map<Integer, Matrix> D = map_hyperexp(mean, scv, 0);
        mu1 = -D.get(0).get(0, 0);
        mu2 = -D.get(0).get(1, 1);
        p = D.get(1).get(0, 0) / mu1;
        HyperExp he = new HyperExp(p, mu1, mu2);
        he.immediate = mean < GlobalConstants.CoarseTol;
        return he;
    }

    /**
     * Fit distribution with given squared coefficient of variation and balanced means i.e.,
     * p/mu1 = (1-p)/mu2
     */
    public static HyperExp fitMeanAndSCVBalanced(double mean, double scv) {
        double p, mu1, mu2;
        mu1 = -(2.0 * (Math.sqrt((scv - 1) / (scv + 1)) / 2.0 - 0.5)) / mean;
        p = 0.5 - Math.sqrt((scv - 1) / (scv + 1)) / 2.0;
        if (mu1 < 0 || p < 0 || p > 1) {
            p = Math.sqrt((scv - 1) / (scv + 1)) / 2.0 + 0.5;
            mu1 = (2 * (Math.sqrt((scv - 1) / (scv + 1)) / 2.0 + 0.5)) / mean;
        }
        mu2 = (1 - p) / p * mu1;
        HyperExp he = new HyperExp(p, mu1, mu2);
        he.immediate = mean < GlobalConstants.CoarseTol;
        return he;
    }
}
