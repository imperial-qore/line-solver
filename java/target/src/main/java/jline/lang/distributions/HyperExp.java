package jline.lang.distributions;

import jline.util.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;


@SuppressWarnings("unchecked")
public class HyperExp extends MarkovianDistribution  implements Serializable {
	
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
	 * @param n - the number of samples
	 * @return - n samples from the distribution
	 */
	@Override
	public List<Double> sample(long n) {
		return this.sample(n,new Random());
	}

	@Override
	public List<Double> sample(long n, Random random) {
		throw new RuntimeException("Not implemented");
	}

	public long getNumberOfPhases() {
        return nPhases;
    }

    public double evalCDF(double t) {
    	if (this.nPhases == 2) {
    		double p = (double) this.getParam(1).getValue();
    		double mu1 = (double) this.getParam(2).getValue();
    		double mu2 = (double) this.getParam(3).getValue();
    		return p*(1-Math.exp(-mu1*t)) + (1-p)*(1-Math.exp(-mu2*t));
    	} else {
    		return super.evalCDF(t);
    	}
    }

	public Map<Integer, Matrix> getPH() {
    	Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
    	
    	double p = (double) this.getParam(1).getValue();
    	double mu1 = (double) this.getParam(2).getValue();
    	double mu2 = (double) this.getParam(3).getValue();
		Matrix D0 = new Matrix(2,2,4);
		Matrix D1 = new Matrix(2,2,4);
		D0.set(0, 0, -mu1); D0.set(1, 1, -mu2);
		D1.set(0, 0, mu1*p); D1.set(0, 1, mu1*(1-p)); D1.set(1, 0, mu2*p); D1.set(1, 1, mu2*(1-p));
		res.put(0, D0);
		res.put(1, D1);
    	
    	return res;
//	NOT USED (FOR N-PHASES)
//    		JLineMatrix D0 = new JLineMatrix(nPhases, nPhases, nPhases);
//    		for(int i = 0; i < nPhases; i++) 
//    			D0.set(i, i, -lambda.get(i));
//    		
//    		JLineMatrix D1 = new JLineMatrix(nPhases, nPhases);
//    		JLineMatrix temP = new JLineMatrix(nPhases, 1, nPhases);
//    		JLineMatrix ones = new JLineMatrix(1, nPhases, nPhases);
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
    		return (2*(p/Math.pow(mu1, 2) + (1-p)/Math.pow(mu2, 2)) - Math.pow(p/mu1 + (1-p)/mu2, 2)) / Math.pow(p/mu1 + (1-p)/mu2, 2);
    	} else {
    		return super.getSCV();
    	}
    }

    public double getRate() {
    	return 1/getMean();
    }

    public double getMean() {
    	if (this.nPhases == 2) {
    		double p = (double) this.getParam(1).getValue();
    		double mu1 = (double) this.getParam(2).getValue();
    		double mu2 = (double) this.getParam(3).getValue();
    		return p/mu1 + (1-p)/mu2;
    	} else {
    		return super.getMean();
    	}
    }

    public double getVar() {
    	return this.getSCV()*Math.pow(this.getMean(), 2);
    }

    public double getSkew() {
        return super.getSkew();
    }

    public String toString() {
        return String.format("jline.HyperExp(%f)", this.getRate());
    }

    public double getRateFromPhase(int phase) {
    	if (phase > this.nPhases)
    		throw new RuntimeException("Exceed the number of phases");

        return ((List<Double>)this.getParam(2)).get(phase - 1);
    }
}
