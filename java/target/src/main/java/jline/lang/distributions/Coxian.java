package jline.lang.distributions;

import java.util.*;


import jline.lang.constant.GlobalConstants;
import jline.lang.processes.MAP;
import jline.util.Matrix;
import static jline.lib.KPCToolbox.*;

@SuppressWarnings("unchecked")
public class Coxian extends MarkovianDistribution{

	List<Double> totalPhaseRate;

	public Coxian(List<Double> mu, List<Double> phi) {
        super("Coxian", 1);
		
        if (!checkParameter(mu, phi))
        	throw new RuntimeException("Parameter Error");

		int nPhases = mu.size();
        this.setParam(1, "mu", mu);
        this.setParam(2, "phi", phi);

		this.totalPhaseRate = new ArrayList<Double>(nPhases);

		for (int i = 0; i < nPhases; i++) {
			double tpr = mu.get(i);
			this.totalPhaseRate.add(tpr);
		}
	}

	private boolean checkParameter(List<Double> p, List<Double> lambda) {	
		if (p.size() == 2 && lambda.size() == 1)
			return true;
		
		if (p.size() != lambda.size()) {
			return false;
		} else {
			return !((Math.abs(lambda.get(lambda.size() - 1) - 1) > GlobalConstants.CoarseTol) && Double.isFinite(lambda.get(lambda.size()-1)));
		}
	}

	@Override
	public List<Double> sample(long n) {
		return this.sample(n,null);
	}
	@Override
	public List<Double> sample(long n, Random random) {
		throw new RuntimeException("Not implemented");
	}

		@Override
	public double getMean() {
		if(this.getNumberOfPhases() == 2) {
			double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
			double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
			double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
			return 1/mu1 + (1-phi1)/mu2;
		} else {
			Map<Integer, Matrix> rep = this.getRepres();
			return map_mean(rep.get(0), rep.get(1));
		}
	}

	@Override
	public double getRate() {
		return 1/getMean();
	}

	@Override
	public double getSCV() {
		if(this.getNumberOfPhases() == 2) {
			double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
			double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
			double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
			double mean = 1/mu1 + (1-phi1)/mu2;
			double var = ((2*mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (2*mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1))/(mu1*mu1*((mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1))) - (1/mu1 - (phi1 - 1)/mu2)*(1/mu1 - (phi1 - 1)/mu2) - (((phi1 - 1)/(mu2*mu2) + (phi1 - 1)/(mu1*mu2))*((2*mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (2*mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1)))/((mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1));
			return var / Math.pow(mean, 2);
		} else {
			Map<Integer, Matrix> rep = this.getRepres();
			Matrix D0 = rep.get(0);
			Matrix D1 = rep.get(1);
			MAP map = new MAP(D0, D1);
			return map.getSCV();
		}
	}

	@Override
	public double getVar() {
		return this.getSCV()*Math.pow(this.getMean(), 2);
	}

	@Override
	public double getSkew() {
		return super.getSkew();
	}

	@Override
	public double evalCDF(double t) {
		return super.evalCDF(t);
	}

	@Override
	public double evalLST(double s) {
		return super.evalLST(s);
	}

	@Override
	public long getNumberOfPhases() {
		return ((List<Double>)this.getParam(1).getValue()).size();
	}

	@Override
	public Map<Integer, Matrix> getPH() {
    	Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
		if(this.getNumberOfPhases() == 2) {
			double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
			double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
			double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
			Matrix matrix1 = new Matrix(2,2,4);
			Matrix matrix2 = new Matrix(2,2,4);
			matrix1.set(0, 0, -mu1); matrix1.set(0, 1, (1-phi1)*mu1); matrix1.set(1, 1, -mu2);
			matrix2.set(0, 0, phi1*mu1); matrix2.set(1, 0, mu2);
			Map<Integer,Matrix> D = map_normalize(matrix1,matrix2);
			res.put(0, D.get(0));
			res.put(1, D.get(1));
		} else {
			Matrix mu = getMu();
			Matrix phi = getPhi();
			//diag(mu(1:end-1).*(1-phi(1:end-1)),1)
			Matrix expression1 = new Matrix(mu.numRows, mu.numRows, mu.numRows - 1);
			for(int i = 0; i < mu.numRows - 1; i++)
				expression1.set(i, i+1, mu.get(i, 0) * (1 - phi.get(i, 0)));
			
			//diag(-mu)
			Matrix tmp = mu.clone();
			tmp.changeSign();
			Matrix expression2 = Matrix.diag(tmp.nz_values);
			
			//diag(-mu)+diag(mu(1:end-1).*(1-phi(1:end-1)),1)
			Matrix expression3 = expression1.add(1, expression2);
			
			//phi.*mu
			Matrix expression4 = new Matrix(mu.numRows, 1);
			for(int i = 0; i < mu.numRows; i++)
				expression4.set(i, 0, mu.get(i, 0) * (phi.get(i, 0)));

			//zeros(length(mu),length(mu)-1)
			Matrix expression5 = new Matrix(mu.numRows, mu.numRows - 1);
			
			//[phi.*mu,zeros(length(mu),length(mu)-1)]
			Matrix expression6 = new Matrix(0,0,0);
			Matrix.concatColumns(expression4, expression5, expression6);

			res = map_normalize(expression3,expression6);
		}
    	return res;
	}
	
	public Matrix getMu() {
		List<Double> mu = (List<Double>) this.getParam(1).getValue();
		Matrix res = new Matrix(mu.size(), 1, mu.size());
		if (this.getNumberOfPhases() == 2) {
			res.set(0, 0, mu.get(0));
			res.set(1, 0, mu.get(1));
		} else {
			for(int i = 0; i < mu.size(); i++) 
				res.set(i, 0, mu.get(i));
		}
		return res;
	}
	
	public Matrix getPhi() {
		List<Double> phi = (List<Double>) this.getParam(2).getValue();
		Matrix res = new Matrix((int)this.getNumberOfPhases(), 1, (int)this.getNumberOfPhases());
		if (this.getNumberOfPhases() == 2) {
			res.set(0, 0, phi.get(0));
			res.set(1, 0, 1);
		} else {
			for(int i = 0; i < phi.size(); i++) 
				res.set(i, 0, phi.get(i));
		}
		return res;
	}

    // Fit a Coxian distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
    public static Coxian fitMeanAndSCV(double mean, double SCV) {

		double n;
		List<Double> mu = new LinkedList<>();
		List<Double> phi = new LinkedList<>();
		double lambda;

		if (SCV >= 1 - GlobalConstants.CoarseTol && SCV <= 1 + GlobalConstants.CoarseTol) {
			n = 1;
			mu.add(1 / mean);
			phi.add(1.0);
		} else if (SCV > 0.5 + GlobalConstants.CoarseTol && SCV < 1 - GlobalConstants.CoarseTol) {
			phi.add(0.0);
			n = 2;
			mu.add(2 / mean / (1 + Math.sqrt(1 + (2 * (SCV - 1)))));
			mu.add(2 / mean / (1 - Math.sqrt(1 + (2 * (SCV - 1)))));
		} else if (SCV <= 0.5 + GlobalConstants.CoarseTol) {
			n = Math.ceil(1 / SCV);
			lambda = n / mean;
			for (int i = 0; i < n; i++) {
				mu.add(lambda);
				phi.add(0.0);
			}
		} else { // SCV > 1 + GlobalConstants.CoarseTol
			n = 2;
			// transform hyperexp into coxian
			mu.add(2 / mean);
			mu.add((2 / mean) / (2 * SCV));
			phi.add(1 - (mu.get(1) / mu.get(0)));
			phi.add(1.0);
		}

		phi.set((int) n - 1, 1.0);
		Coxian cx = new Coxian(mu, phi);
		cx.immediate = mean < GlobalConstants.CoarseTol;
		return cx;
	}

	public double getTotalPhaseRate(int i) {
		return totalPhaseRate.get(i);
	}

}
