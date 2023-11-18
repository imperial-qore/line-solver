package jline.lang.distributions;

import java.util.*;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;
import static jline.lib.KPCToolbox.*;

@SuppressWarnings("unchecked")
public class APH extends MarkovianDistribution {

	List<Double> totalPhaseRate;

	public APH(List<Double> alpha, Matrix T) {
        super("APH", 3);

		int nPhases = alpha.size();

		this.setParam(1, "n", nPhases);
        this.setParam(2, "alpha", alpha);
        this.setParam(3, "T", T);
		this.totalPhaseRate = new ArrayList<Double>(nPhases);

		for (int i = 0; i < nPhases; i++) {
			double tpr = 0.0;
			tpr -= T.get(i,i);
			this.totalPhaseRate.add(tpr);
		}
	}

	public Matrix getInitProb() {
		List<Double> param1 = (List<Double>) this.getParam(2).getValue();
		Matrix alpha = new Matrix(1, param1.size(), param1.size());
		for(int i = 0; i < param1.size(); i++)
			alpha.set(0, i, param1.get(i));
		return alpha;
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
		return super.getMean();
	}

	@Override
	public double getRate() {
		return 1/getMean();
	}

	@Override
	public double getSCV() {
		return super.getSCV();
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
		//Since currently no function to support calculating eigen values, thus calculating expm. This method is now not implemented.
		throw new RuntimeException("Not implemented");
	}

	@Override
	public double evalLST(double s) {
		return super.evalLST(s);
	}

	@Override
	public Map<Integer, Matrix> getPH() {
		Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
		Matrix T = getSubgenerator();

		Matrix ones = new Matrix(T.numCols,1,T.numCols);
		Matrix Te = new Matrix(0,0,0);
		Matrix Tepie = new Matrix(0,0,0);
		ones.fill(1.0);
		T.mult(ones, Te);
		Te.mult(this.getInitProb(), Tepie);
		//Tepie.removeZeros(0);
		Tepie.changeSign();

		res = map_normalize(T,Tepie);
		return res;
	}

	@Override
	public long getNumberOfPhases() {
		return (long) this.getParam(1).getValue();
	}

	public Matrix getSubgenerator() {
		return (Matrix) this.getParam(3).getValue();
	}

	public static APH fitMeanAndSCV(double mean, double scv) {
		APH ex;
		if (mean <= GlobalConstants.FineTol) {
			Matrix T = new Matrix(1,1);
			T.set(0,0,-1/mean);
			List<Double> alpha = Collections.singletonList(1.0);
			ex = new APH(alpha, T);
		} else {
			double e1 = mean;
			double e2 = (1+scv)*e1*e1;
			double cv2 = e2/e1/e1 - 1.0;
			double lambda = 1.0 / e1;
			int N = Math.max((int)Math.ceil(1.0/cv2), 2);
			double p = 1.0 / (cv2 + 1.0 + (cv2-1.0)/(N-1));
			Matrix T = new Matrix(N,N);
			Matrix.eye(N).scale(-lambda*p*N,T);
			for (int i=0;i<N-1;i++) {
				T.set(i,i+1,-1.0*T.get(i,i));
			}
			T.set(N-1,N-1,-lambda*N);
			Matrix alpha = new Matrix(1,N);
			alpha.set(0,0,p);
			alpha.set(0,N-1,1.0 - p);
			ex = new APH(alpha.toList1D(), T);
			ex.immediate = false;
		}
		return ex;
	}

	public double getTotalPhaseRate(int i) {
		return totalPhaseRate.get(i);
	}
}
