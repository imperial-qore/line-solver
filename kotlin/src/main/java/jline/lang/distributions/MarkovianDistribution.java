package jline.lang.distributions;

import java.io.Serializable;
import java.util.List;
import java.util.Map;


import jline.lang.processes.MAP;
import jline.util.Matrix;

import jline.util.Pair;

import static jline.lib.KPCToolbox.*;

abstract public class MarkovianDistribution extends Distribution implements Serializable {
	
	protected Map<Integer, Matrix> representation; // <0, D0>, <1, D1>, <2, D2>
	public Matrix initProb;
	public Matrix invSubgenerator;

    public MarkovianDistribution(String name, int numParam) {
        super(name, numParam, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
    }

	public long getNumberOfPhases() {
		Map<Integer, Matrix> PH = this.getRepres();
		return PH.get(1).numCols;
	}

	public Matrix getD0() {
		return getRepres().get(0);
	}

	public Matrix getD1() {
		return getRepres().get(1);
	}

	public abstract Map<Integer, Matrix> getPH();

    public Matrix getMu() {
    	Matrix aph_1 = getD0();
    	int size = Math.min(aph_1.numCols, aph_1.numRows);
    	Matrix res = new Matrix(size, 1, size);
    	for(int i = 0; i < size; i++) {
    		res.set(i, 0, -aph_1.get(i, i));
    	}
    	return res;
    }
    
    public Matrix getPhi() {
    	Map<Integer, Matrix> aph = this.getRepres();
    	Matrix ones = new Matrix(aph.get(0).numRows, 1, aph.get(0).numRows);
    	Matrix res = new Matrix(aph.get(1).numRows, 1);
    	Matrix mu = getMu();
    	
    	ones.fill(1.0);
    	aph.get(1).mult(ones, res);
    	res.divideRows(mu.nz_values, 0);
    	return res;
    }
    
    public Map<Integer, Matrix> getRepres() {
		try {
			if (this.representation == null)
				this.representation = getPH();
		} catch(Exception e) {
			e.printStackTrace(); //This function might not be implemented
		}
    	return this.representation;
    }

    public double getMean() {
    	Map<Integer, Matrix> rep = getRepres();
    	if (rep.get(0).hasNaN()) {
    		return Double.NaN;
    	} else {
    		Matrix D0 = getD0();
    		Matrix D1 = getD1();
    		return map_mean(D0,D1);
    	}
    }
    
    public double getSCV() {
    	Map<Integer, Matrix> rep = getRepres();
    	if (rep.get(0).hasNaN()) {
    		return Double.NaN;
    	} else {
			Matrix D0 = getD0();
			Matrix D1 = getD1();
			return map_scv(D0,D1);
    	}
    }
    
    public double evalCDF(double t) {
    	Map<Integer, Matrix> rep = getRepres();
		Matrix D0 = getD0();
		Matrix D1 = getD1();
		//return map.evalCDF(t); Not implemented in MAP
		throw new RuntimeException("Not implemented"); // TODO: not implemented
    }
    
    public double evalLST(double t) {
    	Map<Integer, Matrix> rep = getRepres();
		Matrix D0 = rep.get(0);
		Matrix D1 = rep.get(1);
		
		//Below is the function map_pie
		Matrix PIE = map_pie(D0,D1);
		
//		Matrix A = D0.copy();
//		Matrix e = new Matrix(PIE.numCols, 1);
//		CommonOps_DDRM.fill(e, 1);
//
//		//pie*inv(s*eye(size(A))-A)*(-A)*e
//		double[] diagEl = new double[Math.min(A.numRows, A.numCols)];
//		Arrays.fill(diagEl, 1);
//		Matrix eye = CommonOps_DDRM.diagR(A.numRows, A.numCols, diagEl);
//		CommonOps_DDRM.scale(t, eye);
//		Matrix res = new Matrix();
//		CommonOps_DDRM.subtract(eye, A, res);
//		CommonOps_DDRM.invert(res);
//		CommonOps_DDRM.mult(PIE, res.copy(), res);
//		CommonOps_DDRM.changeSign(A);
//		CommonOps_DDRM.mult(res.copy(), A, res);
//		CommonOps_DDRM.mult(res.copy(), e, res);
//    	return res.get(0,0);
		return 0.0;
    }
    
    public double getSkew() {
    	Map<Integer, Matrix> rep = getRepres();
    	if (rep.get(0).hasNaN()) {
    		return Double.NaN;
    	} else {
			Matrix D0 = rep.get(0);
			Matrix D1 = rep.get(1);
			MAP map = new MAP(D0, D1);
			List<Double> m = map.getMoments();
    		double M3 = m.get(2) - 3*m.get(1)*m.get(0) + 2*Math.pow(m.get(0), 3);
    		return M3 / Math.pow(Math.sqrt(map.getSCV())*m.get(0), 3);
    	}
    }
}
