package jline.api;

import jline.util.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
* APIs for Discrete-Time Markov Chains (DTMCs).
*/    
public class DTMC {

	/**
	 * Returns the stochastic complement of a DTMC
	 *
	 * @param P Transition matrix of the DTMC
	 * @param I Indexes of states to be kept in the stochastic complement
	 * @return Transition matrix of the stochastic complement
	 */
    public static Matrix dtmc_stochcomp(Matrix P, List<Integer> I) {
    	//Note that in this function, List is used instead of Matrix for performance consideration
    	if (P == null)
    		throw new RuntimeException("The first parameter of dtmc_stochcomp cannot be null");
    	
    	int lengthP = Math.max(P.getNumCols(), P.getNumRows());
    	if (I == null || I.size() == 0) {
    		I = new ArrayList<Integer>();
    		for(int i = 0; i < (int) Math.ceil(lengthP/2.0); i++)
    			I.add(i);
    	}
    	
    	List<Integer> Ic = IntStream.rangeClosed(0, lengthP-1).boxed().collect(Collectors.toList());
    	Ic.removeAll(I);

		Matrix P11 = new Matrix(I.size(), I.size());
		Matrix P12 = new Matrix(I.size(), Ic.size());
		Matrix P21 = new Matrix(Ic.size(), I.size());
		Matrix P22 = new Matrix(Ic.size(), Ic.size());
		
		for(int colIdx = 0; colIdx < P.getNumCols(); colIdx++) {
			for(int rowIdx = 0; rowIdx < P.getNumRows(); rowIdx++) {
				double value = P.get(rowIdx,colIdx);
				if (value > 0.0) {
					if (I.contains(colIdx)) {
						if (I.contains(rowIdx))
							P11.set(I.indexOf(rowIdx), I.indexOf(colIdx), value);
						else
							P21.set(Ic.indexOf(rowIdx), I.indexOf(colIdx), value);
					} else {
						if (I.contains(rowIdx))
							P12.set(I.indexOf(rowIdx), Ic.indexOf(colIdx), value);
						else
							P22.set(Ic.indexOf(rowIdx), Ic.indexOf(colIdx), value);
					}
				}
			}
		}
		
		double[] values = new double[Ic.size()];
		Arrays.fill(values, 1.0);
		Matrix S2 = Matrix.diagMatrix(null, values, 0, values.length).sub(1, P22);
		
		// S=P11+P12*(S2 \ P21);
		Matrix s2_p21 = new Matrix(0,0);
		Matrix.solve(S2, P21, s2_p21);
		Matrix S = P11.add(1, P12.mult(s2_p21, null));
    	return S;
    }
    
	/**
	 * Returns the steady-state solution of a DTMC
	 *
	 * @param P Transition matrix of the DTMC
	 * @return Steady-state solution vector of the DTMC
	 */
    public static Matrix dtmc_solve(Matrix P) {
    	
    	//P-eye(size(P))
    	for(int i = 0; i < P.getNumRows(); i++) {
    		P.set(i, i, P.get(i,i) - 1.0);
    	}
    	return CTMC.ctmc_solve(P);
    }

	/**
	 * Compute the infinitesimal generator of the time-reversed DTMC
	 *
	 * @param Q Infinitesimal generator of the DTMC
	 * @return Infinitesimal generator of the time-reversed DTMC
	 */
	public static Matrix dtmc_timereverse(Matrix P) {
		Matrix pie = DTMC.dtmc_solve(P);
		Matrix Prev = new Matrix(P.getNumCols(),P.getNumRows());
		//P-eye(size(P))
		for(int i = 0; i < P.getNumRows(); i++) {
			for(int j = 0; j < P.getNumCols(); j++) {
				Prev.set(i,j,P.get(i,j)*pie.get(i)/pie.get(j));
			}
		}
		return Prev.transpose();
	}

}
