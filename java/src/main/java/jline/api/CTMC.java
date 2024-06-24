package jline.api;

import java.util.Set;

import jline.util.Matrix;

/**
* APIs for Continuous-Time Markov Chains (CTMCs).
*/    
public class CTMC {
	/**
	 * Return the steady-state probability of a CTMC
	 *
	 * @param Q Infinitesimal generator of the CTMC
	 * @return Steady-state probability vector
	 */
	public static Matrix ctmc_solve(Matrix Q) {
		Matrix Qmat = Q.clone();
		if (Qmat.length() == 1) {
			Matrix p = new Matrix(1, 1, 1);
			p.set(0, 0, 1);
			return p;
		}

		Qmat = ctmc_makeinfgen(Qmat);
		int n = Qmat.length();

		// B = abs(Q+Q')>0
		Matrix B = Qmat.sub(1, Qmat.transpose());
		B.abs();
		for (int colIdx = 0; colIdx < B.getNumCols(); colIdx++) {
			int col1 = B.getColIndexes()[colIdx];
			int col2 = B.getColIndexes()[colIdx + 1];

			for (int i = col1; i < col2; i++) {
				if (B.getNonZeroValues()[i] > 0)
					B.getNonZeroValues()[i] = 1;
			}
		}

		// [nConnComp, connComp] = weaklyconncomp(B);
		Set<Set<Integer>> sets = Matrix.weaklyConnect(B, null);
		if (sets.size() > 1) {
			Matrix p = new Matrix(1, n);

			for (Set<Integer> set_c : sets) {
				// Qc = Q(connComp==c,connComp==c);
				Matrix Qc = new Matrix(set_c.size(), set_c.size());
				int Qc_row = 0, Qc_col = 0;
				for (Integer q_row : set_c) {
					for (Integer q_col : set_c) {
						Qc.set(Qc_row, Qc_col++, Qmat.get(q_row, q_col));
					}
					Qc_row++;
					Qc_col = 0;
				}
				// Qc = ctmc_makeinfgen(Qc);
				Qc = ctmc_makeinfgen(Qc);
				// p(connComp==c) = ctmc_solve(Qc);
				Matrix ctmc_solve_Qc = ctmc_solve(Qc);
				int idx = 0;
				for (Integer i : set_c) {
					p.set(0, i, ctmc_solve_Qc.get(0, idx++));
				}
			}
			p.divide(p.sumRows(0), p, true);
			return p;
		}

		if (Qmat.getNonZeroLength() == 0) {
			Matrix p = new Matrix(1, n);
			p.fill(1 / n);
			return p;
		}

		Matrix p = new Matrix(1, n);
		Matrix b = new Matrix(n, 1);
		Matrix nnzel = new Matrix(1, n);
		for (int i = 0; i < n; i++)
			nnzel.set(0, i, i);
		Matrix Qnnz = Qmat.clone();
		Matrix bnnz = b.clone();
		Matrix Qnnz_1 = Qmat.clone();
		Matrix bnnz_1 = bnnz.clone();

		boolean isReducible = false;
		boolean goon = true;
		while (goon) {
			// nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0);
			Matrix Qnnz_abs = Qnnz.clone();
			Qnnz_abs.abs();
			Matrix Qnnz_abs_sum_col = Qnnz_abs.sumCols();
			Matrix Qnnz_abs_sum_rows = Qnnz_abs.sumRows();
			Matrix find_res = new Matrix(1, Qnnz_abs_sum_col.getNumCols());
			for (int i = 0; i < Qnnz_abs_sum_col.getNumCols(); i++) {
				if (Qnnz_abs_sum_col.get(i) != 0 && Qnnz_abs_sum_rows.get(i) != 0)
					find_res.set(0, i, 1);
			}
			nnzel = find_res.find().transpose();

			if (nnzel.length() < n && !isReducible) {
				isReducible = true;
				// if (nargin > 1 && options.verbose == 2) % debug
				// fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
				// end
			}

			// Qnnz = Qnnz(nnzel, nnzel);
			Matrix new_Qnnz = new Matrix(nnzel.getNumCols(), nnzel.getNumCols());
			for (int i = 0; i < nnzel.getNumCols(); i++) {
				for (int j = 0; j < nnzel.getNumCols(); j++) {
					new_Qnnz.set(i, j, Qnnz.get((int) nnzel.get(0, i), (int) nnzel.get(0, j)));
				}
			}
			Qnnz = new_Qnnz;

			// bnnz = bnnz(nnzel);
			Matrix new_bnnz = new Matrix(nnzel.getNumCols(), 1);
			for (int i = 0; i < nnzel.getNumCols(); i++) {
				new_bnnz.set(i, 0, bnnz.get((int) nnzel.get(0, i), 0));
			}
			bnnz = new_bnnz;

			// Qnnz = ctmc_makeinfgen(Qnnz);
			Qnnz = ctmc_makeinfgen(Qnnz);

			if ((Qnnz.getNumCols() * Qnnz.getNumRows() == Qnnz_1.getNumCols() * Qnnz_1.getNumRows())
					&& (bnnz.getNumCols() * bnnz.getNumRows() == bnnz_1.getNumCols() * bnnz_1.getNumRows())) {
				goon = false;
			} else {
				Qnnz_1 = Qnnz.clone();
				bnnz_1 = bnnz.clone();
				nnzel = new Matrix(1, Qnnz.length());
				for (int i = 0; i < Qnnz.length(); i++)
					nnzel.set(0, i, i);
			}
		}

		if ((Qnnz == null) || (Qnnz.isEmpty())) {
			p.fill(1 / n);
			return p;
		}

		// Qnnz(:,end) = 1;
		for (int i = 0; i < n; i++)
			Qnnz.set(i, n - 1, 1.0);

		// bnnz(end) = 1;
		bnnz.set(n - 1, 0, 1.0);

		// p(nnzel)=Qnnz'\ bnnz;
		Matrix solve_res = new Matrix(0, 0);
		Matrix.solve(Qnnz.transpose(), bnnz, solve_res);
		solve_res.removeZeros(1e-10);
		for (int i = 0; i < nnzel.getNumCols(); i++)
			p.set(0, (int) nnzel.get(0, i), solve_res.get(i, 0));

		if (p.hasNaN()) {
			// B = abs(Qnnz+Qnnz')>0;
			B = Qnnz.add(1, Qnnz.transpose());
			B.abs();
			for (int colIdx = 0; colIdx < B.getNumCols(); colIdx++) {
				int col1 = B.getColIndexes()[colIdx];
				int col2 = B.getColIndexes()[colIdx + 1];

				for (int i = col1; i < col2; i++) {
					if (B.getNonZeroValues()[i] > 0)
						B.getNonZeroValues()[i] = 1;
				}
			}

			// [nConnComp, connComp] = weaklyconncomp(B);
			sets = Matrix.weaklyConnect(B, null);
			if (sets.size() > 1) {
				// p(nnzel) = zeros(1,n);
				for (int i = 0; i < nnzel.getNumCols(); i++)
					p.remove(0, (int) nnzel.get(0, i));

				for (Set<Integer> set_c : sets) {
					// Qc = Q(connComp==c,connComp==c);
					Matrix Qc = new Matrix(set_c.size(), set_c.size());
					int Qc_row = 0, Qc_col = 0;
					for (Integer q_row : set_c) {
						for (Integer q_col : set_c) {
							Qc.set(Qc_row, Qc_col++, Qmat.get(q_row, q_col));
						}
						Qc_row++;
						Qc_col = 0;
					}
					Qc = ctmc_makeinfgen(Qc);
					// p(intersect(find(connComp==c),nnzel)) = ctmc_solve(Qc);
					Matrix ctmc_solve_Qc = ctmc_solve(Qc);
					int idx = 0;
					for (int i = 0; i < nnzel.getNumCols(); i++) {
						int val = (int) nnzel.get(0, i);
						if (set_c.contains(val))
							p.set(0, val, ctmc_solve_Qc.get(0, idx++));
					}
				}
				p.divide(p.sumRows(0), p, true);
				return p;
			}
		}
		return p;
	}

	/**
	 * Return the transient probability distribution of the CTMC via the
	 * uniformaization method.
	 *
	 * @param pi0 Initial state the CTMC
	 * @param Q   Infinitesimal generator of the CTMC
	 * @param t   Transient analysis period boundary [0,t]
	 * @return Transient probability vector at time t
	 */
	public static Matrix ctmc_uniformization(Matrix pi0, Matrix Q, double t) {
		double tol = 1e-12;
		int maxiter = 100;
		double q = 0;
		int n = Q.getNumCols();
		for (int i = 0; i < n; i++) {
			q = Math.max(q, 1.1 * Math.abs(Q.get(i, i)));
		}
		Matrix Qs = Matrix.eye(n);
		Qs = Qs.add(1 / q, Q);
		int k = 0;
		double s = 1;
		double r = 1;
		int iter = 0;
		int kmax = 1;
		while (iter < maxiter) {
			iter++;
			k++;
			r = r * (q * t) / k;
			s = s + r;
			if (1 - Math.exp(-q * t) * s <= tol) {
				kmax = k;
				break;
			}
		}

		Matrix pi = new Matrix(1, n);
		Matrix riP = new Matrix(1, n);
		Matrix tmp = new Matrix(1, n);
		pi0.scale(Math.exp(-q * t), pi);
		Matrix P = new Matrix(pi0);
		double ri = Math.exp(-q * t);
		for (int j = 0; j < kmax; j++) {
			P.multEq(Qs);
			ri = ri * (q * t / (j + 1));
			pi = pi.add(ri, P);
		}
		return pi;
	}

	/**
	 * Normalize the input matrix diagonal to be an infinitesimal generator.
	 *
	 * @param Q Candidate infinitesimal generator of the CTMC
	 * @return Infinitesimal generator
	 */

	public static Matrix ctmc_makeinfgen(Matrix Q) {
		double[] val = new double[Q.length()];
		for (int i = 0; i < Q.length(); i++)
			val[i] = Q.get(i, i);

		Matrix diag_diag_Q = Matrix.diagMatrix(null, val, 0, val.length);
//		Matrix.diagMatrix(diag_diag_Q, val, 0, val.length);

		Matrix A = Q.sub(1, diag_diag_Q);

		Matrix diag_sum_A_row = new Matrix(0, 0);
		Matrix sum_A_row = A.sumRows();
		Matrix.diagMatrix(diag_sum_A_row, sum_A_row.getNonZeroValues(), 0, sum_A_row.getNonZeroLength());

		Q = A.sub(1, diag_sum_A_row);
		Q.removeZeros(0);
		return Q;
	}

	/**
	 * Compute the infinitesimal generator of the time-reserved CTMC
	 *
	 * @param Q Infinitesimal generator of the CTMC
	 * @return Infinitesimal generator of the time-reversed CTMC
	 */
	public static Matrix ctmc_timereverse(Matrix Q) {
		Matrix piq = CTMC.ctmc_solve(Q);
		Matrix Qrev = new Matrix(Q.getNumCols(), Q.getNumRows());
		for (int i = 0; i < Q.getNumRows(); i++) {
			for (int j = 0; j < Q.getNumCols(); j++) {
				Qrev.set(i, j, Q.get(i, j) * piq.get(i) / piq.get(j));
			}
		}
		return Qrev.transpose();
	}

	/**
	 * Form a random infinitesimal generator of CTMC
	 *
	 * @param length size of random matrix
	 * @return Infinitesimal generator of CTMC
	 */
	public static Matrix ctmc_rand(int length) {
		Matrix rand_matrix = new Matrix(length, length);
		rand_matrix.randMatrix(length);
		return ctmc_makeinfgen(rand_matrix);
	}

	/**
	 * Form a random infinitesimal generator of CTMC
	 *
	 * @param Q infinitesimal generator of CTMC
	 * @param pi0 initial state distribution vector
	 * @param n times of simulations
	 * @return Infinitesimal generator of CTMC
	 */
	public static CtmcSimulationResult ctmc_simulate(Matrix Q, double[] pi0, int n) {
		int numStates = Q.length();
		Random random = new Random();
		if (pi0 == null || pi0.length == 0){
			pi0 = new double[numStates];
			double sum = 0.0;

			for (int i = 0; i < numStates; i++){
				pi0[i] = random.nextDouble();
				sum += pi0[i];
			}

			for (int i = 0; i < numStates; i++){
				pi0[i] /= sum;
			}
		}
		double cumulative = 0.0;
		double r = random.nextDouble();
		int st = 0;
		for (int i = 0; i < pi0.length; i++){
			cumulative += pi0[i];
			if (r < cumulative){
				st = i;
				break;
			}
		}

    	Matrix F = new Matrix(numStates, numStates);
		for (int i = 0; i < numStates; i++) {
			double rowSum = 0.0;
			for (int j = 0; j < numStates; j++) {
				if (i != j) {
					rowSum += Q.get(i, j);
					F.set(i, j, rowSum);
				}
			}
			for (int j = 0; j < numStates; j++) {
				F.set(i, j, F.get(i, j)/rowSum);
			}
		}

		CtmcSimulationResult result = new CtmcSimulationResult();
		result.states = new int[n];
		result.sojournTimes = new double[n];

		for (int i = 0; i < n; i++) {
			result.states[i] = st;
			ExponentialDistribution expDist = new ExponentialDistribution(-1.0 / Q.get(st, st));
			result.sojournTimes[i] = expDist.sample();

			r = random.nextDouble();
			for (int j = 0; j < numStates; j++) {
				if (r < F.get(st, j)) {
					st = j;
					break;
				}
			}
		}
		return result;
	}

	public static class CtmcSimulationResult{
		public int [] states;
		public double[] sojournTimes;
	}

}
