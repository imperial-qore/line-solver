package jline.api;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Matrix;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.util.*;

import static jline.util.Matrix.oner;
import static jline.util.PopulationLattice.hashpop;
import static jline.util.PopulationLattice.pprod;

/**
* APIs for evaluating Product-Form Queueing Networks.
*/    
public class PFQN {

	/**
	 * Main method to compute the normalizing constant of a load-independent model
	 * @param lambda - arrival rate of open classes
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param options - solver options
	 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
	 */
	// TODO: Change the option parameter to varargin once Solver.parseOptions() is implemented.
	public static pfqnNcXQReturn pfqn_nc(Matrix lambda, Matrix L, Matrix N, Matrix Z, SolverOptions options) {
    /*
		SolverOptions options = Solver.parseOptions();
     */
		String method = options.method;
		int Rin = N.length();

		Double lG = null;
		Matrix X = new Matrix(0, 0);
		Matrix Q = new Matrix(0, 0);
		Matrix lambda_new = lambda.clone();
		Matrix L_new = L.clone();
		Matrix N_new = N.clone();
		Matrix Z_new = Z.clone();

		if (N_new.elementSum() < GlobalConstants.FineTol || N_new.isEmpty()) {
			lG = 0.0;
			return new pfqnNcXQReturn(lG, X, Q, method);
		}

		if (lambda_new.isEmpty()) {
			lambda_new = N_new.clone();
			lambda_new.fill(0.0);
		}

		Matrix Qopen = new Matrix(0, lambda_new.length());
		Matrix Ut = new Matrix(1, L_new.getNumRows());
		double lGopen = 0.0;
		for (int i = 0; i < L_new.getNumRows(); i++) {
			Matrix L_row_i = new Matrix(1, L_new.getNumCols());
			Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
			Ut.set(i, 1-lambda_new.mult(L_row_i.transpose()).get(0));
			if (Double.isNaN(Ut.get(i))) {
				Ut.set(i, 0);
			}
			for (int j = 0; j < L_row_i.length(); j++) {
				L_new.set(i, j, L_new.get(i, j)/Ut.get(i));
			}
			Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
			for (int j = 0; j < L_row_i.length(); j++) {
				L_row_i.set(j, L_row_i.get(j)/Ut.get(i));
			}
			Matrix tmp = lambda_new.elementMult(L_row_i, null);
			Qopen = Matrix.concatRows(Qopen, tmp, null);
		}
		Qopen.removeNaN();
		for (int i = 0; i < N_new.length(); i++) {
			if (Double.isInfinite(N_new.get(i))) {
				N_new.set(i, 0);
			}
		}

		Matrix L_tmp = new Matrix(L_new.getNumRows(), 0);
		Matrix N_tmp = new Matrix(N_new.getNumRows(), 0);
		Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
		Matrix lambda_tmp = new Matrix(lambda_new.getNumRows(), 0);
		for (int i = 0; i < N_new.length(); i++) {
			if (Math.abs(N_new.get(i)) >= GlobalConstants.FineTol) {
				Matrix L_col = Matrix.extractColumn(L_new, i, null);
				Matrix N_col = Matrix.extractColumn(N_new, i, null);
				Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
				Matrix lambda_col = Matrix.extractColumn(lambda_new, i, null);
				L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
				N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
				Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
				lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null);
			}
		}
		L_new = L_tmp;
		N_new = N_tmp;
		Z_new = Z_tmp;
		lambda_new = lambda_tmp;

		List<Integer> ocl = new ArrayList<>();
		for (int i = 0; i < N_new.length(); i++) {
			if (Double.isInfinite(N_new.get(i))) {
				ocl.add(i);
			}
		}
		int R = N_new.length();
		Matrix scalevec = new Matrix(1, R);
		scalevec.fill(1.0);
		for (int r = 0; r < R; r++) {
			Matrix L_col_r = new Matrix(L_new.getNumRows(), 1);
			Matrix.extract(L_new, 0, L_new.getNumRows(), r, r+1, L_col_r, 0, 0);
			Matrix Z_col_r = new Matrix(Z_new.getNumRows(), 1);
			Matrix.extract(Z_new, 0, Z_new.getNumRows(), r, r+1, Z_col_r, 0, 0);
			scalevec.set(r, Math.max(L_col_r.elementMax(), Z_col_r.elementMax()));
		}

		for (int i = 0; i < L_new.getNumRows(); i++) {
			for (int j = 0; j < L_new.getNumCols(); j++) {
				L_new.set(i, j, L_new.get(i, j)/scalevec.get(j));
			}
		}

		for (int i = 0; i < Z_new.length(); i++) {
			Z_new.set(i, Z_new.get(i)/scalevec.get(i));
		}

		Matrix Lsum = new Matrix(L_new.getNumRows(), 1);
		Matrix Lmax = new Matrix(L_new.getNumRows(), 1);

		for (int i = 0; i < L_new.getNumRows(); i++) {
			Matrix L_row_i = new Matrix(1, L_new.getNumCols());
			Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
			Lsum.set(i, L_new.sumRows(i));
			Lmax.set(i, L_row_i.elementMax());
		}

		List<Integer> demStations = new ArrayList<>();
		List<Integer> noDemStations = new ArrayList<>();
		L_tmp = new Matrix(0, L_new.getNumCols());

		for (int i = 0; i < L_new.getNumRows(); i++) {
			if (!Double.isNaN(Lmax.get(i)/Lsum.get(i)) && Lmax.get(i)/Lsum.get(i)> GlobalConstants.FineTol) {
				demStations.add(i);
				Matrix L_row_i = new Matrix(1, L_new.getNumCols());
				Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
				L_tmp = Matrix.concatRows(L_tmp, L_row_i,null);
			} else {
				noDemStations.add(i);
			}
		}
		L_new = L_tmp;

		boolean flag = false;
		for (int i = 0; i < N_new.length(); i++) {
			if (Math.abs(L_new.sumCols(i)+Z_new.sumCols(i))<GlobalConstants.FineTol && N_new.get(i) > GlobalConstants.FineTol) {
				flag = true;
				break;
			}
		}

		if (flag) {
			System.out.println("pfqn_nc warning: The model has no positive demands in any class.");
			if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
				lG = 0.0;
			} else {
				Matrix tmp1 = Z_new.sumCols();
				Matrix tmp2 = scalevec.clone();
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));
				}
				lG = -Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum() + N_new.mult(tmp2.transpose()).get(0);
			}
			return new pfqnNcXQReturn(lG, X, Q, method);
		}

		int M = L_new.getNumRows();
		R = L_new.getNumCols();

		if (L_new.isEmpty() || L_new.elementSum() < options.tol) {
			if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
				lG = lGopen;
			} else {
				Matrix tmp1 = Z_new.sumCols();
				Matrix tmp2 = scalevec.clone();
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));
				}
				lG = lGopen - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum()
						+ N_new.mult(tmp2.transpose()).get(0);
			}
			return new pfqnNcXQReturn(lG, X, Q, method);
		} else if (M == 1 && (Z_new.isEmpty() || Z_new.elementSum() < options.tol)) {
			Matrix tmp1 = L_new.sumCols();
			Matrix tmp2 = scalevec.clone();
			for (int i = 0; i < tmp1.length(); i++) {
				tmp1.set(i, Math.log(tmp1.get(i)));
				tmp2.set(i, Math.log(tmp2.get(i)));
			}
			lG = Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum()
					+ N_new.mult(tmp2.transpose()).get(0);
			return new pfqnNcXQReturn(lG, X, Q, method);
		}

		List<Integer> zeroDemandClasses = new ArrayList<>();
		List<Integer> nonzeroDemandClasses = new ArrayList<>();

		for (int i = 0; i < R; i++) {
			if (L_new.sumCols(i) < options.tol) {
				zeroDemandClasses.add(i);
			} else {
				nonzeroDemandClasses.add(i);
			}
		}

		double lGzdem;
		Matrix Nz;
		Matrix Zz = new Matrix(Z_new.getNumRows(), 0);
		for (int i: zeroDemandClasses) {
			Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
			Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i+1, Z_col_i, 0, 0);
			Zz = Matrix.concatColumns(Zz, Z_col_i, null);
		}

		flag = true;
		for (int i = 0; i < Zz.getNumCols(); i++) {
			if (Zz.sumCols(i) >= options.tol) {
				flag = false;
				break;
			}
		}
		if (Z_new.isEmpty() || flag) {
			lGzdem = 0.0;
			Nz = new Matrix(1, 1);
			Nz.fill(0.0);
		} else {
			if (zeroDemandClasses.isEmpty()) {
				lGzdem = 0.0;
				Nz = new Matrix(1, 1);
				Nz.fill(0.0);
			} else {
				Nz = new Matrix(1, 0);
				for (int i: zeroDemandClasses) {
					Matrix N_col_i = new Matrix(1, 1);
					Matrix.extract(N_new, 0, 1, i, i+1, N_col_i, 0, 0);
					Nz = Matrix.concatColumns(Nz, N_col_i, null);
				}

				Matrix tmp1 = Zz.sumCols();
				Matrix tmp2 = new Matrix(1, 0);
				for (int i: zeroDemandClasses) {
					Matrix scalevec_col_i = new Matrix(1, 1);
					Matrix.extract(scalevec, 0, 1, i, i+1, scalevec_col_i, 0, 0);
					tmp2 = Matrix.concatColumns(tmp2, scalevec_col_i, null);
				}
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));

				}

				lGzdem = -Matrix.factln(Nz).elementSum() + Nz.elementMult(tmp1, null).elementSum()
						+ Nz.mult(tmp2.transpose()).get(0);
			}
		}

		L_tmp = new Matrix(L_new.getNumRows(), 0);
		N_tmp = new Matrix(1, 0);
		Z_tmp = new Matrix(Z_new.getNumRows(), 0);
		Matrix scalevecz = new Matrix(1, 0);

		for (int i: nonzeroDemandClasses) {
			Matrix L_col_i = new Matrix(L_new.getNumRows(), 1);
			Matrix N_col_i = new Matrix(1, 1);
			Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
			Matrix scalevec_col_i = new Matrix(1, 1);
			Matrix.extract(L_new, 0, L_new.getNumRows(), i, i+1, L_col_i, 0, 0);
			Matrix.extract(N_new, 0, 1, i, i+1, N_col_i, 0, 0);
			Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i+1, Z_col_i, 0, 0);
			Matrix.extract(scalevec, 0, 1, i, i+1, scalevec_col_i, 0, 0);
			L_tmp = Matrix.concatColumns(L_tmp, L_col_i, null);
			N_tmp = Matrix.concatColumns(N_tmp, N_col_i, null);
			Z_tmp = Matrix.concatColumns(Z_tmp, Z_col_i, null);
			scalevecz = Matrix.concatColumns(scalevecz, scalevec_col_i, null);
		}
		L_new = L_tmp;
		N_new = N_tmp;
		Z_new = Z_tmp;

		pfqnNcXQReturn ret = compute_norm_const(L_new, N_new, Z_new, options);
		Double lGnzdem = ret.lG;
		Matrix Xnnzdem = ret.X;
		Matrix Qnnzdem = ret.Q;
		method = ret.method;

		if (Xnnzdem.isEmpty()) {
			X = new Matrix(0, 0);
			Q = new Matrix(0, 0);
		}

		Matrix tmp = scalevecz.clone();
		for (int i = 0; i < tmp.length(); i++) {
			tmp.set(i, Math.log(tmp.get(i)));
		}
		lG = lGopen + lGnzdem + lGzdem + N_new.mult(tmp.transpose()).get(0);
		return new pfqnNcXQReturn(lG, X, Q, method);
	}

	/**
	 * Main method to compute the normalizing constant of a load-dependent model
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param mu - load-dependent scalings
	 * @param options - solver options
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_ncld(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
		Double lG = Double.NaN;
		Double G = Double.NaN;
		String method = options.method;

		Matrix mu_new;
		if ((int) N.elementSum() >= mu.getNumCols()) {
			mu_new = mu.clone();
		} else {
			mu_new = new Matrix(mu.getNumRows(), 0);
			for (int i = 0; i < N.elementSum(); i++) {
				Matrix mu_col_i = new Matrix(mu.getNumRows(), 1);
				Matrix.extract(mu, 0, mu.getNumRows(), i, i+1, mu_col_i, 0, 0);
				mu_new = Matrix.concatColumns(mu_new, mu_col_i, null);
			}
		}

		Matrix L_new = new Matrix(L.getNumRows(), 0);
		Matrix N_new = new Matrix(N.getNumRows(), 0);
		Matrix Z_new = new Matrix(Z.getNumRows(), 0);
		for (int i = 0; i < N.length(); i++) {
			if (Math.abs(N.get(i)) >= GlobalConstants.FineTol) {
				Matrix L_col_i = Matrix.extractColumn(L, i, null);
				Matrix N_col_i = Matrix.extractColumn(N, i, null);
				Matrix Z_col_i = Matrix.extractColumn(Z, i, null);
				L_new = Matrix.concatColumns(L_new, L_col_i, null);
				N_new = Matrix.concatColumns(N_new, N_col_i, null);
				Z_new = Matrix.concatColumns(Z_new, Z_col_i, null);
			}
		}

		int R = N_new.getNumCols();
		Matrix scalevec = new Matrix(1,R);
		scalevec.fill(1.0);
		for (int r = 0; r < R; r++) {
			Matrix L_col_r = Matrix.extractColumn(L_new, r, null);
			Matrix Z_col_r = Matrix.extractColumn(Z_new, r, null);
			scalevec.set(r, Math.max(L_col_r.elementMax(), Z_col_r.elementMax()));
		}

		for (int i = 0; i < L_new.getNumRows(); i++) {
			for (int j = 0; j < L_new.getNumCols(); j++) {
				L_new.set(i, j, L_new.get(i, j)/scalevec.get(j));
			}
		}

		for (int j = 0; j < Z_new.getNumCols(); j++) {
			Z_new.set(j, Z_new.get(j)/scalevec.get(j));
		}

		Matrix Lsum = new Matrix(L_new.getNumRows(), 1);
		Matrix Lmax = new Matrix(L_new.getNumRows(), 1);

		for (int i = 0; i < L_new.getNumRows(); i++) {
			Matrix L_row_i = new Matrix(1, L_new.getNumCols());
			Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
			Lsum.set(i, L_new.sumRows(i));
			Lmax.set(i, L_row_i.elementMax());
		}

		List<Integer> demStations = new ArrayList<>();
		Matrix L_tmp = new Matrix(0, L_new.getNumCols());
		Matrix mu_tmp = new Matrix(0, mu_new.getNumCols());

		for (int i = 0; i < L_new.getNumRows(); i++) {
			if (!Double.isNaN(Lmax.get(i)/Lsum.get(i)) && Lmax.get(i)/Lsum.get(i)> GlobalConstants.FineTol) {
				demStations.add(i);
				Matrix L_row_i = new Matrix(1, L_new.getNumCols());
				Matrix.extract(L_new, i, i+1, 0, L_new.getNumCols(), L_row_i, 0, 0);
				L_tmp = Matrix.concatRows(L_tmp, L_row_i,null);
				Matrix mu_row_i = Matrix.extractRows(mu_new, i, i+1, null);
				mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null);
			}
		}
		L_new = L_tmp.clone();
		mu_new = mu_tmp.clone();

		boolean flag = false;
		for (int i = 0; i < N_new.getNumCols(); i++) {
			if (Math.abs(L_new.sumCols(i)+Z_new.sumCols(i))<GlobalConstants.FineTol && N_new.get(i) > GlobalConstants.FineTol) {
				flag = true;
				break;
			}
		}

		if (flag) {
			System.out.println("pfqn_ncld warning: The model has no positive demands in any class.");
			if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
				lG = 0.0;
			} else {
				Matrix tmp1 = Z_new.sumCols();
				Matrix tmp2 = scalevec.clone();
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));
				}
				lG = -Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum() + N_new.mult(tmp2.transpose()).get(0);
			}
			G = Double.NaN;
			return new pfqnNcReturn(G, lG, method);
		}

		int M = L_new.getNumRows();
		R = L_new.getNumCols();

		if (L_new.isEmpty() || L_new.elementSum() < options.tol) {
			if (Z_new.isEmpty() || Z_new.elementSum() < options.tol) {
				lG = 0.0;
			} else {
				Matrix tmp1 = Z_new.sumCols();
				Matrix tmp2 = scalevec.clone();
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));
				}
				lG = -Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum()
						+ N_new.mult(tmp2.transpose()).get(0);
			}
			return new pfqnNcReturn(G, lG, method);
		} else if (M == 1 && (Z_new.isEmpty() || Z_new.elementSum() < options.tol)) {
			Matrix tmp1 = L_new.sumCols();
			Matrix tmp2 = scalevec.clone();
			for (int i = 0; i < tmp1.length(); i++) {
				tmp1.set(i, Math.log(tmp1.get(i)));
				tmp2.set(i, Math.log(tmp2.get(i)));
			}
			Matrix tmp3 = mu_new.clone();
			for (int i = 0; i < tmp3.length(); i++) {
				tmp3.set(i, Math.log(tmp3.get(i)));
			}
			lG = Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum()
					+ N_new.mult(tmp2.transpose()).get(0) - tmp3.elementSum();
			return new pfqnNcReturn(G, lG, method);
		}

		List<Integer> zeroDemandClasses = new ArrayList<>();
		List<Integer> nonzeroDemandClasses = new ArrayList<>();

		for (int i = 0; i < R; i++) {
			if (L_new.sumCols(i) < options.tol) {
				zeroDemandClasses.add(i);
			} else {
				nonzeroDemandClasses.add(i);
			}
		}

		double lGzdem;
		Matrix Nz;
		Matrix Zz = new Matrix(Z_new.getNumRows(), 0);
		for (int i: zeroDemandClasses) {
			Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
			Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i+1, Z_col_i, 0, 0);
			Zz = Matrix.concatColumns(Zz, Z_col_i, null);
		}

		flag = true;
		for (int i = 0; i < Zz.getNumCols(); i++) {
			if (Zz.sumCols(i) >= options.tol) {
				flag = false;
				break;
			}
		}
		if (Z_new.isEmpty() || flag) {
			lGzdem = 0.0;
			Nz = new Matrix(1, 1);
			Nz.fill(0.0);
		} else {
			if (zeroDemandClasses.isEmpty()) {
				lGzdem = 0.0;
				Nz = new Matrix(1, 1);
				Nz.fill(0.0);
			} else {
				Nz = new Matrix(1, 0);
				for (int i: zeroDemandClasses) {
					Matrix N_col_i = new Matrix(1, 1);
					Matrix.extract(N_new, 0, 1, i, i+1, N_col_i, 0, 0);
					Nz = Matrix.concatColumns(Nz, N_col_i, null);
				}

				Matrix tmp1 = Zz.sumCols();
				Matrix tmp2 = new Matrix(1, 0);
				for (int i: zeroDemandClasses) {
					Matrix scalevec_col_i = new Matrix(1, 1);
					Matrix.extract(scalevec, 0, 1, i, i+1, scalevec_col_i, 0, 0);
					tmp2 = Matrix.concatColumns(tmp2, scalevec_col_i, null);
				}
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, Math.log(tmp1.get(i)));
					tmp2.set(i, Math.log(tmp2.get(i)));

				}

				lGzdem = -Matrix.factln(Nz).elementSum() + Nz.elementMult(tmp1, null).elementSum()
						+ Nz.mult(tmp2.transpose()).get(0);
			}
		}

		L_tmp = new Matrix(L_new.getNumRows(), 0);
		Matrix N_tmp = new Matrix(1, 0);
		Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
		Matrix scalevecz = new Matrix(1, 0);

		for (int i: nonzeroDemandClasses) {
			Matrix L_col_i = new Matrix(L_new.getNumRows(), 1);
			Matrix N_col_i = new Matrix(1, 1);
			Matrix Z_col_i = new Matrix(Z_new.getNumRows(), 1);
			Matrix scalevec_col_i = new Matrix(1, 1);
			Matrix.extract(L_new, 0, L_new.getNumRows(), i, i+1, L_col_i, 0, 0);
			Matrix.extract(N_new, 0, 1, i, i+1, N_col_i, 0, 0);
			Matrix.extract(Z_new, 0, Z_new.getNumRows(), i, i+1, Z_col_i, 0, 0);
			Matrix.extract(scalevec, 0, 1, i, i+1, scalevec_col_i, 0, 0);
			L_tmp = Matrix.concatColumns(L_tmp, L_col_i, null);
			N_tmp = Matrix.concatColumns(N_tmp, N_col_i, null);
			Z_tmp = Matrix.concatColumns(Z_tmp, Z_col_i, null);
			scalevecz = Matrix.concatColumns(scalevecz, scalevec_col_i, null);
		}
		L_new = L_tmp;
		N_new = N_tmp;
		Z_new = Z_tmp;

		double lGnnzdem;
		if (N_new.elementMin() < 0.0) {
			lGnnzdem = 0.0;
		} else {
			pfqnNcReturn ret = compute_norm_const_ld(L_new, N_new, Z_new, mu_new, options);
			lGnnzdem = ret.lG;
			method = ret.method;
		}

		Matrix tmp = scalevecz.clone();
		for (int i = 0; i < tmp.length(); i++) {
			tmp.set(i, Math.log(tmp.get(i)));
		}
		lG = lGnnzdem + lGzdem + N_new.mult(tmp.transpose()).get(0);
		G = Math.exp(lG);
		return new pfqnNcReturn(G, lG, method);
	}

	/**
	 * Compute the normalizing constant of a single-class load-dependent model
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param mu - load-depedent scalings
	 * @param options - solver options
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_gld(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		Matrix lambda = new Matrix(1, R);
		Double G;
		Double lG;

		if (M == 1) {
			Matrix N_tmp = new Matrix(1, 0);
			Matrix L_tmp = new Matrix(1, 0);
			for (int i = 0; i < R; i++) {
				if (L.get(i) > GlobalConstants.FineTol) {
					Matrix N_tmp2 = new Matrix(1, 1);
					N_tmp2.fill(N.get(i));
					Matrix L_tmp2 = new Matrix(1, 1);
					L_tmp2.fill(Math.log(L.get(0, i)));
					N_tmp = Matrix.concatColumns(N_tmp, N_tmp2, null);
					L_tmp = Matrix.concatColumns(L_tmp, L_tmp2, null);
				}
			}
			Matrix mu_new;
			if ((int) N.elementSum() >= mu.getNumCols()) {
				mu_new = Matrix.extractRows(mu, 0, 1, null);
			} else {
				mu_new = new Matrix(1, 0);
				for (int i = 0; i < N.elementSum(); i++) {
					Matrix mu_col_i = new Matrix(1, 1);
					Matrix.extract(mu, 0, 1, i, i+1, mu_col_i, 0, 0);
					mu_new = Matrix.concatColumns(mu_new, mu_col_i, null);
				}
			}

			for (int i = 0; i < mu_new.length(); i++) {
				mu_new.set(i, Math.log(mu_new.get(i)));
			}

			lG = Maths.factln(N.elementSum()) - Matrix.factln(N).elementSum()
					+ N_tmp.mult(L_tmp.transpose()).get(0) - mu_new.elementSum();
			G = Math.exp(lG);
			return new pfqnNcReturn(G, lG);
		}

		if (R == 1) {
			pfqnNcReturn ret = pfqn_gldsingle(L, N, mu, null);
			lG = ret.lG;
			G = ret.G;
			return new pfqnNcReturn(G, lG);
		}

		if (L == null || L.isEmpty()) {
			G = 0.0;
			lG = Double.NEGATIVE_INFINITY;
			return new pfqnNcReturn(G, lG);
		}

		Matrix mu_new;
		if (mu == null) {
			mu_new = new Matrix(M, (int)N.elementSum());
			mu_new.fill(1.0);
		} else {
			mu_new = mu.clone();
		}

		SolverOptions options_new;
		if (options == null) {
			options_new = SolverNC.defaultOptions();
		} else {
			options_new = options;
		}

		boolean isLoadDep = false;
		boolean[] isInfServer = new boolean[M];
		for (int i = 0; i < M; i++) {
			Matrix mu_row_i = new Matrix(1, (int)N.elementSum());
			Matrix.extract(mu_new, i, i+1, 0, (int)N.elementSum(), mu_row_i, 0, 0);
			boolean flag = true;
			for (int j = 0; j < mu_row_i.getNumCols(); j++) {
				if (Math.abs(mu_row_i.get(j)-(j+1))>GlobalConstants.FineTol) {
					flag = false;
					break;
				}
			}

			if (Math.abs(mu_row_i.elementMin()-1)<GlobalConstants.FineTol && Math.abs(mu_row_i.elementMax()-1)<GlobalConstants.FineTol) {
				isInfServer[i] = false;
			} else if (flag) {
				isInfServer[i] = true;
			} else {
				isInfServer[i] = false;
				isLoadDep = true;
			}
		}

		if (!isLoadDep) {
			Matrix Lli = new Matrix(0, L.getNumCols());
			Matrix Zli = new Matrix(0, L.getNumCols());
			for (int i = 0; i < M; i++) {
				Matrix L_row_i = Matrix.extractRows(L, i, i+1, null);
				if (isInfServer[i]) {
					Zli = Matrix.concatRows(Zli, L_row_i, null);
				} else {
					Lli = Matrix.concatRows(Lli, L_row_i, null);
				}
			}
			if (Lli.isEmpty()) {
				Lli = N.clone();
				Lli.fill(0.0);
			}
			if (Zli.isEmpty()) {
				Zli = N.clone();
				Zli.fill(0.0);
			}
			options_new.method = "exact";
			lG = pfqn_nc(lambda, Lli, N, Zli.sumCols(), options_new).lG;
			G = Math.exp(lG);
			return new pfqnNcReturn(G, lG);
		}

		G = 0.0;
		if (M == 0) {
			lG = Math.log(G);
			return new pfqnNcReturn(G, lG);
		}

		if (Math.abs(N.elementMax())<GlobalConstants.FineTol && Math.abs(N.elementMin())<GlobalConstants.FineTol) {
			G = 1.0;
			lG = Math.log(G);
			return new pfqnNcReturn(G, lG);
		}

		if (R == 1) {
			G = pfqn_gldsingle(L, N, mu_new, null).G;
			lG = Math.log(G);
			return new pfqnNcReturn(G, lG);
		}

		G += pfqn_gld(Matrix.extractRows(L, 0, M-1, null), N,
				Matrix.extractRows(mu_new, 0, M-1, null), options_new).G;

		for (int r = 0; r < R; r++) {
			if (N.get(r) > GlobalConstants.FineTol) {
				Matrix N_1 = N.clone();
				if (R > 1) {
					N_1.set(r, N_1.get(r)-1);
				} else {
					for (int i = 0; i < N_1.length(); i++) {
						N_1.set(i, N_1.get(i)-1);
					}
				}
				G += L.get(M-1, r)/mu_new.get(M-1, 0)
						*pfqn_gld(L, N_1, pfqn_mushift(mu, M-1), options_new).G;
			}
		}
		lG = Math.log(G);
		return new pfqnNcReturn(G, lG);
	}

	/**
	 * Shifst the a load-dependent scaling vector by one position
	 * @param mu - load-dependent scalings
	 * @return normalizing constant and its logarithm
	 */
	public static Matrix pfqn_mushift(Matrix mu, int k) {
		int M = mu.getNumRows();
		int N = mu.getNumCols();
		Matrix mushift = new Matrix(M, N-1);
		Matrix.extract(mu, 0, M, 0, N-1, mushift, 0, 0);
		for (int j = 0; j < N-1; j++) {
			mushift.set(k, j, mu.get(k, j+1));
		}
		return mushift;
	}

	/**
	 * Run a normalizing constant solution method in a load-independent model
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param options - solver options
	 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
	 */
	public static pfqnNcXQReturn compute_norm_const(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		Matrix X = new Matrix(0, 0);
		Matrix Q = new Matrix(0, 0);
		String method = options.method;
		Double lG = null;

		// TODO: several methods missing
		switch (options.method) {
			case "ca": {
				pfqnNcReturn ret = pfqn_ca(L,N,Z.sumCols());
				lG = ret.lG;
				break;
			}
			case "adaptive":
			case "default": {
				Matrix Z_colSum = Z.sumCols();
				if (M > 1) {
					if (R == 1 || (R <= 3 && N.elementSum() < 50)) {
						pfqnNcReturn ret = pfqn_ca(L,N,Z.sumCols());
						lG = ret.lG;
						method = "ca";
					} else {
						if (M <= R) {
							pfqnNcReturn ret = pfqn_le(L,N,Z.sumCols());
							lG = ret.lG;
							method = "le";
						}
					}
				} else if (Z_colSum.getNumCols() == 1 && Math.abs(Z_colSum.get(0)) < GlobalConstants.FineTol) {
					Matrix tmp = L.clone();
					for (int i = 0; i < tmp.length(); i++) {
						tmp.set(i, Math.log(tmp.get(i)));
					}
					lG = -N.mult(tmp.transpose()).get(0);
					method = "exact";
				} else {
					if (N.elementMax() < 10000) {
						pfqnNcReturn ret = pfqn_mmint2_gausslegendre(L,N,Z.sumCols(),null);
						lG = ret.lG;
						method = "gleint";
					} else {
						pfqnNcReturn ret = pfqn_le(L,N,Z.sumCols());
						lG = ret.lG;
						method = "le";
					}
				}
				break;
			}
			case "sampling": {
				pfqnNcReturn ret = pfqn_ls(L,N,Z.sumCols(),options.samples);
				lG = ret.lG;
				method = "ls";
				break;
			}
			case "mmint2":
			case "gleint": {
				if (L.getNumRows() > 1) {
					throw new RuntimeException("The " + options.method + " method requires a model with a delay and a single queueing station.");
				} else {
					pfqnNcReturn ret = pfqn_mmint2_gausslegendre(L,N,Z.sumCols(),null);
					lG = ret.lG;
				}
				break;
			}
			case "le": {
				pfqnNcReturn ret = pfqn_le(L,N,Z.sumCols());
				lG = ret.lG;
				break;
			}
			case "ls": {
				pfqnNcReturn ret = pfqn_ls(L,N,Z.sumCols(),options.samples);
				lG = ret.lG;
				break;
			}
			case "mom": {
				// TODO: incomplete
				if (N.length() <= 1) {
					pfqnNcReturn ret = pfqn_ca(L,N,Z.sumCols());
					lG = ret.lG;
				}
				break;
			}
			case "exact": {
				if (M >= R || N.elementSum() > 10 || Z.sumCols().elementMin() > 0) {
					pfqnNcReturn ret = pfqn_ca(L,N,Z.sumCols());
					lG = ret.lG;
					method = "ca";
				}
				break;
			}
			case "comom": {
				if (R > 1) {
					try {
						if (M <= 1) {
							pfqnComomrmReturn ret = pfqn_comomrm(L,N,Z,1,options.tol);
							lG = ret.lG;
						}
					} catch (Exception e) {
						e.printStackTrace();
						lG = Double.NaN;
					}
				} else {
					pfqnNcReturn ret = pfqn_ca(L,N,Z.sumCols());
					lG = ret.lG;
					method = "ca";
				}
				break;
			}
			default:
				throw new RuntimeException("Unrecognized method: " + options.method);
		}
		return new pfqnNcXQReturn(lG, X, Q, method);
	}

	/**
	 * Run a normalizing constant solution method in a load-dependent model
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param mu - load-depedent scalings
	 * @param options - solver options
	 * @return normalizing constant and its logarithm
	 */
	//TODO: other cases to be implemented
	public static pfqnNcReturn compute_norm_const_ld(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		String method = options.method;
		Double lG = null;

		switch (options.method) {
			case "default":
			case "exact": {
				int D = Z.getNumRows();
				Matrix Lz = Matrix.concatRows(L, Z, null);
				Matrix tmp = new Matrix(1, mu.getNumCols());
				for (int i = 0; i < tmp.length(); i++) {
					tmp.set(i, i+1);
				}
				Matrix muz = Matrix.concatRows(mu, tmp.repmat(D, 1), null);
				if (R == 1) {
					lG = pfqn_gldsingle(Lz, N, muz, options).lG;
				} else if (!(M == 1 && Z.elementMax() > GlobalConstants.FineTol)) {
					pfqnNcReturn ret = pfqn_gld(Lz, N, muz, options);
					lG = ret.lG;
				}
				method = "exact";
				break;
			}
			default:
				throw new RuntimeException("Unrecognized method: " + options.method);
		}
		Double G = Math.exp(lG);
		return new pfqnNcReturn(G, lG, method);
	}

	public static pfqnNcReturn pfqn_panacea(Matrix L, Matrix N, Matrix Z) {
		return pfqn_panacea(L, N, Z, new SolverOptions());
	}

	/**
	 * Compute the PANACEA approximation
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param options - solver options
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_panacea(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
		String method = options.method;
		int M = L.getNumRows();
		int R = L.getNumCols();
		double lG = Double.NaN;
		double G = Double.NaN;
		if (Z.isEmpty() || Z.elementSum() < options.tol) {
			Z = Z.clone();
			Z = Z.add(GlobalConstants.FineTol, Z);
		}


		if (L.isEmpty() || L.elementSum() < options.tol) {
			{
				Matrix tmp1 = Z.sumCols();
				for (int i = 0; i < tmp1.length(); i++) {
					tmp1.set(i, N.get(i) * Math.log(tmp1.get(i)));
				}
				lG = -Matrix.factln(N).elementSum() + tmp1.elementSum();
			}
			G = Math.exp(lG);
			return new pfqnNcReturn(G, lG, method);
		}

		// Compute r = L./repmat(Z,q,1)
		Matrix r = new Matrix(M, R);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < R; j++) {
				r.set(i, j, L.get(i, j) / Z.get(0, j));
			}
		}

		// Find Nt = max(1./r)
		double Nt = Double.MIN_VALUE;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < R; j++) {
				Nt = Math.max(Nt, 1.0 / r.get(i, j));
			}
		}

		// Compute beta = N / Nt
		Matrix beta = new Matrix(1, R);
		for (int j = 0; j < R; j++) {
			beta.set(0, j, (double) N.get(0, j) / Nt);
		}

		// Compute gamma = r * Nt
		Matrix gamma = new Matrix(M, R);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < R; j++) {
				gamma.set(i, j, r.get(i, j) * Nt);
			}
		}

		// Compute alpha = 1 - N * r'
		Matrix alpha = new Matrix(1, M);
		for (int i = 0; i < M; i++) {
			double sum = 0;
			for (int j = 0; j < R; j++) {
				sum += N.get(j) * r.get(i, j);
			}
			alpha.set(0, i, 1 - sum);
		}

		// Compute gammatilde = gamma ./ repmat(alpha', 1, p)
		Matrix gammatilde = new Matrix(M, R);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < R; j++) {
				gammatilde.set(i, j, gamma.get(i, j % R) / alpha.get(0, i));
			}
		}

		// Check if min(alpha) < 0
		double minAlpha = Double.MAX_VALUE;
		for (int i = 0; i < M; i++) {
			minAlpha = Math.min(minAlpha, alpha.get(0, i));
		}
		if (minAlpha < 0) {
			System.out.println("Warning: Model is not in normal usage");
			lG = Double.NaN;
			G = Double.NaN;
			return new pfqnNcReturn(G, lG, method);
		}

		double A0 = 1;
		double A1 = 0;
		for (int j = 0; j < R; j++) {
			Matrix m = new Matrix(1,R);
			m.set(0,j,2.0);
			A1 -= beta.get(j) * pfqn_ca(gammatilde, m).G;
		}

		double A2 = 0;
		for (int j = 0; j < R; j++) {
			Matrix m = new Matrix(1,R);
			m.set(0,j,3.0);
			A2 += 2 * beta.get(j) * pfqn_ca(gammatilde, m).G;

			m.set(0,j,4.0);
			A2 += 3 * Math.pow(beta.get(j), 2) * pfqn_ca(gammatilde, m).G;

			for (int k = 0; k < R; k++) {
				if (k == j) continue;
				m.zero();
				m.set(0,j,2.0);
				m.set(0,k,2.0);
				A2 = A2 + 0.5 * beta.get(j) * beta.get(k) * pfqn_ca(gammatilde, m).G;
			}
		}

		//TODO: to be completed
		Matrix tmp1 = Z.sumCols();
		for (int i = 0; i < tmp1.length(); i++) {
			tmp1.set(i, N.get(i) * Math.log(tmp1.get(i)));
		}
		lG = -Matrix.factln(N).elementSum() + tmp1.elementSum();
		lG += Math.log(A0 + A1/Nt + A2/Nt/Nt);
		for (int s = 0; s < R; s++) {
			lG -= Math.log(alpha.get(0,s));
		}
		G = Math.exp(lG);

		// if ~isfinite(lGn)
		if (!Double.isFinite(lG)) {
			G = Double.NaN;
			lG = Double.NaN;
		}
		return new pfqnNcReturn(G, lG, method);
	}

	/**
	 * Compute the product-form factor relatively to a Delay station
	 * @param Z - think times at the Delay station
	 * @param n - number of jobs for each class
	 * @return product of terms Z[k]^n[k]/n[k]! for all classes k
	 */
	public static double pfqn_pff_delay(Matrix Z, Matrix n) {
		int R = n.length();
		if (n.sumRows().sumCols().get(0, 0) == 0) {
			return 1.0;
		}

		double f = 0;
		for (int r = 0; r < R; r++) {
			if (Z.get(r) > 0) {
				f += Math.log(Z.get(r)) * n.get(r);
				f -= Maths.factln(n.get(r));
			} else if (n.get(r) > 0) {
				return 0.0;
			}
		}
		return Math.exp(f);
	}

	/**
	 * Compute the normalizing constant using the convolution algorithm
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_ca(Matrix L, Matrix N) {
		Matrix Z = N.clone();
		Z.zero();
		return pfqn_ca(L, N, Z);
	}

		/**
         * Compute the normalizing constant using the convolution algorithm
         * @param L - demands at all stations
         * @param N - number of jobs for each class
         * @param Z - think times
         * @return normalizing constant and its logarithm
         */
	public static pfqnNcReturn pfqn_ca(Matrix L, Matrix N, Matrix Z) {
		int M = L.getNumRows();
		int R = L.getNumCols();

		if (M == 0) {
			// Calculate the logarithm of the factorial of every element in N
			Matrix tmp = new Matrix(1, N.length());
			for (int i = 0; i < N.length(); i++) {
				tmp.set(0, i, -Maths.factln(N.get(i)));
			}
			double lGn = tmp.sumRows().sumCols().get(0, 0);

			Matrix tmp2 = Z.sumCols();
			for (int i = 0; i < tmp2.length(); i++) {
				tmp2.set(i, Math.log(tmp2.get(i)));
			}
			if (N.length() == 1) {
				lGn += (N.get(0)*tmp2.sumRows().get(0));
			} else if (tmp2.length() == 1) {
				lGn += (tmp2.get(0)*N.sumRows().sumCols().get(0, 0));
			} else {
				Matrix tmp3 = new Matrix(1, N.length());
				for (int i = 0; i < N.length(); i++) {
					tmp3.set(i, N.get(i));
				}
				lGn += tmp3.elementMult(tmp2, null).sumRows().get(0);
			}
			double Gn = Math.exp(lGn);

			return new pfqnNcReturn(Gn, lGn);
		}

		if (N.elementMin() < 0) {
			return new pfqnNcReturn(0.0, Double.NEGATIVE_INFINITY);
		}

		if (N.sumRows().sumCols().get(0) == 0) {
			return new pfqnNcReturn(1.0, 0.0);
		}

		if (Z == null || Z.isEmpty()) {
			Z = new Matrix(1, R);
			Z.fill(0.0);
		}

		int product_N_plus_one = 1;
		for (int i = 0; i < N.length(); i++) {
			product_N_plus_one *= (N.get(i)+1);
		}
		Matrix G = new Matrix(M+1, product_N_plus_one);
		G.fill(1.0);
		Matrix n = pprod(N);

		while (Math.abs(n.sumRows().sumCols().get(0)+1) > GlobalConstants.FineTol) {
			int idxn = hashpop(n, N);
			G.set(0, idxn, pfqn_pff_delay(Z, n));
			for (int m = 1; m < M+1; m++) {
				G.set(m, idxn, G.get(m-1, idxn));
				for (int r = 0; r < R; r++) {
					if (n.get(r) >= 1) {
						n.set(r, n.get(r)-1);
						int idxn_1r = hashpop(n, N);
						n.set(r, n.get(r)+1);
						double tmp_res = G.get(m, idxn) + L.get(m-1, r)*G.get(m,idxn_1r);
						G.set(m, idxn, tmp_res);
					}
				}
			}
			n = pprod(n, N);
		}

		double Gn = G.get(M, G.getNumCols()-1);
		double lGn = Math.log(Gn);
		return new pfqnNcReturn(Gn, lGn);
	}

	/**
	 * Compute the normalizing constant of a repairmen model using Gauss-Legendre integration
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think times
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_mmint2_gausslegendre(Matrix L, Matrix N, Matrix Z, Integer m) {
		if (m == null) {
			m = 1;
		}

		List<Double> gausslegendreNodes = new ArrayList<>();
		List<Double> gausslegendreWeights = new ArrayList<>();

		if (gausslegendreNodes.isEmpty()) {
			Scanner scan;
			File nodeFile;
			File weightFile;
			try {
				nodeFile = new File(PFQN.class.getResource("/gausslegendre-nodes.txt").toURI());
				weightFile = new File(PFQN.class.getResource("/gausslegendre-weights.txt").toURI());
				try {
					scan = new Scanner(nodeFile);
					while (scan.hasNextDouble()) {
						gausslegendreNodes.add(scan.nextDouble());
					}

					scan = new Scanner(weightFile);
					while (scan.hasNextDouble()) {
						gausslegendreWeights.add(scan.nextDouble());
					}
				} catch (FileNotFoundException e1) {
					e1.printStackTrace();
				}
			} catch (URISyntaxException e2) {
				e2.printStackTrace();
			}

		}

		int n = (int) Math.max(300, Math.min(gausslegendreNodes.size(),
						2*(N.sumRows().sumCols().get(0)+m-1)-1));
		Matrix y = new Matrix(1, n);
		y.fill(0.0);

		if (!(Z.getNumRows() == L.getNumRows() && Z.getNumCols() == L.getNumCols())) {
			throw new RuntimeException("The dimensions of Z and L are not the same.");
		}
		for (int i = 0; i < n; i++) {
			Matrix tmp = L.clone();

			for (int j = 0; j < tmp.getNumRows(); j++) {
				for (int k = 0; k < tmp.getNumCols(); k++) {
					tmp.set(j, k, Math.log(Z.get(j, k) + gausslegendreNodes.get(i)*tmp.get(j, k)));
				}
			}
			y.set(i, (N.mult(tmp.transpose())).get(0, 0));
		}

		Matrix g = y.clone();
		Matrix nodes = new Matrix(1, n);
		Matrix logNodes = new Matrix(1, n);
		Matrix logWeights = new Matrix(1, n);
		for (int i = 0; i < n; i++) {
			nodes.set(i, gausslegendreNodes.get(i));
			logNodes.set(i, Math.log(gausslegendreNodes.get(i)));
			logWeights.set(i, Math.log(gausslegendreWeights.get(i)));
		}

		for (int i = 0; i < n; i++) {
			g.set(i, g.get(i)+logWeights.get(i)-nodes.get(i));
		}

		double coeff = 0;
		for (int i = 0; i < N.length(); i++) {
			coeff -= Maths.factln(N.get(i));
		}
		coeff -= Maths.factln(m-1);
		coeff += (m-1) * logNodes.elementSum();

		double lG = 0.0;
		for (int i = 0; i < g.length(); i++) {
			lG += Math.exp(g.get(i));
		}
		lG = Math.log(lG) + coeff;
		if (!Double.isFinite(lG)) {
			lG = Matrix.logsumexp(g) + coeff;
		}
		double G = Math.exp(lG);
		return new pfqnNcReturn(G, lG);
	}

	/**
	 * Auxiliary function to compute the Hessian used in the logistic expansion method in models with delays
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think times
	 * @param u - term appearing in the integrand
	 * @param v - term appearing in the integrand
	 * @return Hessian matrix
	 */
	protected static Matrix pfqn_le_hessianZ(Matrix L, Matrix N, Matrix Z, Matrix u, double v) {
		int K = L.getNumRows();
		int R = L.getNumCols();
		double Ntot = N.elementSum();
		Matrix A = new Matrix(K, K);
		A.fill(0.0);
		Matrix csi = new Matrix(1, R);
		for (int r = 0; r < R; r++) {
			csi.set(r, N.get(r)/(Z.get(r)+v*u.mult(Matrix.extractColumn(L, r, null)).get(0)));
		}
		Matrix Lhat = new Matrix(K, R);
		Lhat.fill(0.0);
		for (int k = 0; k < K; k++) {
			for (int r = 0; r < R; r++) {
				Lhat.set(k, r, Z.get(r) + v*L.get(k,r));
			}
		}
		double eta = Ntot + K;
		for (int i = 0; i < K; i++) {
			for (int j = 0; j < K; j++) {
				if (i != j) {
					A.set(i, j, -eta*u.get(i)*u.get(j));
					for (int r = 0; r < R; r++) {
						A.set(i, j, A.get(i, j)
										+ csi.get(r)*csi.get(r)*Lhat.get(i,r)*Lhat.get(j,r)*u.get(i)*u.get(j)/N.get(r));
					}
				}
			}
		}
		for (int i = 0; i < K; i++) {
			A.set(i, i, -Matrix.allbut(Matrix.extractRows(A, i, i+1, null), i).elementSum());
		}
		Matrix tmp_A = new Matrix(K, K);
		tmp_A.fill(0.0);
		Matrix.extract(A, 0, K-1, 0, K-1, tmp_A, 0, 0);
		A = tmp_A;
		A.set(K-1, K-1, 1.0);

		for (int r = 0; r < R; r++) {
			Matrix L_col_r = new Matrix(L.getNumRows(), 1);
			Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
			A.set(K-1, K-1,
							A.get(K-1, K-1)-(csi.get(r)*csi.get(r)/N.get(r))*Z.get(r)*u.mult(L_col_r).get(0));
		}

		A.set(K-1, K-1, v*A.get(K-1, K-1));

		for (int i = 0; i < K-1; i++) {
			for (int r = 0; r < R; r++) {
				Matrix L_col_r = new Matrix(L.getNumRows(), 1);
				Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
				A.set(i, K-1, A.get(i, K-1)+
								v*u.get(i)*((csi.get(r)*csi.get(r)/N.get(r))*Lhat.get(i,r)
												*(u.mult(L_col_r).get(0))-csi.get(r)*L.get(i,r)));
			}
			A.set(K-1, i, A.get(i, K-1));
		}

		return A;
	}

	/**
	 * Auxiliary function to compute the Hessian used in the logistic expansion method
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param u0 - term appearing in the integrand
	 * @return normalizing constant and its logarithm
	 */
	protected static Matrix pfqn_le_hessian(Matrix L, Matrix N, Matrix u0) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		double Ntot = N.elementSum();
		Matrix hu = new Matrix(M-1, M-1);
		hu.fill(0.0);

		for (int i = 0; i < M-1; i++) {
			for (int j = 0; j < M-1; j++) {
				if (i != j) {
					hu.set(i, j, -(Ntot+M)*u0.get(i)*u0.get(j));
					for (int r = 0; r < R; r++) {
						Matrix L_col_r = new Matrix(L.getNumRows(), 1);
						Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
						hu.set(i, j, hu.get(i,j)
										+N.get(r)*L.get(i,r)*L.get(j,r)*u0.get(i)*u0.get(j)
										/(u0.mult(L_col_r).get(0)*u0.mult(L_col_r).get(0)));
					}
				} else {
					hu.set(i, j, (Ntot+M)*u0.get(i)*(Matrix.allbut(u0,i).elementSum()));
					for (int r = 0; r < R; r++) {
						Matrix L_col_r = new Matrix(L.getNumRows(), 1);
						Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
						Matrix tmp_L = Matrix.allbut(L_col_r, i).transpose();

						hu.set(i, j, hu.get(i,j)
										-N.get(r)*L.get(i,r)*u0.get(i)
										*(Matrix.allbut(u0,i).mult(tmp_L).get(0))
										/(u0.mult(L_col_r).get(0)*u0.mult(L_col_r).get(0)));
					}
				}
			}
		}
		return hu;
	}

	/**
	 * Fixed-point iteration used in the logistic expansion method
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @return fixed point
	 */
	protected static pfqnLeFpiReturn pfqn_le_fpi(Matrix L, Matrix N) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		Matrix u = new Matrix(M, 1);
		for (int i = 0; i < M; i++) {
			u.set(i, 1.0/M);
		}
		Matrix u_1 = new Matrix(M, 1);
		u_1.fill(Double.POSITIVE_INFINITY);
		Matrix d = new Matrix(0, M);
		Matrix u_abs_diff = new Matrix(M, 1);

		for (int i = 0; i < M; i++) {
			u_abs_diff.set(i, Math.abs(u.get(i)-u_1.get(i)));
		}

		while (u_abs_diff.elementSum() > 1e-10) {
			u_1 = u.clone();
			for (int i = 0; i < M; i++) {
				u.set(i, 1/(N.elementSum()+M));
				for (int r = 0; r < R; r++) {
					Matrix L_col_r = new Matrix(L.getNumRows(), 1);
					Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
					u.set(i, u.get(i)
						+N.get(r)/(N.elementSum()+M)*L.get(i,r)*u_1.get(i)
						/(u_1.transpose().mult(L_col_r).get(0)));
				}
			}

			for (int i = 0; i < M; i++) {
				u_abs_diff.set(i, Math.abs(u.get(i)-u_1.get(i)));
			}


			d = Matrix.concatRows(d, u_abs_diff.transpose(), null);
		}

		return new pfqnLeFpiReturn(u, d);
	}

	/**
	 * Fixed-point iteration used in the logistic expansion method in models with delays
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @return fixed point
	 */
	protected static pfqnLeFpiZReturn pfqn_le_fpiZ(Matrix L, Matrix N, Matrix Z) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		double eta = N.elementSum() + M;
		Matrix u = new Matrix(M, 1);
		for (int i = 0; i < M; i++) {
			u.set(i, 1.0/M);
		}
		double v = eta + 1;
		Matrix u_1 = new Matrix(M, 1);
		u_1.fill(Double.POSITIVE_INFINITY);
		double v_1 = Double.POSITIVE_INFINITY;
		Matrix d = new Matrix(0, M);
		Matrix u_abs_diff = new Matrix(M, 1);

		for (int i = 0; i < M; i++) {
			u_abs_diff.set(i, Math.abs(u.get(i)-u_1.get(i)));
		}

		while (u_abs_diff.elementSum() > 1e-10) {
			u_1 = u.clone();
			v_1 = v;
			for (int i = 0; i < M; i++) {
				u.set(i, 1/eta);
				for (int r = 0; r < R; r++) {
					Matrix L_col_r = new Matrix(L.getNumRows(), 1);
					Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
					u.set(i, u.get(i)
									+N.get(r)/eta*(Z.get(r)+v*L.get(i,r))*u_1.get(i)
									/(Z.get(r)+v*u_1.transpose().mult(L_col_r).get(0)));
				}
			}

			Matrix xi = new Matrix(1, R);
			for (int r = 0; r < R; r++) {
				Matrix L_col_r = new Matrix(L.getNumRows(), 1);
				Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
				xi.set(r, N.get(r)/(Z.get(r)+v*u_1.transpose().mult(L_col_r).get(0)));
			}
			v = eta + 1;
			for (int r = 0; r < R; r++) {
				v -= xi.get(r) * Z.get(r);
			}

			for (int i = 0; i < M; i++) {
				u_abs_diff.set(i, Math.abs(u.get(i)-u_1.get(i)));
			}

			d = Matrix.concatRows(d, u_abs_diff.transpose().elementIncrease(Math.abs(v-v_1)), null);
		}

		return new pfqnLeFpiZReturn(u, v, d);
	}

	// TODO: javadoc
	protected static double simplex_fun(double[] x, Matrix L, Matrix N) {
		int M = x.length + 1;
		Matrix v = new Matrix(1, M);
		for (int i = 0; i < M-1; i++) {
			v.set(i, Math.exp(x[i]));
		}
		v.set(M-1, 1.0);

		Matrix tmp = (v.mult(L)).transpose();
		for (int i = 0; i < tmp.length(); i++) {
			tmp.set(i, Math.log(tmp.get(i)));
		}

		double x_sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			x_sum += x[i];
		}

		double f = Math.exp(N.mult(tmp).elementSum() + x_sum
							- (N.elementSum()+M)*Math.log(v.elementSum()));

		return f;
	}

	/**
	 * Logistic expansion method to compute the normalizing constant
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_le(Matrix L, Matrix N, Matrix Z) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		double lGn;
		double Gn;

		if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0 || L.elementSum() < GlobalConstants.CoarseTol) {
			Matrix tmp = new Matrix(1, Z.getNumCols());
			for (int i = 0; i < tmp.length(); i++) {
				tmp.set(i, Math.log(Z.sumCols(i)));
			}
			lGn = - Matrix.factln(N).elementSum()
							+ N.elementMult(tmp, null).elementSum();
			Gn = Math.exp(lGn);
		} else if (Z == null || Z.isEmpty()) {
			pfqnLeFpiReturn ret = pfqn_le_fpi(L,N);
			Matrix umax = ret.u;
			Matrix A = pfqn_le_hessian(L,N,umax.transpose());
			double S = 0.0;
			for (int r = 0; r < R; r++) {
				Matrix L_col_r = new Matrix(L.getNumRows(), 1);
				Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
				S += N.get(r) * Math.log(umax.transpose().mult(L_col_r).get(0));
			}

			Matrix tmp = new Matrix(1, N.length()+1);
			Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0);
			tmp.set(N.length(), M-1);
			Matrix log_umax = umax.clone();
			for (int i = 0; i < log_umax.length(); i++) {
				log_umax.set(i, Math.log(log_umax.get(i)));
			}
			lGn = Maths.multinomialln(tmp) + Maths.factln(M-1) + (M-1)*Math.log(Math.sqrt(2*Math.PI))
							- Math.log(Math.sqrt(A.det()))
							+ log_umax.elementSum() + S;
			Gn = Math.exp(lGn);
		} else {
			pfqnLeFpiZReturn ret = pfqn_le_fpiZ(L, N, Z);
			Matrix umax = ret.u;
			double vmax = ret.v;
			Matrix A = pfqn_le_hessianZ(L, N, Z, umax.transpose(), vmax);
			double S = 0;
			for (int r = 0; r < R; r++) {
				Matrix L_col_r = new Matrix(L.getNumRows(), 1);
				Matrix.extract(L, 0, L.getNumRows(), r, r+1, L_col_r, 0, 0);
				S += N.get(r) * Math.log(Z.get(r)+vmax*umax.transpose().mult(L_col_r).get(0));
			}
			Matrix log_umax = umax.clone();
			for (int i = 0; i < log_umax.length(); i++) {
				log_umax.set(i, Math.log(log_umax.get(i)));
			}
			lGn = -Matrix.factln(N).elementSum() -vmax + M*Math.log(vmax) + M*Math.log(Math.sqrt(2*Math.PI))
							- Math.log(Math.sqrt(A.det()))
							+ log_umax.elementSum() + S;
			Gn = Math.exp(lGn);
		}
		return new pfqnNcReturn(Gn, lGn);
	}

	/**
	 * Logistic sampling method to compute the normalizing constant
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param I - number of samples
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_ls(Matrix L, Matrix N, Matrix Z, int I) {
		return pfqn_ls(L,N,Z,I,23000);
	}

	/**
	 * Logistic sampling method to compute the normalizing constant
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param I - number of samples
	 * @param seed - random number generation seed
	 * @return normalizing constant and its logarithm
	 */
	public static pfqnNcReturn pfqn_ls(Matrix L, Matrix N, Matrix Z, int I, long seed) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		Matrix Lsum = new Matrix(M, 1);
		for (int i = 0; i < M; i++) {
			Lsum.set(i, L.sumRows(i));
		}
		Matrix L_new = new Matrix(0, R);
		for (int i = 0; i < M; i++) {
			Matrix L_row_i = new Matrix(1, R);
			Matrix.extract(L, i, i+1, 0, R, L_row_i, 0, 0);
			if (Lsum.get(i) > GlobalConstants.CoarseTol) {
				L_new = Matrix.concatRows(L_new, L_row_i, null);
			}
		}
		M = L_new.getNumRows();
		R = L_new.getNumCols();
		double[][] sample = null;
		double lGn;
		double Gn;

		if (L_new.isEmpty() || N.isEmpty() || N.elementSum() == 0 || L_new.elementSum() < GlobalConstants.CoarseTol) {
			Matrix tmp = new Matrix(1, Z.getNumCols());
			for (int i = 0; i < tmp.length(); i++) {
				tmp.set(i, Math.log(Z.sumCols(i)));
			}
			lGn = - Matrix.factln(N).elementSum()
							+ N.elementMult(tmp, null).elementSum();
		} else if (Z == null || Z.isEmpty()) {
			pfqnLeFpiReturn ret = pfqn_le_fpi(L_new,N);
			Matrix umax = ret.u;
			Matrix A = pfqn_le_hessian(L_new,N,umax.transpose());
			Matrix A_t = A.transpose();
			for (int i = 0; i < A.getNumRows(); i++) {
				for (int j = 0; j < A.getNumCols(); j++) {
					A.set(i, j, (A.get(i,j)+A_t.get(i,j))/2.0);
				}
			}
			Matrix iA = A.inv();
			Matrix x0 = new Matrix(1, M-1);
			for (int i = 0; i < M-1; i++) {
				x0.set(i, Math.log(umax.get(i)/umax.get(M-1)));
			}

			double[] x0_array = new double[M-1];
			for (int i = 0; i < M-1; i++) {
				x0_array[i] = x0.get(i);
			}
			double[][] iA_array = new double[iA.getNumRows()][];
			for (int i = 0; i < iA.getNumRows(); i++) {
				double[] tmp_row = new double[iA.getNumCols()];
				for (int j = 0; j < iA.getNumCols(); j++) {
					tmp_row[j] = iA.get(i, j);
				}
				iA_array[i] = tmp_row;
			}
			MultivariateNormalDistribution mvd = new MultivariateNormalDistribution(x0_array, iA_array);
			mvd.reseedRandomGenerator(seed);

			if (sample == null) {
				sample = mvd.sample(I);
			}

			Matrix T = new Matrix(I, 1);
			Matrix dpdf = new Matrix(1, I);
			for (int i = 0; i < I; i++) {
				double[] sample_i = sample[i];
				T.set(i, simplex_fun(sample_i, L_new, N));
				dpdf.set(i, mvd.density(sample_i));
			}

			Matrix tmp = new Matrix(1, N.length()+1);
			Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0);
			tmp.set(N.length(), M-1);
			double div_sum = 0.0;
			for (int i = 0; i < I; i++) {
				div_sum += T.get(i)/dpdf.get(i);
			}
			lGn = Maths.multinomialln(tmp) + Maths.factln(M-1) + Math.log(div_sum/I);
		} else {
			pfqnLeFpiZReturn ret = pfqn_le_fpiZ(L_new, N, Z);
			Matrix umax = ret.u;
			double vmax = ret.v;
			Matrix A = pfqn_le_hessianZ(L_new, N, Z, umax.transpose(), vmax);
			Matrix A_t = A.transpose();
			for (int i = 0; i < A.getNumRows(); i++) {
				for (int j = 0; j < A.getNumCols(); j++) {
					A.set(i, j, (A.get(i,j)+A_t.get(i,j))/2.0);
				}
			}
			Matrix iA = A.inv();
			Matrix x0 = new Matrix(1, M);
			for (int i = 0; i < M-1; i++) {
				x0.set(i, Math.log(umax.get(i)/umax.get(M-1)));
			}
			x0.set(M-1, Math.log(vmax));

			double[] x0_array = new double[M];
			for (int i = 0; i < M; i++) {
				x0_array[i] = x0.get(i);
			}
			double[][] iA_array = new double[iA.getNumRows()][];
			for (int i = 0; i < iA.getNumRows(); i++) {
				double[] tmp_row = new double[iA.getNumCols()];
				for (int j = 0; j < iA.getNumCols(); j++) {
					tmp_row[j] = iA.get(i, j);
				}
				iA_array[i] = tmp_row;
			}
			MultivariateNormalDistribution mvd = new MultivariateNormalDistribution(x0_array, iA_array);

			if (sample == null) {
				sample = mvd.sample(I);
			}

			Matrix T = new Matrix(I, 1);
			double epsilon = 1e-10;
			double eN = epsilon * N.elementSum();
			double eta = N.elementSum() + M*(1+eN);
			int K = M;
			Matrix dpdf = new Matrix(1, I);

			for (int i = 0; i < I; i++) {
				double[] sample_i = sample[i];
				T.set(i, pfqn_ls_helper(sample_i, K, M, eta, eN, L_new, N, Z));
				dpdf.set(i, mvd.density(sample_i));
			}

			double div_sum = 0.0;
			for (int i = 0; i < I; i++) {
				div_sum += T.get(i)/dpdf.get(i);
			}
			lGn = Math.log(Math.exp(-Matrix.factln(N).elementSum())*div_sum/I);
		}
		Gn = Math.exp(lGn);
		return new pfqnNcReturn(Gn, lGn);
	}

	/**
	 * Auxiliary function used in the logistic sampling method
	 */
	protected static double pfqn_ls_helper(double[] x, int K, int M, double eta, double eN, Matrix L, Matrix N, Matrix Z) {
		double res = -Math.exp(x[K-1])+K*(1+eN)*x[M-1];
		double tmp1 = 0.0;
		double tmp2 = 0.0;
		for (int i = 0; i < K-1; i++) {
			tmp1 += x[i];
			tmp2 += Math.exp(x[i]);
		}
		tmp2 = Math.log(tmp2+1);
		tmp2 *= -eta;
		res += (tmp1 + tmp2);

		Matrix L_row_K = new Matrix(1, L.getNumCols());
		Matrix L_first_K_minus_one_rows = new Matrix(K-1, L.getNumCols());
		Matrix x_first_K_minus_one_elements = new Matrix(1, K-1);
		Matrix.extract(L, K-1, K, 0, L.getNumCols(), L_row_K, 0, 0);
		Matrix.extract(L, 0, K-1, 0, L.getNumCols(), L_first_K_minus_one_rows, 0, 0);
		for (int i = 0; i < K-1; i++) {
			x_first_K_minus_one_elements.set(i, Math.exp(x[i]));
		}
		for (int i = 0; i < K-1; i++) {
			for (int j = 0; j < L.getNumCols(); j++) {
				L_first_K_minus_one_rows.set(i, j,
								L_first_K_minus_one_rows.get(i, j)*Math.exp(x[K-1])+Z.get(j));
			}
		}
		x_first_K_minus_one_elements = x_first_K_minus_one_elements.mult(L_first_K_minus_one_rows);
		for (int i = 0; i < L_row_K.length(); i++) {
			L_row_K.set(i, Math.log(L_row_K.get(i)*Math.exp(x[K-1])+Z.get(i)+x_first_K_minus_one_elements.get(i)));
		}
		res += N.mult(L_row_K.transpose()).elementSum();

		res = Math.exp(res);
		return res;
	}

	/**
	 * Sanitizes product-form model parameters to avoid degeneracies
	 * @param lambda - arrival rates for open classes
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param atol - absolute numerical tolerance
	 * @return sanitized parameters
	 */
	public static pfqnNcSanitizeReturn pfqn_nc_sanitize(Matrix lambda, Matrix L, Matrix N, Matrix Z, double atol) {
		Matrix L_new = L.clone();
		Matrix Z_new = Z.clone();
		L_new.removeNaN();
		Z_new.removeNaN();
		Matrix L_tmp = new Matrix(L_new.getNumRows(), 0);
		Matrix N_tmp = new Matrix(N.getNumRows(), 0);
		Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
		Matrix lambda_tmp = new Matrix(lambda.getNumRows(), 0);
		for (int i = 0; i < N.length(); i++) {
			if (Math.abs(N.get(i)) >= GlobalConstants.FineTol && !(L_new.sumCols(i)+Z_new.sumCols(i) < atol)) {
				Matrix L_col = Matrix.extractColumn(L_new, i, null);
				Matrix N_col = Matrix.extractColumn(N, i, null);
				Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
				Matrix lambda_col = Matrix.extractColumn(lambda, i, null);
				L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
				N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
				Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
				lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null);
			}
		}
		L_new = L_tmp;
		Z_new = Z_tmp;
		Matrix N_new = N_tmp;
		Matrix lambda_new = lambda_tmp;

		double lGremaind= 0.0;

		Matrix L_zeroDemand = new Matrix(L_new.getNumRows(), 0);
		Matrix Z_zeroDemand = new Matrix(Z_new.getNumRows(), 0);
		Matrix N_zeroDemand = new Matrix(N_new.getNumRows(), 0);
		L_tmp = new Matrix(L_new.getNumRows(), 0);
		N_tmp = new Matrix(N_new.getNumRows(), 0);
		Z_tmp = new Matrix(Z_new.getNumRows(), 0);

		for (int i = 0; i < L_new.getNumCols(); i++) {
			Matrix L_col = Matrix.extractColumn(L_new, i, null);
			Matrix N_col = Matrix.extractColumn(N_new, i, null);
			Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
			if (L_new.get(i) < atol) {
				L_zeroDemand = Matrix.concatColumns(L_zeroDemand, L_col, null);
				Z_zeroDemand = Matrix.concatColumns(Z_zeroDemand, Z_col, null);
				N_zeroDemand = Matrix.concatColumns(N_zeroDemand, N_col, null);
			} else {
				L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
				N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
				Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
			}
		}

		Matrix log_Z_zeroDemand = Z_zeroDemand.clone();
		Matrix log_N_zeroDemand = N_zeroDemand.clone();

		for (int i = 0; i < Z_zeroDemand.getNumRows(); i++) {
			for (int j = 0; j < Z_zeroDemand.getNumCols(); j++) {
				log_Z_zeroDemand.set(i, j, Math.log(Z_zeroDemand.get(i, j)));
			}
		}

		for (int i = 0; i < N_zeroDemand.getNumRows(); i++) {
			for (int j = 0; j < N_zeroDemand.getNumCols(); j++) {
				log_N_zeroDemand.set(i, j, Math.log(N_zeroDemand.get(i, j)));
			}
		}

		if (!L_zeroDemand.isEmpty()) {
			lGremaind += (N_zeroDemand.mult(log_Z_zeroDemand.transpose()).get(0)-log_N_zeroDemand.elementSum());
		}
		L_new = L_tmp;
		Z_new = Z_tmp;
		N_new = N_tmp;

		if (!L_new.isEmpty()) {
			Matrix Lmax = new Matrix(1, L_new.getNumCols());
				for (int i = 0; i < Lmax.length(); i++) {
					Matrix L_col_i = Matrix.extractColumn(L_new, i, null);
					Lmax.set(i, L_col_i.elementMax());
				}
			if (Lmax.isEmpty()) {
				Lmax = new Matrix(1, Z_new.getNumCols());
				Lmax.ones();
			}
			Matrix repmat_Lmax_L = Lmax.repmat(L_new.getNumRows(), 1);
			Matrix repmat_Lmax_Z = Lmax.repmat(Z_new.getNumRows(), 1);
			for (int i = 0; i < L_new.getNumRows(); i++) {
				for (int j = 0; j < L_new.getNumCols(); j++) {
					L_new.set(i, j, L_new.get(i, j) / repmat_Lmax_L.get(i, j));
				}
			}

			for (int i = 0; i < Z_new.getNumRows(); i++) {
				for (int j = 0; j < Z_new.getNumCols(); j++) {
					Z_new.set(i, j, Z_new.get(i, j) / repmat_Lmax_Z.get(i, j));
				}
			}

			Matrix Lmax_log = Lmax.clone().transpose();
			for (int i = 0; i < Lmax_log.length(); i++) {
				Lmax_log.set(i, Math.log(Lmax_log.get(i)));
			}
			lGremaind += N_new.mult(Lmax_log).get(0);

			if (!Z_new.isEmpty()) {
				Integer[] index = new Integer[Z_new.getNumCols()];
				for (int i = 0; i < index.length; i++) {
					index[i] = i;
				}

				Matrix Z_sum_col = new Matrix(1, Z_new.getNumCols());
				for (int i = 0; i < Z_sum_col.length(); i++) {
					Z_sum_col.set(i, Z_new.sumCols(i));
				}

				Arrays.sort(index, new Comparator<Integer>() {
					public int compare(Integer i1, Integer i2) {
						return ((Double)Z_sum_col.get(i1)).compareTo(Z_sum_col.get(i2));
					}
				});

				L_tmp = new Matrix(L_new.getNumRows(), 0);
				Z_tmp = new Matrix(Z_new.getNumRows(), 0);
				N_tmp = new Matrix(N_new.getNumRows(), 0);

				for (int i = 0; i < index.length; i++) {
					if (!L_new.isEmpty()) {
						L_tmp = Matrix.concatColumns(L_tmp, Matrix.extractColumn(L_new, index[i], null), null);
					}
					Z_tmp = Matrix.concatColumns(Z_tmp, Matrix.extractColumn(Z_new, index[i], null), null);
					N_tmp = Matrix.concatColumns(N_tmp, Matrix.extractColumn(N_new, index[i], null), null);
				}
				L_new = L_tmp;
				Z_new = Z_tmp;
				N_new = N_tmp;
			}
		}

		return new pfqnNcSanitizeReturn(lambda_new, L_new, N_new, Z_new, lGremaind);
	}

	/**
	 * Compute the normalizing constant of a repairmen model using COMOM
	 * @param L - demands at all stations
	 * @param N - number of jobs for each class
	 * @param Z - think time for each class
	 * @param m - multiplicy of queueing station
	 * @param atol - absolute numerical tolerance
	 * @return sanitized parameters
	 */
	public static pfqnComomrmReturn pfqn_comomrm(Matrix L, Matrix N, Matrix Z, Integer m, double atol) {
		int M = L.getNumRows();
		int R = L.getNumCols();
		if (M != 1) {
			throw new RuntimeException("pfqn_comomrm: The solver accepts at most a single queueing station.");
		}
		if (m == null) {
			m = 1;
		}
		Matrix lambda = N.clone();
		lambda.fill(0.0);
		pfqnNcSanitizeReturn ret = pfqn_nc_sanitize(lambda, L, N, Z, atol);
		Matrix L_new = ret.L;
		Matrix N_new = ret.N;
		Matrix Z_new = ret.Z;
		double lG0 = ret.lGremaind;
		List<Integer> zerothinktimes = new ArrayList<>();

		for (int i = 0; i < Z_new.length(); i++) {
			if (Z_new.get(i) < GlobalConstants.FineTol) {
				zerothinktimes.add(i);
			}
		}
		Matrix nvec = new Matrix(1, R);
		nvec.fill(0.0);

		Matrix lh;

		if (!zerothinktimes.isEmpty()) {
			for (int i = 0; i < zerothinktimes.size(); i++) {
				nvec.set(zerothinktimes.get(i), N_new.get(zerothinktimes.get(i)));
			}
			lh = new Matrix(0, 1);
			Matrix tmp = new Matrix(1, 1);
			tmp.set(0, Maths.factln(nvec.elementSum()+m)- Matrix.factln(nvec).elementSum());
			lh = Matrix.concatRows(lh, tmp, null);

			for (int i = 0; i < zerothinktimes.size(); i++) {
				Matrix nvec_s = nvec.clone();
				nvec_s.set(i, nvec_s.get(i)-1);
				tmp.set(0, Maths.factln(nvec_s.elementSum()+m)- Matrix.factln(nvec_s).elementSum());
				lh = Matrix.concatRows(lh, tmp, null);
			}
			tmp.set(0, Maths.factln(nvec.elementSum()+m-1)- Matrix.factln(nvec).elementSum());
			lh = Matrix.concatRows(lh, tmp, null);
			for (int i = 0; i < zerothinktimes.size(); i++) {
				Matrix nvec_s = nvec.clone();
				nvec_s.set(i, nvec_s.get(i)-1);
				tmp.set(0, Maths.factln(nvec_s.elementSum()+m-1)- Matrix.factln(nvec_s).elementSum());
				lh = Matrix.concatRows(lh, tmp, null);
			}
		} else {
			lh = new Matrix(2, 1);
			lh.fill(0.0);
		}
		Matrix h = lh.clone();
		for (int i = 0; i < h.length(); i++) {
			h.set(i, Math.exp(h.get(i)));
		}

		double lG;
		Matrix lGbasis;

		if (zerothinktimes.size() == R) {
			lGbasis = h.clone();
			for (int i = 0; i < lGbasis.length(); i++) {
				lGbasis.set(i, Math.log(lGbasis.get(i)));
			}
			lG = lG0 + Math.log(h.get(h.length()-1-R));
		} else {
			Matrix scale = new Matrix(1, (int) N_new.elementSum());
			scale.fill(1.0);
			double nt = nvec.elementSum();
			Matrix h_1 = h.clone();
			for (int r = zerothinktimes.size()+1; r <= R; r++) {
				Matrix F1r = null;
				Matrix F2r = null;
				for (int Nr = 1; Nr <= N_new.get(r-1); Nr++) {
					nvec.set(r-1, nvec.get(r-1)+1);
					if (Nr == 1) {
						if (r > zerothinktimes.size()+1) {
							Matrix hr = new Matrix(2*r, 1);
							hr.fill(0.0);
							for (int i = 0; i < r-1; i++) {
								hr.set(i, h.get(i));
							}
							for (int i = r; i < 2*r-1; i++) {
								hr.set(i, h.get(i-1));
							}
							h = hr;
							h.set(r-1, h_1.get(0)/scale.get((int)nt-1));
							h.set(h.length()-1, h_1.get(r-1)/scale.get((int)nt-1));
						}

						Matrix A12 = new Matrix(r, r);
						A12.fill(0.0);
						A12.set(0, 0, -1);
						for (int s = 1; s < r; s++) {
							A12.set(s, 0, N_new.get(s-1));
							A12.set(s, s, -Z_new.get(s-1));
						}

						Matrix B2r = Matrix.eye(r);
						Matrix B2r_tmp = Matrix.eye(r);
						for (int i = 0; i < r; i++) {
							B2r.set(i, i, m*L_new.get(0, r-1));
							B2r_tmp.set(i, i, Z_new.get(r-1));
						}
						B2r = Matrix.concatColumns(B2r, B2r_tmp, null);

						Matrix iC = Matrix.eye(r);
						for (int i = 0; i < r; i++) {
							iC.set(i, i, 1.0/m);
							iC.set(0, i, 1.0/m);
						}
						iC.set(0, 0, -1);

						F1r = new Matrix(2*r, 2*r);
						F1r.fill(0.0);
						F1r.set(0, 0, 1);

						F2r = iC.mult(A12).mult(B2r);
						F2r = Matrix.concatRows(F2r, B2r, null);
					}
					h_1 = h;
					Matrix tmp_mat = F1r.clone();
					for (int i = 0; i < tmp_mat.getNumRows(); i++) {
						for (int j = 0; j < tmp_mat.getNumCols(); j++) {
							tmp_mat.set(i, j, F1r.get(i, j)+F2r.get(i, j)/nvec.get(r-1));
						}
					}
					h = tmp_mat.mult(h_1);
					nt = nvec.elementSum();
					scale.set((int) nt-1, Math.abs(h.elementSum()));
					for (int i = 0; i < h.length(); i++) {
						h.set(i, Math.abs(h.get(i))/scale.get((int) nt-1));
					}
				}
			}
			Matrix log_scale = scale.clone();
			for (int i = 0; i < log_scale.length(); i++) {
				log_scale.set(i, Math.log(log_scale.get(i)));
			}
			lG = lG0 + Math.log(h.get(h.length()-1-(R-1))) + log_scale.elementSum();
			lGbasis = h.clone();
			for (int i = 0; i < lGbasis.length(); i++) {
				lGbasis.set(i, Math.log(h.get(i))+log_scale.elementSum());
			}
		}
		return new pfqnComomrmReturn(lG, lGbasis);
	}

	/** Auxiliary function used by pfqn_gld to computer the normalizing constant in a single-class load-dependent model.*/
	public static pfqnNcReturn pfqn_gldsingle(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
		int M = L.getNumRows();
		int R = L.getNumCols();

		if (R > 1) {
			throw new RuntimeException("pfqn_gldsingle: multiclass model detected. pfqn_gldsingle is for single class models.");
		}

		Map<pfqnGldIndex, Double> g = new HashMap<>();
		g.put(new pfqnGldIndex(1, 1, 1), 0.0);
		for (int n = 1; n <= N.get(0); n++) {
			g.put(new pfqnGldIndex(1,n+1,2), 0.0);
		}

		for (int m = 1; m <= M; m++) {
			for (int tm = 1; tm <= N.get(0)+1; tm++) {
				g.put(new pfqnGldIndex(m+1, 1, tm+1), 1.0);
			}
			for (int n = 1; n <= N.get(0); n++) {
				for (int tm = 1; tm <= N.get(0)-n+1; tm++) {
					g.put(new pfqnGldIndex(m+1,n+1,tm+1),
									g.get(new pfqnGldIndex(m,n+1,2))
													+ L.get(m-1)*g.get(new pfqnGldIndex(m+1,n,tm+2))/mu.get(m-1,tm-1));
				}
			}
		}
		double G = g.get(new pfqnGldIndex(M+1, (int)N.get(0)+1, 2));
		double lG = Math.log(G);
		return new pfqnNcReturn(G, lG);
	}

	/** Auxiliary class used to index interim results in pfqn_gld */
	protected static class pfqnGldIndex {
		int a;
		int b;
		int c;

		pfqnGldIndex(int a, int b, int c) {
			this.a = a;
			this.b = b;
			this.c = c;
		}

		@Override
		public int hashCode() {
			return (a+b+c)*a*b*c;
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof pfqnGldIndex)) {
				return false;
			}
			pfqnGldIndex other = (pfqnGldIndex) obj;
			return (this.a == other.a) && (this.b == other.b)
							&& (this.c == other.c);
		}
	}

	/**
	 * Evaluate limited-load dependent (LLD) function
	 *
	 * @param n Queue-length values. The values can be continuous.
	 * @param lldscaling If not null, then the LLD function uses lldscaling to interpolate continuous values of n  
	 * @param nservers If not null, then the LLD function is set to be for a multi-server with nserver stations
	 * @return Interpolated LLD function values
	 */
	public static Matrix pfqn_lldfun(Matrix n, Matrix lldscaling, Matrix nservers) {
		int M = n.length();
		Matrix r = new Matrix(M, 1);
		r.fill(1.0);
		int smax = lldscaling.getNumCols();
		double alpha = 20;

		for(int i = 0; i < M; i++) {
			if (!(nservers == null || nservers.isEmpty())) {
				if (Double.isInfinite(nservers.get(i))) {
					r.set(i, 0, 1);
				} else {
					double val = r.get(i, 0) / Maths.softmin(n.get(i), nservers.get(i), alpha);
					if (Double.isNaN(val))
						r.set(i, 0, 1.0 / Math.min(n.get(i), nservers.get(i)));
					else
						r.set(i, 0, val);
				}
			}

			if (lldscaling != null && !lldscaling.isEmpty()) {
				Matrix lldscaling_i = new Matrix(1,smax);
				Matrix.extract(lldscaling, i, i+1, 0, smax, lldscaling_i, 0, 0);
				if(lldscaling_i.elementMax() != lldscaling_i.elementMin()) {
					double[] X = new double[smax];
					double[] V = new double[smax];
					for(int j = 0; j < smax; j++) {
						X[j] = j+1;
						V[j] = lldscaling.get(i,j);
					}
					r.set(i, 0, r.get(i,0) / (new SplineInterpolator().interpolate(X, V)).value(n.get(i)));
				}
			}
		}
		return r;
	}

	/**
	 * Evaluate class-dependent (CD) scaling function
	 *
	 * @param nvec Per-class queue-length values. The values can be continuous.
	 * @param cdscaling CD function
	 * @param stations Station objects
	 * @return 
	 */
	public static Matrix pfqn_cdfun(Matrix nvec, Map<Station, SerializableFunction<Matrix, Double>> cdscaling, List<Station> stations) {
		int M = nvec.getNumRows();
		Matrix r = new Matrix(M, 1);
		r.fill(1.0);
		if (!(cdscaling == null || cdscaling.size() == 0)) {
			for(int i = 0; i < M; i++)
				r.set(i, 0, 1.0/cdscaling.get(stations.get(i)).apply(Matrix.extractRows(nvec, i, i+1, null)));
		}
		return r;
	}

	/**
	 * Mean Value Analysis (MVA) Algorithm for closed Product-Form Queueing Networks. Exact solution is computed for
	 * several performance measures.
	 * @param L - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param mi - multiplicity vector
	 * @return - mean value of the computed performance measures, including residence times, throughputs, number
	 * of jobs at a specific node, and utilizations.
	 */
    public static pfqnMVAReturn pfqn_mva(Matrix L, Matrix N, Matrix Z, Matrix mi){
		Matrix XN; // throughputs
		Matrix QN; // queue lengths
		Matrix UN; // utilizations
		Matrix CN; // residence times
		double lGN = 0; // log of the normalizing constant
		double InfServ = 1;
		if(Z == null && mi == null){
			InfServ = 0;
		}
		N = N.ceil();
		int M = L.getNumRows(); // M Stations
		int R = L.getNumCols(); // R Classes
		N = N.columnMajorOrder().transpose();
		if(mi == null){
			mi = Matrix.ones(1, M);
		} else {
			/*
			 * If mi is not passed in, then it is initialised to a row vector. If mi is passed in, it is sometimes a
			 * column vector. In order to deal with this, mi is transposed if it is passed in and it has more than one
			 * row.
			 */
			if(mi.getNumRows() > 1){
				mi = mi.transpose();
			}
		}
		if(Z == null || Z.isEmpty()){
			Z = new Matrix(1, R);
		}
		if(!N.any()){
			return new pfqnMVAReturn(new Matrix(0,0), new Matrix(0,0),
					new Matrix(0,0), new Matrix(0,0), 0);
		}
		int NR = N.length();
		if(R != NR){
			throw new RuntimeException("pfqn_mva: Demand matrix and population vector have different number of classes");
		}

       	XN = new Matrix(1, R);
		QN = new Matrix(M, R);
		CN = new Matrix (M, R);
		if(InfServ == 1){
			Z = Z.columnMajorOrder().transpose();
		} else {
			Z = new Matrix(1, R);
		}
		Matrix prods = new Matrix(1, R - 1);
		for(int w = 0; w < R - 1; w++){
			// ones(1, R-(w+1)+1)
			// (w + 2) instead of (w + 1) because Java indexing starts at 0
			Matrix o = Matrix.ones(1, R - (w + 2) + 1);
			// Addition: ones(1,R-(w+1)+1) + N(1, w+1:R)
			for(int i = 0; i < R - (w + 2) + 1; i++){
				o.set(0, i, o.get(0, i) + N.get(0, w + 1 + i));
			}
			// Now take prod(o)
			prods.set(0, w, o.prodVector());
		}

		int firstNonEmpty = R - 1;
		while(firstNonEmpty >= 0 && N.get(0, firstNonEmpty) == 0){
			firstNonEmpty--;
		}
		double totpop =  Matrix.ones(1, N.getNumCols()).add(1, N).prodVector();
		double ctr = totpop;
		Matrix Q = new Matrix((int) totpop, M);
		int currentpop = 1;

		Matrix n = new Matrix(1, R);
		n.set(0, firstNonEmpty, 1);
		while(ctr > 0){
			int s = 0;
			while(s < R){
				int pos_n_1s = 0;
				if(n.get(0, s) > 0){
					n.set(0,s, n.get(0, s) - 1);
					pos_n_1s = (int) n.get(0, R-1);
					// for w=1:R-1
					int w = 0;
					while(w < R - 1){
						pos_n_1s += n.get(0, w) * prods.get(0, w);
						w++;
					}
					n.set(0,s, n.get(0, s) + 1);
				}
				double CNtot = 0;
				int i = 0;
				// Compute the residence times. Compute the total residence time as well to avoid another iteration
				// through all the stations.
				while(i < M){
					double Lis = L.get(i, s);
					CN.set(i, s, Lis * (mi.get(0, i) + Q.get(pos_n_1s, i)));
					CNtot += CN.get(i, s);
					i++;
				}
				// Compute the throughput for class s
				XN.set(0, s, n.get(0,s) / (Z.get(0, s) + CNtot));
				i = 0;
				// Compute the queue lengths
				while(i < M){
					QN.set(i, s, XN.get(0, s) * CN.get(i, s));
					Q.set(currentpop, i, Q.get(currentpop, i) + QN.get(i,s));
					i++;
				}
				s++;
			}
			s = R - 1;
			while(s >= 0 && (n.get(0, s) == N.get(0, s)) || s > firstNonEmpty){
				s--;
			}
			Matrix nonZero = n.find();
			if(!nonZero.isEmpty()){
				int nonZeroIdx = nonZero.getNumRows() - 1;
				while(nonZeroIdx >= 0 && n.get((int) nonZero.get(nonZeroIdx, 0)) <= 0){
					nonZeroIdx--;
				}
				int last_nnz = (int) nonZero.get(nonZeroIdx, 0);
				double sumn = 0, sumN = 0, sumnprime = 0;
				for(int i = 0; i < last_nnz; i++){
					sumn += n.get(0, i);
					sumN += N.get(0, i);
				}
				for(int i = last_nnz + 1; i < R; i++){
					sumnprime += n.get(0, i);
				}
				if(sumn == sumN && sumnprime == 0){
					double logX = Math.log(XN.get(0, last_nnz));
					lGN -= logX;
				}
			}
			if(s == -1){
				break;
			}
			n.set(0, s, n.get(0, s) + 1);
			s++;
			while(s < R){
				n.set(0, s, 0);
				s++;
			}
			ctr--;
			currentpop++;
		}

		UN = new Matrix(M, R); // Utilizations
		for(int m = 0; m < M; m++){
			for(int r = 0; r < R; r++){
				UN.set(m,r, XN.get(0, r) * L.get(m, r));
			}
		}
		return new pfqnMVAReturn(XN, QN, UN, CN, lGN);
    }

	/**
	 * Computes the upper Geometric Square-Root Bound (GSB) for the throughput of the given closed single-class
	 * queueing networks
	 * @param L - service demand matrix
	 * @param N - population
	 * @param Z - think time
	 * @return - the upper GSB for the throughput
	 */
	public static double pfqn_xzgsbup(Matrix L, double N, double Z){
		int M = L.length();
		double maxL = L.elementMax();
		double R = Z + L.elementSum() + maxL * (N - 1);
		for(int i = 0; i < M; i++){
			if(L.get(i) < maxL){
				R += (L.get(i) - maxL) * pfqn_qzgbup(L, N-1, Z, i);
			}
		}
		return 2 * N / (R + Math.sqrt(R*R - 4 * Z * maxL * N));
	}

	/**
	 * Computes the upper Geometric Bound (GB) for the queue length of the given closed single-class
	 * queueing networks
	 * @param L - service demand matrix
	 * @param N - population
	 * @param Z - think time
	 * @param i - station index
	 * @return - the upper GB for the queue length
	 */
	public static double pfqn_qzgbup(Matrix L, double N, double Z, int i){
		double sumL = L.elementSum();
		double sumLSq = L.elementMult(L, null).elementSum();
		double sigma = sumLSq / sumL;
		double Yi = L.get(i) * Maths.min(1.0/L.elementMax(), N/(Z + sumL + sigma * (N - 1 - Z * pfqn_xzabaup(L, N-1, Z))));
		double Qgb;
		if(Yi < 1){
			Qgb = Yi / (1 - Yi) - (Math.pow(Yi, (N + 1)) / (1 - Yi));
		} else {
			Qgb = N;
		}
		return Qgb;
	}

	/**
	 * Computes the upper ABA for the throughput of the given closed single-class queueing networks
	 * @param L - service demand matrix
	 * @param N - population
	 * @param Z - think time
	 * @return - the upper ABA for the throughput
	 */
	public static double pfqn_xzabaup(Matrix L, double N, double Z){
		double e1 = 1/L.elementMax();
		double e2 = N/(L.elementSum() + Z);
		return Maths.min(e1, e2);
	}

	/**
	 * Computes the lower Geometric Square-Root Bound (GSB) for the throughput of the given closed single-class
	 * queueing networks
	 * @param L - service demand matrix
	 * @param N - population
	 * @param Z - think time
	 * @return - the lower GSB for the throughput
	 */
	public static double pfqn_xzgsblow(Matrix L, double N, double Z){
		int M = L.length();
		double maxL = L.elementMax();
		double R = Z + L.elementSum() + maxL * (N - 1);
		for(int i = 0; i < M; i++){
			if(L.get(i) < maxL){
				R += (L.get(i) - maxL) * pfqn_qzgblow(L, N-1, Z, i);
			}
		}
		return 2 * N / (R + Math.sqrt(R*R - 4 * Z * maxL * (N-1)));
	}

	/**
	 * Computes the lower Geometric Bound (GB) for the queue length of the given closed single-class
	 * queueing networks
	 * @param L - service demand matrix
	 * @param N - population
	 * @param Z - think time
	 * @param i - station index
	 * @return - the lower GB for the queue length
	 */
	public static double pfqn_qzgblow(Matrix L, double N, double Z, int i){
		double Yi = N * L.get(i) / (Z + L.elementSum() + L.elementMax() * N);
		double Qgb = Yi / (1 - Yi) - (Math.pow(Yi, N+1) / (1-Yi));
		return Qgb;
	}

	/**
	 * General purpose script to handle mixed Query Networks with multiserver nodes.
	 * @param lambda - arrival rates
	 * @param L - service demand matrix
	 * @param N - population vector
	 * @param Z - think time matrix
	 * @param mi - multiplicity vector
	 * @param S - server count for each station
	 * @return - performance measures of the mixed query network.
	 */
	public static pfqnMVAReturn pfqn_mvams(Matrix lambda, Matrix L, Matrix N, Matrix Z, Matrix mi, Matrix S){
		int M = L.getNumRows();
		int R = L.getNumCols();
		double Ntot = 0;
		boolean NhasInf = false;
		for(int i = 0; i < N.getNumRows(); i++){
			for(int j = 0; j < N.getNumCols(); j++){
				double num = N.get(i,j);
				if(Double.isFinite(num)){
					Ntot += num;
				} else if (Double.isInfinite(num)){
					NhasInf = true;
				}
			}
		}
		Matrix mu = Matrix.ones(M, (int) Ntot);
		if(S == null){
			S = Matrix.ones(M,1);
		}
		if(mi == null){
			mi = Matrix.ones(M, 1);
		}
		if(Z.isEmpty()){
			Z = new Matrix(1, R);
		}
		for(int i = 0; i < M; i++){
			for(int j = 0; j < Ntot; j++){
				mu.set(i, j, Maths.min(j+1, S.get(i)));
			}
		}
		double maxS = 0;
		boolean initialisedMax = false;
		for(int i = 0; i < S.getNumRows(); i++){
			for(int j = 0; j < S.getNumCols(); j++){
				double num = S.get(i,j);
				if(Double.isFinite(num)){
					if(!initialisedMax || num > maxS){
						initialisedMax = true;
						maxS = num;
					}
				}
			}
		}
		pfqnMVAReturn returnObject;
		if(maxS == 1){
			// No multi-server nodes
			if(NhasInf){
				// open or mixed model
				pfqnMVAReturn retMVAMX = pfqn_mvamx(lambda, L, N, Z, mi);
				returnObject = new pfqnMVAReturn(retMVAMX.XN, retMVAMX.QN, retMVAMX.UN, retMVAMX.CN, retMVAMX.lGN);
			} else {
				// closed model
				pfqnMVAReturn retMVA = pfqn_mva(L,N,Z,mi);
				returnObject = new pfqnMVAReturn(retMVA.XN, retMVA.QN, retMVA.UN, retMVA.CN, retMVA.lGN);
			}
		} else {
			// Multi-server nodes
			if(NhasInf){
				// open or mixed model
				if(mi.elementMax() == 1){
					double lG = Double.NaN;
					pfqnMVAReturn retMVALDMS = pfqn_mvaldms(lambda, L, N, Z, S);
					returnObject = new pfqnMVAReturn(retMVALDMS.XN, retMVALDMS.QN, retMVALDMS.UN, retMVALDMS.CN, lG);
				} else {
					throw new RuntimeException("pfqn_mvams: Queue replicas not available in exact MVA for mixed models.");
				}
			} else {
				pfqnMVALDReturn retMVALD = pfqn_mvald(L,N,Z,mu);
				double lG = retMVALD.lGN.get(retMVALD.lGN.size()-1);
				returnObject = new pfqnMVAReturn(retMVALD.XN, retMVALD.QN, retMVALD.UN, retMVALD.CN, lG);
			}
		}
		return returnObject;
	}

	/**
	 * Mean Value Analysis (MVA) method for open and mixed queueing networks with no multi-server nodes.
	 * @param lambda - arrival rates
	 * @param D - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param mi - multiplicity vector
	 * @return - the performance measures for the given network.
	 */
	public static pfqnMVAReturn pfqn_mvamx(Matrix lambda, Matrix D, Matrix N, Matrix Z, Matrix mi){
		for(int i = 0; i < lambda.getNumCols(); i++){
			if(lambda.get(i) > 0 && N.get(i) > 0 && Double.isFinite(N.get(i))){
				throw new RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.");
			}
		}
		int M = D.getNumRows();
		int R = D.getNumCols();
		if(mi == null){
			mi = Matrix.ones(M,1);
		}
		ArrayList<Integer> openClasses = new ArrayList<>();
		ArrayList<Integer> closedClasses = new ArrayList<>();
		for(int i = 0; i < N.length(); i++){
			if(Double.isInfinite(N.get(i))){
				openClasses.add(i);
			} else {
				closedClasses.add(i);
			}
		}
		Matrix XN = new Matrix(1, R); // Throughputs
		Matrix UN = new Matrix(M,R); // Utilizations
		Matrix CN = new Matrix(M,R); // Residence times
		Matrix QN = new Matrix(M,R); // Queue lengths

		// Compute utilizations and throughputs for the open classes
		for(Integer r : openClasses){
			for(int i = 0; i < M; i++){
				UN.set(i,r,lambda.get(r) * D.get(i,r));
			}
			XN.set(0,r,lambda.get(r));
		}

		Matrix UNt = UN.sumRows();
		if(Z.isEmpty()){
			Z = new Matrix(1,R);
		}
		Matrix Dc = new Matrix(D.getNumRows(), closedClasses.size());
		Matrix rep = UNt.repmat(1, closedClasses.size());
		for(int i = 0; i < Dc.getNumRows(); i++){
			int j = 0;
			for(int closedClass : closedClasses){
				Dc.set(i, j, D.get(i, closedClass) / (1 - rep.get(i,j)));
				j++;
			}
		}
		Matrix Nclosed = new Matrix(1, closedClasses.size());
		Matrix Zclosed = new Matrix(1, closedClasses.size());
		int idx = 0;
		for(int closedClass : closedClasses){
			Nclosed.set(0, idx, N.get(closedClass));
			Zclosed.set(0, idx, Z.get(closedClass));
			idx++;
		}
		// Solve the closed model consisting of M centers and just the closed classes with teh inflated service demands.
		pfqnMVAReturn ret1 = pfqn_mva(Dc, Nclosed, Zclosed, mi);
		for(int i = 0; i < closedClasses.size(); i++){
			XN.set(closedClasses.get(i), ret1.XN.get(i));
			for(int j = 0; j < QN.getNumRows(); j++){
				QN.set(j, closedClasses.get(i), ret1.QN.get(j, i));
			}
			for(int j = 0; j < CN.getNumRows(); j++){
				CN.set(j, closedClasses.get(i), ret1.CN.get(j, i));
			}
		}
		for(int i = 0; i < M; i++){
			for(int r : closedClasses){
				UN.set(i, r, XN.get(r) * D.get(i, r));
			}
		}
		// Compute the residence times and the queue lengths for the open classes.
		for(int i = 0; i < M; i++){
			for(int r : openClasses){
				if(ret1.QN.isEmpty()){
					CN.set(i,r,D.get(i,r) / (1 - UNt.get(i)));
				} else {
					CN.set(i, r, D.get(i, r) * (1 + ret1.QN.sumRows(i)) / (1 - UNt.get(i)));
				}
				QN.set(i, r, CN.get(i, r) * XN.get(r));
			}
		}
		return new pfqnMVAReturn(XN, QN, UN, CN, ret1.lGN);
	}

	/**
	 * Wrapper for pfqn_mvaldmx that adjusts utilizations to account for multiservers
	 * @param lambda - arrival rates
	 * @param D - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param S - servers at each station
	 * @return - the performance measures for the given network.
	 */
	public static pfqnMVAReturn pfqn_mvaldms(Matrix lambda, Matrix D, Matrix N, Matrix Z, Matrix S){
		int M = D.getNumRows();
		int R = D.getNumCols();
		double Nct = 0;
		for(int i = 0; i < N.getNumRows(); i++){
			for(int j = 0; j < N.getNumCols(); j++){
				double num = N.get(i, j);
				if(Double.isFinite(num)){
					Nct += num;
				}
			}
		}
		Matrix mu = Matrix.ones(M, (int) Nct);
		for(int i = 0; i < M; i++){
			for(int j = 0; j < mu.getNumCols(); j++){
				mu.set(i, j, Maths.min(j + 1, S.get(i)));
			}
		}
		if(Z.isEmpty()){
			Z = new Matrix(1, R);
		}
		pfqnMVALDMXReturn ret1 = pfqn_mvaldmx(lambda, D, N, Z, mu, S);
		ArrayList<Integer> openClasses = new ArrayList<>();
		ArrayList<Integer> closedClasses = new ArrayList<>();
		for(int i = 0; i < N.length(); i++){
			if(Double.isInfinite(N.get(i))){
				openClasses.add(i);
			} else {
				closedClasses.add(i);
			}
		}
		Matrix XN = ret1.XN;
		Matrix QN = ret1.QN;
		Matrix CN = ret1.CN;
		double lGN = ret1.lGN;
		Matrix UN = new Matrix(M, R);
		for(int r : closedClasses){
			for(int i = 0; i < M; i++){
				UN.set(i, r, XN.get(r) * D.get(i, r) / S.get(i));
			}
		}
		for(int r : openClasses){
			for(int i = 0; i < M; i++){
				UN.set(i, r, lambda.get(r) * D.get(i,r) / S.get(i));
			}
		}
		return new pfqnMVAReturn(XN, QN, UN, CN, lGN);
	}

	/**
	 * Mean Value Analysis (MVA) method for mixed queueing networks with load-dependent nodes.
	 * @param lambda - arrival rates
	 * @param D - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param mu - load dependent rates
	 * @param S - number of servers at each station
	 * @return - the performance measures for the given network.
	 */
	public static pfqnMVALDMXReturn pfqn_mvaldmx(Matrix lambda, Matrix D, Matrix N, Matrix Z, Matrix mu, Matrix S){
		/*
		* Parameter S is not used at all
		* */
		double NfiniteSum = 0;
		for(int i = 0; i < N.getNumRows(); i++){
			for(int j = 0; j < N.getNumCols(); j++){
				double num = N.get(i,j);
				if(Double.isFinite(num)){
					NfiniteSum += num;
				}
			}
		}
		if(mu == null && S == null){
			mu = Matrix.ones(D.getNumRows(), (int) NfiniteSum);
			S = Matrix.ones(D.getNumRows(), 1);
		}
		if(mu.getNumCols() < NfiniteSum){
			throw new RuntimeException("pfqn_mvaldmx: MVALDMX requires to specify the load-dependent rates with one job more than the maximum closed population.");
		}
		Matrix f = lambda.find();
		for(int i = 0; i < f.getNumRows(); i++){
			double num = N.get((int) f.get(i));
			if(num > 0 && Double.isFinite(num)){
				throw new RuntimeException("pfqn_mvaldmx: Arrival rate cannot be specified on closed classes.");
			}
		}
		int M = D.getNumRows();
		int R = D.getNumCols();
		ArrayList<Integer> openClasses = new ArrayList<>();
		ArrayList<Integer> closedClasses = new ArrayList<>();
		for(int i = 0; i < N.length(); i++){
			if(Double.isInfinite(N.get(i))){
				openClasses.add(i);
			} else {
				closedClasses.add(i);
			}
		}
		Matrix XN = new Matrix(1, R);
		Matrix UN = new Matrix(M, R);
		Matrix CN = new Matrix(M, R);
		Matrix QN = new Matrix(M, R);
		double lGN = 0;
		Matrix newMu = new Matrix(mu.getNumRows(), mu.getNumCols() + 1);
		for(int i = 0; i < newMu.getNumRows(); i++){
			for(int j = 0; j < newMu.getNumCols(); j++){
				if(j < mu.getNumCols()){
					newMu.set(i, j, mu.get(i, j));
				} else {
					newMu.set(i, j, mu.get(i, j - 1));
				}
			}
		}
		mu = newMu; // we need up to sum(N)+1, but there is limited load dep
		pfqnMVALDMXECReturn ret1 = pfqn_mvaldmx_ec(lambda, D, new Matrix(mu));
		Matrix EC = ret1.EC;
		Matrix E = ret1.E;
		Matrix Eprime = ret1.Eprime;
		int C = closedClasses.size(); // number of closed classes
		Matrix Dc = new Matrix(D.getNumRows(), C);
		Matrix Nc = new Matrix(1, C);
		Matrix Zc = new Matrix(1, C);
		for(int i = 0; i < C; i++){
			int c = closedClasses.get(i);
			for(int j = 0; j < Dc.getNumRows(); j++){
				Dc.set(j, i, D.get(j, c));
			}
			Nc.set(0, i, N.get(c));
			Zc.set(0, i, Z.get(c));
		}
		Matrix prods = new Matrix(1, C); // needed for fast hashing
		for(int r = 0; r < C; r++){
			double prod = 1;
			for(int i = 0; i < r; i++){
				prod *= (Nc.get(i) + 1);
			}
			prods.set(0, r, prod);
		}
		// Start at nc=(0,...,0)
		Matrix nvec = pprod(Nc);
		// Initialize Pc
		Matrix[] Pc = new Matrix[M];
		double ncProd = 1;
		for(int i = 0; i < C; i++){
			ncProd *= (1 + Nc.get(i));
		}
		for(int i = 0; i < M; i++){
			Pc[i] = new Matrix((int) (1 + Nc.elementSum()), (int) ncProd);
		}
		Matrix x = new Matrix(C, (int) ncProd);
		Matrix[] w = new Matrix[M];
		for(int i = 0; i < M; i++){
			w[i] = new Matrix(C, (int) ncProd);
		}
		for(int i = 0; i < M; i++){
			Pc[i].set(0, hashpop(nvec, Nc, C, prods), 1);
		}
		Matrix u = new Matrix(M, C);
		// Population recursion
		while(!(nvec.getNumRows() == 1 && nvec.getNumCols() == 1 && nvec.get(0,0) == -1)){
			int hnvec = hashpop(nvec, Nc, C, prods);
			double nc = nvec.elementSum();
			for(int i = 0; i < M; i++){
				for(int c = 0; c < C; c++){
					if(nvec.get(c) > 0){
						int hnvec_c = hashpop(oner(nvec,new ArrayList<>(Collections.singletonList(c))), Nc, C, prods);
						// Compute mean residence times
						for(int n = 0; n < nc; n++){
							w[i].set(c, hnvec, w[i].get(c, hnvec) +
								Dc.get(i, c) * (n+1) * EC.get(i, n) * Pc[i].get(n, hnvec_c));
						}
					}
				}
			}
			// Compute tput
			for(int c = 0; c < C; c++){
				double sumw = 0;
				for(int i = 0; i < M; i++){
					sumw += w[i].get(c, hnvec);
				}
				x.set(c, hnvec, nvec.get(c) / (Zc.get(c) + sumw));
			}
			for(int i = 0; i < M; i++){
				for(int n = 0; n < nc; n++){
					for(int c = 0; c < C; c++){
						if(nvec.get(c) > 0){
							int hnvec_c = hashpop(oner(nvec, new ArrayList<>(Collections.singletonList(c))), Nc, C, prods);
							Pc[i].set(1 + n, hnvec, Pc[i].get(1 + n, hnvec) +
								Dc.get(i, c) * EC.get(i, n) * x.get(c, hnvec) * Pc[i].get(n, hnvec_c));
						}
					}
				}
				double sumpc = 0;
				for(int k = 0; k < nc; k++){
					sumpc += Pc[i].get(1 + k, hnvec);
				}
				Pc[i].set(0, hnvec, Maths.max(Math.ulp(1.0), 1 - sumpc));
			}
			Matrix nvecFind = nvec.find();
			if(!nvecFind.isEmpty()){
				int idx = nvecFind.getNumRows() - 1;
				while(idx >= 0 && nvec.get((int) nvecFind.get(idx)) <= 0){
					idx--;
				}

				// now compute the normalizing constant
				int last_nnz = (int) nvecFind.get(idx);
				double sumnvec, sumnc, sumnvecp;
				sumnvec = sumnc = sumnvecp = 0;
				for(int i = 0; i < C; i++){
					if(i < last_nnz){
						sumnvec += nvec.get(i);
						sumnc += Nc.get(i);
					} else if (i > last_nnz){
						sumnvecp += nvec.get(i);
					}
				}
				if(sumnvec == sumnc && sumnvecp == 0){
					double logX = Math.log(XN.get(last_nnz));
					lGN -= logX;
				}
			}
			nvec = pprod(nvec, Nc);
		}
		// compute performance indexes at Nc for closed classes
		int hnvec = hashpop(Nc, Nc, C, prods);
		for(int c = 0; c < C; c++){
			int hnvec_c = hashpop(oner(Nc, new ArrayList<>(Collections.singletonList(c))), Nc, C, prods);
			for(int i = 0; i < M; i++){
				u.set(i, c, 0);
				double sumNc = Nc.elementSum();
				for(int n = 0; n < sumNc; n++){
					u.set(i, c, u.get(i, c) + Dc.get(i, c) * x.get(c, hnvec)
						* Eprime.get(i, n) / E.get(i, n) * Pc[i].get(n, hnvec_c));
				}
			}
		}
		for(int i = 0; i < C; i++){
			// Throughput
			XN.set(closedClasses.get(i), x.get(i, hnvec));
			for(int j = 0; j < M; j++){
				// Utilization
				UN.set(j, closedClasses.get(i), u.get(j, i));
				// Response time
				CN.set(j, closedClasses.get(i), w[j].get(i, hnvec));
				// Queue length
				QN.set(j, closedClasses.get(i), XN.get(closedClasses.get(i)) * CN.get(j, closedClasses.get(i)));
			}
		}
		// Compute performance indexes at Nc for open classes
		for(int r : openClasses){
			XN.set(r, lambda.get(r));
			for(int i = 0; i < M; i++){
				// Queue-length
				QN.set(i, r, 0);
				double sumNc = Nc.elementSum();
				for(int n = 0; n <= sumNc; n++){
					QN.set(i, r, QN.get(i, r) + lambda.get(r) * D.get(i,r) * (n+1) *
							EC.get(i, n) * Pc[i].get(n, hnvec));
				}
				// Response time
				CN.set(i, r, QN.get(i, r) / lambda.get(r));
				// Utilization - the formula from Bruell-Balbo-Ashfari does not
				// match simulation, this appears to be simply lambda_r*D_{ir}
				UN.set(i, r, 0);
				for(int n = -1; n < sumNc; n++){
					UN.set(i, r, UN.get(i, r) + lambda.get(r) * Eprime.get(i, n+2) / E.get(i, n+2)
						* Pc[i].get(n + 1, hnvec));
				}
			}
		}
		Matrix newPc = new Matrix(M, (int) (1 + Nc.elementSum()));
		for(int i = 0; i < M; i++){
			for(int j = 0; j < newPc.getNumCols(); j++){
				newPc.set(i, j, Pc[i].get(j, hnvec));
			}
		}
		return new pfqnMVALDMXReturn(XN, QN, UN, CN, lGN, newPc);
	}

	/** Auxiliary function used by pfqn_mvaldmx to compute the EC terms */
	protected static pfqnMVALDMXECReturn pfqn_mvaldmx_ec(Matrix lambda, Matrix D, Matrix mu){
		int M = mu.getNumRows();
		int R = mu.getNumCols();
		Matrix Lo = new Matrix(M, 1);
		for(int i = 0; i < M; i++){
			Lo.set(i, lambda.mult(Matrix.extractRows(D, i, i+1, null).transpose()).get(0));
		}
		Matrix b = new Matrix(M, 1); // Limited load dependence level
		for(int i = 0; i < M; i++){
			int idx = 0;
			while(idx < mu.getNumCols() && mu.get(i, idx) != mu.get(i, mu.getNumCols() - 1)){
				idx++;
			}
			b.set(i, idx);
		}
		int Nt = mu.getNumCols(); // Compute extra elements if present
		int oldEnd = mu.getNumCols() - 1;
		mu.expandMatrix(mu.getNumRows(), mu.getNumCols() + 2 + (int) b.elementMax(), mu.getNonZeroLength() +
				(2 + (int) b.elementMax()) * Matrix.extractColumn(mu, oldEnd, null).getNonZeroLength());
		for(int i = 0; i < mu.getNumRows(); i++){
			for(int j = oldEnd + 1; j < mu.getNumCols(); j++){
				mu.set(i, j, mu.get(i, oldEnd));
			}
		}
		Matrix C = new Matrix(mu.getNumRows(), mu.getNumCols());
		mu.divide(1, C, false);
		Matrix EC = new Matrix(M, Nt);
		Matrix E = new Matrix(M, 1 + Nt);
		Matrix Eprime = new Matrix(M, 1 + Nt);
		for(int i = 0; i < M; i++){
			Matrix E1 = new Matrix(1 + Nt, 1 + Nt);
			Matrix E2 = new Matrix(1 + Nt, 1 + Nt);
			Matrix E3 = new Matrix(1 + Nt, 1 + Nt);
			Matrix F2 = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);
			Matrix F3 = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);

			Matrix E2prime = new Matrix(1 + Nt, 1 + Nt);
			Matrix F2prime = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);
			for(int n = 0; n <= Nt; n++){
				if(n >= b.get(i) + 1){
					E.set(i, n, 1 / Math.pow(1 - Lo.get(i) * C.get(i, (int) b.get(i)), n+1));
					Eprime.set(i, n, C.get(i, (int) b.get(i)) * E.get(i, n));
				} else {
					// Compute E1
					if(n == 0){
						E1.set(n, 1/(1 - Lo.get(i) * C.get(i, (int) b.get(i))));
						for(int j = 0; j < b.get(i) - 1 + 1; j++){
							E1.set(n, E1.get(n) * C.get(i, j) / C.get(i, (int) b.get(i)));
						}
					} else {
						E1.set(n, 1 / (1 - Lo.get(i) * C.get(i, (int) b.get(i))) *
								C.get(i, (int) b.get(i)) / C.get(i, n - 1) * E1.get(n - 1));
					}

					// Compute F2
					for(int n0 = 0; n0 <= b.get(i) - 2 + 1; n0++){
						if(n0 == 0){
							F2.set(n, n0,1);
						} else {
							F2.set(n, n0, ((double) n+n0)/n0 * Lo.get(i) * C.get(i, n + n0 - 1)
									* F2.get(n, n0 - 1));
						}
					}

					// Compute E2
					double sumf2 = 0;
					for(int k = -1; k < b.get(i) - 2 + 1; k++){
						sumf2 += F2.get(n, 1 + k);
					}
					E2.set(n, sumf2);

					// Compute F3
					for(int n0 = 0; n0 <= b.get(i) - 2 + 1; n0++){
						if(n == 0 && n0 == 0){
							F3.set(n, n0, 1);
							for(int j = 0; j < b.get(i) - 1 + 1; j++){
								F3.set(n, n0, F3.get(n,n0) * C.get(i, j) /
										C.get(i, (int) b.get(i)));
							}
						} else if (n > 0 && n0 == 0){
							F3.set(n, n0, C.get(i, (int) b.get(i)) / C.get(i, n - 1) *
									F3.get(n-1, 0));
						} else {
							F3.set(n, n0, ((double) n+n0)/n0 * Lo.get(i) *
									C.get(i, (int) b.get(i)) * F3.get(n, n0-1));
						}
					}

					// Compute E3
					double sumf3 = 0;
					for(int k = -1; k < b.get(i) - 2 + 1; k++){
						sumf3 += F3.get(n, 1 + k);
					}
					E3.set(n, sumf3);

					// Compute F2prime
					for(int n0 = 0; n0 <= b.get(i) - 2 + 1; n0++){
						if(n0 == 0){
							F2prime.set(n, n0, C.get(i, n));
						} else {
							F2prime.set(n,n0, ((double) n+n0)/n0 * Lo.get(i) * C.get(i, n+n0) *
									F2prime.get(n, n0-1));
						}
					}

					// Compute E2prime
					double sumf2p = 0;
					for(int k = -1; k < b.get(i) - 2 + 1; k++){
						sumf2p += F2prime.get(n, 1+k);
					}
					E2prime.set(n, sumf2p);

					// Compute E, Eprime and EC
					E.set(i, n, E1.get(n) + E2.get(n) - E3.get(n));
					if(n < b.get(i) - 1 + 1){
						Eprime.set(i, n, C.get(i, (int) b.get(i)) * E1.get(n) +
									E2prime.get(n) - C.get(i, (int) b.get(i)) * E3.get(n));
					} else {
						Eprime.set(i, n, C.get(i, (int) b.get(i)) * E.get(i, n));
					}
				}
			}
			for(int n = 0; n < Nt; n++){
				EC.set(i, n, C.get(i, n) * E.get(i, n + 1) / E.get(i, n));
			}
		}
		return new pfqnMVALDMXECReturn(EC, E, Eprime, Lo);
	}

	/**
	 * Mean Value Analysis (MVA) method for closed networks with load dependent service and stabilization
	 */
	public static pfqnMVALDReturn pfqn_mvald(Matrix L, Matrix N, Matrix Z, Matrix mu){
		// Stabilize ensures that probabilities do not become negative
		return pfqn_mvald(L, N, Z, mu, true);
	}

	/**
	 * Mean Value Analysis (MVA) method for closed networks with load dependent service
	 * @param L - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param mu - load dependent service rate. mu[i][j] - load dependent service rate of the ith node when there are
	 *           j jobs in it
	 * @param stabilize - whether to stabilize the probabilities or not (ensures that probabilities do not become negative)
	 * @return - performance measures for the closed network.
	 */
	public static pfqnMVALDReturn pfqn_mvald(Matrix L, Matrix N, Matrix Z, Matrix mu, boolean stabilize){
		boolean warn = true;
		boolean isNumStable = true;
		int M = L.getNumRows(); // number of queues
		int R = L.getNumCols(); // number of classes
		double prodN = 1;
		for(int i = 0; i < N.getNumRows(); i++){
			for(int j = 0; j < N.getNumCols(); j++){
				prodN *= (1 + N.get(i, j));
			}
		}
		Matrix Xs = new Matrix(R, (int) prodN); // throughput for a model with station i less
		Matrix[] pi = new Matrix[M]; // marginal queue-length probabilities pi(k)
		for(int i = 0; i < M; i++){
			pi[i] = Matrix.ones((int) N.elementSum() + 1, (int) prodN);
		}
		Matrix WN = new Matrix(M, R);
		Matrix n = pprod(N);
		List<Double> lGN = new ArrayList<>();
		lGN.add(0.0);
		while(!(n.getNumRows() == 1 && n.getNumCols() == 1 && n.get(0,0) == -1)){
			WN = new Matrix(WN.getNumRows(), WN.getNumCols());
			for(int s = 0; s < R; s++){
				if(n.get(s) > 0){
					for(int i = 0; i < M; i++){
						WN.set(i, s, 0);
						for(int k = 0; k < n.elementSum(); k++){
							WN.set(i, s, WN.get(i, s) + (L.get(i, s)/mu.get(i,k))
								*(k+1)*pi[i].get(k, hashpop(oner(n, new ArrayList<>(Collections.singletonList(s))), N)));
						}
					}
					Xs.set(s, hashpop(n, N), n.get(s)/(Z.get(s) +
							Matrix.extractColumn(WN, s, null).elementSum()));
				}
			}
			// Compute pi(k|n)
			for(int k = 0; k < n.elementSum(); k++){
				for(int i = 0; i < M; i++){
					pi[i].set(k + 1, hashpop(n, N), 0);
				}
				for(int s = 0; s < R; s++){
					if(n.get(s) > 0){
						for(int i = 0; i < M; i++){
							pi[i].set(k+1, hashpop(n, N), pi[i].get(k+1, hashpop(n,N))
									+ (L.get(i, s) / mu.get(i, k)) * Xs.get(s, hashpop(n, N))
									* pi[i].get(k, hashpop(oner(n, new ArrayList<>(Collections.singletonList(s))), N)));
						}
					}
				}
			}
			// Compute pi(0|n)
			for(int i = 0; i < M; i++){
				double sumpi = 0;
				for(int k = 0; k < n.elementSum(); k++){
					sumpi += pi[i].get(k + 1, hashpop(n, N));
				}
				double p0 = 1 - sumpi;
				if(p0 < 0){
					if(warn){
						System.err.println("pfqn_mvald: MVA-LD is numerically unstable on this model, " +
								"forcing all probabilities to be non-negative.");
						warn = false;
						isNumStable = false;
					}
					if(stabilize){
						pi[i].set(0, hashpop(n, N), Math.ulp(1.0));
					} else {
						pi[i].set(0, hashpop(n, N), p0);
					}
				} else {
					pi[i].set(0, hashpop(n, N), p0);
				}
			}
			Matrix nfind = n.find();
			if(!nfind.isEmpty()){
				int idx = nfind.getNumRows() - 1;
				while(idx >= 0 && n.get((int) nfind.get(idx)) < 0){
					idx--;
				}
				int last_nnz = (int) nfind.get(idx);
				double sumn = 0, sumN = 0, sumnp = 0;
				for(int i = 0; i < R; i++){
					if(i < last_nnz){
						sumn += n.get(i);
						sumN += N.get(i);
					} else if(i > last_nnz){
						sumnp += n.get(i);
					}
				}
				if(sumn == sumN && sumnp == 0){
					double logX = Math.log(Xs.get(last_nnz, hashpop(n, N)));
					lGN.add(lGN.get(lGN.size()-1) - logX);
				}
			}
			n = pprod(n, N); // get the next population
		}
		Matrix XN = Matrix.extractColumn(Xs, hashpop(N,N), null).transpose();
		Matrix newpi = new Matrix(M, (int) N.elementSum() + 1);
		for(int i = 0; i < M; i++){
			for(int j = 0; j < (int) N.elementSum() + 1; j++){
				newpi.set(i, j, pi[i].get(j, hashpop(N, N)));
			}
		}
		Matrix QN;
		if(WN.isEmpty()){
			QN = new Matrix(0,0);
		} else {
			QN = WN.elementMult(XN.repmat(M, 1), null);
		}
		Matrix UN = Matrix.ones(newpi.getNumRows(), 1).sub(1,
					Matrix.extractColumn(newpi, 0, null));
		Matrix CN = new Matrix(N.getNumRows(), N.getNumCols());
		for(int i = 0; i < CN.getNumRows(); i++){
			for(int j = 0; j < CN.getNumCols(); j++){
				CN.set(i, j, N.get(i, j) / XN.get(i, j) - Z.get(i, j));
			}
		}
		return new pfqnMVALDReturn(XN, QN, UN, CN, lGN, isNumStable, newpi);
	}

	/** Multiserver version of Krzesinski's Linearizer */
	public static pfqnLinearizerMSReturn pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers){
		int M = L.getNumRows();
		Matrix type = new Matrix(M, 1);
		type.fill(SchedStrategy.toID(SchedStrategy.PS));
		return pfqn_linearizerms(L, N, Z, nservers, type);
	}

	/** Multiserver version of Krzesinski's Linearizer */
	public static pfqnLinearizerMSReturn pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers, Matrix type){
		return pfqn_linearizerms(L, N, Z, nservers, type, GlobalConstants.FineTol);
	}

	/** Multiserver version of Krzesinski's Linearizer */
	public static pfqnLinearizerMSReturn pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers, Matrix type, double tol){
		return pfqn_linearizerms(L, N, Z, nservers, type, tol, 1000);
	}

	/**
	 *  Multiserver version of Krzesinski's Linearizer as described in Conway 1989, Fast Approximate Solution of
	 *  Queueing Networks with Multi-Server Chain- Dependent FCFS Queues. Minor adjustments based on De Souza-Muntz's
	 *  description of the algorithm.
	 * @param L - service demand matrix
	 * @param N - population vector
	 * @param Z - think times
	 * @param nservers - number of servers at each station
	 * @param type - scheduling discipline at each station
	 * @param tol - max tolerance admitted between successive iterations
	 * @param maxiter - maximum number of iterations
	 * @return - the performance measures of the network.
	 */
	public static pfqnLinearizerMSReturn pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers, Matrix type, double tol, int maxiter){
		int M = L.getNumRows();
		int R = L.getNumCols();

		if(Z.isEmpty()){
			Z = new Matrix(1, R);
		}

		// Initialise
		Matrix[] Q = new Matrix[M];
		Matrix PB = new Matrix(M, 1 + R);
		Matrix[] P = new Matrix[M];
		Matrix[] Delta = new Matrix[M];
		for(int i = 0; i < M; i++){
			Q[i] = new Matrix(R, 1 + R);
			P[i] = new Matrix((int) nservers.elementMax(), 1 + R);
			Delta[i] = new Matrix(R, R);
		}

		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				for(int s = 0; s < R; s++){
					Delta[i].set(r, s, 0);
				}
				for(int s = -1; s < R; s++){
					Matrix N_1 = oner(N, new ArrayList<>(Collections.singletonList(s)));
					Q[i].set(r, 1 + s, N_1.get(r) / M);
				}
			}
		}

		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				for(int s = -1; s < R; s++){
					Matrix N_1 = oner(N, new ArrayList<>(Collections.singletonList(s)));
					double pop = N_1.elementSum();
					if(nservers.get(i) > 1){
						double sumQ = 0;
						for(int k = 0; k < R; k++){
							sumQ += Q[i].get(k, 1 + s);
						}
						for(int j = 0; j < nservers.get(i) - 1; j++){
							P[i].set(1+j, 1+s, 2 * sumQ/(pop * (pop+1)));
						}
						PB.set(i, 1 + s, 2 * sumQ / (pop + 1 - nservers.get(i)) / (pop * (pop+1)));
						double sumP = 0;
						for(int k = 0; k < nservers.get(i) - 1; k++){
							sumP += P[i].get(k + 1, 1 + s);
						}
						P[i].set(0, 1+s, 1 - PB.get(i, 1 + s) - sumP);
					}
				}
			}
		}

		int totiter = 0;

		// Main loop
		for(int I = 0; I < 2; I++){
			for(int s = -1; s < R; s++){
				Matrix N_1 = oner(N, new ArrayList<>(Collections.singletonList(s)));
				Matrix Q1 = new Matrix(M, R);
				Matrix P1 = new Matrix(M, (int) nservers.elementMax());
				Matrix PB1 = new Matrix(M,1);
				for(int i = 0; i < M; i++){
					for(int j = 0; j < R; j++){
						Q1.set(i, j, Q[i].get(j, 1+s));
					}
					for(int j = 0; j < P1.getNumCols(); j++){
						P1.set(i, j, P[i].get(j, 1+s));
					}
					PB1.set(i, 0, PB.get(i, 1+s));
				}
				pfqnLinearizerMSCoreReturn ret1 = linearizerms_core(L, M, R, N_1, Z, nservers, Q1,
						P1, PB1, Delta, type, tol, maxiter - totiter);
				for(int i = 0; i < M; i++){
					for(int j = 0; j < R; j++){
						Q[i].set(j, 1 + s, ret1.Q.get(i, j));
					}
					for(int j = 0; j < nservers.elementMax(); j++){
						P[i].set(j, 1 + s, ret1.P.get(i, j));
					}
					PB.set(i, 1+s, ret1.PB.get(i));
				}
				totiter += ret1.iter;
			}
			// Upgrade delta
			for(int i = 0; i < M; i++){
				for(int r = 0; r < R; r++){
					for(int s = 0; s < R; s++){
						Matrix Ns = oner(N, new ArrayList<>(Collections.singletonList(s)));
						Delta[i].set(r, s, Q[i].get(r, 1 + s)/Ns.get(r) - Q[i].get(r, 0)/N.get(r));
					}
				}
			}
		}

		Matrix Q1 = new Matrix(M, R);
		Matrix P1 = new Matrix(M, (int) nservers.elementMax());
		Matrix PB1 = new Matrix(M,1);
		for(int i = 0; i < M; i++){
			for(int j = 0; j < R; j++){
				Q1.set(i, j, Q[i].get(j, 0));
			}
			for(int j = 0; j < P1.getNumCols(); j++){
				P1.set(i, j, P[i].get(j, 0));
			}
			PB1.set(i, 0, PB.get(i, 0));
		}
		pfqnLinearizerMSCoreReturn ret1 = linearizerms_core(L, M, R, N, Z, nservers, Q1, P1,
												PB1, Delta, type, tol, maxiter - totiter);
		totiter += ret1.iter;
		Matrix newQ = ret1.Q;
		Matrix W = ret1.W;
		Matrix X = ret1.T;
		Matrix U = new Matrix(M, R);

		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				if(nservers.get(i) == 1){
					U.set(i, r, X.get(r) * L.get(i, r));
				} else {
					U.set(i, r, X.get(r) * L.get(i, r) / nservers.get(i));
				}
			}
		}
		Matrix C = new Matrix(1, R);
		for(int i = 0; i < R; i++){
			C.set(0, i, N.get(i) / X.get(i) - Z.get(i));
		}
		return new pfqnLinearizerMSReturn(newQ, U, W, C, X, totiter);
	}

	protected static pfqnLinearizerMSCoreReturn linearizerms_core(Matrix L, int M, int R, Matrix N_1, Matrix Z, Matrix nservers,
                                                                Matrix Q, Matrix P, Matrix PB, Matrix[] Delta, Matrix type,
                                                                double tol, int maxiter){
		int iter = 0;
		boolean hasConverged = false;
		Matrix W = null, T = null;
		while(!hasConverged){
			iter++;
			Matrix Qlast = new Matrix(Q);
			// Estimate population
			pfqnLinearizerMSEstimateReturn ret1 = linearizerms_estimate(M, R, N_1, nservers, Q, P, PB, Delta);
			// Forward MVA
			pfqnLinearizerMSForwardMVAReturn ret2 = linearizerms_forwardMVA(L, M, R, N_1, Z, nservers, type, ret1.Q_1, ret1.P_1, ret1.PB_1);
			Q = ret2.Q;
			W = ret2.W;
			T = ret2.T;
			P = ret2.P;
			PB = ret2.PB;
			if(Q.sub(1, Qlast).norm() < tol || iter > maxiter){
				hasConverged = true;
			}
		}
		return new pfqnLinearizerMSCoreReturn(Q, W, T, P, PB, iter);
	}

	protected static pfqnLinearizerMSEstimateReturn linearizerms_estimate(int M, int R, Matrix N_1, Matrix nservers, Matrix Q, Matrix P, Matrix PB, Matrix[] Delta) {
		Matrix[] P_1 = new Matrix[M];
		Matrix[] Q_1 = new Matrix[M];
		for(int i = 0; i < M; i++){
			P_1[i] = new Matrix((int) nservers.elementMax(), 1 + R);
			Q_1[i] = new Matrix(R, 1 + R);
		}
		Matrix PB_1 = new Matrix(M, 1 + R);
		for(int i = 0; i < M; i++){
			if(nservers.get(i) > 1){
				for(int j = -1; j < nservers.get(i) - 1; j++){
					for(int s = -1; s < R; s++){
						P_1[i].set(1+j, 1+s, P.get(i, 1 + j));
					}
				}
				for(int s = -1; s < R; s++){
					PB_1.set(i, 1 + s, PB.get(i, 0));
				}
			}
			for(int r = 0; r < R; r++){
				for(int s = 0; s < R; s++){
					Matrix Ns = oner(N_1, new ArrayList<>(Collections.singletonList(s)));
					Q_1[i].set(r, 1 + s, Ns.get(r)*(Q.get(i, r)/N_1.get(r) + Delta[i].get(r,s)));
				}
			}
		}
		return new pfqnLinearizerMSEstimateReturn(Q_1, P_1, PB_1);
	}

	protected static pfqnLinearizerMSForwardMVAReturn linearizerms_forwardMVA(Matrix L, int M, int R, Matrix N_1, Matrix Z,
                                                                           Matrix nservers, Matrix type, Matrix[] Q_1,
                                                                           Matrix[] P_1, Matrix PB_1) {
		Matrix W = new Matrix(M, R);
		Matrix T = new Matrix(1, R);
		Matrix Q = new Matrix(M, R);
		Matrix P = new Matrix(M, (int) nservers.elementMax());
		Matrix PB = new Matrix(M, 1);
		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				W.set(i, r, L.get(i, r) / nservers.get(i));
				if(L.get(i, r) == 0){
					// 0 service demand at this station => this class does not visit the current node
					continue;
				}
				boolean flag = true; // flag = whether all nodes have the FCFS scheduling strategy
				for(int k = 0; k < M; k++){
					if(type.get(k) != SchedStrategy.toID(SchedStrategy.FCFS)){
						flag = false;
					}
				}
				if(flag){
					for(int s = 0; s < R; s++){
						W.set(i, r, W.get(i, r) + (L.get(i, s)/nservers.get(i)) * Q_1[i].get(s, 1 + r));
					}
				} else {
					for(int s = 0; s < R; s++){
						W.set(i, r, W.get(i, r) + (L.get(i, r)/nservers.get(i)) * Q_1[i].get(s, 1 + r));
					}
				}
				if(nservers.get(i) > 1){
					for(int j = 0; j <= nservers.get(i) - 2; j++){
						if(flag){
							for(int s = 0; s < R; s++){
								W.set(i, r, W.get(i, r) + L.get(i, s) * (nservers.get(i) - 1 - j) * P_1[i].get(j, 1 + r));
							}
						} else {
							for(int s = 0; s < R; s++){
								W.set(i, r, W.get(i, r) + L.get(i, r) * (nservers.get(i) - 1 - j) * P_1[i].get(j, 1 + r));
							}
						}
					}
				}
			}
		}
		for(int r = 0; r < R; r++){
			T.set(r, N_1.get(r) / (Z.get(r) + Matrix.extractColumn(W, r, null).elementSum()));
			for(int i = 0; i < M; i++){
				Q.set(i, r, T.get(r) * W.get(i, r));
			}
		}
		for(int i = 0; i < M; i++){
			if(nservers.get(i) > 1){
				for(int k = 0; k < P.getNumCols(); k++){
					P.set(i, k, 0);
				}
				for(int j = 1; j <= nservers.get(i) - 1; j++){
					for(int s = 0; s < R; s++){
						P.set(i, j, P.get(i, j) + L.get(i, s) * T.get(s) * P_1[i].get(j-1, 1+s) / j);
					}
				}
			}
		}
		for(int i = 0; i < M; i++){
			if(nservers.get(i) > 1){
				PB.set(i, 0, 0);
				for(int s = 0; s < R; s++){
					PB.set(i, 0, PB.get(i, 0) + L.get(i, s) * T.get(s) *
							(PB_1.get(i, 1 + s) + P_1[i].get((int) nservers.get(i) - 1, 1 + s))/nservers.get(i));
				}
			}
		}
		for(int i = 0; i < M; i++){
			if(nservers.get(i) > 1){
				P.set(i, 0, 1 - PB.get(i));
				for(int j = 0; j < nservers.get(i) - 1; j++){
					P.set(i, 0, P.get(i, 0) - P.get(i, 1 + j));
				}
			}
		}
		return new pfqnLinearizerMSForwardMVAReturn(Q, W, T, P, PB);
	}

	/** Bard-Schweitzer approximate mean value analysis algorithm */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N){
		Matrix Z = new Matrix(N.getNumRows(), N.getNumCols());
		return pfqn_bs(L, N, Z);
	}

	/** Bard-Schweitzer approximate mean value analysis algorithm */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N, Matrix Z){
		return pfqn_bs(L, N, Z, 1.0e-6);
	}

	/** Bard-Schweitzer approximate mean value analysis algorithm */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol){
		return pfqn_bs(L, N, Z, tol, 1000);
	}

	/** Bard-Schweitzer approximate mean value analysis algorithm */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter){
		return pfqn_bs(L, N, Z, tol, maxiter, null);
	}

	/** Bard-Schweitzer approximate mean value analysis algorithm */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0){
		int M = L.getNumRows();
		if(QN0 == null || QN0.isEmpty()){
			QN0 = N.repmat(M, 1);
			for(int i = 0; i < QN0.getNumRows(); i++){
				for(int j = 0; j < QN0.getNumCols(); j++){
					QN0.set(i, j, QN0.get(i, j) / M);
				}
			}
		}
		SchedStrategy[] type = new SchedStrategy[M];
		Arrays.fill(type, SchedStrategy.PS);
		return pfqn_bs(L, N, Z, tol, maxiter, QN0, type);
	}

	/**
	 * Bard-Schweitzer approximate mean value analysis algorithm
	 * @param L - the service demand matrix
	 * @param N - the population vector
	 * @param Z - the think times vector
	 * @param tol - max tolerance admitted between successive iterations
	 * @param maxiter - maximum number of iterations
	 * @param QN0 - original queue lengths
	 * @param type - scheduling disciplines at each station
	 * @return - the performance metrics for this network.
	 */
	public static pfqnBSReturn pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0, SchedStrategy[] type){
		int M = L.getNumRows();
		int R = L.getNumCols();
		if(QN0 == null || QN0.isEmpty()){
			QN0 = N.repmat(M, 1);
			for(int i = 0; i < QN0.getNumRows(); i++){
				for(int j = 0; j < QN0.getNumCols(); j++){
					QN0.set(i, j, QN0.get(i, j) / M);
				}
			}
		}
		Matrix CN = new Matrix(M, R);
		Matrix QN = QN0;
		Matrix XN = new Matrix(1, R);
		Matrix UN = new Matrix(M, R);
		int it = 1;
		for(; it <= maxiter; it++){
			Matrix QN_1 = new Matrix(QN);
			for(int r = 0; r < R; r++){
				for(int i = 0; i < M; i++){
					CN.set(i, r, L.get(i, r));
					if(L.get(i, r) == 0){
						// 0 service demand at this station => this class does not visit the current node
						continue;
					}
					for(int s = 0; s < R; s++){
						if(s != r){
							if(type[i] == SchedStrategy.FCFS){
								CN.set(i, r, CN.get(i, r) + L.get(i, s) * QN.get(i, s));
							} else {
								CN.set(i, r, CN.get(i, r) + L.get(i, r) * QN.get(i, s));
							}
						} else {
							CN.set(i, r, CN.get(i, r) + L.get(i, r) * QN.get(i, r) * (N.get(r) - 1)/N.get(r));
						}
					}
				}
				XN.set(r, N.get(r) / (Z.get(r) + Matrix.extractColumn(CN, r, null).elementSum()));
			}
			for(int r = 0; r < R; r++){
				for(int i = 0; i < M; i++){
					QN.set(i, r, XN.get(r) * CN.get(i, r));
				}
			}
			for(int r = 0; r < R; r++){
				for(int i = 0; i < M; i++){
					UN.set(i, r, XN.get(r) * L.get(i, r));
				}
			}
			double maxabs = Double.MIN_VALUE;
			for(int i = 0; i < QN.getNumRows(); i++){
				for(int j = 0; j < QN.getNumCols(); j++){
					double val = Math.abs(1 - QN.get(i, j) / QN_1.get(i, j));
					maxabs = Maths.max(maxabs, val);
				}
			}
			if(maxabs < tol){
				break;
			}
		}
		Matrix RN = XN.repmat(M, 1);
		for(int i = 0; i < RN.getNumRows(); i++){
			for(int j = 0; j < RN.getNumCols(); j++){
				RN.set(i, j, QN.get(i, j) / RN.get(i, j));
			}
		}
		return new pfqnBSReturn(XN, QN, UN, RN, it);
	}

	/**
	 * Linearizer approximate mean value analysis algorithm
	 */
	public static pfqnLinearizerReturn pfqn_linearizer(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type){
		return pfqn_linearizer(L, N, Z, type, 1.0e-8);
	}

	/**
	 * Linearizer approximate mean value analysis algorithm
	 */
	public static pfqnLinearizerReturn pfqn_linearizer(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type, double tol){
		return pfqn_linearizer(L, N, Z, type, tol, 1000);
	}

	/**
	 * Linearizer approximate mean value analysis algorithm
	 * @param L - the service demand matrix
	 * @param N - the population vector
	 * @param Z - the think times
	 * @param type - the types of the scheduling disciplines at each station
	 * @param tol - max tolerance admitted between successive iterations
	 * @param maxiter - maximum number of iterations
	 * @return - the performance measures for the given network
	 */
	public static pfqnLinearizerReturn pfqn_linearizer(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type, double tol, int maxiter){
		return pfqn_gflinearizer(L, N, Z, type, tol, maxiter, 1);
	}

	/* General-form linearizer algorithm */
	public static pfqnLinearizerReturn pfqn_gflinearizer(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type, double tol,
                                                         int maxiter, double alpha){
		Matrix alphaM = new Matrix(1, N.getNumCols());
		alphaM.fill(alpha);
		return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alphaM);
	}

	/**
	 * Extended general form linearizer approximate mean value analysis algorithm
	 * @param L - the service demand matrix
	 * @param N - the population vector
	 * @param Z - the think times
	 * @param type - the types of the scheduling disciplines at each station
	 * @param tol - max tolerance admitted between successive iterations
	 * @param maxiter - maximum number of iterations
	 * @param alpha - matrix of alphas. There is one alpha for each class
	 * @return - the performance measures for the given network
	 */
	public static pfqnLinearizerReturn pfqn_egflinearizer(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type, double tol,
                                                          int maxiter, Matrix alpha){
		int M = L.getNumRows();
		int R = L.getNumCols();

		if(Z.isEmpty()){
			Z = new Matrix(1, R);
		}
		Z = Z.sumCols();
		boolean Lmaxcols = true;
		for(int j = 0; j < L.getNumCols(); j++){
			double maxcol = 0;
			for(int i = 0; i < L.getNumRows(); i++){
				if(L.get(i, j) > maxcol)
					maxcol = L.get(i, j);
			}
			if(maxcol != 0){
				Lmaxcols = false;
				break;
			}
		}
		if(L.isEmpty() || Lmaxcols){
			Matrix X = new Matrix(N.getNumRows(), N.getNumCols());
			for(int i = 0; i < X.getNumRows(); i++){
				for(int j = 0; j < X.getNumCols(); j++){
					X.set(i, j, N.get(i, j) / Z.get(i, j));
				}
			}
			Matrix Q = new Matrix(M, R);
			Matrix U = new Matrix(M, R);
			Matrix W = new Matrix(M, R);
			Matrix C = new Matrix(1, R);
			for(int r = 0; r < R; r++){
				for(int i = 0; i < M; i++){
					U.set(i, r, X.get(r) * L.get(i, r));
				}
			}
			int totiter = 0;
			return new pfqnLinearizerReturn(Q, U, W, C, X, totiter);
		}
		// Initialise
		Matrix[] Q = new Matrix[M];
		Matrix[] Delta = new Matrix[M];
		for(int i = 0; i < M; i++){
			Q[i] = new Matrix(R, 1 + R);
			Delta[i] = new Matrix(R, R);
		}

		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				for(int s = -1; s < R; s++){
					Matrix N_1 = oner(N, new ArrayList<>(Collections.singletonList(s)));
					Q[i].set(r, 1 + s, N_1.get(r) / M);
				}
			}
		}

		int totiter = 0;

		// Main loop
		for(int I = 0; I < 3; I++){
			for(int s = -1; s < R; s++){
				Matrix N_1 = oner(N, new ArrayList<>(Collections.singletonList(s)));
				Matrix Q1 = new Matrix(M, R);
				for(int i = 0; i < M; i++){
					for(int j = 0; j < R; j++){
						Q1.set(i, j, Q[i].get(j, 1+s));
					}
				}
				pfqnLinearizerCoreReturn ret1 = egflinearizer_core(L, M, R, N_1, Z, Q1,
						Delta, type, tol, maxiter - totiter, alpha);
				for(int i = 0; i < M; i++){
					for(int j = 0; j < R; j++){
						Q[i].set(j, 1 + s, ret1.Q.get(i, j));
					}
				}

				totiter += ret1.iter;
			}
			// Upgrade delta
			for(int i = 0; i < M; i++){
				for(int r = 0; r < R; r++){
					if(N.get(r) == 1){
						for(int j = 0; j < R; j++){
							Q[i].set(j, 1+r, 0);
						}
					}
					for(int s = 0; s < R; s++){
						if(N.get(s) > 1){
							Matrix Ns = oner(N, new ArrayList<>(Collections.singletonList(s)));
							Delta[i].set(r, s, Q[i].get(r, 1 + s)/Math.pow(Ns.get(r), alpha.get(r)) - Q[i].get(r, 0)/Math.pow(N.get(r), alpha.get(r)));
						}
					}
				}
			}
		}

		// Core(N)
		Matrix Q1 = new Matrix(M, R);
		for(int i = 0; i < M; i++){
			for(int j = 0; j < R; j++){
				Q1.set(i, j, Q[i].get(j, 0));
			}
		}
		pfqnLinearizerCoreReturn ret1 = egflinearizer_core(L, M, R, N, Z, Q1, Delta, type, tol, maxiter - totiter, alpha);
		Matrix newQ = ret1.Q;
		Matrix W = ret1.W;
		Matrix X = ret1.T;
		totiter += ret1.iter;
		// Compute performance metrics
		Matrix U = new Matrix(M, R);
		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				U.set(i, r, X.get(r) * L.get(i, r));
			}
		}
		Matrix C = new Matrix(1, R);
		for(int i = 0; i < R; i++){
			C.set(0, i, N.get(i) / X.get(i) - Z.get(i));
		}
		return new pfqnLinearizerReturn(newQ, U, W, C, X, totiter);
	}

	protected static pfqnLinearizerCoreReturn egflinearizer_core(Matrix L, int M, int R, Matrix N_1, Matrix Z, Matrix Q,
                                                               Matrix[] Delta, SchedStrategy[] type, double tol,
                                                               int maxiter, Matrix alpha){
		boolean hasConverged = false;
		Matrix W = new Matrix(L);
		int iter = 0;
		Matrix T = null;
		while(!hasConverged){
			Matrix Qlast = new Matrix(Q);
			// Estimate population
			pfqnLinearizerEstimateReturn ret1 = egflinearizer_estimate(L, M, R, N_1, Z, Q, Delta, W, alpha);
			Matrix[] Q_1 = ret1.Q_1;
			// Forward MVA
			pfqnLinearizerForwardMVAReturn ret2 = egflinearizer_forwardMVA(L, M, R, type, N_1, Z, Q_1);
			Q = ret2.Q;
			W = ret2.W;
			T = ret2.T;
			if(Q.sub(1, Qlast).norm() < tol || iter > maxiter){
				hasConverged = true;
			}
			iter++;
		}
		return new pfqnLinearizerCoreReturn(Q, W, T, iter);
	}

	protected static pfqnLinearizerEstimateReturn egflinearizer_estimate(Matrix L, int M, int R, Matrix N_1, Matrix Z,
                                                                      Matrix Q, Matrix[] Delta, Matrix W, Matrix alpha){
		Matrix[] Q_1 = new Matrix[M];
		for(int i = 0; i < M; i++){
			Q_1[i] = new Matrix(R, 1 + R);
		}
		Matrix T_1 = new Matrix(R, 1 + R);
		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				for(int s = 0; s < R; s++){
					Matrix Ns = oner(N_1, new ArrayList<>(Collections.singletonList(s)));
					Q_1[i].set(r, 1 + s, Math.pow(Ns.get(r), alpha.get(r))*(Q.get(i, r)/Math.pow(N_1.get(r), alpha.get(r)) + Delta[i].get(r,s)));
				}
			}
		}
		return new pfqnLinearizerEstimateReturn(Q_1, T_1);
	}

	protected static pfqnLinearizerForwardMVAReturn egflinearizer_forwardMVA(Matrix L, int M, int R, SchedStrategy[] type,
                                                                          Matrix N_1, Matrix Z, Matrix[] Q_1){
		Matrix W = new Matrix(M, R);
		Matrix T = new Matrix(1, R);
		Matrix Q = new Matrix(M, R);

		// Compute residence time
		for(int i = 0; i < M; i++){
			for(int r = 0; r < R; r++){
				boolean flag = true; // flag = whether all nodes have the FCFS scheduling strategy
				for(int k = 0; k < M; k++){
					if (type[k] != SchedStrategy.FCFS) {
						flag = false;
						break;
					}
				}
				if(flag){
					W.set(i, r, L.get(i, r));
					if(L.get(i, r) != 0){
						for(int s = 0; s < R; s++){
							W.set(i, r, W.get(i, r) + L.get(i, s) * Q_1[i].get(s, 1 + r));
						}
					}
				} else {
					double sumQ = 0;
					for(int k = 0; k < Q_1[i].getNumRows(); k++){
						sumQ += Q_1[i].get(k, 1 + r);
					}
					W.set(i, r, L.get(i, r) * (1 + sumQ));
				}
			}
		}

		// Compute throughputs and queue lengths
		for(int r = 0; r < R; r++){
			T.set(r, N_1.get(r) / (Z.get(r) + Matrix.extractColumn(W, r, null).elementSum()));
			for(int i = 0; i < M; i++){
				Q.set(i, r, T.get(r) * W.get(i, r));
			}
		}
		return new pfqnLinearizerForwardMVAReturn(Q, W, T);
	}


	/**
	 * Linearizer method for mixed models with multi-server stations
	 * @param lambda - arrival rate of open classes
	 * @param L - the service demand matrix
	 * @param N - the population vector
	 * @param Z - the think times vector
	 * @param nservers - number of servers at the stations
	 * @param type - scheduling strategy type
	 * @param tol - max tolerance admitted between successive iterations
	 * @param maxiter - maximum number of iterations
	 * @param method - solution method
	 * @return - the performance metrics for this network
	 */
	public static pfqnLinearizerReturn pfqn_linearizermx(Matrix lambda, Matrix L, Matrix N, Matrix Z, Matrix nservers,
                                                         SchedStrategy[] type, double tol, int maxiter, String method){
		for(int i = 0; i < lambda.getNumCols(); i++){
			if(lambda.get(i) > 0 && N.get(i) > 0 && Double.isFinite(N.get(i))){
				throw new RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.");
			}
		}
		int M = L.getNumRows();
		int R = L.getNumCols();

		for(int i = 0; i < lambda.getNumRows(); i++){
			for(int j = 0; j < lambda.getNumCols(); j++){
				if(Double.isNaN(lambda.get(i, j))){
					lambda.set(i, j, 0);
				}
			}
		}
		for(int i = 0; i < L.getNumRows(); i++){
			for(int j = 0; j < L.getNumCols(); j++){
				if(Double.isNaN(L.get(i, j))){
					L.set(i, j, 0);
				}
			}
		}
		for(int i = 0; i < Z.getNumRows(); i++){
			for(int j = 0; j < Z.getNumCols(); j++){
				if(Double.isNaN(Z.get(i, j))){
					Z.set(i, j, 0);
				}
			}
		}

		ArrayList<Integer> openClasses = new ArrayList<>();
		ArrayList<Integer> closedClasses = new ArrayList<>();
		for(int i = 0; i < N.length(); i++){
			if(Double.isInfinite(N.get(i))){
				openClasses.add(i);
			} else {
				closedClasses.add(i);
			}
		}

		Matrix XN = new Matrix(1, R);
		Matrix UN = new Matrix(M,R);
		Matrix WN = new Matrix(M, R);
		Matrix QN = new Matrix(M,R);
		Matrix CN = new Matrix(1,R);

		for(Integer r : openClasses){
			for(int i = 0; i < M; i++){
				UN.set(i,r,lambda.get(r) * L.get(i,r));
			}
			XN.set(0,r,lambda.get(r));
		}

		Matrix UNt = UN.sumRows();

		if(Z.isEmpty()){
			Z = new Matrix(1,R);
		} else {
			Z = Z.sumCols();
		}
		Matrix Dc = new Matrix(L.getNumRows(), closedClasses.size());
		Matrix rep = UNt.repmat(1, closedClasses.size());
		for(int i = 0; i < Dc.getNumRows(); i++){
			int j = 0;
			for(int closedClass : closedClasses){
				Dc.set(i, j, L.get(i, closedClass) / (1 - rep.get(i,j)));
				j++;
			}
		}

		Matrix Nclosed = new Matrix(1, closedClasses.size());
		Matrix Zclosed = new Matrix(1, closedClasses.size());
		int idx = 0;
		for(int closedClass : closedClasses){
			Nclosed.set(0, idx, N.get(closedClass));
			Zclosed.set(0, idx, Z.get(closedClass));
			idx++;
		}

		Matrix QNc, UNc, WNc, CNc, XNc;
		int totiter;

		if(nservers.elementMax() == 1){
			pfqnLinearizerReturn res = null;
			if(method.equals("lin")){
				res = pfqn_linearizer(Dc, Nclosed, Zclosed, type, tol, maxiter);
			} else if(method.equals("gflin")){
				double linAlpha = 2.0;
				res = pfqn_gflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, linAlpha);
			} else {
				// Extended General Form Linearizer
				Matrix alphaM = new Matrix(1, Nclosed.getNumCols());
				/**
				 * Gompertz function to approximate alpha as a function of the network population
				 */
				for(int i = 0; i < Nclosed.getNumCols(); i++){
					alphaM.set(i, 0.6 + 1.4 * Math.exp(-8 * Math.exp(-0.8 * Nclosed.get(i))));
				}
				res = pfqn_egflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, alphaM);
			}
			QNc = res.Q;
			UNc = res.U;
			WNc = res.W;
			CNc = res.C;
			XNc = res.X;
			totiter = res.totiter;
		} else {
			Matrix typeMatrix = new Matrix(type.length, 1);
			for(int i = 0; i < type.length; i++){
				typeMatrix.set(i, SchedStrategy.toID(type[i]));
			}
			pfqnLinearizerMSReturn res = pfqn_linearizerms(Dc, Nclosed, Zclosed, nservers,
					typeMatrix, tol, maxiter);
			QNc = res.Q;
			UNc = res.U;
			WNc = res.R;
			CNc = res.C;
			XNc = res.X;
			totiter = res.totiter;
		}

		for(int i = 0; i < closedClasses.size(); i++){
			XN.set(closedClasses.get(i), XNc.get(i));
			for(int j = 0; j < QN.getNumRows(); j++){
				QN.set(j, closedClasses.get(i), QNc.get(j, i));
			}
			for(int j = 0; j < WN.getNumRows(); j++){
				WN.set(j, closedClasses.get(i), WNc.get(j, i));
			}
			for(int j = 0; j < UN.getNumRows(); j++){
				UN.set(j, closedClasses.get(i), UNc.get(j, i));
			}
			CN.set(closedClasses.get(i), CNc.get(i));
		}

		for(int i = 0; i < M; i++){
			for(int r : closedClasses){
				UN.set(i, r, XN.get(r) * L.get(i, r));
			}
		}

		for(int i = 0; i < M; i++){
			for(int r : openClasses){
				if(QNc.isEmpty()){
					WN.set(i,r,L.get(i,r) / (1 - UNt.get(i)));
				} else {
					WN.set(i, r, L.get(i, r) * (1 + QNc.sumRows(i)) / (1 - UNt.get(i)));
				}
				QN.set(i, r, WN.get(i, r) * XN.get(r));
			}
		}
		for(int r : openClasses){
			CN.set(r, WN.sumCols(r));
		}
		return new pfqnLinearizerReturn(QN, UN, WN, CN, XN, totiter);
	}

    public static class pfqnMVAReturn {
        public Matrix XN;
        public Matrix QN;
        public Matrix UN;
        public Matrix CN;
        public double lGN;

        public pfqnMVAReturn(Matrix XN, Matrix QN, Matrix UN, Matrix CN, double lGN) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.CN = CN;
            this.lGN = lGN;
        }
    }

    public static class pfqnBSReturn {
        public Matrix XN;
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public int it;

        public pfqnBSReturn(Matrix XN, Matrix QN, Matrix UN, Matrix RN, int it) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.it = it;
        }
    }

    public static class pfqnMVALDMXReturn {
        public Matrix XN;
        public Matrix QN;
        public Matrix UN;
        public Matrix CN;
        public double lGN;
        public Matrix Pc;

        public pfqnMVALDMXReturn(Matrix XN, Matrix QN, Matrix UN, Matrix CN, double lGN, Matrix Pc) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.CN = CN;
            this.lGN = lGN;
            this.Pc = Pc;
        }
    }

    public static class pfqnMVALDMXECReturn {
        public Matrix EC;
        public Matrix E;
        public Matrix Eprime;
        public Matrix Lo;

        public pfqnMVALDMXECReturn(Matrix EC, Matrix E, Matrix Eprime, Matrix Lo) {
            this.EC = EC;
            this.E = E;
            this.Eprime = Eprime;
            this.Lo = Lo;
        }
    }

    public static class pfqnMVALDReturn {
        public Matrix XN;
        public Matrix QN;
        public Matrix UN;
        public Matrix CN;
        public List<Double> lGN;
        public boolean isNumStable;
        public Matrix pi;

        public pfqnMVALDReturn(Matrix XN, Matrix QN, Matrix UN, Matrix CN, List<Double> lGN, boolean isNumStable, Matrix pi) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.CN = CN;
            this.lGN = lGN;
            this.isNumStable = isNumStable;
            this.pi = pi;
        }
    }

    public static class pfqnLinearizerMSReturn {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix C;
        public Matrix X;
        public int totiter;

        public pfqnLinearizerMSReturn(Matrix Q, Matrix U, Matrix R, Matrix C, Matrix X, int totiter) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.C = C;
            this.X = X;
            this.totiter = totiter;
        }
    }

    public static class pfqnLinearizerMSCoreReturn {
        public Matrix Q;
        public Matrix W;
        public Matrix T;
        public Matrix P;
        public Matrix PB;
        public int iter;

        public pfqnLinearizerMSCoreReturn(Matrix Q, Matrix W, Matrix T, Matrix P, Matrix PB, int iter) {
            this.Q = Q;
            this.W = W;
            this.T = T;
            this.P = P;
            this.PB = PB;
            this.iter = iter;
        }
    }

    public static class pfqnLinearizerMSEstimateReturn {
        public Matrix[] Q_1;
        public Matrix[] P_1;
        public Matrix PB_1;

        public pfqnLinearizerMSEstimateReturn(Matrix[] Q_1, Matrix[] P_1, Matrix PB_1) {
            this.Q_1 = Q_1;
            this.P_1 = P_1;
            this.PB_1 = PB_1;
        }
    }

    public static class pfqnLinearizerMSForwardMVAReturn {
        public Matrix Q;
        public Matrix W;
        public Matrix T;
        public Matrix P;
        public Matrix PB;

        public pfqnLinearizerMSForwardMVAReturn(Matrix Q, Matrix W, Matrix T, Matrix P, Matrix PB) {
            this.Q = Q;
            this.W = W;
            this.T = T;
            this.P = P;
            this.PB = PB;
        }
    }

    public static class pfqnLinearizerReturn{
        public Matrix Q;
        public Matrix U;
        public Matrix W;
        public Matrix C;
        public Matrix X;
        public int totiter;

        public pfqnLinearizerReturn(Matrix Q, Matrix U, Matrix W, Matrix C, Matrix X, int totiter) {
            this.Q = Q;
            this.U = U;
            this.W = W;
            this.C = C;
            this.X = X;
            this.totiter = totiter;
        }
    }

    public static class pfqnLinearizerCoreReturn {
        public Matrix Q;
        public Matrix W;
        public Matrix T;
        public int iter;

        public pfqnLinearizerCoreReturn(Matrix Q, Matrix W, Matrix T, int iter) {
            this.Q = Q;
            this.W = W;
            this.T = T;
            this.iter = iter;
        }
    }

    public static class pfqnLinearizerEstimateReturn{
        public Matrix[] Q_1;
        public Matrix T_1;

        public pfqnLinearizerEstimateReturn(Matrix[] Q_1, Matrix T_1) {
            this.Q_1 = Q_1;
            this.T_1 = T_1;
        }
    }

    public static class pfqnLinearizerForwardMVAReturn {
        public Matrix Q;
        public Matrix W;
        public Matrix T;

        public pfqnLinearizerForwardMVAReturn(Matrix Q, Matrix W, Matrix T) {
            this.Q = Q;
            this.W = W;
            this.T = T;
        }
    }

    static class pfqnLeFpiReturn {
        public final Matrix u;
        public final Matrix d;

        public pfqnLeFpiReturn(Matrix u, Matrix d) {
            this.u = u;
            this.d = d;
        }
    }

    static class pfqnLeFpiZReturn {
        public final Matrix u;
        public final double v;
        public final Matrix d;

        public pfqnLeFpiZReturn(Matrix u, double v, Matrix d) {
            this.u = u;
            this.v = v;
            this.d = d;
        }
    }

    public static class pfqnComomrmReturn {
        public double lG;
        public Matrix lGbasis;

        public pfqnComomrmReturn(double lG, Matrix lGbasis) {
            this.lG = lG;
            this.lGbasis = lGbasis;
        }
    }

    public static class pfqnNcSanitizeReturn {
        public Matrix lambda;
        public Matrix L;
        public Matrix N;
        public Matrix Z;
        public double lGremaind;

        public pfqnNcSanitizeReturn(Matrix lambda, Matrix L, Matrix N, Matrix Z, double lGremaind) {
            this.lambda = lambda;
            this.L = L;
            this.N = N;
            this.Z = Z;
            this.lGremaind = lGremaind;
        }
    }

    public static class pfqnNcXQReturn {
        public Double G;
        public Double lG;
        public Matrix X;
        public Matrix Q;
        public String method;

        public pfqnNcXQReturn(Double G, Double lG, Matrix X, Matrix Q, String method) {
            this.G = G;
            this.lG = lG;
            this.X = X;
            this.Q = Q;
            this.method = method;
        }

        public pfqnNcXQReturn(Double lG, Matrix X, Matrix Q, String method) {
            this.G = Math.exp(lG);
            this.lG = lG;
            this.X = X;
            this.Q = Q;
            this.method = method;
        }
    }

    public static class pfqnNcReturn {
        public Double G;
        public Double lG;
        public String method;

        public pfqnNcReturn(Double G, Double lG) {
            this.G = G;
            this.lG = lG;
            this.method = null;
        }
        public pfqnNcReturn(Double G, Double lG, String method) {
            this.lG = lG;
            this.G = G;
            this.method = method;
        }
    }
}
