package jline.util;

import java.util.*;

import jline.lang.constant.GlobalConstants;
import org.ejml.data.*;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.interfaces.decomposition.QRSparseDecomposition;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.CommonOps_DSCC;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.MatrixFeatures_DDRM;

import org.ejml.sparse.csc.decomposition.qr.QrLeftLookingDecomposition_DSCC;
import org.jblas.DoubleMatrix;

/**
 * A sparse matrix data structure supporting linear algebra functions similar to those available in MATLAB.
 */
public class Matrix extends DMatrixSparseCSC {

	public Matrix(double scalar) {
		super(1, 1, 1);
		this.set(0,0,scalar);
	}

	public Matrix(int numRows, int numCols, int arrayLength) {
		super(numRows, numCols, arrayLength);
	}

	public Matrix(int numRows, int numCols) {
		super(numRows, numCols, 0);
	}

	public Matrix(Matrix matrix) {
		super(matrix.copy());
	}

	public Matrix(DMatrixSparseCSC matrix) {
		super(matrix);
	}

	public Matrix(List<Double> array) {
		super(array.size(), 1, array.size());
		for(int i = 0; i < array.size(); i++)
			this.set(i, 0, (double) array.get(i));
	}

	public Matrix(SimpleMatrix matrix) {
		super(matrix.numRows(), matrix.numCols());
		for (int i = 0; i < matrix.numRows(); i++) {
			for (int j = 0; j < matrix.numCols(); j++) {
				this.set(i, j, matrix.get(i, j));
			}
		}
	}

	public static Matrix allbut(Matrix y, int xset) {
		int index = -1;
		if (xset >= 0 && xset < y.length()) {
			index = xset;
		}
		Matrix ret;
		if (index == -1) {
			ret = new Matrix(1, y.length());
		} else {
			ret = new Matrix(1, y.length()-1);
		}
		int i = 0;
		for (int j = 0; j < y.length(); j++) {
			if (j != index) {
				ret.set(i, y.get(j));
				i++;
			}
		}
		return ret;
	}

	/**
	 * Returns the position of the given row in the corresponding matrix
	 * @param matrix - the matrix to be searched
	 * @param row - the row
	 * @return - Position of the given row in the matrix, or -1 otherwise
	 */
	public static int matchrow(Matrix matrix, Matrix row){
		if(matrix.getNumCols() != row.getNumCols())
			return -1;
		for(int i = 0; i < matrix.getNumRows(); i++){
			boolean rowsEqual = true;
			for(int j = 0; j < matrix.getNumCols(); j++){
				if(matrix.get(i, j) != row.get(j)){
					rowsEqual = false;
					break;
				}
			}
			if(rowsEqual)
				return i;
		}
		return -1;
	}

	/**
	 * Weakly-connected components of a sub-matrix.
	 *
	 * @param param Input matrix
	 * @param colsToIgnore Indexes to be ignored
	 * @return Weighted connected components
	 */
	public static Set<Set<Integer>> weaklyConnect(Matrix param, Set<Integer> colsToIgnore) {
		UndirectedGraph graph = new UndirectedGraph(param, colsToIgnore);
		graph.computeWeaklyConnectedComponents();
		return graph.getWCC();
	}

	/**
	 * Decrease by one an element of an integer vector.
	 *
	 * @param N integer vector
	 * @param r dimension to decrease
	 * @return Decreased vector
	 */
	public static Matrix oner(Matrix N, List<Integer> r) {
		Matrix res = N.clone();
		for (Integer s : r) {
			if (s >= 0)
				res.set(s, res.get(s) - 1);
		}
		return res;
	}

	public static Matrix factln(Matrix n) {
		Matrix ret = n.clone();
		for (int i = 0; i < n.length(); i++) {
			ret.set(i, Maths.factln(ret.get(i)));
		}
		return ret;
	}

	public static double logsumexp(Matrix x) {
		int n = x.length();
		double a = x.elementMax();
		int k = -1;
		for (int i = 0; i < n; i++) {
			if (x.get(i) == a) {
				k = i;
				break;
			}
		}
		Matrix w = new Matrix(1, n);
		w.fill(0.0);
		double s = 0;

		for (int i = 0; i < n; i++) {
			w.set(i, Math.exp(x.get(i)-a));
			if (i != k) {
				s += w.get(i);
			}
		}
		return (a + Math.log1p(s));
	}

	public static Matrix decorate(Matrix inSpace1, Matrix inSpace2) {

	  // TODO: upfront if clause for 1 parameter, lines 7 to 14

	  if (inSpace1.isEmpty()) {
		inSpace1 = inSpace2.clone();
		return inSpace1;
	  }

	  if (inSpace2.isEmpty()) {
		return inSpace1;
	  }

	  int n1 = inSpace1.getNumRows();
	  int m1 = inSpace1.getNumCols();
	  int n2 = inSpace2.getNumRows();
	  int m2 = inSpace2.getNumCols();

	  inSpace1.repmat(n2, 1);
	  int curStatesStart = 0;
	  int curStatesEnd = n1;

	  for (int s = 0; s < n2; s++) {
		Matrix tmp = new Matrix(1, inSpace2.getNumCols());
		extractRows(inSpace2, s, s + 1, tmp);
		tmp.repmat(curStatesEnd - 1, 1);

		inSpace1.expandMatrix(
			inSpace1.getNumRows() + curStatesEnd - 1,
			inSpace1.getNumCols() + m2,
			(inSpace1.getNumRows() + curStatesEnd - 1) * (inSpace1.getNumCols() + m2));
		for (int i = curStatesStart; i < curStatesEnd; i++) {
		  for (int j = m1; j < m1 + m2; j++) {
			inSpace1.set(i, j, tmp.get(i - curStatesStart, j - m1));
		  }
		}
		curStatesStart += n1;
		curStatesEnd += n1;
	  }

	  return inSpace1;
	}

	public Matrix inv() {
		DMatrixRMaj inverse = new DMatrixRMaj(this.numRows, this.numCols);
		Matrix thisinverse = new Matrix(this.numRows, this.numCols);
		CommonOps_DSCC.invert(this,inverse);
		DConvertMatrixStruct.convert(inverse,thisinverse);
		return thisinverse;
	}

	public static List<Double> intersect(Matrix matrixA, Matrix matrixB) {

		Set<Double> matrixAValues = new HashSet<>();
		Set<Double> matrixBValues = new HashSet<>();
		List<Double> outputValues = new LinkedList<>();

		int rows = matrixA.getNumRows();
		int cols = matrixA.getNumCols();
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				matrixAValues.add(matrixA.get(row, col));
			}
		}

		rows = matrixB.getNumRows();
		cols = matrixB.getNumCols();
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				matrixBValues.add(matrixB.get(row, col));
			}
		}

		for (double value : matrixAValues) {
			if (matrixBValues.contains(value)) {
				outputValues.add(value);
			}
		}

		return outputValues;
	}

	public DMatrixSparseCSC toDMatrixSparseCSC() {
		return this.copy();
	}

	public DMatrixSparseCSC toDMatrixSparseCSC(Matrix matrix) {
		return matrix.copy();
	}

	public void expandMatrix(int rows, int cols, int nz_length) {
		if (rows < this.getNumRows() || cols < this.getNumCols()) {
			return;
		}

		DMatrixSparseTriplet nodeRouting = new DMatrixSparseTriplet(rows, cols, nz_length);
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int col1 = this.col_idx[colIdx];
			int col2 = this.col_idx[colIdx+1];

			for(int i = col1; i < col2; i++) {
				int rowIdx = this.nz_rows[i];
				double value = this.nz_values[i];
				nodeRouting.addItem(rowIdx, colIdx, value);
			}
		}
		this.setTo(DConvertMatrixStruct.convert(nodeRouting, (DMatrixSparseCSC)null));
	}

	public boolean isDiag() {
		if (this.numCols != this.numRows)
			return false;

		if (this.getNonZeroLength() != this.numCols)
			return false;

		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int col1 = this.col_idx[colIdx];
			int col2 = this.col_idx[colIdx+1];

			for(int i = col1; i < col2; i++) {
				int rowIdx = this.nz_rows[i];
				if (rowIdx != colIdx)
					return false;
			}
		}

		return true;
	}

	public Matrix clone() {
		return new Matrix(this);
	}

	public boolean hasNaN() {
		for(int i = 0; i < this.nz_length; i++) {
			if (Double.isNaN(this.nz_values[i]))
				return true;
		}
		return false;
	}

	public double get(int idx) {
		if (idx >= this.numCols * this.numRows)
			throw new RuntimeException("Index out of matrix");

		int row = idx % this.getNumRows();
		int col = idx / this.getNumRows();

		return super.get(row, col);
	}

	public void set(int idx, double val) {
		if (idx >= this.numCols * this.numRows)
			throw new RuntimeException("Index out of matrix");

		int row = idx % this.getNumRows();
		int col = idx / this.getNumRows();

		this.set(row, col, val);
	}

	public void set(int row, int col, double val) {
		super.set(row, col, val);

		//This to ensure the value NaN is replaced
		if (val == 0)
			super.remove(row, col); //Remove to ensure that getNonZeroElement not contains the value with 0
	}

	public Matrix cumsumViaRow() {
		Matrix res = new Matrix(this.numRows, this.numCols, this.numRows * this.numCols);
		for(int i = 0; i < this.numRows; i++)
			res.set(i, 0, this.get(i, 0));

		for(int i = 0; i < this.numRows; i++) {
			for(int j = 1; j < this.numCols; j++) {
				res.set(i, j, this.get(i, j) + res.get(i, j-1));
			}
		}
		return res;
	}

	public Matrix cumsumViaCol() {
		Matrix res = new Matrix(this.numRows, this.numCols, this.numRows * this.numCols);
		for(int i = 0; i < this.numCols; i++)
			res.set(0, i, this.get(0, i));

		for(int i = 0; i < this.numCols; i++) {
			for(int j = 1; j < this.numRows; j++) {
				res.set(j, i, this.get(j, i) + res.get(j-1, i));
			}
		}
		return res;
	}

	public double sumRows(int row) {
		double sum = 0;
		for(int i = 0; i < this.numCols; i++) {
			sum += this.get(row, i);
		}
		return sum;
	}

	public Matrix sumRows() {
		DMatrixRMaj sumrows = CommonOps_DSCC.sumRows(this, null);
		DMatrixSparseCSC tmp = new DMatrixSparseCSC(0,0);
		DConvertMatrixStruct.convert(sumrows, tmp);
		return new Matrix(tmp);
	}
	public double sumCols(int col) {
		double sum = 0;
		for(int i = 0; i < this.numRows; i++) {
			sum += this.get(i, col);
		}
		return sum;
	}

	public Matrix sumCols() {
		DMatrixRMaj sumcols = CommonOps_DSCC.sumCols(this, null);
		DMatrixSparseCSC tmp = new DMatrixSparseCSC(0,0);
		DConvertMatrixStruct.convert(sumcols, tmp);
		return new Matrix(tmp);
	}

	public double sumAbsCols(int col) {
		double sum = 0;
		for(int i = 0; i < this.numRows; i++) {
			sum += Math.abs(this.get(i, col));
		}
		return sum;
	}

	public double sumAbsRows(int row) {
		double sum = 0;
		for(int i = 0; i < this.numRows; i++) {
			sum += Math.abs(this.get(row,i));
		}
		return sum;
	}

	public Matrix repmat(int rows, int cols) {
		Matrix res = this.clone();
		for(int i = 1; i < rows; i++) {
			Matrix tmp = new Matrix(0,0,0);
			CommonOps_DSCC.concatRows(res, this, tmp);
			res = tmp;
		}
		for(int i = 1; i < cols; i++) {
			Matrix tmp = new Matrix(0,0,0);
			CommonOps_DSCC.concatColumns(res, res, tmp);
			res = tmp;
		}
		return res;
	}

	public Matrix find() {
		Matrix res = new Matrix(this.nz_length, 1, this.nz_length);
		int count = 0;
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int col1 = this.col_idx[colIdx];
			int col2 = this.col_idx[colIdx+1];

			for(int i = col1; i < col2; i++) {
				int rowIdx = this.nz_rows[i];
				res.set(count++, 0, colIdx * this.numRows + rowIdx);
			}
		}
		return res;
	}

	public Matrix findNonNegative() {
		List<Integer> array = new ArrayList<Integer>();
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			for(int rowIdx = 0; rowIdx < this.numRows; rowIdx++) {
				if (this.get(rowIdx, colIdx) >= 0)
					array.add(colIdx * this.numRows + rowIdx);
			}
		}

		Matrix res = new Matrix(array.size(), 1, array.size());
		for(int i = 0; i < array.size(); i++)
			res.set(i, 0, array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueInRow(int rowIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			array.add((int)this.get(rowIdx,colIdx));
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueNonZerosInRow(int rowIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int val = (int) this.get(rowIdx,colIdx);
			if (val != 0) {
				array.add(val);
			}
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueNonNegativeInRow(int rowIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int val = (int) this.get(rowIdx,colIdx);
			if (val >0) {
				array.add(val);
			}
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueInCol(int colIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int rowIdx = 0; rowIdx < this.numRows; rowIdx++) {
			array.add((int) this.get(rowIdx,colIdx));
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueNonZerosInCol(int colIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int rowIdx = 0; rowIdx < this.numRows; rowIdx++) {
			int val = (int) this.get(rowIdx,colIdx);
			if (val != 0) {
				array.add(val);
			}
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	// find unique elements in the given row
	public Matrix uniqueNonNegativeInCol(int colIdx) {
		List<Integer> array = new ArrayList<Integer>();
		for(int rowIdx = 0; rowIdx < this.numRows; rowIdx++) {
			int val = (int) this.get(rowIdx,colIdx);
			if (val >0) {
				array.add(val);
			}
		}

		List<Integer> unique_array = Utils.unique(array);
		Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
		for(int i = 0; i < unique_array.size(); i++)
			res.set(i, 0, unique_array.get(i));
		return res;
	}

	public int count(double val) {
		if (val == 0) {
			return this.getNumCols() * this.getNumRows() - this.getNonZeroLength();
		} else {
			int res = 0;
			for(int i = 0; i < this.nz_length; i++) {
				if (this.nz_values[i] == val)
					res++;
			}
			return res;
		}
	}

	public Matrix countEachRow(double val) {
		Matrix res = new Matrix(this.numRows, 1);
		for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
			int col1 = this.col_idx[colIdx];
			int col2 = this.col_idx[colIdx+1];

			for(int i = col1; i < col2; i++) {
				if (this.nz_values[i] == val) {
					int rowIdx = this.nz_rows[i];
					res.set(rowIdx, 0, res.get(rowIdx,0) + 1);
				}
			}
		}
		if (val == 0) {
			for(int i = 0; i < this.getNumRows(); i++) {
				for(int j = 0; j < this.getNumCols(); j++) {
					if (this.get(i,j) == 0) {
						res.set(i, 0, res.get(i,0) + 1);
					}
				}
			}
		}
		return res;
	}

	public int length() {
		return Math.max(numRows, numCols);
	}

	public void abs() {
		for(int i = 0; i < nz_length; i++) {
			this.nz_values[i] = Math.abs(this.nz_values[i]);
		}
	}

	public void removeNegative() {
		int offset = 0;
		for (int i = 0; i < this.numCols; i++) {
			int idx0 = this.col_idx[i] + offset;
			int idx1 = this.col_idx[i + 1];

			for (int j = idx0; j < idx1; j++) {
				double val = this.nz_values[j];
				if (val > 0) {
					this.nz_rows[j - offset] = this.nz_rows[j];
					this.nz_values[j - offset] = val;
				} else {
					offset++;
				}
			}
			this.col_idx[i + 1] -= offset;
		}
		this.nz_length -= offset;
	}

	public void removeInfinity() {
		int offset = 0;
		for (int i = 0; i < this.numCols; i++) {
			int idx0 = this.col_idx[i] + offset;
			int idx1 = this.col_idx[i + 1];

			for (int j = idx0; j < idx1; j++) {
				double val = this.nz_values[j];
				if (Double.isFinite(val)) {
					this.nz_rows[j - offset] = this.nz_rows[j];
					this.nz_values[j - offset] = val;
				} else {
					offset++;
				}
			}
			this.col_idx[i + 1] -= offset;
		}
		this.nz_length -= offset;
	}


	public void removeNaN() {
		if (!hasNaN())
			return;

		int offset = 0;
		for (int i = 0; i < this.numCols; i++) {
			int idx0 = this.col_idx[i] + offset;
			int idx1 = this.col_idx[i + 1];

			for (int j = idx0; j < idx1; j++) {
				double val = this.nz_values[j];
				if (!Double.isNaN(val)) {
					this.nz_rows[j - offset] = this.nz_rows[j];
					this.nz_values[j - offset] = val;
				} else {
					offset++;
				}
			}
			this.col_idx[i + 1] -= offset;
		}
		this.nz_length -= offset;
	}

	public boolean isEmpty() {
		return (this.numCols == 0 || this.numRows == 0);
	}

	public void apply(double source, double target, String op) {
		double tol = GlobalConstants.Zero;
		switch (op) {
			case "equal":
				if (Math.abs(source - 0) < tol) {
					if (Math.abs(target - 0) < tol) return;
					for(int i = 0; i < this.numRows; i++) {
						for(int j = 0; j < this.numCols; j++) {
							if (Math.abs(this.get(i,j) - 0) < tol)
								super.set(i, j, target);
						}
					}
				} else if (Double.isNaN(source)) {
					for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
						int col1 = this.col_idx[colIdx];
						int col2 = this.col_idx[colIdx+1];

						for(int i = col1; i < col2; i++) {
							if (Double.isNaN(this.nz_values[i])) {
								super.set(this.nz_rows[i], colIdx, target);
							}
						}
					}
				} else if (Double.isInfinite(source)) {
					for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
						int col1 = this.col_idx[colIdx];
						int col2 = this.col_idx[colIdx+1];

						for(int i = col1; i < col2; i++) {
							if (Double.isInfinite(this.nz_values[i])) {
								super.set(this.nz_rows[i], colIdx, target);
							}
						}
					}
				} else {
					for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
						int col1 = this.col_idx[colIdx];
						int col2 = this.col_idx[colIdx+1];

						for(int i = col1; i < col2; i++) {
							if (Math.abs(this.nz_values[i] - source) < tol){
								super.set(this.nz_rows[i], colIdx, target);
							}
						}
					}
				}
				break;
			case "notequal":
				if (Math.abs(source - 0) < tol) {
					if (Math.abs(target - 0) < tol) this.zero();
					for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
						int col1 = this.col_idx[colIdx];
						int col2 = this.col_idx[colIdx+1];

						for(int i = col1; i < col2; i++) {
							if ((Math.abs(this.nz_values[i] - 0) >= tol) || (Double.isNaN(this.nz_values[i]))) {
								super.set(this.nz_rows[i], colIdx, target);
							}
						}
					}
				} else if (Double.isNaN(source)) {
					for(int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if (!Double.isNaN(this.get(row, col))) {
								super.set(row, col, target);
							}
						}
					}
				} else if (Double.isInfinite(source)) {
					for(int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if (!Double.isInfinite(this.get(row, col))) {
								super.set(row, col, target);
							}
						}
					}
				} else {
					for(int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if ((Math.abs(this.get(row, col) - source) >= tol) || (Double.isNaN(this.get(row, col)))) {
								super.set(row, col, target);
							}
						}
					}
				}
				break;
			case "great":
				if (Math.abs(source - 0) < tol) {
					for(int i = 0; i < this.numRows; i++) {
						for(int j = 0; j < this.numCols; j++) {
							if (Math.abs(this.get(i,j) - 0) >= tol && Double.compare(this.get(i,j), 0) > 0) {
								super.set(i, j, target);
							}
						}
					}
				} else if (Double.isNaN(source)) {
					throw new RuntimeException("Cannot compare with NaN");
				} else if (Double.isInfinite(source)) {
					throw new RuntimeException("Cannot compare with Infinite");
				} else {
					for (int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if (Math.abs(this.get(row, col) - source) >= tol && Double.compare(this.get(row, col), source) > 0) {
								super.set(row, col, target);
							}
						}
					}
				}
				break;
			case "greatequal":
				if (Math.abs(source - 0) < tol) {
					for(int i = 0; i < this.numRows; i++) {
						for(int j = 0; j < this.numCols; j++) {
							if ((Math.abs(this.get(i,j) - 0) < tol) ||
											(Math.abs(this.get(i,j) - 0) >= tol && Double.compare(this.get(i,j), 0) > 0)) {
								super.set(i, j, target);
							}
						}
					}
				} else if (Double.isNaN(source)) {
					throw new RuntimeException("Cannot compare with NaN");
				} else if (Double.isInfinite(source)) {
					throw new RuntimeException("Cannot compare with Infinite");
				} else {
					for (int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if ((Math.abs(this.get(row, col) - source) < tol) ||
											(Math.abs(this.get(row, col) - source) >= tol && Double.compare(this.get(row, col), source) > 0)) {
								super.set(row, col, target);
							}
						}
					}
				}
				break;
			case "less":
				if (Math.abs(source - 0) < tol) {
					for(int i = 0; i < this.numRows; i++) {
						for(int j = 0; j < this.numCols; j++) {
							if (Math.abs(this.get(i,j) - 0) >= tol && Double.compare(this.get(i,j), 0) < 0) {
								super.set(i, j, target);
							}
						}
					}
				} else if (Double.isNaN(source)) {
					throw new RuntimeException("Cannot compare with NaN");
				} else if (Double.isInfinite(source)) {
					throw new RuntimeException("Cannot compare with Infinite");
				} else {
					for (int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if (Math.abs(this.get(row, col) - source) >= tol && Double.compare(this.get(row, col), source) < 0) {
								super.set(row, col, target);
							}
						}
					}
				}
				break;
			case "lessequal":
				if (Math.abs(source - 0) < tol) {
					for(int i = 0; i < this.numRows; i++) {
						for(int j = 0; j < this.numCols; j++) {
							if ((Math.abs(this.get(i,j) - 0) < tol) ||
											(Math.abs(this.get(i,j) - 0) >= tol && Double.compare(this.get(i,j), 0) < 0)) {
								super.set(i, j, target);
							}
						}
					}
				} else if (Double.isNaN(source)) {
					throw new RuntimeException("Cannot compare with NaN");
				} else if (Double.isInfinite(source)) {
					throw new RuntimeException("Cannot compare with Infinite");
				} else {
					for (int row = 0; row < this.numRows; row++) {
						for(int col = 0; col < this.numCols; col++) {
							if ((Math.abs(this.get(row, col) - source) < tol) ||
											(Math.abs(this.get(row, col) - source) >= tol && Double.compare(this.get(row, col), source) < 0)) {
								super.set(row, col, target);
							}
						}
					}
				}
				break;
			default:
				throw new RuntimeException("Operation is not supproted");
		}

		if (target == 0)
			CommonOps_DSCC.removeZeros(this, 0);
	}

	public Matrix elementIncrease(double val) {
		Matrix res = this.clone();
		for(int row = 0; row < this.numRows; row++) {
			for(int col = 0; col < this.numCols; col++) {
				res.set(row, col, res.get(row, col) + val);
			}
		}
		return res;
	}

	public Matrix meanCol() {
		Matrix res = new Matrix(1, this.numCols);
		for(int col = 0; col < this.numCols; col++) {
			res.set(0, col, this.sumCols(col) / this.numRows);
		}
		return res;
	}

	public Matrix meanRow() {
		Matrix res = new Matrix(this.numRows, 1);
		for(int row = 0; row < this.numRows; row++) {
			res.set(row, 0, this.sumRows(row) / this.numCols);
		}
		return res;
	}

	public Matrix power(double t) {
		Matrix res = this.clone();
		if (t == 0) {
			CommonOps_DSCC.fill(res, 1);
		} else if (t != 1) {
			for(int colIdx = 0; colIdx < this.numCols; colIdx++) {
				int col1 = this.col_idx[colIdx];
				int col2 = this.col_idx[colIdx+1];

				for(int i = col1; i < col2; i++) {
					int rowIdx = this.nz_rows[i];
					res.set(rowIdx, colIdx, Math.pow(res.get(rowIdx,colIdx), t));
				}
			}
		}
		return res;
	}

	public Matrix fromArray2D(int[][] matrix){
		for (int i = 0; i< matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				set(i,j,matrix[i][j]);
			}
		}
		return this;
	}

	public Matrix fromArray2D(double[][] matrix){
		for (int i = 0; i< matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				set(i,j,matrix[i][j]);
			}
		}
		return this;
	}

	public void fill(double val) {
		CommonOps_DSCC.fill(this, val);
	}

	public Matrix transpose() {
		Matrix res = new Matrix(0,0);
		CommonOps_DSCC.transpose(this, res, null);
		return res;
	}

	public Matrix sub(double alpha, Matrix matrix) {
		return new Matrix(CommonOps_DSCC.add(1, this, -alpha, matrix, null, null, null));
	}

	public Matrix add(double alpha, Matrix matrix) {
		return new Matrix(CommonOps_DSCC.add(1, this, alpha, matrix, null, null, null));
	}


	public void divide(double scalar, Matrix outputB, boolean flag) {
		if(flag)
			CommonOps_DSCC.divide(this, scalar, outputB);
		else
			CommonOps_DSCC.divide(scalar, this, outputB);
	}

	public void divideRows(double[] diag, int offset) {
		CommonOps_DSCC.divideRows(diag, offset, this);
	}

	public void multEq(Matrix B) {
		Matrix output = new Matrix(CommonOps_DSCC.mult(this, B, null));
		this.setTo(output);
	}

	public Matrix mult(Matrix B) {
		return new Matrix(CommonOps_DSCC.mult(this, B, null));
	}

	public Matrix mult(Matrix B, Matrix out) {
		return new Matrix(CommonOps_DSCC.mult(this, B, out));
	}

	public double elementSum() {
		return CommonOps_DSCC.elementSum(this);
	}
	public double elementMin() {
		return CommonOps_DSCC.elementMin(this);
	}
	public double elementMax() {
		return CommonOps_DSCC.elementMax(this);
	}

	public Matrix elementMult(Matrix B, Matrix output) {
		return new Matrix(CommonOps_DSCC.elementMult(this, B, output, null, null));
	}

	/**
	 * Performs element-wise division
	 * @param B - the other matrix
	 * @return - A ./ B
	 */
	public Matrix elementDiv(Matrix B){
		if(this.getNumRows() != B.getNumRows() || this.getNumCols() != B.getNumCols()){
			throw new IllegalArgumentException("Matrix dimensions should match for element division!");
		}
		Matrix res = new Matrix(this.getNumRows(), this.getNumCols());
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				res.set(i, j, this.get(i, j) / B.get(i, j));
			}
		}
		return res;
	}

	/**
	 * Performs element-wise multiplication
	 * Note that B is a row vector, and the result is
	 * A_{ij} = \sum_{j=1}^n A_{ij} * B(i)
	 * @param B - the other row vector
	 * @return A_i * B
	 */
	public Matrix elementMultWithVector(Matrix B){
		if(this.getNumCols() != B.getNumCols()){
			throw new IllegalArgumentException("Matrix dimensions should match for element multiplication!");
		}
		Matrix res = new Matrix(this.getNumRows(), this.getNumCols());
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				res.set(i, j, this.get(i, j) * B.get(i));
			}
		}
		return res;
	}

	public void removeZeros(double val) {
		CommonOps_DSCC.removeZeros(this, val);
	}

	public void changeSign() {
		CommonOps_DSCC.changeSign(this, this);
	}

	public static void extract(Matrix src, int srcX0, int srcX1, int srcY0, int srcY1,
														 Matrix dst, int dstY0, int dstX0) {
		CommonOps_DSCC.extract(src, srcX0, srcX1, srcY0, srcY1, dst, dstY0, dstX0);
	}

	public static Matrix extractRows(Matrix A, int row0, int row1, Matrix out ) {
		if (out == null) {
			return new Matrix(CommonOps_DSCC.extractRows(A, row0, row1, out));
		} else {
			CommonOps_DSCC.extractRows(A, row0, row1, out);
			return out;
		}
	}

	public static Matrix extractColumn(Matrix A, int column, Matrix out ) {
		if (out == null) {
			return new Matrix(CommonOps_DSCC.extractColumn(A, column, out));
		} else {
			CommonOps_DSCC.extractColumn(A, column, out);
			return out;
		}
	}

	public static void extractDiag(Matrix A, Matrix outputB ) {
		CommonOps_DSCC.extractDiag(A, outputB);
	}

	public static Matrix concatColumns(Matrix left, Matrix right, Matrix out ) {
		if (out == null) {
			return new Matrix(CommonOps_DSCC.concatColumns(left, right, out));
		} else {
			CommonOps_DSCC.concatColumns(left, right, out);
			return out;
		}
	}

	public static Matrix concatRows(Matrix top, Matrix bottom, Matrix out) {
		if (out == null) {
			return new Matrix(CommonOps_DSCC.concatRows(top, bottom, out));
		} else {
			CommonOps_DSCC.concatRows(top, bottom, out);
			return out;
		}
	}

	public static Matrix diagMatrix(Matrix A, double[] values, int offset, int length) {
		return new Matrix(CommonOps_DSCC.diag(A, values, offset, length));
	}

	public static Matrix diag(double... values ) {
		return new Matrix(CommonOps_DSCC.diag(values));
	}

	public double[] toArray1D(){
		double[] array = new double[numRows*numCols];
		int k = 0;
		for (int i=0;i<numRows;i++){
			for(int j=0;j<numCols;j++){
				array[k] = this.get(i,j);
				k++;
			}
		}
		return array;
	}

	public List<Double> toList1D(){
		List<Double> list= new ArrayList<Double>();
		for (int i=0;i<numRows;i++){
			for(int j=0;j<numCols;j++){
				list.add(this.get(i,j));
			}
		}
		return list;
	}

	public double[][] toArray2D(){
		double[][] array = new double[numRows][numCols];
		for (int i=0;i<numRows;i++){
			for(int j=0;j<numCols;j++){
				array[i][j] = this.get(i,j);
			}
		}
		return array;
	}

	public Double toDouble() {
		return this.get(0, 0);
	}

	public List<List<Double>> toDoubleList(){
		List<List<Double>> array = new ArrayList<>();
		for (int i=0;i<numRows;i++){
			List<Double> row = new ArrayList<>();
			for(int j=0;j<numCols;j++){
				row.add(this.get(i,j));
			}
			array.add(row);
		}
		return array;
	}

	public static boolean solve(Matrix a, Matrix b, Matrix x ) {
		return CommonOps_DSCC.solve(a, b, x);
	}

	public double sumSubMatrix(int startRow, int endRow, int startCol, int endCol) {
		// endRow and endCol are EXCLUSIVE i.e. sums up to but not including that row/col
		double sum = 0;
		for (int i = startRow; i < endRow; i++) {
			for (int j = startCol; j < endCol; j++) {
				sum += this.get(i, j);
			}
		}
		return sum;
	}

	public Matrix sumRows(int startCol, int endCol) {
		// endCol is EXCLUSIVE i.e. sums up to but not including that col

		Matrix result = new Matrix(this.getNumRows(), 1);

		for (int i = 0; i < this.getNumRows(); i++) {
			double sum = 0;
			for (int j = startCol; j < endCol; j++) {
				sum += this.get(i, j);
			}
			result.set(i, 0, sum);
		}
		return result;
	}

	public Matrix sumCols(int startRow, int endRow) {
		// endRow is EXCLUSIVE i.e. sums up to but not including that row

		Matrix result = new Matrix(1, this.getNumCols());

		for (int i = 0; i < this.getNumCols(); i++) {
			double sum = 0;
			for (int j = startRow; j < endRow; j++) {
				sum += this.get(i, j);
			}
			result.set(0, i, sum);
		}
		return result;
	}

	public void ones() {
		int rows = this.getNumRows();
		int cols = this.getNumCols();

		if (this.getNumElements() != rows * cols) {
			throw new RuntimeException("Matrix is too small to fill with ones");
		}

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				this.set(i, j, 1);
			}
		}
	}

	// Kronecker sum of matrices A and B
	public Matrix krons(Matrix other) {

		DMatrixRMaj A = new DMatrixRMaj(this);
		DMatrixRMaj B = new DMatrixRMaj(other);
		DMatrixRMaj C = CommonOps_DDRM.kron(A, CommonOps_DDRM.identity(B.numRows), null);
		DMatrixRMaj D = CommonOps_DDRM.kron(CommonOps_DDRM.identity(A.numRows), B, null);
		DMatrixRMaj output = CommonOps_DDRM.add(C, D, null);
		return new Matrix(SimpleMatrix.wrap(output));
	}

	public Matrix createBlockDiagonal(Matrix matrix2) {
		int m1rows = this.getNumRows();
		int m2rows = 0;
		int m1cols = this.getNumCols();
		int m2cols = 0;

		if (matrix2 != null) {
			m2rows = matrix2.getNumRows();
			m2cols = matrix2.getNumCols();
		}

		Matrix output = new Matrix(m1rows + m2rows, m1cols + m2cols);
		for (int i = 0; i < m1rows; i++) {
			for (int j = 0; j < m1cols; j++) {
				output.set(i, j, this.get(i, j));
			}
		}
		for (int i = 0; i < m2rows; i++) {
			for (int j = 0; j < m2cols; j++) {
				output.set(i + m1rows, j + m1cols, matrix2.get(i, j));
			}
		}

		return output;
	}

	public static Matrix eye(int length) {
		return new Matrix(CommonOps_DSCC.identity(length));
	}

	public static Matrix ones(int rows, int cols) {
		Matrix e = new Matrix(rows,cols);
		e.ones();
		return e;
	}

	public Matrix expm() {

		DoubleMatrix inputToEXPM = new DoubleMatrix(this.numRows, this.numCols);
		for (int i = 0; i < this.numRows; i++) {
			for (int j = 0; j < this.numCols; j++) {
				inputToEXPM.put(i, j, this.get(i, j));
			}
		}

		DoubleMatrix outputFromEXPM = org.jblas.MatrixFunctions.expm(inputToEXPM);
		Matrix output = new Matrix(outputFromEXPM.rows, outputFromEXPM.columns);
		for (int i = 0; i < output.numRows; i++) {
			for (int j = 0; j < output.numCols; j++) {
				output.set(i, j, outputFromEXPM.get(i, j));
			}
		}

		return output;
	}

	public double elementMaxAbs() {
		return CommonOps_DSCC.elementMaxAbs(this);
	}

	public void scale(double scalar) {
		Matrix output = new Matrix(this);
		CommonOps_DSCC.scale(scalar, this, output);
		this.setTo(output);
	}

	public void scale(double scalar, Matrix output) {
		CommonOps_DSCC.scale(scalar, this, output);
	}

	/**
	 * Computes the matrix resulted from ceiling every member of the current matrix
	 * @return - the matrix obtained by ceiling the current matrix
	 */
	public Matrix ceil(){
		Matrix output = new Matrix(this);
		for(int i = 0; i < output.getNumRows(); i++){
			for(int j = 0; j < output.getNumCols(); j++){
				output.set(i, j, Math.ceil(output.get(i, j)));
			}
		}
		return output;
	}

	/**
	 * Equivalent to the colon operator in MATLAB: (:).
	 * @return - the Matrix in a column-major order, flattened as a column vector.
	 */
	public Matrix columnMajorOrder(){
		Matrix output = new Matrix(this.getNumRows() * this.getNumCols(), 1);
		int idx = 0;
		for(int j = 0; j < this.getNumCols(); j++){
			for(int i = 0; i < this.getNumRows(); i++){
				output.set(idx, 0, this.get(i, j));
				idx++;
			}
		}
		return output;
	}

	/**
	 * Checks if two matrices are equal
	 * @param m - the other matrix
	 * @return - true if the current matrix and the other matrix are equal, false otherwise
	 */
	public boolean isEqualTo(Matrix m){
		if(this.getNumRows() != m.getNumRows() || this.getNumCols() != m.getNumCols())
			return false;
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				if(this.get(i, j) != m.get(i, j))
					return false;
			}
		}
		return true;
	}

	/**
	 * Checks if the current matrix has a non-zero element
	 * @return - true if there is a non-zero element in the matrix, false otherwise
	 */
	public boolean any(){
		return this.nz_length > 0;
	}

	/**
	 * Computes the product of the elements of a row/column vector
	 * @return - the product of all the elements of the given vector
	 */
	public double prodVector(){
		if(this.getNumRows() != 1 && this.getNumCols() != 1){
			// Not a row/column vector
			throw new IllegalArgumentException("Argument should be a row/column vector");
		}
		if(this.nz_length < this.getNumRows() * this.getNumCols()){
			// Some elements are 0, so the product will automatically be 0
			return 0;
		}
		double product = 1;
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				product *= this.get(i, j);
			}
		}
		return product;
	}

	/**
	 * Computes the Euclidean norm of a matrix
	 * @return - the Euclidean norm of the given matrix
	 */
	public double norm(){
		double sum = 0;
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				double num = this.get(i, j);
				sum += num * num;
			}
		}
		return Math.sqrt(sum);
	}

	public static Matrix pow(Matrix a, int b){
		Matrix result = a.clone();
		for(int i=1;i<b;i++){
			result = result.mult(a);
		}
		return result;
	}


	public static double first_norm(Matrix a){
		Matrix b = a.clone();
		b.abs();
		double norm = 0;
		for(int i=0;i<a.numCols;i++){
			norm = Math.max(norm,a.sumCols(i));
		}
		return norm;
	}

	public void insert_sub_matrix(int start_row, int start_col, int end_row, int end_col, Matrix matrix_to_be_inserted){
		if(end_col-start_col!=matrix_to_be_inserted.numCols || end_row-start_row!=matrix_to_be_inserted.numRows){
			throw new RuntimeException("matrix_to_be_inserted doesn't fit");
		}
		for(int i=0;i<matrix_to_be_inserted.numRows;i++){
			for (int j=0;j<matrix_to_be_inserted.numCols;j++){
				this.set(start_row+i,start_col+j,matrix_to_be_inserted.get(i,j));
			}
		}
	}

	public static Matrix negative(Matrix a){
		Matrix b = a.clone();
		for(int i=0;i<b.nz_length;i++){
			b.nz_values[i] = -b.nz_values[i];
		}
		return b;
	}

	public static Matrix scale_mult(Matrix a, double n){
		Matrix b = a.clone();
		b.scale(n);
		return b;
	}

	public Matrix eigenvalue(){

		if(this.numCols!=this.numRows){
			throw new RuntimeException("Only square matrix can be eigen decomposited");
		}
		DMatrixRMaj matrix = new DMatrixRMaj(this.numRows ,this.numCols);
		DConvertMatrixStruct.convert(this,matrix);


		EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(matrix.numCols, true);
		eig.decompose(matrix);

		Matrix result = new Matrix(1,numCols,numCols);
		for(int i=0;i<numCols;i++){
			if(eig.getEigenvalue(i).isReal()) {
				result.set(i, eig.getEigenvalue(i).getReal());
			}else {
				result.set(i, Double.NaN);
			}
		}

		if(result.hasNaN()){
			System.out.println("Complex eigenvalues detected, return NAN");
		}
		return result;
	}

	public Matrix eigenvector(){
		if(this.numCols!=this.numRows){
			throw new RuntimeException("Only square matrix can be eigen decomposited");
		}
		DMatrixRMaj matrix = new DMatrixRMaj(this.numRows ,this.numCols);
		DConvertMatrixStruct.convert(this,matrix);


		EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(matrix.numCols, true);
		eig.decompose(matrix);

		Matrix result = new Matrix(numRows,numCols,numRows);
		for(int i=0;i<numCols;i++){
			if(eig.getEigenvalue(i).isReal()) {
				DMatrixRMaj vector = eig.getEigenVector(i);
				for (int j=0;j<numRows;j++){
					result.set(i,j,vector.get(j));
				}
			}else {
				for (int j=0;j<numRows;j++){
					result.set(i,j,Double.NaN);
				}
			}
		}

		if(result.hasNaN()){
			System.out.println("Complex eigenvector detected, return NAN");
		}
		return result;
	}



	public Matrix safe_mult(Matrix B){
		if(B.length()==1&&this.length()!=1){
			return Matrix.scale_mult(this,B.get(0));
		}else if(this.numCols==B.numRows){
			return this.mult(B);
		}else {
			throw new RuntimeException("matrix product of X and Y failed");
		}
	}

	public int rank(){
		DMatrixRMaj matrix = new DMatrixRMaj(this.numRows ,this.numCols);
		DConvertMatrixStruct.convert(this,matrix);
		return MatrixFeatures_DDRM.rank(matrix);
	}

	public Matrix element_power(double a){
		Matrix b = this.clone();
		for(int i=0;i<this.nz_length;i++){
			b.nz_values[i] = Math.pow(nz_values[i],a);
		}
		return b;
	}

	public static double inf_norm(Matrix a){
		Matrix b = a.clone();
		b.abs();
		return b.sumRows().elementMax();
	}

	public double det() {
		return CommonOps_DSCC.det(this);
	}

	public static Matrix cellsum(Map<Integer,Matrix> cell_array){
		if(cell_array.isEmpty()){
			throw new RuntimeException("Hashmap is empty");
		}

		Matrix result = cell_array.get(0).clone();
		for(int i=1;i<cell_array.size();i++){
			if(cell_array.get(i).numRows!=result.numRows||cell_array.get(i).numCols!=result.numCols){
				throw new RuntimeException("Matrix have incompatible sizes for add operation");
			}
			result = result.add(1,cell_array.get(i));
		}
		return result;
	}

	public void removeINF() {
		int offset = 0;
		for (int i = 0; i < this.numCols; i++) {
			for (int j=0;j<this.numRows;j++){
				if(Double.isInfinite(get(j,i))){
					set(j,i,0);
				}
			}
			this.col_idx[i + 1] -= offset;
		}
		this.nz_length -= offset;
	}

	public Matrix element_divide(Matrix b){
		if(numCols!=b.numCols || numRows!=b.numRows){
			throw new RuntimeException("Element divide function requires two Matrix have the same size");
		}
		Matrix result = new Matrix(numRows,numCols,numRows*numCols);
		for(int i=0;i<numRows;i++){
			for (int j=0;j<numCols;j++){
				result.set(i,j,get(i,j)/b.get(i,j));
			}
		}
		return result;
	}

	public static Matrix createLike (Matrix B){
		return new Matrix(B.createLike());
	}

	public boolean hasDuplicates(){
		HashSet<Double> values = new HashSet<>();
		for(int i=0;i<numRows;i++){
			for (int j=0;j<numCols;j++){
				if(values.contains(get(i,j))){
					return true;
				}
				values.add(get(i,j));
			}
		}
		return false;
	}

	public Matrix right_matrix_divide (Matrix b){
		return transpose().left_matrix_divide(b.transpose());
	}

	public Matrix left_matrix_divide(Matrix b){
		if(numRows==numCols){
			Matrix x = new Matrix(0,0,0);
			solve(this,b,x);
			return x;
		}else {
			throw new RuntimeException("left matrix divide of rectangular matrix haven't been implemented!");
		}
	}

	public Map<String,Matrix> QR_decomposition(){
		if(numRows!= numCols){
			throw new RuntimeException("Only square matrix can be decomposited");
		}
		QRSparseDecomposition<DMatrixSparseCSC> decomposition = new QrLeftLookingDecomposition_DSCC(null);
		decomposition.decompose(this.toDMatrixSparseCSC());
		Map<String,Matrix> result = new HashMap<>();
		result.put("Q",new Matrix(decomposition.getQ(null,false)));
		result.put("R",new Matrix(decomposition.getR(null,false)));
		return result;
	}

	public Map<String,Matrix> Schur_decomposition(String method,Integer it_){
		if(numRows!= numCols){
			throw new RuntimeException("Only square matrix can be decomposited");
		}
		if(this.nz_length==0){
			Map<String,Matrix> result = new HashMap<>();
			result.put("U",Matrix.eye(numRows));
			result.put("T",this.clone());
			return result;
		}
		int it = 1;
		if(it_ != null){
			it = it_;
		}
		Map<String,Matrix> result = new HashMap<>();
		if(method.equals("default")){
			Matrix A = clone();
			Matrix U = Matrix.eye(numRows);
			for (int i=0;i<it;i++){
				Matrix Q = A.QR_decomposition().get("Q");
				A = A.QR_decomposition().get("R").mult(Q);
				U = U.mult(Q);
			}
			result.put("T",A);
			result.put("U",U);
		}else {
			throw new RuntimeException("This method hasn't been implemented, use default instead");
		}
		return result;
	}

	public Map<String,Matrix> Schur_decomposition(){
		return Schur_decomposition("default",null);
	}

	public static Matrix lyap(Matrix A, Matrix B, Matrix C,Matrix D){
		return sylv(A,B,C);
	}

	public static Matrix sylv(Matrix A, Matrix B, Matrix C){
		int n = C.numCols;
		Map<String,Matrix> schur_decomposition_A = A.Schur_decomposition();
		Matrix ZA = schur_decomposition_A.get("U");
		Matrix TA = schur_decomposition_A.get("T");
		Matrix ZB;
		Matrix TB;
		String solver_direction;
		if(A.transpose().isEqualTo(B)){
			ZB = ZA;
			TB = TA.transpose();
			solver_direction = "backward";
		}else {
			Map<String,Matrix> schur_decomposition_B = B.Schur_decomposition();
			ZB = schur_decomposition_B.get("U");
			TB = schur_decomposition_B.get("T");
			solver_direction = "forward";
		}

		Matrix F = ZA.transpose().mult(C).mult(ZB);

//		if(TA.isDiag_withintol() && TB.isDiag_withintol()){
//			Matrix TA_diag = new Matrix(0,0,0);
//			Matrix TB_diag = new Matrix(0,0,0);
//			Matrix.extractDiag(TA,TA_diag);
//			Matrix.extractDiag(TB,TB_diag);
//			Matrix L = TA.compatible_sizes_add(TB.transpose());
//			for (int i=0;i<L.numRows;i++){
//				for (int j=0;j<L.numCols;j++){
//					L.set(i,j,-1/L.get(i,j));
//				}
//			}
//			return ZA.mult(L.elementMult(F,null)).mult(ZB.transpose());
//		}

		Matrix Y = new Matrix(C.numRows,C.numCols,C.numRows*C.numCols);
		Matrix P = new Matrix(0,0,0);
		Matrix.extractDiag(TA,P);

		if(solver_direction.equals( "backward")){
			for(int k=n-1;k>0;k--){
				Matrix rhs = Matrix.extractColumn(F,k,null).add(1,Y.mult(Matrix.extractColumn(TB,k,null)));
				for (int i=0;i<TA.numRows;i++){
					TA.set(i,i, P.get(i)+TB.get(k,k));
				}
				Y.insert_sub_matrix(0,k,Y.numRows,k+1,TA.left_matrix_divide(Matrix.scale_mult(rhs,-1)));
			}
		}else{
			for (int k=0;k<n;k++){
				Matrix rhs = Matrix.extractColumn(F,k,null).add(1,Y.mult(Matrix.extractColumn(TB,k,null)));
				for (int i=0;i<TA.numRows;i++){
					TA.set(i,i, P.get(i)+TB.get(k,k));
				}
				Matrix c = TA.left_matrix_divide(Matrix.scale_mult(rhs,-1));
				Y.insert_sub_matrix(0,k,Y.numRows,k+1,TA.left_matrix_divide(Matrix.scale_mult(rhs,-1)));
			}
		}

		Matrix result = ZA.mult(Y).mult(ZB.transpose());
		return result;
	}

	public boolean isDiag_withintol() {
		for(int i=0;i<numRows;i++){
			for (int j=0;j<numCols;j++){
				if(i!=j&&get(i,j)>GlobalConstants.FineTol){
					return false;
				}
			}
		}
		return true;
	}

	public Matrix compatible_sizes_add(Matrix b){
		if(numRows == b.numRows && numCols==b.numCols){
			return this.add(1,b);
		}

		if(this.length()==1){
			return b.elementIncrease(get(0));
		}

		if(b.length()==1){
			return elementIncrease(b.get(0));
		}

		if((numCols==1||numRows==1) && b.numCols!=1 && b.numRows!=1){
			return matrix_add_vector(b,this);
		}

		if(numCols!=1 && numRows!=1 && (b.numRows==1||b.numCols==1) ){
			return matrix_add_vector(this,b);
		}

		if(numCols == 1 && b.numRows==1 ){
			return col_vector_add_row_vector(this,b);
		}

		if(numRows==1 && b.numCols==1){
			return col_vector_add_row_vector(b,this);
		}

		throw new RuntimeException("The size of matrix are not compatible");
	}

	public static Matrix matrix_add_vector(Matrix matrix, Matrix vector){
		Matrix result = matrix.clone();
		if(vector.numCols==1 && vector.numRows == matrix.numRows){
			for(int i=0; i< vector.numRows;i++){
				result.row_increase(i, vector.get(i));
			}
		}else if (vector.numRows==1 && vector.numCols==matrix.numCols ){
			for (int i=0;i<vector.numCols;i++){
				result.col_increase(i, vector.get(i));
			}
		}else {
			throw new RuntimeException("The size of matrix and vector are not compatible");
		}
		return result;
	}

	public void row_increase(int row, double a){
		if(row>=0 && row<numRows) {
			for (int i = 0; i < numCols; i++) {
				set(row, i, get(row, i) + a);
			}
		}else {
			throw new RuntimeException("Row index out of range");
		}
	}

	public void col_increase(int col, double a){
		if(col>=0 && col<numCols) {
			for (int i = 0; i < numRows; i++) {
				set(i,col, get(i,col) + a);
			}
		}else {
			throw new RuntimeException("Column index out of range");
		}
	}

	public static Matrix col_vector_add_row_vector(Matrix col_vector, Matrix row_vector){
		Matrix result = new Matrix(row_vector.numRows,col_vector.numCols, row_vector.numRows* col_vector.numCols);
		for(int i=0;i<col_vector.numCols;i++){
			for(int j=0;j< row_vector.numRows;j++){
				result.set(i,j, col_vector.get(i)+row_vector.get(j));
			}
		}
		return result;
	}

	public Matrix kron(Matrix b){
		DMatrixRMaj A = new DMatrixRMaj(this);
		DMatrixRMaj B = new DMatrixRMaj(b);
		DMatrixRMaj C = CommonOps_DDRM.kron(A, B, null);
		return new Matrix(SimpleMatrix.wrap(C));
	}

	public void print() {
		System.out.println(this);
	}

	@Override
	public String toString(){
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < this.getNumRows(); i++){
			sb.append("|");
			for(int j = 0; j < this.getNumCols(); j++){
				sb.append(" " + this.get(i, j) + " ");
			}
			sb.append("|\n");
		}
		return sb.toString();
	}


	/**
	 * Removes the specified rows from the given matrix. If possible, use a HashSet for the collection of rows to
	 * improve performance.
	 * @param rows - the indices of the rows to be removed
	 */
	public void removeRows(Collection<Integer> rows){
		for(int r : rows){
			if(r < 0 || r >= this.getNumRows()){
				throw new IllegalArgumentException("Cannot remove rows that are outside of the current matrix");
			}
		}
		Matrix newMatrix = new Matrix(this.getNumRows() - rows.size(), this.getNumCols());
		int factor = 0;
		for(int i = 0; i < this.getNumRows(); i++){
			if(rows.contains(i)){
				factor++;
			} else {
				for(int j = 0; j < this.getNumCols(); j++){
					newMatrix.set(i - factor, j, this.get(i, j));
				}
			}
		}
		this.setTo(newMatrix);
	}

	/**
	 * Removes the specified columns from the given matrix. If possible, use a HashSet for the collection of columns to
	 * improve performance.
	 * @param cols - the indices of the columns to be removed
	 */
	public void removeCols(Collection<Integer> cols){
		for(int c : cols){
			if(c < 0 || c >= this.getNumCols()){
				throw new IllegalArgumentException("Cannot remove cols that are outside of the current matrix");
			}
		}
		Matrix newMatrix = new Matrix(this.getNumRows(), this.getNumCols() - cols.size());
		for(int i = 0; i < this.getNumRows(); i++){
			int factor = 0;
			for(int j = 0; j < this.getNumCols(); j++){
				if(cols.contains(j)){
					factor++;
				} else {
					newMatrix.set(i, j - factor, this.get(i, j));
				}
			}
		}
		this.setTo(newMatrix);
	}

	/**
	 * Concatenates the columns of two matrices
	 * @param other - the other matrix which will be concatenated to the current one
	 * @return - a new matrix containing the columns of the two matrices, concatenated
	 */
	public Matrix concatCols(Matrix other){
		if(this.getNumRows() != other.getNumRows()){
			throw new IllegalArgumentException("The two matrices must have the same number of rows in order to have their columns concatenated");
		}
		Matrix res = new Matrix(this.getNumRows(), this.getNumCols() + other.getNumCols());
		for(int i = 0; i < this.getNumRows(); i++){
			for(int j = 0; j < this.getNumCols(); j++){
				res.set(i, j, this.get(i, j));
			}
		}
		for(int i = 0; i < other.getNumRows(); i++){
			for(int j = 0; j < other.getNumCols(); j++){
				res.set(i, this.getNumCols() + j, other.get(i, j));
			}
		}
		return res;
	}

	/**
	 * @param row
	 * @param alpha
	 * @return sum_{j=0}^n |A{ij}|^alpha
	 */
	public double powerSumRows(int row, double alpha){
		double sum = 0;
		for(int i = 0; i < this.numCols; i++) {
			sum += Math.pow(Math.abs(this.get(row, i)),alpha);
		}
		return sum;
	}

	/**
	 * @param col
	 * @param alpha
	 * @return sum_{j=0}^n |A{ji}|^alpha
	 */
	public double powerSumCols(int col, double alpha){
		double sum = 0;
		for(int i = 0; i < this.numRows; i++) {
			sum += Math.pow(Math.abs(this.get(i, col)),alpha);
		}
		return sum;
	}
}
