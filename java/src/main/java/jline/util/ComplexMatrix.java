package jline.util;

import org.apache.commons.math3.complex.Complex;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class ComplexMatrix {
  public Matrix real;
  public Matrix im;

  public ComplexMatrix(int i, int j) {
    this.real = new Matrix(i, j);
    this.im = new Matrix(i, j);
  }
  public ComplexMatrix(Matrix real, Matrix im) {
    assert(real.getNumRows() == im.getNumRows() && real.getNumCols() == im.getNumCols());
    this.real = real;
    this.im = im;
  }

  public ComplexMatrix(DMatrixSparseCSC matrix_real, DMatrixSparseCSC matrix_im) {
    this.real = new Matrix(matrix_real);
    this.im = new Matrix(matrix_im);

  }

  public ComplexMatrix(Matrix real) {
    this.real = real;
    this.im = new Matrix(real.getNumRows(), real.getNumCols());
    this.im.zero();
  }

  public void set(int i, int j, Complex val) {
    this.real.set(i, j, val.getReal());
    this.im.set(i, j, val.getImaginary());
  }

  public void set(int idx, Complex val) {
    this.real.set(idx, val.getReal());
    this.im.set(idx, val.getImaginary());
  }

  public Complex get(int i, int j) {
    return new Complex(this.real.get(i, j), this.im.get(i, j));
  }

  public Complex get(int idx) {
    return new Complex(this.real.get(idx), this.im.get(idx));
  }

  public int getNumRows() {
    return this.real.getNumRows();
  }

  public int getNumCols() {
    return this.real.getNumCols();
  }

  public int getNumElements() {
    return this.real.getNumElements();
  }

  public boolean isEmpty() {
    return this.real.isEmpty() && this.im.isEmpty();
  }

  public ComplexMatrix sumRows() {
    return new ComplexMatrix(this.real.sumRows(), this.im.sumRows());
  }

  public static ComplexMatrix extractRows(ComplexMatrix A, int row0, int row1, ComplexMatrix out) {
    if (out == null) {
      return new ComplexMatrix(CommonOps_DSCC.extractRows(A.real.data, row0, row1, null), CommonOps_DSCC.extractRows(A.im.data, row0, row1, null));
    } else {
      CommonOps_DSCC.extractRows(A.real.data, row0, row1, out.real.data);
      CommonOps_DSCC.extractRows(A.im.data, row0, row1, out.im.data);
      return out;
    }

  }




}
