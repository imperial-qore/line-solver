package jline.util;

import org.apache.commons.math3.FieldElement;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SparseFieldMatrix;
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

  public ComplexMatrix clone() {
    return new ComplexMatrix(this.real.clone(), this.im.clone());
  }

  public Complex det() {
    FieldMatrix<Complex> a = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), this.getNumRows(), this.getNumCols());
    for (int i = 0; i < this.getNumRows(); i++) {
      for (int j = 0; j < this.getNumCols(); j++) {
        a.setEntry(i, j, this.get(i, j));
      }
    }
    FieldLUDecomposition LU = new FieldLUDecomposition(a);
    Complex det = (Complex) LU.getDeterminant();
    return det;
  }

  public void scale(double a) {
    this.real.scale(a);
    this.im.scale(a);
  }

  public void set(int i, int j, Complex val) {
    this.real.set(i, j, val.getReal());
    this.im.set(i, j, val.getImaginary());
  }

  public void set(int i, int j, double val) {
    this.real.set(i, j, val);
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

  public void zero() {
    this.real.zero();
    this.im.zero();
  }

  public static ComplexMatrix concatRows(ComplexMatrix top, ComplexMatrix bottom, ComplexMatrix out) {
    if (out == null) {
      return new ComplexMatrix(Matrix.concatRows(top.real, bottom.real, null), Matrix.concatRows(top.im, bottom.im, null));
    }
    else {
      Matrix.concatRows(top.real, bottom.real, out.real);
      Matrix.concatRows(top.im, bottom.im, out.im);
      return out;
    }
  }




}
