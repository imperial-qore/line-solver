package jline.lang;

import jline.util.Matrix;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;

import java.util.List;


/**
 * Class representing a probabilistic routing matrix
 */
public class ClassSwitchMatrix extends Matrix {
    public ClassSwitchMatrix(double scalar) {
        super(scalar);
    }

    public ClassSwitchMatrix(int numRows, int numCols, int arrayLength) {
        super(numRows, numCols, arrayLength);
    }

    public ClassSwitchMatrix(int numRows, int numCols) {
        super(numRows, numCols);
    }

    public ClassSwitchMatrix(Matrix matrix) {
        super(matrix);
    }

    public ClassSwitchMatrix(DMatrixSparseCSC matrix) {
        super(matrix);
    }

    public ClassSwitchMatrix(List<Double> array) {
        super(array);
    }

    public ClassSwitchMatrix(SimpleMatrix matrix) {
        super(matrix);
    }

    public void set(JobClass jobclass1, JobClass jobclass2, double value) {
        this.set(jobclass1.getIndex()-1,jobclass2.getIndex()-1, value);
    }

}
