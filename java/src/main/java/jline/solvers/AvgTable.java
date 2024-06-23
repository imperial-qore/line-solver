package jline.solvers;

import jline.util.Matrix;

import java.util.ArrayList;
import java.util.List;

public abstract class AvgTable {

    SolverOptions options;
    ArrayList<List<Double>> T;

    public AvgTable(ArrayList<List<Double>> table) {
        this.T = table;
    }
    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public List<Double> get(int col) {
        return this.T.get(col);
    }

    public abstract void print();

    public Matrix toMatrix() {
        return new Matrix(T);
    }

}
