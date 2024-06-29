package jline.util;

import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

/**
 * A class that extends EJML's Equation to use jline.util.Matrix
 */
public class MatrixEquation {
    Equation eq;

    public MatrixEquation() {
        eq = new Equation();
    }

    public MatrixEquation(Equation e) {
        eq = e;
    }

    public MatrixEquation(Object... args) {
        this();
        eq = new Equation();
        eq.alias(args);
    }

    public void alias(Object... args) {
        if (args.length % 2 == 1) {
            throw new RuntimeException("Even number of arguments expected");
        } else {
            for(int i = 0; i < args.length; i += 2) {
                if (args[i] instanceof Matrix) {
                    args[i] = ((Matrix)args[i]).toDMatrixSparseCSC();
                }
            }
            eq.alias(args);
        }
    }

    public MatrixEquation process(String equation) {
        return new MatrixEquation(eq.process(equation));
    }

    public Matrix lookupSimple(String token) {
        SimpleMatrix S = eq.lookupSimple(token);
        return new Matrix(S);
    }
}
