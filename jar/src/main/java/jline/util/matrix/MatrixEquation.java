/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

/**
 * A wrapper class that extends EJML's Equation functionality to work seamlessly with jline.util.matrix.Matrix objects.
 *
 * <p>This class provides a bridge between the EJML equation parsing and evaluation system and the
 * LINE matrix system. It automatically handles conversion between Matrix objects and EJML's internal
 * representations, allowing for intuitive mathematical expressions to be evaluated using LINE matrices.</p>
 *
 * <p>Usage example:</p>
 * <pre>
 * Matrix A = new Matrix(3, 3);
 * Matrix B = new Matrix(3, 3);
 * MatrixEquation eq = new MatrixEquation();
 * eq.alias("A", A, "B", B);
 * MatrixEquation result = eq.process("C = A * B + I");
 * Matrix C = result.lookupSimple("C");
 * </pre>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public class MatrixEquation {
    /**
     * The underlying EJML equation processor
     */
    Equation eq;

    /**
     * Creates a new MatrixEquation with an empty equation context.
     */
    public MatrixEquation() {
        eq = new Equation();
    }

    /**
     * Creates a MatrixEquation wrapping an existing EJML Equation.
     *
     * @param e the EJML Equation to wrap
     */
    public MatrixEquation(Equation e) {
        eq = e;
    }

    /**
     * Creates a new MatrixEquation and immediately aliases the provided variables.
     *
     * @param args variable name and value pairs (must be even number of arguments)
     * @throws RuntimeException if an odd number of arguments is provided
     */
    public MatrixEquation(Object... args) {
        this();
        eq = new Equation();
        eq.alias(args);
    }

    /**
     * Creates aliases for variables that can be used in matrix equations.
     * Arguments should be provided in pairs: variable name followed by variable value.
     * Matrix objects are automatically converted to EJML's sparse matrix format.
     *
     * @param args variable name and value pairs (name1, value1, name2, value2, ...)
     * @throws RuntimeException if an odd number of arguments is provided
     */
    public void alias(Object... args) {
        if (args.length % 2 == 1) {
            throw new RuntimeException("Even number of arguments expected");
        } else {
            for (int i = 0; i < args.length; i += 2) {
                if (args[i] instanceof Matrix) {
                    args[i] = ((Matrix) args[i]).toDMatrixSparseCSC();
                }
            }
            eq.alias(args);
        }
    }

    /**
     * Retrieves a matrix result from the equation context by variable name.
     * The returned matrix is automatically converted from EJML's SimpleMatrix format
     * to the LINE Matrix format.
     *
     * @param token the variable name to look up
     * @return the matrix associated with the given variable name
     */
    public Matrix lookupSimple(String token) {
        SimpleMatrix S = eq.lookupSimple(token);
        return new Matrix(S);
    }

    /**
     * Processes a mathematical equation string and returns a new MatrixEquation containing the results.
     * The equation can reference any variables that have been aliased in this context.
     *
     * @param equation the mathematical equation to process (e.g., "C = A * B + I")
     * @return a new MatrixEquation containing the processed equation results
     */
    public MatrixEquation process(String equation) {
        return new MatrixEquation(eq.process(equation));
    }
}
