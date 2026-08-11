package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Abstract base class for permanent computation solvers.
 *
 * This abstract class provides the common interface for all permanent computation algorithms.
 * All solvers inherit from this class and implement the compute() method for their specific algorithm.
 * The base class provides the solve() method that measures time and memory usage during computation.
 */
public abstract class PermSolver {

    public final Matrix matrix;
    public final int n;
    public double value = 0.0;
    public long time = 0L;
    public long memory = 0L;

    public PermSolver(Matrix matrix) {
        this.matrix = matrix;
        this.n = matrix.getNumRows();
    }

    public Matrix getMatrix() {
        return matrix;
    }

    public int getN() {
        return n;
    }

    public double getValue() {
        return value;
    }

    public long getTime() {
        return time;
    }

    public long getMemory() {
        return memory;
    }

    /**
     * Compute the permanent or approximation for the given matrix.
     * This method must be implemented by all concrete solver classes.
     * The result should be stored in the value field.
     */
    public abstract void compute();

    /**
     * Measure time and memory usage of the compute method.
     * This method wraps the compute() call with performance monitoring.
     */
    public void solve() {
        // Get initial memory usage
        Runtime runtime = Runtime.getRuntime();
        runtime.gc(); // Suggest garbage collection before measurement
        long memoryBefore = runtime.totalMemory() - runtime.freeMemory();

        // Measure computation time
        long startTime = System.currentTimeMillis();

        compute();

        long endTime = System.currentTimeMillis();
        this.time = endTime - startTime;

        // Measure memory usage after computation
        long memoryAfter = runtime.totalMemory() - runtime.freeMemory();
        this.memory = memoryAfter - memoryBefore;
    }
}
