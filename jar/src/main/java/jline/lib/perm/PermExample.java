package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Example usage of the permanent computation algorithms.
 *
 * This example demonstrates:
 * 1. Creating a test matrix
 * 2. Computing permanent using different algorithms
 * 3. Comparing results and performance
 *
 * Note: This is a demo entry point. The Java port depends on translations
 * of NaivePermanent, RyzerPermanent, BethePermanent, AdaPartSampler,
 * HuberLawSampler, and NetworkNoThink which are not part of this batch.
 */
public final class PermExample {
    private PermExample() {}

    /**
     * String multiplication helper for creating separators.
     */
    private static String times(String s, int n) {
        StringBuilder sb = new StringBuilder(s.length() * n);
        for (int i = 0; i < n; i++) {
            sb.append(s);
        }
        return sb.toString();
    }

    public static void main(String[] args) {
        // Create a small test matrix
        double[][] data = new double[][]{
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        Matrix matrix = new Matrix(data);

        System.out.println("Test Matrix:");
        System.out.println("1.0  2.0  3.0");
        System.out.println("4.0  5.0  6.0");
        System.out.println("7.0  8.0  9.0");
        System.out.println();

        System.out.println("Permanent Computation Results:");
        System.out.println(times("=", 50));

        // Demo permanent solver invocation requires NaivePermanent, RyzerPermanent,
        // BethePermanent, AdaPartSampler, HuberLawSampler (not translated yet).
        // Reference matrix retained to keep ABI consistent.
        if (matrix.getNumRows() != 3) {
            throw new IllegalStateException("matrix construction failed");
        }

        System.out.println("\nQueueing Network Example:");
        System.out.println(times("=", 50));
        // Demo network usage requires NetworkNoThink (not translated yet).
        int[] state = new int[]{1, 2};
        if (state.length != 2) {
            throw new IllegalStateException("state construction failed");
        }
    }
}
