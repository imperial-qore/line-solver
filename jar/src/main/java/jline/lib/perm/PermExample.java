package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Example usage of the permanent computation algorithms.
 *
 * This example demonstrates:
 * 1. Creating a test matrix
 * 2. Computing permanent using different algorithms
 * 3. Comparing results and performance
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

    /**
     * Solve with one algorithm and print its value and elapsed time.
     */
    private static void report(String name, PermSolver solver) {
        solver.solve();
        System.out.println(String.format("%-16s %14.6g  %5d ms", name, solver.getValue(), solver.getTime()));
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

        System.out.println("Permanent Computation Results (exact value 450):");
        System.out.println(times("=", 50));

        report("Permanent", new Permanent(matrix));
        report("Naive", new NaivePermanent(matrix));
        report("Ryzer graycode", new RyzerPermanent(matrix, "graycode"));
        report("Ryzer naive", new RyzerPermanent(matrix, "naive"));
        report("Bethe", new BethePermanent(matrix));
        report("Heuristic", new HeuristicPermanent(matrix));
        report("Huber-Law", new HuberLawSampler(matrix, 0.1, 0.000001, 0.1, "sample", 1000, 30000L, false));
        report("AdaPart", new AdaPartSampler(matrix, 100, 30000L, 450, "classic", false));

        System.out.println("\nQueueing Network Example:");
        System.out.println(times("=", 50));

        // Two queues, one class, three jobs: marginal probability of state (1,2)
        double[][] demands = new double[][]{
            {0.4, 0.6},
            {0.9, 0.2}
        };
        NetworkNoThink network = new NetworkNoThink(2, 1, new int[]{3}, demands, false);
        int[] state = new int[]{1, 2};
        NetworkNoThink.MarginalResult exact = network.marginal(new RyzerPermanent(matrix), state, false);
        NetworkNoThink.MarginalResult bethe = network.marginal(new BethePermanent(matrix), state, false);
        System.out.println(String.format("state (1,2) exact marginal %14.6g", exact.probability));
        System.out.println(String.format("state (1,2) Bethe marginal %14.6g", bethe.probability));
    }
}
