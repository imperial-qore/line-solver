package jline.lib.mom;

import jline.lib.mom.solver.LinearSolver;
import jline.lib.mom.solver.MomSolver;
import jline.lib.mom.solver.MomSolverResult;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

public final class Main {
    private Main() {}

    public static void main(String[] args) {
        System.out.println("Testing MOM Solver...");

        // Test 1: Simple 1 station, 2 class network
        System.out.println("\n=== Test 1: Single Station, Two Classes ===");
        RealMatrix L1 = MatrixUtils.createRealMatrix(new double[][]{
                {1.0, 2.0}  // Service rates for class 1 and 2
        });
        int[] N1 = new int[]{5, 3};  // 5 customers of class 1, 3 of class 2
        double[] Z1 = new double[]{0.0, 0.0};  // No think time

        try {
            MomSolver solver1 = new MomSolver();
            MomSolverResult result1 = solver1.solve(L1, N1, Z1);

            System.out.println("MomSolver Results:");
            System.out.println("Throughput: " + Arrays.toString(result1.X.getRow(0)));
            System.out.println("Queue lengths: " + Arrays.toString(result1.Q.getRow(0)));
            System.out.println("Normalizing constants size: " + result1.G.length);
        } catch (Exception e) {
            System.out.println("MomSolver failed: " + e.getMessage());
            e.printStackTrace();
        }

        // Test 2: Two stations, two classes with LinearSolver
        System.out.println("\n=== Test 2: Two Stations, Two Classes (LinearSolver) ===");
        RealMatrix L2 = MatrixUtils.createRealMatrix(new double[][]{
                {10.0, 5.0},
                {8.0, 12.0}
        });
        int[] N2 = new int[]{10, 10};
        double[] Z2 = new double[]{1.0, 2.0};

        try {
            LinearSolver solver2 = new LinearSolver();
            MomSolverResult result2 = solver2.solve(L2, N2, Z2);

            System.out.println("LinearSolver Results:");
            for (int i = 0; i < result2.X.getRowDimension(); i++) {
                System.out.println("Station " + i + " throughput: " + Arrays.toString(result2.X.getRow(i)));
            }
            for (int i = 0; i < result2.Q.getRowDimension(); i++) {
                System.out.println("Station " + i + " queue length: " + Arrays.toString(result2.Q.getRow(i)));
            }
        } catch (Exception e) {
            System.out.println("LinearSolver failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Tests Complete ===");
    }
}
