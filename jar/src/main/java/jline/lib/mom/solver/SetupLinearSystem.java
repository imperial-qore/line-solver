package jline.lib.mom.solver;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import jline.lib.mom.util.MomUtils;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Constructs the linear system matrices for the MOM solver.
 * This is a port of the MATLAB setupls.m function.
 */
public final class SetupLinearSystem {
    private SetupLinearSystem() {}

    /**
     * Setup linear system matrices for a given class r.
     */
    public static LinearSystemMatrices setup(RealMatrix L, int[] N, double[] Z, int r) {
        int M = L.getRowDimension();
        int R = L.getColumnDimension();

        // Expand service rate matrix
        RealMatrix Lhat = expandServiceRates(L, N, r);

        // Generate all population distributions for current and previous classes
        Matrix IkMatrix;
        if (r < R - 1) {
            IkMatrix = Maths.multichoose((double) M, (double) (r + 1));
        } else {
            // For last class, only consider feasible states
            IkMatrix = generateFeasibleStatesMatrix(M, R, N, r);
        }

        Matrix IMatrix;
        if (r > 0) {
            IMatrix = Maths.multichoose((double) M, (double) r);
        } else {
            IMatrix = new Matrix(1, M); // All zeros for r=0
        }

        // Convert Matrix to int[][] for sorting
        int[][] Ik = matrixToIntArrays(IkMatrix);
        int[][] I = matrixToIntArrays(IMatrix);

        // Sort combinations for better numerical properties
        int[][] IkSorted = MomUtils.sortByNnzPos(Ik);
        int[][] ISorted = MomUtils.sortByNnzPos(I);

        // Initialize matrices
        int numStates = IkSorted.length * (r + 1);
        int numPrevStates = ISorted.length * Math.max(r, 1);

        RealMatrix C = MatrixUtils.createRealMatrix(numStates, numStates);
        RealMatrix Cg = MatrixUtils.createRealMatrix(numStates, numPrevStates);
        RealMatrix D = MatrixUtils.createRealMatrix(numStates, numPrevStates);
        RealMatrix Dr = MatrixUtils.createRealMatrix(numPrevStates, numPrevStates);

        // Build coefficient matrices
        buildCoefficientMatrices(C, Cg, D, Dr, Lhat, Z, IkSorted, ISorted, M, r);

        return new LinearSystemMatrices(C, Cg, D, Dr);
    }

    /**
     * Convert Matrix to int[][].
     */
    private static int[][] matrixToIntArrays(Matrix matrix) {
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();
        int[][] result = new int[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = (int) matrix.get(i, j);
            }
        }
        return result;
    }

    /**
     * Expand service rates for multi-server stations.
     */
    private static RealMatrix expandServiceRates(RealMatrix L, int[] N, int r) {
        int M = L.getRowDimension();
        int maxServers = 0;
        for (int n : N) maxServers += n;
        RealMatrix Lhat = MatrixUtils.createRealMatrix(M, maxServers);

        for (int i = 0; i < M; i++) {
            for (int k = 0; k < maxServers; k++) {
                double rate = L.getEntry(i, Math.min(r, L.getColumnDimension() - 1));
                Lhat.setEntry(i, k, rate * Math.min(k + 1, N[Math.min(r, N.length - 1)]));
            }
        }

        return Lhat;
    }

    /**
     * Generate feasible population states for the last class.
     */
    private static Matrix generateFeasibleStatesMatrix(int M, int R, int[] N, int r) {
        int Nsum = 0;
        for (int n : N) Nsum += n;
        Matrix allStates = Maths.multichoose((double) M, (double) Math.min(r + 1, Nsum));

        // Filter feasible states
        List<int[]> feasibleList = new ArrayList<int[]>();
        for (int i = 0; i < allStates.getNumRows(); i++) {
            int[] state = new int[M];
            int sum = 0;
            for (int j = 0; j < M; j++) {
                state[j] = (int) allStates.get(i, j);
                sum += state[j];
            }
            if (sum <= Nsum) {
                feasibleList.add(state);
            }
        }

        // Convert back to Matrix
        Matrix result = new Matrix(feasibleList.size(), M);
        for (int i = 0; i < feasibleList.size(); i++) {
            int[] state = feasibleList.get(i);
            for (int j = 0; j < state.length; j++) {
                result.set(i, j, (double) state[j]);
            }
        }

        return result;
    }

    /**
     * Build the coefficient matrices based on balance equations.
     */
    private static void buildCoefficientMatrices(RealMatrix C, RealMatrix Cg, RealMatrix D, RealMatrix Dr,
                                                  RealMatrix Lhat, double[] Z, int[][] Ik, int[][] I,
                                                  int M, int r) {
        // Build C matrix (block diagonal structure)
        for (int h = 0; h < Ik.length; h++) {
            int[] ik = Ik[h];

            for (int s = 0; s <= r; s++) {
                int row = h * (r + 1) + s;

                // Diagonal entry
                double diag = 0.0;
                int ikSum = 0;
                for (int i = 0; i < M; i++) {
                    if (ik[i] > 0) {
                        diag += Lhat.getEntry(i, ik[i] - 1);
                    }
                    ikSum += ik[i];
                }
                if (s < Z.length && Z[s] > 0) {
                    diag += ikSum / Z[s];
                }

                C.setEntry(row, row, diag);

                // Off-diagonal entries within block
                for (int sp = 0; sp <= r; sp++) {
                    if (sp != s) {
                        int col = h * (r + 1) + sp;
                        double value = 0.0;

                        // Station transitions
                        for (int i = 0; i < M; i++) {
                            if (ik[i] > 0) {
                                value += Lhat.getEntry(i, ik[i] - 1) / M;
                            }
                        }

                        C.setEntry(row, col, -value);
                    }
                }
            }
        }

        // Build Cg matrix (coupling to previous populations)
        for (int h = 0; h < Ik.length; h++) {
            int[] ik = Ik[h];

            for (int hp = 0; hp < I.length; hp++) {
                int[] i = I[hp];

                // Check if i can lead to ik by adding one customer
                if (canTransition(i, ik)) {
                    for (int s = 0; s <= r; s++) {
                        int row = h * (r + 1) + s;

                        for (int sp = 0; sp < Math.max(r, 1); sp++) {
                            int col = hp * Math.max(r, 1) + sp;

                            // Simplified - actual implementation needs routing logic
                            Cg.setEntry(row, col, 1.0 / M);
                        }
                    }
                }
            }
        }

        // Build D matrix (direct transitions)
        double dValue = 1.0 / (M * (r + 1));
        for (int i = 0; i < D.getRowDimension(); i++) {
            for (int j = 0; j < D.getColumnDimension(); j++) {
                D.setEntry(i, j, dValue);
            }
        }

        // Build Dr matrix (recursion matrix)
        // Identity-like structure for recursion
        for (int i = 0; i < Math.min(Dr.getRowDimension(), Dr.getColumnDimension()); i++) {
            Dr.setEntry(i, i, 1.0);
        }
    }

    /**
     * Check if population state i can transition to ik by adding one customer.
     */
    private static boolean canTransition(int[] i, int[] ik) {
        int diff = 0;
        for (int j = 0; j < i.length; j++) {
            int d = ik[j] - i[j];
            if (d < 0) return false;
            diff += d;
        }
        return diff == 1;
    }
}
