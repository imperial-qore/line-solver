/**
 * Method of Moments (MOM) solver for queueing network analysis.
 *
 * @since LINE 3.0
 */
package jline.lib.mom.solver;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import jline.util.Maths;

public class MomSolver {

    public MomSolverResult solve(RealMatrix L, int[] N, double[] Z) {
        int M = L.getRowDimension();
        int R = L.getColumnDimension();

        if (N.length != R) {
            throw new IllegalArgumentException("Population vector N must have length R=" + R);
        }
        if (Z.length != R) {
            throw new IllegalArgumentException("Think time vector Z must have length R=" + R);
        }

        int[] n = new int[R];

        double[] g = null;

        for (int r = 0; r < R; r++) {
            System.out.println("MOM: Processing class R*=" + (r + 1));

            LinearSystemMatrices matrices = setupLinearSystem(L, N, Z, r);

            if (r == 0) {
                g = new double[M + 1];
                for (int i = 0; i < M + 1; i++) g[i] = 1.0;
            } else {
                int prevSize = g.length;
                int currSize = (int) Maths.binomialCoeff(M + r - 1, r) * r;
                double[] newG = new double[currSize];

                int outer = (int) Maths.binomialCoeff(M + r - 2, r - 1);
                for (int i = 0; i < outer; i++) {
                    for (int s = 0; s < r - 1; s++) {
                        newG[i * r + s] = g[i * (r - 1) + s];
                    }
                    newG[i * r + r - 1] = g[prevSize - outer + i];
                }

                RealVector gVector = MatrixUtils.createRealVector(newG);
                double[] gk = blockSolve(M, r, matrices.C,
                        matrices.Cg.operate(gVector).toArray());

                g = new double[gk.length + newG.length];
                System.arraycopy(gk, 0, g, 0, gk.length);
                System.arraycopy(newG, 0, g, gk.length, newG.length);
            }

            for (int nr = 0; nr < N[r]; nr++) {
                n[r]++;
                System.out.println("Population state: " + java.util.Arrays.toString(n));

                RealVector gVec = matrices.Dr.operate(MatrixUtils.createRealVector(g));
                double[] G = new double[gVec.getDimension()];
                for (int i = 0; i < gVec.getDimension(); i++) {
                    G[i] = gVec.getEntry(i) / n[r];
                }

                double[] rhs;
                if (nr == 0) {
                    rhs = matrices.D.operate(MatrixUtils.createRealVector(g)).toArray();
                } else {
                    RealMatrix temp = matrices.D.subtract(
                            matrices.Cg.multiply(matrices.Dr).scalarMultiply(1.0 / n[r]));
                    rhs = temp.operate(MatrixUtils.createRealVector(g)).toArray();
                }

                double[] gk = blockSolve(M, r + 1, matrices.C, rhs);

                g = new double[gk.length + G.length];
                System.arraycopy(gk, 0, g, 0, gk.length);
                System.arraycopy(G, 0, g, gk.length, G.length);
            }
        }

        return computePerformanceMeasures(L, N, Z, g);
    }

    private LinearSystemMatrices setupLinearSystem(RealMatrix L, int[] N, double[] Z, int r) {
        int M = L.getRowDimension();

        int size = (int) Maths.binomialCoeff(M + r, r + 1) * (r + 1);
        int prevSize = (r > 0) ? (int) Maths.binomialCoeff(M + r - 1, r) * r : M + 1;

        RealMatrix C = MatrixUtils.createRealIdentityMatrix(size);
        RealMatrix Cg = MatrixUtils.createRealMatrix(size, prevSize);
        RealMatrix D = MatrixUtils.createRealMatrix(size, prevSize);
        RealMatrix Dr = MatrixUtils.createRealMatrix(prevSize, prevSize);

        return new LinearSystemMatrices(C, Cg, D, Dr);
    }

    private double[] blockSolve(int M, int r, RealMatrix C, double[] rhs) {
        org.apache.commons.math3.linear.DecompositionSolver solver = new LUDecomposition(C).getSolver();
        RealVector solution = solver.solve(MatrixUtils.createRealVector(rhs));
        return solution.toArray();
    }

    private MomSolverResult computePerformanceMeasures(RealMatrix L, int[] N, double[] Z, double[] g) {
        int M = L.getRowDimension();
        int R = L.getColumnDimension();

        RealMatrix X = MatrixUtils.createRealMatrix(M, R);
        RealMatrix Q = MatrixUtils.createRealMatrix(M, R);

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                X.setEntry(i, j, 1.0 / (Z[j] + 1.0 / L.getEntry(i, j)));
                Q.setEntry(i, j, (double) N[j] / M);
            }
        }

        return new MomSolverResult(X, Q, g);
    }
}
