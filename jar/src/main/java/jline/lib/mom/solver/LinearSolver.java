/**
 * Linear solver implementation for MOM using double precision arithmetic.
 *
 * @since LINE 3.0
 */
package jline.lib.mom.solver;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import jline.util.Maths;

public class LinearSolver {

    public MomSolverResult solve(RealMatrix L, int[] N, double[] Z) {
        int M = L.getRowDimension();
        int R = L.getColumnDimension();

        int[] n = new int[R];

        RealVector g = null;
        RealVector gr = null;

        for (int r = 0; r < R; r++) {
            System.out.println("Linear solver: Setting up linear system for R*=" + (r + 1));

            LinearSystemMatrices ls = SetupLinearSystem.setup(L, N, Z, r);

            if (r == 0) {
                double[] arr = new double[M + 1];
                for (int i = 0; i < M + 1; i++) arr[i] = 1.0;
                g = MatrixUtils.createRealVector(arr);
            } else {
                RealVector G = expandG(g, gr, M, r);

                RealVector Gk = solveDense(ls.C, ls.Cg.operate(G).mapMultiply(-1.0));

                g = combineVectors(Gk, G);
            }

            RealMatrix CgDr = ls.Cg.multiply(ls.Dr);

            for (int nr = 0; nr < N[r]; nr++) {
                n[r]++;

                RealVector G = ls.Dr.operate(g).mapDivide((double) n[r]);

                RealVector rhs = ls.D.subtract(CgDr.scalarMultiply(1.0 / n[r])).operate(g);

                RealVector Gk = solveDense(ls.C, rhs);

                g = combineVectors(Gk, G);
            }

            gr = g;
            n[r]++;

            RealVector G = ls.Dr.operate(g).mapDivide((double) n[r]);
            RealVector Gk = solveDense(ls.C, ls.D.operate(g));
            g = combineVectors(Gk, G);
        }

        return computePerformanceMeasures(L, N, Z, g, n);
    }

    private RealVector expandG(RealVector g, RealVector gr, int M, int r) {
        int numNetworks = (int) Maths.binomialCoeff(M + r - 2, r - 1);
        int newSize = numNetworks * r;
        RealVector G = MatrixUtils.createRealVector(new double[newSize]);

        for (int i = 0; i < numNetworks; i++) {
            for (int s = 0; s < r - 1; s++) {
                G.setEntry(i * r + s, g.getEntry(i * (r - 1) + s));
            }
            G.setEntry(i * r + r - 1, gr.getEntry(i));
        }

        return G;
    }

    private RealVector solveDense(RealMatrix C, RealVector b) {
        try {
            DecompositionSolver solver = new LUDecomposition(C).getSolver();
            return solver.solve(b);
        } catch (SingularMatrixException e) {
            SingularValueDecomposition svd = new SingularValueDecomposition(C);
            return svd.getSolver().solve(b);
        }
    }

    private RealVector combineVectors(RealVector v1, RealVector v2) {
        RealVector combined = MatrixUtils.createRealVector(new double[v1.getDimension() + v2.getDimension()]);
        combined.setSubVector(0, v1);
        combined.setSubVector(v1.getDimension(), v2);
        return combined;
    }

    private MomSolverResult computePerformanceMeasures(RealMatrix L, int[] N, double[] Z, RealVector g, int[] n) {
        int M = L.getRowDimension();
        int R = L.getColumnDimension();

        double G0 = g.getEntry(g.getDimension() - 1);

        RealMatrix X = MatrixUtils.createRealMatrix(M, R);
        RealMatrix Q = MatrixUtils.createRealMatrix(M, R);

        double totalPop = 0.0;
        for (int i : n) totalPop += i;

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double lambda = (double) N[r] / (Z[r] + totalPop / L.getEntry(i, r));
                X.setEntry(i, r, lambda);
                Q.setEntry(i, r, lambda / L.getEntry(i, r));
            }
        }

        double norm = X.getFrobeniusNorm();
        if (norm > 0) {
            X.scalarMultiply(G0 / norm);
            Q.scalarMultiply(G0 / norm);
        }

        return new MomSolverResult(X, Q, g.toArray());
    }
}
