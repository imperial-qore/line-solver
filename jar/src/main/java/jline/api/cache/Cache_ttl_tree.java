package jline.api.cache;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.util.FastMath;

import jline.api.mc.Dtmc_solve;
import jline.util.matrix.Matrix;

/**
 * Tree-based TTL cache analysis implementation for the LINE solver framework.
 */
public final class Cache_ttl_tree {
    private Cache_ttl_tree() {}

    public static Matrix cache_ttl_tree(Matrix[] lambda, Matrix[][] R, Matrix m) {
        return cache_ttl_tree(lambda, R, m, null);
    }

    public static Matrix cache_ttl_tree(final Matrix[] lambda, final Matrix[][] R, final Matrix m, Long seed) {
        final Random random = (seed != null) ? new Random(seed) : new Random();
        final int n = lambda[0].getNumRows();
        final int h = lambda[0].getNumCols() - 1;

        double rangeLeft = 0.0;
        double rangeRight = 10.0;
        double[] initialGuess = new double[h];
        for (int i = 0; i < h; i++) initialGuess[i] = rangeLeft + (rangeRight - rangeLeft) * random.nextDouble();

        MultivariateFunction objectiveFunction = new MultivariateFunction() {
            @Override
            public double value(double[] x) {
                Object[] result = ttlTreeTime(x, lambda, R, m, n, h);
                double[] f = (double[]) result[0];
                double sum = 0.0;
                for (double v : f) sum += v * v;
                return sum;
            }
        };

        MultivariateOptimizer optimizer = new PowellOptimizer(1e-6, 1e-6);
        PointValuePair result = optimizer.optimize(
                new MaxEval(100000),
                new ObjectiveFunction(objectiveFunction),
                GoalType.MINIMIZE,
                new InitialGuess(initialGuess));

        double[] optimalListTime = result.getPoint();
        Object[] retT = ttlTreeTime(optimalListTime, lambda, R, m, n, h);
        return (Matrix) retT[1];
    }

    private static Object[] ttlTreeTime(double[] x, Matrix[] lambda, Matrix[][] R, Matrix m, int n, int h) {
        Matrix steadyStateProb = new Matrix(n, h + 1);
        Matrix randProb = new Matrix(n, h + 1);
        Matrix avgTime = new Matrix(n, h + 1);
        double[] cdiff = new double[h];
        double[] capa = new double[h];
        double[] rpDenominator = new double[n];

        for (int i = 0; i < n; i++) {
            Matrix transMatrix = new Matrix(h + 1, h + 1);
            for (int j = 0; j < h + 1; j++) {
                List<Integer> leafNodes = new ArrayList<Integer>();
                for (int k = 0; k < h + 1; k++) {
                    if (R[0][i].get(j, k) != 0.0) leafNodes.add(k);
                }
                for (int k : leafNodes) {
                    if (j == 0) {
                        transMatrix.set(j, k, R[0][i].get(j, k));
                    } else {
                        if (j - 1 < x.length) {
                            transMatrix.set(j, k, (1 - FastMath.exp(-lambda[0].get(i, j) * x[j - 1])) * R[0][i].get(j, k));
                        }
                    }
                    if (j != k && k > 0 && k - 1 < x.length) {
                        transMatrix.set(k, j, FastMath.exp(-lambda[0].get(i, k) * x[k - 1]));
                    }
                }
            }
            List<Integer> missConnection = new ArrayList<Integer>();
            for (int row = 0; row < h + 1; row++) {
                boolean hasConnection = false;
                for (int col = 0; col < h + 1; col++) {
                    if (transMatrix.get(row, col) != 0.0) { hasConnection = true; break; }
                }
                if (!hasConnection) missConnection.add(row);
            }
            List<Integer> dtChain = new ArrayList<Integer>();
            for (int rIdx = 0; rIdx < h + 1; rIdx++) if (!missConnection.contains(rIdx)) dtChain.add(rIdx);

            Matrix reducedTransMatrix = new Matrix(dtChain.size(), dtChain.size());
            for (int a = 0; a < dtChain.size(); a++) {
                for (int b = 0; b < dtChain.size(); b++) {
                    reducedTransMatrix.set(a, b, transMatrix.get(dtChain.get(a), dtChain.get(b)));
                }
            }
            Matrix dtmcProb = Dtmc_solve.dtmc_solve(reducedTransMatrix);

            for (int a = 0; a < dtChain.size(); a++) {
                int originalState = dtChain.get(a);
                steadyStateProb.set(i, originalState, dtmcProb.get(0, a));
                if (originalState > 0 && originalState - 1 < x.length) {
                    avgTime.set(i, originalState, (1 - FastMath.exp(-lambda[0].get(i, originalState) * x[originalState - 1])) / lambda[0].get(i, originalState));
                } else {
                    avgTime.set(i, originalState, 1.0 / lambda[0].get(i, originalState));
                }
                rpDenominator[i] += steadyStateProb.get(i, originalState) * avgTime.get(i, originalState);
            }
            for (int a = 0; a < dtChain.size(); a++) {
                int originalState = dtChain.get(a);
                if (rpDenominator[i] > 0) {
                    randProb.set(i, originalState, steadyStateProb.get(i, originalState) * avgTime.get(i, originalState) / rpDenominator[i]);
                }
            }
        }

        for (int l = 0; l < h; l++) {
            capa[l] = 0.0;
            for (int i = 0; i < n; i++) capa[l] += randProb.get(i, l + 1);
            cdiff[l] = m.get(0, l) - capa[l];
        }
        return new Object[]{cdiff, randProb, cdiff};
    }

    /** Cache ttl tree algorithms. */
    public static final class CacheTtlTreeAlgo {}
}
