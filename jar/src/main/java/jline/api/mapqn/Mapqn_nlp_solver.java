/**
 * NLP Solver Utility for Linearly-Constrained Nonlinear Optimization.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

public final class Mapqn_nlp_solver {
    private Mapqn_nlp_solver() {}

    public interface ObjectiveFn {
        double apply(double[] x);
    }

    public static double[] solve(ObjectiveFn objective, int nVars,
                                 double[][] Aeq, double[] beq,
                                 double[][] Aub, double[] bub,
                                 double[] lb, double[] ub, double[] x0) {
        return solve(objective, nVars, Aeq, beq, Aub, bub, lb, ub, x0, 20, 50000);
    }

    /**
     * Solve a linearly-constrained NLP problem using Augmented Lagrangian + BOBYQA.
     */
    public static double[] solve(final ObjectiveFn objective, final int nVars,
                                 final double[][] Aeq, final double[] beq,
                                 final double[][] Aub, final double[] bub,
                                 final double[] lb, final double[] ub, final double[] x0,
                                 final int maxIter, final int maxEval) {
        final int numEq = (Aeq != null) ? Aeq.length : 0;
        final int numIneq = (Aub != null) ? Aub.length : 0;

        final double[] lambdaEq = new double[numEq];
        final double[] lambdaIneq = new double[numIneq];
        double rho = 1.0;
        double rhoMax = 1e6;
        double rhoMult = 10.0;

        double[] xCurrent = x0.clone();

        for (int iter = 0; iter < maxIter; iter++) {
            final double finalRho = rho;
            MultivariateFunction augLagrangian = new MultivariateFunction() {
                @Override
                public double value(double[] x) {
                    double fval = objective.apply(x);

                    if (Aeq != null && beq != null) {
                        for (int i = 0; i < numEq; i++) {
                            double hi = -beq[i];
                            for (int j = 0; j < nVars; j++) {
                                hi += Aeq[i][j] * x[j];
                            }
                            fval += lambdaEq[i] * hi + (finalRho / 2.0) * hi * hi;
                        }
                    }

                    if (Aub != null && bub != null) {
                        for (int i = 0; i < numIneq; i++) {
                            double gi = -bub[i];
                            for (int j = 0; j < nVars; j++) {
                                gi += Aub[i][j] * x[j];
                            }
                            double shifted = Math.max(0.0, lambdaIneq[i] / finalRho + gi);
                            fval += (finalRho / 2.0) * shifted * shifted;
                        }
                    }

                    return fval;
                }
            };

            int npt = Math.min(2 * nVars + 1, (nVars + 1) * (nVars + 2) / 2);
            double initialRadius = 0.1;
            double stoppingRadius = 1e-8;

            try {
                BOBYQAOptimizer optimizer = new BOBYQAOptimizer(npt, initialRadius, stoppingRadius);
                PointValuePair result = optimizer.optimize(
                        new MaxEval(maxEval),
                        new ObjectiveFunction(augLagrangian),
                        GoalType.MINIMIZE,
                        new SimpleBounds(lb, ub),
                        new InitialGuess(xCurrent)
                );
                xCurrent = result.getPoint();
            } catch (Exception e) {
                // keep current point
            }

            if (Aeq != null && beq != null) {
                for (int i = 0; i < numEq; i++) {
                    double hi = -beq[i];
                    for (int j = 0; j < nVars; j++) {
                        hi += Aeq[i][j] * xCurrent[j];
                    }
                    lambdaEq[i] += rho * hi;
                }
            }

            if (Aub != null && bub != null) {
                for (int i = 0; i < numIneq; i++) {
                    double gi = -bub[i];
                    for (int j = 0; j < nVars; j++) {
                        gi += Aub[i][j] * xCurrent[j];
                    }
                    lambdaIneq[i] = Math.max(0.0, lambdaIneq[i] + rho * gi);
                }
            }

            double maxViolation = 0.0;
            if (Aeq != null && beq != null) {
                for (int i = 0; i < numEq; i++) {
                    double hi = -beq[i];
                    for (int j = 0; j < nVars; j++) {
                        hi += Aeq[i][j] * xCurrent[j];
                    }
                    maxViolation = Math.max(maxViolation, Math.abs(hi));
                }
            }
            if (Aub != null && bub != null) {
                for (int i = 0; i < numIneq; i++) {
                    double gi = -bub[i];
                    for (int j = 0; j < nVars; j++) {
                        gi += Aub[i][j] * xCurrent[j];
                    }
                    maxViolation = Math.max(maxViolation, Math.max(0.0, gi));
                }
            }

            if (maxViolation < 1e-6) break;

            rho = Math.min(rho * rhoMult, rhoMax);
        }

        return xCurrent;
    }
}
