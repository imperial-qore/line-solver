package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;

/**
 * Implementation of bnd_mvaversion.mod linear program.
 */
public final class Mapqn_bnd_lr_mva {
    private Mapqn_bnd_lr_mva() {}

    public static Mapqn_solution solve(MVAVersionParameters params, int objectiveQueue, int objectiveLevel) {
        params.validate();
        if (objectiveQueue < 1 || objectiveQueue > params.M) {
            throw new IllegalArgumentException("Objective queue must be in range 1.." + params.M);
        }
        if (objectiveLevel < 1 || objectiveLevel > params.K) {
            throw new IllegalArgumentException("Objective level must be in range 1.." + params.K);
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();
        registerVariables(model, params);
        addConstraints(model, params);

        String objectiveVarName = "UN_" + objectiveQueue + "_" + objectiveLevel;
        double[] objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName);
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        PointValuePair solution = solver.optimize(objectiveFunction, constraintSet, GoalType.MAXIMIZE);

        return new Mapqn_solution(solution.getValue(), extractVariableValues(model, solution.getPoint()));
    }

    private static void registerVariables(Mapqn_lpmodel model, MVAVersionParameters params) {
        int M = params.M;
        int K = params.K;
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K; k++) model.addVariable("UN_" + i + "_" + k);
        }
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K; k++) model.addVariable("QN_" + i + "_" + k);
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K; k++) {
                for (int i = 1; i <= M; i++) model.addVariable("B_" + j + "_" + k + "_" + i);
            }
        }
    }

    private static void addConstraints(Mapqn_lpmodel model, MVAVersionParameters params) {
        int M = params.M;
        int N = params.N;
        int K = params.K;

        // QNB
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K; k++) {
                for (int j = 1; j <= M; j++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                            .addTerm("QN_" + i + "_" + k, 1.0)
                            .addTerm("B_" + j + "_" + k + "_" + i, -1.0);
                    model.addConstraint(cb.geq(0.0));
                }
            }
        }

        // UMAX
        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K; k++) cb.addTerm("UN_" + i + "_" + k, 1.0);
            model.addConstraint(cb.leq(1.0));
        }

        // POPCONSTR
        Mapqn_lpmodel.LinearConstraintBuilder popConstraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K; k++) popConstraint.addTerm("QN_" + i + "_" + k, 1.0);
        }
        model.addConstraint(popConstraint.eq((double) N));

        // FLOW
        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int w = 1; w <= M; w++) {
                        double qIn = params.q(w - 1, i - 1, k - 1, m - 1);
                        double qOut = params.q(i - 1, w - 1, m - 1, k - 1);
                        cb.addTerm("UN_" + w + "_" + k, qIn);
                        cb.addTerm("UN_" + i + "_" + m, -qOut);
                    }
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        // UBAL
        for (int k = 1; k <= K; k++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int h = 1; h <= K; h++) {
                if (h != k) {
                    for (int w = 1; w <= M; w++) {
                        double qOut = params.q(M - 1, w - 1, k - 1, h - 1);
                        double qIn = params.q(M - 1, w - 1, h - 1, k - 1);
                        cb.addTerm("UN_" + M + "_" + k, qOut);
                        cb.addTerm("UN_" + M + "_" + h, -qIn);
                    }
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        // QBAL
        for (int k = 1; k <= K; k++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int h = 1; h <= K; h++) {
                if (h != k) {
                    for (int w = 1; w <= M; w++) {
                        double q = params.q(M - 1, w - 1, k - 1, h - 1);
                        cb.addTerm("QN_" + M + "_" + k, q);
                    }
                }
            }
            for (int m = 1; m <= K; m++) {
                for (int j = 1; j <= M - 1; j++) {
                    double q = params.q(M - 1, j - 1, m - 1, k - 1);
                    cb.addTerm("UN_" + M + "_" + m, q);
                }
            }
            for (int j = 1; j <= M - 1; j++) {
                double q = params.q(j - 1, M - 1, k - 1, k - 1);
                cb.addTerm("UN_" + j + "_" + k, -q);
            }
            for (int h = 1; h <= K; h++) {
                if (h != k) {
                    for (int w = 1; w <= M; w++) {
                        double q = params.q(M - 1, w - 1, h - 1, k - 1);
                        cb.addTerm("QN_" + M + "_" + h, -q);
                    }
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        // MCC
        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int w = 1; w <= M; w++) {
                        if (w != i) {
                            double q = params.q(i - 1, w - 1, k - 1, m - 1);
                            cb.addTerm("QN_" + i + "_" + k, q);
                        }
                    }
                }
            }
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double q = params.q(j - 1, i - 1, k - 1, m - 1);
                            cb.addTerm("QN_" + j + "_" + k, q);
                        }
                    }
                }
            }
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double q = params.q(j - 1, i - 1, k - 1, m - 1);
                            for (int wp = 1; wp <= M; wp++) {
                                if (wp != i && wp != j) cb.addTerm("B_" + j + "_" + k + "_" + wp, q);
                            }
                        }
                    }
                }
            }
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double q = params.q(j - 1, i - 1, k - 1, m - 1);
                            cb.addTerm("UN_" + j + "_" + k, -(N + 1) * q);
                        }
                    }
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        // MCC2
        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int w = 1; w <= M; w++) {
                        if (w != i) {
                            double q = params.q(i - 1, w - 1, k - 1, m - 1);
                            cb.addTerm("QN_" + i + "_" + k, q);
                        }
                    }
                }
            }
            for (int k = 1; k <= K; k++) {
                for (int m = 1; m <= K; m++) {
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double q = params.q(j - 1, i - 1, k - 1, m - 1);
                            cb.addTerm("B_" + j + "_" + k + "_" + i, -q);
                            cb.addTerm("UN_" + j + "_" + k, -q);
                        }
                    }
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        // QMAX
        for (int w = 1; w <= M; w++) {
            for (int k = 1; k <= K; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                        .addTerm("QN_" + w + "_" + k, 1.0)
                        .addTerm("UN_" + w + "_" + k, -(double) N);
                model.addConstraint(cb.leq(0.0));
            }
        }

        // QMIN
        for (int k = 1; k <= K; k++) {
            for (int j = 1; j <= M; j++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int w = 1; w <= M; w++) cb.addTerm("QN_" + w + "_" + k, 1.0);
                cb.addTerm("UN_" + j + "_" + k, -(double) N);
                model.addConstraint(cb.geq(0.0));
            }
        }
    }

    private static Map<String, Double> extractVariableValues(Mapqn_lpmodel model, double[] point) {
        Map<String, Double> result = new HashMap<String, Double>();
        for (Map.Entry<String, Integer> entry : model.variables.entrySet()) {
            result.put(entry.getKey(), point[entry.getValue()]);
        }
        return result;
    }

    // ---------------------------------------------------------------------
    // Top-level extension functions for Mapqn_solution (MVA version).
    // Merged from former Mapqn_bnd_lr_mva class.
    // ---------------------------------------------------------------------

    public static double getUtilizationMVA(Mapqn_solution sol, int i, int k) {
        return sol.getVariable("UN_" + i + "_" + k);
    }

    public static double getQueueLengthMVA(Mapqn_solution sol, int i, int k) {
        return sol.getVariable("QN_" + i + "_" + k);
    }
}
