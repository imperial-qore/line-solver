package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;

/**
 * Implementation of bnd_mvaversion.mod linear program: the MAP-AMVA
 * optimization program of G. Casale, E. Smirni, "MAP-AMVA: Approximate Mean
 * Value Analysis of Bursty Systems", IEEE/IFIP DSN 2009, pp. 409-418.
 *
 * <p>The constraints assembled below are the paper's: the population constraint
 * (1), the utilization bound (2), the MAP phase balance (3), the flow balance
 * (4), the generalized horizontal cut (12), the vertical-cut MVA relation (13)
 * in the linearized form (18)-(19) whose {@code B(j,k,i)} variables are the
 * E_i^{j,k} of Theorem 4, and the two auxiliary families QN &lt;= N*UN and
 * sum_w QN &gt;= N*UN.
 */
public final class Mapqn_bnd_lr_mva {
    private Mapqn_bnd_lr_mva() {}

    public static Mapqn_solution solve(MVAVersionParameters params, int objectiveQueue, int objectiveLevel) {
        return solve(params, objectiveQueue, objectiveLevel, "max");
    }

    /**
     * As above, optimizing in the requested direction. The LP is a relaxation
     * containing the exact solution, so "max" is a valid upper bound on the
     * utilization at that population level and "min" the matching lower bound.
     *
     * @param sense "min" or "max"
     */
    public static Mapqn_solution solve(MVAVersionParameters params, int objectiveQueue, int objectiveLevel,
                                       String sense) {
        return solve(params, objectiveQueue, objectiveLevel, sense, "UN");
    }

    /**
     * As above, over a chosen variable family and optionally over the SUM of
     * the levels.
     *
     * <p>{@code objectiveLevel == 0} optimizes the aggregate sum_k X(queue, k),
     * which is the quantity the paper's bounds are stated on: U_i(N) = sum_k
     * U_i^k(N) is the utilization of station i, while U_i^k alone is its
     * utilization while the MAP sits in phase k. Optimizing the K terms
     * separately and adding them is also a bound but a strictly looser one,
     * since the phases cannot all peak at once.
     *
     * @param objectiveVar "UN" or "QN", the variable family optimized over
     */
    public static Mapqn_solution solve(MVAVersionParameters params, int objectiveQueue, int objectiveLevel,
                                       String sense, String objectiveVar) {
        params.validate();
        if (!("min".equals(sense) || "max".equals(sense))) {
            throw new IllegalArgumentException("Sense must be 'min' or 'max'");
        }
        if (objectiveQueue < 1 || objectiveQueue > params.M) {
            throw new IllegalArgumentException("Objective queue must be in range 1.." + params.M);
        }
        if (objectiveLevel < 0 || objectiveLevel > params.K) {
            throw new IllegalArgumentException(
                    "Objective level must be in range 0.." + params.K + " (0 aggregates over levels)");
        }
        if (!("UN".equals(objectiveVar) || "QN".equals(objectiveVar))) {
            throw new IllegalArgumentException("objectiveVar must be 'UN' or 'QN'");
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();
        registerVariables(model, params);
        addConstraints(model, params);

        // Objective over one level, or over their sum when objectiveLevel is 0.
        double[] objectiveCoeffs = new double[model.getNumVariables()];
        int firstLevel = objectiveLevel == 0 ? 1 : objectiveLevel;
        int lastLevel = objectiveLevel == 0 ? params.K : objectiveLevel;
        for (int k = firstLevel; k <= lastLevel; k++) {
            objectiveCoeffs[model.getVariableIndex(objectiveVar + "_" + objectiveQueue + "_" + k)] = 1.0;
        }
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        PointValuePair solution = solver.optimize(objectiveFunction, constraintSet,
                "min".equals(sense) ? GoalType.MINIMIZE : GoalType.MAXIMIZE);

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
