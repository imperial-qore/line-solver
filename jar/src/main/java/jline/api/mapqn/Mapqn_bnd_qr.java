/**
 * General Quadratic Reduction Bounds
 *
 * Implements general quadratic reduction methods for computing performance
 * bounds in MAP queueing networks. Provides the core quadratic approximation
 * algorithms used across various MAPQN bound computation techniques.
 */
package jline.api.mapqn;

import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

import java.util.HashMap;
import java.util.Map;

/**
 * Implementation of bnd_quadraticreduction.mod linear program.
 * This model uses quadratic (p2) variables for joint probabilities.
 */
public final class Mapqn_bnd_qr {

    private Mapqn_bnd_qr() {
    }

    /**
     * Create and solve the linear program for bnd_quadraticreduction model.
     */
    public static Mapqn_solution solve(LinearReductionParameters params, int objectiveQueue, int objectivePhase) {
        params.validate();
        if (!(objectiveQueue >= 1 && objectiveQueue <= params.M)) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (!(objectivePhase >= 1 && objectivePhase <= params.K[objectiveQueue - 1])) {
            throw new IllegalArgumentException("Objective phase must be in range 1..K[" + (objectiveQueue - 1) + "]");
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();

        registerVariables(model, params);

        addDefinitionConstraints(model, params);
        addMeanIndicesConstraints(model, params);
        addBalanceConstraints(model, params);
        addBoundConstraints(model, params);
        addQuadraticConstraints(model, params);

        String objectiveVarName = "U_" + objectiveQueue + "_" + objectivePhase;
        double[] objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName);
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        PointValuePair solution = solver.optimize(objectiveFunction, constraintSet, GoalType.MAXIMIZE);

        return new Mapqn_solution(solution.getValue(), extractVariableValues(model, solution.getPoint()));
    }

    private static void registerVariables(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("U_" + i + "_" + k);
            }
        }
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("IT_" + i + "_" + k);
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        model.addVariable("UP_" + j + "_" + k + "_" + i + "_" + h);
                    }
                }
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        model.addVariable("QP_" + j + "_" + k + "_" + i + "_" + h);
                    }
                }
            }
        }
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("Q_" + i + "_" + k);
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    model.addVariable("C_" + j + "_" + k + "_" + i);
                }
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    model.addVariable("I_" + j + "_" + k + "_" + i);
                }
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            model.addVariable("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h);
                        }
                    }
                }
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            model.addVariable("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h);
                        }
                    }
                }
            }
        }
        for (int j = 1; j <= M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int k = 1; k <= K[j - 1]; k++) {
                    for (int i = 1; i <= M; i++) {
                        for (int ni = 0; ni <= N; ni++) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                model.addVariable("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h);
                            }
                        }
                    }
                }
            }
        }
    }

    private static void addQuadraticConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        Mapqn_lpmodel.LinearConstraintBuilder pcl2Constraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            for (int j = 1; j <= M; j++) {
                for (int ni = 1; ni <= N; ni++) {
                    for (int nj = 1; nj <= N; nj++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            for (int k = 1; k <= K[j - 1]; k++) {
                                pcl2Constraint.addTerm("p2_" + i + "_" + ni + "_" + h + "_" + j + "_" + nj + "_" + k,
                                        (double) (nj * ni));
                            }
                        }
                    }
                }
            }
        }
        model.addConstraint(pcl2Constraint.eq((double) (N * N)));

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                            constraint.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                            for (int nj = 1; nj <= N; nj++) {
                                constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                            }
                            model.addConstraint(constraint.eq(0.0));
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            LinearConstraint constraint = model.constraintBuilder()
                                    .addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0)
                                    .addTerm("p2_" + j + "_0_" + k + "_" + i + "_" + ni + "_" + h, -1.0)
                                    .eq(0.0);
                            model.addConstraint(constraint);
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int k = 1; k <= K[j - 1]; k++) {
                    for (int i = 1; i <= M; i++) {
                        for (int ni = 0; ni <= N; ni++) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                LinearConstraint constraint = model.constraintBuilder()
                                        .addTerm("p2_" + i + "_" + ni + "_" + h + "_" + j + "_" + nj + "_" + k, 1.0)
                                        .addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0)
                                        .eq(0.0);
                                model.addConstraint(constraint);
                            }
                        }
                    }
                }
            }
        }
    }

    private static void addDefinitionConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                LinearConstraint constraint = model.constraintBuilder()
                        .addTerm("p1_" + j + "_" + k + "_" + j + "_0_" + k, 1.0)
                        .eq(0.0);
                model.addConstraint(constraint);
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int h = 1; h <= K[j - 1]; h++) {
                        if (h != k) {
                            LinearConstraint constraint = model.constraintBuilder()
                                    .addTerm("p1_" + j + "_" + k + "_" + j + "_" + nj + "_" + h, 1.0)
                                    .eq(0.0);
                            model.addConstraint(constraint);
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        if (j != i) {
                            LinearConstraint constraint = model.constraintBuilder()
                                    .addTerm("p1_" + j + "_" + k + "_" + i + "_" + N + "_" + h, 1.0)
                                    .eq(0.0);
                            model.addConstraint(constraint);
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 1; nj <= N; nj++) {
                    for (int h = 1; h <= K[j - 1]; h++) {
                        LinearConstraint constraint = model.constraintBuilder()
                                .addTerm("p1c_" + j + "_" + k + "_" + j + "_" + nj + "_" + h, 1.0)
                                .eq(0.0);
                        model.addConstraint(constraint);
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                LinearConstraint constraint = model.constraintBuilder()
                        .addTerm("C_" + j + "_" + k + "_" + j, 1.0)
                        .addTerm("Q_" + j + "_" + k, -1.0)
                        .eq(0.0);
                model.addConstraint(constraint);
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int k = 1; k <= K[j - 1]; k++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        for (int ni = 0; ni <= N; ni++) {
                            constraint.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                            constraint.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                        }
                    }
                }
                model.addConstraint(constraint.eq(1.0));
            }
        }
    }

    private static void addMeanIndicesConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                for (int t = 1; t <= M; t++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                    constraint.addTerm("U_" + i + "_" + k, 1.0);
                    for (int nt = 0; nt <= N; nt++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            constraint.addTerm("p1_" + i + "_" + k + "_" + t + "_" + nt + "_" + h, -1.0);
                        }
                    }
                    model.addConstraint(constraint.eq(0.0));
                }
            }
        }

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                for (int t = 1; t <= M; t++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                    constraint.addTerm("IT_" + i + "_" + k, 1.0);
                    for (int nt = 0; nt <= N; nt++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            constraint.addTerm("p1c_" + i + "_" + k + "_" + t + "_" + nt + "_" + h, -1.0);
                        }
                    }
                    model.addConstraint(constraint.eq(0.0));
                }
            }
        }

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                constraint.addTerm("Q_" + i + "_" + k, 1.0);
                for (int ni = 0; ni <= N; ni++) {
                    constraint.addTerm("p1_" + i + "_" + k + "_" + i + "_" + ni + "_" + k, -((double) ni));
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }
    }

    private static void addBalanceConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        double qCoeff1 = params.q(i - 1, j - 1, k - 1, h - 1);
                        double qCoeff2 = params.q(i - 1, j - 1, h - 1, k - 1);
                        constraint.addTerm("U_" + i + "_" + k, qCoeff1);
                        constraint.addTerm("U_" + i + "_" + h, -qCoeff2);
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }

        Mapqn_lpmodel.LinearConstraintBuilder popConstraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                popConstraint.addTerm("Q_" + i + "_" + k, 1.0);
            }
        }
        model.addConstraint(popConstraint.eq((double) N));

        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
            for (int k = 1; k <= K[j - 1]; k++) {
                constraint.addTerm("U_" + j + "_" + k, 1.0);
                constraint.addTerm("IT_" + j + "_" + k, 1.0);
            }
            model.addConstraint(constraint.eq(1.0));
        }
    }

    private static void addBoundConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
            for (int k = 1; k <= K[i - 1]; k++) {
                constraint.addTerm("U_" + i + "_" + k, 1.0);
            }
            model.addConstraint(constraint.leq(1.0));
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                LinearConstraint constraint = model.constraintBuilder()
                        .addTerm("Q_" + j + "_" + k, 1.0)
                        .addTerm("U_" + j + "_" + k, -((double) N))
                        .leq(0.0);
                model.addConstraint(constraint);
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                    constraint.addTerm("C_" + j + "_" + k + "_" + i, 1.0);
                    for (int h = 1; h <= K[i - 1]; h++) {
                        constraint.addTerm("Q_" + i + "_" + h, -1.0);
                    }
                    model.addConstraint(constraint.leq(0.0));
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    LinearConstraint constraint = model.constraintBuilder()
                            .addTerm("C_" + j + "_" + k + "_" + i, 1.0)
                            .addTerm("U_" + j + "_" + k, -((double) N))
                            .leq(0.0);
                    model.addConstraint(constraint);
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();

                    for (int t = 1; t <= M; t++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            for (int nt = 0; nt <= N; nt++) {
                                constraint.addTerm("p1_" + j + "_" + k + "_" + t + "_" + nt + "_" + h, (double) nt);
                                constraint.addTerm("p1c_" + j + "_" + k + "_" + t + "_" + nt + "_" + h, (double) nt);
                            }
                        }
                    }

                    for (int h = 1; h <= K[i - 1]; h++) {
                        for (int ni = 0; ni <= N; ni++) {
                            constraint.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) N));
                            constraint.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) N));
                        }
                    }

                    model.addConstraint(constraint.geq(0.0));
                }
            }
        }
    }

    private static Map<String, Double> extractVariableValues(Mapqn_lpmodel model, double[] point) {
        Map<String, Double> result = new HashMap<String, Double>();
        for (Map.Entry<String, Integer> e : model.variables.entrySet()) {
            result.put(e.getKey(), point[e.getValue()]);
        }
        return result;
    }
}
