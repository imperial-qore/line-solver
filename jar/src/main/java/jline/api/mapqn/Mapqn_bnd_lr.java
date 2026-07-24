/**
 * @file General Linear Reduction Bounds
 *
 * Implements general linear reduction methods for computing performance
 * bounds in MAP queueing networks. Provides the core linear approximation
 * algorithms used across various MAPQN bound computation techniques.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

/**
 * Implementation of bnd_linearreduction_new.mod linear program.
 */
public final class Mapqn_bnd_lr {

    private Mapqn_bnd_lr() {}

    /**
     * Create and solve the linear program for bnd_linearreduction_new model.
     *
     * @param params Model parameters
     * @param objectiveQueue Queue index to maximize utilization (1-based)
     * @param objectivePhase Phase index to maximize utilization (1-based)
     * @return Solution containing the optimal value and variable values
     */
    public static Mapqn_solution solve(LinearReductionParameters params, int objectiveQueue, int objectivePhase) {
        params.validate();
        if (objectiveQueue < 1 || objectiveQueue > params.M) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (objectivePhase < 1 || objectivePhase > params.K[objectiveQueue - 1]) {
            throw new IllegalArgumentException("Objective phase must be in range 1..K[" + (objectiveQueue - 1) + "]");
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();

        // Register all variables
        registerVariables(model, params);

        // Add all constraints
        addDefinitionConstraints(model, params);
        addMeanIndicesConstraints(model, params);
        addBalanceConstraints(model, params);
        addBoundConstraints(model, params);

        // Create objective function - maximize U[objectiveQueue][objectivePhase]
        String objectiveVarName = "U_" + objectiveQueue + "_" + objectivePhase;
        double[] objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName);
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        // Solve the LP
        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        PointValuePair solution = solver.optimize(objectiveFunction, constraintSet, GoalType.MAXIMIZE);

        return new Mapqn_solution(solution.getValue(), extractVariableValues(model, solution.getPoint()));
    }

    private static void registerVariables(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // U variables
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("U_" + i + "_" + k);
            }
        }

        // IT variables
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("IT_" + i + "_" + k);
            }
        }

        // UP variables
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        model.addVariable("UP_" + j + "_" + k + "_" + i + "_" + h);
                    }
                }
            }
        }

        // QP variables
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        model.addVariable("QP_" + j + "_" + k + "_" + i + "_" + h);
                    }
                }
            }
        }

        // Q variables
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                model.addVariable("Q_" + i + "_" + k);
            }
        }

        // C variables
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    model.addVariable("C_" + j + "_" + k + "_" + i);
                }
            }
        }

        // I variables
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    model.addVariable("I_" + j + "_" + k + "_" + i);
                }
            }
        }

        // p1 variables
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

        // p1c variables
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
    }

    private static void addDefinitionConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // ZER1: p1[j,k,j,0,k]=0
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                model.addConstraint(model.constraintBuilder()
                        .addTerm("p1_" + j + "_" + k + "_" + j + "_0_" + k, 1.0)
                        .eq(0.0));
            }
        }

        // ZER2: p1[j,k,j,nj,h]=0 for h<>k
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int h = 1; h <= K[j - 1]; h++) {
                        if (h != k) {
                            model.addConstraint(model.constraintBuilder()
                                    .addTerm("p1_" + j + "_" + k + "_" + j + "_" + nj + "_" + h, 1.0)
                                    .eq(0.0));
                        }
                    }
                }
            }
        }

        // ZER3: p1[j,k,i,N,h]=0 for j<>i
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        if (j != i) {
                            model.addConstraint(model.constraintBuilder()
                                    .addTerm("p1_" + j + "_" + k + "_" + i + "_" + N + "_" + h, 1.0)
                                    .eq(0.0));
                        }
                    }
                }
            }
        }

        // ZER4: p1c[j,k,j,nj,h]=0 for nj>=1
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 1; nj <= N; nj++) {
                    for (int h = 1; h <= K[j - 1]; h++) {
                        model.addConstraint(model.constraintBuilder()
                                .addTerm("p1c_" + j + "_" + k + "_" + j + "_" + nj + "_" + h, 1.0)
                                .eq(0.0));
                    }
                }
            }
        }

        // CEQU: C[j,k,j] = Q[j,k]
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                model.addConstraint(model.constraintBuilder()
                        .addTerm("C_" + j + "_" + k + "_" + j, 1.0)
                        .addTerm("Q_" + j + "_" + k, -1.0)
                        .eq(0.0));
            }
        }

        // ONE1: sum over all p1 and p1c = 1
        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int k = 1; k <= K[j - 1]; k++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        for (int ni = 0; ni <= N; ni++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                            cb.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                        }
                    }
                }
                model.addConstraint(cb.eq(1.0));
            }
        }
    }

    private static void addMeanIndicesConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // UTLB: U[i,k]=sum{nt,h} p1[i,k,t,nt,h]  (for each t)
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                for (int t = 1; t <= M; t++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    cb.addTerm("U_" + i + "_" + k, 1.0);
                    for (int nt = 0; nt <= N; nt++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            cb.addTerm("p1_" + i + "_" + k + "_" + t + "_" + nt + "_" + h, -1.0);
                        }
                    }
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }

        // UTLC: IT[i,k]=sum{nt,h} p1c[i,k,t,nt,h]  (for each t)
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                for (int t = 1; t <= M; t++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    cb.addTerm("IT_" + i + "_" + k, 1.0);
                    for (int nt = 0; nt <= N; nt++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            cb.addTerm("p1c_" + i + "_" + k + "_" + t + "_" + nt + "_" + h, -1.0);
                        }
                    }
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }

        // QLEN: Q[i,k]=sum{ni} ni*p1[i,k,i,ni,k]
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                cb.addTerm("Q_" + i + "_" + k, 1.0);
                for (int ni = 0; ni <= N; ni++) {
                    cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + ni + "_" + k, -((double) ni));
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // UPH1
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                        cb.addTerm("UP_" + j + "_" + k + "_" + i + "_" + h, 1.0);
                        for (int ni = 1; ni <= N; ni++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                            cb.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                        }
                        model.addConstraint(cb.eq(0.0));
                    }
                }
            }
        }

        // UPH2
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                        cb.addTerm("UP_" + j + "_" + k + "_" + i + "_" + h, 1.0);
                        for (int ni = 1; ni <= N; ni++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                        }
                        cb.addTerm("p1_" + i + "_" + h + "_" + j + "_0_" + k, -1.0);
                        model.addConstraint(cb.eq(0.0));
                    }
                }
            }
        }

        // QPH1
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                        cb.addTerm("QP_" + j + "_" + k + "_" + i + "_" + h, 1.0);
                        for (int ni = 1; ni <= N; ni++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) ni));
                            cb.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) ni));
                        }
                        model.addConstraint(cb.eq(0.0));
                    }
                }
            }
        }

        // CLEN
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    cb.addTerm("C_" + j + "_" + k + "_" + i, 1.0);
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) ni));
                        }
                    }
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }

        // ILEN
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    cb.addTerm("I_" + j + "_" + k + "_" + i, 1.0);
                    for (int ni = 0; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            cb.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) ni));
                        }
                    }
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }
    }

    private static void addBalanceConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // SRVB
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        if (h != k) {
                            double qCoeff = params.q(i - 1, j - 1, k - 1, h - 1);
                            cb.addTerm("U_" + i + "_" + k, qCoeff);
                            cb.addTerm("U_" + i + "_" + h, -params.q(i - 1, j - 1, h - 1, k - 1));
                        }
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // MPCB
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int i = 1; i <= M; i++) {
                    cb.addTerm("C_" + j + "_" + k + "_" + i, 1.0);
                }
                cb.addTerm("U_" + j + "_" + k, -((double) N));
                model.addConstraint(cb.eq(0.0));
            }
        }

        // MPCI
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int i = 1; i <= M; i++) {
                    cb.addTerm("I_" + j + "_" + k + "_" + i, 1.0);
                }
                cb.addTerm("IT_" + j + "_" + k, -((double) N));
                model.addConstraint(cb.eq(0.0));
            }
        }

        // ONE
        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K[j - 1]; k++) {
                cb.addTerm("U_" + j + "_" + k, 1.0);
                cb.addTerm("IT_" + j + "_" + k, 1.0);
            }
            model.addConstraint(cb.eq(1.0));
        }

        // POPC
        Mapqn_lpmodel.LinearConstraintBuilder popConstraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                popConstraint.addTerm("Q_" + i + "_" + k, 1.0);
            }
        }
        model.addConstraint(popConstraint.eq((double) N));

        // GFFL0
        for (int i = 1; i <= M; i++) {
            for (int u = 1; u <= K[i - 1]; u++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[j - 1]; k++) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                double qVal = params.q(j - 1, i - 1, k - 1, h - 1);
                                cb.addTerm("p1_" + j + "_" + k + "_" + i + "_0_" + u, qVal);
                            }
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            double qVal = params.q(i - 1, j - 1, k - 1, u - 1);
                            cb.addTerm("p1_" + i + "_" + k + "_" + i + "_1_" + k, -qVal);
                        }
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // GFFL
        for (int i = 1; i <= M; i++) {
            for (int ni = 1; ni <= N - 1; ni++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[j - 1]; k++) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                for (int u = 1; u <= K[i - 1]; u++) {
                                    double qVal = params.q(j - 1, i - 1, k - 1, h - 1);
                                    cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + u, qVal);
                                }
                            }
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                double qVal = params.q(i - 1, j - 1, k - 1, h - 1);
                                cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + (ni + 1) + "_" + k, -qVal);
                            }
                        }
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // UJNT
        for (int i = 1; i <= M; i++) {
            for (int j = 1; j <= M; j++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int ni = 1; ni <= N; ni++) {
                    for (int k = 1; k <= K[j - 1]; k++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                        }
                    }
                }
                for (int nj = 1; nj <= N; nj++) {
                    for (int k = 1; k <= K[j - 1]; k++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            cb.addTerm("p1_" + i + "_" + h + "_" + j + "_" + nj + "_" + k, -1.0);
                        }
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // QBAL
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int h = 1; h <= K[i - 1]; h++) {
                    if (h != k) {
                        for (int j = 1; j <= M; j++) {
                            double qVal = params.q(i - 1, j - 1, k - 1, h - 1);
                            cb.addTerm("Q_" + i + "_" + k, qVal);
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            double qVal = params.q(i - 1, j - 1, h - 1, k - 1);
                            cb.addTerm("U_" + i + "_" + h, qVal);
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int u = 1; u <= K[j - 1]; u++) {
                            for (int w = 1; w <= K[j - 1]; w++) {
                                double qVal = params.q(j - 1, i - 1, u - 1, w - 1);
                                for (int nj = 1; nj <= N; nj++) {
                                    cb.addTerm("p1_" + i + "_" + k + "_" + j + "_" + nj + "_" + u, -qVal);
                                }
                                cb.addTerm("p1_" + j + "_" + u + "_" + i + "_0_" + k, -qVal);
                            }
                        }
                    }
                }
                for (int h = 1; h <= K[i - 1]; h++) {
                    if (h != k) {
                        for (int j = 1; j <= M; j++) {
                            double qVal = params.q(i - 1, j - 1, h - 1, k - 1);
                            cb.addTerm("Q_" + i + "_" + h, -qVal);
                        }
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        // MBH
        for (int i = 1; i <= M; i++) {
            for (int ordphase = 1; ordphase <= K[i - 1]; ordphase++) {
                for (int ni = 0; ni <= N - 2; ni++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k = 1; k <= K[j - 1]; k++) {
                                for (int h = 1; h <= K[j - 1]; h++) {
                                    for (int u = 1; u <= K[i - 1]; u++) {
                                        if (u != ordphase) {
                                            double qVal = params.q(j - 1, i - 1, k - 1, h - 1);
                                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + u, qVal);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k = 1; k <= K[j - 1]; k++) {
                                for (int h = 1; h <= K[j - 1]; h++) {
                                    double qVal = params.q(j - 1, i - 1, k - 1, h - 1);
                                    cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + (ni + 1) + "_" + ordphase, qVal);
                                }
                            }
                        }
                    }
                    for (int k = 1; k <= K[i - 1]; k++) {
                        if (k != ordphase) {
                            double qVal = params.q(i - 1, i - 1, ordphase - 1, k - 1);
                            cb.addTerm("p1_" + i + "_" + ordphase + "_" + i + "_" + (ni + 1) + "_" + ordphase, qVal);
                        }
                    }
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k = 1; k <= K[i - 1]; k++) {
                                if (k != ordphase) {
                                    double qVal = params.q(i - 1, j - 1, k - 1, k - 1);
                                    cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + (ni + 1) + "_" + k, -qVal);
                                }
                            }
                        }
                    }
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k = 1; k <= K[i - 1]; k++) {
                                if (k != ordphase) {
                                    for (int h = 1; h <= K[i - 1]; h++) {
                                        if (h != k) {
                                            double qVal = params.q(i - 1, j - 1, k - 1, h - 1);
                                            cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + (ni + 1) + "_" + k, -qVal);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k = 1; k <= K[i - 1]; k++) {
                                if (k != ordphase) {
                                    double qVal = params.q(i - 1, j - 1, k - 1, ordphase - 1);
                                    cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + (ni + 2) + "_" + k, -qVal);
                                }
                            }
                        }
                    }
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double qVal = params.q(i - 1, j - 1, ordphase - 1, ordphase - 1);
                            cb.addTerm("p1_" + i + "_" + ordphase + "_" + i + "_" + (ni + 2) + "_" + ordphase, -qVal);
                        }
                    }
                    for (int k = 1; k <= K[i - 1]; k++) {
                        if (k != ordphase) {
                            double qVal = params.q(i - 1, i - 1, k - 1, ordphase - 1);
                            cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + (ni + 1) + "_" + k, -qVal);
                        }
                    }
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }

        // MBHN
        for (int i = 1; i <= M; i++) {
            for (int ordphase = 1; ordphase <= K[i - 1]; ordphase++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[j - 1]; k++) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                for (int u = 1; u <= K[i - 1]; u++) {
                                    if (u != ordphase) {
                                        double qVal = params.q(j - 1, i - 1, k - 1, h - 1);
                                        cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + (N - 1) + "_" + u, qVal);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int k = 1; k <= K[i - 1]; k++) {
                    if (k != ordphase) {
                        double qVal = params.q(i - 1, i - 1, ordphase - 1, k - 1);
                        cb.addTerm("p1_" + i + "_" + ordphase + "_" + i + "_" + N + "_" + ordphase, qVal);
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            if (k != ordphase) {
                                double qVal = params.q(i - 1, j - 1, k - 1, k - 1);
                                cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + N + "_" + k, -qVal);
                            }
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            if (k != ordphase) {
                                for (int h = 1; h <= K[i - 1]; h++) {
                                    if (h != k) {
                                        double qVal = params.q(i - 1, j - 1, k - 1, h - 1);
                                        cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + N + "_" + k, -qVal);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int k = 1; k <= K[i - 1]; k++) {
                    if (k != ordphase) {
                        double qVal = params.q(i - 1, i - 1, k - 1, ordphase - 1);
                        cb.addTerm("p1_" + i + "_" + k + "_" + i + "_" + N + "_" + k, -qVal);
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }
    }

    private static void addBoundConstraints(Mapqn_lpmodel model, LinearReductionParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // UUB1
        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int k = 1; k <= K[i - 1]; k++) {
                cb.addTerm("U_" + i + "_" + k, 1.0);
            }
            model.addConstraint(cb.leq(1.0));
        }

        // QUB1
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                        .addTerm("Q_" + j + "_" + k, 1.0)
                        .addTerm("U_" + j + "_" + k, -((double) N));
                model.addConstraint(cb.leq(0.0));
            }
        }

        // CUB1
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    cb.addTerm("C_" + j + "_" + k + "_" + i, 1.0);
                    for (int h = 1; h <= K[i - 1]; h++) {
                        cb.addTerm("Q_" + i + "_" + h, -1.0);
                    }
                    model.addConstraint(cb.leq(0.0));
                }
            }
        }

        // CUB2
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                            .addTerm("C_" + j + "_" + k + "_" + i, 1.0)
                            .addTerm("U_" + j + "_" + k, -((double) N));
                    model.addConstraint(cb.leq(0.0));
                }
            }
        }

        // QMIN
        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                    for (int t = 1; t <= M; t++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            for (int nt = 0; nt <= N; nt++) {
                                cb.addTerm("p1_" + j + "_" + k + "_" + t + "_" + nt + "_" + h, (double) nt);
                                cb.addTerm("p1c_" + j + "_" + k + "_" + t + "_" + nt + "_" + h, (double) nt);
                            }
                        }
                    }
                    for (int h = 1; h <= K[i - 1]; h++) {
                        for (int ni = 0; ni <= N; ni++) {
                            cb.addTerm("p1_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) N));
                            cb.addTerm("p1c_" + j + "_" + k + "_" + i + "_" + ni + "_" + h, -((double) N));
                        }
                    }
                    model.addConstraint(cb.geq(0.0));
                }
            }
        }

        // TEST
        for (int j = 1; j <= M; j++) {
            for (int nj = 1; nj <= N - 1; nj++) {
                for (int k = 1; k <= K[j - 1]; k++) {
                    for (int i = 1; i <= M; i++) {
                        if (i != j) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                                        .addTerm("p1_" + j + "_" + k + "_" + j + "_" + nj + "_" + k, 1.0)
                                        .addTerm("p1_" + i + "_" + h + "_" + j + "_" + nj + "_" + k, -1.0);
                                model.addConstraint(cb.geq(0.0));
                            }
                        }
                    }
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
