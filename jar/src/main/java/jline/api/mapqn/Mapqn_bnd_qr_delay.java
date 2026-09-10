/**
 * @file Quadratic Reduction Bounds for Delay Systems
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;

/**
 * Implementation of bnd_quadraticreduction_delay.mod linear program.
 * This model includes delay constraints and load-dependent service rates.
 */
public final class Mapqn_bnd_qr_delay {

    private Mapqn_bnd_qr_delay() {}

    /** Parameters for the quadratic reduction delay model. */
    public static final class QuadraticDelayParameters {
        public final int M;
        public final int N;
        public final int[] K;
        public final double Z;
        public final double D1;
        public final double[][][] mu;
        public final double[][][] v;
        public final double[][] alpha;
        public final double[][] r;

        public QuadraticDelayParameters(int M, int N, int[] K, double Z, double D1,
                                        double[][][] mu, double[][][] v,
                                        double[][] alpha, double[][] r) {
            this.M = M;
            this.N = N;
            this.K = K;
            this.Z = Z;
            this.D1 = D1;
            this.mu = mu;
            this.v = v;
            this.alpha = alpha;
            this.r = r;
        }

        public void validate() {
            if (M <= 0) throw new IllegalArgumentException("M must be positive");
            if (N <= 0) throw new IllegalArgumentException("N must be positive");
            if (Z <= 0) throw new IllegalArgumentException("Z must be positive");
            if (D1 <= 0) throw new IllegalArgumentException("D1 must be positive");
            if (K.length != M) throw new IllegalArgumentException("K array size must equal M");
            if (mu.length != M) throw new IllegalArgumentException("mu array size must equal M");
            if (v.length != M) throw new IllegalArgumentException("v array size must equal M");
            if (alpha.length != M) throw new IllegalArgumentException("alpha array size must equal M");
            if (r.length != M) throw new IllegalArgumentException("r must be MxM matrix");
            for (int i = 0; i < M; i++) {
                if (r[i].length != M) throw new IllegalArgumentException("r must be MxM matrix");
            }
            for (int i = 0; i < M; i++) {
                if (K[i] <= 0) throw new IllegalArgumentException("K[" + i + "] must be positive");
                if (mu[i].length != K[i]) throw new IllegalArgumentException("mu[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
                for (int j = 0; j < K[i]; j++) {
                    if (mu[i][j].length != K[i]) throw new IllegalArgumentException("mu[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
                }
                if (v[i].length != K[i]) throw new IllegalArgumentException("v[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
                for (int j = 0; j < K[i]; j++) {
                    if (v[i][j].length != K[i]) throw new IllegalArgumentException("v[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
                }
                if (alpha[i].length != N) throw new IllegalArgumentException("alpha[" + i + "] must have size N");
                for (double rij : r[i]) {
                    if (rij < 0) throw new IllegalArgumentException("Routing probabilities must be non-negative");
                }
            }
        }

        /** Compute q parameter with load dependence. */
        public double q(int i, int j, int k, int h, int n) {
            // Matches MATLAB q_func: an empty queue has no completions, and
            // population levels beyond the alpha table use unit scaling
            if (n == 0) {
                return 0.0;
            }
            double alphaVal = (n <= alpha[i].length) ? alpha[i][n - 1] : 1.0;
            if (j != i) {
                return r[i][j] * mu[i][k][h] * alphaVal;
            } else {
                return (v[i][k][h] + r[i][i] * mu[i][k][h]) * alphaVal;
            }
        }
    }

    public static Mapqn_solution solve(
            QuadraticDelayParameters params,
            int objectiveQueue,
            int objectivePhase,
            int objectiveN) {
        return solve(params, objectiveQueue, objectivePhase, objectiveN, "max");
    }

    /**
     * As above, optimizing in the requested direction. The LP is a relaxation
     * containing the exact solution, so "max" is a valid upper bound on the
     * marginal probability and "min" the matching lower bound.
     *
     * @param sense "min" or "max"
     */
    public static Mapqn_solution solve(
            QuadraticDelayParameters params,
            int objectiveQueue,
            int objectivePhase,
            int objectiveN,
            String sense) {
        params.validate();
        if (!("min".equals(sense) || "max".equals(sense))) {
            throw new IllegalArgumentException("Sense must be 'min' or 'max'");
        }
        if (objectiveQueue < 1 || objectiveQueue > params.M) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (objectivePhase < 1 || objectivePhase > params.K[objectiveQueue - 1]) {
            throw new IllegalArgumentException("Objective phase must be in range 1..K[" + (objectiveQueue - 1) + "]");
        }
        if (objectiveN < 0 || objectiveN > params.N) {
            throw new IllegalArgumentException("Objective N must be in range 0..N");
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();

        registerVariables(model, params);
        addVariableBounds(model);
        addDefinitionConstraints(model, params);
        addLittlesLawConstraints(model, params);
        addBalanceConstraints(model, params);
        addCorrelationConstraints(model, params);
        addBoundConstraints(model, params);

        String objectiveVarName = "p2_" + objectiveQueue + "_" + objectiveN + "_" + objectivePhase
                + "_" + objectiveQueue + "_" + objectiveN + "_" + objectivePhase;
        double[] objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName);
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        PointValuePair solution = solver.optimize(
                objectiveFunction,
                constraintSet,
                "min".equals(sense) ? GoalType.MINIMIZE : GoalType.MAXIMIZE);

        return new Mapqn_solution(solution.getValue(), extractVariableValues(model, solution.getPoint()));
    }

    // Box bounds lb=0, ub=1 on every variable (mapqn_bnd_qr_delay.m:112-113).
    private static void addVariableBounds(Mapqn_lpmodel model) {
        for (String name : new java.util.ArrayList<String>(model.variables.keySet())) {
            model.addConstraint(model.constraintBuilder().addTerm(name, 1.0).geq(0.0));
            model.addConstraint(model.constraintBuilder().addTerm(name, 1.0).leq(1.0));
        }
    }

    private static void registerVariables(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;
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

    private static void addDefinitionConstraints(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
            for (int nj = 0; nj <= N; nj++) {
                for (int k = 1; k <= K[j - 1]; k++) {
                    constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + nj + "_" + k, 1.0);
                }
            }
            model.addConstraint(constraint.eq(1.0));
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int h = 1; h <= K[j - 1]; h++) {
                        if (h != k) {
                            Mapqn_lpmodel.LinearConstraintBuilder c = model.constraintBuilder()
                                    .addTerm("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + nj + "_" + h, 1.0);
                            model.addConstraint(c.eq(0.0));
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int ni = 0; ni <= N; ni++) {
                        if (nj != ni) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                Mapqn_lpmodel.LinearConstraintBuilder c = model.constraintBuilder()
                                        .addTerm("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + ni + "_" + h, 1.0);
                                model.addConstraint(c.eq(0.0));
                            }
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int i = 1; i <= M; i++) {
                        if (i != j) {
                            for (int ni = 0; ni <= N; ni++) {
                                if (nj + ni > N) {
                                    for (int h = 1; h <= K[i - 1]; h++) {
                                        Mapqn_lpmodel.LinearConstraintBuilder c = model.constraintBuilder()
                                                .addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, 1.0);
                                        model.addConstraint(c.eq(0.0));
                                    }
                                }
                            }
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
                                Mapqn_lpmodel.LinearConstraintBuilder c = model.constraintBuilder()
                                        .addTerm("p2_" + i + "_" + ni + "_" + h + "_" + j + "_" + nj + "_" + k, 1.0)
                                        .addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                                model.addConstraint(c.eq(0.0));
                            }
                        }
                    }
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int nj = 0; nj <= N; nj++) {
                    for (int i = 1; i <= M; i++) {
                        if (i != j) {
                            Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                            constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + nj + "_" + k, 1.0);
                            for (int ni = 0; ni <= N - nj; ni++) {
                                for (int h = 1; h <= K[i - 1]; h++) {
                                    constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, -1.0);
                                }
                            }
                            model.addConstraint(constraint.eq(0.0));
                        }
                    }
                }
            }
        }
    }

    private static void addLittlesLawConstraints(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;
        double Z = params.Z;
        double D1 = params.D1;

        Mapqn_lpmodel.LinearConstraintBuilder pc2Constraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            for (int j = 1; j <= M; j++) {
                for (int ni = 1; ni <= N; ni++) {
                    for (int nj = 1; nj <= N; nj++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            for (int k = 1; k <= K[j - 1]; k++) {
                                pc2Constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h,
                                        (double) (nj * ni));
                            }
                        }
                    }
                }
            }
        }
        model.addConstraint(pc2Constraint.eq((double) (N * N)));

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int i = 1; i <= M; i++) {
                    for (int nj = 1; nj <= N; nj++) {
                        for (int ni = 1; ni <= N; ni++) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, (double) ni);
                            }
                        }
                    }
                }
                for (int nj = 1; nj <= N; nj++) {
                    constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + nj + "_" + k, -(double) N);
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int i = 1; i <= M; i++) {
                    for (int ni = 1; ni <= N; ni++) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            constraint.addTerm("p2_" + j + "_0_" + k + "_" + i + "_" + ni + "_" + h, (double) ni);
                        }
                    }
                }
                constraint.addTerm("p2_" + j + "_0_" + k + "_" + j + "_0_" + k, -(double) N);
                model.addConstraint(constraint.eq(0.0));
            }
        }

        Mapqn_lpmodel.LinearConstraintBuilder xzConstraint = model.constraintBuilder();
        for (int ni = 1; ni <= N; ni++) {
            for (int k = 1; k <= K[M - 1]; k++) {
                xzConstraint.addTerm("p2_" + M + "_" + ni + "_" + k + "_" + M + "_" + ni + "_" + k, (double) ni);
            }
        }
        for (int k = 1; k <= K[0]; k++) {
            for (int nj = 1; nj <= N; nj++) {
                xzConstraint.addTerm("p2_1_" + nj + "_" + k + "_1_" + nj + "_" + k, -Z / D1);
            }
        }
        model.addConstraint(xzConstraint.eq(0.0));
    }

    private static void addBalanceConstraints(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    for (int h = 1; h <= K[i - 1]; h++) {
                        if (!(h == k && i == j)) {
                            for (int ni = 1; ni <= N; ni++) {
                                double qOut = params.q(i - 1, j - 1, k - 1, h - 1, ni);
                                double qIn = params.q(i - 1, j - 1, h - 1, k - 1, ni);
                                constraint.addTerm("p2_" + i + "_" + ni + "_" + k + "_" + i + "_" + ni + "_" + k, qOut);
                                constraint.addTerm("p2_" + i + "_" + ni + "_" + h + "_" + i + "_" + ni + "_" + h, -qIn);
                            }
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }

        for (int i = 1; i <= M; i++) {
            for (int ni = 1; ni < N; ni++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[j - 1]; k++) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                for (int u = 1; u <= K[i - 1]; u++) {
                                    for (int nj = 1; nj <= N - ni; nj++) {
                                        double q = params.q(j - 1, i - 1, k - 1, h - 1, nj);
                                        constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + u, q);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            for (int h = 1; h <= K[i - 1]; h++) {
                                double q = params.q(i - 1, j - 1, k - 1, h - 1, ni + 1);
                                constraint.addTerm("p2_" + i + "_" + (ni + 1) + "_" + k + "_" + i + "_" + (ni + 1) + "_" + k, -q);
                            }
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }

        for (int i = 1; i <= M; i++) {
            for (int u = 1; u <= K[i - 1]; u++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[j - 1]; k++) {
                            for (int h = 1; h <= K[j - 1]; h++) {
                                for (int nj = 1; nj <= N; nj++) {
                                    double q = params.q(j - 1, i - 1, k - 1, h - 1, nj);
                                    constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_0_" + u, q);
                                }
                            }
                        }
                    }
                }
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k = 1; k <= K[i - 1]; k++) {
                            double q = params.q(i - 1, j - 1, k - 1, u - 1, 1);
                            constraint.addTerm("p2_" + i + "_1_" + k + "_" + i + "_1_" + k, -q);
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }

        // QBAL: Queue balance (mapqn_bnd_qr_delay.m:416-495)
        // LHS1 + LHS2 = RHS1 + RHS2, ported as LHS1 + LHS2 - RHS1 - RHS2 = 0.
        for (int i = 1; i <= M; i++) {
            for (int k = 1; k <= K[i - 1]; k++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                // LHS1: sum{h!=k, j, ni>=1} q(i,j,k,h,ni)*ni*p2(i,ni,k,i,ni,k)
                for (int h = 1; h <= K[i - 1]; h++) {
                    if (h != k) {
                        for (int j = 1; j <= M; j++) {
                            for (int ni = 1; ni <= N; ni++) {
                                double q = params.q(i - 1, j - 1, k - 1, h - 1, ni);
                                constraint.addTerm("p2_" + i + "_" + ni + "_" + k + "_" + i + "_" + ni + "_" + k, q * ni);
                            }
                        }
                    }
                }
                // LHS2: sum{j!=i, h, ni>=1} q(i,j,h,k,ni)*p2(i,ni,h,i,ni,h)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int h = 1; h <= K[i - 1]; h++) {
                            for (int ni = 1; ni <= N; ni++) {
                                double q = params.q(i - 1, j - 1, h - 1, k - 1, ni);
                                constraint.addTerm("p2_" + i + "_" + ni + "_" + h + "_" + i + "_" + ni + "_" + h, q);
                            }
                        }
                    }
                }
                // -RHS1a: sum{j!=i, u in K(j), w in K(j), nj>=1} q(j,i,u,w,nj)*p2(j,nj,u,i,0,k)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int u = 1; u <= K[j - 1]; u++) {
                            for (int w = 1; w <= K[j - 1]; w++) {
                                for (int nj = 1; nj <= N; nj++) {
                                    double q = params.q(j - 1, i - 1, u - 1, w - 1, nj);
                                    constraint.addTerm("p2_" + j + "_" + nj + "_" + u + "_" + i + "_0_" + k, -q);
                                }
                            }
                        }
                    }
                }
                // -RHS1b: sum{j!=i, u, w, nj>=1, ni>=1} q(j,i,u,w,nj)*p2(i,ni,k,j,nj,u)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int u = 1; u <= K[j - 1]; u++) {
                            for (int w = 1; w <= K[j - 1]; w++) {
                                for (int nj = 1; nj <= N; nj++) {
                                    double q = params.q(j - 1, i - 1, u - 1, w - 1, nj);
                                    for (int ni = 1; ni <= N; ni++) {
                                        constraint.addTerm("p2_" + i + "_" + ni + "_" + k + "_" + j + "_" + nj + "_" + u, -q);
                                    }
                                }
                            }
                        }
                    }
                }
                // -RHS2: sum{h!=k, j, ni>=1} q(i,j,h,k,ni)*ni*p2(i,ni,h,i,ni,h)
                for (int h = 1; h <= K[i - 1]; h++) {
                    if (h != k) {
                        for (int j = 1; j <= M; j++) {
                            for (int ni = 1; ni <= N; ni++) {
                                double q = params.q(i - 1, j - 1, h - 1, k - 1, ni);
                                constraint.addTerm("p2_" + i + "_" + ni + "_" + h + "_" + i + "_" + ni + "_" + h, -q * ni);
                            }
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }
    }

    private static void addCorrelationConstraints(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        // COR1a: correlation constraint (mapqn_bnd_qr_delay.m:497-604)
        for (int i = 1; i <= M; i++) {
            for (int kstar = 1; kstar <= K[i - 1]; kstar++) {
                for (int nic = 0; nic <= N - 2; nic++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                    // A: sum{j!=i,kj,hj,u!=kstar,nj=1..N-nic} q(j,i,kj,hj,nj)*p2(j,nj,kj,i,nic,u)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int kj = 1; kj <= K[j - 1]; kj++) {
                                for (int hj = 1; hj <= K[j - 1]; hj++) {
                                    for (int u = 1; u <= K[i - 1]; u++) {
                                        if (u != kstar) {
                                            for (int nj = 1; nj <= N - nic; nj++) {
                                                double q = params.q(j - 1, i - 1, kj - 1, hj - 1, nj);
                                                constraint.addTerm("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + nic + "_" + u, q);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // B: sum{j!=i,kj,hj,nj=1..N-nic} q(j,i,kj,hj,nj)*p2(j,nj,kj,i,nic+1,kstar)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int kj = 1; kj <= K[j - 1]; kj++) {
                                for (int hj = 1; hj <= K[j - 1]; hj++) {
                                    for (int nj = 1; nj <= N - nic; nj++) {
                                        double q = params.q(j - 1, i - 1, kj - 1, hj - 1, nj);
                                        constraint.addTerm("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + (nic + 1) + "_" + kstar, q);
                                    }
                                }
                            }
                        }
                    }
                    // C: sum{k!=kstar} q(i,i,kstar,k,nic+1)*p2(i,nic+1,kstar,i,nic+1,kstar)
                    for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                        if (k2 != kstar) {
                            double q = params.q(i - 1, i - 1, kstar - 1, k2 - 1, nic + 1);
                            constraint.addTerm("p2_" + i + "_" + (nic + 1) + "_" + kstar + "_" + i + "_" + (nic + 1) + "_" + kstar, q);
                        }
                    }
                    // -D: -sum{j!=i,k!=kstar} q(i,j,k,k,nic+1)*p2(i,nic+1,k,i,nic+1,k)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                                if (k2 != kstar) {
                                    double q = params.q(i - 1, j - 1, k2 - 1, k2 - 1, nic + 1);
                                    constraint.addTerm("p2_" + i + "_" + (nic + 1) + "_" + k2 + "_" + i + "_" + (nic + 1) + "_" + k2, -q);
                                }
                            }
                        }
                    }
                    // -E: -sum{j!=i,k!=kstar,h!=k} q(i,j,k,h,nic+1)*p2(i,nic+1,k,i,nic+1,k)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                                if (k2 != kstar) {
                                    for (int h2 = 1; h2 <= K[i - 1]; h2++) {
                                        if (h2 != k2) {
                                            double q = params.q(i - 1, j - 1, k2 - 1, h2 - 1, nic + 1);
                                            constraint.addTerm("p2_" + i + "_" + (nic + 1) + "_" + k2 + "_" + i + "_" + (nic + 1) + "_" + k2, -q);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // -F: -sum{j!=i,k!=kstar} q(i,j,k,kstar,nic+2)*p2(i,nic+2,k,i,nic+2,k)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                                if (k2 != kstar) {
                                    double q = params.q(i - 1, j - 1, k2 - 1, kstar - 1, nic + 2);
                                    constraint.addTerm("p2_" + i + "_" + (nic + 2) + "_" + k2 + "_" + i + "_" + (nic + 2) + "_" + k2, -q);
                                }
                            }
                        }
                    }
                    // -G: -sum{j!=i} q(i,j,kstar,kstar,nic+2)*p2(i,nic+2,kstar,i,nic+2,kstar)
                    for (int j = 1; j <= M; j++) {
                        if (j != i) {
                            double q = params.q(i - 1, j - 1, kstar - 1, kstar - 1, nic + 2);
                            constraint.addTerm("p2_" + i + "_" + (nic + 2) + "_" + kstar + "_" + i + "_" + (nic + 2) + "_" + kstar, -q);
                        }
                    }
                    // -H: -sum{k!=kstar} q(i,i,k,kstar,nic+1)*p2(i,nic+1,k,i,nic+1,k)
                    for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                        if (k2 != kstar) {
                            double q = params.q(i - 1, i - 1, k2 - 1, kstar - 1, nic + 1);
                            constraint.addTerm("p2_" + i + "_" + (nic + 1) + "_" + k2 + "_" + i + "_" + (nic + 1) + "_" + k2, -q);
                        }
                    }
                    model.addConstraint(constraint.eq(0.0));
                }
            }
        }

        // COR1b: correlation constraint boundary (mapqn_bnd_qr_delay.m:606-675)
        for (int i = 1; i <= M; i++) {
            for (int kstar = 1; kstar <= K[i - 1]; kstar++) {
                Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                // A': sum{j!=i,kj,hj,u!=kstar} q(j,i,kj,hj,1)*p2(j,1,kj,i,N-1,u)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int kj = 1; kj <= K[j - 1]; kj++) {
                            for (int hj = 1; hj <= K[j - 1]; hj++) {
                                for (int u = 1; u <= K[i - 1]; u++) {
                                    if (u != kstar) {
                                        double q = params.q(j - 1, i - 1, kj - 1, hj - 1, 1);
                                        constraint.addTerm("p2_" + j + "_1_" + kj + "_" + i + "_" + (N - 1) + "_" + u, q);
                                    }
                                }
                            }
                        }
                    }
                }
                // C': sum{k!=kstar} q(i,i,kstar,k,N)*p2(i,N,kstar,i,N,kstar)
                for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                    if (k2 != kstar) {
                        double q = params.q(i - 1, i - 1, kstar - 1, k2 - 1, N);
                        constraint.addTerm("p2_" + i + "_" + N + "_" + kstar + "_" + i + "_" + N + "_" + kstar, q);
                    }
                }
                // -D': -sum{j!=i,k!=kstar} q(i,j,k,k,N)*p2(i,N,k,i,N,k)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                            if (k2 != kstar) {
                                double q = params.q(i - 1, j - 1, k2 - 1, k2 - 1, N);
                                constraint.addTerm("p2_" + i + "_" + N + "_" + k2 + "_" + i + "_" + N + "_" + k2, -q);
                            }
                        }
                    }
                }
                // -E': -sum{j!=i,k!=kstar,h!=k} q(i,j,k,h,N)*p2(i,N,k,i,N,k)
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                            if (k2 != kstar) {
                                for (int h2 = 1; h2 <= K[i - 1]; h2++) {
                                    if (h2 != k2) {
                                        double q = params.q(i - 1, j - 1, k2 - 1, h2 - 1, N);
                                        constraint.addTerm("p2_" + i + "_" + N + "_" + k2 + "_" + i + "_" + N + "_" + k2, -q);
                                    }
                                }
                            }
                        }
                    }
                }
                // -H': -sum{k!=kstar} q(i,i,k,kstar,N)*p2(i,N,k,i,N,k)
                for (int k2 = 1; k2 <= K[i - 1]; k2++) {
                    if (k2 != kstar) {
                        double q = params.q(i - 1, i - 1, k2 - 1, kstar - 1, N);
                        constraint.addTerm("p2_" + i + "_" + N + "_" + k2 + "_" + i + "_" + N + "_" + k2, -q);
                    }
                }
                model.addConstraint(constraint.eq(0.0));
            }
        }
    }

    private static void addBoundConstraints(Mapqn_lpmodel model, QuadraticDelayParameters params) {
        int M = params.M;
        int N = params.N;
        int[] K = params.K;

        for (int j = 1; j <= M; j++) {
            for (int k = 1; k <= K[j - 1]; k++) {
                for (int i = 1; i <= M; i++) {
                    Mapqn_lpmodel.LinearConstraintBuilder constraint = model.constraintBuilder();
                    for (int t = 1; t <= M; t++) {
                        for (int h = 1; h <= K[t - 1]; h++) {
                            for (int nj = 0; nj <= N; nj++) {
                                for (int nt = 0; nt <= N; nt++) {
                                    constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + t + "_" + nt + "_" + h, (double) nt);
                                }
                            }
                        }
                    }
                    for (int h = 1; h <= K[i - 1]; h++) {
                        for (int nj = 0; nj <= N; nj++) {
                            for (int ni = 0; ni <= N; ni++) {
                                constraint.addTerm("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h, -(double) N);
                            }
                        }
                    }
                    model.addConstraint(constraint.geq(0.0));
                }
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

    /** Convenience accessor (extension-function equivalent). */
    public static double getP2(Mapqn_solution sol, int j, int nj, int k, int i, int ni, int h) {
        return sol.getVariable("p2_" + j + "_" + nj + "_" + k + "_" + i + "_" + ni + "_" + h);
    }

    public static double getMarginalP2(Mapqn_solution sol, int j, int nj, int k) {
        return sol.getVariable("p2_" + j + "_" + nj + "_" + k + "_" + j + "_" + nj + "_" + k);
    }
}
