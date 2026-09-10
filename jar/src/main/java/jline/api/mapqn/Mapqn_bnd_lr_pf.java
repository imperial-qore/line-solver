/**
 * @file Linear Reduction Bounds for Product Form Networks
 */
package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;

/**
 * Implementation of bnd_linearreduction_pf.mod linear program
 * This is the Product Form version without phases
 */
public final class Mapqn_bnd_lr_pf {
    private Mapqn_bnd_lr_pf() {}

    /**
     * Parameters for the Product Form linear reduction model
     */
    public static class PFParameters {
        public final int M;
        public final int N;
        public final double[] mu;
        public final double[][] r;

        public PFParameters(int M, int N, double[] mu, double[][] r) {
            this.M = M;
            this.N = N;
            this.mu = mu;
            this.r = r;
        }

        public void validate() {
            if (!(M > 0)) throw new IllegalArgumentException("M must be positive");
            if (!(N > 0)) throw new IllegalArgumentException("N must be positive");
            if (mu.length != M) throw new IllegalArgumentException("mu array size must equal M");
            if (r.length != M) throw new IllegalArgumentException("r must be MxM matrix");
            for (double[] row : r) {
                if (row.length != M) throw new IllegalArgumentException("r must be MxM matrix");
            }
            for (double v : mu) {
                if (v < 0) throw new IllegalArgumentException("Service rates must be non-negative");
            }
            for (double[] row : r) {
                for (double v : row) {
                    if (v < 0) throw new IllegalArgumentException("Routing probabilities must be non-negative");
                }
            }
        }

        public double q(int i, int j) {
            return r[i][j] * mu[i];
        }
    }

    public static Mapqn_solution solve(PFParameters params) {
        return solve(params, 1);
    }

    public static Mapqn_solution solve(PFParameters params, int objectiveQueue) {
        return solve(params, objectiveQueue, "min");
    }

    /**
     * Solve the product-form linear-reduction LP, optimizing the utilization of
     * the given queue. sense="min" yields the classical LR lower bound,
     * sense="max" the matching upper bound; both are valid bounds because the
     * LP relaxation contains the exact solution.
     *
     * @param params          product-form model parameters
     * @param objectiveQueue  1-based queue index to optimize
     * @param sense           "min" or "max"
     * @return LP solution with objective value and variable values
     */
    public static Mapqn_solution solve(PFParameters params, int objectiveQueue, String sense) {
        params.validate();
        if (!(objectiveQueue >= 1 && objectiveQueue <= params.M)) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (!"min".equals(sense) && !"max".equals(sense)) {
            throw new IllegalArgumentException("Sense must be 'min' or 'max'");
        }

        Mapqn_lpmodel model = new Mapqn_lpmodel();

        registerVariables(model, params);

        addDefinitionConstraints(model, params);
        addMeanIndicesConstraints(model, params);
        addBalanceConstraints(model, params);

        // see _kb/03-api-layer.md for rationale
        addUpperBounds(model, params);

        String objectiveVarName = "U_" + objectiveQueue;
        double[] objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName);
        LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoeffs, 0.0);

        SimplexSolver solver = new SimplexSolver();
        LinearConstraintSet constraintSet = new LinearConstraintSet(model.getConstraints());
        // see _kb/03-api-layer.md for rationale
        PointValuePair solution = solver.optimize(
                objectiveFunction,
                constraintSet,
                new NonNegativeConstraint(true),
                "max".equals(sense) ? GoalType.MAXIMIZE : GoalType.MINIMIZE);

        return new Mapqn_solution(
                solution.getValue(),
                extractVariableValues(model, solution.getPoint()));
    }

    /**
     * Impose the variable upper bounds declared by the MATLAB/Python LP models
     * (U,p1,p1c &lt;= 1; Q,C &lt;= N) as explicit LEQ constraints, since Apache's
     * SimplexSolver assumes only non-negativity.
     */
    private static void addUpperBounds(Mapqn_lpmodel model, PFParameters params) {
        int M = params.M;
        int N = params.N;
        for (int i = 1; i <= M; i++) {
            model.addConstraint(new LinearConstraint(
                    model.createObjectiveCoefficients("U_" + i), Relationship.LEQ, 1.0));
            model.addConstraint(new LinearConstraint(
                    model.createObjectiveCoefficients("Q_" + i), Relationship.LEQ, (double) N));
        }
        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                model.addConstraint(new LinearConstraint(
                        model.createObjectiveCoefficients("C_" + j + "_" + i),
                        Relationship.LEQ, (double) N));
                for (int ni = 0; ni <= N; ni++) {
                    model.addConstraint(new LinearConstraint(
                            model.createObjectiveCoefficients("p1_" + j + "_" + i + "_" + ni),
                            Relationship.LEQ, 1.0));
                    model.addConstraint(new LinearConstraint(
                            model.createObjectiveCoefficients("p1c_" + j + "_" + i + "_" + ni),
                            Relationship.LEQ, 1.0));
                }
            }
        }
    }

    private static void registerVariables(Mapqn_lpmodel model, PFParameters params) {
        int M = params.M;
        int N = params.N;

        for (int i = 1; i <= M; i++) {
            model.addVariable("U_" + i);
        }

        for (int i = 1; i <= M; i++) {
            model.addVariable("Q_" + i);
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                model.addVariable("C_" + j + "_" + i);
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                for (int ni = 0; ni <= N; ni++) {
                    model.addVariable("p1_" + j + "_" + i + "_" + ni);
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                for (int ni = 0; ni <= N; ni++) {
                    model.addVariable("p1c_" + j + "_" + i + "_" + ni);
                }
            }
        }
    }

    private static void addDefinitionConstraints(Mapqn_lpmodel model, PFParameters params) {
        int M = params.M;
        int N = params.N;

        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                    .addTerm("p1_" + j + "_" + j + "_0", 1.0);
            model.addConstraint(cb.eq(0.0));
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                if (j != i) {
                    Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                            .addTerm("p1_" + j + "_" + i + "_" + N, 1.0);
                    model.addConstraint(cb.eq(0.0));
                }
            }
        }

        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder()
                    .addTerm("C_" + j + "_" + j, 1.0)
                    .addTerm("Q_" + j, -1.0);
            model.addConstraint(cb.eq(0.0));
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int ni = 0; ni <= N; ni++) {
                    cb.addTerm("p1_" + j + "_" + i + "_" + ni, 1.0);
                    cb.addTerm("p1c_" + j + "_" + i + "_" + ni, 1.0);
                }
                model.addConstraint(cb.eq(1.0));
            }
        }
    }

    private static void addMeanIndicesConstraints(Mapqn_lpmodel model, PFParameters params) {
        int M = params.M;
        int N = params.N;

        for (int i = 1; i <= M; i++) {
            for (int t = 1; t <= M; t++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                cb.addTerm("U_" + i, 1.0);
                for (int nt = 0; nt <= N; nt++) {
                    cb.addTerm("p1_" + i + "_" + t + "_" + nt, -1.0);
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            cb.addTerm("Q_" + i, 1.0);
            for (int ni = 0; ni <= N; ni++) {
                cb.addTerm("p1_" + i + "_" + i + "_" + ni, -(double) ni);
            }
            model.addConstraint(cb.eq(0.0));
        }

        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= M; i++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                cb.addTerm("C_" + j + "_" + i, 1.0);
                for (int ni = 0; ni <= N; ni++) {
                    cb.addTerm("p1_" + j + "_" + i + "_" + ni, -(double) ni);
                }
                model.addConstraint(cb.eq(0.0));
            }
        }
    }

    private static void addBalanceConstraints(Mapqn_lpmodel model, PFParameters params) {
        int M = params.M;
        int N = params.N;

        for (int j = 1; j <= M; j++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int i = 1; i <= M; i++) {
                cb.addTerm("C_" + j + "_" + i, 1.0);
            }
            cb.addTerm("U_" + j, -(double) N);
            model.addConstraint(cb.eq(0.0));
        }

        Mapqn_lpmodel.LinearConstraintBuilder popConstraint = model.constraintBuilder();
        for (int i = 1; i <= M; i++) {
            popConstraint.addTerm("Q_" + i, 1.0);
        }
        model.addConstraint(popConstraint.eq((double) N));

        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int j = 1; j <= M; j++) {
                if (j != i) {
                    double qji = params.q(j - 1, i - 1);
                    double qij = params.q(i - 1, j - 1);
                    cb.addTerm("p1_" + j + "_" + i + "_0", qji);
                    cb.addTerm("p1_" + i + "_" + i + "_1", -qij);
                }
            }
            model.addConstraint(cb.eq(0.0));
        }

        for (int i = 1; i <= M; i++) {
            for (int ni = 1; ni <= N - 1; ni++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int j = 1; j <= M; j++) {
                    if (j != i) {
                        double qji = params.q(j - 1, i - 1);
                        double qij = params.q(i - 1, j - 1);
                        cb.addTerm("p1_" + j + "_" + i + "_" + ni, qji);
                        cb.addTerm("p1_" + i + "_" + i + "_" + (ni + 1), -qij);
                    }
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        for (int i = 1; i <= M; i++) {
            for (int j = 1; j <= M; j++) {
                Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
                for (int ni = 1; ni <= N; ni++) {
                    cb.addTerm("p1_" + j + "_" + i + "_" + ni, 1.0);
                }
                for (int nj = 1; nj <= N; nj++) {
                    cb.addTerm("p1_" + i + "_" + j + "_" + nj, -1.0);
                }
                model.addConstraint(cb.eq(0.0));
            }
        }

        for (int i = 1; i <= M; i++) {
            Mapqn_lpmodel.LinearConstraintBuilder cb = model.constraintBuilder();
            for (int j = 1; j <= M; j++) {
                if (j != i) {
                    double qij = params.q(i - 1, j - 1);
                    double qji = params.q(j - 1, i - 1);

                    cb.addTerm("U_" + i, qij);

                    for (int nj = 1; nj <= N; nj++) {
                        cb.addTerm("p1_" + i + "_" + j + "_" + nj, -qji);
                    }
                    cb.addTerm("p1_" + j + "_" + i + "_0", -qji);
                }
            }
            model.addConstraint(cb.eq(0.0));
        }
    }

    private static Map<String, Double> extractVariableValues(Mapqn_lpmodel model, double[] point) {
        Map<String, Double> result = new HashMap<String, Double>();
        for (Map.Entry<String, Integer> entry : model.variables.entrySet()) {
            result.put(entry.getKey(), Double.valueOf(point[entry.getValue().intValue()]));
        }
        return result;
    }
}
