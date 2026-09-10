package jline.opt.pareto;

import jline.lang.Network;
import jline.opt.OptimizationProblem;
import jline.opt.objectives.Constraint;
import jline.opt.results.OptimizationResult;
import jline.opt.solver.BisectionSolver;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.DecisionVariable;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Epsilon-constraint sweep for bi-objective tradeoff analysis. Solves the
 * problem once per epsilon, each time adding {@code constraintFactory(epsilon)},
 * and filters to the non-dominated cost frontier. Mirrors native-Python
 * {@code line_solver.opt.pareto.ParetoSweep}.
 */
public class ParetoSweep {

    /** Maps an epsilon value to a constraint. */
    public interface ConstraintFactory {
        Constraint make(double epsilon);
    }

    private final OptimizationProblem problem;
    private final ConstraintFactory constraintFactory;
    private final double[] epsilons;
    private final String solver;   // "de" or "bisection"
    private List<ParetoPoint> points = new ArrayList<ParetoPoint>();

    public ParetoSweep(OptimizationProblem problem, ConstraintFactory constraintFactory,
                       double[] epsilons, String solver) {
        this.problem = problem;
        this.constraintFactory = constraintFactory;
        this.epsilons = epsilons.clone();
        this.solver = solver;
    }

    public ParetoSweep(OptimizationProblem problem, ConstraintFactory constraintFactory,
                       double[] epsilons) {
        this(problem, constraintFactory, epsilons, "de");
    }

    public List<ParetoPoint> getPoints() {
        return new ArrayList<ParetoPoint>(points);
    }

    private OptimizationProblem cloneProblem(double epsilon) {
        OptimizationProblem clone = new OptimizationProblem(problem.getModel());
        for (DecisionVariable var : problem.getVariables()) {
            clone.addVariable(var);
        }
        clone.setObjective(problem.getObjective());
        for (Constraint c : problem.getConstraints()) {
            clone.addConstraint(c);
        }
        clone.setFixedVariables(problem.getFixedVariables());
        for (Pair<Network, Double> s : problem.getScenarios()) {
            clone.addScenario(s.getLeft(), s.getRight());
        }
        clone.addConstraint(constraintFactory.make(epsilon));
        return clone;
    }

    public List<ParetoPoint> solve(LineOptSolverOptions options) {
        points = new ArrayList<ParetoPoint>();
        for (double epsilon : epsilons) {
            OptimizationProblem clone = cloneProblem(epsilon);
            OptimizationResult result;
            if ("bisection".equals(solver)) {
                result = new BisectionSolver(clone).solve();
            } else {
                result = clone.solve(options != null ? options : new LineOptSolverOptions());
            }
            points.add(new ParetoPoint(epsilon, result.objectiveValue, result.feasible, result));
        }
        return new ArrayList<ParetoPoint>(points);
    }

    public List<ParetoPoint> solve() {
        return solve(new LineOptSolverOptions());
    }

    /** Non-dominated frontier (both objective and epsilon minimized), by epsilon. */
    public List<ParetoPoint> getFrontier() {
        List<ParetoPoint> feasible = new ArrayList<ParetoPoint>();
        for (ParetoPoint p : points) {
            if (p.feasible) {
                feasible.add(p);
            }
        }
        List<ParetoPoint> frontier = new ArrayList<ParetoPoint>();
        for (ParetoPoint p : feasible) {
            boolean dominated = false;
            for (ParetoPoint q : feasible) {
                if (q.objectiveValue <= p.objectiveValue && q.epsilon <= p.epsilon
                        && (q.objectiveValue < p.objectiveValue || q.epsilon < p.epsilon)) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) {
                frontier.add(p);
            }
        }
        Collections.sort(frontier, new Comparator<ParetoPoint>() {
            public int compare(ParetoPoint a, ParetoPoint b) {
                return Double.compare(a.epsilon, b.epsilon);
            }
        });
        return frontier;
    }
}
