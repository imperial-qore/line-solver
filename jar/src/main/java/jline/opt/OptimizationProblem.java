package jline.opt;

import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.Objective;
import jline.opt.results.OptimizationResult;
import jline.opt.solver.LineOptSolver;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.DecisionVariable;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Declarative specification of a queueing-network optimization problem: a LINE
 * model, decision variables, an objective, constraints, optional fixed
 * variables (for decomposition) and optional workload scenarios (for robust
 * optimization). Mirrors native-Python
 * {@code line_solver.opt.problem.OptimizationProblem}.
 *
 * <p>The model may be a flat {@link Network} or a {@link LayeredNetwork} (LQN);
 * an LQN takes the SolverLN evaluation path and accepts LQN decision variables
 * only (see {@link LqnAdapter}). Workload scenarios are supported for flat
 * models only.</p>
 */
public class OptimizationProblem {

    private final Network model;
    /** Non-null iff the model is a LayeredNetwork (LQN). */
    private final LayeredNetwork lqnModel;
    private final boolean layered;
    private final List<DecisionVariable> variables = new ArrayList<DecisionVariable>();
    private Objective objective;
    private final List<Constraint> constraints = new ArrayList<Constraint>();
    private List<Pair<DecisionVariable, Object>> fixedVariables =
            new ArrayList<Pair<DecisionVariable, Object>>();
    private final List<Pair<Network, Double>> scenarios =
            new ArrayList<Pair<Network, Double>>();

    // Decision-variable types that operate on an LQN vs a flat Network.
    private static final Set<String> LQN_VAR_TYPES = new HashSet<String>();
    private static final Set<String> FLAT_VAR_TYPES = new HashSet<String>();

    static {
        LQN_VAR_TYPES.add("host_demand");
        LQN_VAR_TYPES.add("think_time");
        LQN_VAR_TYPES.add("task_multiplicity");
        LQN_VAR_TYPES.add("task_replication");
        LQN_VAR_TYPES.add("processor_multiplicity");
        FLAT_VAR_TYPES.add("server_allocation");
        FLAT_VAR_TYPES.add("station_replicas");
        FLAT_VAR_TYPES.add("service_rate");
        FLAT_VAR_TYPES.add("job_population");
        FLAT_VAR_TYPES.add("class_priority");
        FLAT_VAR_TYPES.add("routing");
        FLAT_VAR_TYPES.add("class_mapping");
    }

    public OptimizationProblem(Network model) {
        this.model = model;
        this.lqnModel = null;
        this.layered = false;
    }

    public OptimizationProblem(LayeredNetwork lqnModel) {
        this.model = null;
        this.lqnModel = lqnModel;
        this.layered = true;
    }

    /** True if the model is a LayeredNetwork (LQN). */
    public boolean isLayered() {
        return layered;
    }

    /** The flat network model, or null when the problem is layered. */
    public Network getModel() {
        return model;
    }

    /** The LQN model, or null when the problem is flat. */
    public LayeredNetwork getLqnModel() {
        return lqnModel;
    }

    public List<DecisionVariable> getVariables() {
        return new ArrayList<DecisionVariable>(variables);
    }

    public Objective getObjective() {
        return objective;
    }

    public List<Constraint> getConstraints() {
        return new ArrayList<Constraint>(constraints);
    }

    public OptimizationProblem addVariable(DecisionVariable variable) {
        variables.add(variable);
        return this;
    }

    public OptimizationProblem setObjective(Objective objective) {
        this.objective = objective;
        return this;
    }

    public OptimizationProblem addConstraint(Constraint constraint) {
        constraints.add(constraint);
        return this;
    }

    public OptimizationProblem setFixedVariables(List<Pair<DecisionVariable, Object>> pairs) {
        this.fixedVariables = new ArrayList<Pair<DecisionVariable, Object>>(pairs);
        return this;
    }

    public List<Pair<DecisionVariable, Object>> getFixedVariables() {
        return new ArrayList<Pair<DecisionVariable, Object>>(fixedVariables);
    }

    public OptimizationProblem addScenario(Network scenarioModel, double weight) {
        scenarios.add(new Pair<Network, Double>(scenarioModel, weight));
        return this;
    }

    public List<Pair<Network, Double>> getScenarios() {
        return new ArrayList<Pair<Network, Double>>(scenarios);
    }

    public List<String> validate() {
        List<String> errors = new ArrayList<String>();
        if (model == null && lqnModel == null) {
            errors.add("Model is not set");
        }
        if (variables.isEmpty()) {
            errors.add("No decision variables defined");
        }
        if (objective == null) {
            errors.add("Objective function is not set");
        }
        // Model/variable-kind consistency: LQN models take LQN variables and
        // flat models take flat variables; mixing them silently produces no-ops
        // (a variable whose element is never found in the copy).
        for (DecisionVariable var : variables) {
            String vtype = var.getVariableType();
            if (layered && FLAT_VAR_TYPES.contains(vtype)) {
                errors.add("Variable '" + var.getName() + "' (" + vtype + ") is a flat-network "
                        + "variable but the model is a LayeredNetwork");
            } else if (!layered && LQN_VAR_TYPES.contains(vtype)) {
                errors.add("Variable '" + var.getName() + "' (" + vtype + ") is a LayeredNetwork "
                        + "variable but the model is a flat Network");
            }
        }
        return errors;
    }

    public boolean isValid() {
        return validate().isEmpty();
    }

    public OptimizationResult solve() {
        return solve(new LineOptSolverOptions());
    }

    public OptimizationResult solve(LineOptSolverOptions options) {
        List<String> errors = validate();
        if (!errors.isEmpty()) {
            throw new IllegalStateException("Invalid problem: " + String.join(", ", errors));
        }
        return new LineOptSolver(this, options).solve();
    }

    /** Create a decomposition workflow for this problem. */
    public jline.opt.decomposition.DecompositionWorkflow decompose() {
        return new jline.opt.decomposition.DecompositionWorkflow(this);
    }
}
