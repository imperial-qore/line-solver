/**
 * @file MAPQN Linear Programming Model
 *
 * Base class for representing MAP queueing network linear programming models.
 * Provides the foundation for LP-based optimization methods in MAPQN analysis,
 * including constraint formulation and objective function definition.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 * Base class for representing MAPQN Linear Programming models
 */
public class Mapqn_lpmodel {
    private final List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
    final Map<String, Integer> variables = new HashMap<String, Integer>();
    private int variableCounter = 0;

    /**
     * Register a variable and return its index
     */
    public int addVariable(String name) {
        if (variables.containsKey(name)) {
            return variables.get(name);
        }
        int index = variableCounter++;
        variables.put(name, index);
        return index;
    }

    /**
     * Get variable index by name
     */
    public int getVariableIndex(String name) {
        Integer idx = variables.get(name);
        if (idx == null) {
            throw new IllegalArgumentException("Variable " + name + " not found");
        }
        return idx;
    }

    /**
     * Get total number of variables
     */
    public int getNumVariables() {
        return variableCounter;
    }

    /**
     * Add a constraint to the model
     */
    public void addConstraint(LinearConstraint constraint) {
        constraints.add(constraint);
    }

    /**
     * Get all constraints
     */
    public List<LinearConstraint> getConstraints() {
        return new ArrayList<LinearConstraint>(constraints);
    }

    /**
     * Create a linear constraint builder
     */
    public LinearConstraintBuilder constraintBuilder() {
        return new LinearConstraintBuilder(getNumVariables());
    }

    public LinearConstraintBuilder constraintBuilder(int numVars) {
        return new LinearConstraintBuilder(numVars);
    }

    /**
     * Helper class for building linear constraints
     */
    public class LinearConstraintBuilder {
        private final int numVars;
        private final double[] coefficients;

        public LinearConstraintBuilder(int numVars) {
            this.numVars = numVars;
            this.coefficients = new double[numVars];
        }

        public LinearConstraintBuilder addTerm(String varName, double coefficient) {
            int index = getVariableIndex(varName);
            coefficients[index] += coefficient;
            return this;
        }

        public LinearConstraintBuilder addTerm(int varIndex, double coefficient) {
            coefficients[varIndex] += coefficient;
            return this;
        }

        public LinearConstraint eq(double rhs) {
            return new LinearConstraint(coefficients, Relationship.EQ, rhs);
        }

        public LinearConstraint leq(double rhs) {
            return new LinearConstraint(coefficients, Relationship.LEQ, rhs);
        }

        public LinearConstraint geq(double rhs) {
            return new LinearConstraint(coefficients, Relationship.GEQ, rhs);
        }
    }

    /**
     * Create objective function coefficients array
     */
    public double[] createObjectiveCoefficients(String objectiveVar) {
        double[] coeffs = new double[getNumVariables()];
        int index = getVariableIndex(objectiveVar);
        coeffs[index] = 1.0;
        return coeffs;
    }
}
