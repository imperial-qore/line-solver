/**
 * @file MAPQN Linear Programming Model
 * 
 * Base class for representing MAP queueing network linear programming models.
 * Provides the foundation for LP-based optimization methods in MAPQN analysis,
 * including constraint formulation and objective function definition.
 * 
 * @since LINE 3.0
 */
package jline.api.mapqn

import org.apache.commons.math3.optim.linear.*
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import java.util.*

/**
 * Base class for representing MAPQN Linear Programming models
 */
class Mapqn_lpmodel {
    private val constraints = mutableListOf<LinearConstraint>()
    internal val variables = mutableMapOf<String, Int>()
    private var variableCounter = 0
    
    /**
     * Register a variable and return its index
     */
    fun addVariable(name: String): Int {
        if (variables.containsKey(name)) {
            return variables[name]!!
        }
        val index = variableCounter++
        variables[name] = index
        return index
    }
    
    /**
     * Get variable index by name
     */
    fun getVariableIndex(name: String): Int {
        return variables[name] ?: throw IllegalArgumentException("Variable $name not found")
    }
    
    /**
     * Get total number of variables
     */
    fun getNumVariables(): Int = variableCounter
    
    /**
     * Add a constraint to the model
     */
    fun addConstraint(constraint: LinearConstraint) {
        constraints.add(constraint)
    }
    
    /**
     * Get all constraints
     */
    fun getConstraints(): List<LinearConstraint> = constraints.toList()
    
    /**
     * Create a linear constraint builder
     */
    fun constraintBuilder(numVars: Int = getNumVariables()): LinearConstraintBuilder {
        return LinearConstraintBuilder(numVars)
    }
    
    /**
     * Helper class for building linear constraints
     */
    inner class LinearConstraintBuilder(private val numVars: Int) {
        private val coefficients = DoubleArray(numVars)
        
        fun addTerm(varName: String, coefficient: Double): LinearConstraintBuilder {
            val index = getVariableIndex(varName)
            coefficients[index] += coefficient
            return this
        }
        
        fun addTerm(varIndex: Int, coefficient: Double): LinearConstraintBuilder {
            coefficients[varIndex] += coefficient
            return this
        }
        
        fun eq(rhs: Double): LinearConstraint {
            return LinearConstraint(coefficients, Relationship.EQ, rhs)
        }
        
        fun leq(rhs: Double): LinearConstraint {
            return LinearConstraint(coefficients, Relationship.LEQ, rhs)
        }
        
        fun geq(rhs: Double): LinearConstraint {
            return LinearConstraint(coefficients, Relationship.GEQ, rhs)
        }
    }
    
    /**
     * Create objective function coefficients array
     */
    fun createObjectiveCoefficients(objectiveVar: String): DoubleArray {
        val coeffs = DoubleArray(getNumVariables())
        val index = getVariableIndex(objectiveVar)
        coeffs[index] = 1.0
        return coeffs
    }
}