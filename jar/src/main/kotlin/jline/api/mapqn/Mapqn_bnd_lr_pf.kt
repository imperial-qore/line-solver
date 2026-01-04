/**
 * @file Linear Reduction Bounds for Product Form Networks
 * 
 * Implements linear reduction bounds specialized for product-form MAP
 * queueing networks. Leverages product-form properties to provide
 * efficient bounds computation for networks with separable solutions.
 * 
 * @since LINE 3.0
 */
package jline.api.mapqn

import org.apache.commons.math3.optim.linear.LinearConstraint
import org.apache.commons.math3.optim.linear.LinearConstraintSet
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction
import org.apache.commons.math3.optim.linear.SimplexSolver
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType

/**
 * Implementation of bnd_linearreduction_pf.mod linear program
 * This is the Product Form version without phases
 */
object Mapqn_bnd_lr_pf {
    
    /**
     * Parameters for the Product Form linear reduction model
     */
    class PFParameters(
        val M: Int,  // number of queues
        val N: Int,  // population
        val mu: DoubleArray,  // service rates [i]
        val r: Array<DoubleArray>  // routing probabilities [i][j]
    ) {
        fun validate() {
            require(M > 0) { "M must be positive" }
            require(N > 0) { "N must be positive" }
            require(mu.size == M) { "mu array size must equal M" }
            require(r.size == M && r.all { it.size == M }) { "r must be MxM matrix" }
            require(mu.all { it >= 0 }) { "Service rates must be non-negative" }
            require(r.all { row -> row.all { it >= 0 } }) { "Routing probabilities must be non-negative" }
        }
        
        /**
         * Compute q parameter as defined in AMPL model
         */
        fun q(i: Int, j: Int): Double {
            return r[i][j] * mu[i]
        }
    }
    
    /**
     * Create and solve the linear program for bnd_linearreduction_pf model
     * 
     * @param params Model parameters
     * @param objectiveQueue Queue index to minimize utilization (1-based), default is 1
     * @return Solution containing the optimal value and variable values
     */
    fun solve(params: PFParameters, objectiveQueue: Int = 1): Mapqn_solution {
        params.validate()
        require(objectiveQueue in 1..params.M) { "Objective queue must be in range 1..M" }
        
        val model = Mapqn_lpmodel()
        
        // Register all variables
        registerVariables(model, params)
        
        // Add all constraints
        addDefinitionConstraints(model, params)
        addMeanIndicesConstraints(model, params)
        addBalanceConstraints(model, params)
        
        // Create objective function - minimize U[objectiveQueue]
        val objectiveVarName = "U_$objectiveQueue"
        val objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName)
        val objectiveFunction = LinearObjectiveFunction(objectiveCoeffs, 0.0)
        
        // Solve the LP
        val solver = SimplexSolver()
        val constraintSet = LinearConstraintSet(model.getConstraints())
        val solution = solver.optimize(
            objectiveFunction,
            constraintSet,
            GoalType.MINIMIZE
        )
        
        return Mapqn_solution(
            objectiveValue = solution.value,
            variables = extractVariableValues(model, solution.point)
        )
    }
    
    private fun registerVariables(model: Mapqn_lpmodel, params: PFParameters) {
        val M = params.M
        val N = params.N
        
        // U variables
        for (i in 1..M) {
            model.addVariable("U_$i")
        }
        
        // Q variables
        for (i in 1..M) {
            model.addVariable("Q_$i")
        }
        
        // C variables
        for (j in 1..M) {
            for (i in 1..M) {
                model.addVariable("C_${j}_$i")
            }
        }
        
        // p1 variables
        for (j in 1..M) {
            for (i in 1..M) {
                for (ni in 0..N) {
                    model.addVariable("p1_${j}_${i}_$ni")
                }
            }
        }
        
        // p1c variables
        for (j in 1..M) {
            for (i in 1..M) {
                for (ni in 0..N) {
                    model.addVariable("p1c_${j}_${i}_$ni")
                }
            }
        }
    }
    
    private fun addDefinitionConstraints(model: Mapqn_lpmodel, params: PFParameters) {
        val M = params.M
        val N = params.N
        
        // ZER1: p1[j,j,0]=0
        for (j in 1..M) {
            val constraint = model.constraintBuilder()
                .addTerm("p1_${j}_${j}_0", 1.0)
                .eq(0.0)
            model.addConstraint(constraint)
        }
        
        // ZER3: p1[j,i,N]=0 for j<>i
        for (j in 1..M) {
            for (i in 1..M) {
                if (j != i) {
                    val constraint = model.constraintBuilder()
                        .addTerm("p1_${j}_${i}_$N", 1.0)
                        .eq(0.0)
                    model.addConstraint(constraint)
                }
            }
        }
        
        // CEQU: C[j,j] = Q[j]
        for (j in 1..M) {
            val constraint = model.constraintBuilder()
                .addTerm("C_${j}_$j", 1.0)
                .addTerm("Q_$j", -1.0)
                .eq(0.0)
            model.addConstraint(constraint)
        }
        
        // ONE1: sum{ni} (p1[j,i,ni]+p1c[j,i,ni])=1
        for (j in 1..M) {
            for (i in 1..M) {
                val constraint = model.constraintBuilder()
                for (ni in 0..N) {
                    constraint.addTerm("p1_${j}_${i}_$ni", 1.0)
                    constraint.addTerm("p1c_${j}_${i}_$ni", 1.0)
                }
                model.addConstraint(constraint.eq(1.0))
            }
        }
    }
    
    private fun addMeanIndicesConstraints(model: Mapqn_lpmodel, params: PFParameters) {
        val M = params.M
        val N = params.N
        
        // UTIL: U[i]=sum{t,nt} p1[i,t,nt]
        for (i in 1..M) {
            for (t in 1..M) {
                val constraint = model.constraintBuilder()
                constraint.addTerm("U_$i", 1.0)
                for (nt in 0..N) {
                    constraint.addTerm("p1_${i}_${t}_$nt", -1.0)
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
        
        // QLEN: Q[i]=sum{ni} ni*p1[i,i,ni]
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            constraint.addTerm("Q_$i", 1.0)
            for (ni in 0..N) {
                constraint.addTerm("p1_${i}_${i}_$ni", -ni.toDouble())
            }
            model.addConstraint(constraint.eq(0.0))
        }
        
        // CLEN: C[j,i]=sum{ni} ni*p1[j,i,ni]
        for (j in 1..M) {
            for (i in 1..M) {
                val constraint = model.constraintBuilder()
                constraint.addTerm("C_${j}_$i", 1.0)
                for (ni in 0..N) {
                    constraint.addTerm("p1_${j}_${i}_$ni", -ni.toDouble())
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
    }
    
    private fun addBalanceConstraints(model: Mapqn_lpmodel, params: PFParameters) {
        val M = params.M
        val N = params.N
        
        // MPCB: sum{i} C[j,i]=N*U[j]
        for (j in 1..M) {
            val constraint = model.constraintBuilder()
            for (i in 1..M) {
                constraint.addTerm("C_${j}_$i", 1.0)
            }
            constraint.addTerm("U_$j", -N.toDouble())
            model.addConstraint(constraint.eq(0.0))
        }
        
        // POPC: sum{i} Q[i]=N
        val popConstraint = model.constraintBuilder()
        for (i in 1..M) {
            popConstraint.addTerm("Q_$i", 1.0)
        }
        model.addConstraint(popConstraint.eq(N.toDouble()))
        
        // GFFL0: Global flow for ni=0
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            for (j in 1..M) {
                if (j != i) {
                    val qji = params.q(j - 1, i - 1)
                    val qij = params.q(i - 1, j - 1)
                    constraint.addTerm("p1_${j}_${i}_0", qji)
                    constraint.addTerm("p1_${i}_${i}_1", -qij)
                }
            }
            model.addConstraint(constraint.eq(0.0))
        }
        
        // GFFL: Global flow for ni in 1..N-1
        for (i in 1..M) {
            for (ni in 1..N-1) {
                val constraint = model.constraintBuilder()
                for (j in 1..M) {
                    if (j != i) {
                        val qji = params.q(j - 1, i - 1)
                        val qij = params.q(i - 1, j - 1)
                        constraint.addTerm("p1_${j}_${i}_$ni", qji)
                        constraint.addTerm("p1_${i}_${i}_${ni+1}", -qij)
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
        
        // UJNT: Joint probability symmetry
        for (i in 1..M) {
            for (j in 1..M) {
                val constraint = model.constraintBuilder()
                for (ni in 1..N) {
                    constraint.addTerm("p1_${j}_${i}_$ni", 1.0)
                }
                for (nj in 1..N) {
                    constraint.addTerm("p1_${i}_${j}_$nj", -1.0)
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
        
        // QBAL: Queue balance
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            for (j in 1..M) {
                if (j != i) {
                    val qij = params.q(i - 1, j - 1)
                    val qji = params.q(j - 1, i - 1)
                    
                    constraint.addTerm("U_$i", qij)
                    
                    // sum{nj} p1[i,j,nj]
                    for (nj in 1..N) {
                        constraint.addTerm("p1_${i}_${j}_$nj", -qji)
                    }
                    // p1[j,i,0]
                    constraint.addTerm("p1_${j}_${i}_0", -qji)
                }
            }
            model.addConstraint(constraint.eq(0.0))
        }
    }
    
    private fun extractVariableValues(model: Mapqn_lpmodel, point: DoubleArray): Map<String, Double> {
        val result = mutableMapOf<String, Double>()
        for ((name, index) in model.variables) {
            result[name] = point[index]
        }
        return result
    }
}

// Extension function to allow usage with Mapqn_solution
fun Mapqn_solution.getUtilizationPF(i: Int): Double {
    return getVariable("U_$i")
}

fun Mapqn_solution.getQueueLengthPF(i: Int): Double {
    return getVariable("Q_$i")
}