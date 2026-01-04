/**
 * @file Linear Reduction Bounds via MVA
 * 
 * Implements linear reduction bounds for MAP queueing networks using Mean
 * Value Analysis (MVA) techniques. Combines linear approximation methods
 * with MVA algorithms for efficient MAPQN performance evaluation.
 * 
 * @since LINE 3.0
 */
package jline.api.mapqn

import org.apache.commons.math3.optim.linear.LinearConstraintSet
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction
import org.apache.commons.math3.optim.linear.SimplexSolver
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType

/**
 * Implementation of bnd_mvaversion.mod linear program
 */
object Mapqn_bnd_lr_mva {
    
    /**
     * Create and solve the linear program for bnd_mvaversion model
     * 
     * @param params Model parameters
     * @param objectiveQueue Queue index to maximize utilization (1-based)
     * @param objectiveLevel Level index to maximize utilization (1-based)
     * @return Solution containing the optimal value and variable values
     */
    fun solve(params: MVAVersionParameters, objectiveQueue: Int, objectiveLevel: Int): Mapqn_solution {
        params.validate()
        require(objectiveQueue in 1..params.M) { "Objective queue must be in range 1..M" }
        require(objectiveLevel in 1..params.K) { "Objective level must be in range 1..K" }
        
        val model = Mapqn_lpmodel()
        
        // Register all variables
        registerVariables(model, params)
        
        // Add all constraints
        addConstraints(model, params)
        
        // Create objective function - maximize UN[objectiveQueue][objectiveLevel]
        val objectiveVarName = "UN_${objectiveQueue}_$objectiveLevel"
        val objectiveCoeffs = model.createObjectiveCoefficients(objectiveVarName)
        val objectiveFunction = LinearObjectiveFunction(objectiveCoeffs, 0.0)
        
        // Solve the LP
        val solver = SimplexSolver()
        val constraintSet = LinearConstraintSet(model.getConstraints())
        val solution = solver.optimize(
            objectiveFunction,
            constraintSet,
            GoalType.MAXIMIZE
        )
        
        return Mapqn_solution(
            objectiveValue = solution.value,
            variables = extractVariableValues(model, solution.point)
        )
    }
    
    private fun registerVariables(model: Mapqn_lpmodel, params: MVAVersionParameters) {
        val M = params.M
        val K = params.K
        
        // UN variables (utilization)
        for (i in 1..M) {
            for (k in 1..K) {
                model.addVariable("UN_${i}_$k")
            }
        }
        
        // QN variables (queue length)
        for (i in 1..M) {
            for (k in 1..K) {
                model.addVariable("QN_${i}_$k")
            }
        }
        
        // B variables
        for (j in 1..M) {
            for (k in 1..K) {
                for (i in 1..M) {
                    model.addVariable("B_${j}_${k}_$i")
                }
            }
        }
    }
    
    private fun addConstraints(model: Mapqn_lpmodel, params: MVAVersionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K
        
        // QNB: QN[i,k]>=B[j,k,i]
        for (i in 1..M) {
            for (k in 1..K) {
                for (j in 1..M) {
                    val constraint = model.constraintBuilder()
                        .addTerm("QN_${i}_$k", 1.0)
                        .addTerm("B_${j}_${k}_$i", -1.0)
                        .geq(0.0)
                    model.addConstraint(constraint)
                }
            }
        }
        
        // UMAX: sum{k} UN[i,k] <=1
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            for (k in 1..K) {
                constraint.addTerm("UN_${i}_$k", 1.0)
            }
            model.addConstraint(constraint.leq(1.0))
        }
        
        // POPCONSTR: sum{i,k} QN[i,k]=N
        val popConstraint = model.constraintBuilder()
        for (i in 1..M) {
            for (k in 1..K) {
                popConstraint.addTerm("QN_${i}_$k", 1.0)
            }
        }
        model.addConstraint(popConstraint.eq(N.toDouble()))
        
        // FLOW: Flow balance equations
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            for (k in 1..K) {
                for (m in 1..K) {
                    for (w in 1..M) {
                        val qIn = params.q(w - 1, i - 1, k - 1, m - 1)
                        val qOut = params.q(i - 1, w - 1, m - 1, k - 1)
                        
                        constraint.addTerm("UN_${w}_$k", qIn)
                        constraint.addTerm("UN_${i}_$m", -qOut)
                    }
                }
            }
            model.addConstraint(constraint.eq(0.0))
        }
        
        // UBAL: Balance for MAP queue levels
        for (k in 1..K) {
            val constraint = model.constraintBuilder()
            for (h in 1..K) {
                if (h != k) {
                    for (w in 1..M) {
                        val qOut = params.q(M - 1, w - 1, k - 1, h - 1)
                        val qIn = params.q(M - 1, w - 1, h - 1, k - 1)
                        
                        constraint.addTerm("UN_${M}_$k", qOut)
                        constraint.addTerm("UN_${M}_$h", -qIn)
                    }
                }
            }
            model.addConstraint(constraint.eq(0.0))
        }
        
        // QBAL: Queue balance for MAP queue
        for (k in 1..K) {
            val constraint = model.constraintBuilder()
            
            // Outgoing from level k
            for (h in 1..K) {
                if (h != k) {
                    for (w in 1..M) {
                        val q = params.q(M - 1, w - 1, k - 1, h - 1)
                        constraint.addTerm("QN_${M}_$k", q)
                    }
                }
            }
            
            // Incoming to level k from other queues
            for (m in 1..K) {
                for (j in 1..M-1) {
                    val q = params.q(M - 1, j - 1, m - 1, k - 1)
                    constraint.addTerm("UN_${M}_$m", q)
                }
            }
            
            // Incoming from other queues j to MAP
            for (j in 1..M-1) {
                val q = params.q(j - 1, M - 1, k - 1, k - 1)
                constraint.addTerm("UN_${j}_$k", -q)
            }
            
            // Incoming from other levels
            for (h in 1..K) {
                if (h != k) {
                    for (w in 1..M) {
                        val q = params.q(M - 1, w - 1, h - 1, k - 1)
                        constraint.addTerm("QN_${M}_$h", -q)
                    }
                }
            }
            
            model.addConstraint(constraint.eq(0.0))
        }
        
        // MCC: Mean customer count constraint
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            
            // First sum: outgoing from queue i
            for (k in 1..K) {
                for (m in 1..K) {
                    for (w in 1..M) {
                        if (w != i) {
                            val q = params.q(i - 1, w - 1, k - 1, m - 1)
                            constraint.addTerm("QN_${i}_$k", q)
                        }
                    }
                }
            }
            
            // Second sum: incoming to queue i
            for (k in 1..K) {
                for (m in 1..K) {
                    for (j in 1..M) {
                        if (j != i) {
                            val q = params.q(j - 1, i - 1, k - 1, m - 1)
                            constraint.addTerm("QN_${j}_$k", q)
                        }
                    }
                }
            }
            
            // Third sum: B variables
            for (k in 1..K) {
                for (m in 1..K) {
                    for (j in 1..M) {
                        if (j != i) {
                            val q = params.q(j - 1, i - 1, k - 1, m - 1)
                            for (wp in 1..M) {
                                if (wp != i && wp != j) {
                                    constraint.addTerm("B_${j}_${k}_$wp", q)
                                }
                            }
                        }
                    }
                }
            }
            
            // RHS: (N+1) * incoming utilization
            for (k in 1..K) {
                for (m in 1..K) {
                    for (j in 1..M) {
                        if (j != i) {
                            val q = params.q(j - 1, i - 1, k - 1, m - 1)
                            constraint.addTerm("UN_${j}_$k", -(N + 1) * q)
                        }
                    }
                }
            }
            
            model.addConstraint(constraint.eq(0.0))
        }
        
        // MCC2: Second mean customer count constraint
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            
            // LHS: outgoing
            for (k in 1..K) {
                for (m in 1..K) {
                    for (w in 1..M) {
                        if (w != i) {
                            val q = params.q(i - 1, w - 1, k - 1, m - 1)
                            constraint.addTerm("QN_${i}_$k", q)
                        }
                    }
                }
            }
            
            // RHS: incoming
            for (k in 1..K) {
                for (m in 1..K) {
                    for (j in 1..M) {
                        if (j != i) {
                            val q = params.q(j - 1, i - 1, k - 1, m - 1)
                            constraint.addTerm("B_${j}_${k}_$i", -q)
                            constraint.addTerm("UN_${j}_$k", -q)
                        }
                    }
                }
            }
            
            model.addConstraint(constraint.eq(0.0))
        }
        
        // QMAX: QN[w,k] <= N*UN[w,k]
        for (w in 1..M) {
            for (k in 1..K) {
                val constraint = model.constraintBuilder()
                    .addTerm("QN_${w}_$k", 1.0)
                    .addTerm("UN_${w}_$k", -N.toDouble())
                    .leq(0.0)
                model.addConstraint(constraint)
            }
        }
        
        // QMIN: sum{w} QN[w,k] >= N*UN[j,k]
        for (k in 1..K) {
            for (j in 1..M) {
                val constraint = model.constraintBuilder()
                for (w in 1..M) {
                    constraint.addTerm("QN_${w}_$k", 1.0)
                }
                constraint.addTerm("UN_${j}_$k", -N.toDouble())
                model.addConstraint(constraint.geq(0.0))
            }
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

// Extension functions for MVA version
fun Mapqn_solution.getUtilizationMVA(i: Int, k: Int): Double {
    return getVariable("UN_${i}_$k")
}

fun Mapqn_solution.getQueueLengthMVA(i: Int, k: Int): Double {
    return getVariable("QN_${i}_$k")
}