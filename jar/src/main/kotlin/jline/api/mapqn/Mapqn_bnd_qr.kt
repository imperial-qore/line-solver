/**
 * @file General Quadratic Reduction Bounds
 * 
 * Implements general quadratic reduction methods for computing performance
 * bounds in MAP queueing networks. Provides the core quadratic approximation
 * algorithms used across various MAPQN bound computation techniques.
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
 * Implementation of bnd_quadraticreduction.mod linear program
 * This model uses quadratic (p2) variables for joint probabilities
 */
object Mapqn_bnd_qr {
    
    /**
     * Create and solve the linear program for bnd_quadraticreduction model
     * 
     * @param params Model parameters
     * @param objectiveQueue Queue index to maximize utilization (1-based)
     * @param objectivePhase Phase index to maximize utilization (1-based)
     * @return Solution containing the optimal value and variable values
     */
    fun solve(params: LinearReductionParameters, objectiveQueue: Int, objectivePhase: Int): Mapqn_solution {
        params.validate()
        require(objectiveQueue in 1..params.M) { "Objective queue must be in range 1..M" }
        require(objectivePhase in 1..params.K[objectiveQueue - 1]) { 
            "Objective phase must be in range 1..K[${objectiveQueue - 1}]" 
        }
        
        val model = Mapqn_lpmodel()
        
        // Register all variables
        registerVariables(model, params)
        
        // Add all constraints
        addDefinitionConstraints(model, params)
        addMeanIndicesConstraints(model, params)
        addBalanceConstraints(model, params)
        addBoundConstraints(model, params)
        addQuadraticConstraints(model, params)
        
        // Create objective function - maximize U[objectiveQueue][objectivePhase]
        val objectiveVarName = "U_${objectiveQueue}_$objectivePhase"
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
    
    private fun registerVariables(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K
        
        // U variables
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                model.addVariable("U_${i}_$k")
            }
        }
        
        // IT variables
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                model.addVariable("IT_${i}_$k")
            }
        }
        
        // UP variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        model.addVariable("UP_${j}_${k}_${i}_$h")
                    }
                }
            }
        }
        
        // QP variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        model.addVariable("QP_${j}_${k}_${i}_$h")
                    }
                }
            }
        }
        
        // Q variables
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                model.addVariable("Q_${i}_$k")
            }
        }
        
        // C variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    model.addVariable("C_${j}_${k}_$i")
                }
            }
        }
        
        // I variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    model.addVariable("I_${j}_${k}_$i")
                }
            }
        }
        
        // p1 variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            model.addVariable("p1_${j}_${k}_${i}_${ni}_$h")
                        }
                    }
                }
            }
        }
        
        // p1c variables
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            model.addVariable("p1c_${j}_${k}_${i}_${ni}_$h")
                        }
                    }
                }
            }
        }
        
        // p2 variables (quadratic)
        for (j in 1..M) {
            for (nj in 0..N) {
                for (k in 1..K[j - 1]) {
                    for (i in 1..M) {
                        for (ni in 0..N) {
                            for (h in 1..K[i - 1]) {
                                model.addVariable("p2_${j}_${nj}_${k}_${i}_${ni}_$h")
                            }
                        }
                    }
                }
            }
        }
    }
    
    private fun addQuadraticConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K
        
        // PCL2: sum nj*ni*p2[i,ni,h,j,nj,k]=N*N
        val pcl2Constraint = model.constraintBuilder()
        for (i in 1..M) {
            for (j in 1..M) {
                for (ni in 1..N) {
                    for (nj in 1..N) {
                        for (h in 1..K[i - 1]) {
                            for (k in 1..K[j - 1]) {
                                pcl2Constraint.addTerm("p2_${i}_${ni}_${h}_${j}_${nj}_$k", 
                                    (nj * ni).toDouble())
                            }
                        }
                    }
                }
            }
        }
        model.addConstraint(pcl2Constraint.eq((N * N).toDouble()))
        
        // PI21: p1[j,k,i,ni,h] = sum{nj} p2[j,nj,k,i,ni,h]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            val constraint = model.constraintBuilder()
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", 1.0)
                            for (nj in 1..N) {
                                constraint.addTerm("p2_${j}_${nj}_${k}_${i}_${ni}_$h", -1.0)
                            }
                            model.addConstraint(constraint.eq(0.0))
                        }
                    }
                }
            }
        }
        
        // PI22: p1c[j,k,i,ni,h] = p2[j,0,k,i,ni,h]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            val constraint = model.constraintBuilder()
                                .addTerm("p1c_${j}_${k}_${i}_${ni}_$h", 1.0)
                                .addTerm("p2_${j}_0_${k}_${i}_${ni}_$h", -1.0)
                                .eq(0.0)
                            model.addConstraint(constraint)
                        }
                    }
                }
            }
        }
        
        // PI23: p2[i,ni,h,j,nj,k] = p2[j,nj,k,i,ni,h] (symmetry)
        for (j in 1..M) {
            for (nj in 0..N) {
                for (k in 1..K[j - 1]) {
                    for (i in 1..M) {
                        for (ni in 0..N) {
                            for (h in 1..K[i - 1]) {
                                val constraint = model.constraintBuilder()
                                    .addTerm("p2_${i}_${ni}_${h}_${j}_${nj}_$k", 1.0)
                                    .addTerm("p2_${j}_${nj}_${k}_${i}_${ni}_$h", -1.0)
                                    .eq(0.0)
                                model.addConstraint(constraint)
                            }
                        }
                    }
                }
            }
        }
    }
    
    private fun addDefinitionConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K
        
        // ZER1: p1[j,k,j,0,k]=0
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                val constraint = model.constraintBuilder()
                    .addTerm("p1_${j}_${k}_${j}_0_$k", 1.0)
                    .eq(0.0)
                model.addConstraint(constraint)
            }
        }
        
        // ZER2: p1[j,k,j,nj,h]=0 for h<>k
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (nj in 0..N) {
                    for (h in 1..K[j - 1]) {
                        if (h != k) {
                            val constraint = model.constraintBuilder()
                                .addTerm("p1_${j}_${k}_${j}_${nj}_$h", 1.0)
                                .eq(0.0)
                            model.addConstraint(constraint)
                        }
                    }
                }
            }
        }
        
        // ZER3: p1[j,k,i,N,h]=0 for j<>i
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        if (j != i) {
                            val constraint = model.constraintBuilder()
                                .addTerm("p1_${j}_${k}_${i}_${N}_$h", 1.0)
                                .eq(0.0)
                            model.addConstraint(constraint)
                        }
                    }
                }
            }
        }
        
        // ZER4: p1c[j,k,j,nj,h]=0 for nj>=1
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (nj in 1..N) {
                    for (h in 1..K[j - 1]) {
                        val constraint = model.constraintBuilder()
                            .addTerm("p1c_${j}_${k}_${j}_${nj}_$h", 1.0)
                            .eq(0.0)
                        model.addConstraint(constraint)
                    }
                }
            }
        }
        
        // CEQU: C[j,k,j] = Q[j,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                val constraint = model.constraintBuilder()
                    .addTerm("C_${j}_${k}_$j", 1.0)
                    .addTerm("Q_${j}_$k", -1.0)
                    .eq(0.0)
                model.addConstraint(constraint)
            }
        }
        
        // ONE1: sum over all p1 and p1c = 1
        for (j in 1..M) {
            for (i in 1..M) {
                val constraint = model.constraintBuilder()
                for (k in 1..K[j - 1]) {
                    for (h in 1..K[i - 1]) {
                        for (ni in 0..N) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", 1.0)
                            constraint.addTerm("p1c_${j}_${k}_${i}_${ni}_$h", 1.0)
                        }
                    }
                }
                model.addConstraint(constraint.eq(1.0))
            }
        }
    }
    
    private fun addMeanIndicesConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        // Same as in Mapqn_bnd_lr
        val M = params.M
        val N = params.N
        val K = params.K
        
        // UTLB: U[i,k]=sum{nt,h} p1[i,k,t,nt,h]
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                for (t in 1..M) {
                    val constraint = model.constraintBuilder()
                    constraint.addTerm("U_${i}_$k", 1.0)
                    for (nt in 0..N) {
                        for (h in 1..K[t - 1]) {
                            constraint.addTerm("p1_${i}_${k}_${t}_${nt}_$h", -1.0)
                        }
                    }
                    model.addConstraint(constraint.eq(0.0))
                }
            }
        }
        
        // UTLC: IT[i,k]=sum{nt,h} p1c[i,k,t,nt,h]
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                for (t in 1..M) {
                    val constraint = model.constraintBuilder()
                    constraint.addTerm("IT_${i}_$k", 1.0)
                    for (nt in 0..N) {
                        for (h in 1..K[t - 1]) {
                            constraint.addTerm("p1c_${i}_${k}_${t}_${nt}_$h", -1.0)
                        }
                    }
                    model.addConstraint(constraint.eq(0.0))
                }
            }
        }
        
        // QLEN: Q[i,k]=sum{ni} ni*p1[i,k,i,ni,k]
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                constraint.addTerm("Q_${i}_$k", 1.0)
                for (ni in 0..N) {
                    constraint.addTerm("p1_${i}_${k}_${i}_${ni}_$k", -ni.toDouble())
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
        
        // Additional constraints would follow...
    }
    
    private fun addBalanceConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        // Similar to Mapqn_bnd_lr but with modifications for quadratic model
        val M = params.M
        val N = params.N
        val K = params.K
        
        // SRVB: Service balance
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                for (j in 1..M) {
                    for (h in 1..K[i - 1]) {
                        val qCoeff1 = params.q(i - 1, j - 1, k - 1, h - 1)
                        val qCoeff2 = params.q(i - 1, j - 1, h - 1, k - 1)
                        constraint.addTerm("U_${i}_$k", qCoeff1)
                        constraint.addTerm("U_${i}_$h", -qCoeff2)
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }
        
        // POPC: Population conservation
        val popConstraint = model.constraintBuilder()
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                popConstraint.addTerm("Q_${i}_$k", 1.0)
            }
        }
        model.addConstraint(popConstraint.eq(N.toDouble()))
        
        // ONE: sum{k} (U[j,k]+IT[j,k])=1
        for (j in 1..M) {
            val constraint = model.constraintBuilder()
            for (k in 1..K[j - 1]) {
                constraint.addTerm("U_${j}_$k", 1.0)
                constraint.addTerm("IT_${j}_$k", 1.0)
            }
            model.addConstraint(constraint.eq(1.0))
        }
    }
    
    private fun addBoundConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K
        
        // UUB1: sum{k} U[i,k]<=1
        for (i in 1..M) {
            val constraint = model.constraintBuilder()
            for (k in 1..K[i - 1]) {
                constraint.addTerm("U_${i}_$k", 1.0)
            }
            model.addConstraint(constraint.leq(1.0))
        }
        
        // QUB1: Q[j,k] <= N*U[j,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                val constraint = model.constraintBuilder()
                    .addTerm("Q_${j}_$k", 1.0)
                    .addTerm("U_${j}_$k", -N.toDouble())
                    .leq(0.0)
                model.addConstraint(constraint)
            }
        }
        
        // CUB1: C[j,k,i] <= sum{h} Q[i,h]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                    constraint.addTerm("C_${j}_${k}_$i", 1.0)
                    for (h in 1..K[i - 1]) {
                        constraint.addTerm("Q_${i}_$h", -1.0)
                    }
                    model.addConstraint(constraint.leq(0.0))
                }
            }
        }
        
        // CUB2: C[j,k,i] <= N*U[j,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                        .addTerm("C_${j}_${k}_$i", 1.0)
                        .addTerm("U_${j}_$k", -N.toDouble())
                        .leq(0.0)
                    model.addConstraint(constraint)
                }
            }
        }
        
        // QMIN constraint
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                    
                    // Left side: sum over all queues
                    for (t in 1..M) {
                        for (h in 1..K[t - 1]) {
                            for (nt in 0..N) {
                                constraint.addTerm("p1_${j}_${k}_${t}_${nt}_$h", nt.toDouble())
                                constraint.addTerm("p1c_${j}_${k}_${t}_${nt}_$h", nt.toDouble())
                            }
                        }
                    }
                    
                    // Right side: N * sum for specific queue i
                    for (h in 1..K[i - 1]) {
                        for (ni in 0..N) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", -N.toDouble())
                            constraint.addTerm("p1c_${j}_${k}_${i}_${ni}_$h", -N.toDouble())
                        }
                    }
                    
                    model.addConstraint(constraint.geq(0.0))
                }
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