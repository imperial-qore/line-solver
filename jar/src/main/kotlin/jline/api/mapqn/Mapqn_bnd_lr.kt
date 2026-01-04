/**
 * @file General Linear Reduction Bounds
 * 
 * Implements general linear reduction methods for computing performance
 * bounds in MAP queueing networks. Provides the core linear approximation
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
import org.apache.commons.math3.optim.PointValuePair

/**
 * Implementation of bnd_linearreduction_new.mod linear program
 */
object Mapqn_bnd_lr {
    
    /**
     * Create and solve the linear program for bnd_linearreduction_new model
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
        val M = params.M
        val N = params.N
        val K = params.K

        // UTLB: U[i,k]=sum{nt,h} p1[i,k,t,nt,h]  (for each t)
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

        // UTLC: IT[i,k]=sum{nt,h} p1c[i,k,t,nt,h]  (for each t)
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

        // UPH1: UP[j,k,i,h] = sum{ni in 1..N} (p1[j,k,i,ni,h]+p1c[j,k,i,ni,h])
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        val constraint = model.constraintBuilder()
                        constraint.addTerm("UP_${j}_${k}_${i}_$h", 1.0)
                        for (ni in 1..N) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", -1.0)
                            constraint.addTerm("p1c_${j}_${k}_${i}_${ni}_$h", -1.0)
                        }
                        model.addConstraint(constraint.eq(0.0))
                    }
                }
            }
        }

        // UPH2: UP[j,k,i,h] = (sum{ni in 1..N} p1[j,k,i,ni,h]) + p1[i,h,j,0,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        val constraint = model.constraintBuilder()
                        constraint.addTerm("UP_${j}_${k}_${i}_$h", 1.0)
                        for (ni in 1..N) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", -1.0)
                        }
                        constraint.addTerm("p1_${i}_${h}_${j}_0_$k", -1.0)
                        model.addConstraint(constraint.eq(0.0))
                    }
                }
            }
        }

        // QPH1: QP[j,k,i,h] = sum{ni in 1..N} ni*(p1[j,k,i,ni,h]+p1c[j,k,i,ni,h])
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    for (h in 1..K[i - 1]) {
                        val constraint = model.constraintBuilder()
                        constraint.addTerm("QP_${j}_${k}_${i}_$h", 1.0)
                        for (ni in 1..N) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", -ni.toDouble())
                            constraint.addTerm("p1c_${j}_${k}_${i}_${ni}_$h", -ni.toDouble())
                        }
                        model.addConstraint(constraint.eq(0.0))
                    }
                }
            }
        }

        // CLEN: C[j,k,i] = sum{ni in 0..N} sum{h in 1..K[i]} ni*p1[j,k,i,ni,h]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                    constraint.addTerm("C_${j}_${k}_$i", 1.0)
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", -ni.toDouble())
                        }
                    }
                    model.addConstraint(constraint.eq(0.0))
                }
            }
        }

        // ILEN: I[j,k,i] = sum{ni in 0..N} sum{h in 1..K[i]} ni*p1c[j,k,i,ni,h]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                    constraint.addTerm("I_${j}_${k}_$i", 1.0)
                    for (ni in 0..N) {
                        for (h in 1..K[i - 1]) {
                            constraint.addTerm("p1c_${j}_${k}_${i}_${ni}_$h", -ni.toDouble())
                        }
                    }
                    model.addConstraint(constraint.eq(0.0))
                }
            }
        }
    }
    
    private fun addBalanceConstraints(model: Mapqn_lpmodel, params: LinearReductionParameters) {
        val M = params.M
        val N = params.N
        val K = params.K

        // SRVB: Service balance
        // sum{j in 1..M} sum{h in 1..K[i]:h<>k} q[i,j,k,h]*U[i,k] = sum{j in 1..M} sum{h in 1..K[i]:h<>k} q[i,j,h,k]*U[i,h]
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                for (j in 1..M) {
                    for (h in 1..K[i - 1]) {
                        if (h != k) {
                            val qCoeff = params.q(i - 1, j - 1, k - 1, h - 1)
                            constraint.addTerm("U_${i}_$k", qCoeff)
                            constraint.addTerm("U_${i}_$h", -params.q(i - 1, j - 1, h - 1, k - 1))
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // MPCB: sum{i in 1..M} C[j,k,i] = N*U[j,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                val constraint = model.constraintBuilder()
                for (i in 1..M) {
                    constraint.addTerm("C_${j}_${k}_$i", 1.0)
                }
                constraint.addTerm("U_${j}_$k", -N.toDouble())
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // MPCI: sum{i in 1..M} I[j,k,i] = N*IT[j,k]
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                val constraint = model.constraintBuilder()
                for (i in 1..M) {
                    constraint.addTerm("I_${j}_${k}_$i", 1.0)
                }
                constraint.addTerm("IT_${j}_$k", -N.toDouble())
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // ONE: sum{k} (U[j,k]+IT[j,k])=1
        for (j in 1..M) {
            val constraint = model.constraintBuilder()
            for (k in 1..K[j - 1]) {
                constraint.addTerm("U_${j}_$k", 1.0)
                constraint.addTerm("IT_${j}_$k", 1.0)
            }
            model.addConstraint(constraint.eq(1.0))
        }

        // POPC: Population conservation - sum{i,k} Q[i,k] = N
        val popConstraint = model.constraintBuilder()
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                popConstraint.addTerm("Q_${i}_$k", 1.0)
            }
        }
        model.addConstraint(popConstraint.eq(N.toDouble()))

        // GFFL0: Global flow balance at ni=0
        // sum{j<>i} sum{k in K[j]} sum{h in K[j]} q[j,i,k,h]*p1[j,k,i,0,u] = sum{j<>i} sum{k in K[i]} q[i,j,k,u]*p1[i,k,i,1,k]
        for (i in 1..M) {
            for (u in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                // Left side: arrivals to queue i at level 0
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[j - 1]) {
                            for (h in 1..K[j - 1]) {
                                val qVal = params.q(j - 1, i - 1, k - 1, h - 1)
                                constraint.addTerm("p1_${j}_${k}_${i}_0_$u", qVal)
                            }
                        }
                    }
                }
                // Right side: departures from queue i at level 1
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[i - 1]) {
                            val qVal = params.q(i - 1, j - 1, k - 1, u - 1)
                            constraint.addTerm("p1_${i}_${k}_${i}_1_$k", -qVal)
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // GFFL: Global flow balance at ni=1..N-1
        // sum{j<>i} sum{k in K[j]} sum{h in K[j]} sum{u in K[i]} q[j,i,k,h]*p1[j,k,i,ni,u] =
        //   sum{j<>i} sum{k in K[i]} sum{h in K[i]} q[i,j,k,h]*p1[i,k,i,ni+1,k]
        for (i in 1..M) {
            for (ni in 1..N - 1) {
                val constraint = model.constraintBuilder()
                // Left side
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[j - 1]) {
                            for (h in 1..K[j - 1]) {
                                for (u in 1..K[i - 1]) {
                                    val qVal = params.q(j - 1, i - 1, k - 1, h - 1)
                                    constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$u", qVal)
                                }
                            }
                        }
                    }
                }
                // Right side
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[i - 1]) {
                            for (h in 1..K[i - 1]) {
                                val qVal = params.q(i - 1, j - 1, k - 1, h - 1)
                                constraint.addTerm("p1_${i}_${k}_${i}_${ni + 1}_$k", -qVal)
                            }
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // UJNT: Joint utilization constraint
        // sum{ni in 1..N} sum{k in K[j]} sum{h in K[i]} p1[j,k,i,ni,h] =
        //   sum{nj in 1..N} sum{k in K[j]} sum{h in K[i]} p1[i,h,j,nj,k]
        for (i in 1..M) {
            for (j in 1..M) {
                val constraint = model.constraintBuilder()
                // Left side
                for (ni in 1..N) {
                    for (k in 1..K[j - 1]) {
                        for (h in 1..K[i - 1]) {
                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$h", 1.0)
                        }
                    }
                }
                // Right side
                for (nj in 1..N) {
                    for (k in 1..K[j - 1]) {
                        for (h in 1..K[i - 1]) {
                            constraint.addTerm("p1_${i}_${h}_${j}_${nj}_$k", -1.0)
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // QBAL: Queue balance constraint
        // sum{h<>k} sum{j} q[i,j,k,h]*Q[i,k] + sum{j<>i} sum{h} q[i,j,h,k]*U[i,h] =
        //   sum{j<>i} sum{u in K[j]} sum{w in K[j]} q[j,i,u,w]*(sum{nj in 1..N} p1[i,k,j,nj,u]+p1[j,u,i,0,k])
        //   + sum{h<>k} sum{j} q[i,j,h,k]*Q[i,h]
        for (i in 1..M) {
            for (k in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                // Term 1: sum{h<>k} sum{j} q[i,j,k,h]*Q[i,k]
                for (h in 1..K[i - 1]) {
                    if (h != k) {
                        for (j in 1..M) {
                            val qVal = params.q(i - 1, j - 1, k - 1, h - 1)
                            constraint.addTerm("Q_${i}_$k", qVal)
                        }
                    }
                }
                // Term 2: sum{j<>i} sum{h} q[i,j,h,k]*U[i,h]
                for (j in 1..M) {
                    if (j != i) {
                        for (h in 1..K[i - 1]) {
                            val qVal = params.q(i - 1, j - 1, h - 1, k - 1)
                            constraint.addTerm("U_${i}_$h", qVal)
                        }
                    }
                }
                // Term 3 (RHS, negated): sum{j<>i} sum{u,w} q[j,i,u,w]*(sum{nj} p1[i,k,j,nj,u] + p1[j,u,i,0,k])
                for (j in 1..M) {
                    if (j != i) {
                        for (u in 1..K[j - 1]) {
                            for (w in 1..K[j - 1]) {
                                val qVal = params.q(j - 1, i - 1, u - 1, w - 1)
                                for (nj in 1..N) {
                                    constraint.addTerm("p1_${i}_${k}_${j}_${nj}_$u", -qVal)
                                }
                                constraint.addTerm("p1_${j}_${u}_${i}_0_$k", -qVal)
                            }
                        }
                    }
                }
                // Term 4 (RHS, negated): sum{h<>k} sum{j} q[i,j,h,k]*Q[i,h]
                for (h in 1..K[i - 1]) {
                    if (h != k) {
                        for (j in 1..M) {
                            val qVal = params.q(i - 1, j - 1, h - 1, k - 1)
                            constraint.addTerm("Q_${i}_$h", -qVal)
                        }
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
        }

        // MBH: Marginal balance for ni in 0..N-2
        // This constraint is complex - see AMPL model for reference
        for (i in 1..M) {
            for (ordphase in 1..K[i - 1]) {
                for (ni in 0..N - 2) {
                    val constraint = model.constraintBuilder()
                    // LHS Term 1: sum{j<>i} sum{k in K[j]} sum{h in K[j]} sum{u in K[i]:u<>ordphase} q[j,i,k,h]*p1[j,k,i,ni,u]
                    for (j in 1..M) {
                        if (j != i) {
                            for (k in 1..K[j - 1]) {
                                for (h in 1..K[j - 1]) {
                                    for (u in 1..K[i - 1]) {
                                        if (u != ordphase) {
                                            val qVal = params.q(j - 1, i - 1, k - 1, h - 1)
                                            constraint.addTerm("p1_${j}_${k}_${i}_${ni}_$u", qVal)
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // LHS Term 2: sum{j<>i} sum{k in K[j]} sum{h in K[j]} q[j,i,k,h]*p1[j,k,i,ni+1,ordphase]
                    for (j in 1..M) {
                        if (j != i) {
                            for (k in 1..K[j - 1]) {
                                for (h in 1..K[j - 1]) {
                                    val qVal = params.q(j - 1, i - 1, k - 1, h - 1)
                                    constraint.addTerm("p1_${j}_${k}_${i}_${ni + 1}_$ordphase", qVal)
                                }
                            }
                        }
                    }
                    // LHS Term 3: sum{k in K[i]:k<>ordphase} q[i,i,ordphase,k]*p1[i,ordphase,i,ni+1,ordphase]
                    for (k in 1..K[i - 1]) {
                        if (k != ordphase) {
                            val qVal = params.q(i - 1, i - 1, ordphase - 1, k - 1)
                            constraint.addTerm("p1_${i}_${ordphase}_${i}_${ni + 1}_$ordphase", qVal)
                        }
                    }
                    // RHS Term 1 (negated): sum{j<>i} sum{k in K[i]:k<>ordphase} q[i,j,k,k]*p1[i,k,i,ni+1,k]
                    for (j in 1..M) {
                        if (j != i) {
                            for (k in 1..K[i - 1]) {
                                if (k != ordphase) {
                                    val qVal = params.q(i - 1, j - 1, k - 1, k - 1)
                                    constraint.addTerm("p1_${i}_${k}_${i}_${ni + 1}_$k", -qVal)
                                }
                            }
                        }
                    }
                    // RHS Term 2 (negated): sum{j<>i} sum{k in K[i]:k<>ordphase} sum{h in K[i]:h<>k} q[i,j,k,h]*p1[i,k,i,ni+1,k]
                    for (j in 1..M) {
                        if (j != i) {
                            for (k in 1..K[i - 1]) {
                                if (k != ordphase) {
                                    for (h in 1..K[i - 1]) {
                                        if (h != k) {
                                            val qVal = params.q(i - 1, j - 1, k - 1, h - 1)
                                            constraint.addTerm("p1_${i}_${k}_${i}_${ni + 1}_$k", -qVal)
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // RHS Term 3 (negated): sum{j<>i} sum{k in K[i]:k<>ordphase} q[i,j,k,ordphase]*p1[i,k,i,ni+2,k]
                    for (j in 1..M) {
                        if (j != i) {
                            for (k in 1..K[i - 1]) {
                                if (k != ordphase) {
                                    val qVal = params.q(i - 1, j - 1, k - 1, ordphase - 1)
                                    constraint.addTerm("p1_${i}_${k}_${i}_${ni + 2}_$k", -qVal)
                                }
                            }
                        }
                    }
                    // RHS Term 4 (negated): sum{j<>i} q[i,j,ordphase,ordphase]*p1[i,ordphase,i,ni+2,ordphase]
                    for (j in 1..M) {
                        if (j != i) {
                            val qVal = params.q(i - 1, j - 1, ordphase - 1, ordphase - 1)
                            constraint.addTerm("p1_${i}_${ordphase}_${i}_${ni + 2}_$ordphase", -qVal)
                        }
                    }
                    // RHS Term 5 (negated): sum{k in K[i]:k<>ordphase} q[i,i,k,ordphase]*p1[i,k,i,ni+1,k]
                    for (k in 1..K[i - 1]) {
                        if (k != ordphase) {
                            val qVal = params.q(i - 1, i - 1, k - 1, ordphase - 1)
                            constraint.addTerm("p1_${i}_${k}_${i}_${ni + 1}_$k", -qVal)
                        }
                    }
                    model.addConstraint(constraint.eq(0.0))
                }
            }
        }

        // MBHN: Marginal balance at N (boundary)
        for (i in 1..M) {
            for (ordphase in 1..K[i - 1]) {
                val constraint = model.constraintBuilder()
                // LHS Term 1: sum{j<>i} sum{k in K[j]} sum{h in K[j]} sum{u in K[i]:u<>ordphase} q[j,i,k,h]*p1[j,k,i,N-1,u]
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[j - 1]) {
                            for (h in 1..K[j - 1]) {
                                for (u in 1..K[i - 1]) {
                                    if (u != ordphase) {
                                        val qVal = params.q(j - 1, i - 1, k - 1, h - 1)
                                        constraint.addTerm("p1_${j}_${k}_${i}_${N - 1}_$u", qVal)
                                    }
                                }
                            }
                        }
                    }
                }
                // LHS Term 2: sum{k in K[i]:k<>ordphase} q[i,i,ordphase,k]*p1[i,ordphase,i,N,ordphase]
                for (k in 1..K[i - 1]) {
                    if (k != ordphase) {
                        val qVal = params.q(i - 1, i - 1, ordphase - 1, k - 1)
                        constraint.addTerm("p1_${i}_${ordphase}_${i}_${N}_$ordphase", qVal)
                    }
                }
                // RHS Term 1 (negated): sum{j<>i} sum{k in K[i]:k<>ordphase} q[i,j,k,k]*p1[i,k,i,N,k]
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[i - 1]) {
                            if (k != ordphase) {
                                val qVal = params.q(i - 1, j - 1, k - 1, k - 1)
                                constraint.addTerm("p1_${i}_${k}_${i}_${N}_$k", -qVal)
                            }
                        }
                    }
                }
                // RHS Term 2 (negated): sum{j<>i} sum{k in K[i]:k<>ordphase} sum{h in K[i]:h<>k} q[i,j,k,h]*p1[i,k,i,N,k]
                for (j in 1..M) {
                    if (j != i) {
                        for (k in 1..K[i - 1]) {
                            if (k != ordphase) {
                                for (h in 1..K[i - 1]) {
                                    if (h != k) {
                                        val qVal = params.q(i - 1, j - 1, k - 1, h - 1)
                                        constraint.addTerm("p1_${i}_${k}_${i}_${N}_$k", -qVal)
                                    }
                                }
                            }
                        }
                    }
                }
                // RHS Term 3 (negated): sum{k in K[i]:k<>ordphase} q[i,i,k,ordphase]*p1[i,k,i,N,k]
                for (k in 1..K[i - 1]) {
                    if (k != ordphase) {
                        val qVal = params.q(i - 1, i - 1, k - 1, ordphase - 1)
                        constraint.addTerm("p1_${i}_${k}_${i}_${N}_$k", -qVal)
                    }
                }
                model.addConstraint(constraint.eq(0.0))
            }
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

        // QMIN: sum{t,h,nt} nt*(p1[j,k,t,nt,h]+p1c[j,k,t,nt,h]) >= N*sum{h,ni} (p1[j,k,i,ni,h]+p1c[j,k,i,ni,h])
        for (j in 1..M) {
            for (k in 1..K[j - 1]) {
                for (i in 1..M) {
                    val constraint = model.constraintBuilder()
                    // Left side: sum{t,h,nt} nt*(p1[j,k,t,nt,h]+p1c[j,k,t,nt,h])
                    for (t in 1..M) {
                        for (h in 1..K[t - 1]) {
                            for (nt in 0..N) {
                                constraint.addTerm("p1_${j}_${k}_${t}_${nt}_$h", nt.toDouble())
                                constraint.addTerm("p1c_${j}_${k}_${t}_${nt}_$h", nt.toDouble())
                            }
                        }
                    }
                    // Right side (negated): N*sum{h,ni} (p1[j,k,i,ni,h]+p1c[j,k,i,ni,h])
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

        // TEST: p1[j,k,j,nj,k] >= p1[i,h,j,nj,k] for j<>i, nj in 1..N-1
        for (j in 1..M) {
            for (nj in 1..N - 1) {
                for (k in 1..K[j - 1]) {
                    for (i in 1..M) {
                        if (i != j) {
                            for (h in 1..K[i - 1]) {
                                val constraint = model.constraintBuilder()
                                    .addTerm("p1_${j}_${k}_${j}_${nj}_$k", 1.0)
                                    .addTerm("p1_${i}_${h}_${j}_${nj}_$k", -1.0)
                                    .geq(0.0)
                                model.addConstraint(constraint)
                            }
                        }
                    }
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

/**
 * Solution container for MAPQN linear programs
 */
data class Mapqn_solution(
    val objectiveValue: Double,
    val variables: Map<String, Double>
) {
    /**
     * Get the value of a specific variable
     */
    fun getVariable(name: String): Double {
        return variables[name] ?: 0.0
    }
    
    /**
     * Get utilization for queue i, phase k (1-based indices)
     */
    fun getUtilization(i: Int, k: Int): Double {
        // Try both naming conventions used by different solvers
        val uVar = getVariable("U_${i}_$k")
        if (uVar != 0.0) return uVar
        return getVariable("e_${i}_$k")
    }
    
    /**
     * Get queue length for queue i, phase k (1-based indices)
     */
    fun getQueueLength(i: Int, k: Int): Double {
        return getVariable("Q_${i}_$k")
    }
}