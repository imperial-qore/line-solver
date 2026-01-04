package jline.lib.mom.solver

import org.apache.commons.math3.linear.RealMatrix

/**
 * Container for the linear system matrices used in MOM solver
 * 
 * @property C Main coefficient matrix (block diagonal structure)
 * @property Cg Coupling matrix for previous populations
 * @property D Right-hand side matrix
 * @property Dr Special matrix for recursion
 */
data class LinearSystemMatrices(
    val C: RealMatrix,
    val Cg: RealMatrix,
    val D: RealMatrix,
    val Dr: RealMatrix
)