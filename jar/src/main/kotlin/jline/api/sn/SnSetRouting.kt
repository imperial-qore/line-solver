/**
 * @file Routing Matrix Modification for NetworkStruct
 *
 * Provides functions to directly modify routing matrices (rt, rtnodes)
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for routing modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.util.matrix.Matrix

/**
 * Sets a routing probability between two stateful node-class pairs.
 *
 * Updates the `rt` matrix in the NetworkStruct. The routing matrix is indexed
 * as (fromStateful * nclasses + fromClass, toStateful * nclasses + toClass).
 *
 * @param sn NetworkStruct to modify
 * @param fromStateful Source stateful node index (0-based)
 * @param fromClass Source class index (0-based)
 * @param toStateful Destination stateful node index (0-based)
 * @param toClass Destination class index (0-based)
 * @param probability Routing probability [0, 1]
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh visit ratios after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetRoutingProb(
    sn: NetworkStruct,
    fromStateful: Int,
    fromClass: Int,
    toStateful: Int,
    toClass: Int,
    probability: Double,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        if (fromStateful < 0 || fromStateful >= sn.nstateful) {
            errors.add("fromStateful=$fromStateful is out of bounds [0, ${sn.nstateful - 1}]")
        }
        if (toStateful < 0 || toStateful >= sn.nstateful) {
            errors.add("toStateful=$toStateful is out of bounds [0, ${sn.nstateful - 1}]")
        }
        snValidateClassIndex(sn, fromClass, "fromClass")?.let { errors.add(it) }
        snValidateClassIndex(sn, toClass, "toClass")?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            if (probability.isNaN()) {
                errors.add("probability is NaN")
            } else if (probability < 0 || probability > 1) {
                errors.add("probability=$probability must be in [0, 1]")
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetRoutingProb validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Calculate indices in rt matrix
    val fromIdx = fromStateful * sn.nclasses + fromClass
    val toIdx = toStateful * sn.nclasses + toClass

    // Update rt matrix
    if (snWork.rt != null) {
        snWork.rt.set(fromIdx, toIdx, probability)
    }

    // Auto-refresh visit ratios if requested
    if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
        snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes)
    }

    return snWork
}

/**
 * Sets the entire routing matrix for stateful nodes.
 *
 * Replaces the `rt` matrix in the NetworkStruct.
 * The matrix should be of size (nstateful * nclasses) x (nstateful * nclasses).
 *
 * @param sn NetworkStruct to modify
 * @param rt New routing matrix
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh visit ratios after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetRoutingMatrix(
    sn: NetworkStruct,
    rt: Matrix,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    val expectedSize = sn.nstateful * sn.nclasses

    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        if (rt.numRows != expectedSize || rt.numCols != expectedSize) {
            errors.add("rt matrix dimensions (${rt.numRows}x${rt.numCols}) do not match expected ($expectedSize x $expectedSize)")
        }

        if (validation == ValidationLevel.FULL) {
            validateStochasticMatrix(rt, "rt", errors)
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetRoutingMatrix validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Replace rt matrix
    snWork.rt = rt

    // Auto-refresh visit ratios if requested
    if (autoRefresh && snWork.chains != null && snWork.rtnodes != null) {
        snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes)
    }

    return snWork
}

/**
 * Sets the entire routing matrix for all nodes.
 *
 * Replaces the `rtnodes` matrix in the NetworkStruct.
 * The matrix should be of size (nnodes * nclasses) x (nnodes * nclasses).
 *
 * @param sn NetworkStruct to modify
 * @param rtnodes New node routing matrix
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh visit ratios after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetRoutingNodesMatrix(
    sn: NetworkStruct,
    rtnodes: Matrix,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    val expectedSize = sn.nnodes * sn.nclasses

    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        if (rtnodes.numRows != expectedSize || rtnodes.numCols != expectedSize) {
            errors.add("rtnodes matrix dimensions (${rtnodes.numRows}x${rtnodes.numCols}) do not match expected ($expectedSize x $expectedSize)")
        }

        if (validation == ValidationLevel.FULL) {
            validateStochasticMatrix(rtnodes, "rtnodes", errors)
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetRoutingNodesMatrix validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Replace rtnodes matrix
    snWork.rtnodes = rtnodes

    // Auto-refresh visit ratios if requested
    if (autoRefresh && snWork.chains != null && snWork.rt != null) {
        snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes)
    }

    return snWork
}

/**
 * Validates that a matrix is stochastic (row sums approximately equal 1).
 */
private fun validateStochasticMatrix(m: Matrix, name: String, errors: MutableList<String>) {
    val tolerance = 1e-6

    for (i in 0 until m.numRows) {
        var rowSum = 0.0
        var hasNonZero = false
        var hasNegative = false

        for (j in 0 until m.numCols) {
            val value = m.get(i, j)
            if (value.isNaN()) {
                errors.add("$name[$i,$j] is NaN")
                continue
            }
            if (value < 0) {
                hasNegative = true
            }
            if (value > 0) {
                hasNonZero = true
            }
            rowSum += value
        }

        if (hasNegative) {
            errors.add("$name row $i contains negative values")
        }

        // Only check row sum if row has non-zero entries (active row)
        if (hasNonZero && kotlin.math.abs(rowSum - 1.0) > tolerance) {
            errors.add("$name row $i sum = $rowSum (expected 1.0)")
        }
    }
}

/**
 * Stochastic network SetRouting algorithms
 */
@Suppress("unused")
class SnsetroutingAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
