/**
 * @file Job Population Modification for NetworkStruct
 *
 * Provides functions to directly modify the number of jobs for closed classes
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for population modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.util.matrix.Matrix

/**
 * Sets the number of jobs for a closed class.
 *
 * Updates the `njobs` matrix and recalculates `nclosedjobs`.
 * When `autoRefresh` is true, also calls snRefreshVisits to update visit ratios.
 *
 * @param sn NetworkStruct to modify
 * @param classIdx Class index (0-based)
 * @param nJobs Number of jobs (must be non-negative and finite for closed classes)
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh visit ratios after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetPopulation(
    sn: NetworkStruct,
    classIdx: Int,
    nJobs: Double,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        snValidateClassIndex(sn, classIdx)?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            if (nJobs.isNaN()) {
                errors.add("nJobs is NaN")
            } else if (nJobs < 0) {
                errors.add("nJobs=$nJobs must be non-negative")
            }
            // Check if this is a closed class (should not be Inf)
            val currentNJobs = getPopulationValue(sn.njobs, classIdx)
            if (currentNJobs.isInfinite() && nJobs.isFinite()) {
                errors.add("Class $classIdx is an open class (Inf jobs). Cannot set finite population.")
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetPopulation validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update njobs matrix
    setPopulationValue(snWork.njobs, classIdx, nJobs)

    // Recalculate nclosedjobs
    snWork.nclosedjobs = calculateTotalClosedJobs(snWork.njobs)

    // Auto-refresh visit ratios if requested
    if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
        snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes)
    }

    return snWork
}

/**
 * Sets populations for multiple classes in a single operation.
 *
 * More efficient than calling snSetPopulation multiple times.
 * NaN values in the njobs matrix are skipped (not updated).
 *
 * @param sn NetworkStruct to modify
 * @param njobs Matrix of new job counts (1 x nclasses or nclasses x 1) - NaN values are skipped
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh visit ratios after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetPopulationBatch(
    sn: NetworkStruct,
    njobs: Matrix,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        val totalElements = njobs.numRows * njobs.numCols
        if (totalElements != sn.nclasses) {
            errors.add("njobs matrix total elements ($totalElements) do not match nclasses (${sn.nclasses})")
        }

        if (validation == ValidationLevel.FULL) {
            for (k in 0 until sn.nclasses) {
                val n = getMatrixElement(njobs, k)
                if (!n.isNaN() && n < 0) {
                    errors.add("njobs[$k]=$n must be non-negative")
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetPopulationBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update njobs matrix
    for (k in 0 until sn.nclasses) {
        val n = getMatrixElement(njobs, k)
        if (!n.isNaN()) {
            setPopulationValue(snWork.njobs, k, n)
        }
    }

    // Recalculate nclosedjobs
    snWork.nclosedjobs = calculateTotalClosedJobs(snWork.njobs)

    // Auto-refresh visit ratios if requested
    if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
        snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes)
    }

    return snWork
}

/**
 * Gets population value from njobs matrix.
 * Handles both (1 x K) and (K x 1) layouts.
 */
private fun getPopulationValue(njobs: Matrix?, classIdx: Int): Double {
    if (njobs == null) return Double.POSITIVE_INFINITY

    return if (njobs.numRows == 1) {
        njobs.get(0, classIdx)
    } else {
        njobs.get(classIdx, 0)
    }
}

/**
 * Sets population value in njobs matrix.
 * Handles both (1 x K) and (K x 1) layouts.
 */
private fun setPopulationValue(njobs: Matrix?, classIdx: Int, value: Double) {
    if (njobs == null) return

    if (njobs.numRows == 1) {
        njobs.set(0, classIdx, value)
    } else {
        njobs.set(classIdx, 0, value)
    }
}

/**
 * Gets element from a matrix that may be either (1 x K) or (K x 1).
 */
private fun getMatrixElement(m: Matrix, idx: Int): Double {
    return if (m.numRows == 1) {
        m.get(0, idx)
    } else {
        m.get(idx, 0)
    }
}

/**
 * Calculates total number of closed jobs (sum of finite populations).
 */
private fun calculateTotalClosedJobs(njobs: Matrix?): Int {
    if (njobs == null) return 0

    var total = 0.0
    val totalElements = njobs.numRows * njobs.numCols
    for (i in 0 until totalElements) {
        val value = if (njobs.numRows == 1) {
            njobs.get(0, i)
        } else {
            njobs.get(i, 0)
        }
        if (value.isFinite()) {
            total += value
        }
    }
    return total.toInt()
}

/**
 * Stochastic network SetPopulation algorithms
 */
@Suppress("unused")
class SnsetpopulationAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
