/**
 * @file Class Priority Modification for NetworkStruct
 *
 * Provides functions to directly modify class priorities
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for priority modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.util.matrix.Matrix

/**
 * Sets the priority for a class.
 *
 * Updates the `classprio` matrix in the NetworkStruct.
 * Lower priority values typically indicate higher priority (depends on scheduler).
 *
 * @param sn NetworkStruct to modify
 * @param classIdx Class index (0-based)
 * @param priority Priority value
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetPriority(
    sn: NetworkStruct,
    classIdx: Int,
    priority: Double,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        snValidateClassIndex(sn, classIdx)?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            if (priority.isNaN()) {
                errors.add("priority is NaN")
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetPriority validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update classprio matrix
    setPriorityValue(snWork.classprio, classIdx, priority)

    return snWork
}

/**
 * Sets priorities for multiple classes in a single operation.
 *
 * More efficient than calling snSetPriority multiple times.
 * NaN values in the priorities matrix are skipped (not updated).
 *
 * @param sn NetworkStruct to modify
 * @param priorities Matrix of new priorities (1 x nclasses or nclasses x 1) - NaN values are skipped
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetPriorityBatch(
    sn: NetworkStruct,
    priorities: Matrix,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        val totalElements = priorities.numRows * priorities.numCols
        if (totalElements != sn.nclasses) {
            errors.add("priorities matrix total elements ($totalElements) do not match nclasses (${sn.nclasses})")
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetPriorityBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update classprio matrix
    for (k in 0 until sn.nclasses) {
        val p = getMatrixElement(priorities, k)
        if (!p.isNaN()) {
            setPriorityValue(snWork.classprio, k, p)
        }
    }

    return snWork
}

/**
 * Sets priority value in classprio matrix.
 * Handles both (1 x K) and (K x 1) layouts.
 */
private fun setPriorityValue(classprio: Matrix?, classIdx: Int, value: Double) {
    if (classprio == null) return

    if (classprio.numRows == 1) {
        classprio.set(0, classIdx, value)
    } else {
        classprio.set(classIdx, 0, value)
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
 * Stochastic network SetPriority algorithms
 */
@Suppress("unused")
class SnsetpriorityAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
