/**
 * @file Arrival Rate Modification for NetworkStruct
 *
 * Provides functions to directly modify arrival rates at the Source station
 * in a NetworkStruct without rebuilding the full Network object.
 * Delegates to snSetService with the Source station index.
 *
 * Mirrors MATLAB implementation patterns for arrival parameter modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

/**
 * Sets the arrival rate for a class at the Source station.
 *
 * This is a convenience method that finds the Source station and delegates
 * to snSetService. Arrival rates are stored in the same rates matrix as
 * service rates, at the Source station's row.
 *
 * @param sn NetworkStruct to modify
 * @param classIdx Class index (0-based)
 * @param rate New arrival rate (lambda, must be positive)
 * @param scv Squared coefficient of variation (default 1.0 for Poisson arrivals)
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh process fields after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails or no Source station exists
 */
fun snSetArrival(
    sn: NetworkStruct,
    classIdx: Int,
    rate: Double,
    scv: Double = 1.0,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Find Source station index
    val sourceStationIdx = findSourceStationIndex(sn)

    if (sourceStationIdx < 0) {
        if (validation != ValidationLevel.NONE) {
            throw SnValidationException("snSetArrival: No Source station found in network")
        }
        return sn
    }

    // Validate that this is an open class (has Inf jobs)
    if (validation == ValidationLevel.FULL) {
        val njob = getPopulationValue(sn, classIdx)
        if (njob.isFinite()) {
            throw SnValidationException(
                "snSetArrival: Class $classIdx is a closed class (njobs=$njob). " +
                "Arrival rates can only be set for open classes."
            )
        }
    }

    // Delegate to snSetService
    return snSetService(sn, sourceStationIdx, classIdx, rate, scv, mode, validation, autoRefresh)
}

/**
 * Sets arrival rates for multiple classes in a single operation.
 *
 * More efficient than calling snSetArrival multiple times when updating
 * many classes. NaN values in the rates matrix are skipped.
 *
 * @param sn NetworkStruct to modify
 * @param rates Matrix of new arrival rates (1 x nclasses or nclasses x 1) - NaN values are skipped
 * @param scvs Matrix of new SCVs (same size as rates) - optional, NaN values skipped
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh process fields after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetArrivalBatch(
    sn: NetworkStruct,
    rates: Matrix,
    scvs: Matrix? = null,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Find Source station index
    val sourceStationIdx = findSourceStationIndex(sn)

    if (sourceStationIdx < 0) {
        if (validation != ValidationLevel.NONE) {
            throw SnValidationException("snSetArrivalBatch: No Source station found in network")
        }
        return sn
    }

    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()
        val totalElements = rates.numRows * rates.numCols

        if (totalElements != sn.nclasses) {
            errors.add("rates matrix total elements ($totalElements) do not match nclasses (${sn.nclasses})")
        }

        if (scvs != null) {
            val scvElements = scvs.numRows * scvs.numCols
            if (scvElements != sn.nclasses) {
                errors.add("scvs matrix total elements ($scvElements) do not match nclasses (${sn.nclasses})")
            }
        }

        if (validation == ValidationLevel.FULL) {
            for (k in 0 until sn.nclasses) {
                val r = getMatrixElement(rates, k)
                if (!r.isNaN() && r < 0) {
                    errors.add("rates[$k]=$r must be non-negative")
                }
                if (scvs != null) {
                    val s = getMatrixElement(scvs, k)
                    if (!s.isNaN() && s < 0) {
                        errors.add("scvs[$k]=$s must be non-negative")
                    }
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetArrivalBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update rates at source station
    for (k in 0 until sn.nclasses) {
        val r = getMatrixElement(rates, k)
        if (!r.isNaN()) {
            snWork.rates?.set(sourceStationIdx, k, r)
        }
        if (scvs != null) {
            val s = getMatrixElement(scvs, k)
            if (!s.isNaN()) {
                snWork.scv?.set(sourceStationIdx, k, s)
            }
        }
    }

    // Auto-refresh if requested
    if (autoRefresh) {
        for (k in 0 until sn.nclasses) {
            val r = getMatrixElement(rates, k)
            if (!r.isNaN()) {
                snRefreshProcessFields(snWork, sourceStationIdx, k)
            }
        }
    }

    return snWork
}

/**
 * Finds the index of the Source station in the network.
 *
 * @param sn NetworkStruct to search
 * @return Source station index (0-based), or -1 if not found
 */
private fun findSourceStationIndex(sn: NetworkStruct): Int {
    // First check nodetype for Source
    for (i in sn.nodetype.indices) {
        if (sn.nodetype[i] == NodeType.Source) {
            // Convert node index to station index
            val stationIdx = sn.nodeToStation?.get(0, i)?.toInt() ?: -1
            if (stationIdx >= 0) {
                return stationIdx
            }
        }
    }

    // Fallback: Source is typically the first station (index 0)
    if (sn.nstations > 0 && sn.nodetype.isNotEmpty()) {
        val firstNodeType = sn.nodetype[0]
        if (firstNodeType == NodeType.Source) {
            return 0
        }
    }

    return -1
}

/**
 * Gets population value for a class from njobs matrix.
 * Handles both (1 x K) and (K x 1) layouts.
 */
private fun getPopulationValue(sn: NetworkStruct, classIdx: Int): Double {
    if (sn.njobs == null) return Double.POSITIVE_INFINITY

    return if (sn.njobs.numRows == 1) {
        sn.njobs.get(0, classIdx)
    } else {
        sn.njobs.get(classIdx, 0)
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
 * Stochastic network SetArrival algorithms
 */
@Suppress("unused")
class SnsetarrivalAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
