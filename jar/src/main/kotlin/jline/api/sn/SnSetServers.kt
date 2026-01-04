/**
 * @file Server Count Modification for NetworkStruct
 *
 * Provides functions to directly modify the number of servers at stations
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for server multiplicity modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

/**
 * Sets the number of servers at a station.
 *
 * Updates the `nservers` matrix in the NetworkStruct.
 *
 * @param sn NetworkStruct to modify
 * @param stationIdx Station index (0-based)
 * @param nServers Number of servers (positive, or Double.POSITIVE_INFINITY for infinite server)
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetServers(
    sn: NetworkStruct,
    stationIdx: Int,
    nServers: Double,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        snValidateStationIndex(sn, stationIdx)?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            if (nServers.isNaN()) {
                errors.add("nServers is NaN")
            } else if (nServers <= 0 && !nServers.isInfinite()) {
                // Check node type - some nodes may allow 0 servers
                val nodeIdx = sn.stationToNode?.get(stationIdx, 0)?.toInt() ?: stationIdx
                val nodeType = sn.nodetype.getOrNull(nodeIdx)
                if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                    errors.add("nServers=$nServers must be positive for station type $nodeType")
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetServers validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update nservers matrix
    if (snWork.nservers != null) {
        snWork.nservers.set(stationIdx, 0, nServers)
    }

    return snWork
}

/**
 * Sets the number of servers for multiple stations in a single operation.
 *
 * More efficient than calling snSetServers multiple times.
 * NaN values in the nServers matrix are skipped (not updated).
 *
 * @param sn NetworkStruct to modify
 * @param nServers Matrix of new server counts (nstations x 1) - NaN values are skipped
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetServersBatch(
    sn: NetworkStruct,
    nServers: Matrix,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        val totalElements = nServers.numRows * nServers.numCols
        if (totalElements != sn.nstations) {
            errors.add("nServers matrix total elements ($totalElements) do not match nstations (${sn.nstations})")
        }

        if (validation == ValidationLevel.FULL) {
            for (i in 0 until sn.nstations) {
                val n = getMatrixElement(nServers, i)
                if (!n.isNaN() && n <= 0 && !n.isInfinite()) {
                    val nodeIdx = sn.stationToNode?.get(i, 0)?.toInt() ?: i
                    val nodeType = sn.nodetype.getOrNull(nodeIdx)
                    if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                        errors.add("nServers[$i]=$n must be positive")
                    }
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetServersBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update nservers matrix
    for (i in 0 until sn.nstations) {
        val n = getMatrixElement(nServers, i)
        if (!n.isNaN()) {
            snWork.nservers?.set(i, 0, n)
        }
    }

    return snWork
}

/**
 * Gets element from a matrix that may be either (M x 1) or (1 x M).
 */
private fun getMatrixElement(m: Matrix, idx: Int): Double {
    return if (m.numCols == 1) {
        m.get(idx, 0)
    } else {
        m.get(0, idx)
    }
}

/**
 * Stochastic network SetServers algorithms
 */
@Suppress("unused")
class SnsetserversAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
