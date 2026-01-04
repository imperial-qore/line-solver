/**
 * @file Service Rate Modification for NetworkStruct
 *
 * Provides functions to directly modify service rates and squared coefficient of
 * variation (SCV) in a NetworkStruct without rebuilding the full Network object.
 * Optimized for fast parameter updates in optimization loops.
 *
 * Mirrors MATLAB implementation patterns for service parameter modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.util.matrix.Matrix

/**
 * Sets the service rate at a specific station and class.
 *
 * Updates the `rates` and `scv` matrices in the NetworkStruct.
 * When `autoRefresh` is true, also updates derived process fields
 * (mu, phi, proc, pie, phases, phasessz, phaseshift, procid).
 *
 * @param sn NetworkStruct to modify
 * @param stationIdx Station index (0-based)
 * @param classIdx Class index (0-based)
 * @param rate New service rate (must be positive)
 * @param scv Squared coefficient of variation (default 1.0 for exponential)
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh process fields after modification
 * @return Modified NetworkStruct (same instance if IN_PLACE, new copy if COPY)
 * @throws SnValidationException if validation fails
 */
fun snSetService(
    sn: NetworkStruct,
    stationIdx: Int,
    classIdx: Int,
    rate: Double,
    scv: Double = 1.0,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        snValidateStationIndex(sn, stationIdx)?.let { errors.add(it) }
        snValidateClassIndex(sn, classIdx)?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            if (rate.isNaN()) {
                errors.add("rate is NaN")
            } else if (rate < 0) {
                errors.add("rate=$rate must be non-negative")
            }
            if (scv.isNaN()) {
                errors.add("scv is NaN")
            } else if (scv < 0) {
                errors.add("scv=$scv must be non-negative")
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetService validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Update rates matrix
    if (snWork.rates != null) {
        snWork.rates.set(stationIdx, classIdx, rate)
    }

    // Update scv matrix
    if (snWork.scv != null) {
        snWork.scv.set(stationIdx, classIdx, scv)
    }

    // Auto-refresh process fields if requested
    if (autoRefresh) {
        snRefreshProcessFields(snWork, stationIdx, classIdx)
    }

    return snWork
}

/**
 * Sets service rates for multiple station-class pairs in a single operation.
 *
 * More efficient than calling snSetService multiple times when updating
 * many parameters. NaN values in the rates matrix are skipped (not updated).
 *
 * @param sn NetworkStruct to modify
 * @param rates Matrix of new rates (nstations x nclasses) - NaN values are skipped
 * @param scvs Matrix of new SCVs (nstations x nclasses) - optional, NaN values skipped
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @param autoRefresh If true, refresh process fields after modification
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetServiceBatch(
    sn: NetworkStruct,
    rates: Matrix,
    scvs: Matrix? = null,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL,
    autoRefresh: Boolean = false
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        if (rates.numRows != sn.nstations || rates.numCols != sn.nclasses) {
            errors.add("rates matrix dimensions (${rates.numRows}x${rates.numCols}) do not match (${sn.nstations}x${sn.nclasses})")
        }

        if (scvs != null) {
            if (scvs.numRows != sn.nstations || scvs.numCols != sn.nclasses) {
                errors.add("scvs matrix dimensions (${scvs.numRows}x${scvs.numCols}) do not match (${sn.nstations}x${sn.nclasses})")
            }
        }

        if (validation == ValidationLevel.FULL) {
            for (i in 0 until rates.numRows) {
                for (j in 0 until rates.numCols) {
                    val r = rates.get(i, j)
                    if (!r.isNaN() && r < 0) {
                        errors.add("rates[$i,$j]=$r must be non-negative")
                    }
                    if (scvs != null) {
                        val s = scvs.get(i, j)
                        if (!s.isNaN() && s < 0) {
                            errors.add("scvs[$i,$j]=$s must be non-negative")
                        }
                    }
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetServiceBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Track which station-class pairs were updated for auto-refresh
    val updatedPairs = mutableListOf<Pair<Int, Int>>()

    // Update rates matrix
    for (i in 0 until rates.numRows) {
        for (j in 0 until rates.numCols) {
            val r = rates.get(i, j)
            if (!r.isNaN()) {
                if (snWork.rates != null) {
                    snWork.rates.set(i, j, r)
                }
                updatedPairs.add(Pair(i, j))
            }
        }
    }

    // Update scv matrix
    if (scvs != null && snWork.scv != null) {
        for (i in 0 until scvs.numRows) {
            for (j in 0 until scvs.numCols) {
                val s = scvs.get(i, j)
                if (!s.isNaN()) {
                    snWork.scv.set(i, j, s)
                }
            }
        }
    }

    // Auto-refresh process fields if requested
    if (autoRefresh) {
        for ((stationIdx, classIdx) in updatedPairs) {
            snRefreshProcessFields(snWork, stationIdx, classIdx)
        }
    }

    return snWork
}

/**
 * Stochastic network SetService algorithms
 */
@Suppress("unused")
class SnsetserviceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
