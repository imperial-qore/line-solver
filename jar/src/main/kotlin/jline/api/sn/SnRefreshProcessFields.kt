/**
 * @file Process Fields Refresh for NetworkStruct
 *
 * Provides functions to refresh process-related fields (mu, phi, proc, pie, phases)
 * in a NetworkStruct based on the current rate and SCV values. This allows updating
 * derived process representations after modifying service rates directly.
 *
 * Mirrors MATLAB implementation patterns for process parameter computation.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.api.mam.map_erlang
import jline.api.mam.map_exponential
import jline.api.mam.map_hyperexp
import jline.api.mam.map_pie
import jline.lang.JobClass
import jline.lang.NetworkStruct
import jline.lang.constant.ProcessType
import jline.lang.nodes.Station
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.ceil
import kotlin.math.max

/**
 * Refreshes process fields for a specific station-class pair based on current rate and SCV.
 *
 * Given the rate and SCV in sn.rates and sn.scv, this function computes and updates:
 * - sn.proc: MAP representation {D0, D1}
 * - sn.procid: Process type (EXP, ERLANG, HYPEREXP)
 * - sn.mu: Phase service rates
 * - sn.phi: Phase completion probabilities
 * - sn.pie: Initial phase distribution
 * - sn.phases, sn.phasessz, sn.phaseshift: Phase count information
 *
 * The process type is determined by SCV:
 * - SCV = 1.0: Exponential (1 phase)
 * - SCV < 1.0: Erlang approximation (k phases where k = ceil(1/scv))
 * - SCV > 1.0: Hyperexponential(2) approximation
 *
 * @param sn NetworkStruct to modify (in-place)
 * @param stationIdx Station index (0-based)
 * @param classIdx Class index (0-based)
 * @return The modified NetworkStruct (same instance)
 */
fun snRefreshProcessFields(
    sn: NetworkStruct,
    stationIdx: Int,
    classIdx: Int
): NetworkStruct {
    // Get rate and SCV
    val rate = sn.rates?.get(stationIdx, classIdx) ?: return sn
    val scv = sn.scv?.get(stationIdx, classIdx) ?: 1.0

    // Skip if rate is invalid
    if (rate.isNaN() || rate <= 0 || !rate.isFinite()) {
        return sn
    }

    // Get station and job class objects
    val station = sn.stations.getOrNull(stationIdx) ?: return sn
    val jobClass = sn.jobclasses.getOrNull(classIdx) ?: return sn

    val mean = 1.0 / rate

    // Determine process type and create MAP based on SCV
    val (map, processType, nPhases) = when {
        scv.isNaN() || kotlin.math.abs(scv - 1.0) < 1e-10 -> {
            // Exponential
            Triple(map_exponential(mean), ProcessType.EXP, 1)
        }
        scv < 1.0 -> {
            // Erlang: k = ceil(1/scv)
            val k = max(1, ceil(1.0 / scv).toInt())
            Triple(map_erlang(mean, k), ProcessType.ERLANG, k)
        }
        else -> {
            // Hyperexponential (scv > 1)
            val hyperMap = map_hyperexp(mean, scv, 0.99)
            if (hyperMap != null) {
                Triple(hyperMap, ProcessType.HYPEREXP, 2)
            } else {
                // Fallback to exponential if hyperexp fails
                Triple(map_exponential(mean), ProcessType.EXP, 1)
            }
        }
    }

    // Update process fields
    updateProcessFields(sn, station, jobClass, stationIdx, classIdx, map, processType, nPhases)

    return sn
}

/**
 * Refreshes process fields for all station-class pairs.
 *
 * Iterates through all stations and classes, refreshing process fields
 * for each pair that has a valid rate defined.
 *
 * @param sn NetworkStruct to modify (in-place)
 * @return The modified NetworkStruct (same instance)
 */
fun snRefreshAllProcessFields(sn: NetworkStruct): NetworkStruct {
    for (i in 0 until sn.nstations) {
        for (j in 0 until sn.nclasses) {
            snRefreshProcessFields(sn, i, j)
        }
    }
    return sn
}

/**
 * Updates all process-related fields in NetworkStruct for a given MAP.
 *
 * @param sn NetworkStruct to update
 * @param station Station object
 * @param jobClass JobClass object
 * @param stationIdx Station index
 * @param classIdx Class index
 * @param map MAP representation {D0, D1}
 * @param processType Process type enum
 * @param nPhases Number of phases
 */
private fun updateProcessFields(
    sn: NetworkStruct,
    station: Station,
    jobClass: JobClass,
    stationIdx: Int,
    classIdx: Int,
    map: MatrixCell,
    processType: ProcessType,
    nPhases: Int
) {
    val d0 = map[0]
    val d1 = map[1]

    // Initialize maps if null
    if (sn.proc == null) {
        sn.proc = HashMap()
    }
    if (sn.procid == null) {
        sn.procid = HashMap()
    }
    if (sn.mu == null) {
        sn.mu = HashMap()
    }
    if (sn.phi == null) {
        sn.phi = HashMap()
    }
    if (sn.pie == null) {
        sn.pie = HashMap()
    }

    // Initialize inner maps if null
    if (sn.proc[station] == null) {
        sn.proc[station] = HashMap()
    }
    if (sn.procid[station] == null) {
        sn.procid[station] = HashMap()
    }
    if (sn.mu[station] == null) {
        sn.mu[station] = HashMap()
    }
    if (sn.phi[station] == null) {
        sn.phi[station] = HashMap()
    }
    if (sn.pie[station] == null) {
        sn.pie[station] = HashMap()
    }

    // Update process representation
    sn.proc[station]?.set(jobClass, map)
    sn.procid[station]?.set(jobClass, processType)

    // Update phases
    sn.phases?.set(stationIdx, classIdx, nPhases.toDouble())

    // Update phasessz
    sn.phasessz?.set(stationIdx, classIdx, max(nPhases, 1).toDouble())

    // Recompute phaseshift for this station (cumulative sum across classes)
    if (sn.phaseshift != null) {
        var cumSum = 0.0
        sn.phaseshift.set(stationIdx, 0, 0.0)
        for (c in 0 until sn.nclasses) {
            cumSum += sn.phasessz?.get(stationIdx, c) ?: 1.0
            if (c + 1 < sn.phaseshift.numCols) {
                sn.phaseshift.set(stationIdx, c + 1, cumSum)
            }
        }
    }

    // Update mu (rates from -diag(D0))
    val muMatrix = Matrix(nPhases, 1)
    for (i in 0 until nPhases) {
        muMatrix.set(i, 0, -d0.get(i, i))
    }
    sn.mu[station]?.set(jobClass, muMatrix)

    // Update phi (completion probabilities: sum(D1,2) / -diag(D0))
    val phiMatrix = Matrix(nPhases, 1)
    for (i in 0 until nPhases) {
        var d1RowSum = 0.0
        for (j in 0 until d1.numCols) {
            d1RowSum += d1.get(i, j)
        }
        val d0Diag = -d0.get(i, i)
        phiMatrix.set(i, 0, if (d0Diag != 0.0) d1RowSum / d0Diag else 0.0)
    }
    sn.phi[station]?.set(jobClass, phiMatrix)

    // Update pie (initial phase distribution)
    val pieMatrix = map_pie(map)
    sn.pie[station]?.set(jobClass, pieMatrix)
}

/**
 * Stochastic network RefreshProcessFields algorithms
 */
@Suppress("unused")
class SnrefreshprocessfieldsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
