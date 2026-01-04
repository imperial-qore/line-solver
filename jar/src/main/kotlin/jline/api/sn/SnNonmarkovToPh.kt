/**
 * @file Non-Markovian to Phase-Type Distribution Converter
 *
 * Converts non-Markovian distributions in a queueing network to phase-type (PH)
 * representations using Bernstein polynomial approximation. This enables analytical
 * solvers that require Markovian distributions to handle a broader range of
 * distribution types.
 *
 * Mirrors MATLAB implementation: matlab/src/api/sn/sn_nonmarkov_toph.m
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.api.mam.aph_bernstein
import jline.api.mam.map_erlang
import jline.api.mam.map_pie
import jline.api.mam.map_scale
import jline.io.line_warning
import jline.lang.Event
import jline.lang.JobClass
import jline.lang.NetworkStruct
import jline.lang.Sync
import jline.lang.constant.EventType
import jline.lang.constant.ProcessType
import jline.lang.nodes.Station
import jline.lang.nodes.StatefulNode
import jline.lang.processes.*
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.pow

/**
 * Set of process types that do not require Bernstein conversion.
 * Includes Markovian types (PH, MAP, EXP, etc.) and matrix-exponential
 * types (ME, RAP) which are already in a form usable by matrix-analytic methods.
 */
private val SKIP_BERNSTEIN_CONVERSION = setOf(
    ProcessType.EXP,
    ProcessType.ERLANG,
    ProcessType.HYPEREXP,
    ProcessType.PH,
    ProcessType.APH,
    ProcessType.MAP,
    ProcessType.MMAP,    // Marked MAP - already Markovian
    ProcessType.BMAP,    // Batch MAP - already Markovian
    ProcessType.ME,      // Matrix Exponential - skip conversion, use ME-specific solvers
    ProcessType.RAP,     // Rational Arrival Process - skip conversion, use RAP-specific solvers
    ProcessType.COXIAN,
    ProcessType.COX2,
    ProcessType.MMPP2,
    ProcessType.IMMEDIATE,
    ProcessType.DISABLED
)

/**
 * Converts non-Markovian distributions to PH using Bernstein approximation.
 *
 * This function scans all service and arrival processes in the network
 * structure and converts non-Markovian distributions to Markovian Arrival
 * Processes (MAPs) using the specified approximation method.
 *
 * @param snInput Network structure from getStruct()
 * @param options Solver options with nonmkv and nonmkvorder config
 * @return Network structure with converted processes
 */
fun snNonmarkovToPh(snInput: NetworkStruct, options: SolverOptions): NetworkStruct {
    // Get non-Markovian conversion method from options (default 'bernstein')
    val nonmkvMethod = options.config?.nonmkv ?: "bernstein"

    // If method is 'none', return without any conversion
    if (nonmkvMethod.equals("none", ignoreCase = true)) {
        return snInput
    }

    // Get number of phases from options (default 20)
    val nPhases = options.config?.nonmkvorder ?: 20

    // Iterate through all stations and classes
    for (ist in 0 until snInput.nstations) {
        val station = snInput.stations[ist]

        for (r in 0 until snInput.nclasses) {
            val jobClass = snInput.jobclasses[r]

            // Get process type for this station/class
            val procType = snInput.procid[station]?.get(jobClass) ?: continue

            // Skip types that don't need Bernstein conversion
            if (procType in SKIP_BERNSTEIN_CONVERSION) {
                continue
            }

            // Get target mean from rates
            val rate = snInput.rates.get(ist, r)
            if (rate <= 0 || !rate.isFinite()) {
                continue
            }
            val targetMean = 1.0 / rate

            // Issue warning
            val distName = ProcessType.toText(procType)
            line_warning(
                "snNonmarkovToPh",
                "Distribution $distName at station $ist class $r is non-Markovian and will be converted to PH ($nPhases phases)."
            )

            // Get original process representation
            val origProc = snInput.proc[station]?.get(jobClass)
            if (origProc == null) {
                continue
            }

            // Convert based on distribution type
            val map: MatrixCell = when (procType) {
                ProcessType.GAMMA -> {
                    val shape = origProc[0]?.toDouble() ?: continue
                    val scale = origProc[1]?.toDouble() ?: continue
                    val pdfFunc: (Double) -> Double = { x ->
                        if (x <= 0) 0.0 else Gamma(shape, scale).evalPDF(x)
                    }
                    val (d0, d1) = aph_bernstein(pdfFunc, nPhases)
                    map_scale(pairToMatrixCell(d0, d1), targetMean)
                }

                ProcessType.WEIBULL -> {
                    val shapeParam = origProc[0]?.toDouble() ?: continue  // r (shape)
                    val scaleParam = origProc[1]?.toDouble() ?: continue  // alpha (scale)
                    val pdfFunc: (Double) -> Double = { x ->
                        if (x <= 0) 0.0 else Weibull(shapeParam, scaleParam).evalPDF(x)
                    }
                    val (d0, d1) = aph_bernstein(pdfFunc, nPhases)
                    map_scale(pairToMatrixCell(d0, d1), targetMean)
                }

                ProcessType.LOGNORMAL -> {
                    val mu = origProc[0]?.toDouble() ?: continue
                    val sigma = origProc[1]?.toDouble() ?: continue
                    val pdfFunc: (Double) -> Double = { x ->
                        if (x <= 0) 0.0 else Lognormal(mu, sigma).evalPDF(x)
                    }
                    val (d0, d1) = aph_bernstein(pdfFunc, nPhases)
                    map_scale(pairToMatrixCell(d0, d1), targetMean)
                }

                ProcessType.PARETO -> {
                    val shapeParam = origProc[0]?.toDouble() ?: continue  // alpha (shape)
                    val scaleParam = origProc[1]?.toDouble() ?: continue  // k (scale/minimum)
                    val pdfFunc: (Double) -> Double = { x ->
                        if (x < scaleParam) 0.0
                        else shapeParam * scaleParam.pow(shapeParam) / x.pow(shapeParam + 1)
                    }
                    val (d0, d1) = aph_bernstein(pdfFunc, nPhases)
                    map_scale(pairToMatrixCell(d0, d1), targetMean)
                }

                ProcessType.UNIFORM -> {
                    val minVal = origProc[0]?.toDouble() ?: continue
                    val maxVal = origProc[1]?.toDouble() ?: continue
                    val pdfFunc: (Double) -> Double = { x ->
                        if (x >= minVal && x <= maxVal) 1.0 / (maxVal - minVal) else 0.0
                    }
                    val (d0, d1) = aph_bernstein(pdfFunc, nPhases)
                    map_scale(pairToMatrixCell(d0, d1), targetMean)
                }

                ProcessType.DET -> {
                    // Deterministic: use Erlang approximation directly
                    map_erlang(targetMean, nPhases)
                }

                else -> {
                    // Generic fallback: Erlang approximation
                    map_erlang(targetMean, nPhases)
                }
            }

            // Update the network structure for the converted MAP
            val actualPhases = map[0].numRows
            updateSnForMAP(snInput, station, jobClass, ist, r, map, actualPhases)
        }
    }

    return snInput
}

/**
 * Converts a Pair of matrices (D0, D1) to a MatrixCell.
 */
private fun pairToMatrixCell(d0: Matrix, d1: Matrix): MatrixCell {
    val cell = MatrixCell()
    cell[0] = d0
    cell[1] = d1
    return cell
}

/**
 * Updates all network structure fields for a converted MAP.
 *
 * @param sn Network structure
 * @param station Station object
 * @param jobClass JobClass object
 * @param ist Station index
 * @param r Class index
 * @param map MAP representation {D0, D1}
 * @param nPhases Number of phases in the MAP
 */
private fun updateSnForMAP(
    sn: NetworkStruct,
    station: Station,
    jobClass: JobClass,
    ist: Int,
    r: Int,
    map: MatrixCell,
    nPhases: Int
) {
    val d0 = map[0]
    val d1 = map[1]

    // Save old phasessz before updating (needed for state expansion)
    val oldPhases = sn.phasessz.get(ist, r).toInt()

    // Update process representation
    sn.proc[station]?.set(jobClass, map)
    sn.procid[station]?.set(jobClass, ProcessType.MAP)

    // Update phases
    sn.phases.set(ist, r, nPhases.toDouble())

    // Update phasessz (max of nPhases and 1)
    sn.phasessz.set(ist, r, maxOf(nPhases, 1).toDouble())

    // Recompute phaseshift for this station (cumulative sum across classes)
    var cumSum = 0.0
    sn.phaseshift.set(ist, 0, 0.0)
    for (c in 0 until sn.nclasses) {
        cumSum += sn.phasessz.get(ist, c)
        if (c + 1 < sn.phaseshift.numCols) {
            sn.phaseshift.set(ist, c + 1, cumSum)
        }
    }

    // Update mu (rates from -diag(D0))
    val muMatrix = Matrix(nPhases, 1)
    for (i in 0 until nPhases) {
        muMatrix.set(i, 0, -d0.get(i, i))
    }
    if (sn.mu[station] == null) {
        sn.mu[station] = HashMap()
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
    if (sn.phi[station] == null) {
        sn.phi[station] = HashMap()
    }
    sn.phi[station]?.set(jobClass, phiMatrix)

    // Update pie (initial phase distribution)
    val pieMatrix = map_pie(map)
    if (sn.pie[station] == null) {
        sn.pie[station] = HashMap()
    }
    sn.pie[station]?.set(jobClass, pieMatrix)

    // Update nvars for MAP local variable
    val ind = sn.stationToNode.get(ist).toInt()
    sn.nvars.set(ind, r, sn.nvars.get(ind, r) + 1)

    // Expand state vector for the new phases
    expandStateForMAP(sn, ind, r, oldPhases, nPhases)

    // Add PHASE sync event if phases > 1
    if (nPhases > 1) {
        addPhaseSyncIfNeeded(sn, ind, r)
    }
}

/**
 * Expands the state vector to accommodate additional phases from MAP conversion.
 *
 * When converting a non-Markovian distribution to MAP, the state vector needs to be expanded:
 * 1. The server portion (space_srv) needs additional columns for extra phases
 * 2. The local variable portion (space_var) needs a column for MAP phase memory
 *
 * State format: [space_buf | space_srv | space_var]
 */
private fun expandStateForMAP(sn: NetworkStruct, ind: Int, r: Int, oldPhases: Int, newPhases: Int) {
    val isf = sn.nodeToStateful.get(ind).toInt()
    if (isf < 0 || sn.state.isEmpty()) return

    // Find the stateful node
    var statefulNode: StatefulNode? = null
    for (entry in sn.state.entries) {
        if (entry.key?.getStatefulIndex() == isf) {
            statefulNode = entry.key
            break
        }
    }
    if (statefulNode == null) return

    val stateMatrix = sn.state[statefulNode] ?: return
    if (stateMatrix.isEmpty) return

    val ist = sn.nodeToStation.get(ind).toInt()
    val nRows = stateMatrix.numRows

    // Calculate state vector structure
    val V = sn.nvars.getRow(ind).elementSum().toInt()  // Already incremented
    val V_old = V - 1  // Before we added the new local var

    // Get phases array (already updated)
    val sumK_new = sn.phasessz.getRow(ist).elementSum().toInt()
    var sumK_old = 0
    for (c in 0 until sn.nclasses) {
        sumK_old += if (c == r) oldPhases else sn.phasessz.get(ist, c).toInt()
    }

    // Calculate buffer size
    val currentCols = stateMatrix.numCols
    var bufSize = currentCols - sumK_old - V_old
    if (bufSize < 0) bufSize = 0

    // Extract state portions
    val spaceBuf = if (bufSize > 0) {
        Matrix.extract(stateMatrix, 0, nRows, 0, bufSize)
    } else {
        Matrix(nRows, 0)
    }

    val spaceSrv = if (sumK_old > 0) {
        Matrix.extract(stateMatrix, 0, nRows, bufSize, bufSize + sumK_old)
    } else {
        Matrix(nRows, 0)
    }

    val spaceVar = if (V_old > 0) {
        Matrix.extract(stateMatrix, 0, nRows, bufSize + sumK_old, currentCols)
    } else {
        Matrix(nRows, 0)
    }

    // Expand space_srv: insert zeros for new phases
    val phasesToAdd = newPhases - oldPhases
    val spaceSrvNew = if (phasesToAdd > 0 && spaceSrv.numCols > 0) {
        // Calculate insert position (after class r's old phases)
        var insertPos = 0
        for (c in 0 until r) {
            insertPos += if (c == r) oldPhases else sn.phasessz.get(ist, c).toInt()
        }
        insertPos += oldPhases

        // Create expanded matrix
        val expanded = Matrix(nRows, spaceSrv.numCols + phasesToAdd)
        expanded.zero()
        // Copy before insert
        for (row in 0 until nRows) {
            for (col in 0 until insertPos) {
                if (col < spaceSrv.numCols) {
                    expanded.set(row, col, spaceSrv.get(row, col))
                }
            }
            // Zeros are already there for new phases
            // Copy after insert
            for (col in insertPos until spaceSrv.numCols) {
                expanded.set(row, col + phasesToAdd, spaceSrv.get(row, col))
            }
        }
        expanded
    } else {
        spaceSrv
    }

    // Expand space_var: add 1 column for the MAP local variable
    val spaceVarNew = Matrix(nRows, spaceVar.numCols + 1)
    for (row in 0 until nRows) {
        for (col in 0 until spaceVar.numCols) {
            spaceVarNew.set(row, col, spaceVar.get(row, col))
        }
        spaceVarNew.set(row, spaceVar.numCols, 1.0)  // Initialize to phase 1
    }

    // Reconstruct state
    val newState = Matrix.concatColumns(spaceBuf, spaceSrvNew, null)
    val finalState = Matrix.concatColumns(newState, spaceVarNew, null)
    sn.state[statefulNode] = finalState
}

/**
 * Adds a PHASE sync event for the converted MAP if not already present.
 *
 * When a non-Markovian distribution is converted to a MAP with multiple phases,
 * we need to add a PHASE sync event so that phase transitions can occur during simulation.
 */
private fun addPhaseSyncIfNeeded(sn: NetworkStruct, ind: Int, r: Int) {
    val local = sn.nnodes  // Local event node index

    // Check if PHASE sync already exists for this node/class
    for (syncEntry in sn.sync.values) {
        if (syncEntry?.active?.isNotEmpty() == true) {
            val activeEvent = syncEntry.active[0]
            if (activeEvent?.getEvent() == EventType.PHASE &&
                activeEvent.getNode() == ind &&
                activeEvent.jobClass == r) {
                return  // Already exists
            }
        }
    }

    // Add new PHASE sync
    val newSync = Sync()
    val activeEvent = Event(EventType.PHASE, ind, r)
    val passiveEvent = Event(EventType.LOCAL, local, r)
    passiveEvent.prob = 1.0
    newSync.active[0] = activeEvent
    newSync.passive[0] = passiveEvent

    // Add to sync map
    val nextKey = sn.sync.size
    sn.sync[nextKey] = newSync
}

@Suppress("unused")
class SnNonmarkovTophAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
