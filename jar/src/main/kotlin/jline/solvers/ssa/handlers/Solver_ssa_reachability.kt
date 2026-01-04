package jline.solvers.ssa.handlers

import jline.io.Ret
import jline.lang.NetworkStruct
import jline.lang.constant.EventType
import jline.lang.constant.NodeType
import jline.lang.state.State
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import jline.util.PopulationLattice.pprod
import java.util.*

/**
 * Kotlin migration of solver_ssa_reachability.m
 * Computes the reachable state space for SSA analysis
 * 
 * @param sn Network structure
 * @param options Solver options
 * @return Triple containing SSq (state space matrix), SSh (state hash indices), and updated sn
 */
fun solver_ssa_reachability(sn: NetworkStruct, options: SolverOptions): Triple<Matrix?, Matrix, NetworkStruct> {
    val nstateful = sn.nstateful
    val R = sn.nclasses
    val N = sn.njobs.transpose()
    val sync = sn.sync
    val csmask = sn.csmask
    val stack = mutableListOf<List<Matrix>>()
    val stackIndex = mutableListOf<Int>()
    
    // Initialize with starting state  
    val initialStateCell = mutableListOf<Matrix>()
    for (i in 0 until nstateful) {
        initialStateCell.add(sn.state[sn.stateful[i]]!!)
    }
    stack.add(initialStateCell)
    
    var SSq: Matrix? = null
    val A = sync.size
    val isSimulation = false
    val local = sn.nnodes + 1
    
    // Pre-compute sync node and class information
    val nodeA = IntArray(A) { act -> sync[act]!!.active[0]!!.node }
    val nodeP = IntArray(A) { act -> sync[act]!!.passive[0]!!.node }
    val classA = IntArray(A) { act -> sync[act]!!.active[0]!!.jobClass }
    val classP = IntArray(A) { act -> sync[act]!!.passive[0]!!.jobClass }
    val eventA = Array<EventType>(A) { act -> sync[act]!!.active[0]!!.event }
    val eventP = Array<EventType>(A) { act -> sync[act]!!.passive[0]!!.event }
    
    val space = mutableListOf<Matrix>()
    for (i in 0 until nstateful) {
        space.add(sn.state[sn.stateful[i]]!!.copy())
    }
    
    val SSh = Matrix(1, nstateful)
    SSh.fill(1.0) // Initial state hash (1-based indexing)
    
    var ih = 1
    stackIndex.add(1)
    val maxstatesz = IntArray(nstateful)
    
    while (stack.isNotEmpty()) {
        if (stack.isEmpty()) {
            // Build SSq from space and SSh
            val totalCols = space.sumOf { it.getNumCols() }
            SSq = Matrix(SSh.getNumRows(), totalCols)
            for (i in 0 until SSh.getNumRows()) {
                var colctr = 0
                for (j in 0 until nstateful) {
                    val stateIdx = SSh.get(i, j).toInt() - 1 // Convert to 0-based
                    val stateRow = space[j].getRow(stateIdx)
                    for (k in 0 until stateRow.getNumCols()) {
                        SSq.set(i, colctr + k, stateRow.get(0, k))
                    }
                    colctr += space[j].getNumCols()
                }
            }
            sn.space = space.associateBy { space.indexOf(it) }.mapKeys { sn.stateful[it.key] }.toMutableMap()
            return Triple(SSq, SSh, sn)
        }
        
        // Pop state from stack
        val stateCell = stack.removeAt(stack.size - 1)
        ih = stackIndex.removeAt(stackIndex.size - 1)
        
        val newStateCell = Array(A) { stateCell.toMutableList() }
        
        val enabledSync = mutableListOf<Int>()
        val enabledRates = mutableListOf<Double>()
        
        // Process each synchronization action
        for (act in 0 until A) {
            val updateCondA = true
            
            if (updateCondA) {
                val isf = sn.nodeToStateful.get(nodeA[act]).toInt()
                val activeResult = State.afterEvent(sn, nodeA[act], stateCell[isf], eventA[act], classA[act], isSimulation)
                
                if (activeResult.outspace.isEmpty() || activeResult.outrate.isEmpty()) {
                    continue
                }
                
                newStateCell[act][sn.nodeToStateful.get(nodeA[act]).toInt()] = activeResult.outspace
                val rateA = activeResult.outrate
                val outprobA = activeResult.outprob
                
                for (ia in 0 until activeResult.outspace.getNumRows()) {
                    if (activeResult.outspace.getRow(ia).elementSum() == -1.0) {
                        continue
                    }
                    
                    if (rateA.get(ia, 0) > 0) {
                        var probSyncP = 1.0
                        
                        if (nodeP[act] != local) {
                            if (nodeP[act] == nodeA[act]) {
                                // Self-loop case - simplified implementation
                                val passiveResult = State.afterEvent(sn, nodeP[act], newStateCell[act][sn.nodeToStateful.get(nodeA[act]).toInt()], eventP[act], classP[act], isSimulation)
                                newStateCell[act][sn.nodeToStateful.get(nodeP[act]).toInt()] = passiveResult.outspace
                            } else {
                                // Departure from active
                                val passiveResult = State.afterEvent(sn, nodeP[act], stateCell[sn.nodeToStateful.get(nodeP[act]).toInt()], eventP[act], classP[act], isSimulation)
                                newStateCell[act][sn.nodeToStateful.get(nodeP[act]).toInt()] = passiveResult.outspace
                            }
                            
                            if (!newStateCell[act][sn.nodeToStateful.get(nodeP[act]).toInt()].isEmpty()) {
                                probSyncP = if (sn.isstatedep.get(nodeA[act], 2) == 1.0) {
                                    // State-dependent probability - use default value for now
                                    sync[act]!!.passive[0]!!.prob
                                } else {
                                    sync[act]!!.passive[0]!!.prob
                                }
                            } else {
                                probSyncP = 0.0
                            }
                        }
                        
                        if (!newStateCell[act][sn.nodeToStateful.get(nodeA[act]).toInt()].isEmpty()) {
                            if (nodeP[act] == local) {
                                probSyncP = 1.0
                            }
                            
                            if (!rateA.get(ia, 0).isNaN()) {
                                var allNonEmpty = true
                                for (stateMatrix in newStateCell[act]) {
                                    if (stateMatrix.isEmpty()) {
                                        allNonEmpty = false
                                        break
                                    }
                                }
                                if (allNonEmpty) {
                                    // Check class switching mask
                                    if (nodeP[act] < local && csmask.get(classA[act], classP[act]) == 0.0 && 
                                        sn.nodetype[nodeP[act]] != NodeType.Source && 
                                        (rateA.get(ia, 0) * probSyncP > 0)) {
                                        throw RuntimeException("Error: state-dependent routing at node ${nodeA[act]} violates the class switching mask")
                                    }
                                    
                                    enabledRates.add(rateA.get(ia, 0) * probSyncP)
                                    enabledSync.add(act)
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Process enabled transitions
        for (firingCtr in enabledRates.indices) {
            val firingRate = enabledRates[firingCtr]
            val act = enabledSync[firingCtr]
            val netstates = newStateCell[act]
            
            var allNetStatesNonEmpty = true
            for (stateMatrix in netstates) {
                if (stateMatrix.isEmpty()) {
                    allNetStatesNonEmpty = false
                    break
                }
            }
            if (firingRate > 0 && allNetStatesNonEmpty) {
                val nvec = Matrix(1, netstates.size)
                for (i in netstates.indices) {
                    nvec.set(0, i, (netstates[i].getNumRows() - 1).toDouble())
                }
                var n = pprod(nvec)
                
                while (n != null && !n.isEmpty()) {
                    val newstatec = mutableListOf<Matrix>()
                    
                    for (i in netstates.indices) {
                        val nIdx = n.get(0, i).toInt()
                        maxstatesz[i] = maxOf(maxstatesz[i], netstates[i].getRow(nIdx).getNumCols())
                        val paddedState = Matrix(1, maxstatesz[i])
                        paddedState.fill(0.0)
                        val originalState = netstates[i].getRow(nIdx)
                        for (j in 0 until originalState.getNumCols()) {
                            paddedState.set(0, j + (maxstatesz[i] - originalState.getNumCols()), originalState.get(0, j))
                        }
                        newstatec.add(paddedState)
                    }
                    
                    val hashednewstate = Matrix(1, nstateful)
                    for (i in 0 until nstateful) {
                        val matchIdx = findMatchingRow(space[i], newstatec[i])
                        hashednewstate.set(0, i, (matchIdx + 1).toDouble()) // Convert to 1-based
                    }
                    
                    val jh = findMatchingRow(SSh, hashednewstate)
                    if (jh == -1) {
                        // Add new state to space if needed
                        for (i in 0 until nstateful) {
                            if (hashednewstate.get(0, i) == 0.0) { // matchIdx was -1
                                val newSpace = Matrix(space[i].getNumRows() + 1, space[i].getNumCols())
                                for (r in 0 until space[i].getNumRows()) {
                                    for (c in 0 until space[i].getNumCols()) {
                                        newSpace.set(r, c, space[i].get(r, c))
                                    }
                                }
                                for (c in 0 until newstatec[i].getNumCols()) {
                                    newSpace.set(space[i].getNumRows(), c, newstatec[i].get(0, c))
                                }
                                space[i] = newSpace
                                hashednewstate.set(0, i, space[i].getNumRows().toDouble())
                            }
                        }
                        
                        val newSSh = Matrix(SSh.getNumRows() + 1, SSh.getNumCols())
                        for (r in 0 until SSh.getNumRows()) {
                            for (c in 0 until SSh.getNumCols()) {
                                newSSh.set(r, c, SSh.get(r, c))
                            }
                        }
                        for (c in 0 until hashednewstate.getNumCols()) {
                            newSSh.set(SSh.getNumRows(), c, hashednewstate.get(0, c))
                        }
                        stack.add(newstatec)
                        stackIndex.add(newSSh.getNumRows())
                    }
                    n = pprod(n, nvec)
                }
            }
        }
    }
    
    return Triple(SSq, SSh, sn)
}

/**
 * Helper function to find matching row in matrix (equivalent to MATLAB matchrow)
 */
private fun findMatchingRow(matrix: Matrix, targetRow: Matrix): Int {
    if (targetRow.getNumRows() != 1) {
        throw IllegalArgumentException("Target must be a single row")
    }
    
    for (i in 0 until matrix.getNumRows()) {
        if (matrix.getNumCols() == targetRow.getNumCols()) {
            var match = true
            for (j in 0 until matrix.getNumCols()) {
                if (matrix.get(i, j) != targetRow.get(0, j)) {
                    match = false
                    break
                }
            }
            if (match) {
                return i
            }
        }
    }
    return -1
}