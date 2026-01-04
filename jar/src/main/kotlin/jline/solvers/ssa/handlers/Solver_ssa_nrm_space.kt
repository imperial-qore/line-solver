package jline.solvers.ssa.handlers

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.GlobalConstants.NegInf
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginalAggr
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.RandomManager
import jline.util.matrix.Matrix
import kotlin.math.min

data class SolverSSAResultNRMSpace(val pi: Matrix, val outspace: Matrix, val depRates: Matrix, val sn: NetworkStruct)

fun solver_ssa_nrm_space(sn: NetworkStruct, options: SolverOptions): SolverSSAResultNRMSpace {

    // Set master seed for reproducible SSA simulation
    RandomManager.setMasterSeed(options.seed)

    // ---------------------------------------------------------------------
    // Parameters & shorthands
    // ---------------------------------------------------------------------
    val samples = options.samples
    val R = sn.nclasses.toInt()
    val I = sn.nnodes.toInt()
    val state = sn.state

    // ---------------------------------------------------------------------
    // Stoichiometry & reaction mapping (self-loops included)
    // ---------------------------------------------------------------------
    var S = Matrix(0, I * R)  // will transpose at the end
    val fromIdx = mutableListOf<Int>()
    val toIdx = mutableListOf<MutableList<Int>>()
    val fromIR = mutableListOf<IntArray>()
    val probIR = mutableListOf<MutableList<Double>>()

    // this is currently M^2*R^2, it can be lowered to M*R decoupling the routing
    var k = 0
    for (ind in 0 until I) {
        for (r in 0 until R) {
            k++
            fromIR.add(intArrayOf(ind, r))
            fromIdx.add(ind * R + r)
            probIR.add(mutableListOf())
            toIdx.add(mutableListOf())
            val Srow = DoubleArray(I * R) // build stoichiometry row
            if (sn.isslc.get(r) != 0.0) {
                Srow[fromIdx[k - 1]] = NegInf
            } else {
                Srow[fromIdx[k - 1]] = -1.0
                for (jnd in 0 until I) {
                    for (s in 0 until R) {
                        if (sn.rtnodes.get(ind * R + r, jnd * R + s) > 0) {
                            toIdx[k - 1].add(jnd * R + s)
                            val p = sn.rtnodes.get(ind * R + r, jnd * R + s)
                            probIR[k - 1].add(p)
                            Srow[jnd * R + s] = Srow[jnd * R + s] + p
                        }
                    }
                }
            }
            // Expand S matrix if needed
            if (k > S.numRows) {
                val newS = Matrix(k, I * R)
                for (i in 0 until S.numRows) {
                    for (j in 0 until S.numCols) {
                        newS.set(i, j, S.get(i, j))
                    }
                }
                S = newS
            }
            // Set the row
            for (j in Srow.indices) {
                S.set(k - 1, j, Srow[j])
            }
        }
    }
    S = S.transpose() // states × reactions

    // ---------------------------------------------------------------------
    // Initial state vector
    // ---------------------------------------------------------------------
    val nvec0 = Matrix(I * R, 1) // initial state (aggregate state)
    for (ind in 0 until I) {
        if (sn.isstateful.get(ind) != 0.0) {
            val aggr = toMarginalAggr(sn,
                ind,
                state[sn.stateful.get(sn.nodeToStateful.get(ind).toInt())],
                null,
                null,
                null,
                null,
                null)
            for (r in 0 until R) {
                var nir = aggr.nir[r]
                if (nir.isInfinite()) {
                    if (sn.nodetype.get(ind) == NodeType.Source) {
                        nir = 1.0
                    } else {
                        error("Infinite population error.")
                    }
                }
                nvec0.set(ind * R + r, 0, nir)
            }
        }
    }

    // ---------------------------------------------------------------------
    // Rates and servers setup
    // ---------------------------------------------------------------------
    val mi = DoubleArray(I)
    val rates = Array(I) { DoubleArray(R) }
    for (ind in 0 until I) {
        if (sn.isstation.get(ind) != 0.0) {
            for (r in 0 until R) {
                val ist = sn.nodeToStation.get(ind).toInt()
                val muir = sn.rates.get(ist, r)
                if (!muir.isNaN()) {
                    rates[ind][r] = muir
                }
                mi[ind] = sn.nservers.get(ist)
            }
        } else {
            for (r in 0 until R) {
                rates[ind][r] = GlobalConstants.Immediate
                mi[ind] = GlobalConstants.MaxInt.toDouble()
            }
        }
        if (mi[ind].isInfinite()) {
            mi[ind] = GlobalConstants.MaxInt.toDouble()
        }
    }

    // ---------------------------------------------------------------------
    // Propensity functions
    // ---------------------------------------------------------------------
    val epstol = GlobalConstants.Zero
    val a = Array<(Matrix) -> Double>(fromIdx.size) { j ->
        { X: Matrix ->
            val ind = fromIR[j][0]
            val r = fromIR[j][1]
            val base = rates[ind][r]
            when {
                sn.isstation.get(ind) == 0.0 -> base * min(1.0, X.get(fromIdx[j], 0))
                else -> {
                    when (sn.sched.get(sn.stations.get(sn.nodeToStation.get(ind).toInt()))) {
                        SchedStrategy.EXT -> base
                        SchedStrategy.INF -> base * X.get(fromIdx[j], 0)
                        SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                        SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                            if (R == 1) { // single class
                                base * min(mi[ind], X.get(fromIdx[j], 0))
                            } else {
                                val total = (0 until R).sumOf { X.get(ind * R + it, 0) } + epstol
                                base * (X.get(fromIdx[j], 0) / total) * min(mi[ind], total)
                            }
                        }
                        else -> base * X.get(fromIdx[j], 0)
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Propensity function dependencies
    // ---------------------------------------------------------------------
    val D = Array(S.numCols) { mutableListOf<Int>() }
    for (k in 0 until D.size) {
        val J = mutableListOf<Int>() // set of state variables affected by reaction k
        for (i in 0 until S.numRows) {
            if (S.get(i, k) != 0.0) {
                J.add(i)
            }
        }
        val vecd = mutableListOf<Int>()
        for (j in J.indices) {
            // (ind-1)*R + r -> ind*R + r in 0-based indexing
            val pos = J[j]
            val r = pos % R
            val ind = (pos - r) / R
            for (rr in 0 until R) {
                vecd.add(ind * R + rr)
            }
        }
        // vecd now contains all state variables affected by the firing of
        // reaction k. We now find the propensity functions that depend
        // on those variables
        if (vecd.isNotEmpty()) {
            val vecdUnique = vecd.distinct()
            val vecs = mutableListOf<Int>()
            for (j in vecdUnique.indices) {
                for (kk in 0 until S.numCols) {
                    if (S.get(vecdUnique[j], kk) < 0.0) {
                        vecs.add(kk)
                    }
                }
            }
            D[k].addAll(vecs.distinct())
        }
    }

    // Having accounted for them in D, we can now remove self-loops markings
    for (i in 0 until S.numRows) {
        for (j in 0 until S.numCols) {
            if (S.get(i, j).isInfinite()) {
                S.set(i, j, 0.0)
            }
        }
    }

    // ---------------------------------------------------------------------
    // Run SSA/NRM
    // ---------------------------------------------------------------------
    val reactCache = mutableMapOf<String, DoubleArray>()
    val hashfun = { v: Matrix -> v.toString() }
    val (t, nvecsim) = next_reaction_method_space(S, D, a, nvec0, samples, options, reactCache, hashfun)

    // ---------------------------------------------------------------------
    // Empirical state probabilities
    // ---------------------------------------------------------------------
    val dt = DoubleArray(t.size - 1) { i -> t[i + 1] - t[i] }

    // Extract state columns (equivalent to nvecsim(:,1:end-1).')
    val stateCols = (0 until nvecsim.numCols - 1).map { c ->
        DoubleArray(nvecsim.numRows) { r -> nvecsim.get(r, c) }
    }

    // Find unique states (equivalent to unique(..., 'rows', 'stable'))
    val uniq = linkedMapOf<String, Int>()
    val outspaceRows = mutableListOf<DoubleArray>()
    val ic = IntArray(stateCols.size)
    stateCols.forEachIndexed { i, arr ->
        val key = arr.joinToString(",")
        ic[i] = uniq.getOrPut(key) { outspaceRows += arr; outspaceRows.size - 1 }
    }

    // Accumulate time (equivalent to accumarray(ic, dt))
    val timeAccum = DoubleArray(outspaceRows.size)
    ic.forEachIndexed { i, id -> timeAccum[id] += dt[i] }

    // Normalize to get probabilities
    val total = timeAccum.sum()
    val pi = Matrix(1, timeAccum.size).apply {
        timeAccum.forEachIndexed { i, v -> set(0, i, v / total) }
    }

    val outspace = Matrix(outspaceRows.size, I * R)
    outspaceRows.forEachIndexed { i, row ->
        row.forEachIndexed { j, v -> outspace.set(i, j, v) }
    }

    // ---------------------------------------------------------------------
    // Per-state departure rates (self-loops counted)
    // ---------------------------------------------------------------------
    val numStates = outspace.numRows
    val depRates = Matrix(numStates, I * R)

    for (st in 0 until numStates) {
        val stateVec = Matrix(I * R, 1)
        for (i in 0 until I * R) {
            stateVec.set(i, 0, outspace.get(st, i))
        }
        val a_state = reactCache[hashfun(stateVec)] ?: continue
        for (j in fromIdx.indices) {
            depRates.set(st, fromIdx[j], depRates.get(st, fromIdx[j]) + a_state[j])
        }
    }

    return SolverSSAResultNRMSpace(pi, outspace, depRates, sn)
}

// ======================================================================
// Next-Reaction Method core
// ======================================================================
fun next_reaction_method_space(S: Matrix,
                               D: Array<MutableList<Int>>,
                               a: Array<(Matrix) -> Double>,
                               nvec0: Matrix,
                               samples: Int,
                               options: SolverOptions,
                               reactcache: MutableMap<String, DoubleArray>,
                               hashfun: (Matrix) -> String): Pair<List<Double>, Matrix> {
    val numReactions = S.numCols
    val Ak = DoubleArray(numReactions) { a[it](nvec0) }
    val Pk = DoubleArray(numReactions) { -kotlin.math.ln(Maths.rand()) }
    val Tk = DoubleArray(numReactions)
    val times = mutableListOf(0.0)
    val states = mutableListOf<DoubleArray>()
    var nvec = nvec0.copy()
    reactcache.clear()
    reactcache[hashfun(nvec)] = Ak.copyOf()

    var n = 0
    while (n < samples) {
        val tau = DoubleArray(numReactions) { i ->
            if (Ak[i] > 0) (Pk[i] - Tk[i]) / Ak[i] else Inf
        }
        val kfire = tau.withIndex().minByOrNull { it.value }!!.index
        val dt = tau[kfire]
        times += times.last() + dt

        // Update state
        nvec.addEq(S.getColumn(kfire))

        // Update times and propensities
        for (i in 0 until numReactions) Tk[i] += Ak[i] * dt
        D[kfire].forEach { i -> Ak[i] = a[i](nvec) }
        Pk[kfire] -= kotlin.math.ln(Maths.rand())

        // Cache the propensities for this state
        reactcache[hashfun(nvec)] = Ak.copyOf()
        states += DoubleArray(nvec.numRows) { r -> nvec.get(r, 0) }
        n++
        printProgress(options, n)
    }

    val result = Matrix(S.numRows, states.size + 1)
    // Set initial state
    for (r in 0 until nvec0.numRows) {
        result.set(r, 0, nvec0.get(r, 0))
    }
    // Set subsequent states
    states.forEachIndexed { c, col ->
        col.forEachIndexed { r, v -> result.set(r, c + 1, v) }
    }
    return times to result
}

private fun printProgress(options: SolverOptions, samples_collected: Int) {
    if (options.method != "parallel" && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
        if (samples_collected == 2) {
            System.out.printf("SSA samples: %9d ", samples_collected)
            System.out.flush()
        } else if (options.verbose == VerboseLevel.DEBUG) {
            System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected)
            System.out.flush()
        } else if (samples_collected % 1000 == 0) {
            System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected)
            System.out.flush()
        }
        if (samples_collected == options.samples) {
            println()
        }
    }
}