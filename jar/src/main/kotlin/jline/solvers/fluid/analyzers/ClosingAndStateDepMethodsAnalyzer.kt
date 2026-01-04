/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers

import jline.api.mam.map_mean
import jline.lang.JobClass
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.lang.nodes.Station
import jline.solvers.SolverOptions
import jline.solvers.SolverResult
import jline.solvers.fluid.handlers.PassageTimeODE
import jline.solvers.fluid.handlers.TransientDataHandler
import jline.util.matrix.Matrix
import odesolver.LSODA
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import java.util.*
import kotlin.Array
import kotlin.Boolean
import kotlin.Double.Companion.POSITIVE_INFINITY
import kotlin.DoubleArray
import kotlin.Int
import kotlin.RuntimeException
import kotlin.arrayOfNulls
import kotlin.doubleArrayOf
import kotlin.math.max


class ClosingAndStateDepMethodsAnalyzer : FluidAnalyzer {
    var xvec_t: Matrix? = null
    var xvec_it: Matrix? = null

    private fun solver_fluid_iteration(sn: NetworkStruct,
                                       mu: MutableMap<Station?, MutableMap<JobClass?, Matrix?>?>?,
                                       phi: MutableMap<Station?, MutableMap<JobClass?, Matrix?>?>?,
                                       S: Matrix?,
                                       yDefault: DoubleArray,
                                       slowrate: Matrix,
                                       options: SolverOptions,
                                       result: SolverResult) {
        var iter = 0
        var totalSteps = 0

        // Heuristic to select stiff or non-stiff ODE solver
        val slowrateRows = slowrate.numRows
        val slowrateCols = slowrate.numCols
        var minNonZeroRate = Inf

        var maxNonZeroRate = 0.0
        for (i in 0..<slowrateRows) {
            for (j in 0..<slowrateCols) {
                if ((slowrate.get(i, j) > GlobalConstants.CoarseTol) && Double.isFinite(slowrate.get(i, j))) {
                    if (slowrate.get(i, j) < minNonZeroRate) {
                        minNonZeroRate = slowrate.get(i, j)
                    }
                    if (slowrate.get(i, j) > maxNonZeroRate) {
                        maxNonZeroRate = slowrate.get(i, j)
                    }
                }
            }
        }

        // Initialise ODE
        val ode: FirstOrderDifferentialEquations = PassageTimeODE(sn, mu, phi, sn.proc, sn.rt, S, options)

        //decide whether to use stiff or non-stiff method
        options.stiff = detectStiffnessUsingOstrowski(sn, slowrate)

        var T0 = options.timespan[0]
        var T = 0.0

        val tIterations: MutableList<Matrix> = LinkedList<Matrix>()
        val xVecIterations: MutableList<Matrix> = LinkedList<Matrix>()

        // Matching MATLAB's iterative approach (solver_fluid_iteration.m lines 34-96)
        // MATLAB does not break on convergence in this function - it runs until T >= timespan[1] or iter >= iter_max
        var goon = true

        while ((Double.isFinite(options.timespan[1]) && T < options.timespan[1]) || (goon && iter < options.iter_max)) {
            iter++

            // Determine entry state vector
            val initialState = DoubleArray(yDefault.size)
            val nextState = DoubleArray(yDefault.size)
            for (i in yDefault.indices) {
                initialState[i] = xvec_it!!.get(0, i)
                nextState[i] = 0.0
            }

            // Solve ode until T = time for slowest exit rate (MATLAB lines 44-48)
            // First iteration: T = min(timespan(2), abs(10/min(nonZeroRates)))
            // Subsequent iterations: T = min(timespan(2), abs(10*iter/min(nonZeroRates)))
            T = if (iter == 1) {
                FastMath.min(options.timespan[1], FastMath.abs(10.0 / minNonZeroRate))
            } else {
                FastMath.min(options.timespan[1], FastMath.abs(10.0 * iter / minNonZeroRate))
            }
            val tRange = doubleArrayOf(T0, T)

            if (options.tol > GlobalConstants.CoarseTol && (options.verbose == VerboseLevel.DEBUG)) {
                System.err.println("Fast, non-stiff ODE solver is not yet available in JLINE. Using accurate non-stiff ODE solver instead.")
            }

            var Tmax: Int
            if (options.stiff) {
                val odeSolver: LSODA = if (options.tol > GlobalConstants.CoarseTol) {
                    options.odesolvers.fastStiffODESolver
                } else {
                    options.odesolvers.accurateStiffODESolver
                }

                try {
                    odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)
                } catch (e: RuntimeException) {
                    if (options.verbose != VerboseLevel.SILENT) {
                        println("The initial point is invalid, Fluid solver switching to default initialization.")
                    }
                    odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState)
                }
                // transform the output format
                Tmax = odeSolver.stepsTaken + 1
                val tVec = Matrix(Tmax, 1)
                val xVec = Matrix(Tmax, ode.dimension)
                val tHistory = odeSolver.tvec
                val yHistory = odeSolver.yvec
                for (i in 0..<Tmax) {
                    tVec.set(i, 0, tHistory.get(i)!!)
                    for (j in 0..<ode.dimension) xVec.set(i, j, max(0.0, yHistory.get(i)!![j]!!))
                }
                tIterations.add(tVec)
                xVecIterations.add(xVec)
                this.xvec_it = Matrix.extractRows(xVec, Tmax - 1, Tmax, null)
            } else {
                val odeSolver: FirstOrderIntegrator = options.odesolvers.accurateODESolver
                odeSolver.clearStepHandlers()
                val stepHandler = TransientDataHandler(initialState.size)
                odeSolver.addStepHandler(stepHandler)

                try {
                    odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)
                } catch (e: RuntimeException) {
                    if (options.verbose != VerboseLevel.SILENT) {
                        println("The initial point is invalid, Fluid solver switching to default initialization.")
                    }
                    odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState)
                }

                // Retrieve Transient Data
                tIterations.add(stepHandler.tVec)
                xVecIterations.add(stepHandler.xVec)
                Tmax = stepHandler.tVec.numRows
                this.xvec_it = Matrix.extractRows(stepHandler.xVec, Tmax - 1, Tmax, null)
            }
            totalSteps += Tmax

            // Update T0 for next iteration
            T0 = T

            // Check if we've reached the end time (MATLAB line 93-95)
            if (T >= options.timespan[1]) {
                goon = false
            }
        }

        // Migrating tIterations and xVecIterations from Lists to concatenated JLineMatrix objects
        // Note this is done manually for performance purposes - concatRows is not as efficient
        var nextRow = 0
        val cols = xVecIterations[0].numCols
        this.xvec_t = Matrix(totalSteps, cols)
        result.t = Matrix(totalSteps, 1)
        for (i in 0..<iter) {
            val tIter = tIterations[i]
            val xVecIter = xVecIterations[i]
            val stepsPerIter = tIter.numRows
            for (j in nextRow..<nextRow + stepsPerIter) {
                result.t.set(j, 0, tIter.get(j - nextRow, 0))
                for (k in 0..<cols) {
                    xvec_t!!.set(j, k, xVecIter.get(j - nextRow, k))
                }
            }
            nextRow += stepsPerIter
        }
    }

    private fun solver_fluid(sn: NetworkStruct, options: SolverOptions, result: SolverResult) {
        val M = sn.nstations // Number of stations
        val K = sn.nclasses // Number of classes
        val S = sn.nservers.copy()

        // Making deep copies as going forwards, mu and phi are amended
        val mu: MutableMap<Station?, MutableMap<JobClass?, Matrix?>?> =
            HashMap<Station?, MutableMap<JobClass?, Matrix?>?>()
        val phi: MutableMap<Station?, MutableMap<JobClass?, Matrix?>?> =
            HashMap<Station?, MutableMap<JobClass?, Matrix?>?>()
        for (i in 0..<M) {
            val station = sn.stations.get(i)
            val muCopy: MutableMap<JobClass?, Matrix?> = HashMap<JobClass?, Matrix?>()
            val phiCopy: MutableMap<JobClass?, Matrix?> = HashMap<JobClass?, Matrix?>()
            for (k in 0..<K) {
                val jobClass = sn.jobclasses.get(k)
                muCopy.put(jobClass, sn.mu.get(station)!!.get(jobClass)!!.copy())
                phiCopy.put(jobClass, sn.phi.get(station)!!.get(jobClass)!!.copy())
            }
            mu.put(station, muCopy)
            phi.put(station, phiCopy)
        }

        val match = Matrix(M, K) // Indicates whether a class is served at a station
        val phases = Matrix(M, K)
        val slowrate = Matrix(M, K)
        for (i in 0..<M) {
            val station = sn.stations.get(i)
            for (k in 0..<K) {
                val jobClass = sn.jobclasses.get(k)

                if (mu.get(station)!!.get(jobClass)!!.hasNaN()) {
                    mu.get(station)!!.put(jobClass, Matrix(0, 0))
                    phi.get(station)!!.put(jobClass, Matrix(0, 0))
                }

                if (sn.rt.sumCols((i * K) + k) > 0) {
                    match.set(i, k, 1)
                }

                // Set number of servers in delay station = population
                if (Double.isInfinite(S.get(i, 0))) {
                    S.set(i, 0, sn.nclosedjobs)
                }

                phases.set(i, k, mu.get(station)!!.get(jobClass)!!.length())

                if (mu.get(station)!!.get(jobClass)!!.isEmpty) {
                    slowrate.set(i, k, POSITIVE_INFINITY)
                } else {
                    // service completion (exit) rates in each phase
                    slowrate.set(i, k, mu.get(station)!!.get(jobClass)!!.elementMin())
                }
            }
        }

        val y0 = Matrix(1, 0)
        val assigned = Matrix(1, K) // number of jobs of each class already assigned
        var toAssign: kotlin.Double
        for (i in 0..<M) {
            for (k in 0..<K) {
                if (match.get(i, k) > 0) { // indicates whether a class is served at a station
                    if (Double.isInfinite(sn.njobs.get(0, k))) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                            toAssign = 1.0 // open job pool
                        } else {
                            toAssign = 0.0 // set to zero open jobs everywhere
                        }
                    } else {
                        toAssign = FastMath.floor(sn.njobs.get(0, k) / match.sumCols(k))
                        // If this is the last station for this job class
                        if (match.sumSubMatrix(i + 1, match.numRows, k, k + 1) == 0.0) {
                            toAssign = sn.njobs.get(0, k) - assigned.get(0, k)
                        }
                    }

                    // This is just for PS
                    val originalY0Length = y0.numCols
                    val newY0Length = originalY0Length + 1 + (phases.get(i, k) - 1).toInt()
                    y0.expandMatrix(1, newY0Length, newY0Length)
                    y0.set(0, originalY0Length, toAssign)
                    var col = 0
                    while (col < phases.get(i, k) - 1) {
                        y0.set(0, originalY0Length + 1 + col, 0)
                        col++
                    }
                    assigned.set(0, k, assigned.get(0, k) + toAssign)
                } else {
                    val originalY0Length = y0.numCols
                    val newY0Length = originalY0Length + phases.get(i, k).toInt()
                    y0.expandMatrix(1, newY0Length, newY0Length)
                    var col = 0
                    while (col < phases.get(i, k)) {
                        y0.set(0, originalY0Length + col, 0)
                        col++
                    }
                }
            }
        }

        val y0cols = y0.numCols
        val yDefault = DoubleArray(y0cols)
        if (options.init_sol.isEmpty) {
            xvec_it = y0 // average state embedded at stage change transitions out of e
        } else {
            xvec_it = options.init_sol
            for (i in 0..<y0cols) {
                yDefault[i] = y0.get(0, i) // default solution if init_sol fails
            }
        }

        // Solve ODE
        solver_fluid_iteration(sn, mu, phi, S, yDefault, slowrate, options, result)

        // This part assumes PS, DPS, GPS scheduling
        result.QN = Matrix(M, K)
        result.QNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        val Qt = arrayOfNulls<Matrix>(M)
        result.UNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        val Tmax = xvec_t!!.numRows

        for (i in 0..<M) {
            Qt[i] = Matrix(Tmax, 1)
            for (k in 0..<K) {
                val shift =
                    phases.sumSubMatrix(0, i, 0, phases.numCols).toInt() + phases.sumSubMatrix(i, i + 1, 0, k)
                        .toInt()
                result.QN.set(i, k, xvec_it!!.sumSubMatrix(0, 1, shift, shift + phases.get(i, k).toInt()))

                result.QNt[i][k] = Matrix(Tmax, 1)
                for (step in 0..<Tmax) {
                    result.QNt[i][k].set(step,
                        0,
                        xvec_t!!.sumSubMatrix(step, step + 1, shift, shift + phases.get(i, k).toInt()))
                }

                Qt[i] = Qt[i]!!.add(1.0, result.QNt[i][k])
                // computes queue length in each node and stage, summing the total number in service and
                // waiting in that station results are weighted with the stat prob of the stage
            }
        }

        for (i in 0..<M) {
            if (sn.nservers.get(i, 0) > 0) { // Not INF
                for (k in 0..<K) {
                    result.UNt[i][k] = Matrix(Tmax, 1)
                    for (step in 0..<Tmax) {
                        // If not an infinite server then this is a number between 0 and 1
                        result.UNt[i][k].set(step,
                            0,
                            FastMath.min(result.QNt[i][k].get(step, 0) / S.get(i, 0),
                                result.QNt[i][k].get(step, 0) / Qt[i]!!.get(step, 0)))
                        if (Double.isNaN(result.UNt[i][k].get(step, 0))) {
                            result.UNt[i][k].set(step, 0, 0) // fix cases where QLen is 0
                        }
                    }
                }
            } else {
                // Infinite server
                System.arraycopy(result.QNt[i], 0, result.UNt[i], 0, K)
            }
        }
    }

    override fun analyze(sn: NetworkStruct, options: SolverOptions, result: SolverResult) {
        val M = sn.nstations // Number of stations
        val K = sn.nclasses // Number of classes
        val lambda = sn.mu

        // Inner iteration of fluid analysis
        solver_fluid(sn, options, result)

        // Assumes the existence of a delay node through which all classes pass
        val delayNodes = Matrix(1, M)
        for (i in 0..<M) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                delayNodes.set(0, i, 1)
            }
        }

        // THROUGHPUT - for all classes in each station
        result.TN = Matrix(M, K)
        result.TNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        val Tmax = result.t.numRows
        for (i in 0..<M) {
            for (k in 0..<K) {
                result.TNt[i][k] = Matrix(Tmax, 1)
            }
        }

        val Xservice = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }  // Throughput per class, station and phase

        for (i in 0..<M) {
            if (delayNodes.get(0, i) == 1.0) {
                for (k in 0..<K) {
                    val idx =
                        sn.phases.sumSubMatrix(0, i, 0, K).toInt() + sn.phases.sumSubMatrix(i, i + 1, 0, k).toInt()
                    Xservice[i]!![k] = Matrix(sn.phases.get(i, k).toInt(), 1)
                    for (f in 0..<sn.phases.get(i, k).toInt()) {
                        val lambdaIKF = lambda.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(f, 0)
                        val phiIKF = sn.phi.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(f, 0)

                        result.TN.set(i, k, result.TN.get(i, k) + (xvec_it!!.get(0, idx + f) * lambdaIKF * phiIKF))

                        val tmpTNt = result.QNt[i][k].copy()
                        tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt)
                        result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt)

                        Xservice[i]!![k]!!.set(f, 0, (xvec_it!!.get(0, idx + f) * lambdaIKF))
                    }
                }
            } else {
                val xi = result.QN.sumRows(i) // Number of jobs in the station
                var xi_t = result.QNt[i][0].copy()
                for (r in 1..<K) {
                    xi_t = xi_t.add(1.0, result.QNt[i][r])
                }

                var wni = GlobalConstants.FineTol
                val wi = Matrix(1, K)
                if (xi > 0 || sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                        for (k in 0..<K) {
                            val idx =
                                sn.phases.sumSubMatrix(0, i, 0, K).toInt() + sn.phases.sumSubMatrix(i, i + 1, 0, k)
                                    .toInt()
                            val D0 = sn.proc.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(0)
                            val D1 = sn.proc.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(1)
                            wi.set(0, k, map_mean(D0, D1))

                            wni += wi.get(0, k) * xvec_it!!.sumSubMatrix(0, 1, idx, idx + sn.phases.get(i, k).toInt())
                        }
                    }

                    for (k in 0..<K) {
                        val idx =
                            sn.phases.sumSubMatrix(0, i, 0, K).toInt() + sn.phases.sumSubMatrix(i, i + 1, 0, k).toInt()
                        Xservice[i]!![k] = Matrix(sn.phases.get(i, k).toInt(), 1)

                        for (f in 0..<sn.phases.get(i, k).toInt()) {
                            val lambdaIKF = lambda.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(f, 0)
                            val phiIKF = sn.phi.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(f, 0)

                            when (sn.sched.get(sn.stations.get(i))) {
                                SchedStrategy.EXT -> if (f == 0) {
                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + (lambdaIKF * phiIKF * (1 - xvec_it!!.sumSubMatrix(0,
                                            1,
                                            idx + 1,
                                            idx + sn.phases.get(i, k).toInt()))))

                                    var tmpTNt = xvec_t!!.sumRows(idx + 1, idx + sn.phases.get(i, k).toInt())
                                    val ones = Matrix(tmpTNt.numRows, tmpTNt.numCols)
                                    ones.ones()
                                    tmpTNt = ones.sub(tmpTNt)
                                    tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt)
                                    result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt)

                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        lambdaIKF * (1 - xvec_it!!.sumSubMatrix(0,
                                            1,
                                            idx + 1,
                                            idx + sn.phases.get(i, k).toInt())))
                                } else {
                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + (lambdaIKF * phiIKF * xvec_it!!.get(0, idx + f)))
                                    val tmpTNt = xvec_t!!.sumRows(idx + f, idx + f + 1)
                                    tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt)
                                    result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt)
                                    Xservice[i]!![k]!!.set(f, 0, lambdaIKF * xvec_it!!.get(0, idx + f))
                                }

                                SchedStrategy.INF, SchedStrategy.PS -> {
                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + ((lambdaIKF * phiIKF / xi) * xvec_it!!.get(0,
                                            idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))

                                    val tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT)
                                    val rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / xi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it!!.get(0,
                                            idx + f))
                                }

                                SchedStrategy.DPS -> {
                                    val w = Matrix(1, K)
                                    var p = 0
                                    while (p < K) {
                                        w.set(0, p, sn.schedparam.get(i, p))
                                        p++
                                    }

                                    var tmpQ: Matrix? = Matrix(1, result.QN.numCols)
                                    Matrix.extractRows(result.QN, i, i + 1, tmpQ)
                                    tmpQ = tmpQ!!.transpose()
                                    var wxi = Matrix(1, 1) // number of jobs in the station
                                    wxi = w.mult(tmpQ, wxi)

                                    var wxi_t = result.QNt[i][0].copy()
                                    wxi_t.scaleEq(w.value(), wxi_t)
                                    var r = 1
                                    while (r < result.QNt[i][0].numCols) {
                                        val tmpWxi_t = result.QNt[i][r].copy()
                                        tmpWxi_t.scaleEq(w.get(0, r), tmpWxi_t)
                                        wxi_t = wxi_t.add(1.0, tmpWxi_t)
                                        r++
                                    }

                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + (((lambdaIKF * phiIKF * w.get(0,
                                            k)) / wxi.value()) * xvec_it!!.get(0, idx + f) * FastMath.min(xi,
                                            sn.nservers.get(i, 0))))

                                    var tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF * w.get(0, k), tmpTNT)
                                    var rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / wxi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        ((lambdaIKF * w.get(0, k) / wxi.value()) * FastMath.min(xi,
                                            sn.nservers.get(i, 0)) * xvec_it!!.get(0, idx + f)))
                                }

                                SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFS -> {
                                    var tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT)
                                    var rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / xi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    if (options.method == "default" || options.method == "closing") {
                                        result.TN.set(i,
                                            k,
                                            result.TN.get(i, k) + ((lambdaIKF * phiIKF / xi) * xvec_it!!.get(0,
                                                idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))
                                        Xservice[i]!![k]!!.set(f,
                                            0,
                                            lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it!!.get(0,
                                                idx + f))
                                        break
                                    } else if (options.method == "statedep") {
                                        result.TN.set(i,
                                            k,
                                            result.TN.get(i, k) + (((lambdaIKF * phiIKF * wi.get(0,
                                                k)) / wni) * xvec_it!!.get(0, idx + f) * FastMath.min(xi,
                                                sn.nservers.get(i, 0))))

                                        Xservice[i]!![k]!!.set(f,
                                            0,
                                            ((lambdaIKF * wi.get(0, k) / wni) * FastMath.min(xi,
                                                sn.nservers.get(i, 0)) * xvec_it!!.get(0, idx + f)))
                                    }
                                }

                                SchedStrategy.HOL, SchedStrategy.SJF, SchedStrategy.LJF, 
                                SchedStrategy.SEPT, SchedStrategy.LEPT -> {
                                    // Priority-based and job-size-based scheduling strategies
                                    // Use same logic as FCFS but with priority/size ordering (scheduling parameter dependent)
                                    var tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT)
                                    var rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / xi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    if (options.method == "default" || options.method == "closing") {
                                        result.TN.set(i,
                                            k,
                                            result.TN.get(i, k) + ((lambdaIKF * phiIKF / xi) * xvec_it!!.get(0,
                                                idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))
                                        Xservice[i]!![k]!!.set(f,
                                            0,
                                            lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it!!.get(0,
                                                idx + f))
                                    } else if (options.method == "statedep") {
                                        // Use scheduling parameters for priority-based strategies
                                        val priority = sn.schedparam.get(i, k)
                                        result.TN.set(i,
                                            k,
                                            result.TN.get(i, k) + (((lambdaIKF * phiIKF * priority) / wni) * 
                                                xvec_it!!.get(0, idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))

                                        Xservice[i]!![k]!!.set(f,
                                            0,
                                            ((lambdaIKF * priority / wni) * FastMath.min(xi,
                                                sn.nservers.get(i, 0)) * xvec_it!!.get(0, idx + f)))
                                    }
                                }

                                SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                                    // Generalized processor sharing and priority variants - use DPS-like logic
                                    val w = Matrix(1, K)
                                    var p = 0
                                    while (p < K) {
                                        w.set(0, k, sn.schedparam.get(i, k))
                                        p++
                                    }

                                    var tmpQ: Matrix? = Matrix(1, result.QN.numCols)
                                    Matrix.extractRows(result.QN, i, i + 1, tmpQ)
                                    tmpQ = tmpQ!!.transpose()
                                    var wxi = Matrix(1, 1) // number of jobs in the station
                                    wxi = w.mult(tmpQ, wxi)

                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + (((lambdaIKF * phiIKF * w.get(0,
                                            k)) / wxi.value()) * xvec_it!!.get(0, idx + f) * FastMath.min(xi,
                                            sn.nservers.get(i, 0))))

                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        ((lambdaIKF * w.get(0, k) / wxi.value()) * FastMath.min(xi,
                                            sn.nservers.get(i, 0)) * xvec_it!!.get(0, idx + f)))
                                }

                                SchedStrategy.LCFSPR -> {
                                    // Last Come First Served with Preemptive Resume - similar to LCFS but preemptive
                                    var tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT)
                                    var rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / xi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + ((lambdaIKF * phiIKF / xi) * xvec_it!!.get(0,
                                            idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))
                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it!!.get(0,
                                            idx + f))
                                }

                                SchedStrategy.POLLING, SchedStrategy.FORK, SchedStrategy.REF -> {
                                    // Specialized scheduling strategies - use basic FCFS logic as fallback
                                    var tmpTNT = Matrix(xvec_t!!.numRows, 1)
                                    Matrix.extractColumn(xvec_t, idx + f, tmpTNT)
                                    tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT)
                                    var rows = tmpTNT.numRows
                                    var row = 0
                                    while (row < rows) {
                                        result.TNt[i][k].set(row,
                                            0,
                                            result.TNt[i][k].get(row, 0) + (tmpTNT.get(row, 0) / xi_t.get(row,
                                                0) * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))))
                                        row++
                                    }

                                    result.TN.set(i,
                                        k,
                                        result.TN.get(i, k) + ((lambdaIKF * phiIKF / xi) * xvec_it!!.get(0,
                                            idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))))
                                    Xservice[i]!![k]!!.set(f,
                                        0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it!!.get(0,
                                            idx + f))
                                }

                                else -> 
                                    throw RuntimeException("Unsupported scheduling policy: ${sn.sched.get(sn.stations.get(i))}")
                            }
                        }
                    }
                }
            }
        }

        // Response Times - this is approximate, Little's law does not hold in transient
        result.RN = Matrix(M, K)
        for (i in 0..<M) {
            for (k in 0..<K) {
                if (result.TN.get(i, k) > 0) {
                    result.RN.set(i, k, result.QN.get(i, k) / result.TN.get(i, k))
                }
            }
        }

        // Utilisation
        result.UN = Matrix(M, K)
        for (i in 0..<M) {
            for (k in 0..<K) {
                var sumForUN = 0.0
                var sumForTN = 0.0
                val rows = Xservice[i]!![k]!!.numRows
                for (row in 0..<rows) {
                    if (Xservice[i]!![k]!!.get(row, 0) > 0) {
                        sumForUN += (Xservice[i]!![k]!!.get(row, 0) / lambda.get(sn.stations.get(i))!!
                            .get(sn.jobclasses.get(k))!!.get(row, 0))
                        sumForTN += Xservice[i]!![k]!!.get(row, 0)
                    }
                }
                result.UN.set(i, k, sumForUN)
                if ((sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) && (options.method == "statedep")) {
                    result.TN.set(i, k, sumForTN)
                }
            }
        }

        for (i in 0..<M) {
            for (k in 0..<K) {
                if (delayNodes.get(0, i) == 0.0) {
                    result.UN.set(i, k, result.UN.get(i, k) / sn.nservers.get(i, 0))
                }
            }
        }
        result.WN = Matrix(M, K) // TODO
        result.AN = Matrix(M, K) // TODO
    }

    override fun getXVecIt(): Matrix {
        return this.xvec_it!!
    }

    fun detectStiffnessUsingOstrowski(sn: NetworkStruct, rate: Matrix): Boolean {
        var transitionMatrix = sn.rt
        for (i in 0..<transitionMatrix.numRows) {
            var p = 0.0
            for (j in 0..<transitionMatrix.numCols) {
                if (i != j) p -= transitionMatrix.get(i, j)
            }
            transitionMatrix.set(i, i, p)
        }

        val expandRate = Matrix(1, rate.numCols * rate.numRows)
        for (i in 0..<rate.numCols * rate.numRows) {
            val r = i / rate.numCols
            val c = i % rate.numCols
            expandRate.set(0, i, rate.get(r, c))
        }
        transitionMatrix = transitionMatrix.elementMultWithVector(expandRate)
        val bound = Array<DoubleArray?>(transitionMatrix.numCols) { DoubleArray(2) }
        val n = transitionMatrix.numCols
        val alpha = 0.5
        var stiff = false

        for (i in 0..<n) {
            bound[i]!![0] = transitionMatrix.get(i, i)
            val rSum = transitionMatrix.sumAbsRows(i) - FastMath.abs(transitionMatrix.get(i, i))
            val cSum = transitionMatrix.sumAbsCols(i) - FastMath.abs(transitionMatrix.get(i, i))
            bound[i]!![1] = FastMath.pow(rSum, alpha) * FastMath.pow(cSum, 1 - alpha)
            stiff = stiff || (bound[i]!![0] < 0 && FastMath.abs(bound[i]!![1]) < FastMath.abs(bound[i]!![0]))
        }
        return stiff
    }
}
