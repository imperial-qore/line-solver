/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers

import jline.io.line_warning
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.solvers.SolverOptions
import jline.solvers.SolverResult
import jline.solvers.fluid.handlers.MatrixMethodODE
import jline.solvers.fluid.handlers.TransientDataHandler
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixEquation
import odesolver.LSODA
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import java.util.*
import kotlin.Any
import kotlin.Array
import kotlin.DoubleArray
import kotlin.Int
import kotlin.arrayOfNulls
import kotlin.doubleArrayOf
import kotlin.math.max

class MatrixMethodAnalyzer : FluidAnalyzer {
    var xvec_t: Matrix? = null
    var xvec_it: Matrix? = null

    override fun analyze(sn: NetworkStruct, options: SolverOptions, result: SolverResult) {
        val M = sn.nstations
        val K = sn.nclasses

        val S = sn.nservers.copy()
        val initialPopulation = sn.njobs.elementSum()
        val SRows = S.numRows
        for (i in 0..<SRows) {
            if (Double.isInfinite(S.get(i, 0))) {
                S.set(i, 0, initialPopulation)
            }
        }

        // Extract station-to-station routing matrix from stateful-to-stateful sn.rt
        val stationIndices = ArrayList<Int>()
        for (ist in 0..<M) {
            val isf = sn.stationToStateful.get(ist).toInt()
            for (r in 0..<K) {
                stationIndices.add(isf * K + r)
            }
        }
        val PSize = M * K
        val P = Matrix(PSize, PSize)
        for (row in 0..<PSize) {
            for (col in 0..<PSize) {
                P.set(row, col, sn.rt.get(stationIndices[row], stationIndices[col]))
            }
        }

        // Remove Sink->Source feedback routing for open classes
        for (srcIst in 0..<M) {
            if (sn.sched.get(sn.stations.get(srcIst)) == SchedStrategy.EXT) {
                for (r in 0..<K) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        val srcCol = srcIst * K + r
                        for (fromIst in 0..<M) {
                            if (fromIst != srcIst) {
                                for (fromR in 0..<K) {
                                    P.set(fromIst * K + fromR, srcCol, 0.0)
                                }
                            }
                        }
                    }
                }
            }
        }

        // Build block-diagonal matrices for ODE
        var psi = Matrix(0, 0)
        var A = Matrix(0, 0)
        var B = Matrix(0, 0)

        for (i in 0..<M) {
            val station = sn.stations.get(i)
            for (r in 0..<K) {
                val jobClass = sn.jobclasses.get(r)
                if (sn.phases.get(i, r) == 0.0) {
                    val zeroMatrix = Matrix(1, 1, 1)
                    val nanMatrix = Matrix(1, 1, 1)
                    nanMatrix.set(0, 0, kotlin.Double.NaN)
                    psi = psi.createBlockDiagonal(zeroMatrix)
                    A = A.createBlockDiagonal(nanMatrix)
                    B = B.createBlockDiagonal(zeroMatrix)
                } else {
                    psi = psi.createBlockDiagonal(sn.proc.get(station)!!.get(jobClass)!!.get(0))
                    A = A.createBlockDiagonal(sn.pie.get(station)!!.get(jobClass)!!.transpose())
                    B = B.createBlockDiagonal(sn.proc.get(station)!!.get(jobClass)!!.get(1).sumRows())
                }
            }
        }

        // W = psi + B*P*A' using station-level P
        val W = calculateW(psi, A, B, P)

        // Handle degenerate case where W is empty
        if (W.numRows == 0 || W.numCols == 0) {
            result.QN = Matrix(M, K)
            result.UN = Matrix(M, K)
            result.RN = Matrix(M, K)
            result.TN = Matrix(M, K)
            result.WN = Matrix(M, K)
            result.AN = Matrix(M, K)
            result.t = Matrix(1, 1)
            result.QNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
            result.UNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
            result.TNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
            for (i in 0..<M) {
                for (k in 0..<K) {
                    result.QNt[i][k] = Matrix(1, 1)
                    result.UNt[i][k] = Matrix(1, 1)
                    result.TNt[i][k] = Matrix(1, 1)
                }
            }
            this.xvec_t = Matrix(1, options.init_sol.length())
            this.xvec_it = Matrix(1, options.init_sol.length())
            return
        }

        // Build arrival rate vector ALambda for mixed/open networks
        val sourceArrivals = Matrix(M, K)
        for (srcIst in 0..<M) {
            if (sn.sched.get(sn.stations.get(srcIst)) == SchedStrategy.EXT) {
                for (r in 0..<K) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        sourceArrivals.set(srcIst, r, sn.rates.get(srcIst, r))
                    }
                }
            }
        }

        // Total states including placeholders for disabled classes
        var totalStates = 0
        for (i in 0..<M) {
            for (r in 0..<K) {
                val np = sn.phases.get(i, r).toInt()
                totalStates += if (np == 0) 1 else np
            }
        }

        // Build ALambda_full matching full W matrix structure
        val ALambdaFull = Matrix(totalStates, 1)
        var state = 0
        for (ist in 0..<M) {
            val station = sn.stations.get(ist)
            for (r in 0..<K) {
                val jobClass = sn.jobclasses.get(r)
                val nPhases = sn.phases.get(ist, r).toInt()
                if (nPhases > 0) {
                    if (sn.sched.get(station) == SchedStrategy.EXT) {
                        state += nPhases
                    } else {
                        var arrivalRateToQueue = 0.0
                        for (srcIst in 0..<M) {
                            if (sourceArrivals.get(srcIst, r) > 0) {
                                arrivalRateToQueue += sourceArrivals.get(srcIst, r) * P.get(srcIst * K + r, ist * K + r)
                            }
                        }
                        if (arrivalRateToQueue > 0) {
                            val pieVec = sn.pie.get(station)!!.get(jobClass)!!
                            for (k in 0..<nPhases) {
                                ALambdaFull.set(state, 0, pieVec.get(k) * arrivalRateToQueue)
                                state++
                            }
                        } else {
                            state += nPhases
                        }
                    }
                } else {
                    state++
                }
            }
        }

        // Filter disabled transitions (NaN columns in W)
        val keep = ArrayList<Int>()
        val sumWCols = W.sumCols()
        for (col in 0..<W.numCols) {
            if (!Double.isNaN(sumWCols.get(0, col))) {
                keep.add(col)
            }
        }

        val WFiltered = Matrix(keep.size, keep.size)
        for (i in 0..<keep.size) {
            for (j in 0..<keep.size) {
                WFiltered.set(i, j, W.get(keep[i], keep[j]))
            }
        }

        val ALambda = Matrix(keep.size, 1)
        for (i in 0..<keep.size) {
            ALambda.set(i, 0, ALambdaFull.get(keep[i], 0))
        }

        // Build Qa, SQC, SUC, STC, x0_build with full structure then filter
        val QaFull = Matrix(1, totalStates)
        val SQCFull = Matrix(M * K, totalStates)
        val SUCFull = Matrix(M * K, totalStates)
        val STCFull = Matrix(M * K, totalStates)
        val x0Build = Matrix(totalStates, 1)

        state = 0
        var initSolIdx = 0
        for (ist in 0..<M) {
            val station = sn.stations.get(ist)
            for (r in 0..<K) {
                val jobClass = sn.jobclasses.get(r)
                val nPhases = sn.phases.get(ist, r).toInt()
                if (nPhases == 0) {
                    QaFull.set(0, state, ist.toDouble())
                    state++
                } else {
                    for (k in 0..<nPhases) {
                        QaFull.set(0, state, ist.toDouble())
                        SQCFull.set(ist * K + r, state, 1.0)
                        SUCFull.set(ist * K + r, state, 1.0 / S.get(ist, 0))
                        STCFull.set(ist * K + r, state, sn.proc.get(station)!!.get(jobClass)!!.get(1).sumRows(k))
                        x0Build.set(state, 0, options.init_sol.get(0, initSolIdx))
                        initSolIdx++
                        state++
                    }
                }
            }
        }

        // Apply keep filtering
        val nStatesFiltered = keep.size
        val Qa = Matrix(1, nStatesFiltered)
        val SQC = Matrix(M * K, nStatesFiltered)
        val SUC = Matrix(M * K, nStatesFiltered)
        val STC = Matrix(M * K, nStatesFiltered)
        val x0 = DoubleArray(nStatesFiltered)

        for (i in 0..<nStatesFiltered) {
            val ki = keep[i]
            Qa.set(0, i, QaFull.get(0, ki))
            for (row in 0..<M * K) {
                SQC.set(row, i, SQCFull.get(row, ki))
                SUC.set(row, i, SUCFull.get(row, ki))
                STC.set(row, i, STCFull.get(row, ki))
            }
            x0[i] = x0Build.get(ki, 0)
        }

        // Identify Source station states and zero their x0
        val isSourceState = BooleanArray(nStatesFiltered)
        for (s in 0..<nStatesFiltered) {
            val ist = Qa.get(0, s).toInt()
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                isSourceState[s] = true
                x0[s] = 0.0
            }
        }

        // Build SQ matrix
        val SQ = Matrix(nStatesFiltered, nStatesFiltered)
        for (s in 0..<nStatesFiltered) {
            val ist = Qa.get(0, s).toInt()
            for (col in 0..<nStatesFiltered) {
                if (Qa.get(0, col) == ist.toDouble()) {
                    SQ.set(s, col, 1)
                }
            }
        }

        val Sa = Matrix(nStatesFiltered, 1)
        for (i in 0..<nStatesFiltered) {
            Sa.set(i, 0, S.get(Qa.get(0, i).toInt(), 0))
        }

        var minNonZeroRate = Inf
        for (i in 0..<WFiltered.numRows) {
            for (j in 0..<WFiltered.numCols) {
                val tmpRate = FastMath.abs(WFiltered.get(i, j))
                if (tmpRate < minNonZeroRate && tmpRate > 0) {
                    minNonZeroRate = tmpRate
                }
            }
        }
        val tRange = doubleArrayOf(options.timespan[0],
            FastMath.min(options.timespan[1], FastMath.abs(10 * options.iter_max / minNonZeroRate)))

        val initialState = x0.copyOf()
        val nextState = DoubleArray(nStatesFiltered)

        val ode: FirstOrderDifferentialEquations?
        if (options.config.pstar.size == 0) {
            ode = MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState)
        } else {
            ode = MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState, sn, options.config.pstar)
        }

        val Tmax: Int
        if (options.stiff) {
            val odeSolver: LSODA
            if (options.tol > GlobalConstants.CoarseTol) {
                odeSolver = options.odesolvers.fastStiffODESolver
            } else {
                odeSolver = options.odesolvers.accurateStiffODESolver
            }
            odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)
            Tmax = odeSolver.stepsTaken + 1
            val tVec = Matrix(Tmax, 1)
            val xVec = Matrix(Tmax, ode.dimension)
            val tHistory = odeSolver.tvec
            val yHistory = odeSolver.yvec
            for (i in 0..<Tmax) {
                tVec.set(i, 0, tHistory.get(i)!!)
                for (j in 0..<ode.dimension) xVec.set(i, j, max(0.0, yHistory.get(i)!![j]!!))
            }
            result.t = tVec
            this.xvec_t = xVec
        } else {
            if (options.tol > GlobalConstants.CoarseTol && (options.verbose == VerboseLevel.DEBUG)) {
                line_warning(mfilename(object : Any() {}),
                    "Fast, non-stiff ODE solver is not yet available in JLINE. Using accurate non-stiff ODE solver instead.")
            }
            val odeSolver = options.odesolvers.accurateODESolver
            odeSolver.clearStepHandlers()
            val stepHandler = TransientDataHandler(nStatesFiltered)
            odeSolver.addStepHandler(stepHandler)
            odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)
            result.t = stepHandler.tVec
            Tmax = result.t.numRows
            this.xvec_t = stepHandler.xVec
        }
        this.xvec_it = Matrix.extractRows(xvec_t, Tmax - 1, Tmax, null)

        result.QNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        result.UNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        result.TNt = Array<Array<Matrix?>?>(M) { arrayOfNulls<Matrix>(K) }
        for (i in 0..<M) {
            for (k in 0..<K) {
                result.QNt[i][k] = Matrix(Tmax, 1)
                result.UNt[i][k] = Matrix(Tmax, 1)
                result.TNt[i][k] = Matrix(Tmax, 1)
            }
        }

        for (step in 0..<Tmax) {
            var x = Matrix.extractRows(xvec_t, step, step + 1, null)
            x = x.transpose()
            val theta = x.copy()
            val SQx = SQ.mult(x, null)
            for (phase in 0..<nStatesFiltered) {
                if (isSourceState[phase]) {
                    theta.set(phase, 0, 0.0)
                } else {
                    val valSQx = SQx.get(phase, 0) + GlobalConstants.FineTol
                    val valSa = Sa.get(phase, 0)
                    theta.set(phase, 0, x.get(phase, 0) / valSQx * FastMath.min(valSa, valSQx))
                }
            }

            val QNtmp = SQC.mult(x, null)
            val UNtmp = SUC.mult(theta, null)
            val TNtmp = STC.mult(theta, null)

            for (i in 0..<M) {
                for (k in 0..<K) {
                    result.QNt[i][k].set(step, 0, QNtmp.get(i * K + k, 0))
                    result.UNt[i][k].set(step, 0, UNtmp.get(i * K + k, 0))
                    result.TNt[i][k].set(step, 0, TNtmp.get(i * K + k, 0))
                }
            }
        }

        result.QN = Matrix(M, K)
        result.UN = Matrix(M, K)
        result.RN = Matrix(M, K)
        result.TN = Matrix(M, K)
        result.WN = Matrix(M, K)
        result.AN = Matrix(M, K)
        for (i in 0..<M) {
            for (j in 0..<K) {
                result.QN.set(i, j, result.QNt[i][j].get(Tmax - 1, 0))
                result.UN.set(i, j, result.UNt[i][j].get(Tmax - 1, 0))
                result.TN.set(i, j, result.TNt[i][j].get(Tmax - 1, 0))
                result.RN.set(i, j, result.QN.get(i, j) / result.TN.get(i, j))
            }
        }

        // Set throughput at Source stations to arrival rates for open classes
        for (ist in 0..<M) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                for (r in 0..<K) {
                    if (!Double.isNaN(sn.rates.get(ist, r)) && sn.rates.get(ist, r) > 0) {
                        result.TN.set(ist, r, sn.rates.get(ist, r))
                    }
                }
            }
        }
    }

    override fun getXVecIt(): Matrix? {
        return this.xvec_it
    }

    private fun calculateW(psi: Matrix, A: Matrix?, B: Matrix?, P: Matrix): Matrix {
        val calculateW = MatrixEquation()
        calculateW.alias(psi, "psi", A, "A", P, "P", B, "B")
        calculateW.process("W = psi + B*P*A'")
        var W = Matrix(calculateW.lookupSimple("W"))

        if (W.hasNaN()) {
            val WRows = W.numRows
            val WCols = W.numCols
            val colsWithNans: MutableList<Int?> = LinkedList<Int?>()
            val sumWCols = W.sumCols()
            for (i in 0..<WCols) {
                if (Double.isNaN(sumWCols.get(0, i))) {
                    colsWithNans.add(i)
                }
            }

            val tmpW = Matrix(W.numRows - colsWithNans.size, W.numCols - colsWithNans.size)
            var tmpWRow = 0
            var i = 0
            while (i < WRows) {
                while (i < WRows && colsWithNans.contains(i)) {
                    i++
                }
                if (i >= WRows) break
                var tmpWCol = 0
                var j = 0
                while (j < WCols) {
                    while (j < WCols && colsWithNans.contains(j)) {
                        j++
                    }
                    if (j >= WCols) break
                    tmpW.set(tmpWRow, tmpWCol, W.get(i, j))
                    tmpWCol++
                    j++
                }
                tmpWRow++
                i++
            }

            W = tmpW
        }

        return W
    }
}
