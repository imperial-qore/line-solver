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
import jline.lang.constant.JobClassType
import jline.lang.constant.NodeType
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
        val M = sn.nstations // Number of stations
        val K = sn.nclasses // Number of classes

        val S = sn.nservers.copy()
        val initialPopulation = sn.njobs.elementSum()
        val SRows = S.numRows
        for (i in 0..<SRows) {
            if (Double.isInfinite(S.get(i, 0))) {
                S.set(i, 0, initialPopulation)
            }
        }

        // W (complete graph of transition rates) as per Ruuskanen et al., PEVA 151 (2021)
        var psi = Matrix(0, 0)
        var A = Matrix(0, 0)
        var B = Matrix(0, 0)
        val lambda = Matrix(M * K, 1)

        for (i in 0..<M) {
            if (sn.nodetype.get(sn.stationToNode.get(i).toInt()) == NodeType.Source) {
                for (k in 0..<K) {
                    if (sn.jobclasses.get(k).jobClassType == JobClassType.OPEN) {
                        lambda.set(i * K + k, 0, sn.rates.get(i, k))
                    }
                }
            }
        }

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

        // W (complete graph of transition rates) as per Ruuskanen et al., PEVA 151 (2021)
        val W = calculateW(sn, psi, A, B)

        // Handle degenerate case where W is empty (no valid transitions)
        if (W.numRows == 0 || W.numCols == 0) {
            // Return zero results for this degenerate model
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

        val ALambda = A.mult(lambda, null)

        val totalNumPhases = sn.phases.elementSum().toInt()
        // State mapping to queues (called Q(a) in Ruuskanen et al.)
        val Qa = Matrix(1, totalNumPhases)
        // To compute per-class queue length, utilisation and throughput at the end
        val SQC = Matrix(M * K, totalNumPhases)
        val SUC = Matrix(M * K, totalNumPhases)
        val STC = Matrix(M * K, totalNumPhases)
        // To compute total queue length in ODEs
        val SQ = Matrix(totalNumPhases, totalNumPhases)

        var state = 0
        for (i in 0..<M) {
            val station = sn.stations.get(i)
            for (r in 0..<K) {
                val jobClass = sn.jobclasses.get(r)
                val nPhases = sn.phases.get(i, r).toInt()
                for (k in 0..<nPhases) {
                    Qa.set(0, state, i)
                    SQC.set(i * K + r, state, 1.0)
                    SUC.set(i * K + r, state, 1.0 / S.get(i, 0))
                    STC.set(i * K + r, state, sn.proc.get(station)!!.get(jobClass)!!.get(1).sumRows(k))
                    state++
                }
            }
        }

        var nextSQRow = 0
        for (i in 0..<M) {
            for (r in 0..<K) {
                val nPhases = sn.phases.get(i, r).toInt()
                for (k in 0..<nPhases) {
                    for (col in 0..<totalNumPhases) {
                        if (Qa.get(0, col) == i.toDouble()) {
                            SQ.set(nextSQRow, col, 1) // Setting weights
                        }
                    }
                    nextSQRow++
                }
            }
        }

        val Sa = Matrix(totalNumPhases, 1)
        for (i in 0..<totalNumPhases) {
            Sa.set(i, 0, S.get(Qa.get(0, i).toInt(), 0))
        }

        var minNonZeroRate = Inf
        val WRows = W.numRows
        val WCols = W.numCols
        for (i in 0..<WRows) {
            for (j in 0..<WCols) {
                val tmpRate = FastMath.abs(W.get(i, j))
                if (tmpRate < minNonZeroRate && tmpRate > 0) {
                    minNonZeroRate = tmpRate
                }
            }
        }
        val tRange = doubleArrayOf(options.timespan[0],
            FastMath.min(options.timespan[1], FastMath.abs(10 * options.iter_max / minNonZeroRate)))

        val initSolLength = options.init_sol.length()
        val initialState = DoubleArray(initSolLength)
        val nextState = DoubleArray(initSolLength)
        for (i in 0..<initSolLength) {
            initialState[i] = options.init_sol.get(0, i)
            nextState[i] = 0.0
        }

        // Choose between original compact matrix form representation, and p-norm smoothed
        // representation as per Ruuskanen et al., PEVA 151 (2021).
        val ode: FirstOrderDifferentialEquations?

        if (options.config.pstar.size == 0) { // If pStar values do not exist
            ode = MatrixMethodODE(W, SQ, S, Qa, ALambda, initSolLength)
        } else {
            ode = MatrixMethodODE(W, SQ, S, Qa, ALambda, initSolLength, sn, options.config.pstar)
        }

        val Tmax: Int
        if (options.stiff) {
            val odeSolver: LSODA
            if (options.tol > GlobalConstants.CoarseTol) {
                odeSolver = options.odesolvers.fastStiffODESolver
            } else {
                odeSolver = options.odesolvers.accurateStiffODESolver
            }
            //System.out.print("Starting ODE integration cycle...");
            odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)
            //System.out.println("done.");
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
            val stepHandler = TransientDataHandler(initSolLength)
            odeSolver.addStepHandler(stepHandler)
            //System.out.print("Starting ODE integration cycle...");
            odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState)

            //System.out.println("done.");

            // Retrieve Transient Data
            result.t = stepHandler.tVec
            Tmax = result.t.numRows
            this.xvec_t = stepHandler.xVec
        }
        this.xvec_it = Matrix.extractRows(xvec_t, Tmax - 1, Tmax, null)

        // Use Transient Data to Store Results
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
            val theta = x.copy() // Theta per Ruuskanen et al., PEVA 151 (2021).
            val SQx = SQ.mult(x, null)
            for (phase in 0..<totalNumPhases) {
                val valSQx = SQx.get(phase, 0) + GlobalConstants.FineTol
                val valSa = Sa.get(phase, 0) + GlobalConstants.FineTol
                theta.set(phase, 0, x.get(phase, 0) / valSQx * FastMath.min(valSa, valSQx))
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
        result.WN = Matrix(M, K) // TODO
        result.AN = Matrix(M, K) // TODO
        for (i in 0..<M) {
            for (j in 0..<K) {
                result.QN.set(i, j, result.QNt[i][j].get(Tmax - 1, 0))
                result.UN.set(i, j, result.UNt[i][j].get(Tmax - 1, 0))
                result.TN.set(i, j, result.TNt[i][j].get(Tmax - 1, 0))
                result.RN.set(i, j, result.QN.get(i, j) / result.TN.get(i, j))
            }
        }
    }

    override fun getXVecIt(): Matrix? {
        return this.xvec_it
    }

    private fun calculateW(sn: NetworkStruct, psi: Matrix, A: Matrix?, B: Matrix?): Matrix {
        // ODE building as per Ruuskanen et al., PEVA 151 (2021).

        val calculateW = MatrixEquation()
        calculateW.alias(psi, "psi", A, "A", sn.rt, "P", B, "B")
        calculateW.process("W = psi + B*P*A'")
        var W = Matrix(calculateW.lookupSimple("W"))

        // Remove disabled transitions
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
