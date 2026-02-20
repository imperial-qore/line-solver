package jline.solvers.ctmc.analyzers

import jline.api.mam.map_mean
import jline.api.mc.ctmc_solve
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.lang.nodeparam.CacheNodeParam
import jline.lang.state.ToMarginal
import jline.solvers.SolverOptions
import jline.solvers.ctmc.SolverCTMC
import jline.solvers.ctmc.handlers.Solver_ctmc.Companion.solver_ctmc
import jline.util.Maths
import jline.util.MatFileUtils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.sn.snNonmarkovToPh

class Solver_ctmc_analyzer(private val solverCTMC: SolverCTMC?) {
    companion object {
        @JvmStatic
        fun solver_ctmc_analyzer(snInput: NetworkStruct, options: SolverOptions): SolverCTMC.AnalyzerResult {
            // Convert non-Markovian distributions to PH
            var sn = snNonmarkovToPh(snInput, options)
            val M = sn.nstations
            val K = sn.nclasses
            val S = sn.nservers
            //    NK -> (int) k * 1
            val NK = sn.njobs.transpose()
            val schedid = sn.sched
            val Tstart = System.nanoTime()
            val PH = sn.proc

            // Note: hide_immediate now selectively preserves Cache immediate transitions
            // in solver_ctmc, so we no longer need to disable it entirely for Cache nodes
            val solverCTMCResult = solver_ctmc(sn, options)
            val InfGen = solverCTMCResult.getQ()
            val StateSpace = solverCTMCResult.getStateSpace()
            val StateSpaceAggr = solverCTMCResult.getStateSpaceAggr()
            val EventFiltration = solverCTMCResult.getDfilt()
            val arvRates = solverCTMCResult.getArvRates()
            val depRates = solverCTMCResult.getDepRates()
            sn = solverCTMCResult.getSn()

            for (isf in 0..<sn.nstateful) {
                if (sn.state.get(sn.stateful.get(isf))!!.getNumCols() < sn.space.get(sn.stateful.get(isf))!!
                        .getNumCols()) {
                    val state_matrix = Matrix(1, sn.space.get(sn.stateful.get(isf))!!.getNumCols())
                    state_matrix.zero()

                    val startIdx =
                        (sn.space.get(sn.stateful.get(isf))!!.getNumCols() - sn.state.get(sn.stateful.get(isf))!!
                            .getNumCols())
                    val endIdx = state_matrix.getNumCols()
                    for (col in startIdx..<endIdx) {
                        state_matrix.set(0, col, sn.state.get(sn.stateful.get(isf))!!.get(col - startIdx))
                    }
                    sn.state.replace(sn.stateful.get(isf), state_matrix)
                }
            }

            val sncopy = sn
            var fname = ""
            if (options.keep) {
                try {
                    MatFileUtils.ensureWorkspaceDirectoryExists()
                    fname = MatFileUtils.genFilename("workspace")
                    MatFileUtils.saveCTMCWorkspace(StateSpace, InfGen, null, fname)
                } catch (e: Exception) {
                    line_warning("solver_ctmc_analyzer", "Could not save workspace to .mat file: %s", e.message)
                    fname = ""
                }
            } else {
                fname = ""
            }

            val wset = Matrix(1, InfGen.length())
            for (col in 0..<wset.getNumCols()) {
                wset.set(0, col, col)
            }

            val probSysState = ctmc_solve(InfGen)

            // Check if CTMC solve failed and returned NaN values
            if (probSysState.hasNaN() || probSysState.isEmpty()) {
                throw RuntimeException("CTMC solver failed to compute steady-state probabilities for this cache model. " +
                    "This may indicate numerical instability or an invalid model configuration.")
            }

            for (row in 0..<probSysState.getNumRows()) {
                for (col in 0..<probSysState.getNumCols()) {
                    if (probSysState.get(row, col) < 0) {
                        probSysState.set(row, col, 0)
                    }
                }
            }

            val sum = probSysState.sumSubMatrix(0, probSysState.getNumRows(), 0, probSysState.getNumCols())
            if (sum > 0) {
                probSysState.divide(sum, probSysState, true)
            } else {
                throw RuntimeException("CTMC solver computed zero total probability. This indicates an invalid model configuration.")
            }

            val XN = Matrix(1, K)
            XN.zero()
            val UN = Matrix(M, K)
            UN.zero()
            val QN = Matrix(M, K)
            QN.zero()
            val RN = Matrix(M, K)
            RN.zero()
            val TN = Matrix(M, K)
            TN.zero()
            val CN = Matrix(1, K)
            CN.zero()

            val istSpaceShift = Matrix(1, M)
            istSpaceShift.zero()

            for (i in 0..<M) {
                if (i == 0) {
                    istSpaceShift.set(0, i, 0)
                } else {
                    val temp = istSpaceShift.get(0, i - 1) + sn.space.get(sn.stateful.get(i - 1))!!.getNumCols()
                    istSpaceShift.set(0, i, temp)
                }
            }

            var refsf = 0.0
            for (k in 0..<K) {
                refsf = sn.stationToStateful.get(sn.refstat.get(k).toInt())
                XN.set(0, k, refsf)
                var sumValue = 0.0
                for (i in 0..<wset.getNumCols()) {
                    val index = wset.get(i).toInt()
                    sumValue += probSysState.get(index) * arvRates[index]!![refsf.toInt()]!![k]
                }
                XN.set(0, k, sumValue)
            }

            for (i in 0..<M) {
                val isf = sn.stationToStateful.get(i).toInt()
                var ind = sn.stationToNode.get(i).toInt()
                for (k in 0..<K) {
                    var sumTN = 0.0
                    var sumQN = 0.0
                    for (index in 0..<wset.getNumCols()) {
                        val wstIdx = wset.get(index).toInt()
                        val depRate = depRates[wstIdx]!![isf]!![k]
                        val probState = probSysState.get(wstIdx)
                        sumTN += probState * depRate
                        val ssaValue = StateSpaceAggr.get(wstIdx, i * K + k)
                        val product = probState * ssaValue
                        sumQN += product
                    }
                    TN.set(i, k, sumTN)
                    QN.set(i, k, sumQN)
                }
                if (sn.nodetype.get(ind) != NodeType.Source) {
                    when (schedid.get(sn.stations.get(i))) {
                        SchedStrategy.INF -> {
                            var k = 0
                            while (k < K) {
                                UN.set(i, k, QN.get(i, k))
                                k++
                            }
                        }
                        SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS -> if (sn.lldscaling.isEmpty() && sn.cdscaling.isEmpty()) {
                            var k = 0
                            while (k < K) {
                                if (!PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.isEmpty()) {
                                    val value: MatrixCell = PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!
                                    val mean = map_mean(value) / S.get(i)
                                    var UNarv_ik = 0.0
                                    var idx = 0
                                    while (idx < wset.length()) {
                                        val wsetIdx = wset.get(idx).toInt()
                                        UNarv_ik += probSysState.get(idx) * arvRates[wsetIdx]!![isf]!![k]
                                        idx++
                                    }
                                    UNarv_ik = UNarv_ik * mean
                                    val UNdep_ik = TN.get(i, k) * mean
                                    UN.set(i, k, Maths.max(UNarv_ik, UNdep_ik))
                                }
                                k++
                            }
                        } else {
                            ind = sn.stationToNode.get(i).toInt()
                            var col = 0
                            while (col < K) {
                                UN.set(i, col, 0)
                                col++
                            }
                            var index = 0
                            while (index < wset.getNumCols()) {
                                val st = wset.get(index).toInt()
                                val StateSpaceColStart = istSpaceShift.get(i).toInt()
                                val StateSpaceColEnd =
                                    istSpaceShift.get(i).toInt() + sn.space.get(sn.stateful.get(i))!!.getNumCols()
                                val toMarginalResult = ToMarginal.toMarginal(sn,
                                    ind,
                                    Matrix.extract(StateSpace, st, st + 1, StateSpaceColStart, StateSpaceColEnd),
                                    null,
                                    null,
                                    null,
                                    null,
                                    null)
                                val ni = toMarginalResult.ni
                                val nir = toMarginalResult.nir
                                var checkNi = true
                                var niIndex = 0
                                while (niIndex < ni.length()) {
                                    if (ni.get(niIndex) <= 0) {
                                        checkNi = false
                                    }
                                    niIndex++
                                }
                                if (checkNi) {
                                    var k = 0
                                    while (k < K) {
                                        val v = probSysState.get(st) * nir.get(k) * sn.schedparam.get(i, k)
                                        val dividend = nir.mult(sn.schedparam.getRow(i).transpose())
                                        var addSum = 0.0
                                        var divIndex = 0
                                        while (divIndex < dividend.length()) {
                                            addSum += v / dividend.get(divIndex)
                                            divIndex++
                                        }
                                        UN.set(i, k, addSum + UN.get(i, k))
                                        k++
                                    }
                                }
                                index++
                            }
                        }
                        else -> if ((sn.lldscaling==null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.isEmpty())) {
                            var k = 0
                            while (k < K) {
                                if (!PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.isEmpty()) {
                                    var UNarv_ik = 0.0
                                    val value: MatrixCell = PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!
                                    val mean = map_mean(value)
                                    var idx = 0
                                    while (idx < wset.length()) {
                                        val wsetIdx = wset.get(idx).toInt()
                                        UNarv_ik += probSysState.get(idx) * arvRates[wsetIdx]!![isf]!![k]
                                        idx++
                                    }
                                    UNarv_ik = UNarv_ik * mean / S.get(i)
                                    val UNdep_ik = TN.get(i, k) * map_mean(value) / S.get(i)
                                    UN.set(i, k, Maths.max(UNarv_ik, UNdep_ik))
                                }
                                k++
                            }
                        } else {
                            ind = sn.stationToNode.get(i).toInt()
                            var col = 0
                            while (col < K) {
                                UN.set(i, col, 0)
                                col++
                            }

                            var index = 0
                            while (index < wset.length()) {
                                val st = wset.get(index).toInt()
                                val StateSpaceColStart = istSpaceShift.get(i).toInt()
                                val StateSpaceColEnd =
                                    istSpaceShift.get(i).toInt() + sn.space.get(sn.stateful.get(i))!!.getNumCols()

                                val toMarginalResult = ToMarginal.toMarginal(sn,
                                    ind,
                                    Matrix.extract(StateSpace, st, st + 1, StateSpaceColStart, StateSpaceColEnd),
                                    null,
                                    null,
                                    null,
                                    null,
                                    null)
                                val ni = toMarginalResult.ni
                                val sir = toMarginalResult.sir
                                var checkNir = true
                                var niIndex = 0
                                while (niIndex < ni.length()) {
                                    if (ni.get(niIndex) <= 0) {
                                        checkNir = false
                                    }
                                    niIndex++
                                }
                                if (checkNir) {
                                    var k = 0
                                    while (k < K) {
                                        val v = UN.get(i, k) + probSysState.get(st) * sir.get(k) / S.get(i)
                                        UN.set(i, k, v)
                                        k++
                                    }
                                }
                                index++
                            }
                        }
                    }
                }
            }

            for (k in 0..<K) {
                for (i in 0..<M) {
                    if (TN.get(i, k) > 0) {
                        RN.set(i, k, QN.get(i, k) / TN.get(i, k))
                    } else {
                        RN.set(i, k, 0)
                    }
                }
                CN.set(k, NK.get(k) / XN.get(k))
            }
            QN.setNaNToZero()
            CN.setNaNToZero()
            RN.setNaNToZero()
            UN.setNaNToZero()
            XN.setNaNToZero()
            TN.setNaNToZero()

            val Tstop = System.nanoTime()
            val runtime = (Tstop - Tstart).toDouble() / 1000000000.0

            val TNcache = Matrix(sn.nstateful, K)
            for (k in 0..<K) {
                for (isf in 0..<sn.nstateful) {
                    val ind = sncopy.statefulToNode.get(isf).toInt()
                    if (sncopy.nodetype.get(ind) == NodeType.Cache) {
                        var TNcacheSum = 0.0
                        for (i in 0..<wset.getNumCols()) {
                            val index = wset.get(i).toInt()
                            TNcacheSum += probSysState.get(i) * depRates[index]!![isf]!![k]
                        }
                        TNcache.set(isf, k, TNcacheSum)
                    }
                }
            }

            for (k in 0..<K) {
                for (isf in 0..<sncopy.nstateful) {
                    val ind = sncopy.statefulToNode.get(isf).toInt()
                    if (sncopy.nodetype.get(ind) == NodeType.Cache) {
                        val cacheParam = sncopy.nodeparam.get(sn.nodes.get(ind)) as CacheNodeParam
                        if (cacheParam.hitclass.length() > k) {
                            val h = cacheParam.hitclass.get(k).toInt()
                            val m = cacheParam.missclass.get(k).toInt()
                            if (h == -1 || m == -1) {
                                // Hit or miss class not defined for this class - set NaN
                                if (cacheParam.actualhitprob == null) {
                                    cacheParam.actualhitprob = Matrix(1, K)
                                    cacheParam.actualhitprob.fill(Double.NaN)
                                }
                                cacheParam.actualhitprob.set(k, Double.NaN)
                                if (cacheParam.actualmissprob == null) {
                                    cacheParam.actualmissprob = Matrix(1, K)
                                    cacheParam.actualmissprob.fill(Double.NaN)
                                }
                                cacheParam.actualmissprob.set(k, Double.NaN)
                            } else {
                                val actualhitprobValue =
                                    TNcache.get(isf, h) / (TNcache.get(isf, h) + TNcache.get(isf, m))
                                val actualmissprobValue =
                                    TNcache.get(isf, m) / (TNcache.get(isf, h) + TNcache.get(isf, m))
                                if (cacheParam.actualhitprob == null) {
                                    cacheParam.actualhitprob = Matrix(1, K)
                                    cacheParam.actualhitprob.fill(Double.NaN)
                                }
                                cacheParam.actualhitprob.set(k, actualhitprobValue)
                                if (cacheParam.actualmissprob == null) {
                                    cacheParam.actualmissprob = Matrix(1, K)
                                    cacheParam.actualmissprob.fill(Double.NaN)
                                }
                                cacheParam.actualmissprob.set(k, actualmissprobValue)
                            }
                        }
                    }
                }
            }
            return SolverCTMC.AnalyzerResult(QN,
                UN,
                RN,
                TN,
                CN,
                XN,
                InfGen,
                StateSpace,
                StateSpaceAggr,
                EventFiltration,
                runtime,
                fname,
                sncopy)
        }
    }
}