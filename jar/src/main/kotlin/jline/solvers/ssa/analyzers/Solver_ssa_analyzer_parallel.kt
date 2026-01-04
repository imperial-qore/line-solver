package jline.solvers.ssa.analyzers

import jline.api.mam.map_mean
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.StatefulNode
import jline.lang.state.EventCache
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAResult
import jline.solvers.ssa.SolverSSA
import jline.solvers.ssa.handlers.solver_ssa
import jline.util.Maths
import jline.util.RandomManager
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.TimeUnit

/**
 * Averages a collection of matrices by summing them and dividing by the count.
 * Uses Matrix.addEq() for efficient in-place addition.
 *
 * @param matrices Collection of matrices to average (must all have same dimensions)
 * @param numThreads Number of threads (used as divisor for averaging)
 * @return A new matrix containing the element-wise average
 */
private fun averageMatrices(matrices: Collection<Matrix>, numThreads: Int): Matrix {
    if (matrices.isEmpty()) {
        throw IllegalArgumentException("Cannot average empty collection of matrices")
    }
    val first = matrices.first()
    val result = Matrix(first.numRows, first.numCols)
    result.zero()
    for (m in matrices) {
        result.addEq(m)
    }
    return Matrix.scaleMult(result, 1.0 / numThreads.toDouble())
}

fun solver_ssa_analyzer_parallel(sn: NetworkStruct,
                                 init_state: MutableMap<StatefulNode?, Matrix?>,
                                 options: SolverOptions,
                                 solverSSA: SolverSSA): SSAResult {
    val M = sn.nstations
    val K = sn.nclasses
    val PH = sn.proc
    val S = sn.nservers
    val NK = sn.njobs.transpose()

    val QNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()
    val UNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()
    val RNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()
    val TNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()
    val CNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()
    val XNs: MutableMap<Int, Matrix> = ConcurrentHashMap<Int, Matrix>()


    options.samples = FastMath.ceil(options.samples / solverSSA.numThreads.toDouble()).toInt()
    solverSSA.eventCache = EventCache(true, options.config.eventcache)

    for (t in 0..<solverSSA.numThreads) {
        val threadIndex = t
        solverSSA.threadPool.submit(Runnable {
            // Get thread-specific random generator for reproducible parallel execution
            val threadRandom = RandomManager.getParallelRandom("SSA_PARALLEL", threadIndex)
            var probSysState = Matrix(0, 0)
            var SSq: Matrix? = Matrix(0, 0)
            var arvRates: MutableMap<Int?, Matrix> = HashMap<Int?, Matrix>()
            var depRates: MutableMap<Int?, Matrix> = HashMap<Int?, Matrix>()

            when (options.method) {
                "ssa.parallel", "parallel" -> {
                    val result = solver_ssa(sn, solverSSA.eventCache, init_state, options, solverSSA)
                    probSysState = result.pi
                    SSq = result.SSq
                    arvRates = result.arvRates
                    depRates = result.depRates
                }
            }
            val XN = Matrix(1, K)
            XN.fill(Double.NaN)
            val UN = Matrix(M, K)
            UN.fill(Double.NaN)
            val QN = Matrix(M, K)
            QN.fill(Double.NaN)
            val RN = Matrix(M, K)
            RN.fill(Double.NaN)
            val TN = Matrix(M, K)
            TN.fill(Double.NaN)
            val CN = Matrix(1, K)
            CN.fill(Double.NaN)

            for (k in 0..<K) {
                val refsf = sn.stationToStateful.get(sn.refstat.get(k).toInt())
                val departure = depRates.get(k)
                val dep_wset_refsf = Matrix.extractColumn(departure, refsf.toInt(), null)
                XN.set(k, probSysState.mult(dep_wset_refsf).toDouble())
                for (i in 0..<M) {
                    val isf = sn.stationToStateful.get(i).toInt()

                    val dep_isf = Matrix.extractColumn(departure, isf, null)
                    TN.set(i, k, probSysState.mult(dep_isf).toDouble())

                    val ssq_extracted = Matrix.extractColumn(SSq, i * K + k, null)
                    QN.set(i, k, probSysState.mult(ssq_extracted).toDouble())

                    when (sn.sched.get(sn.stations.get(i))) {
                        SchedStrategy.INF -> UN.set(i, k, QN.get(i, k))
                        else -> if (!PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.isEmpty) {
                            // probSysState*arvRates(:,i,k)*map_mean(PH{i}{k})/S(i);
                            val arrival = arvRates.get(k)
                            val arv_ik = Matrix.extractColumn(arrival, i, null)
                            val map_mean = map_mean(PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(0),
                                PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(1)) / S.get(i)
                            UN.set(i, k, probSysState.mult(arv_ik).toDouble() * map_mean)
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
            QN.apply(Double.NaN, 0.0, "equal")
            CN.apply(Double.NaN, 0.0, "equal")
            RN.apply(Double.NaN, 0.0, "equal")
            UN.apply(Double.NaN, 0.0, "equal")
            XN.apply(Double.NaN, 0.0, "equal")
            TN.apply(Double.NaN, 0.0, "equal")

            QNs.put(threadIndex, QN)
            CNs.put(threadIndex, CN)
            RNs.put(threadIndex, RN)
            UNs.put(threadIndex, UN)
            XNs.put(threadIndex, XN)
            TNs.put(threadIndex, TN)
        })
        //            threads[t].start();
    }
    solverSSA.threadPool.shutdown()
    try {
        solverSSA.threadPool.awaitTermination(600, TimeUnit.SECONDS)
    } catch (e: InterruptedException) {
        throw RuntimeException(e)
    }


    // Average results from all parallel threads using optimized addEq() method
    val XN = averageMatrices(XNs.values, solverSSA.numThreads)
    val UN = averageMatrices(UNs.values, solverSSA.numThreads)
    val QN = averageMatrices(QNs.values, solverSSA.numThreads)
    val RN = averageMatrices(RNs.values, solverSSA.numThreads)
    val TN = averageMatrices(TNs.values, solverSSA.numThreads)
    val CN = averageMatrices(CNs.values, solverSSA.numThreads)

    val tranSysState: MutableMap<Int?, Matrix?>? = null
    val tranSync: Matrix? = null
    val res = SSAResult(QN, UN, RN, TN, CN, XN, tranSysState, tranSync, sn)
    res.method = "parallel"
    return res
}