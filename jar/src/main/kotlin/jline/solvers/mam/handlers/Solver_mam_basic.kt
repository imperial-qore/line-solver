package jline.solvers.mam.handlers

import jline.api.mam.*
import jline.api.qsys.qsys_mapdc
import jline.api.sn.snGetDemandsChain
import jline.lang.JobClass
import jline.lang.NetworkStruct
import jline.lang.nodes.Station
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.lang.constant.ProcessType
import jline.lang.constant.SchedStrategy
import jline.lib.butools.MMAPPH1FCFS
import jline.lib.butools.MMAPPH1NPPR
import jline.lib.butools.MMAPPH1PRPR
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMResult
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

/**
 * Determines if arrivals at a station are RAP (Rational Arrival Process).
 * Checks source nodes for open models or considers inter-arrival in closed models.
 */
private fun hasRAPArrivals(sn: NetworkStruct, stationIdx: Int): Boolean {
    val station = sn.stations[stationIdx]
    // Check if any class at this station has RAP arrival process
    for (jobClass in sn.jobclasses) {
        val procType = sn.procid[station]?.get(jobClass)
        if (procType == ProcessType.RAP) {
            return true
        }
    }
    return false
}

/**
 * Determines if service at a station uses RAP or ME distributions.
 */
private fun hasRAPorMEService(sn: NetworkStruct, stationIdx: Int): Boolean {
    val station = sn.stations[stationIdx]
    for (jobClass in sn.jobclasses) {
        val procType = sn.procid[station]?.get(jobClass)
        if (procType == ProcessType.RAP || procType == ProcessType.ME) {
            return true
        }
    }
    return false
}

/**
 * Routes to appropriate queue solver based on arrival and service process types.
 *
 * Current routing logic:
 * - RAP arrivals + RAP service -> qbd_raprap1 (QBD method for RAP/RAP/1)
 * - All other cases -> MMAPPH1FCFS (BUTools MAP/PH/1 solver)
 *
 * Note: MMAPPH1FCFS can handle ME service if it represents a valid PH distribution
 * (i.e., alpha is a probability vector and A is a sub-generator). For general ME
 * distributions that are not valid PH, a dedicated ME solver would be needed.
 *
 * @param arrivalProc Arrival process (MAP or RAP)
 * @param serviceProc Service process per class (PH, ME, or RAP)
 * @param isRAParrivals True if arrivals are RAP
 * @param isRAPorMEservice True if service uses RAP or ME
 * @return Queue solution with performance metrics
 */
private fun solveQueue(
    arrivalProc: MatrixCell,
    servicePie: Map<Int?, Matrix>,
    serviceD0: Map<Int?, Matrix>,
    isRAParrivals: Boolean,
    isRAPorMEservice: Boolean,
    numQLMoms: Int?,
    numQLProbs: Int?,
    numSTMoms: Int?
): Map<String, MutableMap<Int, Matrix>> {

    // Use RAP/RAP/1 solver if both arrivals and service are RAP
    if (isRAParrivals && isRAPorMEservice) {
        // TODO: Implement proper RAP/RAP/1 or RAP/ME/1 solver routing
        // For now, fall through to MMAPPH1FCFS which works for valid PH representations
    }

    // Default: use MMAPPH1FCFS (handles MAP/PH, and ME if it's a valid PH)
    return MMAPPH1FCFS(
        arrivalProc,
        servicePie,
        serviceD0,
        numQLMoms,
        numQLProbs,
        numSTMoms,
        null,
        false,
        false,
        null,
        null
    )
}

fun solver_mam_basic(sn: NetworkStruct, options: SolverOptions): MAMResult {
    val config = options.config
    val tol = options.tol

    val PH = sn.proc
    val I = sn.nnodes
    val M = sn.nstations
    val K = sn.nclasses
    val C = sn.nchains
    val N = sn.njobs.transpose()
    val V = Matrix.cellsum(sn.visits)
    var S = Matrix.ones(sn.rates.numRows, sn.rates.numCols)
    S = S.elementDivide(sn.rates);
    val Lchain = snGetDemandsChain(sn).Dchain

    val QN = Matrix(M, K, M * K)
    val UN = Matrix(M, K, M * K)
    val RN = Matrix(M, K, M * K)
    val TN = Matrix(M, K, M * K)
    val WN = Matrix(M, K, M * K)
    val AN = Matrix(M, K, M * K)
    var CN = Matrix(1, K, K)
    val XN = Matrix(1, K, K)

    val pie: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    val D0: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()

    val lambda = Matrix(1, C, C)
    val chainSysArrivals: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    var TN_1 = Matrix(M, K, M * K)
    for (i in 0..<M) {
        for (j in 0..<K) {
            TN_1.set(i, j, Double.Companion.POSITIVE_INFINITY)
        }
    }

    var it = 0

    for (ist in 0..<M) {
        if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS || sn.sched.get(
                sn.stations.get(ist)
            ) == SchedStrategy.HOL || sn.sched.get(
                sn.stations.get(ist)
            ) == SchedStrategy.FCFSPRPRIO || sn.sched.get(
                sn.stations.get(ist)
            ) == SchedStrategy.PS
        ) {
            pie.put(ist, MatrixCell())
            D0.put(ist, MatrixCell())
            for (k in 0..<K) {
                PH.get(sn.stations.get(ist))!!.put(
                    sn.jobclasses.get(k),
                    map_scale(
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1),
                        S.get(ist, k) / sn.nservers.get(ist)
                    )
                )
                pie.get(ist)!!.set(
                    k,
                    map_pie(
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1)
                    )
                )
                D0.get(ist)!!.set(k, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0))
                if (D0.get(ist)!!.get(k).hasNaN()) {
                    PH.get(sn.stations.get(ist))!!
                        .put(sn.jobclasses.get(k), map_exponential(GlobalConstants.Immediate))
                    pie.get(ist)!!.set(k, Matrix.singleton(1.0))
                    D0.get(ist)!!.set(k, Matrix.singleton(-GlobalConstants.Immediate))
                }
            }
        }
    }

    var isOpen = false
    var isClosed = false

    for (i in 0..<sn.njobs.numRows) {
        for (j in 0..<sn.njobs.numCols) {
            if (!java.lang.Double.isFinite(sn.njobs.get(i, j))) {
                isOpen = true
                break
            }
        }
    }

    for (i in 0..<sn.njobs.numRows) {
        for (j in 0..<sn.njobs.numCols) {
            if (java.lang.Double.isFinite(sn.njobs.get(i, j))) {
                isClosed = true
                break
            }
        }
    }

    val isMixed = isOpen && isClosed

    if (isMixed) {
        // treat open as having higher priority than closed
        // sn.classprio(~isfinite(sn.njobs)) = 1 + max(sn.classprio(isfinite(sn.njobs)));
        // COMMENTED OUT to match MATLAB behavior - priority adjustment disabled
    }

    val lambdas_inchain = MatrixCell(C)
    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        lambdas_inchain.set(c, Matrix(inchain.numCols, 1, inchain.numCols))
        val ist = sn.refstat.get(inchain.value().toInt(), 0).toInt()
        if (!PH.containsKey(sn.stations.get(ist))) {
            PH.put(sn.stations.get(ist), HashMap<JobClass?, MatrixCell?>())
        }
        for (j in 0..<inchain.numCols) {
            lambdas_inchain.get(c).set(j, 0, sn.rates.get(ist, inchain.get(0, j).toInt()))
        }
        var sum = 0.0
        for (k in 0..<lambdas_inchain.get(c).numRows) {
            if (java.lang.Double.isFinite(lambdas_inchain.get(c).get(k, 0))) {
                sum += lambdas_inchain.get(c).get(k, 0)
            }
        }
        lambda.set(0, c, sum)
        var openChain = false
        for (j in 0..<inchain.numCols) {
            if (Utils.isInf(sn.njobs.get(inchain.get(0, j).toInt()))) {
                openChain = true
                break
            }
        }
        if (openChain) {
            for (k in 0..<K) {
                // if the matrix is NaN, then matrix has only one value and its is Double.NaN
                if (java.lang.Double.isNaN(
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0).get(0)
                    )
                ) {
                    PH.get(sn.stations.get(ist))!!
                        .put(sn.jobclasses.get(k), map_exponential(Double.Companion.POSITIVE_INFINITY))
                }
            }
            var k = inchain.get(0)
            chainSysArrivals.put(c, MatrixCell())
            chainSysArrivals.get(c)!!
                .set(0, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(0))
            chainSysArrivals.get(c)!!
                .set(1, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(1))
            chainSysArrivals.get(c)!!
                .set(2, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(1))
            for (ki in 1..<inchain.length()) {
                k = inchain.get(ki)
                if (java.lang.Double.isNaN(
                        PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(0)
                            .get(0)
                    )
                ) {
                    PH.get(sn.stations.get(ist))!!
                        .put(sn.jobclasses.get(k.toInt()), map_exponential(Double.Companion.POSITIVE_INFINITY))
                }
                val MMAPS: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
                MMAPS.put(0, chainSysArrivals.get(c))
                MMAPS.put(1, MatrixCell())
                MMAPS.get(1)!!.set(0, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(0))
                MMAPS.get(1)!!.set(1, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(1))
                MMAPS.get(1)!!.set(2, PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k.toInt()))!!.get(1))
                chainSysArrivals.put(
                    c,
                    mmap_super_safe(
                        MMAPS.filterValues { it != null }.mapValues { it.value!! },
                        config.space_max,
                        "default"
                    )
                )
            }
            for (i in 0..<inchain.numCols) {
                TN.set(ist, inchain.get(0, i).toInt(), lambdas_inchain.get(c).get(i, 0))
            }
        }
    }

    val sd = Matrix(sn.nservers.numRows, sn.nservers.numCols, sn.nservers.nonZeroLength)
    for (i in 0..<sn.nservers.numRows) {
        if (!Utils.isInf(sn.nservers.get(i, 0))) {
            sd.set(i, 0, 1)
        }
    }

    //TODO: check from here
    var dif_matrix = TN.add(-1.0, TN_1)
    var dif = FastMath.abs(dif_matrix.elementMaxAbs())
    while (dif > tol && it <= options.iter_max) {
        it++
        TN_1 = TN.copy()
        val Umax_matrix = Matrix(M, 1, M)
        for (i in 0..<M) {
            if (sd.get(i, 0) == 0.0) {
                Umax_matrix.set(i, 0, -Double.Companion.POSITIVE_INFINITY)
            } else {
                var sum = 0.0
                for (j in 0..<K) {
                    sum += UN.get(i, j)
                }
                Umax_matrix.set(i, 0, sum)
            }
        }

        val Umax = Umax_matrix.elementMax()
        if (Umax >= 1) { // DO NOT CHANGE TO >1 AS MATLAB
            lambda.divideEq(Umax)
        } else {
            for (c in 0..<C) {
                val inchain: Matrix = sn.inchain.get(c)!!
                var openChain = false
                for (j in 0..<inchain.numCols) {
                    if (Utils.isInf(sn.njobs.get(inchain.get(0, j).toInt()))) {
                        openChain = true
                        break
                    }
                }
                if (!openChain) {
                    var Nc = 0.0
                    for (i in 0..<inchain.length()) {
                        Nc = Nc + sn.njobs.get(inchain.get(i).toInt())
                    }
                    val QN_omitnan = Matrix(QN)
                    QN_omitnan.removeNaN()
                    val QN_col_sum = Matrix(inchain.length(), 1, inchain.length())
                    for (i in 0..<inchain.length()) {
                        QN_col_sum.set(i, 0, QN_omitnan.sumCols(inchain.get(i).toInt()))
                    }
                    var QNc = QN_col_sum.elementSum()
                    QNc = FastMath.max(options.tol, QNc)
                    val TNlb = Nc / Lchain.sumCols(c)
                    if (it == 1) {
                        lambda.set(0, c, TNlb)
                    } else {
                        lambda.set(
                            0,
                            c,
                            lambda.get(0, c) * it / options.iter_max + (Nc / QNc) * lambda.get(
                                0,
                                c
                            ) * (options.iter_max - it) / options.iter_max
                        )
                    }
                }
            }
        }

        for (c in 0..<C) {
            val inchain: Matrix = sn.inchain.get(c)!!
            val lambda_c = Matrix(1, 1, 1)
            lambda_c.set(0, 0, lambda.get(c))
            chainSysArrivals.put(c, mmap_exponential(lambda_c))
            for (m in 0..<M) {
                for (i in 0..<inchain.length()) {
                    TN.set(m, inchain.get(i).toInt(), V.get(m, inchain.get(i).toInt()) * lambda.get(c))
                }
            }
        }

        for (ind in 0..<I) {
            if (sn.isstation.get(ind) == 1.0) {
                val ist = sn.nodeToStation.get(ind).toInt()
                if (sn.nodetype.get(ind) == NodeType.Join) {
                    for (c in 0..<C) {
                        val inchain: Matrix = sn.inchain.get(c)!!
                        for (i in 0..<inchain.length()) {
                            var fanin = 0
                            for (j in 0..<sn.rtnodes.numRows) {
                                if (sn.rtnodes.get(j, ((ind - 1) * K + inchain.get(i)).toInt()) != 0.0) {
                                    fanin++
                                }
                            }
                            val k = inchain.get(i).toInt()
                            TN.set(ist, k, lambda.get(c) * V.get(ist, k) / fanin)
                            UN.set(ist, k, 0.0)
                            QN.set(ist, k, 0.0)
                            RN.set(ist, k, 0.0)
                        }
                    }
                } else if (sn.nodetype.get(ind) == NodeType.Fork) {
                    // Fork nodes: fan-out operation, set metrics to zero like Join nodes
                    // Traffic flow is handled through routing matrix in solver_mam_traffic
                    for (c in 0..<C) {
                        val inchain: Matrix = sn.inchain.get(c)!!
                        for (i in 0..<inchain.length()) {
                            val k = inchain.get(i).toInt()
                            // Fork throughput is based on input traffic
                            TN.set(ist, k, lambda.get(c) * V.get(ist, k))
                            UN.set(ist, k, 0.0)
                            QN.set(ist, k, 0.0)
                            RN.set(ist, k, 0.0)
                        }
                    }
                } else {
                    if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                        for (c in 0..<C) {
                            val inchain: Matrix = sn.inchain.get(c)!!
                            for (i in 0..<inchain.length()) {
                                val k = inchain.get(i).toInt()
                                TN.set(ist, k, lambda.get(c) * V.get(ist, k))
                                UN.set(ist, k, S.get(ist, k) * TN.get(ist, k))
                                QN.set(ist, k, TN.get(ist, k) * S.get(ist, k) * V.get(ist, k))
                                if (V.get(ist, k) <= GlobalConstants.Zero) {
                                    RN.set(ist, k, 0.0)
                                } else {
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
                                }
                            }
                        }
                    } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.PS) {
                        for (c in 0..<C) {
                            val inchain: Matrix = sn.inchain.get(c)!!
                            for (i in 0..<inchain.length()) {
                                val k = inchain.get(i).toInt()
                                TN.set(ist, k, lambda.get(c) * V.get(ist, k))
                                UN.set(ist, k, S.get(ist, k) * TN.get(ist, k))
                            }
                            val Uden = FastMath.min(1 - GlobalConstants.FineTol, UN.sumRows(ist))
                            for (i in 0..<inchain.length()) {
                                val k = inchain.get(i).toInt()
                                QN.set(ist, k, UN.get(ist, k) / (1 - Uden))
                                if (V.get(ist, k) <= GlobalConstants.Zero) {
                                    RN.set(ist, k, 0.0)
                                } else {
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
                                }
                            }
                        }
                    } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL || sn.sched.get(
                            sn.stations.get(
                                ist
                            )
                        ) == SchedStrategy.FCFS || sn.sched.get(
                            sn.stations.get(
                                ist
                            )
                        ) == SchedStrategy.FCFSPRPRIO
                    ) {
                        val chainArrivalAtNode: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
                        val rates: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
                        var aggrArrivalAtNode: MatrixCell? = MatrixCell()
                        for (c in 0..<C) {
                            val a = Matrix.extractRows(V, ist, ist + 1, null)
                            a.scaleEq(lambda.get(c))
                            if (c == 0) {
                                rates.put(ist, MatrixCell())
                            }
                            rates.get(ist)!!.set(c, a)
                            val inchain: Matrix = sn.inchain.get(c)!!
                            val markProb = Matrix(1, inchain.length(), inchain.length())
                            for (i in 0..<inchain.length()) {
                                markProb.set(0, i, rates.get(ist)!!.get(c).get(inchain.get(i).toInt()))
                            }
                            markProb.scaleEq(1.0 / markProb.elementSum())
                            markProb.removeNaN()
                            chainArrivalAtNode.put(c, mmap_mark(chainSysArrivals.get(c)!!, markProb))
                            chainArrivalAtNode.put(c, mmap_normalize(chainArrivalAtNode.get(c)!!))
                            val b = Matrix(
                                rates.get(ist)!!.get(c).numRows,
                                rates.get(ist)!!.get(c).numCols,
                                rates.get(ist)!!.get(c).numRows * rates.get(ist)!!.get(c).numCols
                            )
                            for (i in 0..<b.numRows) {
                                for (j in 0..<b.numCols) {
                                    b.set(i, j, 1 / rates.get(ist)!!.get(c).get(i, j))
                                }
                            }
                            chainArrivalAtNode.put(c, mmap_scale(chainArrivalAtNode.get(c)!!, b))
                            if (c == 0) {
                                val MMAPS: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
                                MMAPS.put(0, chainArrivalAtNode.get(c))
                                MMAPS.put(1, mmap_exponential(Matrix.singleton(0.0), 1))
                                val aggrArrivalAtNode_temp =
                                    mmap_super_safe(
                                        MMAPS.filterValues { it != null }.mapValues { it.value!! },
                                        config.space_max,
                                        "default"
                                    )
                                aggrArrivalAtNode!!.set(0, aggrArrivalAtNode_temp!!.get(0))
                                aggrArrivalAtNode.set(1, aggrArrivalAtNode_temp.get(1))
                                aggrArrivalAtNode.set(2, aggrArrivalAtNode_temp.get(1))
                                val lc =
                                    map_lambda(chainArrivalAtNode.get(c)!!.get(0), chainArrivalAtNode.get(c)!!.get(1))
                                if (lc > 0) {
                                    aggrArrivalAtNode = mmap_scale(aggrArrivalAtNode, Matrix.singleton(1.0 / lc))
                                }
                            } else {
                                val MMAPS: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
                                MMAPS.put(0, aggrArrivalAtNode)
                                MMAPS.put(1, chainArrivalAtNode.get(c))
                                aggrArrivalAtNode =
                                    mmap_super_safe(
                                        MMAPS.filterValues { it != null }.mapValues { it.value!! },
                                        config.space_max,
                                        "default"
                                    )
                            }
                        }
                        var Qret: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
                        val priorityAnalysis = analyzePriorities(sn)

                        val schedStrat = sn.sched.get(sn.stations.get(ist))
                        if ((schedStrat == SchedStrategy.HOL || schedStrat == SchedStrategy.FCFSPRPRIO) && !priorityAnalysis.isIdentical
                        ) {
                            if (priorityAnalysis.isAllDistinct) {
                                // Get unique priorities and their indices (equivalent to MATLAB's [uK,iK] = unique(sn.classprio))
                                val priorityMap = mutableMapOf<Double, Int>()
                                for (k in 0..<K) {
                                    priorityMap[sn.classprio.get(k)] = k
                                }
                                // BUTools convention: D1=lowest priority, DK=highest priority
                                // LINE convention: lower value = higher priority
                                // sortedDescending puts highest LINE priority (lowest value) last, matching BUTools
                                val sortedPriorities = priorityMap.keys.sortedDescending()
                                val iK = IntArray(K) // Priority-sorted indices
                                for (i in sortedPriorities.indices) {
                                    iK[i] = priorityMap[sortedPriorities[i]]!!
                                }

                                // Reorder arrival process: aggrArrivalAtNode{[1;2+iK]}
                                val reorderedArrival = MatrixCell()
                                reorderedArrival.set(
                                    0,
                                    aggrArrivalAtNode!!.get(0)
                                ) // D0 matrix (index 1 in MATLAB, 0 in Java)
                                for (i in iK.indices) {
                                    reorderedArrival.set(
                                        1 + i,
                                        aggrArrivalAtNode.get(1 + iK[i])
                                    ) // D1 matrices (2+iK in MATLAB, 1+iK in Java)
                                }

                                // Reorder pie and D0 by priority: {pie{ist}{iK}}, {D0{ist,iK}}
                                val reorderedPie = mutableMapOf<Int, Matrix>()
                                val reorderedD0 = mutableMapOf<Int, Matrix>()
                                for (i in iK.indices) {
                                    reorderedPie[i] = pie.get(ist)!!.get(iK[i])
                                    reorderedD0[i] = D0.get(ist)!!.get(iK[i])
                                }

                                // Convert maps to MatrixCell objects for priority analysis
                                val reorderedPieCell = MatrixCell()
                                val reorderedD0Cell = MatrixCell()
                                for ((key, value) in reorderedPie) {
                                    reorderedPieCell.set(key, value)
                                }
                                for ((key, value) in reorderedD0) {
                                    reorderedD0Cell.set(key, value)
                                }

                                // Select appropriate priority solver based on scheduling strategy
                                val prioResult = if (schedStrat == SchedStrategy.FCFSPRPRIO) {
                                    // Use MMAPPH1PRPR for preemptive-resume priority
                                    MMAPPH1PRPR(
                                        reorderedArrival,
                                        reorderedPieCell,
                                        reorderedD0Cell,
                                        1,
                                        null, null, null, null, null, null
                                    )
                                } else {
                                    // Use MMAPPH1NPPR for non-preemptive (HOL) priority
                                    MMAPPH1NPPR(
                                        reorderedArrival,
                                        reorderedPieCell,
                                        reorderedD0Cell,
                                        1,
                                        null, null, null, null, null, null
                                    )
                                }

                                val prioResultMoms = prioResult.get("ncMoms")?.let { moms ->
                                    moms.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }.toMutableMap()
                                } ?: mutableMapOf()

                                // Restore original class order: [Qret{iK'}] equivalent
                                Qret = HashMap<Int?, Matrix?>()
                                for (i in iK.indices) {
                                    if (prioResultMoms.containsKey(i)) {
                                        Qret[iK[i]] = prioResultMoms[i]
                                    }
                                }
                            } else {
                                throw RuntimeException(priorityAnalysis.message)
                            }
                        } else {
                            val sn_rates_ist_k = Matrix(1, K, K)
                            for (i in 0..<K) {
                                sn_rates_ist_k.set(i, sn.rates.get(ist, i) * sn.nservers.get(ist))
                            }
                            val aggrUtil =
                                mmap_lambda(aggrArrivalAtNode!!).elementDivide(
                                    sn_rates_ist_k.elementIncrease(
                                        GlobalConstants.FineTol
                                    )
                                ).elementSum()
                            if (aggrUtil < 1 - GlobalConstants.FineTol) {
                                var closed = true
                                for (i in 0..<N.length()) {
                                    if (Utils.isInf(N.get(i))) {
                                        closed = false
                                        break
                                    }
                                }

                                // Check for MAP/D/c: single-class open model with deterministic service
                                var isMapDc = false
                                if (!closed && K == 1) {
                                    val station = sn.stations.get(ist)
                                    val jobClass = sn.jobclasses.get(0)
                                    val procType = sn.procid[station]?.get(jobClass)
                                    isMapDc = (procType == ProcessType.DET)
                                }

                                if (isMapDc) {
                                    // Use exact MAP/D/c solver
                                    val arrivalMAP = mmap_shorten(aggrArrivalAtNode)
                                    val D0_arr = arrivalMAP.get(0)
                                    val D1_arr = arrivalMAP.get(1)
                                    val detServiceTime = S.get(ist, 0)
                                    val numServers = sn.nservers.get(ist).toInt()

                                    try {
                                        val mapdcResult = qsys_mapdc(D0_arr, D1_arr, detServiceTime, numServers)
                                        Qret.put(0, Matrix.singleton(mapdcResult.meanQueueLength))
                                        // Mark this station as using MAP/D/c solver
                                    } catch (e: Exception) {
                                        // Fall back to standard solver if MAP/D/c fails
                                        Qret = MMAPPH1FCFS(
                                            mmap_shorten(aggrArrivalAtNode),
                                            pie.get(ist)!!.toMap(),
                                            D0.get(ist)!!.toMap(),
                                            1,
                                            null,
                                            null,
                                            null,
                                            false,
                                            false,
                                            null,
                                            null
                                        ).get("ncMoms")?.let { moms ->
                                            moms.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }.toMutableMap()
                                        } ?: mutableMapOf()
                                    }
                                } else if (!closed) {
                                    Qret = MMAPPH1FCFS(
                                        mmap_shorten(aggrArrivalAtNode),
                                        pie.get(ist)!!.toMap(),
                                        D0.get(ist)!!.toMap(),
                                        1,
                                        null,
                                        null,
                                        null,
                                        false,
                                        false,
                                        null,
                                        null
                                    ).get("ncMoms")?.let { moms ->
                                        moms.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }.toMutableMap()
                                    } ?: mutableMapOf()
                                } else {
                                    val finite_N = N.copy()
                                    finite_N.removeInfinite()
                                    val maxLevel = finite_N.elementMax() + 1
                                    val D = mmap_shorten(aggrArrivalAtNode)
                                    var pdistr: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
                                    if (map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                                        for (k in 0..<K) {
                                            val pdistrK = Matrix(1, 2, 2)
                                            pdistrK.set(0, 1 - GlobalConstants.FineTol)
                                            pdistrK.set(1, GlobalConstants.FineTol)
                                            pdistr.put(k, pdistrK)
                                            Qret.put(k, Matrix.singleton(GlobalConstants.FineTol / sn.rates.get(ist)))
                                        }
                                    } else {
                                        pdistr = MMAPPH1FCFS(
                                            D,
                                            pie.get(ist)!!.toMap(),
                                            D0.get(ist)!!.toMap(),
                                            null,
                                            maxLevel.toInt(),
                                            null,
                                            null,
                                            false,
                                            false,
                                            null,
                                            null
                                        ).get("ncDistr")?.let { distr ->
                                            distr.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }
                                                .toMutableMap()
                                        } ?: mutableMapOf()
                                        for (k in 0..<K) {
                                            pdistr.put(
                                                k,
                                                Matrix.extractRows(
                                                    pdistr.get(k)!!.transpose(),
                                                    0,
                                                    N.get(k).toInt() + 1,
                                                    null
                                                )
                                            )
                                            pdistr.get(k)!!.absEq()
                                            var sum = 0.0
                                            for (i in 0..<pdistr.get(k)!!.length() - 1) {
                                                sum = sum + pdistr.get(k)!!.get(i)
                                            }
                                            pdistr.get(k)!!.set(pdistr.get(k)!!.length() - 1, abs(1 - sum))
                                            pdistr.get(k)!!.scaleEq(1 / pdistr.get(k)!!.elementSum())
                                            val a = Matrix(1, N.get(k).toInt() + 1, N.get(k).toInt() + 1)
                                            for (i in 0..<a.length()) {
                                                a.set(i, i.toDouble())
                                            }
                                            val b = Matrix(1, N.get(k).toInt() + 1, N.get(k).toInt() + 1)
                                            for (i in 0..<a.length()) {
                                                b.set(i, pdistr.get(k)!!.get(i))
                                            }
                                            Qret.put(
                                                k,
                                                Matrix.singleton(max(0.0, min(N.get(k), a.mult(b.transpose()).get(0))))
                                            )
                                        }
                                    }
                                }
                            } else {
                                // Saturated queue - preserve njobs values (including infinity for open classes)
                                for (k in 0..<K) {
                                    val njobValue = sn.njobs.get(k)
                                    Qret.put(k, Matrix.singleton(njobValue))
                                }
                            }
                        }
                        for (i in 0..<Qret.size) {
                            QN.set(ist, i, Qret.get(i)!!.get(0))
                        }
                        for (k in 0..<K) {
                            var c = 0
                            for (i in 0..<sn.chains.numRows) {
                                if (sn.chains.get(i, k) != 0.0) {
                                    c = i
                                    break
                                }
                            }
                            TN.set(ist, k, rates.get(ist)!!.get(c).get(k))
                            UN.set(ist, k, TN.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist))
                            QN.set(ist, k, Qret.get(k)!!.get(0))
                            QN.set(
                                ist,
                                k,
                                QN.get(ist, k) + TN.get(ist, k) * S.get(
                                    ist,
                                    k
                                ) * (sn.nservers.get(ist) - 1) / sn.nservers.get(ist)
                            )
                            if (V.get(ist, k) <= GlobalConstants.Zero) {
                                RN.set(ist, k, 0.0)
                            } else {
                                RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
                            }
                        }
                    }
                }
            } else {
                // Other node types are handled by default traffic flow calculations
                // Fork nodes are now supported above
            }
        }
        dif_matrix = TN.add(-1.0, TN_1)
        dif = dif_matrix.elementMaxAbs()
    }
    val totiter = it + 2

    CN = RN.sumCols()
    QN.absEq()
    for (i in 1..2) {
        for (c in 0..<C) {
            val inchain: Matrix = sn.inchain.get(c)!!
            var Nc = 0.0
            for (j in 0..<inchain.length()) {
                Nc = Nc + sn.njobs.get(inchain.get(j).toInt())
            }
            if (java.lang.Double.isFinite(Nc)) {
                var QNc = 0.0
                for (j in 0..<inchain.length()) {
                    QNc = QNc + QN.sumCols(inchain.get(j).toInt())
                }
                for (m in 0..<inchain.length()) {
                    for (n in 0..<QN.numRows) {
                        QN.set(n, inchain.get(m).toInt(), QN.get(n, inchain.get(m).toInt()) * (Nc / QNc))
                    }
                }
            }

            for (ind in 0..<I) {
                for (k in 0..<inchain.length()) {
                    if (sn.isstation.get(ind) == 1.0) {
                        val ist = sn.nodeToStation.get(ind).toInt()
                        if (V.get(ist, inchain.get(k).toInt()) > GlobalConstants.Zero) {
                            if (Utils.isInf(sn.nservers.get(ist))) {
                                RN.set(ist, inchain.get(k).toInt(), S.get(ist, inchain.get(k).toInt()))
                            } else {
                                RN.set(
                                    ist,
                                    inchain.get(k).toInt(),
                                    max(
                                        S.get(ist, inchain.get(k).toInt()),
                                        QN.get(ist, inchain.get(k).toInt()) / TN.get(ist, inchain.get(k).toInt())
                                    )
                                )
                            }
                        } else {
                            RN.set(ist, inchain.get(k).toInt(), 0)
                        }
                        QN.set(
                            ist,
                            inchain.get(k).toInt(),
                            RN.get(ist, inchain.get(k).toInt()) * TN.get(ist, inchain.get(k).toInt())
                        )
                    }
                }
            }
            if (Nc == 0.0) {
                for (k in 0..<QN.numRows) {
                    QN.set(k, c, 0)
                }
                for (k in 0..<UN.numRows) {
                    UN.set(k, c, 0)
                }
                for (k in 0..<RN.numRows) {
                    RN.set(k, c, 0)
                }
                for (k in 0..<TN.numRows) {
                    TN.set(k, c, 0)
                }
                CN.set(0, c, 0)
                XN.set(0, c, 0)
            }
        }
    }
    val result = MAMResult()
    result.QN = QN
    result.UN = UN
    result.RN = RN
    result.TN = TN
    result.CN = CN
    result.XN = XN
    result.WN = WN
    result.AN = AN
    result.iter = totiter
    return result
}