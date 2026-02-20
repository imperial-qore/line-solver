package jline.solvers.mam.handlers

import jline.api.mam.*
import jline.io.line_warning
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.lib.butools.MMAPPH1FCFS
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMResult
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.lang.Double
import kotlin.Any
import kotlin.Int

fun solver_mam(sn: NetworkStruct, options: SolverOptions): MAMResult {
    val startTime = System.nanoTime()
    val method = options.method
    val config = options.config
    val PH = sn.proc
    val M = sn.nstations
    val K = sn.nclasses
    val C = sn.nchains
    val V = Matrix(sn.visits.get(0))
    if (sn.nchains > 1) {
        for (i in 1..<sn.nchains) {
            V.add(1.0, sn.visits.get(i))
        }
    }

    val QN = Matrix(M, K, M * K)
    val UN = Matrix(M, K, M * K)
    val RN = Matrix(M, K, M * K)
    val TN = Matrix(M, K, M * K)
    val AN = Matrix(M, K, M * K)
    val WN = Matrix(M, K, M * K)
    val CN = Matrix(1, K, K)
    val XN = Matrix(1, K, K)

    val lambda = Matrix(1, K, K)
    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        val lambdas_inchain = Matrix(1, inchain.length(), inchain.length())
        for (i in 0..<inchain.length()) {
            val rate = sn.rates.get(sn.refstat.get(inchain.get(0).toInt()).toInt(), inchain.get(i).toInt())
            lambdas_inchain.set(0, i, rate)
        }
        var sum = 0.0
        for (j in 0..<lambdas_inchain.length()) {
            if (Double.isFinite(lambdas_inchain.get(j))) {
                sum = sum + lambdas_inchain.get(j)
            }
        }
        for (i in 0..<inchain.length()) {
            lambda.set(0, inchain.get(i).toInt(), sum)
        }
    }

    val chain = Matrix(1, K, K)

    for (k in 0..<K) {
        for (i in 0..<chain.numCols) {
            if (sn.chains.get(i, k) != 0.0) {
                chain.set(0, k, i)
                break
            }
        }
    }
    
    // Check for non-FCFS queues
    for (ist in 0..<sn.nstations) {
        val sched = sn.sched.get(sn.stations.get(ist))
        if (sched != SchedStrategy.EXT && sched != SchedStrategy.FCFS && sched != SchedStrategy.HOL) {
            if (options.verbose != VerboseLevel.SILENT) {
                line_warning(mfilename(object : Any() {}), "The dec.mmap method does not support non-FCFS queues.")
            }
            val emptyResult = MAMResult()
            emptyResult.QN = Matrix(0, 0)
            emptyResult.UN = Matrix(0, 0)
            emptyResult.RN = Matrix(0, 0)
            emptyResult.TN = Matrix(0, 0)
            emptyResult.CN = Matrix(0, 0)
            emptyResult.XN = Matrix(0, 0)
            emptyResult.iter = 0
            emptyResult.method = ""
            emptyResult.runtime = (System.nanoTime() - startTime) / 1000000000.0
            return emptyResult
        }
    }

    var isopen = true
    for (i in 0..<sn.njobs.length()) {
        if (Double.isFinite(sn.njobs.get(i))) {
            isopen = false
            break
        }
    }

    val result = MAMResult()
    if (isopen) {
        var last_it = 0
        val pie: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
        val D0: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
        for (ist in 0..<M) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                for (i in 0..<TN.numCols) {
                    if (!Double.isNaN(sn.rates.get(ist, i))) {
                        TN.set(ist, i, sn.rates.get(ist, i))
                    }
                }
            } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL) {
                for (k in 0..<K) {
                    val D0_ = PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0).copy()
                    D0_.scaleEq(1.0 / sn.nservers.get(ist))
                    val D1 = PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1).copy()
                    D1.scaleEq(1.0 / sn.nservers.get(ist))
                    PH.get(sn.stations.get(ist))!!.put(sn.jobclasses.get(k),
                        map_scale(PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                            PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1),
                            map_mean(D0_, D1)))
                    if (k == 0) {
                        pie.put(ist, MatrixCell())
                        D0.put(ist, MatrixCell())
                    }
                    pie.get(ist)!!.set(k,
                        map_pie(PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                            PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1)))

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
        val it_max = options.iter_max
        val DEP: MutableMap<Int?, MutableMap<Int?, MatrixCell?>?> =
            ph_reindex(sn).mapKeys { it.key as Int? }.mapValues { it.value?.toMutableMap()?.mapKeys { k -> k.key as Int? }?.toMutableMap() }.toMutableMap()
        for (it in 0..<it_max) {
            last_it = it
            if (it == 1) {
                for (ind in 0..<M) {
                    for (r in 0..<K) {
                        val ist = sn.nodeToStation.get(ind).toInt()
                        DEP.get(ind)!!.put(r,
                            map_scale(PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!.get(0),
                                PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!.get(1),
                                1 / (lambda.get(r) * V.get(ind, r))))
                    }
                }
            }

            val ARV = solver_mam_traffic(sn, DEP, config)
            val QN_1 = QN.copy()
            for (ist in 0..<M) {
                val ind = sn.stationToNode.get(ist).toInt()
                if (sn.nodetype.get(ind) == NodeType.Queue) {
                    if (ARV.get(ind)!!.get(0).length() > config.space_max) {
                        println("Arrival process at node " + ind + " is now at " + ARV.get(ind)!!.get(0)
                            .length() + " states. Compressing")
                        val mmapCell = ARV.get(ind)!!
                        val mmapArray = Array(mmapCell.size()) { mmapCell.get(it) }
                        val compressedArray = Mmap_compress.mmap_compress(mmapArray)
                        val compressedCell = MatrixCell(compressedArray.size)
                        for (i in compressedArray.indices) {
                            compressedCell.set(i, compressedArray[i])
                        }
                        ARV.put(ind, compressedCell)
                    }
                    val D = MatrixCell()
                    for (i in 0..<ARV.get(ind)!!.size()) {
                        if (i == 0) {
                            D.set(0, ARV.get(ind)!!.get(0))
                        } else if (i > 1) {
                            D.set(i, ARV.get(ind)!!.get(i - 1))
                        }
                    }
                    val mmapResult: Map<String, MutableMap<Int, Matrix>> = MMAPPH1FCFS(D,
                        pie.get(ist)!!.toMap(),
                        D0.get(ist)!!.toMap(),
                        1,
                        2,
                        null,
                        null,
                        false,
                        false,
                        null,
                        null)
                    val Qret: MutableMap<Int, Matrix> = mmapResult.get("ncMoms")!!
                    for (k in 0..<K) {
                        QN.set(ist, k, Qret[k]!!.elementSum())
                    }
                    TN.insertSubMatrix(ist, 0, ist + 1, TN.numCols, mmap_lambda(ARV.get(ind)!!))
                }
                for (k in 0..<K) {
                    UN.set(ist,
                        k,
                        TN.get(ist, k) * map_mean(PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                            PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1)))
                    // Add number of jobs at the surrogate delay server
                    val meanServiceTime = map_mean(PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                            PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!.get(1))
                    QN.set(ist,
                        k,
                        QN.get(ist, k) + TN.get(ist, k) * meanServiceTime * sn.nservers.get(ist) * (sn.nservers.get(ist) - 1) / sn.nservers.get(ist))

                    if (V.get(ist, k) <= GlobalConstants.Zero) {
                        RN.set(ist, k, 0.0)
                    } else {
                        RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
                    }
                }
            }

            if (it >= 3) {
                val Q_diff = QN.add(-1.0, QN_1)
                Q_diff.absEq()
                if (Q_diff.elementDivide(QN_1).elementMax() < options.iter_tol) {
                    break
                }
            }

            for (ist in 0..<M) {
                val ind = sn.stationToNode.get(ist).toInt()
                if (sn.nodetype.get(ind) == NodeType.Queue) {
                    for (r in 0..<K) {
                        val types = Matrix(1, K - 1, K - 1)
                        var idx = 0
                        for (i in 0..<K) {
                            if (i != r) {
                                types.set(idx, i.toDouble())
                                idx++
                            }
                        }
                        var A = mmap_hide(ARV.get(ind)!!, types)
                        val tA = A!!.get(1).sumRows()
                        val pieA = map_pie(A.get(0), A.get(1))

                        val S: MatrixCell = PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!
                        val pieS = map_pie(S.get(0), S.get(1))
                        val tS = S.get(1).sumRows()

                        val rho = UN.sumRows(ist)
                        
                        // Calculate queue seen upon arrival with PASTA
                        val AQ = QN.sumRows(ist)
                        val Afull = AQ / rho // from AQ = (1-rho)*0 + rho*AQfull
                        val pfullonarrival = 1.0 - (Afull) / (1.0 + Afull)

                        A = mmap_scale(A, Matrix.singleton(map_mean(A.get(0), A.get(1)) - map_mean(S.get(0), S.get(1))))
                        val zAS = Matrix.createLike(tA.mult(pieS))
                        val zSA = Matrix.createLike(tS.mult(pieA))
                        val zA = Matrix.createLike(A!!.get(1))

                        val DEP0ir = Matrix.concatRows(Matrix.concatColumns(S.get(0),
                            Matrix.scaleMult(tS.mult(pieA), rho),
                            null), Matrix.concatColumns(zAS, A.get(0), null), null)
                        val DEP1ir =
                            Matrix.concatRows(Matrix.concatColumns(Matrix.scaleMult(S.get(1), 1 - rho), zSA, null),
                                Matrix.concatColumns(tA.mult(pieS), zA, null),
                                null)
                        DEP.get(ind)!!.put(r, map_normalize(DEP0ir, DEP1ir))
                        DEP.get(ind)!!.put(r,
                                map_scale(DEP.get(ind)!!.get(r)!!.get(0),
                                    DEP.get(ind)!!.get(r)!!.get(1),
                                    1 / (lambda.get(r) * V.get(ind, r))))
                    }
                }
            }
            result.iter = it
        }
        result.iter = last_it
        if (options.verbose != VerboseLevel.SILENT) {
            println("MAM parametric decomposition completed in " + last_it + " iterations")
        }
    } else {
        if (options.verbose != VerboseLevel.SILENT) {
            line_warning(mfilename(object : Any() {}), "This model is not supported by SolverMAM yet. Returning with no result.")
        }
    }
    result.QN = QN
    result.UN = UN
    result.RN = RN
    result.TN = TN
    result.WN = WN
    result.AN = AN
    result.CN = CN
    result.XN = XN
    result.runtime = (System.nanoTime() - startTime) / 1000000000.0
    return result
}