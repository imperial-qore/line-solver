/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.api.pfqn.mva.pfqn_qzgblow
import jline.api.pfqn.mva.pfqn_qzgbup
import jline.api.pfqn.mva.pfqn_xzgsblow
import jline.api.pfqn.mva.pfqn_xzgsbup
import jline.api.pfqn.pfqn_harel_bounds
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * MVA Analyzer class for bounding methods
 */
fun solver_mva_bound_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val res = MVAResult()
    val startTime = System.nanoTime()
    var endTime = startTime
    val iter = 1
    val method = options.method
    var QN = Matrix(0, 0)
    var UN = Matrix(0, 0)
    var RN = Matrix(0, 0)
    var TN = Matrix(0, 0)
    var CN = Matrix(0, 0)
    val WN = Matrix(0, 0)
    var XN = Matrix(0, 0)
    var lG = Double.NaN
    when (method) {
        "aba.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                val Dsum = D.elementSum()
                CN = Matrix(1, 1)
                CN.set(0, 0, Z + N * Dsum)
                XN = Matrix(1, 1)
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum)))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i))
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N)
                        }
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "aba.lower" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dsum = D.elementSum()
                XN = Matrix(1, 1)
                XN.set(0, 0, N / (Z + N * Dsum))
                CN = Matrix(1, 1)
                CN.set(0, 0, Z + Dsum)
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "bjb.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                val Dsum = D.elementSum()

                val Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum))
                val Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum)

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1))
                XN = Matrix(1, 1)
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum + D.meanCol().value() * (N - 1 - Z * Xaba_upper_1))))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA upper
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i))
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N)
                        }
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "bjb.lower" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                val Dsum = D.elementSum()

                val Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum))
                val Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum)

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + Dsum + D.meanCol().value() * (N - 1 - Z * Xaba_upper_1))
                XN = Matrix(1, 1)
                XN.set(0, 0, N / (Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1)))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA lower
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "pb.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                val Dsum = D.elementSum()

                val Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum))
                val Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum)

                val Dpb2 = D.elementPower(2.0).elementSum() / Dsum
                val DpbN = D.elementPower(N).elementSum() / D.elementPower(N - 1).elementSum()

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1))
                XN = Matrix(1, 1)
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1))))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA upper
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i))
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N)
                        }
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "pb.lower" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                val Dsum = D.elementSum()

                val Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum))
                val Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum)

                val Dpb2 = D.elementPower(2.0).elementSum() / Dsum
                val DpbN = D.elementPower(N).elementSum() / D.elementPower(N - 1).elementSum()

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1))
                XN = Matrix(1, 1)
                XN.set(0, 0, N / (Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1)))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA lower
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "sb.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) throw RuntimeException("Unsupported method for a model with infinite-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()

                val A3 = D.elementPower(3.0).elementSum()
                val A2 = D.elementPower(2.0).elementSum()
                val A1 = D.elementPower(1.0).elementSum()

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + A1 + (N - 1) * (A1 * A2 + A3) / (A1 * A1 + A2))
                XN = Matrix(1, 1)
                XN.set(0, 0, Maths.min(1 / Dmax, N / CN.value()))
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA lower
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "sb.lower" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) throw RuntimeException("Unsupported method for a model with infinite-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                D.elementMax()

                val AN = D.elementPower(N).elementSum()
                val A1 = D.elementPower(1.0).elementSum()

                CN = Matrix(1, 1)
                CN.set(0, 0, Z + A1 + (N - 1) * FastMath.pow(AN / A1, 1 / (N - 1)))
                XN = Matrix(1, 1)
                XN.set(0, 0, N / CN.value())
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN undefined in the literature so we use ABA lower
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "gb.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                val Dmax = D.elementMax()
                XN = Matrix(1, 1)
                XN.set(0, 0, Maths.min(1 / Dmax, pfqn_xzgsbup(D, N, Z)))
                val ret = pfqn_xzgsblow(D, N, Z)
                CN = Matrix(1, 1)
                CN.set(0, 0, N / ret)
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                val XNlow = ret
                var k = 0
                RN = Matrix(sn.sched.size, 1)
                QN = Matrix(sn.sched.size, 1)
                run {
                    var i = 0
                    while (i < sn.sched.size) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i))
                            QN.set(i, 0, XN.value() * RN.get(i, 0))
                        } else {
                            QN.set(i, 0, pfqn_qzgbup(D, N, Z, k))
                            RN.set(i, 0, QN.get(i, 0) / XNlow / V.get(i))
                            k++ // increment after setting QN and RN because of 0-indexing
                        }
                        i++
                    }
                }/* RN(sn.schedid == SchedStrategy.ID_INF,1) = 1 ./ sn.rates(sn.schedid == SchedStrategy.ID_INF,1);
                 *   is redundant in LINE, RN is already set in the previous for-loop*/
                UN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < UN.numRows) {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                        i++
                    }
                }
                var i = 0
                while (i < sn.sched.size) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "gb.lower" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var Z = 0.0
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i)
                        } else {
                            nondelays++
                        }
                        i++
                    }
                }
                val D = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs.toDouble()
                XN = Matrix(1, 1)
                XN.set(0, 0, pfqn_xzgsblow(D, N, Z))
                val ret = pfqn_xzgsbup(D, N, Z)
                CN = Matrix(1, 1)
                CN.set(0, 0, N / ret)
                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                val XNup = ret
                var k = 0
                RN = Matrix(sn.sched.size, 1)
                QN = Matrix(sn.sched.size, 1)
                run {
                    var i = 0
                    while (i < sn.sched.size) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i))
                            QN.set(i, 0, XN.value() * RN.get(i, 0))
                        } else {
                            QN.set(i, 0, pfqn_qzgblow(D, N, Z, k))
                            RN.set(i, 0, QN.get(i, 0) / XNup / V.get(i))
                            k++ // increment after setting QN and RN because of 0-indexing
                        }
                        i++
                    }
                }/* RN(sn.schedid == SchedStrategy.ID_INF,1) = 1 ./ sn.rates(sn.schedid == SchedStrategy.ID_INF,1);
                 *   is redundant in LINE, RN is already set in the previous for-loop*/
                UN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < UN.numRows) {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                        i++
                    }
                }
                var i = 0
                while (i < sn.sched.size) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }

        "harel.lower", "harel.upper" -> {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                // Closed single-class queueing network - Harel et al. bounds
                run {
                    var i = 0
                    while (i < sn.nservers.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) throw RuntimeException(
                            "Unsupported method for a model with multi-server stations.")
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) throw RuntimeException(
                            "Harel bounds do not support models with think times (infinite-server stations).")
                        i++
                    }
                }
                val V: Matrix = sn.visits.get(0)!!
                var nondelays = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            nondelays++
                        }
                        i++
                    }
                }
                // Build loading vector rho = D (service demands) for stations
                val rho = Matrix(nondelays, 1)
                var idx = 0
                run {
                    var i = 0
                    while (i < V.numRows) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            rho.set(idx, 0, V.get(i) / sn.rates.get(i))
                            idx++
                        }
                        i++
                    }
                }
                val N = sn.nclosedjobs
                val Dmax = rho.elementMax()

                // Call the Harel bounds API
                val bounds = pfqn_harel_bounds(rho, N, 0.0, minOf(N, 7))

                // Select bound based on method
                val throughputBound = if (method == "harel.lower") {
                    bounds.LB
                } else {
                    // For upper bound, use the tightest available UB(n) which is UB(min(N,7))
                    val maxN = minOf(N, 7)
                    if (maxN >= 2) bounds.UB[maxN] else bounds.TH[1]
                }

                XN = Matrix(1, 1)
                XN.set(0, 0, if (method == "harel.upper") {
                    Maths.min(1 / Dmax, throughputBound)
                } else {
                    throughputBound
                })

                CN = Matrix(1, 1)
                CN.set(0, 0, N.toDouble() / XN.value())

                TN = Matrix(V.numRows, 1)
                run {
                    var i = 0
                    while (i < V.numRows) {
                        TN.set(i, 0, V.get(i, 0) * XN.value())
                        i++
                    }
                }
                // RN uses simple approximation
                RN = Matrix(sn.rates.numRows, 1)
                run {
                    var i = 0
                    while (i < RN.numRows) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i))
                        i++
                    }
                }
                QN = Matrix(TN.numRows, 1)
                run {
                    var i = 0
                    while (i < TN.numRows) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0))
                        i++
                    }
                }
                UN = Matrix(TN.numRows, 1)
                var i = 0
                while (i < TN.numRows) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0))
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i))
                    }
                    i++
                }
                lG = -N * FastMath.log(XN.value())
            }
            endTime = System.nanoTime()
        }
    }
    res.QN = QN
    res.UN = UN
    res.RN = RN
    res.TN = TN
    res.CN = CN
    res.XN = XN
    res.AN = Matrix(0, 0)
    res.WN = WN
    res.logNormConstAggr = lG
    res.runtime = (endTime - startTime) / 1000000000.0
    res.iter = iter
    res.method = method
    return res
}

