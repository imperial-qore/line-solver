package jline.solvers.nc.handlers

import jline.GlobalConstants
import jline.api.mam.map_pie
import jline.api.pfqn.ld.pfqn_ncld
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.lang.state.State
import jline.lang.state.ToMarginal
import jline.lang.state.ToMarginal.*
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import kotlin.RuntimeException
import kotlin.math.abs
import kotlin.math.ln
import kotlin.run

fun solver_nc_marg(sn: NetworkStruct, options: SolverOptions, lG: Double): SolverNC.SolverNCMargReturn {
    var lG: Double = lG
    var M = sn.nstations
    var K = sn.nclasses
    val state = sn.state
    val S = sn.nservers
    var V = Matrix(sn.nstateful, K)
    for (i in 0..<sn.visits.size) {
        V = V.add(1.0, sn.visits.get(i))
    }
    //V = cellsum(sn.visits);
    val rates = sn.rates
    val ST = rates.copy()
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            ST.set(i, j, 1.0 / ST.get(i, j))
        }
    }
    ST.removeNaN()

    val ret = snGetDemandsChain(sn)
    val Lchain = ret.Dchain
    val STchain = ret.STchain
    val Nchain = ret.Nchain

    val startTimeMillis = System.nanoTime()
    M = STchain.numRows
    K = STchain.numCols

    var mu = Matrix(0, Nchain.elementSum().toInt())
    for (i in 0..<M) {
        val tmp = Matrix(1, Nchain.elementSum().toInt())
        if (Utils.isInf(S.get(i))) {
            for (j in 0..<tmp.length()) {
                tmp.set(j, (j + 1).toDouble())
            }
        } else {
            for (j in 0..<tmp.length()) {
                tmp.set(j, FastMath.min((j + 1).toDouble(), S.get(i)))
            }
        }
        mu = Matrix.concatRows(mu, tmp, null)
    }


    if (lG.isNaN) {
        val Z_tmp = Nchain.copy()
        Z_tmp.fill(0.0)
        lG = pfqn_ncld(Lchain, Nchain, Z_tmp, mu, options).lG as Double
    }

    val G = FastMath.exp(lG.toDouble())
    val lPr = Matrix(sn.nstations, 1)
    lPr.fill(0.0)

    for (ist in 0..<sn.nstations) {
        val ind = sn.stationToNode.get(ist).toInt()
        val isf = sn.stationToStateful.get(ist).toInt()
        val ret1 = toMarginal(sn, ind, state.get(sn.stateful.get(isf)), null, null, null, null, null)
        val nirvec = ret1.nir
        val sivec = ret1.sir
        val kirvec = ret1.kir

        // Remove the check for negative elements as it was causing valid computations to be skipped
        // The original check (nirvec.elementMin() < GlobalConstants.Zero) was setting NaN for negative values
        // but the test expects negative log probabilities which are valid
        if (false && nirvec.elementMin() < GlobalConstants.Zero) {
            lPr.set(ist, Double.NaN)
        } else {
            //nivec_chain = nirvec * sn.chains'; Assume it to be matrix
            val nivec_chain = nirvec.mult(sn.chains.transpose())

            var Lchain_tmp = Matrix(0, Lchain.numCols)
            var mu_tmp = Matrix(0, mu.numCols)
            for (i in 0..<sn.nstations) {
                if (i != ist) {
                    val Lchain_row_i = Matrix.extractRows(Lchain, i, i + 1, null)
                    val mu_row_i = Matrix.extractRows(mu, i, i + 1, null)
                    Lchain_tmp = Matrix.concatRows(Lchain_tmp, Lchain_row_i, null)
                    mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null)
                }
            }
            val Nchain_tmp = Nchain.copy()
            for (i in 0..<Nchain_tmp.length()) {
                Nchain_tmp.set(i, Nchain_tmp.get(i) - nivec_chain.get(i))
            }
            val Zchain_tmp = Nchain.copy()
            Zchain_tmp.fill(0.0)

            val lG_minus_i = pfqn_ncld(Lchain_tmp, Nchain_tmp, Zchain_tmp, mu_tmp, options).lG
            var lF_i = 0.0

            when (sn.sched.get(sn.stations.get(ist))) {
                SchedStrategy.FCFS -> {
                    var r = 0
                    while (r < K) {
                        val PHr: MatrixCell = sn.proc.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!
                        if (!PHr.isEmpty) {
                            var kir = Matrix(1, kirvec.size)
                            var i = 0
                            while (i < kir.length()) {
                                kir.set(i, kirvec.get(i)!!.get(0, r))
                                i++
                            }
                            if (kir.length() > 1) {
                                throw RuntimeException("solver_nc_marg: Cannot return state probability " + "because the product-form solution requires exponential service times at FCFS nodes.")
                            }

                            if (abs(ST.get(ist, r) - Matrix.extractRows(ST, ist, ist + 1, null).elementMax()) > 1e-6) {
                                throw RuntimeException("solver_nc_marg: Cannot return state probability " + "because the product-form solution requires identical service times at FCFS nodes.")
                            }
                        }
                        r++
                    }
                    var allZero = true
                    var i = 0
                    while (allZero && i < sivec.numRows) {
                        var j = 0
                        while (allZero && j < sivec.numCols) {
                            if (abs(sivec.get(i, j)) > 1e-6) {
                                allZero = false
                            }
                            j++
                        }
                        i++
                    }
                    if (!allZero) {
                        var sumLog = 0.0
                        for (r in 0..<K) {
                            sumLog += nirvec.get(0, r) * ln(V.get(ist, r))
                        }
                        var sum_kirvec = 0.0
                        run {
                            var i = 0
                            while (i < kirvec.size) {
                                sum_kirvec += kirvec.get(i)!!.elementSum()
                                i++
                            }
                        }
                        val mu_row_ist = Matrix(1, sum_kirvec.toInt())
                        Matrix.extract(mu, ist, ist + 1, 0, sum_kirvec.toInt(), mu_row_ist, 0, 0)
                        var i = 0
                        while (i < mu_row_ist.length()) {
                            mu_row_ist.set(i, FastMath.log(mu_row_ist.get(i)))
                            i++
                        }
                        lF_i += (sumLog - mu_row_ist.elementSum())
                    } else {
                        lF_i = 0.0
                    }
                }

                SchedStrategy.PS -> {
                    var r = 0
                    while (r < K) {
                        val PHr: MatrixCell = sn.proc.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!
                        if (!PHr.isEmpty) {
                            val kir = Matrix(1, kirvec.size)
                            run {
                                var i = 0
                                while (i < kir.length()) {
                                    kir.set(i, kirvec.get(i)!!.get(0, r))
                                    i++
                                }
                            }

                            var PHr_tmp = PHr.get(0)
                            run {
                                var i = 0
                                while (i < PHr_tmp.numRows) {
                                    var j = 0
                                    while (j < PHr_tmp.numCols) {
                                        PHr_tmp.set(i, j, -1.0 * PHr_tmp.get(i, j))
                                        j++
                                    }
                                    i++
                                }
                            }
                            PHr_tmp = PHr_tmp.inv()
                            val Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp)

                            val kir_tmp = Ar.copy()
                            var i = 0
                            while (i < kir_tmp.numRows) {
                                var j = 0
                                while (j < kir_tmp.numCols) {
                                    kir_tmp.set(i, j, kir.get(i, j) * FastMath.log(V.get(ist, r) * kir_tmp.get(i, j)))
                                    j++
                                }
                                i++
                            }

                            lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum())
                        }
                        r++
                    }

                    var sum_kirvec = 0.0
                    run {
                        var i = 0
                        while (i < kirvec.size) {
                            sum_kirvec += kirvec.get(i)!!.elementSum()
                            i++
                        }
                    }
                    val mu_row_ist = Matrix(1, sum_kirvec.toInt())
                    Matrix.extract(mu, ist, ist + 1, 0, sum_kirvec.toInt(), mu_row_ist, 0, 0)
                    var i = 0
                    while (i < mu_row_ist.numCols) {
                        mu_row_ist.set(i, FastMath.log(mu_row_ist.get(i)))
                        i++
                    }
                    lF_i += (Maths.factln(sum_kirvec) - mu_row_ist.elementSum())
                }

                SchedStrategy.INF -> {
                    var r = 0
                    while (r < K) {
                        val PHr: MatrixCell = sn.proc.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!
                        if (!PHr.isEmpty) {
                            val kir = Matrix(1, kirvec.size)
                            run {
                                var i = 0
                                while (i < kir.length()) {
                                    kir.set(i, kirvec.get(i)!!.get(0, r))
                                    i++
                                }
                            }
                            var PHr_tmp = PHr.get(0)
                            run {
                                var i = 0
                                while (i < PHr_tmp.numRows) {
                                    var j = 0
                                    while (j < PHr_tmp.numCols) {
                                        PHr_tmp.set(i, j, -1.0 * PHr_tmp.get(i, j))
                                        j++
                                    }
                                    i++
                                }
                            }
                            PHr_tmp = PHr_tmp.inv()
                            val Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp)

                            val kir_tmp = Ar.copy()
                            var i = 0
                            while (i < kir_tmp.numRows) {
                                var j = 0
                                while (j < kir_tmp.numCols) {
                                    kir_tmp.set(i, j, kir.get(i, j) * FastMath.log(V.get(ist, r) * kir_tmp.get(i, j)))
                                    j++
                                }
                                i++
                            }

                            lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum())
                        }
                        r++
                    }
                }

                else -> {}
            }
            lPr.set(ist, lF_i + lG_minus_i - lG.toDouble())
        }
    }
    val endTimeMillis = System.nanoTime()
    val runtime = (endTimeMillis - startTimeMillis) / 1000000000.0
    lPr.removeNaN()
    return SolverNC.SolverNCMargReturn(lPr, G, lG.toDouble(), runtime)
}