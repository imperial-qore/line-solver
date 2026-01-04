package jline.solvers.ssa.analyzers

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAResult
import jline.solvers.ssa.handlers.solver_ssa_nrm_space
import jline.util.matrix.Matrix

fun solver_ssa_analyzer_nrm_space(sn: NetworkStruct, options: SolverOptions): SSAResult {
    val M = sn.nstations
    val K = sn.nclasses
    val S = sn.nservers
    val NK = sn.njobs.transpose()
    val PH = sn.proc
    val sched = sn.sched

    val XN = Matrix(1, K).apply { fill(Double.NaN) }
    val UN = Matrix(M, K).apply { fill(Double.NaN) }
    val QN = Matrix(M, K).apply { fill(Double.NaN) }
    val RN = Matrix(M, K).apply { fill(Double.NaN) }
    val TN = Matrix(M, K).apply { fill(Double.NaN) }
    val CN = Matrix(1, K).apply { fill(Double.NaN) }

    val result = solver_ssa_nrm_space(sn, options)
    var pi = result.pi
    val space = result.outspace
    val depRates = result.depRates

    if (pi.numRows > 1) {
        pi = pi.transpose()
    }

    for (k in 0 until K) {
        val refnd = sn.stationToNode.get(sn.refstat.get(k).toInt()).toInt()
        val colView = depRates.getColumnView(refnd * K + k)
        XN.set(0, k, pi.multColumnView(colView))
    }

    val rates = sn.rates
    for (ist in 0 until M) {
        val ind = sn.stationToNode.get(ist).toInt()
        for (k in 0 until K) {
            val dcolView = depRates.getColumnView(ind * K + k)
            val scolView = space.getColumnView(ind * K + k)
            TN.set(ist, k, pi.multColumnView(dcolView))
            QN.set(ist, k, pi.multColumnView(scolView))
        }

        when (sched.get(sn.stations.get(ist))) {
            SchedStrategy.INF, SchedStrategy.EXT -> {
                for (k in 0 until K) {
                    UN.set(ist, k, QN.get(ist, k))
                }
            }
            SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
            SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                for (k in 0 until K) {
                    val phEntry = PH.get(sn.stations.get(ist))?.get(sn.jobclasses.get(k))
                    if (phEntry != null && !phEntry.isEmpty) {
                        val dcolView = depRates.getColumnView(ind * K + k)
                        val throughput = pi.multColumnView(dcolView)
                        UN.set(ist, k, throughput / rates.get(ist, k) / S.get(ist))
                    }
                }
            }
            else -> {}
        }
    }

    for (k in 0 until K) {
        for (ist in 0 until M) {
            if (TN.get(ist, k) > 0) {
                RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
            } else {
                RN.set(ist, k, 0.0)
            }
        }
        CN.set(0, k, NK.get(k) / XN.get(0, k))
    }

    QN.apply(Double.NaN, 0.0, "equal")
    UN.apply(Double.NaN, 0.0, "equal")
    RN.apply(Double.NaN, 0.0, "equal")
    XN.apply(Double.NaN, 0.0, "equal")
    TN.apply(Double.NaN, 0.0, "equal")
    CN.apply(Double.NaN, 0.0, "equal")

    return SSAResult(QN, UN, RN, TN, CN, XN, null, null, sn)
}