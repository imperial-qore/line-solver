package jline.solvers.qns.analyzers

import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.NetworkSolver
import jline.solvers.qns.QNSResult
import jline.solvers.qns.SolverQNS
import jline.solvers.qns.handlers.Solver_qns
import jline.util.matrix.Matrix
import jline.api.sn.snGetArvRFromTput

/**
 * Analyzer for the QNS solver
 */
class Solver_qns_analyzer(private val solver: SolverQNS) {

    /**
     * Run the QNS analyzer
     */
    fun runAnalyzer(): QNSResult {
        val startTime = System.nanoTime()

        val sn = solver.sn
        val options = solver.options

        // Validate options
        if (options.method == null) {
            options.method = "default"
        }

        // Map method to multiserver approximation if needed
        when (options.method) {
            "conway" -> options.config.multiserver = "conway"
            "rolia" -> options.config.multiserver = "rolia"
            "zhou" -> options.config.multiserver = "zhou"
            "suri" -> options.config.multiserver = "suri"
            "reiser" -> options.config.multiserver = "reiser"
            "schmidt" -> options.config.multiserver = "schmidt"
            "default" -> options.config.multiserver = "default"
        }

        // Run the core solver
        val handler = Solver_qns(sn, options)
        val result = handler.solve()

        // Calculate arrival rates from throughputs
        val T = solver.getAvgTputHandles()
        val AN = snGetArvRFromTput(sn, result.TN, T)

        val endTime = System.nanoTime()
        val runtime = (endTime - startTime) / 1e9

        // Determine the actual method used
        val actualMethod = when (options.method) {
            "default" -> "default/${result.method}"
            else -> result.method
        }

        return QNSResult(
            result.QN,
            result.UN,
            result.RN,
            result.TN,
            AN,
            result.WN,
            result.CN,
            result.XN,
            runtime,
            actualMethod,
            0  // QNS doesn't have iterations
        )
    }
}