package jline.solvers.qns.handlers

import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.GlobalConstants
import jline.VerboseLevel
import jline.solvers.jmt.SolverJMT
import jline.solvers.qns.QNSResult
import jline.util.matrix.Matrix
import jline.io.lineTempName
import jline.api.sn.snGetDemandsChain
import jline.api.sn.snDeaggregateChainResults
import java.io.*
import java.nio.file.Files

/**
 * Core handler for the QNS solver
 */
class Solver_qns(private val sn: NetworkStruct, private val options: SolverOptions) {

    private var actualMethod: String = options.method ?: "default"

    /**
     * Solve the queueing network using the external qnsolver tool
     */
    fun solve(): QNSResult {
        val M = sn.nstations
        val K = sn.nclasses

        // Initialize result matrices
        val QN = Matrix.zeros(M, K)
        val UN = Matrix.zeros(M, K)
        val RN = Matrix.zeros(M, K)
        val TN = Matrix.zeros(M, K)
        val AN = Matrix.zeros(M, K)
        val WN = Matrix.zeros(M, K)
        val CN = Matrix.zeros(1, K)
        val XN = Matrix.zeros(1, K)

        // Track the actual method that will be used (line 30 in MATLAB)
        actualMethod = options.method ?: "default"

        // Map method to multiserver configuration (lines 32-41 in MATLAB)
        when (options.method) {
            "conway" -> options.config.multiserver = "conway"
            "reiser" -> options.config.multiserver = "reiser"
            "rolia" -> options.config.multiserver = "rolia"
            "zhou" -> options.config.multiserver = "zhou"
            "suri" -> options.config.multiserver = "suri"
            "schmidt" -> options.config.multiserver = "schmidt"
        }

        // Create temporary directory for files
        val tempDirPath = lineTempName("qns")
        val tempDir = File(tempDirPath)

        try {
            // Write model to JMVA format
            val modelFile = File(tempDir, "model.jmva")
            writeJMVAFile(modelFile)

            // Prepare result file
            val resultFile = File(tempDir, "result.jmva")
            val logFile = File(tempDir, "console.out")

            // Build command
            val cmd = buildCommand(modelFile, resultFile, logFile)

            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                println("SolverQNS command: $cmd")
            }

            // Execute command
            // Note: Output/error redirection is handled in the command string itself
            val exitCode = if (isWindows()) {
                val process = ProcessBuilder("cmd", "/c", cmd)
                    .directory(tempDir)
                    .redirectOutput(ProcessBuilder.Redirect.DISCARD)
                    .redirectError(ProcessBuilder.Redirect.DISCARD)
                    .start()
                process.waitFor()
            } else {
                val process = ProcessBuilder("sh", "-c", cmd)
                    .directory(tempDir)
                    .redirectOutput(ProcessBuilder.Redirect.DISCARD)
                    .redirectError(ProcessBuilder.Redirect.DISCARD)
                    .start()
                process.waitFor()
            }

            if (exitCode != 0) {
                val logContent = if (logFile.exists()) logFile.readText() else "No log file"
                throw RuntimeException("QNS solver failed with exit code: $exitCode\nLog: $logContent")
            }

            // Parse results from the output file
            val (Uchain, Qchain, Wchain, Tchain) = parseResults(resultFile, sn.nchains)

            // Convert chain results to station-class results
            val demandResults = snGetDemandsChain(sn)
            val Lchain = demandResults.Dchain
            val STchain = demandResults.STchain
            val Vchain = demandResults.Vchain
            val alpha = demandResults.alpha

            val Xchain = Matrix.zeros(1, sn.nchains)

            // Calculate system throughput for each chain
            // MATLAB: ref= zeros(sn.nchains,1); (line 132, not used)
            // MATLAB: Xchain(c)=Tchain(sn.refstat(c),c); (line 135)
            for (c in 0 until sn.nchains) {
                val refstat = sn.refstat.get(c).toInt() // sn.refstat is 0-based index
                Xchain.set(0, c, Tchain.get(refstat, c))
            }

            // Response times per chain (Rchain = Wchain for dollar output)
            val Rchain = Wchain.copy()
            Rchain.removeNaN()  // Set NaN to 0 like MATLAB: Rchain(isnan(Rchain))=0

            // Adjust utilizations for multi-server stations
            for (i in 0 until sn.nstations) {
                val servers = sn.nservers.get(i)
                if (!servers.isInfinite() && servers > 1) {
                    for (c in 0 until sn.nchains) {
                        Uchain.set(i, c, Uchain.get(i, c) / servers)
                    }
                }
            }

            // Deaggregate chain results to station-class results
            val results = snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha,
                Qchain, Uchain, Rchain, Tchain, null, Xchain)

            // Update result matrices
            QN.setTo(results.Q)
            UN.setTo(results.U)
            RN.setTo(results.R)
            TN.setTo(results.T)
            CN.setTo(results.C)
            XN.setTo(results.X)

            // WN is the same as RN for QNS solver
            WN.setTo(results.R)

            // AN will be calculated later in the analyzer

        } finally {
            // Clean up temporary directory
            deleteDirectory(tempDir)
        }

        return QNSResult(QN, UN, RN, TN, AN, WN, CN, XN, 0.0, actualMethod, 0)
    }

    /**
     * Write the network model to JMVA format
     */
    private fun writeJMVAFile(modelFile: File) {
        // Use the existing JMT solver functionality to write JMVA format
        SolverJMT.writeJMVA(sn, modelFile.absolutePath, options)
    }

    /**
     * Build the command line for the qnsolver
     */
    private fun buildCommand(modelFile: File, resultFile: File, logFile: File): String {
        val cmd = StringBuilder()
        cmd.append("qnsolver")
        cmd.append(" -l ").append(modelFile.absolutePath)

        // Add multiserver method if needed (lines 42-72 in MATLAB)
        if (hasMultiServer() && options.config.multiserver != null) {
            when (options.config.multiserver) {
                "default", "conway" -> {
                    cmd.append(" -mconway")
                    actualMethod = "conway"
                }
                "reiser" -> {
                    cmd.append(" -mreiser")
                    actualMethod = "reiser"
                }
                "rolia" -> {
                    cmd.append(" -mrolia")
                    actualMethod = "rolia"
                }
                "zhou" -> {
                    cmd.append(" -mzhou")
                    actualMethod = "zhou"
                }
                // Note: 'suri' and 'schmidt' are missing from MATLAB command building section
                // even though they are mapped in lines 32-41. Keeping for exact parity.
            }
        }

        cmd.append(" -o ").append(resultFile.absolutePath)
        cmd.append(" > ").append(logFile.absolutePath).append(" 2>&1")

        return cmd.toString()
    }

    /**
     * Check if the model has multi-server stations
     */
    private fun hasMultiServer(): Boolean {
        for (i in 0 until sn.nstations) {
            val servers = sn.nservers.get(i).toInt()
            if (servers > 1 && servers != Int.MAX_VALUE) {
                return true
            }
        }
        return false
    }

    /**
     * Check if running on Windows
     */
    private fun isWindows(): Boolean {
        return System.getProperty("os.name").lowercase().contains("win")
    }

    /**
     * Parse the QNS solver results from the output file
     */
    private fun parseResults(resultFile: File, nchains: Int): ParsedResults {
        // Initialize matrices with zeros
        val Uchain = Matrix(sn.nstations, nchains)
        val Qchain = Matrix(sn.nstations, nchains)
        val Wchain = Matrix(sn.nstations, nchains)
        val Tchain = Matrix(sn.nstations, nchains)

        Uchain.fill(0.0)
        Qchain.fill(0.0)
        Wchain.fill(0.0)
        Tchain.fill(0.0)

        if (!resultFile.exists()) {
            throw RuntimeException("QNS result file not found: ${resultFile.absolutePath}")
        }

        resultFile.forEachLine { line ->
            if (line.contains(",") && !line.contains("$")) {
                // Use nclasses to determine output format - matches MATLAB solver_qns.m line 110
                // The aggregate column presence in qnsolver output depends on nclasses, not nchains
                val parsed = if (sn.nclasses == 1) {
                    parseDollarOutputSingleClass(line, nchains)
                } else {
                    parseDollarOutput(line, nchains)
                }

                if (parsed != null) {
                    // Find which station index matches the parsed station name
                    var stationIdx = -1
                    for (i in 0 until sn.nstations) {
                        val nodeIdx = sn.stationToNode.get(i).toInt()
                        val stationNodeName = sn.nodenames[nodeIdx]
                        if (stationNodeName == parsed.statName) {
                            stationIdx = i
                            break
                        }
                    }

                    // If found, populate the matrices
                    if (stationIdx != -1) {
                        for (c in 0 until nchains) {
                            Qchain.set(stationIdx, c, parsed.Q[c])
                            Wchain.set(stationIdx, c, parsed.W[c])
                            Uchain.set(stationIdx, c, parsed.U[c])
                            Tchain.set(stationIdx, c, parsed.T[c])
                        }
                    }
                }
            }
        }

        return ParsedResults(Uchain, Qchain, Wchain, Tchain)
    }

    /**
     * Parse a single line of output for multi-class models
     * Format: Station, $Q(Chain01), $Q(Chain02), $Q, $R(Chain01), $R(Chain02), $R, $U(Chain01), $U(Chain02), $U, $X(Chain01), $X(Chain02), $X
     * Note: R values map to W (residence time), X values map to T (throughput)
     */
    private fun parseDollarOutput(line: String, nchains: Int): ParsedLine? {
        val parts = line.replace(" ", "").split(",")
        if (parts.size < 1 + 4 * (nchains + 1)) return null

        val statName = parts[0]
        val Q = DoubleArray(nchains)
        val W = DoubleArray(nchains)  // Will hold R values from output
        val U = DoubleArray(nchains)
        val T = DoubleArray(nchains)  // Will hold X values from output

        var ptr = 1
        // $Q values
        for (r in 0 until nchains) {
            Q[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += nchains + 1 // skip aggregate $Q value

        // $R values (map to W - residence times)
        for (r in 0 until nchains) {
            W[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += nchains + 1 // skip aggregate $R value

        // $U values
        for (r in 0 until nchains) {
            U[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += nchains + 1 // skip aggregate $U value

        // $X values (map to T - throughputs)
        for (r in 0 until nchains) {
            T[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }

        return ParsedLine(statName, Q, W, U, T)
    }

    /**
     * Parse a single line of output for single-class models
     * For single class, the output format doesn't have aggregate values
     */
    private fun parseDollarOutputSingleClass(line: String, nchains: Int): ParsedLine? {
        val parts = line.replace(" ", "").split(",")
        if (parts.size < 5) return null

        val statName = parts[0]
        val Q = DoubleArray(nchains)
        val W = DoubleArray(nchains)
        val U = DoubleArray(nchains)
        val T = DoubleArray(nchains)

        var ptr = 1
        // Q values
        for (r in 0 until nchains) {
            Q[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += 1  // No aggregate for single class

        // R values (map to W)
        for (r in 0 until nchains) {
            W[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += 1

        // U values
        for (r in 0 until nchains) {
            U[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }
        ptr += 1

        // X values (map to T)
        for (r in 0 until nchains) {
            T[r] = parts[ptr + r].toDoubleOrNull() ?: 0.0
        }

        return ParsedLine(statName, Q, W, U, T)
    }


    /**
     * Delete a directory and all its contents
     */
    private fun deleteDirectory(dir: File): Boolean {
        if (dir.isDirectory) {
            val children = dir.listFiles()
            children?.forEach { child ->
                deleteDirectory(child)
            }
        }
        return dir.delete()
    }

    // Data classes for results
    private data class ParsedLine(
        val statName: String,
        val Q: DoubleArray,
        val W: DoubleArray,
        val U: DoubleArray,
        val T: DoubleArray
    )

    private data class ParsedResults(
        val Uchain: Matrix,
        val Qchain: Matrix,
        val Wchain: Matrix,
        val Tchain: Matrix
    )

}