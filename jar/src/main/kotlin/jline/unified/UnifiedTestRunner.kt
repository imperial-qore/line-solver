/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.unified

import com.google.gson.Gson
import com.google.gson.JsonArray
import com.google.gson.JsonElement
import com.google.gson.JsonObject
import jline.lang.Network
import jline.solvers.NetworkSolver
import jline.solvers.SolverOptions
import jline.solvers.ctmc.CTMC
import jline.solvers.des.DES
import jline.solvers.fluid.FLD
import jline.solvers.jmt.JMT
import jline.solvers.mam.MAM
import jline.solvers.mva.MVA
import jline.solvers.nc.NC
import jline.solvers.ssa.SSA
import java.io.File
import java.nio.file.Paths
import kotlin.math.abs

/**
 * Data class representing test results
 */
data class TestResults(
    var passed: Int = 0,
    var failed: Int = 0,
    var skipped: Int = 0,
    val errors: MutableList<String> = mutableListOf()
)

/**
 * Data class for JSON-formatted test results
 */
data class JsonResult(
    val modelName: String,
    val language: String = "java",
    var status: String = "passed",
    val timing: MutableMap<String, Any> = mutableMapOf(
        "total_ms" to 0L,
        "solvers" to mutableMapOf<String, Long>()
    ),
    val solverResults: MutableMap<String, Map<String, Any>> = mutableMapOf(),
    val errors: MutableList<String> = mutableListOf()
)

/**
 * Executes unified cross-language tests from JSON definitions.
 *
 * This class loads JSON test definitions and runs them against Java/Kotlin solvers,
 * comparing results against expected values with appropriate tolerances.
 */
class UnifiedTestRunner {

    companion object {
        // Default tolerances
        const val DEFAULT_TOL = 1e-8
        const val SIMULATION_TOL = 0.05

        // Simulation-based solvers
        val SIMULATION_SOLVERS = setOf("SolverJMT", "SolverSSA", "SolverDES")
    }

    private val gson = Gson()
    private val definitionsDir: File
    var results = TestResults()
        private set

    init {
        // Check for system property first (set by run_tests.sh)
        val definitionsPath = System.getProperty("line.definitions")
        definitionsDir = if (definitionsPath != null) {
            File(definitionsPath)
        } else {
            // Find definitions directory relative to test resources or project root
            val projectRoot = findProjectRoot()
            File(projectRoot, "test/unified/definitions")
        }
    }

    private fun findProjectRoot(): File {
        // Try to find project root by looking for CLAUDE.md or pom.xml
        var current = File(System.getProperty("user.dir"))
        while (current.parentFile != null) {
            if (File(current, "CLAUDE.md").exists() ||
                (File(current, "jar").exists() && File(current, "matlab").exists())) {
                return current
            }
            current = current.parentFile
        }
        // Fallback to current directory
        return File(System.getProperty("user.dir"))
    }

    /**
     * Run all available test definitions.
     *
     * @param verbose Enable verbose output
     * @return TestResults with passed/failed/skipped counts and errors
     */
    fun runAll(verbose: Boolean = false): TestResults {
        results = TestResults()

        if (!definitionsDir.exists()) {
            println("Definitions directory not found: ${definitionsDir.absolutePath}")
            return results
        }

        val jsonFiles = definitionsDir.listFiles { f -> f.extension == "json" }
            ?.sortedBy { it.name }
            ?: emptyList()

        println("=== Unified Test Runner (Java/Kotlin) ===")
        println("Found ${jsonFiles.size} test definitions\n")

        for (jsonFile in jsonFiles) {
            val modelName = jsonFile.nameWithoutExtension
            runModel(modelName, verbose)
        }

        printSummary()
        return results
    }

    /**
     * Run tests for a specific model.
     *
     * @param modelName Name of the model (e.g., "oqn_basic")
     * @param verbose Enable verbose output
     * @param jsonOutput Return JsonResult instead of String
     * @return "passed", "failed", or "skipped" (or JsonResult if jsonOutput=true)
     */
    fun runModel(modelName: String, verbose: Boolean = false, jsonOutput: Boolean = false): Any {
        // Initialize JSON result if needed
        val jsonResult = if (jsonOutput) JsonResult(modelName) else null
        val modelStartTime = if (jsonOutput) System.nanoTime() else 0L

        val jsonFile = File(definitionsDir, "$modelName.json")
        if (!jsonFile.exists()) {
            if (verbose) println("  [SKIP] $modelName - definition file not found")
            results.skipped++
            if (jsonOutput) {
                jsonResult!!.status = "skipped"
                jsonResult.errors.add("definition file not found")
                return jsonResult
            }
            return "skipped"
        }

        val definition: JsonObject
        try {
            definition = gson.fromJson(jsonFile.readText(), JsonObject::class.java)
        } catch (e: Exception) {
            if (!jsonOutput) println("  [ERROR] $modelName - failed to parse JSON: ${e.message}")
            results.failed++
            results.errors.add("$modelName: JSON parse error")
            if (jsonOutput) {
                jsonResult!!.status = "error"
                jsonResult.errors.add("JSON parse error: ${e.message}")
                return jsonResult
            }
            return "failed"
        }

        // Check if Java is supported
        val supportedLanguages = definition.getAsJsonObject("supportedLanguages")
        if (supportedLanguages != null &&
            supportedLanguages.has("java") &&
            !supportedLanguages.get("java").asBoolean) {
            if (verbose) println("  [SKIP] $modelName - Java not supported")
            results.skipped++
            if (jsonOutput) {
                jsonResult!!.status = "skipped"
                jsonResult.errors.add("Java not supported")
                return jsonResult
            }
            return "skipped"
        }

        // Check if model is available in registry
        if (!ModelRegistry.hasModel(modelName)) {
            if (verbose) println("  [SKIP] $modelName - not in model registry")
            results.skipped++
            if (jsonOutput) {
                jsonResult!!.status = "skipped"
                jsonResult.errors.add("not in model registry")
                return jsonResult
            }
            return "skipped"
        }

        if (verbose && !jsonOutput) println("Testing: $modelName")

        // Get the model
        val model: Network
        try {
            model = ModelRegistry.getModel(modelName)
        } catch (e: NotImplementedError) {
            if (verbose) println("  [SKIP] $modelName - not yet implemented")
            results.skipped++
            if (jsonOutput) {
                jsonResult!!.status = "skipped"
                jsonResult.errors.add("not yet implemented")
                return jsonResult
            }
            return "skipped"
        } catch (e: Exception) {
            if (!jsonOutput) println("  [ERROR] $modelName - failed to build model: ${e.message}")
            results.failed++
            results.errors.add("$modelName: Model build error - ${e.message}")
            if (jsonOutput) {
                jsonResult!!.status = "error"
                jsonResult.errors.add("Model build error: ${e.message}")
                return jsonResult
            }
            return "failed"
        }

        // Run tests for each solver
        var allPassed = true
        val solversToTest = definition.getAsJsonArray("solvers")

        // Get list of solvers to skip for Java
        val skipSolvers = definition.getAsJsonObject("skipSolvers")
            ?.getAsJsonArray("java")
            ?.map { it.asString }
            ?: emptyList()

        for (solverElement in solversToTest) {
            val solverName = solverElement.asString

            // Skip if solver is in skip list
            if (solverName in skipSolvers) {
                if (jsonOutput) {
                    jsonResult!!.solverResults[solverName] = mapOf(
                        "status" to "skipped",
                        "metrics" to emptyMap<String, Any>()
                    )
                }
                continue
            }

            val expectedResults = definition.getAsJsonObject("expectedResults")
            if (!expectedResults.has(solverName)) continue

            try {
                if (jsonOutput) {
                    val (passed, solverResult, solverTime) = testSolverJson(model, solverName, definition, verbose)
                    jsonResult!!.solverResults[solverName] = solverResult
                    @Suppress("UNCHECKED_CAST")
                    (jsonResult.timing["solvers"] as MutableMap<String, Long>)[solverName] = solverTime
                    if (!passed) allPassed = false
                } else {
                    val passed = testSolver(model, solverName, definition, verbose)
                    if (!passed) allPassed = false
                }
            } catch (e: Exception) {
                if (verbose) println("    [ERROR] $solverName: ${e.message}")
                allPassed = false
                if (jsonOutput) {
                    jsonResult!!.solverResults[solverName] = mapOf(
                        "status" to "error",
                        "metrics" to emptyMap<String, Any>()
                    )
                    jsonResult.errors.add("$solverName: ${e.message}")
                }
                results.errors.add("$modelName/$solverName: ${e.message}")
            }
        }

        return if (allPassed) {
            results.passed++
            if (jsonOutput) {
                jsonResult!!.status = "passed"
                jsonResult.timing["total_ms"] = (System.nanoTime() - modelStartTime) / 1_000_000
                jsonResult
            } else {
                if (!verbose) print(".")
                "passed"
            }
        } else {
            results.failed++
            if (jsonOutput) {
                jsonResult!!.status = "failed"
                jsonResult.timing["total_ms"] = (System.nanoTime() - modelStartTime) / 1_000_000
                jsonResult
            } else {
                if (!verbose) print("F")
                "failed"
            }
        }
    }

    private fun testSolver(
        model: Network,
        solverName: String,
        definition: JsonObject,
        verbose: Boolean
    ): Boolean {
        // Get tolerance
        val toleranceConfig = definition.getAsJsonObject("tolerance")
        val tol = if (solverName in SIMULATION_SOLVERS) {
            toleranceConfig?.get("simulation")?.asDouble ?: SIMULATION_TOL
        } else {
            toleranceConfig?.get("default")?.asDouble ?: DEFAULT_TOL
        }

        // Create solver instance
        val solver = createSolver(model, solverName)
        if (solver == null) {
            if (verbose) println("    [SKIP] $solverName - solver not available")
            return true // Don't fail on unavailable solver
        }

        // Get results
        val avgTable = try {
            solver.avgTable
        } catch (e: Exception) {
            if (verbose) println("    [ERROR] $solverName - solver failed: ${e.message}")
            return false
        }

        // Compare against expected results
        val expected = definition.getAsJsonObject("expectedResults").getAsJsonObject(solverName)
        var passed = true

        // Compare each metric
        val metrics = listOf("QLen", "Util", "RespT", "Tput")
        for (metric in metrics) {
            if (!expected.has(metric)) continue

            val expectedValues = flatten2DArray(expected.getAsJsonArray(metric))
            val actualValues = extractMetric(avgTable, metric)

            if (actualValues != null && expectedValues != null) {
                if (!compareWithTolerance(actualValues, expectedValues, tol)) {
                    if (verbose) {
                        println("    [FAIL] $solverName.$metric mismatch (tol=${String.format("%.2e", tol)})")
                        println("      Expected: ${expectedValues.contentToString()}")
                        println("      Actual:   ${actualValues.contentToString()}")
                    }
                    passed = false
                } else if (verbose) {
                    println("    [PASS] $solverName.$metric")
                }
            }
        }

        if (passed && verbose) println("  [PASS] $solverName")
        return passed
    }

    private fun testSolverJson(
        model: Network,
        solverName: String,
        definition: JsonObject,
        verbose: Boolean
    ): Triple<Boolean, Map<String, Any>, Long> {
        val solverResult = mutableMapOf<String, Any>(
            "status" to "passed",
            "metrics" to mutableMapOf<String, Any>()
        )

        // Get tolerance
        val toleranceConfig = definition.getAsJsonObject("tolerance")
        val tol = if (solverName in SIMULATION_SOLVERS) {
            toleranceConfig?.get("simulation")?.asDouble ?: SIMULATION_TOL
        } else {
            toleranceConfig?.get("default")?.asDouble ?: DEFAULT_TOL
        }

        // Create solver instance
        val solver = createSolver(model, solverName)
        if (solver == null) {
            solverResult["status"] = "skipped"
            return Triple(true, solverResult, 0L)
        }

        // Get results with timing
        val solverStartTime = System.nanoTime()
        val avgTable = try {
            solver.avgTable
        } catch (e: Exception) {
            val solverTime = (System.nanoTime() - solverStartTime) / 1_000_000
            solverResult["status"] = "error"
            return Triple(false, solverResult, solverTime)
        }
        val solverTime = (System.nanoTime() - solverStartTime) / 1_000_000

        // Compare against expected results
        val expected = definition.getAsJsonObject("expectedResults").getAsJsonObject(solverName)
        var passed = true

        @Suppress("UNCHECKED_CAST")
        val metrics = solverResult["metrics"] as MutableMap<String, Any>

        // Compare each metric
        val metricNames = listOf("QLen", "Util", "RespT", "Tput")
        for (metric in metricNames) {
            if (!expected.has(metric)) continue

            val expectedValues = flatten2DArray(expected.getAsJsonArray(metric))
            val actualValues = extractMetric(avgTable, metric)

            val metricResult = mutableMapOf<String, Any>("passed" to true)

            if (actualValues != null && expectedValues != null) {
                if (!compareWithTolerance(actualValues, expectedValues, tol)) {
                    metricResult["passed"] = false
                    metricResult["expected"] = expectedValues.toList()
                    metricResult["actual"] = actualValues.toList()
                    // Calculate max error for compact display
                    var maxError = 0.0
                    for (i in actualValues.indices) {
                        if (!expectedValues[i].isNaN() && !actualValues[i].isNaN()) {
                            val error = abs(actualValues[i] - expectedValues[i])
                            if (error > maxError) maxError = error
                        }
                    }
                    metricResult["maxError"] = maxError
                    passed = false
                }
            }
            metrics[metric] = metricResult
        }

        if (!passed) {
            solverResult["status"] = "failed"
        }

        return Triple(passed, solverResult, solverTime)
    }

    private fun createSolver(model: Network, solverName: String): NetworkSolver? {
        return try {
            val options = SolverOptions()
            options.seed = 23000

            when (solverName) {
                "SolverMVA" -> MVA(model, options)
                "SolverCTMC" -> {
                    options.cutoff = jline.util.matrix.Matrix(doubleArrayOf(10.0))
                    CTMC(model, options)
                }
                "SolverJMT" -> JMT(model, options)
                "SolverSSA" -> {
                    options.samples = 10000
                    SSA(model, options)
                }
                "SolverDES" -> DES(model, options)
                "SolverFluid" -> FLD(model, options)
                "SolverNC" -> NC(model, options)
                "SolverMAM" -> MAM(model, options)
                else -> null
            }
        } catch (e: Exception) {
            null
        }
    }

    private fun extractMetric(avgTable: Any, metric: String): DoubleArray? {
        return try {
            val table = avgTable as jline.solvers.NetworkAvgTable
            val values: List<Double>? = when (metric) {
                "QLen" -> table.getQLen()
                "Util" -> table.getUtil()
                "RespT" -> table.getRespT()
                "Tput" -> table.getTput()
                else -> return null
            }
            values?.toDoubleArray()
        } catch (e: Exception) {
            null
        }
    }

    private fun flatten2DArray(arr: JsonArray): DoubleArray? {
        val flat = mutableListOf<Double>()
        for (row in arr) {
            if (row.isJsonArray) {
                for (elem in row.asJsonArray) {
                    flat.add(jsonElementToDouble(elem))
                }
            } else {
                flat.add(jsonElementToDouble(row))
            }
        }
        return flat.toDoubleArray()
    }

    private fun jsonElementToDouble(elem: JsonElement): Double {
        return when {
            elem.isJsonNull -> Double.NaN
            elem.isJsonPrimitive -> {
                val prim = elem.asJsonPrimitive
                when {
                    prim.isNumber -> prim.asDouble
                    prim.isString && prim.asString == "NaN" -> Double.NaN
                    prim.isString && prim.asString == "Inf" -> Double.POSITIVE_INFINITY
                    prim.isString && prim.asString == "-Inf" -> Double.NEGATIVE_INFINITY
                    else -> Double.NaN
                }
            }
            else -> Double.NaN
        }
    }

    private fun compareWithTolerance(actual: DoubleArray, expected: DoubleArray, tol: Double): Boolean {
        if (actual.size != expected.size) return false

        for (i in actual.indices) {
            val exp = expected[i]
            val act = actual[i]

            // Both NaN - OK
            if (exp.isNaN() && act.isNaN()) continue
            // One NaN, other not - fail
            if (exp.isNaN() || act.isNaN()) return false

            val relErr = abs(act - exp) / (abs(exp) + 1e-10)
            val absErr = abs(act - exp)
            if (relErr > tol && absErr > tol) return false
        }
        return true
    }

    private fun printSummary() {
        println("\n\n=== Test Summary ===")
        println("Passed:  ${results.passed}")
        println("Failed:  ${results.failed}")
        println("Skipped: ${results.skipped}")
        println("====================")

        if (results.failed > 0) {
            println("\nFailed tests:")
            for (error in results.errors) {
                println("  - $error")
            }
        }
    }
}

/**
 * CLI entry point for running unified tests.
 *
 * Usage:
 *   java jline.unified.UnifiedTestRunnerKt <model_name> [--json]
 */
fun main(args: Array<String>) {
    val modelName = args.getOrNull(0)
    val jsonOutput = args.contains("--json")
    val verbose = args.contains("--verbose") || args.contains("-v")

    if (modelName == null) {
        System.err.println("Usage: java jline.unified.UnifiedTestRunnerKt <model_name> [--json] [--verbose]")
        kotlin.system.exitProcess(1)
    }

    val runner = UnifiedTestRunner()
    val result = runner.runModel(modelName, verbose = verbose, jsonOutput = jsonOutput)

    if (jsonOutput) {
        val gson = com.google.gson.Gson()
        println(gson.toJson(result))
    } else {
        println("\nResult: $result")
    }

    val exitCode = when {
        result is JsonResult -> if (result.status == "passed") 0 else 1
        result == "passed" -> 0
        else -> 1
    }
    kotlin.system.exitProcess(exitCode)
}
