/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.unified

import com.google.gson.Gson
import com.google.gson.JsonObject
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.DisplayName
import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.TestFactory
import org.junit.jupiter.api.TestInstance
import org.junit.jupiter.api.BeforeAll
import java.io.File

/**
 * JUnit5 test suite for unified cross-language tests.
 *
 * This class provides JUnit5-compatible test methods that run the unified
 * tests against all available JSON definitions.
 *
 * Usage:
 *   mvn test -Dtest=UnifiedTestSuite
 */
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
@DisplayName("Unified Cross-Language Tests")
class UnifiedTestSuite {

    private lateinit var runner: UnifiedTestRunner
    private lateinit var definitionsDir: File
    private val gson = Gson()

    @BeforeAll
    fun setup() {
        runner = UnifiedTestRunner()
        definitionsDir = findDefinitionsDir()
    }

    private fun findDefinitionsDir(): File {
        var current = File(System.getProperty("user.dir"))
        while (current.parentFile != null) {
            val testDir = File(current, "test/unified/definitions")
            if (testDir.exists()) return testDir
            current = current.parentFile
        }
        // Fallback
        return File("test/unified/definitions")
    }

    /**
     * Dynamic test factory that generates a test for each model definition.
     *
     * This uses JUnit5's dynamic test feature to create one test per model,
     * allowing individual model tests to be run independently.
     */
    @TestFactory
    @DisplayName("Model Tests")
    fun testAllModels(): List<DynamicTest> {
        if (!definitionsDir.exists()) {
            return listOf(
                DynamicTest.dynamicTest("No definitions found") {
                    fail<Nothing>("Definitions directory not found: ${definitionsDir.absolutePath}")
                }
            )
        }

        val jsonFiles = definitionsDir.listFiles { f -> f.extension == "json" }
            ?.sortedBy { it.name }
            ?: emptyList()

        if (jsonFiles.isEmpty()) {
            return listOf(
                DynamicTest.dynamicTest("No definitions found") {
                    fail<Nothing>("No JSON definition files found in ${definitionsDir.absolutePath}")
                }
            )
        }

        return jsonFiles.map { jsonFile ->
            val modelName = jsonFile.nameWithoutExtension
            DynamicTest.dynamicTest("Test model: $modelName") {
                testModel(modelName)
            }
        }
    }

    private fun testModel(modelName: String) {
        val result = runner.runModel(modelName, verbose = true)

        when (result) {
            "skipped" -> {
                // Check reason for skip
                val jsonFile = File(definitionsDir, "$modelName.json")
                if (!jsonFile.exists()) {
                    org.junit.jupiter.api.Assumptions.assumeTrue(false,
                        "Definition file not found for $modelName")
                }

                if (!ModelRegistry.hasModel(modelName)) {
                    org.junit.jupiter.api.Assumptions.assumeTrue(false,
                        "Model $modelName not in registry")
                }

                try {
                    ModelRegistry.getModel(modelName)
                } catch (e: NotImplementedError) {
                    org.junit.jupiter.api.Assumptions.assumeTrue(false,
                        "Model $modelName not yet implemented")
                }
            }
            "failed" -> {
                fail("Model $modelName failed. Errors: ${runner.results.errors}")
            }
            "passed" -> {
                // Test passed
            }
        }
    }

    /**
     * Test that the ModelRegistry has models registered.
     */
    @org.junit.jupiter.api.Test
    @DisplayName("Registry has models registered")
    fun testRegistryHasModels() {
        val models = ModelRegistry.getAvailableModels()
        assertTrue(models.isNotEmpty(), "Registry should have models registered")
    }

    /**
     * Test that basic models can be built from the registry.
     */
    @TestFactory
    @DisplayName("Registry can build models")
    fun testRegistryCanBuildModels(): List<DynamicTest> {
        val basicModels = listOf(
            "oqn_basic", "cqn_repairmen", "mqn_basic",
            "cqn_bcmp_theorem", "oqn_fourqueues"
        )

        return basicModels.map { modelName ->
            DynamicTest.dynamicTest("Build model: $modelName") {
                if (!ModelRegistry.hasModel(modelName)) {
                    org.junit.jupiter.api.Assumptions.assumeTrue(false,
                        "Model $modelName not registered")
                    return@dynamicTest
                }

                try {
                    val model = ModelRegistry.getModel(modelName)
                    assertNotNull(model, "Model $modelName should not be null")
                } catch (e: NotImplementedError) {
                    org.junit.jupiter.api.Assumptions.assumeTrue(false,
                        "Model $modelName not yet implemented")
                }
            }
        }
    }

    /**
     * Test that all definition files are valid JSON.
     */
    @TestFactory
    @DisplayName("Definitions are valid JSON")
    fun testDefinitionsValidJson(): List<DynamicTest> {
        if (!definitionsDir.exists()) {
            return emptyList()
        }

        val jsonFiles = definitionsDir.listFiles { f -> f.extension == "json" }
            ?.sortedBy { it.name }
            ?: emptyList()

        return jsonFiles.map { jsonFile ->
            DynamicTest.dynamicTest("Valid JSON: ${jsonFile.name}") {
                try {
                    val content = jsonFile.readText()
                    val data = gson.fromJson(content, JsonObject::class.java)

                    assertTrue(data.has("version"),
                        "${jsonFile.name} missing 'version'")
                    assertTrue(data.has("modelName"),
                        "${jsonFile.name} missing 'modelName'")
                    assertTrue(data.has("expectedResults"),
                        "${jsonFile.name} missing 'expectedResults'")
                } catch (e: Exception) {
                    fail("Invalid JSON in ${jsonFile.name}: ${e.message}")
                }
            }
        }
    }

    /**
     * Test that definitions conform to schema structure.
     */
    @TestFactory
    @DisplayName("Definitions match schema")
    fun testDefinitionsMatchSchema(): List<DynamicTest> {
        if (!definitionsDir.exists()) {
            return emptyList()
        }

        val jsonFiles = definitionsDir.listFiles { f -> f.extension == "json" }
            ?.sortedBy { it.name }
            ?: emptyList()

        val requiredFields = listOf("version", "modelName", "modelType", "solvers", "expectedResults")
        val validTypes = listOf("open", "closed", "mixed", "layered")

        return jsonFiles.map { jsonFile ->
            DynamicTest.dynamicTest("Schema check: ${jsonFile.name}") {
                val data = gson.fromJson(jsonFile.readText(), JsonObject::class.java)

                for (field in requiredFields) {
                    assertTrue(data.has(field),
                        "${jsonFile.name} missing required field '$field'")
                }

                assertEquals("2.0", data.get("version").asString,
                    "${jsonFile.name} should have version 2.0")

                val modelType = data.get("modelType").asString
                assertTrue(modelType in validTypes,
                    "${jsonFile.name} has invalid modelType: $modelType")
            }
        }
    }
}
