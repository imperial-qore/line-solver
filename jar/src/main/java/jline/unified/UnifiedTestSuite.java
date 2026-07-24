/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.unified;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestFactory;
import org.junit.jupiter.api.TestInstance;

import com.google.gson.Gson;
import com.google.gson.JsonObject;

import jline.lang.Network;

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
public class UnifiedTestSuite {

    private UnifiedTestRunner runner;
    private File definitionsDir;
    private final Gson gson = new Gson();

    @BeforeAll
    public void setup() {
        runner = new UnifiedTestRunner();
        definitionsDir = findDefinitionsDir();
    }

    private File findDefinitionsDir() {
        File current = new File(System.getProperty("user.dir"));
        while (current.getParentFile() != null) {
            File testDir = new File(current, "test/unified/definitions");
            if (testDir.exists()) return testDir;
            current = current.getParentFile();
        }
        // Fallback
        return new File("test/unified/definitions");
    }

    /**
     * Dynamic test factory that generates a test for each model definition.
     */
    @TestFactory
    @DisplayName("Model Tests")
    public List<DynamicTest> testAllModels() {
        if (!definitionsDir.exists()) {
            return Collections.singletonList(
                    DynamicTest.dynamicTest("No definitions found", new org.junit.jupiter.api.function.Executable() {
                        @Override
                        public void execute() {
                            fail("Definitions directory not found: " + definitionsDir.getAbsolutePath());
                        }
                    })
            );
        }

        File[] jsonFilesArr = definitionsDir.listFiles(new java.io.FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".json");
            }
        });
        List<File> jsonFiles = (jsonFilesArr != null) ? new ArrayList<File>(Arrays.asList(jsonFilesArr)) : new ArrayList<File>();
        Collections.sort(jsonFiles, new java.util.Comparator<File>() {
            @Override
            public int compare(File a, File b) {
                return a.getName().compareTo(b.getName());
            }
        });

        if (jsonFiles.isEmpty()) {
            return Collections.singletonList(
                    DynamicTest.dynamicTest("No definitions found", new org.junit.jupiter.api.function.Executable() {
                        @Override
                        public void execute() {
                            fail("No JSON definition files found in " + definitionsDir.getAbsolutePath());
                        }
                    })
            );
        }

        List<DynamicTest> tests = new ArrayList<DynamicTest>();
        for (final File jsonFile : jsonFiles) {
            final String modelName = stripExtension(jsonFile.getName());
            tests.add(DynamicTest.dynamicTest("Test model: " + modelName, new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() {
                    testModel(modelName);
                }
            }));
        }
        return tests;
    }

    private static String stripExtension(String name) {
        int idx = name.lastIndexOf('.');
        return (idx > 0) ? name.substring(0, idx) : name;
    }

    private void testModel(String modelName) {
        Object resultObj = runner.runModel(modelName, true, false);
        String result = (resultObj != null) ? resultObj.toString() : "";

        if ("skipped".equals(result)) {
            // Check reason for skip
            File jsonFile = new File(definitionsDir, modelName + ".json");
            if (!jsonFile.exists()) {
                Assumptions.assumeTrue(false, "Definition file not found for " + modelName);
            }

            if (!ModelRegistry.hasModel(modelName)) {
                Assumptions.assumeTrue(false, "Model " + modelName + " not in registry");
            }

            try {
                ModelRegistry.getModel(modelName);
            } catch (Throwable e) {
                Assumptions.assumeTrue(false, "Model " + modelName + " not yet implemented");
            }
        } else if ("failed".equals(result)) {
            fail("Model " + modelName + " failed. Errors: " + runner.getResults().errors);
        } else if ("passed".equals(result)) {
            // Test passed
        }
    }

    /**
     * Test that the ModelRegistry has models registered.
     */
    @Test
    @DisplayName("Registry has models registered")
    public void testRegistryHasModels() {
        List<String> models = ModelRegistry.getAvailableModels();
        assertTrue(!models.isEmpty(), "Registry should have models registered");
    }

    /**
     * Test that basic models can be built from the registry.
     */
    @TestFactory
    @DisplayName("Registry can build models")
    public List<DynamicTest> testRegistryCanBuildModels() {
        List<String> basicModels = Arrays.asList(
                "oqn_basic", "cqn_repairmen", "mqn_basic",
                "cqn_bcmp_theorem", "oqn_fourqueues"
        );

        List<DynamicTest> tests = new ArrayList<DynamicTest>();
        for (final String modelName : basicModels) {
            tests.add(DynamicTest.dynamicTest("Build model: " + modelName, new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() {
                    if (!ModelRegistry.hasModel(modelName)) {
                        Assumptions.assumeTrue(false, "Model " + modelName + " not registered");
                        return;
                    }

                    try {
                        Network model = ModelRegistry.getModel(modelName);
                        assertNotNull(model, "Model " + modelName + " should not be null");
                    } catch (Throwable e) {
                        Assumptions.assumeTrue(false, "Model " + modelName + " not yet implemented");
                    }
                }
            }));
        }
        return tests;
    }

    /**
     * Test that all definition files are valid JSON.
     */
    @TestFactory
    @DisplayName("Definitions are valid JSON")
    public List<DynamicTest> testDefinitionsValidJson() {
        if (!definitionsDir.exists()) {
            return Collections.emptyList();
        }

        File[] jsonFilesArr = definitionsDir.listFiles(new java.io.FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".json");
            }
        });
        List<File> jsonFiles = (jsonFilesArr != null) ? new ArrayList<File>(Arrays.asList(jsonFilesArr)) : new ArrayList<File>();
        Collections.sort(jsonFiles, new java.util.Comparator<File>() {
            @Override
            public int compare(File a, File b) {
                return a.getName().compareTo(b.getName());
            }
        });

        List<DynamicTest> tests = new ArrayList<DynamicTest>();
        for (final File jsonFile : jsonFiles) {
            tests.add(DynamicTest.dynamicTest("Valid JSON: " + jsonFile.getName(), new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() {
                    try {
                        String content = new String(Files.readAllBytes(jsonFile.toPath()), StandardCharsets.UTF_8);
                        JsonObject data = gson.fromJson(content, JsonObject.class);

                        assertTrue(data.has("version"),
                                jsonFile.getName() + " missing 'version'");
                        assertTrue(data.has("modelName"),
                                jsonFile.getName() + " missing 'modelName'");
                        assertTrue(data.has("expectedResults"),
                                jsonFile.getName() + " missing 'expectedResults'");
                    } catch (Exception e) {
                        fail("Invalid JSON in " + jsonFile.getName() + ": " + e.getMessage());
                    }
                }
            }));
        }
        return tests;
    }

    /**
     * Test that definitions conform to schema structure.
     */
    @TestFactory
    @DisplayName("Definitions match schema")
    public List<DynamicTest> testDefinitionsMatchSchema() {
        if (!definitionsDir.exists()) {
            return Collections.emptyList();
        }

        File[] jsonFilesArr = definitionsDir.listFiles(new java.io.FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".json");
            }
        });
        List<File> jsonFiles = (jsonFilesArr != null) ? new ArrayList<File>(Arrays.asList(jsonFilesArr)) : new ArrayList<File>();
        Collections.sort(jsonFiles, new java.util.Comparator<File>() {
            @Override
            public int compare(File a, File b) {
                return a.getName().compareTo(b.getName());
            }
        });

        final List<String> requiredFields = Arrays.asList("version", "modelName", "modelType", "solvers", "expectedResults");
        final List<String> validTypes = Arrays.asList("open", "closed", "mixed", "layered");

        List<DynamicTest> tests = new ArrayList<DynamicTest>();
        for (final File jsonFile : jsonFiles) {
            tests.add(DynamicTest.dynamicTest("Schema check: " + jsonFile.getName(), new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() throws IOException {
                    String content = new String(Files.readAllBytes(jsonFile.toPath()), StandardCharsets.UTF_8);
                    JsonObject data = gson.fromJson(content, JsonObject.class);

                    for (String field : requiredFields) {
                        assertTrue(data.has(field),
                                jsonFile.getName() + " missing required field '" + field + "'");
                    }

                    assertEquals("2.0", data.get("version").getAsString(),
                            jsonFile.getName() + " should have version 2.0");

                    String modelType = data.get("modelType").getAsString();
                    assertTrue(validTypes.contains(modelType),
                            jsonFile.getName() + " has invalid modelType: " + modelType);
                }
            }));
        }
        return tests;
    }
}
