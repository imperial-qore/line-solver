package jline;

import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test that website examples are structurally correct.
 * <p>
 * This test verifies that all Java examples shown on the LINE website solver*.html
 * pages have the correct structure and can be used by users.
 */
public class WebsiteExamplesTest {

    @Test
    public void testWebsiteExamplesExist() throws Exception {
        // Verify the JSON file exists
        String jsonPath = findJsonFile();
        assertTrue(Files.exists(Paths.get(jsonPath)),
                "website_examples.json not found. Run doc/extract_website_examples.py first.");

        // Read and parse the JSON file manually (simple parsing without dependencies)
        List<String> javaCodeBlocks = extractJavaCodeBlocks(jsonPath);

        int totalTests = javaCodeBlocks.size();
        int passedTests = 0;
        List<String> failedTests = new ArrayList<>();

        for (int i = 0; i < javaCodeBlocks.size(); i++) {
            String code = javaCodeBlocks.get(i);
            String testName = String.format("java_example_%d", i + 1);

            try {
                if (verifyJavaCodeStructure(code)) {
                    passedTests++;
                } else {
                    failedTests.add(testName);
                }
            } catch (Exception e) {
                failedTests.add(testName);
            }
        }

        if (!failedTests.isEmpty()) {
            System.out.println("\nFailed tests:");
            for (String test : failedTests) {
                System.out.printf("  - %s%n", test);
            }
        }


        // Assert all tests passed
        if (passedTests != totalTests) {
            fail(String.format("%d out of %d tests failed.", totalTests - passedTests, totalTests));
        }
    }

    /**
     * Extract Java code blocks from the JSON file using simple string parsing
     */
    private List<String> extractJavaCodeBlocks(String jsonPath) throws IOException {
        List<String> codeBlocks = new ArrayList<>();
        StringBuilder currentCode = new StringBuilder();
        boolean inJavaCode = false;
        boolean nextLineIsCode = false;

        try (BufferedReader reader = new BufferedReader(new FileReader(jsonPath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();

                // Detect Java lang marker
                if (trimmed.contains("\"lang\": \"java\"")) {
                    inJavaCode = true;
                    continue;
                }

                // Detect code field start
                if (inJavaCode && trimmed.startsWith("\"code\":")) {
                    nextLineIsCode = true;
                    currentCode = new StringBuilder();
                    // Handle case where code starts on same line
                    int codeStart = line.indexOf("\"code\": \"") + 9;
                    if (codeStart > 9) {
                        String codePart = line.substring(codeStart);
                        currentCode.append(unescapeJson(codePart));
                    }
                    continue;
                }

                // Collect code lines
                if (nextLineIsCode) {
                    if (trimmed.startsWith("\"output\":") || trimmed.startsWith("\"description\":")) {
                        // End of code block
                        String code = currentCode.toString();
                        if (code.contains("import jline.") && code.contains("Network")) {
                            codeBlocks.add(code);
                        }
                        inJavaCode = false;
                        nextLineIsCode = false;
                        currentCode = new StringBuilder();
                    } else {
                        currentCode.append(unescapeJson(line));
                    }
                }
            }
        }

        return codeBlocks;
    }

    /**
     * Simple JSON string unescape
     */
    private String unescapeJson(String json) {
        return json.replace("\\n", "\n")
                .replace("\\\"", "\"")
                .replace("\\\\", "\\")
                .replace("\\t", "\t");
    }

    /**
     * Verify that Java code has the expected structure
     */
    private boolean verifyJavaCodeStructure(String code) {
        // Check for required imports
        if (!code.contains("import jline.")) {
            return false;
        }

        // Check for Network creation
        if (!code.contains("new Network(")) {
            return false;
        }

        // Check for key model components
        boolean hasNodes = code.contains("Queue") || code.contains("Source") ||
                code.contains("Delay") || code.contains("Sink");
        boolean hasClasses = code.contains("OpenClass") || code.contains("ClosedClass");

        // Check for solver usage (most examples should have this)
        boolean hasSolver = code.contains("Solver") || code.contains("MVA") ||
                code.contains("JMT") || code.contains("SSA") ||
                code.contains("CTMC") || code.contains("MAM") ||
                code.contains("AUTO") || code.contains("NC") ||
                code.contains("DES");

        return hasNodes && hasClasses && hasSolver;
    }

    /**
     * Find the JSON file with website examples
     */
    private String findJsonFile() throws IOException {
        // Try relative path from jar/ directory
        String path1 = "../doc/website_examples.json";
        if (Files.exists(Paths.get(path1))) {
            return path1;
        }

        // Try from repository root
        String path2 = "doc/website_examples.json";
        if (Files.exists(Paths.get(path2))) {
            return path2;
        }

        throw new IOException("Cannot find website_examples.json");
    }
}
