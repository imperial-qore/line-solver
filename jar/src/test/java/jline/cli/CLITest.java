package jline.cli;

import jline.GlobalConstants;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;

/**
 * Unit tests for LineCLI command line interface
 */
public class CLITest {

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();
    private final PrintStream originalOut = System.out;
    private final PrintStream originalErr = System.err;

    @TempDir
    Path tempDir;

    @BeforeEach
    public void setUpStreams() {
        System.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(errContent));
    }

    @AfterEach
    public void restoreStreams() {
        System.setOut(originalOut);
        System.setErr(originalErr);
    }

    /**
     * Test that help is displayed when no arguments are provided
     */
    @Test
    public void testNoArguments() throws IOException {
        String result = LineCLI.parseArgs(new String[]{});
        assertNull(result, "parseArgs should return null when displaying help");
        String output = outContent.toString();
        assertTrue(output.contains("LINE Solver - Command Line Interface"));
        assertTrue(output.contains("USAGE:"));
        assertTrue(output.contains("-h, --help"));
    }

    /**
     * Test help option (-h)
     */
    @Test
    public void testHelpShortOption() throws IOException {
        String result = LineCLI.parseArgs(new String[]{"-h"});
        assertNull(result, "parseArgs should return null when displaying help");
        String output = outContent.toString();
        assertTrue(output.contains("LINE Solver - Command Line Interface"));
        assertTrue(output.contains("USAGE:"));
        assertTrue(output.contains("-h, --help"));
    }

    /**
     * Test help option (--help)
     */
    @Test
    public void testHelpLongOption() throws IOException {
        String result = LineCLI.parseArgs(new String[]{"--help"});
        assertNull(result, "parseArgs should return null when displaying help");
        String output = outContent.toString();
        assertTrue(output.contains("LINE Solver - Command Line Interface"));
        assertTrue(output.contains("USAGE:"));
        assertTrue(output.contains("-h, --help"));
    }

    /**
     * Test version option (-V)
     */
    @Test
    public void testVersionShortOption() throws IOException {
        String result = LineCLI.parseArgs(new String[]{"-V"});
        assertNull(result, "parseArgs should return null when displaying version");
        String output = outContent.toString();
        assertTrue(output.contains(GlobalConstants.getVersion()), "Version output should contain version number");
    }

    /**
     * Test version option (--version)
     */
    @Test
    public void testVersionLongOption() throws IOException {
        String result = LineCLI.parseArgs(new String[]{"--version"});
        assertNull(result, "parseArgs should return null when displaying version");
        String output = outContent.toString();
        assertTrue(output.contains(GlobalConstants.getVersion()), "Version output should contain version number");
    }

    /**
     * Test default output format is readable
     */
    @Test
    public void testDefaultOutputFormat() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        // Copy the working JSIM file from test resources instead of using inline content
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test with file input - default output should be readable
        String result = LineCLI.parseArgs(new String[]{"-f", testFile.toString(), "-s", "fluid"});
        assertNotNull(result, "Result should not be null for valid input");
        // Check that output is in readable format (contains table headers)
        assertTrue(result.contains("Station") || result.contains("JobClass") || 
                  result.contains("---"), "Default output should be in readable format");
        assertFalse(result.startsWith("[{"), "Default output should not be JSON");
    }

    /**
     * Test explicit JSON output format
     */
    @Test
    public void testJsonOutputFormat() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        // Copy the working JSIM file from test resources instead of using inline content
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test with explicit JSON output
        String result = LineCLI.parseArgs(new String[]{"-f", testFile.toString(), "-s", "fluid", "-o", "json"});
        assertNotNull(result, "Result should not be null for valid input");
        assertTrue(result.startsWith("{"), "JSON output should start with {");
        assertTrue(result.contains("\"avg\"") || result.contains("\"sys\"") || result.contains("\"all\""),
            "JSON output should contain analysis type key");
    }

    /**
     * Test readable output format
     */
    @Test
    public void testReadableOutputFormat() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        // Copy the working JSIM file from test resources instead of using inline content
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test with explicit readable output
        String result = LineCLI.parseArgs(new String[]{"-f", testFile.toString(), "-s", "fluid", "-o", "readable"});
        assertNotNull(result, "Result should not be null for valid input");
        assertTrue(result.contains("Station") || result.contains("JobClass") || 
                  result.contains("---"), "Readable output should contain table formatting");
        assertFalse(result.startsWith("[{"), "Readable output should not be JSON");
    }

    /**
     * Test input format options
     */
    @Test
    public void testInputFormatOptions() throws IOException {
        // Test that input format options are accepted
        String[] inputFormats = {"jsim", "jsimg", "jsimw", "lqnx", "xml"};
        
        for (String format : inputFormats) {
            outContent.reset();
            errContent.reset();
            
            // Just test that the option is accepted (will fail on file read)
            try {
                LineCLI.parseArgs(new String[]{"-i", format, "-f", "nonexistent.file"});
            } catch (Exception e) {
                // Expected - file doesn't exist
            }
            
            String errorOutput = errContent.toString();
            assertFalse(errorOutput.contains("Unknown input format"), 
                       "Input format " + format + " should be recognized");
        }
    }

    /**
     * Test solver options
     */
    @Test
    public void testSolverOptions() throws IOException {
        String[] solvers = {"mva", "ctmc", "fluid", "jmt", "nc", "ssa", "ln", "lqns"};
        
        for (String solver : solvers) {
            outContent.reset();
            errContent.reset();
            
            // Just test that the option is accepted
            try {
                LineCLI.parseArgs(new String[]{"-s", solver, "-f", "nonexistent.file"});
            } catch (Exception e) {
                // Expected - file doesn't exist
            }
            
            String errorOutput = errContent.toString();
            assertFalse(errorOutput.contains("Unknown solver"), 
                       "Solver " + solver + " should be recognized");
        }
    }

    /**
     * Test analysis type options
     */
    @Test
    public void testAnalysisOptions() throws IOException {
        String[] analysisTypes = {"all", "avg", "sys"};
        
        for (String analysis : analysisTypes) {
            outContent.reset();
            errContent.reset();
            
            // Just test that the option is accepted
            try {
                LineCLI.parseArgs(new String[]{"-a", analysis, "-f", "nonexistent.file"});
            } catch (Exception e) {
                // Expected - file doesn't exist
            }
            
            String errorOutput = errContent.toString();
            assertFalse(errorOutput.contains("Unknown analysis type"), 
                       "Analysis type " + analysis + " should be recognized");
        }
    }

    /**
     * Test verbosity options
     */
    @Test
    public void testVerbosityOptions() throws IOException {
        String[] verbosityLevels = {"normal", "silent"};
        
        for (String verbosity : verbosityLevels) {
            outContent.reset();
            errContent.reset();
            
            // Just test that the option is accepted
            try {
                LineCLI.parseArgs(new String[]{"-v", verbosity, "-f", "nonexistent.file"});
            } catch (Exception e) {
                // Expected - file doesn't exist
            }
            
            String errorOutput = errContent.toString();
            assertFalse(errorOutput.contains("Unknown verbosity"), 
                       "Verbosity level " + verbosity + " should be recognized");
        }
    }

    /**
     * Test seed option
     */
    @Test
    public void testSeedOption() throws IOException {
        outContent.reset();
        errContent.reset();
        
        // Test that seed option is accepted
        try {
            LineCLI.parseArgs(new String[]{"-d", "12345", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }
        
        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter -d requires"), 
                   "Seed option should be recognized");
    }

    /**
     * Test error handling for missing parameter values
     */
    @Test
    public void testMissingParameterValue() throws IOException {
        outContent.reset();
        errContent.reset();
        
        String result = LineCLI.parseArgs(new String[]{"-f"});
        assertNull(result, "Should return null for missing parameter value");
        
        String errorOutput = errContent.toString();
        assertTrue(errorOutput.contains("Parameter -f requires a value"), 
                  "Should report missing value for -f");
    }

    /**
     * Test error handling for unknown options
     */
    @Test
    public void testUnknownOption() throws IOException {
        outContent.reset();
        errContent.reset();
        
        String result = LineCLI.parseArgs(new String[]{"--unknown", "value"});
        assertNull(result, "Should return null for unknown option");
        
        String errorOutput = errContent.toString();
        assertTrue(errorOutput.contains("Unknown parameter: --unknown"), 
                  "Should report unknown parameter");
    }

    /**
     * Test combination of options
     */
    @Test
    public void testCombinedOptions() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        // Copy the working JSIM file from test resources instead of using inline content
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test combination of multiple options
        String result = LineCLI.parseArgs(new String[]{
            "-f", testFile.toString(),
            "-i", "jsim",
            "-o", "json",
            "-s", "fluid",
            "-a", "avg",
            "-v", "silent"
        });
        
        assertNotNull(result, "Result should not be null for valid combined options");
        assertTrue(result.startsWith("{"), "Output should be JSON as specified");
    }

    /**
     * Test long form options
     */
    @Test
    public void testLongFormOptions() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        // Copy the working JSIM file from test resources instead of using inline content
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test long form options
        String result = LineCLI.parseArgs(new String[]{
            "--file", testFile.toString(),
            "--input", "jsim",
            "--output", "json",
            "--solver", "fluid",
            "--analysis", "avg",
            "--verbosity", "silent",
            "--seed", "42"
        });

        assertNotNull(result, "Result should not be null for valid long form options");
        assertTrue(result.startsWith("{"), "Output should be JSON as specified");
    }

    /**
     * Test extended analysis type options including new types
     */
    @Test
    public void testExtendedAnalysisOptions() throws IOException {
        String[] analysisTypes = {"all", "avg", "sys", "stage", "chain", "node", "nodechain",
            "cdf-respt", "cdf-passt", "perct-respt",
            "tran-avg", "tran-cdf-respt", "tran-cdf-passt",
            "prob", "prob-aggr", "prob-marg", "prob-sys", "prob-sys-aggr",
            "sample", "sample-aggr", "sample-sys", "sample-sys-aggr",
            "reward", "reward-steady", "reward-value"};

        for (String analysis : analysisTypes) {
            outContent.reset();
            errContent.reset();

            // Just test that the option is accepted
            try {
                LineCLI.parseArgs(new String[]{"-a", analysis, "-f", "nonexistent.file"});
            } catch (Exception e) {
                // Expected - file doesn't exist
            }

            String errorOutput = errContent.toString();
            assertFalse(errorOutput.contains("Unknown analysis type"),
                       "Analysis type " + analysis + " should be recognized");
        }
    }

    /**
     * Test comma-separated multi-analysis types
     */
    @Test
    public void testMultiAnalysisOptions() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test comma-separated analysis types
        String result = LineCLI.parseArgs(new String[]{
            "-f", testFile.toString(),
            "-s", "fluid",
            "-a", "avg,sys"
        });

        assertNotNull(result, "Result should not be null for multi-analysis");
    }

    /**
     * Test node index parameter (-n)
     */
    @Test
    public void testNodeIndexOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that node option is accepted
        try {
            LineCLI.parseArgs(new String[]{"-n", "1", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter -n requires"),
                   "Node option should be recognized");
    }

    /**
     * Test class index parameter (-c)
     */
    @Test
    public void testClassIndexOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that class option is accepted
        try {
            LineCLI.parseArgs(new String[]{"-c", "0", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter -c requires"),
                   "Class option should be recognized");
    }

    /**
     * Test events parameter
     */
    @Test
    public void testEventsOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that events option is accepted
        try {
            LineCLI.parseArgs(new String[]{"--events", "5000", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter --events requires"),
                   "Events option should be recognized");
    }

    /**
     * Test percentiles parameter
     */
    @Test
    public void testPercentilesOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that percentiles option is accepted
        try {
            LineCLI.parseArgs(new String[]{"--percentiles", "50,90,95,99", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter --percentiles requires"),
                   "Percentiles option should be recognized");
    }

    /**
     * Test reward-name parameter
     */
    @Test
    public void testRewardNameOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that reward-name option is accepted
        try {
            LineCLI.parseArgs(new String[]{"--reward-name", "QLen", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter --reward-name requires"),
                   "Reward-name option should be recognized");
    }

    /**
     * Test MAM solver option
     */
    @Test
    public void testMamSolverOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that MAM solver option is accepted
        try {
            LineCLI.parseArgs(new String[]{"-s", "mam", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Unknown solver"),
                   "MAM solver should be recognized");
    }

    /**
     * Test state parameter
     */
    @Test
    public void testStateOption() throws IOException {
        outContent.reset();
        errContent.reset();

        // Test that state option is accepted
        try {
            LineCLI.parseArgs(new String[]{"--state", "1,0,0", "-f", "nonexistent.file"});
        } catch (Exception e) {
            // Expected - file doesn't exist
        }

        String errorOutput = errContent.toString();
        assertFalse(errorOutput.contains("Parameter --state requires"),
                   "State option should be recognized");
    }

    /**
     * Test stage analysis with fluid solver
     */
    @Test
    public void testStageAnalysis() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test stage analysis type
        String result = LineCLI.parseArgs(new String[]{
            "-f", testFile.toString(),
            "-s", "fluid",
            "-a", "stage"
        });

        // Result may be null if getStageTable returns null for this model
        // Just verify no exception is thrown
    }

    /**
     * Test help text includes new analysis types
     */
    @Test
    public void testHelpIncludesNewAnalysisTypes() throws IOException {
        LineCLI.parseArgs(new String[]{"-h"});
        String output = outContent.toString();

        // Check that help mentions new analysis types
        assertTrue(output.contains("stage"), "Help should mention stage analysis");
        assertTrue(output.contains("cdf-respt"), "Help should mention cdf-respt analysis");
        assertTrue(output.contains("sample"), "Help should mention sample analysis");
        assertTrue(output.contains("reward"), "Help should mention reward analysis");
    }

    /**
     * Test help text includes new parameters
     */
    @Test
    public void testHelpIncludesNewParameters() throws IOException {
        LineCLI.parseArgs(new String[]{"-h"});
        String output = outContent.toString();

        // Check that help mentions new parameters
        assertTrue(output.contains("-n") || output.contains("--node"), "Help should mention node parameter");
        assertTrue(output.contains("--events"), "Help should mention events parameter");
        assertTrue(output.contains("--percentiles"), "Help should mention percentiles parameter");
        assertTrue(output.contains("--reward-name"), "Help should mention reward-name parameter");
    }

    /**
     * Test combined extended options
     */
    @Test
    public void testCombinedExtendedOptions() throws IOException {
        // Create a simple test JSIM file
        Path testFile = tempDir.resolve("test.jsim");
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        // Test combination of new and existing options
        String result = LineCLI.parseArgs(new String[]{
            "-f", testFile.toString(),
            "-s", "ssa",
            "-a", "avg",
            "-n", "1",
            "--events", "100",
            "-v", "silent"
        });

        assertNotNull(result, "Result should not be null for valid combined options");
    }

    /**
     * Test that invalid analysis type is rejected
     */
    @Test
    public void testInvalidAnalysisType() throws IOException {
        outContent.reset();
        errContent.reset();

        Path testFile = tempDir.resolve("test.jsim");
        Path sourceFile = Paths.get("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
        String jsimContent = new String(Files.readAllBytes(sourceFile));
        Files.write(testFile, jsimContent.getBytes());

        String result = LineCLI.parseArgs(new String[]{
            "-f", testFile.toString(),
            "-s", "mva",
            "-a", "invalid-analysis"
        });

        assertNull(result, "Should return null for invalid analysis type");
        String errorOutput = errContent.toString();
        assertTrue(errorOutput.contains("Invalid analysis type") || errorOutput.contains("Unknown analysis type"),
                  "Should report invalid analysis type");
    }
}