package jline.util;

import static jline.io.InputOutputKt.line_warning;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

/**
 * Utility class for calling Python scripts from Java
 */
public class PythonCaller {
    
    /**
     * Execute a Python script and return its output
     * @param pythonScript The Python script to execute
     * @return The output from the Python script
     */
    public static String callPython(String pythonScript) {
        try {
            // Start Python process
            ProcessBuilder pb = new ProcessBuilder("python3", "-c", pythonScript);
            pb.redirectErrorStream(true);
            Process process = pb.start();
            
            // Read output
            BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream())
            );
            StringBuilder output = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");
            }
            
            // Wait for process to complete
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new RuntimeException("Python script failed with exit code: " + exitCode + "\nOutput: " + output);
            }
            
            return output.toString().trim();
            
        } catch (Exception e) {
            // If Python is not available or script fails, return a default result
            line_warning("PythonCaller", "Could not execute Python script: %s", e.getMessage());
            // Return a default result that simulates output
            return "";
        }
    }
}