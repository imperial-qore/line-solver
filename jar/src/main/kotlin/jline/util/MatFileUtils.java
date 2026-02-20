/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.util.matrix.Matrix;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.MatlabType;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility class for saving matrices and workspaces to MATLAB .mat files
 * using the MFL (MATLAB File Library) for Java.
 */
public class MatFileUtils {
    
    /**
     * Saves a single matrix to a .mat file
     * 
     * @param matrix The matrix to save
     * @param variableName The name of the variable in the .mat file
     * @param filename The output filename
     * @throws IOException if file writing fails
     */
    public static void saveMatrix(Matrix matrix, String variableName, String filename) throws IOException {
        MatFile matFile = Mat5.newMatFile();
        
        // Convert LINE Matrix to MFL Matrix
        us.hebi.matlab.mat.types.Matrix mflMatrix = Mat5.newMatrix(matrix.getNumRows(), matrix.getNumCols(), MatlabType.Double);
        for (int i = 0; i < matrix.getNumRows(); i++) {
            for (int j = 0; j < matrix.getNumCols(); j++) {
                mflMatrix.setDouble(i, j, matrix.get(i, j));
            }
        }
        
        // Add matrix to the .mat file
        matFile.addArray(variableName, mflMatrix);
        
        // Write to file
        Mat5.writeToFile(matFile, filename);
    }
    
    /**
     * Saves multiple matrices to a .mat file as a workspace
     * 
     * @param matrices Map of variable names to matrices
     * @param filename The output filename
     * @throws IOException if file writing fails
     */
    public static void saveWorkspace(Map<String, Matrix> matrices, String filename) throws IOException {
        MatFile matFile = Mat5.newMatFile();
        
        for (Map.Entry<String, Matrix> entry : matrices.entrySet()) {
            String varName = entry.getKey();
            Matrix matrix = entry.getValue();
            
            if (matrix != null && !matrix.isEmpty()) {
                // Convert LINE Matrix to MFL Matrix
                us.hebi.matlab.mat.types.Matrix mflMatrix = Mat5.newMatrix(matrix.getNumRows(), matrix.getNumCols(), MatlabType.Double);
                for (int i = 0; i < matrix.getNumRows(); i++) {
                    for (int j = 0; j < matrix.getNumCols(); j++) {
                        mflMatrix.setDouble(i, j, matrix.get(i, j));
                    }
                }
                
                // Add matrix to the .mat file
                matFile.addArray(varName, mflMatrix);
            }
        }
        
        // Write to file
        Mat5.writeToFile(matFile, filename);
    }
    
    /**
     * Saves CTMC solver workspace to a .mat file
     * 
     * @param stateSpace The state space matrix
     * @param infGen The infinitesimal generator matrix
     * @param pi The steady-state probability vector
     * @param filename The output filename
     * @throws IOException if file writing fails
     */
    public static void saveCTMCWorkspace(Matrix stateSpace, Matrix infGen, Matrix pi, String filename) throws IOException {
        Map<String, Matrix> workspace = new HashMap<>();
        
        if (stateSpace != null) workspace.put("StateSpace", stateSpace);
        if (infGen != null) workspace.put("InfGen", infGen);
        if (pi != null) workspace.put("pi", pi);
        
        saveWorkspace(workspace, filename);
    }
    
    /**
     * Generates a timestamp-based filename for .mat files in the appropriate workspace directory
     * 
     * @param prefix The prefix for the filename
     * @return A filename with timestamp in the appropriate workspace directory
     */
    public static String genFilename(String prefix) {
        long timestamp = System.currentTimeMillis();
        String filename = prefix + "_" + timestamp + ".mat";
        String workspaceDir = getWorkspaceDirectory();
        return workspaceDir + "/" + filename;
    }
    
    /**
     * Determines the appropriate workspace directory based on the runtime context
     * 
     * @return The workspace directory path relative to the JAR location
     */
    public static String getWorkspaceDirectory() {
        // Get the location of the JAR file
        String jarPath = MatFileUtils.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        File jarFile = new File(jarPath);
        File jarDir = jarFile.getParentFile();
        
        // Check if we're running in a Python context by looking for python-specific system properties
        // or by checking the classpath for Python-related components
        String pythonHome = System.getProperty("python.home");
        String pythonPath = System.getProperty("python.path");
        
        // Also check if we're in a Python environment by looking at the current working directory
        String currentDir = System.getProperty("user.dir");
        
        File workspaceDir;
        if (pythonHome != null || pythonPath != null || 
            (currentDir != null && currentDir.contains("python"))) {
            // Return ../python/workspace relative to JAR location
            workspaceDir = new File(jarDir.getParentFile(), "python/workspace");
        } else {
            // Return ../jar/workspace relative to JAR location
            workspaceDir = new File(jarDir.getParentFile(), "jar/workspace");
        }
        
        return workspaceDir.getPath();
    }
    
    /**
     * Ensures the directory exists for the given filename
     * 
     * @param filename The filename to create directory for
     * @throws IOException if directory creation fails
     */
    public static void ensureDirectoryExists(String filename) throws IOException {
        File file = new File(filename);
        File parentDir = file.getParentFile();
        if (parentDir != null && !parentDir.exists()) {
            Files.createDirectories(parentDir.toPath());
        }
    }
    
    /**
     * Ensures the appropriate workspace directory exists
     * 
     * @throws IOException if directory creation fails
     */
    public static void ensureWorkspaceDirectoryExists() throws IOException {
        String workspaceDir = getWorkspaceDirectory();
        File dir = new File(workspaceDir);
        if (!dir.exists()) {
            Files.createDirectories(dir.toPath());
        }
    }
}