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
     * The directory the {@code keep} option writes its {@code .mat} workspaces to.
     *
     * <p>{@code $LINE_WORKSPACE_ROOT/line_workspace/mat}, or the system temp dir
     * under the same name, which is where the MATLAB twin puts them too:
     * {@code solver_ctmc_marg.m} and its siblings name the file with
     * {@code lineTempName}.</p>
     *
     * <p>This used to resolve {@code <installRoot>/jar/workspace}, or
     * {@code <installRoot>/python/workspace} when the working directory
     * happened to contain the substring "python". Both are wrong for the same
     * reason: an installation directory is not a place to write to. It may be
     * read-only, it may be shared between users, and on a release archive the
     * only trace of the scheme was a pair of empty {@code workspace} folders
     * shipped so that the write would have somewhere to land. The cwd
     * substring test was also a guess about the caller that a caller under,
     * say, {@code /home/python-dev} answered wrongly.</p>
     *
     * @return the workspace directory path, which need not exist yet
     */
    public static String getWorkspaceDirectory() {
        String envRoot = System.getenv("LINE_WORKSPACE_ROOT");
        File baseRoot;
        if (envRoot != null && !envRoot.trim().isEmpty()) {
            baseRoot = new File(envRoot.trim());
        } else {
            String tmp = System.getProperty("java.io.tmpdir");
            baseRoot = new File(tmp == null || tmp.isEmpty() ? "/tmp" : tmp);
        }
        return new File(new File(baseRoot, "line_workspace"), "mat").getPath();
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