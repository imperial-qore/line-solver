/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.cli;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import static jline.io.SysUtilsKt.lineTempName;

/**
 * A client for sending layered network models to a Docker-based LINE server.
 */
public class LineDockerClient {
    
    /**
     * Sends a layered network model to a Docker-based LINE server.
     * 
     * @param modelXmlPath Path to the model XML file
     * @param outputPath Path where the server response will be saved
     * @param portNumber Port number of the LINE server
     */
    public static void sendModel(String modelXmlPath, String outputPath, String portNumber) {
        sendModel(modelXmlPath, outputPath, "127.0.0.1", portNumber);
    }
    
    /**
     * Sends a layered network model to a Docker-based LINE server.
     * 
     * @param modelXmlPath Path to the model XML file
     * @param outputPath Path where the server response will be saved
     * @param ipNumber IP address of the LINE server
     * @param portNumber Port number of the LINE server
     */
    public static void sendModel(String modelXmlPath, String outputPath, String ipNumber, String portNumber) {
        String filePath = null;
        try {
            filePath = lineTempName("layered");
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        
        try {
            File modelXML = new File(modelXmlPath);
            File out = new File(outputPath);
            ProcessBuilder processBuilder = new ProcessBuilder("java", "-jar", 
                filePath + "/lineclient.jar", ipNumber, portNumber, "-i", "lqnx");
            processBuilder.redirectErrorStream(true);
            processBuilder.redirectInput(modelXML);
            processBuilder.redirectOutput(out);
            
            Process process = processBuilder.start();
            process.waitFor();
            
            InputStream inputStream = process.getInputStream();
            BufferedReader input = new BufferedReader(new InputStreamReader(inputStream));
            String ss = null;
            while ((ss = input.readLine()) != null) {
                System.out.println(ss);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}