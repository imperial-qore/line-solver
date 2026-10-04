/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.cli;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.CountDownLatch;

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
        try {
            String model = new String(Files.readAllBytes(Paths.get(modelXmlPath)), StandardCharsets.UTF_8);
            // first line carries the CLI arguments, the server overrides the leading IP,PORT pair
            String message = ipNumber + "," + portNumber + ",-i,lqnx\n" + model;
            final StringBuilder reply = new StringBuilder();
            final CountDownLatch replied = new CountDownLatch(1);
            LineWebSocketClient client = new LineWebSocketClient(new URI("ws://" + ipNumber + ":" + portNumber), message) {
                @Override
                public void onMessage(String response) {
                    reply.append(response);
                    replied.countDown();
                }

                @Override
                public void onClose(int code, String reason, boolean remote) {
                    replied.countDown();
                }

                @Override
                public void onError(Exception ex) {
                    ex.printStackTrace();
                    replied.countDown();
                }
            };
            client.connectBlocking();
            replied.await();
            client.closeBlocking();
            Files.write(Paths.get(outputPath), reply.toString().getBytes(StandardCharsets.UTF_8));
        } catch (IOException e) {
            e.printStackTrace();
        } catch (URISyntaxException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }
}
