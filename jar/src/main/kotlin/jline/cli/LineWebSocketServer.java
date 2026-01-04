/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import jline.lang.Model;
import org.java_websocket.WebSocket;
import org.java_websocket.handshake.ClientHandshake;
import org.java_websocket.server.WebSocketServer;

import java.io.IOException;
import java.net.InetSocketAddress;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * LineWebSocketServer is a WebSocket server that receives client connections, processes incoming messages,
 * and interacts with the LINE Solver. It handles both text and binary messages and can be started
 * on a specified port or host and port.
 */
public class LineWebSocketServer extends WebSocketServer {

    public String message;

    /**
     * Constructs a new LineWebSocketServer instance on the specified port.
     *
     * @param port the port on which the server will listen.
     */
    public LineWebSocketServer(int port) {
        super(new InetSocketAddress(port));
    }

    /**
     * Constructs a new LineWebSocketServer instance on the specified host and port.
     *
     * @param host the hostname or IP address on which the server will listen.
     * @param port the port on which the server will listen.
     */
    public LineWebSocketServer(String host, int port) {
        super(new InetSocketAddress(host, port));
    }

    /**
     * Constructs a new LineWebSocketServer instance with a specified address.
     *
     * @param address the address on which the server will listen.
     */
    public LineWebSocketServer(InetSocketAddress address) {
        super(address);
    }

    /**
     * The main method, which serves as the entry point for the server application.
     * It initializes and starts the server on a specified port.
     *
     * @param args command-line arguments (currently unused).
     */
    public static void main(String[] args) {
//        int serverPort = 5863;
//
//        LineWebSocketServer server = new LineWebSocketServer(serverPort);
//        server.run();
        try {
            String ret = LineCLI.parseArgs(args);
            if (ret != null) {
                System.out.println(ret);
            }
        } catch (IOException e) {
            line_error(mfilename(new Object[]{}), "IOException: Failed to read file.");
        }
    }

    /**
     * Called when a connection is closed.
     *
     * @param conn   the WebSocket connection.
     * @param code   the status code indicating the reason for closure.
     * @param reason a textual description of the reason for closure.
     * @param remote indicates whether the connection was closed by the remote host.
     */
    @Override
    public void onClose(WebSocket conn, int code, String reason, boolean remote) {
        System.out.printf("%s%n", reason);
    }

    /**
     * Called when an error occurs.
     *
     * @param conn the WebSocket connection.
     * @param ex   the exception that was thrown.
     */
    @Override
    public void onError(WebSocket conn, Exception ex) {
        ex.printStackTrace();
    }

    /**
     * Called when a message is received from a client.
     * This method processes the message, saves the client data to a temporary file,
     * and invokes the LINE Solver with the specified arguments.
     *
     * @param conn    the WebSocket connection.
     * @param message the message received from the client.
     */
    @Override
    public void onMessage(WebSocket conn, String message) {
        Path tmpfile;
        try {
            tmpfile = Files.createTempFile("linetmp", null);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        String[] msg = message.split("\n", 2);
        String args = msg[0];
        message = msg[1];
        try {
            Files.write(tmpfile, message.getBytes(StandardCharsets.UTF_8));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.printf("Client model saved in: %s%n", tmpfile);
        String[] splitargs = args.split(",");
        splitargs[0] = "--file";
        splitargs[1] = tmpfile.toString();
        String ret;
        try {
            ret = LineCLI.parseArgs(splitargs);
        } catch (IOException e) {
            ret = "Failed to read file.";
        }
        conn.send(ret);
        conn.close();
    }

    /**
     * Called when a binary message is received from a client.
     * This method is currently a placeholder for future binary message processing.
     *
     * @param conn    the WebSocket connection.
     * @param message the binary message received from the client.
     */
    @Override
    public void onMessage(WebSocket conn, ByteBuffer message) {
        // Binary messages are not supported by LINE server
        // All communication should be text-based for model definitions and commands
        line_error(mfilename(new Object[]{}), "Warning: Received unsupported binary message from " + conn.getRemoteSocketAddress());
    }

    /**
     * Called when a new client connection is opened.
     *
     * @param conn      the WebSocket connection.
     * @param handshake the client handshake data.
     */
    @Override
    public void onOpen(WebSocket conn, ClientHandshake handshake) {
        System.out.println("New connection to " + conn.getRemoteSocketAddress());
    }

    /**
     * Called when the server is started.
     * This method prints the server's address and port to the console.
     */
    @Override
    public void onStart() {
        System.out.println("--------------------------------------------------------------------");
        System.out.println("LINE Solver - Command Line Interface");
        System.out.println("Copyright (c) 2012-2026, QORE group, Imperial College London");
        System.out.printf("Version %s. All rights reserved.%n", new Model("").getVersion());
        System.out.println("--------------------------------------------------------------------");
        System.out.printf("Running in server mode on %s.%n", this.getAddress());
        System.out.println("Press Q to stop the server at any time.");
    }
}
