/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import org.apache.commons.cli.*;
import org.java_websocket.client.WebSocketClient;
import org.java_websocket.handshake.ServerHandshake;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.Scanner;

/**
 * LineWebSocketClient is a WebSocket client used to communicate with the LINE Solver server.
 * It extends the WebSocketClient class and provides methods for connecting to the server
 * and sending messages.
 */
public class LineWebSocketClient extends WebSocketClient {

    /**
     * The message to be sent to the server.
     */
    public String message;

    /**
     * Constructs a new LineWebSocketClient instance.
     *
     * @param serverURI the URI of the server to connect to.
     * @param msg       the message to be sent to the server.
     */
    public LineWebSocketClient(URI serverURI, String msg) {
        super(serverURI);
        this.message = msg;
    }

    /**
     * The main method, which serves as the entry point for the program.
     * It processes command-line arguments and establishes a connection with the server.
     *
     * @param args command-line arguments specifying the server IP, port, and other options.
     * @throws URISyntaxException if the URI syntax is incorrect.
     * @throws ParseException     if there is an error in parsing the command-line arguments.
     */
    public static void main(String[] args) throws URISyntaxException, ParseException {
        Options options = new Options();
        options.addOption("h", false, "Help");
        options.addOption("m", false, "Max number of requests");
        options.addOption("s", false, "Solver");
        options.addOption("a", false, "Analysis");
        options.addOption("f", false, "File");
        options.addOption("v", false, "Verbosity");
        options.addOption("i", false, "Input");
        options.addOption("o", false, "Output");
        options.addOption("d", false, "Seed");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption("h")) {
            System.out.println("LINE Solver - Client Interface");
            System.out.println("Copyright (c) 2012-2026, QORE group, Imperial College London");
            System.out.println("-----------------------------------------------------------------------");
            System.out.println("Usage: cat MODEL | java -jar linecli.jar SERVER_IP SERVER_PORT\n");
            System.out.println("Example: cat mymodel.jsimg | java -jar linecli.jar 192.168.0.1 5863 ");
        } else {

            String message = String.join(",", args);
            Scanner sc = new Scanner(System.in);
            for (int i = 1; sc.hasNext(); i++) {
                message = message + "\n" + sc.nextLine();
            }
            LineWebSocketClient c = new LineWebSocketClient(new URI("ws://" + args[0] + ":" + args[1]), message);
            c.connect();
        }
    }

    @Override
    public void onClose(int code, String reason, boolean remote) {
        // Handle the WebSocket connection closing
    }

    @Override
    public void onError(Exception ex) {
        // Handle errors during communication
        ex.printStackTrace();
    }

    @Override
    public void onMessage(String message) {
        // Handle incoming messages from the server
        System.out.println(message);
    }

    @Override
    public void onOpen(ServerHandshake handshakedata) {
        // Handle the WebSocket connection opening
        send(this.message);
    }
}
