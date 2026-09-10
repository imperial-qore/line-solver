/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import jline.solvers.ag.handlers.RCATModel;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_warning;

/**
 * Coordinator-side connection to one ag-worker process.
 *
 * <p>The connection is opened lazily and the agents are assigned once; after
 * that each sweep is one request carrying the reversed rates and one reply
 * carrying the agents' stationary vectors. Any failure marks the client dead
 * rather than propagating: {@link AgExec} then solves those agents locally,
 * which it can do because an agent depends on the rest of the model only through
 * the reversed rates.</p>
 */
public final class AgWorkerClient implements AutoCloseable {

    private final String endpoint;
    private final int timeoutMillis;

    private Socket socket;
    private BufferedReader in;
    private PrintWriter out;
    private boolean live = true;
    private List<Integer> owned = new ArrayList<Integer>();

    public AgWorkerClient(String endpoint, double timeoutSeconds) {
        this.endpoint = endpoint;
        this.timeoutMillis = (int) Math.max(1000.0, timeoutSeconds * 1000.0);
    }

    public String endpoint() {
        return endpoint;
    }

    public boolean isLive() {
        return live;
    }

    /** Mark the worker unusable for the rest of the run and drop its socket. */
    public void kill() {
        live = false;
        close();
    }

    private void connect() throws Exception {
        int colon = endpoint.lastIndexOf(':');
        if (colon <= 0) {
            throw new IllegalArgumentException("malformed endpoint '" + endpoint
                    + "', expected host:port");
        }
        String host = endpoint.substring(0, colon);
        int port = Integer.parseInt(endpoint.substring(colon + 1).trim());
        socket = new Socket();
        socket.connect(new InetSocketAddress(host, port), timeoutMillis);
        socket.setSoTimeout(timeoutMillis);
        out = new PrintWriter(socket.getOutputStream(), true);
        in = new BufferedReader(new InputStreamReader(socket.getInputStream(), "UTF-8"));
    }

    /** Ship the static half of the owned agents. Sent once per solve. */
    public void assign(List<Integer> agents, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                       int[] ACT, int[] PSV, int numActions, int[] N, RCATModel rcat)
            throws Exception {
        if (socket == null) connect();
        this.owned = agents;

        JsonArray arr = new JsonArray();
        for (int t = 0; t < agents.size(); t++) {
            int k = agents.get(t).intValue();
            arr.add(AgWire.agent(k, Aa, Pb, L, ACT, PSV, numActions, N, rcat));
        }
        JsonObject msg = new JsonObject();
        msg.addProperty("op", "assign");
        msg.add("agents", arr);

        JsonObject reply = call(msg);
        if (!"assigned".equals(reply.get("op").getAsString())) {
            throw new IllegalStateException("worker did not acknowledge the assignment");
        }
    }

    /**
     * One sweep: send the reversed rates, receive the owned agents' stationary
     * vectors in the order they were assigned.
     */
    public List<Matrix> sweep(Matrix x, List<Integer> agents) throws Exception {
        JsonObject msg = new JsonObject();
        msg.addProperty("op", "sweep");
        msg.add("x", AgWire.rates(x));

        JsonObject reply = call(msg);
        if (!"swept".equals(reply.get("op").getAsString())) {
            throw new IllegalStateException("unexpected reply '"
                    + reply.get("op").getAsString() + "'");
        }
        JsonArray got = reply.getAsJsonArray("agents");
        if (got.size() != agents.size()) {
            throw new IllegalStateException("worker answered for " + got.size()
                    + " agent(s), " + agents.size() + " were assigned");
        }
        // Index by k rather than by position: the reply order is the worker's,
        // and reading it positionally would silently transpose two agents'
        // stationary vectors if a worker ever reordered them.
        List<Matrix> out = new ArrayList<Matrix>(agents.size());
        for (int t = 0; t < agents.size(); t++) out.add(null);
        for (int t = 0; t < got.size(); t++) {
            JsonObject e = got.get(t).getAsJsonObject();
            int k = e.get("k").getAsInt();
            int slot = agents.indexOf(Integer.valueOf(k));
            if (slot < 0) {
                throw new IllegalStateException("worker answered for agent " + k
                        + ", which it was never assigned");
            }
            out.set(slot, AgWire.toVector(e.getAsJsonArray("pi")));
        }
        for (int t = 0; t < out.size(); t++) {
            if (out.get(t) == null) {
                throw new IllegalStateException("worker did not answer for agent "
                        + agents.get(t));
            }
        }
        return out;
    }

    private JsonObject call(JsonObject msg) throws Exception {
        out.println(msg.toString());
        String line = in.readLine();
        if (line == null) {
            throw new IllegalStateException("worker closed the connection");
        }
        JsonObject reply = new JsonParser().parse(line).getAsJsonObject();
        if (reply.has("error")) {
            throw new IllegalStateException(reply.get("error").getAsString());
        }
        return reply;
    }

    public void close() {
        if (socket != null) {
            try {
                if (out != null) {
                    JsonObject bye = new JsonObject();
                    bye.addProperty("op", "bye");
                    out.println(bye.toString());
                }
                socket.close();
            } catch (Exception e) {
                line_warning("AgWorkerClient", "closing " + endpoint + ": " + e.getMessage());
            }
            socket = null;
            in = null;
            out = null;
        }
    }
}
