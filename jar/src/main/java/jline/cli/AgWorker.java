/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.cli;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import jline.solvers.ag.AgAgent;
import jline.solvers.ag.AgWire;
import jline.util.matrix.Matrix;

/**
 * Worker process for the AG solver's {@code cluster} execution backend.
 *
 * <pre>java -cp jline.jar jline.cli.AgWorker -p 5870</pre>
 *
 * <p>It holds a set of agents and answers sweeps. THE WORKER HOLDS NO MODEL: an
 * agent's generator is its own local rate matrix plus the passive matrices of
 * its actions scaled by the reversed rates, so everything the worker needs
 * arrives in the assignment and everything that changes per sweep is the vector
 * of reversed rates. That is the whole reason the decomposition distributes
 * well -- there is no shared state to keep coherent, and a worker that dies
 * costs its agents' wall clock and nothing else.</p>
 *
 * <p>One connection is one session: assign once, sweep many times, bye. Several
 * coordinators may connect at once and each gets its own agent set, so a worker
 * can be shared between runs without them seeing each other.</p>
 */
public final class AgWorker {

    private AgWorker() {}

    public static void main(String[] args) throws Exception {
        int port = 5870;
        int backlog = 16;
        boolean quiet = false;

        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if (("-p".equals(a) || "--port".equals(a)) && i + 1 < args.length) {
                port = Integer.parseInt(args[++i]);
            } else if (("-b".equals(a) || "--backlog".equals(a)) && i + 1 < args.length) {
                backlog = Integer.parseInt(args[++i]);
            } else if ("-q".equals(a) || "--quiet".equals(a)) {
                quiet = true;
            } else if ("-h".equals(a) || "--help".equals(a)) {
                System.out.println("usage: java -cp jline.jar jline.cli.AgWorker [-p PORT] "
                        + "[-b BACKLOG] [-q]");
                System.out.println("  -p  listen port (default 5870)");
                System.out.println("  -b  accept backlog (default 16)");
                System.out.println("  -q  do not log sessions");
                return;
            } else {
                System.err.println("AgWorker: unknown argument '" + a + "'");
                System.exit(2);
            }
        }

        final boolean log = !quiet;
        ServerSocket server = new ServerSocket(port, backlog);
        if (log) {
            System.out.println("AgWorker listening on port " + port);
        }
        try {
            while (true) {
                final Socket client = server.accept();
                Thread t = new Thread(new Runnable() {
                    public void run() {
                        serve(client, log);
                    }
                });
                t.setDaemon(true);
                t.start();
            }
        } finally {
            server.close();
        }
    }

    /** One coordinator session: its own agents, its own lifetime. */
    private static void serve(Socket client, boolean log) {
        List<Agent> agents = new ArrayList<Agent>();
        try {
            BufferedReader in = new BufferedReader(
                    new InputStreamReader(client.getInputStream(), "UTF-8"));
            PrintWriter out = new PrintWriter(client.getOutputStream(), true);

            String line;
            while ((line = in.readLine()) != null) {
                JsonObject msg;
                try {
                    msg = new JsonParser().parse(line).getAsJsonObject();
                } catch (Exception e) {
                    out.println(error("malformed request: " + e.getMessage()));
                    continue;
                }
                String op = msg.has("op") ? msg.get("op").getAsString() : "";

                if ("assign".equals(op)) {
                    agents.clear();
                    JsonArray arr = msg.getAsJsonArray("agents");
                    JsonArray ks = new JsonArray();
                    for (int i = 0; i < arr.size(); i++) {
                        Agent a = Agent.read(arr.get(i).getAsJsonObject());
                        agents.add(a);
                        ks.add(Integer.valueOf(a.k));
                    }
                    JsonObject reply = new JsonObject();
                    reply.addProperty("op", "assigned");
                    reply.add("k", ks);
                    out.println(reply.toString());
                    if (log) {
                        System.out.println("assigned " + agents.size() + " agent(s) to "
                                + client.getRemoteSocketAddress());
                    }

                } else if ("sweep".equals(op)) {
                    if (agents.isEmpty()) {
                        out.println(error("sweep before assign"));
                        continue;
                    }
                    Matrix x = AgWire.toRates(msg.getAsJsonArray("x"));
                    JsonArray answers = new JsonArray();
                    for (int i = 0; i < agents.size(); i++) {
                        Agent a = agents.get(i);
                        Matrix Qk = a.generator(x);
                        Matrix pik = AgAgent.stationary(Qk, a.mph, a.nlev, a.level);
                        JsonObject e = new JsonObject();
                        e.addProperty("k", Integer.valueOf(a.k));
                        e.add("pi", AgWire.vector(pik));
                        answers.add(e);
                    }
                    JsonObject reply = new JsonObject();
                    reply.addProperty("op", "swept");
                    reply.add("agents", answers);
                    out.println(reply.toString());

                } else if ("bye".equals(op)) {
                    break;

                } else {
                    out.println(error("unknown op '" + op + "'"));
                }
            }
        } catch (Exception e) {
            if (log) {
                System.out.println("session ended: " + e.getMessage());
            }
        } finally {
            try {
                client.close();
            } catch (Exception e) {
                // The coordinator has already gone; nothing left to report to.
            }
        }
    }

    private static String error(String message) {
        JsonObject o = new JsonObject();
        o.addProperty("error", message);
        return o.toString();
    }

    /**
     * One agent as the worker holds it: the static matrices, and the action
     * indices that scale them.
     */
    private static final class Agent {
        int k;
        int n;
        int mph;
        int nlev;
        int[] level;
        Matrix L;
        int[] passiveC;
        Matrix[] passiveM;
        int[] activeC;
        Matrix[] activeM;

        static Agent read(JsonObject o) {
            Agent a = new Agent();
            a.k = o.get("k").getAsInt();
            a.n = o.get("n").getAsInt();
            a.mph = o.get("mph").getAsInt();
            a.nlev = o.get("nlev").getAsInt();

            JsonArray lvl = o.getAsJsonArray("level");
            a.level = new int[lvl.size()];
            for (int i = 0; i < lvl.size(); i++) a.level[i] = lvl.get(i).getAsInt();

            a.L = AgWire.fromTriplets(o.getAsJsonArray("L"), a.n);

            JsonArray pas = o.getAsJsonArray("passive");
            a.passiveC = new int[pas.size()];
            a.passiveM = new Matrix[pas.size()];
            for (int i = 0; i < pas.size(); i++) {
                JsonObject e = pas.get(i).getAsJsonObject();
                a.passiveC[i] = e.get("c").getAsInt();
                a.passiveM[i] = AgWire.fromTriplets(e.getAsJsonArray("M"), a.n);
            }

            JsonArray act = o.getAsJsonArray("active");
            a.activeC = new int[act.size()];
            a.activeM = new Matrix[act.size()];
            for (int i = 0; i < act.size(); i++) {
                JsonObject e = act.get(i).getAsJsonObject();
                a.activeC[i] = e.get("c").getAsInt();
                a.activeM[i] = AgWire.fromTriplets(e.getAsJsonArray("M"), a.n);
            }
            return a;
        }

        /**
         * The agent's generator at the reversed rates x, assembled in the same
         * order as the coordinator's own {@code agentGenerator}: local rates with
         * their row-sum diagonal, then the passive matrices scaled by x, then the
         * active ones, then the generator fix-up. The ORDER matters -- floating
         * point addition is not associative, and a different order would make a
         * remote agent's answer differ from a local one in the last bits, which
         * is exactly what the backend identity tests forbid.
         */
        Matrix generator(Matrix x) {
            int[] N = new int[] {n};
            Matrix[] Ls = new Matrix[] {L};

            int numActions = 0;
            for (int i = 0; i < passiveC.length; i++) {
                numActions = Math.max(numActions, passiveC[i] + 1);
            }
            for (int i = 0; i < activeC.length; i++) {
                numActions = Math.max(numActions, activeC[i] + 1);
            }
            numActions = Math.max(numActions, x.getNumRows());

            Matrix[] Pb = new Matrix[numActions];
            Matrix[] Aa = new Matrix[numActions];
            int[] PSV = new int[numActions];
            int[] ACT = new int[numActions];
            for (int c = 0; c < numActions; c++) {
                PSV[c] = -1;
                ACT[c] = -1;
            }
            for (int i = 0; i < passiveC.length; i++) {
                Pb[passiveC[i]] = passiveM[i];
                PSV[passiveC[i]] = 0;
            }
            for (int i = 0; i < activeC.length; i++) {
                Aa[activeC[i]] = activeM[i];
                ACT[activeC[i]] = 0;
            }
            return AgAgent.generator(0, x, Aa, Pb, Ls, ACT, PSV, numActions, N);
        }
    }
}
