/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

public final class Infer_quick_model {
    private Infer_quick_model() {}

    /**
     * Generate a simple queueing network based on given parameters.
     */
    public static Network infer_quick_model(
            boolean open,
            List<SchedStrategy> stations,
            double[][] classes,
            int[] servers,
            int[] jobs,
            double[][] routing) {
        int numStations = stations.size();
        int numClasses = classes.length;
        int[] actualServers;
        if (servers != null) {
            actualServers = servers;
        } else {
            actualServers = new int[numStations];
            for (int i = 0; i < numStations; i++) actualServers[i] = 1;
        }
        int[] actualJobs;
        if (jobs != null) {
            actualJobs = jobs;
        } else {
            actualJobs = new int[numClasses];
            for (int i = 0; i < numClasses; i++) actualJobs[i] = 1;
        }

        Network model = new Network("quickModel");
        Node[] nodes = new Node[numStations + (open ? 2 : 0)];

        if (open) {
            nodes[0] = new Source(model, "mySource");
            nodes[numStations + 1] = new Sink(model, "mySink");
        }

        for (int i = 0; i < numStations; i++) {
            int idx = open ? i + 1 : i;
            Queue queue = new Queue(model, "QueueStation" + (i + 1), stations.get(i));
            queue.setNumberOfServers(actualServers[i]);
            nodes[idx] = queue;
        }

        RoutingMatrix P = model.initRoutingMatrix();
        JobClass[] jobClasses = new JobClass[numClasses];

        for (int c = 0; c < numClasses; c++) {
            if (!open) {
                Station refNode = (Station) nodes[0];
                jobClasses[c] = new ClosedClass(model, "Class" + (c + 1), actualJobs[c], refNode);
            } else {
                jobClasses[c] = new OpenClass(model, "Class" + (c + 1));
            }

            for (int i = 0; i < numStations; i++) {
                int idx = open ? i + 1 : i;
                Queue queue = (Queue) nodes[idx];
                queue.setService(jobClasses[c], Exp.fitMean(classes[c][i]));
            }
        }

        if (open) {
            List<Node> nodeList = new ArrayList<Node>();
            for (int i = 0; i < nodes.length; i++) {
                if (nodes[i] != null) nodeList.add(nodes[i]);
            }
            for (int c = 0; c < numClasses; c++) {
                P.set(jobClasses[c], Network.serialRouting(nodeList));
            }
        } else {
            for (int c = 0; c < numClasses; c++) {
                if (routing != null) {
                    P.set(jobClasses[c], new Matrix(routing[c]));
                } else {
                    List<Node> nodeList = new ArrayList<Node>();
                    for (int i = 0; i < numStations; i++) {
                        if (nodes[i] != null) nodeList.add(nodes[i]);
                    }
                    P.set(jobClasses[c], Network.serialRouting(nodeList));
                }
            }
        }

        model.link(P);
        return model;
    }

    public static Network infer_quick_model(boolean open, List<SchedStrategy> stations, double[][] classes) {
        return infer_quick_model(open, stations, classes, null, null, null);
    }
}
