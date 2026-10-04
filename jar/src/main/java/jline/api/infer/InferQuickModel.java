/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.infer;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;

/**
 * Convenience factory for the single-layer networks the inference estimators fit.
 *
 * <p>The model is named {@code quickModel} and has one Queue {@code QueueStation<i>}
 * per entry of {@code stations}, with that discipline and {@code servers[i]} servers,
 * and one class {@code Class<c>} per row of {@code demands}, served at station i by
 * {@code Exp.fitMean(demands[c][i])}.
 *
 * <p>OPEN: a Source {@code mySource} and a Sink {@code mySink}, every class routed
 * serially Source, QueueStation1, ..., QueueStationM, Sink. No arrival process is
 * set: the estimators supply it. CLOSED: class c has {@code jobs[c]} jobs referenced
 * at QueueStation1 and follows {@code routing[c]} (an M x M station-to-station
 * matrix) when given, and the cycle 1, 2, ..., M, 1 otherwise.
 *
 * <p>Port of MATLAB infer_quick_model.m and the Python-native infer_quick_model.
 */
public final class InferQuickModel {

    private InferQuickModel() {
    }

    /** Serial routing, one server per station and one job per closed class. */
    public static Network infer_quick_model(boolean isOpen, SchedStrategy[] stations, double[][] demands) {
        return infer_quick_model(isOpen, stations, demands, null, null, null);
    }

    /**
     * @param isOpen   true for an open network, false for a closed one
     * @param stations (M) scheduling discipline of each queue
     * @param demands  (K x M) mean service demand of class c at station i
     * @param servers  (M) server count per station, or null for one each
     * @param jobs     (K) population per closed class, or null for one each
     * @param routing  (K x M x M) per-class station routing of a closed model, or null for serial
     * @return the linked network
     */
    public static Network infer_quick_model(boolean isOpen, SchedStrategy[] stations, double[][] demands,
                                            int[] servers, int[] jobs, double[][][] routing) {
        int M = stations.length;
        int K = demands.length;
        if (M == 0) {
            throw new IllegalArgumentException("infer_quick_model: at least one station is required");
        }
        if (K == 0) {
            throw new IllegalArgumentException("infer_quick_model: at least one class is required");
        }
        for (double[] row : demands) {
            if (row.length != M) {
                throw new IllegalArgumentException("infer_quick_model: demands must have one column per station");
            }
        }
        if (servers != null && servers.length != M) {
            throw new IllegalArgumentException("infer_quick_model: servers must have one entry per station");
        }
        if (jobs != null && jobs.length != K) {
            throw new IllegalArgumentException("infer_quick_model: jobs must have one entry per class");
        }
        if (!isOpen && routing != null) {
            if (routing.length != K) {
                throw new IllegalArgumentException("infer_quick_model: routing must have one matrix per class");
            }
            for (double[][] P : routing) {
                if (P.length != M) {
                    throw new IllegalArgumentException("infer_quick_model: each routing matrix must be M x M");
                }
                for (double[] row : P) {
                    if (row.length != M) {
                        throw new IllegalArgumentException("infer_quick_model: each routing matrix must be M x M");
                    }
                }
            }
        }

        Network model = new Network("quickModel");
        Source source = null;
        Sink sink = null;
        if (isOpen) {
            source = new Source(model, "mySource");
            sink = new Sink(model, "mySink");
        }

        Queue[] queue = new Queue[M];
        for (int i = 0; i < M; i++) {
            queue[i] = new Queue(model, "QueueStation" + (i + 1), stations[i]);
            queue[i].setNumberOfServers(servers == null ? 1 : servers[i]);
        }

        JobClass[] cls = new JobClass[K];
        for (int c = 0; c < K; c++) {
            String nm = "Class" + (c + 1);
            cls[c] = isOpen ? new OpenClass(model, nm)
                    : new ClosedClass(model, nm, jobs == null ? 1 : jobs[c], queue[0]);
            for (int i = 0; i < M; i++) {
                queue[i].setService(cls[c], Exp.fitMean(demands[c][i]));
            }
        }

        RoutingMatrix P = model.initRoutingMatrix();
        for (int c = 0; c < K; c++) {
            JobClass r = cls[c];
            if (isOpen) {
                P.set(r, r, source, queue[0], 1.0);
                for (int i = 0; i + 1 < M; i++) {
                    P.set(r, r, queue[i], queue[i + 1], 1.0);
                }
                P.set(r, r, queue[M - 1], sink, 1.0);
            } else if (routing != null) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        if (routing[c][i][j] > 0) {
                            P.set(r, r, queue[i], queue[j], routing[c][i][j]);
                        }
                    }
                }
            } else {
                for (int i = 0; i < M; i++) {
                    P.set(r, r, queue[i], queue[(i + 1) % M], 1.0);
                }
            }
        }
        model.link(P);
        return model;
    }
}
