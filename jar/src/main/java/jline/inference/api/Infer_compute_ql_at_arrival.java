/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Collections;

public final class Infer_compute_ql_at_arrival {
    private Infer_compute_ql_at_arrival() {}

    private static final class EventEntry {
        final double time;
        final int type;
        final int sortedIdx;
        final int classIdx;
        EventEntry(double time, int type, int sortedIdx, int classIdx) {
            this.time = time; this.type = type; this.sortedIdx = sortedIdx; this.classIdx = classIdx;
        }
    }

    /**
     * Compute per-class queue lengths at arrival.
     */
    public static Matrix infer_compute_ql_at_arrival(double[] at, int[] atJobid, double[] rt, int[] rtJobid,
                                                     int[] classVec, int R) {
        int n = at.length;

        HashMap<Integer, Double> rtMap = new HashMap<Integer, Double>(rtJobid.length);
        for (int i = 0; i < rtJobid.length; i++) {
            rtMap.put(rtJobid[i], rt[i]);
        }
        double[] rtMatched = new double[n];
        for (int i = 0; i < n; i++) {
            Double v = rtMap.get(atJobid[i]);
            if (v == null) {
                throw new IllegalArgumentException(
                    "infer_compute_ql_at_arrival: not all arrival job IDs found in response time job IDs.");
            }
            rtMatched[i] = v.doubleValue();
        }

        // Sort arrivals by time
        Integer[] sortIdxBoxed = new Integer[n];
        for (int i = 0; i < n; i++) sortIdxBoxed[i] = i;
        final double[] atRef = at;
        java.util.Arrays.sort(sortIdxBoxed, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(atRef[a], atRef[b]);
            }
        });
        int[] sortIdx = new int[n];
        double[] atSorted = new double[n];
        int[] classSorted = new int[n];
        double[] rtSorted = new double[n];
        for (int i = 0; i < n; i++) {
            sortIdx[i] = sortIdxBoxed[i];
            atSorted[i] = at[sortIdx[i]];
            classSorted[i] = classVec[sortIdx[i]];
            rtSorted[i] = rtMatched[sortIdx[i]];
        }

        double[] exitTimes = new double[n];
        for (int i = 0; i < n; i++) exitTimes[i] = atSorted[i] + rtSorted[i];

        ArrayList<EventEntry> events = new ArrayList<EventEntry>(2 * n);
        for (int i = 0; i < n; i++) {
            events.add(new EventEntry(atSorted[i], 1, i, classSorted[i]));
            events.add(new EventEntry(exitTimes[i], -1, i, classSorted[i]));
        }

        // Sort by time; departures (-1) before arrivals (+1) at same time
        Collections.sort(events, new Comparator<EventEntry>() {
            @Override
            public int compare(EventEntry a, EventEntry b) {
                int c = Double.compare(a.time, b.time);
                if (c != 0) return c;
                return Integer.compare(a.type, b.type);
            }
        });

        int[] state = new int[R];
        Matrix qlSorted = new Matrix(n, R);

        for (EventEntry ev : events) {
            int c = ev.classIdx;
            if (ev.type == 1) {
                state[c]++;
                for (int r = 0; r < R; r++) {
                    qlSorted.set(ev.sortedIdx, r, state[r]);
                }
            } else {
                state[c]--;
            }
        }

        Matrix ql = new Matrix(n, R);
        for (int si = 0; si < n; si++) {
            int origIdx = sortIdx[si];
            for (int r = 0; r < R; r++) {
                ql.set(origIdx, r, qlSorted.get(si, r));
            }
        }

        return ql;
    }
}
