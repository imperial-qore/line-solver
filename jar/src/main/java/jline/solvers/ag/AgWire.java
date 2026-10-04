/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import jline.solvers.ag.handlers.RCATModel;
import jline.util.matrix.Matrix;

/**
 * The ag-worker wire protocol, encoder and decoder in one place.
 *
 * <p>ONE CLASS ON PURPOSE. The coordinator and the worker are different
 * processes and, in the MATLAB and Python clients, different languages; a
 * protocol written twice is a protocol that drifts, and the failure it produces
 * is a wrong number rather than a parse error. Everything either side puts on
 * the wire is built here.</p>
 *
 * <p>The transport is newline-delimited JSON over a plain TCP socket: one JSON
 * object per line, request then reply, no framing beyond the newline. It is
 * deliberately NOT the existing LineWebSocket protocol, whose unit is a whole
 * model and whose first line is a CLI argument vector -- the unit here is one
 * agent, and the per-sweep message is a vector of reversed rates.</p>
 *
 * <p>Messages:</p>
 * <ul>
 *   <li><b>assign</b> {@code {"op":"assign","agents":[...]}} carries the STATIC
 *       half of each agent -- its local rate matrix and the passive/active
 *       matrices of the actions it takes part in. Sent once per solve, because
 *       none of it changes across the fixed point. Reply
 *       {@code {"op":"assigned","k":[...]}}.</li>
 *   <li><b>sweep</b> {@code {"op":"sweep","x":[...]}} carries the reversed
 *       rates, one double per action -- the entire coupling between agents.
 *       Reply {@code {"op":"swept","agents":[{"k":..,"pi":[..]}]}}.</li>
 *   <li><b>bye</b> {@code {"op":"bye"}} closes the session.</li>
 * </ul>
 *
 * <p>Matrices travel as {@code [row, col, value]} triplets with 0-based indices,
 * because an action's matrices touch one level transition and are nearly empty;
 * a dense N-by-N payload would be almost all zeros.</p>
 */
public final class AgWire {

    private AgWire() {}

    /** Non-zero entries of M as a JSON array of [row, col, value] triplets. */
    public static JsonArray triplets(Matrix m) {
        JsonArray out = new JsonArray();
        if (m == null) return out;
        int rows = m.getNumRows();
        int cols = m.getNumCols();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double v = m.get(i, j);
                if (v != 0.0) {
                    JsonArray t = new JsonArray();
                    t.add(Integer.valueOf(i));
                    t.add(Integer.valueOf(j));
                    t.add(Double.valueOf(v));
                    out.add(t);
                }
            }
        }
        return out;
    }

    /** Rebuild an n-by-n matrix from the triplets written by {@link #triplets}. */
    public static Matrix fromTriplets(JsonArray a, int n) {
        Matrix m = new Matrix(n, n);
        if (a == null) return m;
        for (int t = 0; t < a.size(); t++) {
            JsonArray e = a.get(t).getAsJsonArray();
            m.set(e.get(0).getAsInt(), e.get(1).getAsInt(), e.get(2).getAsDouble());
        }
        return m;
    }

    /** The static description of agent k, as the worker needs it. */
    public static JsonObject agent(int k, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                   int[] ACT, int[] PSV, int numActions, int[] N,
                                   RCATModel rcat) {
        JsonObject a = new JsonObject();
        a.addProperty("k", Integer.valueOf(k));
        a.addProperty("n", Integer.valueOf(N[k]));
        a.addProperty("mph", Integer.valueOf(rcat.mph[k]));
        a.addProperty("nlev", Integer.valueOf(rcat.nlev[k]));

        JsonArray lvl = new JsonArray();
        int[] level = rcat.level[k];
        for (int i = 0; i < level.length; i++) lvl.add(Integer.valueOf(level[i]));
        a.add("level", lvl);

        a.add("L", triplets(L[k]));

        JsonArray passive = new JsonArray();
        JsonArray active = new JsonArray();
        for (int c = 0; c < numActions; c++) {
            if (PSV[c] == k && Pb[c] != null) {
                JsonObject e = new JsonObject();
                e.addProperty("c", Integer.valueOf(c));
                e.add("M", triplets(Pb[c]));
                passive.add(e);
            } else if (ACT[c] == k && Aa[c] != null) {
                JsonObject e = new JsonObject();
                e.addProperty("c", Integer.valueOf(c));
                e.add("M", triplets(Aa[c]));
                active.add(e);
            }
        }
        a.add("passive", passive);
        a.add("active", active);
        return a;
    }

    /** The reversed-rate vector, the whole per-sweep payload. */
    public static JsonArray rates(Matrix x) {
        JsonArray out = new JsonArray();
        for (int c = 0; c < x.getNumRows(); c++) out.add(Double.valueOf(x.get(c, 0)));
        return out;
    }

    /** A stationary vector as a JSON array. */
    public static JsonArray vector(Matrix pi) {
        JsonArray out = new JsonArray();
        for (int j = 0; j < pi.getNumCols(); j++) out.add(Double.valueOf(pi.get(0, j)));
        return out;
    }

    /** The inverse of {@link #vector}. */
    public static Matrix toVector(JsonArray a) {
        Matrix m = new Matrix(1, a.size());
        for (int j = 0; j < a.size(); j++) m.set(0, j, a.get(j).getAsDouble());
        return m;
    }

    /** The inverse of {@link #rates}, as the column vector the solver uses. */
    public static Matrix toRates(JsonArray a) {
        Matrix m = new Matrix(a.size(), 1);
        for (int c = 0; c < a.size(); c++) m.set(c, 0, a.get(c).getAsDouble());
        return m;
    }

    public static List<Integer> toIntList(JsonArray a) {
        List<Integer> out = new ArrayList<Integer>();
        for (int i = 0; i < a.size(); i++) out.add(Integer.valueOf(a.get(i).getAsInt()));
        return out;
    }
}
