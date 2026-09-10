/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * Topology and coefficients of the product-form state-dependent routing of
 * Krzesinski (1987), "Multiclass Queueing Networks with State-Dependent
 * Routing", Performance Evaluation 7(2):125-143, the multiclass generalization
 * of Towsley (1980), J. ACM 27(2):323-337.
 *
 * <p>A network is split into a subnetwork Q(V,V) under SDR and its complement
 * M-V. Q(V,V) has one entry center e and one departure center d, both outside
 * it, and is partitioned into disjoint branches arranged in a hierarchy of
 * nested subnetworks V_1 &gt; V_2 &gt; ... &gt; V_T. Each branch has one entry
 * center, one departure center, and may hold several centers between them.</p>
 *
 * <p>Branch index 1 denotes the complement M-V and is unused, following the
 * paper's own indexing so that the coefficients d_tb transcribe straight from
 * the text. The SDR branches are numbered 2..B, that is indices 1..B-1 of the
 * zero-based arrays here.</p>
 *
 * <p>The same structure is carried twice: once in node indices, read by the
 * per-state routing function, and once in station indices, read by the
 * product-form solver. {@link #toStationIndices} performs the conversion.</p>
 *
 * @see jline.api.pfqn.Pfqn_sdr
 */
public class StateDepRouting implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Entry center e of Q(V,V). */
    public int entry;
    /** Departure center d of Q(V,V); may equal {@link #entry}. */
    public int departure;
    /** branch[b] holds the centers of branch b, b &gt;= 1; branch[0] is unused. */
    public int[][] branch;
    /** entryOf[b] is the entry center e(b) of branch b. */
    public int[] entryOf;
    /** departureOf[b] is the departure center d(b) of branch b. */
    public int[] departureOf;
    /** level[b] is the unique t with B_b in V_t - V_{t+1}; level[0] is unused. */
    public int[] level;
    /** Coefficients C_t of eq. (11), length T. */
    public double[] C;
    /** Coefficients d_tb of eq. (11), T rows by B columns. */
    public double[][] d;

    /** Node objects as declared, kept so the structure can be re-resolved after a refresh. */
    public transient Node entryNode;
    /** Departure node as declared. */
    public transient Node departureNode;
    /** Branch nodes as declared; branchNodes.get(0) is null. */
    public transient List<List<Node>> branchNodes;

    public StateDepRouting() {
    }

    /** Number of branch indices, including the unused complement index 0. */
    public int getNumberOfBranches() {
        return this.branch == null ? 0 : this.branch.length;
    }

    /** Number of levels of subnetwork nesting. */
    public int getNumberOfLevels() {
        return this.C == null ? 0 : this.C.length;
    }

    /**
     * Returns a copy of this structure with every center index mapped through
     * the supplied node-to-station index vector.
     *
     * @param nodeToStation node-to-station index vector, negative where the node is not a station
     * @param nodeNames node names, used only to report a center that is not a station
     * @return the station-indexed twin of this structure
     */
    public StateDepRouting toStationIndices(Matrix nodeToStation, List<String> nodeNames) {
        StateDepRouting out = new StateDepRouting();
        out.entry = mapOne(this.entry, nodeToStation, nodeNames);
        out.departure = mapOne(this.departure, nodeToStation, nodeNames);
        int B = this.branch.length;
        out.branch = new int[B][];
        out.entryOf = new int[B];
        out.departureOf = new int[B];
        for (int b = 1; b < B; b++) {
            int[] src = this.branch[b];
            int[] dst = new int[src.length];
            for (int k = 0; k < src.length; k++) {
                dst[k] = mapOne(src[k], nodeToStation, nodeNames);
            }
            out.branch[b] = dst;
            out.entryOf[b] = mapOne(this.entryOf[b], nodeToStation, nodeNames);
            out.departureOf[b] = mapOne(this.departureOf[b], nodeToStation, nodeNames);
        }
        out.level = this.level.clone();
        out.C = this.C.clone();
        out.d = new double[this.d.length][];
        for (int t = 0; t < this.d.length; t++) {
            out.d[t] = this.d[t].clone();
        }
        return out;
    }

    private static int mapOne(int ind, Matrix nodeToStation, List<String> nodeNames) {
        int ist = (int) nodeToStation.get(0, ind);
        if (ist < 0) {
            String name = (nodeNames != null && ind < nodeNames.size()) ? nodeNames.get(ind) : Integer.toString(ind);
            throw new RuntimeException("Node " + name + " takes part in state-dependent routing but is not a station: "
                    + "the product form is over queue lengths, and a stateless node holds none.");
        }
        return ist;
    }

    /**
     * Structural equality, ignoring the declared node objects. A network admits
     * one subnetwork Q(V,V): the routing probabilities of Krzesinski (1987) are
     * chain independent, so every class routed by it must declare the same
     * branches, nesting and coefficients.
     *
     * @param other the structure to compare against
     * @return true when the two declarations describe the same subnetwork
     */
    public boolean sameAs(StateDepRouting other) {
        if (other == null) {
            return false;
        }
        if (this.entry != other.entry || this.departure != other.departure) {
            return false;
        }
        if (this.branch.length != other.branch.length || this.C.length != other.C.length) {
            return false;
        }
        for (int b = 1; b < this.branch.length; b++) {
            if (!java.util.Arrays.equals(this.branch[b], other.branch[b])) {
                return false;
            }
            if (this.entryOf[b] != other.entryOf[b] || this.departureOf[b] != other.departureOf[b]) {
                return false;
            }
        }
        if (!java.util.Arrays.equals(this.level, other.level)) {
            return false;
        }
        if (!java.util.Arrays.equals(this.C, other.C)) {
            return false;
        }
        for (int t = 0; t < this.d.length; t++) {
            if (!java.util.Arrays.equals(this.d[t], other.d[t])) {
                return false;
            }
        }
        return true;
    }

    /** Deep copy, sharing the declared node objects. */
    public StateDepRouting copy() {
        StateDepRouting out = new StateDepRouting();
        out.entry = this.entry;
        out.departure = this.departure;
        out.branch = new int[this.branch.length][];
        for (int b = 1; b < this.branch.length; b++) {
            out.branch[b] = this.branch[b].clone();
        }
        out.entryOf = this.entryOf.clone();
        out.departureOf = this.departureOf.clone();
        out.level = this.level.clone();
        out.C = this.C.clone();
        out.d = new double[this.d.length][];
        for (int t = 0; t < this.d.length; t++) {
            out.d[t] = this.d[t].clone();
        }
        out.entryNode = this.entryNode;
        out.departureNode = this.departureNode;
        if (this.branchNodes != null) {
            out.branchNodes = new ArrayList<List<Node>>(this.branchNodes);
        }
        return out;
    }
}
