/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * A fork firing synchronization on an FJ tag-augmented struct (see
 * ModelAdapter.fjtag). The firing atomically consumes one parent job of
 * class jobclass held at the (stateful) Fork node and emits weight
 * (tasksPerLink) siblings per branch, in the auxiliary classes of the
 * entry's tag, at the branch head nodes.
 */
public class FJSync implements Serializable {

    public Event active;
    public int fork;       // fork node index (0-based)
    public int join;       // join node index (0-based)
    public int jobclass;   // original class index (0-based)
    public int tag;        // origin-job slot (0-based)
    public int[] branchheads; // branch head node indices (0-based)
    public int[] auxclasses;  // auxiliary class per branch for this tag (0-based)
    public Matrix auxall;     // B x T auxiliary class indices for the tag-occupancy scan
    public double weight;     // tasksPerLink: siblings emitted per branch
    public double prob;

    public FJSync() {
        this.weight = 1.0;
        this.prob = 1.0;
    }
}
