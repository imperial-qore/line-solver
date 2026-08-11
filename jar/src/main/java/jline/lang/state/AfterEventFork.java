/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * Event handler for stateful Fork nodes (FJ tag-augmented structs only,
 * see ModelAdapter.fjtag). The fork state is a per-class count of parent
 * jobs momentarily held before the fork firing. Arrivals are buffered
 * here; the atomic multi-branch emission is not a DEP event but a fork
 * firing synchronization (sn.fjsync) handled by AfterFJEvent.
 */
public class AfterEventFork implements Serializable {

    static Ret.EventResult afterEventFork(NetworkStruct sn, int ind, EventType event, int jobClass, boolean isSimulation, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar) {
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1.0);

        if (event == EventType.ARV) {
            Matrix srv = spaceSrv.copy();
            for (int row = 0; row < srv.getNumRows(); row++) {
                srv.set(row, jobClass, srv.get(row, jobClass) + 1);
            }
            outspace = Matrix.concatColumns(srv, spaceVar, null);
            outrate = new Matrix(outspace.getNumRows(), 1);
            outrate.ones();
            outrate.scaleEq(-1); // passive action, rate is unspecified
            outprob = new Matrix(outspace.getNumRows(), 1);
            outprob.ones();
        }
        // DEP from a Fork occurs only through sn.fjsync firings (see refreshSync)
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}
