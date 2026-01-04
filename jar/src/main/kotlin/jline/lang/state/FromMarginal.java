package jline.lang.state;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.InputOutputKt;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.*;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Station;
import jline.lang.nodes.Transition;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.UniqueRowResult;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.Collections;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;
import static jline.lang.constant.SchedStrategy.FCFS;
import static jline.lang.state.State.spaceClosedSingle;

public class FromMarginal implements Serializable {
    public static Matrix fromMarginal(NetworkStruct sn, int ind, Matrix n) {

        // Generate states such that the marginal queue-lengths are as in vector n
        // n(r): number of jobs at the station in class r
        int R = sn.nclasses;
        Matrix S = sn.nservers;
        Matrix state = new Matrix(0, 0);
        Matrix space = new Matrix(0, 0);

        if (sn.isstation.get(ind) == 1) {
            boolean mapOrMMPP2ProcTypes = false;
            for (int i = 0; i < sn.procid.get((Station) sn.nodes.get(ind)).size(); i++) {
                if (sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(i)) == ProcessType.MAP || sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(i)) == ProcessType.MMPP2) {
                    mapOrMMPP2ProcTypes = true;
                }
            }
            if (mapOrMMPP2ProcTypes) {
                if (sn.nservers.get(ind) > 1) {
                    throw new RuntimeException("Multiserver MAP stations are not supported.");
                }
                if (sn.sched.get(((Station) sn.nodes.get(ind))) != FCFS && !sn.nodetype.contains(NodeType.Source)) {
                    throw new RuntimeException("Non-FCFS MAP stations are not supported.");
                }
            }
        }

        int ist = (int) sn.nodeToStation.get(ind);
        int isf = (int) sn.nodeToStateful.get(ind);

        if (sn.isstateful.get(ind, 0) == 1 && sn.isstation.get(ind, 0) == 0) {
            for (int r = 0; r < R; r++) {
                Matrix init_r = spaceClosedSingle(1, n.get(r));
                state = Matrix.cartesian(state, init_r);
            }
            return Matrix.cartesian(space, state);
        }

        Matrix phases = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
                phases.set(0, r, 0);
            } else {
                phases.set(0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
            }
        }

        boolean anyNGreaterThanClassCap = false;

        for (int nIndex = 0; nIndex < n.length(); nIndex++) {
            if (n.get(nIndex) > sn.classcap.get(ist, nIndex)) {
                anyNGreaterThanClassCap = true;
            }
        }

        if (sn.sched.get(((Station) sn.nodes.get(ind))) != SchedStrategy.EXT && anyNGreaterThanClassCap) {
            return space;
        }

        // Generate local-state space
        switch (sn.nodetype.get(ind)) {
            case Queue:
            case Delay:
            case Source:
                switch (sn.sched.get(sn.stations.get(ist))) {
                    case EXT:
                        // source node case, treated as an infinite pool of jobs in buffer and a server for each class
                        for (int r = 0; r < R; r++) {
                            Matrix init_r = new Matrix(0, 0);
                            if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty() && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) { // if service disabled
                                init_r = new Matrix(1, (int) phases.get(r));
                                init_r.zero();
                            } else {
                                init_r = spaceClosedSingle(phases.get(r), 1);
                            }
                            state = Matrix.cartesian(state, init_r);
                        }
                        space = Matrix.cartesian(space, state); // server part

                        Matrix ones = new Matrix(space.getNumRows(), 1);
                        ones.ones();
                        Matrix infBuffer = ones.mult(new Matrix(1, 1).fromArray2D(new double[][]{{Inf}}));
                        // Attach infinite buffer to all state spaces containing job distribution across servers
                        space = infBuffer.concatCols(space);
                        break;
                    case INF:
                    case PS:
                    case DPS:
                    case GPS:
                    case LPS:
                    case PSPRIO:
                    case DPSPRIO:
                    case GPSPRIO:
                        // In these policies we only track the jobs in the servers
                        for (int r = 0; r < R; r++) {
                            Matrix init_r = spaceClosedSingle(phases.get(r), n.get(r));
                            state = Matrix.cartesian(state, init_r);
                        }
                        space = Matrix.cartesian(space, state);
                        break;
                    case SIRO:
                    case LEPT:
                    case SEPT:
                    case SRPT:
                    case SRPTPRIO:
                        // In these policies we track an un-ordered buffer and the jobs in the servers.
                        // We build list of job classes in the node, with repetition
                        if (n.elementSum() <= S.get(ist)) {
                            // buffer will be empty as we have enough servers to handle all tasks in parallel
                            for (int r = 0; r < R; r++) {
                                Matrix init_r = spaceClosedSingle(phases.get(r), n.get(r));
                                state = Matrix.cartesian(state, init_r);
                            }
                            Matrix newStates = new Matrix(state.getNumRows(), R);
                            newStates.zero();
                            newStates = newStates.concatCols(state);
                            space = Matrix.cartesian(space, newStates);
                        } else {
                            Matrix si = Maths.multiChooseCon(n, S.get(ist)); // jobs of class r that are running
                            Matrix mi_buf = n.repmat(si.getNumRows(), 1).sub(1, si);
                            for (int k = 0; k < si.getNumRows(); k++) {
                                // determine number of class r jobs running in phase j
                                Matrix kstate = new Matrix(0, 0);
                                for (int r = 0; r < R; r++) {
                                    Matrix init_r = spaceClosedSingle(phases.get(r), si.get(k, r));
                                    kstate = Matrix.cartesian(kstate, init_r);
                                }
                                state = Matrix.extractRows(mi_buf, k, k + 1, null).repmat(kstate.getNumRows(), 1).concatCols(kstate);
                                if (space.isEmpty()) {
                                    space = state.copy();
                                } else {
                                    space = Matrix.concatRows(space, state, null);
                                }
                            }
                        }
                        break;
                    case FCFS:
                    case HOL:
                    case FCFSPRIO:
                    case LCFS:
                    case LCFSPRIO:
                        Matrix vi = new Matrix(0, 0);
                        Matrix mi = new Matrix(0, 0);
                        double sizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) + Maths.factln(sn.cap.get(ist));
                        sizeEstimator = FastMath.round(sizeEstimator / FastMath.log(10));
                        if (sizeEstimator > 3) {
                            // Large state space warning - equivalent to MATLAB line_warning
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                InputOutputKt.line_warning(mfilename(new Object() {}), 
                                    "State space size is very large: 1e" + (int) sizeEstimator + " states. " +
                                    "This may cause performance issues. Consider using a smaller model or force=true option.");
                            }
                        }
                        if (n.elementSum() == 0) {
                            // For empty queues, create proper state space
                            space = new Matrix(1, (int) (1 + phases.elementSum()));
                            space.zero(); // Initialize with zeros
                            if (!sn.nodetype.get(ind).equals(NodeType.Source)) {
                                for (int r = 0; r < R; r++) {
                                    switch (sn.procid.get(sn.stations.get((int) sn.nodeToStation.get(ind))).get(sn.jobclasses.get(r))) {
                                        case MAP:
                                        case MMPP2:
                                            List<Double> phasesRange = new ArrayList<Double>();
                                            for (double i = 1; i <= sn.phases.get(ind, r); i++) {
                                                phasesRange.add(i);
                                            }
                                            // no transpose as constructor creates a column vector
                                            space = Matrix.cartesian(space, new Matrix(phasesRange));
                                    }
                                }
                            }
                            return space;
                        }
                        vi = new Matrix(0, 0);
                        for (int r = 0; r < R; r++) {
                            if (n.get(0, r) > 0) {
                                Matrix newVi = new Matrix(1, vi.getNumCols() + (int) n.get(0, r));
                                for (int i = 0; i < vi.getNumCols(); i++) {
                                    newVi.set(0, i, vi.get(0, i));
                                }
                                for (int i = vi.getNumCols(); i < newVi.getNumCols(); i++) {
                                    newVi.set(0, i, r + 1);
                                }
                                vi = newVi.copy();
                            }
                        }
                        // gen permutation of their positions in the waiting buffer
                        mi = Maths.uniquePerms(vi);
                        Matrix mi_buf = new Matrix(0, 0);
                        // now generate server states
                        if (mi.isEmpty()) {
                            mi_buf = new Matrix(1, (int) FastMath.max(0, n.elementSum() - S.get(ist)));
                            mi_buf.zero();
                            state = new Matrix(1, R);
                            state.zero();
                            state = Matrix.cartesian(state, mi_buf.concatCols(state));
                        } else {
                            int numCols = (int) Maths.min(n.elementSum(), sn.cap.get(ist));
                            Matrix miClone = mi.copy();
                            mi = new Matrix(mi.getNumRows(), numCols);

                            for (int row = 0; row < miClone.getNumRows(); row++) {
                                for (int col = miClone.getNumCols() - numCols; col < miClone.getNumCols(); col++) {
                                    mi.set(row, col, miClone.get(row, col));
                                }
                            }
                            mi = Maths.uniqueAndSort(mi);

                            // mi_buf: class of job in buffer position i (0=empty)
                            int numColumnsRight = (int) Maths.max((mi.getNumCols() - S.get(ist)), 0);
                            Matrix right = new Matrix(mi.getNumRows(), numColumnsRight);
                            Matrix.extract(mi, 0, mi.getNumRows(), 0, numColumnsRight, right, 0, 0);
                            double x = sn.cap.get(ist);
                            int numColumnsLeft = (int) Maths.max(0, (Maths.min(n.elementSum(), sn.cap.get(ist)) - S.get(ist) - right.getNumCols()));
                            Matrix left = new Matrix(mi.getNumRows(), numColumnsLeft);
                            mi_buf = left.concatCols(right);
                            if (mi_buf.isEmpty()) {
                                mi_buf = new Matrix(mi.getNumRows(), 1);
                                mi_buf.zero();
                            }
                            // mi_srv: class of job running in server i
                            int numColsSrv = (int) Maths.max(S.get(ist), 1);
                            Matrix miSrv = new Matrix(mi.getNumRows(), numColsSrv);

                            int colForMiSrv = 0;
                            for (int row = 0; row < miSrv.getNumRows(); row++) {
                                colForMiSrv = 0;
                                for (int col = (int) Maths.max(mi.getNumCols() - S.get(ist), 0); col < mi.getNumCols(); col++) {
                                    miSrv.set(row, colForMiSrv, mi.get(row, col));
                                    colForMiSrv++;
                                }
                            }

                            // si: number of class r jobs that are running
                            Matrix si = new Matrix(miSrv.getNumRows(), R);
                            for (int k = 0; k < mi.getNumRows(); k++) {
                                Matrix miSrvKRow = Matrix.extractRows(miSrv, k, k + 1, null);
                                Matrix histRow = Maths.binHist(miSrvKRow, 1, R);
                                for (int j = 0; j < R; j++) {
                                    si.set(k, j, histRow.get(j));
                                }
                            }

                            for (int k = 0; k < si.getNumRows(); k++) {
                                // determine number of class r jobs running in phase
                                // j in server state mi_srv(k,:) and build state
                                Matrix kState = new Matrix(0, 0);
                                kState.zero();
                                List<Integer> map_cols = new ArrayList<Integer>();
                                for (int r = 0; r < R; r++) {
                                    Matrix init_r = spaceClosedSingle(phases.get(r), si.get(k, r));
                                    
                                    // Handle MAP/MMPP2 process types
                                    if (sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MAP || 
                                        sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MMPP2) {
                                        
                                        if (si.get(k, r) == 0) {
                                            // Create phase range [1:phases(r)]
                                            List<Double> phasesRange = new ArrayList<Double>();
                                            for (double i = 1; i <= phases.get(r); i++) {
                                                phasesRange.add(i);
                                            }
                                            init_r = Matrix.cartesian(init_r, new Matrix(phasesRange));
                                        } else {
                                            // Add zero column
                                            init_r = Matrix.cartesian(init_r, new Matrix(1, 1));
                                        }
                                        
                                        // Update init_r: if last element is 0, set it to first non-zero element index
                                        for (int i = 0; i < init_r.getNumRows(); i++) {
                                            if (init_r.get(i, init_r.getNumCols() - 1) == 0) {
                                                for (int j = 0; j < init_r.getNumCols(); j++) {
                                                    if (init_r.get(i, j) != 0) {
                                                        init_r.set(i, init_r.getNumCols() - 1, j + 1);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    kState = Matrix.cartesian(kState, init_r).copy();
                                    
                                    // Track MAP/MMPP2 column positions
                                    if (sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MAP || 
                                        sn.procid.get((Station) sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MMPP2) {
                                        map_cols.add(kState.getNumCols());
                                    }
                                }
                                
                                // Reorder kState columns so MAP/MMPP2 columns come at the end
                                if (!map_cols.isEmpty()) {
                                    List<Integer> otherCols = new ArrayList<Integer>();
                                    for (int i = 1; i <= kState.getNumCols(); i++) {
                                        if (!map_cols.contains(i)) {
                                            otherCols.add(i);
                                        }
                                    }
                                    List<Integer> newColOrder = new ArrayList<Integer>();
                                    newColOrder.addAll(otherCols);
                                    newColOrder.addAll(map_cols);
                                    
                                    if (newColOrder.size() == kState.getNumCols()) {
                                        Matrix reorderedKState = new Matrix(kState.getNumRows(), kState.getNumCols());
                                        for (int i = 0; i < kState.getNumRows(); i++) {
                                            for (int j = 0; j < newColOrder.size(); j++) {
                                                reorderedKState.set(i, j, kState.get(i, newColOrder.get(j) - 1));
                                            }
                                        }
                                        kState = reorderedKState;
                                    }
                                }
                                Matrix miBufreplicated = Matrix.extractRows(mi_buf, k, k + 1, null).repmat(kState.getNumRows(), 1);
                                miBufreplicated = miBufreplicated.concatCols(kState).copy();

                                if (state.isEmpty()) {
                                    state = miBufreplicated;
                                } else {
                                    state = Matrix.concatRows(state, miBufreplicated, null);
                                }
                            }
                        }
                        space = state;
                        break;
                    case LCFSPI:
                    case LCFSPR:
                    case LCFSPIPRIO:
                    case LCFSPRPRIO:
                        Matrix vi_lpr = new Matrix(0, 0);
                        Matrix mi_lpr = new Matrix(0, 0);
                        // sum(n) - 1 due to Maths.factln including + 1
                        double lcfsprSizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) + Maths.factln(sn.cap.get(ist));
                        lcfsprSizeEstimator = FastMath.round(lcfsprSizeEstimator / FastMath.log(10));
                        if (lcfsprSizeEstimator > 3) {
                            // Large state space warning - equivalent to MATLAB line_warning
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                InputOutputKt.line_warning(mfilename(new Object() {}), 
                                    "State space size is very large: 1e" + (int) lcfsprSizeEstimator + " states. " +
                                    "This may cause performance issues. Consider using a smaller model or force=true option.");
                            }
                        }

                        if (n.elementSum() == 0) {
                            Matrix newSpace = new Matrix(1, (int) (1 + phases.elementSum()));
                            newSpace.zero();
                            space = newSpace;
                            return space;
                        }
                        // Similar to FCFS/HOL/LCFS case we track an ordered buffer and the jobs in the servers
                        // but in this case due to pre-emption jobs in buffer can be in not initial phase

                        // build list of job classes in the node, with repetition

                        vi_lpr = new Matrix(0, 0);
                        for (int r = 0; r < R; r++) {
                            if (n.get(0, r) > 0) {
                                Matrix newVi = new Matrix(1, vi_lpr.getNumCols() + (int) n.get(0, r));
                                for (int i = 0; i < vi_lpr.getNumCols(); i++) {
                                    newVi.set(0, i, vi_lpr.get(0, i));
                                }
                                for (int i = vi_lpr.getNumCols(); i < newVi.getNumCols(); i++) {
                                    newVi.set(0, i, r + 1);
                                }
                                vi_lpr = newVi.copy();
                            }
                        }
                        // gen permutation of their positions in the waiting buffer
                        mi_lpr = Maths.uniquePerms(vi_lpr);
                        // now generate server states
                        if (mi_lpr.isEmpty()) {
                            Matrix mi_buf_lpr = new Matrix(1, (int) Maths.max(0, n.elementSum() - S.get(ist)));
                            mi_buf_lpr.zero();
                            state = new Matrix(1, R);
                            state.zero();
                            state = Matrix.cartesian(state, mi_buf_lpr.concatCols(state));
                        } else {
                            int numCols = (int) Maths.min(n.elementSum(), sn.cap.get(ist));
                            Matrix miClone = mi_lpr.copy();
                            mi_lpr = new Matrix(mi_lpr.getNumRows(), numCols);

                            for (int row = 0; row < miClone.getNumRows(); row++) {
                                for (int col = miClone.getNumCols() - numCols; col < miClone.getNumCols(); col++) {
                                    mi_lpr.set(row, col, miClone.get(row, col));
                                }
                            }
                            // mi_buf: class of job in buffer position i (0 = empty)
                            int numColumnsRight = (int) Maths.max((mi_lpr.getNumCols() - S.get(ist)), 0);
                            Matrix right = new Matrix(mi_lpr.getNumRows(), numColumnsRight);
                            Matrix.extract(mi_lpr, 0, mi_lpr.getNumRows(), 0, numColumnsRight, right, 0, 0);
                            int numColumnsLeft = (int) Maths.max(0, (Maths.min(n.elementSum(), sn.cap.get(ist)) - S.get(ist) - right.getNumCols()));
                            Matrix left = new Matrix(mi_lpr.getNumRows(), numColumnsLeft);
                            Matrix mi_buf_lpr = left.concatCols(right);
                            if (mi_buf_lpr.isEmpty()) {
                                mi_buf_lpr = new Matrix(mi_lpr.getNumRows(), 1);
                                mi_buf_lpr.zero();
                            }
                            // miSrv: class of job running in server i
                            int numColsSrv = (int) Maths.max(S.get(ist), 1);
                            Matrix miSrv = new Matrix(mi_lpr.getNumRows(), numColsSrv);

                            int colForMiSrv = 0;
                            for (int row = 0; row < miSrv.getNumRows(); row++) {
                                colForMiSrv = 0;
                                for (int col = mi_lpr.getNumCols() - numColsSrv; col < mi_lpr.getNumCols(); col++) {
                                    miSrv.set(row, colForMiSrv, mi_lpr.get(row, col));
                                    colForMiSrv++;
                                }
                            }

                            // si: number of class r jobs that are running
                            Matrix si = new Matrix(miSrv.getNumRows(), R);
                            for (int k = 0; k < mi_lpr.getNumRows(); k++) {
                                Matrix miSrvKRow = Matrix.extractRows(miSrv, k, k + 1, null);
                                Matrix histRow = Maths.binHist(miSrvKRow, 1, R);
                                for (int j = 0; j < R; j++) {
                                    si.set(k, j, histRow.get(j));
                                }
                            }
                            for (int k = 0; k < si.getNumRows(); k++) {
                                // determine number of class r jobs running in phase j
                                // in server state miSrv(k, :) and build state
                                Matrix kState = new Matrix(0, 0);
                                for (int r = 0; r < R; r++) {
                                    kState = Matrix.cartesian(kState, spaceClosedSingle(phases.get(r), si.get(k, r)));
                                }
                                // generate job phases for all buffer states since we have pre-emption
                                Matrix bkState = new Matrix(0, 0);

                                Matrix jobsInBuffer = Matrix.extractRows(mi_buf_lpr, k, k + 1, null);
                                for (int j = 0; j < jobsInBuffer.length(); j++) {
                                    double job = jobsInBuffer.get(j);
                                    if (job > 0) {
                                        List<Double> phasesJRange = new ArrayList<Double>();
                                        for (double i = 1; i <= phases.get((int) job - 1); i++) {
                                            phasesJRange.add(i);
                                        }
                                        // no transpose as constructor makes column vector
                                        bkState = Matrix.cartesian(bkState, new Matrix(phasesJRange));
                                    }
                                    // Note: when job == 0 (empty buffer position), we skip phase enumeration
                                    // as there's no job in that position (matching MATLAB behavior)
                                }
                                Matrix bufStateTmp = Matrix.cartesian(Matrix.extractRows(mi_buf_lpr, k, k + 1, null), bkState);
                                // here we interleave positions of class and phases in buffer
                                Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                                bufState.zero();

                                // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                                int colForBufStateTmp = 0;
                                for (int row = 0; row < bufState.getNumRows(); row++) {
                                    for (int col = 0; col < bufState.getNumCols(); col += 2) {
                                        if (colForBufStateTmp < mi_buf_lpr.getNumCols()) {
                                            bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                            colForBufStateTmp++;
                                        }
                                    }
                                    colForBufStateTmp = 0;
                                }
                                colForBufStateTmp = mi_buf_lpr.getNumCols();
                                for (int row = 0; row < bufState.getNumRows(); row++) {
                                    for (int col = 1; col < bufState.getNumCols(); col += 2) {
                                        if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                                            bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                            colForBufStateTmp++;
                                        }
                                    }
                                    colForBufStateTmp = mi_buf_lpr.getNumCols();
                                }
                                if (state.isEmpty()) {
                                    state = Matrix.cartesian(bufState, kState);
                                } else {
                                    state = Matrix.concatRows(state, Matrix.cartesian(bufState, kState), null);
                                }
                            }
                        }
                        space = state;
                        break;
                    case SJF:
                    case LJF:
                        // in these policies the state space includes continuous
                        // random variables for the service times
                        throw new RuntimeException("The scheduling policy does not admit a discrete state space.");
                }

                for (int r = 0; r < R; r++) {
                    RoutingStrategy routingStrategy = sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r));
                    if (routingStrategy == RoutingStrategy.RROBIN || routingStrategy == RoutingStrategy.WRROBIN) {
                        // Get the outlinks for round-robin routing
                        Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(r));
                        if (outlinks != null) {
                            // Convert outlinks to column vector if needed
                            Matrix outlinkVector = new Matrix(outlinks.getNumRows() * outlinks.getNumCols(), 1);
                            int idx = 0;
                            for (int i = 0; i < outlinks.getNumRows(); i++) {
                                for (int j = 0; j < outlinks.getNumCols(); j++) {
                                    outlinkVector.set(idx++, 0, outlinks.get(i, j));
                                }
                            }
                            space = Matrix.cartesian(space, outlinkVector);
                        }
                    }
                }
                break;
            case Cache:
                switch (sn.sched.get(ist)) {
                    case INF:
                        // in this policy we only track the jobs in the servers
                        for (int r = 0; r < R; r++) {
                            Matrix init_r = State.spaceClosedSingle(phases.get(r), (int) n.get(r));
                            state = Matrix.cartesian(state, init_r);
                        }
                        space = Matrix.cartesian(space, state);
                        break;
                    default:
                        // For other scheduling policies in cache nodes
                        throw new RuntimeException("Unsupported scheduling policy for Cache node: " + sn.sched.get(ist));
                }
                
                // Handle round-robin routing for cache nodes
                for (int r = 0; r < R; r++) {
                    RoutingStrategy routingStrategy = sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r));
                    if (routingStrategy == RoutingStrategy.RROBIN) {
                        // Get the outlinks for round-robin routing
                        Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(r));
                        if (outlinks != null) {
                            // Convert outlinks to column vector if needed
                            Matrix outlinkVector = new Matrix(outlinks.getNumRows() * outlinks.getNumCols(), 1);
                            int idx = 0;
                            for (int i = 0; i < outlinks.getNumRows(); i++) {
                                for (int j = 0; j < outlinks.getNumCols(); j++) {
                                    outlinkVector.set(idx++, 0, outlinks.get(i, j));
                                }
                            }
                            space = Matrix.cartesian(space, outlinkVector);
                        }
                    }
                }
                break;
        }

        // Required to sort empty state as first
        List<Matrix> uniqueRows = new ArrayList<Matrix>();
        for (int i = 0; i < space.getNumRows(); i++) {
            Matrix tmp = new Matrix(1, space.getNumCols());
            Matrix tmp2 = new Matrix(1, space.getNumCols());
            Matrix.extractRows(space, i, i + 1, tmp);
            boolean unique = true;
            for (int j = i + 1; j < space.getNumRows(); j++) {
                Matrix.extractRows(space, j, j + 1, tmp2);
                if (tmp.isEqualTo(tmp2)) {
                    unique = false;
                }
            }
            if (unique) {
                uniqueRows.add(tmp);
            }
        }

        Comparator<Matrix> lexico = (mat1, mat2) -> {
            for (int col = 0; col < mat1.getNumCols(); col++) {
                int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
                if (result != 0) {
                    return result;
                }
            }
            return 0;
        };

        uniqueRows.sort(lexico);
        Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
        // So that states with jobs in phase 1 comes earlier
        int row = 0;
        for (int i = uniqueRows.size() - 1; i >= 0; i--) {
            for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
                newSpace.set(row, j, uniqueRows.get(i).get(0, j));
            }
            row++;
        }

        return newSpace;
    }

    public static Matrix fromMarginalAndRunning(NetworkStruct sn, int ind, Matrix n, Matrix s) {
        return fromMarginalAndRunning(sn, ind, n, s, true);
    }

    public static Matrix fromMarginalAndRunning(Network sn, int ind, Matrix n, Matrix s) {
        return fromMarginalAndRunning(sn.getStruct(true), ind, n, s, true);
    }

    public static Matrix fromMarginalAndRunning(NetworkStruct sn, int ind, Matrix n, Matrix s, boolean optionsForce) {
        int ist = (int) sn.nodeToStation.get(ind);
        int isf = (int) sn.nodeToStateful.get(ind);

        // generate one initial state such that the marginal queue-lengths are as in vector n
        // n(r): number of jobs at the station in class r
        // s(r): jobs of class r that are running
        int R = sn.nclasses;
        Matrix S = sn.nservers;
        Matrix K = new Matrix(1, R);

        for (int r = 0; r < R; r++) {
            if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
                K.set(0, r, 0);
            } else {
                K.set(0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
            }
        }
        Matrix state = new Matrix(0, 0);
        Matrix space = new Matrix(0, 0);
        LinkedList<Integer> exceeded = new LinkedList<Integer>();
        for (int i = 0; i < sn.classcap.getNumCols(); i++) {
            if (n.get(0, i) > sn.classcap.get(ist, i)) {
                exceeded.add(i);
            }
        }
        if (!exceeded.isEmpty()) {
            for (Integer r : exceeded) {
                if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty() && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
                    InputOutputKt.line_warning(mfilename(new Object() {
                    }), "State vector at station " + ist + " exceeds the class capacity. Some service classes are disabled.\n");
                } else {
                    InputOutputKt.line_warning(mfilename(new Object() {
                    }), "State vector at station " + ist + " exceeds the class capacity.\n");
                }
            }
            return space;
        }

        if (sn.nservers.get(ist, 0) > 0 && s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols()) > sn.nservers.get(ist, 0)) {
            return space;
        }

        if ((sn.nodetype.get(ind) == NodeType.Queue) || (sn.nodetype.get(ind) == NodeType.Delay) || (sn.nodetype.get(ind) == NodeType.Source)) {
            switch (sn.sched.get(sn.stations.get(ist))) {
                case EXT:
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(r), 0);
                        if (Utils.isInf(sn.njobs.get(r))) {
                            if (Double.isNaN(sn.rates.get(ist, r))) {
                                init.set(0, 0, 0); // class is not processed at this source
                            } else {
                                init.set(0, 0, 1);
                            }
                        }
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    Matrix ones = new Matrix(space.getNumRows(), 1);
                    ones.ones();
                    Matrix infBuffer = ones.mult(new Matrix(1, 1).fromArray2D(new double[][]{{Inf}}));
                    space = infBuffer.concatCols(space);
                    break;
                case INF:
                case PS:
                case DPS:
                case GPS:
                    // in these policies we only track the jobs in the servers
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(r), n.get(r));
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    break;
                case SIRO:
                case LEPT:
                case SEPT:
                case SRPT:
                case SRPTPRIO:
                    // in these policies we track an un-ordered buffer and the jobs in the servers
                    // we build a list of job classes in the node with repetition
                    if (n.elementSum() <= S.get(ist)) {
                        for (int r = 0; r < R; r++) {
                            Matrix init = spaceClosedSingle(K.get(r), n.get(r));
                            state = Matrix.cartesian(state, init);
                        }
                        Matrix newStates = new Matrix(state.getNumRows(), R);
                        newStates.zero();
                        newStates = newStates.concatCols(state);
                        space = Matrix.cartesian(space, newStates);
                    } else {
                        Matrix si = s.copy();
                        Matrix mi_buf = n.repmat(si.getNumRows(), 1).sub(1, si); // jobs of class r in buffer
                        for (int k = 0; k < si.getNumRows(); k++) {
                            Matrix kstate = new Matrix(0, 0);
                            for (int r = 0; r < R; r++) {
                                Matrix init = spaceClosedSingle(K.get(r), si.get(k, r));
                                kstate = Matrix.cartesian(kstate, init);
                            }
                            state = Matrix.extractRows(mi_buf, k, k + 1, null).repmat(kstate.getNumRows(), 1).concatCols(kstate);
                            if (space.isEmpty()) {
                                space = state.copy();
                            } else {
                                space = Matrix.concatRows(space, state, null);
                            }
                        }
                    }
                    break;
                case FCFS:
                case HOL:
                case FCFSPRIO:
                case LCFS:
                case LCFSPRIO:
                    double sizeEstimator = Maths.multinomialln(n.sub(1, s));
                    sizeEstimator = FastMath.round(sizeEstimator / FastMath.log(10));
                    if (sizeEstimator > 2) {
                        if (!optionsForce) {
                            InputOutputKt.line_error(mfilename(new Object() {
                            }), "State space size is very large: 1e" + FastMath.round(sizeEstimator / FastMath.log(10)) + " states. Stopping execution. " + "Set options.force=true to bypass this control.\n");
                        }
                    }
                    if (n.elementSum() == 0) {
                        space = new Matrix(1, (int) (1 + K.elementSum()));
                        space.zero();
                        return space;
                    }
                    // in these policies we track an ordered buffer and
                    // the jobs in the servers

                    // build list of job classes in the buffer, with repetition
                    Matrix inbuf = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(0, r) > s.get(r)) {
                            int numNewCols = (int) (n.get(0, r) - s.get(0, r));
                            Matrix newInBuf = new Matrix(1, inbuf.getNumCols() + numNewCols);
                            for (int i = 0; i < inbuf.getNumCols(); i++) {
                                newInBuf.set(0, i, inbuf.get(0, i));
                            }
                            for (int i = inbuf.getNumCols(); i < newInBuf.getNumCols(); i++) {
                                newInBuf.set(0, i, r + 1);
                            }
                            inbuf = newInBuf.copy();
                        }
                    }

                    // gen permutation of their positions in the fcfs buffer
                    Matrix mi = Maths.uniquePerms(inbuf);
                    if (mi.isEmpty()) {
                        Matrix mi_buf = new Matrix(1, (int) Maths.max(0, n.elementSum() - S.get(ist)));
                        state = new Matrix(1, R);
                        state.zero();
                        state = Matrix.cartesian(state, mi_buf.concatCols(state));
                    } else {
                        // mi_buf: class of job in buffer position i (0=empty)
                        Matrix mi_buf = new Matrix(0, 0);
                        double sumN = n.elementSum();
                        double sums = s.elementSum();
                        if (sumN > sums) {
                            mi_buf = new Matrix(mi.getNumRows(), (int) sumN - (int) sums);
                            for (int row = 0; row < mi_buf.getNumRows(); row++) {
                                for (int col = 0; col < sumN - sums; col++) {
                                    mi_buf.set(row, col, mi.get(row, col));
                                }
                            }
                        } else {
                            mi_buf = new Matrix(1, 1);
                            mi_buf.set(0, 0, 0);
                        }

                        // si: number of class r jobs that are running
                        Matrix si = s.copy();
                        for (int b = 0; b < mi_buf.getNumRows(); b++) {
                            for (int k = 0; k < si.getNumRows(); k++) {
                                Matrix kstate = new Matrix(0, 0);
                                for (int r = 0; r < R; r++) {
                                    Matrix init = spaceClosedSingle(K.get(r), si.get(k, r));
                                    kstate = Matrix.cartesian(kstate, init);
                                }
                                Matrix miBufRep = Matrix.extractRows(mi_buf, b, b + 1, null).repmat(kstate.getNumRows(), 1);
                                miBufRep = miBufRep.concatCols(kstate);
                                if (state.isEmpty()) {
                                    state = miBufRep;
                                } else {
                                    state = Matrix.concatRows(state, miBufRep, null);
                                }
                            }
                        }
                    }
                    space = state;
                    break;

                case LCFSPI:
                case LCFSPR:
                case LCFSPIPRIO:
                case LCFSPRPRIO:
                    double sizeEstimatorLPR = Maths.multinomialln(n.sub(1, s));
                    sizeEstimatorLPR = FastMath.round(sizeEstimatorLPR / FastMath.log(10));
                    if (sizeEstimatorLPR > 2) {
                        if (!optionsForce) {
                            System.err.format("State space size is very large: 1e%d states. Stopping execution. Set options = true," + "to bypass this control.\n", FastMath.round(sizeEstimatorLPR / FastMath.log(10)));
                        }
                    }
                    if (n.elementSum() == 0) {
                        space = new Matrix(1, (int) (1 + K.elementSum()));
                        space.zero();
                        return space;
                    }
                    // in these policies we track an ordered buffer and the jobs in the servers

                    // build list of job classes in the buffer with repetition
                    Matrix inbufLpr = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(r) > s.get(r)) {
                            int numNewCols = (int) (n.get(r) - s.get(r));
                            Matrix newInBuf = new Matrix(1, inbufLpr.getNumCols() + numNewCols);
                            for (int i = 0; i < inbufLpr.getNumCols(); i++) {
                                newInBuf.set(0, i, inbufLpr.get(0, i));
                            }
                            for (int i = inbufLpr.getNumCols(); i < newInBuf.getNumCols(); i++) {
                                newInBuf.set(0, i, r + 1);
                            }
                            inbufLpr = newInBuf.copy();
                        }
                    }

                    // gen permutation of their positions in the FCFS buffer
                    Matrix miLpr = Maths.uniquePerms(inbufLpr);
                    if (miLpr.isEmpty()) {
                        Matrix mi_buf = new Matrix(1, (int) FastMath.max(0, n.elementSum() - S.get(ist)));
                        mi_buf.zero();
                        state = new Matrix(1, R);
                        state.zero();
                        state = Matrix.cartesian(state, mi_buf.concatCols(state));
                    } else {
                        // mi_buf: class of job in buffer position i (0=empty)
                        Matrix mi_buf = new Matrix(0, 0);
                        if (n.elementSum() > s.elementSum()) {
                            double sumN = n.elementSum();
                            double sums = s.elementSum();
                            mi_buf = new Matrix(miLpr.getNumRows(), (int) sumN - (int) sums);
                            for (int row = 0; row < mi_buf.getNumRows(); row++) {
                                for (int col = 0; col < sumN - sums; col++) {
                                    mi_buf.set(row, col, miLpr.get(row, col));
                                }
                            }
                        } else {
                            mi_buf = new Matrix(1, 1);
                            mi_buf.set(0, 0, 0);
                        }

                        // si: number of class r jobs that are running
                        Matrix si = s.copy();
                        for (int b = 0; b < mi_buf.getNumRows(); b++) {
                            for (int k = 0; k < si.getNumRows(); k++) {
                                Matrix kState = new Matrix(0, 0);
                                for (int r = 0; r < R; r++) {
                                    Matrix init = spaceClosedSingle(K.get(r), si.get(k, r));
                                    kState = Matrix.cartesian(kState, init);
                                }
                                Matrix bkState = new Matrix(0, 0);
                                Matrix jobsInBuffer = Matrix.extractRows(mi_buf, b, b + 1, null);
                                for (int j = 0; j < jobsInBuffer.length(); j++) {
                                    double job = jobsInBuffer.get(j);
                                    if (job > 0) {
                                        List<Double> phasesJRange = new ArrayList<Double>();
                                        for (double i = 1; i <= K.get((int) job - 1); i++) {
                                            phasesJRange.add(i);
                                        }
                                        // no transpose as constructor makes column vector
                                        bkState = Matrix.cartesian(bkState, new Matrix(phasesJRange));
                                    } else {
                                        bkState = new Matrix(1, 1);
                                        bkState.zero();
                                    }
                                }
                                Matrix bufStateTmp = Matrix.cartesian(Matrix.extractRows(mi_buf, b, b + 1, null), bkState);
                                // here we interleave positions of class and phases in buffer
                                Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                                bufState.zero();

                                // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                                int colForBufStateTmp = 0;
                                for (int row = 0; row < bufState.getNumRows(); row++) {
                                    for (int col = 0; col < bufState.getNumCols(); col += 2) {
                                        if (colForBufStateTmp < mi_buf.getNumCols()) {
                                            bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                            colForBufStateTmp++;
                                        }
                                    }
                                    colForBufStateTmp = 0;
                                }
                                colForBufStateTmp = mi_buf.getNumCols();
                                for (int row = 0; row < bufState.getNumRows(); row++) {
                                    for (int col = 1; col < bufState.getNumCols(); col += 2) {
                                        if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                                            bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                            colForBufStateTmp++;
                                        }
                                    }
                                    colForBufStateTmp = mi_buf.getNumCols();
                                }
                                if (state.isEmpty()) {
                                    state = Matrix.cartesian(bufState, kState);
                                } else {
                                    state = Matrix.concatRows(state, Matrix.cartesian(bufState, kState), null);
                                }
                            }
                        }
                    }
                    space = state;
                    break;

                case SJF:
                case LJF:
                    // in these policies the state space includes continuous random variables
                    // for the service times
                    System.err.format("The scheduling policy does not admit a discrete state space");

            }


        } else if (sn.nodetype.get(ind) == NodeType.Cache) {
            // Handle cache node - only supports INF scheduling strategy
            switch (sn.sched.get(sn.stations.get(ist))) {
                case INF:
                    // In this policy we only track the jobs in the servers
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(r), n.get(r));
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    break;
                default:
                    throw new RuntimeException("Cache nodes only support INF scheduling strategy");
            }
        }

        //Required to sort empty state as first
        List<Matrix> uniqueRows = new ArrayList<Matrix>();
        for (int i = 0; i < space.getNumRows(); i++) {
            Matrix tmp = new Matrix(1, space.getNumCols());
            Matrix tmp2 = new Matrix(1, space.getNumCols());
            Matrix.extractRows(space, i, i + 1, tmp);
            boolean unique = true;
            for (int j = i + 1; j < space.getNumRows(); j++) {
                Matrix.extractRows(space, j, j + 1, tmp2);
                if (tmp.isEqualTo(tmp2)) {
                    unique = false;
                }
            }
            if (unique) {
                uniqueRows.add(tmp);
            }
        }

        Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
        // this ensures that states where jobs start in phase 1 are first, which is used eg
        // in SSA
        int row = 0;
        Comparator<Matrix> lexico = (mat1, mat2) -> {
            for (int col = 0; col < mat1.getNumCols(); col++) {
                int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
                if (result != 0) {
                    return result;
                }
            }
            return 0;
        };

        uniqueRows.sort(lexico);

        for (int i = uniqueRows.size() - 1; i >= 0; i--) {
            for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
                newSpace.set(row, j, uniqueRows.get(i).get(0, j));
            }
            row++;
        }

        return newSpace;
    }

    public static Matrix fromMarginalAndStarted(NetworkStruct sn, int ind, Matrix n, Matrix s) {
        return fromMarginalAndStarted(sn, ind, n, s, true);
    }

    public static Matrix fromMarginalAndStarted(Network network, int ind, Matrix n, Matrix s) {
        return fromMarginalAndStarted(network.getStruct(true), ind, n, s, true);
    }

    public static Matrix fromMarginalAndStarted(NetworkStruct sn, int ind, Matrix n, Matrix s, Boolean optionsForce) {
        // generate one initial state such that the marginal queue-lengths are as in vector n
        // n(r): number of jobs at the station in class r
        // s(r): jobs of class r that are running
        int R = sn.nclasses;
        Matrix S = sn.nservers;
        int ist = (int) sn.nodeToStation.get(ind);
        Matrix K = new Matrix(1, R);

        if (!Double.isNaN(sn.nodeToStation.get(ind))) {
            for (int r = 0; r < R; r++) {
                if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
                    K.set(0, r, 0);
                } else {
                    K.set(0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
                }
            }
        } else {
            if (sn.nodes.get(ind) instanceof Transition) {
                K = Matrix.zeros(1, ((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).nmodes);
                for (int m = 0; m < ((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).nmodes; m++) {
                    if (((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).firingproc.isEmpty()) {
                        K.set(0, m, 0);
                    } else {
                        K.set(0, m, ((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).firingproc.get(m).get(0).length());
                    }
                }
            }
        }

        Matrix state = new Matrix(0, 0);
        Matrix space = new Matrix(0, 0);
        LinkedList<Integer> exceeded = new LinkedList<Integer>();
        for (int i = 0; i < sn.classcap.getNumCols(); i++) {
            if ((sn.nodes.get(ist) instanceof Station) && n.get(0, i) > sn.classcap.get(ist, i)) {
                exceeded.add(i);
            }
        }
        if (!exceeded.isEmpty()) {
            for (Integer r : exceeded) {
                if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty() && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
                    System.err.format("State vector at station %d exceeds the class capacity. Some service classes are disabled.\n", ist);
                } else {
                    System.err.format("State vector at station %d exceeds the class capacity.\n", ist);
                }
            }
            return space;
        }

        if ((sn.nodes.get(ist) instanceof Station) && (sn.nservers.get(ist, 0) > 0 && s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols()) > sn.nservers.get(ist, 0))) {
            return space;
        }

        // Generate local-state space
        if ((sn.nodetype.get(ind) == NodeType.Queue) || (sn.nodetype.get(ind) == NodeType.Delay) || (sn.nodetype.get(ind) == NodeType.Source)) {
            switch (sn.sched.get(sn.stations.get(ist))) {
                case EXT:
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(0, r), 0);
                        if (Utils.isInf(sn.njobs.get(0, r))) {
                            if ((!sn.proc.isEmpty()) && (!sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
                                init.set(0, 0, 0); // class is not processed at this source
                            } else {
                                // init the job generation
                                init.set(0, 0, 1);
                            }
                        }
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    Matrix ones = new Matrix(space.getNumRows(), 1);
                    ones.ones();
                    Matrix infBuffer = ones.mult(new Matrix(1, 1).fromArray2D(new double[][]{{Inf}}));
                    space = infBuffer.concatCols(space);
                    break;
                case INF:
                case PS:
                case DPS:
                case GPS:
                case PSPRIO:
                case DPSPRIO:
                case GPSPRIO:
                    // In these policies we only track the jobs in the servers
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(0, r), 0);
                        init.set(0, 0, n.get(0, r));
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    break;

                case SIRO:
                case LEPT:
                case SEPT:
                case SRPT:
                case SRPTPRIO:
                    // In these policies we track an un-ordered buffer and the jobs in the servers build list
                    // of job classes in the node, with repetition
                    if (n.elementSum() <= S.get(ist)) {
                        for (int r = 0; r < R; r++) {
                            Matrix init = spaceClosedSingle(K.get(0, r), 0);
                            init.set(0, 0, n.get(0, r));
                            state = Matrix.cartesian(state, init);
                        }
                        Matrix newStates = new Matrix(state.getNumRows(), R);
                        newStates.zero();
                        newStates = newStates.concatCols(state);
                        space = Matrix.cartesian(space, newStates);
                    } else {
                        Matrix si = s.copy();
                        Matrix mi_buf = n.repmat(si.getNumRows(), 1).sub(1, si); // jobs of class r in buffer
                        for (int k = 0; k < si.getNumRows(); k++) {
                            Matrix kstate = new Matrix(0, 0);
                            for (int r = 0; r < R; r++) {
                                Matrix init = spaceClosedSingle(K.get(0, r), 0);
                                init.set(0, 0, si.get(k, r));
                                kstate = Matrix.cartesian(kstate, init);
                            }
                            state = Matrix.extractRows(mi_buf, k, k + 1, null).repmat(kstate.getNumRows(), 1).concatCols(kstate);
                            if (space.isEmpty()) {
                                space = state.copy();
                            } else {
                                space = Matrix.concatRows(space, state, null);
                            }
                        }
                    }
                    break;
                case FCFS:
                case HOL:
                case FCFSPRIO:
                case LCFS:
                case LCFSPRIO:
                    Matrix inbuf = new Matrix(0, 0);
                    double sizeEstimator = 0;
                    Matrix mi = new Matrix(0, 0);
                    Matrix mi_buf = new Matrix(0, 0);
                    Matrix mi_srv = new Matrix(0, 0);
                    if (n.elementSum() == 0) {
                        // For empty queues, return proper initial state with zeros
                        // Buffer part: single zero (empty buffer)
                        state = new Matrix(1, 1);
                        state.set(0, 0, 0);
                        // Server part: K.elementSum() zeros (idle servers)
                        for (int i = 0; i < K.elementSum(); i++) {
                            state = state.concatCols(Matrix.singleton(0));
                        }
                        return state;
                    }

                    // In these policies we track an ordered buffer and the jobs in the servers
                    // build list of job classes in the buffer, with repetition
                    inbuf = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(0, r) > 0) {
                            int numNewCols = (int) (n.get(0, r) - s.get(0, r));
                            Matrix newInBuf = new Matrix(1, inbuf.getNumCols() + numNewCols);
                            for (int i = 0; i < inbuf.getNumCols(); i++) {
                                newInBuf.set(0, i, inbuf.get(0, i));
                            }
                            for (int i = inbuf.getNumCols(); i < newInBuf.getNumCols(); i++) {
                                newInBuf.set(0, i, r + 1);
                            }
                            inbuf = newInBuf.copy();
                        }
                    }

                    sizeEstimator = Maths.multinomialln(n);
                    sizeEstimator = FastMath.round(sizeEstimator / FastMath.log(10));
                    if (sizeEstimator > 2) {
                        if (!optionsForce) {
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                System.err.format("State space size is large: 1e%d states. Cannot generate valid state space. Initializing station %d from a default state.\n", (int) sizeEstimator, ind);
                            }
                            state = inbuf.copy();
                            return state;
                        }
                    }

                    // Gen permutation of their positions in the FCFS buffer
                    mi = Maths.uniquePerms(inbuf);
                    double sumN = n.sumSubMatrix(0, n.getNumRows(), 0, n.getNumCols());
                    double sumS = s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols());
                    if (mi.isEmpty()) {
                        mi_buf = new Matrix(1, (int) FastMath.max(1, sumN - S.get(ist, 0)));
                        mi_buf.zero();
                        state = new Matrix(1, (int) K.sumSubMatrix(0, K.getNumRows(), 0, K.getNumCols()));
                        state.zero();
                        Matrix newState = new Matrix(1, mi_buf.getNumCols() + state.getNumCols());
                        for (int i = 0; i < mi_buf.getNumCols(); i++) {
                            newState.set(0, i, mi_buf.get(0, i));
                        }
                        for (int i = mi_buf.getNumCols(); i < newState.getNumCols(); i++) {
                            newState.set(0, i, state.get(0, i - mi_buf.getNumCols()));
                        }
                        state = newState.copy();
                    } else {
                        // mi_buf: class of job in buffer position i (0 = empty)
                        if (sumN > sumS) {
                            mi_buf = new Matrix(mi.getNumRows(), (int) sumN - (int) sumS);
                            for (int row = 0; row < mi_buf.getNumRows(); row++) {
                                for (int col = 0; col < sumN - sumS; col++) {
                                    mi_buf.set(row, col, mi.get(row, col));
                                }
                            }
                        } else {
                            mi_buf = new Matrix(1, 1);
                            mi_buf.set(0, 0, 0);
                        }
                    }

                    // mi_srv: class of jobs running in the server of i
                    mi_srv = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        Matrix new_mi_srv = new Matrix(1, mi_srv.getNumCols() + (int) s.get(0, r));
                        for (int i = 0; i < mi_srv.getNumCols(); i++) {
                            new_mi_srv.set(0, i, mi_srv.get(0, i));
                        }
                        for (int i = mi_srv.getNumCols(); i < new_mi_srv.getNumCols(); i++) {
                            new_mi_srv.set(0, i, r);
                        }
                        mi_srv = new_mi_srv.copy();
                    }

                    // si: number of class r jobs that are running
                    Matrix si = s.copy();
                    for (int b = 0; b < mi_buf.getNumRows(); b++) {
                        for (int k = 0; k < si.getNumRows(); k++) {
                            Matrix kState = new Matrix(0, 0);
                            for (int r = 0; r < R; r++) {
                                Matrix init = spaceClosedSingle(K.get(0, r), 0);
                                init.set(0, 0, si.get(k, r));
                                kState = Matrix.cartesian(kState, init);
                            }
                            Matrix miBufRep = Matrix.extractRows(mi_buf, b, b + 1, null).repmat(kState.getNumRows(), 1);
                            miBufRep = miBufRep.concatCols(kState);
                            if (state.isEmpty()) {
                                state = miBufRep;
                            } else {
                                state = Matrix.concatRows(state, miBufRep, null);
                            }

                        }
                    }
                    space = state;
                    break;
                case LCFSPI:
                case LCFSPR:
                case LCFSPIPRIO:
                case LCFSPRPRIO:
                    Matrix inbuf_lpr = new Matrix(0, 0);
                    double sizeEstimator_lpr = 0;
                    Matrix mi_lpr = new Matrix(0, 0);
                    Matrix mi_buf_lpr = new Matrix(0, 0);
                    Matrix mi_srv_lpr = new Matrix(0, 0);
                    if (n.elementSum() == 0) {
                        space = new Matrix(1, (int) (1 + K.elementSum()));
                        return space;
                    }
                    // in this policy we track an ordered buffer and the jobs in the servers
                    // build list of job classes in the buffer, with repetition

                    inbuf_lpr = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(0, r) > 0) {
                            int numNewCols = (int) (n.get(0, r) - s.get(0, r));
                            Matrix newInBuf = new Matrix(1, inbuf_lpr.getNumCols() + numNewCols);
                            for (int i = 0; i < inbuf_lpr.getNumCols(); i++) {
                                newInBuf.set(0, i, inbuf_lpr.get(0, i));
                            }
                            for (int i = inbuf_lpr.getNumCols(); i < newInBuf.getNumCols(); i++) {
                                newInBuf.set(0, i, r + 1);
                            }
                            inbuf_lpr = newInBuf.copy();
                        }
                    }
                    sizeEstimator_lpr = Maths.multinomialln(n);
                    sizeEstimator_lpr = FastMath.round(sizeEstimator_lpr / FastMath.log(10));
                    if (sizeEstimator_lpr > 2) {
                        if (!optionsForce) {
                            System.err.format("State space size is very large: 1e%f states. Cannot generate valid state space. Initializing station %d from a default state.\n", sizeEstimator_lpr, ind);
                            state = inbuf_lpr.copy();
                            return state;
                        }
                    }

                    // gen permutation of their positions in the buffer
                    mi_lpr = Maths.uniquePerms(inbuf_lpr);
                    mi_buf_lpr = new Matrix(0, 0);
                    if (mi_lpr.isEmpty()) {
                        mi_buf_lpr = new Matrix(1, (int) Maths.max(1, n.elementSum() - S.get(ist)));
                        state = new Matrix(1, (int) K.elementSum());
                        state.zero();
                        state = mi_buf_lpr.concatCols(state);
                    } else {
                        // mi_buf: class of job in buffer position i (0 = empty)
                        if (n.elementSum() > s.elementSum()) {
                            mi_buf_lpr = new Matrix(mi_lpr.getNumRows(), (int) (n.elementSum() - s.elementSum()));
                            for (int row = 0; row < mi_buf_lpr.getNumRows(); row++) {
                                for (int col = 0; col < mi_buf_lpr.getNumCols(); col++) {
                                    mi_buf_lpr.set(row, col, mi_lpr.get(row, col));
                                }
                            }
                        } else {
                            mi_buf_lpr = new Matrix(1, 1);
                            mi_buf_lpr.zero();
                        }
                    }

                    // mi_srv: class of jobs running in the server of i
                    mi_srv_lpr = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(0, r) > 0) {
                            Matrix newMiSrv = new Matrix(1, mi_srv_lpr.getNumCols() + (int) s.get(0, r));
                            for (int i = 0; i < mi_srv_lpr.getNumCols(); i++) {
                                newMiSrv.set(0, i, mi_srv_lpr.get(0, i));
                            }
                            for (int i = mi_srv_lpr.getNumCols(); i < newMiSrv.getNumCols(); i++) {
                                newMiSrv.set(0, i, r + 1);
                            }
                            mi_srv_lpr = newMiSrv.copy();
                        }
                    }


                    // si: number of class r jobs that are running
                    Matrix si_lpr = s.copy();

                    for (int b = 0; b < mi_buf_lpr.getNumRows(); b++) {
                        for (int k = 0; k < si_lpr.getNumRows(); k++) {
                            Matrix kState = new Matrix(0, 0);
                            for (int r = 0; r < R; r++) {
                                Matrix init = spaceClosedSingle(K.get(r), 0);
                                init.set(0, 0, si_lpr.get(k, r));
                                kState = Matrix.cartesian(kState, init);
                            }
                            Matrix bkState = new Matrix(0, 0);
                            Matrix jobsInBuffer = Matrix.extractRows(mi_buf_lpr, b, b + 1, null);
                            for (int j = 0; j < jobsInBuffer.length(); j++) {
                                double job = jobsInBuffer.get(j);
                                if (job > 0) {
                                    List<Double> phasesJRange = new ArrayList<Double>();
                                    for (double i = 1; i <= K.get((int) job - 1); i++) {
                                        phasesJRange.add(i);
                                    }
                                    // no transpose as constructor makes column vector
                                    bkState = Matrix.cartesian(bkState, new Matrix(phasesJRange));
                                } else {
                                    bkState = new Matrix(1, 1);
                                    bkState.zero();
                                }
                            }
                            Matrix bufStateTmp = Matrix.cartesian(Matrix.extractRows(mi_buf_lpr, b, b + 1, null), bkState);
                            // here we interleave positions of class and phases in buffer
                            Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                            bufState.zero();

                            // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                            int colForBufStateTmp = 0;
                            for (int row = 0; row < bufState.getNumRows(); row++) {
                                for (int col = 0; col < bufState.getNumCols(); col += 2) {
                                    if (colForBufStateTmp < mi_buf_lpr.getNumCols()) {
                                        bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                        colForBufStateTmp++;
                                    }
                                }
                                colForBufStateTmp = 0;
                            }
                            colForBufStateTmp = mi_buf_lpr.getNumCols();
                            for (int row = 0; row < bufState.getNumRows(); row++) {
                                for (int col = 1; col < bufState.getNumCols(); col += 2) {
                                    if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                                        bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                                        colForBufStateTmp++;
                                    }
                                }
                                colForBufStateTmp = mi_buf_lpr.getNumCols();
                            }
                            if (state.isEmpty()) {
                                state = Matrix.cartesian(bufState, kState);
                            } else {
                                state = Matrix.concatRows(state, Matrix.cartesian(bufState, kState), null);
                            }
                        }
                    }
                    space = state;
                    break;
                case POLLING:
                    // POLLING servers track an ordered buffer and jobs in the servers, similar to FCFS
                    // but may also need to track polling state/position
                    Matrix inbuf_poll = new Matrix(0, 0);
                    double sizeEstimator_poll = 0;
                    Matrix mi_poll = new Matrix(0, 0);
                    Matrix mi_buf_poll = new Matrix(0, 0);
                    Matrix mi_srv_poll = new Matrix(0, 0);
                    if (n.elementSum() == 0) {
                        space = new Matrix(1, (int) (1 + Maths.max(R, K.elementSum())));
                        return space;
                    }

                    // Build list of job classes in the buffer, with repetition
                    inbuf_poll = new Matrix(0, 0);
                    for (int r = 0; r < R; r++) {
                        if (n.get(0, r) > 0) {
                            int numNewCols = (int) (n.get(0, r) - s.get(0, r));
                            Matrix newInBuf = new Matrix(1, inbuf_poll.getNumCols() + numNewCols);
                            for (int i = 0; i < inbuf_poll.getNumCols(); i++) {
                                newInBuf.set(0, i, inbuf_poll.get(0, i));
                            }
                            for (int i = inbuf_poll.getNumCols(); i < newInBuf.getNumCols(); i++) {
                                newInBuf.set(0, i, r + 1);
                            }
                            inbuf_poll = newInBuf.copy();
                        }
                    }

                    sizeEstimator_poll = Maths.multinomialln(n);
                    sizeEstimator_poll = FastMath.round(sizeEstimator_poll / FastMath.log(10));
                    if (sizeEstimator_poll > 2) {
                        if (!optionsForce) {
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                System.err.format("State space size is large: 1e%d states. Cannot generate valid state space. Initializing station %d from a default state.\n", (int) sizeEstimator_poll, ind);
                            }
                            state = inbuf_poll.copy();
                            return state;
                        }
                    }

                    // Generate permutation of their positions in the polling buffer
                    mi_poll = Maths.uniquePerms(inbuf_poll);
                    double sumN_poll = n.sumSubMatrix(0, n.getNumRows(), 0, n.getNumCols());
                    double sumS_poll = s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols());
                    if (mi_poll.isEmpty()) {
                        mi_buf_poll = new Matrix(1, (int) FastMath.max(1, sumN_poll - S.get(ist, 0)));
                        mi_buf_poll.zero();
                        state = new Matrix(1, (int) K.sumSubMatrix(0, K.getNumRows(), 0, K.getNumCols()));
                        state.zero();
                        Matrix newState = new Matrix(1, mi_buf_poll.getNumCols() + state.getNumCols());
                        for (int i = 0; i < mi_buf_poll.getNumCols(); i++) {
                            newState.set(0, i, mi_buf_poll.get(0, i));
                        }
                        for (int i = mi_buf_poll.getNumCols(); i < newState.getNumCols(); i++) {
                            newState.set(0, i, state.get(0, i - mi_buf_poll.getNumCols()));
                        }
                        state = newState.copy();
                    } else {
                        // mi_buf_poll: class of job in buffer position i (0 = empty)
                        if (sumN_poll > sumS_poll) {
                            mi_buf_poll = new Matrix(mi_poll.getNumRows(), (int) sumN_poll - (int) sumS_poll);
                            for (int row = 0; row < mi_buf_poll.getNumRows(); row++) {
                                for (int col = 0; col < sumN_poll - sumS_poll; col++) {
                                    mi_buf_poll.set(row, col, mi_poll.get(row, col));
                                }
                            }
                        } else {
                            mi_buf_poll = new Matrix(1, 1);
                            mi_buf_poll.zero();
                        }

                        // mi_srv_poll: class of jobs running in the server
                        mi_srv_poll = new Matrix(0, 0);
                        for (int r = 0; r < R; r++) {
                            if (n.get(0, r) > 0) {
                                Matrix newMiSrv = new Matrix(1, mi_srv_poll.getNumCols() + (int) s.get(0, r));
                                for (int i = 0; i < mi_srv_poll.getNumCols(); i++) {
                                    newMiSrv.set(0, i, mi_srv_poll.get(0, i));
                                }
                                for (int i = mi_srv_poll.getNumCols(); i < newMiSrv.getNumCols(); i++) {
                                    newMiSrv.set(0, i, r + 1);
                                }
                                mi_srv_poll = newMiSrv.copy();
                            }
                        }

                        // si_poll: number of class r jobs that are running
                        Matrix si_poll = s.copy();

                        for (int b = 0; b < mi_buf_poll.getNumRows(); b++) {
                            for (int k = 0; k < si_poll.getNumRows(); k++) {
                                Matrix kState = new Matrix(0, 0);
                                for (int r = 0; r < R; r++) {
                                    Matrix init = spaceClosedSingle(K.get(r), 0);
                                    init.set(0, 0, si_poll.get(k, r));
                                    kState = Matrix.cartesian(kState, init);
                                }
                                Matrix bkState = new Matrix(0, 0);
                                Matrix jobsInBuffer = Matrix.extractRows(mi_buf_poll, b, b + 1, null);
                                for (int j = 0; j < jobsInBuffer.length(); j++) {
                                    double job = jobsInBuffer.get(j);
                                    if (job > 0) {
                                        List<Double> phasesJRange = new ArrayList<Double>();
                                        for (double i = 1; i <= K.get((int) job - 1); i++) {
                                            phasesJRange.add(i);
                                        }
                                        // no transpose as constructor makes column vector
                                        bkState = Matrix.cartesian(bkState, new Matrix(phasesJRange));
                                    } else {
                                        bkState = Matrix.cartesian(bkState, new Matrix(1, 1));
                                    }
                                }
                                if (state.isEmpty()) {
                                    state = Matrix.cartesian(bkState, kState);
                                } else {
                                    state = Matrix.concatRows(state, Matrix.cartesian(bkState, kState), null);
                                }
                            }
                        }
                    }
                    space = state;
                    break;
                case SJF:
                case LJF:
                    // In these policies the state space includes continuous random variables for the service
                    // times in these policies we only track the jobs in the servers
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(0, r), 0);
                        init.set(0, 0, n.get(0, r));
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    System.err.format("The schedyling policy does not admit a discrete state space");
                    break;
            }

        } else if (sn.nodetype.get(ind) == NodeType.Cache) {
            // Handle cache node - only supports INF scheduling strategy
            switch (sn.sched.get(sn.stations.get(ist))) {
                case INF:
                    // In this policy we only track the jobs in the servers
                    for (int r = 0; r < R; r++) {
                        Matrix init = spaceClosedSingle(K.get(0, r), 0);
                        init.set(0, 0, n.get(0, r));
                        state = Matrix.cartesian(state, init);
                    }
                    space = Matrix.cartesian(space, state);
                    break;
                default:
                    throw new RuntimeException("Cache nodes only support INF scheduling strategy");
            }
        } else if (sn.nodetype.get(ind) == NodeType.Transition) {
            InputOutputKt.line_error(mfilename(new Object() {
            }), "fromMarginalAndStarted cannot be used on Petri net elements");
        } else if (sn.nodetype.get(ind) == NodeType.Place) {
            InputOutputKt.line_error(mfilename(new Object() {
            }), "fromMarginalAndStarted cannot be used on Petri net elements");
        } else if (sn.nodetype.get(ind) == NodeType.Join) {
            space = new Matrix(1, 1);
            space.zero();
        }

        // Required to sort empty state as first
        List<Matrix> uniqueRows = new ArrayList<Matrix>();
        for (int i = 0; i < space.getNumRows(); i++) {
            Matrix tmp = new Matrix(1, space.getNumCols());
            Matrix tmp2 = new Matrix(1, space.getNumCols());
            Matrix.extractRows(space, i, i + 1, tmp);
            boolean unique = true;
            for (int j = i + 1; j < space.getNumRows(); j++) {
                Matrix.extractRows(space, j, j + 1, tmp2);
                if (tmp.isEqualTo(tmp2)) {
                    unique = false;
                }
            }
            if (unique) {
                uniqueRows.add(tmp);
            }
        }

        Comparator<Matrix> lexico = (mat1, mat2) -> {
            for (int col = 0; col < mat1.getNumCols(); col++) {
                int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
                if (result != 0) {
                    return result;
                }
            }
            return 0;
        };

        uniqueRows.sort(lexico);
        Matrix newSpace;
        if (uniqueRows.isEmpty()) {
            // Return empty matrix with appropriate dimensions
            newSpace = new Matrix(0, space.getNumCols());
        } else {
            newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
            // this ensures that states where jobs start in phase 1 are first, which is used eg
            // in SSA
            int row = 0;
            for (int i = uniqueRows.size() - 1; i >= 0; i--) {
                for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
                    newSpace.set(row, j, uniqueRows.get(i).get(0, j));
                }
                row++;
            }
        }

        return newSpace;
    }

    public static Matrix fromMarginalBounds(
            NetworkStruct sn, int ind, Matrix ub, double cap, SolverOptions options) {
        if (options == null) {
            options = Solver.defaultOptions();
        }

        double ist = sn.nodeToStation.get(ind);
        Matrix space = new Matrix(0, 0);
        //    space.set(0, 0, Inf);
        Matrix lb = ub.copy();
        lb.zero();
        int R = sn.nclasses;

        boolean isVectorLB = lb.length() != 1;
        boolean isVectorUB = ub.length() != 1;

        if (isVectorLB != isVectorUB) {
            line_error(mfilename(new Object[]{}), "Bounds must either be both vectors or both scalars");
        }

        if (isVectorUB && isVectorLB) {
            Matrix nmax = fromMarginal(sn, ind, ub);

            if (nmax.isEmpty()) {
                nmax = fromMarginal(sn, ind, ub);
            }
            // Generate population combinations following MATLAB lines 29-34
            Matrix n = PopulationLattice.pprodcon(null, lb, ub);
            boolean checkN = true;
            for (int row = 0; row < n.getNumRows(); row++) {
                for (int col = 0; col < n.getNumCols(); col++) {
                    if (n.get(row, col) == -1) {
                        checkN = false;
                    }
                }
            }
            while (checkN) {
                Matrix state = fromMarginal(sn, ind, n);

                int colNum = nmax.getNumCols();
                int originalRows = space.getNumRows();
                int cols = state.getNumCols();
                int originalCols = nmax.getNumCols() - state.getNumCols();
                space.expandMatrix(space.getNumRows() + state.getNumRows(), colNum, space.getNumNonZeros());


                for (int row = 0; row < state.getNumRows(); row++) {
                    for (int col = 0; col < cols; col++) {
//            System.out.println("rowIdx: " + (row + originalRows) + ", col: " + (col + originalCols) + " space size: " + space.getNumRows() + " * " + space.getNumCols());
//            System.out.println("space size: row - " + space.getNumRows() + "   col - " + space.getNumCols());
//            System.out.println("current row: " + originalRows + row + "   current col: " + col + originalCols);
//            System.out.println("nmax cols: " + nmax.getNumCols() + "   state cols: " + state.getNumCols());
                        space.set(originalRows + row, col + originalCols, state.get(row, col));
                    }
                }
                n = PopulationLattice.pprodcon(n, lb, ub);
                for (int row = 0; row < n.getNumRows(); row++) {
                    for (int col = 0; col < n.getNumCols(); col++) {
                        if (n.get(row, col) == -1) {
                            checkN = false;
                        }
                    }
                }
            }
        } else {
            if (ub.get(0, 0) >= lb.get(0, 0)) {
                for (int bi = (int) ub.get(0, 0); bi >= (int) lb.get(0, 0); bi--) {
                    Matrix nset = Maths.multichoose((double)R, (double)bi);
                    for (int j = 0; j < nset.getNumRows(); j++) {
                        Matrix state = fromMarginal(sn, ind, nset.getRow(j));
                        int originalRows = space.getNumRows();
                        int originalCols = space.getNumCols();
                        int stateRows = state.getNumRows();
                        int stateCols = state.getNumCols();

                        if (bi == (int) ub.get(0, 0) && j == 0) {
                            //              if (Double.isInfinite(space.get(0, 0))) {
                            //                space.expandMatrix(
                            //                    space.getNumRows() + state.getNumRows() - 1,
                            //                    state.getNumCols(),
                            //                    space.getNumNonZeros());
                            //                for (int row = 0; row < state.getNumRows(); row++) {
                            //                  for (int col = 0; col < stateCols; col++) {
                            //                    space.set(originalRows + row - state.getNumRows(), col,
                            // state.get(row, col));
                            //                  }
                            //                }
                            //              } else {
                            //                space.expandMatrix(
                            //                    space.getNumRows() + state.getNumRows(),
                            //                    space.getNumCols(),
                            //                    space.getNumNonZeros());
                            //                for (int row = 0; row < state.getNumRows(); row++) {
                            //                  for (int col = 0; col < originalRows; col++) {
                            //                    space.set(originalRows + row, col, state.get(row, col));
                            //                  }
                            //                }
                            //              }
                            space.expandMatrix(
                                    space.getNumRows() + state.getNumRows(),
                                    state.getNumCols(),
                                    space.getNumNonZeros());
                            for (int row = 0; row < state.getNumRows(); row++) {
                                for (int col = 0; col < stateCols; col++) {
                                    space.set(originalRows + row, col, state.get(row, col));
                                }
                            }

                        } else {
                            double newCols = Maths.max(originalRows, originalCols + stateCols);
                            space.expandMatrix(
                                    space.getNumRows() + state.getNumRows(),
                                    space.getNumCols(),
                                    space.getNumNonZeros());
                            for (int row = 0; row < stateRows; row++) {
                                for (int col = space.getNumCols() - state.getNumCols();
                                     col < space.getNumCols();
                                     col++) {
                                    space.set(
                                            space.getNumRows() + row - state.getNumRows(),
                                            col,
                                            state.get(row, col - (space.getNumCols() - state.getNumCols())));
                                }
                            }
                        }
                    }
                }
            }
        }

        UniqueRowResult uniqueResult = Matrix.uniqueRows(space);
        space = uniqueResult.sortedMatrix;
        if (sn.isstateful.get(ind, 0) == 1) {
            List<Integer> keep = new ArrayList<Integer>();
            for (int s = 0; s < space.getNumRows(); s++) {
                State.StateMarginalStatistics result =
                        ToMarginal.toMarginal(sn, ind, space.getRow(s), null, null, null, null, null);
                Matrix ni = result.ni;
                Matrix nir = result.nir;
                if (sn.isstation.get(ind, 0) == 1) {
                    boolean check = true;
                    Matrix compareMatrix = sn.classcap.getRow((int) ist);
                    for (int col = 0; col < compareMatrix.getNumCols(); col++) {
                        if (compareMatrix.get(0, col) > 1000000000) {
                            compareMatrix.set(0, col, Inf);
                        }
                    }
                    for (int col = 0; col < nir.getNumCols(); col++) {
                        if (nir.get(col) > compareMatrix.get(col)) {
                            check = false;
                            break;
                        }
                    }
                    for (int row = 0; row < ni.getNumRows(); row++) {
                        for (int col = 0; col < ni.getNumCols(); col++) {
                            if (ni.get(row, col) > cap) {
                                check = false;
                                break;
                            }
                        }
                    }
                    if (check) {
                        keep.add(s);
                    }
                } else {
                    if (ni.get(0, 0) <= cap) {
                        keep.add(s);
                    }
                }
            }
            Matrix newSpace = new Matrix(keep.size(), space.getNumCols());
            for (int i = 0; i < keep.size(); i++) {
                int rowIndex = keep.get(i);

                for (int col = 0; col < space.getNumCols(); col++) {
                    newSpace.set(i, col, space.get(rowIndex, col));
                }
            }
            space = newSpace.copy();
        }
        Matrix reverse = space.reverseRows();
        Matrix space_cp = space;
        return space.reverseRows();
    }
}