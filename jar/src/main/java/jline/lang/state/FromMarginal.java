package jline.lang.state;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.InputOutput;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.*;
import jline.lang.NodeParam;
import jline.lang.Mode;
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

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;
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

        if (sn.isstation.get(ind) == 1 && sn.nodetype.get(ind) != NodeType.Place) {
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

        // Synchronous call (REPLY signal): the node holds one server per job that has
        // left for its callee and is waiting for the reply. Those servers are not
        // derivable from the marginal n, so enumerate the held-server counts here and
        // build the rest of the state with the REMAINING servers -- with b servers held,
        // only S-b jobs can be in service, a configuration the plain enumeration never
        // produces. Recurse on a struct with the block cleared, then append its columns,
        // which are the last ones in the local-variable layout (see ReplyBlock).
        if (ReplyBlock.holds(sn, ind)) {
            List<Integer> rclasses = new ArrayList<Integer>();
            for (int r = 0; r < R; r++) {
                if (sn.replyblock.get(ind, r) > 0) {
                    rclasses.add(r);
                }
            }
            int Sist = (int) S.get(ist);
            NetworkStruct snb = sn.shallowCopy();
            snb.replyblock = sn.replyblock.copy();
            snb.nvars = sn.nvars.copy();
            snb.nservers = sn.nservers.copy();
            for (int r = 0; r < R; r++) {
                snb.replyblock.set(ind, r, 0);
                snb.nvars.set(ind, 2 * R + 1 + r, 0);
            }
            Matrix bspace = new Matrix(0, 0);
            for (int i = 0; i < rclasses.size(); i++) {
                Matrix bcol = new Matrix(Sist + 1, 1);
                for (int v = 0; v <= Sist; v++) {
                    bcol.set(v, 0, v);
                }
                bspace = Matrix.cartesian(bspace, bcol);
            }
            List<Matrix> subspaces = new ArrayList<Matrix>();
            List<Matrix> bkept = new ArrayList<Matrix>();
            int maxw = 0;
            for (int bi = 0; bi < bspace.getNumRows(); bi++) {
                Matrix b = Matrix.extractRows(bspace, bi, bi + 1, null);
                if (b.elementSum() > Sist) {
                    continue;
                }
                snb.nservers.set(ist, 0, Sist - b.elementSum());
                Matrix subspace = fromMarginal(snb, ind, n);
                if (subspace == null || subspace.getNumRows() == 0) {
                    continue;
                }
                subspaces.add(subspace);
                bkept.add(b);
                maxw = Math.max(maxw, subspace.getNumCols());
            }
            // Held servers push jobs into the buffer, so the sub-spaces have different
            // buffer widths. The buffer is RIGHT-aligned (empty slots pad the left), so
            // widen the narrow rows on the left before stacking them.
            Matrix stacked = new Matrix(0, 0);
            for (int bi = 0; bi < subspaces.size(); bi++) {
                Matrix subspace = subspaces.get(bi);
                if (subspace.getNumCols() < maxw) {
                    Matrix pad = new Matrix(subspace.getNumRows(), maxw - subspace.getNumCols());
                    pad.zero();
                    subspace = pad.concatCols(subspace);
                }
                Matrix rows = subspace.concatCols(bkept.get(bi).repmat(subspace.getNumRows(), 1));
                stacked = stacked.isEmpty() ? rows : Matrix.concatRows(stacked, rows, null);
            }
            if (stacked.getNumRows() == 0) {
                return new Matrix(0, maxw + rclasses.size());
            }
            return reverseRows(Maths.uniqueAndSort(stacked));
        }

        if (sn.isstateful.get(ind, 0) == 1 && sn.isstation.get(ind, 0) == 0) {
            for (int r = 0; r < R; r++) {
                Matrix init_r = spaceClosedSingle(1, n.get(r));
                state = Matrix.cartesian(state, init_r);
            }
            return Matrix.cartesian(space, state);
        }

        Matrix phases = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (sn.proc == null || sn.proc.get(sn.stations.get(ist)) == null || sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)) == null || sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
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

        SchedStrategy schedAtInd = sn.sched.get(((Station) sn.nodes.get(ind)));
        if (schedAtInd != null && schedAtInd != SchedStrategy.EXT && anyNGreaterThanClassCap) {
            return space;
        }

        // Generate local-state space
        switch (sn.nodetype.get(ind)) {
            case Place:
                // Place nodes hold tokens without service - track count per class
                // Similar to INF scheduling: just track jobs (tokens) per class
                for (int r = 0; r < R; r++) {
                    Matrix init_r = spaceClosedSingle(Math.max(1, phases.get(r)), n.get(r));
                    state = Matrix.cartesian(state, init_r);
                }
                space = Matrix.cartesian(space, state);
                break;
            case Queue:
            case Delay:
            case Source:
                switch (sn.sched.get(sn.stations.get(ist))) {
                    case EXT:
                        // source node case, treated as an infinite pool of jobs in buffer and a server for each class
                        for (int r = 0; r < R; r++) {
                            Matrix init_r = new Matrix(0, 0);
                            if (sn.markidx != null && ist < sn.markidx.getNumRows()
                                    && sn.markidx.get(ist, r) > 1) {
                                // marked (MMAP) non-carrier class: the modulating
                                // chain lives in the carrier's phase block; single
                                // always-zero column (mirrors sn.phasessz)
                                init_r = new Matrix(1, 1);
                                init_r.zero();
                            } else if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty() && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) { // if service disabled
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
                    case POLLING: {
                        // Un-ordered per-class buffers and a single server, as in SIRO,
                        // but NOT work-conserving over the station: the server may be
                        // idle while jobs wait, because it is walking towards a buffer
                        // (a switchover) or is parked. Both the empty-facility and the
                        // one-job-in-service configurations are therefore enumerated;
                        // the controller columns appended after the routing variables
                        // discard the combinations the discipline cannot occupy.
                        if (S.get(ist) != 1) {
                            throw new RuntimeException("Polling stations must have a single server.");
                        }
                        int sumKpoll = 0;
                        for (int r = 0; r < R; r++) {
                            sumKpoll += (int) phases.get(r);
                        }
                        // service facility empty
                        Matrix emptyRow = new Matrix(1, R + sumKpoll);
                        emptyRow.zero();
                        for (int r = 0; r < R; r++) {
                            emptyRow.set(0, r, n.get(r));
                        }
                        space = emptyRow;
                        for (int p = 0; p < R; p++) {
                            if (n.get(p) <= 0) {
                                continue;
                            }
                            Matrix srvp = new Matrix(0, 0);
                            for (int cls = 0; cls < R; cls++) {
                                Matrix init_r = spaceClosedSingle(phases.get(cls), (cls == p) ? 1 : 0);
                                srvp = Matrix.cartesian(srvp, init_r);
                            }
                            Matrix bufp = new Matrix(1, R);
                            for (int r = 0; r < R; r++) {
                                bufp.set(0, r, n.get(r));
                            }
                            bufp.set(0, p, n.get(p) - 1);
                            Matrix rowsp = bufp.repmat(srvp.getNumRows(), 1).concatCols(srvp);
                            space = Matrix.concatRows(space, rowsp, null);
                        }
                        break;
                    }
                    case SIRO:
                    case LEPT:
                    case SEPT:
                    case SRPT:
                    case SRPTPRIO:
                    case SETF:
                    case FSP:
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
                    case EDD:
                        // Retrial station: enumerate every (in-service, orbit) split,
                        // not just the work-conserving one. The server holds csrv in
                        // 0..min(n,S) jobs and the rest orbit in the buffer (class-id,
                        // right-aligned). The idle-server states are essential: an
                        // orbiting job re-enters service only through a RETRY event.
                        // Scoped to a single populated class (analyzer guard enforces).
                        boolean isRetrialStationFM = false;
                        if (sn.retrialProc != null) {
                            java.util.Map<jline.lang.JobClass, jline.util.matrix.MatrixCell> rpFM = sn.retrialProc.get(sn.stations.get(ist));
                            if (rpFM != null) {
                                for (jline.lang.JobClass jcFM : rpFM.keySet()) {
                                    jline.util.matrix.MatrixCell mcFM = rpFM.get(jcFM);
                                    if (mcFM != null && mcFM.size() > 0) { isRetrialStationFM = true; break; }
                                }
                            }
                        }
                        if (isRetrialStationFM) {
                            int rrFM = -1;
                            for (int r = 0; r < R; r++) if (n.get(0, r) > 0) { rrFM = r; break; }
                            if (rrFM < 0) {
                                space = new Matrix(1, (int) phases.elementSum());
                                space.zero();
                            } else {
                                int maxorbit = (int) n.get(0, rrFM);
                                int SistFM = (int) S.get(ist);
                                Matrix accFM = new Matrix(0, 0);
                                for (int csrv = 0; csrv <= Math.min(maxorbit, SistFM); csrv++) {
                                    int orbit = maxorbit - csrv;
                                    Matrix bufFM = new Matrix(1, maxorbit);
                                    bufFM.zero();
                                    for (int b = maxorbit - orbit; b < maxorbit; b++) bufFM.set(0, b, rrFM + 1);
                                    Matrix srvFM = new Matrix(0, 0);
                                    for (int cls = 0; cls < R; cls++) {
                                        Matrix init_c = spaceClosedSingle(phases.get(cls), (cls == rrFM) ? csrv : 0);
                                        srvFM = Matrix.cartesian(srvFM, init_c);
                                    }
                                    Matrix rowsFM = bufFM.repmat(srvFM.getNumRows(), 1).concatCols(srvFM);
                                    accFM = accFM.isEmpty() ? rowsFM : Matrix.concatRows(accFM, rowsFM, null);
                                }
                                space = accFM;
                            }
                            break;
                        }
                        Matrix vi = new Matrix(0, 0);
                        Matrix mi = new Matrix(0, 0);
                        double sizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) + Maths.factln(sn.cap.get(ist));
                        sizeEstimator = FastMath.round(sizeEstimator / FastMath.log(10));
                        if (sizeEstimator > 3) {
                            // Large state space warning - equivalent to MATLAB line_warning
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                InputOutput.line_warning(mfilename(new Object() {}), 
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
                            // Routing vars precede the node block in the nvars layout.
                            space = appendRouteVars(sn, ind, R, space);
                            // True BAS blocked marker: keep the state width consistent for an
                            // empty BAS station (blocked=0). Otherwise the n=0 state is one
                            // column short of the n>0 states, corrupting solvers (e.g. SSA)
                            // that build station states directly from fromMarginal.
                            // Gate on the dedicated sn.isbasblocking field (set for the
                            // blocking station under BOTH declaration forms) rather than the
                            // shared marker column, which a width-1 polling controller also
                            // sets. See BUG-83.
                            if (sn.isbasblocking != null && ind < sn.isbasblocking.length()
                                    && sn.isbasblocking.get(ind) == 1) {
                                Matrix zeroBlocked = new Matrix(1, 1);
                                zeroBlocked.set(0, 0, 0);
                                space = Matrix.cartesian(space, zeroBlocked);
                            }
                            // The empty-station early return bypasses the general appends
                            // at the end of the method, so the breakdown status must be
                            // enumerated here too. Emitting only status 0 would drop the
                            // empty-and-up state, and the station would then behave as if
                            // it could never empty, inflating its queue length by about
                            // one job -- the same failure mode appendRouteVars documents.
                            space = appendBreakdownStatus(sn, ind, space);
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
                                int startCol = miClone.getNumCols() - numCols;
                                for (int col = startCol; col < miClone.getNumCols(); col++) {
                                    mi.set(row, col - startCol, miClone.get(row, col));
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
                    case FCFSPR:
                    case FCFSPI:
                    case FCFSPRPRIO:
                    case FCFSPIPRIO:
                    case LCFSPI:
                    case LCFSPR:
                    case LCFSPIPRIO:
                    case LCFSPRPRIO:
                    case EDF:
                        Matrix vi_lpr = new Matrix(0, 0);
                        Matrix mi_lpr = new Matrix(0, 0);
                        // sum(n) - 1 due to Maths.factln including + 1
                        double lcfsprSizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) + Maths.factln(sn.cap.get(ist));
                        lcfsprSizeEstimator = FastMath.round(lcfsprSizeEstimator / FastMath.log(10));
                        if (lcfsprSizeEstimator > 3) {
                            // Large state space warning - equivalent to MATLAB line_warning
                            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                                InputOutput.line_warning(mfilename(new Object() {}), 
                                    "State space size is very large: 1e" + (int) lcfsprSizeEstimator + " states. " +
                                    "This may cause performance issues. Consider using a smaller model or force=true option.");
                            }
                        }

                        if (n.elementSum() == 0) {
                            Matrix newSpace = new Matrix(1, (int) (1 + phases.elementSum()));
                            newSpace.zero();
                            space = newSpace;
                            space = appendRouteVars(sn, ind, R, space);
                            space = appendBreakdownStatus(sn, ind, space);
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
                                int startCol = miClone.getNumCols() - numCols;
                                for (int col = startCol; col < miClone.getNumCols(); col++) {
                                    mi_lpr.set(row, col - startCol, miClone.get(row, col));
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
                    case OI:
                    case PAS: {
                        // Pass-and-swap / order-independent: the local state is the
                        // ordered list of 1-based class indices, left-aligned and
                        // right zero-padded to the station capacity (no server split).
                        if (Double.isInfinite(sn.cap.get(ist))) {
                            throw new RuntimeException("PAS stations require finite capacity for state-space generation.");
                        }
                        int Wpas = (int) sn.cap.get(ist);
                        int totalPas = (int) n.elementSum();
                        if (totalPas == 0) {
                            space = new Matrix(1, Wpas);
                            space.zero();
                        } else if (totalPas > Wpas) {
                            space = new Matrix(0, Wpas);
                        } else {
                            Matrix viPas = new Matrix(0, 0);
                            for (int r = 0; r < R; r++) {
                                if (n.get(0, r) > 0) {
                                    Matrix newVi = new Matrix(1, viPas.getNumCols() + (int) n.get(0, r));
                                    for (int i = 0; i < viPas.getNumCols(); i++) newVi.set(0, i, viPas.get(0, i));
                                    for (int i = viPas.getNumCols(); i < newVi.getNumCols(); i++) newVi.set(0, i, r + 1);
                                    viPas = newVi.copy();
                                }
                            }
                            Matrix miPas = Maths.uniquePerms(viPas);
                            space = new Matrix(miPas.getNumRows(), Wpas);
                            space.zero();
                            for (int row = 0; row < miPas.getNumRows(); row++)
                                for (int col = 0; col < miPas.getNumCols(); col++)
                                    space.set(row, col, miPas.get(row, col));
                        }
                        // Return directly: the ordered-list states are already
                        // unique and fixed-width (mirrors Python/MATLAB).
                        return space;
                    }
                    case SJF:
                    case LJF:
                        // in these policies the state space includes continuous
                        // random variables for the service times
                        throw new RuntimeException("The scheduling policy does not admit a discrete state space.");
                }

                space = appendRouteVars(sn, ind, R, space);
                // The polling controller trails the routing variables, matching the
                // nvars column order (modulation, routing, node block).
                space = Polling.space(sn, ind, space, phases);
                // True BAS blocked marker (nvars col 2R == 1): a completed job held at the
                // server awaiting room downstream. Meaningful only when the station holds
                // >=1 job; enumerate {0,1} then, else {0}.
                if (sn.isbasblocking != null && ind < sn.isbasblocking.length()
                        && sn.isbasblocking.get(ind) == 1) {
                    Matrix blockedCol;
                    if (n.elementSum() > 0) {
                        blockedCol = new Matrix(2, 1);
                        blockedCol.set(0, 0, 0);
                        blockedCol.set(1, 0, 1);
                    } else {
                        blockedCol = new Matrix(1, 1);
                        blockedCol.set(0, 0, 0);
                    }
                    space = Matrix.cartesian(space, blockedCol);
                }
                // Server breakdown status (nvars col 2R == 1): 0 = down, 1 = up.
                // Enumerated for EVERY marginal including the empty one, because the
                // failure clock runs whenever the server is up, idle or busy. The
                // column is exclusive with the BAS marker and the polling controller,
                // which Network.refreshLocalVars enforces, so it is always the
                // trailing column here.
                space = appendBreakdownStatus(sn, ind, space);
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

        if (uniqueRows.isEmpty()) {
            // No states for this (marginal, started) combination, e.g. an
            // infeasible started vector for a preempt-priority station (a
            // lower-priority job running while a higher-priority one waits).
            // Return an empty (0-row) space rather than indexing uniqueRows.get(0);
            // the caller concatenates per-(n,s) spaces, so an empty one is a no-op.
            return new Matrix(0, space.getNumCols());
        }
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
                    InputOutput.line_warning(mfilename(new Object() {
                    }), "State vector at station " + ist + " exceeds the class capacity. Some service classes are disabled.\n");
                } else {
                    InputOutput.line_warning(mfilename(new Object() {
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
                            } else if (sn.markidx != null && ist < sn.markidx.getNumRows()
                                    && sn.markidx.get(ist, r) > 1) {
                                init.set(0, 0, 0); // marked non-carrier: chain lives in the carrier block
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
                case LPS:
                case PSPRIO:
                case DPSPRIO:
                case GPSPRIO:
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
                case SETF:
                case FSP:
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
                case EDD:
                    double sizeEstimator = Maths.multinomialln(n.sub(1, s));
                    sizeEstimator = FastMath.round(sizeEstimator / FastMath.log(10));
                    if (sizeEstimator > 2) {
                        if (!optionsForce) {
                            InputOutput.line_error(mfilename(new Object() {
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

                case FCFSPR:
                case FCFSPI:
                case FCFSPRPRIO:
                case FCFSPIPRIO:
                case LCFSPI:
                case LCFSPR:
                case LCFSPIPRIO:
                case LCFSPRPRIO:
                case EDF:
                    double sizeEstimatorLPR = Maths.multinomialln(n.sub(1, s));
                    sizeEstimatorLPR = FastMath.round(sizeEstimatorLPR / FastMath.log(10));
                    if (sizeEstimatorLPR > 2) {
                        if (!optionsForce) {
                            System.err.format("State space size is very large: 1e%d states. Stopping execution. Set options = true," + "to bypass this control.\n", FastMath.round(sizeEstimatorLPR / FastMath.log(10)));
                        }
                    }
                    if (n.elementSum() == 0) {
                        // Even-width empty preempt buffer: one empty (class,phase)
                        // pair (width 2), see fromMarginalAndStarted.
                        space = new Matrix(1, (int) (2 + K.elementSum()));
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

    /**
     * Wrapper: the discipline branches below return early from several places, so the
     * synchronous-call (REPLY) counter columns are appended here, once, for every exit
     * path. An initial or user-supplied state has no call outstanding, so the counters
     * are zero -- but the columns must be present, otherwise the row is narrower than
     * the enumerated local space, matchrow fails, and Solver_ctmc silently skips its
     * unreachable-state pruning, leaving the enumerated-but-unreachable "counter set
     * while every job is here" states as a second absorbing class.
     *
     * @param sn           network structure
     * @param ind          node index
     * @param n            per-class marginal queue lengths
     * @param s            per-class jobs in service
     * @param optionsForce force flag
     * @return the local state space rows
     */
    public static Matrix fromMarginalAndStarted(NetworkStruct sn, int ind, Matrix n, Matrix s, Boolean optionsForce) {
        Matrix space = subFromMarginalAndStarted(sn, ind, n, s, optionsForce);
        if (ReplyBlock.holds(sn, ind) && space != null && space.getNumRows() > 0) {
            int width = 0;
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.replyblock.get(ind, r) > 0) {
                    width++;
                }
            }
            Matrix pad = new Matrix(space.getNumRows(), width);
            pad.zero();
            space = space.concatCols(pad);
        }
        return space;
    }

    private static Matrix subFromMarginalAndStarted(NetworkStruct sn, int ind, Matrix n, Matrix s, Boolean optionsForce) {
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
                        // firingproc is keyed by Mode, not by mode index: passing
                        // the int autoboxes to an Integer that never matches, so
                        // get returns null and the chained get(0) throws.
                        Mode modeObjK = ((Transition) sn.nodes.get(ind)).getModes().get(m);
                        K.set(0, m, ((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).firingproc.get(modeObjK).get(0).length());
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
                            } else if (sn.markidx != null && ist < sn.markidx.getNumRows()
                                    && sn.markidx.get(ist, r) > 1) {
                                init.set(0, 0, 0); // marked non-carrier: chain lives in the carrier block
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
                case LPS:
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

                case POLLING:
                case SIRO:
                case LEPT:
                case SEPT:
                case SRPT:
                case SRPTPRIO:
                case SETF:
                case FSP:
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
                case EDD:
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
                case FCFSPR:
                case FCFSPI:
                case FCFSPRPRIO:
                case FCFSPIPRIO:
                case LCFSPI:
                case LCFSPR:
                case LCFSPIPRIO:
                case LCFSPRPRIO:
                case EDF:
                    Matrix inbuf_lpr = new Matrix(0, 0);
                    double sizeEstimator_lpr = 0;
                    Matrix mi_lpr = new Matrix(0, 0);
                    Matrix mi_buf_lpr = new Matrix(0, 0);
                    Matrix mi_srv_lpr = new Matrix(0, 0);
                    if (n.elementSum() == 0) {
                        // Even-width empty preempt buffer: these policies track the
                        // buffer as (class,phase) PAIRS, so the buffer must be even
                        // (an odd/lone column is half a pair and makes the simulator
                        // drop the first preempted job). Use one empty pair (width 2),
                        // matching the native-Python convention.
                        space = new Matrix(1, (int) (2 + K.elementSum()));
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
                case OI:
                case PAS: {
                    // Pass-and-swap / order-independent: ordered list of 1-based
                    // class indices, left-aligned and right zero-padded to cap.
                    // The started counts s are immaterial (no server/buffer split).
                    if (Double.isInfinite(sn.cap.get(ist))) {
                        throw new RuntimeException("PAS stations require finite capacity for state-space generation.");
                    }
                    int Wpas = (int) sn.cap.get(ist);
                    int totalPas = (int) n.elementSum();
                    if (totalPas == 0) {
                        space = new Matrix(1, Wpas);
                        space.zero();
                    } else if (totalPas > Wpas) {
                        space = new Matrix(0, Wpas);
                    } else {
                        Matrix viPas = new Matrix(0, 0);
                        for (int r = 0; r < R; r++) {
                            if (n.get(0, r) > 0) {
                                Matrix newVi = new Matrix(1, viPas.getNumCols() + (int) n.get(0, r));
                                for (int i = 0; i < viPas.getNumCols(); i++) newVi.set(0, i, viPas.get(0, i));
                                for (int i = viPas.getNumCols(); i < newVi.getNumCols(); i++) newVi.set(0, i, r + 1);
                                viPas = newVi.copy();
                            }
                        }
                        Matrix miPas = Maths.uniquePerms(viPas);
                        space = new Matrix(miPas.getNumRows(), Wpas);
                        space.zero();
                        for (int row = 0; row < miPas.getNumRows(); row++)
                            for (int col = 0; col < miPas.getNumCols(); col++)
                                space.set(row, col, miPas.get(row, col));
                    }
                    return space;
                }
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
                    System.err.format("The scheduling policy does not admit a discrete state space");
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
            InputOutput.line_error(mfilename(new Object() {
            }), "fromMarginalAndStarted cannot be used on Petri net elements");
        } else if (sn.nodetype.get(ind) == NodeType.Place) {
            InputOutput.line_error(mfilename(new Object() {
            }), "fromMarginalAndStarted cannot be used on Petri net elements");
        } else if (sn.nodetype.get(ind) == NodeType.Join) {
            if (sn.isfjaugmented) {
                // FJ tag-augmented struct: the join state is the per-class
                // count vector of buffered jobs/siblings, deterministic given
                // the marginals (the started counts are immaterial, no service)
                space = new Matrix(1, sn.nclasses);
                for (int r = 0; r < sn.nclasses; r++) {
                    space.set(0, r, n.get(r));
                }
            } else {
                space = new Matrix(1, 1);
                space.zero();
            }
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
            line_error(mfilename(new Object(){}), "Bounds must either be both vectors or both scalars");
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
            int nvarsSum = (int) sn.nvars.sumRows(ind);
            for (int s = 0; s < space.getNumRows(); s++) {
                // toMarginal expects state_i to include trailing nvars columns
                // and strips them off; fromMarginalBounds builds buffer-only
                // rows, so pad with zeros for non-station stateful nodes to
                // match the expected layout.
                Matrix stateRow;
                if (sn.isstation.get(ind, 0) == 0 && nvarsSum > 0) {
                    Matrix pad = new Matrix(1, nvarsSum);
                    pad.zero();
                    stateRow = Matrix.concatColumns(space.getRow(s), pad, null);
                } else {
                    stateRow = space.getRow(s);
                }
                State.StateMarginalStatistics result =
                        ToMarginal.toMarginal(sn, ind, stateRow, null, null, null, null, null);
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

    /**
     * Appends the round-robin routing-variable columns for node ind to space.
     * Every exit of fromMarginal must produce rows carrying these columns,
     * including the empty-station early returns: an empty state emitted
     * without the pointer column misaligns in fromMarginalBounds and silently
     * drops the empty configuration from the state space, which inflates the
     * station QLen by about one job (the chain can then never empty the
     * station).
     *
     * @param sn    network structure
     * @param ind   node index
     * @param R     number of classes
     * @param space the state rows built so far
     * @return the rows with the routing-variable columns appended
     */
    /**
     * Appends the server-breakdown status column (nvars col 2R) for node ind.
     *
     * <p>The status is 0 = down, 1 = up and is enumerated for EVERY marginal,
     * including the empty station, because the failure clock runs whenever the
     * server is up, idle or busy. Exactly one feature owns the shared trailing
     * local-variable column per station, which Network.refreshLocalVars
     * enforces.</p>
     *
     * @param sn    network structure
     * @param ind   node index
     * @param space the state rows built so far
     * @return the rows with the status column appended, unchanged when the station
     *         is not subject to breakdowns
     */
    private static Matrix appendBreakdownStatus(NetworkStruct sn, int ind, Matrix space) {
        if (sn.hasbreakdown == null || ind >= sn.hasbreakdown.length() || sn.hasbreakdown.get(ind) != 1) {
            return space;
        }
        Matrix statusCol = new Matrix(2, 1);
        statusCol.set(0, 0, 0);
        statusCol.set(1, 0, 1);
        return Matrix.cartesian(space, statusCol);
    }

    /**
     * Row-reversal, so that states where jobs start in phase 1 come first (the
     * {@code space(end:-1:1,:)} convention of the MATLAB enumeration).
     *
     * @param m matrix to reverse
     * @return a copy of {@code m} with the row order reversed
     */
    private static Matrix reverseRows(Matrix m) {
        Matrix out = new Matrix(m.getNumRows(), m.getNumCols());
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                out.set(i, j, m.get(m.getNumRows() - 1 - i, j));
            }
        }
        return out;
    }

    private static Matrix appendRouteVars(NetworkStruct sn, int ind, int R, Matrix space) {
        for (int r = 0; r < R; r++) {
            RoutingStrategy routingStrategy = sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r));
            if (routingStrategy == RoutingStrategy.WRROBIN) {
                // The WRROBIN slot holds a POSITION (1..cycle length) in the
                // weighted-outlink cycle, not a node index: a destination may
                // repeat in the cycle, and both the pointer advance and the
                // routing function (sub_rr_wrr) interpret the slot as a
                // position. Enumerating outlink values here instead produced
                // states the position-based dynamics can never reach or leave,
                // degenerating the chain (matches MATLAB spaceLocalVars).
                Matrix wol = sn.nodeparam.get(sn.nodes.get(ind)).weightedOutlinks.get(sn.jobclasses.get(r));
                int cyc = (wol == null) ? 0 : (int) wol.length();
                if (cyc == 0) {
                    Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(r));
                    cyc = (outlinks == null) ? 0 : (int) outlinks.length();
                }
                if (cyc > 0) {
                    Matrix positions = new Matrix(cyc, 1);
                    for (int i = 0; i < cyc; i++) {
                        positions.set(i, 0, i + 1);
                    }
                    space = Matrix.cartesian(space, positions);
                }
            } else if (routingStrategy == RoutingStrategy.RROBIN) {
                Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(r));
                if (outlinks != null) {
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
        return space;
    }

}