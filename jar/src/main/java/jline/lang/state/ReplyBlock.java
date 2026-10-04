package jline.lang.state;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Layout of, and accessors for, the synchronous-call (REPLY signal) blocked-server
 * block that a node carries in its local-variable state.
 *
 * <p>A job of a class r with {@code sn.syncreply(r) >= 0} makes a synchronous call:
 * it leaves the station for the callee but KEEPS its server, which stays held until
 * the matching REPLY signal class comes back. The held servers are not derivable
 * from the marginal state (the job is at the callee, not here), so they are counted
 * per calling class in this block. LDES keys the same information by job id
 * (Solver_ssj.pendingReplyMap); a CTMC has no job identity, so it carries counts.
 *
 * <p>The block trails the modulation, routing and node blocks of {@code sn.nvars}
 * (columns 0..R-1, R..2R-1 and 2R), occupying columns 2R+1+r. Appending keeps every
 * existing nvars reader valid, and the columns stay zero-width for models without
 * reply signals, so no other model changes state width.
 *
 * <p>Port of MATLAB {@code State.replyBlockInfo} and {@code State.replyBlocked}.
 */
public class ReplyBlock implements Serializable {

    private static final long serialVersionUID = 1L;

    /** Layout of the reply block at one node. */
    public static class Info implements Serializable {
        private static final long serialVersionUID = 1L;
        /** Calling classes (0-based) that can hold a server at this node, in class order. */
        public List<Integer> classes = new ArrayList<Integer>();
        /** Offset of the block inside the node's local-variable vector. */
        public int off = 0;
        /** slot[r] is the 0-based index of class r inside the local-variable vector, -1 when absent. */
        public int[] slot;
        /** Number of columns in the block. */
        public int width = 0;
    }

    private ReplyBlock() {
    }

    /**
     * Layout of the reply block that node {@code ind} carries.
     *
     * @param sn  network structure
     * @param ind node index
     * @return the block layout; width 0 when the node holds no synchronous call
     */
    public static Info info(NetworkStruct sn, int ind) {
        int R = sn.nclasses;
        Info rinfo = new Info();
        rinfo.slot = new int[R];
        for (int r = 0; r < R; r++) {
            rinfo.slot[r] = -1;
        }
        if (sn.nvars == null || sn.nvars.getNumCols() < 3 * R + 1 || ind >= sn.nvars.getNumRows()) {
            return rinfo;
        }
        int off = 0;
        for (int c = 0; c <= 2 * R; c++) {
            off += (int) sn.nvars.get(ind, c);
        }
        rinfo.off = off;
        int pos = off;
        for (int r = 0; r < R; r++) {
            if (sn.nvars.get(ind, 2 * R + 1 + r) > 0) {
                rinfo.slot[r] = pos;
                pos++;
                rinfo.classes.add(r);
                rinfo.width++;
            }
        }
        return rinfo;
    }

    /**
     * True iff node {@code ind} holds servers for any synchronous call.
     *
     * @param sn  network structure
     * @param ind node index
     * @return true when the node carries a reply block
     */
    public static boolean holds(NetworkStruct sn, int ind) {
        if (sn.replyblock == null || sn.replyblock.isEmpty() || ind >= sn.replyblock.getNumRows()) {
            return false;
        }
        for (int r = 0; r < sn.replyblock.getNumCols(); r++) {
            if (sn.replyblock.get(ind, r) > 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * Per-class counts of servers held at node {@code ind} by jobs that made a
     * synchronous call and are waiting for their REPLY signal.
     *
     * @param sn       network structure
     * @param ind      node index
     * @param spaceVar local-variable part of the state, one row per state
     * @return (rows x nclasses) counts; all zero when the node carries no block
     */
    public static Matrix blocked(NetworkStruct sn, int ind, Matrix spaceVar) {
        int R = sn.nclasses;
        int nrows = (spaceVar == null || spaceVar.getNumRows() == 0) ? 1 : spaceVar.getNumRows();
        Matrix b = new Matrix(nrows, R);
        b.zero();
        if (spaceVar == null || spaceVar.isEmpty()) {
            return b;
        }
        Info rinfo = info(sn, ind);
        if (rinfo.width == 0) {
            return b;
        }
        for (int i = 0; i < rinfo.classes.size(); i++) {
            int r = rinfo.classes.get(i).intValue();
            if (rinfo.slot[r] < spaceVar.getNumCols()) {
                for (int row = 0; row < spaceVar.getNumRows(); row++) {
                    b.set(row, r, spaceVar.get(row, rinfo.slot[r]));
                }
            }
        }
        return b;
    }

    /**
     * Total number of servers held for pending replies, one entry per row of
     * {@code spaceVar}. Zero when the node carries no synchronous-call block, so
     * callers can subtract it from the server count unconditionally.
     *
     * @param sn       network structure
     * @param ind      node index
     * @param spaceVar local-variable part of the state, one row per state
     * @return (rows x 1) totals
     */
    public static Matrix blockedTotal(NetworkStruct sn, int ind, Matrix spaceVar) {
        Matrix b = blocked(sn, ind, spaceVar);
        Matrix nb = new Matrix(b.getNumRows(), 1);
        nb.zero();
        for (int row = 0; row < b.getNumRows(); row++) {
            double s = 0;
            for (int r = 0; r < b.getNumCols(); r++) {
                s += b.get(row, r);
            }
            nb.set(row, 0, s);
        }
        return nb;
    }
}
