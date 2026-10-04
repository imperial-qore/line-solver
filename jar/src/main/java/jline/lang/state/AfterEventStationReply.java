package jline.lang.state;

import java.io.Serializable;
import java.util.Map;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Arrival of a REPLY signal class at the FCFS station that is holding a server for
 * the matching synchronous call. The reply completes the call:
 *
 * <ol>
 *   <li>one held server of the calling class is released (the reply block counter is
 *       decremented);</li>
 *   <li>the reply itself takes that server, never a waiting job, and is served here
 *       (typically Immediate) before being routed on.</li>
 * </ol>
 *
 * <p>Unlike a NEGATIVE or CATASTROPHE signal the reply is NOT annihilated: it is a job
 * that carries the call result onward. This mirrors LDES, where a REPLY arrival frees
 * the blocked server and then continues routing.
 *
 * <p>The event is passive: the rate is set by the active departure at the callee.
 *
 * <p>Port of MATLAB {@code State.afterEventStationReply}.
 */
public class AfterEventStationReply implements Serializable {

    private static final long serialVersionUID = 1L;

    private AfterEventStationReply() {
    }

    /**
     * Applies a REPLY arrival at node {@code ind}.
     *
     * @param sn       network structure
     * @param ind      node index of the holding station
     * @param ist      station index of the holding station
     * @param jobClass 0-based index of the arriving REPLY class
     * @param K        per-class number of service phases
     * @param Ks       per-class phase offset into the server block
     * @param S        per-station number of servers
     * @param pie      per-station, per-class service entry-phase distribution
     * @param spaceBuf buffer part of the state
     * @param spaceSrv server part of the state
     * @param spaceVar local-variable part of the state
     * @return successor states, rates (-1, passive) and entry-phase probabilities
     */
    public static Ret.EventResult apply(NetworkStruct sn, int ind, int ist, int jobClass,
                                        Matrix K, Matrix Ks, Matrix S,
                                        Map<Station, Map<JobClass, Matrix>> pie,
                                        Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar) {
        ReplyBlock.Info rinfo = ReplyBlock.info(sn, ind);

        // The calling class this reply releases: the class whose expected reply is
        // jobClass and which holds a block here. sn.syncreply is stored 0-based.
        int callclass = -1;
        for (int i = 0; i < rinfo.classes.size(); i++) {
            int r = rinfo.classes.get(i).intValue();
            if ((int) sn.syncreply.get(r, 0) == jobClass) {
                callclass = r;
                break;
            }
        }

        Matrix outspace = new Matrix(0, 0);
        Matrix outprob = new Matrix(0, 0);
        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
            Matrix buf = Matrix.extractRows(spaceBuf, row, row + 1, null);
            Matrix srv = Matrix.extractRows(spaceSrv, row, row + 1, null);
            Matrix var = Matrix.extractRows(spaceVar, row, row + 1, null);

            if (callclass >= 0 && rinfo.slot[callclass] >= 0
                    && var.get(0, rinfo.slot[callclass]) > 0) {
                // PASS-THROUGH: the released server is taken by the REPLY itself, never
                // by a waiting job. The reply is the released work of a call this station
                // already paid for, so it does not queue behind the residents (queueing
                // it gave QLen 0.125 and residence 0.0996 against 0 in LDES, and stole
                // capacity, costing 5% of throughput). Its service is typically
                // Immediate, so the server is handed back at once and the ordinary FCFS
                // departure path then promotes the head of line -- which also keeps
                // sum(srv) <= S, unlike admitting the reply on top of a promoted job.
                var.set(0, rinfo.slot[callclass], var.get(0, rinfo.slot[callclass]) - 1);
            }

            // The reply job joins the station: into a free server (enumerating its entry
            // phase) or, if all remaining servers are busy, at the tail of the buffer.
            double nb = ReplyBlock.blockedTotal(sn, ind, var).get(0, 0);
            double Seff = S.get(ist) - nb;
            if (srv.elementSum() < Seff) {
                Matrix pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
                for (int kentry = 0; kentry < (int) K.get(jobClass); kentry++) {
                    if (pentry.get(kentry) <= 0) {
                        continue;
                    }
                    Matrix srvK = srv.copy();
                    int col = (int) Ks.get(jobClass) + kentry;
                    srvK.set(0, col, srvK.get(0, col) + 1);
                    Matrix rowOut = buf.concatCols(srvK).concatCols(var);
                    outspace = outspace.isEmpty() ? rowOut : Matrix.concatRows(outspace, rowOut, null);
                    Matrix p = new Matrix(1, 1);
                    p.set(0, 0, pentry.get(kentry));
                    outprob = outprob.isEmpty() ? p : Matrix.concatRows(outprob, p, null);
                }
                continue;
            }
            // All available servers busy: queue at the tail. The FCFS buffer is
            // right-aligned with the newest job leftmost, so the tail is the last empty
            // slot; grow the buffer by one column when it is full, as the ordinary
            // arrival branch of AfterEventStation does.
            int emptypos = -1;
            for (int c = buf.getNumCols() - 1; c >= 0; c--) {
                if (buf.get(0, c) == 0) {
                    emptypos = c;
                    break;
                }
            }
            if (emptypos < 0) {
                Matrix grown = new Matrix(1, buf.getNumCols() + 1);
                grown.zero();
                for (int c = 0; c < buf.getNumCols(); c++) {
                    grown.set(0, c + 1, buf.get(0, c));
                }
                buf = grown;
                emptypos = 0;
            }
            buf.set(0, emptypos, jobClass + 1);
            Matrix rowOut = buf.concatCols(srv).concatCols(var);
            outspace = outspace.isEmpty() ? rowOut : Matrix.concatRows(outspace, rowOut, null);
            Matrix p = new Matrix(1, 1);
            p.set(0, 0, 1.0);
            outprob = outprob.isEmpty() ? p : Matrix.concatRows(outprob, p, null);
        }
        Matrix outrate = new Matrix(outspace.getNumRows(), 1);
        outrate.fill(-1.0); // passive action
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}
