/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.NodeParam;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.processes.Distribution;
import jline.lang.processes.Markovian;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static jline.api.mam.Map_pie.map_pie;

/**
 * State-space helpers for polling stations. Port of the MATLAB
 * {@code State.polling*} functions.
 *
 * <p>The polling controller is stored in the trailing local-variable block of
 * the station state, whose width is {@code sn.nvars(ind, 2*R)}. The block holds,
 * in order, the columns {@code [pos, swk, ctr]}; each is materialized only when
 * the discipline actually needs it, so that a polling station never carries
 * state that its dynamics cannot distinguish.</p>
 *
 * <ul>
 * <li>{@code pos} index of the buffer the server is currently at (serving) or
 * heading to (switching). Materialized only when at least one switchover is
 * non-immediate: while a job is in service pos always equals the class of that
 * job, and while the station is empty and every switchover is immediate the
 * server position is unobservable.</li>
 * <li>{@code swk} 0 when the server sits at pos, otherwise the phase of the
 * switchover PH into buffer pos. Materialized only when some switchover is
 * non-immediate.</li>
 * <li>{@code ctr} the visit budget, materialized for every discipline except
 * EXHAUSTIVE. GATED: jobs of class pos admitted at the polling instant that have
 * not completed yet (the job in service counts as one of them), so the visit ends
 * when ctr reaches 0. KLIMITED: services still permitted in this visit, the one
 * in progress included. DECREMENTING: the target class-pos population, i.e. one
 * below the level found at the polling instant (semi-exhaustive).</li>
 * </ul>
 *
 * <p>The three tangible controller configurations are SERVING(p) (swk=0, one
 * class-p job in the service facility), SWITCHING(p) (swk&gt;0, service facility
 * empty) and PARKED (swk=0, facility empty, station empty; reachable only when
 * every switchover is immediate, in which case a server that completes a full lap
 * without finding work would otherwise cycle in zero time forever).</p>
 *
 * <p>Immediate switchovers are NOT represented as states: taking their rate
 * literally would put a ~1e8 rate in the generator (stiff, and a spurious state
 * per buffer). They are instead folded into the enclosing transition by
 * {@link #next}, which walks the cyclic order until a tangible state.</p>
 */
public class Polling implements Serializable {

    private static final long serialVersionUID = 1L;

    /** Landing mode: start a visit at the resolved buffer. */
    public static final int MODE_VISIT = 1;
    /** Landing mode: enter the switchover into the resolved buffer. */
    public static final int MODE_SWITCH = 2;
    /** Landing mode: park (no work anywhere, every switchover immediate). */
    public static final int MODE_PARK = 0;

    /**
     * Derived description of the polling controller at a station.
     */
    public static class Info implements Serializable {
        private static final long serialVersionUID = 1L;
        /** Polling discipline, identical across the class buffers. */
        public PollingType ptype;
        /** K parameter of KLIMITED polling. */
        public int pk;
        /** Buffers the server visits (a class disabled at the station is dropped). */
        public boolean[] polled;
        /** True when the leg entering buffer q takes a strictly positive time. */
        public boolean[] hasSw;
        /** Number of phases of the switchover entering buffer q. */
        public int[] Ksw;
        /** D0 of the switchover entering buffer q. */
        public Matrix[] swD0;
        /** D1 of the switchover entering buffer q. */
        public Matrix[] swD1;
        /** Entry probability vector of the switchover entering buffer q. */
        public Matrix[] swPie;
        /** Whether the pos/swk/ctr columns are materialized. */
        public boolean wpos, wswk, wctr;
        /** Total width of the polling block. */
        public int width;
    }

    /**
     * Derived description of the polling controller at the station of node ind,
     * or null when the node is not a polling station. Memoized on the node
     * parameters: this is read once per state per synchronization while the
     * generator is built, and rebuilding it (map_pie and all) per call is orders
     * of magnitude slower than the analysis itself.
     *
     * @param sn network structure
     * @param ind node index
     * @return the controller description, or null
     */
    public static Info info(NetworkStruct sn, int ind) {
        if (sn.isstation.get(ind, 0) == 0) {
            return null;
        }
        int ist = (int) sn.nodeToStation.get(0, ind);
        if (ist < 0) {
            return null;
        }
        Station station = sn.stations.get(ist);
        if (sn.sched.get(station) != SchedStrategy.POLLING) {
            return null;
        }
        NodeParam np = (sn.nodeparam == null) ? null : sn.nodeparam.get(sn.nodes.get(ind));
        QueueNodeParam qnp = (np instanceof QueueNodeParam) ? (QueueNodeParam) np : null;
        if (qnp != null && qnp.pollinfo != null) {
            return qnp.pollinfo;
        }

        int R = sn.nclasses;
        Info pinfo = new Info();
        pinfo.ptype = PollingType.EXHAUSTIVE;
        pinfo.pk = 1;
        pinfo.polled = new boolean[R];
        pinfo.hasSw = new boolean[R];
        pinfo.Ksw = new int[R];
        pinfo.swD0 = new Matrix[R];
        pinfo.swD1 = new Matrix[R];
        pinfo.swPie = new Matrix[R];

        // Buffers the server actually visits. A class disabled at this station can
        // never hold a job, and State.afterEvent short-circuits every event carrying
        // it (K(class)==0), so it cannot be given a switchover to walk through: such
        // a buffer is dropped from the cyclic order rather than polled forever.
        List<Integer> polledList = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            MatrixCell proc_ir = (sn.proc == null || sn.proc.get(station) == null)
                    ? null : sn.proc.get(station).get(sn.jobclasses.get(r));
            boolean enabled = proc_ir != null && proc_ir.size() > 0 && !proc_ir.get(0).hasNaN();
            pinfo.polled[r] = enabled;
            if (enabled) {
                polledList.add(Integer.valueOf(r));
            }
        }
        if (polledList.isEmpty()) {
            throw new RuntimeException("Polling station " + ist + " has no class with an enabled service.");
        }

        if (qnp != null) {
            if (qnp.pollingType != null) {
                pinfo.ptype = qnp.pollingType;
            }
            if (qnp.pollingPar != null) {
                pinfo.pk = qnp.pollingPar.intValue();
            }
        }

        // Switchover indexing. Queue.setSwitchover(class_i, distrib) is the time to
        // switch FROM buffer i to the next one, i.e. Takagi's r_i, and that is the
        // convention of the exact formulas in api/polling. The controller indexes its
        // states by DESTINATION, so the lookup is shifted onto the previous polled
        // buffer here, once, rather than at every use.
        for (int jj = 0; jj < polledList.size(); jj++) {
            int q = polledList.get(jj).intValue();
            int src = polledList.get(Math.floorMod(jj - 1, polledList.size())).intValue();
            Distribution sw = (qnp == null || qnp.switchoverTime == null)
                    ? null : qnp.switchoverTime.get(sn.jobclasses.get(src));
            if (sw == null || sw.isImmediate()) {
                continue; // zero-time switchover: folded, never a state
            }
            if (!(sw instanceof Markovian)) {
                throw new RuntimeException("Polling switchover distributions must be phase-type or Markovian.");
            }
            MatrixCell proc = ((Markovian) sw).getProcess();
            pinfo.hasSw[q] = true;
            pinfo.Ksw[q] = proc.get(0).getNumRows();
            pinfo.swD0[q] = proc.get(0);
            pinfo.swD1[q] = proc.get(1);
            pinfo.swPie[q] = map_pie(proc.get(0), proc.get(1));
        }

        // Column widths. pos and swk exist only to encode SWITCHING(p); with no
        // non-immediate switchover the server is either serving (pos = class in
        // service) or parked (pos unobservable), so neither column carries
        // information. ctr exists for every discipline that bounds a visit.
        boolean anySw = false;
        for (int r = 0; r < R; r++) {
            if (pinfo.hasSw[r]) { anySw = true; break; }
        }
        pinfo.wpos = anySw;
        pinfo.wswk = anySw;
        pinfo.wctr = pinfo.ptype != PollingType.EXHAUSTIVE;
        pinfo.width = (pinfo.wpos ? 1 : 0) + (pinfo.wswk ? 1 : 0) + (pinfo.wctr ? 1 : 0);

        if (qnp != null) {
            qnp.pollinfo = pinfo;
        }
        return pinfo;
    }

    /**
     * Initial value of the ctr column of a visit that starts at a buffer holding
     * nbufq waiting jobs. nbufq is the class-q population at the polling instant:
     * the service facility is empty when a visit starts, so it is the whole
     * class-q population at the station.
     *
     * @param pinfo controller description
     * @param nbufq waiting jobs of the visited class at the polling instant
     * @return the initial visit budget
     */
    public static int budget(Info pinfo, int nbufq) {
        switch (pinfo.ptype) {
            case EXHAUSTIVE:
                return 0; // unused: the visit ends when the buffer drains
            case GATED:
                return nbufq; // serve exactly the jobs found at the polling instant
            case KLIMITED:
                return pinfo.pk; // serve at most K, fewer if the buffer drains first
            case DECREMENTING:
                return nbufq - 1; // serve until the population drops one below the level found
            default:
                throw new RuntimeException("Unsupported polling type: " + pinfo.ptype);
        }
    }

    /**
     * Resolve the tangible controller state a polling server reaches once it stops
     * serving buffer pos. Port of MATLAB State.pollingNext.
     *
     * <p>The walk is what makes an Immediate switchover a zero-time leg rather than
     * a state: the server passes straight through such a buffer when it has no
     * work, and only stops at a buffer that either has work or costs time to reach.
     * This is the classical cyclic-polling discipline, in which the server visits
     * the buffers in strict cyclic order and pays the switchover of every buffer it
     * moves to, whether or not that buffer turns out to have work.</p>
     *
     * @param pinfo   controller description
     * @param pos     buffer the server is leaving, or has arrived at
     * @param nbuf    per-class waiting counts (the service facility is empty)
     * @param R       number of classes
     * @param arrived when true the server has just ARRIVED at pos and pos itself is
     *                examined first, without charging its switchover a second time;
     *                when false the server is LEAVING pos and the walk starts at pos+1
     * @return {q, mode, budget}; mode is MODE_VISIT, MODE_SWITCH or MODE_PARK
     */
    public static int[] next(Info pinfo, int pos, int[] nbuf, int R, boolean arrived) {
        if (arrived && pinfo.polled[pos] && nbuf[pos] > 0) {
            // the switchover into pos has already been paid, so a visit starts here
            return new int[]{pos, MODE_VISIT, budget(pinfo, nbuf[pos])};
        }
        int p = pos;
        for (int step = 1; step <= R; step++) { // a full lap, so the last buffer examined is pos
            p = (p + 1) % R;
            if (!pinfo.polled[p]) {
                continue;
            }
            if (pinfo.hasSw[p]) {
                // moving to p costs a strictly positive switchover: the server dwells
                // in it regardless of whether p has work to offer on arrival
                return new int[]{p, MODE_SWITCH, 0};
            }
            if (nbuf[p] > 0) {
                return new int[]{p, MODE_VISIT, budget(pinfo, nbuf[p])};
            }
        }
        return new int[]{pos, MODE_PARK, 0};
    }

    /**
     * Every controller configuration compatible with a service facility holding a
     * class-srvclass job (-1 when empty) and buffers holding nbuf. Returns one
     * {pos, swk, ctr} triple per configuration, in full regardless of which columns
     * are materialized (project with {@link #project}), and an empty list when the
     * combination is unoccupiable.
     *
     * <p>PARKED is pruned even when no column is materialized: leaving an idle
     * server with a non-empty station in the state space would make the generator
     * reducible.</p>
     *
     * @param pinfo    controller description
     * @param srvclass class in the service facility, -1 when empty
     * @param nbuf     per-class waiting counts
     * @param R        number of classes
     * @return the admissible {pos, swk, ctr} triples
     */
    public static List<int[]> blocks(Info pinfo, int srvclass, int[] nbuf, int R) {
        List<int[]> trips = new ArrayList<int[]>();
        if (srvclass >= 0) {
            if (!pinfo.polled[srvclass]) {
                return trips; // a buffer outside the cyclic order can hold no job
            }
            // SERVING(p): the server stands at the buffer of the job it is serving, so
            // pos is pinned to srvclass. The visit budget is free within the bounds its
            // discipline can have left it in.
            int lo, hi;
            switch (pinfo.ptype) {
                case EXHAUSTIVE:
                    lo = 0; hi = 0;
                    break;
                case GATED:
                    // ctr counts the admitted jobs not yet completed, the one in service
                    // included, so ctr >= 1; the ctr-1 uncompleted others all wait in the buffer
                    lo = 1; hi = nbuf[srvclass] + 1;
                    break;
                case KLIMITED:
                    lo = 1; hi = pinfo.pk;
                    break;
                case DECREMENTING:
                    // ctr is the population the visit drives the class down to; it started
                    // one below the level found and the visit is still running
                    lo = 0; hi = nbuf[srvclass];
                    break;
                default:
                    throw new RuntimeException("Unsupported polling type: " + pinfo.ptype);
            }
            for (int ctr = lo; ctr <= hi; ctr++) {
                trips.add(new int[]{srvclass, 0, ctr});
            }
        } else {
            // SWITCHING(q): the facility must be empty and q must cost time to reach.
            // Jobs may wait meanwhile: this is exactly what makes a polling station
            // non-work-conserving.
            for (int q = 0; q < R; q++) {
                if (!pinfo.hasSw[q]) {
                    continue;
                }
                for (int swk = 1; swk <= pinfo.Ksw[q]; swk++) {
                    trips.add(new int[]{q, swk, 0});
                }
            }
            boolean anySw = false;
            for (int r = 0; r < R; r++) {
                if (pinfo.hasSw[r]) { anySw = true; break; }
            }
            int tot = 0;
            for (int r = 0; r < R; r++) {
                tot += nbuf[r];
            }
            if (!anySw && tot == 0) {
                int first = 0;
                for (int r = 0; r < R; r++) {
                    if (pinfo.polled[r]) { first = r; break; }
                }
                trips.add(new int[]{first, 0, 0}); // parked, pos canonical
            }
        }
        return trips;
    }

    /**
     * Result of a landing: the successor rows and the probability of each.
     */
    public static class Landing implements Serializable {
        private static final long serialVersionUID = 1L;
        /** Successor [buffer, server, local-var] rows. */
        public List<Matrix> rows = new ArrayList<Matrix>();
        /** Probability of each successor row. */
        public List<Double> probs = new ArrayList<Double>();
    }

    /**
     * Materialize the state rows a polling server lands in after {@link #next} has
     * resolved (q, mode, budget). The three space rows describe the station at the
     * instant the decision is taken, i.e. with the completed job (if any) already
     * removed from the service facility.
     *
     * <p>The probabilities split a landing across the entry phases of a phase-type:
     * which phase a service or a switchover starts in is a random choice, so one
     * decision yields one row per entry phase, weighted by the corresponding entry
     * probability. Callers fold them into the transition rate rather than into
     * outprob, since the branching happens at the instant the active event fires.</p>
     *
     * @param pinfo    controller description
     * @param q        resolved buffer
     * @param mode     MODE_VISIT, MODE_SWITCH or MODE_PARK
     * @param budget   visit budget for MODE_VISIT
     * @param spaceBuf single buffer row
     * @param spaceSrv single server row
     * @param spaceVar single local-variable row
     * @param K        per-class phase counts
     * @param Ks       per-class phase offsets
     * @param pieist   per-class entry probability vectors of the station
     * @param R        number of classes
     * @return the successor rows and their probabilities
     */
    public static Landing land(Info pinfo, int q, int mode, int budget,
                               Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar,
                               Matrix K, Matrix Ks, Map<JobClass, Matrix> pieist,
                               List<JobClass> jobclasses, int R) {
        Landing res = new Landing();
        if (mode == MODE_VISIT) {
            // start or continue a visit at q: pull a waiting class-q job into service
            Matrix buf = spaceBuf.copy();
            buf.set(0, q, buf.get(0, q) - 1);
            Matrix pentry = pieist.get(jobclasses.get(q));
            for (int kentry = 0; kentry < (int) K.get(q); kentry++) {
                double p = pentry.get(kentry);
                if (p <= 0) {
                    continue;
                }
                Matrix srv = spaceSrv.copy();
                int col = (int) Ks.get(q) + kentry;
                srv.set(0, col, srv.get(0, col) + 1);
                Matrix var = set(pinfo, spaceVar, q, 0, budget);
                res.rows.add(catRow(buf, srv, var));
                res.probs.add(Double.valueOf(p));
            }
        } else if (mode == MODE_SWITCH) {
            // enter the switchover into q: the facility stays empty while walking
            Matrix swpie = pinfo.swPie[q];
            for (int kentry = 0; kentry < pinfo.Ksw[q]; kentry++) {
                double p = swpie.get(kentry);
                if (p <= 0) {
                    continue;
                }
                Matrix var = set(pinfo, spaceVar, q, kentry + 1, 0);
                res.rows.add(catRow(spaceBuf, spaceSrv, var));
                res.probs.add(Double.valueOf(p));
            }
        } else {
            // park: held until the next arrival, see next()
            Matrix var = set(pinfo, spaceVar, q, 0, 0);
            res.rows.add(catRow(spaceBuf, spaceSrv, var));
            res.probs.add(Double.valueOf(1.0));
        }
        return res;
    }

    private static Matrix catRow(Matrix buf, Matrix srv, Matrix var) {
        Matrix out = buf;
        if (srv != null && srv.getNumCols() > 0) {
            out = Matrix.concatColumns(out, srv, null);
        }
        if (var != null && var.getNumCols() > 0) {
            out = Matrix.concatColumns(out, var, null);
        }
        return out;
    }

    /**
     * Append the polling controller columns to the rows of space, which must hold
     * the [buffer, server, routing-variable] layout of a polling station. Rows are
     * expanded into one row per controller configuration the discipline can occupy,
     * and rows for which no configuration exists are dropped.
     *
     * <p>Enumerating the controller per row rather than as a blind cartesian
     * product is what keeps the state space tight and the chain irreducible: pos is
     * pinned to the class in service, a switchover excludes a busy service facility,
     * and a park excludes a non-empty station. A cartesian product would instead
     * admit states such as "serving buffer 1 while a class-2 job occupies the
     * server", which no transition can reach or leave consistently with the
     * marginals.</p>
     *
     * @param sn    network structure
     * @param ind   node index
     * @param space the [buffer, server, routing-var] rows
     * @param K     per-class phase counts
     * @return the rows with the controller appended
     */
    public static Matrix space(NetworkStruct sn, int ind, Matrix space, Matrix K) {
        Info pinfo = info(sn, ind);
        if (pinfo == null || space == null || space.isEmpty() || pinfo.width == 0) {
            return space;
        }
        int R = sn.nclasses;
        int sumK = 0;
        int[] Ks = new int[R];
        for (int r = 0; r < R; r++) {
            Ks[r] = sumK;
            sumK += (int) K.get(r);
        }
        List<Matrix> out = new ArrayList<Matrix>();
        for (int row = 0; row < space.getNumRows(); row++) {
            int[] nbuf = new int[R];
            for (int r = 0; r < R; r++) {
                nbuf[r] = (int) Math.round(space.get(row, r));
            }
            // class occupying the single service facility, -1 when it is empty
            int srvclass = -1;
            for (int r = 0; r < R; r++) {
                double tot = 0;
                for (int k = 0; k < (int) K.get(r); k++) {
                    tot += space.get(row, R + Ks[r] + k);
                }
                if (tot > 0) {
                    srvclass = r;
                    break;
                }
            }
            List<int[]> trips = blocks(pinfo, srvclass, nbuf, R);
            for (int b = 0; b < trips.size(); b++) {
                Matrix rowMat = Matrix.extractRows(space, row, row + 1, null);
                out.add(rowMat.concatCols(project(pinfo, trips.get(b))));
            }
        }
        if (out.isEmpty()) {
            return new Matrix(0, space.getNumCols() + pinfo.width);
        }
        Matrix res = out.get(0);
        for (int i = 1; i < out.size(); i++) {
            res = Matrix.concatRows(res, out.get(i), null);
        }
        return res;
    }

    /**
     * Project a full {pos, swk, ctr} triple onto the materialized columns. The
     * elided columns are reconstructible from the rest of the state (see
     * {@link #get}), so keeping them would split each state into copies no
     * observation can tell apart.
     *
     * @param pinfo controller description
     * @param trip  the full triple
     * @return a 1-by-width row, empty when nothing is materialized
     */
    public static Matrix project(Info pinfo, int[] trip) {
        Matrix row = new Matrix(1, pinfo.width);
        int c = 0;
        if (pinfo.wpos) { row.set(0, c, trip[0]); c++; }
        if (pinfo.wswk) { row.set(0, c, trip[1]); c++; }
        if (pinfo.wctr) { row.set(0, c, trip[2]); c++; }
        return row;
    }

    /**
     * Read the polling controller out of the local-variable columns of a single
     * state row. The polling block occupies the trailing width columns of
     * spaceVar, so it is located by counting back from the end.
     *
     * <p>The elided columns are reconstructed here, so callers always see a
     * complete controller. pos: with no materialized column every switchover is
     * immediate, so the server is either serving (and then it stands at the buffer
     * of the job in service) or parked (and then its position is unobservable, so
     * the first polled buffer is canonical). swk: no switchover takes time, so the
     * server is never inside one. ctr: the discipline is EXHAUSTIVE, which bounds a
     * visit by the buffer draining rather than by a budget.</p>
     *
     * @param pinfo    controller description
     * @param spaceVar the local-variable columns of one state row
     * @param srvclass class in the service facility, -1 when empty
     * @return {pos, swk, ctr}
     */
    public static int[] get(Info pinfo, Matrix spaceVar, int srvclass) {
        int pos;
        int swk = 0;
        int ctr = 0;
        int base = ((spaceVar == null) ? 0 : spaceVar.getNumCols()) - pinfo.width;
        int c = base;
        if (pinfo.wpos) {
            pos = (int) Math.round(spaceVar.get(0, c));
            c++;
        } else if (srvclass >= 0) {
            pos = srvclass;
        } else {
            pos = 0;
            for (int r = 0; r < pinfo.polled.length; r++) {
                if (pinfo.polled[r]) { pos = r; break; }
            }
        }
        if (pinfo.wswk) {
            swk = (int) Math.round(spaceVar.get(0, c));
            c++;
        }
        if (pinfo.wctr) {
            ctr = (int) Math.round(spaceVar.get(0, c));
        }
        return new int[]{pos, swk, ctr};
    }

    /**
     * Write the polling controller into the local-variable columns of a state row.
     * Columns that are not materialized are dropped: they are reconstructible from
     * the rest of the state (see {@link #get}).
     *
     * @param pinfo    controller description
     * @param spaceVar the local-variable columns of one state row (copied)
     * @param pos      server position
     * @param swk      switchover phase, 0 when the server sits at pos
     * @param ctr      visit budget
     * @return the updated local-variable row
     */
    public static Matrix set(Info pinfo, Matrix spaceVar, int pos, int swk, int ctr) {
        Matrix out = spaceVar.copy();
        int base = out.getNumCols() - pinfo.width;
        int c = base;
        if (pinfo.wpos) { out.set(0, c, pos); c++; }
        if (pinfo.wswk) { out.set(0, c, swk); c++; }
        if (pinfo.wctr) { out.set(0, c, ctr); }
        return out;
    }
}
