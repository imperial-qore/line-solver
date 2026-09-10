/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.processes.DiscreteDistribution;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Passive arrival of a G-network signal class at a station. Port of the MATLAB
 * {@code State.afterEventStationSignal} and {@code State.signalBatchPMF}.
 *
 * <p>A signal never joins the station: it acts on the jobs already there and is
 * annihilated. CATASTROPHE empties the station, ignoring the removal count
 * distribution (a catastrophe removes all jobs by definition). NEGATIVE removes
 * a batch of jobs:</p>
 *
 * <ul>
 *   <li>Victims. When the signal declares a target class ({@code sn.signaltarget >= 0},
 *       set via forJobClass) only that class is eligible; otherwise every non-signal
 *       class is eligible. The untargeted case is the classic Gelenbe negative
 *       customer and agrees with SolverMAM and SolverLDES, which are both
 *       class-agnostic.</li>
 *   <li>Count. {@code sn.signalremdist} gives the batch-size pmf. An oversized batch
 *       empties the eligible jobs rather than driving the queue negative, so the pmf
 *       tail P(B >= n) lumps onto "remove all n" (the min(B,n) clipping used by LDES
 *       and the tail term used by MAM).</li>
 *   <li>Policy. {@code sn.signalrempolicy} selects the victim: FCFS removes the oldest
 *       waiting job, LCFS the newest, RANDOM draws uniformly over waiting and
 *       in-service jobs. FCFS/LCFS only touch an in-service job when no job is
 *       waiting (two-tier, as in LDES).</li>
 * </ul>
 *
 * <p>The station state carries no arrival-time order for in-service jobs, so when a
 * policy has to reach into the servers the victim is drawn uniformly across the
 * occupied phases (exact whenever at most one job of the class is in service, which
 * covers every single-server station).</p>
 */
public class SignalRemoval implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * True if the signal empties a station on arrival. The signaltype is consulted
     * alongside the iscatastrophe flag so the two encodings of a catastrophe cannot
     * disagree (SolverMAM applies the same test).
     */
    public static boolean isCatastropheSignal(NetworkStruct sn, int jobClass) {
        if (sn.iscatastrophe != null && jobClass < sn.iscatastrophe.length()
                && sn.iscatastrophe.get(jobClass) > 0) {
            return true;
        }
        if (sn.signaltype != null && jobClass < sn.signaltype.size()) {
            SignalType st = sn.signaltype.get(jobClass);
            return st == SignalType.CATASTROPHE;
        }
        return false;
    }

    /**
     * Batch-size distribution of the jobs removed by a signal class when ntot
     * eligible jobs are present. Without a removal distribution a signal removes
     * exactly one job. With one, the pmf is clipped at ntot: a batch larger than
     * the eligible population empties it instead of driving the queue negative, so
     * the tail P(B >= ntot) lumps onto "remove ntot".
     *
     * @return a two-row matrix: row 0 holds the batch sizes, row 1 their probabilities.
     */
    public static Matrix signalBatchPMF(NetworkStruct sn, int jobClass, int ntot) {
        if (sn.signalremdist == null || jobClass >= sn.signalremdist.size()
                || sn.signalremdist.get(jobClass) == null) {
            return singleRemoval();
        }
        DiscreteDistribution dist = sn.signalremdist.get(jobClass);

        // head holds B = 0 .. ntot-1, the tail lumps P(B >= ntot) onto "remove ntot"
        List<Double> kvals = new ArrayList<Double>();
        List<Double> kprobs = new ArrayList<Double>();
        double headsum = 0.0;
        for (int b = 0; b < ntot; b++) {
            double p = dist.evalPMF(b);
            headsum += p;
            kvals.add((double) b);
            kprobs.add(p);
        }
        kvals.add((double) ntot);
        kprobs.add(Math.max(0.0, 1.0 - headsum));

        double total = 0.0;
        List<Double> keptVals = new ArrayList<Double>();
        List<Double> keptProbs = new ArrayList<Double>();
        for (int i = 0; i < kvals.size(); i++) {
            if (kprobs.get(i) > 0) {
                keptVals.add(kvals.get(i));
                keptProbs.add(kprobs.get(i));
                total += kprobs.get(i);
            }
        }
        if (keptVals.isEmpty()) {
            return singleRemoval();
        }
        Matrix out = new Matrix(2, keptVals.size());
        for (int i = 0; i < keptVals.size(); i++) {
            out.set(0, i, keptVals.get(i));
            out.set(1, i, keptProbs.get(i) / total);
        }
        return out;
    }

    private static Matrix singleRemoval() {
        Matrix out = new Matrix(2, 1);
        out.set(0, 0, 1.0);
        out.set(1, 0, 1.0);
        return out;
    }

    /**
     * Passive arrival of a G-network signal class at station ist.
     *
     * <p>Removal is in general a random choice, so the generator needs every
     * destination state and its probability. A simulation instead advances along one
     * sample path, so under isSimulation the set is collapsed to a single successor
     * drawn from outprob (the rate stays -1, as for the other passive actions).</p>
     *
     * @return the weighted successor states; the rate is -1 throughout (passive action).
     */
    public static Ret.EventResult handleSignalArrival(NetworkStruct sn, int ind, int ist, Matrix inspace, int jobClass,
                                                      boolean isSimulation,
                                                      Matrix phasessz, Matrix phaseshift, Matrix K, Matrix Ks, Matrix S,
                                                      Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar) {
        State.StateMarginalStatistics sgs =
                ToMarginal.toMarginal(sn, ind, inspace, phasessz, phaseshift, spaceBuf, spaceSrv, spaceVar);

        double[] buf = rowOf(spaceBuf);
        double[] srv = rowOf(spaceSrv);
        double[] var = rowOf(spaceVar);

        if (isCatastropheSignal(sn, jobClass)) {
            double[] out = concat(new double[buf.length], new double[srv.length], var);
            return single(out);
        }

        // Eligible victim classes.
        int tgt = -1;
        if (sn.signaltarget != null && jobClass < sn.signaltarget.length()) {
            tgt = (int) Math.round(sn.signaltarget.get(jobClass));
        }
        List<Integer> tgtclasses = new ArrayList<Integer>();
        if (tgt >= 0) {
            tgtclasses.add(tgt);
        } else {
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.issignal == null || sn.issignal.get(r) == 0) {
                    tgtclasses.add(r);
                }
            }
        }
        List<Integer> eligible = new ArrayList<Integer>();
        double ntotd = 0.0;
        for (int i = 0; i < tgtclasses.size(); i++) {
            int r = tgtclasses.get(i);
            double nr = sgs.nir.get(0, r);
            if (nr > 0) {
                eligible.add(r);
                ntotd += nr;
            }
        }
        int ntot = (int) Math.round(ntotd);
        if (eligible.isEmpty() || ntot <= 0) {
            return single(concat(buf, srv, var)); // no victim: the signal vanishes
        }

        Matrix pmf = signalBatchPMF(sn, jobClass, ntot);

        RemovalPolicy policy = RemovalPolicy.RANDOM;
        if (sn.signalrempolicy != null && jobClass < sn.signalrempolicy.size()
                && sn.signalrempolicy.get(jobClass) != null) {
            policy = sn.signalrempolicy.get(jobClass);
        }

        SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
        double nservers = S.get(ist);

        Outcome acc = new Outcome();
        for (int ik = 0; ik < pmf.getNumCols(); ik++) {
            double kprob = pmf.get(1, ik);
            if (kprob <= 0) {
                continue;
            }
            int k = (int) Math.round(pmf.get(0, ik));
            if (k <= 0) {
                acc.add(concat(buf, srv, var), kprob);
                continue;
            }
            Outcome batch = removeBatch(sched, buf, srv, var, k, eligible, policy, K, Ks, nservers);
            for (int j = 0; j < batch.space.size(); j++) {
                acc.add(batch.space.get(j), kprob * batch.prob.get(j));
            }
        }

        // Merge duplicate destination states so the generator sees one entry each.
        acc = acc.merged();
        if (isSimulation && acc.space.size() > 1) {
            acc = acc.sampleOne();
        }
        int n = acc.space.size();
        int width = acc.space.get(0).length;
        Matrix outspace = new Matrix(n, width);
        Matrix outrate = new Matrix(n, 1);
        Matrix outprob = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double[] row = acc.space.get(i);
            for (int c = 0; c < width; c++) {
                outspace.set(i, c, row[c]);
            }
            outrate.set(i, 0, -1.0); // passive action
            outprob.set(i, 0, acc.prob.get(i));
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    /**
     * Remove k jobs one at a time; sequential uniform draws without replacement
     * reproduce a uniform choice of the removed subset.
     */
    private static Outcome removeBatch(SchedStrategy sched, double[] buf, double[] srv, double[] var, int k,
                                       List<Integer> tgtclasses, RemovalPolicy policy, Matrix K, Matrix Ks,
                                       double nservers) {
        int nb = buf.length;
        int ns = srv.length;
        Outcome cur = new Outcome();
        cur.add(concat(buf, srv, var), 1.0);
        for (int step = 0; step < k; step++) {
            Outcome next = new Outcome();
            for (int row = 0; row < cur.space.size(); row++) {
                double[] st = cur.space.get(row);
                double[] b = slice(st, 0, nb);
                double[] s = slice(st, nb, nb + ns);
                double[] v = slice(st, nb + ns, st.length);
                Outcome one = removeOne(sched, b, s, v, tgtclasses, policy, K, Ks, nservers);
                if (one.space.isEmpty()) {
                    // nothing left to remove: the state is already drained
                    next.add(st, cur.prob.get(row));
                } else {
                    for (int j = 0; j < one.space.size(); j++) {
                        next.add(one.space.get(j), cur.prob.get(row) * one.prob.get(j));
                    }
                }
            }
            cur = next.merged();
        }
        return cur;
    }

    /** Enumerate the single-victim outcomes and their probabilities. */
    private static Outcome removeOne(SchedStrategy sched, double[] buf, double[] srv, double[] var,
                                     List<Integer> tgtclasses, RemovalPolicy policy, Matrix K, Matrix Ks,
                                     double nservers) {
        Outcome out = new Outcome();
        Waiting w = waitingVictims(sched, buf, tgtclasses);
        InService sv = inServiceVictims(srv, tgtclasses, K, Ks);
        double nwait = 0.0;
        for (int i = 0; i < w.weight.size(); i++) {
            nwait += w.weight.get(i);
        }
        double nsrv = 0.0;
        for (int i = 0; i < sv.count.size(); i++) {
            nsrv += sv.count.get(i);
        }
        if (nwait == 0 && nsrv == 0) {
            return out;
        }

        // FCFS/LCFS rank the waiting line by age, which only an ordered buffer
        // records; a per-class count buffer (SIRO/SEPT/LEPT) carries no age, so an
        // age-based policy degenerates to a uniform draw over the waiting jobs.
        boolean ageOrdered = w.isOrdered && (policy == RemovalPolicy.FCFS || policy == RemovalPolicy.LCFS);
        if (ageOrdered && nwait > 0) {
            int pick = 0;
            for (int i = 1; i < w.pos.size(); i++) {
                if (policy == RemovalPolicy.FCFS) {
                    if (w.pos.get(i) > w.pos.get(pick)) pick = i; // head of line: the last occupied slot
                } else {
                    if (w.pos.get(i) < w.pos.get(pick)) pick = i; // most recent arrival: the first occupied slot
                }
            }
            double[] b2 = dropWaiting(buf, w.pos.get(pick), w.isOrdered, w.isPairBuf, w.cls.get(pick));
            out.add(concat(b2, srv, var), 1.0);
            return out;
        }

        double total;
        if (policy == RemovalPolicy.RANDOM) {
            total = nwait + nsrv; // uniform over waiting and in-service jobs alike
        } else {
            total = nwait > 0 ? nwait : nsrv; // FCFS/LCFS drain the waiting line before the servers
        }

        if (nwait > 0) {
            for (int i = 0; i < w.pos.size(); i++) {
                double[] b2 = dropWaiting(buf, w.pos.get(i), w.isOrdered, w.isPairBuf, w.cls.get(i));
                out.add(concat(b2, srv, var), w.weight.get(i) / total);
            }
        }
        if (policy == RemovalPolicy.RANDOM || nwait == 0) {
            for (int j = 0; j < sv.cls.size(); j++) {
                if (sv.count.get(j) <= 0) {
                    continue;
                }
                double[][] bs = dropInService(sched, buf, srv, sv.cls.get(j), sv.phase.get(j), Ks, nservers);
                out.add(concat(bs[0], bs[1], var), sv.count.get(j) / total);
            }
        }
        return out.merged();
    }

    /**
     * Eligible waiting jobs, as (position, class, multiplicity). Three buffer layouts
     * are in use: an ordered list of class ids (FCFS family), an ordered list of
     * (class, phase) pairs (preemptive family), and a per-class count vector
     * (SIRO/SEPT/LEPT). Stations that admit every job into service have no waiting
     * line at all. Only the ordered layouts record arrival order.
     */
    private static Waiting waitingVictims(SchedStrategy sched, double[] buf, List<Integer> tgtclasses) {
        Waiting w = new Waiting();
        if (buf.length == 0) {
            return w;
        }
        switch (sched) {
            case FCFS:
            case HOL:
            case LCFS:
            case LCFSPRIO:
                w.isOrdered = true;
                for (int c = 0; c < buf.length; c++) {
                    int v = (int) Math.round(buf[c]);
                    if (v > 0 && tgtclasses.contains(v - 1)) {
                        w.pos.add(c);
                        w.cls.add(v - 1);
                        w.weight.add(1.0);
                    }
                }
                break;
            case FCFSPR:
            case FCFSPI:
            case FCFSPRPRIO:
            case FCFSPIPRIO:
            case LCFSPR:
            case LCFSPI:
            case LCFSPRPRIO:
            case LCFSPIPRIO:
                w.isOrdered = true;
                w.isPairBuf = true;
                for (int c = 0; c + 1 < buf.length; c += 2) {
                    int v = (int) Math.round(buf[c]);
                    if (v > 0 && tgtclasses.contains(v - 1)) {
                        w.pos.add(c);
                        w.cls.add(v - 1);
                        w.weight.add(1.0);
                    }
                }
                break;
            case SIRO:
            case SEPT:
            case LEPT:
                // per-class counts: each eligible class is one victim kind carrying as
                // many interchangeable jobs as its count
                for (int i = 0; i < tgtclasses.size(); i++) {
                    int r = tgtclasses.get(i);
                    if (r < buf.length && buf[r] > 0) {
                        w.pos.add(r);
                        w.cls.add(r);
                        w.weight.add(buf[r]);
                    }
                }
                break;
            default:
                // PS/INF/DPS/GPS/LPS and the source hold no waiting jobs
                break;
        }
        return w;
    }

    private static InService inServiceVictims(double[] srv, List<Integer> tgtclasses, Matrix K, Matrix Ks) {
        InService sv = new InService();
        for (int i = 0; i < tgtclasses.size(); i++) {
            int r = tgtclasses.get(i);
            for (int p = 0; p < (int) K.get(r); p++) {
                double n = srv[(int) Ks.get(r) + p];
                if (n > 0) {
                    sv.cls.add(r);
                    sv.phase.add(p);
                    sv.count.add(n);
                }
            }
        }
        return sv;
    }

    /**
     * Remove a waiting job, keeping the layout invariant (empty slots pad the left,
     * the head of line stays rightmost).
     */
    private static double[] dropWaiting(double[] buf, int pos, boolean isOrdered, boolean isPairBuf, int cls) {
        if (isPairBuf) {
            return deletePrepend(buf, pos, 2);
        } else if (isOrdered) {
            return deletePrepend(buf, pos, 1);
        }
        double[] out = buf.clone();
        out[cls] -= 1; // per-class count buffer
        return out;
    }

    /**
     * Remove an in-service job and, at a station that keeps a waiting line, pull the
     * head of line into the freed server.
     *
     * @return a two-element array holding the new buffer and the new server vector.
     */
    private static double[][] dropInService(SchedStrategy sched, double[] buf, double[] srv, int cls, int phase,
                                            Matrix Ks, double nservers) {
        double[] b = buf.clone();
        double[] s = srv.clone();
        s[(int) Ks.get(cls) + phase] -= 1;
        double busy = 0.0;
        for (int i = 0; i < s.length; i++) {
            busy += s[i];
        }
        if (b.length == 0 || busy >= nservers) {
            return new double[][]{b, s};
        }
        switch (sched) {
            case FCFS:
            case HOL:
            case LCFS:
            case LCFSPRIO: {
                int headpos = -1;
                for (int c = b.length - 1; c >= 0; c--) {
                    if (b[c] > 0) {
                        headpos = c;
                        break;
                    }
                }
                if (headpos >= 0) {
                    int promo = (int) Math.round(b[headpos]);
                    // the slot is vacated, not blanked in place: the waiting line stays
                    // right-aligned with the empty slots padding the left
                    b = deletePrepend(b, headpos, 1);
                    s[(int) Ks.get(promo - 1)] += 1;
                }
                break;
            }
            case FCFSPR:
            case FCFSPI:
            case FCFSPRPRIO:
            case FCFSPIPRIO:
            case LCFSPR:
            case LCFSPI:
            case LCFSPRPRIO:
            case LCFSPIPRIO: {
                int headpos = -1;
                for (int c = 0; c + 1 < b.length; c += 2) {
                    if (b[c] > 0) {
                        headpos = c;
                    }
                }
                if (headpos >= 0) {
                    int promo = (int) Math.round(b[headpos]);
                    int promophase = (int) Math.round(b[headpos + 1]);
                    if (promophase < 1) {
                        promophase = 1;
                    }
                    b = deletePrepend(b, headpos, 2);
                    s[(int) Ks.get(promo - 1) + promophase - 1] += 1; // resumes at its stored phase
                }
                break;
            }
            case SIRO:
            case SEPT:
            case LEPT: {
                // the freed server takes a waiting job; the count buffer carries no order,
                // so any waiting class is equally eligible. Promote the lowest-indexed
                // waiting class to keep the map single-valued: the SIRO service order is
                // itself resolved by the rate function.
                int promo = -1;
                for (int c = 0; c < b.length; c++) {
                    if (b[c] > 0) {
                        promo = c;
                        break;
                    }
                }
                if (promo >= 0) {
                    b[promo] -= 1;
                    s[(int) Ks.get(promo)] += 1;
                }
                break;
            }
            default:
                // no waiting line to promote from
                break;
        }
        return new double[][]{b, s};
    }

    // Delete `width` entries at `pos` and re-pad the same number of zeros on the left,
    // so the buffer keeps its canonical right-aligned form.
    private static double[] deletePrepend(double[] buf, int pos, int width) {
        double[] out = new double[buf.length];
        int w = width; // leave the leading columns as the padding zeros
        for (int c = 0; c < buf.length; c++) {
            if (c >= pos && c < pos + width) {
                continue;
            }
            out[w] = buf[c];
            w++;
        }
        return out;
    }

    private static double[] rowOf(Matrix m) {
        if (m == null || m.getNumCols() == 0 || m.getNumRows() == 0) {
            return new double[0];
        }
        double[] out = new double[m.getNumCols()];
        for (int c = 0; c < m.getNumCols(); c++) {
            out[c] = m.get(0, c);
        }
        return out;
    }

    private static double[] slice(double[] a, int from, int to) {
        double[] out = new double[to - from];
        System.arraycopy(a, from, out, 0, to - from);
        return out;
    }

    private static double[] concat(double[] a, double[] b, double[] c) {
        double[] out = new double[a.length + b.length + c.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        System.arraycopy(c, 0, out, a.length + b.length, c.length);
        return out;
    }

    private static Ret.EventResult single(double[] row) {
        Matrix outspace = new Matrix(1, row.length);
        for (int c = 0; c < row.length; c++) {
            outspace.set(0, c, row[c]);
        }
        Matrix outrate = new Matrix(1, 1);
        outrate.set(0, 0, -1.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1.0);
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    /** A weighted set of successor states. */
    private static class Outcome {
        List<double[]> space = new ArrayList<double[]>();
        List<Double> prob = new ArrayList<Double>();

        void add(double[] s, double p) {
            space.add(s);
            prob.add(p);
        }

        /**
         * Draw a single successor with probability proportional to its weight, as the
         * balking branch of MATLAB afterEventStation does for its passive actions.
         */
        Outcome sampleOne() {
            double tot = 0.0;
            for (int i = 0; i < prob.size(); i++) {
                tot += prob.get(i);
            }
            double rnd = Maths.rand();
            double cum = 0.0;
            int pick = space.size() - 1;
            for (int i = 0; i < space.size(); i++) {
                cum += prob.get(i) / tot;
                if (rnd <= cum) {
                    pick = i;
                    break;
                }
            }
            Outcome out = new Outcome();
            out.add(space.get(pick), 1.0);
            return out;
        }

        /** Merge duplicate rows, summing their probabilities and keeping first-seen order. */
        Outcome merged() {
            Outcome out = new Outcome();
            Map<String, Integer> seen = new LinkedHashMap<String, Integer>();
            for (int i = 0; i < space.size(); i++) {
                String key = keyOf(space.get(i));
                Integer at = seen.get(key);
                if (at == null) {
                    seen.put(key, out.space.size());
                    out.add(space.get(i), prob.get(i));
                } else {
                    out.prob.set(at, out.prob.get(at) + prob.get(i));
                }
            }
            return out;
        }

        private static String keyOf(double[] row) {
            StringBuilder sb = new StringBuilder();
            for (int c = 0; c < row.length; c++) {
                sb.append((long) Math.round(row[c])).append(',');
            }
            return sb.toString();
        }
    }

    private static class Waiting {
        List<Integer> pos = new ArrayList<Integer>();
        List<Integer> cls = new ArrayList<Integer>();
        List<Double> weight = new ArrayList<Double>();
        boolean isOrdered = false;
        boolean isPairBuf = false;
    }

    private static class InService {
        List<Integer> cls = new ArrayList<Integer>();
        List<Integer> phase = new ArrayList<Integer>();
        List<Double> count = new ArrayList<Double>();
    }
}
