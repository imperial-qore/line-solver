package jline.solvers.ssa.handlers;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.BalkingStrategy;
import jline.lang.constant.BalkingThreshold;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.PollingType;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.NodeParam;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.state.AfterEventStation;
import jline.lang.state.Polling;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.util.IndexedMinHeap;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

public final class Solver_ssa_nrm {
    private Solver_ssa_nrm() {}

    /**
     * Class-dependence factor for a class-r completion at station ist: the
     * class-r component of the 1xR scaling vector beta(n) returned by the
     * sn.cdscaling handle, evaluated on the per-class population of node ind
     * read from the NRM state vector X (see FesBetaFunction and
     * State.cdclassfactor / MATLAB solver_ssa_nrm cdfac). Returns 1 when the
     * station declares no class dependence.
     */
    static double cdfac(NetworkStruct sn, int ist, Matrix X, int ind, int R, int r) {
        return cdfacPop(sn, ist, stationPop(X, ind, R), R, r);
    }

    /**
     * Class-dependence factor evaluated at an explicit per-class population
     * vector. The DPSPRIO/GPSPRIO rate laws evaluate it at the
     * priority-restricted population rather than the full one, so the vector
     * cannot always be read straight out of the state.
     */
    static double cdfacPop(NetworkStruct sn, int ist, double[] npop, int R, int r) {
        jline.util.SerializableFunction<Matrix, Matrix> f =
                (sn.cdscaling != null) ? sn.cdscaling.get(sn.stations.get(ist)) : null;
        // Joint-dependence eta_i (non-product-form) is evaluated identically and
        // multiplied in, so a jd-only or cd+jd station is handled uniformly.
        jline.util.SerializableFunction<Matrix, Matrix> g =
                (sn.jdscaling != null) ? sn.jdscaling.get(sn.stations.get(ist)) : null;
        if (f == null && g == null) {
            return 1.0;
        }
        Matrix nvec = new Matrix(1, R);
        for (int rr = 0; rr < R; rr++) {
            nvec.set(0, rr, npop[rr]);
        }
        double fac = 1.0;
        if (f != null) {
            Matrix v = f.apply(nvec);
            if (v != null && v.length() > 0) {
                fac *= v.get(Math.min(r, v.length() - 1));
            }
        }
        if (g != null) {
            Matrix v = g.apply(nvec);
            if (v != null && v.length() > 0) {
                fac *= v.get(Math.min(r, v.length() - 1));
            }
        }
        return fac;
    }

    /**
     * Per-class population vector of the station sitting at node ind.
     */
    static double[] stationPop(Matrix X, int ind, int R) {
        double[] npop = new double[R];
        for (int r = 0; r < R; r++) {
            npop[r] = X.get(ind * R + r, 0);
        }
        return npop;
    }

    static double sumPop(double[] npop) {
        double s = 0.0;
        for (int r = 0; r < npop.length; r++) {
            s += npop[r];
        }
        return s;
    }

    /**
     * Layout of the NRM state vector: one slot per (node, class, service PHASE),
     * so phase-type service is represented exactly instead of being collapsed
     * onto its mean rate. Phase counts differ per (station, class) via
     * sn.phasessz, hence an explicit offset map rather than arithmetic on R.
     * <p>
     * The layout is backward compatible BY CONSTRUCTION: with nph == 1 everywhere
     * phOff[ind][r] == ind*R + r, so slot(ind,r,0) is exactly the flat class index
     * the engine used before phase expansion, and an exponential model reproduces
     * its previous state vector, reactions and rates unchanged. buildSmap asserts
     * that identity.
     * </p><p>
     * Every consumer of a state index must decode through {@link #node} /
     * {@link #cls} rather than through pos/R and pos%R: once any class carries
     * more than one phase the arithmetic decode names the WRONG node, and the
     * resulting stale propensities are invisible to an exponential model.
     * </p>
     */
    static final class Smap {
        int R;
        int NS;              // total number of state slots
        int[][] nph;         // nnodes x R, phases of each (node,class)
        int[][] phOff;       // nnodes x R, first slot of each (node,class)
        int[] node;          // NS, owning node of each slot
        int[] cls;           // NS, owning class of each slot
        int[] phase;         // NS, phase index within its (node,class)
        boolean expanded;    // NS > nnodes*R, i.e. some class has several phases

        int slot(int ind, int r, int k) {
            return phOff[ind][r] + k;
        }
    }

    /**
     * Build the slot map. Non-station nodes, and station-classes whose service is
     * exponential or disabled, carry a single phase.
     */
    static Smap buildSmap(NetworkStruct sn, int I, int R) {
        Smap sm = new Smap();
        sm.R = R;
        sm.nph = new int[I][R];
        sm.phOff = new int[I][R];
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) {
                int np = 1;
                if (sn.isstation.get(ind) != 0.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (ist >= 0) {
                        np = Math.max(1, (int) sn.phasessz.get(ist, r));
                    }
                }
                sm.nph[ind][r] = np;
            }
        }
        int ns = 0;
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) {
                sm.phOff[ind][r] = ns;
                ns += sm.nph[ind][r];
            }
        }
        sm.NS = ns;
        sm.expanded = ns > I * R;
        sm.node = new int[ns];
        sm.cls = new int[ns];
        sm.phase = new int[ns];
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) {
                for (int k = 0; k < sm.nph[ind][r]; k++) {
                    int s = sm.phOff[ind][r] + k;
                    sm.node[s] = ind;
                    sm.cls[s] = r;
                    sm.phase[s] = k;
                }
            }
        }
        if (!sm.expanded) {
            // Self-check on the generalization: a single-phase model must land on
            // the pre-expansion flat layout, so every exponential model is
            // guaranteed to reproduce its previous results.
            for (int ind = 0; ind < I; ind++) {
                for (int r = 0; r < R; r++) {
                    if (sm.phOff[ind][r] != ind * R + r) {
                        throw new RuntimeException(
                                "NRM slot map: a single-phase model must reduce to the flat (node,class) layout.");
                    }
                }
            }
        }
        return sm;
    }

    /**
     * Per-class populations at node ind, each class summed over its phases. The
     * scheduling rate laws are class-level and are unchanged by phase expansion;
     * only the per-phase share (see {@link #kirFrac}) is layered on top.
     */
    static double[] classCounts(Matrix X, Smap sm, int ind) {
        double[] npop = new double[sm.R];
        for (int r = 0; r < sm.R; r++) {
            npop[r] = classPop(X, sm, ind, r);
        }
        return npop;
    }

    /** Integer per-class populations at node ind (rounded), used by the polling
     * controller walk, which reasons over waiting counts. */
    static int[] classCountsInt(Matrix X, Smap sm, int ind, int R) {
        int[] out = new int[R];
        for (int r = 0; r < R; r++) {
            out[r] = (int) Math.round(classPop(X, sm, ind, r));
        }
        return out;
    }

    /** Index drawn from the (unnormalized, nonnegative) weight vector w. */
    static int pollDraw(double[] w) {
        double tot = 0.0;
        for (int i = 0; i < w.length; i++) tot += w[i];
        if (tot <= 0.0) return 0;
        double u = Maths.rand() * tot;
        double c = 0.0;
        for (int i = 0; i < w.length; i++) {
            c += w[i];
            if (u < c) return i;
        }
        return w.length - 1;
    }

    /**
     * Controller row {mode,pos,swk,ctr} the server lands in after Polling.next
     * resolves (q, mode, budget): SERVING q with the visit budget, SWITCHING into
     * q with the entry phase drawn from the switchover PH, or PARKED at q.
     * Specialized to the single-server polling station the NRM carries (exponential
     * service, so no in-service phase).
     */
    static int[] pollLand(Polling.Info pinf, int q, int mode, int budget) {
        if (mode == Polling.MODE_VISIT) {
            return new int[]{Polling.MODE_VISIT, q, 0, budget};
        }
        if (mode == Polling.MODE_SWITCH) {
            Matrix pie = pinf.swPie[q];
            double[] w = new double[pie.getNumCols()];
            for (int i = 0; i < w.length; i++) w[i] = pie.get(0, i);
            return new int[]{Polling.MODE_SWITCH, q, pollDraw(w), 0};
        }
        return new int[]{Polling.MODE_PARK, q, 0, 0};
    }

    /**
     * Advance the polling controller of any node the firing of reaction kfire
     * touched, mirroring State.afterEventStation for POLLING: a service completion
     * ends the visit unless the discipline still admits another job of the served
     * class and on ending walks the cyclic order; a switchover reaction advances
     * the switchover PH one phase or on absorption arrives at the target buffer;
     * an arrival to a PARKED server wakes it and the walk resolves at once.
     * Returns true when a controller moved, forcing a full propensity refresh.
     */
    static boolean pollAdvance(int kfire, int destPos, Matrix nvec, Smap sm, int R,
            Polling.Info[] pinfo, boolean[] isPollNode, int[][] pollCtrl,
            boolean[] isPollSwRx, int[] pollSwNode, List<int[]> fromIR, int nDepRx, boolean[] isPhaseRx) {
        boolean changed = false;
        int srcNode = fromIR.get(kfire)[0];
        if (kfire < isPollSwRx.length && isPollSwRx[kfire]) {
            int pind = pollSwNode[kfire];
            Polling.Info pinf = pinfo[pind];
            int[] ctrl = pollCtrl[pind];
            int pos = ctrl[1];
            int swk = ctrl[2];
            Matrix D0 = pinf.swD0[pos];
            int Ksw = pinf.Ksw[pos];
            double[] w = new double[Ksw + 1];
            for (int kd = 0; kd < Ksw; kd++) {
                if (kd != swk && D0.get(swk, kd) > 0.0) w[kd] = D0.get(swk, kd);
            }
            double rowsum = 0.0;
            for (int c = 0; c < Ksw; c++) rowsum += D0.get(swk, c);
            w[Ksw] = Math.max(0.0, -rowsum);   // absorption (D1 row sum)
            int pick = pollDraw(w);
            if (pick < Ksw) {
                ctrl[2] = pick;                 // internal phase advance
            } else {
                int[] nb = classCountsInt(nvec, sm, pind, R);
                int[] nx = Polling.next(pinf, pos, nb, R, true);
                pollCtrl[pind] = pollLand(pinf, nx[0], nx[1], nx[2]);
            }
            changed = true;
        } else if (isPollNode[srcNode] && kfire < nDepRx && !isPhaseRx[kfire]) {
            Polling.Info pinf = pinfo[srcNode];
            int[] ctrl = pollCtrl[srcNode];
            int pos = ctrl[1];
            int ctr = ctrl[3];
            int[] nb = classCountsInt(nvec, sm, srcNode, R);   // after the departure
            int ctrnext;
            boolean goon;
            switch (pinf.ptype) {
                case EXHAUSTIVE:
                    ctrnext = 0;
                    goon = nb[pos] > 0;
                    break;
                case GATED:
                    ctrnext = ctr - 1;
                    goon = ctrnext > 0;
                    break;
                case KLIMITED:
                    ctrnext = ctr - 1;
                    goon = ctrnext > 0 && nb[pos] > 0;
                    break;
                case DECREMENTING:
                    ctrnext = ctr;
                    goon = nb[pos] > ctr;
                    break;
                default:
                    ctrnext = 0;
                    goon = false;
                    break;
            }
            if (goon) {
                pollCtrl[srcNode] = new int[]{Polling.MODE_VISIT, pos, 0, ctrnext};
            } else {
                int[] nx = Polling.next(pinf, pos, nb, R, false);
                pollCtrl[srcNode] = pollLand(pinf, nx[0], nx[1], nx[2]);
            }
            changed = true;
        }
        if (destPos >= 0) {
            int jnd = sm.node[destPos];
            if (isPollNode[jnd] && pollCtrl[jnd] != null && pollCtrl[jnd][0] == Polling.MODE_PARK) {
                Polling.Info pinf = pinfo[jnd];
                int[] nb = classCountsInt(nvec, sm, jnd, R);   // includes the arrival
                int[] nx = Polling.next(pinf, pollCtrl[jnd][1], nb, R, true);
                pollCtrl[jnd] = pollLand(pinf, nx[0], nx[1], nx[2]);
                changed = true;
            }
        }
        return changed;
    }

    /** Population of class r at node ind, summed over its phases. */
    static double classPop(Matrix X, Smap sm, int ind, int r) {
        double n = 0.0;
        int off = sm.phOff[ind][r];
        for (int k = 0; k < sm.nph[ind][r]; k++) {
            n += X.get(off + k, 0);
        }
        return n;
    }

    /**
     * Share of its class that the jobs in one phase represent: kir/nir. The
     * class-level rate law is split across the class's phases in this ratio,
     * which is how State.afterEventStation writes every phase-aware case (DPS,
     * for instance, uses (kir/nir) * [class share]). For a single-phase class
     * this is 1 whenever the class is present, so an exponential model is
     * unaffected.
     */
    static double kirFrac(Matrix X, Smap sm, int slot, int ind, int r) {
        double nir = classPop(X, sm, ind, r);
        if (nir <= 0.0) {
            return 0.0;
        }
        return X.get(slot, 0) / nir;
    }

    /**
     * Entry-phase distribution of a class-s job arriving at node jnd: pie of its
     * service process there. A non-station node, or a station whose process is
     * absent (a disabled class), has a single phase entered with probability 1.
     */
    static double[] entryProbs(NetworkStruct sn, Smap sm, int jnd, int s) {
        int np = sm.nph[jnd][s];
        double[] pe = new double[np];
        if (np <= 1) {
            pe[0] = 1.0;
            return pe;
        }
        int ist = (int) sn.nodeToStation.get(jnd);
        Matrix pie = null;
        if (sn.pie != null && sn.pie.get(sn.stations.get(ist)) != null) {
            pie = sn.pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(s));
        }
        double tot = 0.0;
        if (pie != null) {
            for (int k = 0; k < Math.min(np, pie.length()); k++) {
                double v = pie.get(k);
                if (!Double.isNaN(v) && v > 0.0) {
                    pe[k] = v;
                    tot += v;
                }
            }
        }
        if (tot <= 0.0) {
            // no entry distribution declared: enter the first phase
            for (int k = 0; k < np; k++) pe[k] = 0.0;
            pe[0] = 1.0;
            return pe;
        }
        for (int k = 0; k < np; k++) {
            pe[k] /= tot;
        }
        return pe;
    }

    /**
     * Normalized per-class DPS/GPS weights of station ist, from sn.schedparam.
     */
    static double[] schedWeights(NetworkStruct sn, int ist, int R) {
        double[] w = new double[R];
        double tot = 0.0;
        for (int r = 0; r < R; r++) {
            w[r] = sn.schedparam.get(ist, r);
            tot += w[r];
        }
        if (tot <= 0.0) {
            throw new RuntimeException("Station " + ist
                    + " has weighted (DPS/GPS) scheduling with non-positive total weight.");
        }
        for (int r = 0; r < R; r++) {
            w[r] /= tot;
        }
        return w;
    }

    /**
     * DPS sharing factor: class r receives w_r*n_r/(w.n) of the single server.
     */
    static double dpsshare(double[] w, double[] npop, int r) {
        double den = 0.0;
        for (int s = 0; s < npop.length; s++) {
            den += w[s] * npop[s];
        }
        if (den <= 0.0) {
            return 0.0;
        }
        return w[r] * npop[r] / den;
    }

    /**
     * GPS sharing factor: w_r/(w.c) with c_s = 1{n_s>0}, so the weights are
     * split across the active classes rather than across the jobs.
     */
    static double gpsshare(double[] w, double[] npop, int r) {
        if (npop[r] <= 0.0) {
            return 0.0;
        }
        double den = 0.0;
        for (int s = 0; s < npop.length; s++) {
            if (npop[s] > 0.0) {
                den += w[s];
            }
        }
        if (den <= 0.0) {
            return 0.0;
        }
        return w[r] / den;
    }

    /**
     * True when class r belongs to the most urgent non-empty priority group
     * (LINE orders priorities with lower value = more urgent).
     */
    static boolean isUrgent(NetworkStruct sn, double[] npop, int r) {
        double minPrio = Double.MAX_VALUE;
        boolean any = false;
        for (int s = 0; s < npop.length; s++) {
            if (npop[s] > 0.0) {
                any = true;
                double p = sn.classprio.get(s);
                if (p < minPrio) {
                    minPrio = p;
                }
            }
        }
        return any && sn.classprio.get(r) == minPrio;
    }

    /**
     * Population vector restricted to the priority group of class r.
     */
    static double[] prioGroup(NetworkStruct sn, double[] npop, int r) {
        double[] act = new double[npop.length];
        for (int s = 0; s < npop.length; s++) {
            if (sn.classprio.get(s) == sn.classprio.get(r)) {
                act[s] = npop[s];
            }
        }
        return act;
    }

    /**
     * Population the lld factor is evaluated at: the full station population
     * below capacity, the priority-group population above it.
     */
    static double prioPop(NetworkStruct sn, double[] npop, int r, double c) {
        double ni = sumPop(npop);
        if (ni <= c || !isUrgent(sn, npop, r)) {
            return ni;
        }
        return sumPop(prioGroup(sn, npop, r));
    }

    /**
     * Population vector the cd factor is evaluated at for DPSPRIO/GPSPRIO.
     * PSPRIO instead uses the full vector in both branches; that asymmetry is
     * inherited from State.afterEventStation and reproduced here.
     */
    static double[] prioVec(NetworkStruct sn, double[] npop, int r, double c) {
        double ni = sumPop(npop);
        if (ni <= c || !isUrgent(sn, npop, r)) {
            return npop;
        }
        return prioGroup(sn, npop, r);
    }

    /**
     * PSPRIO: PS below capacity; above it only the most urgent non-empty group
     * shares the servers and every other class is frozen.
     */
    static double psprioshare(NetworkStruct sn, double[] npop, int r, double c) {
        double ni = sumPop(npop);
        if (ni <= 0.0) {
            return 0.0;
        }
        if (ni <= c) {
            return (npop[r] / ni) * Math.min(ni, c);
        }
        if (!isUrgent(sn, npop, r)) {
            return 0.0;
        }
        double niprio = sumPop(prioGroup(sn, npop, r));
        return (npop[r] / niprio) * Math.min(niprio, c);
    }

    /**
     * DPSPRIO: DPS below capacity, DPS restricted to the urgent group above it.
     */
    static double dpsprioshare(NetworkStruct sn, double[] w, double[] npop, int r, double c) {
        double ni = sumPop(npop);
        if (ni <= 0.0) {
            return 0.0;
        }
        if (ni <= c) {
            return dpsshare(w, npop, r);
        }
        if (!isUrgent(sn, npop, r)) {
            return 0.0;
        }
        return dpsshare(w, prioGroup(sn, npop, r), r);
    }

    /**
     * Finite capacity regions (DROP rule).
     * <p>
     * A region constrains an aggregate of the per-class populations of its
     * member stations, which is a linear function of the NRM state vector, so
     * admission is a multiplicative 0/1 gate on the routing draw. The DROP rule
     * censors the refused transition, and censoring an exponential transition is
     * exactly what zeroing its share of the propensity does. WAITQ instead parks
     * a refused job in a per-region FIFO and admits it head-of-line as capacity
     * frees; the FIFO is carried explicitly here (see {@link #fcrReleaseCascade})
     * so both rules run in this engine.
     * </p>
     */
    static final class Fcr {
        boolean on = false;
        int F = 0;
        boolean[][] memberNode;   // F x nnodes
        double[][] classCap;      // F x K
        double[] globalCap;
        double[] memCap;
        double[][] sz;            // F x K
        Matrix[] A;
        Matrix[] b;
        // Per-(region,class) admission rule: DROP destroys a refused job, WAITQ
        // parks it in the region FIFO. Mirrors MATLAB fcr.waitq
        // (regionrule ~= DropStrategy.DROP). Compared by NAME through
        // DropStrategy.fromID, never by raw id, because the DropStrategy numeric
        // encodings are not portable across codebases.
        boolean[][] waitq;        // F x K
        boolean anyWaitq = false;
    }

    /**
     * Per-region member nodes and admission caps, mirroring the FCR precompute
     * of the serial engine (Solver_ssa) and of MATLAB solver_ssa_nrm.
     */
    static Fcr fcrPrecompute(NetworkStruct sn) {
        Fcr fcr = new Fcr();
        if (sn.nregions == 0) {
            return fcr;
        }
        int K = sn.nclasses;
        int M = sn.nstations;
        int F = sn.nregions;
        fcr.on = true;
        fcr.F = F;
        fcr.memberNode = new boolean[F][sn.nnodes];
        fcr.classCap = new double[F][];
        fcr.globalCap = new double[F];
        fcr.memCap = new double[F];
        fcr.sz = new double[F][];
        fcr.A = new Matrix[F];
        fcr.b = new Matrix[F];
        // Per-(region,class) DROP/WAITQ rule, compared BY NAME. A region with no
        // WAITQ class carries no FIFO, so pure-DROP models pay nothing extra.
        fcr.waitq = new boolean[F][K];
        if (sn.regionrule != null && !sn.regionrule.isEmpty()) {
            for (int f = 0; f < F; f++) {
                for (int r = 0; r < K; r++) {
                    DropStrategy rule = DropStrategy.fromID((int) sn.regionrule.get(f, r));
                    if (rule != DropStrategy.Drop) {
                        fcr.waitq[f][r] = true;
                        fcr.anyWaitq = true;
                    }
                }
            }
        }
        for (int f = 0; f < F; f++) {
            Matrix Rmat = sn.region.get(f);            // M x (K+1)
            double[] memvec = new double[M];
            for (int i = 0; i < M; i++) memvec[i] = -1.0;
            if (sn.regionmaxmem != null && sn.regionmaxmem.size() > f && sn.regionmaxmem.get(f) != null) {
                Matrix mm = sn.regionmaxmem.get(f);
                for (int i = 0; i < Math.min(M, mm.length()); i++) memvec[i] = mm.get(i);
            }
            List<Integer> members = new ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                boolean anyCap = false;
                for (int c = 0; c < Rmat.getNumCols(); c++) {
                    if (Rmat.get(i, c) != -1.0) { anyCap = true; break; }
                }
                if (anyCap || memvec[i] != -1.0) members.add(i);
            }
            double[] ccap = new double[K];
            for (int r = 0; r < K; r++) {
                double best = Double.POSITIVE_INFINITY;
                for (Integer i : members) {
                    double v = Rmat.get(i, r);
                    if (v != -1.0 && v < best) best = v;
                }
                ccap[r] = best;
            }
            fcr.classCap[f] = ccap;
            double gbest = Double.POSITIVE_INFINITY;
            for (Integer i : members) {
                double v = Rmat.get(i, K);
                if (v != -1.0 && v < gbest) gbest = v;
            }
            fcr.globalCap[f] = gbest;
            double mbest = Double.POSITIVE_INFINITY;
            for (Integer i : members) {
                if (memvec[i] != -1.0 && memvec[i] < mbest) mbest = memvec[i];
            }
            fcr.memCap[f] = mbest;
            double[] szf = new double[K];
            for (int r = 0; r < K; r++) szf[r] = sn.regionsz.get(f, r);
            fcr.sz[f] = szf;
            if (sn.regionlincon != null && sn.regionlincon.containsKey(f)
                    && sn.regionlincon.get(f) != null && sn.regionlincon.get(f).size() >= 2
                    && sn.regionlincon.get(f).get(0) != null
                    && sn.regionlincon.get(f).get(0).length() > 0) {
                fcr.A[f] = sn.regionlincon.get(f).get(0);
                fcr.b[f] = sn.regionlincon.get(f).get(1);
            }
            for (Integer i : members) {
                fcr.memberNode[f][(int) sn.stationToNode.get(i)] = true;
            }
        }
        return fcr;
    }

    /** True if per-class population vector xn breaks any admission constraint. */
    static boolean fcrViolates(double[] xn, double[] ccap, double gcap, double memcap,
                               double[] sz, Matrix A, Matrix b) {
        double tot = 0.0;
        double mem = 0.0;
        for (int r = 0; r < xn.length; r++) {
            if (xn[r] > ccap[r]) return true;
            tot += xn[r];
            mem += xn[r] * sz[r];
        }
        if (tot > gcap || mem > memcap) return true;
        if (A != null) {
            for (int i = 0; i < A.getNumRows(); i++) {
                double lhs = 0.0;
                for (int r = 0; r < Math.min(xn.length, A.getNumCols()); r++) {
                    lhs += A.get(i, r) * xn[r];
                }
                if (lhs > b.get(i)) return true;
            }
        }
        return false;
    }

    /** Per-class population of a region, read off the NRM state vector. */
    static double[] fcrRegionPop(Matrix nvec, boolean[] memberRow, int R, Smap sm) {
        double[] x = new double[R];
        for (int jnd = 0; jnd < memberRow.length; jnd++) {
            if (!memberRow[jnd]) continue;
            for (int r = 0; r < R; r++) x[r] += classPop(nvec, sm, jnd, r);
        }
        return x;
    }

    /**
     * True if a class-dstClass job may enter dstNode, having just left srcNode
     * as srcClass. Only regions containing the destination can refuse the move;
     * a move whose source is in the same region frees a slot first, so the
     * departure is accounted for before the arrival is tested.
     */
    static boolean fcrAdmits(Fcr fcr, Matrix nvec, int srcNode, int srcClass,
                             int dstNode, int dstClass, int R, Smap sm) {
        return fcrRefusingRegion(fcr, nvec, srcNode, srcClass, dstNode, dstClass, R, sm) < 0;
    }

    /**
     * Index of the FIRST region that refuses a class-dstClass job entering
     * dstNode, having just left srcNode as srcClass; -1 if every region admits
     * it. Same admission test as {@link #fcrAdmits}, but it names the refusing
     * region so the caller can consult that region's DROP/WAITQ rule. Mirrors
     * MATLAB fcrRefusingRegion. A srcNode &lt; 0 means the mover has no live
     * source in the state (a WAITQ release, whose job already left its source
     * when it was parked), so no source slot is freed.
     */
    static int fcrRefusingRegion(Fcr fcr, Matrix nvec, int srcNode, int srcClass,
                                 int dstNode, int dstClass, int R, Smap sm) {
        if (!fcr.on) return -1;
        for (int f = 0; f < fcr.F; f++) {
            if (!fcr.memberNode[f][dstNode]) continue; // does not constrain the destination
            double[] x = fcrRegionPop(nvec, fcr.memberNode[f], R, sm);
            if (srcNode >= 0 && fcr.memberNode[f][srcNode]) x[srcClass] -= 1.0;
            x[dstClass] += 1.0;
            if (fcrViolates(x, fcr.classCap[f], fcr.globalCap[f], fcr.memCap[f],
                    fcr.sz[f], fcr.A[f], fcr.b[f])) {
                return f;
            }
        }
        return -1;
    }

    /**
     * Strict-FIFO head-of-line release of parked WAITQ tokens: admit each
     * region's FIFO head while the admission constraints permit, applying the
     * arrival to the destination station (entry-phase slot plus buffer join).
     * Mirrors MATLAB fcrReleaseCascade. A token encodes (dstNode, dstClass) as
     * dstNode*R + dstClass; the entry phase is drawn at release, as a routed
     * arrival draws it. Loops until a full pass frees nothing, so a release that
     * frees capacity elsewhere cascades. Returns the number of jobs admitted.
     */
    static int fcrReleaseCascade(Fcr fcr, Matrix nvec, ArrayDeque<Integer>[] buffers,
                                 ArrayDeque<Integer>[] fcrBuf, double[] mi, int R,
                                 NetworkStruct sn, Smap sm,
                                 double[][][] svcph, boolean[] bufPHNode) {
        int released = 0;
        boolean progress = true;
        while (progress) {
            progress = false;
            for (int f = 0; f < fcrBuf.length; f++) {
                if (fcrBuf[f].isEmpty()) continue;
                int tok = fcrBuf[f].peekFirst();
                int dstNode = tok / R;
                int dstClass = tok % R;
                // The parked job already left its source, so admission is tested
                // with the source term absent (srcNode = -1).
                if (fcrRefusingRegion(fcr, nvec, -1, dstClass, dstNode, dstClass, R, sm) >= 0) {
                    continue; // head-of-line: this FIFO stays blocked
                }
                if (bufPHNode[dstNode]) {
                    // Buffered-PH destination: the released job lands in the class
                    // total slot; whether it enters service (and its entry phase) is
                    // decided in applyArrivalBuffer against the server occupancy.
                    int dslot = sm.phOff[dstNode][dstClass];
                    nvec.set(dslot, 0, nvec.get(dslot, 0) + 1.0);
                } else {
                    double[] pentry = entryProbs(sn, sm, dstNode, dstClass);
                    int ke = drawFromDist(pentry);
                    int dslot = sm.phOff[dstNode][dstClass] + ke;
                    nvec.set(dslot, 0, nvec.get(dslot, 0) + 1.0);
                }
                applyArrivalBuffer(dstNode, dstClass, nvec, buffers, mi, R, sn, sm, svcph, bufPHNode);
                fcrBuf[f].pollFirst();
                released++;
                progress = true;
            }
        }
        return released;
    }

    /** Index drawn from the (unnormalized, nonnegative) weight vector p. */
    static int drawFromDist(double[] p) {
        double tot = 0.0;
        for (int i = 0; i < p.length; i++) tot += p[i];
        if (tot <= 0.0) return 0;
        double u = Maths.rand();
        double c = 0.0;
        for (int i = 0; i < p.length; i++) {
            c += p[i] / tot;
            if (c > u) return i;
        }
        return p.length - 1;
    }

    /**
     * Share of the routing draw the regions currently admit. A reaction with no
     * destination (self-looping class) moves no job across a region boundary and
     * is never gated.
     */

    /**
     * Balking (QUEUE_LENGTH strategy).
     * <p>
     * An arrival that balks is lost: it has left its source but never joins the
     * destination, so the departure rate is unchanged and only the arrival
     * outcome differs (State.afterEventStation scales the admitted branches by
     * 1-balkProb and adds a balked branch of probability balkProb that leaves
     * the destination state untouched). Only QUEUE_LENGTH is a pure function of
     * the state vector; EXPECTED_WAIT / COMBINED depend on the mean wait and are
     * rejected by Solver_ssa_analyzer.
     * </p>
     */
    static final class Balk {
        boolean on = false;
        int[] nodeToStation;                    // node -> station, -1 for non-stations
        BalkingStrategy[][] strategy;           // station x class, null = none
        List<BalkingThreshold>[][] thresholds;  // station x class
    }

    /**
     * Per (station,class) balking threshold table, indexed by node so the draw
     * can be evaluated straight off the NRM state vector.
     */
    @SuppressWarnings("unchecked")
    static Balk balkPrecompute(NetworkStruct sn) {
        Balk balk = new Balk();
        if (sn.balkingStrategy == null || sn.balkingStrategy.isEmpty()) {
            return balk;
        }
        int M = sn.nstations;
        int K = sn.nclasses;
        BalkingStrategy[][] strategy = new BalkingStrategy[M][K];
        boolean any = false;
        for (int ist = 0; ist < M; ist++) {
            Map<JobClass, BalkingStrategy> smap = sn.balkingStrategy.get(sn.stations.get(ist));
            if (smap == null) continue;
            for (int r = 0; r < K; r++) {
                BalkingStrategy s = smap.get(sn.jobclasses.get(r));
                strategy[ist][r] = s;
                if (s == BalkingStrategy.QUEUE_LENGTH) any = true;
            }
        }
        if (!any) {
            return balk;
        }
        balk.on = true;
        balk.strategy = strategy;
        balk.thresholds = (List<BalkingThreshold>[][]) new List[M][K];
        for (int ist = 0; ist < M; ist++) {
            Map<JobClass, List<BalkingThreshold>> tmap =
                    (sn.balkingThresholds != null) ? sn.balkingThresholds.get(sn.stations.get(ist)) : null;
            if (tmap == null) continue;
            for (int r = 0; r < K; r++) {
                balk.thresholds[ist][r] = tmap.get(sn.jobclasses.get(r));
            }
        }
        balk.nodeToStation = new int[sn.nnodes];
        for (int ind = 0; ind < sn.nnodes; ind++) {
            balk.nodeToStation[ind] = (sn.isstation.get(ind) != 0.0) ? (int) sn.nodeToStation.get(ind) : -1;
        }
        return balk;
    }

    /**
     * True if the job routed to state slot destPos balks. The threshold table is
     * scanned in order and the FIRST interval containing the pre-arrival total
     * station population wins, matching State.afterEventStation.
     */
    static boolean balkDraw(Balk balk, Matrix nvec, int destPos, int R, Smap sm) {
        int jnd = sm.node[destPos];
        int s = sm.cls[destPos];
        int ist = balk.nodeToStation[jnd];
        if (ist < 0 || balk.strategy[ist][s] != BalkingStrategy.QUEUE_LENGTH) {
            return false;
        }
        List<BalkingThreshold> th = balk.thresholds[ist][s];
        if (th == null) {
            return false;
        }
        double qlen = 0.0; // pre-arrival total population of the destination station
        for (int r = 0; r < R; r++) qlen += classPop(nvec, sm, jnd, r);
        int q = (int) Math.round(qlen);
        double balkProb = 0.0;
        for (BalkingThreshold t : th) {
            if (t.matches(q)) { balkProb = t.getProbability(); break; }
        }
        return balkProb > 0.0 && Maths.rand() < balkProb;
    }

    /**
     * Finite station capacity (loss). An open-class arrival at a physically
     * finite-capacity destination station is lost once the station is full,
     * mirroring the reference producer's hasRoom gate + State.arrivalIsLost in
     * AfterEventStation.handleArv: a physical cap refuses the placement, and an
     * open (njobs == Inf) class self-loops the job away (offered, not carried).
     * Closed classes are NOT dropped here -- a refused closed job must block
     * upstream, not vanish from the conserved population, so they fall through
     * unchanged (the NRM handler's pre-existing behavior). Buffered-station total
     * occupancy (buffer + in service) is capped at sn.cap, per-class at
     * sn.classcap (0 = no per-class bound, MATLAB convention). Inert unless the
     * destination declares a physical drop rule.
     */
    static boolean capacityLoss(NetworkStruct sn, Matrix nvec, int destPos, int R, Smap sm) {
        int jnd = sm.node[destPos];
        int dstC = sm.cls[destPos];
        int ist = (int) sn.nodeToStation.get(jnd);
        if (ist < 0) {
            return false;
        }
        if (!State.isPhysicalCapacity(sn, ist, dstC) || !State.arrivalIsLost(sn, ist, dstC)) {
            return false;
        }
        double capLimit = sn.cap.get(ist);
        double total = 0.0; // pre-arrival total population of the destination station
        for (int r = 0; r < R; r++) total += classPop(nvec, sm, jnd, r);
        if (Double.isFinite(capLimit) && total >= capLimit) {
            return true;
        }
        double classCapLimit = (sn.classcap != null) ? sn.classcap.get(ist, dstC) : 0.0;
        if (classCapLimit > 0 && classPop(nvec, sm, jnd, dstC) >= classCapLimit) {
            return true;
        }
        return false;
    }

    /**
     * Round-robin routing pointers.
     * <p>
     * The pointer that RROBIN/WRROBIN walk is a per-(node,class) local variable,
     * not a population, and no rate depends on it: it only decides where a
     * departure goes. In a generator that makes it a genuine extra state
     * dimension, but a simulator can carry it as auxiliary state alongside the
     * buffers, which is what happens here -- so it must NOT enter the reaction
     * network. State.afterEventRouter advances the pointer on the departure and
     * the routing closure then reads state_AFTER, so the destination used is the
     * one the pointer lands on: advance first, then select.
     * </p>
     */
    static final class Rr {
        boolean on = false;
        boolean[][] isrr;   // node x class
        int[][] cycle;      // ordered destination node list walked per dispatch
        int[] pos;          // current position in that list, flattened node*R + class
    }

    /**
     * Per-(node,class) round-robin pointers, seeded from the initial state so a
     * warm start is honored. RROBIN's slot holds a destination NODE INDEX,
     * WRROBIN's a POSITION in the weighted cycle (each outlink repeated by its
     * weight), matching State.spaceLocalVars and AfterEventRouter.
     */
    static Rr rrPrecompute(NetworkStruct sn) {
        Rr rr = new Rr();
        if (sn.routing == null) {
            return rr;
        }
        int R = sn.nclasses;
        rr.isrr = new boolean[sn.nnodes][R];
        rr.cycle = new int[sn.nnodes * R][];
        rr.pos = new int[sn.nnodes * R];
        for (int ind = 0; ind < sn.nnodes; ind++) {
            Map<JobClass, jline.lang.constant.RoutingStrategy> rmap = sn.routing.get(sn.nodes.get(ind));
            if (rmap == null) continue;
            for (int r = 0; r < R; r++) {
                jline.lang.constant.RoutingStrategy rs = rmap.get(sn.jobclasses.get(r));
                boolean isWRR = rs == jline.lang.constant.RoutingStrategy.WRROBIN;
                if (rs != jline.lang.constant.RoutingStrategy.RROBIN && !isWRR) {
                    continue;
                }
                NodeParam np = sn.nodeparam.get(sn.nodes.get(ind));
                Matrix cyc = null;
                if (isWRR && np.weightedOutlinks != null) {
                    cyc = np.weightedOutlinks.get(sn.jobclasses.get(r));
                }
                if (cyc == null || cyc.length() == 0) {
                    cyc = (np.outlinks != null) ? np.outlinks.get(sn.jobclasses.get(r)) : null;
                }
                if (cyc == null || cyc.length() == 0) {
                    continue;
                }
                rr.on = true;
                rr.isrr[ind][r] = true;
                int len = (int) cyc.length();
                int[] cv = new int[len];
                for (int i = 0; i < len; i++) cv[i] = (int) cyc.get(i);
                rr.cycle[ind * R + r] = cv;
                rr.pos[ind * R + r] = rrSeed(sn, ind, r, cv, isWRR);
            }
        }
        return rr;
    }

    /**
     * Position the pointer of (ind,class r) starts at, read off the initial state
     * slot. The routing variables sit in the trailing sum(nvars(ind,:)) columns
     * of the state row and the (R+r)-th variable is this class's pointer, which
     * is the slot AfterEventRouter advances. Returns 0 (cycle start) when the
     * node carries no state row or the slot holds nothing recognizable.
     */
    private static int rrSeed(NetworkStruct sn, int ind, int r, int[] cyc, boolean isWRR) {
        if (sn.isstateful.get(ind) == 0.0) {
            return 0;
        }
        Matrix st = sn.state.get(sn.stateful.get((int) sn.nodeToStateful.get(ind)));
        if (st == null || st.getNumRows() == 0) {
            return 0;
        }
        int R = sn.nclasses;
        int nvarSum = 0;
        for (int c = 0; c <= R + r; c++) nvarSum += (int) sn.nvars.get(ind, c);
        // nvarSum indexes the VAR BLOCK, which is right-aligned in the state row
        // behind the buffer/server columns, so it has to be offset past them.
        // AfterEventRouter gets away with using nvarSum directly only because it
        // is handed spaceVar already split off from the row; MATLAB's
        // rrPrecompute applies it to the whole row, which coincides with this
        // expression only for a node that has no buffer/server columns (a
        // Router). A station that routes with RROBIN would land inside its own
        // server block there, so the offset is kept explicit here.
        int slot = st.getNumCols() - (int) sn.nvars.sumRows(ind) + nvarSum - 1;
        if (slot < 0 || slot >= st.getNumCols()) {
            return 0;
        }
        int v = (int) st.get(0, slot);
        if (isWRR) {
            return (v >= 1 && v <= cyc.length) ? (v - 1) : 0;
        }
        for (int i = 0; i < cyc.length; i++) {
            if (cyc[i] == v) return i;
        }
        return 0;
    }

    /** Advance the pointer cyclically and return the destination node it lands on. */
    static int rrNext(Rr rr, int ind, int r, int R) {
        int[] cyc = rr.cycle[ind * R + r];
        int p = rr.pos[ind * R + r];
        p = (p >= cyc.length - 1) ? 0 : (p + 1);
        rr.pos[ind * R + r] = p;
        return cyc[p];
    }

    /**
     * G-network signals.
     *
     * A signal class never joins the station it reaches: it removes jobs there
     * and is annihilated (State.afterEventStationSignal). Like balking, this is
     * an arrival-side effect only -- the departure that released the signal
     * keeps its rate, just its arrival outcome differs -- so no propensity is
     * gated on it.
     *
     * The reference enumerates every victim subset with its probability because
     * it builds a generator; a simulator only needs one trajectory, so the batch
     * size and the victims are sampled instead. The two agree in distribution.
     */
    static final class Sig {
        boolean on;
        NetworkStruct sn;
        boolean[] issignal;   // by class
        int[] nonsignal;      // classes eligible as victims of an untargeted signal
    }

    static Sig sigPrecompute(NetworkStruct sn) {
        Sig sig = new Sig();
        if (sn.issignal == null || sn.issignal.isEmpty()) {
            return sig;
        }
        int R = sn.nclasses;
        boolean[] issignal = new boolean[R];
        boolean any = false;
        int nnon = 0;
        for (int r = 0; r < R; r++) {
            issignal[r] = sn.issignal.get(r) != 0.0;
            if (issignal[r]) any = true; else nnon++;
        }
        if (!any) {
            return sig;
        }
        sig.on = true;
        sig.sn = sn;
        sig.issignal = issignal;
        sig.nonsignal = new int[nnon];
        int x = 0;
        for (int r = 0; r < R; r++) {
            if (!issignal[r]) sig.nonsignal[x++] = r;
        }
        return sig;
    }

    /** True when the state slot destPos is a signal class at a station. */
    static boolean sigIsSignalArrival(Sig sig, int destPos, int R, Smap sm) {
        int s = sm.cls[destPos];
        int jnd = sm.node[destPos];
        return sig.issignal[s] && sig.sn.isstation.get(jnd) != 0.0;
    }

    /**
     * True if signal class cls empties a station on arrival. The signaltype is
     * consulted alongside the iscatastrophe flag so that the two encodings of a
     * catastrophe cannot disagree (State.isCatastropheSignal, SolverMAM and
     * Solver_ssj all apply the same test).
     */
    static boolean sigIsCatastrophe(NetworkStruct sn, int cls) {
        if (sn.iscatastrophe != null && !sn.iscatastrophe.isEmpty()
                && sn.iscatastrophe.length() > cls && sn.iscatastrophe.get(cls) != 0.0) {
            return true;
        }
        return sn.signaltype != null && sn.signaltype.size() > cls
                && sn.signaltype.get(cls) == SignalType.CATASTROPHE;
    }

    /**
     * Apply the arrival of a signal class at a station: pick the victims and
     * remove them.
     */
    static void sigApply(Sig sig, Matrix nvec, ArrayDeque<Integer>[] buffers, int destPos,
                         int R, double[] mi, Smap sm) {
        NetworkStruct sn = sig.sn;
        int jnd = sm.node[destPos];
        int cls = sm.cls[destPos];
        int ist = (int) sn.nodeToStation.get(jnd);

        // CATASTROPHE empties the station of every job, ignoring the batch-size
        // distribution: a catastrophe removes all jobs by definition.
        if (sigIsCatastrophe(sn, cls)) {
            for (int r = 0; r < R; r++) {
                for (int k = 0; k < sm.nph[jnd][r]; k++) nvec.set(sm.slot(jnd, r, k), 0, 0.0);
            }
            buffers[jnd].clear();
            return;
        }

        // Eligible victim classes. A signal that declares a target (forJobClass,
        // sn.signaltarget >= 0 since the JAR stores it 0-based) only removes that
        // class; otherwise every non-signal class is eligible, which is the
        // classic Gelenbe negative customer and is what SolverMAM and SolverLDES
        // both do.
        int tgt = -1;
        if (sn.signaltarget != null && !sn.signaltarget.isEmpty() && sn.signaltarget.length() > cls) {
            tgt = (int) sn.signaltarget.get(cls);
        }
        int[] cand = (tgt >= 0) ? new int[]{tgt} : sig.nonsignal;
        int neligible = 0;
        double ntot = 0.0;
        for (int x = 0; x < cand.length; x++) {
            if (classPop(nvec, sm, jnd, cand[x]) > 0.0) {
                neligible++;
                ntot += classPop(nvec, sm, jnd, cand[x]);
            }
        }
        if (neligible == 0 || ntot <= 0.0) {
            return; // no victim: the signal simply vanishes
        }
        int[] tgtclasses = new int[neligible];
        int y = 0;
        for (int x = 0; x < cand.length; x++) {
            if (classPop(nvec, sm, jnd, cand[x]) > 0.0) tgtclasses[y++] = cand[x];
        }

        // Batch size. Without a removal distribution a signal removes exactly one
        // job. With one, the draw is clipped at the eligible population, so an
        // oversized batch empties it rather than driving the queue negative; this
        // is the min(B, n) clipping that Solver_ssj and SolverMAM also apply, and
        // is distributionally the pmf that the reference enumerates (the tail
        // P(B >= ntot) lumps onto "remove ntot").
        int k = 1;
        DiscreteDistribution dist = (sn.signalremdist != null && sn.signalremdist.size() > cls)
                ? sn.signalremdist.get(cls) : null;
        if (dist != null) {
            int sampled = (int) dist.sample(1)[0];
            k = (int) Math.min(sampled, ntot);
            if (k < 0) k = 0;
        }

        RemovalPolicy policy = RemovalPolicy.RANDOM;
        if (sn.signalrempolicy != null && sn.signalrempolicy.size() > cls
                && sn.signalrempolicy.get(cls) != null) {
            policy = sn.signalrempolicy.get(cls);
        }

        for (int step = 0; step < k; step++) {
            if (!sigRemoveOne(sn, nvec, buffers, jnd, ist, tgtclasses, policy, R, mi, sm)) {
                break; // already drained
            }
        }
    }

    /**
     * Remove one class-r job from node jnd, drawing its service phase in
     * proportion to the phase occupancies. The victim policies below select a
     * JOB, not a phase, and the class-r jobs in service are exchangeable given
     * the phase counts, so the phase of the victim is drawn from those counts. A
     * single-phase class is a plain decrement and consumes no random number, so
     * an exponential model keeps its previous random stream.
     */
    static void removeOneJob(Matrix nvec, Smap sm, int jnd, int r) {
        int np = sm.nph[jnd][r];
        int off = sm.phOff[jnd][r];
        if (np <= 1) {
            nvec.set(off, 0, nvec.get(off, 0) - 1.0);
            return;
        }
        double tot = classPop(nvec, sm, jnd, r);
        if (tot <= 0.0) {
            nvec.set(off, 0, nvec.get(off, 0) - 1.0);
            return;
        }
        double u = Maths.rand() * tot;
        double acc = 0.0;
        for (int k = 0; k < np; k++) {
            acc += nvec.get(off + k, 0);
            if (u < acc) {
                nvec.set(off + k, 0, nvec.get(off + k, 0) - 1.0);
                return;
            }
        }
        nvec.set(off + np - 1, 0, nvec.get(off + np - 1, 0) - 1.0);
    }

    /**
     * Remove one victim under the signal's removal policy. Waiting jobs live in
     * the buffer; the rest of each class population is in service. Returns false
     * when nothing eligible is left.
     */
    static boolean sigRemoveOne(NetworkStruct sn, Matrix nvec, ArrayDeque<Integer>[] buffers,
                                int jnd, int ist, int[] tgtclasses,
                                RemovalPolicy policy, int R, double[] mi, Smap sm) {
        Integer[] buf = buffers[jnd].toArray(new Integer[0]);
        List<Integer> waitIdx = new ArrayList<Integer>();   // eligible waiting positions
        for (int i = 0; i < buf.length; i++) {
            for (int x = 0; x < tgtclasses.length; x++) {
                if (buf[i] - 1 == tgtclasses[x]) { waitIdx.add(i); break; }
            }
        }
        int nwait = waitIdx.size();
        double nsrv = 0.0;
        for (int x = 0; x < tgtclasses.length; x++) {
            nsrv += Math.max(0.0, classPop(nvec, sm, jnd, tgtclasses[x]) - bufCount(buf, tgtclasses[x]));
        }
        if (nwait == 0 && nsrv == 0.0) {
            return false;
        }

        // FCFS/LCFS rank the waiting line by age, which only an ordered buffer
        // records. The NRM buffer is newest-first / oldest-last, so the head of
        // line (oldest) is the last eligible position and the most recent arrival
        // the first. A per-class count buffer (SIRO/SEPT/LEPT) carries no age, so
        // an age-based policy degenerates to a uniform draw there, exactly as in
        // the reference.
        SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
        boolean isOrdered = sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL
                || sched == SchedStrategy.LCFS;
        boolean ageOrdered = isOrdered
                && (policy == RemovalPolicy.FCFS || policy == RemovalPolicy.LCFS);
        if (ageOrdered && nwait > 0) {
            int pick = (policy == RemovalPolicy.FCFS)
                    ? waitIdx.get(nwait - 1)   // head of line: the oldest waiting job
                    : waitIdx.get(0);          // the most recent arrival
            int victim = buf[pick] - 1;
            removeAt(buffers[jnd], pick);
            removeOneJob(nvec, sm, jnd, victim);
            return true;
        }

        // RANDOM draws uniformly over waiting and in-service alike; FCFS/LCFS
        // drain the waiting line before reaching into the servers.
        double total;
        if (policy == RemovalPolicy.RANDOM) {
            total = nwait + nsrv;
        } else {
            total = nwait;
            if (total == 0.0) total = nsrv;
        }
        double u = Maths.rand() * total;
        if (nwait > 0 && (policy != RemovalPolicy.RANDOM || u < nwait)) {
            // a waiting victim, uniform over the eligible positions
            int pick = waitIdx.get((int) Math.floor(Maths.rand() * nwait));
            int victim = buf[pick] - 1;
            removeAt(buffers[jnd], pick);
            removeOneJob(nvec, sm, jnd, victim);
            return true;
        }

        // an in-service victim, uniform over the eligible in-service jobs
        double acc = 0.0;
        double target = Maths.rand() * nsrv;
        for (int x = 0; x < tgtclasses.length; x++) {
            int r = tgtclasses[x];
            double cnt = Math.max(0.0, classPop(nvec, sm, jnd, r) - bufCount(buf, r));
            acc += cnt;
            if (cnt > 0.0 && target < acc) {
                removeOneJob(nvec, sm, jnd, r);
                // the freed server pulls the head of line in, which in the NRM is
                // just the waiting job leaving the buffer (in-service is derived
                // as population minus buffer occupancy)
                double totalNew = 0.0;
                for (int s = 0; s < R; s++) totalNew += classPop(nvec, sm, jnd, s);
                if (buffers[jnd].size() > Math.max(0.0, totalNew - mi[jnd])) {
                    buffers[jnd].pollLast();   // head of line: the oldest waiting job
                }
                return true;
            }
        }
        return false;
    }

    /** Occupancy of class r (0-based) in a newest-first buffer snapshot. */
    private static int bufCount(Integer[] buf, int r) {
        int n = 0;
        for (int i = 0; i < buf.length; i++) {
            if (buf[i] - 1 == r) n++;
        }
        return n;
    }

    /**
     * Abandonment rate of a single waiting class-r job at station ist
     * (sn.impatienceMu), or 0 when the station-class does not renege. Only
     * exponential patience reaches here: Solver_ssa_analyzer routes phase-type
     * patience to the serial engine, since the aggregate rate law below is
     * correct only when patience is memoryless.
     */
    static double renegeRateOf(NetworkStruct sn, int ist, int r) {
        if (sn.impatienceClass == null || sn.impatienceMu == null) {
            return 0.0;
        }
        Map<JobClass, ImpatienceType> cmap = sn.impatienceClass.get(sn.stations.get(ist));
        if (cmap == null || cmap.get(sn.jobclasses.get(r)) != ImpatienceType.RENEGING) {
            return 0.0;
        }
        Map<JobClass, Matrix> mumap = sn.impatienceMu.get(sn.stations.get(ist));
        if (mumap == null) {
            return 0.0;
        }
        Matrix mu = mumap.get(sn.jobclasses.get(r));
        if (mu == null || mu.length() == 0) {
            return 0.0;
        }
        return mu.get(0);
    }

    /**
     * Normalized DPS/GPS-family weights of every station, indexed by station,
     * with a null entry for stations whose policy does not read weights.
     * State.afterEventStation rejects multi-server DPS/GPS, so this fails on
     * those rather than silently simulating a different station.
     */
    static double[][] buildSchedWeights(NetworkStruct sn, int R) {
        double[][] wnorm = new double[sn.nstations][];
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.DPS || s == SchedStrategy.GPS
                    || s == SchedStrategy.DPSPRIO || s == SchedStrategy.GPSPRIO) {
                wnorm[ist] = schedWeights(sn, ist, R);
                if (sn.nservers.get(ist) > 1) {
                    throw new RuntimeException("Multi-server " + s + " stations are not supported yet.");
                }
            }
        }
        return wnorm;
    }

    /**
     * True for the processor-sharing family, whose utilization is the share of
     * service capacity a class receives rather than a server-occupancy count.
     */
    static boolean isPsFamily(SchedStrategy sched) {
        return sched == SchedStrategy.PS || sched == SchedStrategy.LPS
                || sched == SchedStrategy.DPS || sched == SchedStrategy.GPS
                || sched == SchedStrategy.PSPRIO || sched == SchedStrategy.DPSPRIO
                || sched == SchedStrategy.GPSPRIO;
    }

    /**
     * Class-k utilization of a PS-family station at population npop: the share
     * of service capacity class k receives, divided by the server count. The
     * lld/cd scalings are excluded on purpose, as they rescale work rather than
     * server occupancy, matching the PS branch of the MATLAB solver_ssa_nrm.
     */
    static double psFamilyUtil(NetworkStruct sn, double[] w, double[] npop, int k,
                               double servers, SchedStrategy sched) {
        if (sched == SchedStrategy.PS || sched == SchedStrategy.LPS) {
            double totalPop = sumPop(npop);
            if (totalPop <= 0.0) {
                return 0.0;
            }
            return (npop[k] / totalPop) * Math.min(servers, totalPop) / servers;
        }
        if (sched == SchedStrategy.DPS) {
            return dpsshare(w, npop, k) / servers;
        }
        if (sched == SchedStrategy.GPS) {
            return gpsshare(w, npop, k) / servers;
        }
        if (sched == SchedStrategy.PSPRIO) {
            return psprioshare(sn, npop, k, servers) / servers;
        }
        if (sched == SchedStrategy.DPSPRIO) {
            return dpsprioshare(sn, w, npop, k, servers) / servers;
        }
        if (sched == SchedStrategy.GPSPRIO) {
            return gpsprioshare(sn, w, npop, k, servers) / servers;
        }
        return 0.0;
    }

    /**
     * GPSPRIO: GPS below capacity, GPS restricted to the urgent group above it.
     */
    static double gpsprioshare(NetworkStruct sn, double[] w, double[] npop, int r, double c) {
        double ni = sumPop(npop);
        if (ni <= 0.0) {
            return 0.0;
        }
        if (ni <= c) {
            return gpsshare(w, npop, r);
        }
        if (!isUrgent(sn, npop, r)) {
            return 0.0;
        }
        return gpsshare(w, prioGroup(sn, npop, r), r);
    }

    /**
     * Limited load-dependent scaling factor at total station population ntot:
     * the tabulated sn.lldscaling entry of station ist (as in
     * State.afterEventStation and the MATLAB solver_ssa_nrm lldfac / python
     * nrm _lldfac helpers). Returns 1 when the station has no load dependence
     * or is empty; beyond the tabulated limit the last entry is used.
     */
    static double lldfac(NetworkStruct sn, int ist, double ntot) {
        if (sn.lldscaling == null || sn.lldscaling.isEmpty() || ntot < 1.0) {
            return 1.0;
        }
        int limit = sn.lldscaling.getNumCols();
        int idx = (int) Math.min(Math.round(ntot), (long) limit) - 1;
        if (ist >= 0 && ist < sn.lldscaling.getNumRows() && idx >= 0 && idx < limit) {
            return sn.lldscaling.get(ist, idx);
        }
        return 1.0;
    }

    public static SolverSSAResultNRM solver_ssa_nrm(NetworkStruct sn_in, SolverOptions options) {
        RandomManager.setMasterSeed(options.seed);

        // Immediate feedback (sn.immfeed) is not modeled by the NRM reaction
        // network: self-loops are handled as ordinary class-switching with
        // re-queueing. Warn so an immfeed model is not silently given the
        // re-queueing result; the serial SSA engine (Solver_ssa) honors it.
        if (sn_in.immfeed != null && sn_in.immfeed.elementSum() > 0) {
            line_warning(mfilename(new Object() {
            }), "SolverSSA (method=nrm) does not model immediate feedback (immfeed); self-loops are treated as class-switching with re-queueing. Use method=serial for immediate feedback.");
        }

        final NetworkStruct sn = sn_in.copy();
        final int samples = options.samples;
        final int R = sn.nclasses;
        final int I = sn.nnodes;

        // The state vector counts jobs per (node, class, PHASE) rather than per
        // (node, class), so that phase-type service is represented exactly instead
        // of being collapsed onto its mean rate. With nph == 1 everywhere the slot
        // map is the old flat (node,class) layout, so an exponential model is
        // unaffected -- buildSmap asserts that reduction.
        final Smap sm = buildSmap(sn, I, R);
        final int NS = sm.NS;

        // Buffered phase-type service: see _kb/06-solver-catalog.md (SSA/NRM
        // section) for the svcph auxiliary-structure rationale.
        final Set<SchedStrategy> bufPHSchedSet = new java.util.HashSet<SchedStrategy>(
                java.util.Arrays.asList(SchedStrategy.FCFS, SchedStrategy.LCFS,
                        SchedStrategy.SIRO, SchedStrategy.HOL, SchedStrategy.SEPT,
                        SchedStrategy.LEPT));
        final boolean[][] bufPHClass = new boolean[I][R];
        final boolean[] bufPHNode = new boolean[I];
        int maxnphTmp = 1;
        for (int ind = 0; ind < I; ind++) {
            boolean isBufSt = false;
            if (sn.isstation.get(ind) != 0.0) {
                int istB = (int) sn.nodeToStation.get(ind);
                if (istB >= 0) {
                    isBufSt = bufPHSchedSet.contains(sn.sched.get(sn.stations.get(istB)));
                }
            }
            for (int r = 0; r < R; r++) {
                if (sm.nph[ind][r] > maxnphTmp) maxnphTmp = sm.nph[ind][r];
                if (isBufSt && sm.nph[ind][r] > 1) {
                    bufPHClass[ind][r] = true;
                    bufPHNode[ind] = true;
                }
            }
        }
        final int maxnph = maxnphTmp;
        // In-service phase multiset of each buffered-PH node; null for every other
        // node. Captured by the propensity closures (which read it in place) and
        // by the engine (which mutates it on arrivals, departures and phase moves).
        final double[][][] svcph = new double[I][][];
        for (int ind = 0; ind < I; ind++) {
            if (bufPHNode[ind]) {
                svcph[ind] = new double[R][maxnph];
            }
        }

        // Cache nodes. A Cache is an immediate class switch: a job arrives in a
        // READ class, reads an item drawn from pread, and leaves in the hit or
        // miss class depending on whether the item is cached, after which the
        // replacement policy rewrites the contents. The NRM models this as a
        // state-dependent class-switch reaction at the cache node (consume
        // [cache,readClass], produce [cache,hitClass] or [cache,missClass],
        // chosen at firing by the cache access), mirroring State.afterEventCache.
        // The cache CONTENTS ride alongside the engine in cacheContents0 (like
        // pollCtrl/svcph); the hit/miss draw and replacement read and rewrite them
        // at firing. A class r is a READ class of cache ind iff its pread entry is
        // a non-empty probability row and its hitclass is set.
        final boolean[] isCacheNode = new boolean[I];
        final boolean[][] isCacheReadClass = new boolean[I][R];
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Cache) {
                isCacheNode[ind] = true;
                jline.lang.NodeParam npc = sn.nodeparam.get(sn.nodes.get(ind));
                if (npc instanceof jline.lang.nodeparam.CacheNodeParam) {
                    jline.lang.nodeparam.CacheNodeParam cnp = (jline.lang.nodeparam.CacheNodeParam) npc;
                    for (int r = 0; r < R; r++) {
                        boolean readr = cnp.pread != null && cnp.pread.containsKey(r)
                                && cnp.pread.get(r) != null && !cnp.pread.get(r).isEmpty();
                        boolean hasHit = cnp.hitclass != null && r < cnp.hitclass.length()
                                && cnp.hitclass.get(r) > 0;
                        // A retrieval class (created by setRetrievalSystem) reads its
                        // own one-hot item to COMPLETE a miss; it has hitclass 0 but a
                        // miss class, so it is a read class too.
                        boolean isRetr = cnp.retrievalClassIndices != null
                                && cnp.retrievalClassIndices.contains(r);
                        if (readr && (hasHit || isRetr)) isCacheReadClass[ind][r] = true;
                    }
                }
            }
        }

        // Stochastic Petri net path: a model with Transition nodes is a Petri net,
        // not a queueing network. Its dynamics are firings of transition modes over
        // a place marking, not job departures routed by the rt matrix, so the
        // generic (node,class) departure grid below does not apply (a firing
        // produces to several places deterministically, never a routing draw).
        // Route it to the dedicated builder/runner, which shares the Gibson & Bruck
        // clocks but its own firing application and vanishing-marking collapse.
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Transition) {
                return nrm_spn(sn, options, sm, samples);
            }
        }

        // Polling controller descriptions, one per polling station. The NRM tracks
        // only the controller [mode,pos,swk,ctr] of a polling server, not the
        // service phase of the single job in service, so phase-type service at a
        // polling station is rejected here (the switchover itself may be
        // phase-type: its phase rides in the controller, not in nvec).
        final Polling.Info[] pinfo = new Polling.Info[I];
        final boolean[] isPollNode = new boolean[I];
        boolean anyPollTmp = false;
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstation.get(ind) == 0.0) continue;
            int istP = (int) sn.nodeToStation.get(ind);
            if (sn.sched.get(sn.stations.get(istP)) != SchedStrategy.POLLING) continue;
            pinfo[ind] = Polling.info(sn, ind);
            isPollNode[ind] = true;
            anyPollTmp = true;
            for (int r = 0; r < R; r++) {
                if (sm.nph[ind][r] > 1) {
                    throw new RuntimeException("NRM polling supports exponential service only; station "
                            + istP + " class " + r + " has phase-type service. Use method=serial.");
                }
            }
        }
        final boolean anyPoll = anyPollTmp;
        // Shared polling controller store, captured by the propensity closures and
        // mutated in place by the direct engine; null at every non-polling node.
        final int[][] pollCtrl = new int[I][];

        Matrix S = new Matrix(0, NS);
        final List<Integer> fromIdx = new ArrayList<Integer>();
        final List<List<Integer>> toIdx = new ArrayList<List<Integer>>();
        final List<int[]> fromIR = new ArrayList<int[]>();
        final List<List<Double>> probIR = new ArrayList<List<Double>>();
        final List<Integer> depPhase = new ArrayList<Integer>();   // source phase each dep/phase reaction reads
        final List<Boolean> phaseFlag = new ArrayList<Boolean>();  // true for phase-transition reactions
        final List<Double> phaseRate = new ArrayList<Double>();    // D0(k,k') of a phase transition
        final List<Integer> phaseTo = new ArrayList<Integer>();    // target phase kb of a phase transition (-1 if none)
        final List<Boolean> bufSvcFlag = new ArrayList<Boolean>(); // dep reaction of a buffered-PH class (reads svcph)
        final List<Boolean> cacheFlag = new ArrayList<Boolean>();  // cache-access reaction (read -> hit/miss at a cache node)
        final List<Integer> cacheHitSlot = new ArrayList<Integer>();  // nvec slot the hit-class job is produced into
        final List<Integer> cacheMissSlot = new ArrayList<Integer>(); // nvec slot the miss-class job is produced into
        // Policies whose waiting jobs (PAS/OI: whose whole job list) are carried
        // in a per-node buffer alongside the reaction network, and whose initial
        // contents therefore have to be read out of the initial state.
        final Set<SchedStrategy> bufferedSched = new java.util.HashSet<SchedStrategy>(
                java.util.Arrays.asList(SchedStrategy.FCFS, SchedStrategy.LCFS,
                        SchedStrategy.SIRO, SchedStrategy.HOL, SchedStrategy.SEPT,
                        SchedStrategy.LEPT, SchedStrategy.LCFSPR, SchedStrategy.PAS));

        // Departure reactions, one per (node, class, PHASE). A departure is the
        // absorption of the phase-type service process, so it fires at mu(k)*phi(k)
        // and the job re-enters its destination in an entry phase drawn from pie:
        // the destination draw therefore carries the product of the routing
        // probability and the entry-phase probability. The weighted-destination
        // sampler takes that product unchanged.
        int k = 0;
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) {
                for (int ph = 0; ph < sm.nph[ind][r]; ph++) {
                    k++;
                    fromIR.add(new int[]{ind, r});
                    if (isCacheReadClass[ind][r]) {
                        // Cache access: consume the read-class job at the cache; its
                        // production (hit or miss class, at the SAME cache node) and
                        // the contents update are resolved at firing by cacheAccess.
                        // No static routing: rtnodes has no out-edge for the read
                        // class. Caches have nph==1, so this fires once (ph==0).
                        jline.lang.nodeparam.CacheNodeParam cnp =
                                (jline.lang.nodeparam.CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                        fromIdx.add(sm.slot(ind, r, 0));
                        depPhase.add(0);
                        phaseFlag.add(Boolean.FALSE);
                        phaseRate.add(0.0);
                        phaseTo.add(-1);
                        bufSvcFlag.add(Boolean.FALSE);
                        cacheFlag.add(Boolean.TRUE);
                        // The outcome class (hit/miss/retrieval) is resolved at
                        // firing, so these slots are informational only; a retrieval
                        // class has hitclass 0 (and hitclass may not be sized to R),
                        // so guard the lookup.
                        int hcR = (r < cnp.hitclass.length() && cnp.hitclass.get(r) > 0)
                                ? sm.slot(ind, (int) cnp.hitclass.get(r), 0) : -1;
                        int mcR = (r < cnp.missclass.length() && cnp.missclass.get(r) > 0)
                                ? sm.slot(ind, (int) cnp.missclass.get(r), 0) : -1;
                        cacheHitSlot.add(hcR);
                        cacheMissSlot.add(mcR);
                        probIR.add(new ArrayList<Double>());
                        toIdx.add(new ArrayList<Integer>());
                        double[] SrowC = new double[NS];
                        SrowC[sm.slot(ind, r, 0)] = -1.0;
                        if (k > S.getNumRows()) {
                            Matrix newS = new Matrix(k, NS);
                            for (int i = 0; i < S.getNumRows(); i++)
                                for (int j = 0; j < S.getNumCols(); j++) newS.set(i, j, S.get(i, j));
                            S = newS;
                        }
                        for (int j = 0; j < SrowC.length; j++) S.set(k - 1, j, SrowC[j]);
                        continue;
                    }
                    cacheFlag.add(Boolean.FALSE);
                    cacheHitSlot.add(-1);
                    cacheMissSlot.add(-1);
                    // At a buffered-PH source only the jobs in service carry a
                    // phase, tracked in svcph; nvec holds the whole class population
                    // in its first phase slot. A departure therefore removes one job
                    // from that total slot regardless of which service phase
                    // completed -- the completing phase ph is carried in depPhase and
                    // consumed from svcph at firing.
                    fromIdx.add(bufPHClass[ind][r] ? sm.slot(ind, r, 0) : sm.slot(ind, r, ph));
                    depPhase.add(ph);
                    phaseFlag.add(Boolean.FALSE);
                    phaseRate.add(0.0);
                    phaseTo.add(-1);
                    bufSvcFlag.add(bufPHClass[ind][r] ? Boolean.TRUE : Boolean.FALSE);
                    probIR.add(new ArrayList<Double>());
                    toIdx.add(new ArrayList<Integer>());
                    double[] Srow = new double[NS];
                    if (sn.isslc.get(r) != 0.0) {
                        Srow[fromIdx.get(k - 1)] = GlobalConstants.NegInf;
                    } else {
                        Srow[fromIdx.get(k - 1)] = -1.0;
                        for (int jnd = 0; jnd < I; jnd++) {
                            for (int s = 0; s < R; s++) {
                                double p = sn.rtnodes.get(ind * R + r, jnd * R + s);
                                if (p > 0) {
                                    if (bufPHClass[jnd][s]) {
                                        // A job arriving at a buffered-PH destination
                                        // lands in the total-population slot; whether
                                        // it enters service (and its entry phase) or
                                        // waits is decided at firing from the server
                                        // occupancy and pie, not the routing draw. So
                                        // the destination collapses to the single
                                        // total slot with weight p.
                                        int dslot = sm.slot(jnd, s, 0);
                                        toIdx.get(k - 1).add(dslot);
                                        probIR.get(k - 1).add(p);
                                        Srow[dslot] = Srow[dslot] + p;
                                    } else {
                                        double[] pentry = entryProbs(sn, sm, jnd, s);
                                        for (int ke = 0; ke < sm.nph[jnd][s]; ke++) {
                                            if (pentry[ke] <= 0.0) continue;
                                            int dslot = sm.slot(jnd, s, ke);
                                            toIdx.get(k - 1).add(dslot);
                                            probIR.get(k - 1).add(p * pentry[ke]);
                                            Srow[dslot] = Srow[dslot] + p * pentry[ke];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (k > S.getNumRows()) {
                        Matrix newS = new Matrix(k, NS);
                        for (int i = 0; i < S.getNumRows(); i++) {
                            for (int j = 0; j < S.getNumCols(); j++) {
                                newS.set(i, j, S.get(i, j));
                            }
                        }
                        S = newS;
                    }
                    for (int j = 0; j < Srow.length; j++) {
                        S.set(k - 1, j, Srow[j]);
                    }
                }
            }
        }
        final int nDepRx = k;   // departure reactions occupy 0..nDepRx-1

        // Phase-transition reactions, one per (node, class, ka -> kb) with
        // D0(ka,kb) > 0. These move a job between the phases of its own service
        // process and so never leave the node; D0's off-diagonal carries their
        // rates (State.afterEventStation, EventType.PHASE).
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstation.get(ind) == 0.0) continue;
            int ist = (int) sn.nodeToStation.get(ind);
            for (int r = 0; r < R; r++) {
                if (sm.nph[ind][r] <= 1) continue;
                MatrixCell proc = phEntryOf(sn.proc, sn, ist, r);
                if (proc == null || proc.isEmpty()) continue;
                Matrix D0 = proc.get(0);
                if (D0 == null) continue;
                for (int ka = 0; ka < sm.nph[ind][r]; ka++) {
                    for (int kb = 0; kb < sm.nph[ind][r]; kb++) {
                        if (ka == kb || D0.get(ka, kb) <= 0.0) continue;
                        k++;
                        fromIR.add(new int[]{ind, r});
                        fromIdx.add(sm.slot(ind, r, ka));
                        depPhase.add(ka);
                        phaseFlag.add(Boolean.TRUE);
                        phaseRate.add(D0.get(ka, kb));
                        phaseTo.add(kb);
                        bufSvcFlag.add(Boolean.FALSE);
                        probIR.add(new ArrayList<Double>());
                        toIdx.add(new ArrayList<Integer>());
                        if (k > S.getNumRows()) {
                            Matrix newS = new Matrix(k, NS);
                            for (int i = 0; i < S.getNumRows(); i++) {
                                for (int j = 0; j < S.getNumCols(); j++) {
                                    newS.set(i, j, S.get(i, j));
                                }
                            }
                            S = newS;
                        }
                        // A buffered-PH class keeps its in-service phase counts in
                        // svcph, not in nvec: a phase transition moves a job between
                        // phases of the same in-service composition, so it leaves
                        // nvec (the class total) unchanged. The stoichiometry column
                        // is therefore all zeros; the move is applied to svcph at
                        // firing and, like a retry/switchover, its dependency set
                        // must be supplied through a forced refresh (D cannot derive
                        // it from an all-zero column).
                        if (!bufPHClass[ind][r]) {
                            S.set(k - 1, sm.slot(ind, r, ka), -1.0);
                            S.set(k - 1, sm.slot(ind, r, kb), 1.0);
                        }
                    }
                }
            }
        }

        // Reneging: each waiting (queued, not-in-service) class-r job abandons at
        // the memoryless rate sn.impatienceMu, so the aggregate rate out of the
        // state is (waiting count)*mu and the job leaves the system. This is a
        // reaction the (node,class) departure grid above cannot express -- it
        // consumes a job without producing one -- so it is appended as an extra
        // column whose stoichiometry is a bare -1 at the source slot. A renege is
        // not a departure and must not count towards throughput; the TN
        // accumulator sums the per-phase DEPARTURE reactions only, so the appended
        // columns and the phase-transition columns both stay out of it.
        final int nDep = k;   // departure and phase reactions occupy 0..nDep-1
        final List<Double> renegeMu = new ArrayList<Double>();
        final List<Boolean> renegeFlag = new ArrayList<Boolean>();
        for (int j = 0; j < nDep; j++) {
            renegeMu.add(0.0);
            renegeFlag.add(Boolean.FALSE);
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int r = 0; r < R; r++) {
                double mu = renegeRateOf(sn, ist, r);
                if (mu <= 0.0) continue;
                k++;
                fromIR.add(new int[]{ind, r});
                // Reneging needs memoryless patience and so only reaches a
                // single-phase class, whose sole slot is the class's first one.
                fromIdx.add(sm.slot(ind, r, 0));
                depPhase.add(0);
                phaseTo.add(-1);
                bufSvcFlag.add(Boolean.FALSE);
phaseFlag.add(Boolean.FALSE);
                phaseRate.add(0.0);
                probIR.add(new ArrayList<Double>());
                toIdx.add(new ArrayList<Integer>());
                if (k > S.getNumRows()) {
                    Matrix newS = new Matrix(k, NS);
                    for (int i = 0; i < S.getNumRows(); i++) {
                        for (int j = 0; j < S.getNumCols(); j++) {
                            newS.set(i, j, S.get(i, j));
                        }
                    }
                    S = newS;
                }
                S.set(k - 1, sm.slot(ind, r, 0), -1.0);   // job abandons and leaves the system
                renegeMu.add(mu);
                renegeFlag.add(Boolean.TRUE);
            }
        }
        // Retrial: an orbiting class-r job retries entry at the memoryless rate
        // sn.retrialMu and succeeds only when a server is free; otherwise the
        // event is a no-op and is simply not generated. A retry moves a job from
        // the orbit into service WITHOUT changing any population, so its
        // stoichiometry column is ALL ZEROS -- which is why its dependency set,
        // and its membership of every other reaction's, both have to be supplied
        // by hand below: D is derived from the SIGN of the S entries, so an
        // all-zero column would otherwise leave every rate at the node stale
        // after a retry fires and never appear in any other reaction's set.
        final List<Double> retryMu = new ArrayList<Double>();
        final List<Boolean> retryFlag = new ArrayList<Boolean>();
        for (int j = 0; j < k; j++) {
            retryMu.add(0.0);
            retryFlag.add(Boolean.FALSE);
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            if (!isRetrialStation(ind, sn)) continue;
            for (int r = 0; r < R; r++) {
                double mu = retrialRateOf(sn, ist, r);
                if (mu <= 0.0) continue;
                k++;
                fromIR.add(new int[]{ind, r});
                // A retrial station buffers its orbit, so it only reaches a
                // single-phase class: the class's first slot is its sole slot.
                fromIdx.add(sm.slot(ind, r, 0));
                depPhase.add(0);
                phaseTo.add(-1);
                bufSvcFlag.add(Boolean.FALSE);
phaseFlag.add(Boolean.FALSE);
                phaseRate.add(0.0);
                probIR.add(new ArrayList<Double>());
                toIdx.add(new ArrayList<Integer>());
                if (k > S.getNumRows()) {
                    Matrix newS = new Matrix(k, NS);
                    for (int i = 0; i < S.getNumRows(); i++) {
                        for (int j = 0; j < S.getNumCols(); j++) {
                            newS.set(i, j, S.get(i, j));
                        }
                    }
                    S = newS;
                }
                // row k-1 stays all zero: a retry moves no job between nodes
                retryMu.add(mu);
                retryFlag.add(Boolean.TRUE);
                renegeMu.add(0.0);
                renegeFlag.add(Boolean.FALSE);
            }
        }
        // Polling switchover reactions: one per polling station that has a timed
        // switchover. A server dwelling in a switchover is a genuine timed event
        // that moves no job, so it is appended as an all-zero stoichiometry column
        // whose propensity reads the controller and whose firing samples the
        // switchover PH -- exactly as a retry is appended. A station whose
        // switchovers are all immediate never dwells and gets no reaction.
        final List<Integer> pollSwNodeList = new ArrayList<Integer>();
        final List<Integer> pollSwRxIdx = new ArrayList<Integer>();
        for (int ind = 0; ind < I; ind++) {
            if (!isPollNode[ind]) continue;
            boolean hasAnySw = false;
            for (int r = 0; r < R; r++) {
                if (pinfo[ind].hasSw[r]) { hasAnySw = true; break; }
            }
            if (!hasAnySw) continue;
            k++;
            fromIR.add(new int[]{ind, 0});    // class field is a sentinel, never read as a class
            fromIdx.add(sm.slot(ind, 0, 0));  // unused slot: a switchover consumes no job
            depPhase.add(0);
            phaseTo.add(-1);
            bufSvcFlag.add(Boolean.FALSE);
phaseFlag.add(Boolean.FALSE);
            phaseRate.add(0.0);
            probIR.add(new ArrayList<Double>());
            toIdx.add(new ArrayList<Integer>());
            if (k > S.getNumRows()) {
                Matrix newS = new Matrix(k, NS);
                for (int i = 0; i < S.getNumRows(); i++) {
                    for (int j = 0; j < S.getNumCols(); j++) {
                        newS.set(i, j, S.get(i, j));
                    }
                }
                S = newS;
            }
            // row k-1 stays all zero: a switchover moves no job between nodes
            renegeMu.add(0.0);
            renegeFlag.add(Boolean.FALSE);
            retryMu.add(0.0);
            retryFlag.add(Boolean.FALSE);
            pollSwNodeList.add(ind);
            pollSwRxIdx.add(k - 1);
        }
        final boolean[] isPollSwRx = new boolean[fromIdx.size()];
        final int[] pollSwNode = new int[fromIdx.size()];
        for (int j = 0; j < pollSwNode.length; j++) pollSwNode[j] = -1;
        for (int t = 0; t < pollSwRxIdx.size(); t++) {
            isPollSwRx[pollSwRxIdx.get(t)] = true;
            pollSwNode[pollSwRxIdx.get(t)] = pollSwNodeList.get(t);
        }

        final boolean[] isRenegeRx = new boolean[renegeFlag.size()];
        for (int j = 0; j < isRenegeRx.length; j++) {
            isRenegeRx[j] = renegeFlag.get(j);
        }
        final boolean[] isRetryRx = new boolean[retryFlag.size()];
        for (int j = 0; j < isRetryRx.length; j++) {
            isRetryRx[j] = retryFlag.get(j);
        }
        final boolean[] isPhaseRx = new boolean[phaseFlag.size()];
        for (int j = 0; j < isPhaseRx.length; j++) {
            isPhaseRx[j] = phaseFlag.get(j);
        }
        final boolean[] isBufSvcRx = new boolean[bufSvcFlag.size()];
        for (int j = 0; j < isBufSvcRx.length; j++) {
            isBufSvcRx[j] = bufSvcFlag.get(j);
        }
        // Cache-reaction flags/slots: only the departure reactions add to these
        // lists (via the cache branch above); the phase/renege/retry/switchover
        // reactions appended afterwards do not, so pad to the final count.
        while (cacheFlag.size() < fromIdx.size()) { cacheFlag.add(Boolean.FALSE); cacheHitSlot.add(-1); cacheMissSlot.add(-1); }
        final boolean[] isCacheRx = new boolean[cacheFlag.size()];
        final int[] cacheHitSlotArr = new int[cacheFlag.size()];
        final int[] cacheMissSlotArr = new int[cacheFlag.size()];
        for (int j = 0; j < isCacheRx.length; j++) {
            isCacheRx[j] = cacheFlag.get(j);
            cacheHitSlotArr[j] = cacheHitSlot.get(j);
            cacheMissSlotArr[j] = cacheMissSlot.get(j);
        }
        // Cache contents: item held in each of the totalCacheCapacity slots, laid
        // out list by list, initialised to items 1..totalCacheCapacity (any valid
        // ordered placement is a correct warm start since the chain is ergodic).
        final int[][] cacheContents0 = new int[I][];
        for (int ind = 0; ind < I; ind++) {
            if (isCacheNode[ind]) {
                jline.lang.nodeparam.CacheNodeParam cnp =
                        (jline.lang.nodeparam.CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                int tcc = cnp.totalCacheCapacity;
                // With a retrieval system the contents are followed by a per-item
                // occupancy bitmap: column tcc+item-1 is 1 iff item is being fetched.
                int width = (cnp.retrievalSystemCapacity > 0) ? tcc + cnp.nitems : tcc;
                cacheContents0[ind] = new int[width];
                for (int c = 0; c < tcc; c++) cacheContents0[ind][c] = c + 1;
            }
        }
        // Destination slot a BEGUN retrieval routes to (the fetch queue). A begin
        // must place the job at the retrieval class's routed destination, NOT at
        // [cache,retrievalClass] -- otherwise the cache-access reaction fires again
        // and completes the miss instantly, collapsing the fetch. Only queue-returns
        // occupy [cache,retrievalClass] and trigger completion. Resolved from rtnodes.
        final int[][] cacheRetrDest = new int[I][R];
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) cacheRetrDest[ind][r] = -1;
            if (isCacheNode[ind]) {
                jline.lang.nodeparam.CacheNodeParam cnp =
                        (jline.lang.nodeparam.CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                if (cnp.retrievalClassIndices != null) {
                    for (Integer rc : cnp.retrievalClassIndices) {
                        for (int b = 0; b < I * R; b++) {
                            if (sn.rtnodes.get(ind * R + rc, b) > 0) {
                                int jnd = b / R; int s = b % R;
                                cacheRetrDest[ind][rc] = sm.slot(jnd, s, 0);
                                break;
                            }
                        }
                    }
                }
            }
        }
        final int[] phaseToArr = new int[phaseTo.size()];
        for (int j = 0; j < phaseToArr.length; j++) {
            phaseToArr[j] = phaseTo.get(j);
        }

        S = S.transpose();

        final Matrix nvec0 = new Matrix(NS, 1);
        @SuppressWarnings("unchecked")
        final ArrayDeque<Integer>[] buffers0 = new ArrayDeque[I];
        for (int ind = 0; ind < I; ind++) {
            buffers0[ind] = new ArrayDeque<Integer>();
        }
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                Matrix state_i = sn.state.get(sn.stateful.get((int) sn.nodeToStateful.get(ind)));
                if (state_i == null) {
                    throw new RuntimeException("State matrix for stateful node " + ind + " is null");
                }
                State.StateMarginalStatistics aggr = ToMarginal.toMarginalAggr(sn, ind,
                        state_i, null, null, null, null, null);
                for (int r = 0; r < R; r++) {
                    double nir = aggr.nir.get(r);
                    if (Double.isInfinite(nir)) {
                        if (sn.nodetype.get(ind) == NodeType.Source) {
                            nir = 1.0;
                        } else {
                            throw new RuntimeException("Infinite population error.");
                        }
                    }
                    // Spread the class population across its phases. The marginal
                    // the initial state carries is per class, not per phase, so the
                    // entry distribution pie is the natural allocation: it is the
                    // phase a job starts service in. A single-phase class puts
                    // everything in its sole slot, reproducing the old flat layout.
                    if (sm.nph[ind][r] <= 1 || bufPHClass[ind][r]) {
                        // A buffered-PH class keeps its whole population in slot 0
                        // (nvec is the class total there); the in-service phase
                        // composition lives in svcph, seeded below.
                        nvec0.set(sm.slot(ind, r, 0), 0, nir);
                    } else {
                        double[] pe = entryProbs(sn, sm, ind, r);
                        double left = nir;
                        for (int ke = 0; ke < sm.nph[ind][r]; ke++) {
                            double take;
                            if (ke == sm.nph[ind][r] - 1) {
                                take = left;
                            } else {
                                take = Math.min(left, Math.round(nir * pe[ke]));
                            }
                            nvec0.set(sm.slot(ind, r, ke), 0, take);
                            left -= take;
                        }
                    }
                }

                // Populate buffers for buffered nodes from the raw state vector
                // (only stations have FCFS/LCFS scheduling; skip non-station
                // stateful nodes such as RROBIN dispatchers/Routers and Caches)
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist >= 0 && bufferedSched.contains(sn.sched.get(sn.stations.get(ist)))) {
                    Matrix Kmat = new Matrix(1, sn.phasessz.getNumCols());
                    Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), Kmat, 0, 0);
                    int sumK = (int) Kmat.elementSum();
                    int sumNvars = (int) sn.nvars.sumRows(ind);
                    int bufCols = state_i.getNumCols() - sumK - sumNvars;

                    if (isListSched(sn.sched.get(sn.stations.get(ist)))) {
                        // PAS/OI has no server/phase split at all: the full ordered
                        // list occupies every column ahead of the routing variables,
                        // c(1) oldest, zero-padded on the right (ToMarginal and
                        // afterEventStationPas both read it as W = nCols - V). That
                        // is already the order the NRM needs, so it is copied
                        // verbatim rather than reversed.
                        for (int pos = 0; pos < state_i.getNumCols() - sumNvars; pos++) {
                            int classId = (int) state_i.get(0, pos);
                            if (classId >= 1 && classId <= R) {
                                buffers0[ind].addLast(classId);
                            }
                        }
                    } else if (isCountBufferedSched(sn.sched.get(sn.stations.get(ist)))) {
                        // SIRO/SEPT/LEPT keep an UN-ordered buffer: the first R
                        // columns hold the per-class counts of waiting jobs, not
                        // class ids (see State.fromMarginalAndRunning). Expand
                        // them into the NRM's ordered list; the order within it
                        // is immaterial for these disciplines, which select by
                        // class and never by position.
                        for (int r2 = 0; r2 < Math.min(R, bufCols); r2++) {
                            int cnt = (int) state_i.get(0, r2);
                            for (int c = 0; c < cnt; c++) {
                                buffers0[ind].addLast(r2 + 1);
                            }
                        }
                    } else {
                        // FCFS/HOL/LCFS keep an ordered list of class ids
                        for (int pos = 0; pos < bufCols; pos++) {
                            int classId = (int) state_i.get(0, pos);
                            if (classId >= 1 && classId <= R) {
                                buffers0[ind].addLast(classId);
                            }
                            // classId == 0 means empty position, skip
                        }
                    }
                }

                // Seed the in-service phase multiset of a buffered-PH node. The
                // jobs in service are the class total minus the ones waiting in the
                // buffer just built; their starting phases are drawn from the entry
                // distribution pie, the same allocation the INF/PS init uses. Only
                // in-service jobs get a phase -- waiting jobs carry none.
                if (bufPHNode[ind]) {
                    for (int r = 0; r < R; r++) {
                        int waiting = bufCount(buffers0[ind], r);
                        double insvc = Math.max(0.0, aggr.nir.get(r) - waiting);
                        if (sm.nph[ind][r] <= 1) {
                            svcph[ind][r][0] = insvc;
                        } else {
                            double[] pe = entryProbs(sn, sm, ind, r);
                            double left = insvc;
                            for (int ke = 0; ke < sm.nph[ind][r]; ke++) {
                                double take = (ke == sm.nph[ind][r] - 1)
                                        ? left : Math.min(left, Math.round(insvc * pe[ke]));
                                svcph[ind][r][ke] = take;
                                left -= take;
                            }
                        }
                    }
                }
            }
        }

        final double[] mi = new double[I];
        final double[][] rates = new double[I][R];
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstation.get(ind) != 0.0) {
                int ist = (int) sn.nodeToStation.get(ind);
                for (int r = 0; r < R; r++) {
                    double muir = sn.rates.get(ist, r);
                    if (!Double.isNaN(muir)) {
                        rates[ind][r] = muir;
                    }
                }
                mi[ind] = sn.nservers.get(ist);
            } else {
                for (int r = 0; r < R; r++) {
                    rates[ind][r] = GlobalConstants.Immediate;
                    mi[ind] = GlobalConstants.MaxInt;
                }
            }
            if (Double.isInfinite(mi[ind])) {
                mi[ind] = GlobalConstants.MaxInt;
            }
        }

        final double epstol = GlobalConstants.Zero;

        // Normalized per-class weights of the DPS/GPS-family stations.
        final double[][] wnorm = buildSchedWeights(sn, R);

        // Finite capacity regions (DROP rule); inert when the model has none.
        final Fcr fcr = fcrPrecompute(sn);

        // Balking thresholds; inert when no station declares QUEUE_LENGTH balking.
        final Balk balk = balkPrecompute(sn);

        // G-network signal classes; inert when no class is a signal.
        final Sig sig = sigPrecompute(sn);

        // Round-robin pointers, seeded from the initial state; inert when no
        // (node,class) routes with RROBIN/WRROBIN. No propensity reads them, so
        // they stay out of the reaction network.
        final Rr rrp = rrPrecompute(sn);

        // Rate of the service-process event each reaction carries: the absorption
        // mu(k)*phi(k) for a departure, the off-diagonal D0(k,k') for a phase
        // change. For a single-phase class this is just the exponential rate, so
        // an exponential model sees exactly the rates it saw before.
        final double[] rateOf = new double[fromIdx.size()];
        for (int j = 0; j < nDep; j++) {
            int ind = fromIR.get(j)[0];
            int r = fromIR.get(j)[1];
            if (isPhaseRx[j]) {
                rateOf[j] = phaseRate.get(j);
            } else if (sn.isstation.get(ind) != 0.0) {
                int ist = (int) sn.nodeToStation.get(ind);
                int ph = depPhase.get(j);
                MatrixCell proc = phEntryOf(sn.proc, sn, ist, r);
                Matrix muM = (sn.mu.get(sn.stations.get(ist)) != null)
                        ? sn.mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)) : null;
                Matrix phiM = (sn.phi.get(sn.stations.get(ist)) != null)
                        ? sn.phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)) : null;
                if (sm.nph[ind][r] > 1 && proc != null && !proc.isEmpty()
                        && muM != null && phiM != null) {
                    rateOf[j] = muM.get(ph) * phiM.get(ph);
                } else {
                    rateOf[j] = rates[ind][r];
                }
            } else {
                rateOf[j] = rates[ind][r];
            }
        }

        @SuppressWarnings("unchecked")
        final BiFunction<Matrix, ArrayDeque<Integer>[], Double>[] a = new BiFunction[fromIdx.size()];
        for (int j = 0; j < nDep; j++) {
            final int jj = j;
            a[j] = new BiFunction<Matrix, ArrayDeque<Integer>[], Double>() {
                @Override
                public Double apply(Matrix X, ArrayDeque<Integer>[] bufs) {
                    int ind = fromIR.get(jj)[0];
                    int rr = fromIR.get(jj)[1];
                    // Buffered phase-type service. Only the jobs in service carry a
                    // phase, and their per-phase counts live in svcph[ind][r][k],
                    // not in nvec. Both the departure (absorption of phase kk) and
                    // the internal phase transition (kk -> kb) therefore fire at
                    // rateOf(j) times the number of class-rr jobs in service in the
                    // source phase kk -- exactly the INF-family law rateOf*kir, but
                    // with kir read from the in-service multiset svcph rather than
                    // from nvec (whose class total also counts the waiting jobs). The
                    // load-/class-dependent factors still read the total population.
                    if (sn.isstation.get(ind) != 0.0 && bufPHClass[ind][rr]) {
                        int istB = (int) sn.nodeToStation.get(ind);
                        int kk = depPhase.get(jj);
                        double[] npopB = classCounts(X, sm, ind);
                        return rateOf[jj] * svcph[ind][rr][kk]
                                * lldfac(sn, istB, sumPop(npopB))
                                * cdfacPop(sn, istB, npopB, R, rr);
                    }
                    // rate(class r, phase k) = (kir_k/nir) * classShare_r(class
                    // counts) * rate_k. Phase expansion touches neither the class
                    // share, which the helpers below still read off the class-level
                    // counts, nor the class-level counts themselves: it only splits
                    // the share across a class's phases in the ratio kir/nir. For a
                    // single-phase class kirFrac is 1 and classPop is the old slot
                    // value, so every rate is exactly the pre-expansion one.
                    double base = rateOf[jj] * kirFrac(X, sm, fromIdx.get(jj), ind, rr);
                    if (sn.isstation.get(ind) == 0.0) {
                        return base * Math.min(1.0, classPop(X, sm, ind, rr));
                    }
                    int ist = (int) sn.nodeToStation.get(ind);
                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    if (sched == SchedStrategy.EXT) {
                        // A Source's arrival rate is state-independent: its slot
                        // holds a fictitious token that the departure consumes and
                        // nothing replenishes, so scaling by the token count (as
                        // the population-based laws below do) would latch the rate
                        // to zero the first time the count crosses 0 and silence
                        // the Source forever. Phase-type arrivals never reach here:
                        // the source token is not a population of jobs in phases,
                        // so phaseNrmOK keeps non-exponential EXT on the serial
                        // engine and nph is always 1 at a Source.
                        return rateOf[jj];
                    }
                    if (sched == SchedStrategy.INF) {
                        return base * classPop(X, sm, ind, rr)
                                * cdfacPop(sn, ist, classCounts(X, sm, ind), R, rr);
                    }
                    if (sched == SchedStrategy.PS || sched == SchedStrategy.LPS) {
                        // LPS shares the PS rate law in State.afterEventStation:
                        // the sharing limit is the server count, so min(ni,c)
                        // already covers it.
                        if (R == 1) {
                            double nir = classPop(X, sm, ind, rr);
                            return base * Math.min(mi[ind], nir)
                                    * lldfac(sn, ist, nir)
                                    * cdfacPop(sn, ist, classCounts(X, sm, ind), R, rr);
                        }
                        double[] npop = classCounts(X, sm, ind);
                        double total = epstol + sumPop(npop);
                        return base * (npop[rr] / total) * Math.min(mi[ind], total)
                                * lldfac(sn, ist, sumPop(npop))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    if (sched == SchedStrategy.DPS) {
                        double[] npop = classCounts(X, sm, ind);
                        return base * dpsshare(wnorm[ist], npop, rr)
                                * lldfac(sn, ist, sumPop(npop))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    if (sched == SchedStrategy.GPS) {
                        double[] npop = classCounts(X, sm, ind);
                        return base * gpsshare(wnorm[ist], npop, rr)
                                * lldfac(sn, ist, sumPop(npop))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    if (sched == SchedStrategy.PSPRIO) {
                        double[] npop = classCounts(X, sm, ind);
                        return base * psprioshare(sn, npop, rr, mi[ind])
                                * lldfac(sn, ist, prioPop(sn, npop, rr, mi[ind]))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    if (sched == SchedStrategy.DPSPRIO) {
                        double[] npop = classCounts(X, sm, ind);
                        return base * dpsprioshare(sn, wnorm[ist], npop, rr, mi[ind])
                                * lldfac(sn, ist, prioPop(sn, npop, rr, mi[ind]))
                                * cdfacPop(sn, ist, prioVec(sn, npop, rr, mi[ind]), R, rr);
                    }
                    if (sched == SchedStrategy.GPSPRIO) {
                        double[] npop = classCounts(X, sm, ind);
                        return base * gpsprioshare(sn, wnorm[ist], npop, rr, mi[ind])
                                * lldfac(sn, ist, prioPop(sn, npop, rr, mi[ind]))
                                * cdfacPop(sn, ist, prioVec(sn, npop, rr, mi[ind]), R, rr);
                    }
                    if (isListSched(sched)) {
                        // Position p of the ordered list is served at
                        // Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}), and
                        // pass-and-swap decides which class that completion
                        // ejects. The class-r departure rate is therefore the
                        // total Delta_mu over the positions whose pass-and-swap
                        // ejects a class-r job, which is exactly what
                        // afterEventStationPas enumerates. mu(c) already carries
                        // any load/class dependence, so no lld/cd factor applies.
                        return oirate(sn, ind, bufs[ind], rr);
                    }
                    if (isBufferedSched(sched)) {
                        // Invariant: buffers[ind].size == max(0, total - mi[ind]),
                        // except at a retrial station, where a departure does not
                        // promote and the orbit can be occupied while servers idle.
                        // State update must ensure that buffer only fills when all servers busy.
                        int waiting = 0;
                        for (Integer cl : bufs[ind]) {
                            if (cl == rr + 1) waiting++;
                        }
                        double[] npop = classCounts(X, sm, ind);
                        double inService = npop[rr] - waiting;
                        if (inService <= 0) {
                            return 0.0;
                        }
                        // rate proportional to jobs actually being served, scaled by
                        // the load-dependent factor at the total station population
                        return base * inService * lldfac(sn, ist, sumPop(npop))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    if (sched == SchedStrategy.POLLING) {
                        // A polling station has a single server that serves exactly
                        // one job, of the class the controller currently attends. The
                        // class-r departure therefore fires only while the controller
                        // is SERVING class r, at the plain service rate of the one job
                        // in service -- never scaled by the class population, since the
                        // other class-r jobs wait in the buffer for the server to come
                        // back to them.
                        int[] ctrl = pollCtrl[ind];
                        if (ctrl == null || ctrl[0] != Polling.MODE_VISIT || ctrl[1] != rr) {
                            return 0.0;
                        }
                        double[] npop = classCounts(X, sm, ind);
                        return base * lldfac(sn, ist, sumPop(npop))
                                * cdfacPop(sn, ist, npop, R, rr);
                    }
                    return base * classPop(X, sm, ind, rr);
                }
            };
        }

        // Reneging propensities: only the jobs actually waiting can abandon, so
        // the rate is the class-r buffer occupancy, exactly the quantity the
        // FCFS-family rate law already relies on.
        for (int j = nDep; j < a.length; j++) {
            final int ind = fromIR.get(j)[0];
            final int rr = fromIR.get(j)[1];
            if (isRetryRx[j]) {
                // Retrial propensity: only jobs actually in orbit retry, and only
                // a free server admits them. Because a departure does not promote
                // at a retrial station, the buffer invariant of the other buffered
                // policies is FALSE here -- the orbit can be full while servers
                // idle -- so "total > mi" no longer means "the servers are busy"
                // and the test must consult the in-service count directly.
                final double mu = retryMu.get(j);
                a[j] = new BiFunction<Matrix, ArrayDeque<Integer>[], Double>() {
                    @Override
                    public Double apply(Matrix X, ArrayDeque<Integer>[] bufs) {
                        double total = sumPop(classCounts(X, sm, ind));
                        if (total - bufs[ind].size() >= mi[ind]) {
                            return 0.0;
                        }
                        return mu * bufCount(bufs[ind], rr);
                    }
                };
                continue;
            }
            final double mu = renegeMu.get(j);
            a[j] = new BiFunction<Matrix, ArrayDeque<Integer>[], Double>() {
                @Override
                public Double apply(Matrix X, ArrayDeque<Integer>[] bufs) {
                    return mu * bufCount(bufs[ind], rr);
                }
            };
        }

        // Polling switchover propensities. The renege/retry loop above set these
        // appended columns to a zero-rate renege closure; restore the switchover
        // law here. A switchover fires only while the controller is walking (mode
        // SWITCH) at the total leaving rate -D0(swk,swk) of the current phase swk
        // of the switchover PH into buffer pos. The competition between advancing
        // to another phase and absorbing (arriving at pos) is resolved at firing
        // time, exactly as a routed departure resolves its destination.
        for (int j = 0; j < isPollSwRx.length; j++) {
            if (!isPollSwRx[j]) continue;
            final int pind = pollSwNode[j];
            final Polling.Info pinf = pinfo[pind];
            a[j] = new BiFunction<Matrix, ArrayDeque<Integer>[], Double>() {
                @Override
                public Double apply(Matrix X, ArrayDeque<Integer>[] bufs) {
                    int[] ctrl = pollCtrl[pind];
                    if (ctrl == null || ctrl[0] != Polling.MODE_SWITCH) {
                        return 0.0;
                    }
                    int pos = ctrl[1];
                    int swk = ctrl[2];
                    Matrix D0 = pinf.swD0[pos];
                    return -D0.get(swk, swk);
                }
            };
        }

        // Finite capacity regions do NOT gate the propensities. Under the DROP
        // rule the refused job is DESTROYED, not held back: the departure fires
        // at its full rate and the job simply never reaches the destination.
        // Scaling the propensity by the admitted share instead censors the
        // transition, which keeps the job at its SOURCE -- a different model,
        // and one that diverges as soon as the source is a real queue rather
        // than a Source node (an interior region makes the upstream queue grow
        // without bound while nothing is ever lost). The two coincide only at a
        // Source, whose population is fictitious, which is why every FCR fixture
        // placed a region on a Source-fed station and never saw the difference.
        // Refusal is applied at firing time on the drawn destination instead,
        // exactly as a balk is (see fireReaction), matching SOLVER_SSA and the
        // exact CTMC.

        @SuppressWarnings("unchecked")
        final List<Integer>[] D = new List[S.getNumCols()];
        for (int kk = 0; kk < D.length; kk++) {
            D[kk] = new ArrayList<Integer>();
            if (isRetryRx[kk]) {
                // A retry has an all-zero stoichiometry column, so the generic
                // derivation below would return an EMPTY dependency set and leave
                // every rate at the node stale. A retry does change the in-service
                // composition, hence every reaction whose source is this node.
                int indK = fromIR.get(kk)[0];
                for (int j = 0; j < fromIdx.size(); j++) {
                    if (sm.node[fromIdx.get(j)] == indK) {
                        D[kk].add(j);
                    }
                }
                continue;
            }
            List<Integer> J = new ArrayList<Integer>();
            for (int i = 0; i < S.getNumRows(); i++) {
                if (S.get(i, kk) != 0.0) J.add(i);
            }
            List<Integer> vecd = new ArrayList<Integer>();
            for (int j = 0; j < J.size(); j++) {
                // Decode through the slot map, never arithmetically: with phase
                // expansion a state index is a (node,class,PHASE) slot, so pos%R
                // names the WRONG node as soon as any class has more than one
                // phase, and D then refreshes the wrong propensities and leaves
                // rates stale. The defect is invisible to an exponential model,
                // where the arithmetic decode still happens to be correct, so a
                // green all-exponential test proves nothing about it. Collect
                // EVERY slot of each affected node, because a rate law reads its
                // node's whole class-count vector (classCounts sums each class
                // over its phases) and the per-phase share reads the sibling
                // phases of its own class.
                int ind = sm.node[J.get(j)];
                for (int rcls = 0; rcls < R; rcls++) {
                    for (int ph = 0; ph < sm.nph[ind][rcls]; ph++) {
                        vecd.add(sm.slot(ind, rcls, ph));
                    }
                }
            }
            if (!vecd.isEmpty()) {
                Set<Integer> vecdUniqueSet = new java.util.LinkedHashSet<Integer>(vecd);
                List<Integer> vecdUnique = new ArrayList<Integer>(vecdUniqueSet);
                Set<Integer> vecsUnique = new java.util.LinkedHashSet<Integer>();
                for (Integer u : vecdUnique) {
                    for (int kk2 = 0; kk2 < S.getNumCols(); kk2++) {
                        // No fcr.on widening: regions no longer gate the
                        // propensities (see the FCR note above), so a
                        // departure's rate depends only on its own station's
                        // populations, as in the unregulated case. Admission is
                        // resolved at firing time on the drawn destination and
                        // changes no rate.
                        boolean touches = S.get(u, kk2) < 0.0;
                        if (touches) {
                            vecsUnique.add(kk2);
                        }
                    }
                }
                D[kk].addAll(vecsUnique);
            }
        }

        // The derivation above collects reactions by the SIGN of their S entries,
        // so it can never place an all-zero retry column in any OTHER reaction's
        // dependency set. A retry still has to be refreshed whenever the node it
        // serves changes, because its rate reads both the orbit occupancy and
        // whether a server is free: without this patch a retry blocked at a busy
        // server keeps its zero rate after the server frees, so the orbit never
        // drains and the station grows without bound.
        for (int j = 0; j < isRetryRx.length; j++) {
            if (!isRetryRx[j]) continue;
            int indj = fromIR.get(j)[0];
            for (int kk = 0; kk < S.getNumCols(); kk++) {
                boolean touches = fromIR.get(kk)[0] == indj;
                // Every slot of the node, not just R of them: with phase expansion
                // a class owns several consecutive slots.
                for (int s = 0; s < R && !touches; s++) {
                    for (int ph = 0; ph < sm.nph[indj][s] && !touches; ph++) {
                        if (S.get(sm.slot(indj, s, ph), kk) != 0.0) touches = true;
                    }
                }
                if (touches && !D[kk].contains(Integer.valueOf(j))) {
                    D[kk].add(j);
                }
            }
        }

        for (int i = 0; i < S.getNumRows(); i++) {
            for (int j = 0; j < S.getNumCols(); j++) {
                if (Double.isInfinite(S.get(i, j))) {
                    S.set(i, j, 0.0);
                }
            }
        }

        int M = sn.nstations;
        int K = sn.nclasses;

        Matrix QN = new Matrix(M, K); QN.fill(0.0);
        Matrix UN = new Matrix(M, K); UN.fill(0.0);
        Matrix RN = new Matrix(M, K); RN.fill(0.0);
        Matrix TN = new Matrix(M, K); TN.fill(0.0);
        Matrix CN = new Matrix(1, K); CN.fill(0.0);
        Matrix XN = new Matrix(1, K); XN.fill(0.0);

        long[][] cacheProd = new long[I][R];   // per (cache node, PRODUCED class) count
        nrm_direct(S, D, a, nvec0, buffers0, samples, options, QN, UN, RN, TN, CN, XN, sn, fromIdx, fromIR, mi, fcr, balk, sig, rrp, isRenegeRx, isRetryRx, sm, nDepRx, isPhaseRx, anyPoll, pinfo, isPollNode, pollCtrl, isPollSwRx, pollSwNode, svcph, bufPHNode, isBufSvcRx, depPhase, phaseToArr, isCacheRx, cacheHitSlotArr, cacheMissSlotArr, isCacheNode, cacheContents0, cacheProd, cacheRetrDest);

        // Write the measured hit/miss probabilities into each Cache node's param
        // (afterEventCache convention: actualhitprob(r) = hit throughput /
        // (hit+miss) throughput at the cache, per read class r), so the node-table
        // reconstruction reports the cache hit/miss class throughputs.
        for (int ind = 0; ind < I; ind++) {
            if (isCacheNode[ind]) {
                // Write onto the ORIGINAL struct sn_in (not the internal copy):
                // its nodeparam is keyed by the model's node objects, which is the
                // struct SolverSSA reads to push the hit prob onto the cache node.
                jline.lang.nodeparam.CacheNodeParam cnp =
                        (jline.lang.nodeparam.CacheNodeParam) sn_in.nodeparam.get(sn_in.nodes.get(ind));
                if (cnp.actualhitprob == null || cnp.actualhitprob.length() < K) {
                    cnp.actualhitprob = new Matrix(1, K); cnp.actualhitprob.fill(Double.NaN);
                }
                if (cnp.actualmissprob == null || cnp.actualmissprob.length() < K) {
                    cnp.actualmissprob = new Matrix(1, K); cnp.actualmissprob.fill(Double.NaN);
                }
                // The hit/miss probability of a read class is the throughput of its
                // hit class over hit+miss at the cache, from the per-produced-class
                // counts. A retrieval completion produces the miss class (counted
                // here); a delayed hit produces nothing (excluded). Retrieval classes
                // (hitclass 0) are internal and get no probability of their own.
                for (int r = 0; r < R; r++) {
                    if (isCacheReadClass[ind][r] && r < cnp.hitclass.length() && cnp.hitclass.get(r) > 0) {
                        int hc = (int) cnp.hitclass.get(r);
                        int mc = (int) cnp.missclass.get(r);
                        long tot = cacheProd[ind][hc] + cacheProd[ind][mc];
                        if (tot > 0) {
                            cnp.actualhitprob.set(r, (double) cacheProd[ind][hc] / tot);
                            cnp.actualmissprob.set(r, (double) cacheProd[ind][mc] / tot);
                        }
                    }
                }
            }
        }

        return new SolverSSAResultNRM(QN, UN, RN, TN, CN, XN, sn);
    }

    // ======================================================================
    // Stochastic Petri net (Place / Transition) via the Next-Reaction Method
    //
    // A Place holds a per-class token count (a population slot of the state
    // vector); a timed Transition mode is a reaction whose stoichiometry column
    // is the arc incidence -- input (enabling) arcs consume, output (firing)
    // arcs produce. Enabling is a propensity gate (all input places at or above
    // their arc weight, every inhibitor place strictly below its threshold); a
    // single-server mode fires at its exponential rate, an infinite/k-server
    // mode at that rate times its enabling degree. Each firing applies the
    // mode's stoichiometry once, the atomic GSPN firing shared by the exact CTMC
    // (single server), JMT and GreatSPN. IMMEDIATE modes fire in zero time and
    // are resolved by vanishing-marking elimination (see spnCollapse). A
    // non-exponential firing distribution needs per-mode in-flight phase state
    // the reaction network does not carry and is rejected by the featset.
    // ======================================================================

    /** One transition-mode reaction of a stochastic Petri net. */
    private static final class SpnRx {
        int node;
        int mode;
        double[] Svec;       // NS net stoichiometry (firing minus enabling)
        int[] enSlot;        // input-place slots
        double[] enW;        // input-arc weights
        int[] inhSlot;       // inhibitor-place slots
        double[] inhThr;     // inhibitor thresholds
        double baseRate;     // exponential firing rate (sum of D1)
        double nservers;     // mode servers (MaxInt for infinite)
        double weight;       // firing weight (immediate conflict resolution)
        double prio;         // firing priority (immediate conflict resolution)
    }

    public static SolverSSAResultNRM nrm_spn(NetworkStruct sn, SolverOptions options, Smap sm, int samples) {
        final int R = sn.nclasses;
        final int I = sn.nnodes;
        final int M = sn.nstations;
        final int K = sn.nclasses;
        final int NS = sm.NS;

        List<SpnRx> rx = new ArrayList<SpnRx>();     // timed modes
        List<SpnRx> imm = new ArrayList<SpnRx>();    // immediate modes
        // consumers[ind*R+c]: indices into rx of timed modes consuming from place
        // ind, class c. Place throughput is their aggregate firing rate (once per
        // firing, unweighted -- the depRates the CTMC accumulates from PRE events).
        Map<Integer, List<Integer>> consumers = new java.util.HashMap<Integer, List<Integer>>();
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Transition) continue;
            TransitionNodeParam tp = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
            for (int m = 0; m < tp.nmodes; m++) {
                SpnRx rec = spnBuildMode(sn, sm, ind, m, tp, NS);
                if (tp.timing.get(m) == TimingStrategy.IMMEDIATE) {
                    imm.add(rec);
                } else {
                    rx.add(rec);
                    int ridx = rx.size() - 1;
                    for (int a = 0; a < rec.enSlot.length; a++) {
                        int p = sm.node[rec.enSlot[a]];
                        int c = sm.cls[rec.enSlot[a]];
                        int key = p * R + c;
                        List<Integer> lst = consumers.get(key);
                        if (lst == null) { lst = new ArrayList<Integer>(); consumers.put(key, lst); }
                        lst.add(ridx);
                    }
                }
            }
        }
        // Source arrivals. A Source is not a Transition, so its Poisson arrival is
        // not one of the transition modes above; it needs its own reaction or the
        // fed Place stays empty and the net deadlocks. Add one arrival reaction per
        // (Source node, open class, routed Place-class edge). Splitting a Poisson
        // stream by the independent routing probabilities yields independent Poisson
        // streams, so an edge of probability p carries rate lambda*p exactly. The
        // reaction has an EMPTY enabling set (always enabled, state-independent
        // propensity = lambda*p) and deposits +1 token into the routed Place slot.
        // producers[node*R+class] indexes these so the Source station reports its
        // arrival rate as throughput (the open class's reference-station throughput).
        Map<Integer, List<Integer>> producers = new java.util.HashMap<Integer, List<Integer>>();
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Source) continue;
            int ist = (int) sn.nodeToStation.get(ind);
            for (int r = 0; r < R; r++) {
                double lambda = sn.rates.get(ist, r);
                if (Double.isNaN(lambda) || lambda <= 0.0) continue;
                Map<JobClass, ProcessType> pmap = sn.procid.get(sn.stations.get(ist));
                ProcessType pt = (pmap != null) ? pmap.get(sn.jobclasses.get(r)) : null;
                if (pt != null && pt != ProcessType.EXP) {
                    throw new RuntimeException("Source node " + ind + " class " + r
                            + " has a non-exponential arrival, which the NRM SPN path does not support.");
                }
                boolean foundPlace = false;
                for (int jnd = 0; jnd < I; jnd++) {
                    if (sn.nodetype.get(jnd) != NodeType.Place) continue;
                    for (int s = 0; s < R; s++) {
                        double p = sn.rtnodes.get(ind * R + r, jnd * R + s);
                        if (p <= 0.0) continue;
                        foundPlace = true;
                        SpnRx rec = new SpnRx();
                        rec.node = ind;
                        rec.mode = -1;   // arrival, not a transition mode
                        double[] Svec = new double[NS];
                        Svec[sm.slot(jnd, s, 0)] += 1.0;
                        rec.Svec = Svec;
                        rec.enSlot = new int[0];
                        rec.enW = new double[0];
                        rec.inhSlot = new int[0];
                        rec.inhThr = new double[0];
                        rec.baseRate = lambda * p;   // Poisson thinning by the routing prob
                        rec.nservers = 1.0;          // constant propensity = baseRate
                        rec.weight = 1.0;
                        rec.prio = 1.0;
                        rx.add(rec);
                        int ridx = rx.size() - 1;
                        int key = ind * R + r;
                        List<Integer> lst = producers.get(key);
                        if (lst == null) { lst = new ArrayList<Integer>(); producers.put(key, lst); }
                        lst.add(ridx);
                    }
                }
                if (!foundPlace) {
                    throw new RuntimeException("Source node " + ind + " class " + r
                            + " does not route to any Place; the NRM SPN path needs a Source->Place arc.");
                }
            }
        }

        final int nR = rx.size();
        if (nR == 0) {
            throw new RuntimeException("Stochastic Petri net has no timed reaction; nothing to simulate.");
        }

        // Initial marking: token counts per (place, class) from the initial state.
        Matrix nvec = new Matrix(NS, 1);
        nvec.fill(0.0);
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Place || sn.isstateful.get(ind) == 0.0) continue;
            Matrix state_i = sn.state.get(sn.stateful.get((int) sn.nodeToStateful.get(ind)));
            State.StateMarginalStatistics aggr = ToMarginal.toMarginalAggr(sn, ind, state_i, null, null, null, null, null);
            for (int c = 0; c < R; c++) {
                double nir = aggr.nir.get(c);
                if (Double.isInfinite(nir)) {
                    throw new RuntimeException("Infinite marking at a Place is not supported.");
                }
                nvec.set(sm.slot(ind, c, 0), 0, nir);
            }
        }

        final int maxImmSteps = 100000;

        // Finite-capacity Place DROP enforcement. A Place with a finite per-class
        // capacity (sn.classcap) or total capacity (sn.cap) loses any arriving token
        // that would exceed it (JMT/CTMC loss semantics: an M/M/1/1 Place with cap 1
        // holds mean 0.333 at rho=0.5, not the unbounded-M/M/1 value 1.0). Without
        // this the deposit nvec + Svec accumulates tokens past capacity. Precompute
        // the per-slot per-class caps, the per-place total caps, and each reaction's
        // deposited slots so the clamp in the loop touches only what just grew.
        // Mirrors the MATLAB solver_ssa_nrm applyPlaceCaps and the Python native
        // _apply_place_caps.
        double[] pcapSlot = new double[NS];                                  // per-(place,class) slot cap
        java.util.Arrays.fill(pcapSlot, Double.POSITIVE_INFINITY);
        List<Double> placeTotalCap = new ArrayList<Double>();                // total cap per capped place
        List<int[]> placeTotalSlots = new ArrayList<int[]>();                // its slots
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Place || sn.isstateful.get(ind) == 0.0) continue;
            int ist = (int) sn.nodeToStation.get(ind);
            int[] slotsHere = new int[R];
            for (int c = 0; c < R; c++) {
                int slot = sm.slot(ind, c, 0);
                slotsHere[c] = slot;
                if (sn.classcap != null && ist < sn.classcap.getNumRows() && c < sn.classcap.getNumCols()) {
                    double cc = sn.classcap.get(ist, c);
                    if (isFiniteCap(cc)) pcapSlot[slot] = cc;
                }
            }
            double tcap = Double.POSITIVE_INFINITY;
            if (sn.cap != null && ist < sn.cap.getNumRows() && sn.cap.getNumCols() > 0) {
                tcap = sn.cap.get(ist, 0);
            }
            if (isFiniteCap(tcap)) {
                placeTotalCap.add(tcap);
                placeTotalSlots.add(slotsHere);
            }
        }
        boolean hasPlaceCaps = !placeTotalCap.isEmpty();
        for (int s = 0; s < NS && !hasPlaceCaps; s++) {
            if (isFiniteCap(pcapSlot[s])) hasPlaceCaps = true;
        }
        int[][] depSlots = new int[nR][];
        for (int k = 0; k < nR; k++) {
            double[] Svec = rx.get(k).Svec;
            int cnt = 0;
            for (int s = 0; s < NS; s++) if (Svec[s] > 0.0) cnt++;
            int[] dep = new int[cnt];
            int at = 0;
            for (int s = 0; s < NS; s++) if (Svec[s] > 0.0) dep[at++] = s;
            depSlots[k] = dep;
        }
        int[] allSlots = new int[NS];
        for (int s = 0; s < NS; s++) allSlots[s] = s;

        // ------------------------------------------------------------------
        // Next-Reaction Method run loop
        // ------------------------------------------------------------------
        spnCollapse(nvec, imm, maxImmSteps);
        if (hasPlaceCaps) applyPlaceCaps(nvec, allSlots, pcapSlot, placeTotalCap, placeTotalSlots);
        double[] Ak = new double[nR];
        for (int k = 0; k < nR; k++) Ak[k] = spnProp(nvec, rx.get(k));
        double[] Pk = new double[nR];
        for (int k = 0; k < nR; k++) Pk[k] = -Math.log(Maths.rand());
        double[] Tk = new double[nR];
        double[] tau = new double[nR];
        for (int k = 0; k < nR; k++) tau[k] = (Ak[k] == 0.0) ? Double.POSITIVE_INFINITY : (Pk[k] - Tk[k]) / Ak[k];

        Matrix QN = new Matrix(M, K); QN.fill(0.0);
        Matrix UN = new Matrix(M, K); UN.fill(0.0);
        Matrix RN = new Matrix(M, K); RN.fill(0.0);
        Matrix TN = new Matrix(M, K); TN.fill(0.0);
        Matrix CN = new Matrix(1, K); CN.fill(0.0);
        Matrix XN = new Matrix(1, K); XN.fill(0.0);
        Matrix NK = sn.njobs.transpose();
        double totalTime = 0.0;

        int n = 0;
        while (n < samples) {
            int kfire = -1;
            double dt = Double.POSITIVE_INFINITY;
            for (int k = 0; k < nR; k++) {
                if (tau[k] < dt) { dt = tau[k]; kfire = k; }
            }
            if (Double.isInfinite(dt)) {
                throw new RuntimeException("Deadlock: no transition is enabled. Quitting nrm method.");
            }
            totalTime += dt;

            for (int ist = 0; ist < M; ist++) {
                int ind = (int) sn.stationToNode.get(ist);
                for (int c = 0; c < K; c++) {
                    double tokens = classPop(nvec, sm, ind, c);
                    QN.set(ist, c, QN.get(ist, c) + tokens * dt);
                    UN.set(ist, c, UN.get(ist, c) + tokens * dt);
                    double depr = 0.0;
                    List<Integer> cons = consumers.get(ind * R + c);
                    if (cons != null) {
                        for (int a = 0; a < cons.size(); a++) depr += Ak[cons.get(a)];
                    }
                    // A Source station has no consuming transition; its throughput
                    // is the aggregate arrival rate it injects (producers).
                    List<Integer> prod = producers.get(ind * R + c);
                    if (prod != null) {
                        for (int a = 0; a < prod.size(); a++) depr += Ak[prod.get(a)];
                    }
                    TN.set(ist, c, TN.get(ist, c) + depr * dt);
                }
            }

            // Fire the selected timed mode (single atomic firing), then collapse any
            // immediate transitions the new marking enabled. A finite-capacity DROP
            // Place loses any token the firing pushed above its capacity, before the
            // immediate cascade sees the new marking.
            SpnRx rf = rx.get(kfire);
            for (int s = 0; s < NS; s++) {
                if (rf.Svec[s] != 0.0) nvec.set(s, 0, nvec.get(s, 0) + rf.Svec[s]);
            }
            if (hasPlaceCaps) applyPlaceCaps(nvec, depSlots[kfire], pcapSlot, placeTotalCap, placeTotalSlots);
            spnCollapse(nvec, imm, maxImmSteps);

            for (int k = 0; k < nR; k++) Tk[k] += Ak[k] * dt;
            for (int k = 0; k < nR; k++) Ak[k] = spnProp(nvec, rx.get(k));
            Pk[kfire] -= Math.log(Maths.rand());
            for (int k = 0; k < nR; k++) tau[k] = (Ak[k] == 0.0) ? Double.POSITIVE_INFINITY : (Pk[k] - Tk[k]) / Ak[k];

            n++;
        }

        if (totalTime > 0.0) {
            for (int ist = 0; ist < M; ist++) {
                for (int c = 0; c < K; c++) {
                    QN.set(ist, c, QN.get(ist, c) / totalTime);
                    UN.set(ist, c, UN.get(ist, c) / totalTime);
                    TN.set(ist, c, TN.get(ist, c) / totalTime);
                }
            }
        }
        for (int c = 0; c < K; c++) {
            XN.set(0, c, TN.get((int) sn.refstat.get(c), c));
            for (int ist = 0; ist < M; ist++) {
                if (TN.get(ist, c) > 0.0) RN.set(ist, c, QN.get(ist, c) / TN.get(ist, c));
            }
            if (XN.get(0, c) > 0.0) CN.set(0, c, NK.get(c) / XN.get(0, c));
        }
        return new SolverSSAResultNRM(QN, UN, RN, TN, CN, XN, sn);
    }

    private static SpnRx spnBuildMode(NetworkStruct sn, Smap sm, int ind, int m, TransitionNodeParam tp, int NS) {
        final int R = sn.nclasses;
        final int I = sn.nnodes;
        SpnRx rec = new SpnRx();
        rec.node = ind;
        rec.mode = m;
        double[] Svec = new double[NS];
        List<Integer> enSlot = new ArrayList<Integer>();
        List<Double> enW = new ArrayList<Double>();
        Matrix en = tp.enabling.get(m);
        Matrix fir = tp.firing.get(m);
        Matrix inh = tp.inhibiting.get(m);
        for (int p = 0; p < I; p++) {
            for (int c = 0; c < R; c++) {
                double w = en.get(p, c);
                if (w > 0.0) {
                    int slot = sm.slot(p, c, 0);
                    enSlot.add(slot);
                    enW.add(w);
                    Svec[slot] -= w;
                }
                double fw = fir.get(p, c);
                if (fw > 0.0) {
                    Svec[sm.slot(p, c, 0)] += fw;
                }
            }
        }
        List<Integer> inhSlot = new ArrayList<Integer>();
        List<Double> inhThr = new ArrayList<Double>();
        for (int p = 0; p < I; p++) {
            for (int c = 0; c < R; c++) {
                double thr = inh.get(p, c);
                if (!Double.isInfinite(thr)) {
                    inhSlot.add(sm.slot(p, c, 0));
                    inhThr.add(thr);
                }
            }
        }
        rec.Svec = Svec;
        rec.enSlot = toIntArray(enSlot);
        rec.enW = toDoubleArray(enW);
        rec.inhSlot = toIntArray(inhSlot);
        rec.inhThr = toDoubleArray(inhThr);
        if (tp.timing.get(m) != TimingStrategy.IMMEDIATE) {
            double fK = tp.firingphases.get(0, m);
            // firingproc / firingpie are keyed by the Mode object (not the mode
            // index); fetch the mode from the transition node to look it up.
            jline.lang.Mode modeObj = ((jline.lang.nodes.Transition) sn.nodes.get(ind)).getModes().get(m);
            MatrixCell proc = tp.firingproc.get(modeObj);
            if (Double.isNaN(fK) || fK != 1.0 || proc == null || proc.isEmpty()) {
                throw new RuntimeException("Transition node " + ind + " mode " + m
                        + " has non-exponential firing, which the NRM SPN path does not support.");
            }
            rec.baseRate = proc.get(1).elementSum();
        }
        double ns = tp.nmodeservers.get(0, m);
        if (Double.isInfinite(ns)) ns = GlobalConstants.MaxInt;
        rec.nservers = ns;
        rec.weight = tp.fireweight.get(0, m);
        rec.prio = tp.firingprio.get(0, m);
        return rec;
    }

    private static int[] toIntArray(List<Integer> l) {
        int[] a = new int[l.size()];
        for (int i = 0; i < l.size(); i++) a[i] = l.get(i);
        return a;
    }

    private static double[] toDoubleArray(List<Double> l) {
        double[] a = new double[l.size()];
        for (int i = 0; i < l.size(); i++) a[i] = l.get(i);
        return a;
    }

    /**
     * Enabling degree of a mode: min over input arcs of floor(tokens/weight),
     * zeroed by any active inhibitor arc. A mode with no input arc is treated as
     * single-degree.
     */
    private static double spnEnDegree(Matrix nvec, SpnRx rx) {
        for (int i = 0; i < rx.inhSlot.length; i++) {
            if (nvec.get(rx.inhSlot[i], 0) >= rx.inhThr[i]) return 0.0;
        }
        if (rx.enSlot.length == 0) return 1.0;
        double d = Double.POSITIVE_INFINITY;
        for (int i = 0; i < rx.enSlot.length; i++) {
            d = Math.min(d, Math.floor(nvec.get(rx.enSlot[i], 0) / rx.enW[i]));
        }
        return d;
    }

    /**
     * Propensity of a timed mode: exponential rate times the effective server
     * count, min(enabling degree, mode servers). Single-server modes therefore
     * fire at their rate whenever enabled; infinite/k-server modes scale by the
     * enabling degree.
     */
    private static double spnProp(Matrix nvec, SpnRx rx) {
        double d = spnEnDegree(nvec, rx);
        double eff = Math.min(d, rx.nservers);
        if (eff <= 0.0) return 0.0;
        return rx.baseRate * eff;
    }

    /** A capacity is enforceable only if it is a finite number (Inf or NaN means none). */
    private static boolean isFiniteCap(double c) {
        return !Double.isInfinite(c) && !Double.isNaN(c);
    }

    /**
     * Drop tokens a firing pushed above a Place per-class or total capacity. Only
     * the just-deposited slots (Svec &gt; 0) can overflow, so the clamp is local.
     * Mirrors the MATLAB applyPlaceCaps and the Python native _apply_place_caps.
     */
    private static void applyPlaceCaps(Matrix nvec, int[] deposited, double[] pcapSlot,
                                       List<Double> placeTotalCap, List<int[]> placeTotalSlots) {
        for (int a = 0; a < deposited.length; a++) {
            int j = deposited[a];
            if (nvec.get(j, 0) > pcapSlot[j]) {
                nvec.set(j, 0, pcapSlot[j]);
            }
        }
        for (int p = 0; p < placeTotalCap.size(); p++) {
            double tcap = placeTotalCap.get(p);
            int[] slots = placeTotalSlots.get(p);
            double sum = 0.0;
            for (int a = 0; a < slots.length; a++) sum += nvec.get(slots[a], 0);
            double excess = sum - tcap;
            if (excess > 0.0) {
                for (int a = 0; a < deposited.length && excess > 0.0; a++) {
                    int j = deposited[a];
                    boolean inPlace = false;
                    for (int b = 0; b < slots.length; b++) {
                        if (slots[b] == j) { inPlace = true; break; }
                    }
                    if (inPlace && nvec.get(j, 0) > 0.0) {
                        double d = Math.min(excess, nvec.get(j, 0));
                        nvec.set(j, 0, nvec.get(j, 0) - d);
                        excess -= d;
                    }
                }
            }
        }
    }

    /**
     * Vanishing-marking elimination: fire enabled immediate transitions until the
     * marking is tangible, highest firing priority first and, among equal
     * priority, in proportion to firing weight. Immediate firings take zero time.
     */
    private static void spnCollapse(Matrix nvec, List<SpnRx> imm, int maxsteps) {
        if (imm.isEmpty()) return;
        int steps = 0;
        while (true) {
            List<Integer> enabled = new ArrayList<Integer>();
            for (int m = 0; m < imm.size(); m++) {
                if (spnEnDegree(nvec, imm.get(m)) >= 1.0) enabled.add(m);
            }
            if (enabled.isEmpty()) return;
            double maxprio = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < enabled.size(); i++) {
                maxprio = Math.max(maxprio, imm.get(enabled.get(i)).prio);
            }
            List<Integer> top = new ArrayList<Integer>();
            for (int i = 0; i < enabled.size(); i++) {
                if (imm.get(enabled.get(i)).prio == maxprio) top.add(enabled.get(i));
            }
            int pick;
            if (top.size() == 1) {
                pick = top.get(0);
            } else {
                double tot = 0.0;
                for (int i = 0; i < top.size(); i++) tot += imm.get(top.get(i)).weight;
                double u = Maths.rand() * tot;
                double acc = 0.0;
                pick = top.get(top.size() - 1);
                for (int i = 0; i < top.size(); i++) {
                    acc += imm.get(top.get(i)).weight;
                    if (u < acc) { pick = top.get(i); break; }
                }
            }
            SpnRx rf = imm.get(pick);
            for (int s = 0; s < rf.Svec.length; s++) {
                if (rf.Svec[s] != 0.0) nvec.set(s, 0, nvec.get(s, 0) + rf.Svec[s]);
            }
            steps++;
            if (steps > maxsteps) {
                throw new RuntimeException("Immediate-transition livelock: the vanishing-marking collapse did not reach a tangible marking.");
            }
        }
    }

    // ======================================================================
    // Next-Reaction Method (Gibson & Bruck original) with priority queue
    // ======================================================================
    public static void nrm_direct(
            Matrix S,
            List<Integer>[] D,
            BiFunction<Matrix, ArrayDeque<Integer>[], Double>[] a,
            Matrix nvec0,
            ArrayDeque<Integer>[] buffers0,
            int samples,
            SolverOptions options,
            Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
            NetworkStruct sn,
            List<Integer> fromIdx,
            List<int[]> fromIR,
            double[] mi,
            Fcr fcr,
            Balk balk,
            Sig sig,
            Rr rrp,
            boolean[] isRenegeRx,
            boolean[] isRetryRx,
            Smap sm,
            int nDepRx,
            boolean[] isPhaseRx,
            boolean anyPoll,
            Polling.Info[] pinfo,
            boolean[] isPollNode,
            int[][] pollCtrl,
            boolean[] isPollSwRx,
            int[] pollSwNode,
            double[][][] svcph,
            boolean[] bufPHNode,
            boolean[] isBufSvcRx,
            List<Integer> depPhase,
            int[] phaseToArr,
            boolean[] isCacheRx,
            int[] cacheHitSlotArr,
            int[] cacheMissSlotArr,
            boolean[] isCacheNode,
            int[][] cacheContents0,
            long[][] cacheProd,
            int[][] cacheRetrDest) {
        int numReactions = S.getNumCols();
        Matrix nvec = nvec0.copy();
        // Working copy of the per-cache contents (item in each slot), mutated by
        // the cache access at firing.
        int[][] cacheContents = new int[cacheContents0.length][];
        for (int i = 0; i < cacheContents0.length; i++) {
            if (cacheContents0[i] != null) cacheContents[i] = cacheContents0[i].clone();
        }
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] buffers = new ArrayDeque[buffers0.length];
        for (int i = 0; i < buffers0.length; i++) buffers[i] = new ArrayDeque<Integer>(buffers0[i]);

        // Seed each polling controller from the initial per-class populations, so
        // the first propensity evaluation sees a controller in the reachable
        // controller space (Polling.next's landing rule), not an unset one.
        if (anyPoll) {
            for (int ind = 0; ind < pollCtrl.length; ind++) {
                if (!isPollNode[ind]) continue;
                int[] nb = classCountsInt(nvec, sm, ind, sn.nclasses);
                int[] nx = Polling.next(pinfo[ind], 0, nb, sn.nclasses, true);
                pollCtrl[ind] = pollLand(pinfo[ind], nx[0], nx[1], nx[2]);
            }
        }

        int M = sn.nstations;
        int K = sn.nclasses;
        int R = sn.nclasses;
        // Per-region WAITQ FIFO of parked (dstNode, dstClass) tokens, encoded as
        // dstNode*R + dstClass. Empty and untouched unless a region uses WAITQ.
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] fcrBuf = new ArrayDeque[fcr.on ? fcr.F : 0];
        for (int i = 0; i < fcrBuf.length; i++) fcrBuf[i] = new ArrayDeque<Integer>();
        // Reaction-index list used to refresh EVERY reaction after a WAITQ
        // release, which can change populations at arbitrary destination nodes.
        List<Integer> allReactionsIdx = new ArrayList<Integer>(numReactions);
        for (int i = 0; i < numReactions; i++) allReactionsIdx.add(i);
        // Same normalized weights the propensities use, so the PS-family
        // utilization accumulator below applies identical sharing factors.
        double[][] wnorm = buildSchedWeights(sn, R);
        Matrix NK = sn.njobs.transpose();

        java.util.Map<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, MatrixCell>> PH = sn.proc;
        Matrix S_servers = sn.nservers;

        final Routing rt = buildRouting(S);
        final boolean[] jsqReaction = buildJsqFlags(sn, rt, fromIdx, R, sm);
        final int[] sqD = buildSqK(sn, rt, fromIdx, R, sm);

        // Initialise propensities and absolute firing times: tau_k = -ln(U) / a_k
        double[] Ak = new double[numReactions];
        for (int i = 0; i < numReactions; i++) Ak[i] = a[i].apply(nvec, buffers);
        IndexedMinHeap pq = new IndexedMinHeap(numReactions);
        for (int kk = 0; kk < numReactions; kk++) {
            pq.key[kk] = (Ak[kk] > 0.0) ? -Math.log(Maths.rand()) / Ak[kk] : Double.POSITIVE_INFINITY;
        }
        pq.buildHeap();

        double simTime = 0.0;
        double totalTime = 0.0;
        int n = 0;

        long tmoStart = System.nanoTime();
        while (n < samples && !Solver.timeExceeded(tmoStart, options.timeout)) {
            int kfire = pq.peekMin();
            double tauFire = pq.peekMinKey();
            if (tauFire == Double.POSITIVE_INFINITY) break; // absorbing state

            double dt = tauFire - simTime;
            totalTime += dt;

            for (int ist = 0; ist < M; ist++) {
                int ind = (int) sn.stationToNode.get(ist);
                for (int kk = 0; kk < K; kk++) {
                    // nvec counts jobs per phase now, so the class population is
                    // the sum over that class's phases.
                    double currentPop = classPop(nvec, sm, ind, kk);
                    QN.set(ist, kk, QN.get(ist, kk) + currentPop * dt);

                    // Throughput is the total absorption rate of the class: with
                    // phase expansion each phase owns its own departure reaction,
                    // so they are summed. Phase-change reactions move no job and
                    // are excluded, as are the appended renege/retry columns.
                    double depRate = 0.0;
                    for (int idx = 0; idx < nDepRx; idx++) {
                        if (fromIR.get(idx)[0] == ind && fromIR.get(idx)[1] == kk && !isPhaseRx[idx]) {
                            depRate += Ak[idx];
                        }
                    }
                    TN.set(ist, kk, TN.get(ist, kk) + depRate * dt);

                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    if (sched == SchedStrategy.INF || sched == SchedStrategy.EXT) {
                        UN.set(ist, kk, UN.get(ist, kk) + currentPop * dt);
                    } else if (isPsFamily(sched)) {
                        MatrixCell phEntry = phEntryOf(PH, sn, ist, kk);
                        if (phEntry != null && !phEntry.isEmpty()) {
                            UN.set(ist, kk, UN.get(ist, kk) + psFamilyUtil(sn, wnorm[ist],
                                    classCounts(nvec, sm, ind), kk, S_servers.get(ist), sched) * dt);
                        }
                    } else if (isListSched(sched)) {
                        // PAS/OI has no server/buffer split, so its utilization is
                        // the time-average count of in-service jobs per class over
                        // the servers, not a shared-capacity fraction.
                        UN.set(ist, kk, UN.get(ist, kk)
                                + pasUtil(sn, ind, buffers[ind], kk, S_servers.get(ist)) * dt);
                    } else if (isBufferedSched(sched)) {
                        MatrixCell phEntry = phEntryOf(PH, sn, ist, kk);
                        if (phEntry != null && !phEntry.isEmpty()) {
                            int waiting = 0;
                            for (Integer cl : buffers[ind]) {
                                if (cl == kk + 1) waiting++;
                            }
                            double inService = currentPop - waiting;
                            double servers = S_servers.get(ist);
                            UN.set(ist, kk, UN.get(ist, kk) + (inService / servers) * dt);
                        }
                    } else if (sched == SchedStrategy.POLLING) {
                        // The single server is busy on exactly one class-kk job
                        // while the controller serves class kk, and idle (switching
                        // or parked) otherwise; so class-kk utilization is the
                        // fraction of time the controller is SERVING class kk.
                        int[] ctrl = pollCtrl[ind];
                        if (ctrl != null && ctrl[0] == Polling.MODE_VISIT && ctrl[1] == kk) {
                            UN.set(ist, kk, UN.get(ist, kk) + dt / S_servers.get(ist));
                        }
                    }
                }
            }

            simTime = tauFire;
            int destPos;
            boolean cacheChanged = false;
            if (isCacheRx[kfire]) {
                // Cache access. The read-class job reads an item drawn from pread;
                // the cache contents decide a hit or a miss and the replacement
                // policy rewrites them (State.afterEventCache, READ, isSimulation).
                // The job leaves in the hit or miss class at the SAME cache node,
                // and the existing immediate forwarding routes it downstream.
                int cn = fromIR.get(kfire)[0];
                int rdc = fromIR.get(kfire)[1];
                int[] res = cacheAccessNrm(sn, cn, rdc, cacheContents[cn]);
                int outClass = res[0];
                int cat = res[1];
                nvec.set(fromIdx.get(kfire), 0, nvec.get(fromIdx.get(kfire), 0) - 1.0);
                destPos = -1;
                if (outClass >= 0) {
                    if (cat == 4) {
                        // BEGIN retrieval: the job travels to the fetch queue and
                        // returns before the miss completes, so it is placed at the
                        // retrieval class's routed destination, NOT left at the cache.
                        destPos = cacheRetrDest[cn][outClass];
                    } else {
                        // Hit or miss/completion: the job leaves in the hit or miss
                        // class at the SAME cache node; count the production per
                        // produced class (State.afterEventCache convention).
                        destPos = sm.slot(cn, outClass, 0);
                        cacheProd[cn][outClass]++;
                    }
                    nvec.set(destPos, 0, nvec.get(destPos, 0) + 1.0);
                }
                // outClass < 0 is a delayed hit: the request is absorbed (produces
                // nothing), coalescing onto the in-flight retrieval.
                cacheChanged = true;
            } else {
                destPos = fireReaction(kfire, nvec, S, rt, fromIdx, jsqReaction, sqD, R, fcr, fcrBuf, fromIR, balk, sig, rrp, buffers, mi, sm, isPhaseRx, sn);
            }

            if (!isCacheRx[kfire]) {
                maintainBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, isRenegeRx, isRetryRx, sm, svcph, bufPHNode, isBufSvcRx, depPhase, isPhaseRx, phaseToArr);
            }

            // A firing that touches a buffered-PH node changes svcph, which the
            // static dependency set does not cover (svcph is not in the
            // stoichiometry), so force a full refresh -- exactly as a polling
            // controller move or a WAITQ release does.
            boolean svcChanged = bufPHNode[fromIR.get(kfire)[0]] || (destPos >= 0 && bufPHNode[sm.node[destPos]]);

            // Advance any polling controller the firing moved. A controller move
            // flips service gates or the switchover rate, i.e. changes rates
            // outside the static dependency set of the fired reaction.
            boolean pollChanged = anyPoll && pollAdvance(kfire, destPos, nvec, sm, R,
                    pinfo, isPollNode, pollCtrl, isPollSwRx, pollSwNode, fromIR, nDepRx, isPhaseRx);

            // WAITQ: admit parked jobs whose regions this firing may have
            // relieved. A release changes populations at arbitrary destination
            // nodes.
            int nReleased = 0;
            if (fcr.on && fcr.anyWaitq) {
                nReleased = fcrReleaseCascade(fcr, nvec, buffers, fcrBuf, mi, R, sn, sm, svcph, bufPHNode);
            }

            // Refresh the affected reactions' propensities and putative times.
            // The fired reaction always draws a fresh exponential; every other
            // affected reaction rescales its residual. Either a polling controller
            // move or a WAITQ release changes rates outside the fired reaction's
            // static dependency set, so either widens the refresh set to all
            // reactions rather than only D[kfire].
            Iterable<Integer> refreshSet = (pollChanged || nReleased > 0 || svcChanged || cacheChanged) ? allReactionsIdx : D[kfire];
            for (Integer i : refreshSet) {
                double aOld = Ak[i];
                Ak[i] = a[i].apply(nvec, buffers);
                double aNew = Ak[i];
                double newTau;
                if (i == kfire) {
                    newTau = (aNew > 0.0) ? simTime - Math.log(Maths.rand()) / aNew : Double.POSITIVE_INFINITY;
                } else if (aNew <= 0.0) {
                    newTau = Double.POSITIVE_INFINITY;
                } else if (aOld <= 0.0) {
                    // Reaction was dormant, now active - no residual time to rescale
                    newTau = simTime - Math.log(Maths.rand()) / aNew;
                } else {
                    newTau = (aOld / aNew) * (pq.getKey(i) - simTime) + simTime;
                }
                pq.update(i, newTau);
            }

            n++;
            printProgress(options, n);
        }

        if (totalTime > 0) {
            for (int ist = 0; ist < M; ist++) {
                for (int kk = 0; kk < K; kk++) {
                    QN.set(ist, kk, QN.get(ist, kk) / totalTime);
                    UN.set(ist, kk, UN.get(ist, kk) / totalTime);
                    TN.set(ist, kk, TN.get(ist, kk) / totalTime);
                }
            }
        }

        // Class-dependent stations report utilization as T*S/peak, using the
        // declared per-class peak rate (sn.cdscalingpeak), matching the analytic
        // solvers and serial SSA; the in-service accumulation above divides by
        // the server count (1 for a cd station), a different quantity.
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
                jline.lang.nodes.Station cdStat = sn.stations.get(ist);
                if (sn.cdscaling.get(cdStat) == null) continue;
                Matrix cdPeakVec = (sn.cdscalingpeak != null) ? sn.cdscalingpeak.get(cdStat) : null;
                for (int kk = 0; kk < K; kk++) {
                    double cdRate = sn.rates.get(ist, kk);
                    double cdPk = (cdPeakVec != null) ? cdPeakVec.get(0, kk) : 1.0;
                    if (Double.isFinite(cdRate) && cdRate > 0 && cdPk > 0) {
                        UN.set(ist, kk, TN.get(ist, kk) / cdRate / cdPk);
                    } else {
                        UN.set(ist, kk, 0.0);
                    }
                }
            }
        }

        // Joint-dependent stations report utilization as T*S/peak, using the
        // declared per-class peak rate (sn.jdscalingpeak), mirroring the
        // class-dependence block.
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
                jline.lang.nodes.Station jdStat = sn.stations.get(ist);
                if (sn.jdscaling.get(jdStat) == null) continue;
                Matrix jdPeakVec = (sn.jdscalingpeak != null) ? sn.jdscalingpeak.get(jdStat) : null;
                for (int kk = 0; kk < K; kk++) {
                    double jdRate = sn.rates.get(ist, kk);
                    double jdPk = (jdPeakVec != null) ? jdPeakVec.get(0, kk) : 1.0;
                    if (Double.isFinite(jdRate) && jdRate > 0 && jdPk > 0) {
                        UN.set(ist, kk, TN.get(ist, kk) / jdRate / jdPk);
                    } else {
                        UN.set(ist, kk, 0.0);
                    }
                }
            }
        }

        for (int kk = 0; kk < K; kk++) {
            XN.set(0, kk, TN.get((int) sn.refstat.get(kk), kk));
            for (int ist = 0; ist < M; ist++) {
                if (TN.get(ist, kk) > 0) {
                    RN.set(ist, kk, QN.get(ist, kk) / TN.get(ist, kk));
                } else {
                    RN.set(ist, kk, 0.0);
                }
            }
            if (XN.get(0, kk) > 0) {
                CN.set(0, kk, NK.get(kk) / XN.get(0, kk));
            }
        }

        QN.apply(Double.NaN, 0.0, "equal");
        UN.apply(Double.NaN, 0.0, "equal");
        RN.apply(Double.NaN, 0.0, "equal");
        XN.apply(Double.NaN, 0.0, "equal");
        TN.apply(Double.NaN, 0.0, "equal");
        CN.apply(Double.NaN, 0.0, "equal");
    }

    // 0-indexed offset of 1-indexed position J in 1-indexed list I of a cache.
    private static int cposC(int[] m, int i, int j) {
        int off = 0;
        for (int x = 0; x < i - 1; x++) off += m[x];
        return off + (j - 1);
    }

    // Index (1..len) drawn from an unnormalised nonnegative weight list/array,
    // mirroring MATLAB drawFromDist (normalised cumsum > U).
    private static int drawCumList(List<Double> p) {
        double tot = 0.0;
        for (int i = 0; i < p.size(); i++) tot += p.get(i);
        if (tot <= 0.0) return 1;
        double r = Maths.rand() * tot, c = 0.0;
        for (int i = 0; i < p.size(); i++) { c += p.get(i); if (c > r) return i + 1; }
        return p.size();
    }

    private static int drawCumArr(double[] w) {
        double tot = 0.0;
        for (int i = 0; i < w.length; i++) tot += w[i];
        if (tot <= 0.0) return 1;
        double r = Maths.rand() * tot, c = 0.0;
        for (int i = 0; i < w.length; i++) { c += w[i]; if (c > r) return i + 1; }
        return w.length;
    }

    private static double[] cacheRow(Matrix mat, int row) {
        double[] out = new double[mat.getNumCols()];
        for (int j = 0; j < out.length; j++) out[j] = mat.get(row, j);
        return out;
    }

    /**
     * Simulate one cache READ at cache node IND by class CLS over the cache
     * CONTENTS (item id in each of totalCacheCapacity slots, laid out list by
     * list). Draws the requested item from pread and the miss-insertion list from
     * the access-cost row, rewrites CONTENTS per the replacement policy, and
     * returns whether it was a hit. A faithful port of State.afterEventCache
     * (READ, isSimulation) for the non-retrieval case; reads from a snapshot and
     * writes to CONTENTS to mirror MATLAB's varp/var copy semantics exactly.
     */
    private static int[] cacheAccessNrm(NetworkStruct sn, int ind, int cls, int[] contents) {
        jline.lang.nodeparam.CacheNodeParam np =
                (jline.lang.nodeparam.CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
        Matrix mm = np.itemcap;
        int h = mm.getNumCols();
        int[] m = new int[h];
        for (int i = 0; i < h; i++) m[i] = (int) mm.get(i);
        int tcc = np.totalCacheCapacity;
        jline.lang.constant.ReplacementStrategy rep = np.replacestrat;
        double qadm = np.qlru;
        int hitClass = (cls < np.hitclass.length()) ? (int) np.hitclass.get(cls) : 0;
        int missClass = (cls < np.missclass.length()) ? (int) np.missclass.get(cls) : 0;
        boolean isFromRetrieval = np.retrievalClassIndices != null
                && np.retrievalClassIndices.contains(cls);
        boolean hasRetrieval = np.retrievalSystemCapacity > 0;
        List<Double> p = np.pread.get(cls);
        int k = drawCumList(p);                                  // requested item id 1..n
        int l = drawCumArr(cacheRow(np.accost[cls][k - 1], 0));  // target list for a miss (1 => reject)
        int posk = 0;
        for (int c = 0; c < tcc; c++) { if (contents[c] == k) { posk = c + 1; break; } }
        if (isFromRetrieval) posk = 0;   // a returning retrieval always COMPLETES its miss
        int[] orig = contents.clone();

        if (posk == 0) {
            // CACHE MISS / retrieval. Consult the retrieval system first: an item
            // with a retrieval class is fetched rather than admitted directly.
            if (hasRetrieval && !isFromRetrieval) {
                int rClass = -1;
                if (np.retrievalClasses != null && (k - 1) < np.retrievalClasses.getNumRows()
                        && cls < np.retrievalClasses.getNumCols()) {
                    rClass = (int) np.retrievalClasses.get(k - 1, cls);
                }
                if (rClass != -1) {
                    boolean inRetrieval = (tcc + (k - 1) < contents.length) && contents[tcc + (k - 1)] != 0;
                    if (inRetrieval) {
                        // DELAYED HIT: served by the in-flight retrieval, absorbed
                        // (produces nothing). outClass -1 is the absorb sentinel
                        // (class 0 is a valid class in the JAR's 0-based indexing).
                        return new int[]{-1, 3};
                    } else {
                        // BEGIN retrieval: switch to the item's retrieval class and
                        // mark the item as being fetched.
                        contents[tcc + (k - 1)] = 1;
                        return new int[]{rClass, 4};
                    }
                }
            }
            // COMPLETE the miss: returning retrieval, or a plain miss with no
            // retrieval class. Clear the retrieval bit (if any) and admit item k.
            if (isFromRetrieval && (tcc + (k - 1) < contents.length)) contents[tcc + (k - 1)] = 0;
            int listidx = l - 1;
            if (rep == jline.lang.constant.ReplacementStrategy.FIFO
                    || rep == jline.lang.constant.ReplacementStrategy.LRU
                    || rep == jline.lang.constant.ReplacementStrategy.SFIFO
                    || rep == jline.lang.constant.ReplacementStrategy.HLRU) {
                if (listidx > 0) {
                    for (int x = 2; x <= m[listidx - 1]; x++) contents[cposC(m, listidx, x)] = orig[cposC(m, listidx, x - 1)];
                    contents[cposC(m, listidx, 1)] = k;
                }
            } else if (rep == jline.lang.constant.ReplacementStrategy.RR) {
                if (listidx > 0) {
                    int rpos = 1 + (int) (Maths.rand() * m[listidx - 1]);
                    contents[cposC(m, listidx, rpos)] = k;
                }
            } else if (rep == jline.lang.constant.ReplacementStrategy.QLRU) {
                if (listidx > 0 && Maths.rand() <= qadm) {
                    for (int x = 2; x <= m[listidx - 1]; x++) contents[cposC(m, listidx, x)] = orig[cposC(m, listidx, x - 1)];
                    contents[cposC(m, listidx, 1)] = k;
                }
            }
            return new int[]{missClass, 2};
        }

        // Locate the 1-indexed list i and position j of posk.
        int cum = 0, i = 1;
        for (; i <= h; i++) { if (posk <= cum + m[i - 1]) break; cum += m[i - 1]; }
        int j = posk - cum;

        if (i < h) {
            // CACHE HIT in list i < h: promote toward the last list.
            double[] accr = cacheRow(np.accost[cls][k - 1], i);      // MATLAB row 1+i
            double[] sub = java.util.Arrays.copyOfRange(accr, i, accr.length); // cols (1+i):end
            int inew = i + drawCumArr(sub) - 1;                     // target list, can stay at i
            if (rep == jline.lang.constant.ReplacementStrategy.FIFO) {
                if (inew != i) {
                    contents[cposC(m, i, j)] = orig[cposC(m, inew, m[inew - 1])];
                    for (int x = 2; x <= m[inew - 1]; x++) contents[cposC(m, inew, x)] = orig[cposC(m, inew, x - 1)];
                    contents[cposC(m, inew, 1)] = k;
                }
            } else if (rep == jline.lang.constant.ReplacementStrategy.RR) {
                int rpos = 1 + (int) (Maths.rand() * m[inew - 1]);
                contents[cposC(m, i, j)] = orig[cposC(m, inew, rpos)];
                contents[cposC(m, inew, rpos)] = k;
            } else if (rep == jline.lang.constant.ReplacementStrategy.LRU
                    || rep == jline.lang.constant.ReplacementStrategy.SFIFO
                    || rep == jline.lang.constant.ReplacementStrategy.HLRU
                    || rep == jline.lang.constant.ReplacementStrategy.QLRU) {
                for (int x = 2; x <= j; x++) contents[cposC(m, i, x)] = orig[cposC(m, i, x - 1)];
                contents[cposC(m, i, 1)] = orig[cposC(m, inew, m[inew - 1])];
                for (int x = 2; x <= m[inew - 1]; x++) contents[cposC(m, inew, x)] = orig[cposC(m, inew, x - 1)];
                contents[cposC(m, inew, 1)] = k;
            }
        } else {
            // CACHE HIT in the last list h.
            if (rep == jline.lang.constant.ReplacementStrategy.LRU
                    || rep == jline.lang.constant.ReplacementStrategy.HLRU
                    || rep == jline.lang.constant.ReplacementStrategy.QLRU) {
                for (int x = 2; x <= j; x++) contents[cposC(m, h, x)] = orig[cposC(m, h, x - 1)];
                contents[cposC(m, h, 1)] = orig[cposC(m, h, j)];
            }
            // RR/FIFO/SFIFO: no reordering.
        }
        return new int[]{hitClass, 1};
    }

    // ======================================================================
    // Anderson's Modified Next-Reaction Method with direct metric computation
    // ======================================================================
    public static void modified_nrm_direct(
            Matrix S,
            List<Integer>[] D,
            BiFunction<Matrix, ArrayDeque<Integer>[], Double>[] a,
            Matrix nvec0,
            ArrayDeque<Integer>[] buffers0,
            int samples,
            SolverOptions options,
            Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
            NetworkStruct sn,
            List<Integer> fromIdx,
            List<int[]> fromIR,
            double[] mi,
            Fcr fcr,
            Balk balk,
            Sig sig,
            Rr rrp,
            boolean[] isRenegeRx,
            boolean[] isRetryRx,
            Smap sm,
            int nDepRx,
            boolean[] isPhaseRx,
            double[][][] svcph,
            boolean[] bufPHNode,
            boolean[] isBufSvcRx,
            List<Integer> depPhase,
            int[] phaseToArr) {
        int numReactions = S.getNumCols();
        double[] Ak = new double[numReactions];
        for (int i = 0; i < numReactions; i++) Ak[i] = a[i].apply(nvec0, buffers0);
        double[] Pk = new double[numReactions];
        for (int i = 0; i < numReactions; i++) Pk[i] = -Math.log(Maths.rand());
        double[] Tk = new double[numReactions];

        Matrix nvec = nvec0.copy();
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] buffers = new ArrayDeque[buffers0.length];
        for (int i = 0; i < buffers0.length; i++) buffers[i] = new ArrayDeque<Integer>(buffers0[i]);
        // Per-region WAITQ FIFO of parked (dstNode, dstClass) tokens.
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] fcrBuf = new ArrayDeque[fcr.on ? fcr.F : 0];
        for (int i = 0; i < fcrBuf.length; i++) fcrBuf[i] = new ArrayDeque<Integer>();

        int M = sn.nstations;
        int K = sn.nclasses;
        int R = sn.nclasses;
        // Same normalized weights the propensities use, so the PS-family
        // utilization accumulator below applies identical sharing factors.
        double[][] wnorm = buildSchedWeights(sn, R);
        Matrix NK = sn.njobs.transpose();

        double totalTime = 0.0;
        java.util.Map<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, MatrixCell>> PH = sn.proc;
        Matrix S_servers = sn.nservers;

        final Routing rt = buildRouting(S);
        final boolean[] jsqReaction = buildJsqFlags(sn, rt, fromIdx, R, sm);
        final int[] sqD = buildSqK(sn, rt, fromIdx, R, sm);

        int n = 0;
        long tmoStart = System.nanoTime();
        while (n < samples && !Solver.timeExceeded(tmoStart, options.timeout)) {
            double[] tau = new double[numReactions];
            for (int i = 0; i < numReactions; i++) {
                tau[i] = (Ak[i] > 0) ? (Pk[i] - Tk[i]) / Ak[i] : GlobalConstants.Inf;
            }
            int kfire = 0;
            double minTau = tau[0];
            for (int i = 1; i < numReactions; i++) {
                if (tau[i] < minTau) {
                    minTau = tau[i];
                    kfire = i;
                }
            }
            double dt = tau[kfire];
            totalTime += dt;

            for (int ist = 0; ist < M; ist++) {
                int ind = (int) sn.stationToNode.get(ist);
                for (int kk = 0; kk < K; kk++) {
                    // nvec counts jobs per phase now, so the class population is
                    // the sum over that class's phases.
                    double currentPop = classPop(nvec, sm, ind, kk);
                    QN.set(ist, kk, QN.get(ist, kk) + currentPop * dt);

                    // Throughput is the total absorption rate of the class: with
                    // phase expansion each phase owns its own departure reaction,
                    // so they are summed. Phase-change reactions move no job and
                    // are excluded, as are the appended renege/retry columns.
                    double depRate = 0.0;
                    for (int idx = 0; idx < nDepRx; idx++) {
                        if (fromIR.get(idx)[0] == ind && fromIR.get(idx)[1] == kk && !isPhaseRx[idx]) {
                            depRate += Ak[idx];
                        }
                    }
                    TN.set(ist, kk, TN.get(ist, kk) + depRate * dt);

                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    if (sched == SchedStrategy.INF || sched == SchedStrategy.EXT) {
                        UN.set(ist, kk, UN.get(ist, kk) + currentPop * dt);
                    } else if (isPsFamily(sched)) {
                        MatrixCell phEntry = phEntryOf(PH, sn, ist, kk);
                        if (phEntry != null && !phEntry.isEmpty()) {
                            UN.set(ist, kk, UN.get(ist, kk) + psFamilyUtil(sn, wnorm[ist],
                                    classCounts(nvec, sm, ind), kk, S_servers.get(ist), sched) * dt);
                        }
                    } else if (isListSched(sched)) {
                        // PAS/OI has no server/buffer split, so its utilization is
                        // the time-average count of in-service jobs per class over
                        // the servers, not a shared-capacity fraction.
                        UN.set(ist, kk, UN.get(ist, kk)
                                + pasUtil(sn, ind, buffers[ind], kk, S_servers.get(ist)) * dt);
                    } else if (isBufferedSched(sched)) {
                        MatrixCell phEntry = phEntryOf(PH, sn, ist, kk);
                        if (phEntry != null && !phEntry.isEmpty()) {
                            int waiting = 0;
                            for (Integer cl : buffers[ind]) {
                                if (cl == kk + 1) waiting++;
                            }
                            double inService = currentPop - waiting;
                            double servers = S_servers.get(ist);
                            UN.set(ist, kk, UN.get(ist, kk) + (inService / servers) * dt);
                        }
                    }
                }
            }

            int destPos = fireReaction(kfire, nvec, S, rt, fromIdx, jsqReaction, sqD, R, fcr, fcrBuf, fromIR, balk, sig, rrp, buffers, mi, sm, isPhaseRx, sn);

            maintainBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, isRenegeRx, isRetryRx, sm, svcph, bufPHNode, isBufSvcRx, depPhase, isPhaseRx, phaseToArr);

            boolean svcChanged = bufPHNode[fromIR.get(kfire)[0]] || (destPos >= 0 && bufPHNode[sm.node[destPos]]);

            int nReleased = 0;
            if (fcr.on && fcr.anyWaitq) {
                nReleased = fcrReleaseCascade(fcr, nvec, buffers, fcrBuf, mi, R, sn, sm, svcph, bufPHNode);
            }

            for (int i = 0; i < numReactions; i++) Tk[i] += Ak[i] * dt;
            if (nReleased > 0 || svcChanged) {
                for (int i = 0; i < numReactions; i++) Ak[i] = a[i].apply(nvec, buffers);
            } else {
                for (Integer i : D[kfire]) Ak[i] = a[i].apply(nvec, buffers);
            }
            Pk[kfire] -= Math.log(Maths.rand());

            n++;
            printProgress(options, n);
        }

        if (totalTime > 0) {
            for (int ist = 0; ist < M; ist++) {
                for (int kk = 0; kk < K; kk++) {
                    QN.set(ist, kk, QN.get(ist, kk) / totalTime);
                    UN.set(ist, kk, UN.get(ist, kk) / totalTime);
                    TN.set(ist, kk, TN.get(ist, kk) / totalTime);
                }
            }
        }

        // Class-dependent stations report utilization as T*S/peak, using the
        // declared per-class peak rate (sn.cdscalingpeak), matching the analytic
        // solvers and serial SSA; the in-service accumulation above divides by
        // the server count (1 for a cd station), a different quantity.
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
                jline.lang.nodes.Station cdStat = sn.stations.get(ist);
                if (sn.cdscaling.get(cdStat) == null) continue;
                Matrix cdPeakVec = (sn.cdscalingpeak != null) ? sn.cdscalingpeak.get(cdStat) : null;
                for (int kk = 0; kk < K; kk++) {
                    double cdRate = sn.rates.get(ist, kk);
                    double cdPk = (cdPeakVec != null) ? cdPeakVec.get(0, kk) : 1.0;
                    if (Double.isFinite(cdRate) && cdRate > 0 && cdPk > 0) {
                        UN.set(ist, kk, TN.get(ist, kk) / cdRate / cdPk);
                    } else {
                        UN.set(ist, kk, 0.0);
                    }
                }
            }
        }

        // Joint-dependent stations report utilization as T*S/peak, using the
        // declared per-class peak rate (sn.jdscalingpeak), mirroring the
        // class-dependence block.
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
                jline.lang.nodes.Station jdStat = sn.stations.get(ist);
                if (sn.jdscaling.get(jdStat) == null) continue;
                Matrix jdPeakVec = (sn.jdscalingpeak != null) ? sn.jdscalingpeak.get(jdStat) : null;
                for (int kk = 0; kk < K; kk++) {
                    double jdRate = sn.rates.get(ist, kk);
                    double jdPk = (jdPeakVec != null) ? jdPeakVec.get(0, kk) : 1.0;
                    if (Double.isFinite(jdRate) && jdRate > 0 && jdPk > 0) {
                        UN.set(ist, kk, TN.get(ist, kk) / jdRate / jdPk);
                    } else {
                        UN.set(ist, kk, 0.0);
                    }
                }
            }
        }

        for (int kk = 0; kk < K; kk++) {
            XN.set(0, kk, TN.get((int) sn.refstat.get(kk), kk));
            for (int ist = 0; ist < M; ist++) {
                if (TN.get(ist, kk) > 0) {
                    RN.set(ist, kk, QN.get(ist, kk) / TN.get(ist, kk));
                } else {
                    RN.set(ist, kk, 0.0);
                }
            }
            if (XN.get(0, kk) > 0) {
                CN.set(0, kk, NK.get(kk) / XN.get(0, kk));
            }
        }

        QN.apply(Double.NaN, 0.0, "equal");
        UN.apply(Double.NaN, 0.0, "equal");
        RN.apply(Double.NaN, 0.0, "equal");
        XN.apply(Double.NaN, 0.0, "equal");
        TN.apply(Double.NaN, 0.0, "equal");
        CN.apply(Double.NaN, 0.0, "equal");
    }

    private static MatrixCell phEntryOf(
            java.util.Map<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, MatrixCell>> PH,
            NetworkStruct sn, int ist, int kk) {
        java.util.Map<jline.lang.JobClass, MatrixCell> stProc = PH.get(sn.stations.get(ist));
        return (stProc != null) ? stProc.get(sn.jobclasses.get(kk)) : null;
    }

    /**
     * Buffer maintenance for the firing of reaction kfire. A renege takes a
     * separate path: it removes a job that was WAITING, so no server is freed and
     * no queued job is promoted; the abandoning job simply leaves the buffer.
     * State.afterEventStation drops the newest waiting job of the class and notes
     * that for memoryless patience all waiting jobs are exchangeable, so the
     * choice cannot affect the marginal distribution.
     */
    private static void maintainBuffers(
            int kfire,
            Matrix nvec,
            ArrayDeque<Integer>[] buffers,
            List<int[]> fromIR,
            int destPos,
            double[] mi,
            int R,
            NetworkStruct sn,
            boolean[] isRenegeRx,
            boolean[] isRetryRx,
            Smap sm,
            double[][][] svcph,
            boolean[] bufPHNode,
            boolean[] isBufSvcRx,
            List<Integer> depPhase,
            boolean[] isPhaseRx,
            int[] phaseToArr) {
        if (kfire < isRetryRx.length && isRetryRx[kfire]) {
            // A successful retry moves one orbiting job into the free server. The
            // population is unchanged (it was already counted at the station), so
            // only the orbit shrinks; in-service is read back as population minus
            // orbit occupancy. For memoryless retrials the orbiting jobs of a
            // class are exchangeable, so which one leaves cannot matter.
            buffers[fromIR.get(kfire)[0]].removeFirstOccurrence(Integer.valueOf(fromIR.get(kfire)[1] + 1));
            return;
        }
        if (kfire < isRenegeRx.length && isRenegeRx[kfire]) {
            // buffers are ordered newest-first, so the first occurrence is the
            // newest waiting job of the class
            buffers[fromIR.get(kfire)[0]].removeFirstOccurrence(Integer.valueOf(fromIR.get(kfire)[1] + 1));
            return;
        }
        if (kfire < isPhaseRx.length && isPhaseRx[kfire] && bufPHNode[fromIR.get(kfire)[0]]) {
            // A buffered-PH phase transition moves one in-service job between phases
            // of its own service process. It frees no server and adds no arrival, so
            // the buffer is untouched and only svcph changes (INF/PS phase moves are
            // already applied to nvec via the stoichiometry and fall to updateBuffers
            // below as a no-op, as before).
            int ind = fromIR.get(kfire)[0];
            int r = fromIR.get(kfire)[1];
            svcph[ind][r][depPhase.get(kfire)] -= 1.0;
            svcph[ind][r][phaseToArr[kfire]] += 1.0;
            return;
        }
        updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, sm, svcph, bufPHNode, isBufSvcRx, depPhase);
    }

    private static void updateBuffers(
            int kfire,
            Matrix nvec,
            ArrayDeque<Integer>[] buffers,
            List<int[]> fromIR,
            int destPos,
            double[] mi,
            int R,
            NetworkStruct sn,
            Smap sm,
            double[][][] svcph,
            boolean[] bufPHNode,
            boolean[] isBufSvcRx,
            List<Integer> depPhase) {
        int ind = fromIR.get(kfire)[0];

        // Buffered-PH departure: the completing job leaves service, so drop it from
        // the in-service phase it occupied (carried in depPhase). The promotion
        // below refills the freed server from the buffer at a fresh entry phase.
        if (bufPHNode[ind] && kfire < isBufSvcRx.length && isBufSvcRx[kfire]) {
            svcph[ind][fromIR.get(kfire)[1]][depPhase.get(kfire)] -= 1.0;
        }

        // A pass-and-swap station keeps the full ordered list, so a departure is
        // not a promotion but a rewrite: the completing position's chain shifts
        // classes along and removes one slot.
        if (isListSched(ind, sn) && !buffers[ind].isEmpty()) {
            oiDepart(sn, ind, buffers[ind], fromIR.get(kfire)[1]);
            return;
        }

        // Handle departure from a buffered source node: the freed server takes
        // the waiting job selected by the station's discipline. A retrial station
        // is the exception -- the freed server is NOT filled from the orbit,
        // orbiting jobs re-enter only through RETRY events at the memoryless
        // retrial rate (State.afterEventStation suppresses promotion likewise).
        if (isBuffered(ind, sn) && !buffers[ind].isEmpty() && !isRetrialStation(ind, sn)) {
            int pos = pickFromBuffer(buffers[ind], sn, (int) sn.nodeToStation.get(ind));
            int promoted = elemAt(buffers[ind], pos) - 1;   // class id -> 0-based class
            removeAt(buffers[ind], pos);
            if (bufPHNode[ind] && promoted >= 0) {
                // The promoted waiting job starts service now, entering a phase drawn
                // from its entry distribution pie (the same allocation the init uses).
                int ke = drawEntryPhase(sn, sm, ind, promoted);
                svcph[ind][promoted][ke] += 1.0;
            }
        }

        // Destination is the (state-row) position selected by the routing draw
        if (destPos < 0) {
            return;
        }
        applyArrivalBuffer(sm.node[destPos], sm.cls[destPos], nvec, buffers, mi, R, sn, sm, svcph, bufPHNode);
    }

    /**
     * Join a just-arrived class-r job to the ordered buffer of destination node
     * jnd, if that node is buffered. nvec already includes the arrival. Shared by
     * {@link #updateBuffers} (routed arrivals) and {@link #fcrReleaseCascade}
     * (WAITQ releases), so the two paths cannot drift. Mirrors MATLAB
     * applyArrivalBuffer.
     */
    private static void applyArrivalBuffer(int jnd, int r, Matrix nvec,
                                           ArrayDeque<Integer>[] buffers, double[] mi,
                                           int R, NetworkStruct sn, Smap sm,
                                           double[][][] svcph, boolean[] bufPHNode) {
        if (isListSched(jnd, sn)) {
            // PAS/OI: the arrival simply joins the back of the ordered list; there
            // is no server/buffer split, so no capacity test against mi. Capacity
            // is the station's own cap, and an arrival past it is lost.
            if (buffers[jnd].size() < sn.cap.get((int) sn.nodeToStation.get(jnd))) {
                buffers[jnd].addLast(r + 1); // append at the back (newest last)
            }
            return;
        }

        // Handle arrival at buffered destination node
        if (isBuffered(jnd, sn)) {
            double totalAtDest = 0.0;
            for (int i = 0; i < R; i++) totalAtDest += classPop(nvec, sm, jnd, i);
            boolean enteredService = false;
            if (isRetrialStation(jnd, sn)) {
                // A retrial station breaks the buffer invariant the other policies
                // share: because a departure does not promote, the orbit can be
                // occupied while servers sit idle, so "total > mi" no longer means
                // "the servers are busy". An arrival must consult the servers
                // directly and only join the orbit when none is free.
                double inSvc = (totalAtDest - 1.0) - buffers[jnd].size();
                if (inSvc >= mi[jnd]) {
                    buffers[jnd].addFirst(r + 1);
                } else {
                    enteredService = true;
                }
            } else if (totalAtDest > mi[jnd]) {
                if (isPreemptive(jnd, sn)) {
                    // Preempt-resume: the arrival seizes a server and the incumbent
                    // it displaces is the one that joins the buffer. Buffering the
                    // incumbent rather than the arrival is what leaves the new job
                    // in service, since in-service is read back as population minus
                    // buffer occupancy. Promotion is newest-first, so the most
                    // recently preempted job is the one that resumes.
                    int c = pickPreempted(nvec, buffers[jnd], jnd, r, R, sm);
                    if (c >= 0) {
                        buffers[jnd].addFirst(c + 1); // addFirst
                    }
                    enteredService = true;
                } else {
                    // All servers busy - arriving job joins the buffer
                    buffers[jnd].addFirst(r + 1);
                }
            } else {
                // A server is free: the job goes straight into service.
                enteredService = true;
            }
            if (enteredService && bufPHNode[jnd]) {
                int ke = drawEntryPhase(sn, sm, jnd, r);
                svcph[jnd][r][ke] += 1.0;
            }
        }
    }

    /**
     * Sample the service phase a class-r job starts in at node jnd from its entry
     * distribution pie. A single-phase class always enters phase 0.
     */
    private static int drawEntryPhase(NetworkStruct sn, Smap sm, int jnd, int r) {
        int np = sm.nph[jnd][r];
        if (np <= 1) return 0;
        double[] pe = entryProbs(sn, sm, jnd, r);
        double tot = 0.0;
        for (int i = 0; i < np; i++) tot += pe[i];
        if (tot <= 0.0) return 0;
        double u = Maths.rand() * tot;
        double acc = 0.0;
        for (int i = 0; i < np; i++) {
            acc += pe[i];
            if (u <= acc) return i;
        }
        return np - 1;
    }

    /** Class id at ordered position idx (0 = newest / head) of a buffer deque. */
    private static int elemAt(ArrayDeque<Integer> dq, int idx) {
        int i = 0;
        for (Integer v : dq) {
            if (i == idx) return v;
            i++;
        }
        return -1;
    }

    // Per-reaction routing for inverse-CDF destination sampling. Mirrors MATLAB
    // next_reaction_method_direct: P = S with negative entries incremented by 1,
    // so each firing moves exactly one job to a single stochastically chosen
    // destination (keeping the marginal state integer). Without this, a fractional
    // routing column would be applied in full, draining immediate pass-through
    // nodes (e.g. ClassSwitch) below zero and biasing station queue lengths.
    private static final class Routing {
        final int[] nnzP;       // number of candidate destinations per reaction
        final int[][] destRow;  // destination state-row indices per reaction
        final double[][] cdf;   // cumulative routing probabilities per reaction
        Routing(int n) { nnzP = new int[n]; destRow = new int[n][]; cdf = new double[n][]; }
    }

    private static Routing buildRouting(Matrix S) {
        int numReactions = S.getNumCols();
        int rows = S.getNumRows();
        Routing rt = new Routing(numReactions);
        for (int k = 0; k < numReactions; k++) {
            List<Integer> dest = new ArrayList<Integer>();
            List<Double> pr = new ArrayList<Double>();
            for (int row = 0; row < rows; row++) {
                double v = S.get(row, k);
                double p = (v < 0.0) ? v + 1.0 : v; // P = S; P(P<0) += 1
                if (p > 0.0) { dest.add(row); pr.add(p); }
            }
            int nd = dest.size();
            rt.nnzP[k] = nd;
            rt.destRow[k] = new int[nd];
            rt.cdf[k] = new double[nd];
            double cum = 0.0;
            for (int x = 0; x < nd; x++) {
                rt.destRow[k][x] = dest.get(x);
                cum += pr.get(x);
                rt.cdf[k][x] = cum;
            }
        }
        return rt;
    }

    // Apply one reaction firing to the marginal state, sampling a single
    // destination via inverse-CDF when the reaction has multiple candidates.
    // Returns the destination state-row index (-1 for self-loop / no move).
    /**
     * Flags reactions whose source class routes with JSQ. Such reactions join
     * the candidate whose destination node holds the smallest total population
     * at firing time (each candidate evaluated on its own queue, never the
     * routing node's; ties split uniformly) instead of inverse-CDF sampling.
     */
    private static boolean[] buildJsqFlags(NetworkStruct sn, Routing rt, List<Integer> fromIdx, int R, Smap sm) {
        boolean[] flags = new boolean[rt.nnzP.length];
        if (sn.routing == null) return flags;
        for (int k = 0; k < rt.nnzP.length; k++) {
            if (rt.nnzP[k] > 1) {
                int ind = sm.node[fromIdx.get(k)];
                int r = sm.cls[fromIdx.get(k)];
                Map<jline.lang.JobClass, jline.lang.constant.RoutingStrategy> rmap = sn.routing.get(sn.nodes.get(ind));
                if (rmap != null && rmap.get(sn.jobclasses.get(r)) == jline.lang.constant.RoutingStrategy.JSQ) {
                    flags[k] = true;
                }
            }
        }
        return flags;
    }

    /**
     * Number of candidates each SQ-routed reaction draws, or 0 when the
     * reaction does not route with SQ.
     * <p>
     * Power-of-k choices samples k candidates uniformly WITH replacement and
     * joins the one holding the smallest total population, ties broken by first
     * occurrence in the sampled tuple. This is the sampled form of the marginal
     * enumerated by sub_sq in MNetwork.refreshRoutingMatrix and of LDES's
     * selectSQDestination; drawing directly is equivalent for a
     * simulator and avoids enumerating the ndest^d tuples. The m>0 variant
     * forces the previous pick as the last candidate, which needs a
     * per-(node,class) memory the reaction network does not carry, so those
     * models are routed to the serial engine by Solver_ssa_analyzer.
     * </p>
     */
    private static int[] buildSqK(NetworkStruct sn, Routing rt, List<Integer> fromIdx, int R, Smap sm) {
        int[] kk = new int[rt.nnzP.length];
        if (sn.routing == null) return kk;
        for (int k = 0; k < rt.nnzP.length; k++) {
            if (rt.nnzP[k] <= 1) continue;
            int ind = sm.node[fromIdx.get(k)];
            int r = sm.cls[fromIdx.get(k)];
            jline.lang.nodes.Node node = sn.nodes.get(ind);
            Map<jline.lang.JobClass, jline.lang.constant.RoutingStrategy> rmap = sn.routing.get(node);
            if (rmap == null
                    || rmap.get(sn.jobclasses.get(r)) != jline.lang.constant.RoutingStrategy.SQ) {
                continue;
            }
            int kval = 2; // sub_sq default when nodeparam carries no d
            if (sn.nodeparam != null && sn.nodeparam.get(node) != null
                    && sn.nodeparam.get(node).d != null) {
                Integer kp = sn.nodeparam.get(node).d.get(sn.jobclasses.get(r));
                if (kp != null) kval = kp;
            }
            kk[k] = Math.max(1, Math.min(kval, rt.destRow[k].length));
        }
        return kk;
    }

    private static int fireReaction(int kfire, Matrix nvec, Matrix S, Routing rt,
                                    List<Integer> fromIdx, boolean[] jsqReaction,
                                    int[] sqD, int R,
                                    Fcr fcr, ArrayDeque<Integer>[] fcrBuf,
                                    List<int[]> fromIR, Balk balk,
                                    Sig sig, Rr rrp, ArrayDeque<Integer>[] buffers, double[] mi,
                                    Smap sm, boolean[] isPhaseRx, NetworkStruct sn) {
        // A phase transition moves a job between the phases of its own service
        // process. It never leaves the node, so it bypasses the routing/balking/
        // signal machinery entirely and just applies its stoichiometry column.
        if (kfire < isPhaseRx.length && isPhaseRx[kfire]) {
            nvec.addEq(S.getColumn(kfire));
            return -1;
        }
        if (rt.nnzP[kfire] > 1) {
            int sel;
            int nCandAll = rt.destRow[kfire].length;
            // A finite capacity region does NOT filter the routing draw.
            // Routing picks the destination first and the region decides
            // admission at the destination's entry afterwards, dropping the job
            // on refusal; a routing strategy that steered around full regions
            // would be a different (and better-behaved) model than the one
            // SOLVER_SSA and the CTMC implement. The refusal check is applied to
            // the drawn destination below.
            if (jsqReaction[kfire]) {
                // JSQ: join the destination node with the smallest total
                // population (ties split uniformly)
                int nCand = rt.destRow[kfire].length;
                double[] npop = new double[nCand];
                double minPop = Double.POSITIVE_INFINITY;
                for (int x = 0; x < nCand; x++) {
                    int jnd = sm.node[rt.destRow[kfire][x]];
                    double pop = 0.0;
                    for (int rr = 0; rr < R; rr++) {
                        pop += classPop(nvec, sm, jnd, rr);
                    }
                    npop[x] = pop;
                    if (pop < minPop) minPop = pop;
                }
                int nMin = 0;
                for (int x = 0; x < nCand; x++) {
                    if (npop[x] == minPop) nMin++;
                }
                int pick = (nMin > 1) ? (int) (Maths.rand() * nMin) : 0;
                sel = 0;
                for (int x = 0; x < nCand; x++) {
                    if (npop[x] == minPop) {
                        if (pick == 0) { sel = x; break; }
                        pick--;
                    }
                }
            } else if (rrp.on && rrp.isrr[fromIR.get(kfire)[0]][fromIR.get(kfire)[1]]) {
                // Round-robin: ADVANCE the pointer, THEN take the destination it
                // lands on (AfterEventRouter advances on DEP and the routing
                // closure reads state_after). The pointer steers the draw only,
                // so the region gate is applied by the propensity as usual and
                // the pick is not re-tested here.
                int jnd = rrNext(rrp, fromIR.get(kfire)[0], fromIR.get(kfire)[1], R);
                // A phase-type destination contributes ONE candidate per entry
                // phase, each weighted by pentry in the routing matrix. The pointer
                // fixes the NODE; the entry PHASE must still be drawn from pentry
                // among that node's candidates. Taking the first match (phase 0)
                // would enter every job in a fixed phase and bias the service time
                // -- the RROBIN + phase-type residence bug (RUN-10). Sample among
                // the matching candidates in proportion to their routing weights.
                double[] cdRr = rt.cdf[kfire];
                int firstMatch = -1, nMatch = 0;
                double wsum = 0.0;
                for (int x = 0; x < nCandAll; x++) {
                    if (sm.node[rt.destRow[kfire][x]] == jnd) {
                        if (firstMatch < 0) firstMatch = x;
                        nMatch++;
                        wsum += cdRr[x] - (x > 0 ? cdRr[x - 1] : 0.0);
                    }
                }
                if (firstMatch < 0) {
                    throw new RuntimeException("Round-robin selected node " + jnd
                            + ", which is not a routing destination of node " + fromIR.get(kfire)[0] + ".");
                }
                if (nMatch == 1 || wsum <= 0.0) {
                    sel = firstMatch;
                } else {
                    double uRr = Maths.rand() * wsum;
                    double accRr = 0.0;
                    sel = firstMatch;
                    for (int x = 0; x < nCandAll; x++) {
                        if (sm.node[rt.destRow[kfire][x]] == jnd) {
                            accRr += cdRr[x] - (x > 0 ? cdRr[x - 1] : 0.0);
                            if (accRr > uRr) { sel = x; break; }
                        }
                    }
                }
            } else if (sqD[kfire] > 0) {
                // SQ: draw k candidates uniformly with replacement from
                // the admitted destinations and keep the least loaded; the
                // strict comparison retains the first occurrence, which is the
                // tie rule of sub_sq.
                int nCand = rt.destRow[kfire].length;
                double[] npop = new double[nCand];
                for (int x = 0; x < nCand; x++) {
                    int jnd = sm.node[rt.destRow[kfire][x]];
                    double pop = 0.0;
                    for (int rr = 0; rr < R; rr++) pop += classPop(nvec, sm, jnd, rr);
                    npop[x] = pop;
                }
                sel = 0;
                double bestPop = Double.POSITIVE_INFINITY;
                int draws = Math.min(sqD[kfire], nCand);
                for (int t = 0; t < draws; t++) {
                    int x = (int) (Maths.rand() * nCand);
                    if (npop[x] < bestPop) {
                        bestPop = npop[x];
                        sel = x;
                    }
                }
            } else {
                double u = Maths.rand();
                sel = rt.cdf[kfire].length - 1;
                for (int x = 0; x < rt.cdf[kfire].length; x++) {
                    if (rt.cdf[kfire][x] > u) { sel = x; break; }
                }
            }
            int src = fromIdx.get(kfire);
            int destPos = rt.destRow[kfire][sel];
            // Balking is decided on the pre-arrival population, so it is drawn
            // before the state is updated. A balked job is lost: the source still
            // releases it, the destination never receives it.
            boolean balked = balk.on && balkDraw(balk, nvec, destPos, R, sm);
            // An open arrival at a full physically-capped destination is lost,
            // exactly as a balked one is: the source releases it, the destination
            // never receives it. Mirrors AfterEventStation.handleArv.
            if (!balked && capacityLoss(sn, nvec, destPos, R, sm)) {
                balked = true;
            }
            // A region refuses the drawn destination on the same pre-arrival
            // population. Under DROP the refused job is lost, exactly as a balked
            // one is; under WAITQ it is parked in the refusing region's FIFO and
            // admitted later, head-of-line. Either way the source still departs
            // and the destination is not entered now.
            if (!balked && fcr.on) {
                int dstN = sm.node[destPos];
                int dstC = sm.cls[destPos];
                int fref = fcrRefusingRegion(fcr, nvec, fromIR.get(kfire)[0], fromIR.get(kfire)[1], dstN, dstC, R, sm);
                if (fref >= 0) {
                    balked = true;
                    if (fcr.waitq[fref][dstC]) {
                        fcrBuf[fref].addLast(dstN * R + dstC);
                    }
                }
            }
            nvec.set(src, 0, nvec.get(src, 0) - 1.0);
            if (balked) {
                return -1;
            }
            if (sig.on && sigIsSignalArrival(sig, destPos, R, sm)) {
                // the signal is annihilated on arrival: it never joins the station
                sigApply(sig, nvec, buffers, destPos, R, mi, sm);
                return -1;
            }
            nvec.set(destPos, 0, nvec.get(destPos, 0) + 1.0);
            return destPos;
        }
        if (rt.nnzP[kfire] == 1) {
            int destPos = rt.destRow[kfire][0];
            boolean lost = balk.on && balkDraw(balk, nvec, destPos, R, sm);
            if (!lost && capacityLoss(sn, nvec, destPos, R, sm)) {
                lost = true;
            }
            // Single-destination departures cross region boundaries too, so the
            // region gate applies here exactly as it does to a drawn destination.
            // Renege and retry columns carry no destination and never reach here.
            if (!lost && fcr.on) {
                int dstN = sm.node[destPos];
                int dstC = sm.cls[destPos];
                int fref = fcrRefusingRegion(fcr, nvec, fromIR.get(kfire)[0], fromIR.get(kfire)[1], dstN, dstC, R, sm);
                if (fref >= 0) {
                    lost = true;
                    if (fcr.waitq[fref][dstC]) {
                        fcrBuf[fref].addLast(dstN * R + dstC);
                    }
                }
            }
            if (lost) {
                // lost or parked on arrival: apply the source departure only
                int src = fromIdx.get(kfire);
                nvec.set(src, 0, nvec.get(src, 0) - 1.0);
                return -1;
            }
            if (sig.on && sigIsSignalArrival(sig, destPos, R, sm)) {
                // the signal is annihilated on arrival: it never joins the station
                int src = fromIdx.get(kfire);
                nvec.set(src, 0, nvec.get(src, 0) - 1.0);
                sigApply(sig, nvec, buffers, destPos, R, mi, sm);
                return -1;
            }
            nvec.addEq(S.getColumn(kfire));
            return destPos;
        }
        // no destination: a self-loop, or a renege whose column is a bare -1
        nvec.addEq(S.getColumn(kfire));
        return -1;
    }

    /**
     * Non-preemptive policies that hold waiting jobs in an ordered buffer. They
     * share the rate law (a class-r completion fires at mu_r times the class-r
     * jobs actually in service) and differ only in which waiting job is promoted
     * on a departure; see {@link #pickFromBuffer}.
     */
    static boolean isCountBufferedSched(SchedStrategy sched) {
        // Of the buffered policies, those whose state keeps the buffer as
        // per-class counts rather than an ordered list of class ids
        // (State.fromMarginalAndRunning).
        return sched == SchedStrategy.SIRO || sched == SchedStrategy.SEPT
                || sched == SchedStrategy.LEPT;
    }

    static boolean isBufferedSched(SchedStrategy sched) {
        return sched == SchedStrategy.FCFS || sched == SchedStrategy.LCFS
                || sched == SchedStrategy.SIRO || sched == SchedStrategy.HOL
                || sched == SchedStrategy.SEPT || sched == SchedStrategy.LEPT
                || sched == SchedStrategy.LCFSPR;
    }

    /**
     * Preempt-resume policies: an arrival at a fully busy station takes a server
     * and pushes back into the buffer the incumbent it displaced, rather than
     * queueing itself (State.afterEventStation, the LCFSPR arrival group). The
     * rate law is unchanged -- it still counts the jobs actually in service -- so
     * no extra state is needed: in-service is population minus buffer occupancy,
     * and that automatically names the new arrival as the one being served. With
     * exponential service preempt-resume needs no stored phase, because a resumed
     * job has the same memoryless residual as a fresh one. LCFSPI is deliberately
     * absent: SolverSSA.getFeatureSet does not advertise it (nor does SolverCTMC),
     * so the NRM must not claim it either.
     */
    static boolean isPreemptiveSched(SchedStrategy sched) {
        return sched == SchedStrategy.LCFSPR;
    }

    private static boolean isPreemptive(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0) return false;
        return isPreemptiveSched(sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ind))));
    }

    /**
     * Class of the incumbent displaced by an arrival of class arrClass (0-based)
     * at node jnd, drawn in proportion to the servers' class occupancies, as
     * State.afterEventStation weights its preemption branches by
     * si_preempt/sum(space_srv). Returns -1 when no server is occupied. nvec
     * already counts the arrival, so it is discounted here to recover the
     * pre-arrival in-service composition (in-service = population minus buffer).
     */
    private static int pickPreempted(Matrix nvec, ArrayDeque<Integer> buf, int jnd, int arrClass, int R, Smap sm) {
        double[] insvc = new double[R];
        double tot = 0.0;
        for (int r = 0; r < R; r++) {
            int waiting = 0;
            for (Integer cl : buf) {
                if (cl == r + 1) waiting++;
            }
            insvc[r] = classPop(nvec, sm, jnd, r) - waiting;
            if (r == arrClass) {
                insvc[r] -= 1.0; // discount the job that just arrived
            }
            if (insvc[r] < 0.0) insvc[r] = 0.0;
            tot += insvc[r];
        }
        if (tot <= 0.0) {
            return -1;
        }
        double u = Maths.rand() * tot;
        double acc = 0.0;
        int last = -1;
        for (int r = 0; r < R; r++) {
            if (insvc[r] <= 0.0) continue;
            last = r;
            acc += insvc[r];
            if (u < acc) return r;
        }
        return last;
    }

    /**
     * Stations whose buffer holds the FULL ordered job list rather than only the
     * waiting jobs: pass-and-swap / order-independent. There is no server/buffer
     * split at all and the rate is a function mu(c) of the whole list, so these
     * carry a different buffer invariant -- size == total, rather than
     * max(0, total - mi) -- and the list runs OLDEST-FIRST, the reverse of every
     * other buffered policy here.
     * <p>
     * sn.sched carries PAS for both PAS and OI stations: OI is canonicalized to
     * pass-and-swap with an all-zero swap graph (Network.refreshLocalVars), so
     * reading the graph covers both and OI needs no separate case.
     * </p>
     */
    static boolean isListSched(SchedStrategy sched) {
        return sched == SchedStrategy.PAS;
    }

    private static boolean isListSched(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0) return false;
        return isListSched(sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ind))));
    }

    private static QueueNodeParam pasParam(NetworkStruct sn, int ind) {
        NodeParam np = sn.nodeparam.get(sn.nodes.get(ind));
        QueueNodeParam qp = (np instanceof QueueNodeParam) ? (QueueNodeParam) np : null;
        if (qp == null || qp.svcRateFun == null) {
            throw new RuntimeException(
                    "PAS/OI station has no service rate function mu(c); set it via setService(c -> ...).");
        }
        return qp;
    }

    /** The ordered list as a plain 1-based array, oldest first. */
    private static int[] listOf(ArrayDeque<Integer> buf) {
        int[] c = new int[buf.size()];
        int i = 0;
        for (Integer v : buf) {
            c[i++] = v;
        }
        return c;
    }

    /**
     * Service rate increment of position p (0-based) of the ordered list c:
     * Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}). svcRateFun is called with
     * a row Matrix of 0-based class indices, matching afterEventStationPas.
     */
    private static double[] oiDeltaMu(QueueNodeParam qp, int[] c) {
        double[] d = new double[c.length];
        double muPrev = 0.0; // mu of the empty prefix is 0
        for (int p = 0; p < c.length; p++) {
            Matrix prefix = new Matrix(1, p + 1);
            for (int i = 0; i <= p; i++) prefix.set(0, i, c[i] - 1);
            double muCur = qp.svcRateFun.apply(prefix);
            d[p] = muCur - muPrev;
            muPrev = muCur;
        }
        return d;
    }

    /**
     * Aggregate class-r (0-based) departure rate of a pass-and-swap station
     * holding the ordered list buf. Mirrors the DEP branch of
     * AfterEventStation.afterEventStationPas: every position contributes its own
     * service token at Delta_mu, and pass-and-swap decides which class actually
     * leaves, so the class-r rate is the total Delta_mu over the positions whose
     * pass-and-swap ejects a class-r job.
     */
    static double oirate(NetworkStruct sn, int ind, ArrayDeque<Integer> buf, int r) {
        if (buf.isEmpty()) {
            return 0.0;
        }
        QueueNodeParam qp = pasParam(sn, ind);
        int[] c = listOf(buf);
        double[] d = oiDeltaMu(qp, c);
        double rt = 0.0;
        for (int p = 0; p < c.length; p++) {
            if (d[p] <= 0.0) continue; // position p receives no service
            if (AfterEventStation.passAndSwap(c, p, qp.swapGraph).depClass == r + 1) {
                rt += d[p];
            }
        }
        return rt;
    }

    /**
     * Class-k (0-based) utilization contribution of a pass-and-swap station
     * holding the ordered list buf: the number of class-k jobs in service over
     * the server count. "In service" means the positions whose marginal rate
     * increment Delta_mu is POSITIVE, so a job served by several server types
     * still counts once rather than 1/rate. This is the sir that the PAS branch
     * of ToMarginal.toMarginalAggr reports, which is what solver_ctmc_analyzer
     * divides by the server count -- the PS-family share formula does not apply
     * here, because a PAS station has no server/buffer split to share out.
     */
    static double pasUtil(NetworkStruct sn, int ind, ArrayDeque<Integer> buf, int k, double servers) {
        return pasInSvc(sn, ind, buf, k) / servers;
    }

    /**
     * Number of class-r (0-based) jobs in service at a PAS/OI station holding the
     * ordered list buf: the positions with a positive Delta_mu.
     */
    static int pasInSvc(NetworkStruct sn, int ind, ArrayDeque<Integer> buf, int r) {
        if (buf.isEmpty()) {
            return 0;
        }
        int[] c = listOf(buf);
        double[] d = oiDeltaMu(pasParam(sn, ind), c);
        int n = 0;
        for (int p = 0; p < c.length; p++) {
            if (d[p] > 0.0 && c[p] == r + 1) n++;
        }
        return n;
    }

    /**
     * Apply the pass-and-swap rewrite for a class-r (0-based) departure at the
     * PAS station ind. Which position completed is redrawn here among those whose
     * pass-and-swap ejects class r, weighted by that position's own Delta_mu,
     * which is the same split afterEventStationPas enumerates.
     */
    private static void oiDepart(NetworkStruct sn, int ind, ArrayDeque<Integer> buf, int r) {
        QueueNodeParam qp = pasParam(sn, ind);
        int[] c = listOf(buf);
        double[] d = oiDeltaMu(qp, c);
        List<Integer> pos = new ArrayList<Integer>();
        List<Double> w = new ArrayList<Double>();
        double tot = 0.0;
        for (int p = 0; p < c.length; p++) {
            if (d[p] <= 0.0) continue;
            if (AfterEventStation.passAndSwap(c, p, qp.swapGraph).depClass == r + 1) {
                pos.add(p);
                w.add(d[p]);
                tot += d[p];
            }
        }
        if (pos.isEmpty()) {
            return; // this class cannot depart from the current list
        }
        double u = Maths.rand() * tot;
        double acc = 0.0;
        int pick = pos.get(pos.size() - 1);
        for (int x = 0; x < pos.size(); x++) {
            acc += w.get(x);
            if (u < acc) { pick = pos.get(x); break; }
        }
        int[] cnew = AfterEventStation.passAndSwap(c, pick, qp.swapGraph).cnew;
        buf.clear();
        for (int i = 0; i < cnew.length; i++) {
            buf.addLast(cnew[i]);
        }
    }

    /**
     * Retrial orbits. An orbiting class-r job retries entry at the memoryless
     * rate sn.retrialMu and succeeds only when a server is free; otherwise the
     * event is a no-op and is simply not generated (State.afterEventStation,
     * EventType.RETRY). The orbit needs no new state: orbiting jobs are already
     * counted in the station population and held in the buffer, so orbit_r is
     * exactly the buffer occupancy the FCFS-family rate law already reads.
     */
    static boolean isRetrialStation(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0 || sn.retrialProc == null) return false;
        Map<JobClass, MatrixCell> pmap =
                sn.retrialProc.get(sn.stations.get((int) sn.nodeToStation.get(ind)));
        if (pmap == null) return false;
        for (MatrixCell pc : pmap.values()) {
            if (pc != null) return true;
        }
        return false;
    }

    /**
     * Retrial rate of a single orbiting class-r job at station ist
     * (sn.retrialMu), or 0 when the station-class has no orbit. retrialType=0
     * means EXP in MATLAB but "none" in Python, so the presence test goes through
     * retrialProc rather than through the type.
     */
    static double retrialRateOf(NetworkStruct sn, int ist, int r) {
        if (sn.retrialMu == null) {
            return 0.0;
        }
        Map<JobClass, Matrix> mumap = sn.retrialMu.get(sn.stations.get(ist));
        if (mumap == null) {
            return 0.0;
        }
        Matrix mu = mumap.get(sn.jobclasses.get(r));
        if (mu == null || mu.length() == 0) {
            return 0.0;
        }
        return mu.get(0);
    }

    /** Occupancy of class r (0-based) in a buffer. */
    private static int bufCount(ArrayDeque<Integer> buf, int r) {
        int n = 0;
        for (Integer cl : buf) {
            if (cl == r + 1) n++;
        }
        return n;
    }

    private static boolean isBuffered(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0) return false;
        int ist = (int) sn.nodeToStation.get(ind);
        return isBufferedSched(sn.sched.get(sn.stations.get(ist)));
    }

    /**
     * Index of the waiting job that the discipline at station ist promotes into
     * service. The buffer is ordered newest-first / oldest-last, matching the
     * convention of State.afterEventStation's space_buf.
     */
    private static int pickFromBuffer(ArrayDeque<Integer> buf, NetworkStruct sn, int ist) {
        SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
        int n = buf.size();
        if (sched == SchedStrategy.FCFS) {
            return n - 1;                       // oldest
        }
        if (sched == SchedStrategy.LCFS || sched == SchedStrategy.LCFSPR) {
            return 0;                           // newest / most recently preempted
        }
        if (sched == SchedStrategy.SIRO) {
            // Uniform over the waiting jobs: State.afterEventStation promotes a
            // class-r job with probability (nir(r)-sir(r))/(ni-sum(sir)), i.e.
            // the waiting class-r fraction, which is exactly a uniform draw.
            return (int) Math.floor(Maths.rand() * n);
        }
        Integer[] arr = buf.toArray(new Integer[0]);
        if (sched == SchedStrategy.HOL) {
            // Highest priority (lowest classprio value); FCFS within the group,
            // so the oldest = the last matching position.
            double best = Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                double p = sn.classprio.get(arr[i] - 1);
                if (p < best) best = p;
            }
            for (int i = n - 1; i >= 0; i--) {
                if (sn.classprio.get(arr[i] - 1) == best) return i;
            }
        }
        if (sched == SchedStrategy.SEPT || sched == SchedStrategy.LEPT) {
            // Shortest (SEPT) or longest (LEPT) expected processing time, i.e.
            // the highest or lowest service rate. Oldest first within a class.
            double best = (sched == SchedStrategy.SEPT) ? -Double.MAX_VALUE : Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                double rate = sn.rates.get(ist, arr[i] - 1);
                if (Double.isNaN(rate)) continue;
                if (sched == SchedStrategy.SEPT ? rate > best : rate < best) best = rate;
            }
            for (int i = n - 1; i >= 0; i--) {
                double rate = sn.rates.get(ist, arr[i] - 1);
                if (!Double.isNaN(rate) && rate == best) return i;
            }
        }
        throw new RuntimeException("pickFromBuffer: unsupported buffered policy " + sched);
    }

    /**
     * Remove the element at index idx (0 = newest). ArrayDeque has no positional
     * removal, so rotate the head off, drop the target, and restore order.
     */
    private static void removeAt(ArrayDeque<Integer> dq, int idx) {
        ArrayDeque<Integer> head = new ArrayDeque<Integer>();
        for (int i = 0; i < idx; i++) {
            head.addLast(dq.pollFirst());
        }
        dq.pollFirst();
        while (!head.isEmpty()) {
            dq.addFirst(head.pollLast());
        }
    }

    private static void printProgress(SolverOptions options, int samples_collected) {
        if (System.console() != null && !"parallel".equals(options.method)
                && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
            if (samples_collected == 2) {
                System.out.printf("\nSSA samples: %9d ", samples_collected);
                System.out.flush();
            } else if (samples_collected % 1000 == 0) {
                System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected);
                System.out.flush();
            }
            if (samples_collected == options.samples) {
                System.out.println();
            }
        }
    }
}
