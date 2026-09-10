package jline.solvers.ag.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.SignalType;
import jline.lang.processes.DiscreteDistribution;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Port of {@code build_rcat} in solver_ag.m: the network to its RCAT
 * components and actions.
 *
 * Each component is a QBD over (queue length, phase), the phase being the pair
 * (arrival phase, service phase). Every scalar rate of the exponential
 * construction becomes a matrix block and every unit passive entry an identity
 * block, so with exponential processes the matrices are entry for entry the
 * ones this builder produced before.
 */
public final class Solver_ag_build {
    private Solver_ag_build() {}

    /** Row/column offset of level N (0-based) in a component with MPH phases. */
    static int blk(int n, int mph) {
        return n * mph;
    }

    /** Add BLOCK, scaled by SCALE, into the (LI, LJ) level block of TARGET. */
    static void addBlock(Matrix target, int li, int lj, int mph, Matrix block, double scale) {
        int r0 = blk(li, mph);
        int c0 = blk(lj, mph);
        for (int i = 0; i < mph; i++) {
            for (int j = 0; j < mph; j++) {
                double v = block.get(i, j);
                if (v != 0.0) {
                    target.set(r0 + i, c0 + j, target.get(r0 + i, c0 + j) + scale * v);
                }
            }
        }
    }

    /** Add SCALE times the identity into the (LI, LJ) level block of TARGET. */
    static void addIdentity(Matrix target, int li, int lj, int mph, double scale) {
        int r0 = blk(li, mph);
        int c0 = blk(lj, mph);
        for (int i = 0; i < mph; i++) {
            target.set(r0 + i, c0 + i, target.get(r0 + i, c0 + i) + scale);
        }
    }

    /** Set the (LI, LJ) level block of TARGET to the identity. */
    static void setIdentity(Matrix target, int li, int lj, int mph) {
        int r0 = blk(li, mph);
        int c0 = blk(lj, mph);
        for (int i = 0; i < mph; i++) {
            target.set(r0 + i, c0 + i, 1.0);
        }
    }

    /**
     * True when (D0, D1) is a genuine MAP rather than a RAP or an ME process:
     * non-negative off-diagonal rates in D0, non-negative rates in D1, and
     * (D0 + D1) an infinitesimal generator. A CTMC assembled from anything else
     * is a rational generator whose stationary solution is a signed vector.
     */
    static boolean isMarkovianMap(Matrix D0, Matrix D1) {
        if (D0 == null || D1 == null) return false;
        int ns = D0.getNumRows();
        if (D0.getNumCols() != ns || D1.getNumRows() != ns || D1.getNumCols() != ns) return false;
        double scale = 1.0;
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                if (Double.isNaN(D0.get(i, j)) || Double.isInfinite(D0.get(i, j))) return false;
                if (Double.isNaN(D1.get(i, j)) || Double.isInfinite(D1.get(i, j))) return false;
                scale = Math.max(scale, Math.abs(D0.get(i, j)));
                scale = Math.max(scale, Math.abs(D1.get(i, j)));
            }
        }
        double tol = 1e-9 * scale;
        for (int i = 0; i < ns; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < ns; j++) {
                if (i != j && D0.get(i, j) < -tol) return false;
                if (D1.get(i, j) < -tol) return false;
                rowSum += D0.get(i, j) + D1.get(i, j);
            }
            if (Math.abs(rowSum) > tol) return false;
        }
        return true;
    }

    /**
     * (D0,D1) of the process at (IST,R), or the exponential pair built from
     * sn.rates when the struct carries no usable matrix representation.
     *
     * A non-Markovian pair is refused here and answered as its mean rate; the
     * solver-level gate rejects those models before they reach this point.
     */
    static Matrix[] procMap(NetworkStruct sn, int ist, int r) {
        Matrix D0 = null;
        Matrix D1 = null;
        if (sn.proc != null && ist < sn.stations.size() && r < sn.jobclasses.size()) {
            Station station = sn.stations.get(ist);
            JobClass jobClass = sn.jobclasses.get(r);
            java.util.Map<JobClass, MatrixCell> pm = sn.proc.get(station);
            if (pm != null) {
                MatrixCell cell = pm.get(jobClass);
                if (cell != null && cell.size() >= 2) {
                    Matrix c0 = cell.get(0);
                    Matrix c1 = cell.get(1);
                    if (isMarkovianMap(c0, c1)) {
                        D0 = c0;
                        D1 = c1;
                    }
                }
            }
        }
        if (D0 == null) {
            double rate = sn.rates.get(ist, r);
            if (Double.isNaN(rate) || rate <= 0) rate = 0.0;
            D0 = new Matrix(1, 1);
            D0.set(0, 0, -rate);
            D1 = new Matrix(1, 1);
            D1.set(0, 0, rate);
        }
        return new Matrix[]{D0, D1};
    }

    /**
     * Accumulate the level n to level m block of a batch removal into B, scaled
     * by RATE. Landing on the empty level absorbs the whole upper tail of the
     * pmf, which keeps the block stochastic once the batch exceeds the queue
     * length. The phase is untouched: a removal takes a waiting job, not the one
     * in service.
     */
    static void addBatchRemoval(Matrix B, DiscreteDistribution dist, double rate, int nlev, int mph) {
        for (int n = 1; n < nlev; n++) {
            for (int m = 1; m <= n; m++) {
                double pk = dist.evalPMF((double) (n - m));
                if (pk > 0) addIdentity(B, n, m, mph, rate * pk);
            }
            double cdf = 0.0;
            for (int j = 0; j < n; j++) cdf += dist.evalPMF((double) j);
            double tail = 1.0 - cdf;
            if (tail > 0) addIdentity(B, n, 0, mph, rate * tail);
        }
    }

    public static RCATModel solver_ag_build(NetworkStruct sn, int maxStates) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix rt = sn.rt;

        List<Integer> sourceStations = new ArrayList<Integer>();
        List<Integer> queueStations = new ArrayList<Integer>();
        List<Integer> sinkNodes = new ArrayList<Integer>();
        for (int nodeIdx = 0; nodeIdx < sn.nodetype.size(); nodeIdx++) {
            if (sn.nodetype.get(nodeIdx) == NodeType.Sink) sinkNodes.add(nodeIdx);
        }
        for (int ist = 0; ist < M; ist++) {
            int nodeIdx = (int) sn.stationToNode.get(ist);
            NodeType nt = sn.nodetype.get(nodeIdx);
            if (nt == NodeType.Source) sourceStations.add(ist);
            else queueStations.add(ist);
        }

        int processIdx = 0;
        Matrix processMap = new Matrix(M, K);
        processMap.fill(-1.0);
        for (int ist : queueStations) {
            for (int r = 0; r < K; r++) {
                if (sn.issignal != null && sn.issignal.get(r, 0) > 0) continue;
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && rate > 0) {
                    processMap.set(ist, r, processIdx);
                    processIdx++;
                }
            }
        }
        int numProcesses = processIdx;
        if (numProcesses == 0) {
            return new RCATModel(new Matrix[1][1], new Matrix(0, 2), processMap,
                    new ArrayList<ActionInfo>(), new int[0], new int[0], new int[0],
                    new int[0][], new double[0][], new double[0][]);
        }

        // The QBD shape of every component: both the service MAP of its station
        // and the arrival MAP of the external streams reaching it.
        RcatComponent[] comp = new RcatComponent[numProcesses];
        int[] N = new int[numProcesses];
        int[] nlev = new int[numProcesses];
        int[] mph = new int[numProcesses];
        int[][] level = new int[numProcesses][];
        double[][] svcrate = new double[numProcesses][];
        double[][] svcdown = new double[numProcesses][];
        for (int p = 0; p < numProcesses; p++) {
            RcatComponent c = new RcatComponent();
            outerLoop:
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    if ((int) processMap.get(ist, r) == p) {
                        c.ist = ist;
                        c.r = r;
                        break outerLoop;
                    }
                }
            }
            Matrix[] svc = procMap(sn, c.ist, c.r);
            c.Ds0 = svc[0];
            c.Ds1 = svc[1];
            arrivalMap(sn, c, rt, sourceStations, K);
            c.ns = c.Ds0.getNumRows();
            c.na = c.Da0.getNumRows();
            c.mph = c.na * c.ns;
            double njobsR = sn.njobs.get(c.r);
            c.nlev = Double.isInfinite(njobsR) ? maxStates : ((int) njobsR + 1);
            c.Dsvc = Matrix.eye(c.na).kron(c.Ds1);
            c.N = c.nlev * c.mph;
            comp[p] = c;

            N[p] = c.N;
            nlev[p] = c.nlev;
            mph[p] = c.mph;
            level[p] = new int[c.N];
            for (int n = 0; n < c.nlev; n++) {
                for (int j = 0; j < c.mph; j++) level[p][n * c.mph + j] = n;
            }
            double[] row = new double[c.mph];
            for (int i = 0; i < c.mph; i++) {
                double s = 0.0;
                for (int j = 0; j < c.mph; j++) s += c.Dsvc.get(i, j);
                row[i] = s;
            }
            svcdown[p] = row;
            svcrate[p] = new double[c.N];
            for (int n = 1; n < c.nlev; n++) {
                for (int j = 0; j < c.mph; j++) svcrate[p][n * c.mph + j] = row[j];
            }
        }

        List<ActionInfo> actionMap = new ArrayList<ActionInfo>();
        for (int ist : queueStations) {
            for (int r = 0; r < K; r++) {
                if ((int) processMap.get(ist, r) >= 0) {
                    // Removal signal: NEGATIVE or CATASTROPHE (the two are
                    // distinct SignalType values, so both must be tested).
                    boolean isNegativeClass = sn.issignal != null
                            && sn.issignal.get(r, 0) > 0
                            && sn.signaltype != null
                            && (sn.signaltype.get(r) == SignalType.NEGATIVE
                                || sn.signaltype.get(r) == SignalType.CATASTROPHE);
                    boolean isCatastropheClass = (sn.iscatastrophe != null
                            && sn.iscatastrophe.get(r, 0) > 0)
                            || (sn.signaltype != null && sn.signaltype.get(r) == SignalType.CATASTROPHE);
                    DiscreteDistribution removalDist = null;
                    if (sn.signalremdist != null && r < sn.signalremdist.size()) {
                        removalDist = sn.signalremdist.get(r);
                    }
                    for (int jst : queueStations) {
                        for (int s = 0; s < K; s++) {
                            if ((int) processMap.get(jst, s) >= 0) {
                                double probIJRS = rt.get(ist * K + r, jst * K + s);
                                if (probIJRS > 0 && (ist != jst || r != s)) {
                                    actionMap.add(new ActionInfo(ist, r, jst, s, probIJRS, isNegativeClass, isCatastropheClass, removalDist));
                                }
                            }
                        }
                    }
                }
            }
        }
        int numActions = actionMap.size();

        Matrix[][] R = new Matrix[numActions + 1][Math.max(numProcesses, 2)];
        Matrix AP = new Matrix(Math.max(numActions, 1), 2);

        for (int p = 0; p < numProcesses; p++) {
            R[numActions][p] = buildLocalRates(sn, comp[p], rt, sinkNodes, K);
        }

        for (int a = 0; a < numActions; a++) {
            ActionInfo am = actionMap.get(a);
            int pActive = (int) processMap.get(am.fromStation, am.fromClass);
            AP.set(a, 0, (double) pActive);
            RcatComponent pa = comp[pActive];
            double prob = am.prob;

            // Active matrix: level n -> n-1 carrying the service completion
            // block kron(I, D1^s), scaled by this action's routing probability.
            Matrix Aa = new Matrix(pa.N, pa.N);
            for (int n = 1; n < pa.nlev; n++) addBlock(Aa, n, n - 1, pa.mph, pa.Dsvc, prob);
            // see _kb/06-solver-catalog.md for rationale. The boundary self-loop
            // is written on the DIAGONAL so it stays inert in the generator while
            // still contributing the pi(i)/pi(i) = 1 ratio the INAP estimators
            // read off the top level.
            if (!Double.isInfinite(sn.njobs.get(am.fromClass))) {
                int off = blk(pa.nlev - 1, pa.mph);
                for (int i = 0; i < pa.mph; i++) {
                    Aa.set(off + i, off + i, Aa.get(off + i, off + i) + svcdown[pActive][i] * prob);
                }
            }
            R[a][0] = Aa;

            int pPassive = (int) processMap.get(am.toStation, am.toClass);
            AP.set(a, 1, (double) pPassive);
            RcatComponent pp = comp[pPassive];

            Matrix Pb = new Matrix(pp.N, pp.N);
            if (am.isNegative) {
                if (am.isCatastrophe) {
                    // CATASTROPHE: every level drops to level 0
                    for (int n = 0; n < pp.nlev; n++) setIdentity(Pb, n, 0, pp.mph);
                } else if (am.removalDistribution != null) {
                    addBatchRemoval(Pb, am.removalDistribution, 1.0, pp.nlev, pp.mph);
                    setIdentity(Pb, 0, 0, pp.mph);  // an empty queue absorbs the signal
                } else {
                    setIdentity(Pb, 0, 0, pp.mph);
                    for (int n = 1; n < pp.nlev; n++) setIdentity(Pb, n, n - 1, pp.mph);
                }
            } else {
                // POSITIVE: a normal arrival. The phase is untouched: a job
                // joining does not restart the server, and the service phase
                // frozen at level 0 is the one the last completion left behind,
                // which for a phase-type is already its entry distribution.
                for (int n = 0; n < pp.nlev - 1; n++) setIdentity(Pb, n, n + 1, pp.mph);
                setIdentity(Pb, pp.nlev - 1, pp.nlev - 1, pp.mph);
            }
            R[a][1] = Pb;
        }
        return new RCATModel(R, AP, processMap, actionMap, N, nlev, mph, level, svcrate, svcdown);
    }

    /**
     * External (Source) streams reaching the component, as one arrival MAP for
     * the positive customers plus the scalar rates of the removal signals.
     *
     * Each stream is thinned by its routing probability -- a MAP thinned with
     * probability p is (D0 + (1-p) D1, p D1) -- and the streams are superposed
     * by the Kronecker sum, so several Poisson sources still collapse to the
     * single rate sum this analyzer used before. Removal signals stay scalar: a
     * signal is a trigger with no service, and its arrival process is required
     * exponential.
     */
    private static void arrivalMap(NetworkStruct sn, RcatComponent c, Matrix rt,
                                   List<Integer> sourceStations, int K) {
        boolean haveArrival = false;
        Matrix Da0 = new Matrix(1, 1);
        Matrix Da1 = new Matrix(1, 1);

        for (int isrc : sourceStations) {
            for (int sSrc = 0; sSrc < K; sSrc++) {
                boolean isSignalSrc = sn.issignal != null && sn.issignal.get(sSrc, 0) > 0;
                double probSrc;
                if (isSignalSrc) {
                    double pSum = 0.0;
                    for (int sDst = 0; sDst < K; sDst++) pSum += rt.get(isrc * K + sSrc, c.ist * K + sDst);
                    probSrc = pSum;
                } else {
                    probSrc = rt.get(isrc * K + sSrc, c.ist * K + c.r);
                }
                double srcRate = sn.rates.get(isrc, sSrc);
                if (!(probSrc > 0) || Double.isNaN(srcRate)) continue;

                // see _kb/06-solver-catalog.md for rationale
                boolean isRemovalSrc = isSignalSrc && sn.signaltype != null
                        && (sn.signaltype.get(sSrc) == SignalType.NEGATIVE
                            || sn.signaltype.get(sSrc) == SignalType.CATASTROPHE);
                if (isRemovalSrc) {
                    boolean isCat = (sn.iscatastrophe != null && sn.iscatastrophe.get(sSrc, 0) > 0)
                            || sn.signaltype.get(sSrc) == SignalType.CATASTROPHE;
                    if (isCat) {
                        c.lamCat += srcRate * probSrc;
                        continue;
                    }
                    DiscreteDistribution removalDist = null;
                    if (sn.signalremdist != null && sSrc < sn.signalremdist.size()) {
                        removalDist = sn.signalremdist.get(sSrc);
                    }
                    if (removalDist != null) {
                        c.batchRates.add(srcRate * probSrc);
                        c.batchDists.add(removalDist);
                    } else {
                        c.lamNeg += srcRate * probSrc;
                    }
                    continue;
                }

                if (!(srcRate > 0)) continue;
                Matrix[] src = procMap(sn, isrc, sSrc);
                Matrix S0 = src[0];
                Matrix S1 = src[1];
                if (probSrc < 1) {
                    S0 = S0.add(1.0 - probSrc, S1);
                    S1 = S1.scale(probSrc);
                }
                if (haveArrival) {
                    Da0 = Da0.krons(S0);
                    Da1 = Da1.krons(S1);
                } else {
                    Da0 = S0;
                    Da1 = S1;
                    haveArrival = true;
                }
            }
        }
        c.Da0 = Da0;
        c.Da1 = Da1;
    }

    private static Matrix buildLocalRates(NetworkStruct sn, RcatComponent c, Matrix rt,
                                          List<Integer> sinkNodes, int K) {
        int mph = c.mph;
        int nlev = c.nlev;
        Matrix L = new Matrix(c.N, c.N);

        // Level-local blocks: the arrival phase always runs, the service phase
        // only while the server is busy (qbd_mapmap1's Lbar = kron(D0^a, I) at
        // level 0 and L = krons(D0^a, D0^s) above it). With one phase each these
        // are pure diagonals, which ctmc_makeinfgen discards and rebuilds from
        // the row sums.
        addBlock(L, 0, 0, mph, c.Da0.kron(Matrix.eye(c.ns)), 1.0);
        Matrix Lbusy = c.Da0.krons(c.Ds0);
        for (int n = 1; n < nlev; n++) addBlock(L, n, n, mph, Lbusy, 1.0);

        // Positive arrivals: level n -> n+1, carrying kron(D1^a, I).
        Matrix Aup = c.Da1.kron(Matrix.eye(c.ns));
        for (int n = 0; n < nlev - 1; n++) addBlock(L, n, n + 1, mph, Aup, 1.0);
        // At the truncation the job is lost but the arrival process still moves
        // on, so the block stays on the top level. With a single arrival phase
        // this is a pure diagonal and is discarded, exactly as before.
        addBlock(L, nlev - 1, nlev - 1, mph, Aup, 1.0);

        // Catastrophe arrivals: every busy level drops to level 0
        if (c.lamCat > 0) {
            for (int n = 1; n < nlev; n++) addIdentity(L, n, 0, mph, c.lamCat);
        }
        for (int b = 0; b < c.batchRates.size(); b++) {
            addBatchRemoval(L, c.batchDists.get(b), c.batchRates.get(b), nlev, mph);
        }
        // Single-removal negative arrivals: level n -> n-1 (busy levels only)
        if (c.lamNeg > 0) {
            for (int n = 1; n < nlev; n++) addIdentity(L, n, n - 1, mph, c.lamNeg);
        }

        // Service completions that are not synchronizing actions: departures to
        // a Sink (level down) and self-routing (level unchanged, service
        // restarted).
        double muIr = sn.rates.get(c.ist, c.r);
        if (!Double.isNaN(muIr) && muIr > 0) {
            int nodeIdx = (int) sn.stationToNode.get(c.ist);
            double probSink = 0.0;
            if (sn.rtnodes != null && sn.rtnodes.getNumRows() > 0) {
                for (int jsnk : sinkNodes) {
                    for (int s = 0; s < K; s++) {
                        int fromIdx = nodeIdx * K + c.r;
                        int toIdx = jsnk * K + s;
                        if (fromIdx < sn.rtnodes.getNumRows() && toIdx < sn.rtnodes.getNumCols()) {
                            probSink += sn.rtnodes.get(fromIdx, toIdx);
                        }
                    }
                }
            }
            double probSelf = rt.get(c.ist * K + c.r, c.ist * K + c.r);
            if (probSink > 0) {
                for (int n = 1; n < nlev; n++) addBlock(L, n, n - 1, mph, c.Dsvc, probSink);
            }
            if (probSelf > 0) {
                for (int n = 1; n < nlev; n++) addBlock(L, n, n, mph, c.Dsvc, probSelf);
            }
        }
        return L;
    }
}
