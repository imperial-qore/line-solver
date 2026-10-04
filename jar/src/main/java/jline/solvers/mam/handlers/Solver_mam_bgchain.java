/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.api.mam.Mmap_super_safe;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnRtStations;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Mixed-network solver that treats the CLOSED classes as a background
 * modulating chain and the OPEN classes as matrix-analytic queues driven by it.
 *
 * <p>The closed population vector of a mixed network is a finite continuous-time
 * Markov chain in its own right: it is the only part of the model whose state
 * space is bounded. This method solves it exactly (given the capacity the open
 * work leaves free) and hands the open classes a station-local Markovian
 * environment read off that chain, so each open station becomes a level-dependent
 * QBD whose phase carries the number of closed jobs competing for its server.
 * The two halves meet at a fixed point on the closed capacity share.</p>
 *
 * <ol>
 * <li>background chain: the closed population vector over the stations the
 *     closed classes visit ({@link Mam_bgchain_ctmc})</li>
 * <li>environment: that chain lumped onto the closed occupancy of one station
 *     ({@link Mam_bgchain_env})</li>
 * <li>open station: a MAP/PH/c queue modulated by that environment, solved as a
 *     level-dependent QBD ({@link Mam_bgchain_station})</li>
 * <li>fixed point: the capacity share feeds step 1 and closes</li>
 * </ol>
 *
 * <p><b>Tagged-class iteration.</b> Step 1 is a population process of dimension
 * (closed chains) x (stations), so its state space is exponential in the number
 * of closed chains. With R &gt; 1 closed chains the method keeps ONE chain free
 * at a time: the tagged chain r is carried exactly, the other R-1 are replaced
 * by flow-equivalent aggregate classes whose population is their total and whose
 * service time and routing at each station are their throughput-weighted means
 * (Chandy-Herzog-Woo aggregation). Every chain takes its turn as the tagged one
 * and reads its own metrics off the chain it is exact in; the open-class results
 * are averaged over the passes.</p>
 *
 * <p><b>How much to aggregate</b> is {@code config.bgaggr}, the number G of
 * aggregate classes; the background chain then carries 1 + G. G = 1 is the
 * classic tagged/aggregate pair and the default, so the chain stays two-class
 * whatever R is; G &gt;= R-1 aggregates nothing, carries every closed chain
 * exactly, and answers in ONE pass instead of solving the same chain R times.
 * Passing R reaches that, so asking for no aggregation needs no magic value. The
 * cost is the state space, the product over the 1 + G classes of
 * nchoosek(N_b + Mc - 1, Mc - 1), capped by {@code bgstates_max}.</p>
 *
 * <p><b>Which chains share a group</b> is decided by similarity of per-station
 * SERVICE DEMAND, in {@link Mam_bgchain_groups}. An aggregate carries the
 * flow-weighted mean of its members' service times and routing, so it is exact
 * when they place the same demand at every station and distorts in proportion to
 * how far apart they are; grouping the demand-similar chains together keeps the
 * aggregation where it is harmless and away from the chains it would
 * misrepresent.</p>
 *
 * <p><b>Exactness</b>, as measured against SolverCTMC and exact MVA on mixed
 * models of two to four stations: PS or INF with ANY service law (exponential,
 * Erlang, HyperExp, Coxian), any number of servers, Poisson or MAP arrivals and
 * one to four closed chains agree to 4-5 significant digits, as does FCFS with
 * class-independent rates. FCFS with class-DEPENDENT rates keeps the closed
 * queue lengths within ~1% while the open queue length reads 14-20% low, because
 * the server is held here in random order rather than head-of-line.</p>
 *
 * <p>PS is INSENSITIVE to the service law beyond its mean, and the method
 * honours that rather than approximating it: at a PS station the open service is
 * replaced by the exponential of the same mean before the QBD is built. This QBD
 * tracks ONE service phase for the whole station, so carrying the phase-type
 * there makes the open queue length inherit the SCV-sensitivity of an M/PH/1
 * FCFS queue -- measured, a HyperExp of SCV 4 read 21% high where the exact
 * answer is the exponential one to five digits. At an FCFS station the service
 * law IS carried, and the background chain reads only the MEAN closed service
 * time, exact under PS by the same insensitivity and a first-moment surrogate
 * under FCFS.</p>
 */
public final class Solver_mam_bgchain {
    private Solver_mam_bgchain() {}

    /**
     * Number of states of the background-chain CTMC this solver would build on
     * the model, WITHOUT building it. Mirrors {@code mam_bgchain_states.m}.
     *
     * <p>The size is what decides whether bgchain is affordable and
     * {@link Mam_bgchain_ctmc} only discovers it after the partition is fixed,
     * so the default-method chooser needs it up front. The count follows the
     * partition used below: a pass carries the tagged closed chain as
     * background class 0 and the demand-similar groups of the other closed
     * chains as classes 1..G, each enumerating the compositions of its
     * population over the stations its members visit. Merging two chains onto
     * the UNION of their supports can raise the count as easily as lower it, so
     * the passes are enumerated rather than bounded and the largest returned:
     * that is the one Mam_bgchain_ctmc would refuse.</p>
     *
     * @return the state count, 0 when bgchain does not apply to the model at
     *         all, or {@link Double#POSITIVE_INFINITY} when it overflows
     */
    public static double bgchainStates(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int C = sn.nchains;

        int bgaggrOpt = 1;
        if (options != null && options.config != null) {
            Object cfgAggr = options.config.get("bgaggr");
            if (cfgAggr instanceof Number) bgaggrOpt = ((Number) cfgAggr).intValue();
        }

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);

        List<Integer> closedChains = new ArrayList<Integer>();
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            boolean open = false;
            for (int idx = 0; idx < inchain.length(); idx++) {
                int k = (int) inchain.get(idx);
                if (Double.isInfinite(sn.njobs.get(0, k))) open = true;
            }
            if (!open && dem.Nchain.get(0, c) > 0) {
                closedChains.add(Integer.valueOf(c));
            }
        }
        int R = closedChains.size();
        if (R == 0) {
            return 0.0;
        }

        List<Integer> cstList = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            for (int ci = 0; ci < R; ci++) {
                if (dem.Vchain.get(i, closedChains.get(ci).intValue()) > GlobalConstants.Zero) {
                    cstList.add(Integer.valueOf(i));
                    break;
                }
            }
        }
        int Mc = cstList.size();
        if (Mc == 0) {
            return 0.0;
        }
        int[] cst = new int[Mc];
        for (int ii = 0; ii < Mc; ii++) cst[ii] = cstList.get(ii).intValue();

        int naggr = Math.min(Math.max(bgaggrOpt, 1), Math.max(R - 1, 1));
        boolean noAggr = (R == 1) || (naggr >= R - 1);
        int npass = noAggr ? 1 : R;

        double worst = 0.0;
        for (int pidx = 0; pidx < npass; pidx++) {
            List<List<Integer>> members = new ArrayList<List<Integer>>();
            if (noAggr) {
                for (int oi = 0; oi < R; oi++) {
                    List<Integer> one = new ArrayList<Integer>();
                    one.add(closedChains.get(oi));
                    members.add(one);
                }
            } else {
                List<Integer> others = new ArrayList<Integer>();
                for (int oi = 0; oi < R; oi++) {
                    if (oi != pidx) others.add(closedChains.get(oi));
                }
                double[][] D = new double[Mc][others.size()];
                for (int ii = 0; ii < Mc; ii++) {
                    for (int oi = 0; oi < others.size(); oi++) {
                        D[ii][oi] = dem.Dchain.get(cst[ii], others.get(oi).intValue());
                    }
                }
                int[] grp = Mam_bgchain_groups.mam_bgchain_groups(D, naggr);
                List<Integer> tag = new ArrayList<Integer>();
                tag.add(closedChains.get(pidx));
                members.add(tag);
                for (int g = 0; g < naggr; g++) {
                    List<Integer> grpMembers = new ArrayList<Integer>();
                    for (int oi = 0; oi < others.size(); oi++) {
                        if (grp[oi] == g) grpMembers.add(others.get(oi));
                    }
                    members.add(grpMembers);
                }
            }

            double n = 1.0;
            for (int b = 0; b < members.size(); b++) {
                List<Integer> mem = members.get(b);
                double Nb = 0;
                for (int oi = 0; oi < mem.size(); oi++) Nb += dem.Nchain.get(0, mem.get(oi).intValue());
                int m = 0;
                for (int ii = 0; ii < Mc; ii++) {
                    for (int oi = 0; oi < mem.size(); oi++) {
                        if (dem.Vchain.get(cst[ii], mem.get(oi).intValue()) > GlobalConstants.Zero) {
                            m++;
                            break;
                        }
                    }
                }
                if (m == 0) m = 1;   // an empty class still needs one slot to be indexed by
                n *= binomial((int) Math.round(Nb) + m - 1, m - 1);
                if (!Double.isFinite(n)) {
                    return Double.POSITIVE_INFINITY;
                }
            }
            worst = Math.max(worst, n);
        }
        return worst;
    }

    /** nchoosek in floating point, so a chain far above any usable size still compares. */
    private static double binomial(int n, int k) {
        if (k < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 1.0;
        int kk = Math.min(k, n - k);
        double acc = 1.0;
        for (int i = 1; i <= kk; i++) {
            acc = acc * (n - kk + i) / i;
        }
        return acc;
    }

    private static final String MFILENAME = "solver_mam_bgchain";

    public static MAMResult solver_mam_bgchain(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        double tol = options.tol;

        int spaceMax = 128;
        Object cfgSpace = options.config.get("space_max");
        if (cfgSpace instanceof Number) spaceMax = ((Number) cfgSpace).intValue();
        int qbdPhasesMax = 500;
        Object cfgPhases = options.config.get("qbdphases_max");
        if (cfgPhases instanceof Number) qbdPhasesMax = ((Number) cfgPhases).intValue();
        int bgaggrOpt = 1;
        Object cfgAggr = options.config.get("bgaggr");
        if (cfgAggr instanceof Number) bgaggrOpt = ((Number) cfgAggr).intValue();

        Matrix S = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                double r = sn.rates.get(i, k);
                S.set(i, k, (Double.isFinite(r) && r > 0) ? 1.0 / r : 0.0);
            }
        }

        Pair<Matrix, Matrix> rtv = SnRtStations.snRtStations(sn);
        Matrix rtst = rtv.getLeft();
        Matrix V = rtv.getRight();
        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);

        boolean[] isopenchain = new boolean[C];
        List<Integer> openChains = new ArrayList<Integer>();
        List<Integer> closedChains = new ArrayList<Integer>();
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            boolean open = false;
            for (int idx = 0; idx < inchain.length(); idx++) {
                int k = (int) inchain.get(idx);
                if (Double.isInfinite(sn.njobs.get(0, k))) open = true;
            }
            isopenchain[c] = open;
            if (open) {
                openChains.add(Integer.valueOf(c));
            } else if (dem.Nchain.get(0, c) > 0) {
                closedChains.add(Integer.valueOf(c));
            }
        }
        int R = closedChains.size();
        if (R == 0) {
            InputOutput.line_error(MFILENAME, "The bgchain method requires at least one closed class: the "
                    + "background chain IS the closed population vector, so a purely open model has nothing "
                    + "to build it from. Use dec.source.");
        }

        // Stations the closed chains visit: the support of the background chain
        List<Integer> cstList = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            boolean visited = false;
            for (int ci = 0; ci < R; ci++) {
                if (dem.Vchain.get(i, closedChains.get(ci).intValue()) > GlobalConstants.Zero) visited = true;
            }
            if (visited) cstList.add(Integer.valueOf(i));
        }
        int Mc = cstList.size();
        if (Mc == 0) {
            InputOutput.line_error(MFILENAME, "The closed classes of this model visit no station.");
        }
        int[] cst = new int[Mc];
        for (int ii = 0; ii < Mc; ii++) cst[ii] = cstList.get(ii).intValue();

        // Chain-level station routing, folding the class axis of sn.rt
        List<Matrix> Pchain = new ArrayList<Matrix>();
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            Matrix P = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int ki = 0; ki < inchain.length(); ki++) {
                    int k = (int) inchain.get(ki);
                    double a = dem.alpha.get(i, k);
                    if (a <= 0) continue;
                    for (int j = 0; j < M; j++) {
                        double acc = 0;
                        for (int kpi = 0; kpi < inchain.length(); kpi++) {
                            int kp = (int) inchain.get(kpi);
                            acc += rtst.get(i * K + k, j * K + kp);
                        }
                        if (acc != 0) P.set(i, j, P.get(i, j) + a * acc);
                    }
                }
            }
            Pchain.add(P);
        }

        // Open arrival streams
        double[] lambdaChain = new double[C];
        Map<Integer, MatrixCell> chainArrival = new HashMap<Integer, MatrixCell>();
        for (int ci = 0; ci < openChains.size(); ci++) {
            int c = openChains.get(ci).intValue();
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            int isrc = (int) sn.refstat.get((int) inchain.get(0), 0);
            double lam = 0;
            MatrixCell acc = null;
            Station srcStation = sn.stations.get(isrc);
            for (int ki = 0; ki < inchain.length(); ki++) {
                int k = (int) inchain.get(ki);
                double rk = sn.rates.get(isrc, k);
                if (!Double.isFinite(rk) || rk <= 0) continue;
                lam += rk;
                JobClass jc = sn.jobclasses.get(k);
                MatrixCell proc = sn.proc.get(srcStation).get(jc);
                if (proc == null || proc.get(0) == null || proc.get(0).hasNaN()) continue;
                MatrixCell mk = new MatrixCell(3);
                mk.set(0, proc.get(0));
                mk.set(1, proc.get(1));
                mk.set(2, proc.get(1));
                if (acc == null) {
                    acc = mk;
                } else {
                    Map<Integer, MatrixCell> pair = new HashMap<Integer, MatrixCell>();
                    pair.put(Integer.valueOf(0), acc);
                    pair.put(Integer.valueOf(1), mk);
                    acc = Mmap_super_safe.mmap_super_safe(pair, spaceMax, "default");
                }
            }
            lambdaChain[c] = lam;
            if (acc != null) {
                MatrixCell map = new MatrixCell(2);
                map.set(0, acc.get(0));
                map.set(1, acc.get(1));
                chainArrival.put(Integer.valueOf(c), map);
            }
        }

        Matrix lambdaOpen = new Matrix(M, K);
        boolean[] isopenclass = new boolean[K];
        for (int ci = 0; ci < openChains.size(); ci++) {
            int c = openChains.get(ci).intValue();
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            for (int ki = 0; ki < inchain.length(); ki++) {
                int k = (int) inchain.get(ki);
                isopenclass[k] = true;
                for (int i = 0; i < M; i++) {
                    lambdaOpen.set(i, k, lambdaChain[c] * V.get(i, k));
                }
            }
        }

        MAMResult result = new MAMResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);

        int Ntot = 0;
        for (int ci = 0; ci < R; ci++) Ntot += (int) dem.Nchain.get(0, closedChains.get(ci).intValue());

        // cshare[i][e]: mean number of servers of station i that its e closed jobs
        // hold once the open work has taken its share. Starts at min(e,c).
        double[][] cshare = new double[M][Ntot + 1];
        for (int i = 0; i < M; i++) {
            for (int e = 0; e <= Ntot; e++) {
                cshare[i][e] = Math.min(e, sn.nservers.get(i, 0));
            }
        }
        double[] Xclosed = new double[C];
        for (int ci = 0; ci < R; ci++) {
            int c = closedChains.get(ci).intValue();
            double denom = 0;
            for (int i = 0; i < M; i++) denom += dem.Dchain.get(i, c);
            if (denom > 0) Xclosed[c] = dem.Nchain.get(0, c) / denom;
        }

        // How many aggregate classes the background chain carries, and which chains
        // share each of them. config.bgaggr is the number of AGGREGATE classes G:
        // G = 1 is the classic tagged/aggregate pair, G >= R-1 aggregates nothing.
        int naggr = Math.min(Math.max(bgaggrOpt, 1), Math.max(R - 1, 1));
        // With nothing left to aggregate ONE background chain carries every closed
        // chain exactly, so the tagged loop would solve the same chain R times over.
        boolean noAggr = (R == 1) || (naggr >= R - 1);

        // The grouping is a property of the demands, not of the iterate, so it is
        // fixed once here rather than recomputed inside the fixed point.
        List<List<Integer>> othersOf = new ArrayList<List<Integer>>();
        List<int[]> grpOf = new ArrayList<int[]>();
        for (int ridx = 0; ridx < R; ridx++) {
            List<Integer> others = new ArrayList<Integer>();
            for (int oi = 0; oi < R; oi++) {
                if (oi != ridx) others.add(closedChains.get(oi));
            }
            othersOf.add(others);
            if (noAggr || others.isEmpty()) {
                grpOf.add(new int[0]);
            } else {
                double[][] D = new double[Mc][others.size()];
                for (int ii = 0; ii < Mc; ii++) {
                    for (int oi = 0; oi < others.size(); oi++) {
                        D[ii][oi] = dem.Dchain.get(cst[ii], others.get(oi).intValue());
                    }
                }
                grpOf.add(Mam_bgchain_groups.mam_bgchain_groups(D, naggr));
            }
        }
        final int npass = noAggr ? 1 : R;

        Matrix TNprev = new Matrix(M, K);
        TNprev.fill(Double.POSITIVE_INFINITY);
        int totiter = 0;
        double relax = 0.5;

        while (maxAbsDiff(result.TN, TNprev) > tol && totiter < options.iter_max) {
            totiter++;
            TNprev = result.TN.copy();

            double[] QopenAcc = new double[M];
            double[] UopenAcc = new double[M];
            double[][] cshareAcc = new double[M][Ntot + 1];

            for (int pidx = 0; pidx < npass; pidx++) {
                int r = closedChains.get(pidx).intValue();

                // Background classes: class 0 is the tagged chain, classes 1..G the
                // flow-equivalent aggregates of the demand-similar groups. With no
                // aggregation every closed chain is a class of its own, in chain order.
                List<List<Integer>> members = new ArrayList<List<Integer>>();
                if (noAggr) {
                    for (int oi = 0; oi < R; oi++) {
                        List<Integer> one = new ArrayList<Integer>();
                        one.add(closedChains.get(oi));
                        members.add(one);
                    }
                } else {
                    List<Integer> others = othersOf.get(pidx);
                    int[] grp = grpOf.get(pidx);
                    List<Integer> tag = new ArrayList<Integer>();
                    tag.add(Integer.valueOf(r));
                    members.add(tag);
                    for (int g = 0; g < naggr; g++) {
                        List<Integer> mem = new ArrayList<Integer>();
                        for (int oi = 0; oi < others.size(); oi++) {
                            if (grp[oi] == g) mem.add(others.get(oi));
                        }
                        members.add(mem);
                    }
                }
                int B = members.size();
                int[] Nb = new int[B];
                double[][] STb = new double[Mc][B];
                List<Matrix> Pb = new ArrayList<Matrix>();
                // A class can only hold jobs at the stations its members visit; see
                // Mam_bgchain_ctmc on why the union instead makes the chain reducible.
                boolean[][] suppb = new boolean[Mc][B];
                for (int b = 0; b < B; b++) {
                    List<Integer> mem = members.get(b);
                    int nmem = mem.size();
                    for (int oi = 0; oi < nmem; oi++) {
                        int o = mem.get(oi).intValue();
                        Nb[b] += (int) dem.Nchain.get(0, o);
                        for (int ii = 0; ii < Mc; ii++) {
                            if (dem.Vchain.get(cst[ii], o) > GlobalConstants.Zero) suppb[ii][b] = true;
                        }
                    }
                    if (nmem == 1) {
                        // a group of one is carried exactly: no mean to take
                        int o = mem.get(0).intValue();
                        for (int ii = 0; ii < Mc; ii++) STb[ii][b] = dem.STchain.get(cst[ii], o);
                        Pb.add(subMatrix(Pchain.get(o), cst));
                    } else if (nmem == 0) {
                        Pb.add(new Matrix(Mc, Mc));
                    } else {
                        double[][] w = new double[Mc][nmem];
                        for (int ii = 0; ii < Mc; ii++) {
                            double rowsum = 0;
                            for (int oi = 0; oi < nmem; oi++) {
                                int o = mem.get(oi).intValue();
                                w[ii][oi] = Xclosed[o] * dem.Vchain.get(cst[ii], o);
                                rowsum += w[ii][oi];
                            }
                            for (int oi = 0; oi < nmem; oi++) {
                                w[ii][oi] = (rowsum > 0) ? w[ii][oi] / rowsum : 1.0 / nmem;
                            }
                        }
                        Matrix Pagg = new Matrix(Mc, Mc);
                        for (int ii = 0; ii < Mc; ii++) {
                            double st = 0;
                            for (int oi = 0; oi < nmem; oi++) {
                                int o = mem.get(oi).intValue();
                                st += w[ii][oi] * dem.STchain.get(cst[ii], o);
                                Matrix Po = Pchain.get(o);
                                for (int jj = 0; jj < Mc; jj++) {
                                    double v = w[ii][oi] * Po.get(cst[ii], cst[jj]);
                                    if (v != 0) Pagg.set(ii, jj, Pagg.get(ii, jj) + v);
                                }
                            }
                            STb[ii][b] = st;
                        }
                        Pb.add(Pagg);
                    }
                }
                for (int b = 0; b < B; b++) rowNormalize(Pb.get(b));

                SchedStrategy[] schedC = new SchedStrategy[Mc];
                double[] nserversC = new double[Mc];
                double[][] cshareC = new double[Mc][];
                for (int ii = 0; ii < Mc; ii++) {
                    schedC[ii] = sn.sched.get(sn.stations.get(cst[ii]));
                    nserversC[ii] = sn.nservers.get(cst[ii], 0);
                    cshareC[ii] = cshare[cst[ii]];
                }

                Mam_bgchain_ctmc.Result bg = Mam_bgchain_ctmc.mam_bgchain_ctmc(
                        Nb, STb, Pb, schedC, nserversC, cshareC, suppb, options);

                // Closed-class metrics of every chain this pass carries EXACTLY: the
                // tagged one always, and every chain when nothing was aggregated.
                int bmax = noAggr ? B : 1;
                for (int b = 0; b < bmax; b++) {
                    if (members.get(b).isEmpty()) continue;
                    int rb = members.get(b).get(0).intValue();
                    Matrix inchain = sn.inchain.get(Integer.valueOf(rb));
                    for (int ki = 0; ki < inchain.length(); ki++) {
                        int k = (int) inchain.get(ki);
                        for (int i = 0; i < M; i++) {
                            result.QN.set(i, k, 0);
                            result.UN.set(i, k, 0);
                            result.RN.set(i, k, 0);
                            result.TN.set(i, k, 0);
                        }
                    }
                    for (int ii = 0; ii < Mc; ii++) {
                        int i = cst[ii];
                        for (int ki = 0; ki < inchain.length(); ki++) {
                            int k = (int) inchain.get(ki);
                            double a = dem.alpha.get(i, k);
                            if (a <= 0) continue;
                            // THROUGHPUT splits by VISIT share, OCCUPANCY by
                            // DEMAND share. A chain queue divided by alpha alone
                            // gives every class of a station the same response
                            // time, impossible at a Delay where R must be the
                            // class service time; the weight is alpha*ST/STchain,
                            // the rule snDeaggregateChainResults applies.
                            double stc = dem.STchain.get(i, rb);
                            double w = (stc > GlobalConstants.Zero) ? a * S.get(i, k) / stc : a;
                            double q = bg.QLen[ii][b] * w;
                            double x = bg.Tput[ii][b] * a;
                            result.QN.set(i, k, q);
                            result.TN.set(i, k, x);
                            if (schedC[ii] == SchedStrategy.INF) {
                                result.UN.set(i, k, q);
                            } else {
                                result.UN.set(i, k, bg.Ubusy[ii][b] * w);
                            }
                            result.RN.set(i, k, (x > GlobalConstants.Zero) ? q / x : 0.0);
                        }
                    }
                    int iref = (int) sn.refstat.get((int) inchain.get(0), 0);
                    double tputref = 0;
                    for (int ki = 0; ki < inchain.length(); ki++) {
                        tputref += result.TN.get(iref, (int) inchain.get(ki));
                    }
                    double vref = dem.Vchain.get(iref, rb);
                    Xclosed[rb] = (vref > GlobalConstants.Zero) ? tputref / vref : tputref;
                }

                double[] Uclosed = new double[M];
                for (int ii = 0; ii < Mc; ii++) {
                    double u = 0;
                    for (int b = 0; b < B; b++) u += bg.Ubusy[ii][b];
                    Uclosed[cst[ii]] = u;
                }

                openPass(sn, options, bg, cst, S, V, lambdaOpen, lambdaChain, chainArrival,
                        openChains, Uclosed, cshare, qbdPhasesMax, spaceMax,
                        QopenAcc, UopenAcc, cshareAcc);
            }

            double[] Qopen = new double[M];
            for (int i = 0; i < M; i++) {
                Qopen[i] = QopenAcc[i] / npass;
                for (int e = 0; e <= Ntot; e++) {
                    cshare[i][e] = (1 - relax) * cshare[i][e] + relax * (cshareAcc[i][e] / npass);
                }
            }

            // Open-class metrics from the aggregate station results
            for (int i = 0; i < M; i++) {
                List<Integer> kopen = new ArrayList<Integer>();
                for (int k = 0; k < K; k++) {
                    if (isopenclass[k] && lambdaOpen.get(i, k) > GlobalConstants.Zero) {
                        kopen.add(Integer.valueOf(k));
                    }
                }
                if (kopen.isEmpty()) {
                    for (int k = 0; k < K; k++) {
                        if (!isopenclass[k]) continue;
                        result.TN.set(i, k, lambdaOpen.get(i, k));
                        result.QN.set(i, k, 0);
                        result.UN.set(i, k, 0);
                        result.RN.set(i, k, 0);
                    }
                    continue;
                }
                double lamtot = 0;
                double work = 0;
                for (int kk = 0; kk < kopen.size(); kk++) {
                    int k = kopen.get(kk).intValue();
                    lamtot += lambdaOpen.get(i, k);
                    work += lambdaOpen.get(i, k) * S.get(i, k);
                }
                double Smix = work / lamtot;
                SchedStrategy sched = sn.sched.get(sn.stations.get(i));
                for (int kk = 0; kk < kopen.size(); kk++) {
                    int k = kopen.get(kk).intValue();
                    result.TN.set(i, k, lambdaOpen.get(i, k));
                    if (sched == SchedStrategy.EXT) {
                        result.QN.set(i, k, 0);
                        result.UN.set(i, k, 0);
                        result.RN.set(i, k, 0);
                    } else if (sched == SchedStrategy.INF) {
                        result.RN.set(i, k, S.get(i, k));
                        result.QN.set(i, k, lambdaOpen.get(i, k) * S.get(i, k));
                        result.UN.set(i, k, lambdaOpen.get(i, k) * S.get(i, k));
                    } else {
                        double Rtot = Qopen[i] / lamtot;
                        double rk;
                        if (sched == SchedStrategy.PS) {
                            // processor sharing: residence scales with the demand
                            rk = Rtot * S.get(i, k) / Smix;
                        } else {
                            // FCFS and its variants: the wait is class-blind, the
                            // service time is not
                            rk = Math.max(S.get(i, k), Rtot - Smix + S.get(i, k));
                        }
                        result.RN.set(i, k, rk);
                        result.QN.set(i, k, lambdaOpen.get(i, k) * rk);
                        // Utilization Law: a c-server station holds TN*S/c
                        result.UN.set(i, k, lambdaOpen.get(i, k) * S.get(i, k) / sn.nservers.get(i, 0));
                    }
                }
            }
        }

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            for (int ki = 0; ki < inchain.length(); ki++) {
                int k = (int) inchain.get(ki);
                result.XN.set(0, k, isopenchain[c] ? lambdaChain[c] : Xclosed[c]);
            }
        }
        for (int k = 0; k < K; k++) {
            double acc = 0;
            for (int i = 0; i < M; i++) acc += result.RN.get(i, k);
            result.CN.set(0, k, acc);
        }
        sanitize(result.QN);
        sanitize(result.UN);
        sanitize(result.RN);
        sanitize(result.TN);
        sanitize(result.CN);
        sanitize(result.XN);

        result.iter = totiter;
        result.method = "bgchain";
        return result;
    }

    /**
     * One pass of the open side: for each station, lump the background chain onto
     * its closed occupancy and solve the resulting modulated level-dependent QBD.
     * Accumulates into QopenAcc/UopenAcc/cshareAcc.
     */
    private static void openPass(NetworkStruct sn, SolverOptions options, Mam_bgchain_ctmc.Result bg,
                                 int[] cst, Matrix S, Matrix V, Matrix lambdaOpen, double[] lambdaChain,
                                 Map<Integer, MatrixCell> chainArrival, List<Integer> openChains,
                                 double[] Uclosed, double[][] cshare, int qbdPhasesMax, int spaceMax,
                                 double[] QopenAcc, double[] UopenAcc, double[][] cshareAcc) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int Ngrid = cshare[0].length;

        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            boolean solved = false;
            if (sched != SchedStrategy.EXT && sched != SchedStrategy.INF) {
                List<Integer> kopen = new ArrayList<Integer>();
                for (int k = 0; k < K; k++) {
                    if (lambdaOpen.get(i, k) > GlobalConstants.Zero) kopen.add(Integer.valueOf(k));
                }
                if (!kopen.isEmpty()) {
                    MatrixCell Da = null;
                    for (int ci = 0; ci < openChains.size(); ci++) {
                        int c = openChains.get(ci).intValue();
                        if (lambdaChain[c] <= GlobalConstants.Zero) continue;
                        Matrix inchain = sn.inchain.get(Integer.valueOf(c));
                        double vsum = 0;
                        for (int ki = 0; ki < inchain.length(); ki++) {
                            vsum += V.get(i, (int) inchain.get(ki));
                        }
                        double rate_ic = lambdaChain[c] * vsum;
                        if (rate_ic <= GlobalConstants.Zero) continue;
                        MatrixCell base = chainArrival.get(Integer.valueOf(c));
                        if (base == null) continue;
                        MatrixCell scaled = Map_scale.map_scale(base, 1.0 / rate_ic);
                        if (Da == null) {
                            Da = scaled;
                        } else {
                            Map<Integer, MatrixCell> pair = new HashMap<Integer, MatrixCell>();
                            pair.put(Integer.valueOf(0), asMarked(Da));
                            pair.put(Integer.valueOf(1), asMarked(scaled));
                            MatrixCell sup = Mmap_super_safe.mmap_super_safe(pair, spaceMax, "default");
                            MatrixCell merged = new MatrixCell(2);
                            merged.set(0, sup.get(0));
                            merged.set(1, sup.get(1));
                            Da = merged;
                        }
                    }
                    if (Da != null) {
                        // Arrival-weighted phase-type mixture of the open service laws
                        double lamtot = 0;
                        for (int kk = 0; kk < kopen.size(); kk++) {
                            lamtot += lambdaOpen.get(i, kopen.get(kk).intValue());
                        }
                        List<Matrix> pies = new ArrayList<Matrix>();
                        List<Matrix> subgens = new ArrayList<Matrix>();
                        int msTotal = 0;
                        // PROCESSOR SHARING IS INSENSITIVE to the service law beyond its mean, so
                        // carrying the phase-type representation at a PS station is not merely unnecessary,
                        // it is WRONG. This QBD tracks ONE service phase for the whole station, which makes
                        // the open queue length inherit the SCV-sensitivity of an M/PH/1 FCFS queue;
                        // measured against SolverCTMC, a HyperExp of SCV 4 then read 21% high where the
                        // exact answer is the exponential one to five digits. The exponential of the same
                        // mean is exact here, and it shrinks the QBD's phase count as a side effect.
                        final boolean isPSstation = (sched == SchedStrategy.PS);
                        for (int kk = 0; kk < kopen.size(); kk++) {
                            int k = kopen.get(kk).intValue();
                            MatrixCell phk = sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(k));
                            MatrixCell scaledPh;
                            if (isPSstation || phk == null || phk.get(0) == null
                                    || phk.get(0).hasNaN()) {
                                Matrix d0 = new Matrix(1, 1);
                                Matrix d1 = new Matrix(1, 1);
                                double rate = (S.get(i, k) > 0) ? 1.0 / S.get(i, k) : GlobalConstants.Immediate;
                                d0.set(0, 0, -rate);
                                d1.set(0, 0, rate);
                                scaledPh = new MatrixCell(2);
                                scaledPh.set(0, d0);
                                scaledPh.set(1, d1);
                            } else {
                                scaledPh = Map_scale.map_scale(phk, S.get(i, k));
                            }
                            Matrix pik = Map_pie.map_pie(scaledPh);
                            pies.add(pik.scale(lambdaOpen.get(i, k) / lamtot));
                            subgens.add(scaledPh.get(0));
                            msTotal += scaledPh.get(0).getNumRows();
                        }
                        Matrix alphaS = new Matrix(1, msTotal);
                        Matrix Tblk = new Matrix(msTotal, msTotal);
                        int off = 0;
                        for (int kk = 0; kk < pies.size(); kk++) {
                            Matrix pik = pies.get(kk);
                            Matrix sub = subgens.get(kk);
                            int n = sub.getNumRows();
                            for (int a = 0; a < n; a++) {
                                alphaS.set(0, off + a, pik.get(0, a));
                                for (int b = 0; b < n; b++) {
                                    double v = sub.get(a, b);
                                    if (v != 0) Tblk.set(off + a, off + b, v);
                                }
                            }
                            off += n;
                        }

                        Matrix A;
                        int[] esup;
                        int ipos = indexOf(cst, i);
                        if (ipos >= 0) {
                            Mam_bgchain_env.Env env = Mam_bgchain_env.mam_bgchain_env(bg, ipos);
                            A = env.A;
                            esup = env.esup;
                        } else {
                            A = new Matrix(1, 1);
                            esup = new int[] { 0 };
                        }

                        int nphases = Da.get(0).getNumRows() * esup.length * msTotal;
                        if (nphases > qbdPhasesMax) {
                            InputOutput.line_error(MFILENAME, "The modulated QBD of station " + i + " needs "
                                    + nphases + " phases (" + Da.get(0).getNumRows() + " arrival x " + esup.length
                                    + " environment x " + msTotal + " service), above the limit of " + qbdPhasesMax
                                    + ". The environment axis is the closed population held by the station, so it "
                                    + "grows with the closed population. Raise options.config.qbdphases_max, or "
                                    + "reduce the closed population or the order of the arrival and service "
                                    + "processes.");
                        }

                        double[] gref = new double[esup.length];
                        for (int e = 0; e < esup.length; e++) {
                            gref[e] = cshare[i][Math.min(esup[e], Ngrid - 1)];
                        }
                        double svcWork = 0;
                        for (int kk = 0; kk < kopen.size(); kk++) {
                            int k = kopen.get(kk).intValue();
                            svcWork += lambdaOpen.get(i, k) * S.get(i, k);
                        }
                        int Kmax = cutoff(options, lamtot, svcWork / lamtot, sn.nservers.get(i, 0), Uclosed[i]);

                        Mam_bgchain_station.StationResult res = Mam_bgchain_station.mam_bgchain_station(
                                Da.get(0), Da.get(1), alphaS, Tblk, A, esup, sn.nservers.get(i, 0),
                                gref, Kmax, options);

                        QopenAcc[i] += res.QLen;
                        UopenAcc[i] += res.Util;
                        // The QBD only saw the environment states the background chain
                        // reaches; interpolate the rest so the next background chain has
                        // a share wherever it may go, clipped to what a closed job can hold.
                        for (int e = 0; e < Ngrid; e++) {
                            double g = interpolate(res.esup, res.cshare, e);
                            cshareAcc[i][e] += Math.min(Math.max(g, 0.0), Math.min(e, sn.nservers.get(i, 0)));
                        }
                        solved = true;
                    }
                }
            }
            if (!solved) {
                // a station with no open queue keeps the share it had
                for (int e = 0; e < Ngrid; e++) cshareAcc[i][e] += cshare[i][e];
            }
        }
    }

    private static MatrixCell asMarked(MatrixCell map) {
        MatrixCell mk = new MatrixCell(3);
        mk.set(0, map.get(0));
        mk.set(1, map.get(1));
        mk.set(2, map.get(1));
        return mk;
    }

    /** Linear interpolation with flat-slope extrapolation, over an ascending support. */
    private static double interpolate(int[] xs, double[] ys, double x) {
        int n = xs.length;
        if (n == 1) return ys[0];
        if (x <= xs[0]) {
            double slope = (ys[1] - ys[0]) / (xs[1] - xs[0]);
            return ys[0] + slope * (x - xs[0]);
        }
        if (x >= xs[n - 1]) {
            double slope = (ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]);
            return ys[n - 1] + slope * (x - xs[n - 1]);
        }
        for (int idx = 0; idx < n - 1; idx++) {
            if (x >= xs[idx] && x <= xs[idx + 1]) {
                double w = (x - xs[idx]) / (xs[idx + 1] - xs[idx]);
                return ys[idx] + w * (ys[idx + 1] - ys[idx]);
            }
        }
        return ys[n - 1];
    }

    /**
     * Truncation level of the open queue: the explicit cutoff when the user set
     * one, else enough levels for the geometric tail left by the closed traffic
     * to be negligible.
     */
    private static int cutoff(SolverOptions options, double lambda, double Smix, double nservers, double Uclosed) {
        if (options.cutoff != null && options.cutoff.length() > 0) {
            double cutoffVal = options.cutoff.get(0, 0);
            if (Double.isFinite(cutoffVal) && cutoffVal > 0) {
                return Math.max(2, (int) Math.round(cutoffVal));
            }
        }
        double free = Math.max(GlobalConstants.FineTol, 1.0 - Uclosed);
        double rho = lambda * Smix / (nservers * free);
        rho = Math.min(Math.max(rho, 1e-3), 1 - 1e-3);
        int Kmax = (int) Math.ceil(Math.log(1e-8) / Math.log(rho));
        return Math.min(Math.max(Kmax, 20), 200);
    }

    private static int indexOf(int[] arr, int v) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == v) return i;
        }
        return -1;
    }

    private static Matrix subMatrix(Matrix P, int[] idx) {
        Matrix out = new Matrix(idx.length, idx.length);
        for (int i = 0; i < idx.length; i++) {
            for (int j = 0; j < idx.length; j++) {
                double v = P.get(idx[i], idx[j]);
                if (v != 0) out.set(i, j, v);
            }
        }
        return out;
    }

    /**
     * Row-normalizes a routing matrix, leaving an all-zero row as a self-loop so
     * the background chain stays a proper Markov chain on its support.
     */
    private static void rowNormalize(Matrix P) {
        int n = P.getNumRows();
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) s += P.get(i, j);
            if (s > GlobalConstants.Zero) {
                for (int j = 0; j < n; j++) {
                    double v = P.get(i, j);
                    if (v != 0) P.set(i, j, v / s);
                }
            } else {
                for (int j = 0; j < n; j++) P.set(i, j, 0);
                P.set(i, i, 1.0);
            }
        }
    }

    private static double maxAbsDiff(Matrix a, Matrix b) {
        double m = 0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                double d = Math.abs(a.get(i, j) - b.get(i, j));
                if (d > m) m = d;
            }
        }
        return m;
    }

    private static void sanitize(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                double v = m.get(i, j);
                if (!Double.isFinite(v)) m.set(i, j, 0);
            }
        }
    }
}
