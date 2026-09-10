package jline.solvers.mam.analyzers;

import java.util.Map;

import jline.api.qsys.Qsys_is_retrial;
import jline.api.qsys.RetrialInfo;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnHasForkJoin;
import jline.api.sn.SnHasLoadDependence;
import jline.api.sn.SnHasMixedClasses;
import jline.api.sn.SnIsClosedModel;
import jline.api.sn.SnIsDiscreteTime;
import jline.api.sn.SnIsOpenModel;
import jline.api.sn.SnNonmarkovToPh;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.solvers.mam.handlers.Solver_mam;
import jline.solvers.mam.handlers.Solver_mam_basic;
import jline.solvers.mam.handlers.Solver_mam_basic_mmap;
import jline.solvers.mam.handlers.Solver_mam_bgchain;
import jline.solvers.mam.handlers.Solver_mam_ldqbd;
import jline.solvers.mam.handlers.Solver_mam_retrial;
import jline.solvers.mam.handlers.Solver_mna_closed;
import jline.solvers.mam.handlers.Solver_mna_open;

public final class Solver_mam_analyzer {
    private Solver_mam_analyzer() {}

    /**
     * Check if the model is a single-class closed Delay+Queue, the exact regime
     * of Solver_mam_ldqbd: one class, finite population, exactly two stations,
     * one INF (Delay) and one FCFS (Queue). Mirrors the isClosedDelayQueue
     * subfunction of the MATLAB solver_mam_analyzer.
     */
    /**
     * Whether any station declares a setup/delay-off server. Mirrors {@code hasSetup}
     * in the MATLAB solver_mam_analyzer.
     *
     * <p>Setup/delay-off is read by {@code Solver_mam_basic} and, since 2026-09, by
     * {@code Solver_mam_ldqbd}. The bgchain and mna analyzers still DROP it, so the
     * default dispatch must not hand them such a model. The ldqbd branch does not
     * consult this: {@code isClosedDelayQueue} admits exactly the regime the exact
     * analysis covers and refuses the rest.
     *
     * @param sn the network structure
     * @return true when some station carries setup/delay-off
     */
    private static boolean hasSetup(NetworkStruct sn) {
        if (sn.hassetup == null) {
            return false;
        }
        for (int i = 0; i < sn.hassetup.getNumRows(); i++) {
            if (sn.hassetup.get(i, 0) == 1.0) {
                return true;
            }
        }
        return false;
    }

    private static boolean isClosedDelayQueue(NetworkStruct sn) {
        if (sn.nclasses != 1 || sn.nstations != 2) {
            return false;
        }
        for (int r = 0; r < sn.njobs.getNumCols(); r++) {
            if (!Double.isFinite(sn.njobs.get(0, r))) {
                return false;
            }
        }
        int nDelay = 0;
        int nQueue = 0;
        int queueIdx = -1;
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.INF) {
                nDelay++;
            } else if (sched == SchedStrategy.FCFS) {
                nQueue++;
                queueIdx = i;
            }
        }
        if (nDelay != 1 || nQueue != 1) {
            return false;
        }
        // A SETUP/DELAY-OFF QUEUE is exact here only at a single server with
        // exponential service and no load dependence, which is what
        // Qbd_setupdelayoff_closed models. Outside that, LDQBD refuses by name,
        // so the default has to fall through to the decomposition rather than
        // reach the refusal; dec.source carries the same analysis approximately.
        if (sn.hassetup != null && sn.hassetup.getNumRows() > queueIdx
                && sn.hassetup.get(queueIdx, 0) == 1.0) {
            if (sn.nservers.get(queueIdx) > 1 || SnHasLoadDependence.snHasLoadDependence(sn)) {
                return false;
            }
            MatrixCell PHq = sn.proc.get(sn.stations.get(queueIdx)).get(sn.jobclasses.get(0));
            Matrix D0q = PHq.get(0);
            if (D0q.getNumRows() != 1 || D0q.getNumCols() != 1) {
                return false;   // exponential service only
            }
        }
        return true;
    }

    /**
     * Whether the closed MNA analyzer covers this model. Mirrors the mnaApplies
     * subfunction of the MATLAB solver_mam_analyzer: the round-robin split is carried
     * by the open traffic equations only, a self-looping class has no inter-station
     * flow to decompose, and a station whose discipline the flow sweep does not update
     * would keep a zero queue length. Any of those routes the default to dec.source.
     */
    private static boolean mnaApplies(NetworkStruct sn) {
        // Solver_mna_closed drives its bisection over CLASSES but stores the throughput
        // in the CHAIN-indexed lambda, and renormalizes chain c's queue lengths with the
        // class-indexed njobs(c). Both are only correct when each chain holds exactly one
        // class, so a class-switching model (fewer chains than classes) is outside what
        // the analyzer computes and belongs to dec.source.
        if (sn.nchains != sn.nclasses) {
            return false;
        }
        if (sn.routing != null) {
            for (Map<JobClass, RoutingStrategy> perClass : sn.routing.values()) {
                if (perClass == null) {
                    continue;
                }
                for (RoutingStrategy rs : perClass.values()) {
                    if (rs == RoutingStrategy.RROBIN) {
                        return false;
                    }
                }
            }
        }

        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched != SchedStrategy.INF && sched != SchedStrategy.PS
                    && sched != SchedStrategy.FCFS && sched != SchedStrategy.EXT) {
                return false;
            }
            // The PS branch of Solver_mna_closed forms U = S*T and the geometric bound
            // from it WITHOUT dividing by the number of servers, so a multiserver PS
            // station is misrepresented; dec.source is exact on the non-queueing regime
            // (c >> N) that shape usually stands for.
            if (sched == SchedStrategy.PS && sn.nservers.get(i) > 1) {
                return false;
            }
            // A multiclass FCFS station makes the flow sweep superpose one MMAP per
            // class and then solve MMAPPH1FCFS at level sum(N)+1: measured against an
            // exact CTMC that costs 12-61s where dec.source costs 0.03s and is not more
            // accurate (mean relative error 0.11-0.38 against 0.04-0.19). The
            // single-class case is both cheap and better, so keep only that one.
            if (sched == SchedStrategy.FCFS && sn.nclasses > 1) {
                return false;
            }
        }

        Matrix V = new Matrix(sn.nstations, sn.nclasses);
        for (Matrix vc : sn.visits.values()) {
            if (vc == null) {
                continue;
            }
            for (int i = 0; i < sn.nstations; i++) {
                for (int k = 0; k < sn.nclasses; k++) {
                    V.set(i, k, V.get(i, k) + vc.get(i, k));
                }
            }
        }
        for (int k = 0; k < sn.nclasses; k++) {
            if (!Double.isFinite(sn.njobs.get(0, k))) {
                continue;
            }
            int seen = 0;
            int at = -1;
            for (int i = 0; i < sn.nstations; i++) {
                if (V.get(i, k) > GlobalConstants.FineTol) {
                    seen++;
                    at = i;
                }
            }
            if (seen == 1) {
                SchedStrategy sched = sn.sched.get(sn.stations.get(at));
                if (sched != SchedStrategy.INF && sched != SchedStrategy.EXT) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Whether the closed MNA outer bisection closed on N. Mirrors the mnaConserves
     * subfunction of the MATLAB solver_mam_analyzer: Solver_mna_closed rescales each
     * chain onto its population as a last step, so a diverged bisection still returns
     * queue lengths that sum to N and the failure is invisible in QN. R and T are NOT
     * rescaled, so Little's law over the whole network, sum_i T(i,k)*R(i,k) = N_k,
     * still reads the raw iterate: a converged run lands within 5e-4 of N and a
     * diverged one is orders of magnitude out, or negative.
     */
    private static boolean mnaConserves(NetworkStruct sn, MAMResult result) {
        if (result == null || result.QN == null || result.RN == null || result.TN == null) {
            return false;
        }
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (!Double.isFinite(result.QN.get(i, k)) || !Double.isFinite(result.RN.get(i, k))
                        || !Double.isFinite(result.TN.get(i, k))) {
                    return false;
                }
            }
        }
        for (int k = 0; k < sn.nclasses; k++) {
            double nk = sn.njobs.get(0, k);
            if (!Double.isFinite(nk) || nk <= 0) {
                continue;
            }
            double npred = 0;
            for (int i = 0; i < sn.nstations; i++) {
                npred += result.TN.get(i, k) * result.RN.get(i, k);
            }
            if (Math.abs(npred - nk) > 0.01 * nk) {
                return false;
            }
        }
        return true;
    }

    /**
     * Whether the background-chain method covers this model. Mirrors the
     * bgchainApplies subfunction of the MATLAB solver_mam_analyzer: no class
     * priorities, no fork-join, and at least one station visited by the closed
     * classes, which is what the background chain is built over.
     */
    private static boolean bgchainApplies(NetworkStruct sn, SolverOptions options) {
        boolean prioSched = false;
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.HOL || sched == SchedStrategy.FCFSPRPRIO) {
                prioSched = true;
                break;
            }
        }
        if (prioSched && sn.classprio != null && sn.classprio.length() > 1) {
            double first = sn.classprio.get(0);
            for (int k = 1; k < sn.classprio.length(); k++) {
                if (sn.classprio.get(k) != first) {
                    return false;
                }
            }
        }

        if (SnHasForkJoin.snHasForkJoin(sn)) {
            return false;
        }

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        boolean visited = false;
        for (int c = 0; c < sn.nchains; c++) {
            if (!Double.isFinite(dem.Nchain.get(c)) || dem.Nchain.get(c) <= 0) {
                continue;
            }
            for (int i = 0; i < sn.nstations; i++) {
                if (dem.Vchain.get(i, c) > GlobalConstants.Zero) {
                    visited = true;
                    break;
                }
            }
            if (visited) {
                break;
            }
        }
        if (!visited) {
            return false;
        }

        // The background chain enumerates the closed population vector over the
        // stations the closed classes visit, so a large closed population makes
        // Mam_bgchain_ctmc refuse the model outright. That refusal is right when
        // the user asked for bgchain by name and wrong as a default, which must
        // land on a method that answers: size the chain first and leave those
        // models to dec.source. The limit is the one Mam_bgchain_ctmc enforces.
        int bgstatesMax = 20000;
        if (options != null && options.config != null) {
            Object cfg = options.config.get("bgstates_max");
            if (cfg instanceof Number) bgstatesMax = ((Number) cfg).intValue();
        }
        return Solver_mam_bgchain.bgchainStates(sn, options) <= bgstatesMax;
    }

    /**
     * Whether the background chain represents this closed model's service laws
     * exactly. Mirrors the bgchainClosedExact subfunction of the MATLAB
     * solver_mam_analyzer.
     *
     * <p>A station that is not PS or INF must satisfy BOTH conditions below,
     * because the background chain makes two separate first-moment substitutions
     * there.</p>
     *
     * <ol>
     * <li>Mam_bgchain_ctmc builds its generator from the MEAN service time alone.
     * That is exact at a PS or INF station, which is insensitive to the service
     * law beyond its first moment, and exact under any discipline when the law IS
     * exponential. Measured on a closed Delay+FCFS cycle with Erlang-3 service,
     * the mean-only chain reads 2.8% off SolverCTMC.</li>
     * <li>The capacity a station's closed jobs hold is split over the background
     * classes in proportion to their COUNTS, which is service in random order.
     * That is exact under PS, and exact under FCFS only when the classes are
     * served at the SAME rate -- an FCFS station with class-dependent rates reads
     * 25.2% off SolverCTMC on a two-chain closed cycle, against 3.5e-16 when the
     * two rates are made equal.</li>
     * </ol>
     *
     * <p>Solver_mna_closed carries the phase-type representation instead, so
     * neither surrogate may be chosen as the DEFAULT. Asking for bgchain by name
     * still gets it, with both approximations documented in
     * Solver_mam_bgchain.</p>
     */
    private static boolean bgchainClosedExact(NetworkStruct sn) {
        for (int ist = 0; ist < sn.nstations; ist++) {
            Station station = sn.stations.get(ist);
            SchedStrategy sched = sn.sched.get(station);
            if (sched == SchedStrategy.INF || sched == SchedStrategy.PS || sched == SchedStrategy.EXT) {
                continue;
            }
            Map<JobClass, ProcessType> procMap = sn.procid == null ? null : sn.procid.get(station);
            double rateHere = Double.NaN;
            for (int r = 0; r < sn.nclasses; r++) {
                double rate = sn.rates.get(ist, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    continue;
                }
                ProcessType pt = procMap == null ? null : procMap.get(sn.jobclasses.get(r));
                if (pt != ProcessType.EXP) {
                    return false;
                }
                if (Double.isNaN(rateHere)) {
                    rateHere = rate;
                } else if (Math.abs(rate - rateHere) > GlobalConstants.CoarseTol * rateHere) {
                    return false;
                }
            }
        }
        return true;
    }

    public static MAMResult solver_mam_analyzer(NetworkStruct snInput, SolverOptions options) {
        long start = System.nanoTime();

        // A finite buffer a CLOSED class can fill: every MAM analyzer models a
        // buffer as a LOSS buffer, and a closed job that finds no room blocks
        // instead. Same predicate SolverMAM.supportsModelMethod reports, so the
        // report and the run cannot answer differently.
        String bufferWhy = jline.solvers.mam.SolverMAM.bufferRefusal(snInput);
        if (!bufferWhy.isEmpty()) {
            throw new RuntimeException(bufferWhy);
        }

        // Discrete-time (slotted) models are recognized from the distributions
        // and routed to the Q-MAM discrete-time algorithms. The test must run
        // BEFORE SnNonmarkovToPh, which would fit a continuous surrogate to a
        // Geometric and erase the lattice; see _kb/06-solver-catalog.md for the
        // LAS-DA convention.
        SnIsDiscreteTime.Result dt = SnIsDiscreteTime.snIsDiscreteTime(snInput, options);
        if (dt.isDiscreteTime) {
            MAMResult dtResult = Solver_mam_dt.solver_mam_dt(snInput, options, dt.slotLength);
            dtResult.runtime = (System.nanoTime() - start) / 1e9;
            return dtResult;
        }

        // see _kb/06-solver-catalog.md for rationale
        MAMResult exactMapMap1 = Solver_mam_mapmap1_exact.solver_mam_mapmap1_exact(snInput);
        if (exactMapMap1 != null) {
            exactMapMap1.runtime = (System.nanoTime() - start) / 1e9;
            return exactMapMap1;
        }

        // see _kb/06-solver-catalog.md for rationale
        MAMResult exactMmck = Solver_mam_mmck_exact.solver_mam_mmck_exact(snInput);
        if (exactMmck != null) {
            exactMmck.runtime = (System.nanoTime() - start) / 1e9;
            return exactMmck;
        }

        NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(snInput, options);

        // Check if the model is mixed (has both open and closed classes)
        boolean isOpen = SnIsOpenModel.snIsOpenModel(sn);
        boolean isClosed = SnIsClosedModel.snIsClosedModel(sn);
        // isOpen means EVERY class is open and isClosed that every class is closed, so a
        // mixed model is neither, not both
        boolean isMixed = SnHasMixedClasses.snHasMixedClasses(sn);

        // Mixed models are supported by the dec.source method

        options.config.merge = "super";
        options.config.compress = "mixture.order1";
        options.config.space_max = 128;
        Object etaqaTrunc = options.config.get("etaqa_trunc");
        int etaqaTruncVal = 0;
        if (etaqaTrunc instanceof Integer) {
            etaqaTruncVal = ((Integer) etaqaTrunc).intValue();
        }
        if (!options.config.containsKey("etaqa_trunc") || etaqaTruncVal == 0) {
            options.config.put("etaqa_trunc", Integer.valueOf(8));
        }
        MAMResult result = new MAMResult();

        if ("dec.mmap".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using dec.mmap method, calling solver_mam");
            result = Solver_mam.solver_mam(sn, options);
            result.method = "dec.mmap";
        } else if ("dec.source.mmap".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using dec.source.mmap method, calling solver_mam_basic_mmap");
            result = Solver_mam_basic_mmap.solver_mam_basic_mmap(sn, options);
            result.method = "dec.source.mmap";
        } else if ("default".equals(options.method) || "dec.source".equals(options.method)) {
            if (SnHasForkJoin.snHasForkJoin(sn) && SnIsOpenModel.snIsOpenModel(sn)) {
                InputOutput.line_debug(options.verbose, "Detected open fork-join topology, using dec.source.mmap");
                result = Solver_mam_basic_mmap.solver_mam_basic_mmap(sn, options);
                result.method = "dec.source.mmap";
            } else {

                // Check if network is a valid BMAP/PH/N/N bufferless retrial queue
                RetrialInfo retInfo;
                try {
                    retInfo = Qsys_is_retrial.qsys_is_retrial(sn);
                } catch (Exception e) {
                    retInfo = new RetrialInfo();
                }

                if (retInfo.isRetrial()) {
                    // Use BMAP/PH/N/N retrial solver
                    if ("default".equals(options.method)) {
                        InputOutput.line_debug(options.verbose, "Default method: using retrial for BMAP/PH/N/N bufferless topology");
                    }
                    InputOutput.line_debug(options.verbose, "Detected BMAP/PH/N/N retrial topology, using retrial method");
                    result = Solver_mam_retrial.solver_mam_retrial(sn, options);
                    result.method = "retrial";
                } else if ("default".equals(options.method) && isClosedDelayQueue(sn)) {
                    // see _kb/06-solver-catalog.md for rationale
                    InputOutput.line_debug(options.verbose, "Default method: using LDQBD for single-class closed Delay/Queue");
                    result = Solver_mam_ldqbd.solver_mam_ldqbd(sn, options);
                    result.method = "ldqbd";
                // bgchain DROPS setup/delay-off: it would answer with the always-warm
                    // chain and nothing would say so. Only Solver_mam_basic and
                    // Solver_mam_ldqbd read it, so a setup model must not reach here.
                } else if ("default".equals(options.method) && isClosed && !hasSetup(sn)
                        && bgchainApplies(sn, options)
                        && bgchainClosedExact(sn)) {
                    // A closed model is the degenerate case of the background chain: with no
                    // open work to take a share of the servers the chain is the EXACT closed
                    // CTMC at chain granularity, so it dominates the mna fixed point wherever
                    // its state space fits. bgchainApplies sizes that state space and
                    // bgchainClosedExact checks the chain is built from a service law it
                    // represents exactly.
                    InputOutput.line_debug(options.verbose, "Default method: using bgchain for a closed model");
                    result = jline.solvers.mam.handlers.Solver_mam_bgchain.solver_mam_bgchain(sn, options);
                    result.method = "bgchain";
                // mna DROPS setup/delay-off for the same reason as bgchain.
                } else if ("default".equals(options.method) && isClosed && !hasSetup(sn)
                        && mnaApplies(sn)) {
                    // A closed model has no arrival stream for dec.source to build its
                    // Poisson surrogate from: it replaces each closed chain by a source at
                    // the current throughput iterate and never enforces the population, so
                    // the answer neither conserves N nor separates the classes. MNA closes
                    // the same traffic equations by bisecting the per-class throughput
                    // against N; see _kb/06-solver-catalog.md.
                    InputOutput.line_debug(options.verbose, "Default method: using mna for a closed model");
                    result = Solver_mna_closed.solver_mna_closed(sn, options);
                    result.method = "mna";
                    if (!mnaConserves(sn, result)) {
                        // The outer bisection did not close on the population: the last
                        // step rescales each chain onto N regardless, so the failure is
                        // invisible in QN alone and only Little's law on the unrescaled
                        // R and T still shows it.
                        InputOutput.line_debug(options.verbose,
                                "mna did not close on the population, falling back to dec.source");
                        result = Solver_mam_basic.solver_mam_basic(sn, options);
                        result.method = "dec.source";
                    }
                } else if ("default".equals(options.method) && isMixed && !hasSetup(sn)
                        && bgchainApplies(sn, options)) {
                    // Mixed: the closed classes are solved exactly as a background chain and
                    // the open ones as QBDs driven by it, which is 4-5 significant digits
                    // against CTMC where dec.source is 10-24% out.
                    InputOutput.line_debug(options.verbose, "Default method: using bgchain for a mixed model");
                    result = jline.solvers.mam.handlers.Solver_mam_bgchain.solver_mam_bgchain(sn, options);
                    result.method = "bgchain";
                } else {
                    // arrival process per chain rescaled by visits at each node
                    if ("default".equals(options.method)) {
                        InputOutput.line_debug(options.verbose, "Default method: using dec.source");
                    }
                    InputOutput.line_debug(options.verbose, "Using default/dec.source method, calling solver_mam_basic");
                    result = Solver_mam_basic.solver_mam_basic(sn, options);
                    result.method = "dec.source";
                }
            }
        } else if ("dec.poisson".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using dec.poisson method with space_max=1, calling solver_mam_basic");
            options.config.space_max = 1;
            result = Solver_mam_basic.solver_mam_basic(sn, options);
            result.method = "dec.poisson";
        } else if ("mna".equals(options.method)) {
            if (SnIsClosedModel.snIsClosedModel(sn)) {
                InputOutput.line_debug(options.verbose, "MNA method (closed)");
                result = Solver_mna_closed.solver_mna_closed(sn, options);
                result.method = "mna";
            } else if (SnIsOpenModel.snIsOpenModel(sn)) {
                InputOutput.line_debug(options.verbose, "MNA method (open)");
                result = (MAMResult) Solver_mna_open.solver_mna_open(sn, options);
                result.method = "mna";
            } else {
                throw new RuntimeException("The mna method in SolverMAM does not support mixed models.");
            }
        } else if ("bgchain".equals(options.method)) {
            // Mixed networks: the closed classes are a background modulating
            // chain, the open classes are QBDs driven by it. With several closed
            // chains an outer iteration tags one chain at a time and aggregates
            // the rest, so the background chain always carries two classes.
            InputOutput.line_debug(options.verbose, "Using bgchain method, calling solver_mam_bgchain");
            result = jline.solvers.mam.handlers.Solver_mam_bgchain.solver_mam_bgchain(sn, options);
            result.method = "bgchain";
        } else if ("ldqbd".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using LDQBD method for single-class closed network");
            result = Solver_mam_ldqbd.solver_mam_ldqbd(sn, options);
            result.method = "ldqbd";
        } else if ("retrial".equals(options.method)) {
            // The predicate SolverMAM.supportsModelMethod asks, so the gate that
            // decides whether to OFFER 'retrial' and this run cannot drift apart.
            String retrialWhy = jline.solvers.mam.SolverMAM.mamRetrialRefusal(sn);
            if (!retrialWhy.isEmpty()) {
                throw new RuntimeException(retrialWhy);
            }
            InputOutput.line_debug(options.verbose, "Using retrial method for BMAP/PH/N/N bufferless topology");
            result = Solver_mam_retrial.solver_mam_retrial(sn, options);
            result.method = "retrial";
        } else {
            throw new RuntimeException("Unknown method");
        }

        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                for (int j = 0; j < result.TN.getNumCols(); j++) {
                    result.TN.set(i, j, sn.rates.get(i, j));
                }
            }
        }

        // Handle self-looping classes: override metrics if the method produced
        // incorrect values (Inf/0). Mirrors MATLAB solver_mam_ag.m lines 675-725.
        int M = sn.nstations;
        int K = sn.nclasses;
        if (sn.isslc != null) {
            boolean hasSlc = false;
            for (int r = 0; r < K; r++) {
                if (sn.isslc.get(r) == 1.0) {
                    hasSlc = true;
                    break;
                }
            }
            if (hasSlc) {
                for (int r = 0; r < K; r++) {
                    if (sn.isslc.get(r) == 1.0) {
                        int refst = (int) sn.refstat.get(r);
                        if (refst >= 0 && refst < M) {
                            double qVal = result.QN.get(refst, r);
                            // Override if QN is 0, Inf, or NaN (method didn't compute correctly)
                            if (qVal == 0.0 || Double.isInfinite(qVal) || Double.isNaN(qVal)) {
                                // Clear SLC metrics at all stations first
                                for (int i = 0; i < M; i++) {
                                    result.QN.set(i, r, 0.0);
                                    result.UN.set(i, r, 0.0);
                                    result.TN.set(i, r, 0.0);
                                    result.RN.set(i, r, 0.0);
                                }
                                // All jobs stay at reference station
                                result.QN.set(refst, r, sn.njobs.get(r));
                                double muIr = sn.rates.get(refst, r);
                                if (!Double.isNaN(muIr) && muIr > 0) {
                                    double nservers = sn.nservers.get(refst);
                                    if (Double.isInfinite(nservers)) {
                                        // Delay (infinite server)
                                        result.UN.set(refst, r, result.QN.get(refst, r));
                                        result.TN.set(refst, r, muIr * result.QN.get(refst, r));
                                    } else {
                                        // Queue (finite server): use remaining capacity
                                        double otherUtil = 0.0;
                                        for (int s = 0; s < K; s++) {
                                            if (s != r && sn.isslc.get(s) != 1.0) {
                                                otherUtil += result.UN.get(refst, s);
                                            }
                                        }
                                        double remainingCapacity = Math.max(0.0, 1.0 - otherUtil);
                                        double slcDemand = result.QN.get(refst, r) / muIr;
                                        result.UN.set(refst, r, Math.min(slcDemand, remainingCapacity));
                                        result.TN.set(refst, r, muIr * result.UN.get(refst, r));
                                    }
                                }
                            }
                        }
                    }
                }
                // Recompute response times from Little's law for all entries
                for (int i = 0; i < M; i++) {
                    for (int r = 0; r < K; r++) {
                        if (result.TN.get(i, r) > 0) {
                            result.RN.set(i, r, result.QN.get(i, r) / result.TN.get(i, r));
                        } else {
                            result.RN.set(i, r, 0.0);
                        }
                    }
                }
            }
        }

        result.QN.setNaNToZero();
        result.CN.setNaNToZero();
        result.RN.setNaNToZero();
        result.UN.setNaNToZero();
        result.XN.setNaNToZero();
        result.TN.setNaNToZero();
        long finish = System.nanoTime();
        result.runtime = (finish - start) / 1000000000.0;

        return result;
    }
}
