package jline.solvers.mam.analyzers;

import jline.api.qsys.Qsys_is_retrial;
import jline.api.qsys.RetrialInfo;
import jline.api.sn.SnHasForkJoin;
import jline.api.sn.SnIsClosedModel;
import jline.api.sn.SnIsOpenModel;
import jline.api.sn.SnNonmarkovToPh;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.solvers.mam.handlers.Solver_mam;
import jline.solvers.mam.handlers.Solver_mam_ag;
import jline.solvers.mam.handlers.Solver_mam_basic;
import jline.solvers.mam.handlers.Solver_mam_basic_mmap;
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
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.INF) {
                nDelay++;
            } else if (sched == SchedStrategy.FCFS) {
                nQueue++;
            }
        }
        return nDelay == 1 && nQueue == 1;
    }

    public static MAMResult solver_mam_analyzer(NetworkStruct snInput, SolverOptions options) {
        long start = System.nanoTime();

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

        // Convert non-Markovian distributions to PH
        NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(snInput, options);

        // Check if the model is mixed (has both open and closed classes)
        boolean isOpen = SnIsOpenModel.snIsOpenModel(sn);
        boolean isClosed = SnIsClosedModel.snIsClosedModel(sn);
        boolean isMixed = isOpen && isClosed;

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
        } else if ("inap".equals(options.method) || "inapplus".equals(options.method)
                || "inapinf".equals(options.method) || "exact".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using RCAT method: " + options.method);
            result = Solver_mam_ag.solver_mam_ag(sn, options);
            result.method = options.method;
        } else if ("ldqbd".equals(options.method)) {
            InputOutput.line_debug(options.verbose, "Using LDQBD method for single-class closed network");
            result = Solver_mam_ldqbd.solver_mam_ldqbd(sn, options);
            result.method = "ldqbd";
        } else if ("retrial".equals(options.method)) {
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
