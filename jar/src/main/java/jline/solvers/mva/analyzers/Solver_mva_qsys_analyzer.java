/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import java.util.Map;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.apache.commons.math3.util.FastMath;

import jline.api.qsys.Qsys_gig1_approx_allencunneen;
import jline.api.qsys.Qsys_gig1_approx_heyman;
import jline.api.qsys.Qsys_gig1_approx_klb;
import jline.api.qsys.Qsys_gig1_approx_kobayashi;
import jline.api.qsys.Qsys_gig1_approx_marchal;
import jline.api.qsys.Qsys_gig1_ubnd_kingman;
import jline.api.qsys.Qsys_erlanga;
import jline.api.qsys.Qsys_gig1_bnds_extremal;
import jline.api.qsys.Qsys_ggnm_diffusion;
import jline.api.qsys.Qsys_gigk_approx;
import jline.api.qsys.Qsys_gigk_approx_kingman;
import jline.api.qsys.Qsys_gigk_approx_whitt;
import jline.api.qsys.Qsys_mg1k_loss_mgs;
import jline.api.qsys.Qsys_mgisrgi_whitt;
import jline.api.qsys.Qsys_mmk_qed;
import jline.api.qsys.QsysAbandonResult;
import jline.api.sn.SnIsMm1kLoss;
import jline.api.sn.SnPatienceHandles;
import jline.api.qsys.Qsys_gm1;
import jline.api.qsys.Qsys_mg1;
import jline.api.qsys.Qsys_mg1_prio;
import jline.api.qsys.Qsys_mm1;
import jline.api.qsys.Qsys_mmk;
import jline.api.qsys.Qsys_mxm1;
import jline.api.qsys.Qsys_phm1;
import jline.api.mam.Map_pie;
import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mva_qsys_analyzer {
    private Solver_mva_qsys_analyzer() {}

    /**
     * MVA Query System analyzer
     */
    public static MVAResult solver_mva_qsys_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        String method = options.method;
        Matrix QN = new Matrix(sn.nstations, sn.nclasses);
        Matrix UN = new Matrix(sn.nstations, sn.nclasses);
        Matrix RN = new Matrix(sn.nstations, sn.nclasses);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix CN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        Matrix WN = new Matrix(sn.nstations, sn.nclasses);
        Matrix XN = new Matrix(sn.nstations, sn.nclasses);
        double lG = Double.NaN;
        int it = 1;

        int source_ist = -1;
        int queue_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }

        // For single-chain qsys models, use chain index 0
        int chainIdx = 0;

        if (sn.visits == null || sn.visits.isEmpty() || sn.visits.get(chainIdx) == null) {
            throw new RuntimeException("Invalid visits matrix: chain " + chainIdx + " not found in visits matrix");
        }
        if (sn.stationToStateful == null || sn.stationToStateful.length() <= queue_ist) {
            throw new RuntimeException("Invalid stationToStateful mapping: queue station " + queue_ist + " not found");
        }

        int statefulIndex = (int) sn.stationToStateful.get(queue_ist);
        Matrix visitsMatrix = sn.visits.get(chainIdx);
        if (visitsMatrix == null || visitsMatrix.length() <= statefulIndex) {
            throw new RuntimeException("Invalid visits matrix: stateful index " + statefulIndex + " not found in visits for chain " + chainIdx);
        }

        double lambda = sn.rates.get(source_ist) * visitsMatrix.get(statefulIndex);
        int k = (int) sn.nservers.get(queue_ist);
        double mu = sn.rates.get(queue_ist);
        double ca = FastMath.sqrt(sn.scv.get(source_ist));
        double cs = FastMath.sqrt(sn.scv.get(queue_ist));
        // null unless the queue reneges
        SnPatienceHandles.Handles hpat = SnPatienceHandles.snPatienceHandles(sn, queue_ist, 0);

        // Finite-capacity loss branch (M/M/1/K with tail drop), the ONE finite
        // buffer this solver honours. The loss probability is the moment-based
        // (MacGregor Smith) Qsys_mg1k_loss_mgs, exact only at scv=1; the queue
        // length comes from the truncated M/M/1/K distribution. Being an
        // approximation in general it is not offered under method='exact', and
        // finiteCapacityReason exempts exactly this shape for every other name,
        // so the gate and the branch have to agree on which models arrive here.
        if (SnIsMm1kLoss.snIsMm1kLoss(sn)) {
            if ("exact".equals(method)) {
                throw new RuntimeException("M/M/1/K tail-drop is solved by the approximate "
                        + "'mg1k.mgs' method (MacGregor Smith); it is not available under "
                        + "method='exact'. Use the default method, or SolverCTMC/SolverNC for "
                        + "an exact result.");
            }
            double Kcap = sn.cap.get(queue_ist);
            double rhoK = lambda / mu;
            double Ploss = (Double) Qsys_mg1k_loss_mgs
                    .qsys_mg1k_loss_mgs(lambda, mu, cs * cs, (int) Math.round(Kcap)).get("lossprob");
            double Tq = lambda * (1.0 - Ploss);   // carried throughput
            double Uq = Tq / mu;                  // single-server utilization
            double Lsys;
            if (Math.abs(rhoK - 1.0) < 1e-10) {
                Lsys = Kcap / 2.0;                // L'Hopital limit at rho=1
            } else {
                double rKp1 = FastMath.pow(rhoK, Kcap + 1.0);
                Lsys = rhoK / (1.0 - rhoK) - (Kcap + 1.0) * rKp1 / (1.0 - rKp1);
            }
            double visitsK = visitsMatrix.get(statefulIndex);
            double Rq = Lsys / Tq;                // per-visit response time, by Little
            for (int r = 0; r < sn.nclasses; r++) {
                RN.set(queue_ist, r, Rq);
                QN.set(queue_ist, r, Lsys);
                UN.set(queue_ist, r, Uq);
                TN.set(queue_ist, r, Tq);         // carried (effective) rate
                TN.set(source_ist, r, lambda);    // offered arrival rate
                XN.set(queue_ist, r, Tq);         // system throughput = carried rate
                CN.set(queue_ist, r, Rq * visitsK);
            }
            long endLoss = System.nanoTime();
            res.QN = QN;
            res.UN = UN;
            res.RN = RN;
            res.TN = TN;
            res.CN = CN;
            res.XN = XN;
            res.AN = AN;
            res.WN = WN;
            res.logNormConstAggr = 0.0;
            res.runtime = (endLoss - startTime) / 1000000000.0;
            res.iter = it;
            res.method = "mg1k.mgs";
            return res;
        }

        // Check if this is a BMAP arrival process (MX/M/1)
        Station sourceStation = sn.stations.get(source_ist);
        JobClass jobClass = sn.jobclasses.get(0);
        ProcessType sourceProcType = null;
        if (sn.procid != null) {
            Map<JobClass, ProcessType> innerMap = sn.procid.get(sourceStation);
            if (innerMap != null) {
                sourceProcType = innerMap.get(jobClass);
            }
        }
        boolean isBMAP = sourceProcType == ProcessType.BMAP;

        if ("exact".equals(method)) {
            if (isBMAP && cs == 1.0 && k == 1) {
                method = "mxm1";
            } else if (ca == 1.0 && cs == 1.0 && k == 1) {
                method = "mm1";
            } else if (ca == 1.0 && cs == 1.0 && k > 1) {
                method = "mmk";
            } else if (ca == 1.0 && k == 1) {
                method = "mg1";
            } else if (cs == 1.0 && k == 1) {
                method = "gm1";
            } else {
                throw new RuntimeException("MVA exact method unavailable for this model.");
            }
        }
        if ("default".equals(method)) {
            if (hpat != null) {
                // A station customers walk away from is a different model, not
                // a correction to one: nothing in the G/G/k family below carries
                // an abandonment rate, so the choice is made here and not by
                // ca/cs.
                method = hpat.isExponential ? "erlanga" : "mgisrgi";
            } else if (isBMAP && cs == 1.0 && k == 1) {
                method = "mxm1";
            } else if (ca == 1.0 && cs == 1.0 && k == 1) {
                method = "mm1";
            } else if (ca == 1.0 && cs == 1.0 && k > 1) {
                method = "mmk";
            } else if (ca == 1.0 && k == 1) {
                method = "mg1";
            } else if (cs == 1.0 && k == 1) {
                method = "gm1";
            } else if (k > 1) {
                method = "gigk";
            } else {
                method = "gig1.klb";
            }
        }
        double R = 0.0;
        double lambdaEffective = lambda;

        // Whitt family, full metric set. These methods answer a station whose
        // CARRIED throughput is below the offered rate -- customers abandon, or
        // are blocked -- so Little's law on lambda would silently overstate the
        // queue and the common tail below cannot be used.
        if ("erlanga".equals(method) || "mgisrgi".equals(method)
                || "gigk.diffusion".equals(method)) {
            double cap = sn.cap != null ? sn.cap.get(queue_ist) : Double.POSITIVE_INFINITY;
            // An uncapped station carries Integer.MAX_VALUE here, not Inf as in
            // MATLAB and Python; taking it literally asks for an array of that
            // length. See _kb/04-networkstruct.md (cap).
            boolean unbounded = !Double.isFinite(cap) || cap >= Integer.MAX_VALUE;
            // waiting spaces, servers excluded
            double room = unbounded ? Double.POSITIVE_INFINITY : Math.max(0.0, cap - k);
            double Lsys;
            double Tq;
            double Uq;
            if ("gigk.diffusion".equals(method)) {
                Map<String, Double> dif = Qsys_ggnm_diffusion.qsys_ggnm_diffusion(lambda, mu, k, room, ca, cs);
                Lsys = dif.get("meanNumber");
                Tq = dif.get("throughput");
                Uq = dif.get("utilization");
            } else {
                if (hpat == null) {
                    throw new RuntimeException("method '" + method
                            + "' needs a reneging patience law on the queue.");
                }
                QsysAbandonResult ab;
                if ("erlanga".equals(method) || hpat.isExponential) {
                    // Exponential patience makes the state-dependent
                    // approximation exact, so take the exact chain either way.
                    ab = Qsys_erlanga.qsys_erlanga(lambda, mu, hpat.rate, k, room);
                } else {
                    ab = Qsys_mgisrgi_whitt.qsys_mgisrgi_whitt(lambda, mu, k, room, hpat.asPatience());
                }
                Lsys = ab.meanNumber;
                Tq = ab.throughput;
                Uq = ab.utilization;
            }
            double visitsQ = visitsMatrix.get(statefulIndex);
            // Little's law on the CARRIED rate, as SolverCTMC reports it.
            double Rq = Tq > 0 ? Lsys / Tq : 0.0;
            for (int r = 0; r < sn.nclasses; r++) {
                RN.set(queue_ist, r, Rq);
                QN.set(queue_ist, r, Lsys);
                UN.set(queue_ist, r, Uq);
                TN.set(queue_ist, r, Tq);
                TN.set(source_ist, r, lambda / visitsQ);
                XN.set(queue_ist, r, Tq);
                AN.set(queue_ist, r, lambda);
                CN.set(queue_ist, r, Rq * visitsQ);
            }
            long endAbandon = System.nanoTime();
            res.QN = QN;
            res.UN = UN;
            res.RN = RN;
            res.TN = TN;
            res.CN = CN;
            res.XN = XN;
            res.AN = AN;
            res.WN = WN;
            res.logNormConstAggr = 0.0;
            res.runtime = (endAbandon - startTime) / 1000000000.0;
            res.iter = it;
            res.method = method;
            return res;
        }

        if ("mm1".equals(method)) {
            Qsys_mm1.qsys_mm1(lambda, mu);
            R = Ret.qsys.W;
        } else if ("rqt".equals(method)) {
            // Robust Queueing Theory single-queue solution: the arrival and
            // service flows enter as polyhedral uncertainty sets
            double rho1 = lambda / (k * mu);
            double gammaA = ca / lambda;
            double gammaS = jline.api.qsys.Qsys_gigk_rqt_gamma.qsys_gigk_rqt_gamma(rho1, mu, gammaA, cs / mu, k);
            R = jline.api.qsys.Qsys_gigk_rqt.qsys_gigk_rqt(lambda, mu, gammaA, gammaS, k, 2.0, 2.0)[0];
        } else if ("rqna".equals(method)) {
            // Robust Queueing (RQ) single-queue solution: characterize the
            // arrival flow by its index of dispersion for counts (IDC).
            MatrixCell arvMAP = null;
            if (sn.proc != null) {
                Map<JobClass, MatrixCell> innerMap = sn.proc.get(sourceStation);
                if (innerMap != null) {
                    arvMAP = innerMap.get(jobClass);
                }
            }
            if (arvMAP == null) {
                throw new RuntimeException("RQNA arrival process not found for source station");
            }
            final MatrixCell arvMAPf = arvMAP;
            double rho1 = lambda / mu;
            java.util.function.DoubleUnaryOperator IaFun1 = new java.util.function.DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double x) {
                    return jline.api.mam.Map_count_idc.map_count_idc(arvMAPf, x);
                }
            };
            double[] zwqx = jline.api.qsys.Qsys_gig1_rq.qsys_gig1_rq(rho1, mu, cs * cs, IaFun1);
            R = zwqx[1] + 1.0 / mu;
        } else if ("mxm1".equals(method)) {
            MatrixCell proc = null;
            if (sn.proc != null) {
                Map<JobClass, MatrixCell> innerMap = sn.proc.get(sourceStation);
                if (innerMap != null) {
                    proc = innerMap.get(jobClass);
                }
            }
            if (proc == null) {
                throw new RuntimeException("BMAP process not found for source station");
            }

            Matrix D0 = proc.get(0);
            int maxBatchSize = proc.size() - 2;

            double lambdaBatch = -D0.get(0, 0);

            double meanBatchSize = 0.0;
            double secondMomentBatchSize = 0.0;
            double totalRate = 0.0;

            for (int idx = 2; idx < proc.size(); idx++) {
                int batchSize = idx - 1;
                Matrix Dk = proc.get(idx);
                double rate_k = Dk.get(0, 0);
                totalRate += rate_k;
                meanBatchSize += batchSize * rate_k;
                secondMomentBatchSize += batchSize * batchSize * rate_k;
            }

            if (totalRate > 0) {
                meanBatchSize /= totalRate;
                secondMomentBatchSize /= totalRate;
            }

            Qsys_mxm1.qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize);
            R = Ret.qsys.W;

            lambdaEffective = lambdaBatch * meanBatchSize;

            method = "mxm1";
        } else if ("mmk".equals(method)) {
            Qsys_mmk.qsys_mmk(lambda, mu, k);
            R = Ret.qsys.W;
        } else if ("mg1".equals(method) || "mgi1".equals(method)) {
            Qsys_mg1.qsys_mg1(lambda, mu, cs);
            R = Ret.qsys.W;
        } else if ("gigk".equals(method)) {
            Qsys_gigk_approx.qsys_gigk_approx(lambda, mu, ca, cs, k);
            R = Ret.qsys.W;
        } else if ("gigk.kingman_approx".equals(method)) {
            Qsys_gigk_approx_kingman.qsys_gigk_approx_kingman(lambda, mu, ca, cs, k);
            R = Ret.qsys.W;
        } else if ("gig1.kingman".equals(method)) {
            Qsys_gig1_ubnd_kingman.qsys_gig1_ubnd_kingman(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gig1.gelenbe".equals(method)) {
            R = (Double) jline.api.qsys.Qsys_gig1_approx_gelenbe
                    .qsys_gig1_approx_gelenbe(lambda, mu, ca, cs).get("W");
        } else if ("gig1.kimura".equals(method)) {
            R = (Double) jline.api.qsys.Qsys_gig1_approx_kimura
                    .qsys_gig1_approx_kimura(lambda, mu, ca, cs).get("W");
        } else if ("gig1.heyman".equals(method)) {
            Qsys_gig1_approx_heyman.qsys_gig1_approx_heyman(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gig1".equals(method) || "gig1.allen".equals(method)) {
            Qsys_gig1_approx_allencunneen.qsys_gig1_approx_allencunneen(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gig1.kobayashi".equals(method)) {
            Qsys_gig1_approx_kobayashi.qsys_gig1_approx_kobayashi(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gig1.klb".equals(method)) {
            Qsys_gig1_approx_klb.qsys_gig1_approx_klb(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gig1.marchal".equals(method)) {
            Qsys_gig1_approx_marchal.qsys_gig1_approx_marchal(lambda, mu, ca, cs);
            R = Ret.qsys.W;
        } else if ("gigk.whitt".equals(method)) {
            R = (Double) Qsys_gigk_approx_whitt.qsys_gigk_approx_whitt(lambda, mu, ca, cs, k).get("W");
        } else if ("qed".equals(method)) {
            R = Qsys_mmk_qed.qsys_mmk_qed(lambda, mu, k).get("meanWait") + 1.0 / mu;
        } else if ("gig1.extremal".equals(method)) {
            // The upper end, gig1.kingman already reporting a bound.
            R = Qsys_gig1_bnds_extremal.qsys_gig1_bnds_extremal(lambda, mu, ca, cs).get("upperBound")
                    + 1.0 / mu;
        } else if ("gm1".equals(method) || "gim1".equals(method)) {
            // see _kb/06-solver-catalog.md for rationale
            Double Rgm1 = null;
            try {
                Map<JobClass, MatrixCell> procMap = (sn.proc != null) ? sn.proc.get(sourceStation) : null;
                MatrixCell arrCell = (procMap != null) ? procMap.get(jobClass) : null;
                // see _kb/06-solver-catalog.md for rationale
                if (ProcessType.isMarkovian(sourceProcType)
                        && arrCell != null && arrCell.size() >= 2
                        && arrCell.get(0) != null && arrCell.get(1) != null
                        && arrCell.get(0).getNumRows() == arrCell.get(0).getNumCols()
                        && arrCell.get(0).getNumRows() > 0
                        && !arrCell.get(0).hasNaN()) {
                    Matrix D0 = arrCell.get(0);
                    Matrix pieM = Map_pie.map_pie(D0, arrCell.get(1));
                    int n = D0.getNumRows();
                    double[] alpha = new double[n];
                    double[][] Tm = new double[n][n];
                    for (int a = 0; a < n; a++) {
                        alpha[a] = pieM.get(a);
                        for (int b = 0; b < n; b++) {
                            Tm[a][b] = D0.get(a, b);
                        }
                    }
                    Rgm1 = Qsys_phm1.qsys_phm1(alpha, Tm, mu).getMeanSojournTime();
                }
            } catch (Exception phEx) {
                Rgm1 = null;
            }
            // see _kb/06-solver-catalog.md for rationale
            if (Rgm1 == null && sourceProcType == ProcessType.DET && lambda > 0.0) {
                final double meanIa = 1.0 / lambda;
                UnivariateFunction fdet = new UnivariateFunction() {
                    @Override
                    public double value(double x) {
                        return FastMath.exp(-(mu * (1.0 - x)) * meanIa) - x;
                    }
                };
                try {
                    double sigma = gm1CaudalSigma(fdet);
                    if (!Double.isNaN(sigma)) {
                        Qsys_gm1.qsys_gm1(sigma, mu);
                        Rgm1 = Ret.qsys.W;
                    }
                } catch (Exception detEx) {
                    Rgm1 = null;
                }
            }
            if (Rgm1 == null && sn.lst != null && sn.lst.get(sourceStation) != null
                    && sn.lst.get(sourceStation).get(jobClass) != null) {
                try {
                    // sn.lst is COMPLEX, since transform inversion and root
                    // location need arguments off the real axis; the sigma-root
                    // below walks the real line, so it passes Complex(s, 0).
                    final SerializableFunction<org.apache.commons.math3.complex.Complex,
                            org.apache.commons.math3.complex.Complex> F =
                            sn.lst.get(sourceStation).get(jobClass);
                    final UnivariateFunction LA = new UnivariateFunction() {
                        @Override
                        public double value(double s) {
                            return F.apply(new org.apache.commons.math3.complex.Complex(s, 0.0)).getReal();
                        }
                    };
                    final double muFinal = mu;
                    UnivariateFunction func = new UnivariateFunction() {
                        @Override
                        public double value(double x) {
                            return LA.value(muFinal - muFinal * x) - x;
                        }
                    };
                    // see _kb/06-solver-catalog.md for rationale
                    double sigma = gm1CaudalSigma(func);
                    if (!Double.isNaN(sigma)) {
                        Qsys_gm1.qsys_gm1(sigma, mu);
                        Rgm1 = Ret.qsys.W;
                    }
                } catch (Exception lstEx) {
                    Rgm1 = null;
                }
            }
            if (Rgm1 == null) {
                // No usable PH or LST (e.g. Det/Uniform): G/G/1 KLB approximation.
                Qsys_gig1_approx_klb.qsys_gig1_approx_klb(lambda, mu, ca, cs);
                Rgm1 = Ret.qsys.W;
            }
            R = Rgm1;
        } else {
            throw new RuntimeException("Unsupported method for a model with 1 station and 1 class.");
        }

        double visits = visitsMatrix.get(statefulIndex);
        for (int r = 0; r < sn.nclasses; r++) {
            // see _kb/06-solver-catalog.md for rationale
            RN.set(queue_ist, r, R);
            CN.set(queue_ist, r, R * visits);
            XN.set(queue_ist, r, lambdaEffective);
            UN.set(queue_ist, r, lambdaEffective / mu / k);
            // see _kb/06-solver-catalog.md for rationale
            TN.set(source_ist, r, lambdaEffective / visits);
            TN.set(queue_ist, r, lambdaEffective);
            // Station queue length by Little's law at the station: offered
            // rate x per-visit residence time.
            QN.set(queue_ist, r, XN.get(queue_ist, r) * RN.get(queue_ist, r));
        }
        lG = 0.0;
        long endTime = System.nanoTime();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000000000.0;
        res.iter = it;
        res.method = method;
        return res;
    }

    // First (caudal / non-trivial) root in (0,1) of the G/M/1 fixed-point residual
    // f(sigma) = LST(mu*(1-sigma)) - sigma. The trivial root is sigma=1, and f can
    // be positive at both ends of (0,1), so bracketing the whole interval misses
    // the caudal root; scan for the FIRST sign change instead (matches MATLAB
    // fzero(@(x) LA(mu-mu*x)-x, 0.5)). Returns NaN if no interior root is found.
    private static double gm1CaudalSigma(UnivariateFunction f) {
        int n = 400;
        double loB = 1e-9;
        double hiB = 1.0 - 1e-6;
        double prevX = loB;
        double prevF = f.value(prevX);
        for (int i = 1; i < n; i++) {
            double x = loB + (hiB - loB) * i / (n - 1);
            double fx = f.value(x);
            if (prevF * fx < 0.0) {
                // MATLAB's fzero converges to machine precision, and the parity
                // tolerance on the resulting W is itself 1e-6 relative, so the
                // BrentSolver DEFAULT absolute accuracy (1e-6) is the same order
                // as the quantity being asserted. Ask for the root properly.
                return new BrentSolver(1e-14, 1e-14).solve(1000, f, prevX, x);
            }
            prevX = x;
            prevF = fx;
        }
        return Double.NaN;
    }

    /**
     * Exact non-preemptive priority (HOL) analyzer for a single open M/G/1
     * queue with Poisson per-class arrivals: dispatches to the Cobham formula
     * (Qsys_mg1_prio) instead of the AMVA preemptive shadow-server
     * approximation, which underestimates the waiting time of every class.
     */
    public static MVAResult solver_mva_qsys_prio_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        int R = sn.nclasses;
        Matrix QN = new Matrix(sn.nstations, R);
        Matrix UN = new Matrix(sn.nstations, R);
        Matrix RN = new Matrix(sn.nstations, R);
        Matrix TN = new Matrix(sn.nstations, R);
        Matrix CN = new Matrix(sn.nstations, R);
        Matrix AN = new Matrix(sn.nstations, R);
        Matrix WN = new Matrix(sn.nstations, R);
        Matrix XN = new Matrix(sn.nstations, R);

        int source_ist = -1;
        int queue_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }

        // Order classes by priority (0 = highest priority first, stable)
        Integer[] order = new Integer[R];
        for (int r = 0; r < R; r++) order[r] = r;
        final Matrix prios = sn.classprio;
        java.util.Arrays.sort(order, (a, b) -> Double.compare(prios.get(a), prios.get(b)));

        Matrix lambdaOrd = new Matrix(1, R);
        Matrix muOrd = new Matrix(1, R);
        Matrix csOrd = new Matrix(1, R);
        for (int j = 0; j < R; j++) {
            int r = order[j];
            lambdaOrd.set(0, j, sn.rates.get(source_ist, r));
            muOrd.set(0, j, sn.rates.get(queue_ist, r));
            double scv = sn.scv.get(queue_ist, r);
            csOrd.set(0, j, (Double.isFinite(scv) && scv > 0) ? FastMath.sqrt(scv) : 1.0);
        }

        Ret.qsys_prio prio = Qsys_mg1_prio.qsys_mg1_prio(lambdaOrd, muOrd, csOrd);
        Matrix W = Ret.qsys_prio.W;

        for (int j = 0; j < R; j++) {
            int r = order[j];
            double lam = lambdaOrd.get(0, j);
            double mu_r = muOrd.get(0, j);
            double W_r = W.get(j);
            RN.set(queue_ist, r, W_r);
            CN.set(queue_ist, r, W_r);
            XN.set(queue_ist, r, lam);
            UN.set(queue_ist, r, lam / mu_r);
            TN.set(queue_ist, r, lam);
            AN.set(queue_ist, r, lam);
            QN.set(queue_ist, r, lam * W_r);
            TN.set(source_ist, r, lam);
            AN.set(source_ist, r, lam);
        }

        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = 0.0;
        res.iter = 0;
        res.method = "mg1.prio";
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return res;
    }

    /**
     * Numerically exact analyzer for a single open M/M/1-DPS queue: dispatches
     * to the truncated-CTMC DPS solver (Qsys_mm1_dps). The AMVA-DPS cross-term
     * correction violates the equal-rate conservation law (the total count must
     * equal the M/M/1 value when all service rates are equal).
     */
    public static MVAResult solver_mva_qsys_dps_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        int R = sn.nclasses;
        Matrix QN = new Matrix(sn.nstations, R);
        Matrix UN = new Matrix(sn.nstations, R);
        Matrix RN = new Matrix(sn.nstations, R);
        Matrix TN = new Matrix(sn.nstations, R);
        Matrix CN = new Matrix(sn.nstations, R);
        Matrix AN = new Matrix(sn.nstations, R);
        Matrix WN = new Matrix(sn.nstations, R);
        Matrix XN = new Matrix(1, R);

        int source_ist = -1;
        int queue_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }

        Matrix lambda = new Matrix(1, R);
        Matrix mu = new Matrix(1, R);
        Matrix w = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            lambda.set(0, r, sn.rates.get(source_ist, r));
            mu.set(0, r, sn.rates.get(queue_ist, r));
            double wr = (sn.schedparam != null && !sn.schedparam.isEmpty())
                    ? sn.schedparam.get(queue_ist, r) : 1.0;
            w.set(0, r, wr > 0 ? wr : 1.0);
        }

        Matrix T = jline.api.qsys.Qsys_mm1_dps.qsys_mm1_dps(lambda, mu, w);

        for (int r = 0; r < R; r++) {
            double lam = lambda.get(0, r);
            if (lam <= 0) continue;
            double T_r = T.get(0, r);
            RN.set(queue_ist, r, T_r);
            CN.set(queue_ist, r, T_r);
            XN.set(0, r, lam);
            UN.set(queue_ist, r, lam / mu.get(0, r));
            TN.set(queue_ist, r, lam);
            AN.set(queue_ist, r, lam);
            QN.set(queue_ist, r, lam * T_r);
            TN.set(source_ist, r, lam);
            AN.set(source_ist, r, lam);
        }

        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = 0.0;
        res.iter = 0;
        res.method = "mm1.dps";
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return res;
    }
}
