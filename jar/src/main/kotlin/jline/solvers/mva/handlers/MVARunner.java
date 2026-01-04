/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva.handlers;

import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.AvgHandle;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.SolverAvgHandles;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.SolverMVA;
import jline.lang.FeatureSet;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.api.sn.SnGetResidTFromRespTKt.snGetResidTFromRespT;
import static jline.api.sn.SnHasPollingKt.snHasPolling;
import static jline.io.InputOutputKt.*;
import static jline.solvers.mva.analyzers.Solver_mva_analyzerKt.solver_mva_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_bound_analyzerKt.solver_mva_bound_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_cache_analyzerKt.solver_mva_cache_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_cacheqn_analyzerKt.solver_mva_cacheqn_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_polling_analyzerKt.solver_mva_polling_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_qsys_analyzerKt.solver_mva_qsys_analyzer;
import static jline.solvers.mva.analyzers.Solver_mvald_analyzerKt.solver_mvald_analyzer;
import static jline.lang.ModelAdapter.*;
import static jline.util.Utils.isInf;

public class MVARunner {

    private final Network model;
    private final boolean enableChecks;
    protected NetworkStruct sn;
    protected SolverOptions options;
    protected MVAResult res;

    public MVARunner(Network model, SolverOptions options, boolean enableChecks) {
        this.model = model;
        this.options = options;
        this.enableChecks = enableChecks;
        this.sn = this.model.getStruct(false);
        this.res = null;
    }

    /**
     * runAnalyzer() method from LINE.
     *
     * @param avgHandles - the average handles for the model
     * @return - the performance measures corresponding to the given network
     */
    public SolverResult runAnalyzer(SolverAvgHandles avgHandles) {
        long T0 = System.nanoTime();
        int iter = 0;

        if (this.enableChecks && !FeatureSet.supports(SolverMVA.getFeatureSet(), this.model.getUsedLangFeatures())) {
            line_error(mfilename(new Object() {
            }), "This model contains features not supported by the solver.");
        }

        // Case 'java' or 'jline.amva' can be ignored, so we can remove the switch
        this.sn = this.model.getStruct(false);
        boolean forkLoop = true;
        int forkIter = 0;
        int forkNodes = 0;
        for (NodeType n : this.sn.nodetype) {
            if (n == NodeType.Fork)
                forkNodes++;
        }
        Matrix forkLambda = new Matrix(1, 2 * this.sn.nclasses * forkNodes);
        forkLambda.fill(GlobalConstants.FineTol);
        Matrix QN = new Matrix(1, this.sn.nclasses);
        QN.fill(GlobalConstants.Immediate);
        Matrix QN_1 = new Matrix(1, this.sn.nclasses);
        Matrix UN = new Matrix(1, this.sn.nclasses);
        boolean forceOneMoreIteration = false;
        MVAResult ret = new MVAResult();
        Network nonfjmodel = null;
        Matrix fjclassmap, fjforkmap;
        fjclassmap = fjforkmap = null;
        Map<Integer, Integer> fj_auxiliary_delays, fanout;
        fj_auxiliary_delays = null;
        fanout = null;
        Matrix outerForks = null;
        Matrix parentForks = null;
        while (forkLoop && forkIter < this.options.iter_max) {
            // Clear out any fields from the past returns
            ret = new MVAResult();
            if (this.model.hasFork()) {
                forkIter += 1;
                if (forkIter == 1) {
                    Ret.FJApprox approxReturn = null;
                    switch (options.config.fork_join) {
                        case "heidelberger-trivedi":
                        case "ht":
                            approxReturn = ht(this.model);
                            break;
                        case "mmt":
                        case "default":
                        case "fjt":
                            approxReturn = mmt(this.model, forkLambda);
                            Ret.FJsortForks sortForksReturn = sort_forks(sn, approxReturn.nonfjmodel.getStruct(false),
                                    approxReturn.fjforkmap, approxReturn.fjclassmap, approxReturn.nonfjmodel);
                            outerForks = sortForksReturn.outerForks;
                            parentForks = sortForksReturn.parentForks;
                            break;
                    }
                    nonfjmodel = approxReturn.nonfjmodel;
                    fjclassmap = approxReturn.fjclassmap;
                    fjforkmap = approxReturn.fjforkmap;
                    fj_auxiliary_delays = approxReturn.fj_auxiliary_delays;
                    fanout = approxReturn.fanout;
                } else if (!options.config.fork_join.equals("heidelberger-trivedi") && !options.config.fork_join.equals("ht")) {
                    for (int r = 0; r < fjclassmap.length(); r++) {
                        int s = (int) fjclassmap.get(r);
                        if (s > -1) {
                            Source nonfjSource = nonfjmodel.getSource();
                            if (fanout.get(r) > 0) {
                                if (!nonfjSource.getArrivalProcess(nonfjmodel.getClasses().get(r)).isDisabled()) {
                                    Exp e = (Exp) nonfjSource.getArrivalProcess(nonfjmodel.getClasses().get(r));
                                    e.updateRate((fanout.get(r) - 1) * forkLambda.get(r));
                                }
                            }
                        }
						/* Changed from refreshStruct(true) to refreshRates(null, null) for performance optimization.
						   This change has been validated and produces equivalent results while being more efficient.
						   refreshRates only updates rate/service information rather than rebuilding the entire structure.
						   */
                        //nonfjmodel.refreshStruct(true);
                        nonfjmodel.refreshRates(null, null);
                    }
                }
                this.sn = nonfjmodel.getStruct(false);
                if (forkIter > 2) {
                    Matrix convCheck = Matrix.ones(QN_1.getNumRows(), QN_1.getNumCols()).sub(1, QN_1.elementDiv(QN));
                    double maxAbs = -1;
                    for (int i = 0; i < convCheck.getNumRows(); i++) {
                        for (int j = 0; j < convCheck.getNumCols(); j++) {
                            if (!Double.isNaN(convCheck.get(i, j)) && FastMath.abs(convCheck.get(i, j)) > maxAbs) {
                                maxAbs = FastMath.abs(convCheck.get(i, j));
                            }
                        }
                    }
                    if (maxAbs < GlobalConstants.CoarseTol) {
                        forkLoop = false;
                    } else {
                        QN_1 = QN;
                    }
                } else {
                    if (this.model.hasOpenClasses()) {
                        int sourceIndex = this.model.getIndexSourceNode();
                        Matrix UNnosource = new Matrix(UN);
                        for (int i = 0; i < UN.getNumCols(); i++) {
                            UN.set(sourceIndex, i, 0);
                        }
                        Matrix util = UNnosource.sumRows();
                        for (int i = 0; i < util.getNumRows(); i++) {
                            double Uiopen = 0.0;
                            for (int j = 0; j < UN.getNumCols(); j++) {
                                // sum util of open classes
                                if (isInf(sn.njobs.get(j))) {
                                    Uiopen += UN.get(i, j);
                                }
                            }
                            if (Uiopen > 0.99 && this.sn.nservers.get(i) != Integer.MAX_VALUE) {
                                System.out.println("The model may be unstable: the utilization of station " + i + " for open classes exceeds 99 percent.\n");
                            }
                        }
                    }
                    QN_1 = QN;
                }
            } else {
                forkLoop = false;
            }

            if (this.options.method.equals("exact") && !this.model.hasProductFormSolution()) {
                line_error(mfilename(new Object() {
                }), "The exact method requires the model to have a product-form solution. This model does not have one.");
            }
            if (this.options.method.equals("mva") && !this.model.hasProductFormSolution()) {
                line_warning(mfilename(new Object() {
                }), "The exact method requires the model to have a product-form solution. This model does not have one. SolverMVA will return an approximation generated by an exact MVA algorithm.");
            }
            String method = this.options.method;
            line_debug(this.options.verbose, String.format("MVA analyzer starting: method=%s, nclasses=%d, nclosed=%d, nnodes=%d",
                method, this.sn.nclasses, this.sn.nclosedjobs, this.sn.nnodes));

            if (this.sn.nclasses == 1 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink)) {
                // Single-class open queueing system
                line_debug(this.options.verbose, "Detected open queueing system (Source-Queue-Sink), calling solver_mva_qsys_analyzer");
                ret = solver_mva_qsys_analyzer(this.sn, this.options.copy());
            } else if (this.sn.nclasses > 1 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink) && snHasPolling(this.sn)) {
                // Multi-class open polling system
                line_debug(this.options.verbose, "Detected multi-class polling system (Source-Queue-Sink with POLLING scheduling), calling solver_mva_polling_analyzer");
                ret = solver_mva_polling_analyzer(this.sn, this.options.copy());
            } else if (this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Cache, NodeType.Sink)) {
                // Non-rentrant cache
                // Random initialisation
                for (int ind = 0; ind < this.sn.nnodes; ind++) {
                    if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                        Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                        Matrix prob = new Matrix(cacheNode.getHitClass());
                        for (int i = 0; i < prob.getNumRows(); i++) {
                            for (int j = 0; j < prob.getNumCols(); j++) {
                                if (prob.get(i, j) > 0) {
                                    prob.set(i, j, 0.5);
                                }
                            }
                        }
                        cacheNode.setResultHitProb(prob);
                        Matrix missProb = new Matrix(prob.getNumRows(), prob.getNumCols());
                        for (int i = 0; i < prob.getNumRows(); i++) {
                            for (int j = 0; j < prob.getNumCols(); j++) {
                                missProb.set(i, j, 1 - prob.get(i, j));
                            }
                        }
                        cacheNode.setResultMissProb(missProb);
                    }
                }
                this.model.refreshChains(true);
                // Start iteration
                line_debug(this.options.verbose, "Detected cache system (Source-Cache-Sink), calling solver_mva_cache_analyzer");
                ret = solver_mva_cache_analyzer(this.sn, this.options.copy());

                for (int ind = 0; ind < this.sn.nnodes; ind++) {
                    if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                        Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                        Matrix hitClass = cacheNode.getHitClass();
                        Matrix missClass = cacheNode.getMissClass();
                        Matrix hitProb = new Matrix(1, hitClass.length());
                        for (int k = 0; k < hitClass.length(); k++) {
                            int chain_k = 0;
                            for (; chain_k < this.sn.chains.getNumRows(); chain_k++) {
                                if (this.sn.chains.get(chain_k, k) > 0)
                                    break;
                            }
                            Matrix inchain = new Matrix(1, this.sn.chains.getNumCols());
                            for (int i = 0; i < inchain.getNumCols(); i++) {
                                inchain.set(0, i, this.sn.chains.get(chain_k, i) > 0 ? 1 : 0);
                            }
                            int h = (int) hitClass.get(k);
                            int m = (int) missClass.get(k);
                            if (h > -1 && m > -1) {
                                double sumXN = 0;
                                for (int i = 0; i < inchain.getNumCols(); i++) {
                                    if (inchain.get(i) > 0 && !Double.isNaN(ret.XN.get(i))) {
                                        sumXN += ret.XN.get(i);
                                    }
                                }
                                hitProb.set(k, ret.XN.get(h) / sumXN);
                            }
                        }
                        Matrix missProb = new Matrix(1, hitClass.length());
                        for (int i = 0; i < hitClass.length(); i++) {
                            missProb.set(i, 1 - hitProb.get(i));
                        }
                        cacheNode.setResultHitProb(hitProb);
                        cacheNode.setResultMissProb(missProb);
                    }
                }
                this.model.refreshStruct(true);
            } else {
                // Queueing network
                boolean cachePresent = false;
                for (NodeType t : this.sn.nodetype) {
                    if (t == NodeType.Cache) {
                        cachePresent = true;
                        break;
                    }
                }
                if (cachePresent) {
                    // Integrated Cache Queueing
                    line_debug(this.options.verbose, "Detected cache+queueing network, calling solver_mva_cacheqn_analyzer");
                    ret = solver_mva_cacheqn_analyzer(this.sn, this.options.copy());
                    for (int ind = 0; ind < this.sn.nnodes; ind++) {
                        if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                            Cache cache = (Cache) this.model.getNodes().get(ind);
                            cache.setResultHitProb(Matrix.extractRows(ret.hitProb, ind, ind + 1, null));
                            cache.setResultMissProb(Matrix.extractRows(ret.missProb, ind, ind + 1, null));
                        }
                    }
                    this.model.refreshStruct(true);
                } else {
                    // Ordinary queueing network
                    switch (method) {
                        case "aba.upper":
                        case "aba.lower":
                        case "bjb.upper":
                        case "bjb.lower":
                        case "pb.upper":
                        case "pb.lower":
                        case "gb.upper":
                        case "gb.lower":
                        case "sb.upper":
                        case "sb.lower":
                        case "harel.lower":
                        case "harel.upper":
                            line_debug(this.options.verbose, "Using bound analyzer method: " + method + ", calling solver_mva_bound_analyzer");
                            ret = solver_mva_bound_analyzer(this.sn, this.options.copy());
                            break;
                        default: // this is the main solution block for standard queueing networks
                            if ((this.sn.lldscaling != null && !this.sn.lldscaling.isEmpty()) ||
                                    (this.sn.cdscaling != null && !this.sn.cdscaling.isEmpty())) {
                                line_debug(this.options.verbose, "Detected load-dependent scaling, calling solver_mvald_analyzer");
                                ret = solver_mvald_analyzer(this.sn, this.options.copy());
                            } else {
                                line_debug(this.options.verbose, "Using standard MVA method, calling solver_mva_analyzer");
                                ret = solver_mva_analyzer(this.sn, this.options.copy());
                            }
                    }
                }
            }
            if (this.model.hasFork()) {
                NetworkStruct nonfjstruct = this.sn;
                this.sn = this.model.getStruct(false);
                for (int f = 0; f < this.sn.nodetype.size(); f++) {
                    if (this.sn.nodetype.get(f) != NodeType.Fork) {
                        continue;
                    }
                    // Declare variables outside switch to avoid scope issues
                    int joinIdx;

                    switch (options.config.fork_join) {
                        case "mmt":
                        case "default":
                        case "fjt":
                            Matrix TNfork = new Matrix(1, this.sn.nclasses);
                            for (int c = 0; c < this.sn.nchains; c++) {
                                Matrix inchain = Matrix.extractRows(this.sn.chains, c, c + 1, null).find();
                                Matrix chainVisits = this.sn.visits.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int r = (int) inchain.get(i);
                                    double visitSum = 0, TNsum = 0;
                                    for (int j = 0; j < inchain.length(); j++) {
                                        int k = (int) inchain.get(j);
                                        visitSum += chainVisits.get((int) this.sn.stationToStateful.get((int) this.sn.refstat.get(r)), k);
                                        TNsum += ret.TN.get((int) this.sn.refstat.get(r), k);
                                    }
                                    TNfork.set(r, (this.sn.nodevisits.get(c).get((int) parentForks.get(f), r) / visitSum) * TNsum);
                                }
                            }
                            // Find the join associated to the fork node f
                            joinIdx = -1;
                            for (int i = 0; i < this.sn.fj.getNumCols(); i++) {
                                if (this.sn.fj.get(f, i) != 0) {
                                    joinIdx = i;
                                    break;
                                }
                            }
                            ArrayList<Integer> forkauxclasses = new ArrayList<>();
                            for (int i = 0; i < fjforkmap.length(); i++) {
                                if (fjforkmap.get(i) == f) {
                                    forkauxclasses.add(i);
                                }
                            }
                            for (int s : forkauxclasses) {
                                int r = (int) fjclassmap.get(s);
                                if (joinIdx == -1) {
                                    forkLambda.set(s, (forkLambda.get(s) + TNfork.get(r)) / 2.0);
                                } else {
                                    double tnSum = 0;
                                    for (int i = 0; i < fjclassmap.length(); i++) {
                                        if (fjclassmap.get(i) == r) {
                                            tnSum += ret.TN.get((int) sn.nodeToStation.get(joinIdx), i);
                                        }
                                    }
                                    ret.TN.set((int) sn.nodeToStation.get(joinIdx), r, ret.TN.get((int) sn.nodeToStation.get(joinIdx), r) + tnSum - ret.TN.get((int) sn.nodeToStation.get(joinIdx), s));
                                    forkLambda.set(s, (forkLambda.get(s) + ret.TN.get((int) sn.nodeToStation.get(joinIdx), r)) / 2.0);
                                }
                                if (joinIdx == -1) {
                                    // No join nodes for this fork
                                    continue;
                                }
                                // Find the parallel paths coming out of the fork
                                ArrayList<Integer> toMerge = new ArrayList<>();
                                toMerge.add(r);
                                toMerge.add(s);
                                JobClass classr = nonfjmodel.getJobClassFromIndex(r);
//                                Matrix P = nonfjstruct.rtorig.get(classr).get(classr);
//                                Matrix ri = ModelAdapter.findPaths(this.sn, P, f, joinIdx, r,
//                                        toMerge, ret.QN, ret.TN, 0,
//                                        fjclassmap, fjforkmap, nonfjmodel);
                                int Psz = nonfjstruct.rtorig.get(classr).get(classr).getNumRows();
                                int Pcs_sz = nonfjmodel.getNumberOfClasses() * Psz;
                                Matrix Pcs = new Matrix(Pcs_sz, Pcs_sz);
                                for (int rcs = 0; rcs < nonfjstruct.nclasses; rcs++) {
                                    for (int scs = 0; scs < nonfjstruct.nclasses; scs++) {
                                        Matrix block_rs = nonfjstruct.rtorig.get(nonfjmodel.getClasses().get(rcs)).get(nonfjmodel.getClasses().get(scs));
                                        Pcs.setSliceEq(rcs * Psz, (rcs + 1) * Psz, scs * Psz, (scs + 1) * Psz, block_rs);
                                    }
                                }

                                Matrix ri = findPathsCS(this.sn, Pcs, f, joinIdx, r,
                                        toMerge, ret.QN, ret.TN, 0,
                                        fjclassmap, fjforkmap, nonfjmodel);
                                Matrix lambdai = Matrix.ones(ri.getNumRows(), ri.getNumCols()).elementDiv(ri);
                                double d0 = 0;
                                int parallel_branches = ri.getNonZeroLength();
                                for (int pow = 0; pow < parallel_branches; pow++) {
                                    Matrix nk = Maths.nCk(lambdai, pow + 1);
                                    nk = nk.sumRows();
                                    double currentSum = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
                                    d0 += FastMath.pow(-1, pow) * currentSum;
                                }
                                // Set the synchronisation delays
                                // Always set for auxiliary class s to ensure proper queue length accumulation
                                // Only set for original class r if this is the outer fork for that class
                                double syncDelay = d0 * ((ForkNodeParam) sn.nodeparam.get(this.model.getNodes().get(f))).fanOut - ri.elementSum() / ri.length();
                                ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(s), Exp.fitMean(syncDelay));
                                if (outerForks.get(f, r) != 0) {
                                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(r), Exp.fitMean(syncDelay));
                                }
                                nonfjmodel.refreshRates(null, null);
//                                if (forkIter==11) {
//                                    new SolverMVA(nonfjmodel).getAvgTable().print();
//                                    nonfjmodel.jsimwView();
//                                }
                            }
                            break;
                        case "heidelberger-trivedi":
                        case "ht":
                            // Find the join associated to the fork node f
                            joinIdx = -1;
                            for (int i = 0; i < this.sn.fj.getNumCols(); i++) {
                                if (this.sn.fj.get(f, i) != 0) {
                                    joinIdx = i;
                                    break;
                                }
                            }
                            for (int c = 0; c < this.sn.nchains; c++) {
                                Matrix inchain = Matrix.extractRows(this.sn.chains, c, c + 1, null).find();
                                for (int i = 0; i < inchain.length(); i++) {
                                    int r = (int) inchain.get(i);
                                    if (this.sn.nodevisits.get(c).get(f, r) == 0) {
                                        continue;
                                    }
                                    // Obtain the response times at the parallel branches
                                    int artificialClasses = 0;
                                    for (int j = 0; j < fjclassmap.length(); j++) {
                                        if (fjclassmap.get(j) == r) {
                                            artificialClasses++;
                                        }
                                    }
                                    Matrix ri = new Matrix(1, artificialClasses);
                                    int idx = 0;
                                    for (int j = 0; j < fjclassmap.length(); j++) {
                                        if (fjclassmap.get(j) == r) {
                                            // Compute the response time at all stations
                                            double colSum = 0;
                                            for (int k = 0; k < ret.RN.getNumRows(); k++) {
                                                if (k == nonfjstruct.nodeToStation.get(fj_auxiliary_delays.get(joinIdx)) ||
                                                        k == nonfjstruct.nodeToStation.get(joinIdx)) {
                                                    // Do not add the times at the auxiliary delay station or at the join delay
                                                    continue;
                                                }
                                                if (!Double.isNaN(ret.RN.get(k, j)) && !isInf(ret.RN.get(k, j))) {
                                                    colSum += ret.RN.get(k, j);
                                                }
                                            }
                                            ri.set(0, idx, colSum);
                                            idx++;
                                        }
                                    }
                                    Matrix lambdai = Matrix.ones(1, artificialClasses).elementDiv(ri);
                                    double d0 = 0;
                                    int parallel_branches = artificialClasses;
                                    for (int pow = 0; pow < parallel_branches; pow++) {
                                        Matrix nk = Maths.nCk(lambdai, pow + 1);
                                        nk = nk.sumRows();
                                        double currentSum = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
                                        d0 += FastMath.pow(-1, pow) * currentSum;
                                    }
                                    Matrix di = new Matrix(1, artificialClasses);
                                    for (int j = 0; j < artificialClasses; j++) {
                                        di.set(j, d0 * ((ForkNodeParam) sn.nodeparam.get(this.model.getNodes().get(f))).fanOut - ri.get(j));
                                    }
                                    double r0 = 0;
                                    for (int j = 0; j < inchain.length(); j++) {
                                        int k = (int) inchain.get(j);
                                        for (int l = 0; l < ret.RN.getNumRows(); l++) {
                                            if (l == nonfjstruct.nodeToStation.get(joinIdx)) {
                                                // Do not count the times at the join delay
                                                continue;
                                            }
                                            if (!Double.isNaN(ret.RN.get(l, k)) && !isInf(ret.RN.get(l, k))) {
                                                r0 += ret.RN.get(l, k);
                                            }
                                        }
                                    }
                                    // Update the delays at the join node and at the auxiliary delay
                                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(r), Exp.fitMean(d0 * ((ForkNodeParam) sn.nodeparam.get(this.model.getNodes().get(f))).fanOut));
                                    idx = 0;
                                    for (int s = 0; s < fjclassmap.length(); s++) {
                                        if (fjclassmap.get(s) == r) {
                                            ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(s), Exp.fitMean(di.get(idx)));
                                            idx++;
                                            ((Delay) nonfjmodel.getNodes().get(fj_auxiliary_delays.get(joinIdx))).setService(nonfjmodel.getJobClassFromIndex(s), Exp.fitMean(r0));
                                        }
                                    }
                                }
                            }
                    }
                }
                // Merge the results
                // Declare variables outside switch to avoid scope issues
                HashSet<Integer> originalClasses;
                Integer[] cls;
                Matrix TN_orig;

                switch (options.config.fork_join) {
                    case "heidelberger-trivedi":
                    case "ht":
                        nonfjmodel.refreshStruct(true);

                        // Save the throughputs of the original classes at the join node
                        ArrayList<Integer> joinIdx = new ArrayList<>();
                        for (int i = 0; i < this.sn.nodetype.size(); i++) {
                            if (this.sn.nodetype.get(i) == NodeType.Join) {
                                joinIdx.add(i);
                            }
                        }
                        originalClasses = new HashSet<>();
                        for (int i = 0; i < fjclassmap.length(); i++) {
                            if (fjclassmap.get(i) > -1) {
                                originalClasses.add((int) fjclassmap.get(i));
                            }
                        }
                        cls = originalClasses.toArray(new Integer[0]);
                        TN_orig = new Matrix(joinIdx.size(), cls.length);
                        for (int i = 0; i < joinIdx.size(); i++) {
                            for (int j = 0; j < cls.length; j++) {
                                TN_orig.set(i, j, ret.TN.get((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j]));
                            }
                        }

                        // Delete the queue lengths, response times, throughputs and utilizations of the original classes at the join nodes
                        for (int i = 0; i < joinIdx.size(); i++) {
                            for (int j = 0; j < cls.length; j++) {
                                ret.QN.set((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j], 0);
                                ret.RN.set((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j], 0);
                                ret.TN.set((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j], 0);
                                ret.UN.set((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j], 0);
                            }
                        }

                        // Remove the performance measures at the auxiliary delays
                        Collection<Integer> auxDelayValues = fj_auxiliary_delays.values();
                        HashSet<Integer> auxDelayStationIdx = new HashSet<>();
                        for (int i : auxDelayValues) {
                            auxDelayStationIdx.add((int) nonfjstruct.nodeToStation.get(i));
                        }
                        ret.QN.removeRows(auxDelayStationIdx);
                        ret.UN.removeRows(auxDelayStationIdx);
                        ret.RN.removeRows(auxDelayStationIdx);
                        ret.TN.removeRows(auxDelayStationIdx);

                        // Merge back artificial classes into their original classes
                        for (int r = 0; r < fjclassmap.length(); r++) {
                            int s = (int) fjclassmap.get(r);
                            if (s > -1) {
                                for (int i = 0; i < ret.QN.getNumRows(); i++) {
                                    ret.QN.set(i, s, ret.QN.get(i, s) + ret.QN.get(i, r));
                                    ret.UN.set(i, s, ret.UN.get(i, s) + ret.UN.get(i, r));
                                    // Add all throughputs of the auxiliary classes to facilitate the computation of the response times
                                    ret.TN.set(i, s, ret.TN.get(i, s) + ret.TN.get(i, r));
                                    ret.RN.set(i, s, ret.QN.get(i, s) / ret.TN.get(i, s));
                                }
                            }
                        }
                        // Set the throughputs back to their initial values for the original classes
                        for (int i = 0; i < joinIdx.size(); i++) {
                            for (int j = 0; j < cls.length; j++) {
                                ret.TN.set((int) nonfjstruct.nodeToStation.get(joinIdx.get(i)), cls[j], TN_orig.get(i, j));
                            }
                        }
                        break;
                    case "mmt":
                    case "default":
                    case "fjt":
                        // Save the throughputs of the original classes at the join node and at the source node
                        ArrayList<Integer> joinSourceIdx = new ArrayList<>();
                        for (int i = 0; i < this.sn.nodetype.size(); i++) {
                            if (this.sn.nodetype.get(i) == NodeType.Join || this.sn.nodetype.get(i) == NodeType.Source) {
                                joinSourceIdx.add(i);
                            }
                        }
                        originalClasses = new HashSet<>();
                        for (int i = 0; i < fjclassmap.length(); i++) {
                            if (fjclassmap.get(i) > -1) {
                                originalClasses.add((int) fjclassmap.get(i));
                            }
                        }
                        cls = originalClasses.toArray(new Integer[0]);
                        TN_orig = new Matrix(joinSourceIdx.size(), cls.length);
                        for (int i = 0; i < joinSourceIdx.size(); i++) {
                            for (int j = 0; j < cls.length; j++) {
                                TN_orig.set(i, j, ret.TN.get((int) nonfjstruct.nodeToStation.get(joinSourceIdx.get(i)), cls[j]));
                            }
                        }

                        // Merge back artificial classes into their original classes
                        for (int r = 0; r < fjclassmap.length(); r++) {
                            int s = (int) fjclassmap.get(r);
                            if (s > -1) {
                                for (int i = 0; i < ret.QN.getNumRows(); i++) {
                                    if (r < ret.QN.getNumCols()) {
                                        ret.QN.set(i, s, ret.QN.get(i, s) + ret.QN.get(i, r));
                                        ret.UN.set(i, s, ret.UN.get(i, s) + ret.UN.get(i, r));
                                        // Add all throughputs of the auxiliary classes to facilitate the computation of the response times
                                        ret.TN.set(i, s, ret.TN.get(i, s) + ret.TN.get(i, r));
                                        ret.RN.set(i, s, ret.QN.get(i, s) / ret.TN.get(i, s));
                                    }
                                }
                            }
                        }
                        // Set the throughputs back to their initial values for the original classes
                        for (int i = 0; i < joinSourceIdx.size(); i++) {
                            for (int j = 0; j < cls.length; j++) {
                                ret.TN.set((int) nonfjstruct.nodeToStation.get(joinSourceIdx.get(i)), cls[j], TN_orig.get(i, j));
                            }
                        }
                        break;
                }
                HashSet<Integer> artificialClasses = new HashSet<>();
                for (int i = 0; i < fjclassmap.length(); i++) {
                    if (fjclassmap.get(i) > -1) {
                        artificialClasses.add(i);
                    }
                }
                ret.QN.removeCols(artificialClasses);
                ret.UN.removeCols(artificialClasses);
                ret.RN.removeCols(artificialClasses);
                ret.TN.removeCols(artificialClasses);
                ret.CN.removeCols(artificialClasses);
                ret.XN.removeCols(artificialClasses);
            }
            QN = ret.QN;
            iter += ret.iter;
        }
        this.sn = this.model.getStruct(true);

        // Compute average arrival rate at steady-state
        AvgHandle T = avgHandles.getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(sn, ret.TN, T);

        AvgHandle W = avgHandles.getAvgResidTHandles();
        Matrix WN;
        if (ret.RN != null && !ret.RN.isEmpty()) {
            WN = snGetResidTFromRespT(sn, ret.RN, W);
        } else {
            // Handle case where RN is null or empty by creating a zero matrix
            WN = new Matrix(sn.nstations, sn.nclasses);
        }

        this.res = new MVAResult();
        this.res.method = ret.method;
        this.res.QN = ret.QN; //
        this.res.RN = ret.RN; //
        this.res.XN = ret.XN;
        this.res.UN = ret.UN; //
        this.res.TN = ret.TN; //
        this.res.CN = ret.CN;
        this.res.AN = AN;
        this.res.WN = WN;
        this.res.runtime = (System.nanoTime() - T0) / 1000000000.0;
        this.res.iter = iter;
        this.res.logNormConstAggr = ret.logNormConstAggr;

        return this.res;
    }

    /**
     * Check if the node types match exactly the expected types (order independent)
     */
    private boolean checkNodeTypes(List<NodeType> nodeTypes, NodeType... expectedTypes) {
        if (nodeTypes.size() != expectedTypes.length) {
            return false;
        }
        
        // Count occurrences of each node type
        Map<NodeType, Integer> actualCounts = new HashMap<>();
        Map<NodeType, Integer> expectedCounts = new HashMap<>();
        
        for (NodeType type : nodeTypes) {
            actualCounts.put(type, actualCounts.getOrDefault(type, 0) + 1);
        }
        
        for (NodeType type : expectedTypes) {
            expectedCounts.put(type, expectedCounts.getOrDefault(type, 0) + 1);
        }
        
        return actualCounts.equals(expectedCounts);
    }

}