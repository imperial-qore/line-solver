/**
 * @file Solver-agnostic driver of the fork-join fixed point
 *
 * @since LINE 3.0
 */
package jline.solvers.fj;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Source;
import jline.io.Ret;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import jline.api.fj.FJ_ordstat_exp;
import jline.api.sn.SnJoinQuorum;
import jline.lang.ModelAdapter;
import jline.lang.nodes.Delay;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.nodes.Node;
import jline.lang.processes.Exp;

import static jline.io.InputOutput.line_debug;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;
import static jline.util.Utils.isInf;
import static jline.lang.ModelAdapter.ht;
import static jline.lang.ModelAdapter.mmt;
import static jline.lang.ModelAdapter.sort_forks;
import static jline.lang.ModelAdapter.findPathsCS;

/**
 * Solver-agnostic driver of the fork-join fixed point.
 *
 * <p>The transformation (ModelAdapter.mmt, or ht under
 * options.config.fork_join) turns the model into a plain network in which every
 * fork is a router, every join a zero-service delay, and the parallelism is
 * carried by auxiliary open classes of arrival rate (fanout-1)*forkLambda. Each
 * pass solves that network, recomputes the synchronisation delays from the
 * resulting metrics and updates forkLambda; the loop ends when the queue
 * lengths stop moving.
 *
 * <p>{@link InnerSolve} is the inner solve. Nothing in the loop reads a solver
 * internal, so any solver whose analyzer honours that contract can be driven
 * here: SolverMVA passes its analyzer dispatch, SolverNC the
 * normalizing-constant analyzers. Port of
 * matlab/src/solvers/@NetworkSolver/fjFixedPoint.m.
 *
 * <p>The carrier is {@link MVAResult} for historical reasons: it is the shape
 * the loop already consumed when it lived inside MVARunner (QN, UN, RN, TN, CN,
 * XN, iter, logNormConstAggr, method), and reusing it keeps the loop body
 * identical to the one that produced today's numbers.
 */
public final class FJFixedPoint {
    private FJFixedPoint() {}

    /** The inner solve of one pass of the fixed point. */
    public interface InnerSolve {
        /**
         * Solves the transformed network.
         *
         * <p>The transformed MODEL is passed alongside its struct because a
         * state-based inner solve needs it: SolverFluid's initSol keys sn.state
         * by the node objects of the model it was built from, so handing it the
         * original model with the transformed struct would look up stations that
         * are not in it. A struct-only solver (MVA, NC) can ignore the argument.
         *
         * @param net the transformed network, or the original one when it has no fork
         * @param sn the struct of the transformed network
         * @param options the solver options
         * @return the metrics of that solve
         */
        MVAResult solve(Network net, NetworkStruct sn, SolverOptions options);
    }

    /** State the fixed point retains across calls, as MATLAB keeps on the solver. */
    public static class FJState {
        /** Cached MMT transformation, reused while the struct version matches. */
        public Ret.FJApprox mmtCache;
        /** Auxiliary-class arrival rates, the warm start of the next call. */
        public Matrix fjForkLambda;

        public FJState() {
        }

        public FJState(Ret.FJApprox mmtCache, Matrix fjForkLambda) {
            this.mmtCache = mmtCache;
            this.fjForkLambda = fjForkLambda;
        }
    }

    /** Outcome of the fixed point: the merged metrics and the retained state. */
    public static class FJOutcome {
        public MVAResult ret;
        public Matrix QN;
        public int iter;
        public NetworkStruct sn;
        public FJState state;
    }

    /**
     * Runs the fork-join fixed point, or the inner solve once when the model has
     * no fork.
     *
     * @param model the model being solved
     * @param sn its struct, already compiled
     * @param options the solver options
     * @param state the retained transformation cache and warm start
     * @param fcn the inner solve
     * @param T0 the wall-clock marker of the enclosing analyzer
     * @return the merged metrics, the working struct and the updated state
     */
    public static FJOutcome run(Network model, NetworkStruct sn, SolverOptions options,
                                FJState state, InnerSolve fcn, long T0) {
        int iter = 0;
        boolean quorumClamped = false;
        if (state == null) {
            state = new FJState();
        }
        boolean forkLoop = true;
        int forkIter = 0;
        int forkNodes = 0;
        for (NodeType n : sn.nodetype) {
            if (n == NodeType.Fork)
                forkNodes++;
        }
        Matrix forkLambda = new Matrix(1, 2 * sn.nclasses * forkNodes);
        forkLambda.fill(GlobalConstants.FineTol);
        // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint warm-start reuse gate)
        if (options.config.fj_warmstart && state.fjForkLambda != null
                && state.fjForkLambda.getNumRows() == forkLambda.getNumRows()
                && state.fjForkLambda.getNumCols() == forkLambda.getNumCols()) {
            forkLambda = state.fjForkLambda.copy();
        }
        Matrix QN = new Matrix(1, sn.nclasses);
        QN.fill(GlobalConstants.Immediate);
        Matrix QN_1 = new Matrix(1, sn.nclasses);
        Matrix UN = new Matrix(1, sn.nclasses);
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
        while (forkLoop && forkIter < options.iter_max && !Solver.timeExceeded(T0, options.timeout)) {
            // Clear out any fields from the past returns
            ret = new MVAResult();
            if (model.hasFork()) {
                forkIter += 1;
                if (forkIter == 1) {
                    Ret.FJApprox approxReturn = null;
                    switch (options.config.fork_join) {
                        case "heidelberger-trivedi":
                        case "ht":
                            approxReturn = ht(model);
                            break;
                        case "mmt":
                        case "default":
                        case "fjt":
                            // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint warm-start reuse gate)
                            boolean cacheUsable = state.mmtCache != null
                                    && state.mmtCache.baseSn != null
                                    && state.mmtCache.baseSn == model.getStruct(false)
                                    && ModelAdapter.refreshServicesFromBase(state.mmtCache);
                            if (cacheUsable) {
                                approxReturn = state.mmtCache;
                                outerForks = state.mmtCache.outerForks;
                                parentForks = state.mmtCache.parentForks;
                                approxReturn.nonfjmodel.refreshRates(null, null);
                            } else {
                                approxReturn = mmt(model, forkLambda);
                                Ret.FJsortForks sortForksReturn = sort_forks(sn, approxReturn.nonfjmodel.getStruct(false),
                                        approxReturn.fjforkmap, approxReturn.fjclassmap, approxReturn.nonfjmodel);
                                outerForks = sortForksReturn.outerForks;
                                parentForks = sortForksReturn.parentForks;
                                approxReturn.baseSn = model.getStruct(false);
                                approxReturn.outerForks = outerForks;
                                approxReturn.parentForks = parentForks;
                                state.mmtCache = approxReturn;
                            }
                            break;
                    }
                    nonfjmodel = approxReturn.nonfjmodel;
                    fjclassmap = approxReturn.fjclassmap;
                    fjforkmap = approxReturn.fjforkmap;
                    fj_auxiliary_delays = approxReturn.fj_auxiliary_delays;
                    fanout = approxReturn.fanout;
                } else if (!options.config.fork_join.equals("heidelberger-trivedi") && !options.config.fork_join.equals("ht")) {
                    Source nonfjSource = nonfjmodel.getSource();
                    for (int r = 0; r < fjclassmap.length(); r++) {
                        int s = (int) fjclassmap.get(r);
                        if (s > -1) {
                            if (fanout.get(r) > 0) {
                                if (!nonfjSource.getArrivalProcess(nonfjmodel.getClasses().get(r)).isDisabled()) {
                                    // REBIND rather than Exp.updateRate. The rate
                                    // reaches sn.rates either way, which is all an
                                    // MVA or NC inner solve reads, but the fluid
                                    // drift reads the PHASE representation, and
                                    // mutating the distribution in place left that
                                    // at the initial GlobalConstants.FineTol rate:
                                    // the auxiliary classes then carried no traffic
                                    // and the fixed point converged on the untouched
                                    // transformed model.
                                    nonfjSource.setArrival(nonfjmodel.getClasses().get(r),
                                            new Exp((fanout.get(r) - 1) * forkLambda.get(r)));
                                }
                            }
                        }
                    }
                    nonfjmodel.refreshRates(null, null);
                    // refreshRates writes sn.rates and sn.scv only. An MVA or NC
                    // inner solve reads the rate, but the fluid drift reads
                    // sn.mu/sn.phi/sn.proc, and updateRate on an Exp leaves the
                    // SCV at 1, so refreshProcesses would skip the phase refresh
                    // and the auxiliary source would keep integrating at its
                    // initial GlobalConstants.FineTol rate on every pass.
                    nonfjmodel.refreshProcessPhases(null, null);
                    nonfjmodel.refreshProcessRepresentations();
                }
                sn = nonfjmodel.getStruct(false);
                if (jline.io.LineConsole.ownsLog()) {
                    if (QN_1.getNumRows() == QN.getNumRows() && QN_1.getNumCols() == QN.getNumCols()) {
                        jline.io.LineConsole.step(
                                "fork-join iteration %d: queue lengths moved by at most %.3e",
                                forkIter, QN_1.sub(QN).elementMaxAbs());
                    } else {
                        jline.io.LineConsole.step(
                                "fork-join iteration %d: transformed model rebuilt", forkIter);
                    }
                }
                if (forkIter > 2 && QN_1.getNumRows() == QN.getNumRows() && QN_1.getNumCols() == QN.getNumCols()) {
                    // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint mixed absolute/relative convergence test)
                    boolean converged = true;
                    for (int i = 0; i < QN.getNumRows() && converged; i++) {
                        for (int j = 0; j < QN.getNumCols(); j++) {
                            double qn = QN.get(i, j);
                            double qn1 = QN_1.get(i, j);
                            double delta = FastMath.abs(qn1 - qn);
                            // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint iter_tol vs CoarseTol)
                            double thresh = GlobalConstants.Zero + options.iter_tol * FastMath.abs(qn);
                            if (Double.isNaN(delta) || delta > thresh) {
                                converged = false;
                                break;
                            }
                        }
                    }
                    if (converged) {
                        forkLoop = false;
                    } else {
                        QN_1 = QN;
                    }
                } else {
                    if (model.hasOpenClasses()) {
                        int sourceIndex = model.getIndexSourceNode();
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
                            if (Uiopen > 0.99 && sn.nservers.get(i) != Integer.MAX_VALUE) {
                                System.out.println("The model may be unstable: the utilization of station " + i + " for open classes exceeds 99 percent.\n");
                            }
                        }
                    }
                    QN_1 = QN;
                }
            } else {
                forkLoop = false;
            }

            ret = fcn.solve(nonfjmodel != null ? nonfjmodel : model, sn, options);
            if (model.hasFork()) {
                NetworkStruct nonfjstruct = sn;
                sn = model.getStruct(false);
                // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint Pcs precompute and sync delay setting)
                Matrix Pcs = null;
                if (options.config.fork_join.equals("mmt") || options.config.fork_join.equals("default") || options.config.fork_join.equals("fjt")) {
                    Map<JobClass, Map<JobClass, Matrix>> nonfjRtorig = nonfjmodel.getLinkedRoutingMatrix();
                    JobClass anyClass = nonfjmodel.getClasses().get(0);
                    int Psz = nonfjRtorig.get(anyClass).get(anyClass).getNumRows();
                    int nclassesNonfj = nonfjmodel.getNumberOfClasses();
                    int Pcs_sz = nclassesNonfj * Psz;
                    Pcs = new Matrix(Pcs_sz, Pcs_sz);
                    for (int rcs = 0; rcs < nclassesNonfj; rcs++) {
                        for (int scs = 0; scs < nclassesNonfj; scs++) {
                            Matrix block_rs = nonfjRtorig.get(nonfjmodel.getClasses().get(rcs)).get(nonfjmodel.getClasses().get(scs));
                            Pcs.setSliceEq(rcs * Psz, (rcs + 1) * Psz, scs * Psz, (scs + 1) * Psz, block_rs);
                        }
                    }
                }
                for (int f = 0; f < sn.nodetype.size(); f++) {
                    if (sn.nodetype.get(f) != NodeType.Fork) {
                        continue;
                    }
                    // Declare variables outside switch to avoid scope issues
                    int joinIdx;

                    switch (options.config.fork_join) {
                        case "mmt":
                        case "default":
                        case "fjt":
                            Matrix TNfork = new Matrix(1, sn.nclasses);
                            for (int c = 0; c < sn.nchains; c++) {
                                Matrix inchain = Matrix.extractRows(sn.chains, c, c + 1, null).find();
                                Matrix chainVisits = sn.visits.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int r = (int) inchain.get(i);
                                    double visitSum = 0, TNsum = 0;
                                    for (int j = 0; j < inchain.length(); j++) {
                                        int k = (int) inchain.get(j);
                                        visitSum += chainVisits.get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), k);
                                        TNsum += ret.TN.get((int) sn.refstat.get(r), k);
                                    }
                                    TNfork.set(r, (sn.nodevisits.get(c).get((int) parentForks.get(f), r) / visitSum) * TNsum);
                                }
                            }
                            // Find the join associated to the fork node f
                            joinIdx = -1;
                            for (int i = 0; i < sn.fj.getNumCols(); i++) {
                                if (sn.fj.get(f, i) != 0) {
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
                                    int joinStat = (int) sn.nodeToStation.get(joinIdx);
                                    double tnJoinR = ret.TN.get(joinStat, r);
                                    double tnJoinS = ret.TN.get(joinStat, s);
                                    double tnSum = 0;
                                    for (int i = 0; i < fjclassmap.length(); i++) {
                                        if (fjclassmap.get(i) == r) {
                                            tnSum += ret.TN.get(joinStat, i);
                                        }
                                    }
                                    ret.TN.set(joinStat, r, ret.TN.get(joinStat, r) + tnSum - ret.TN.get(joinStat, s));
                                    forkLambda.set(s, (forkLambda.get(s) + ret.TN.get(joinStat, r)) / 2.0);
                                }
                                if (joinIdx == -1 || outerForks.get(f, r) == 0) {
                                    // No join nodes for this fork, or not the outer fork for this class
                                    continue;
                                }
                                // Find the parallel paths coming out of the fork
                                ArrayList<Integer> toMerge = new ArrayList<>();
                                toMerge.add(r);
                                toMerge.add(s);
                                Matrix ri = findPathsCS(sn, Pcs, f, joinIdx, r,
                                        toMerge, ret.QN, ret.TN, 0,
                                        fjclassmap, fjforkmap, nonfjmodel);
                                // tasksPerLink = w sends w IDENTICAL tasks down each link,
                                // so the join synchronises on w*B siblings and not on B:
                                // the sibling set is each branch's completion time
                                // REPLICATED w times, and the order statistic is taken
                                // over that multiset. Scaling E[X_(k)] by w instead (what
                                // this did before) is w*H_B/mu where the answer is
                                // H_(w*B)/mu, which OVER-states the delay by more the
                                // larger w is. w = 1 replicates to itself, so nothing
                                // moves there.
                                int w = forkTasksPerLink(sn, model.getNodes().get(f));
                                ri = replicate(ri, w);
                                // The join fires on the k-th sibling completion, k = the
                                // sibling count on a standard join and the declared quorum
                                // on a PARTIAL one. The quorum is declared against the
                                // SIBLING count w*B, which is the replicated length.
                                int kreq = SnJoinQuorum.snJoinQuorum(sn, model.getNodes().get(joinIdx),
                                        model.getClasses().get(r), ri.length());
                                double d0 = FJ_ordstat_exp.fj_ordstat_exp(ri, kreq);
                                // see _kb/06-solver-catalog.md (JAR-only implementation notes: FJFixedPoint Pcs precompute and sync delay setting)
                                double syncDelay = d0 - ri.elementSum() / ri.length();
                                if (syncDelay < 0) {
                                    // The quorum is met BEFORE the branch the transform's own
                                    // token walks, so the parent ought to leave ahead of it. The
                                    // MMT cannot express that: its token is a job of the closed
                                    // chain and must finish its branch, and that closed token is
                                    // what keeps the branch stable, so it cannot be made open
                                    // either. The delay floors at zero, which OVER-states the
                                    // cycle time.
                                    quorumClamped = true;
                                    syncDelay = 0;
                                }
                                ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(s), Exp.fitMean(syncDelay));
                                if (outerForks.get(f, r) != 0) {
                                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(r), Exp.fitMean(syncDelay));
                                }
                            }
                            break;
                        case "heidelberger-trivedi":
                        case "ht":
                            // Find the join associated to the fork node f
                            joinIdx = -1;
                            for (int i = 0; i < sn.fj.getNumCols(); i++) {
                                if (sn.fj.get(f, i) != 0) {
                                    joinIdx = i;
                                    break;
                                }
                            }
                            for (int c = 0; c < sn.nchains; c++) {
                                Matrix inchain = Matrix.extractRows(sn.chains, c, c + 1, null).find();
                                for (int i = 0; i < inchain.length(); i++) {
                                    int r = (int) inchain.get(i);
                                    if (sn.nodevisits.get(c).get(f, r) == 0) {
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
                                    // The join fires on the k-th branch completion, k = the
                                    // branch count on a standard join and the declared quorum
                                    // on a PARTIAL one.
                                    int kreq = SnJoinQuorum.snJoinQuorum(sn, model.getNodes().get(joinIdx),
                                            model.getClasses().get(r), artificialClasses);
                                    double d0 = FJ_ordstat_exp.fj_ordstat_exp(ri, kreq);
                                    Matrix di = new Matrix(1, artificialClasses);
                                    for (int j = 0; j < artificialClasses; j++) {
                                        // Under a quorum the k-th completion can precede a
                                        // branch's own, and then that branch waits no further.
                                        // see the mmt branch for why the floor is a boundary of
                                        // the transform and not a choice.
                                        // No fanOut factor here: ModelAdapter.ht REFUSES
                                        // tasksPerLink > 1 by name, so w is 1 on every model
                                        // that reaches this branch.
                                        double dij = d0 - ri.get(j);
                                        if (dij < 0) {
                                            quorumClamped = true;
                                            dij = 0;
                                        }
                                        di.set(j, dij);
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
                                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getClasses().get(r), Exp.fitMean(d0));
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
                // Batch refreshRates after all sync delay updates (moved out of inner loops for performance)
                if (!options.config.fork_join.equals("heidelberger-trivedi") && !options.config.fork_join.equals("ht")) {
                    nonfjmodel.refreshRates(null, null);
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
                        for (int i = 0; i < sn.nodetype.size(); i++) {
                            if (sn.nodetype.get(i) == NodeType.Join) {
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
                        for (int i = 0; i < sn.nodetype.size(); i++) {
                            if (sn.nodetype.get(i) == NodeType.Join || sn.nodetype.get(i) == NodeType.Source) {
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
        // The fork-join loop previously exhausted options.iter_max silently, so a
        // non-converged MMT fixed point was returned as a normal result.
        if (model.hasFork()) {
            if (quorumClamped) {
                line_warning(mfilename(new Object() {
                }), "A quorum join fires before the branch the fork-join transformation follows, which it cannot represent: the synchronisation delay is floored at zero, which OVER-states the cycle time and so under-states the throughput. SolverLDES and SolverJMT simulate the quorum on their sample path; SolverCTMC and SolverSSA refuse it, because a quorum fork-join has an unbounded state space.");
            }
            if (forkLoop && forkIter >= options.iter_max) {
                line_warning(mfilename(new Object() {
                }), "The fork-join (mmt) fixed point did not converge in options.iter_max=" + options.iter_max + " iterations; returning the interim solution.");
            }
            // Retain the MMT iterate so that a subsequent runAnalyzer call on this
            // solver (an outer LN iteration) resumes the fixed point here.
            state.fjForkLambda = forkLambda;
        }

        FJOutcome out = new FJOutcome();
        out.ret = ret;
        out.QN = QN;
        out.iter = iter;
        out.sn = sn;
        out.state = state;
        return out;
    }

    /**
     * The tasksPerLink of a Fork node, at least 1. A Fork with no nodeparam entry, or
     * one whose fanOut was never set, emits one task per link.
     */
    private static int forkTasksPerLink(NetworkStruct sn, Node forkNode) {
        if (sn == null || sn.nodeparam == null || forkNode == null) {
            return 1;
        }
        NodeParam param = sn.nodeparam.get(forkNode);
        if (!(param instanceof ForkNodeParam)) {
            return 1;
        }
        double w = ((ForkNodeParam) param).fanOut;
        if (Double.isNaN(w) || w < 1) {
            return 1;
        }
        return (int) Math.max(1, Math.round(w));
    }

    /**
     * The row vector {@code v} concatenated with itself {@code w} times, which is the
     * sibling multiset of a fork emitting w identical tasks per link. w = 1 returns a
     * copy of the input, so nothing moves on an ordinary fork.
     */
    private static Matrix replicate(Matrix v, int w) {
        if (w <= 1 || v == null || v.length() == 0) {
            return v;
        }
        Matrix out = new Matrix(1, v.length() * w);
        for (int k = 0; k < w; k++) {
            for (int i = 0; i < v.length(); i++) {
                out.set(0, k * v.length() + i, v.get(i));
            }
        }
        return out;
    }
}
