package jline.lang;

import jline.io.Ret;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.processes.Det;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Immediate;
import jline.lang.sections.Forker;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.io.*;
import java.util.*;

import static jline.io.InputOutput.*;
import static jline.api.sn.SnGetDemandsChain.snGetDemandsChain;

import jline.api.fj.FJ_ordstat_exp;
import jline.api.sn.SnJoinQuorum;
import jline.api.fes.FESAggregator;
import jline.api.fes.FESResult;
import jline.api.fes.FESOptions;

/**
 * Static class to transform and adapt models, providing functionality for:
 * - Creating tagged job models for response time analysis
 * - Fork-join network transformations (formerly from FJ.java)
 * - Model preprocessing and adaptation operations
 */
public class ModelAdapter {
    /**
     * Routing block of a class pair, materialised if absent. MATLAB keeps the
     * routing as a dense cell array in which an unused class pair is an
     * all-zero matrix, whereas the Java map omits it, so a model whose classes
     * are only sparsely linked (an LQN2QN conversion, for instance) would
     * dereference null here. The block is sized from the ones already present,
     * which give the node count of the model this routing belongs to.
     */
    private static Matrix rtBlock(Map<JobClass, Map<JobClass, Matrix>> P, JobClass from, JobClass to) {
        Map<JobClass, Matrix> row = P.get(from);
        if (row == null) {
            row = new HashMap<>();
            P.put(from, row);
        }
        Matrix block = row.get(to);
        if (block != null) {
            return block;
        }
        int n = 0;
        outer:
        for (Map<JobClass, Matrix> other : P.values()) {
            for (Matrix m : other.values()) {
                if (m != null) {
                    n = m.getNumRows();
                    break outer;
                }
            }
        }
        block = new Matrix(n, n);
        row.put(to, block);
        return block;
    }

    
    /**
     * Result of tagging a chain in a model
     */
    public static class TaggedChainResult {
        private final Network taggedModel;
        private final JobClass taggedJob;
        
        public TaggedChainResult(Network taggedModel, JobClass taggedJob) {
            this.taggedModel = taggedModel;
            this.taggedJob = taggedJob;
        }
        
        public Network getTaggedModel() {
            return taggedModel;
        }
        
        public JobClass getTaggedJob() {
            return taggedJob;
        }
    }
    
    /**
     * Create a tagged version of a job chain for response time analysis
     * 
     * @param model The original model
     * @param chain The chain to tag
     * @param jobclass The specific job class to tag (optional, defaults to first class in chain)
     * @param suffix The suffix to add to tagged class names (optional, defaults to ".tagged")
     * @return TaggedChainResult containing the tagged model and tagged job class
     */
    /**
     * Build a tagged copy of {@code model}: one job of {@code jobclass} is moved
     * out of its own class into a new class of population 1, so that a solver
     * can follow that single job.
     *
     * <p>Port of matlab/src/io/@ModelAdapter/tagChain.m, which is the reference.
     *
     * <p>THIS REPLACED A STUB. The previous body created the tagged class and
     * stopped: it gave the class NO SERVICE at any station and NO ROUTING, so
     * the tagged model had a job that could neither be served nor move, and
     * SolverCTMC produced an EMPTY (0x0) generator from it. It also called
     * {@code addJobClass} on top of the ClosedClass constructor, which already
     * registers the class, and its population decrement was guarded by
     * {@code > 1} so a two-job class kept both jobs and the tagged model carried
     * one more job than the original. Everything downstream degraded silently
     * because {@code getCdfRespT} wrapped the lot in a catch that returned
     * zeros.
     *
     * <p>The three things that make a tagged class real, and that the stub
     * omitted: a service process at EVERY station cloned from the source class,
     * the source class's routing replicated for the new class over the linked
     * routing matrix, and one job actually MOVED rather than added.
     */
    public static TaggedChainResult tagChain(Network model, Chain chain, JobClass jobclass, String suffix) {
        if (suffix == null || suffix.isEmpty()) {
            suffix = ".tagged";
        }
        if (jobclass == null && !chain.getClasses().isEmpty()) {
            jobclass = chain.getClasses().get(0);
        }
        if (jobclass == null) {
            throw new RuntimeException("tagChain: the chain carries no class to tag");
        }

        Network taggedModel = model.copy();

        // Resolve the chain's classes inside the COPY: the Chain object holds
        // references into the original model, which the copy does not share.
        List<JobClass> chainInTagged = new ArrayList<JobClass>();
        for (JobClass src : chain.getClasses()) {
            for (JobClass cls : taggedModel.getClasses()) {
                if (cls.getName().equals(src.getName())) {
                    chainInTagged.add(cls);
                    break;
                }
            }
        }
        JobClass sourceInTagged = null;
        for (JobClass cls : taggedModel.getClasses()) {
            if (cls.getName().equals(jobclass.getName())) {
                sourceInTagged = cls;
                break;
            }
        }
        if (sourceInTagged == null || chainInTagged.isEmpty()) {
            throw new RuntimeException("tagChain: the class to tag is not present in the model copy");
        }

        // The linked routing matrix must be read BEFORE the new classes exist,
        // so that it is the original class-pair map and not a half-extended one.
        Map<JobClass, Map<JobClass, Matrix>> P = taggedModel.getLinkedRoutingMatrix();
        if (P == null) {
            throw new RuntimeException("tagChain: the model has no linked routing matrix; "
                    + "link() it before asking for a tagged copy");
        }

        int nnodes = taggedModel.getNumberOfNodes();
        List<JobClass> taggedClasses = new ArrayList<JobClass>();
        for (JobClass src : chainInTagged) {
            double pop = (src == sourceInTagged) ? 1.0 : 0.0;
            Station refstat = (src instanceof ClosedClass)
                    ? ((ClosedClass) src).getReferenceStation() : null;
            if (refstat == null && !taggedModel.getStations().isEmpty()) {
                refstat = taggedModel.getStations().get(0);
            }
            // The constructor registers the class with the model; calling
            // addJobClass on top of it adds the same class twice.
            ClosedClass tagged = new ClosedClass(taggedModel, src.getName() + suffix, pop,
                    refstat, src.getPriority());

            // A class with no service process is not served anywhere, and the
            // state-space generator then has no transition to build.
            for (Station st : taggedModel.getStations()) {
                if (st instanceof ServiceStation) {
                    ServiceStation ss = (ServiceStation) st;
                    Distribution d = ss.getServiceProcess(src);
                    if (d != null) {
                        ss.setService(tagged, d);
                    }
                }
            }
            taggedClasses.add(tagged);
        }

        // MOVE the job: the tagged class gained one, so the source loses one.
        if (sourceInTagged instanceof ClosedClass) {
            ClosedClass src = (ClosedClass) sourceInTagged;
            src.setPopulation(src.getPopulation() - 1);
        }

        // Replicate the chain's routing for the tagged classes, pair by pair.
        for (int ir = 0; ir < chainInTagged.size(); ir++) {
            JobClass fromOld = chainInTagged.get(ir);
            JobClass fromNew = taggedClasses.get(ir);
            Map<JobClass, Matrix> rowOld = P.get(fromOld);
            Map<JobClass, Matrix> rowNew = P.get(fromNew);
            if (rowNew == null) {
                rowNew = new HashMap<JobClass, Matrix>();
                P.put(fromNew, rowNew);
            }
            for (int is = 0; is < chainInTagged.size(); is++) {
                JobClass toOld = chainInTagged.get(is);
                JobClass toNew = taggedClasses.get(is);
                Matrix block = (rowOld == null) ? null : rowOld.get(toOld);
                rowNew.put(toNew, block == null ? new Matrix(nnodes, nnodes) : block.copy());
            }
        }
        // Every remaining class pair needs a block, or the relink sees a hole.
        for (JobClass r : taggedModel.getClasses()) {
            Map<JobClass, Matrix> row = P.get(r);
            if (row == null) {
                row = new HashMap<JobClass, Matrix>();
                P.put(r, row);
            }
            for (JobClass sCls : taggedModel.getClasses()) {
                if (!row.containsKey(sCls) || row.get(sCls) == null) {
                    row.put(sCls, new Matrix(nnodes, nnodes));
                }
            }
        }

        taggedModel.relinkFromRtorig(P);
        taggedModel.refreshStruct(true);

        return new TaggedChainResult(taggedModel,
                taggedClasses.isEmpty() ? jobclass : taggedClasses.get(taggedClasses.size() - 1));
    }

    /**
     * Convenience method with default parameters
     */
    public static TaggedChainResult tagChain(Network model, Chain chain, JobClass jobclass) {
        return tagChain(model, chain, jobclass, ".tagged");
    }
    
    /**
     * Convenience method with default parameters
     */
    public static TaggedChainResult tagChain(Network model, Chain chain) {
        return tagChain(model, chain, null, ".tagged");
    }

    // ========== Fork-Join Methods (formerly from FJ.java) ==========

    /**
     * Finds the response times along each path leading out of startNode up to (and not including) endNode
     */
    public static Matrix findPaths(NetworkStruct sn, Matrix P, int startNode, int endNode, int r, ArrayList<Integer> toMerge,
                                   Matrix QN, Matrix TN, double currentTime, Matrix fjclassmap, Matrix fjforkmap,
                                   Network nonfjmodel) {
        return findPaths(sn, P, startNode, endNode, r, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap,
                nonfjmodel, new java.util.HashSet<Integer>());
    }

    /** Simple paths only; see the ONPATH note on findPathsCS. */
    public static Matrix findPaths(NetworkStruct sn, Matrix P, int startNode, int endNode, int r, ArrayList<Integer> toMerge,
                                   Matrix QN, Matrix TN, double currentTime, Matrix fjclassmap, Matrix fjforkmap,
                                   Network nonfjmodel, java.util.Set<Integer> onPath) {
        if (startNode == endNode) {
            double qLen = 0;
            double tput = 0;
            for (int s : toMerge) {
                qLen += QN.get((int) sn.nodeToStation.get(startNode), s);
                tput += TN.get((int) sn.nodeToStation.get(startNode), s);
            }
            Matrix ri = new Matrix(1, 1);
            ri.set(0, 0, currentTime - qLen / tput);
            return ri;
        }
        Matrix ri = new Matrix(1, 0);
        onPath.add(startNode);
        for (int i = 0; i < P.getNumCols(); i++) {
            if (P.get(startNode, i) == 0) {
                continue;
            }
            if (i != endNode && onPath.contains(i)) {
                continue;
            }
            double qLen = 0;
            double tput = 1;
            if (sn.nodeToStation.get(i) > -1) {
                tput = 0;
                for (int s : toMerge) {
                    qLen += QN.get((int) sn.nodeToStation.get(i), s);
                    tput += TN.get((int) sn.nodeToStation.get(i), s);
                }
            }
            if (sn.nodetype.get(i) == NodeType.Fork) {
                int joinIdx = 0;
                while (joinIdx < sn.fj.getNumCols() && sn.fj.get(i, joinIdx) == 0) {
                    joinIdx++;
                }
                // Check if joinIdx is out of bounds (no join found for this fork)
                if (joinIdx >= sn.fj.getNumCols()) {
                    // No join found - this is a pure fork, handle gracefully
                    Matrix emptyResult = new Matrix(1, 1);
                    emptyResult.set(0, 0, currentTime);
                    return emptyResult;
                }
                int s = 0;
                while (s < fjforkmap.length() && (fjforkmap.get(s) != i || fjclassmap.get(s) != r)) {
                    s++;
                }
                toMerge.add(s);
                Matrix paths = findPaths(sn, P, i, joinIdx, r, toMerge, QN,
                        TN, 0, fjclassmap, fjforkmap, nonfjmodel, new java.util.HashSet<Integer>(onPath));
                // The inner join fires on the k-th branch completion, k = the branch
                // count on a standard join and the declared quorum on a PARTIAL one.
                int kreq = SnJoinQuorum.snJoinQuorum(sn, nonfjmodel.getNodes().get(joinIdx),
                        nonfjmodel.getJobClasses().get(r), paths.length());
                double d0 = FJ_ordstat_exp.fj_ordstat_exp(paths, kreq);
                double innerSync = Math.max(d0 - paths.elementSum() / paths.length(), 0);
                for (int cls : toMerge) {
                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getJobClasses().get(cls),
                            Exp.fitMean(innerSync));
                }
                toMerge.remove(toMerge.size() - 1);
                ri = ri.concatCols(findPaths(sn, P, joinIdx, endNode, r,
                        toMerge, QN, TN, currentTime + d0, fjclassmap,
                        fjforkmap, nonfjmodel, new java.util.HashSet<Integer>(onPath)));
            } else {
                ri = ri.concatCols(findPaths(sn, P, i, endNode, r, toMerge,
                        QN, TN, currentTime + qLen / tput, fjclassmap, fjforkmap,
                        nonfjmodel, new java.util.HashSet<Integer>(onPath)));
            }
        }
        return ri;
    }

    /**
     * Finds the response times along each path leading out of curNode up to (and not including) endNode
     * Variant for models with class switching
     */
    public static Matrix findPathsCS(NetworkStruct sn, Matrix P, int curNode, int endNode, int curClass, ArrayList<Integer> toMerge,
                                     Matrix QN, Matrix TN, double currentTime, Matrix fjclassmap, Matrix fjforkmap,
                                     Network nonfjmodel) {
        return findPathsCS(sn, P, curNode, endNode, curClass, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap,
                nonfjmodel, new java.util.HashSet<Long>());
    }

    /**
     * Enumerates the SIMPLE paths only: the call classes carry a geometric loop
     * (server -> Aux -> server) whenever a call mean exceeds one, so the routing
     * graph between a fork and its join is cyclic and the path set would be
     * infinite without ONPATH. A repeated visit adds no new branch, its residence
     * time is already carried by QN/TN at the station.
     */
    public static Matrix findPathsCS(NetworkStruct sn, Matrix P, int curNode, int endNode, int curClass, ArrayList<Integer> toMerge,
                                     Matrix QN, Matrix TN, double currentTime, Matrix fjclassmap, Matrix fjforkmap,
                                     Network nonfjmodel, java.util.Set<Long> onPath) {
        if (curNode == endNode) {
            double qLen = 0;
            double tput = 0;
            for (int s : toMerge) {
                qLen += QN.get((int) sn.nodeToStation.get(curNode), s);
                tput += TN.get((int) sn.nodeToStation.get(curNode), s);
            }
            Matrix ri = new Matrix(1, 1);
            ri.set(0, 0, currentTime - qLen / tput);
            return ri;
        }
        int orignodes = nonfjmodel.getStruct(false).rtorig.get(nonfjmodel.getClasses().get(0)).get(nonfjmodel.getClasses().get(0)).getNumCols();
        Matrix ri = new Matrix(1, 0);
        long here = (long) curClass * orignodes + curNode;
        onPath.add(here);
        for (int transition = 0; transition < P.getNumCols(); transition++) {
            if (P.get(curClass * orignodes + curNode, transition) == 0) {
                continue;
            }
            ArrayList<Integer> curMerge = new ArrayList<>(toMerge);
            int nextClass = (int) FastMath.floor((transition) / (double) orignodes);
            int nextNode = transition - (nextClass) * orignodes;
            if (nextNode != endNode && onPath.contains((long) transition)) {
                continue;
            }
            curMerge.set(0, nextClass);
            double qLen = 0;
            double tput = 1;
            if (sn.nodeToStation.get(nextNode) > -1) {
                tput = 0;
                for (int s : curMerge) {
                    qLen += QN.get((int) sn.nodeToStation.get(nextNode), s);
                    tput += TN.get((int) sn.nodeToStation.get(nextNode), s);
                }
            }
            if (sn.nodetype.get(nextNode) == NodeType.Fork) {
                int joinIdx = 0;
                while (joinIdx < sn.fj.getNumCols() && sn.fj.get(nextNode, joinIdx) == 0) {
                    joinIdx++;
                }
                // Check if joinIdx is out of bounds (no join found for this fork)
                if (joinIdx >= sn.fj.getNumCols()) {
                    // No join found - this is a pure fork, handle gracefully
                    Matrix emptyResult = new Matrix(1, 1);
                    emptyResult.set(0, 0, currentTime);
                    return emptyResult;
                }
                int s = 0;
                while (s < fjforkmap.length() && (fjforkmap.get(s) != nextNode || fjclassmap.get(s) != curClass)) {
                    s++;
                }
                toMerge.add(s);
                Matrix paths = findPathsCS(sn, P, nextNode, joinIdx, curClass, toMerge, QN,
                        TN, 0, fjclassmap, fjforkmap, nonfjmodel, new java.util.HashSet<Long>(onPath));
                // The inner join fires on the k-th branch completion, k = the branch
                // count on a standard join and the declared quorum on a PARTIAL one.
                int kreq = SnJoinQuorum.snJoinQuorum(sn, nonfjmodel.getNodes().get(joinIdx),
                        nonfjmodel.getJobClasses().get(curClass), paths.length());
                double d0 = FJ_ordstat_exp.fj_ordstat_exp(paths, kreq);
                double innerSync = Math.max(d0 - paths.elementSum() / paths.length(), 0);
                // Match MATLAB: loop over [curMerge, s] to include inner auxiliary class
                ArrayList<Integer> mergeWithS = new ArrayList<>(curMerge);
                mergeWithS.add(s);
                for (int cls : mergeWithS) {
                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getJobClasses().get(cls),
                            Exp.fitMean(innerSync));
                }
                toMerge.remove(toMerge.size() - 1);
                ri = ri.concatCols(findPathsCS(sn, P, joinIdx, endNode, nextClass,
                        curMerge, QN, TN, currentTime + d0, fjclassmap,
                        fjforkmap, nonfjmodel, new java.util.HashSet<Long>(onPath)));
            } else {
                ri = ri.concatCols(findPathsCS(sn, P, nextNode, endNode, nextClass, curMerge,
                        QN, TN, currentTime + qLen / tput, fjclassmap, fjforkmap,
                        nonfjmodel, new java.util.HashSet<Long>(onPath)));
            }
        }
        return ri;
    }

    /**
     * Heidelberger-Trivedi fork-join queueing network transformation.
     * Transforms the queueing network containing a FJ subsystem into a queueing network without one.
     * Fork nodes changed to Router nodes. Join nodes changed to Delay nodes.
     * One artificial class is created for each parallel branch and for each class.
     * Another delay is added to model the sojourn time of the original classes.
     * --
     * This approach is derived by PHILIP HEIDELBERGER and KISHOR S. TRIVEDI in
     * "Analytic Queueing Models for Programs with Internal Concurrency"
     *
     * @param model - the original network
     * @return - queueing network with no FJ system, the class and the fork maps for the artificial classes, and
     * the auxiliary delay map (each join node is mapped to a corresponding auxiliary delay).
     */
    public static Ret.FJApprox ht(Network model) {
        NetworkStruct sn = model.getStruct(true);
        HashMap<Integer, Integer> fjclassmap = new HashMap<>(); // s = fjclassmap(r) for auxiliary class r gives the index s of the original class
        int fjclassmapSize = 0;
        HashMap<Integer, Integer> fjforkmap = new HashMap<>(); // f = fjforkmap(r) for auxiliary class r gives the associated fork node f
        int fjforkmapSize = 0;
        Network nonfjmodel = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(model);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            nonfjmodel = (Network) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            line_error(mfilename(new Object() {
            }), "Could not create a copy of the model in the Heidelberger-Trivedi method");
            return null;
        }
        nonfjmodel.setAllowReplace(true);
        Map<JobClass, Map<JobClass, Matrix>> P = nonfjmodel.getStruct(false).rtorig;
        nonfjmodel.resetNetwork(true);
        nonfjmodel.resetStruct();
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);
        HashMap<Integer, ArrayList<Integer>> forkedClasses = new HashMap<>();
        ArrayList<Integer> forkIndexes = new ArrayList<>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Fork) {
                forkIndexes.add(i);
            }
        }

        // d = fj_auxiliary_delays{j} for join j gives the additional delay d to mimic the sojourn time of the original classes
        HashMap<Integer, Integer> fj_auxiliary_delays = new HashMap<>();

        // Replace each fork with a router
        for (int f : forkIndexes) {
            Node oldNode = nonfjmodel.getNodes().get(f);
            Router newRouter = new Router(nonfjmodel, oldNode.getName());
            ArrayList<Integer> forked = new ArrayList<>();
            for (int i = 0; i < Vnodes.getNumCols(); i++) {
                if (Vnodes.get(f, i) > GlobalConstants.Zero) {
                    forked.add(i);
                }
            }
            forkedClasses.put(f, forked);
        }

        // Replace each join with a delay
        for (int j = 0; j < sn.nodetype.size(); j++) {
            if (sn.nodetype.get(j) != NodeType.Join) {
                continue;
            }
            Node oldNode = nonfjmodel.getNodes().get(j);
            Delay auxiliaryDelay = new Delay(nonfjmodel, oldNode.getName());
            nonfjmodel.getStations().set(model.getNodes().get(j).getStationIdx(), (Station) nonfjmodel.getNodes().get(j));
            for (int c = 0; c < nonfjmodel.getClasses().size(); c++) {
                ((Delay) nonfjmodel.getNodes().get(j)).setService(nonfjmodel.getClasses().get(c), Immediate.getInstance());
            }

            // Add another delay to mimic the sojourn time of the original classes for the artificial classes
            Delay newDelay = new Delay(nonfjmodel, "Auxiliary Delay - " + nonfjmodel.getNodes().get(j).getName());
            fj_auxiliary_delays.put(j, newDelay.getNodeIndex());
            for (int r = 0; r < nonfjmodel.getClasses().size(); r++) {
                newDelay.setService(nonfjmodel.getClasses().get(r), Immediate.getInstance());
            }
            for (JobClass r : P.keySet()) {
                for (JobClass s : P.get(r).keySet()) {
                    Matrix Prs = P.get(r).get(s);
                    Matrix newPrs = new Matrix(Prs.getNumRows() + 1, Prs.getNumCols() + 1);
                    for (int i = 0; i < Prs.getNumRows(); i++) {
                        for (int k = 0; k < Prs.getNumCols(); k++) {
                            newPrs.set(i, k, Prs.get(i, k));
                        }
                    }
                    for (int k = 0; k < Prs.getNumCols(); k++) {
                        newPrs.set(newDelay.getNodeIndex(), k, Prs.get(j, k));
                    }
                    P.get(r).put(s, newPrs);
                }
                for (int i = 0; i < P.get(r).get(r).getNumCols(); i++) {
                    P.get(r).get(r).set(j, i, 0);
                }
                P.get(r).get(r).set(j, newDelay.getNodeIndex(), 1);
            }
        }

        nonfjmodel.setConnectionMatrix(new Matrix(nonfjmodel.getNodes().size(), nonfjmodel.getNodes().size()));

        // Create the artificial classes
        for (int f : forkIndexes) {
            List<OutputStrategy> outputStrategies = model.getNodes().get(f).getOutput().getOutputStrategies();
            int joinIdx = -1;
            for (int i = 0; i < sn.fj.getNumCols(); i++) {
                if (sn.fj.get(f, i) != 0) {
                    if (joinIdx == -1) {
                        joinIdx = i;
                    } else {
                        line_error(mfilename(new Object() {
                        }), "LINE supports at present only a single join station per fork node.");
                    }
                }
            }
            
            // Safety check: if no join found for this fork, skip processing to avoid matrix bounds errors
            if (joinIdx == -1) {
                line_warning(mfilename(new Object() {}), String.format("No join node found for fork node %d in ht method. This fork-join structure is not supported.",
                        f));
            }
            
            ArrayList<Integer> forkedChains = new ArrayList<>();
            for (int i = 0; i < sn.chains.getNumRows(); i++) {
                double sum = 0;
                for (int c : forkedClasses.get(f)) {
                    sum += sn.chains.get(i, c);
                }
                if (sum != 0) {
                    forkedChains.add(i);
                }
            }
            for (int fc : forkedChains) {
                HashMap<Ret.FJAuxClassKey, JobClass> auxClasses = new HashMap<>();
                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0 || sn.nodevisits.get(fc).get(f, r) == 0) {
                        continue;
                    }
                    // Assumption: every class forks into exactly the same parallel branches
                    int parallelBranches = 0;
                    for (OutputStrategy o : outputStrategies) {
                        if (o.getJobClass() == model.getClasses().get(r) && o.getDestination() != null) {
                            parallelBranches++;
                        }
                    }
                    for (int par = 0; par < parallelBranches; par++) {
                        // One artificial class for each parallel branch
                        if (((Forker) model.getNodes().get(f).getOutput()).tasksPerLink > 1) {
                            line_error(mfilename(new Object() {
                            }), "Multiple tasks per link are not supported in H-T.");
                        }
                        int auxPopulation = (int) (((Forker) model.getNodes().get(f).getOutput()).tasksPerLink *
                                ((ClosedClass) model.getClasses().get(r)).getPopulation());
                        auxClasses.put(new Ret.FJAuxClassKey(r, par), new ClosedClass(nonfjmodel,
                                nonfjmodel.getClasses().get(r).getName() + "." + nonfjmodel.getNodes().get(f).getName()
                                        + ".B" + par, auxPopulation, (Station) nonfjmodel.getNodes().get(fj_auxiliary_delays.get(joinIdx)), 0));
                        fjclassmap.put(auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex() - 1, nonfjmodel.getClasses().get(r).getIndex() - 1);
                        if (auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex() - 1 >= fjclassmapSize) {
                            fjclassmapSize = auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex();
                        }
                        fjforkmap.put(auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex() - 1, f);
                        if (auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex() - 1 >= fjforkmapSize) {
                            fjforkmapSize = auxClasses.get(new Ret.FJAuxClassKey(r, par)).getIndex();
                        }

                        // Set the service rates at the join node and at the stations
                        for (int i = 0; i < sn.nnodes; i++) {
                            if (sn.isstation.get(i) != 0) {
                                switch (sn.nodetype.get(i)) {
                                    case Join:
                                        ((jline.lang.nodes.Queue) nonfjmodel.getNodes().get(i)).setService(auxClasses.get(new Ret.FJAuxClassKey(r, par)),
                                                Immediate.getInstance());
                                        break;
                                    case Source:
                                    case Fork:
                                        // No-op
                                        break;
                                    default:
                                        Distribution distributionCopy = null;
                                        try {
                                            ByteArrayOutputStream bos = new ByteArrayOutputStream();
                                            ObjectOutputStream out = new ObjectOutputStream(bos);
                                            out.writeObject(((jline.lang.nodes.Queue) model.getNodes().get(i)).getService(model.getClasses().get(r)));
                                            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                                            ObjectInputStream in = new ObjectInputStream(bis);
                                            distributionCopy = (Distribution) in.readObject();
                                        } catch (IOException | ClassNotFoundException e) {
                                            line_error(mfilename(new Object() {
                                            }), "Could not copy the distribution of the original class in H-T");
                                        }
                                        ((jline.lang.nodes.Queue) nonfjmodel.getNodes().get(i)).setService(auxClasses.get(new
                                                Ret.FJAuxClassKey(r, par)), distributionCopy);
                                }
                            }
                        }
                        ((jline.lang.nodes.Queue) nonfjmodel.getNodes().get(fj_auxiliary_delays.get(joinIdx))).setService(auxClasses.get(new Ret.FJAuxClassKey(r, par)),
                                Immediate.getInstance());
                    }
                }

                // Set the routing of the artificial classes
                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0 || sn.nodevisits.get(fc).get(f, r) == 0) {
                        continue;
                    }
                    for (int s = 0; s < sn.chains.getNumCols(); s++) {
                        if (sn.chains.get(fc, s) == 0 || sn.nodevisits.get(fc).get(f, s) == 0) {
                            continue;
                        }
                        int par = 0;
                        for (OutputStrategy o : outputStrategies) {
                            if (o.getJobClass() != model.getClasses().get(r)) {
                                continue;
                            }
                            JobClass rpar = auxClasses.get(new Ret.FJAuxClassKey(r, par));
                            JobClass spar = auxClasses.get(new Ret.FJAuxClassKey(s, par));
                            if (!P.containsKey(rpar)) {
                                P.put(rpar, new HashMap<>());
                            }
                            Map<JobClass, Matrix> m = P.get(rpar);
                            m.put(spar, new Matrix(P.get(nonfjmodel.getJobClassFromIndex(r)).get(nonfjmodel.getJobClassFromIndex(s))));
                            for (int i = 0; i < m.get(spar).getNumCols(); i++) {
                                m.get(spar).set(f, i, 0);
                            }
                            m.get(spar).set(f, o.getDestination().getNodeIndex(), 1);
                            m.get(spar).set(joinIdx, fj_auxiliary_delays.get(joinIdx), 1);
                            for (int i = 0; i < m.get(spar).getNumCols(); i++) {
                                m.get(spar).set(fj_auxiliary_delays.get(joinIdx), i, 0);
                            }
                            m.get(spar).set(fj_auxiliary_delays.get(joinIdx), f, 1);
                            par++;
                        }
                        // Route the original classes straight to the join to avoid the interference with the artificial classes
                        Matrix Prs = P.get(nonfjmodel.getJobClassFromIndex(r)).get(nonfjmodel.getJobClassFromIndex(s));
                        for (int i = 0; i < Prs.getNumCols(); i++) {
                            Prs.set(f, i, 0);
                        }
                        Prs.set(f, joinIdx, 1);
                    }
                }
            }
        }
        RoutingMatrix routingMatrix = new RoutingMatrix(nonfjmodel, nonfjmodel.getClasses(), nonfjmodel.getNodes());
        int numNodesHT = nonfjmodel.getNumberOfNodes();
        for (JobClass r : P.keySet()) {
            for (JobClass s : P.get(r).keySet()) {
                Matrix Prs = P.get(r).get(s);
                // Iterate only up to nonfjmodel's node count since resetNetwork may have removed ClassSwitch/Logger nodes
                int maxRows = Math.min(Prs.getNumRows(), numNodesHT);
                int maxCols = Math.min(Prs.getNumCols(), numNodesHT);
                for (int i = 0; i < maxRows; i++) {
                    for (int j = 0; j < maxCols; j++) {
                        if (Prs.get(i, j) != 0) {
                            routingMatrix.set(r, s, nonfjmodel.getNodes().get(i), nonfjmodel.getNodes().get(j), Prs.get(i, j));
                        }
                    }
                }
            }
        }
        nonfjmodel.relink(routingMatrix);
        Matrix fjclassmapMatrix = new Matrix(1, fjclassmapSize);
        fjclassmapMatrix.fill(-1);
        for (int r : fjclassmap.keySet()) {
            fjclassmapMatrix.set(r, fjclassmap.get(r));
        }
        Matrix fjforkmapMatrix = new Matrix(1, fjforkmapSize);
        fjforkmapMatrix.fill(-1);
        for (int r : fjforkmap.keySet()) {
            fjforkmapMatrix.set(r, fjforkmap.get(r));
        }
        return new Ret.FJApprox(nonfjmodel, fjclassmapMatrix, fjforkmapMatrix, fj_auxiliary_delays, null);
    }

    /**
     * Fork-Join Transform approach with default forkLambda parameter
     */
    public static Ret.FJApprox mmt(Network model) {
        return mmt(model, new Matrix(1, model.getNumberOfClasses()).add(GlobalConstants.FineTol));
    }

    /**
     * Fork-Join Transform approach to evaluate queueing networks including fork-join systems. An equivalent network is
     * created where the fork nodes are replaced by routers, the join nodes are replaced by delays, and the parallelism
     * induced by a fork-join system is simulated through the addition of artificial open customer classes.
     *
     * @param model      - the original queueing network
     * @param forkLambda - the arrival rates of the artificial classes
     * @return - the equivalent queueing network with the fork-join systems replaced with other nodes, a mapping of the
     * artificial classes and their corresponding original classes, a mapping of the artificial classes and their FJ
     * systems, and the fanout of each artificial class
     */
    public static Ret.FJApprox mmt(Network model, Matrix forkLambda) {
        NetworkStruct sn = model.getStruct(true);
        HashMap<Integer, Integer> fjclassmap = new HashMap<>(); // s = fjclassmap(r) for auxiliary class r gives the index s of the original class
        int fjclassmapSize = 0;
        HashMap<Integer, Integer> fjforkmap = new HashMap<>(); // f = fjforkmap(r) for auxiliary class r gives the associated fork node f
        int fjforkmapSize = 0;
        HashMap<Integer, Integer> fanout = new HashMap<>(); // fo = fanout(r) is the number of output jobs across all links for the (fork f,class s) pair modelled by auxiliary class r
        int fanoutSize = 0;
        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
        Map<String, Integer> serviceSrc = new HashMap<String, Integer>();
        Map<String, Object> serviceSrcObj = new HashMap<String, Object>();
        Set<String> immediateSlots = new HashSet<String>();
        Map<Integer, Integer> auxArrivalSrc = new HashMap<Integer, Integer>();
        Set<Integer> auxDisabled = new HashSet<Integer>();
        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.isstation.get(i) == 0) {
                continue;
            }
            NodeType nt = sn.nodetype.get(i);
            if (nt == NodeType.Source || nt == NodeType.Fork || nt == NodeType.Join) {
                continue;
            }
            if (!(model.getNodes().get(i) instanceof jline.lang.nodes.Queue)) {
                continue;
            }
            for (int r = 0; r < model.getClasses().size(); r++) {
                Object svc = ((jline.lang.nodes.Queue) model.getNodes().get(i)).getService(model.getClasses().get(r));
                if (svc != null) {
                    serviceSrc.put(Ret.FJApprox.slot(i, r), r);
                    serviceSrcObj.put(Ret.FJApprox.slot(i, r), svc);
                }
            }
        }
        Network nonfjmodel = null;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(model);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            nonfjmodel = (Network) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            line_error(mfilename(new Object() {
            }), "Could not create a copy of the model in the Heidelberger-Trivedi method");
            return null;
        }
        nonfjmodel.setAllowReplace(true);
        NetworkStruct copySn = nonfjmodel.getStruct(false);
        Map<JobClass, Map<JobClass, Matrix>> P = copySn.rtorig;
        // Capture original nodes AFTER getStruct() (which may add ClassSwitch nodes) but BEFORE resetNetwork
        // This ensures originalNodes indices match the routing matrix P indices
        List<Node> originalNodes = new ArrayList<>(nonfjmodel.getNodes());
        nonfjmodel.resetNetwork(true);
        nonfjmodel.resetStruct();
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);
        HashMap<Integer, ArrayList<Integer>> forkedClasses = new HashMap<>();
        ArrayList<Integer> forkIndexes = new ArrayList<>();
        int maxForkIdx = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Fork) {
                forkIndexes.add(i);
                if (i > maxForkIdx) {
                    maxForkIdx = i;
                }
            }
        }

        Matrix origfanout = new Matrix(maxForkIdx + 1, nonfjmodel.getNumberOfClasses());

        // Replace each fork with a router
        for (int f : forkIndexes) {
            // Read output strategies from the original model's Fork (matching MATLAB: model.nodes{f}).
            // Use index comparison for JobClass (not identity) to handle deserialized objects.
            List<OutputStrategy> forkOutputStrategies = model.getNodes().get(f).getOutputStrategies();
            for (JobClass r : P.keySet()) {
                int parallelBranches = 0;
                for (OutputStrategy o : forkOutputStrategies) {
                    if (o.getJobClass().getIndex() == r.getIndex() && o.getDestination() != null) {
                        parallelBranches++;
                    }
                }
                if (parallelBranches > 0) {
                    origfanout.set(f, r.getIndex() - 1, parallelBranches);
                    for (JobClass s : P.get(r).keySet()) {
                        Matrix Prs = P.get(r).get(s);
                        for (int j = 0; j < Prs.getNumCols(); j++) {
                            Prs.set(f, j, Prs.get(f, j) / parallelBranches);
                        }
                    }
                } else {
                    // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                    for (JobClass s : P.get(r).keySet()) {
                        Matrix Prs = P.get(r).get(s);
                        for (int j = 0; j < Prs.getNumCols(); j++) {
                            Prs.set(f, j, 0.0);
                        }
                    }
                }
            }
            Node oldNode = nonfjmodel.getNodes().get(f);
            Router newRouter = new Router(nonfjmodel, oldNode.getName());
            ArrayList<Integer> forked = new ArrayList<>();
            for (int i = 0; i < Vnodes.getNumCols(); i++) {
                if (Vnodes.get(f, i) > GlobalConstants.Zero) {
                    forked.add(i);
                }
            }
            forkedClasses.put(f, forked);
        }
        // Replace each join with a delay
        for (int j = 0; j < sn.nodetype.size(); j++) {
            if (sn.nodetype.get(j) != NodeType.Join) {
                continue;
            }
            Node oldNode = nonfjmodel.getNodes().get(j);
            Delay auxiliaryDelay = new Delay(nonfjmodel, oldNode.getName());
            nonfjmodel.getStations().set(model.getNodes().get(j).getStationIdx(), (Station) nonfjmodel.getNodes().get(j));
            for (int c = 0; c < nonfjmodel.getClasses().size(); c++) {
                ((Delay) nonfjmodel.getNodes().get(j)).setService(nonfjmodel.getClasses().get(c), Immediate.getInstance());
                immediateSlots.add(Ret.FJApprox.slot(j, c));
            }
        }
        Source source = null;
        Sink sink = null;
        if (nonfjmodel.hasOpenClasses()) {
            source = nonfjmodel.getSource();
            sink = nonfjmodel.getSink();
        } else {
            source = new Source(nonfjmodel, "Source");
            sink = new Sink(nonfjmodel, "Sink");
            for (JobClass r : P.keySet()) {
                for (JobClass s : P.get(r).keySet()) {
                    Matrix Prs = P.get(r).get(s);
                    Matrix newPrs = new Matrix(Prs.getNumRows() + 2, Prs.getNumCols() + 2);
                    for (int i = 0; i < Prs.getNumRows(); i++) {
                        for (int j = 0; j < Prs.getNumCols(); j++) {
                            newPrs.set(i, j, Prs.get(i, j));
                        }
                    }
                    P.get(r).put(s, newPrs);
                }
            }
        }
        nonfjmodel.setConnectionMatrix(new Matrix(nonfjmodel.getNumberOfNodes(), nonfjmodel.getNumberOfNodes()));
        ArrayList<Integer> allAuxClassIndices = new ArrayList<>(); // 0-based aux class indices for post-relink routing fix

        for (int f : forkIndexes) {
            int joinIdx = -1;
            for (int i = 0; i < sn.fj.getNumCols(); i++) {
                if (sn.fj.get(f, i) != 0) {
                    if (joinIdx == -1) {
                        joinIdx = i;
                    } else {
                        line_error(mfilename(new Object() {
                        }), "LINE supports at present only a single join station per fork node.");
                    }
                }
            }

            // Safety check: if no join found for this fork, report error with helpful message
            if (joinIdx == -1) {
                line_warning(mfilename(new Object() {}), String.format("No join node found for fork node %d. ",
                        f, sn.fj.getNumRows(), sn.fj.getNumCols(), sn.fj.isEmpty() || sn.fj.elementSum() == 0));
            }

            ArrayList<Integer> forkedChains = new ArrayList<>();
            for (int i = 0; i < sn.chains.getNumRows(); i++) {
                double sum = 0;
                for (int c : forkedClasses.get(f)) {
                    sum += sn.chains.get(i, c);
                }
                if (sum != 0) {
                    forkedChains.add(i);
                }
            }
            for (int fc : forkedChains) {
                ArrayList<JobClass> oclass = new ArrayList<>();

                // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                ArrayList<Integer> inchainBfs = new ArrayList<>();
                for (int ci = 0; ci < sn.chains.getNumCols(); ci++) {
                    if (sn.chains.get(fc, ci) != 0) {
                        inchainBfs.add(ci);
                    }
                }
                int refStationIdx = (int) sn.refstat.get(inchainBfs.get(0));
                int refStatefulIdx = (int) sn.stationToStateful.get(refStationIdx);
                int refNodeIdx = (int) sn.statefulToNode.get(refStatefulIdx);
                int K_bfs = sn.nclasses;
                boolean[][] reachableFromRef = new boolean[sn.nnodes][K_bfs];
                ArrayList<int[]> bfsQueue = new ArrayList<>();
                for (int ci : inchainBfs) {
                    reachableFromRef[refNodeIdx][ci] = true;
                    bfsQueue.add(new int[]{refNodeIdx, ci});
                }
                while (!bfsQueue.isEmpty()) {
                    int[] pair = bfsQueue.remove(0);
                    int cn = pair[0]; int cc = pair[1];
                    for (int destNode = 0; destNode < sn.nnodes; destNode++) {
                        for (int dci : inchainBfs) {
                            int srcIdx = cn * K_bfs + cc;
                            int dstIdx = destNode * K_bfs + dci;
                            if (!reachableFromRef[destNode][dci] && sn.rtnodes.get(srcIdx, dstIdx) > 0) {
                                reachableFromRef[destNode][dci] = true;
                                bfsQueue.add(new int[]{destNode, dci});
                            }
                        }
                    }
                }

                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0) {
                        continue;
                    }
                    oclass.add(new OpenClass(nonfjmodel, nonfjmodel.getJobClasses().get(r).getName() + "." +
                            nonfjmodel.getNodes().get(f).getName()));
                    // Store the 0-indexed original class index directly to ensure correct mapping
                    fjclassmap.put(oclass.get(oclass.size() - 1).getIndex() - 1, r);
                    if (oclass.get(oclass.size() - 1).getIndex() - 1 >= fjclassmapSize) {
                        fjclassmapSize = oclass.get(oclass.size() - 1).getIndex();
                    }
                    fjforkmap.put(oclass.get(oclass.size() - 1).getIndex() - 1, f);
                    if (oclass.get(oclass.size() - 1).getIndex() - 1 >= fjforkmapSize) {
                        fjforkmapSize = oclass.get(oclass.size() - 1).getIndex();
                    }
                    int s = fjclassmap.get(oclass.get(oclass.size() - 1).getIndex() - 1);
                    // fanout is the SIBLING count, links times tasksPerLink, which is what
                    // the auxiliary open class's rate (fanout-1)*forkLambda must carry: one
                    // sibling is the closed token's own, the other fanout-1 are open
                    // traffic. FJFixedPoint synchronises on the same count by replicating
                    // each branch time tasksPerLink times.
                    int fanoutValue = (int) (origfanout.get(f, r) * ((Forker) model.getNodes().get(f).getOutput()).tasksPerLink);
                    fanout.put(oclass.get(oclass.size() - 1).getIndex() - 1, fanoutValue);
                    allAuxClassIndices.add(oclass.get(oclass.size() - 1).getIndex() - 1); // 0-based
                    boolean disableAux = origfanout.get(f, r) == 0 || !reachableFromRef[f][r];
                    int auxIdx = oclass.get(oclass.size() - 1).getIndex() - 1;
                    if (disableAux) {
                        source.setArrival(oclass.get(oclass.size() - 1), Disabled.getInstance());
                        auxDisabled.add(auxIdx);
                    } else {
                        source.setArrival(oclass.get(oclass.size() - 1), new Exp(forkLambda.get(r)));
                    }
                    // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                    auxArrivalSrc.put(auxIdx, r);
                    // joins are now Delays, let us set their service time
                    for (int i = 0; i < sn.nnodes; i++) {
                        if (sn.isstation.get(i) != 0) {
                            switch (sn.nodetype.get(i)) {
                                case Join:
                                    ((Delay) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), Immediate.getInstance());
                                    // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                                    immediateSlots.add(Ret.FJApprox.slot(i, auxIdx));
                                    break;
                                case Source:
                                case Fork:
                                    // no-op
                                    break;
                                default:
                                    Distribution distributionCopy = null;
                                    try {
                                        ByteArrayOutputStream bos = new ByteArrayOutputStream();
                                        ObjectOutputStream out = new ObjectOutputStream(bos);
                                        out.writeObject(((jline.lang.nodes.Queue) model.getNodes().get(i)).getService(model.getClasses().get(r)));
                                        ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                                        ObjectInputStream in = new ObjectInputStream(bis);
                                        distributionCopy = (Distribution) in.readObject();
                                    } catch (IOException | ClassNotFoundException e) {
                                        line_error(mfilename(new Object(){}), "Could not copy the distribution of the original class in FJT");
                                    }
                                    ((jline.lang.nodes.Queue) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), distributionCopy);
                                    serviceSrc.put(Ret.FJApprox.slot(i, auxIdx), r);
                                    serviceSrcObj.put(Ret.FJApprox.slot(i, auxIdx),
                                            ((jline.lang.nodes.Queue) model.getNodes().get(i)).getService(model.getClasses().get(r)));
                            }
                        }
                    }
                }
                int rIdx = 0;
                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0) {
                        continue;
                    }
                    int sIdx = 0;
                    for (int s = 0; s < sn.chains.getNumCols(); s++) {
                        if (sn.chains.get(fc, s) == 0) {
                            continue;
                        }
                        if (!P.containsKey(oclass.get(rIdx))) {
                            P.put(oclass.get(rIdx), new HashMap<>());
                        }
                        P.get(oclass.get(rIdx)).put(oclass.get(sIdx),
                                new Matrix(rtBlock(P, nonfjmodel.getJobClasses().get(r),
                                        nonfjmodel.getJobClasses().get(s))));
                        sIdx++;
                    }
                    rIdx++;
                }
                rIdx = 0;
                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0) {
                        continue;
                    }
                    int sIdx = 0;
                    for (int s = 0; s < sn.chains.getNumCols(); s++) {
                        if (sn.chains.get(fc, s) == 0) {
                            continue;
                        }
                        Matrix Prs = rtBlock(P, oclass.get(rIdx), oclass.get(sIdx));
                        for (int i = 0; i < Prs.getNumCols(); i++) {
                            Prs.set(source.getNodeIndex(), i, 0.0);
                            if (joinIdx >= 0)
                                Prs.set(joinIdx, i, 0.0);
                        }
                        sIdx++;
                    }
                    if (origfanout.get(f, r) > 0) {
                        Matrix Prr = rtBlock(P, oclass.get(rIdx), oclass.get(rIdx));
                        Prr.set(source.getNodeIndex(), f, 1);
                        if (joinIdx >= 0) {
                            Prr.set(joinIdx, sink.getNodeIndex(), 1);
                        }
                    }
                    rIdx++;
                }

                // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                ArrayList<Integer> inchain = new ArrayList<>();
                for (int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) != 0) {
                        inchain.add(r);
                    }
                }
                boolean allForked = true;
                for (int r : inchain) {
                    if (origfanout.get(f, r) == 0) {
                        allForked = false;
                        break;
                    }
                }
                if (allForked) {
                    // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                    int pnnodes = rtBlock(P, nonfjmodel.getJobClasses().get(inchain.get(0)),
                            nonfjmodel.getJobClasses().get(inchain.get(0))).getNumRows();
                    boolean[] fjNodeMask = new boolean[pnnodes];
                    fjNodeMask[f] = true;
                    if (joinIdx >= 0) fjNodeMask[joinIdx] = true;
                    // Source and Sink are mmt infrastructure, always in scope
                    for (int nd = 0; nd < pnnodes && nd < nonfjmodel.getNodes().size(); nd++) {
                        if (nonfjmodel.getNodes().get(nd) instanceof jline.lang.nodes.Source
                                || nonfjmodel.getNodes().get(nd) instanceof jline.lang.nodes.Sink) {
                            fjNodeMask[nd] = true;
                        }
                    }
                    // Class-aware BFS: find all (node, class) pairs reachable from Fork
                    int maxclass = 0;
                    for (int c : inchain) { if (c > maxclass) maxclass = c; }
                    maxclass++;
                    boolean[][] visited = new boolean[pnnodes][maxclass];
                    ArrayList<int[]> bfsQ = new ArrayList<>();
                    // Seed: classes entering the fork
                    for (int r : inchain) {
                        for (int s : inchain) {
                            JobClass rClass = nonfjmodel.getJobClasses().get(r);
                            JobClass sClass = nonfjmodel.getJobClasses().get(s);
                            if (!P.containsKey(rClass) || !P.get(rClass).containsKey(sClass)) continue;
                            Matrix Prs = P.get(rClass).get(sClass);
                            boolean hasEntry = false;
                            for (int j = 0; j < Prs.getNumCols(); j++) {
                                if (Prs.get(f, j) > 0) { hasEntry = true; break; }
                            }
                            if (hasEntry && !visited[f][r]) {
                                visited[f][r] = true;
                                bfsQ.add(new int[]{f, r});
                            }
                        }
                    }
                    while (!bfsQ.isEmpty()) {
                        int[] pair = bfsQ.remove(0);
                        int cn = pair[0]; int cc = pair[1];
                        for (int s : inchain) {
                            JobClass ccClass = nonfjmodel.getJobClasses().get(cc);
                            JobClass sClass = nonfjmodel.getJobClasses().get(s);
                            if (!P.containsKey(ccClass) || !P.get(ccClass).containsKey(sClass)) continue;
                            Matrix Prs = P.get(ccClass).get(sClass);
                            for (int nd = 0; nd < pnnodes; nd++) {
                                if (Prs.get(cn, nd) > 0 && !visited[nd][s]) {
                                    visited[nd][s] = true;
                                    fjNodeMask[nd] = true;
                                    if (joinIdx < 0 || nd != joinIdx) {
                                        bfsQ.add(new int[]{nd, s});
                                    }
                                }
                            }
                        }
                    }
                    // Clear outgoing aux routing at nodes outside fork-join scope.
                    for (int nd = 0; nd < pnnodes; nd++) {
                        if (fjNodeMask[nd]) continue;
                        for (int ri = 0; ri < oclass.size(); ri++) {
                            for (int si = 0; si < oclass.size(); si++) {
                                if (!P.containsKey(oclass.get(ri))) continue;
                                if (!P.get(oclass.get(ri)).containsKey(oclass.get(si))) continue;
                                Matrix Prs = P.get(oclass.get(ri)).get(oclass.get(si));
                                for (int j = 0; j < Prs.getNumCols(); j++) {
                                    Prs.set(nd, j, 0.0);
                                }
                            }
                        }
                    }
                }
            }
        }
        RoutingMatrix routingMatrix = new RoutingMatrix(nonfjmodel, nonfjmodel.getClasses(), nonfjmodel.getNodes());
        for (JobClass r : P.keySet()) {
            for (JobClass s : P.get(r).keySet()) {
                Matrix Prs = P.get(r).get(s);
                // Iterate only up to nonfjmodel's node count since resetNetwork may have removed ClassSwitch/Logger nodes
                int numNodesMMT = nonfjmodel.getNumberOfNodes();
                int maxRows = Math.min(Prs.getNumRows(), numNodesMMT);
                int maxCols = Math.min(Prs.getNumCols(), numNodesMMT);
                for (int i = 0; i < maxRows; i++) {
                    for (int j = 0; j < maxCols; j++) {
                        if (Prs.get(i, j) != 0) {
                            routingMatrix.set(r, s, nonfjmodel.getNodes().get(i), nonfjmodel.getNodes().get(j), Prs.get(i, j));
                        }
                    }
                }
            }
        }
        nonfjmodel.relink(routingMatrix);
        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
        for (int nd = 0; nd < nonfjmodel.getNodes().size(); nd++) {
            List<OutputStrategy> os = nonfjmodel.getNodes().get(nd).getOutputStrategies();
            for (int ci : allAuxClassIndices) {
                for (OutputStrategy o : os) {
                    if (o.getJobClass().getIndex() - 1 == ci && o.getRoutingStrategy() != RoutingStrategy.PROB) {
                        o.setRoutingStrategy(RoutingStrategy.DISABLED);
                    }
                }
            }
        }
        for (int f : forkIndexes) {
            for (OutputStrategy o : nonfjmodel.getNodeByIndex(f).getOutputStrategies()) {
                if (o.getRoutingStrategy() == RoutingStrategy.RAND) {
                    o.setRoutingStrategy(RoutingStrategy.DISABLED);
                }
            }
        }

        Matrix fjclassmapMatrix = new Matrix(1, fjclassmapSize);
        fjclassmapMatrix.fill(-1);
        for (int r : fjclassmap.keySet()) {
            fjclassmapMatrix.set(r, fjclassmap.get(r));
        }
        Matrix fjforkmapMatrix = new Matrix(1, fjforkmapSize);
        fjforkmapMatrix.fill(-1);
        for (int r : fjforkmap.keySet()) {
            fjforkmapMatrix.set(r, fjforkmap.get(r));
        }
        // The transformed model's nodes were copied from a model with FEWER
        // classes, so each carries a state row of the old width. Network.initDefault
        // PRESERVES an existing state, so that stale row would survive and decode
        // to a population of zero: the fluid inner solve then integrates an empty
        // network and every metric comes back 0. Clear the states here, where the
        // class count changes, rather than in each consumer.
        for (jline.lang.nodes.Node nd : nonfjmodel.getNodes()) {
            if (nd.isStateful()) {
                ((jline.lang.nodes.StatefulNode) nd).setState(new Matrix(0, 0));
            }
        }

        Ret.FJApprox mmtReturn = new Ret.FJApprox(nonfjmodel, fjclassmapMatrix, fjforkmapMatrix, null, fanout);
        mmtReturn.serviceSrc = serviceSrc;
        mmtReturn.serviceSrcObj = serviceSrcObj;
        mmtReturn.immediateSlots = immediateSlots;
        mmtReturn.auxArrivalSrc = auxArrivalSrc;
        mmtReturn.auxDisabled = auxDisabled;
        mmtReturn.baseModel = model;
        mmtReturn.forkLambdaInit = forkLambda == null ? null : forkLambda.copy();
        return mmtReturn;
    }

    /**
     * Re-feed a transformed model produced by {@link #mmt} from the current service
     * parameters of the base model it was derived from, so that the transformation
     * can be reused across the iterations of an outer fixed point instead of being
     * rebuilt. SolverLN re-solves each layer once per iteration and only the rates
     * change between iterations; the fork topology the transformation encodes does
     * not. Rebuilding costs a full serialisation deep copy of the model each time.
     * <p>
     * Three kinds of slot, each handled differently:
     * <ul>
     *   <li>{@code serviceSrc}: base-derived, re-read from the base model.</li>
     *   <li>{@code immediateSlots}: owned by the transformation, and RESET rather
     *       than skipped -- the fork loop overwrites every join with its current
     *       synchronisation delay on each pass, so a reused model still holds the
     *       previous outer iteration's converged value. Leaving them alone silently
     *       warm-starts the fork loop.</li>
     *   <li>auxiliary arrivals: owned by the forkLambda fixed point, reset to what
     *       a cold call would have set.</li>
     * </ul>
     * Reading a transformation-owned slot from the base model is the converse error
     * and corrupts the transform outright. Both failure modes are silent.
     *
     * @return false when the provenance cannot be applied, in which case the caller
     *         must fall back to a cold {@link #mmt}.
     */
    public static boolean refreshServicesFromBase(Ret.FJApprox mmtResult) {
        if (mmtResult == null || mmtResult.serviceSrc == null || mmtResult.baseModel == null
                || mmtResult.nonfjmodel == null) {
            return false;
        }
        Network base = mmtResult.baseModel;
        Network nonfj = mmtResult.nonfjmodel;
        int nclasses = base.getClasses().size();

        for (Map.Entry<String, Integer> e : mmtResult.serviceSrc.entrySet()) {
            String[] parts = e.getKey().split(",");
            int i = Integer.parseInt(parts[0]);
            int c = Integer.parseInt(parts[1]);
            int r = e.getValue();
            if (i >= nonfj.getNodes().size() || c >= nonfj.getClasses().size() || r >= nclasses) {
                return false;
            }
            if (!(base.getNodes().get(i) instanceof jline.lang.nodes.Queue)
                    || !(nonfj.getNodes().get(i) instanceof jline.lang.nodes.Queue)) {
                return false;
            }
            Distribution svc = (Distribution) ((jline.lang.nodes.Queue) base.getNodes().get(i))
                    .getService(base.getClasses().get(r));
            if (svc == null) {
                return false;
            }
            // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
            if (mmtResult.serviceSrcObj != null && mmtResult.serviceSrcObj.get(e.getKey()) == svc) {
                continue;
            }
            // Copy mirrors the cold path, which hands every slot its own
            // distribution object rather than aliasing one across classes.
            Distribution svcCopy = null;
            try {
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);
                out.writeObject(svc);
                ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                ObjectInputStream in = new ObjectInputStream(bis);
                svcCopy = (Distribution) in.readObject();
            } catch (IOException | ClassNotFoundException ex) {
                return false;
            }
            ((jline.lang.nodes.Queue) nonfj.getNodes().get(i)).setService(nonfj.getClasses().get(c), svcCopy);
            if (mmtResult.serviceSrcObj != null) {
                mmtResult.serviceSrcObj.put(e.getKey(), svc);
            }
        }

        if (mmtResult.immediateSlots != null) {
            for (String key : mmtResult.immediateSlots) {
                String[] parts = key.split(",");
                int i = Integer.parseInt(parts[0]);
                int c = Integer.parseInt(parts[1]);
                if (i >= nonfj.getNodes().size() || c >= nonfj.getClasses().size()) {
                    return false;
                }
                if (!(nonfj.getNodes().get(i) instanceof ServiceStation)) {
                    return false;
                }
                ((ServiceStation) nonfj.getNodes().get(i)).setService(nonfj.getClasses().get(c),
                        Immediate.getInstance());
            }
        }

        if (mmtResult.auxArrivalSrc != null && !mmtResult.auxArrivalSrc.isEmpty()) {
            Source src = nonfj.getSource();
            if (src == null) {
                return false;
            }
            for (Map.Entry<Integer, Integer> e : mmtResult.auxArrivalSrc.entrySet()) {
                int c = e.getKey();
                int r = e.getValue();
                if (c >= nonfj.getClasses().size()) {
                    return false;
                }
                if (mmtResult.auxDisabled != null && mmtResult.auxDisabled.contains(c)) {
                    src.setArrival(nonfj.getClasses().get(c), Disabled.getInstance());
                } else if (mmtResult.forkLambdaInit != null && r < mmtResult.forkLambdaInit.length()) {
                    src.setArrival(nonfj.getClasses().get(c), new Exp(mmtResult.forkLambdaInit.get(r)));
                }
            }
        }
        return true;
    }

    private static void nested_forks(int startNode, int endNode, Matrix conn, Matrix forks, NetworkStruct sn) {
        if (startNode == endNode) {
            return;
        }
        for (int i = 0; i < conn.getNumCols(); i++) {
            if (conn.get(startNode, i) == 0) {
                continue;
            }
            if (sn.nodetype.get(i) == NodeType.Fork) {
                forks.set(i, 0);
            }
            nested_forks(i, endNode, conn, forks, sn);
        }
    }

    /**
     * Determines a directed acyclic graph of relationships among fork nodes.
     */
    public static Ret.FJsortForks sort_forks(NetworkStruct sn, NetworkStruct nonfjstruct, Matrix fjforkmap, Matrix fjclassmap, Network nonfjmodel) {
        Matrix forks = new Matrix(sn.nodetype.size(), (int) fjclassmap.elementMax() + 1);
        // Initialize forks(f, :) = 1 for all Fork nodes and all original classes
        // This matches MATLAB: forks(find(sn.nodetype == NodeType.Fork), :) = 1
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Fork) {
                for (int r = 0; r < forks.getNumCols(); r++) {
                    forks.set(i, r, 1);
                }
            }
        }
        Matrix parents = new Matrix(1, sn.nodetype.size());
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Fork) {
                parents.set(i, i);
            }
        }
        for (int f = 0; f < sn.nodetype.size(); f++) {
            if (sn.nodetype.get(f) == NodeType.Fork) {
                int joinIdx = 0;
                while (joinIdx < sn.fj.getNumCols() && sn.fj.get(f, joinIdx) == 0) {
                    joinIdx++;
                }
                // Check if joinIdx is out of bounds (no join found for this fork)
                if (joinIdx >= sn.fj.getNumCols()) {
                    line_warning(mfilename(new Object() {}), String.format("No join node found for fork node %d in nested forks processing.", f));
                    continue; // Skip this fork
                }
                for (int s = 0; s < fjforkmap.length(); s++) {
                    if (fjforkmap.get(s) == f) {
                        int r = (int) fjclassmap.get(s);
                        Matrix nested = new Matrix(sn.nodetype.size(), 1);
                        for (int i = 0; i < sn.nodetype.size(); i++) {
                            if (sn.nodetype.get(i) == NodeType.Fork) {
                                nested.set(i, 1);
                            }
                        }
                        nested_forks(f, joinIdx, nonfjstruct.rtorig.get(nonfjmodel.getJobClasses().get(r)).get(nonfjmodel.getJobClasses().get(r)), nested, sn);
                        for (int i = 0; i < forks.getNumRows(); i++) {
                            forks.set(i, r, (int) forks.get(i, r) & (int) nested.get(i));
                        }
                        for (int i = 0; i < nested.length(); i++) {
                            if (nested.get(i) == 0) {
                                parents.set(i, parents.get(f));
                            }
                        }
                    }
                }
            }
        }
        return new Ret.FJsortForks(forks, parents);
    }

    // ========== Chain Aggregation Methods ==========

    /**
     * Result of aggregating chains in a model
     */
    public static class AggregateChainResult {
        private final Network chainModel;
        private final Matrix alpha;
        private final DeaggInfo deaggInfo;

        public AggregateChainResult(Network chainModel, Matrix alpha, DeaggInfo deaggInfo) {
            this.chainModel = chainModel;
            this.alpha = alpha;
            this.deaggInfo = deaggInfo;
        }

        public Network getChainModel() {
            return chainModel;
        }

        public Matrix getAlpha() {
            return alpha;
        }

        public DeaggInfo getDeaggInfo() {
            return deaggInfo;
        }
    }

    /**
     * Deaggregation information for converting chain-level results back to class-level
     */
    public static class DeaggInfo {
        public Matrix alpha;
        public Matrix Vchain;
        public Matrix STchain;
        public Matrix Lchain;
        public Matrix SCVchain;
        public Matrix Nchain;
        public Matrix lambdaChain;
        public boolean[] isOpenChain;
        public Map<Integer, Matrix> inchain;
        public Matrix refstat;
        public Matrix refstatchain;
        public NetworkStruct originalSn;
        public boolean isAggregated;
        public int nclasses;
        public int nchains;
    }

    /**
     * Transform a multi-class model into an equivalent chain-aggregated model
     *
     * This function transforms a queueing network model with multiple classes
     * into a stochastically equivalent model where each chain becomes a single
     * class. Classes belonging to the same chain (i.e., classes that can switch
     * into each other) are merged into one aggregate class.
     *
     * The aggregated model preserves:
     * - Total chain population (closed chains)
     * - Total arrival rate (open chains)
     * - Service demands at each station
     * - Routing structure at the chain level
     *
     * @param model Source Network model with potentially multiple classes per chain
     * @return AggregateChainResult containing the aggregated model, alpha, and deaggInfo
     */
    /**
     * Returns a copy of the model with the specified job class removed,
     * updating all station configurations, routing matrices and class-dependent
     * parameters accordingly. The original model is left unchanged.
     *
     * @param model    the source network model
     * @param jobclass the job class to remove
     * @return a new model without the specified class
     */
    public static Network removeClass(Network model, JobClass jobclass) {
        Network newmodel = model.copy();
        newmodel.removeClass(jobclass);
        return newmodel;
    }

    public static AggregateChainResult aggregateChains(Network model) {
        return aggregateChains(model, "");
    }

    /**
     * Transform a multi-class model into an equivalent chain-aggregated model
     *
     * @param model Source Network model with potentially multiple classes per chain
     * @param suffix Optional suffix for chain class names
     * @return AggregateChainResult containing the aggregated model, alpha, and deaggInfo
     */
    public static AggregateChainResult aggregateChains(Network model, String suffix) {
        if (suffix == null) {
            suffix = "";
        }

        // Get network structure
        NetworkStruct sn = model.getStruct(true);

        // Extract dimensions
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;

        // If each class is its own chain, just return a copy
        if (C == K) {
            Network chainModel = model.copy();
            Matrix alpha = Matrix.eye(M);
            if (alpha.getNumCols() < K) {
                // Expand to M x K
                Matrix newAlpha = new Matrix(M, K);
                for (int i = 0; i < M && i < K; i++) {
                    newAlpha.set(i, i, 1.0);
                }
                alpha = newAlpha;
            }
            DeaggInfo deaggInfo = new DeaggInfo();
            deaggInfo.alpha = alpha;
            deaggInfo.inchain = sn.inchain;
            deaggInfo.originalSn = sn;
            deaggInfo.isAggregated = false;
            return new AggregateChainResult(chainModel, alpha, deaggInfo);
        }

        // Get aggregation parameters from existing API
        Ret.snGetDemands demands = snGetDemandsChain(sn);
        Matrix Lchain = demands.Dchain;
        Matrix STchain = demands.STchain;
        Matrix Vchain = demands.Vchain;
        Matrix alpha = demands.alpha;
        Matrix Nchain = demands.Nchain;
        Matrix SCVchain = demands.SCVchain;
        Matrix refstatchain = demands.refstatchain;

        // Determine which chains are open vs closed
        boolean[] isOpenChain = new boolean[C];
        Matrix lambdaChain = new Matrix(1, C);
        int sourceIdx = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                sourceIdx = i;
                break;
            }
        }
        int sourceStationIdx = -1;
        if (sourceIdx >= 0) {
            sourceStationIdx = (int) sn.nodeToStation.get(sourceIdx);
        }

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            boolean isOpen = false;
            for (int col = 0; col < inchain.getNumCols(); col++) {
                int classIdx = (int) inchain.get(0, col);
                if (Double.isInfinite(sn.njobs.get(classIdx))) {
                    isOpen = true;
                    break;
                }
            }
            isOpenChain[c] = isOpen;
            if (isOpen && sourceStationIdx >= 0) {
                // Sum arrival rates for all classes in this chain
                double totalLambda = 0.0;
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int classIdx = (int) inchain.get(0, col);
                    double rate = sn.rates.get(sourceStationIdx, classIdx);
                    if (!Double.isNaN(rate) && Double.isFinite(rate)) {
                        totalLambda += rate;
                    }
                }
                lambdaChain.set(0, c, totalLambda);
            }
        }

        // Create the new aggregated network
        Network chainModel = new Network(model.getName() + "_aggregated");

        // Map from original node index to new node
        Map<Integer, Node> nodeMap = new HashMap<>();
        Map<Integer, Station> stationMap = new HashMap<>();

        // First, identify which nodes to copy (skip auto-added ClassSwitch nodes)
        List<Integer> nodesToCopy = new ArrayList<>();
        for (int i = 0; i < model.getNodes().size(); i++) {
            Node node = model.getNodes().get(i);
            // Skip auto-added class switch nodes
            if (node instanceof ClassSwitch && ((ClassSwitch) node).autoAdded) {
                continue;
            }
            nodesToCopy.add(i);
        }

        // Create nodes in the new model
        for (int i : nodesToCopy) {
            Node node = model.getNodes().get(i);

            if (node instanceof Source) {
                nodeMap.put(i, new Source(chainModel, node.getName()));
            } else if (node instanceof Sink) {
                nodeMap.put(i, new Sink(chainModel, node.getName()));
            } else if (node instanceof Delay) {
                nodeMap.put(i, new Delay(chainModel, node.getName()));
            } else if (node instanceof jline.lang.nodes.Queue) {
                jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                jline.lang.nodes.Queue newQueue = new jline.lang.nodes.Queue(chainModel, node.getName(), queue.getSchedStrategy());
                if (!Double.isInfinite(queue.getNumberOfServers())) {
                    newQueue.setNumberOfServers(queue.getNumberOfServers());
                }
                nodeMap.put(i, newQueue);
            } else if (node instanceof Router) {
                nodeMap.put(i, new Router(chainModel, node.getName()));
            } else if (node instanceof ClassSwitch) {
                // User-defined ClassSwitch nodes should not exist in the
                // aggregated model since class switching is eliminated
                continue;
            } else {
                line_warning(mfilename(new Object() {}),
                    String.format("Node type %s not fully supported in chain aggregation.", node.getClass().getSimpleName()));
                continue;
            }

            // Store station mapping
            Node newNode = nodeMap.get(i);
            if (newNode instanceof Station) {
                int stationIdx = (int) sn.nodeToStation.get(i);
                if (stationIdx >= 0) {
                    stationMap.put(stationIdx, (Station) newNode);
                }
            }
        }

        // Create chain classes
        List<JobClass> chainClass = new ArrayList<>();
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);

            // Build chain class name from original class names
            String chainClassName;
            if (inchain.getNumCols() == 1) {
                int classIdx = (int) inchain.get(0, 0);
                chainClassName = sn.classnames.get(classIdx);
            } else {
                chainClassName = "Chain" + c;
            }
            if (!suffix.isEmpty()) {
                chainClassName = chainClassName + suffix;
            }

            if (isOpenChain[c]) {
                // Open chain
                chainClass.add(new OpenClass(chainModel, chainClassName));
            } else {
                // Closed chain - find reference station
                int refStationIdx = (int) refstatchain.get(c);
                Station refStation = stationMap.get(refStationIdx);
                if (refStation == null) {
                    line_error(mfilename(new Object() {}),
                        String.format("Reference station %d for chain %d not found in aggregated model.", refStationIdx, c));
                    return null;
                }
                chainClass.add(new ClosedClass(chainModel, chainClassName, Nchain.get(0, c), refStation));
            }
        }

        // Set arrival rates for open chains
        boolean hasOpenChains = false;
        for (int c = 0; c < C; c++) {
            if (isOpenChain[c]) {
                hasOpenChains = true;
                break;
            }
        }
        if (hasOpenChains) {
            Source chainSource = chainModel.getSource();
            for (int c = 0; c < C; c++) {
                if (isOpenChain[c] && lambdaChain.get(0, c) > 0) {
                    chainSource.setArrival(chainClass.get(c), Exp.fitMean(1.0 / lambdaChain.get(0, c)));
                }
            }
        }

        // Set service times at each station
        for (int i = 0; i < M; i++) {
            Station station = stationMap.get(i);
            if (station == null || station instanceof Source) {
                continue;
            }

            for (int c = 0; c < C; c++) {
                double stchain = STchain.get(i, c);
                if (stchain > 0 && Double.isFinite(stchain)) {
                    double scv = SCVchain.get(i, c);
                    if (Double.isNaN(scv) || scv <= 0) {
                        scv = 1.0; // Default to exponential
                    }

                    // Choose distribution based on SCV
                    Distribution dist;
                    if (Math.abs(scv - 1.0) < GlobalConstants.FineTol) {
                        // Exponential (SCV = 1)
                        dist = Exp.fitMean(stchain);
                    } else if (scv < 1.0) {
                        // SCV < 1: use Erlang or deterministic
                        if (scv < GlobalConstants.FineTol) {
                            dist = new Det(stchain);
                        } else {
                            // Erlang: SCV = 1/k, so k = 1/SCV
                            int k = (int) Math.round(1.0 / scv);
                            if (k < 1) {
                                k = 1;
                            }
                            dist = Erlang.fitMeanAndOrder(stchain, k);
                        }
                    } else {
                        // SCV > 1: use HyperExp
                        dist = HyperExp.fitMeanAndSCV(stchain, scv);
                    }

                    if (station instanceof ServiceStation) {
                        ((ServiceStation) station).setService(chainClass.get(c), dist);
                    }
                } else {
                    // No service for this chain at this station
                    if (station instanceof ServiceStation) {
                        ((ServiceStation) station).setService(chainClass.get(c), Disabled.getInstance());
                    }
                }
            }
        }

        // Build routing matrix from aggregated routing probabilities
        RoutingMatrix P = chainModel.initRoutingMatrix();
        int I_new = chainModel.getNodes().size();

        // Create mapping from original station to new model node index
        int[] stationToNewNode = new int[M];
        Arrays.fill(stationToNewNode, -1);
        for (int i = 0; i < M; i++) {
            int iNode = (int) sn.stationToNode.get(i);
            Node mappedNode = nodeMap.get(iNode);
            if (mappedNode != null) {
                // Find the index of mappedNode in chainModel.nodes
                for (int n = 0; n < I_new; n++) {
                    if (chainModel.getNodes().get(n) == mappedNode) {
                        stationToNewNode[i] = n;
                        break;
                    }
                }
            }
        }

        // For each chain, compute routing probabilities between stations
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);

            // Build routing for this chain
            for (int i = 0; i < M; i++) {
                int iNode = (int) sn.stationToNode.get(i);
                if (!nodeMap.containsKey(iNode) || stationToNewNode[i] < 0) {
                    continue;
                }

                // Skip source/sink in routing calculations for closed chains
                if (!isOpenChain[c] &&
                    (sn.nodetype.get(iNode) == NodeType.Source || sn.nodetype.get(iNode) == NodeType.Sink)) {
                    continue;
                }

                // Get stateful index for station i (sn.rt is indexed by stateful nodes)
                int isf_i = (int) sn.stationToStateful.get(i);

                for (int j = 0; j < M; j++) {
                    int jNode = (int) sn.stationToNode.get(j);
                    if (!nodeMap.containsKey(jNode) || stationToNewNode[j] < 0) {
                        continue;
                    }

                    // Get stateful index for station j
                    int isf_j = (int) sn.stationToStateful.get(j);

                    // Compute aggregated routing probability from station i to station j
                    double pij = 0.0;
                    for (int kCol = 0; kCol < inchain.getNumCols(); kCol++) {
                        int k = (int) inchain.get(0, kCol);
                        for (int sCol = 0; sCol < inchain.getNumCols(); sCol++) {
                            int s = (int) inchain.get(0, sCol);
                            // Use stateful indices to access sn.rt
                            int fromIdx = (isf_i) * K + k;
                            int toIdx = (isf_j) * K + s;
                            if (fromIdx < sn.rt.getNumRows() && toIdx < sn.rt.getNumCols()) {
                                double p_ks = sn.rt.get(fromIdx, toIdx);
                                if (alpha.get(i, k) > 0 && p_ks > 0) {
                                    pij += alpha.get(i, k) * p_ks;
                                }
                            }
                        }
                    }

                    if (pij > GlobalConstants.FineTol) {
                        int iNodeNew = stationToNewNode[i];
                        int jNodeNew = stationToNewNode[j];
                        P.set(chainClass.get(c), chainClass.get(c),
                              chainModel.getNodes().get(iNodeNew),
                              chainModel.getNodes().get(jNodeNew), pij);
                    }
                }
            }

            // Normalize routing probabilities (should be close to 1 already if indexing is correct)
            // Get routing matrix for this chain-to-chain and normalize rows
            Matrix Pcc = P.get(chainClass.get(c), chainClass.get(c));
            if (Pcc != null) {
                for (int iNew = 0; iNew < Pcc.getNumRows(); iNew++) {
                    double rowSum = 0.0;
                    for (int jNew = 0; jNew < Pcc.getNumCols(); jNew++) {
                        rowSum += Pcc.get(iNew, jNew);
                    }
                    if (rowSum > GlobalConstants.FineTol) {
                        if (Math.abs(rowSum - 1.0) > 0.01) {
                            line_warning(mfilename(new Object() {}),
                                String.format("Large normalization correction at node %d for chain %d: rowSum=%.4f", iNew, c, rowSum));
                        }
                        for (int jNew = 0; jNew < Pcc.getNumCols(); jNew++) {
                            Pcc.set(iNew, jNew, Pcc.get(iNew, jNew) / rowSum);
                        }
                    }
                }
            }
        }

        // Link the model with the routing matrix
        chainModel.link(P);

        // Build deaggregation information
        DeaggInfo deaggInfo = new DeaggInfo();
        deaggInfo.alpha = alpha;
        deaggInfo.Vchain = Vchain;
        deaggInfo.STchain = STchain;
        deaggInfo.Lchain = Lchain;
        deaggInfo.SCVchain = SCVchain;
        deaggInfo.Nchain = Nchain;
        deaggInfo.lambdaChain = lambdaChain;
        deaggInfo.isOpenChain = isOpenChain;
        deaggInfo.inchain = sn.inchain;
        deaggInfo.refstat = sn.refstat;
        deaggInfo.refstatchain = refstatchain;
        deaggInfo.originalSn = sn;
        deaggInfo.isAggregated = true;
        deaggInfo.nclasses = K;
        deaggInfo.nchains = C;

        return new AggregateChainResult(chainModel, alpha, deaggInfo);
    }

    // ========== Flow-Equivalent Server (FES) Aggregation ==========

    /**
     * Replace a station subset with a Flow-Equivalent Server (FES)
     *
     * This function replaces a subset of stations in a closed product-form
     * queueing network with a single Flow-Equivalent Server (FES). The FES has
     * Limited Joint Dependence (LJD) service rates where the rate for class-c
     * in state (n1,...,nK) equals the throughput of class-c in an isolated
     * subnetwork consisting only of the subset stations.
     *
     * @param model Closed product-form Network model
     * @param stationSubset List of Station objects to aggregate
     * @return FESResult containing the FES model, FES station, and deaggregation info
     */
    public static FESResult aggregateFES(Network model, List<Station> stationSubset) {
        return FESAggregator.aggregateFES(model, stationSubset);
    }

    /**
     * Replace a station subset with a Flow-Equivalent Server (FES) with options
     *
     * @param model Closed product-form Network model
     * @param stationSubset List of Station objects to aggregate
     * @param options FES aggregation options
     * @return FESResult containing the FES model, FES station, and deaggregation info
     */
    public static FESResult aggregateFES(Network model, List<Station> stationSubset, FESOptions options) {
        return FESAggregator.aggregateFES(model, stationSubset, options);
    }

    /**
     * Result of the FJ tag augmentation (fjtag).
     */
    public static class FJTagResult {
        public final Network fjmodel;
        public final NetworkStruct fjsn;
        public final Matrix fjclassmap; // (1,Kaug) original class of each auxiliary class, -1 for originals

        public FJTagResult(Network fjmodel, NetworkStruct fjsn, Matrix fjclassmap) {
            this.fjmodel = fjmodel;
            this.fjsn = fjsn;
            this.fjclassmap = fjclassmap;
        }
    }

    /**
     * Validate that a fork-join model is within the supported feature set of
     * the native CTMC/SSA fork-join implementation: closed chains only,
     * non-nested fork-join pairs, standard join strategy, integer
     * tasksPerLink. Mirrors matlab/src/api/fj/sn_fj_validate.m.
     */
    /**
     * Can the exact fork-join construction be asked for this model?
     *
     * <p>The fork-join model class {@link #fjValidate} admits, asked as a
     * predicate rather than thrown. {@code SolverCTMC.supportsModelMethod} and
     * {@code SolverSSA.supportsModelMethod} call it so that a caller (help,
     * findSolver, SolverAUTO) sees the verdict before paying for a run, and BOTH
     * analyzers reach the SAME rules through {@link #fjtag}. The sentence the
     * validator throws names "the native CTMC/SSA fork-join implementation",
     * which is why this predicate lives beside it rather than in either
     * solver.</p>
     *
     * <p>IT WRAPS THE VALIDATOR RATHER THAN RESTATING IT, and that is the point:
     * the rules are eight and they move (pairing, join strategy,
     * tasks-per-link, branch probability, open classes through a fork), so a
     * second copy would be a second thing to keep in step. There is exactly one
     * body of rules and two ways in -- one that throws, for the run, and this
     * one, which answers.</p>
     *
     * <p>WHAT IT REFUSES AND WHY THE ANALYZER IS RIGHT TO. The fork-join PAIRING
     * is a declaration carried by the Join ({@code new Join(model, name, fork)}
     * in all four codebases), not a derivation from the routing: a nested model
     * such as fj_basic_nesting has two forks and two joins whose pairing the
     * routing alone does not determine. So a Join built without naming its fork
     * leaves {@code sn.fj} empty, and "Fork nodes without a matched Join" is the
     * honest answer to a model that declares none.</p>
     *
     * @param sn the network structure
     * @return empty string when the fork-join construction may run, else the refusal
     */
    public static String fjSupportsReason(NetworkStruct sn) {
        boolean anyFJ = false;
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Fork || nt == NodeType.Join) {
                anyFJ = true;
                break;
            }
        }
        if (!anyFJ) {
            return "";
        }
        try {
            fjValidate(sn);
        } catch (RuntimeException e) {
            // the validator's own refusal, turned into an answer
            return e.getMessage() == null ? "This fork-join model is not supported." : e.getMessage();
        }
        return "";
    }

    public static void fjValidate(NetworkStruct sn) {
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);
        for (int f = 0; f < sn.nodetype.size(); f++) {
            if (sn.nodetype.get(f) != NodeType.Fork) continue;
            List<Integer> joins = new ArrayList<Integer>();
            for (int j = 0; j < sn.fj.getNumCols(); j++) {
                if (sn.fj.get(f, j) > 0) joins.add(j);
            }
            if (joins.isEmpty()) {
                throw new RuntimeException("Fork nodes without a matched Join are not supported by the native CTMC/SSA fork-join implementation.");
            }
            if (joins.size() > 1) {
                throw new RuntimeException("Multiple Join nodes per Fork are not supported by the native CTMC/SSA fork-join implementation.");
            }
            for (int r = 0; r < sn.nclasses; r++) {
                if (Vnodes.get(f, r) <= 0) continue;
                int c = -1;
                for (int ci = 0; ci < sn.chains.getNumRows(); ci++) {
                    if (sn.chains.get(ci, r) > 0) { c = ci; break; }
                }
                double chainPop = 0;
                if (c >= 0) {
                    for (int kk = 0; kk < sn.nclasses; kk++) {
                        if (sn.chains.get(c, kk) > 0) chainPop += sn.njobs.get(kk);
                    }
                } else {
                    chainPop = sn.njobs.get(r);
                }
                if (Double.isInfinite(sn.njobs.get(r)) || Double.isInfinite(chainPop)) {
                    throw new RuntimeException("Open classes routed through a Fork are not supported by the native CTMC/SSA fork-join implementation.");
                }
            }
        }
        for (int j = 0; j < sn.nodetype.size(); j++) {
            if (sn.nodetype.get(j) != NodeType.Join) continue;
            boolean matched = false;
            for (int f = 0; f < sn.fj.getNumRows(); f++) {
                if (sn.fj.get(f, j) > 0) { matched = true; break; }
            }
            if (!matched) {
                throw new RuntimeException("Join nodes without a matched Fork are not supported by the native CTMC/SSA fork-join implementation.");
            }
            Object npObj = sn.nodeparam.get(sn.nodes.get(j));
            if (npObj instanceof jline.lang.nodeparam.JoinNodeParam) {
                jline.lang.nodeparam.JoinNodeParam jnp = (jline.lang.nodeparam.JoinNodeParam) npObj;
                if (jnp.joinStrategy != null) {
                    for (jline.lang.constant.JoinStrategy js : jnp.joinStrategy.values()) {
                        if (js != null && js != jline.lang.constant.JoinStrategy.STD) {
                            throw new RuntimeException("Only JoinStrategy.STD is supported by the native CTMC/SSA fork-join implementation.");
                        }
                    }
                }
            }
        }
    }

    /**
     * Fold the auxiliary-class metric columns computed on an FJ tag-augmented
     * struct back into the original classes: queue lengths, utilizations and
     * throughputs of the sibling classes are exact aggregates of the class
     * they were forked from; response times are recomputed by Little's law
     * after folding. Mirrors matlab/src/api/fj/sn_fj_foldback.m. The input
     * matrices are modified in place up to column Korig; callers should
     * extract the first Korig columns afterwards.
     */
    public static void fjFoldback(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix fjclassmap, int Korig) {
        int Kaug = fjclassmap.getNumCols();
        for (int a = 0; a < Kaug; a++) {
            int r = (int) fjclassmap.get(0, a);
            if (r >= 0) {
                for (int ist = 0; ist < QN.getNumRows(); ist++) {
                    QN.set(ist, r, QN.get(ist, r) + QN.get(ist, a));
                    UN.set(ist, r, UN.get(ist, r) + UN.get(ist, a));
                    TN.set(ist, r, TN.get(ist, r) + TN.get(ist, a));
                }
            }
        }
        for (int r = 0; r < Korig; r++) {
            for (int ist = 0; ist < QN.getNumRows(); ist++) {
                if (TN.get(ist, r) > 0) {
                    RN.set(ist, r, QN.get(ist, r) / TN.get(ist, r));
                } else {
                    RN.set(ist, r, 0.0);
                }
            }
        }
    }

    /**
     * Build a tag-augmented copy of a closed fork-join model for exact native
     * analysis by SolverCTMC/SolverSSA. For each (fork f, class r) pair with
     * matched join j, branch b and tag t (one per chain population slot), an
     * auxiliary transient closed class with population 0 is created. The tag
     * identifies the origin job: the fork firing (AfterFJEvent) emits weight
     * (tasksPerLink) siblings per branch in the tag's auxiliary classes using
     * the lowest free tag, and the join fires only when all siblings OF THE
     * SAME TAG are buffered, releasing one class-r job (AfterEventJoin).
     * Mirrors matlab/src/io/@ModelAdapter/fjtag.m.
     */
    public static FJTagResult fjtag(Network model) {
        NetworkStruct sn = model.getStruct(true);
        fjValidate(sn);

        int K = sn.nclasses;
        int I = sn.nnodes;
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);

        Network fjmodel;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(model);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            fjmodel = (Network) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("fjtag: could not create a copy of the model.");
        }
        fjmodel.setAllowReplace(true);
        fjmodel.isFJAugmented = true;
        NetworkStruct copySn = fjmodel.getStruct(false);
        Map<JobClass, Map<JobClass, Matrix>> P = copySn.rtorig;
        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
        int Ip = 0;
        for (Map<JobClass, Matrix> row : P.values()) {
            for (Matrix Pblock : row.values()) {
                if (Pblock != null && Pblock.getNumRows() > Ip) Ip = Pblock.getNumRows();
            }
        }
        fjmodel.resetNetwork(true);
        fjmodel.resetStruct();

        // replace each Fork with a StatefulFork (positional replacement by name)
        Map<Integer, Double> forkTasksPerLink = new HashMap<Integer, Double>();
        for (int f = 0; f < sn.nodetype.size(); f++) {
            if (sn.nodetype.get(f) != NodeType.Fork) continue;
            double tpl = ((Forker) model.getNodes().get(f).getOutput()).tasksPerLink;
            forkTasksPerLink.put(f, tpl);
            Node oldNode = fjmodel.getNodes().get(f);
            StatefulFork newFork = new StatefulFork(fjmodel, oldNode.getName());
            ((Forker) newFork.getOutput()).tasksPerLink = tpl;
        }

        Matrix njobs = sn.njobs;
        List<int[]> forkinfoFJR = new ArrayList<int[]>();          // {f, j, r, w, T}
        List<int[]> forkinfoHeads = new ArrayList<int[]>();        // branch heads per row
        List<List<List<Integer>>> forkinfoSets = new ArrayList<List<List<Integer>>>(); // branch sets per row
        List<Matrix> forkinfoAux = new ArrayList<Matrix>();        // B x T aux class indices per row
        List<Integer> auxOfClass = new ArrayList<Integer>();       // fjclassmap entries for aux classes (by creation order)

        for (int f = 0; f < sn.nodetype.size(); f++) {
            if (sn.nodetype.get(f) != NodeType.Fork) continue;
            int j = -1;
            for (int jj = 0; jj < sn.fj.getNumCols(); jj++) {
                if (sn.fj.get(f, jj) > 0) { j = jj; break; }
            }
            double tplD = forkTasksPerLink.get(f);
            if (tplD != Math.rint(tplD)) {
                throw new RuntimeException("Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.");
            }
            int w = Math.max(1, (int) tplD);
            for (int r = 0; r < K; r++) {
                if (Vnodes.get(f, r) <= 0) continue;
                // branch heads: nodes receiving class r directly from the fork
                List<Integer> heads = new ArrayList<Integer>();
                for (int jnd = 0; jnd < I; jnd++) {
                    if (sn.rtnodes.get(f * K + r, jnd * K + r) > 0) heads.add(jnd);
                }
                int B = heads.size();
                if (B < 2) {
                    throw new RuntimeException("Degenerate forks with a single output link are not supported by the native CTMC/SSA fork-join implementation.");
                }
                // branch discovery: class-r BFS closure from each head up to the join
                List<List<Integer>> branchsets = new ArrayList<List<Integer>>();
                for (int b = 0; b < B; b++) {
                    List<Integer> visitset = new ArrayList<Integer>();
                    visitset.add(heads.get(b));
                    List<Integer> frontier = new ArrayList<Integer>();
                    frontier.add(heads.get(b));
                    while (!frontier.isEmpty()) {
                        int cn = frontier.remove(0);
                        if (sn.nodetype.get(cn) == NodeType.Fork) {
                            throw new RuntimeException("Nested fork-join is not supported by the native CTMC/SSA fork-join implementation.");
                        }
                        if (sn.nodetype.get(cn) == NodeType.Join && cn != j) {
                            throw new RuntimeException("Overlapping fork-join pairs are not supported by the native CTMC/SSA fork-join implementation.");
                        }
                        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
                        if (cn >= Ip) {
                            throw new RuntimeException("Class switching between fork and join is not supported by the native CTMC/SSA fork-join implementation.");
                        }
                        for (int jnd = 0; jnd < I; jnd++) {
                            for (int s = 0; s < K; s++) {
                                if (sn.rtnodes.get(cn * K + r, jnd * K + s) > 0) {
                                    if (s != r) {
                                        throw new RuntimeException("Class switching between fork and join is not supported by the native CTMC/SSA fork-join implementation.");
                                    }
                                    if (jnd != j && !visitset.contains(jnd)) {
                                        visitset.add(jnd);
                                        frontier.add(jnd);
                                    }
                                }
                            }
                        }
                    }
                    // trap check: every branch node must reach the join within the branch
                    List<Integer> canreach = new ArrayList<Integer>();
                    canreach.add(j);
                    boolean changed = true;
                    while (changed) {
                        changed = false;
                        for (int cn : visitset) {
                            if (!canreach.contains(cn)) {
                                for (int jnd : new ArrayList<Integer>(canreach)) {
                                    if (sn.rtnodes.get(cn * K + r, jnd * K + r) > 0) {
                                        canreach.add(cn);
                                        changed = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    for (int cn : visitset) {
                        if (!canreach.contains(cn)) {
                            throw new RuntimeException("Fork branches from which the Join is unreachable are not supported by the native CTMC/SSA fork-join implementation.");
                        }
                    }
                    branchsets.add(visitset);
                }
                // tag pool size = chain population (class switching outside the
                // fork-join section can concentrate the whole chain in class r)
                int c = -1;
                for (int ci = 0; ci < sn.chains.getNumRows(); ci++) {
                    if (sn.chains.get(ci, r) > 0) { c = ci; break; }
                }
                double chainPop = 0;
                for (int kk = 0; kk < K; kk++) {
                    if (c >= 0 && sn.chains.get(c, kk) > 0) chainPop += njobs.get(kk);
                }
                int T = (int) Math.rint(chainPop);
                Matrix auxmatrix = new Matrix(B, T);
                Station refstatStation = fjmodel.getStations().get((int) sn.refstat.get(r, 0));
                int prio = (sn.classprio != null && !sn.classprio.isEmpty()) ? (int) sn.classprio.get(r) : 0;
                jline.lang.nodes.Join joinNode = (jline.lang.nodes.Join) fjmodel.getNodes().get(j);
                for (int t = 0; t < T; t++) {
                    for (int b = 0; b < B; b++) {
                        String auxname = sn.classnames.get(r) + "_f" + (f + 1) + "_b" + (b + 1) + "_t" + (t + 1);
                        ClosedClass auxclass = new ClosedClass(fjmodel, auxname, 0.0, refstatStation, prio);
                        int a = auxclass.getIndex() - 1; // 0-based
                        auxmatrix.set(b, t, a);
                        while (auxOfClass.size() <= a - K) auxOfClass.add(-1);
                        auxOfClass.set(a - K, r);
                        // sibling service on the branch copies the original class
                        for (int cn : branchsets.get(b)) {
                            if (sn.isstation.get(cn, 0) == 1 && sn.nodetype.get(cn) != NodeType.Join) {
                                Distribution distributionCopy;
                                try {
                                    ByteArrayOutputStream bos = new ByteArrayOutputStream();
                                    ObjectOutputStream out = new ObjectOutputStream(bos);
                                    out.writeObject(((jline.lang.nodes.Queue) model.getNodes().get(cn)).getService(model.getClasses().get(r)));
                                    ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                                    ObjectInputStream in = new ObjectInputStream(bis);
                                    distributionCopy = (Distribution) in.readObject();
                                } catch (IOException | ClassNotFoundException e) {
                                    throw new RuntimeException("fjtag: could not copy the service distribution of the original class.");
                                }
                                ((jline.lang.nodes.Queue) fjmodel.getNodes().get(cn)).setService(auxclass, distributionCopy);
                            }
                        }
                        // register the auxiliary class at the join input section
                        joinNode.setStrategy(auxclass, jline.lang.constant.JoinStrategy.STD);
                        joinNode.setRequired(auxclass, -1);
                        // sibling routing: copy the class-r branch routing; the
                        // auxiliary class terminates at the join (no outgoing row)
                        JobClass rClassNew = fjmodel.getJobClasses().get(r);
                        Matrix Prr = (P.containsKey(rClassNew) && P.get(rClassNew).containsKey(rClassNew))
                                ? P.get(rClassNew).get(rClassNew) : null;
                        if (Prr != null) {
                            Matrix auxP = new Matrix(Prr.getNumRows(), Prr.getNumCols());
                            for (int cn : branchsets.get(b)) {
                                for (int jnd = 0; jnd < Prr.getNumCols(); jnd++) {
                                    if (Prr.get(cn, jnd) != 0) auxP.set(cn, jnd, Prr.get(cn, jnd));
                                }
                            }
                            if (!P.containsKey(auxclass)) P.put(auxclass, new HashMap<JobClass, Matrix>());
                            P.get(auxclass).put(auxclass, auxP);
                        }
                    }
                }
                int[] fjr = new int[]{f, j, r, w, T};
                forkinfoFJR.add(fjr);
                int[] headsArr = new int[B];
                for (int b = 0; b < B; b++) headsArr[b] = heads.get(b);
                forkinfoHeads.add(headsArr);
                forkinfoSets.add(branchsets);
                forkinfoAux.add(auxmatrix);
            }
        }

        int Kaug = fjmodel.getJobClasses().size();

        // rebuild the routing matrix and relink
        RoutingMatrix routingMatrix = new RoutingMatrix(fjmodel, fjmodel.getClasses(), fjmodel.getNodes());
        int numNodesFj = fjmodel.getNumberOfNodes();
        for (JobClass r : P.keySet()) {
            for (JobClass s : P.get(r).keySet()) {
                Matrix Prs = P.get(r).get(s);
                int maxRows = Math.min(Prs.getNumRows(), numNodesFj);
                int maxCols = Math.min(Prs.getNumCols(), numNodesFj);
                for (int i = 0; i < maxRows; i++) {
                    for (int jj = 0; jj < maxCols; jj++) {
                        if (Prs.get(i, jj) != 0) {
                            routingMatrix.set(r, s, fjmodel.getNodes().get(i), fjmodel.getNodes().get(jj), Prs.get(i, jj));
                        }
                    }
                }
            }
        }
        fjmodel.relink(routingMatrix);

        // see _kb/04-networkstruct.md (refreshStruct.m: fjclassmap pairing) for rationale
        Matrix csm = fjmodel.getCsMatrix();
        Matrix csmNew = new Matrix(Kaug, Kaug);
        csmNew.zero();
        if (csm != null && !csm.isEmpty()) {
            for (int i = 0; i < Math.min(Kaug, csm.getNumRows()); i++) {
                for (int jj = 0; jj < Math.min(Kaug, csm.getNumCols()); jj++) {
                    csmNew.set(i, jj, csm.get(i, jj));
                }
            }
        }
        for (int a = 0; a < Kaug; a++) {
            csmNew.set(a, a, 1.0);
            if (a >= K) {
                int r = auxOfClass.get(a - K);
                if (r >= 0) {
                    csmNew.set(r, a, 1.0);
                    csmNew.set(a, r, 1.0);
                }
            }
        }
        fjmodel.setCsMatrix(csmNew);

        // re-initialize the default state: the copy inherits the original
        // model's states, which have pre-augmentation class widths
        for (Node nodeI : fjmodel.getNodes()) {
            if (nodeI instanceof StatefulNode) {
                ((StatefulNode) nodeI).setState(new Matrix(0, 0));
            }
        }
        fjmodel.resetStruct();
        try {
            fjmodel.initDefault();
        } catch (Exception e) {
            throw new RuntimeException("fjtag: default state initialization failed on the tag-augmented copy.", e);
        }
        NetworkStruct fjsn = fjmodel.getStruct(true);

        // ---- sn post-edits ----
        Matrix fjclassmap = new Matrix(1, Kaug);
        fjclassmap.fill(-1);
        for (int a = K; a < Kaug; a++) {
            fjclassmap.set(0, a, auxOfClass.get(a - K));
        }

        // see _kb/04-networkstruct.md (Python io/model_adapter.py: MMT cache-reuse invariants, JAR notes) for rationale
        for (int r = 0; r < K; r++) {
            int corig = -1, cnew = -1;
            for (int ci = 0; ci < sn.chains.getNumRows(); ci++) {
                if (sn.chains.get(ci, r) > 0) { corig = ci; break; }
            }
            for (int ci = 0; ci < fjsn.chains.getNumRows(); ci++) {
                if (fjsn.chains.get(ci, r) > 0) { cnew = ci; break; }
            }
            if (corig < 0 || cnew < 0) continue;
            Matrix vOld = sn.visits.get(corig);
            Matrix vNew = fjsn.visits.get(cnew);
            for (int isfnew = 0; isfnew < fjsn.nstateful && isfnew < vNew.getNumRows(); isfnew++) {
                int ind = (int) fjsn.statefulToNode.get(isfnew);
                if (ind < sn.nnodes && sn.isstateful.get(ind, 0) == 1) {
                    int isfold = (int) sn.nodeToStateful.get(ind);
                    if (isfold < vOld.getNumRows()) {
                        vNew.set(isfnew, r, vOld.get(isfold, r));
                    }
                } else { // stateful Fork: nonzero marker for capacity gating
                    vNew.set(isfnew, r, Vnodes.get(ind, r) > 0 ? 1.0 : 0.0);
                }
            }
            Matrix nvOld = sn.nodevisits.get(corig);
            Matrix nvNew = fjsn.nodevisits.get(cnew);
            for (int ind = 0; ind < Math.min(nvOld.getNumRows(), nvNew.getNumRows()); ind++) {
                nvNew.set(ind, r, nvOld.get(ind, r));
            }
        }

        // auxiliary-class visits (zero-versus-nonzero capacity gate), classcap,
        // nodeparam fj blocks and the fjsync entries
        Map<Integer, FJSync> fjsync = new LinkedHashMap<Integer, FJSync>();
        for (int row = 0; row < forkinfoFJR.size(); row++) {
            int f = forkinfoFJR.get(row)[0];
            int j = forkinfoFJR.get(row)[1];
            int r = forkinfoFJR.get(row)[2];
            int w = forkinfoFJR.get(row)[3];
            int T = forkinfoFJR.get(row)[4];
            int[] heads = forkinfoHeads.get(row);
            List<List<Integer>> branchsets = forkinfoSets.get(row);
            Matrix auxmatrix = forkinfoAux.get(row);
            int B = heads.length;
            int cnew = -1;
            for (int ci = 0; ci < fjsn.chains.getNumRows(); ci++) {
                if (fjsn.chains.get(ci, r) > 0) { cnew = ci; break; }
            }
            for (int b = 0; b < B; b++) {
                List<Integer> support = new ArrayList<Integer>(branchsets.get(b));
                support.add(j);
                for (int t = 0; t < T; t++) {
                    int a = (int) auxmatrix.get(b, t);
                    Matrix vNew = fjsn.visits.get(cnew);
                    Matrix nvNew = fjsn.nodevisits.get(cnew);
                    for (int isf = 0; isf < vNew.getNumRows(); isf++) vNew.set(isf, a, 0.0);
                    for (int ind = 0; ind < nvNew.getNumRows(); ind++) nvNew.set(ind, a, 0.0);
                    for (int ist = 0; ist < fjsn.classcap.getNumRows(); ist++) fjsn.classcap.set(ist, a, 0.0);
                    for (int cn : support) {
                        if (cn < nvNew.getNumRows()) nvNew.set(cn, a, 1.0);
                        if (fjsn.isstateful.get(cn, 0) == 1) {
                            int isf = (int) fjsn.nodeToStateful.get(cn);
                            if (isf < vNew.getNumRows()) vNew.set(isf, a, 1.0);
                        }
                        if (fjsn.isstation.get(cn, 0) == 1) {
                            int ist = (int) fjsn.nodeToStation.get(cn);
                            // each auxiliary class holds at most tasksPerLink
                            // siblings network-wide (STD join, one tag at a time)
                            fjsn.classcap.set(ist, a, w);
                        }
                    }
                }
            }
            // nodeparam fj blocks read by AfterEventJoin/AfterFJEvent
            Node forkNode = fjsn.nodes.get(f);
            Object fnpObj = fjsn.nodeparam.get(forkNode);
            jline.lang.nodeparam.ForkNodeParam fnp;
            if (fnpObj instanceof jline.lang.nodeparam.ForkNodeParam) {
                fnp = (jline.lang.nodeparam.ForkNodeParam) fnpObj;
            } else {
                fnp = new jline.lang.nodeparam.ForkNodeParam();
                fnp.fanOut = w;
                fjsn.nodeparam.put(forkNode, fnp);
            }
            if (fnp.fjClasses == null) {
                fnp.fjClasses = new ArrayList<Integer>();
                fnp.fjJoins = new ArrayList<Integer>();
                fnp.fjAuxmatrix = new HashMap<Integer, Matrix>();
                fnp.fjBranchheads = new HashMap<Integer, int[]>();
            }
            fnp.fjClasses.add(r);
            fnp.fjJoins.add(j);
            fnp.fjAuxmatrix.put(r, auxmatrix);
            fnp.fjBranchheads.put(r, heads);
            Node joinNodeSn = fjsn.nodes.get(j);
            Object jnpObj = fjsn.nodeparam.get(joinNodeSn);
            jline.lang.nodeparam.JoinNodeParam jnp;
            if (jnpObj instanceof jline.lang.nodeparam.JoinNodeParam) {
                jnp = (jline.lang.nodeparam.JoinNodeParam) jnpObj;
            } else {
                jnp = new jline.lang.nodeparam.JoinNodeParam();
                fjsn.nodeparam.put(joinNodeSn, jnp);
            }
            jnp.fjFork = f;
            if (jnp.fjOrigclasses == null) {
                jnp.fjOrigclasses = new ArrayList<Integer>();
                jnp.fjAuxmatrix = new HashMap<Integer, Matrix>();
                jnp.fjRequired = new HashMap<Integer, Matrix>();
            }
            jnp.fjOrigclasses.add(r);
            jnp.fjAuxmatrix.put(r, auxmatrix);
            Matrix req = new Matrix(B, 1);
            req.fill(w); // STD, tasksPerLink siblings per branch; slot for PARTIAL/fanIn
            jnp.fjRequired.put(r, req);
            // fork firing synchronizations: one entry per (fork, class, tag)
            for (int t = 0; t < T; t++) {
                FJSync entry = new FJSync();
                entry.active = new Event(jline.lang.constant.EventType.FIRE, f, r);
                entry.fork = f;
                entry.join = j;
                entry.jobclass = r;
                entry.tag = t;
                entry.branchheads = heads;
                int[] auxcl = new int[B];
                for (int b = 0; b < B; b++) auxcl[b] = (int) auxmatrix.get(b, t);
                entry.auxclasses = auxcl;
                entry.auxall = auxmatrix;
                entry.weight = w;
                entry.prob = 1.0;
                fjsync.put(fjsync.size(), entry);
            }
        }

        fjsn.fjsync = fjsync;
        fjsn.fjclassmap = fjclassmap;
        fjsn.isfjaugmented = true;

        return new FJTagResult(fjmodel, fjsn, fjclassmap);
    }
}