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

import static jline.io.InputOutputKt.*;
import static jline.api.sn.SnGetDemandsChainKt.snGetDemandsChain;

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
    public static TaggedChainResult tagChain(Network model, Chain chain, JobClass jobclass, String suffix) {
        if (suffix == null || suffix.isEmpty()) {
            suffix = ".tagged";
        }
        if (jobclass == null && !chain.getClasses().isEmpty()) {
            jobclass = chain.getClasses().get(0);
        }
        
        // Create a copy of the model
        Network taggedModel = model.copy();
        
        // For simplicity, create a new tagged job class based on the original
        // This is a simplified implementation that focuses on the core functionality
        JobClass taggedJob = null;
        
        try {
            // Find the original job class in the tagged model
            JobClass originalInTagged = null;
            for (JobClass cls : taggedModel.getClasses()) {
                if (cls.getName().equals(jobclass.getName())) {
                    originalInTagged = cls;
                    break;
                }
            }
            
            if (originalInTagged != null) {
                // Create a tagged version of the job class
                if (originalInTagged instanceof ClosedClass) {
                    ClosedClass original = (ClosedClass) originalInTagged;
                    ClosedClass tagged = new ClosedClass(taggedModel, 
                        original.getName() + suffix, 
                        1, // Tagged job has population 1
                        original.getReferenceStation(), 
                        original.getPriority());
                    
                    taggedModel.addJobClass(tagged);
                    taggedJob = tagged;
                    
                    // Reduce original population by 1
                    if (original.getPopulation() > 1) {
                        original.setPopulation(original.getPopulation() - 1);
                    }
                } else if (originalInTagged instanceof OpenClass) {
                    // For open classes, create a closed tagged version with population 1
                    OpenClass original = (OpenClass) originalInTagged;
                    ClosedClass tagged = new ClosedClass(taggedModel,
                        original.getName() + suffix,
                        1, // Tagged job has population 1
                        null, // No specific reference station
                        original.getPriority());
                    
                    taggedModel.addJobClass(tagged);
                    taggedJob = tagged;
                }
            }
            
            // Reset and refresh the model
            taggedModel.reset(true);
            taggedModel.refreshStruct(true);
            
        } catch (Exception e) {
            // If tagging fails, return a simple copy with the original job
            taggedJob = jobclass;
        }
        
        return new TaggedChainResult(taggedModel, taggedJob);
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
        for (int i = 0; i < P.getNumCols(); i++) {
            if (P.get(startNode, i) == 0) {
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
                        TN, 0, fjclassmap, fjforkmap, nonfjmodel);
                Matrix lambdai = Matrix.ones(paths.getNumRows(), paths.getNumCols()).elementDiv(paths);
                double d0 = 0;
                int parallel_branches = paths.length();
                for (int pow = 0; pow < parallel_branches; pow++) {
                    Matrix nk = Maths.nCk(lambdai, pow + 1);
                    nk = nk.sumRows();
                    double currentSum = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
                    d0 += FastMath.pow(-1, pow) * currentSum;
                }
                for (int cls : toMerge) {
                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getJobClasses().get(cls),
                            Exp.fitMean(d0 - paths.elementSum() / paths.length()));
                }
                toMerge.remove(toMerge.size() - 1);
                ri = ri.concatCols(findPaths(sn, P, joinIdx, endNode, r,
                        toMerge, QN, TN, currentTime + d0, fjclassmap,
                        fjforkmap, nonfjmodel));
            } else {
                ri = ri.concatCols(findPaths(sn, P, i, endNode, r, toMerge,
                        QN, TN, currentTime + qLen / tput, fjclassmap, fjforkmap,
                        nonfjmodel));
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
        for (int transition = 0; transition < P.getNumCols(); transition++) {
            if (P.get(curClass * orignodes + curNode, transition) == 0) {
                continue;
            }
            ArrayList<Integer> curMerge = new ArrayList<>(toMerge);
            int nextClass = (int) FastMath.floor((transition) / (double) orignodes);
            int nextNode = transition - (nextClass) * orignodes;
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
                        TN, 0, fjclassmap, fjforkmap, nonfjmodel);
                Matrix lambdai = Matrix.ones(paths.getNumRows(), paths.getNumCols()).elementDiv(paths);
                double d0 = 0;
                int parallel_branches = paths.length();
                for (int pow = 0; pow < parallel_branches; pow++) {
                    Matrix nk = Maths.nCk(lambdai, pow + 1);
                    nk = nk.sumRows();
                    double currentSum = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
                    d0 += FastMath.pow(-1, pow) * currentSum;
                }
                // Match MATLAB: loop over [curMerge, s] to include inner auxiliary class
                ArrayList<Integer> mergeWithS = new ArrayList<>(curMerge);
                mergeWithS.add(s);
                for (int cls : mergeWithS) {
                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getJobClasses().get(cls),
                            Exp.fitMean(d0 - paths.elementSum() / paths.length()));
                }
                toMerge.remove(toMerge.size() - 1);
                ri = ri.concatCols(findPathsCS(sn, P, joinIdx, endNode, nextClass,
                        curMerge, QN, TN, currentTime + d0, fjclassmap,
                        fjforkmap, nonfjmodel));
            } else {
                ri = ri.concatCols(findPathsCS(sn, P, nextNode, endNode, nextClass, curMerge,
                        QN, TN, currentTime + qLen / tput, fjclassmap, fjforkmap,
                        nonfjmodel));
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
            for (JobClass r : P.keySet()) {
                List<OutputStrategy> outputStrategies = nonfjmodel.getNodes().get(f).getOutputStrategies();
                int parallelBranches = 0;
                for (OutputStrategy o : outputStrategies) {
                    if (o.getJobClass() == r && o.getDestination() != null) {
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
                    if (((Forker) model.getNodes().get(f).getOutput()).tasksPerLink > 1) {
                        line_warning(mfilename(new Object() {
                        }), "There are no synchronisation delays implemented in MMT for multiple tasks per link. Results may be inaccurate.");
                    }
                    int fanoutValue = (int) (origfanout.get(f, r) * ((Forker) model.getNodes().get(f).getOutput()).tasksPerLink);
                    fanout.put(oclass.get(oclass.size() - 1).getIndex() - 1, fanoutValue);
                    allAuxClassIndices.add(oclass.get(oclass.size() - 1).getIndex() - 1); // 0-based
                    // Use tolerance for floating point comparison (values < FineTol are effectively zero)
                    // Also disable classes with zero fanout: they pass through the fork via
                    // class-switching but aren't actually forked.
                    if (origfanout.get(f, r) == 0 || Math.abs(sn.nodevisits.get(fc).get(f, r)) < GlobalConstants.FineTol) {
                        source.setArrival(oclass.get(oclass.size() - 1), Disabled.getInstance());
                    } else {
                        source.setArrival(oclass.get(oclass.size() - 1), new Exp(forkLambda.get(r)));
                    }
                    // joins are now Delays, let us set their service time
                    for (int i = 0; i < sn.nnodes; i++) {
                        if (sn.isstation.get(i) != 0) {
                            switch (sn.nodetype.get(i)) {
                                case Join:
                                    ((Delay) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), Immediate.getInstance());
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
                                        line_error(mfilename(new Object[]{}), "Could not copy the distribution of the original class in FJT");
                                    }
                                    ((jline.lang.nodes.Queue) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), distributionCopy);
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
                                new Matrix(P.get(nonfjmodel.getJobClasses().get(r)).get(nonfjmodel.getJobClasses().get(s))));
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
                        Matrix Prs = P.get(oclass.get(rIdx)).get(oclass.get(sIdx));
                        for (int i = 0; i < Prs.getNumCols(); i++) {
                            Prs.set(source.getNodeIndex(), i, 0.0);
                            if (joinIdx >= 0)
                                Prs.set(joinIdx, i, 0.0);
                        }
                        sIdx++;
                    }
                    if (origfanout.get(f, r) > 0) {
                        P.get(oclass.get(rIdx)).get(oclass.get(rIdx)).set(source.getNodeIndex(), f, 1);
                        if (joinIdx >= 0) {
                            P.get(oclass.get(rIdx)).get(oclass.get(rIdx)).set(joinIdx, sink.getNodeIndex(), 1);
                        }
                    }
                    rIdx++;
                }

                // Check if all classes in this chain have non-zero fanout at this fork.
                // BFS scope clearing only applies when all classes are actually forked;
                // when some classes pass through the fork via class-switching (origfanout=0),
                // clearing would disconnect their aux chain and break MVA convergence.
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
                    // Determine fork-join scope via class-aware BFS from Fork, stopping
                    // at Join. Tracks (node, class) pairs so that shared stations (e.g.
                    // processors serving both fork-branch and return-path classes) are
                    // correctly handled: only fork-branch class routing is followed.
                    int pnnodes = P.get(nonfjmodel.getJobClasses().get(inchain.get(0)))
                            .get(nonfjmodel.getJobClasses().get(inchain.get(0))).getNumRows();
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
        // Fix spurious routing for aux classes at non-scope nodes and CS nodes.
        // relink() only calls setProbRouting for non-zero P entries, so nodes
        // where P was zeroed retain default RAND routing that inherits physical
        // connections from original classes. CS nodes created by link() for
        // cross-class routing also get default RAND for aux classes. Override
        // all non-PROB routing for aux classes to DISABLED.
        for (int nd = 0; nd < nonfjmodel.getNodes().size(); nd++) {
            List<OutputStrategy> os = nonfjmodel.getNodes().get(nd).getOutputStrategies();
            for (int ci : allAuxClassIndices) {
                if (ci < os.size()) {
                    if (os.get(ci).getRoutingStrategy() != RoutingStrategy.PROB) {
                        os.get(ci).setRoutingStrategy(RoutingStrategy.DISABLED);
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
        return new Ret.FJApprox(nonfjmodel, fjclassmapMatrix, fjforkmapMatrix, null, fanout);
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
}