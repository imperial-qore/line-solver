package jline.api;

import jline.lang.*;
import jline.lang.constant.NodeType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.distributions.Disabled;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;
import jline.lang.nodes.Queue;
import jline.lang.nodes.*;
import jline.lang.sections.Forker;
import jline.util.Maths;
import jline.util.Matrix;

import java.io.*;
import java.util.*;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * API for the methods used to deal with fork-join (FJ) systems in the MVA solver.
 */
public class FJ {
    /**
     * Heidelberger-Trivedi fork-join queueing network transformation.
     * Transforms the queueing network containing a FJ subsystem into a queueing network without one.
     * Fork nodes changed to Router nodes. Join nodes changed to Delay nodes.
     * One artificial class is created for each parallel branch and for each class.
     * Another delay is added to model the sojourn time of the original classes.
     * --
     * This approach is derived by PHILIP HEIDELBERGER and KISHOR S. TRIVEDI in
     * "Analytic Queueing Models for Programs with Internal Concurrency"
     * @param model - the original network
     * @return - queueing network with no FJ system, the class and the fork maps for the artificial classes, and
     * the auxiliary delay map (each join node is mapped to a corresponding auxiliary delay).
     */
    public static FJApproxReturn ht(Network model){
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
            line_error(mfilename(new Object(){}),"Could not create a copy of the model in the Heidelberger-Trivedi method");
            return null;
        }
//        nonfjmodel.allowReplace = true; -- no allowReplace implemented in JLINE
        Map<JobClass, Map<JobClass, Matrix>> P = nonfjmodel.getStruct(false).rtorig;
        nonfjmodel.resetNetwork(true);
        nonfjmodel.resetStruct();
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);
        HashMap<Integer, ArrayList<Integer>> forkedClasses = new HashMap<>();
        ArrayList<Integer> forkIndexes = new ArrayList<>();
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Fork){
                forkIndexes.add(i);
            }
        }

        // d = fj_auxiliary_delays{j} for join j gives the additional delay d to mimic the sojourn time of the original classes
        HashMap<Integer, Integer> fj_auxiliary_delays = new HashMap<>();

        // Replace each fork with a router
        for(int f : forkIndexes){
            nonfjmodel.getNodes().set(f, new Router(nonfjmodel, nonfjmodel.getNodes().get(f).getName()));
            /*
             * Done because the Router is added twice, once in the constructor, and once in the replacement.
             * Hence, we need to remove the reference added by the constructor from the nodes list.
             */
            nonfjmodel.getNodes().remove(nonfjmodel.getNodes().size() - 1);
            ArrayList<Integer> forked = new ArrayList<>();
            for(int i = 0; i < Vnodes.getNumCols(); i++){
                if(Vnodes.get(f, i) > 0){
                    forked.add(i);
                }
            }
            forkedClasses.put(f, forked);
        }

        // Replace each join with a delay
        for(int j = 0; j < sn.nodetypes.size(); j++){
            if(sn.nodetypes.get(j) != NodeType.Join){
                continue;
            }
            nonfjmodel.getNodes().set(j, new Delay(nonfjmodel, nonfjmodel.getNodes().get(j).getName()));
            nonfjmodel.getStations().set(model.getNodes().get(j).getStationIdx(), (Station) nonfjmodel.getNodes().get(j));
            /*
             * Same problem as above, the delay is added twice to the stations and to the node list. Must remove the
             * last reference from the two lists.
             */
            nonfjmodel.getNodes().remove(nonfjmodel.getNodes().size() - 1);
            nonfjmodel.getStations().remove(nonfjmodel.getStations().size() - 1);
            for(int c = 0; c < nonfjmodel.getClasses().size(); c++){
                ((Delay) nonfjmodel.getNodes().get(j)).setService(nonfjmodel.getClasses().get(c), new Immediate());
            }

            // Add another delay to mimic the sojourn time of the original classes for the artificial classes
            Delay newDelay = new Delay(nonfjmodel, "Auxiliary Delay - " + nonfjmodel.getNodes().get(j).getName());
            fj_auxiliary_delays.put(j, newDelay.getNodeIdx());
            for(int r = 0; r < nonfjmodel.getClasses().size(); r++){
                newDelay.setService(nonfjmodel.getClasses().get(r), new Immediate());
            }
            for(JobClass r : P.keySet()){
                for(JobClass s : P.get(r).keySet()){
                    Matrix Prs = P.get(r).get(s);
                    Matrix newPrs = new Matrix(Prs.getNumRows() + 1, Prs.getNumCols() + 1);
                    for(int i = 0; i < Prs.getNumRows(); i++){
                        for(int k = 0; k < Prs.getNumCols(); k++){
                            newPrs.set(i, k, Prs.get(i, k));
                        }
                    }
                    for(int k = 0; k < Prs.getNumCols(); k++){
                        newPrs.set(newDelay.getNodeIdx(), k, Prs.get(j, k));
                    }
                    P.get(r).put(s, newPrs);
                }
                for(int i = 0; i < P.get(r).get(r).getNumCols(); i++){
                    P.get(r).get(r).set(j, i, 0);
                }
                P.get(r).get(r).set(j, newDelay.getNodeIdx(), 1);
            }
        }

        nonfjmodel.setConnectionMatrix(new Matrix(nonfjmodel.getNodes().size(), nonfjmodel.getNodes().size()));

        // Create the artificial classes
        for(int f : forkIndexes){
            List<OutputStrategy> outputStrategies = model.getNodes().get(f).getOutput().getOutputStrategies();
            int joinIdx = -1;
            for(int i = 0; i < sn.fj.getNumCols(); i++){
                if(sn.fj.get(f, i) != 0){
                    if(joinIdx == -1){
                        joinIdx = i;
                    } else {
                        throw new RuntimeException("SolverMVA supports at present only a single join station per fork node.");
                    }
                }
            }
            ArrayList<Integer> forkedChains = new ArrayList<>();
            for(int i = 0; i < sn.chains.getNumRows(); i++){
                double sum = 0;
                for(int c : forkedClasses.get(f)){
                    sum += sn.chains.get(i, c);
                }
                if(sum != 0){
                    forkedChains.add(i);
                }
            }
            for(int fc : forkedChains){
                HashMap<AuxClassKey, JobClass> auxClasses = new HashMap<>();
                for(int r = 0; r < sn.chains.getNumCols(); r++){
                    if(sn.chains.get(fc, r) == 0 || sn.nodevisits.get(fc).get(f, r) == 0){
                        continue;
                    }
                    // Assumption: every class forks into exactly the same parallel branches
                    int parallelBranches = 0;
                    for(OutputStrategy o : outputStrategies){
                        if(o.getJobClass() == model.getClasses().get(r) && o.getDestination() != null){
                            parallelBranches++;
                        }
                    }
                    for(int par = 0; par < parallelBranches; par++){
                        // One artificial class for each parallel branch
                        if(((Forker) model.getNodes().get(f).getOutput()).taskPerLink > 1){
                            line_error(mfilename(new Object(){}),"Multiple tasks per link are not supported in H-T.");
                        }
                        int auxPopulation = (int) (((Forker) model.getNodes().get(f).getOutput()).taskPerLink *
                                ((ClosedClass) model.getClasses().get(r)).getPopulation());
                        auxClasses.put(new AuxClassKey(r, par), new ClosedClass(nonfjmodel,
                                nonfjmodel.getClasses().get(r).getName() + "." + nonfjmodel.getNodes().get(f).getName()
                                + ".B" + par, auxPopulation, (Station) nonfjmodel.getNodes().get(fj_auxiliary_delays.get(joinIdx)), 0));
                        fjclassmap.put(auxClasses.get(new AuxClassKey(r, par)).getIndex() - 1, nonfjmodel.getClasses().get(r).getIndex() - 1);
                        if(auxClasses.get(new AuxClassKey(r, par)).getIndex() - 1 >= fjclassmapSize){
                            fjclassmapSize = auxClasses.get(new AuxClassKey(r, par)).getIndex();
                        }
                        fjforkmap.put(auxClasses.get(new AuxClassKey(r, par)).getIndex() - 1, f);
                        if(auxClasses.get(new AuxClassKey(r, par)).getIndex() - 1 >= fjforkmapSize){
                            fjforkmapSize = auxClasses.get(new AuxClassKey(r, par)).getIndex();
                        }

                        // Set the service rates at the join node and at the stations
                        for(int i = 0; i < sn.nnodes; i++){
                            if(sn.isstation.get(i) != 0){
                                switch (sn.nodetypes.get(i)){
                                    case Join:
                                        ((Queue) nonfjmodel.getNodes().get(i)).setService(auxClasses.get(new AuxClassKey(r, par)),
                                                new Immediate());
                                        break;
                                    case Source: case Fork:
                                        // No-op
                                        break;
                                    default:
                                        Distribution distributionCopy = null;
                                        try{
                                            ByteArrayOutputStream bos = new ByteArrayOutputStream();
                                            ObjectOutputStream out = new ObjectOutputStream(bos);
                                            out.writeObject(((Queue) model.getNodes().get(i)).getService(model.getClasses().get(r)));
                                            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                                            ObjectInputStream in = new ObjectInputStream(bis);
                                            distributionCopy = (Distribution) in.readObject();
                                        } catch (IOException | ClassNotFoundException e) {
                                            line_error(mfilename(new Object(){}),"Could not copy the distribution of the original class in H-T");
                                        }
                                        ((Queue) nonfjmodel.getNodes().get(i)).setService(auxClasses.get(new
                                                        AuxClassKey(r, par)), distributionCopy);
                                }
                            }
                        }
                        ((Queue) nonfjmodel.getNodes().get(fj_auxiliary_delays.get(joinIdx))).setService(auxClasses.get(new AuxClassKey(r, par)),
                                new Immediate());
                    }
                }

                // Set the routing of the artificial classes
                for(int r = 0; r < sn.chains.getNumCols(); r++) {
                    if (sn.chains.get(fc, r) == 0 || sn.nodevisits.get(fc).get(f, r) == 0) {
                        continue;
                    }
                    for(int s = 0; s < sn.chains.getNumCols(); s++) {
                        if (sn.chains.get(fc, s) == 0 || sn.nodevisits.get(fc).get(f, s) == 0) {
                            continue;
                        }
                        int par = 0;
                        for(OutputStrategy o : outputStrategies){
                            if(o.getJobClass() != model.getClasses().get(r)){
                                continue;
                            }
                            JobClass rpar = auxClasses.get(new AuxClassKey(r, par));
                            JobClass spar = auxClasses.get(new AuxClassKey(s, par));
                            if(!P.containsKey(rpar)){
                                P.put(rpar, new HashMap<>());
                            }
                            Map<JobClass, Matrix> m = P.get(rpar);
                            m.put(spar, new Matrix(P.get(nonfjmodel.getJobClassFromIndex(r)).get(nonfjmodel.getJobClassFromIndex(s))));
                            for(int i = 0; i < m.get(spar).getNumCols(); i++){
                                m.get(spar).set(f, i, 0);
                            }
                            m.get(spar).set(f, o.getDestination().getNodeIdx(), 1);
                            m.get(spar).set(joinIdx, fj_auxiliary_delays.get(joinIdx), 1);
                            for(int i = 0; i < m.get(spar).getNumCols(); i++){
                                m.get(spar).set(fj_auxiliary_delays.get(joinIdx), i, 0);
                            }
                            m.get(spar).set(fj_auxiliary_delays.get(joinIdx), f, 1);
                            par++;
                        }
                        // Route the original classes straight to the join to avoid the interference with the artificial classes
                        Matrix Prs = P.get(nonfjmodel.getJobClassFromIndex(r)).get(nonfjmodel.getJobClassFromIndex(s));
                        for(int i = 0; i < Prs.getNumCols(); i++){
                            Prs.set(f, i, 0);
                        }
                        Prs.set(f, joinIdx, 1);
                    }
                }
            }
        }
        RoutingMatrix routingMatrix = new RoutingMatrix(nonfjmodel, nonfjmodel.getClasses(), nonfjmodel.getNodes());
        for(JobClass r : P.keySet()){
            for(JobClass s : P.get(r).keySet()){
                Matrix Prs = P.get(r).get(s);
                for(int i = 0; i < Prs.getNumRows(); i++){
                    for(int j = 0; j < Prs.getNumCols(); j++){
                        if(Prs.get(i, j) != 0){
                            routingMatrix.set(r, s, nonfjmodel.getNodes().get(i), nonfjmodel.getNodes().get(j), Prs.get(i, j));
                        }
                    }
                }
            }
        }
        nonfjmodel.link(routingMatrix);
        Matrix fjclassmapMatrix = new Matrix(1, fjclassmapSize);
        fjclassmapMatrix.fill(-1);
        for(int r : fjclassmap.keySet()){
            fjclassmapMatrix.set(r, fjclassmap.get(r));
        }
        Matrix fjforkmapMatrix = new Matrix(1, fjforkmapSize);
        fjforkmapMatrix.fill(-1);
        for(int r : fjforkmap.keySet()){
            fjforkmapMatrix.set(r, fjforkmap.get(r));
        }
        return new FJApproxReturn(nonfjmodel, fjclassmapMatrix, fjforkmapMatrix, fj_auxiliary_delays, null);
    }

    /**
     * Fork-Join Transform approach to evaluate queueing networks including fork-join systems. An equivalent network is
     * created where the fork nodes are replaced by routers, the join nodes are replaced by delays, and the parallelism
     * induced by a fork-join system is simulated through the addition of artificial open customer classes.
     * @param model - the original queueing network
     * @param forkLambda - the arrival rates of the artificial classes
     * @return - the equivalent queueing network with the fork-join systems replaced with other nodes, a mapping of the
     * artificial classes and their corresponding original classes, a mapping of the artificial classes and their FJ
     * systems, and the fanout of each artificial class
     */
    public static FJApproxReturn mmt(Network model, Matrix forkLambda){
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
            line_error(mfilename(new Object(){}),"Could not create a copy of the model in the Heidelberger-Trivedi method");
            return null;
        }
//        nonfjmodel.allowReplace = true; -- no allowReplace implemented in JLINE
        Map<JobClass, Map<JobClass, Matrix>> P = nonfjmodel.getStruct(false).rtorig;
        nonfjmodel.resetNetwork(true);
        nonfjmodel.resetStruct();
        Matrix Vnodes = Matrix.cellsum(sn.nodevisits);
        HashMap<Integer, ArrayList<Integer>> forkedClasses = new HashMap<>();
        ArrayList<Integer> forkIndexes = new ArrayList<>();
        int maxForkIdx = -1;
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Fork){
                forkIndexes.add(i);
                if(i > maxForkIdx){
                    maxForkIdx = i;
                }
            }
        }

        Matrix origfanout = new Matrix(maxForkIdx + 1, nonfjmodel.getNumberOfClasses());

        // Replace each fork with a router
        for(int f : forkIndexes){
            for(JobClass r : P.keySet()){
                List<OutputStrategy> outputStrategies = nonfjmodel.getNodes().get(f).getOutputStrategies();
                int parallelBranches = 0;
                for(OutputStrategy o : outputStrategies){
                    if(o.getJobClass() == r && o.getDestination() != null){
                        parallelBranches++;
                    }
                }
                if(parallelBranches > 0){
                    origfanout.set(f, r.getIndex() - 1, parallelBranches);
                    for(JobClass s : P.get(r).keySet()){
                        Matrix Prs = P.get(r).get(s);
                        for(int j = 0; j < Prs.getNumCols(); j++){
                            Prs.set(f, j, Prs.get(f, j) / parallelBranches);
                        }
                    }
                }
            }
            nonfjmodel.getNodes().set(f, new Router(nonfjmodel, nonfjmodel.getNodes().get(f).getName()));
            /*
             * Done because the Router is added twice, once in the constructor, and once in the replacement.
             * Hence, we need to remove the reference added by the constructor from the nodes list.
             */
            nonfjmodel.getNodes().remove(nonfjmodel.getNodes().size() - 1);
            ArrayList<Integer> forked = new ArrayList<>();
            for(int i = 0; i < Vnodes.getNumCols(); i++){
                if(Vnodes.get(f, i) > 0){
                    forked.add(i);
                }
            }
            forkedClasses.put(f, forked);
        }
        // Replace each join with a delay
        for(int j = 0; j < sn.nodetypes.size(); j++) {
            if (sn.nodetypes.get(j) != NodeType.Join) {
                continue;
            }
            nonfjmodel.getNodes().set(j, new Delay(nonfjmodel, nonfjmodel.getNodes().get(j).getName()));
            nonfjmodel.getStations().set(model.getNodes().get(j).getStationIdx(), (Station) nonfjmodel.getNodes().get(j));
            /*
             * Same problem as above, the delay is added twice to the stations and to the node list. Must remove the
             * last reference from the two lists.
             */
            nonfjmodel.getNodes().remove(nonfjmodel.getNodes().size() - 1);
            nonfjmodel.getStations().remove(nonfjmodel.getStations().size() - 1);
            for (int c = 0; c < nonfjmodel.getClasses().size(); c++) {
                ((Delay) nonfjmodel.getNodes().get(j)).setService(nonfjmodel.getClasses().get(c), new Immediate());
            }
        }
        Source source = null;
        Sink sink = null;
        if(nonfjmodel.hasOpenClasses()){
            source = nonfjmodel.getSource();
            sink = nonfjmodel.getSink();
        } else {
            source = new Source(nonfjmodel, "Source");
            sink = new Sink(nonfjmodel, "Sink");
            for(JobClass r : P.keySet()){
                for(JobClass s : P.get(r).keySet()){
                    Matrix Prs = P.get(r).get(s);
                    Matrix newPrs = new Matrix(Prs.getNumRows() + 2, Prs.getNumCols() + 2);
                    for(int i = 0; i < Prs.getNumRows(); i++){
                        for(int j = 0; j < Prs.getNumCols(); j++){
                            newPrs.set(i, j, Prs.get(i, j));
                        }
                    }
                    P.get(r).put(s, newPrs);
                }
            }
        }
        nonfjmodel.setConnectionMatrix(new Matrix(nonfjmodel.getNumberOfNodes(), nonfjmodel.getNumberOfNodes()));

        for(int f : forkIndexes) {
            int joinIdx = -1;
            for (int i = 0; i < sn.fj.getNumCols(); i++) {
                if (sn.fj.get(f, i) != 0) {
                    if (joinIdx == -1) {
                        joinIdx = i;
                    } else {
                        throw new RuntimeException("SolverMVA supports at present only a single join station per fork node.");
                    }
                }
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
                    oclass.add(new OpenClass(nonfjmodel, nonfjmodel.getJobClass().get(r).getName() + "." +
                            nonfjmodel.getNodes().get(f).getName()));
                    fjclassmap.put(oclass.get(oclass.size() - 1).getIndex() - 1, nonfjmodel.getJobClass().get(r).getIndex() - 1);
                    if(oclass.get(oclass.size() - 1).getIndex() - 1 >= fjclassmapSize){
                        fjclassmapSize = oclass.get(oclass.size() - 1).getIndex();
                    }
                    fjforkmap.put(oclass.get(oclass.size() - 1).getIndex() - 1, f);
                    if(oclass.get(oclass.size() - 1).getIndex() - 1 >= fjforkmapSize){
                        fjforkmapSize = oclass.get(oclass.size() - 1).getIndex();
                    }
                    int s = fjclassmap.get(oclass.get(oclass.size() - 1).getIndex() - 1);
                    if(((Forker) model.getNodes().get(f).getOutput()).taskPerLink > 1){
                        line_error(mfilename(new Object(){}),"There are no synchronisation delays implemented in FJT for multiple tasks per link.");
                    }
                    fanout.put(oclass.get(oclass.size() - 1).getIndex() - 1, (int) (origfanout.get(f, r) * ((Forker) model.getNodes().get(f).getOutput()).taskPerLink));
                    if(sn.nodevisits.get(fc).get(f, r) == 0){
                        source.setArrival(oclass.get(oclass.size() - 1), Disabled.getInstance());
                    } else {
                        source.setArrival(oclass.get(oclass.size() - 1), new Exp(forkLambda.get(r)));
                    }
                    // joins are now Delays, let us set their service time
                    for(int i = 0; i < sn.nnodes; i++){
                        if(sn.isstation.get(i) != 0){
                            switch (sn.nodetypes.get(i)){
                                case Join:
                                    ((Delay) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), new Immediate());
                                    break;
                                case Source: case Fork:
                                    // no-op
                                    break;
                                default:
                                    Distribution distributionCopy = null;
                                    try{
                                        ByteArrayOutputStream bos = new ByteArrayOutputStream();
                                        ObjectOutputStream out = new ObjectOutputStream(bos);
                                        out.writeObject(((Queue) model.getNodes().get(i)).getService(model.getClasses().get(r)));
                                        ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                                        ObjectInputStream in = new ObjectInputStream(bis);
                                        distributionCopy = (Distribution) in.readObject();
                                    } catch (IOException | ClassNotFoundException e) {
                                        System.err.println("Could not copy the distribution of the original class in FJT");
                                    }
                                    ((Queue) nonfjmodel.getNodes().get(i)).setService(oclass.get(oclass.size() - 1), distributionCopy);
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
                    for(int s = 0; s < sn.chains.getNumCols(); s++){
                        if(sn.chains.get(fc, s) == 0){
                            continue;
                        }
                        if(!P.containsKey(oclass.get(rIdx))){
                            P.put(oclass.get(rIdx), new HashMap<>());
                        }
                        P.get(oclass.get(rIdx)).put(oclass.get(sIdx),
                                new Matrix(P.get(nonfjmodel.getJobClass().get(r)).get(nonfjmodel.getJobClass().get(s))));
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
                    for(int s = 0; s < sn.chains.getNumCols(); s++){
                        if(sn.chains.get(fc, s) == 0){
                            continue;
                        }
                        Matrix Prs = P.get(oclass.get(rIdx)).get(oclass.get(sIdx));
                        for(int i = 0; i < Prs.getNumCols(); i++){
                            Prs.set(source.getNodeIdx(), i, 0);
                            Prs.set(joinIdx, i, 0);
                        }
                        sIdx++;
                    }
                    P.get(oclass.get(rIdx)).get(oclass.get(rIdx)).set(source.getNodeIdx(), f, 1);
                    P.get(oclass.get(rIdx)).get(oclass.get(rIdx)).set(joinIdx, sink.getNodeIdx(), 1);
                    rIdx++;
                }
            }
        }
        RoutingMatrix routingMatrix = new RoutingMatrix(nonfjmodel, nonfjmodel.getClasses(), nonfjmodel.getNodes());
        for(JobClass r : P.keySet()){
            for(JobClass s : P.get(r).keySet()){
                Matrix Prs = P.get(r).get(s);
                for(int i = 0; i < Prs.getNumRows(); i++){
                    for(int j = 0; j < Prs.getNumCols(); j++){
                        if(Prs.get(i, j) != 0){
                            routingMatrix.set(r, s, nonfjmodel.getNodes().get(i), nonfjmodel.getNodes().get(j), Prs.get(i, j));
                        }
                    }
                }
            }
        }
        nonfjmodel.link(routingMatrix);
        for(int f : forkIndexes){
            for(OutputStrategy o : nonfjmodel.getNodes().get(f).getOutputStrategies()){
                if(o.getRoutingStrategy() == RoutingStrategy.RAND){
                    o.setRoutingStrategy(RoutingStrategy.DISABLED);
                }
            }
        }
        Matrix fjclassmapMatrix = new Matrix(1, fjclassmapSize);
        fjclassmapMatrix.fill(-1);
        for(int r : fjclassmap.keySet()){
            fjclassmapMatrix.set(r, fjclassmap.get(r));
        }
        Matrix fjforkmapMatrix = new Matrix(1, fjforkmapSize);
        fjforkmapMatrix.fill(-1);
        for(int r : fjforkmap.keySet()){
            fjforkmapMatrix.set(r, fjforkmap.get(r));
        }
        return new FJApproxReturn(nonfjmodel, fjclassmapMatrix, fjforkmapMatrix, null, fanout);
    }

    /**
     * Finds the response times along each path leading out of startNode up to (and not including) endNode
     */
    public static Matrix findPaths(NetworkStruct sn, Matrix P, int startNode, int endNode, int r, ArrayList<Integer> toMerge,
                                   Matrix QN, Matrix TN, double currentTime, Matrix fjclassmap, Matrix fjforkmap,
                                   Network nonfjmodel){
        if(startNode == endNode){
            double qLen = 0;
            double tput = 0;
            for(int s : toMerge){
                qLen += QN.get((int) sn.nodeToStation.get(startNode), s);
                tput += TN.get((int) sn.nodeToStation.get(startNode), s);
            }
            Matrix ri = new Matrix(1,1);
            ri.set(0,0,currentTime - qLen/tput);
            return ri;
        }
        Matrix ri = new Matrix(1,0);
        for(int i = 0; i < P.getNumCols(); i++){
            if(P.get(startNode, i) == 0){
                continue;
            }
            double qLen = 0;
            double tput = 1;
            if(sn.nodeToStation.get(i) > -1){
                tput = 0;
                for(int s : toMerge){
                    qLen += QN.get((int) sn.nodeToStation.get(i), s);
                    tput += TN.get((int) sn.nodeToStation.get(i), s);
                }
            }
            if (sn.nodetypes.get(i) == NodeType.Fork){
                int joinIdx = 0;
                while(joinIdx < sn.fj.getNumCols() && sn.fj.get(i, joinIdx) == 0){
                    joinIdx++;
                }
                int s = 0;
                while(s < fjforkmap.length() && (fjforkmap.get(s) != i || fjclassmap.get(s) != r)){
                    s++;
                }
                toMerge.add(s);
                Matrix paths = findPaths(sn, P, i, joinIdx, r, toMerge, QN,
                        TN, 0, fjclassmap, fjforkmap, nonfjmodel);
                Matrix lambdai = Matrix.ones(paths.getNumRows(), paths.getNumCols()).elementDiv(paths);
                double d0 = 0;
                int parallel_branches = paths.length();
                for(int pow = 0; pow < parallel_branches; pow++){
                    Matrix nk = Maths.nchoosek(lambdai, pow + 1);
                    nk = nk.sumRows();
                    double currentSum = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
                    d0 += Math.pow(-1, pow) * currentSum;
                }
                for(int cls : toMerge){
                    ((Delay) nonfjmodel.getNodes().get(joinIdx)).setService(nonfjmodel.getJobClass().get(cls),
                            Exp.fitMean(d0 - paths.elementSum()/paths.length()));
                }
                toMerge.remove(toMerge.size() - 1);
                ri = ri.concatCols(findPaths(sn, P, joinIdx, endNode, r,
                        toMerge, QN, TN, currentTime + d0, fjclassmap,
                        fjforkmap, nonfjmodel));
            } else {
                ri = ri.concatCols(findPaths(sn, P, i, endNode, r, toMerge,
                        QN, TN, currentTime + qLen/tput, fjclassmap, fjforkmap,
                        nonfjmodel));
            }
        }
        return ri;
    }

    /**
     * Determines a directed acyclic graph of relationships among fork nodes.
     */
    public static FJsortForksReturn sortForks(NetworkStruct sn, NetworkStruct nonfjstruct, Matrix fjforkmap, Matrix fjclassmap, Network nonfjmodel){
        Matrix forks = new Matrix(sn.nodetypes.size(), (int) fjclassmap.elementMax() + 1);
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Fork){
                for(int j = 0; j < forks.getNumCols(); j++){
                    forks.set(i, j, 1);
                }
            }
        }
        Matrix parents = new Matrix(1, sn.nodetypes.size());
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Fork){
                parents.set(i, i);
            }
        }
        for(int f = 0; f < sn.nodetypes.size(); f++){
            if(sn.nodetypes.get(f) == NodeType.Fork){
                int joinIdx = 0;
                while(joinIdx < sn.fj.getNumCols() && sn.fj.get(f, joinIdx) == 0){
                    joinIdx++;
                }
                for(int s = 0; s < fjforkmap.length(); s++){
                    if(fjforkmap.get(s) == f){
                        int r = (int) fjclassmap.get(s);
                        Matrix nested = new Matrix(sn.nodetypes.size(), 1);
                        for(int i = 0; i < sn.nodetypes.size(); i++){
                            if(sn.nodetypes.get(i) == NodeType.Fork){
                                nested.set(i, 1);
                            }
                        }
                        nestedForks(f, joinIdx, nonfjstruct.rtorig.get(nonfjmodel.getJobClass().get(r)).get(nonfjmodel.getJobClass().get(r)), nested, sn);
                        for(int i = 0; i < forks.getNumRows(); i++){
                            forks.set(i, r, (int) forks.get(i, r) & (int) nested.get(i));
                        }
                        for(int i = 0; i < nested.length(); i++){
                            if(nested.get(i) == 0){
                                parents.set(i, parents.get(f));
                            }
                        }
                    }
                }
            }
        }
        return new FJsortForksReturn(forks, parents);
    }

    private static void nestedForks(int startNode, int endNode, Matrix conn, Matrix forks, NetworkStruct sn){
        if(startNode == endNode){
            return;
        }
        for(int i = 0; i < conn.getNumCols(); i++){
            if(conn.get(startNode, i) == 0){
                continue;
            }
            if(sn.nodetypes.get(i) == NodeType.Fork){
                forks.set(i, 0);
            }
            nestedForks(i, endNode, conn, forks, sn);
        }
    }

    private static class AuxClassKey {

        private final int x;
        private final int y;

        public AuxClassKey(int x, int y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof AuxClassKey)) return false;
            AuxClassKey key = (AuxClassKey) o;
            return x == key.x && y == key.y;
        }

        @Override
        public int hashCode() {
            return Objects.hash(x, y);
        }
    }

    public static class FJsortForksReturn{
        public Matrix outerForks;
        public Matrix parentForks;

        public FJsortForksReturn(Matrix outerForks, Matrix parentForks) {
            this.outerForks = outerForks;
            this.parentForks = parentForks;
        }
    }

    public static class FJApproxReturn{
        public Network nonfjmodel;
        public Matrix fjclassmap;
        public Matrix fjforkmap;
        public Map<Integer, Integer> fj_auxiliary_delays;
        public Map<Integer, Integer> fanout;

        public FJApproxReturn(Network nonfjmodel, Matrix fjclassmap, Matrix fjforkmap, Map<Integer, Integer> fj_auxiliary_delays, Map<Integer, Integer> fanout) {
            this.nonfjmodel = nonfjmodel;
            this.fjclassmap = fjclassmap;
            this.fjforkmap = fjforkmap;
            this.fj_auxiliary_delays = fj_auxiliary_delays;
            this.fanout = fanout;
        }
    }
}
