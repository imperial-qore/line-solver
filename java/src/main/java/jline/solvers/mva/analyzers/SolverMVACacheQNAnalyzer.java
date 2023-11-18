package jline.solvers.mva.analyzers;

import jline.api.CACHE;
import jline.api.SN;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;
import jline.lang.NodeParam;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.api.DTMC.dtmc_stochcomp;

/**
 * MVA Analyzer class for solver_mva_cacheqn_analyzer
 */
public class SolverMVACacheQNAnalyzer implements MVAAnalyzer{
    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        NetworkStruct snorig = sn;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(sn);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            snorig = (NetworkStruct) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            System.err.println("Could not create a copy of the NetworkStruct in SolverMVACacheQNAnalyzer");
        }
        sn = snorig;
        Random random = new Random(options.seed);
        int I = sn.nnodes;
        int K = sn.nclasses;
        ArrayList<Integer> statefulNodes = new ArrayList<>();
        for(int i = 0; i < sn.isstateful.length(); i++){
            if(sn.isstateful.get(i) == 1){
                statefulNodes.add(i);
            }
        }
        Matrix statefulNodeClasses = new Matrix(1, statefulNodes.size() * K);
        int idx = 0;
        for(int ind : statefulNodes){
            for(int i = 0; i < K; i++){
                statefulNodeClasses.set(idx, ind * K + i);
                idx++;
            }
        }
        Matrix lambda = new Matrix(1, K);
        Matrix lambda_1 = new Matrix(1, K);
        ArrayList<Integer> caches = new ArrayList<>();
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Cache){
                caches.add(i);
            }
        }
        Matrix hitprob = new Matrix(sn.nodetypes.size(), K);
        Matrix missprob = new Matrix(sn.nodetypes.size(), K);

        for(int it = 1; it <= options.iter_max; it++){
            ArrayList<Integer> inputClass = new ArrayList<>();
            for(int ind : caches){
                NodeParam ch = sn.nodeparam.get(sn.nodes.get(ind));
                Matrix hitClass = ch.hitclass;
                Matrix missClass = ch.missclass;
                for(int i = 0; i < hitClass.length(); i++){
                    if(hitClass.get(i) != -1){
                        inputClass.add(i);
                    }
                }
                if(it == 1){
                    for(int i = 0; i < inputClass.size(); i++){
                        lambda_1.set(inputClass.get(i), random.nextDouble());
                    }
                    lambda = new Matrix(lambda_1);
                    sn.nodetypes.set(ind, NodeType.ClassSwitch);
                }

                // Solution of the isolated cache
                Matrix m = ch.itemcap;
                int n = ch.nitems;
                int h = m.length();
                int u = lambda.length();

                Matrix missrate = new Matrix(sn.nodetypes.size(), u);

                Matrix[] lambda_cache = new Matrix[u];
                for(int i = 0; i < u; i++){
                    lambda_cache[i] = new Matrix(n, h + 1);
                }

                for(int v = 0; v < u; v++){
                    for(int k = 0; k < n; k++){
                        for(int l = 0; l < h + 1; l++){
                            if(ch.pread.getOrDefault(v, null) != null){
                                lambda_cache[v].set(k, l, lambda.get(v) * ch.pread.get(v).get(k));
                            }
                        }
                    }
                }

                Matrix[][] R = ch.accost;
                Matrix gamma = CACHE.cache_gamma_lp(lambda_cache, R).gamma;
                Matrix pij;
                if(options.method.equals("exact")){
                    pij = CACHE.cache_mva(gamma, m).pij;
                    Matrix newpij = new Matrix(pij.getNumRows(), pij.getNumCols() + 1);
                    for(int i = 0; i < newpij.getNumRows(); i++){
                        for(int j = 0; j < newpij.getNumCols(); j++){
                            if(j == 0){
                                newpij.set(i, j, Math.abs(1 - pij.sumRows(i)));
                            } else {
                                newpij.set(i, j, pij.get(i, j - 1));
                            }
                        }
                    }
                    pij = newpij;
                } else {
                    pij = CACHE.cache_prob_asy(gamma, m); // FPI method
                }
                for(int i = 0; i < missprob.getNumCols(); i++){
                    missprob.set(ind, i, 0);
                }
                for(int v = 0; v < u; v++){
                    Matrix A = Matrix.extractColumn(lambda_cache[v], 0, null).transpose();
                    Matrix B = Matrix.extractColumn(pij, 0, null);
                    missrate.set(ind, v, A.mult(B).get(0));
                }
                for(int i = 0; i < missprob.getNumCols(); i++){
                    double val = missrate.get(ind, i)/lambda.get(i);
                    if(Double.isNaN(val)){
                        missprob.set(ind, i, 0);
                    } else {
                        missprob.set(ind, i, val);
                    }
                }
                for(int i = 0; i < hitprob.getNumCols(); i++){
                    double val = 1 - missrate.get(ind, i)/lambda.get(i);
                    if(Double.isNaN(val)){
                        hitprob.set(ind, i, 0);
                    } else {
                        hitprob.set(ind, i, val);
                    }
                }
                // Bring back isolated model results into the queueing model
                for(int i = 0; i < inputClass.size(); i++){
                    int r = inputClass.get(i);
                    for(int j = 0; j < sn.rtnodes.getNumCols(); j++){
                        sn.rtnodes.set(ind * K + r, j, 0);
                    }
                    for(int jnd = 0; jnd < I; jnd++){
                        if(sn.connmatrix.get(ind, jnd) == 1){
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + hitClass.get(r)), hitprob.get(ind, r));
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + missClass.get(r)), missprob.get(ind, r));
                        }
                    }
                }
                List<Integer> statefulNodeClassesList = new ArrayList<>(statefulNodeClasses.length());
                for(int i = 0; i < statefulNodeClasses.length(); i++){
                    statefulNodeClassesList.add((int) statefulNodeClasses.get(i));
                }
                sn.rt = dtmc_stochcomp(sn.rtnodes, statefulNodeClassesList);
            }
            sn = SN.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);

            switch(options.method){
                case "aba.upper": case "aba.lower": case "bjb.upper": case "bjb.lower": case "pb.upper": case "pb.lower":
                case "gb.upper": case "gb.lower": case "sb.upper": case "sb.lower":
                    new SolverMVABoundAnalyzer().analyze(sn, options, res);
                    break;
                default:
                    if(!(sn.lldscaling == null || sn.lldscaling.isEmpty()) || !(sn.cdscaling == null || sn.cdscaling.isEmpty())){
                        new SolverMVALDAnalyzer().analyze(sn, options, res);
                    } else {
                        new SolverMVAAnalyzer().analyze(sn, options, res);
                    }
                    break;
            }

            Matrix nodevisits = null;
            for(Integer key : sn.nodevisits.keySet()){
                if(nodevisits == null){
                    nodevisits = sn.nodevisits.get(key);
                } else {
                    nodevisits = nodevisits.add(1, sn.nodevisits.get(key));
                }
            }
            for(int ind : caches) {
                for (int i = 0; i < inputClass.size(); i++) {
                    int r = inputClass.get(i);
                    int c = 0;
                    while (c < sn.chains.getNumRows() && sn.chains.get(c, r) == 0) {
                        c++;
                    }
                    ArrayList<Integer> inchain = new ArrayList<>();
                    for (int j = 0; j < sn.chains.getNumCols(); j++) {
                        if (sn.chains.get(c, j) != 0) {
                            inchain.add(j);
                        }
                    }
                    double sumXN = 0;
                    for (int ix : inchain) {
                        sumXN += res.XN.get(ix);
                    }
                    if (sn.refclass.get(c) > -1) {
                        lambda.set(r, sumXN * nodevisits.get(ind, r) / nodevisits.get((int) sn.stationToNode.get((int) sn.refstat.get(r)), (int) sn.refclass.get(c)));
                    } else {
                        lambda.set(r, sumXN * nodevisits.get(ind, r) / nodevisits.get((int) sn.stationToNode.get((int) sn.refstat.get(r)), r));
                    }
                }
            }
            if(lambda.sub(1, lambda_1).elementSum() < options.iter_tol){
                res.iter = it;
                break;
            }
            lambda_1 = lambda;
        }
        res.hitProb = hitprob;
        res.missProb = missprob;
    }
}
