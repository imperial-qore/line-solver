package jline.api.sn;

import java.util.ArrayList;
import java.util.HashSet;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

public final class SnGetProductFormChainParams {
    private SnGetProductFormChainParams() {}

    /**
     * Calculate the parameters at class and chain level for a queueing network model.
     *
     * @param sn NetworkStruct object for the queueing network model.
     * @return queueing network parameters
     */
    public static Ret.snGetProductFormParams snGetProductFormChainParams(NetworkStruct sn) {
        Ret.snGetProductFormParams ret1 = SnGetProductFormParams.snGetProductFormParams(sn);
        Matrix lambda = ret1.lambda;
        Matrix mu = ret1.mu;
        ArrayList<Integer> queueIndex = new ArrayList<Integer>();
        ArrayList<Integer> delayIndex = new ArrayList<Integer>();
        HashSet<Integer> ignoreIndex = new HashSet<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            NodeType nt = sn.nodetype.get(i);
            if (nt == NodeType.Queue) {
                queueIndex.add(i);
            } else if (nt == NodeType.Delay) {
                delayIndex.add(i);
            } else if (nt == NodeType.Source || nt == NodeType.Join) {
                ignoreIndex.add((int) sn.nodeToStation.get(i));
            }
        }
        Ret.snGetDemands ret2 = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Dchain = ret2.Dchain;
        Matrix Vchain = ret2.Vchain.copy();
        Matrix Nchain = ret2.Nchain;
        Matrix lambda_chains = new Matrix(1, sn.nchains);

        Matrix D_chains = new Matrix(queueIndex.size(), sn.nchains);
        Matrix Z_chains = new Matrix(delayIndex.size(), sn.nchains);

        for (int c = 0; c < sn.nchains; c++) {
            double lambdaSum = 0.0;
            Matrix inc = sn.inchain.get(c);
            for (int i = 0; i < inc.getNumCols(); i++) {
                int idx = (int) inc.get(i);
                if (!Double.isNaN(lambda.get(idx))) {
                    lambdaSum += lambda.get(idx);
                }
            }
            lambda_chains.set(0, c, lambdaSum);
            for (int i = 0; i < queueIndex.size(); i++) {
                D_chains.set(i, c, Dchain.get((int) sn.nodeToStation.get(queueIndex.get(i)), c));
            }
            for (int i = 0; i < delayIndex.size(); i++) {
                Z_chains.set(i, c, Dchain.get((int) sn.nodeToStation.get(delayIndex.get(i)), c));
            }
        }
        Matrix S = new Matrix(queueIndex.size(), 1);
        for (int i = 0; i < queueIndex.size(); i++) {
            S.set(i, 0, sn.nservers.get((int) sn.nodeToStation.get(queueIndex.get(i))));
        }
        Ret.snGetProductFormParams ret = new Ret.snGetProductFormParams();
        ret.S = S;
        ret.lambda = lambda_chains;
        ret.N = Nchain;
        Vchain.removeRows(ignoreIndex);
        ret.D = D_chains;
        ret.Z = Z_chains;
        ret.mu = mu;
        ret.V = Vchain;
        if (ret.Z.isEmpty()) {
            ret.Z = new Matrix(ret.N.getNumRows(), ret.N.getNumCols());
        }
        return ret;
    }
}
