/**
 * @file Product-form queueing network parameter extraction
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class SnGetProductFormParams {
    private SnGetProductFormParams() {}

    /**
     * Calculate the parameters at class level for a queueing network model.
     */
    public static Ret.snGetProductFormParams snGetProductFormParams(NetworkStruct sn) {
        int R = sn.nclasses;
        Matrix N = sn.njobs;
        List<Integer> queueIndices = new ArrayList<Integer>();
        List<Integer> delayIndices = new ArrayList<Integer>();
        int sourceIndex = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            NodeType nt = sn.nodetype.get(i);
            if (nt == NodeType.Queue) {
                queueIndices.add(i);
            } else if (nt == NodeType.Delay) {
                delayIndices.add(i);
            } else if (nt == NodeType.Source) {
                sourceIndex = i;
            }
        }
        int Mq = queueIndices.size();
        int Mz = delayIndices.size();

        Matrix lambda = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (Utils.isInf(N.get(r))) {
                lambda.set(0, r, sn.rates.get((int) sn.nodeToStation.get(sourceIndex), r));
            }
        }
        Matrix S = new Matrix(queueIndices.size(), 1);
        double Smax = Double.MIN_VALUE;
        for (int i = 0; i < queueIndices.size(); i++) {
            S.set(i, 0, sn.nservers.get((int) sn.nodeToStation.get(queueIndices.get(i))));
            if (Double.isFinite(S.get(i, 0)) && S.get(i, 0) > Smax) Smax = S.get(i, 0);
        }
        Matrix D = new Matrix(Mq, R);
        double Nct = 0.0;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                if (Double.isFinite(N.get(i, j))) Nct += N.get(i, j);
            }
        }
        Matrix mu = Matrix.ones(Mq, (int) (Math.ceil(Nct) + Smax));
        for (int ist = 0; ist < Mq; ist++) {
            for (int r = 0; r < R; r++) {
                int c = 0;
                while (c < sn.chains.getNumRows()) {
                    if (sn.chains.get(c, r) != 0.0) break;
                    c++;
                }
                if (sn.refclass.get(c) > 0) {
                    D.set(ist, r,
                        sn.visits.get(c).get((int) sn.nodeToStateful.get(queueIndices.get(ist)), r)
                            / sn.rates.get((int) sn.nodeToStation.get(queueIndices.get(ist)), r)
                            / sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), (int) sn.refclass.get(c)));
                } else {
                    D.set(ist, r,
                        sn.visits.get(c).get((int) sn.nodeToStateful.get(queueIndices.get(ist)), r)
                            / sn.rates.get((int) sn.nodeToStation.get(queueIndices.get(ist)), r));
                }
            }
            for (int j = 0; j < mu.getNumCols(); j++) {
                mu.set(ist, j, Maths.min((double) (j + 1), sn.nservers.get((int) sn.nodeToStation.get(queueIndices.get(ist)))));
            }
        }
        Matrix Z = new Matrix((int) Maths.max(1.0, (double) Mz), R);
        for (int ist = 0; ist < Mz; ist++) {
            for (int r = 0; r < R; r++) {
                int c = 0;
                while (c < sn.chains.getNumRows()) {
                    if (sn.chains.get(c, r) != 0.0) break;
                    c++;
                }
                if (sn.refclass.get(c) > 0) {
                    Z.set(ist, r,
                        sn.visits.get(c).get((int) sn.nodeToStateful.get(delayIndices.get(ist)), r)
                            / sn.rates.get((int) sn.nodeToStation.get(delayIndices.get(ist)), r)
                            / sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), (int) sn.refclass.get(c)));
                } else {
                    Z.set(ist, r,
                        sn.visits.get(c).get((int) sn.nodeToStateful.get(delayIndices.get(ist)), r)
                            / sn.rates.get((int) sn.nodeToStation.get(delayIndices.get(ist)), r));
                }
            }
        }
        Matrix V = null;
        for (Integer i : sn.visits.keySet()) {
            Matrix visitsValue = sn.visits.get(i);
            if (V == null) {
                V = visitsValue;
            } else {
                V = V.add(1.0, visitsValue);
            }
        }
        D.apply(Double.NaN, 0.0, "equal");
        Z.apply(Double.NaN, 0.0, "equal");
        Ret.snGetProductFormParams ret = new Ret.snGetProductFormParams();
        ret.lambda = lambda;
        ret.D = D;
        ret.N = N;
        ret.Z = Z;
        ret.mu = mu;
        ret.S = S;
        ret.V = V;
        return ret;
    }
}
