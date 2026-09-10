package jline.opt.sensitivity;

import jline.api.pfqn.sens.Pfqn_sens;
import jline.io.Ret;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.opt.results.SensitivityData;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Analytic d(metric)/d(rate) for closed single-server product-form networks
 * with unit visit ratios, via the differentiated-MVA primitive
 * {@link Pfqn_sens}. Mirrors native-Python {@code _closed_sensitivities}: the
 * full cross-station Jacobian is assembled and the rate chain rule
 * d/d(rate) = -(L/rate) * d/dL is applied. Returns {@code null} for open,
 * multiserver, or non-unit-visit networks.
 */
final class ClosedSensitivity {

    private ClosedSensitivity() {
    }

    static SensitivityData compute(Network model) {
        NetworkStruct sn = model.getStruct();
        int R = sn.nclasses;
        Matrix N = sn.njobs;
        if (N == null || N.length() == 0) {
            return null;
        }
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(r)) || Double.isNaN(N.get(r))) {
                return null;   // closed only
            }
        }

        Ret.snGetProductFormParams pf = model.getProductFormParameters();
        Matrix D = pf.D;                       // Mq x R
        Matrix S = pf.S;                       // Mq x 1 server counts
        if (S != null) {
            for (int i = 0; i < S.length(); i++) {
                double s = S.get(i);
                if (!Double.isNaN(s) && !Double.isInfinite(s) && s > 1) {
                    return null;               // single-server only
                }
            }
        }
        int Mq = D.getNumRows();

        // Z (1 x R): sum delay demands over delay stations
        Matrix Z = new Matrix(1, R);
        if (pf.Z != null && pf.Z.getNumRows() > 0) {
            for (int r = 0; r < R; r++) {
                double sum = 0.0;
                for (int z = 0; z < pf.Z.getNumRows(); z++) {
                    sum += pf.Z.get(z, r);
                }
                Z.set(0, r, sum);
            }
        }

        Matrix Nrow = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Nrow.set(0, r, N.get(r));
        }

        Ret.pfqnSens sens = Pfqn_sens.pfqn_sens(D, Nrow, Z);

        Matrix rates = sn.rates;
        Matrix nodeToStation = sn.nodeToStation;
        List<Integer> queueNodes = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Queue) {
                queueNodes.add(i);
            }
        }
        // queue-station order must line up with the product-form D rows
        if (queueNodes.size() != Mq) {
            return null;
        }

        // precompute rate params: (queue j, class s) -> (paramIndex, chain, key)
        List<int[]> paramIdx = new ArrayList<int[]>();   // {j, s, p}
        List<Double> paramChain = new ArrayList<Double>();
        List<String> paramKey = new ArrayList<String>();
        for (int j = 0; j < queueNodes.size(); j++) {
            int nodeJ = queueNodes.get(j);
            int sj = (int) nodeToStation.get(nodeJ);
            for (int s = 0; s < R; s++) {
                double ratej = rates.get(sj, s);
                if (Double.isNaN(ratej) || Double.isInfinite(ratej) || ratej <= 0 || D.get(j, s) <= 0) {
                    continue;
                }
                paramIdx.add(new int[]{j, s, j * R + s});
                paramChain.add(-D.get(j, s) / ratej);
                paramKey.add(SensitivityData.paramKey(name(sn, nodeJ), className(sn, s)));
            }
        }

        SensitivityData out = new SensitivityData();
        for (int i = 0; i < queueNodes.size(); i++) {
            int nodeI = queueNodes.get(i);
            String nameI = name(sn, nodeI);
            for (int r = 0; r < R; r++) {
                if (D.get(i, r) <= 0) {
                    continue;
                }
                String mkey = SensitivityData.metricKey(nameI, className(sn, r));
                for (int pk = 0; pk < paramIdx.size(); pk++) {
                    int p = paramIdx.get(pk)[2];
                    double chain = paramChain.get(pk);
                    String pkey = paramKey.get(pk);
                    out.add("RespT", mkey, pkey, sens.dR[p].get(i, r) * chain);
                    out.add("QLen", mkey, pkey, sens.dQ[p].get(i, r) * chain);
                    out.add("Tput", mkey, pkey, sens.dX.get(r, p) * chain);   // unit visits
                    out.add("Util", nameI, pkey, sens.dU[p].get(i, r) * chain);
                }
            }
        }
        return out;
    }

    private static String name(NetworkStruct sn, int nodeIndex) {
        return sn.nodenames.get(nodeIndex);
    }

    private static String className(NetworkStruct sn, int r) {
        return sn.classnames.get(r);
    }
}
