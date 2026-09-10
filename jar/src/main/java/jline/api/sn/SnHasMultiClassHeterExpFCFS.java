package jline.api.sn;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class SnHasMultiClassHeterExpFCFS {
    private SnHasMultiClassHeterExpFCFS() {}

    /**
     * Checks if the network has one or more stations with multiclass heterogeneous FCFS
     * and exponential service.
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasMultiClassHeterExpFCFS(NetworkStruct sn) {
        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.FCFS) {
                continue;
            }
            Matrix row = Matrix.extractRows(sn.rates, i, i + 1, null);
            // Compute max and min ignoring NaN values to match MATLAB's max/min behavior
            double rateMax = Double.NEGATIVE_INFINITY;
            double rateMin = Double.POSITIVE_INFINITY;
            boolean hasRate = false;
            for (int j = 0; j < row.getNumElements(); j++) {
                double v = row.get(j);
                if (!Double.isNaN(v)) {
                    if (v > rateMax) rateMax = v;
                    if (v < rateMin) rateMin = v;
                    hasRate = true;
                }
            }
            if (hasRate && rateMax - rateMin > 0) {
                Matrix scvs = Matrix.extractRows(sn.scv, i, i + 1, null);
                double scvMax = Double.NEGATIVE_INFINITY;
                double scvMin = Double.POSITIVE_INFINITY;
                for (int j = 0; j < scvs.getNumElements(); j++) {
                    double v = scvs.get(j);
                    if (!Double.isNaN(v)) {
                        if (v > scvMax) scvMax = v;
                        if (v < scvMin) scvMin = v;
                    }
                }
                if (scvMax < 1 + GlobalConstants.FineTol && scvMin > 1 - GlobalConstants.FineTol) {
                    return true;
                }
            }
        }
        return false;
    }
}
