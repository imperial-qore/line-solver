package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

public final class SnHasMultiClassHeterFCFS {
    private SnHasMultiClassHeterFCFS() {}

    /**
     * Checks if the network has one or more stations with multiclass heterogeneous FCFS.
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return true if there is at least one heterogeneous multiclass FCFS station
     */
    public static boolean snHasMultiClassHeterFCFS(NetworkStruct sn) {
        for (int i = 0; i < sn.sched.size(); i++) {
            Station station = sn.stations.get(i);
            if (sn.sched.get(station) != SchedStrategy.FCFS) {
                continue;
            }
            Matrix row = Matrix.extractRows(sn.rates, i, i + 1, null);
            // see _kb/03-api-layer.md for rationale
            double max = Double.NEGATIVE_INFINITY;
            double min = Double.POSITIVE_INFINITY;
            boolean hasValue = false;
            for (int j = 0; j < row.getNumElements(); j++) {
                double v = row.get(j);
                if (!Double.isNaN(v)) {
                    if (v > max) max = v;
                    if (v < min) min = v;
                    hasValue = true;
                }
            }
            if (hasValue && max - min > 0) {
                return true;
            }
        }
        return false;
    }
}
