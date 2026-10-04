package jline.api.sn;

import jline.api.mam.Map_pie;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.JobClass;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Map;

/**
 * Detects whether a network has a bursty (non-renewal) external arrival process.
 *
 * <p>Returns true if any Source station has an arrival process with
 * autocorrelated inter-arrival times (a non-renewal Markovian arrival process
 * such as an MMPP/MAP), as opposed to a renewal process (Poisson, or any
 * i.i.d. renewal process such as Erlang/HyperExp/Coxian/APH). Detection is
 * exact: a MAP with matrices (D0,D1) is renewal iff D1 equals its rank-one
 * renewal form t0*pie, where t0 = -D0*e and pie is the embedded stationary
 * vector; any departure signals correlation between inter-arrival times.
 *
 * @since LINE 3.0
 */
public final class SnHasBurstyArrival {
    private SnHasBurstyArrival() {}

    /**
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if some external arrival process is non-renewal (bursty)
     */
    public static boolean snHasBurstyArrival(NetworkStruct sn) {
        if (sn.proc == null) {
            return false;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            int nd = (int) sn.stationToNode.get(ist);
            if (sn.nodetype.get(nd) != NodeType.Source) {
                continue;
            }
            Station st = sn.stations.get(ist);
            Map<JobClass, MatrixCell> procMap = sn.proc.get(st);
            if (procMap == null) {
                continue;
            }
            for (int r = 0; r < sn.nclasses; r++) {
                MatrixCell mapproc = procMap.get(sn.jobclasses.get(r));
                if (mapproc == null || mapproc.size() < 2 || mapproc.get(0) == null || mapproc.get(1) == null) {
                    continue;
                }
                Matrix D0 = mapproc.get(0);
                Matrix D1 = mapproc.get(1);
                int n = D1.getNumRows();
                if (n <= 1) {
                    continue;   // single-phase arrival is Poisson, hence renewal
                }
                Matrix pie = Map_pie.map_pie(D0, D1);           // 1 x n
                Matrix e = Matrix.ones(n, 1);
                Matrix D1ren = D1.mult(e).mult(pie);            // rank-one renewal form t0*pie
                double diff = D1.sub(D1ren).norm();
                double base = Math.max(1.0, D1.norm());
                if (diff > 1e-8 * base) {
                    return true;
                }
            }
        }
        return false;
    }
}
