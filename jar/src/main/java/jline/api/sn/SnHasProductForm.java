/**
 * @file Stochastic network product-form solution checker
 *
 * Determines if a queueing network has a known product-form solution by validating
 * scheduling disciplines, class structures, and network topology. Product-form networks
 * enable efficient exact analysis using algorithms like MVA and convolution.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasProductForm {
    private SnHasProductForm() {}

    /**
     * Checks if the network has a known product-form solution.
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasProductForm(NetworkStruct sn) {
        boolean ret = true;
        boolean hasLcfs = false;
        boolean hasLcfspr = false;
        int lcfsCount = 0;
        int lcfsprCount = 0;

        for (int i = 0; i < sn.sched.size(); i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.LCFS) {
                hasLcfs = true;
                lcfsCount++;
            } else if (sched == SchedStrategy.LCFSPR) {
                hasLcfspr = true;
                lcfsprCount++;
            } else if (sched == SchedStrategy.INF
                    || sched == SchedStrategy.PS
                    || sched == SchedStrategy.FCFS
                    || sched == SchedStrategy.EXT) {
                // ok
            } else {
                ret = false;
            }
        }

        // see _kb/03-api-layer.md for rationale
        if (hasLcfs && hasLcfspr && lcfsCount == 1 && lcfsprCount == 1) {
            // This is a valid product-form configuration
        } else if (hasLcfs) {
            // LCFS alone (without LCFSPR pairing) is not product-form
            ret = false;
        }

        ret = ret && !SnHasMultiClassHeterFCFS.snHasMultiClassHeterFCFS(sn);
        ret = ret && !SnHasPriorities.snHasPriorities(sn);
        ret = ret && !SnHasForkJoin.snHasForkJoin(sn);
        ret = ret && !SnHasSDRouting.snHasSDRouting(sn);
        // BCMP asks for infinite buffers. Nothing here read sn.cap/sn.classcap/
        // sn.droprule, so a BAS-blocked station or any binding finite buffer passed
        // the gate and the network read as product form while its truncation
        // couples the station occupancies.
        ret = ret && !SnHasBlocking.snHasBlocking(sn);

        // BCMP type 1 asks the FCFS service to be exponential. SnHasMultiClassHeterFCFS
        // compares the class MEANS only, so a class-homogeneous Erlang, hyper-exponential
        // or deterministic FCFS station passed this gate and was dispatched to exact MVA,
        // which reads the means alone and returns the exponential answer with no warning.
        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.FCFS) {
                continue;
            }
            for (int r = 0; r < sn.scv.getNumCols(); r++) {
                double scv = sn.scv.get(i, r);
                if (Double.isFinite(scv) && scv > 0) {
                    ret = ret && (scv > 1 - GlobalConstants.FineTol) && (scv < 1 + GlobalConstants.FineTol);
                }
            }
        }
        return ret;
    }
}
