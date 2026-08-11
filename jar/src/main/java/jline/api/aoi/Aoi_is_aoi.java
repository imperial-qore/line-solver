/**
 * @file AoI topology validation
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.ArrayList;
import java.util.List;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aoi_is_aoi {
    private Aoi_is_aoi() {}

    /**
     * Check if network is a valid AoI topology for analytical AoI analysis.
     */
    public static AoiValidationResult aoi_is_aoi(NetworkStruct sn) {
        // Check 1: Must have at least one open class
        Matrix njobs = sn.njobs;
        boolean hasOpen = false;
        List<Integer> openClasses = new ArrayList<Integer>();
        for (int r = 0; r < sn.nclasses; r++) {
            if (Double.isInfinite(njobs.get(r))) {
                hasOpen = true;
                openClasses.add(r);
            }
        }
        if (!hasOpen) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Not an open model - all classes are closed");
        }

        if (openClasses.size() > 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Multiple open classes found (" + openClasses.size()
                            + ") - AoI analysis requires single class");
        }

        int sourceNodeIdx = -1;
        int queueNodeIdx = -1;
        int sinkNodeIdx = -1;
        int sourceCount = 0;
        int queueCount = 0;
        int sinkCount = 0;

        List<NodeType> nodetype = sn.nodetype;
        for (int i = 0; i < nodetype.size(); i++) {
            NodeType nt = nodetype.get(i);
            if (nt == NodeType.Source) {
                sourceNodeIdx = i;
                sourceCount++;
            } else if (nt == NodeType.Queue) {
                queueNodeIdx = i;
                queueCount++;
            } else if (nt == NodeType.Sink) {
                sinkNodeIdx = i;
                sinkCount++;
            }
        }

        if (sourceCount == 0) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "", "No source node found");
        }
        if (sourceCount > 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Multiple source nodes found (" + sourceCount + ")");
        }

        if (sinkCount == 0) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "", "No sink node found");
        }
        if (sinkCount > 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Multiple sink nodes found (" + sinkCount + ")");
        }

        if (queueCount == 0) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "", "No queue node found");
        }
        if (queueCount > 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Multiple queue nodes found (" + queueCount
                            + ") - AoI analysis supports single queue only");
        }

        int sourceStation = (int) sn.nodeToStation.get(sourceNodeIdx);
        int queueStation = (int) sn.nodeToStation.get(queueNodeIdx);

        int nservers = (int) sn.nservers.get(queueStation);
        if (nservers != 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Queue has " + nservers + " servers - AoI analysis requires single server (c=1)");
        }

        double cap = sn.cap.get(queueStation);
        if (Double.isInfinite(cap) || cap > 2 || cap < 1) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Queue capacity is " + cap
                            + " - AoI analysis requires capacity 1 (bufferless) or 2 (single-buffer)");
        }

        int capacity = (int) cap;
        String systemType = (capacity == 1) ? "bufferless" : "singlebuffer";

        SchedStrategy schedStrategy = sn.sched.get(sn.stations.get(queueStation));
        String schedName = (schedStrategy != null) ? schedStrategy.toString() : "unknown";

        if (schedStrategy != SchedStrategy.FCFS
                && schedStrategy != SchedStrategy.LCFS
                && schedStrategy != SchedStrategy.LCFSPR) {
            return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                    "Unsupported scheduling strategy (" + schedName
                            + ") - AoI analysis supports FCFS, LCFS, or LCFSPR only");
        }

        if (capacity == 2) {
            int classIdx = openClasses.get(0);
            Station sourceStationObj = sn.stations.get(sourceStation);
            JobClass jobClassObj = sn.jobclasses.get(classIdx);
            MatrixCell arrivalProc = (sn.proc != null && sn.proc.get(sourceStationObj) != null)
                    ? sn.proc.get(sourceStationObj).get(jobClassObj) : null;
            if (arrivalProc != null && arrivalProc.size() >= 1) {
                Matrix D0 = arrivalProc.get(0);
                if (D0 != null && D0.getNumRows() > 1) {
                    return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                            "Single-buffer (capacity=2) requires exponential arrivals (Poisson process)");
                }
            }
        }

        int K = sn.nclasses;
        int classIdx2 = openClasses.get(0);
        int rtIdx = queueStation * K + classIdx2;
        if (sn.rt != null && rtIdx < sn.rt.getNumRows() && rtIdx < sn.rt.getNumCols()) {
            if (sn.rt.get(rtIdx, rtIdx) > 0) {
                return new AoiValidationResult(false, -1, -1, -1, -1, -1, -1, "", "",
                        "Self-loop detected at queue - violates AoI model assumptions");
            }
        }

        return new AoiValidationResult(true, sourceNodeIdx, queueNodeIdx, sinkNodeIdx,
                sourceStation, queueStation, capacity, schedName, systemType, "");
    }
}
