package jline.lang.state;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.GlobalConstants;
import jline.util.matrix.Matrix;

import java.io.Serializable;

public class AfterEventRouter implements Serializable {
    static Ret.EventResult afterEventRouter(NetworkStruct sn, int ind, EventType event, int jobClass, boolean isSimulation, Object eventCache, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1.0);

        switch (event) {
            case ARV:
                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                    spaceSrv.set(row, jobClass, spaceSrv.get(row, jobClass) + 1);
                }
                // buf is empty
                outspace = Matrix.concatColumns(spaceSrv, spaceVar, null);
                // passive action, rate is unspecified
                outrate = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                outrate.ones();
                outrate.scaleEq(-1);
                break;
            case DEP:
                if (spaceSrv.get(jobClass) > 0) {
                    for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        spaceSrv.set(row, jobClass, spaceSrv.get(row, jobClass) - 1);
                    }
                    switch (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass))) {
                        case RROBIN:
                            // Calculate sum of nvars for this node up to current class
                            // MATLAB: sum(sn.nvars(ind,1:(sn.nclasses+class))) sums columns 1 to (nclasses+class) (1-based)
                            // Java: columns 0 to (nclasses+jobClass) (0-based), same logical data
                            // Since MATLAB class = jobClass + 1, extract columns 0 to nclasses+jobClass
                            int nvar_sum = 0;
                            for (int c = 0; c <= sn.nclasses + jobClass; c++) {
                                nvar_sum += (int) sn.nvars.get(ind, c);
                            }

                            Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(jobClass));
                            int idx = -1;

                            // Convert from MATLAB 1-based position to Java 0-based index
                            // (same pattern as AfterEventStation line 1064)
                            int spaceVarIdx = nvar_sum - 1;

                            // outlinks is a row vector (1×N), use length() not getNumRows()
                            // (same pattern as AfterEventStation line 1071)
                            int numOutlinks = (int) outlinks.length();

                            // Find current outlink index
                            for (int i = 0; i < numOutlinks; i++) {
                                if (spaceVar.get(spaceVarIdx) == outlinks.get(i)) {
                                    idx = i;
                                    break;
                                }
                            }

                            // Update to next outlink (round-robin)
                            if (idx < numOutlinks - 1) {
                                spaceVar.set(spaceVarIdx, outlinks.get(idx + 1));
                            } else {
                                spaceVar.set(spaceVarIdx, outlinks.get(0));
                            }
                            break;
                    }
                    // buf is empty
                    outspace = Matrix.concatColumns(spaceSrv, spaceVar, null);
                    // immediate action
                    outrate = new Matrix(outspace.getNumRows(), 1);
                    outrate.ones();
                    outrate.scaleEq(GlobalConstants.Immediate);
                }
                break;
        }

        // Ensure outprob has the same number of rows as outspace
        // In MATLAB, outprob = 1 (scalar) which gets expanded to match outspace rows
        if (!outspace.isEmpty() && outspace.getNumRows() > 1) {
            outprob = new Matrix(outspace.getNumRows(), 1);
            outprob.ones();
        }

        // Handle simulation logic
        if (isSimulation) {
            if (eventCache != null && key != null) {
                // Store in cache - simplified for Java (cache implementation depends on specific cache type)
            }

            if (outspace.getNumRows() > 1) {
                double totRate = outrate.elementSum();
                Matrix cumRate = outrate.copy();
                double cumSum = 0;
                for (int i = 0; i < cumRate.getNumRows(); i++) {
                    cumSum += cumRate.get(i, 0);
                    cumRate.set(i, 0, cumSum / totRate);
                }

                // Select action based on random number - MATLAB: firing_ctr = 1 + max([0,find( rand > cum_rate' )])
                double rand = Math.random();
                int firingCtr = 0; // Default to first action (0-indexed)
                for (int i = 0; i < cumRate.getNumRows(); i++) {
                    if (rand > cumRate.get(i, 0)) {
                        firingCtr = i + 1; // Found last index where rand > cumRate, so firingCtr is next
                    } else {
                        break; // First cumRate element >= rand, stop searching
                    }
                }

                // Extract selected row
                Matrix selectedSpace = new Matrix(1, outspace.getNumCols());
                Matrix.extract(outspace, firingCtr, firingCtr + 1, 0, outspace.getNumCols(), selectedSpace, 0, 0);
                outspace = selectedSpace;

                // outrate becomes sum of all rates (total rate)
                Matrix totalRate = new Matrix(1, 1);
                totalRate.set(0, 0, totRate);
                outrate = totalRate;

                // outprob becomes scalar 1.0 for single selected action
                outprob = new Matrix(1, 1);
                outprob.set(0, 0, 1.0);
            }
        }

        return new Ret.EventResult(outspace, outrate, outprob);
    }
}