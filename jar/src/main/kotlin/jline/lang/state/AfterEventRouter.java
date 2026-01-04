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
                            // MATLAB: sum(sn.nvars(ind,1:(sn.nclasses+class))) - 1-indexed in MATLAB
                            int nvar_sum = 0;
                            for (int c = 0; c < sn.nclasses + jobClass; c++) {
                                nvar_sum += (int) sn.nvars.get(ind, c);
                            }
                            
                            Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(jobClass));
                            int idx = -1;
                            
                            // Find current outlink index
                            for (int i = 0; i < outlinks.getNumRows(); i++) {
                                if (spaceVar.get(nvar_sum) == outlinks.get(i, 0)) {
                                    idx = i;
                                    break;
                                }
                            }
                            
                            // Update to next outlink (round-robin)
                            if (idx < outlinks.getNumRows() - 1) {
                                spaceVar.set(nvar_sum, outlinks.get(idx + 1, 0));
                            } else {
                                spaceVar.set(nvar_sum, outlinks.get(0, 0));
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
                
                // outprob remains 1.0 for single selected action (scalar)
            }
        }
        
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}