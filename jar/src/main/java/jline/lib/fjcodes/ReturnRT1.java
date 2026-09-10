/**
 * Compute response time percentiles for K=1 (single queue case).
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.HashMap;
import java.util.Map;

import jline.io.InputOutput;
import jline.lib.butools.MMAPPH1FCFS;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class ReturnRT1 {
    private ReturnRT1() {}

    /**
     * Compute response time percentiles for K=1 Fork-Join queue.
     */
    public static Matrix returnRT1(FJArrival arrival, FJService service, double[] pers) {
        // For K=1, the Fork-Join queue reduces to a single MAP/PH/1 queue.
        MatrixCell D = new MatrixCell(2);
        D.set(0, arrival.lambda0);
        D.set(1, arrival.lambda1);

        Map<Integer, Matrix> sigma = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> S = new HashMap<Integer, Matrix>();

        sigma.put(0, service.tau_st.transpose());
        S.put(0, service.ST);

        try {
            Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                    D, sigma, S,
                    null, null, null, null,
                    false, true, 1e-14, null);

            Map<Integer, Matrix> stDistrPH_alpha = result.get("stDistrPH_alpha");
            Map<Integer, Matrix> stDistrPH_A = result.get("stDistrPH_A");

            if (stDistrPH_alpha == null || stDistrPH_A == null) {
                throw new RuntimeException("MMAPPH1FCFS did not return sojourn time PH distribution");
            }

            Matrix res_alpha = stDistrPH_alpha.get(0);
            Matrix Smat = stDistrPH_A.get(0);

            if (res_alpha == null || Smat == null) {
                throw new RuntimeException("MMAPPH1FCFS returned null sojourn time distribution for class 0");
            }

            return ReturnPer.returnPer(res_alpha, Smat, pers);

        } catch (Exception e) {
            InputOutput.line_warning("returnRT1", "MMAPPH1FCFS failed, using simplified calculation: %s", e.getMessage());

            Matrix res_alpha = service.tau_st.transpose();
            Matrix Smat = service.ST;
            return ReturnPer.returnPer(res_alpha, Smat, pers);
        }
    }
}
