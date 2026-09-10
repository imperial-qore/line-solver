package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * Build PH representation for 2-node Fork-Join job.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class BuildServiceH {
    private BuildServiceH() {}

    /**
     * Build service representation for a 2-node FJ job.
     *
     * This function constructs the Phase-Type (PH) representation for the service
     * time of a 2-node Fork-Join job, which is the maximum of two independent
     * service times.
     *
     * @param service Single subtask service process
     * @return FJServiceH with PH representation for 2-node job
     */
    public static FJServiceH build_Service_h(FJService service) {
        int dim_single = service.tau_st.length();
        Matrix phases_single = FJUtils.build_index(dim_single, 1);

        // Possible service phases for the 2-node FJ job
        // Each row represents a state: [longest queue phases, shortest queue phases]
        int totalStates = dim_single * dim_single;
        Matrix service_phases = new Matrix(totalStates, 2 * dim_single);

        int k = 0;
        for (int i = 0; i < dim_single; i++) {
            for (int j = 0; j < dim_single; j++) {
                // [longest, shortest]
                for (int col = 0; col < dim_single; col++) {
                    service_phases.set(k, col, phases_single.get(i, col));
                    service_phases.set(k, dim_single + col, phases_single.get(j, col));
                }
                k++;
            }
        }

        // PH representation for the service time of a 2-node FJ job
        // beta = kron(tau_st, tau_st)
        Matrix beta = service.tau_st.kron(service.tau_st).transpose();  // Row vector

        // S = kronsum(ST, ST) = ST \otimes I + I \otimes ST
        Matrix S = FJUtils.kronsum(service.ST, service.ST);

        return new FJServiceH(service_phases, beta, S);
    }
}
