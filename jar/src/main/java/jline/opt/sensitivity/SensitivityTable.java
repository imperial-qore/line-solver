package jline.opt.sensitivity;

import jline.lang.Network;
import jline.opt.results.SensitivityData;

/**
 * Analytic performance sensitivities for a product-form model, mirroring
 * native-Python {@code compute_model_sensitivities}. Computes exact
 * d(metric)/d(rate) for open product-form networks (closed-form single-server /
 * delay derivatives) and closed single-server unit-visit networks (via the
 * differentiated-MVA primitive {@code Pfqn_sens}); returns {@code null} for
 * mixed, multiserver or non-unit-visit topologies, so the optimizer falls back
 * to finite differences.
 */
public class SensitivityTable {

    private SensitivityTable() {
    }

    /**
     * Analytic sensitivities for the model, or {@code null} when unavailable.
     * The open and closed product-form branches are implemented in
     * {@link OpenSensitivity} and {@link ClosedSensitivity}.
     */
    public static SensitivityData compute(Network model) {
        return compute(model, false);
    }

    /**
     * Analytic sensitivities for the model, or {@code null} when unavailable.
     * The open and closed product-form branches are fast but narrow: they return
     * {@code null} for open, multiserver, or non-unit-visit networks. When they do
     * and {@code useCTMC} is true, the exact generator-derivative fallback
     * {@link CtmcSensitivity} is tried instead, which is limited by the state-space
     * size rather than by product-form assumptions and produces only QLen. Mirrors
     * MATLAB {@code opt.sens.computeModelSensitivities}.
     */
    public static SensitivityData compute(Network model, boolean useCTMC) {
        SensitivityData sens = null;
        try {
            SensitivityData open = OpenSensitivity.compute(model);
            if (open != null) {
                return open;
            }
            sens = ClosedSensitivity.compute(model);
        } catch (RuntimeException e) {
            sens = null;
        }
        if (sens == null && useCTMC) {
            try {
                sens = CtmcSensitivity.compute(model);
            } catch (RuntimeException e) {
                sens = null;
            }
        }
        return sens;
    }
}
