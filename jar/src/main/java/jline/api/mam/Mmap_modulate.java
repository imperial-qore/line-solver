package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import java.util.function.Function;

public final class Mmap_modulate {
    private Mmap_modulate() {}

    /**
     * Modulates an MMAP by another MMAP, creating a compound arrival process.
     * The modulating process controls the rate of the modulated process.
     *
     * @param baseMmap          The base MMAP to be modulated
     * @param modulatingMmap    The MMAP that provides modulation
     * @param modulationFactor  Strength of modulation
     * @return Modulated MMAP
     */
    public static MatrixCell mmap_modulate(MatrixCell baseMmap, MatrixCell modulatingMmap,
                                           double modulationFactor) {
        int baseOrder = baseMmap.get(0).getNumRows();
        int modulatingOrder = modulatingMmap.get(0).getNumRows();
        int baseClasses = baseMmap.size() - 2;
        int modulatingClasses = modulatingMmap.size() - 2;

        // Create compound state space (Kronecker product structure)
        int compoundOrder = baseOrder * modulatingOrder;
        int totalClasses = baseClasses; // Keep base classes, modulation affects rates

        MatrixCell modulated = new MatrixCell(totalClasses + 2);

        // Initialize matrices
        Matrix D0 = new Matrix(compoundOrder, compoundOrder);
        Matrix D1 = new Matrix(compoundOrder, compoundOrder);
        modulated.set(0, D0);
        modulated.set(1, D1);

        for (int c = 0; c < totalClasses; c++) {
            modulated.set(c + 2, new Matrix(compoundOrder, compoundOrder));
        }

        // Extract matrices from input MMAPs
        Matrix baseD0 = baseMmap.get(0);
        Matrix modulatingD0 = modulatingMmap.get(0);

        // Compute modulating process steady-state probabilities
        double[] modulatingPi = computeSteadyState(modulatingMmap);

        // Compute modulation rates for each state of modulating process
        double[] modulationRates = new double[modulatingOrder];
        for (int i = 0; i < modulatingOrder; i++) {
            double totalRate = 0.0;
            for (int c = 1; c < modulatingMmap.size(); c++) {
                Matrix Dc = modulatingMmap.get(c);
                for (int j = 0; j < modulatingOrder; j++) {
                    totalRate += Dc.get(i, j);
                }
            }
            modulationRates[i] = Math.max(0.1, totalRate * modulationFactor);
        }

        // Build compound generator matrix
        for (int i = 0; i < baseOrder; i++) {
            for (int j = 0; j < baseOrder; j++) {
                for (int k = 0; k < modulatingOrder; k++) {
                    for (int l = 0; l < modulatingOrder; l++) {
                        int compoundI = i * modulatingOrder + k;
                        int compoundJ = j * modulatingOrder + l;

                        if (compoundI < compoundOrder && compoundJ < compoundOrder) {
                            if (i == j) {
                                // Modulating process transitions (base state unchanged)
                                D0.set(compoundI, compoundJ, D0.get(compoundI, compoundJ) + modulatingD0.get(k, l));
                            } else if (k == l) {
                                // Base process transitions (modulating state unchanged)
                                D0.set(compoundI, compoundJ, D0.get(compoundI, compoundJ) + baseD0.get(i, j) * modulationRates[k]);
                            }
                        }
                    }
                }
            }
        }

        // Build class arrival matrices with modulation
        for (int c = 0; c < totalClasses; c++) {
            Matrix baseDc = baseMmap.get(c + 2);
            Matrix modulatedDc = modulated.get(c + 2);

            for (int i = 0; i < baseOrder; i++) {
                for (int j = 0; j < baseOrder; j++) {
                    for (int k = 0; k < modulatingOrder; k++) {
                        for (int l = 0; l < modulatingOrder; l++) {
                            int compoundI = i * modulatingOrder + k;
                            int compoundJ = j * modulatingOrder + l;

                            if (compoundI < compoundOrder && compoundJ < compoundOrder && k == l) {
                                // Base arrivals modulated by current modulating state
                                double modulatedRate = baseDc.get(i, j) * modulationRates[k];
                                modulatedDc.set(compoundI, compoundJ, modulatedRate);
                                D1.set(compoundI, compoundJ, D1.get(compoundI, compoundJ) + modulatedRate);
                            }
                        }
                    }
                }
            }
        }

        return modulated;
    }

    /**
     * Modulates an MMAP by another MMAP using the default modulation factor 1.0.
     */
    public static MatrixCell mmap_modulate(MatrixCell baseMmap, MatrixCell modulatingMmap) {
        return mmap_modulate(baseMmap, modulatingMmap, 1.0);
    }

    /**
     * Time-varying modulation of an MMAP.
     * The modulation factor varies according to a specified pattern.
     *
     * @param baseMmap          The base MMAP to modulate
     * @param modulationPattern Function that maps time to modulation factor
     * @param timeHorizon       Time horizon for discretization
     * @param numSteps          Number of time steps
     * @return Time-modulated MMAP (approximated as mixture)
     */
    public static MatrixCell mmap_modulate_time_varying(MatrixCell baseMmap,
                                                        Function<Double, Double> modulationPattern,
                                                        double timeHorizon,
                                                        int numSteps) {
        // Create mixture of MMAPs for different time periods
        List<MatrixCell> timeStepMmaps = new ArrayList<MatrixCell>();
        double[] weights = new double[numSteps];
        for (int i = 0; i < numSteps; i++) weights[i] = 1.0 / numSteps;

        for (int step = 0; step < numSteps; step++) {
            double time = timeHorizon * step / numSteps;
            double modulationFactor = modulationPattern.apply(time);

            // Create modulated MMAP for this time period
            Matrix scaleFactor = Matrix.singleton(1.0 / modulationFactor); // Scale mean inter-arrival time
            MatrixCell modulated = Mmap_scale.mmap_scale(baseMmap, scaleFactor);
            if (modulated != null) {
                timeStepMmaps.add(modulated);
            }
        }

        // Return mixture representing time-varying modulation
        return Mmap_mixture_order2.mmap_mixture_order2(timeStepMmaps, weights);
    }

    /**
     * Time-varying modulation with default time horizon (10.0) and number of steps (10).
     */
    public static MatrixCell mmap_modulate_time_varying(MatrixCell baseMmap,
                                                        Function<Double, Double> modulationPattern) {
        return mmap_modulate_time_varying(baseMmap, modulationPattern, 10.0, 10);
    }

    /**
     * Cross-modulation between two MMAPs.
     * Each MMAP modulates the other, creating mutual influence.
     *
     * @param mmap1                     First MMAP
     * @param mmap2                     Second MMAP
     * @param crossModulationStrength   Strength of cross-modulation
     * @return Cross-modulated MMAP system
     */
    public static MatrixCell mmap_cross_modulate(MatrixCell mmap1, MatrixCell mmap2,
                                                 double crossModulationStrength) {
        // First, mmap1 modulates mmap2
        MatrixCell mmap2ModulatedBy1 = mmap_modulate(mmap2, mmap1, crossModulationStrength);

        // Then, mmap2 modulates mmap1
        MatrixCell mmap1ModulatedBy2 = mmap_modulate(mmap1, mmap2, crossModulationStrength);

        // Combine using superposition with equal weights
        MatrixCell result = Mmap_super.mmap_super(mmap1ModulatedBy2, mmap2ModulatedBy1);
        if (result == null) {
            throw new RuntimeException("Failed to create superposition");
        }
        return result;
    }

    /**
     * Cross-modulation between two MMAPs with default strength 0.5.
     */
    public static MatrixCell mmap_cross_modulate(MatrixCell mmap1, MatrixCell mmap2) {
        return mmap_cross_modulate(mmap1, mmap2, 0.5);
    }

    /**
     * Compute steady-state probabilities for an MMAP (simplified)
     */
    private static double[] computeSteadyState(MatrixCell mmap) {
        int order = mmap.get(0).getNumRows();
        double[] pi = new double[order];
        for (int i = 0; i < order; i++) pi[i] = 1.0 / order; // Uniform distribution as approximation

        // In practice, would solve pi * Q = 0 where Q is the generator matrix
        // This is a simplified version
        return pi;
    }
}
