/**
 * @file M3PP interleaved multi-process fitting
 *
 * Implements fitting and interleaving of k second-order M3PP processes with varying
 * class counts. Constructs higher-order M3PP from multiple independent processes
 * through state space composition and parameter aggregation.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;

public final class M3pp_interleave_fitc {
    private M3pp_interleave_fitc() {}

    /**
     * Fits k second-order M3PP[m_j] and interleaves them into a M3PP[m] of order k+1.
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc(
            double[] av,
            double[] btv,
            double[] binfv,
            double[][] acc,
            double[][] gtcc,
            double t,
            double tinf,
            int[] mapping,
            boolean reorder) {

        int k = av.length;
        if (btv.length != k) throw new IllegalArgumentException("btv must have length k");
        if (binfv.length != k) throw new IllegalArgumentException("binfv must have length k");
        if (acc.length != k) throw new IllegalArgumentException("acc must have length k");
        if (gtcc.length != k) throw new IllegalArgumentException("gtcc must have length k");

        List<MatrixCell> componentM3pps = new ArrayList<MatrixCell>();
        List<ComponentSpec> componentSpecs = new ArrayList<ComponentSpec>();

        for (int i = 0; i < k; i++) {
            double[] classRates = acc[i];
            double[] classVarCov = gtcc[i];
            double processRate = av[i];

            ComponentSpec spec = new ComponentSpec(i, processRate, classRates, classVarCov, btv[i], binfv[i]);
            componentSpecs.add(spec);

            MatrixCell component = fitM3ppWithClassCharacteristics(
                    processRate, btv[i], binfv[i], classRates, classVarCov, t, tinf);

            componentM3pps.add(component);
        }

        MatrixCell interleaved = interleaveM3pps(componentSpecs, componentM3pps, mapping, reorder);
        return new Pair<MatrixCell, List<MatrixCell>>(interleaved, componentM3pps);
    }

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc(
            double[] av, double[] btv, double[] binfv,
            double[][] acc, double[][] gtcc, double t, double tinf) {
        return m3pp_interleave_fitc(av, btv, binfv, acc, gtcc, t, tinf, null, false);
    }

    private static class ComponentSpec {
        final int processIndex;
        final double processRate;
        final double[] classRates;
        final double[] classVarCov;
        final double idc_t;
        final double idc_inf;

        ComponentSpec(int processIndex, double processRate, double[] classRates,
                      double[] classVarCov, double idc_t, double idc_inf) {
            this.processIndex = processIndex;
            this.processRate = processRate;
            this.classRates = classRates;
            this.classVarCov = classVarCov;
            this.idc_t = idc_t;
            this.idc_inf = idc_inf;
        }
    }

    private static MatrixCell fitM3ppWithClassCharacteristics(
            double processRate, double bt, double binf,
            double[] classRates, double[] classVarCov, double t, double tinf) {

        int numClasses = classRates.length;
        Matrix[] mmpp = new Matrix[2];

        double lambda1 = processRate * (1 + bt) / 2;
        double lambda2 = processRate * (1 - bt) / 2;
        double mu1 = lambda1 / 10;
        double mu2 = lambda2 / 10;

        mmpp[0] = new Matrix(2, 2);
        mmpp[0].set(0, 0, -(lambda1 + mu1));
        mmpp[0].set(0, 1, mu1);
        mmpp[0].set(1, 0, mu2);
        mmpp[0].set(1, 1, -(lambda2 + mu2));

        mmpp[1] = new Matrix(2, 2);
        mmpp[1].set(0, 0, lambda1);
        mmpp[1].set(0, 1, 0.0);
        mmpp[1].set(1, 0, 0.0);
        mmpp[1].set(1, 1, lambda2);

        MatrixCell m3pp = new MatrixCell(2 + numClasses);
        m3pp.set(0, mmpp[0]);
        m3pp.set(1, mmpp[1]);

        double totalClassRate = 0.0;
        for (double cr : classRates) totalClassRate += cr;

        for (int i = 0; i < numClasses; i++) {
            double proportion = (totalClassRate > 0) ? classRates[i] / totalClassRate : 1.0 / numClasses;
            double varianceWeight = (classVarCov.length > 0) ? classVarCov[i] : 1.0;
            m3pp.set(2 + i, mmpp[1].scale(proportion * Math.max(0.1, varianceWeight)));
        }
        return m3pp;
    }

    private static MatrixCell interleaveM3pps(
            List<ComponentSpec> specs, List<MatrixCell> components,
            int[] mapping, boolean reorder) {
        if (components.isEmpty()) {
            throw new IllegalArgumentException("Cannot interleave empty list of M3PPs");
        }
        int totalClasses = 0;
        for (ComponentSpec spec : specs) {
            totalClasses += spec.classRates.length;
        }
        int interleavedOrder = Math.min(components.size() + 1, 10);
        MatrixCell interleaved = new MatrixCell(totalClasses + 2);
        Matrix D0 = new Matrix(interleavedOrder, interleavedOrder);
        Matrix D1 = new Matrix(interleavedOrder, interleavedOrder);
        interleaved.set(0, D0);
        interleaved.set(1, D1);
        for (int c = 0; c < totalClasses; c++) {
            interleaved.set(2 + c, new Matrix(interleavedOrder, interleavedOrder));
        }

        int globalClassIndex = 0;
        for (int compIdx = 0; compIdx < components.size(); compIdx++) {
            MatrixCell component = components.get(compIdx);
            ComponentSpec spec = specs.get(compIdx);
            int compClasses = spec.classRates.length;

            int[] stateMapping = createStateMapping(compIdx, components.size(), interleavedOrder);

            Matrix compD0 = component.get(0);
            for (int i = 0; i < Math.min(compD0.getNumRows(), interleavedOrder); i++) {
                for (int j = 0; j < Math.min(compD0.getNumCols(), interleavedOrder); j++) {
                    int globalI = stateMapping[i % stateMapping.length];
                    int globalJ = stateMapping[j % stateMapping.length];
                    if (globalI < interleavedOrder && globalJ < interleavedOrder) {
                        D0.set(globalI, globalJ, D0.get(globalI, globalJ) + compD0.get(i, j) / components.size());
                    }
                }
            }

            for (int c = 0; c < compClasses; c++) {
                Matrix compClass = component.get(2 + c);
                int targetClassIndex = (reorder && mapping != null)
                        ? mapping[globalClassIndex % mapping.length]
                        : globalClassIndex;

                if (targetClassIndex < totalClasses) {
                    Matrix interleavedClass = interleaved.get(2 + targetClassIndex);
                    for (int i = 0; i < Math.min(compClass.getNumRows(), interleavedOrder); i++) {
                        for (int j = 0; j < Math.min(compClass.getNumCols(), interleavedOrder); j++) {
                            int globalI = stateMapping[i % stateMapping.length];
                            int globalJ = stateMapping[j % stateMapping.length];
                            if (globalI < interleavedOrder && globalJ < interleavedOrder) {
                                double coordinationFactor = computeCoordinationFactor(spec, compIdx, components.size());
                                double arrivalRate = compClass.get(i, j) * coordinationFactor;
                                interleavedClass.set(globalI, globalJ, interleavedClass.get(globalI, globalJ) + arrivalRate);
                                D1.set(globalI, globalJ, D1.get(globalI, globalJ) + arrivalRate);
                            }
                        }
                    }
                }
                globalClassIndex++;
            }
        }
        addCoordinationTransitions(interleaved, specs);
        ensureGeneratorProperty(D0);
        return interleaved;
    }

    private static int[] createStateMapping(int componentIndex, int totalComponents, int interleavedOrder) {
        int[] mapping = new int[interleavedOrder];
        int offset = componentIndex;
        for (int i = 0; i < mapping.length; i++) {
            mapping[i] = (offset + i * totalComponents) % interleavedOrder;
        }
        return mapping;
    }

    private static double computeCoordinationFactor(ComponentSpec spec, int compIndex, int totalComponents) {
        double sumClass = 0.0;
        for (double v : spec.classRates) sumClass += v;
        double relativerate = spec.processRate / (sumClass + 1e-6);
        double positionFactor = (compIndex + 1.0) / totalComponents;
        return Math.min(1.0, relativerate * positionFactor);
    }

    private static void addCoordinationTransitions(MatrixCell interleaved, List<ComponentSpec> specs) {
        Matrix D0 = interleaved.get(0);
        int order = D0.getNumRows();
        double sum = 0.0;
        for (ComponentSpec s : specs) sum += s.processRate;
        double avg = specs.isEmpty() ? 0.0 : sum / specs.size();
        double coordinationRate = avg * 0.1;
        for (int i = 0; i < order - 1; i++) {
            for (int j = i + 1; j < order; j++) {
                D0.set(i, j, D0.get(i, j) + coordinationRate / order);
                D0.set(j, i, D0.get(j, i) + coordinationRate / order);
            }
        }
    }

    private static void ensureGeneratorProperty(Matrix D0) {
        for (int i = 0; i < D0.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (i != j) rowSum += D0.get(i, j);
            }
            D0.set(i, i, -rowSum);
        }
    }
}
