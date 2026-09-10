package jline.api.sn;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Collects the (D0,D1) modulation records of every non-renewal arrival or
 * service process declared in the network.
 *
 * <p>A MAP with matrices (D0,D1) is a Poisson-like point process modulated by
 * the CTMC with generator Q = D0 + D1 (the phase process), whose conditional
 * intensity in phase k is lambda(k) = sum_j D1(k,j). This returns one record
 * per modulating process, so that a solver-agnostic transformation can replace
 * each of them by a random-environment stage set (see jline.io.MAPQN2RENV).
 *
 * <p>Only processes declared as MAP, MMPP2 or MMAP are reported: every other
 * distribution is stored in sn.proc in (D0,D1) form as well (Erlang, Coxian,
 * APH, ...), but those are renewal processes that carry no modulation and are
 * supported natively by the phase-type solvers.
 *
 * <p>Marked processes (MMAP) at a Source are reported as a single record whose
 * classes list holds every marked class, since all marks share one phase
 * process; the per-class intensity comes from the mark-specific D1 matrices.
 *
 * <p>Mirrors matlab/src/api/sn/sn_map_modulation.m.
 *
 * @since LINE 3.0
 */
public final class SnMapModulation {
    private SnMapModulation() {}

    /** One modulating process of the network. */
    public static final class MapModulation {
        /** Station index of the modulated process. */
        public final int ist;
        /** Node index of that station. */
        public final int node;
        /** True when the process is an arrival process, false for a service process. */
        public final boolean arrival;
        /** Class indices governed by this process (several only for a marked MAP). */
        public final List<Integer> classes;
        /** Hidden-phase generator part D0. */
        public final Matrix D0;
        /** Per-class D1 blocks, aligned with {@link #classes}. */
        public final List<Matrix> D1;
        /** Number of phases of the modulating chain. */
        public final int order;
        /** True when every D1 block is diagonal, i.e. the process is an MMPP. */
        public final boolean isMMPP;

        MapModulation(int ist, int node, boolean arrival, List<Integer> classes,
                      Matrix D0, List<Matrix> D1, int order, boolean isMMPP) {
            this.ist = ist;
            this.node = node;
            this.arrival = arrival;
            this.classes = classes;
            this.D0 = D0;
            this.D1 = D1;
            this.order = order;
            this.isMMPP = isMMPP;
        }

        /** @return the phase-conditional intensity of class index c (position in classes) in phase k */
        public double intensity(int c, int k) {
            Matrix D1c = D1.get(c);
            double rate = 0;
            for (int j = 0; j < D1c.getNumCols(); j++) {
                rate += D1c.get(k, j);
            }
            return rate;
        }

        /** @return the generator D0 + sum_c D1_c of the phase process */
        public Matrix phaseGenerator() {
            Matrix Q = D0.copy();
            for (Matrix D1c : D1) {
                Q = Q.add(1.0, D1c);
            }
            return Q;
        }
    }

    /**
     * @param sn the NetworkStruct object for the queueing network model
     * @return one record per MAP/MMPP2/MMAP process, empty when the model has none
     */
    public static List<MapModulation> snMapModulation(NetworkStruct sn) {
        List<MapModulation> mods = new LinkedList<>();
        if (sn == null || sn.procid == null || sn.proc == null) {
            return mods;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            Station st = sn.stations.get(ist);
            int nd = (int) sn.stationToNode.get(ist);
            boolean arrival = sn.nodetype.get(nd) == NodeType.Source;
            Map<JobClass, ProcessType> pidMap = sn.procid.get(st);
            Map<JobClass, MatrixCell> procMap = sn.proc.get(st);
            if (pidMap == null || procMap == null) {
                continue;
            }
            boolean[] done = new boolean[sn.nclasses];
            for (int r = 0; r < sn.nclasses; r++) {
                if (done[r]) {
                    continue;
                }
                ProcessType pt = pidMap.get(sn.jobclasses.get(r));
                if (pt != ProcessType.MAP && pt != ProcessType.MMPP2 && pt != ProcessType.MMAP) {
                    continue;
                }
                MatrixCell mapIr = mapOf(procMap, sn.jobclasses.get(r));
                if (mapIr == null) {
                    continue;
                }
                int markOfR = markIndex(sn, ist, r);
                if (markOfR > 0) {
                    // MMAP: one phase process shared by every marked class of the
                    // station, one D1 block per mark
                    List<Integer> marked = new ArrayList<>();
                    int carrier = -1;
                    int minMark = Integer.MAX_VALUE;
                    for (int rr = 0; rr < sn.nclasses; rr++) {
                        int mk = markIndex(sn, ist, rr);
                        if (mk > 0) {
                            marked.add(rr);
                            if (mk < minMark) {
                                minMark = mk;
                                carrier = rr;
                            }
                        }
                    }
                    MatrixCell mapC = mapOf(procMap, sn.jobclasses.get(carrier));
                    if (mapC == null || mapC.size() < 2 + marked.size()) {
                        throw new RuntimeException(String.format(
                                "The marked arrival process at station %d carries %d mark matrices for %d marked "
                                        + "classes; the (D0,D1,D1^(1),...,D1^(C)) form is required.",
                                ist + 1, mapC == null ? 0 : Math.max(0, mapC.size() - 2), marked.size()));
                    }
                    List<Matrix> D1c = new ArrayList<>();
                    boolean isMMPP = true;
                    for (int k = 0; k < marked.size(); k++) {
                        Matrix D1k = mapC.get(1 + markIndex(sn, ist, marked.get(k)));
                        D1c.add(D1k);
                        isMMPP = isMMPP && isDiagonal(D1k);
                    }
                    mods.add(new MapModulation(ist, nd, arrival, marked, mapC.get(0), D1c,
                            mapC.get(0).getNumRows(), isMMPP));
                    for (int rr : marked) {
                        done[rr] = true;
                    }
                } else {
                    List<Integer> classes = new ArrayList<>();
                    classes.add(r);
                    List<Matrix> D1c = new ArrayList<>();
                    D1c.add(mapIr.get(1));
                    mods.add(new MapModulation(ist, nd, arrival, classes, mapIr.get(0), D1c,
                            mapIr.get(0).getNumRows(), isDiagonal(mapIr.get(1))));
                    done[r] = true;
                }
            }
        }
        return mods;
    }

    /** Per station-class (D0,D1,...) representation, null when absent or disabled. */
    private static MatrixCell mapOf(Map<JobClass, MatrixCell> procMap, JobClass jobclass) {
        MatrixCell mapIr = procMap.get(jobclass);
        if (mapIr == null || mapIr.size() < 2 || mapIr.get(0) == null || mapIr.get(1) == null) {
            return null;
        }
        if (mapIr.get(0).hasNaN()) {
            return null; // disabled
        }
        return mapIr;
    }

    /** 1-based mark index of class r at station ist, or -1 when the class is not marked. */
    private static int markIndex(NetworkStruct sn, int ist, int r) {
        if (sn.markidx == null || sn.markidx.isEmpty()) {
            return -1;
        }
        if (ist >= sn.markidx.getNumRows() || r >= sn.markidx.getNumCols()) {
            return -1;
        }
        return (int) sn.markidx.get(ist, r);
    }

    /** True when D1 has no off-diagonal mass, i.e. the process is an MMPP. */
    private static boolean isDiagonal(Matrix D1) {
        double off = 0;
        for (int i = 0; i < D1.getNumRows(); i++) {
            for (int j = 0; j < D1.getNumCols(); j++) {
                if (i != j) {
                    off += Math.abs(D1.get(i, j));
                }
            }
        }
        return off <= 1e-14 * Math.max(1.0, D1.norm());
    }
}
