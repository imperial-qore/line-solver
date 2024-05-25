package jline.api;

import jline.lang.NetworkStruct;
import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;

import java.util.HashMap;
import java.util.Map;

import static jline.lib.M3A.*;

/**
 * APIs for approximating non-product-form queueing networks
 */
public class NPFQN {
    /**
     * Merges MMAP traffic flows
     */
    // TODO: different style used here to specify input configuration
    public static Map<Integer, Matrix> npfqn_traffic_merge(Map<Integer, Map<Integer, Matrix>> MMAPa, String config_merge_, String config_compress_) {
        String merge = "default";
        String compress = "default";
        if (config_merge_.equals("super") || config_merge_.equals("mixture") || config_merge_.equals("interpos")) {
            merge = config_merge_;
        } else {
            throw new RuntimeException("Unsupported configuration for merge.");
        }

        int n = MMAPa.size();
        Map<Integer, Matrix> SMMAP = new HashMap<>();
        //TODO: incomplete implementation
        return mmap_normalize(SMMAP);
    }

    /**
     * Merges MMAP traffic flows with class switching with default parameters
     */
    public static Map<Integer, Matrix> npfqn_traffic_merge_cs(Map<Integer, Map<Integer, Matrix>> MMAPs, Matrix prob) {
        return npfqn_traffic_merge_cs(MMAPs, prob, "default");
    }

    /**
     * Merges MMAP traffic flows with class switching
     */
    public static Map<Integer, Matrix> npfqn_traffic_merge_cs(Map<Integer, Map<Integer, Matrix>> MMAPs, Matrix prob, String config) {
        int n = MMAPs.size();
        int R = prob.getNumCols();
        Map<Integer, Map<Integer, Matrix>> MMAP_copy = new HashMap<>();
        for (Map.Entry<Integer, Map<Integer, Matrix>> entry : MMAPs.entrySet()) {
            Integer key = entry.getKey();
            Map<Integer, Matrix> innerMap = entry.getValue();

            Map<Integer, Matrix> deepCopyInnerMap = new HashMap<>();

            for (Map.Entry<Integer, Matrix> innerEntry : innerMap.entrySet()) {
                Integer innerKey = innerEntry.getKey();
                Matrix innerValue = innerEntry.getValue();

                // Perform a deep copy of innerValue if needed

                deepCopyInnerMap.put(innerKey, innerValue);
            }

            MMAP_copy.put(key, deepCopyInnerMap);
        }
        for (int i = 0; i < n; i++) {
            Matrix P = new Matrix(R, R, R * R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    P.set(r, s, prob.get((i - 1) * R + r, s));
                }
            }
            MMAP_copy.put(i, mmap_mark(MMAP_copy.get(i), P));
        }
        Map<Integer, Matrix> SMMAP = new HashMap<>();
        if (n == 1) {
            SMMAP = MMAPs.get(0);
        } else {
            if (config.equals("default") || config.equals("super")) {
                SMMAP = MMAPs.get(0);
                for (int j = 1; j < n; j++) {
                    SMMAP = mmap_super(SMMAP, MMAP_copy.get(j), "match");
                }
            }
        }
        return SMMAP;
    }

    /**
     * Splits MMAP traffic flows with class switching
     * @param MMAP - a MMAP to be split
     * @param P - a class switching matrix
     * @return a split MMAP
     */
    public static Map<Integer, Map<Integer, Matrix>> npfqn_traffic_split_cs(Map<Integer, Matrix> MMAP, Matrix P) {
        int n = MMAP.size();
        int R = P.getNumRows();
        int J = P.getNumCols();
        int M = Math.round(J / (float) R);
        Map<Integer, Map<Integer, Matrix>> SMMAP = new HashMap<>();
        for (int j = 0; j < M; j++) {
            SMMAP.put(j, new HashMap<>());
            SMMAP.get(j).put(0, MMAP.get(0).add(1, MMAP.get(1)));
            SMMAP.get(j).put(1, new Matrix(MMAP.get(1).getNumRows(), MMAP.get(1).getNumCols()));
            for (int s = 0; s < R; s++) {
                SMMAP.get(j).put(2 + s, new Matrix(SMMAP.get(j).get(0).getNumRows(), SMMAP.get(j).get(0).getNumCols()));
                for (int r = 0; r < R; r++) {
                    Matrix a = MMAP.get(2 + r).clone();
                    a.scale(P.get(r, (j - 1) * R + s));
                    SMMAP.get(j).put(2 + s, SMMAP.get(j).get(2 + s).add(1, a));
                    SMMAP.get(j).put(1, SMMAP.get(j).get(1).add(1, a));
                    SMMAP.get(j).put(0, SMMAP.get(j).get(0).add(1, a));
                }
            }
            SMMAP.put(j, mmap_normalize(SMMAP.get(j)));
        }
        return SMMAP;
    }

    public static npfqnNonexpApproxReturn npfqn_nonexp_approx(String method, NetworkStruct sn, Matrix ST, Matrix V, Matrix SCV, Matrix T, Matrix U, Matrix gamma, Matrix nservers) {
        int M = sn.nstations;
        Matrix rho = new Matrix(M, 1);
        rho.zero();
        Matrix scva = rho.clone();
        Matrix scvs = rho.clone();
        Matrix eta = rho.clone();
        switch (method) {
            case "default": case "none": case "hmva":
                return new npfqnNonexpApproxReturn(ST, gamma, nservers, rho, scva, scvs, eta);
            case "interp":
                for (int i=0; i<M; i++) {
                    Matrix nnzClasses = new Matrix(1, ST.getNumCols());
                    nnzClasses.zero();
                    for (int j=0; j<ST.getNumCols(); j++) {
                        if (Double.isFinite(ST.get(i, j)) && Double.isFinite(SCV.get(i, j))) {
                            nnzClasses.set(i, j, 1);
                        }
                    }
                    for (int j=0; j<nnzClasses.getNumElements(); j++) {
                        if (nnzClasses.get(j) > 0) {
                            rho.set(i, rho.get(i) + U.get(i, j));
                        }
                    }
                    if (nnzClasses.elementSum() != 0) {
                        switch (sn.sched.get(sn.stations.get(i))) {
                            case FCFS:
                                Matrix STinnz = new Matrix(1, (int) nnzClasses.elementSum());
                                Matrix SCVinnz = STinnz.clone();
                                Matrix Tinnz = STinnz.clone();
                                int tempj = 0;
                                for (int j=0; j<nnzClasses.getNumElements(); j++) {
                                    if (nnzClasses.get(j) > 0) {
                                        STinnz.set(tempj, ST.get(i, j));
                                        SCVinnz.set(tempj, ST.get(i, j));
                                        T.set(tempj, ST.get(i, j));
                                        tempj++;
                                    }
                                }
                                if (STinnz.elementMax() - STinnz.elementMin() > 0 || SCVinnz.elementMax() > 1 + GlobalConstants.FineTol || SCVinnz.elementMin() < 1 - GlobalConstants.FineTol) {
                                    scva.set(i, 1);
                                    scvs.set(i, SCVinnz.mult(Tinnz.transpose()).toDouble() / Tinnz.elementSum());
                                    gamma.set(i, (Math.pow(rho.get(i), nservers.get(i)) + rho.get(i)) / 2);
                                    if (scvs.get(i) > 1-1e-6 && scvs.get(i) < 1+1e-6 && nservers.get(i) == 1) {
                                        eta.set(i, rho.get(i));
                                    } else {
                                        eta.set(i, Math.exp(-2 * (1 - rho.get(i)) / (scvs.get(i) + scva.get(i) * rho.get(i))));
                                    }
                                    int order = 8;
                                    double ai = Math.pow(rho.get(i), order);
                                    double bi = Math.pow(rho.get(i), order);
                                    for (int k=0; k<nnzClasses.getNumElements(); k++) {
                                        if (nnzClasses.get(k) > 0 && sn.rates.get(i, k) > 0) {
                                            ST.set(i, k, Math.max(0, 1-ai) * ST.get(i, k) + ai * (bi * eta.get(i) + Math.max(0, 1 - bi) * gamma.get(i)) * (nservers.get(i) / Tinnz.elementSum()));
                                        }
                                    }
                                    nservers.set(i, 1);
                                }
                                // Continue from line 20
                        }
                    }
                }
        }
        return new npfqnNonexpApproxReturn(ST, gamma, nservers, rho, scva, scvs, eta);
    }

    public static class npfqnNonexpApproxReturn {
        public Matrix ST;
        public Matrix gamma;
        public Matrix nservers;
        public Matrix rho;
        public Matrix scva;
        public Matrix scvs;
        public Matrix eta;

        public npfqnNonexpApproxReturn(Matrix ST, Matrix gamma, Matrix nservers, Matrix rho, Matrix scva, Matrix scvs, Matrix eta) {
            this.ST = ST;
            this.gamma = gamma;
            this.nservers = nservers;
            this.rho = rho;
            this.scva = scva;
            this.scvs = scvs;
            this.eta = eta;
        }
    }

}
