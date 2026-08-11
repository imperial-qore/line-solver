package jline.opt.sensitivity;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.opt.results.SensitivityData;
import jline.util.matrix.Matrix;

import java.util.Map;

/**
 * Analytic d(metric)/d(rate) for open product-form networks (single-server
 * queueing and infinite-server delay stations, which decouple). Mirrors
 * native-Python {@code _open_sensitivities}. Returns {@code null} for closed,
 * mixed, or multiserver networks (finite-difference fallback).
 *
 * <p>References: Z. Liu and P. Nain, INRIA RR-1144 (1989), Thm 3.2 for open
 * BCMP networks.</p>
 */
final class OpenSensitivity {

    private OpenSensitivity() {
    }

    static SensitivityData compute(Network model) {
        NetworkStruct sn = model.getStruct();
        int R = sn.nclasses;
        Matrix njobs = sn.njobs;
        if (njobs == null || njobs.length() == 0) {
            return null;
        }
        for (int r = 0; r < R; r++) {
            if (!Double.isInfinite(njobs.get(r))) {
                return null;   // open only: every class must be infinite-population
            }
        }

        int nstations = sn.nstations;
        Matrix rates = sn.rates;
        Matrix nservers = sn.nservers;

        // visits per (station, class): sum over chains
        double[][] visits = new double[nstations][R];
        boolean anyVisit = false;
        if (sn.visits != null && !sn.visits.isEmpty()) {
            for (Map.Entry<Integer, Matrix> e : sn.visits.entrySet()) {
                Matrix vm = e.getValue();
                for (int i = 0; i < nstations; i++) {
                    for (int r = 0; r < R; r++) {
                        visits[i][r] += vm.get(i, r);
                    }
                }
            }
            anyVisit = true;
        }
        if (!anyVisit) {
            for (int i = 0; i < nstations; i++) {
                for (int r = 0; r < R; r++) {
                    visits[i][r] = 1.0;
                }
            }
        }

        // arrival rate per class from EXT (Source) stations
        double[] lam = new double[R];
        for (int i = 0; i < nstations; i++) {
            if (schedOf(sn, i) == SchedStrategy.EXT) {
                for (int r = 0; r < R; r++) {
                    double rate = rates.get(i, r);
                    if (!Double.isNaN(rate) && !Double.isInfinite(rate)) {
                        lam[r] += rate;
                    }
                }
            }
        }

        SensitivityData sens = new SensitivityData();

        for (int i = 0; i < nstations; i++) {
            SchedStrategy sc = schedOf(sn, i);
            if (sc == SchedStrategy.EXT) {
                continue;
            }
            boolean isDelay = (sc == SchedStrategy.INF);
            if (!isDelay && nservers.get(i) > 1) {
                return null;   // multiserver M/M/c not handled here
            }
            String stName = stationNodeName(sn, i);

            double[] D = new double[R];
            double[] rho = new double[R];
            for (int r = 0; r < R; r++) {
                double mu = rates.get(i, r);
                if (!Double.isNaN(mu) && !Double.isInfinite(mu) && mu > 0 && visits[i][r] > 0) {
                    D[r] = visits[i][r] / mu;
                    rho[r] = lam[r] * D[r];
                }
            }
            double U = 0.0;
            if (!isDelay) {
                for (int r = 0; r < R; r++) {
                    U += rho[r];
                }
            }
            double denom = 1.0 - U;
            if (!isDelay && denom <= 0) {
                return null;   // unstable: sensitivity diverges
            }

            // Util_i = sum_s rho_s ; d U_i / d mu_s = -rho_s / mu_s
            for (int s = 0; s < R; s++) {
                double muS = rates.get(i, s);
                if (!Double.isNaN(muS) && !Double.isInfinite(muS) && muS > 0 && visits[i][s] > 0) {
                    sens.add("Util", stName, SensitivityData.paramKey(stName, className(sn, s)),
                            -rho[s] / muS);
                }
            }

            for (int r = 0; r < R; r++) {
                double mu = rates.get(i, r);
                if (visits[i][r] <= 0 || Double.isNaN(mu) || Double.isInfinite(mu) || mu <= 0) {
                    continue;
                }
                String cl = className(sn, r);
                String mkey = SensitivityData.metricKey(stName, cl);
                for (int s = 0; s < R; s++) {
                    double muS = rates.get(i, s);
                    if (Double.isNaN(muS) || Double.isInfinite(muS) || muS <= 0 || visits[i][s] <= 0) {
                        continue;
                    }
                    String pkey = SensitivityData.paramKey(stName, className(sn, s));
                    double dU = isDelay ? 0.0 : (-rho[s] / muS);          // d U_i / d mu_s
                    double dDr = (s == r) ? (-D[r] / muS) : 0.0;          // d D_r / d mu_s
                    double dRhoR = (s == r) ? (-rho[r] / muS) : 0.0;      // d rho_r / d mu_s
                    double dResp;
                    double dQ;
                    if (isDelay) {
                        dResp = dDr;
                        dQ = dRhoR;
                    } else {
                        dResp = (dDr * denom + D[r] * dU) / (denom * denom);
                        dQ = (dRhoR * denom + rho[r] * dU) / (denom * denom);
                    }
                    sens.add("RespT", mkey, pkey, dResp);
                    sens.add("QLen", mkey, pkey, dQ);
                }
                // Tput sensitivity to its own rate is 0 for open networks
                sens.add("Tput", mkey, SensitivityData.paramKey(stName, cl), 0.0);
            }
        }
        return sens;
    }

    private static SchedStrategy schedOf(NetworkStruct sn, int stationIndex) {
        if (sn.stations == null || stationIndex >= sn.stations.size()) {
            return null;
        }
        Station st = sn.stations.get(stationIndex);
        return sn.sched != null ? sn.sched.get(st) : null;
    }

    private static String stationNodeName(NetworkStruct sn, int stationIndex) {
        try {
            int nodeIdx = (int) sn.stationToNode.get(stationIndex);
            if (nodeIdx >= 0 && nodeIdx < sn.nodenames.size()) {
                return sn.nodenames.get(nodeIdx);
            }
        } catch (RuntimeException e) {
            // fall through
        }
        return Integer.toString(stationIndex);
    }

    private static String className(NetworkStruct sn, int r) {
        if (sn.classnames != null && r < sn.classnames.size()) {
            return sn.classnames.get(r);
        }
        return Integer.toString(r);
    }
}
