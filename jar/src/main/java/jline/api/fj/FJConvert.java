/**
 * @file Convert LINE distributions to FJ_codes format
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Map_lambda;
import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.lib.fjcodes.FJArrival;
import jline.lib.fjcodes.FJService;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class FJConvert {
    private FJConvert() {}

    /**
     * Extract Fork-Join parameters from network structure
     */
    public static Pair<List<FJArrival>, List<FJService>> extractFJParams(NetworkStruct sn, FJInfo fjInfo) {
        List<FJArrival> arrivals = new ArrayList<FJArrival>();
        List<FJService> services = new ArrayList<FJService>();

        for (int r = 0; r < sn.nclasses; r++) {
            Station sourceStation = sn.stations.get((int) sn.nodeToStation.get(fjInfo.getSourceIdx()));
            JobClass jobClassR = sn.jobclasses.get(r);
            MatrixCell arrivalProc = (sn.proc != null && sn.proc.get(sourceStation) != null)
                    ? sn.proc.get(sourceStation).get(jobClassR) : null;
            double lambda = sn.rates.get(fjInfo.getSourceIdx(), r);

            FJArrival arrival;
            if (arrivalProc != null && arrivalProc.size() >= 2 && arrivalProc.get(0).getNumRows() > 1) {
                arrival = convertToFJArrival(arrivalProc.get(0), arrivalProc.get(1));
            } else {
                Matrix l0 = new Matrix(1, 1);
                l0.set(0, 0, -lambda);
                Matrix l1 = new Matrix(1, 1);
                l1.set(0, 0, lambda);
                arrival = new FJArrival(lambda, l0, l1, 1);
            }
            arrivals.add(arrival);

            int firstQueueIdx = fjInfo.getQueueIndices()[0];
            Station queueStation = sn.stations.get((int) sn.nodeToStation.get(firstQueueIdx));
            MatrixCell serviceProc = (sn.proc != null && sn.proc.get(queueStation) != null)
                    ? sn.proc.get(queueStation).get(jobClassR) : null;
            ProcessType procType = (sn.procid != null && sn.procid.get(queueStation) != null)
                    ? sn.procid.get(queueStation).get(jobClassR) : null;
            double mu_val = 1.0;
            if (sn.mu != null && sn.mu.get(queueStation) != null && sn.mu.get(queueStation).get(jobClassR) != null) {
                mu_val = sn.mu.get(queueStation).get(jobClassR).get(0, 0);
            }

            FJService service;
            if (serviceProc != null && serviceProc.size() >= 2
                    && procType != null && procType != ProcessType.EXP
                    && serviceProc.get(0).getNumRows() > 1) {
                service = convertToFJServiceFromProc(serviceProc);
            } else {
                Matrix ST = new Matrix(1, 1);
                ST.set(0, 0, -mu_val);
                Matrix St = new Matrix(1, 1);
                St.set(0, 0, mu_val);
                Matrix tau = new Matrix(1, 1);
                tau.set(0, 0, 1.0);
                service = new FJService(mu_val, ST, St, tau, 1);
            }
            services.add(service);

            double meanServiceRate;
            if (serviceProc != null && serviceProc.size() >= 2 && !serviceProc.get(0).hasNaN()) {
                meanServiceRate = Map_lambda.map_lambda(serviceProc.get(0), serviceProc.get(1));
            } else {
                meanServiceRate = mu_val;
            }
            if (lambda >= meanServiceRate) {
                throw new IllegalStateException("Class " + r + " unstable: arrival rate " + lambda
                        + " >= service rate " + meanServiceRate);
            }
        }

        return new Pair<List<FJArrival>, List<FJService>>(arrivals, services);
    }

    /**
     * Convert LINE MAP to FJ arrival format
     */
    public static FJArrival convertToFJArrival(Matrix D0, Matrix D1) {
        int n = D0.getNumRows();

        double lambda = Map_lambda.map_lambda(D0, D1);

        return new FJArrival(lambda, D0, D1, (n == 1) ? 1 : 2);
    }

    /**
     * Convert LINE PH to FJ service format
     */
    public static FJService convertToFJService(Matrix S, Matrix s, Matrix tau) {
        int n = S.getNumRows();

        Matrix negSinv = S.scale(-1.0).inv();
        Matrix ones = Matrix.ones(n, 1);
        double meanServiceTime = tau.mult(negSinv).mult(ones).get(0, 0);
        double mu = 1.0 / meanServiceTime;

        return new FJService(mu, S, s, tau, (n == 1) ? 1 : 2);
    }

    /**
     * Convert LINE process cell {D0, D1} to FJ service format.
     */
    public static FJService convertToFJServiceFromProc(MatrixCell proc) {
        Matrix D0 = proc.get(0);
        Matrix D1 = proc.get(1);
        int n = D0.getNumRows();

        Matrix ones = Matrix.ones(n, 1);
        Matrix s = D0.scale(-1.0).mult(ones);

        Matrix tau = Map_pie.map_pie(D0, D1);

        Matrix negSinv = D0.scale(-1.0).inv();
        double meanServiceTime = tau.mult(negSinv).mult(ones).get(0, 0);
        double mu = 1.0 / meanServiceTime;

        return new FJService(mu, D0, s, tau, (n == 1) ? 1 : 2);
    }
}
