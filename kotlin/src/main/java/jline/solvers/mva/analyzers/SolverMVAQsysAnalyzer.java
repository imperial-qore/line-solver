package jline.solvers.mva.analyzers;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;
import org.apache.commons.lang3.NotImplementedException;

/**
 * MVA Query System analyzer
 */
public class SolverMVAQsysAnalyzer implements MVAAnalyzer{
    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        long startTime = System.currentTimeMillis();
        Matrix QN = new Matrix(3,1);
        Matrix UN = new Matrix(3,1);
        Matrix RN = new Matrix(3,1);
        Matrix TN = new Matrix(3,1);
        Matrix CN = new Matrix(3,1);
        Matrix AN = new Matrix(3,1);
        Matrix WN = new Matrix(3,1);
        Matrix XN = new Matrix(3,1);
        double lG = Double.NaN;
        int it = 1;

        String method = options.method;
        int source_ist = -1, queue_ist = -1;
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Source){
                source_ist = (int) sn.nodeToStation.get(i);
            } else if(sn.nodetypes.get(i) == NodeType.Queue){
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }
        double lambda = sn.rates.get(source_ist) * sn.visits.get(0).get((int) sn.stationToStateful.get(queue_ist));
        int k = (int) sn.nservers.get(queue_ist);
        double mu = sn.rates.get(queue_ist);
        double ca = Math.sqrt(sn.scv.get(source_ist));
        double cs = Math.sqrt(sn.scv.get(queue_ist));
        if(method.equals("exact")){
            if(ca == 1 && cs == 1 && k == 1){
                method = "mm1";
            } else if(ca == 1 && cs == 1 && k > 1){
                method = "mmk";
            } else if(ca == 1 && k == 1){
                method = "mg1";
            } else if(cs == 1 && k == 1){
                method = "gm1";
            } else {
                throw new RuntimeException("MVA exact method unavailable for this model.");
            }
        }
        if (method.equals("default")) {
            if (k > 1) {
                method = "gigk";
            } else {
                method = "gig1.klb";
            }
        }
        double R = 0;
        switch(method){
            case "mm1":
                R = qsys_mm1(lambda, mu).W;
                break;
            case "mmk":
                R = qsys_mmk(lambda,mu,k).W;
                break;
            case "mg1": case "mgi1":  // verified
                R = qsys_mg1(lambda,mu,cs).W;
                break;
            case "gigk":
                R = qsys_gigk_approx(lambda,mu,ca,cs,k).W;
                break;
            case "gigk.kingman_approx":
                R = qsys_gigk_approx_kingman(lambda,mu,ca,cs,k).W;
                break;
            case "gig1": case "gig1.kingman":  // verified
                R = qsys_gig1_ubnd_kingman(lambda,mu,ca,cs).W;
                break;
            case "gig1.heyman":
                R = qsys_gig1_approx_heyman(lambda,mu,ca,cs).W;
                break;
            case "gig1.allen":
                R = qsys_gig1_approx_allencunneen(lambda,mu,ca,cs).W;
                break;
            case "gig1.kobayashi":
                R = qsys_gig1_approx_kobayashi(lambda,mu,ca,cs).W;
                break;
            case "gig1.klb":
                R = qsys_gig1_approx_klb(lambda,mu,ca,cs).W;
                if(options.method.equals("default")){
                    method = "default [gig1.klb]";
                }
                break;
            case "gig1.marchal": // verified
                R = qsys_gig1_approx_marchal(lambda,mu,ca,cs).W;
                break;
            case "gm1": case "gim1":
                // sigma = Load at arrival instants (Laplace transform of the inter-arrival times)
                // TODO: implement Fzero
                throw new NotImplementedException("Fzero is not currently implemented in Java");
            default:
                throw new RuntimeException("Unsupported method for a model with 1 station and 1 class.");
        }
        RN.set(queue_ist, 0, R * sn.visits.get(0).get((int) sn.stationToStateful.get(queue_ist)));
        CN.set(queue_ist, 0, RN.get(0,0));
        XN.set(queue_ist, 0, lambda);
        UN.set(queue_ist, 0, lambda/mu/k);
        TN.set(source_ist, 0, lambda);
        TN.set(queue_ist, 0, lambda);
        QN.set(queue_ist, 0, XN.get(queue_ist, 0) * RN.get(queue_ist, 0));
        lG = 0;
        long endTime = System.currentTimeMillis();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000.0;
        res.iter = it;
    }


    public qsysReturn qsys_mm1(double lambda, double mu){
        double rho = lambda/mu;
        double W = rho/(1-rho) / lambda;
        return new qsysReturn(W, rho);
    }

    public qsysReturn qsys_mmk(double lambda, double mu, int k){
        double rho = lambda/mu/k;
        double Q = rho/(1 - rho)*ErlangC(k, rho) + k * rho;
        double W = Q / lambda;
        return new qsysReturn(W, rho);
    }

    public qsysReturn qsys_mg1(double lambda, double mu, double cs){
        double rho = lambda/mu;
        double Q = rho + rho * rho / (2 * (1 - rho)) + lambda * lambda * cs * cs / (mu * mu)/(2 * (1 - rho));
        double W = Q/lambda;
        double rhohat = Q/(1 + Q);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gigk_approx(double lambda, double mu, double ca, double cs, int k){
        double rho = lambda / (mu * k);
        double alpha = 0;
        if(rho > 0.7){
            alpha = (Math.pow(rho, k) + rho) / 2;
        } else {
            alpha = Math.pow(rho, (k+1)/2.0);
        }
        double W = (alpha/mu) * (1/(1 - rho)) * (ca * ca + cs * cs)/(2*k) + 1/mu;
        double rhohat = W * lambda/(1 + W * lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gigk_approx_kingman(double lambda, double mu, double ca, double cs, int k){
        double W = (ca * ca + cs * cs) / 2 * (qsys_mmk(lambda, mu, k).W - 1/mu) + 1/mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_ubnd_kingman(double lambda, double mu, double ca, double cs){
        double rho = lambda/mu;
        double W = rho/(1-rho)*(ca * ca + cs * cs)/2*(1/mu) + (1/mu);
        double rhohat = W*lambda/(1+W*lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_approx_heyman(double lambda, double mu, double ca, double cs){
        double rho = lambda/mu;
        double W = rho/(1 - rho) / mu * (ca * ca + cs * cs) / 2 + 1 / mu;
        double rhohat = W*lambda/(1+W*lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_approx_allencunneen(double lambda, double mu, double ca, double cs){
        double rho = lambda/mu;
        double W = (rho / (1 - rho)) / mu * ((cs * cs + ca * ca)/2) + 1/mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_approx_kobayashi(double lambda, double mu, double ca, double cs){
        double rho = lambda/mu;
        double rhohat = Math.exp(-2 * (1 - rho) / (rho * (ca * ca + cs * cs / rho)));
        double W = rhohat / (1 - rhohat) / lambda;
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_approx_klb(double lambda, double mu, double ca, double cs){
        // Kramer-Langenbach-Belz formula
        double rho = lambda/mu;
        double g;
        if(ca <= 1){
            g = Math.exp(-2 * (1 - rho) * Math.pow(1 - ca * ca, 2) / (3 * rho * (ca * ca + cs * cs)));
        } else {
            g = Math.exp(-(1 - rho) * (ca * ca - 1) / (ca * ca + 4 * cs * cs));
        }
        double W = 1 / mu * ((rho / (1 - rho)) * ((cs * cs + ca * ca) / 2) * g + 1);
        double rhohat = W * lambda / (1 + W * lambda);
        return new qsysReturn(W, rhohat);
    }

    public qsysReturn qsys_gig1_approx_marchal(double lambda, double mu, double ca, double cs){
        double rho = lambda/mu;
        double Wmm1 = rho/(1 - rho);
        double W = Wmm1 * (1 + cs * cs) / 2 / mu * (ca + rho * rho * cs * cs)/(1 + rho * rho * cs * cs) + 1/mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new qsysReturn(W, rhohat);
    }


    /**
     * The probability that an arriving customer is forced to join the queue (i.e. all servers are occupied)
     * @param k - the number of servers
     * @param rho - utilization
     * @return - probability that an arriving customer is forced to join the queue (i.e. all servers are occupied)
     */
    public double ErlangC(int k, double rho){
        double S = 0;
        int factj = 1;
        int factj_1 = 1;
        for(int j = 0; j <= k-1; j++){
            if(j == 0){
                S += Math.pow(k * rho, j);
            } else {
                factj = j * factj_1;
                S += Math.pow(k * rho, j)/factj;
                factj_1 = factj;
            }
        }
        return 1 / (1 + (1 - rho) * (k * factj_1) / Math.pow(k * rho, k) * S);
    }

    public static class qsysReturn {
        public double W;
        public double rho;

        public qsysReturn(double W, double rho) {
            this.W = W;
            this.rho = rho;
        }
    }
}
