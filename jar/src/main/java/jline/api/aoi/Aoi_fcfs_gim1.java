/**
 * @file GI/M/1 FCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_fcfs_gim1 {
    private Aoi_fcfs_gim1() {}

    /**
     * Mean AoI and LST for GI/M/1 FCFS queue.
     *
     * <p>Exact mean: E[A] = lambda*E[Y^2]/2 + 1/mu + lambda*(-Y*'(eta))/eta
     * with eta = mu*(1-sigma), sigma the root of Y*(mu - mu*sigma) = sigma;
     * E[Apeak] = E[Y] + 1/(mu*(1-sigma)).</p>
     */
    public static AoiLstResult aoi_fcfs_gim1(final LstFunction Y_lst, final double mu,
                                             final double E_Y, final double E_Y2) {
        if (mu <= 0) throw new IllegalArgumentException("Service rate mu must be positive");
        if (E_Y <= 0) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");
        if (E_Y2 < E_Y * E_Y) throw new IllegalArgumentException("Second moment E_Y2 must be >= E_Y^2");

        double lambda = 1.0 / E_Y;
        double rho = lambda / mu;
        if (rho >= 1) {
            throw new IllegalArgumentException("System unstable: rho = 1/(E_Y*mu) = "
                    + String.format("%.4f", rho) + " >= 1");
        }

        final double sigma = findSigmaGIM1(Y_lst, mu);

        double E_D = 1.0 / (mu * (1.0 - sigma));
        // see _kb/03-api-layer.md for rationale
        double eta = mu * (1.0 - sigma);
        double hstep = 1e-6 * Math.max(1.0, eta);
        double dYstar = (Y_lst.evaluate(eta + hstep) - Y_lst.evaluate(eta - hstep)) / (2.0 * hstep);
        double meanAoI = lambda * E_Y2 / 2.0 + 1.0 / mu + lambda * (-dYstar) / eta;
        double peakAoI = E_Y + E_D;

        // LST of AoI (Inoue et al. 2019, Theorem 3), via the general age formula
        //   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
        // the cycle average of exp(-s*age) over a departure interval. In GI/M/1
        // the system time is EXPONENTIAL at rate eta = mu*(1-sigma), so
        // T*(s) = eta/(s+eta), and Lindley gives W' + Y = max(Y, T) with
        // T ~ Exp(eta) independent of the next interarrival Y, so
        //   E[exp(-s*max(Y,T))] = Y*(s) - (s/(s+eta)) * Y*(s+eta),
        // and the peak adds one fresh Exp(mu) service. A*(0) = 1 follows from
        // the defining relation Y*(eta) = sigma.
        //
        // THE PREVIOUS FORM WAS NOT AN LST: (mu*sigma(s))/(s+mu-mu*sigma(s))*D*(s)
        // gives sigma/(1-sigma) at s = 0 rather than 1, and it re-solved
        // sigma(s) by bisection at every point with a SILENT fallback to
        // sigma(0). Checked against simulation on E2/M/1: at s = 0.2 the old
        // form gave 0.23056, the form below 0.59319, and the sample path 0.59332.
        final double etaF = eta;
        final double lambdaF = lambda;
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) {
                    return 1.0; // A*(0) = 1 for any proper LST
                }
                double T_s = etaF / (s + etaF);
                double peak_s = (mu / (s + mu))
                        * (Y_lst.evaluate(s) - (s / (s + etaF)) * Y_lst.evaluate(s + etaF));
                return (lambdaF / s) * (T_s - peak_s);
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }

    private static double findSigmaGIM1(LstFunction Y_lst, double mu) {
        double lo = 0.001;
        double hi = 0.999;
        int maxIter = 100;
        double tol = 1e-12;

        double fLo = Y_lst.evaluate(mu - mu * lo) - lo;
        double fHi = Y_lst.evaluate(mu - mu * hi) - hi;

        if (fLo * fHi < 0) {
            for (int iter = 0; iter < maxIter; iter++) {
                double mid = (lo + hi) / 2.0;
                double fMid = Y_lst.evaluate(mu - mu * mid) - mid;
                if (Math.abs(fMid) < tol) {
                    return mid;
                }
                double fLoNew = Y_lst.evaluate(mu - mu * lo) - lo;
                if (fLoNew * fMid < 0) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return (lo + hi) / 2.0;
        } else {
            double sigma = 0.5;
            for (int iter = 0; iter < maxIter; iter++) {
                double sigmaNew = Y_lst.evaluate(mu - mu * sigma);
                if (Math.abs(sigmaNew - sigma) < tol) {
                    return sigmaNew;
                }
                sigma = sigmaNew;
            }
            return sigma;
        }
    }

    private static double findSigmaGIM1Shifted(LstFunction Y_lst, double mu, double s, double sigma0) {
        double lo = 0.001;
        double hi = 0.999;
        int maxIter = 100;
        double tol = 1e-12;

        double fLo = Y_lst.evaluate(s + mu - mu * lo) - lo;
        double fHi = Y_lst.evaluate(s + mu - mu * hi) - hi;

        if (fLo * fHi < 0) {
            for (int iter = 0; iter < maxIter; iter++) {
                double mid = (lo + hi) / 2.0;
                double fMid = Y_lst.evaluate(s + mu - mu * mid) - mid;
                if (Math.abs(fMid) < tol) {
                    return mid;
                }
                double fLoNew = Y_lst.evaluate(s + mu - mu * lo) - lo;
                if (fLoNew * fMid < 0) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return (lo + hi) / 2.0;
        } else {
            double sigma = sigma0;
            for (int iter = 0; iter < maxIter; iter++) {
                double sigmaNew = Y_lst.evaluate(s + mu - mu * sigma);
                if (Math.abs(sigmaNew - sigma) < tol) {
                    return sigmaNew;
                }
                sigma = sigmaNew;
            }
            return sigma;
        }
    }
}
