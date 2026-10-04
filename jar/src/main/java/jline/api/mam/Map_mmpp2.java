/**
 * @file MMPP(2) fitting from moments and autocorrelation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_mmpp2 {
    private Map_mmpp2() {}

    /**
     * Fits an MMPP(2) as a MAP from four descriptive parameters.
     *
     * Matches the requested mean, SCV, skewness and lag-1 autocorrelation
     * exactly. Throws when the request lies outside the MMPP(2) feasible set,
     * rather than returning a non-MAP.
     *
     * @param MEAN the mean inter-arrival time
     * @param SCV_param the squared coefficient of variation, at least 1
     * @param SKEW the skewness, or -1 for the minimum-skewness fit
     * @param ACF1 the lag-1 autocorrelation, or -1 for the maximum feasible one
     * @return the fitted MMPP(2), as {D0,D1}
     */
    public static MatrixCell map_mmpp2(double MEAN, double SCV_param, double SKEW, double ACF1) {
        double FEASTOL = Math.pow(10.0, (double) (-Map_feastol.map_feastol()));
        double E1 = MEAN;
        double E2 = (1.0 + SCV_param) * E1 * E1;
        double E1_2 = E1 * E1;
        double E1_3 = E1_2 * E1;
        double E3 = -(2.0 * E1_3 - 3.0 * E1 * E2 - SKEW * Math.pow(E2 - E1_2, 1.5));

        // The closed form below solves the moment-matching equations, but nothing
        // in it constrains the solution to be a MAP: outside the MMPP(2) feasible
        // set it returns negative rates, i.e. a D1 with negative entries and a D0
        // with a positive diagonal. Rowsums stay zero, so the usual generator
        // check does not catch it. Reject the request instead of returning a
        // non-MAP.
        if (SCV_param < 1.0 - FEASTOL) {
            throw new IllegalArgumentException(String.format(
                    "map_mmpp2: SCV=%g is infeasible, the inter-arrival times of an MMPP(2) "
                            + "are over-dispersed (SCV>=1).", SCV_param));
        }
        if (Math.abs(SCV_param - 1.0) <= FEASTOL) {
            throw new IllegalArgumentException(String.format(
                    "map_mmpp2: SCV=1 is the Poisson boundary, where the MMPP(2) fit is "
                            + "degenerate: the decay rate G2=ACF1/(1-1/SCV)/0.5 divides by zero and "
                            + "every rate comes back NaN. Use map_exponential(%g) for a Poisson process.",
                    MEAN));
        }

        // ACF1=RHO0MAX is attained in the limit of a decay rate G2->1; no MMPP(2)
        // exceeds it, and none is negatively autocorrelated
        double RHO0MAX = 0.5 * (1.0 - 1.0 / SCV_param);
        if (ACF1 != -1.0) {
            if (ACF1 < -FEASTOL) {
                throw new IllegalArgumentException(String.format(
                        "map_mmpp2: ACF1=%g is infeasible, an MMPP(2) cannot be negatively "
                                + "autocorrelated. Pass ACF1=-1 to request the maximum feasible "
                                + "autocorrelation.", ACF1));
            }
            if (ACF1 > RHO0MAX + FEASTOL) {
                throw new IllegalArgumentException(String.format(
                        "map_mmpp2: ACF1=%g exceeds the maximum lag-1 autocorrelation %g feasible "
                                + "at SCV=%g. Pass ACF1=-1 to request it.", ACF1, RHO0MAX, SCV_param));
            }
        }

        double G2;
        if (ACF1 == -1.0) {
            G2 = 1.0 - 10.0 * Math.pow(10.0, (double) (-Map_feastol.map_feastol()));
        } else {
            G2 = ACF1 / (1.0 - 1.0 / SCV_param) / 0.5;
        }

        double SCV = SCV_param;
        if (SKEW == -1.0 && SCV > 1.0) {
            E3 = (3.0 / 2.0 + 0.001) * E2 * E2 / E1;
        }

        // the SKEW==-1 branch above sits just above E3MIN, which is the infimum of
        // the third moment over the class
        double E3MIN = (3.0 / 2.0) * E2 * E2 / E1;
        if (E3 < E3MIN - FEASTOL) {
            throw new IllegalArgumentException(String.format(
                    "map_mmpp2: SKEW=%g gives E3=%g, below the minimum third moment %g feasible "
                            + "at SCV=%g. Pass SKEW=-1 for the minimum-skewness fit.",
                    SKEW, E3, E3MIN, SCV_param));
        }

        SCV = (E2 - E1_2) / E1_2;

        double mu00;
        double mu11;
        double q01;
        double q10;

        if (G2 < 1e-6) {
            double SCV2 = SCV * SCV;
            mu00 = 2.0 * (6.0 * E1_3 * SCV - E3) / E1 / (6.0 * E1_3 * SCV + 3.0 * E1_3 * SCV2 + 3.0 * E1_3 - 2.0 * E3);
            mu11 = 0.0;
            double E1_5 = E1_3 * E1_2;
            q01 = 9.0 * E1_5 * (SCV - 1.0) * (SCV2 - 2.0 * SCV + 1.0) / (6.0 * E1_3 * SCV - E3) / (6.0 * E1_3 * SCV + 3.0 * E1_3 * SCV2 + 3.0 * E1_3 - 2.0 * E3);
            q10 = -3.0 * (SCV - 1.0) * E1_2 / (6.0 * E1_3 * SCV - E3);
        } else {
            double E1_5 = E1_3 * E1_2;
            double E1_6 = E1_3 * E1_3;
            double SCV2 = SCV * SCV;
            double SCV3 = SCV2 * SCV;
            double G2_2 = G2 * G2;
            double G2_3 = G2_2 * G2;

            double DISC = E3 * E3 - 12.0 * E1_3 * SCV * E3 + 6.0 * E1_3 * G2 * E3 - 6.0 * G2 * SCV * E1_3 * E3
                    + 18.0 * G2 * SCV3 * E1_6 - 18.0 * E1_6 * G2 * SCV2 + 9.0 * E1_6 * G2_2 + 36.0 * E1_6 * SCV2
                    + 18.0 * E1_6 * G2 * SCV - 18.0 * E1_6 * SCV * G2_2 + 9.0 * E1_6 * SCV2 * G2_2 - 18.0 * E1_6 * G2;

            double A = -3.0 * E1_3 * G2 + 3.0 * E1_3 * G2 * SCV - 6.0 * E1_3 * SCV + E3 + Math.sqrt(DISC);
            double B = -3.0 * E1_3 * SCV2 - 6.0 * E1_3 * SCV - 3.0 * E1_3 + 2.0 * E3;
            double F = A / B;

            mu11 = F / E1;

            double mu00_numer = -4.0 * E3 * G2 + 4.0 * F * E3 * G2 - 18.0 * E1_3 * F * G2 - 18.0 * E1_3 * F * G2 * SCV2
                    - 12.0 * E1_3 * G2_2 - 12.0 * E1_3 * F * G2_2 * SCV + 12.0 * E1_3 * F * G2 * SCV
                    + 12.0 * E1_3 * G2 * SCV2 - 9.0 * E1_3 * F * SCV + 3.0 * E1_3 * F + 12.0 * E1_3 * G2_2 * SCV
                    + 9.0 * E1_3 * F * SCV2 + 12.0 * E1_3 * G2 + 12.0 * E1_3 * F * G2_2 - 3.0 * E1_3 * F * SCV3;

            double mu00_denom = 12.0 * E1_3 * G2_3 * SCV + 3.0 * E1_3 * SCV3 * G2 - 12.0 * E1_3 * G2_3
                    + 18.0 * E1_3 * G2_2 * SCV2 - 3.0 * E1_3 * G2 + 27.0 * E1_3 * F * G2 * SCV2
                    - 9.0 * E1_3 * G2 * SCV2 + 18.0 * E1_3 * G2_2 - 12.0 * E1_3 * G2_2 * SCV
                    + 9.0 * E1_3 * G2 * SCV - 12.0 * E1_3 * F * G2_3 * SCV - 9.0 * E1_3 * F * SCV3 * G2
                    - 24.0 * E1_3 * F * G2_2 * SCV2 - F * E3 * SCV2 + 4.0 * F * E3 * G2_2 + 12.0 * E1_3 * F * G2_3
                    - F * E3 + 2.0 * F * E3 * SCV + 9.0 * E1_3 * F * G2 + 24.0 * E1_3 * F * G2_2 * SCV
                    - 27.0 * E1_3 * F * G2 * SCV + 6.0 * E1_3 * F * SCV - 12.0 * E1_3 * F * SCV2
                    - 24.0 * E1_3 * F * G2_2 + 6.0 * E1_3 * F * SCV3 - 4.0 * E3 * G2_2;

            mu00 = G2 * mu00_numer / mu00_denom / E1;

            double q01_numer = -6.0 * F * E1_2 * SCV + 12.0 * F * E1_2 * G2 * SCV - 6.0 * G2 * SCV * E1_2
                    - 3.0 * F * E1_2 * G2 + mu11 * E3 + 3.0 * E1_2 * G2 + 6.0 * F * E1_2 * SCV2
                    - 9.0 * F * E1_2 * SCV2 * G2 + 3.0 * E1_2 * G2 * SCV2 - E3 * mu11 * SCV
                    - 6.0 * F * E1_2 * G2_2 * SCV + 6.0 * E1_2 * G2_2 * SCV + 3.0 * F * E1_2 * G2_2
                    - G2 * mu11 * E3 - 3.0 * E1_2 * G2_2 + 3.0 * F * E1_2 * SCV2 * G2_2
                    - 3.0 * E1_2 * SCV2 * G2_2 + G2 * SCV * mu11 * E3;

            double q01_denom = -45.0 * F * E1_5 * G2 * SCV2 + 18.0 * G2_2 * E1_5 * SCV + 18.0 * E1_5 * G2_3
                    - 27.0 * E1_5 * G2_2 * SCV2 + 6.0 * E1_2 * G2_2 * E3 - 27.0 * E1_5 * G2_2
                    - 18.0 * E1_5 * G2_3 * SCV - 18.0 * E1_5 * G2 * SCV + 18.0 * E1_5 * G2 * SCV2
                    + 3.0 * E1_2 * G2 * E3 - 3.0 * E1_2 * G2 * E3 * SCV + mu11 * E3 * E3
                    + 3.0 * F * E1_2 * G2 * SCV * E3 - 36.0 * F * E1_5 * G2_2 * SCV
                    + 36.0 * F * E1_5 * G2_2 + 36.0 * F * E1_5 * SCV2 + 45.0 * F * E1_5 * G2 * SCV
                    - 12.0 * F * E1_2 * SCV * E3 - 3.0 * F * E1_2 * G2 * E3
                    + 9.0 * F * E1_5 * G2 * SCV3 + 36.0 * F * E1_5 * G2_2 * SCV2
                    - 6.0 * F * E1_2 * G2_2 * E3 + 18.0 * F * E1_5 * G2_3 * SCV
                    - 18.0 * F * E1_5 * G2_3 - 9.0 * F * E1_5 * G2;

            q01 = -3.0 * E1_2 * q01_numer / q01_denom;

            double q10_inner = -3.0 * E1_3 * F * SCV3 - 3.0 * E1_3 * F * G2 * SCV2 + 6.0 * E1_3 * SCV2
                    + 3.0 * E1_3 * G2 * SCV2 + 3.0 * E1_3 * F * SCV2 + 6.0 * E1_3 * F * G2 * SCV
                    - E3 * SCV - 6.0 * E1_3 * SCV + F * E3 * SCV - 6.0 * E1_3 * G2 * SCV
                    - 3.0 * E1_3 * F * SCV - F * E3 + 3.0 * E1_3 * G2 - 3.0 * E1_3 * F * G2
                    + 3.0 * E1_3 * F + E3;

            q10 = 3.0 * q10_inner * E1_2 * (-1.0 + G2) / DISC;
        }

        // Catch-all: the checks above cover the known infeasible directions, but
        // the authoritative test is the solution itself. An MMPP(2) has
        // non-negative rates by construction, so anything else is not a MAP and
        // must not be returned.
        double[] rates = new double[]{mu00, mu11, q01, q10};
        for (int i = 0; i < rates.length; i++) {
            if (Double.isNaN(rates[i]) || rates[i] < -FEASTOL) {
                throw new IllegalArgumentException(String.format(
                        "map_mmpp2: (MEAN=%g,SCV=%g,SKEW=%g,ACF1=%g) is not MMPP(2)-feasible: the "
                                + "fit gives [mu00 mu11 q01 q10]=[%g %g %g %g], which is not a MAP.",
                        MEAN, SCV_param, SKEW, ACF1, mu00, mu11, q01, q10));
            }
        }
        // a request on the feasibility boundary lands on a zero rate up to
        // roundoff: clear that sign flip, but leave small POSITIVE rates alone --
        // ACF1=-1 asks for G2=1-1e-7, whose near-uncoupled chain has legitimate
        // rates around 1e-9
        if (mu00 < 0) mu00 = 0;
        if (mu11 < 0) mu11 = 0;
        if (q01 < 0) q01 = 0;
        if (q10 < 0) q10 = 0;

        Matrix D0 = new Matrix(2, 2, 2);
        D0.set(0, 0, -mu00 - q01);
        D0.set(0, 1, q01);
        D0.set(1, 0, q10);
        D0.set(1, 1, -mu11 - q10);

        Matrix D1 = new Matrix(2, 2, 2);
        D1.set(0, 0, mu00);
        D1.set(0, 1, 0.0);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, mu11);

        MatrixCell MAP = new MatrixCell();
        MAP.set(0, D0);
        MAP.set(1, D1);
        return MAP;
    }
}
