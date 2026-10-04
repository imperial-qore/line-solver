/**
 * Markov Modulated Poisson Process two-state fitting algorithms.
 *
 * The general branch is transcribed from the MATLAB closed form in
 * matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m. That closed form is exact: the
 * fitted MAP reproduces E1, E2, E3 and the autocorrelation decay rate
 * G2 = rho(k+1)/rho(k) identically, as verified symbolically in SageMath by
 * reducing each moment identity modulo SQ^2 - DISC. The expressions are
 * machine-generated from the MATLAB text; do not rewrite them by hand.
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

public final class Mmpp2_fit {

    private Mmpp2_fit() {}

    /**
     * Closed-form mu00 (arrival rate in the first environment state) of the MMPP(2) matching (E1, E2, E3, G2).
     *
     * Transcribed from matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m, whose closed
     * form is exact: the fitted MAP reproduces E1, E2, E3 and the decay rate
     * G2 = rho(k+1)/rho(k) identically (verified symbolically in SageMath).
     * SQ is the square root of the discriminant shared by the four rates.
     */
    public static double mmpp2_fit_mu00(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / (E1 * E1);
        double SQ = FastMath.sqrt(FastMath.pow(E3, 2)-12*FastMath.pow(E1, 3)*SCV*E3+6*FastMath.pow(E1, 3)*G2*E3-6*G2*SCV*FastMath.pow(E1, 3)*E3+18*G2*FastMath.pow(SCV, 3)*
            FastMath.pow(E1, 6)-18*FastMath.pow(E1, 6)*G2*FastMath.pow(SCV, 2)+9*FastMath.pow(E1, 6)*FastMath.pow(G2, 2)+36*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 6)*G2*SCV-18*FastMath.pow(E1, 6)*SCV*FastMath.pow(G2, 2)+9*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 6)*G2);
        return G2*(-4*E3*G2+4*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3*G2-18*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*G2-18*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*FastMath.pow(SCV, 2)-
            12*FastMath.pow(E1, 3)*FastMath.pow(G2, 2)-12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-
            6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(G2, 2)*SCV+12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*SCV+12*
            FastMath.pow(E1, 3)*G2*FastMath.pow(SCV, 2)-9*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-
            6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*SCV+3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)+12*FastMath.pow(E1, 3)*
            FastMath.pow(G2, 2)*SCV+9*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 2)+
            12*FastMath.pow(E1, 3)*G2+12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)-
            3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 3))/(12*FastMath.pow(E1, 3)*
            FastMath.pow(G2, 3)*SCV+3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 3)*G2-12*FastMath.pow(E1, 3)*FastMath.pow(G2, 3)+18*FastMath.pow(E1, 3)*
            FastMath.pow(G2, 2)*FastMath.pow(SCV, 2)-3*FastMath.pow(E1, 3)*G2+27*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2*FastMath.pow(SCV, 2)-9*FastMath.pow(E1, 3)*G2*FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 3)*FastMath.pow(G2, 2)-
            12*FastMath.pow(E1, 3)*FastMath.pow(G2, 2)*SCV+9*FastMath.pow(E1, 3)*G2*SCV-12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 3)*SCV-9*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(SCV, 3)*G2-24*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*
            FastMath.pow(SCV, 2)-(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3*FastMath.pow(SCV, 2)+4*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*E3*FastMath.pow(G2, 2)+12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(G2, 3)-(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*
            FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3+2*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*E3*SCV+9*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*
            FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2+24*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*
            SCV-27*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*
            FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*SCV+6*FastMath.pow(E1, 3)*
            (-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*SCV-12*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(SCV, 2)-24*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)+
            6*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 3)-4*E3*FastMath.pow(G2, 2))/
            E1;
    }

    /**
     * Closed-form mu11 (arrival rate in the second environment state) of the MMPP(2) matching (E1, E2, E3, G2).
     *
     * Transcribed from matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m, whose closed
     * form is exact: the fitted MAP reproduces E1, E2, E3 and the decay rate
     * G2 = rho(k+1)/rho(k) identically (verified symbolically in SageMath).
     * SQ is the square root of the discriminant shared by the four rates.
     */
    public static double mmpp2_fit_mu11(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / (E1 * E1);
        double SQ = FastMath.sqrt(FastMath.pow(E3, 2)-12*FastMath.pow(E1, 3)*SCV*E3+6*FastMath.pow(E1, 3)*G2*E3-6*G2*SCV*FastMath.pow(E1, 3)*E3+18*G2*FastMath.pow(SCV, 3)*
            FastMath.pow(E1, 6)-18*FastMath.pow(E1, 6)*G2*FastMath.pow(SCV, 2)+9*FastMath.pow(E1, 6)*FastMath.pow(G2, 2)+36*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 6)*G2*SCV-18*FastMath.pow(E1, 6)*SCV*FastMath.pow(G2, 2)+9*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 6)*G2);
        return (-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3);
    }

    /**
     * Closed-form q01 (environment transition rate out of the first state) of the MMPP(2) matching (E1, E2, E3, G2).
     *
     * Transcribed from matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m, whose closed
     * form is exact: the fitted MAP reproduces E1, E2, E3 and the decay rate
     * G2 = rho(k+1)/rho(k) identically (verified symbolically in SageMath).
     * SQ is the square root of the discriminant shared by the four rates.
     */
    public static double mmpp2_fit_q01(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / (E1 * E1);
        double SQ = FastMath.sqrt(FastMath.pow(E3, 2)-12*FastMath.pow(E1, 3)*SCV*E3+6*FastMath.pow(E1, 3)*G2*E3-6*G2*SCV*FastMath.pow(E1, 3)*E3+18*G2*FastMath.pow(SCV, 3)*
            FastMath.pow(E1, 6)-18*FastMath.pow(E1, 6)*G2*FastMath.pow(SCV, 2)+9*FastMath.pow(E1, 6)*FastMath.pow(G2, 2)+36*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 6)*G2*SCV-18*FastMath.pow(E1, 6)*SCV*FastMath.pow(G2, 2)+9*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 6)*G2);
        return -3*FastMath.pow(E1, 2)*(-6*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*SCV+12*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*SCV-6*G2*SCV*FastMath.pow(E1, 2)-3*(-3*FastMath.pow(E1, 3)*G2+
            3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2+(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3+3*FastMath.pow(E1, 2)*
            G2+6*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 2)-9*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 2)*G2+3*FastMath.pow(E1, 2)*G2*FastMath.pow(SCV, 2)-
            E3*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*SCV-6*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*SCV+6*FastMath.pow(E1, 2)*FastMath.pow(G2, 2)*SCV+3*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)-G2*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*E3-3*FastMath.pow(E1, 2)*FastMath.pow(G2, 2)+3*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-3*FastMath.pow(E1, 2)*FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)+G2*SCV*(-
            3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3)/(-45*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-
            6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*FastMath.pow(SCV, 2)+18*FastMath.pow(G2, 2)*FastMath.pow(E1, 5)*SCV+18*FastMath.pow(E1, 5)*
            FastMath.pow(G2, 3)-27*FastMath.pow(E1, 5)*FastMath.pow(G2, 2)*FastMath.pow(SCV, 2)+6*FastMath.pow(E1, 2)*FastMath.pow(G2, 2)*
            E3-27*FastMath.pow(E1, 5)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 5)*FastMath.pow(G2, 3)*SCV-18*FastMath.pow(E1, 5)*G2*
            SCV+18*FastMath.pow(E1, 5)*G2*FastMath.pow(SCV, 2)+3*FastMath.pow(E1, 2)*G2*E3-3*FastMath.pow(E1, 2)*G2*E3*SCV+(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/E1/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(E3, 2)+3*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2*SCV*E3-36*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*
            SCV+36*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*
            FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)+36*(-
            3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 2)+45*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*SCV-12*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-
            6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*SCV*E3-3*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2*E3+9*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*FastMath.pow(SCV, 3)+
            36*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*FastMath.pow(SCV, 2)-6*(-
            3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 2)/(-3*FastMath.pow(E1, 3)*
            FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 2)*E3+18*(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 3)*SCV-18*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(G2, 3)-9*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)*FastMath.pow(E1, 5)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2);
    }

    /**
     * Closed-form q10 (environment transition rate out of the second state) of the MMPP(2) matching (E1, E2, E3, G2).
     *
     * Transcribed from matlab/lib/kpctoolbox/mmpp/mmpp2_fit3.m, whose closed
     * form is exact: the fitted MAP reproduces E1, E2, E3 and the decay rate
     * G2 = rho(k+1)/rho(k) identically (verified symbolically in SageMath).
     * SQ is the square root of the discriminant shared by the four rates.
     */
    public static double mmpp2_fit_q10(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / (E1 * E1);
        double SQ = FastMath.sqrt(FastMath.pow(E3, 2)-12*FastMath.pow(E1, 3)*SCV*E3+6*FastMath.pow(E1, 3)*G2*E3-6*G2*SCV*FastMath.pow(E1, 3)*E3+18*G2*FastMath.pow(SCV, 3)*
            FastMath.pow(E1, 6)-18*FastMath.pow(E1, 6)*G2*FastMath.pow(SCV, 2)+9*FastMath.pow(E1, 6)*FastMath.pow(G2, 2)+36*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 6)*G2*SCV-18*FastMath.pow(E1, 6)*SCV*FastMath.pow(G2, 2)+9*FastMath.pow(E1, 6)*
            FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 6)*G2);
        return 3*(-3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*
            FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*FastMath.pow(SCV, 3)-3*FastMath.pow(E1, 3)*
            (-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-
            6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*FastMath.pow(SCV, 2)+6*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)+
            3*FastMath.pow(E1, 3)*G2*FastMath.pow(SCV, 2)+3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*
            SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*FastMath.pow(SCV, 2)+6*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*
            SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*G2*SCV-E3*
            SCV-6*FastMath.pow(E1, 3)*SCV+(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*E3*SCV-6*FastMath.pow(E1, 3)*
            G2*SCV-3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-
            3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)*SCV-(-3*FastMath.pow(E1, 3)*
            G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*
            SCV-3*FastMath.pow(E1, 3)+2*E3)*E3+3*FastMath.pow(E1, 3)*G2-3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*
            G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/(-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+
            2*E3)*G2+3*FastMath.pow(E1, 3)*(-3*FastMath.pow(E1, 3)*G2+3*FastMath.pow(E1, 3)*G2*SCV-6*FastMath.pow(E1, 3)*SCV+E3+SQ)/
            (-3*FastMath.pow(E1, 3)*FastMath.pow(SCV, 2)-6*FastMath.pow(E1, 3)*SCV-3*FastMath.pow(E1, 3)+2*E3)+E3)*FastMath.pow(E1, 2)*
            (-1+G2)/(FastMath.pow(E3, 2)-12*FastMath.pow(E1, 3)*SCV*E3+6*FastMath.pow(E1, 3)*G2*E3-6*G2*SCV*FastMath.pow(E1, 3)*E3+
            18*G2*FastMath.pow(SCV, 3)*FastMath.pow(E1, 6)-18*FastMath.pow(E1, 6)*G2*FastMath.pow(SCV, 2)+9*FastMath.pow(E1, 6)*FastMath.pow(G2, 2)+
            36*FastMath.pow(E1, 6)*FastMath.pow(SCV, 2)+18*FastMath.pow(E1, 6)*G2*SCV-18*FastMath.pow(E1, 6)*SCV*FastMath.pow(G2, 2)+
            9*FastMath.pow(E1, 6)*FastMath.pow(SCV, 2)*FastMath.pow(G2, 2)-18*FastMath.pow(E1, 6)*G2);
    }
    /**
     * Fits a 2-phase Markov modulated Poisson process (MMPP(2)) to the first three
     * moments E1, E2, E3 and the autocorrelation decay rate G2 = rho(k+1)/rho(k).
     */
    public static MatrixCell mmpp2_fit(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / E1 / E1;
        double mu00;
        double mu11;
        double q01;
        double q10;
        if (G2 < 1e-6) {
            mu00 = 2 * (6 * FastMath.pow(E1, 3) * SCV - E3) / E1 / (6 * FastMath.pow(E1, 3) * SCV + 3 * FastMath.pow(E1,
                3) * SCV * SCV + 3 * FastMath.pow(E1, 3) - 2 * E3);
            mu11 = 0.0;
            q01 = 9 * FastMath.pow(E1, 5) * (SCV - 1) * (SCV * SCV - 2 * SCV + 1) / (6 * FastMath.pow(E1,
                3) * SCV - E3) / (6 * FastMath.pow(E1, 3) * SCV + 3 * FastMath.pow(E1, 3) * SCV * SCV + 3 * FastMath.pow(E1,
                3) - 2 * E3);
            q10 = -3 * (SCV - 1) * FastMath.pow(E1, 2) / (6 * FastMath.pow(E1, 3) * SCV - E3);
        } else {
            mu00 = mmpp2_fit_mu00(E1, E2, E3, G2);
            mu11 = mmpp2_fit_mu11(E1, E2, E3, G2);
            q01 = mmpp2_fit_q01(E1, E2, E3, G2);
            q10 = mmpp2_fit_q10(E1, E2, E3, G2);
        }

        Matrix d0 = new Matrix(new double[][]{{-mu00 - q01, q01}, {q10, -mu11 - q10}});
        Matrix d1 = new Matrix(new double[][]{{mu00, 0.0}, {0.0, mu11}});
        MatrixCell map = new MatrixCell();
        map.set(0, d0);
        map.set(1, d1);
        return map;
    }
}
