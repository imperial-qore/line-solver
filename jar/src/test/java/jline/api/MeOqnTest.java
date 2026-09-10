package jline.api;

import jline.api.nc.MeOqnResult;
import jline.api.nc.Me_oqn;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Validates the Maximum Entropy open-network algorithm (Kouvatsos 1994)
 * against closed-form results: M/M/1, M/M/1 with Bernoulli feedback
 * (Jackson-exact), GE/M/1 (paper closed form), tandem M/M/1, M/M/c
 * (Erlang-C exact) and the GE/GE/inf building block.
 */
public class MeOqnTest {

    private static final double TOL = 1e-8;

    private static Matrix scalar(double v) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, v);
        return m;
    }

    private static Matrix[][] noRouting(int M, int R) {
        Matrix[][] P = new Matrix[M][M];
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                P[j][i] = new Matrix(R, 1);
            }
        }
        return P;
    }

    @Test
    public void testMM1() {
        MeOqnResult res = Me_oqn.me_oqn(1, 1, scalar(0.6), scalar(1.0),
                scalar(1.0), scalar(1.0), noRouting(1, 1));
        assertEquals(1.5, res.getL().get(0, 0), TOL);
        assertEquals(1.0, res.getCd().get(0, 0), TOL);
        assertEquals(2.5, res.getW().get(0, 0), TOL);
    }

    @Test
    public void testMM1Feedback() {
        // M/M/1 with 50% Bernoulli feedback: composite service is exponential,
        // hence Jackson-exact L = rho/(1-rho) = 1 at rho = 0.5
        Matrix[][] P = noRouting(1, 1);
        P[0][0].set(0, 0, 0.5);
        MeOqnResult res = Me_oqn.me_oqn(1, 1, scalar(1.0), scalar(1.0),
                scalar(4.0), scalar(1.0), P);
        assertEquals(1.0, res.getL().get(0, 0), TOL);
        assertEquals(2.0, res.getLambda().get(0, 0), TOL); // visit-inclusive
        assertEquals(0.5, res.getRho().get(0, 0), TOL);
    }

    @Test
    public void testGEM1() {
        // GE/M/1 with Ca=5, rho=0.5: paper closed form
        // L = rho/2*(Ca+1) + rho^2*(Ca+Cs)/(2*(1-rho)) = 3.0
        MeOqnResult res = Me_oqn.me_oqn(1, 1, scalar(1.0), scalar(5.0),
                scalar(2.0), scalar(1.0), noRouting(1, 1));
        assertEquals(3.0, res.getL().get(0, 0), TOL);
    }

    @Test
    public void testTandemMM1() {
        Matrix lambda0 = new Matrix(2, 1);
        lambda0.set(0, 0, 0.5);
        Matrix Ca0 = new Matrix(2, 1);
        Ca0.set(0, 0, 1.0);
        Matrix mu = new Matrix(2, 1);
        mu.set(0, 0, 1.0);
        mu.set(1, 0, 0.8);
        Matrix Cs = Matrix.ones(2, 1);
        Matrix[][] P = noRouting(2, 1);
        P[0][1].set(0, 0, 1.0);
        MeOqnResult res = Me_oqn.me_oqn(2, 1, lambda0, Ca0, mu, Cs, P);
        assertEquals(1.0, res.getL().get(0, 0), TOL);
        assertEquals(5.0 / 3.0, res.getL().get(1, 0), TOL);
    }

    @Test
    public void testMMcErlangC() {
        // M/M/3 with lambda=2, mu=1: Erlang-C exact L = 26/9
        MeOqnResult res = Me_oqn.me_oqn(1, 1, scalar(2.0), scalar(1.0),
                scalar(1.0), scalar(1.0), noRouting(1, 1), scalar(3.0),
                new jline.api.nc.MeOqnOptions());
        assertEquals(26.0 / 9.0, res.getL().get(0, 0), TOL);
        assertEquals(1.0, res.getCd().get(0, 0), TOL);
        assertEquals(2.0 / 3.0, res.getRho().get(0, 0), TOL);
    }

    @Test
    public void testInfiniteServer() {
        // GE/GE/inf: L = lambda/mu, departures inherit the arrival scv
        MeOqnResult res = Me_oqn.me_oqn(1, 1, scalar(3.0), scalar(1.0),
                scalar(2.0), scalar(1.0), noRouting(1, 1),
                scalar(Double.POSITIVE_INFINITY), new jline.api.nc.MeOqnOptions());
        assertEquals(1.5, res.getL().get(0, 0), TOL);
        assertEquals(1.0, res.getCd().get(0, 0), TOL);
    }

    @Test
    public void testTwoClassTandemParity() {
        // Two-class tandem; reference values from the MATLAB implementation
        Matrix lambda0 = new Matrix(2, 2);
        lambda0.set(0, 0, 0.3);
        lambda0.set(0, 1, 0.2);
        Matrix Ca0 = Matrix.ones(2, 2);
        Matrix mu = new Matrix(2, 2);
        mu.set(0, 0, 1.0);
        mu.set(0, 1, 1.0);
        mu.set(1, 0, 0.8);
        mu.set(1, 1, 0.8);
        Matrix Cs = Matrix.ones(2, 2);
        Matrix[][] P = noRouting(2, 2);
        P[0][1].set(0, 0, 1.0);
        P[0][1].set(1, 0, 1.0);
        MeOqnResult res = Me_oqn.me_oqn(2, 2, lambda0, Ca0, mu, Cs, P);
        assertEquals(0.6, res.getL().get(0, 0), TOL);
        assertEquals(0.4, res.getL().get(0, 1), TOL);
        assertEquals(1.24, res.getL().get(1, 0), TOL);
        assertEquals(0.8266666666666667, res.getL().get(1, 1), TOL);
    }
}
