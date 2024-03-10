package jline.api;

import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static jline.api.PFQN.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test class for the PFQN API. These tests were determined by using release 2.0.27 of the LINE
 * software as a baseline. These tests use as inputs and outputs the values observed through the MATLAB
 * workspace whenever the corresponding functions were called/executed.
 */
class PFQNTest {
    private final double tolerance = 1.0e-15;

    @Test
    void pfqn_caTest1() {
        Matrix L = new Matrix(0, 0);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(2, 3);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.3);
        Z.set(0, 1, 0.4);
        Z.set(0, 2, 0.5);
        Z.set(1, 0, 0.2);
        Z.set(1, 1, 0.1);
        Z.set(1, 2, 0.0);
        pfqnNcReturn ret = PFQN.pfqn_ca(L, N, Z);
        assertEquals(-6.643789733147672, ret.lG, tolerance);
        assertEquals(0.001302083333333, ret.G, tolerance);
    }

    @Test
    void pfqn_caTest2() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(2, 3);
        L.set(0, 0, 0.3);
        L.set(0, 1, 0.22);
        L.set(0, 2, 0.155);
        N.set(0, 0, 1);
        N.set(0, 1, -2);
        N.set(0, 2, -6);
        Z.set(0, 0, 0.3);
        Z.set(0, 1, 0.4);
        Z.set(0, 2, 0.5);
        Z.set(1, 0, 0.2);
        Z.set(1, 1, 0.1);
        Z.set(1, 2, 0.0);
        pfqnNcReturn ret = PFQN.pfqn_ca(L, N, Z);
        assertEquals(Double.NEGATIVE_INFINITY, ret.lG, tolerance);
        assertEquals(0.0, ret.G, tolerance);
    }

    @Test
    void pfqn_caTest3() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(2, 3);
        L.set(0, 0, 0.3);
        L.set(0, 1, 0.22);
        L.set(0, 2, 0.155);
        N.set(0, 0, 0);
        N.set(0, 1, 0);
        N.set(0, 2, 0);
        Z.set(0, 0, 0.3);
        Z.set(0, 1, 0.4);
        Z.set(0, 2, 0.5);
        Z.set(1, 0, 0.2);
        Z.set(1, 1, 0.1);
        Z.set(1, 2, 0.0);
        pfqnNcReturn ret = PFQN.pfqn_ca(L, N, Z);
        assertEquals(0.0, ret.lG, tolerance);
        assertEquals(1.0, ret.G, tolerance);
    }

    @Test
    void pfqn_caTest4() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(2, 3);
        L.set(0, 0, 0.3);
        L.set(0, 1, 0.22);
        L.set(0, 2, 0.155);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.3);
        Z.set(0, 1, 0.4);
        Z.set(0, 2, 0.5);
        Z.set(1, 0, 0.2);
        Z.set(1, 1, 0.1);
        Z.set(1, 2, 0.0);
        pfqnNcReturn ret = PFQN.pfqn_ca(L, N, Z);
        assertEquals(-4.037788879321289, ret.lG, tolerance);
        assertEquals(0.017636425600000, ret.G, tolerance);
    }

    @Test
    void pfqn_mmint2_gausslegendreTest1() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Integer m = null;
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 5);
        N.set(0, 0, 87);
        N.set(0, 1, 42);
        N.set(0, 2, 100);
        Z.set(0, 0, 3.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 1.2);
        pfqnNcReturn ret = PFQN.pfqn_mmint2_gausslegendre(L, N, Z, m);
        assertEquals(4.000255068564353e+191, ret.G, 4.000255068564353e+191*tolerance);
        assertEquals(4.411801108880907e+02, ret.lG, 4.411801108880907e+02*tolerance);
    }

    @Test
    void pfqn_mmint2_gausslegendreTest2() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Integer m = 3;
        L.set(0, 0, 1);
        L.set(0, 1, 2);
        L.set(0, 2, 3);
        N.set(0, 0, 87);
        N.set(0, 1, 42);
        N.set(0, 2, 100);
        Z.set(0, 0, 3.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 1.2);
        pfqnNcReturn ret = PFQN.pfqn_mmint2_gausslegendre(L, N, Z, m);
        assertEquals(Double.POSITIVE_INFINITY, ret.G, tolerance);
        assertEquals(5.155963658778453e+03, ret.lG, 5.155963658778453e+03*tolerance);
    }

    @Test
    void pfqn_mmint2_gausslegendreTest3() {
        Matrix L = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Integer m = 5;
        L.set(0, 0, 3.24);
        N.set(0, 0, 32);
        Z.set(0, 0, 0.4215);
        pfqnNcReturn ret = PFQN.pfqn_mmint2_gausslegendre(L, N, Z, m);
        assertEquals(Double.POSITIVE_INFINITY, ret.G, tolerance);
        assertEquals(5.234927321273997e+03, ret.lG, 5.234927321273997e+03*tolerance);
    }

    @Test
    void pfqn_mmint2_gausslegendreTest4() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Integer m = 1;
        L.set(0, 0, 0.211000000000000);
        L.set(0, 1, 0.021340000000000);
        L.set(0, 2, 2.45199840000000e-05);
        N.set(0, 0, 21);
        N.set(0, 1, 8);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.00948300000000000);
        Z.set(0, 1, 0.249100000000000);
        Z.set(0, 2, 0.00567810000000000);
        pfqnNcReturn ret = PFQN.pfqn_mmint2_gausslegendre(L, N, Z, m);
        assertEquals(8.641953895633185e-28, ret.G, 8.641953895633185e-28*tolerance);
        assertEquals( -62.315753901256400, ret.lG, tolerance);
    }

    @Test
    void pfqn_leTest1() {
        Matrix L = new Matrix(0, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        Z.set(0, 0, 0.144);
        Z.set(0, 1, 3.214);
        Z.set(0, 2, 55.221);
        pfqnNcReturn ret = PFQN.pfqn_le(L, N, Z);
        assertEquals(7.265624088932120e-08, ret.G, 7.265624088932120e-08*tolerance);
        assertEquals(  -16.437526547118978, ret.lG, tolerance);
    }

    @Test
    void pfqn_leTest2() {
        Matrix L = new Matrix(3, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = null;
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 2);
        L.set(1, 0, 4);
        L.set(1, 1, 2);
        L.set(1, 2, 1);
        L.set(2, 0, 0.24);
        L.set(2, 1, 2);
        L.set(2, 2, 1.44);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        pfqnNcReturn ret = PFQN.pfqn_le(L, N, Z);
        assertEquals(7.005301586504859e+07, ret.G, 7.005301586504859e+07*tolerance);
        assertEquals(  18.064762882854776, ret.lG, tolerance);
    }

    @Test
    void pfqn_leTest3() {
        Matrix L = new Matrix(3, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 2);
        L.set(1, 0, 4);
        L.set(1, 1, 2);
        L.set(1, 2, 1);
        L.set(2, 0, 0.24);
        L.set(2, 1, 2);
        L.set(2, 2, 1.44);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        Z.set(0, 0, 5.533);
        Z.set(0, 1, 2.1567);
        Z.set(0, 2, 9.20195);
        pfqnNcReturn ret = PFQN.pfqn_le(L, N, Z);
        assertEquals(3.598124832428049e+08, ret.G, 3.598124832428049e+08*tolerance);
        assertEquals(    19.701093573828281, ret.lG, tolerance);
    }

    @Test
    void pfqn_leTest4() {
        Matrix L = new Matrix(3, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        L.set(0, 0, 1);
        L.set(1, 0, 4);
        L.set(2, 0, 0.24);
        N.set(0, 0, 7);
        Z.set(0, 0, 5.533);
        pfqnNcReturn ret = PFQN.pfqn_le(L, N, Z);
        assertEquals(6.532800848885124e+04, ret.G, 6.532800848885124e+04*tolerance);
        assertEquals(    11.087176143501352, ret.lG, tolerance);
    }

    @Test
    void pfqn_lsTest1() {
        Matrix L = new Matrix(0, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        Z.set(0, 0, 0.144);
        Z.set(0, 1, 3.214);
        Z.set(0, 2, 55.221);
        int I = 10000000;
        pfqnNcReturn ret = PFQN.pfqn_ls(L, N, Z, I);
        assertEquals(7.265624088932120e-08, ret.G, 7.265624088932120e-08*tolerance);
        assertEquals(  -16.437526547118978, ret.lG, tolerance);
    }

    @Test
    void pfqn_lsTest2() {
        Matrix L = new Matrix(3, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = null;
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 2);
        L.set(1, 0, 4);
        L.set(1, 1, 2);
        L.set(1, 2, 1);
        L.set(2, 0, 0.24);
        L.set(2, 1, 2);
        L.set(2, 2, 1.44);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        int I = 10000;
        pfqnNcReturn ret = PFQN.pfqn_ls(L, N, Z, I);
        assertEquals(8.433706150183025e+07, ret.G, 8.433706150183025e+07*0.1);
        assertEquals(  18.250331964578145, ret.lG, 18.250331964578145*0.01);
    }

    @Test
    void pfqn_lsTest3() {
        Matrix L = new Matrix(4, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = null;
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 2);
        L.set(1, 0, 4);
        L.set(1, 1, 2);
        L.set(1, 2, 1);
        L.set(2, 0, 0.24);
        L.set(2, 1, 2);
        L.set(2, 2, 1.44);
        L.set(3, 0, 0.000001);
        L.set(3, 1, 0.000001);
        L.set(3, 2, 0.000001);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        int I = 10000;
        pfqnNcReturn ret = PFQN.pfqn_ls(L, N, Z, I);
        assertEquals(8.309387711087850e+07, ret.G, 8.433706150183025e+07*0.1);
        assertEquals(  18.235481576134909, ret.lG, 18.250331964578145*0.01);
    }

    @Test
    void pfqn_lsTest4() {
        Matrix L = new Matrix(3, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 1);
        L.set(0, 1, 3);
        L.set(0, 2, 2);
        L.set(1, 0, 4);
        L.set(1, 1, 2);
        L.set(1, 2, 1);
        L.set(2, 0, 0.24);
        L.set(2, 1, 2);
        L.set(2, 2, 1.44);
        N.set(0, 0, 7);
        N.set(0, 1, 2);
        N.set(0, 2, 1);
        Z.set(0, 0, 0.144);
        Z.set(0, 1, 3.214);
        Z.set(0, 2, 55.221);
        int I = 10000;
        pfqnNcReturn ret = PFQN.pfqn_ls(L, N, Z, I);
        assertEquals(5.368090955068173e+08, ret.G, 5.368090955068173e+08*0.1);
        assertEquals(20.101153087417874, ret.lG, 20.101153087417874*0.01);
    }

    @Test
    void pfqn_nc_sanitizeTest1() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix Z = new Matrix(3, 4);
        Matrix lambda = new Matrix(1, 4);
        L.set(0, 0, 0.1);
        L.set(0, 1, 0.2);
        L.set(0, 2, 0.3);
        L.set(0, 3, 0.4);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        N.set(0, 3, 4);
        Z.set(0, 0, 0.0001);
        Z.set(0, 1, 0.153);
        Z.set(0, 2, 0.129);
        Z.set(0, 3, 0.8482);
        Z.set(1, 0, 0.0004);
        Z.set(1, 1, 0.122);
        Z.set(1, 2, 1.23);
        Z.set(1, 3, 3.219);
        Z.set(2, 0, 0.00001);
        Z.set(2, 1, 3.21);
        Z.set(2, 2, 4);
        Z.set(2, 3, 5);
        lambda.set(0, 0, 0.553);
        lambda.set(0, 1, 0.142);
        lambda.set(0, 2, 0.251);
        lambda.set(0, 3, 0.984);
        double atol = 0.15;
        pfqnNcSanitizeReturn ret = PFQN.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix lambda_new = ret.lambda;
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lGremaind = ret.lGremaind;
        assertEquals(0.142, lambda_new.get(0), 0.142*tolerance);
        assertEquals(0.251, lambda_new.get(1), 0.251*tolerance);
        assertEquals(0.984, lambda_new.get(2), 0.984*tolerance);
        assertEquals(1, L_new.get(0), tolerance);
        assertEquals(1, L_new.get(1), tolerance);
        assertEquals(1, L_new.get(2), tolerance);
        assertEquals(2, N_new.get(0), tolerance);
        assertEquals(3, N_new.get(1), tolerance);
        assertEquals(4, N_new.get(2), tolerance);
        assertEquals(0.765, Z_new.get(0, 0), tolerance);
        assertEquals(0.43, Z_new.get(0, 1), tolerance);
        assertEquals(2.1205, Z_new.get(0, 2), tolerance);
        assertEquals(0.61, Z_new.get(1, 0), tolerance);
        assertEquals(4.100000000000001, Z_new.get(1, 1), tolerance);
        assertEquals(8.047499999999999, Z_new.get(1, 2), tolerance);
        assertEquals(16.049999999999997, Z_new.get(2, 0), tolerance);
        assertEquals(13.333333333333334, Z_new.get(2, 1), tolerance);
        assertEquals(12.5, Z_new.get(2, 2), tolerance);


        assertEquals(-10.495957165342629, lGremaind, tolerance);
    }

    @Test
    void pfqn_nc_sanitizeTest2() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix Z = new Matrix(1, 4);
        Matrix lambda = new Matrix(1, 4);
        L.set(0, 0, 0.1);
        L.set(0, 1, 0.2);
        L.set(0, 2, 0.3);
        L.set(0, 3, 0.4);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        N.set(0, 3, 4);
        Z.set(0, 0, 0.00001);
        Z.set(0, 1, 3.21);
        Z.set(0, 2, 4);
        Z.set(0, 3, 5);
        lambda.set(0, 0, 0.553);
        lambda.set(0, 1, 0.142);
        lambda.set(0, 2, 0.251);
        lambda.set(0, 3, 0.984);
        double atol = 0.001;
        pfqnNcSanitizeReturn ret = PFQN.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix lambda_new = ret.lambda;
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lGremaind = ret.lGremaind;
        assertEquals(0.553, lambda_new.get(0), 0.553*tolerance);
        assertEquals(0.142, lambda_new.get(1), 0.142*tolerance);
        assertEquals(0.251, lambda_new.get(2), 0.251*tolerance);
        assertEquals(0.984, lambda_new.get(3), 0.984*tolerance);
        assertEquals(1, L_new.get(0), tolerance);
        assertEquals(1, L_new.get(1), tolerance);
        assertEquals(1, L_new.get(2), tolerance);
        assertEquals(1, L_new.get(3), tolerance);
        assertEquals(1, N_new.get(0), tolerance);
        assertEquals(4, N_new.get(1), tolerance);
        assertEquals(3, N_new.get(2), tolerance);
        assertEquals(2, N_new.get(3), tolerance);
        assertEquals(0.0001, Z_new.get(0, 0), tolerance);
        assertEquals(12.5, Z_new.get(0, 1), tolerance);
        assertEquals(13.333333333333334, Z_new.get(0, 2), tolerance);
        assertEquals(16.049999999999997, Z_new.get(0, 3), tolerance);

        assertEquals(-12.798542258336674, lGremaind, tolerance);
    }

    @Test
    void pfqn_nc_sanitizeTest3() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix Z = new Matrix(1, 4);
        Matrix lambda = new Matrix(1, 4);
        L.set(0, 0, 0.1);
        L.set(0, 1, 0.2);
        L.set(0, 2, 0.3);
        L.set(0, 3, 0.4);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        N.set(0, 3, 4);
        Z.set(0, 0, 6.32);
        Z.set(0, 1, 1.25);
        Z.set(0, 2, 5.322);
        Z.set(0, 3, 0.221);
        lambda.set(0, 0, 0.553);
        lambda.set(0, 1, 0.142);
        lambda.set(0, 2, 0.251);
        lambda.set(0, 3, 0.984);
        double atol = 0.001;
        pfqnNcSanitizeReturn ret = PFQN.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix lambda_new = ret.lambda;
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lGremaind = ret.lGremaind;
        assertEquals(0.553, lambda_new.get(0), 0.553*tolerance);
        assertEquals(0.142, lambda_new.get(1), 0.142*tolerance);
        assertEquals(0.251, lambda_new.get(2), 0.251*tolerance);
        assertEquals(0.984, lambda_new.get(3), 0.984*tolerance);
        assertEquals(1, L_new.get(0), tolerance);
        assertEquals(1, L_new.get(1), tolerance);
        assertEquals(1, L_new.get(2), tolerance);
        assertEquals(1, L_new.get(3), tolerance);
        assertEquals(4, N_new.get(0), tolerance);
        assertEquals(2, N_new.get(1), tolerance);
        assertEquals(3, N_new.get(2), tolerance);
        assertEquals(1, N_new.get(3), tolerance);
        assertEquals(0.5525, Z_new.get(0, 0), tolerance);
        assertEquals(6.25, Z_new.get(0, 1), tolerance);
        assertEquals(17.740000000000002, Z_new.get(0, 2), tolerance);
        assertEquals(63.200000000000003, Z_new.get(0, 3), tolerance);

        assertEquals(-12.798542258336674, lGremaind, tolerance);
    }

    @Test
    void pfqn_nc_sanitizeTest4() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix Z = new Matrix(1, 4);
        Matrix lambda = new Matrix(1, 4);
        L.set(0, 0, 0.1);
        L.set(0, 1, 0.2);
        L.set(0, 2, 0.3);
        L.set(0, 3, 0.4);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        N.set(0, 3, 4);
        Z.set(0, 0, 6.32);
        Z.set(0, 1, 1.25);
        Z.set(0, 2, 5.322);
        Z.set(0, 3, 0.221);
        lambda.set(0, 0, 0.553);
        lambda.set(0, 1, 0.142);
        lambda.set(0, 2, 0.251);
        lambda.set(0, 3, 0.984);
        double atol = 1;
        pfqnNcSanitizeReturn ret = PFQN.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix lambda_new = ret.lambda;
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lGremaind = ret.lGremaind;
        assertEquals(0.553, lambda_new.get(0), 0.553*tolerance);
        assertEquals(0.142, lambda_new.get(1), 0.142*tolerance);
        assertEquals(0.251, lambda_new.get(2), 0.251*tolerance);
        assertTrue(L_new.isEmpty());
        assertTrue(N_new.isEmpty());
        assertTrue(Z_new.isEmpty());
        assertEquals(5.513794359225622, lGremaind, 1e5*tolerance);
    }

    @Test
    void pfqn_nc_sanitizeTest5() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix Z = new Matrix(1, 4);
        Matrix lambda = new Matrix(1, 4);
        L.set(0, 0, 0.1);
        L.set(0, 1, 0.2);
        L.set(0, 2, 0.3);
        L.set(0, 3, 0.4);
        N.set(0, 0, 1);
        N.set(0, 1, 2);
        N.set(0, 2, 3);
        N.set(0, 3, 4);
        Z.set(0, 0, 6.32);
        Z.set(0, 1, 1.25);
        Z.set(0, 2, 5.322);
        Z.set(0, 3, 0.221);
        lambda.set(0, 0, 0.553);
        lambda.set(0, 1, 0.142);
        lambda.set(0, 2, 0.251);
        lambda.set(0, 3, 0.984);
        double atol = 0.3;
        pfqnNcSanitizeReturn ret = PFQN.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix lambda_new = ret.lambda;
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lGremaind = ret.lGremaind;
        assertEquals(0.553, lambda_new.get(0), 0.553*tolerance);
        assertEquals(0.142, lambda_new.get(1), 0.142*tolerance);
        assertEquals(0.251, lambda_new.get(2), 0.251*tolerance);
        assertEquals(0.984, lambda_new.get(3), 0.984*tolerance);
        assertEquals(1, L_new.get(0), tolerance);
        assertEquals(1, L_new.get(1), tolerance);
        assertEquals(4, N_new.get(0), tolerance);
        assertEquals(3, N_new.get(1), tolerance);
        assertEquals(0.5525, Z_new.get(0, 0), tolerance);
        assertEquals(17.740000000000002, Z_new.get(0, 1), tolerance);

        assertEquals(-5.680222210247187, lGremaind, tolerance);
    }

    @Test
    void pfqn_comomrmTest1() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 2.1);
        L.set(0, 1, 5.66);
        L.set(0, 2, 9.224);
        N.set(0, 0, 9);
        N.set(0, 1, 4);
        N.set(0, 2, 2);
        Z.set(0, 0, 1.224);
        Z.set(0, 1, 0.251);
        Z.set(0, 2, 4.2211);
        int m = 1;
        double atol = 0.2;
        pfqnComomrmReturn ret = pfqn_comomrm(L, N, Z, m, atol);
        double lG = ret.lG;
        Matrix lGbasis = ret.lGbasis;
        assertEquals(6, lGbasis.getNumRows());
        assertEquals(1, lGbasis.getNumCols());
        assertEquals(14.392932838406633, lGbasis.get(0), 1e5*tolerance);
        assertEquals(13.030110328778971, lGbasis.get(1), 1e5*tolerance);
        assertEquals(12.309110365151165, lGbasis.get(2), 1e5*tolerance);
        assertEquals(11.646858331869151, lGbasis.get(3), 1e5*tolerance);
        assertEquals(10.352237763479343, lGbasis.get(4), 1e5*tolerance);
        assertEquals(9.629154581674436, lGbasis.get(5), 1e5*tolerance);
        assertEquals(29.701607569227431, lG, 1e5*tolerance);
    }

    @Test
    void pfqn_comomrmTest2() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 2.1);
        L.set(0, 1, 5.66);
        L.set(0, 2, 9.224);
        N.set(0, 0, 9);
        N.set(0, 1, 4);
        N.set(0, 2, 2);
        Z.set(0, 0, 1.224);
        Z.set(0, 1, 0.251);
        Z.set(0, 2, 4.2211);
        int m = 2;
        double atol = 0.5;
        pfqnComomrmReturn ret = pfqn_comomrm(L, N, Z, m, atol);
        double lG = ret.lG;
        Matrix lGbasis = ret.lGbasis;
        assertEquals(6, lGbasis.getNumRows());
        assertEquals(1, lGbasis.getNumCols());
        assertEquals(16.509667835109035, lGbasis.get(0), 1e5*tolerance);
        assertEquals(15.083238517521064, lGbasis.get(1), 1e5*tolerance);
        assertEquals(14.364043210388438, lGbasis.get(2), 1e5*tolerance);
        assertEquals(14.392932838406631, lGbasis.get(3), 1e5*tolerance);
        assertEquals(13.030110328778971, lGbasis.get(4), 1e5*tolerance);
        assertEquals(12.309110365151165, lGbasis.get(5), 1e5*tolerance);
        assertEquals(32.447682075764909, lG, 1e5*tolerance);
    }

    @Test
    void pfqn_comomrmTest3() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 2.1);
        L.set(0, 1, 5.66);
        L.set(0, 2, 9.224);
        N.set(0, 0, 9);
        N.set(0, 1, 4);
        N.set(0, 2, 2);
        Z.set(0, 0, 1.224);
        Z.set(0, 1, 0.251);
        Z.set(0, 2, 1e-9);
        int m = 2;
        double atol = 0.5;
        pfqnComomrmReturn ret = pfqn_comomrm(L, N, Z, m, atol);
        double lG = ret.lG;
        Matrix lGbasis = ret.lGbasis;
        assertEquals(6, lGbasis.getNumRows());
        assertEquals(1, lGbasis.getNumCols());
        assertEquals(17.144180857969538, lGbasis.get(0), 1e5*tolerance);
        assertEquals(15.023945666602936, lGbasis.get(1), 1e5*tolerance);
        assertEquals(15.714262659486160, lGbasis.get(2), 1e5*tolerance);
        assertEquals(15.023945666609865, lGbasis.get(3), 1e5*tolerance);
        assertEquals(12.967020460050561, lGbasis.get(4), 1e5*tolerance);
        assertEquals(13.657139249487017, lGbasis.get(5), 1e5*tolerance);
        assertEquals(33.078694903968140, lG, 1e5*tolerance);
    }

    @Test
    void pfqn_comomrmTest4() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        L.set(0, 0, 2.1);
        L.set(0, 1, 5.66);
        L.set(0, 2, 9.224);
        N.set(0, 0, 9);
        N.set(0, 1, 4);
        N.set(0, 2, 2);
        Z.set(0, 0, 1e-9);
        Z.set(0, 1, 1e-9);
        Z.set(0, 2, 1e-9);
        int m = 2;
        double atol = 0.5;
        pfqnComomrmReturn ret = pfqn_comomrm(L, N, Z, m, atol);
        double lG = ret.lG;
        Matrix lGbasis = ret.lGbasis;
        assertEquals(8, lGbasis.getNumRows());
        assertEquals(1, lGbasis.getNumCols());
        assertEquals(16.832044959147531, lGbasis.get(0), 1e5*tolerance);
        assertEquals(14.691978795651258, lGbasis.get(1), 1e5*tolerance);
        assertEquals(15.385125976211203, lGbasis.get(2), 1e5*tolerance);
        assertEquals(16.196056192427530, lGbasis.get(3), 1e5*tolerance);
        assertEquals(13.998831615091312, lGbasis.get(4), 1e5*tolerance);
        assertEquals(11.919390073411476, lGbasis.get(5), 1e5*tolerance);
        assertEquals(12.612537253971421, lGbasis.get(6), 1e5*tolerance);
        assertEquals(13.423467470187749, lGbasis.get(7), 1e5*tolerance);
        assertEquals(32.053580852449592, lG, 1e5*tolerance);
    }

    @Test
    void compute_norm_constTest1() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "ca";
        L.set(0, 0, 0.2);
        L.set(0, 1, 0.338);
        L.set(0, 2, 1.1948);
        L.set(1, 0, 7.291);
        L.set(1, 1, 2);
        L.set(1, 2, 3.9908);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.221);
        Z.set(0, 1, 0.333);
        Z.set(0, 2, 0.1986);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(26.052104159913611, ret.lG, tolerance);
        assertEquals("ca", ret.method);
    }

    @Test
    void compute_norm_constTest2() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "adaptive";
        L.set(0, 0, 0.2);
        L.set(0, 1, 0.338);
        L.set(0, 2, 1.1948);
        L.set(1, 0, 7.291);
        L.set(1, 1, 2);
        L.set(1, 2, 3.9908);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.221);
        Z.set(0, 1, 0.333);
        Z.set(0, 2, 0.1986);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(26.052104159913611, ret.lG, tolerance);
        assertEquals("ca", ret.method);
    }

    @Test
    void compute_norm_constTest3() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "adaptive";
        L.set(0, 0, 0.2);
        L.set(0, 1, 0.338);
        L.set(0, 2, 1.1948);
        L.set(1, 0, 7.291);
        L.set(1, 1, 2);
        L.set(1, 2, 3.9908);
        N.set(0, 0, 40);
        N.set(0, 1, 5);
        N.set(0, 2, 23);
        Z.set(0, 0, 0.221);
        Z.set(0, 1, 0.333);
        Z.set(0, 2, 0.1986);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(1.700939228986259e+02, ret.lG, 1.700939228986259e+02*tolerance);
        assertEquals("le", ret.method);
    }

    @Test
    void compute_norm_constTest4() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 0.2);
        L.set(0, 1, 0.338);
        L.set(0, 2, 1.1948);
        N.set(0, 0, 40);
        N.set(0, 1, 5);
        N.set(0, 2, 23);
        Z.set(0, 0, 0.221);
        Z.set(0, 1, 0.333);
        Z.set(0, 2, 0.1986);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(-9.688694712607486, ret.lG, tolerance);
        assertEquals("gleint", ret.method);
    }

    @Test
    void compute_norm_constTest5() {
        Matrix L = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 0.221);
        N.set(0, 0, 76);
        Z.set(0, 0, 0.0);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(1.147290358872932e+02, ret.lG, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void compute_norm_constTest6() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "sampling";
        L.set(0, 0, 0.221);
        L.set(0, 1, 1.9847);
        L.set(0, 2, 3.8799);
        L.set(1, 0, 0.2);
        L.set(1, 1, 0.338);
        L.set(1, 2, 1.1948);
        N.set(0, 0, 3);
        N.set(0, 1, 8);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.2441);
        Z.set(0, 1, 2.14325);
        Z.set(0, 2, 0.00921);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(14.456204385514516, ret.lG, 0.01*14.456204385514516);
        assertEquals("ls", ret.method);
    }

    @Test
    void compute_norm_constTest7() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "exact";
        L.set(0, 0, 0.221);
        L.set(0, 1, 1.9847);
        L.set(0, 2, 3.8799);
        L.set(1, 0, 0.2);
        L.set(1, 1, 0.338);
        L.set(1, 2, 1.1948);
        N.set(0, 0, 3);
        N.set(0, 1, 8);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.2441);
        Z.set(0, 1, 2.14325);
        Z.set(0, 2, 0.00921);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(14.472899169591592, ret.lG, tolerance);
        assertEquals("ca", ret.method);
    }

    @Test
    void compute_norm_constTest8() {
        Matrix L = new Matrix(1, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        SolverOptions options = new SolverOptions();
        options.method = "comom";
        L.set(0, 0, 0.221);
        L.set(0, 1, 1.9847);
        L.set(0, 2, 3.8799);
        N.set(0, 0, 3);
        N.set(0, 1, 8);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.2441);
        Z.set(0, 1, 2.14325);
        Z.set(0, 2, 0.00921);
        pfqnNcXQReturn ret = compute_norm_const(L,N,Z,options);
        assertEquals(14.042355293715579, ret.lG, tolerance);
        assertEquals("comom", ret.method);
    }

    @Test
    void pfqn_ncTest1() {
        SolverOptions options = new SolverOptions();
        options.method = "";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 2);
        lambda.set(0, 1, 3);
        lambda.set(0, 2, 1);
        L.set(0, 0, 9);
        L.set(0, 1, 2);
        L.set(0, 2, 1);
        L.set(1, 0, 1);
        L.set(1, 1, 2);
        L.set(1, 2, 3);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 0.1);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(-29.122675992706263, ret.lG, tolerance);
    }

    @Test
    void pfqn_ncTest2() {
        SolverOptions options = new SolverOptions();
        options.method = "le";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 0);
        lambda.set(0, 1, 0);
        lambda.set(0, 2, 0);
        L.set(0, 0, 9);
        L.set(0, 1, 2);
        L.set(0, 2, 1);
        L.set(1, 0, 1);
        L.set(1, 1, 2);
        L.set(1, 2, 3);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 0.1);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(24.556348056490499, ret.lG, tolerance);
    }

    @Test
    void pfqn_ncTest3() {
        SolverOptions options = new SolverOptions();
        options.method = "";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 0);
        lambda.set(0, 1, 0);
        lambda.set(0, 2, 0.2);
        L.set(0, 0, 9);
        L.set(0, 1, 2);
        L.set(0, 2, 0);
        L.set(1, 0, 1);
        L.set(1, 1, 2);
        L.set(1, 2, 0);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 0);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(Double.NaN, ret.lG, tolerance);
    }

    @Test
    void pfqn_ncTest4() {
        SolverOptions options = new SolverOptions();
        options.method = "";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 0);
        lambda.set(0, 1, 0);
        lambda.set(0, 2, 0);
        L.set(0, 0, 9);
        L.set(0, 1, 2);
        L.set(0, 2, 4);
        L.set(1, 0, 0);
        L.set(1, 1, 0);
        L.set(1, 2, 0);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0);
        Z.set(0, 1, 0);
        Z.set(0, 2, 0);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(26.643426748808118, ret.lG, tolerance);
    }

    @Test
    void pfqn_ncTest5() {
        SolverOptions options = new SolverOptions();
        options.method = "le";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 0);
        lambda.set(0, 1, 0);
        lambda.set(0, 2, 0);
        L.set(0, 0, 0);
        L.set(0, 1, 2);
        L.set(0, 2, 4);
        L.set(1, 0, 0);
        L.set(1, 1, 9.21);
        L.set(1, 2, 18.2134);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 1.392);
        Z.set(0, 1, 2.473);
        Z.set(0, 2, 1.094);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(22.217278817021842, ret.lG, tolerance);
    }

    @Test
    void pfqn_ncTest6() {
        SolverOptions options = new SolverOptions();
        options.method = "ca";
        Matrix lambda = new Matrix(1, 3);
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        lambda.set(0, 0, 0);
        lambda.set(0, 1, 0);
        lambda.set(0, 2, 0);
        L.set(0, 0, 0);
        L.set(0, 1, 2);
        L.set(0, 2, 4);
        L.set(1, 0, 8);
        L.set(1, 1, 9.21);
        L.set(1, 2, 18.2134);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 1.392);
        Z.set(0, 1, 2.473);
        Z.set(0, 2, 1.094);
        pfqnNcXQReturn ret = pfqn_nc(lambda, L, N, Z, options);
        assertEquals(38.695882901953304, ret.lG, tolerance);
    }

    @Test
    void pfqn_gldsingleTest1() {
        Matrix L = new Matrix(3, 1);
        Matrix N = new Matrix(1, 1);
        Matrix mu = new Matrix(3, 5);
        L.set(0, 0.123);
        L.set(1, 0.456);
        L.set(2, 0.789);
        N.set(0, 5);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        pfqnNcReturn ret = pfqn_gldsingle(L,N,mu,null);
        assertEquals(1.232841907438090, ret.lG, tolerance);
        assertEquals(3.430966182727118, ret.G, tolerance);
    }

    @Test
    void pfqn_gldsingleTest2() {
        Matrix L = new Matrix(4, 1);
        Matrix N = new Matrix(1, 1);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0.123);
        L.set(1, 0.456);
        L.set(2, 0.789);
        L.set(3, 0.0012);
        N.set(0, 6);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(0, 5, 0.098);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(1, 5, 8.12985);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1.6);
        mu.set(3, 1, 1.7);
        mu.set(3, 2, 1.8);
        mu.set(3, 3, 1.9);
        mu.set(3, 4, 2.0);
        mu.set(3, 5, 4.198765);
        pfqnNcReturn ret = pfqn_gldsingle(L,N,mu,null);
        assertEquals(0.914829050814071, ret.lG, tolerance);
        assertEquals(2.496348466688459, ret.G, tolerance);
    }

    @Test
    void pfqn_mushiftTest() {
        Matrix mu = new Matrix(3, 5);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        Matrix ret = pfqn_mushift(mu, 2);
        assertEquals(3, ret.getNumRows());
        assertEquals(4, ret.getNumCols());
        assertEquals(0.1, ret.get(0, 0), tolerance);
        assertEquals(0.2, ret.get(0, 1), tolerance);
        assertEquals(0.3, ret.get(0, 2), tolerance);
        assertEquals(0.4, ret.get(0, 3), tolerance);
        assertEquals(0.6, ret.get(1, 0), tolerance);
        assertEquals(0.7, ret.get(1, 1), tolerance);
        assertEquals(0.8, ret.get(1, 2), tolerance);
        assertEquals(0.9, ret.get(1, 3), tolerance);
        assertEquals(1.2, ret.get(2, 0), tolerance);
        assertEquals(1.3, ret.get(2, 1), tolerance);
        assertEquals(1.4, ret.get(2, 2), tolerance);
        assertEquals(1.5, ret.get(2, 3), tolerance);
    }

    @Test
    void pfqn_gldTest1() {
        Matrix L = new Matrix(1, 4);
        Matrix N = new Matrix(1, 4);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0.123);
        L.set(1, 0.456);
        L.set(2, 0.789);
        L.set(3, -0.0012);
        N.set(0, 1);
        N.set(1, 0);
        N.set(2, 2);
        N.set(3, 3);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(0, 5, 0.098);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(1, 5, 8.12985);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1.6);
        mu.set(3, 1, 1.7);
        mu.set(3, 2, 1.8);
        mu.set(3, 3, 1.9);
        mu.set(3, 4, 2.0);
        mu.set(3, 5, 4.198765);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(10.573017244839603, ret.lG, tolerance);
        assertEquals(3.906636887755102e+04, ret.G, 3.906636887755102e+04*tolerance);
    }

    @Test
    void pfqn_gldTest2() {
        Matrix L = new Matrix(4, 1);
        Matrix N = new Matrix(1, 1);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0.123);
        L.set(1, 0.456);
        L.set(2, 0.789);
        L.set(3, 0.0012);
        N.set(0, 6);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(0, 5, 0.098);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(1, 5, 8.12985);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1.6);
        mu.set(3, 1, 1.7);
        mu.set(3, 2, 1.8);
        mu.set(3, 3, 1.9);
        mu.set(3, 4, 2.0);
        mu.set(3, 5, 4.198765);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(0.914829050814071, ret.lG, tolerance);
        assertEquals(2.496348466688459, ret.G, tolerance);
    }

    @Test
    void pfqn_gldTest3() {
        Matrix L = new Matrix(4, 0);
        Matrix N = new Matrix(1, 1);
        Matrix mu = new Matrix(4, 6);
        N.set(0, 6);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(0, 5, 0.098);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(1, 5, 8.12985);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1.6);
        mu.set(3, 1, 1.7);
        mu.set(3, 2, 1.8);
        mu.set(3, 3, 1.9);
        mu.set(3, 4, 2.0);
        mu.set(3, 5, 4.198765);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(Double.NEGATIVE_INFINITY, ret.lG, tolerance);
        assertEquals(0.0, ret.G, tolerance);
    }

    @Test
    void pfqn_gldTest4() {
        Matrix L = new Matrix(4, 4);
        Matrix N = new Matrix(1, 4);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0, 0.123);
        L.set(0, 1, 0.456);
        L.set(0, 2, 0.789);
        L.set(0, 3, 0.0012);
        L.set(1, 0, 9.2019);
        L.set(1, 1, 0.2958);
        L.set(1, 2, 8.9487);
        L.set(1, 3, 0.2911);
        L.set(2, 0, 0.1232);
        L.set(2, 1, 0.4344);
        L.set(2, 2, 8.789);
        L.set(2, 3, 0.3222);
        L.set(3, 0, 9.123);
        L.set(3, 1, 5.456);
        L.set(3, 2, 8.789);
        L.set(3, 3, 0.9987);
        N.set(0, 1);
        N.set(1, 0);
        N.set(2, 2);
        N.set(3, 3);
        mu.set(0, 0, 0.1);
        mu.set(0, 1, 0.2);
        mu.set(0, 2, 0.3);
        mu.set(0, 3, 0.4);
        mu.set(0, 4, 0.5);
        mu.set(0, 5, 0.098);
        mu.set(1, 0, 0.6);
        mu.set(1, 1, 0.7);
        mu.set(1, 2, 0.8);
        mu.set(1, 3, 0.9);
        mu.set(1, 4, 1.0);
        mu.set(1, 5, 8.12985);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1.6);
        mu.set(3, 1, 1.7);
        mu.set(3, 2, 1.8);
        mu.set(3, 3, 1.9);
        mu.set(3, 4, 2.0);
        mu.set(3, 5, 4.198765);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(11.221814851185433, ret.lG, tolerance);
        assertEquals(7.474329970526634e+04, ret.G, 7.474329970526634e+04*tolerance);
    }

    @Test
    void pfqn_gldTest5() {
        Matrix L = new Matrix(4, 4);
        Matrix N = new Matrix(1, 4);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0, 0.123);
        L.set(0, 1, 0.456);
        L.set(0, 2, 0.789);
        L.set(0, 3, 0.0012);
        L.set(1, 0, 9.2019);
        L.set(1, 1, 0.2958);
        L.set(1, 2, 8.9487);
        L.set(1, 3, 0.2911);
        L.set(2, 0, 0.1232);
        L.set(2, 1, 0.4344);
        L.set(2, 2, 8.789);
        L.set(2, 3, 0.3222);
        L.set(3, 0, 9.123);
        L.set(3, 1, 5.456);
        L.set(3, 2, 8.789);
        L.set(3, 3, 0.9987);
        N.set(0, 1);
        N.set(1, 0);
        N.set(2, 2);
        N.set(3, 3);
        mu.set(0, 0, 1);
        mu.set(0, 1, 1);
        mu.set(0, 2, 1);
        mu.set(0, 3, 1);
        mu.set(0, 4, 1);
        mu.set(0, 5, 1);
        mu.set(1, 0, 1);
        mu.set(1, 1, 2);
        mu.set(1, 2, 3);
        mu.set(1, 3, 4);
        mu.set(1, 4, 5);
        mu.set(1, 5, 6);
        mu.set(2, 0, 1.1);
        mu.set(2, 1, 1.2);
        mu.set(2, 2, 1.3);
        mu.set(2, 3, 1.4);
        mu.set(2, 4, 1.5);
        mu.set(2, 5, 3.19485);
        mu.set(3, 0, 1);
        mu.set(3, 1, 2);
        mu.set(3, 2, 3);
        mu.set(3, 3, 4);
        mu.set(3, 4, 5);
        mu.set(3, 5, 6);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(8.872892797319089, ret.lG, tolerance);
        assertEquals(7.135893838940448e+03, ret.G, 7.135893838940448e+03*tolerance);
    }

    @Test
    void pfqn_gldTest6() {
        Matrix L = new Matrix(4, 4);
        Matrix N = new Matrix(1, 4);
        Matrix mu = new Matrix(4, 6);
        L.set(0, 0, 0.123);
        L.set(0, 1, 0.456);
        L.set(0, 2, 0.789);
        L.set(0, 3, 0.0012);
        L.set(1, 0, 9.2019);
        L.set(1, 1, 0.2958);
        L.set(1, 2, 8.9487);
        L.set(1, 3, 0.2911);
        L.set(2, 0, 0.1232);
        L.set(2, 1, 0.4344);
        L.set(2, 2, 8.789);
        L.set(2, 3, 0.3222);
        L.set(3, 0, 9.123);
        L.set(3, 1, 5.456);
        L.set(3, 2, 8.789);
        L.set(3, 3, 0.9987);
        N.set(0, 1);
        N.set(1, 0);
        N.set(2, 2);
        N.set(3, 3);
        mu.set(0, 0, 1);
        mu.set(0, 1, 1);
        mu.set(0, 2, 1);
        mu.set(0, 3, 1);
        mu.set(0, 4, 1);
        mu.set(0, 5, 1);
        mu.set(1, 0, 1);
        mu.set(1, 1, 2);
        mu.set(1, 2, 3);
        mu.set(1, 3, 4);
        mu.set(1, 4, 5);
        mu.set(1, 5, 6);
        mu.set(2, 0, 1);
        mu.set(2, 1, 1);
        mu.set(2, 2, 1);
        mu.set(2, 3, 1);
        mu.set(2, 4, 1);
        mu.set(2, 5, 1);
        mu.set(3, 0, 1);
        mu.set(3, 1, 2);
        mu.set(3, 2, 3);
        mu.set(3, 3, 4);
        mu.set(3, 4, 5);
        mu.set(3, 5, 6);
        pfqnNcReturn ret = pfqn_gld(L,N,mu,null);
        assertEquals(9.242452277384288, ret.lG, tolerance);
        assertEquals(1.032633056063182e+04, ret.G, 1.032633056063182e+04*tolerance);
    }

    @Test
    void compute_norm_const_ldTest1() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 9.87);
        L.set(0, 1, 8.47);
        L.set(0, 2, 1.224);
        L.set(1, 0, 0.224);
        L.set(1, 1, 0.9874);
        L.set(1, 2, 0.008);
        N.set(0, 0, 1);
        N.set(0, 1, 0);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.442);
        Z.set(0, 1, 0.1948);
        Z.set(0, 2, 0.99983);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = compute_norm_const_ld(L,N,Z,mu,options);
        assertEquals(3.267460334816642, ret.lG, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void compute_norm_const_ldTest2() {
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 9.87);
        L.set(1, 0, 0.224);
        N.set(0, 0, 3);
        Z.set(0, 0, 0.442);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = compute_norm_const_ld(L,N,Z,mu,options);
        assertEquals(5.274008718004701, ret.lG, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void pfqn_ncldTest1() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 9.87);
        L.set(0, 1, 8.47);
        L.set(0, 2, 1.224);
        L.set(1, 0, 0.224);
        L.set(1, 1, 0.9874);
        L.set(1, 2, 0.008);
        N.set(0, 0, 1);
        N.set(0, 1, 0);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.442);
        Z.set(0, 1, 0.1948);
        Z.set(0, 2, 0.99983);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = pfqn_ncld(L,N,Z,mu,options);
        assertEquals(3.267460334816641, ret.lG, tolerance);
        assertEquals(26.244602131765181, ret.G, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void pfqn_ncldTest2() {
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 9.87);
        L.set(1, 0, 0.224);
        N.set(0, 0, 3);
        Z.set(0, 0, 0.442);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = pfqn_ncld(L,N,Z,mu,options);
        assertEquals(5.274008718004701, ret.lG, tolerance);
        assertEquals(1.951968854186666e+02, ret.G, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void pfqn_ncldTest3() {
        SolverOptions options = new SolverOptions();
        options.method = "default";
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Matrix mu = new Matrix(2, 12);
        L.set(0, 0, 9);
        L.set(0, 1, 2);
        L.set(0, 2, 0);
        L.set(1, 0, 1);
        L.set(1, 1, 2);
        L.set(1, 2, 0);
        N.set(0, 0, 4);
        N.set(0, 1, 5);
        N.set(0, 2, 3);
        Z.set(0, 0, 0.2);
        Z.set(0, 1, 0.3);
        Z.set(0, 2, 0.1);
        pfqnNcReturn ret = pfqn_ncld(L,N,Z,mu,options);
        assertEquals(Double.POSITIVE_INFINITY, ret.lG, tolerance);
        assertEquals(Double.POSITIVE_INFINITY, ret.G, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void pfqn_ncldTest4() {
        Matrix L = new Matrix(2, 3);
        Matrix N = new Matrix(1, 3);
        Matrix Z = new Matrix(1, 3);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 9.87);
        L.set(0, 1, 8.47);
        L.set(0, 2, 0);
        L.set(1, 0, 0.224);
        L.set(1, 1, 0.9874);
        L.set(1, 2, 0);
        N.set(0, 0, 1);
        N.set(0, 1, 0);
        N.set(0, 2, 2);
        Z.set(0, 0, 0.442);
        Z.set(0, 1, 0.1948);
        Z.set(0, 2, 0.99983);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = pfqn_ncld(L,N,Z,mu,options);
        assertEquals(1.661310754977759, ret.lG, tolerance);
        assertEquals(5.266209032245200, ret.G, tolerance);
        assertEquals("exact", ret.method);
    }

    @Test
    void pfqn_ncldTest5() {
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 3);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 0);
        L.set(1, 0, 0);
        N.set(0, 0, 3);
        Z.set(0, 0, 0.442);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(1, 0, 1);
        mu.set(1, 1, 1);
        mu.set(1, 2, 1);
        pfqnNcReturn ret = pfqn_ncld(L,N,Z,mu,options);
        assertEquals(-4.241095659941372, ret.lG, tolerance);
        assertEquals(Double.NaN, ret.G, tolerance);
        assertEquals("default", ret.method);
    }

    @Test
    void pfqn_mvaTest1() {
        Matrix L = new Matrix(1,1);
        Matrix N = new Matrix(1,1);
        Matrix Z = new Matrix(1,1);
        Matrix mi = new Matrix(1,1);
        L.set(0,0,0.416071428571429);
        mi.set(0,0,1);
        N.set(0,0,4);
        Z.set(0,0,0.484523809523810);
        pfqnMVAReturn ret = PFQN.pfqn_mva(L, N, Z, mi);
        assertEquals(-2.349956418185016, ret.lGN, tolerance);
        assertEquals(1.220823439816656, ret.CN.get(0,0), tolerance);
        assertEquals(0.975921892112806, ret.UN.get(0,0), tolerance);
        assertEquals(2.863518712189093, ret.QN.get(0,0), tolerance);
        assertEquals(2.345563345850518, ret.XN.get(0,0), tolerance);
    }

    @Test
    void pfqn_mvaTest2(){
        Matrix L =  new Matrix(2,2);
        Matrix mi = new Matrix(2,1);
        Matrix N = new Matrix(1, 2);
        Matrix Z = new Matrix(2,2);
        L.set(0,0,10);
        L.set(0,1,5);
        L.set(1,0,5);
        L.set(1,1,9);
        mi.set(0,0,1);
        mi.set(1,0,1);
        N.set(0,0,1);
        N.set(0,1,2);
        Z.set(0,0,91);
        Z.set(0,1,92);
        Z.set(1,0,93);
        Z.set(1,1,94);
        pfqnMVAReturn ret = PFQN.pfqn_mva(L, N, Z, mi);
        assertEquals(13.342678933707491, ret.lGN, tolerance);

        assertEquals(10.969277369104285, ret.CN.get(0,0), tolerance);
        assertEquals(5.743201888607152, ret.CN.get(0,1), tolerance);
        assertEquals(5.903504976200779, ret.CN.get(1,0), tolerance);
        assertEquals(10.242546122234852, ret.CN.get(1,1), tolerance);

        assertEquals(0.092701789854549, ret.UN.get(0,0), tolerance);
        assertEquals(0.091755116448852, ret.UN.get(0,1), tolerance);
        assertEquals(0.046350894927275, ret.UN.get(1,0), tolerance);
        assertEquals(0.165159209607933, ret.UN.get(1,1), tolerance);

        assertEquals(0.101687164552697, ret.QN.get(0,0), tolerance);
        assertEquals(0.105393631615683, ret.QN.get(0,1), tolerance);
        assertEquals(0.054726547770905, ret.QN.get(1,0), tolerance);
        assertEquals(0.187961202435678, ret.QN.get(1,1), tolerance);

        assertEquals(0.009270178985455, ret.XN.get(0,0), tolerance);
        assertEquals(0.018351023289770, ret.XN.get(0,1), tolerance);
    }

    @Test
    void pfqn_mvaTest3(){
        Matrix L =  new Matrix(2,1);
        Matrix mi = null;
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1,1);
        L.set(0,0,0.0417);
        L.set(1,0,0.5);
        N.set(0,0,16);
        pfqnMVAReturn ret = PFQN.pfqn_mva(L, N, Z, mi);
        assertEquals(-11.003270782052011, ret.lGN, tolerance);

        assertEquals(0.045494217761292, ret.CN.get(0,0), tolerance);
        assertEquals(7.954505782238709, ret.CN.get(1,0), tolerance);

        assertEquals(0.083400000000000, ret.UN.get(0,0), tolerance);
        assertEquals(1.000000000000000, ret.UN.get(1,0), tolerance);

        assertEquals(0.090988435522583, ret.QN.get(0,0), tolerance);
        assertEquals(15.909011564477417, ret.QN.get(1,0), tolerance);

        assertEquals(2, ret.XN.get(0,0), tolerance);
    }

    @Test
    void pfqn_mvaTest4(){
        Matrix L =  new Matrix(2,1);
        Matrix mi = null;
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1,1);
        L.set(0,0,0.0444);
        L.set(1,0,0.500000000000000);
        N.set(0,0,15);
        pfqnMVAReturn ret = PFQN.pfqn_mva(L, N, Z, mi);

        assertEquals(-10.304214841550015, ret.lGN, tolerance);

        assertEquals(0.048726953467954, ret.CN.get(0,0), tolerance);
        assertEquals(7.451273046532047, ret.CN.get(1,0), tolerance);

        assertEquals(0.088800000000000, ret.UN.get(0,0), tolerance);
        assertEquals(1.000000000000000, ret.UN.get(1,0), tolerance);

        assertEquals(0.097453906935908, ret.QN.get(0,0), tolerance);
        assertEquals(14.902546093064091, ret.QN.get(1,0), tolerance);

        assertEquals(2.000000000000000, ret.XN.get(0,0), tolerance);
    }

    @Test
    void pfqn_xzgsbupTest1(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,1.599585062240664);
        L.set(1,0,1.931535269709544);
        L.set(2,0,0.827800829875519);
        L.set(3,0,1.931535269709540);
        double N = 100;
        double Z = 4.041493775933610;
        double x = PFQN.pfqn_xzgsbup(L, N, Z);
        assertEquals(0.512489589584437, x, tolerance);
    }

    @Test
    void pfqn_xzgsbupTest2(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.814102564102564);
        L.set(1,0,0.722934472934473);
        double N = 10;
        double Z = 1.632478632478633;
        double x = PFQN.pfqn_xzgsbup(L, N, Z);
        assertEquals(1.165952105623232, x, tolerance);
    }

    @Test
    void pfqn_xzgsbupTest3(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.8141);
        L.set(1,0,Double.NaN);
        double N = 10;
        double Z = 1.2;
        double x = PFQN.pfqn_xzgsbup(L, N, Z);
        assertEquals(Double.NaN, x);
    }

    @Test
    void pfqn_xzgsbupTest4(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.450000000000000);
        double N = 10;
        double Z = 1.200000000000000;
        double x = PFQN.pfqn_xzgsbup(L, N, Z);
        assertEquals(2.222222222222223, x, tolerance);
    }

    @Test
    void pfqn_xzgsbupTest5(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,1.599585062240664);
        L.set(1,0,1.931535269709544);
        L.set(2,0,0.827800829875519);
        L.set(3,0,1.931535269709540);
        double N = 100;
        double Z = Double.NaN;
        double x = PFQN.pfqn_xzgsbup(L, N, Z);
        assertEquals(Double.NaN, x, tolerance);
    }

    @Test
    void pfqn_qzgbupTest1(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,1.5996);
        L.set(1,0,1.9315);
        L.set(2,0,0.8278);
        L.set(3,0,1.9315);
        double N = 99;
        double Z = 4.0415;
        int i = 0;
        assertEquals(4.819523915279272, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance);
    }
    @Test
    void pfqn_qzgbupTest2(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,1.5996);
        L.set(1,0,1.9315);
        L.set(2,0,0.8278);
        L.set(3,0,1.9315);
        double N = 99;
        double Z = 4.0415;
        int i = 2;
        assertEquals(0.750022651082722, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance);
    }

    @Test
    void pfqn_qzgbupTest3(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.8141);
        L.set(1,0,0.7229);
        double N = 9;
        double Z = 1.6325;
        int i = 1;
        assertEquals(4.060176885695963, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance);
    }

    @Test
    void pfqn_qzgbupTest4(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.8141);
        L.set(1,0,0.7229);
        double N = 10;
        double Z = 1.6325;
        int i = 0;
        assertEquals(7.649396164358148, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance);
    }

    @Test
    void pfqn_qzgbupTest5(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.45);
        L.set(1,0,Double.NaN);
        double N = 10;
        double Z = 1.2;
        int i = 0;
        assertEquals(10, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance);
    }

    @Test
    void pfqn_qzgbupTest6(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.45);
        L.set(1,0,Double.NaN);
        double N = 10;
        double Z = 1.2;
        int i = 1;
        assertEquals(10, PFQN.pfqn_qzgbup(L, N, Z, i), tolerance * 10);
    }

    @Test
    void pfqn_xzabaupTest1(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,1.5996);
        L.set(1,0,1.9315);
        L.set(2,0,0.8278);
        L.set(3,0,1.9315);
        double N = 98;
        double Z = 4.0415;
        assertEquals(0.517732332384157, PFQN.pfqn_xzabaup(L, N, Z), tolerance);
    }

    @Test
    void pfqn_xzabaupTest2(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.8141);
        L.set(1,0,0.7229);
        double N = 8;
        double Z = 1.6325;
        assertEquals(1.228350325512836, PFQN.pfqn_xzabaup(L, N, Z), tolerance);
    }

    @Test
    void pfqn_xzabaupTest3(){
        Matrix L = new Matrix(2,1);
        L.set(0,0,0.45);
        L.set(1,0,Double.NaN);
        double N = 10;
        double Z = 1.2;
        assertEquals(2.222222222222222, PFQN.pfqn_xzabaup(L, N, Z), tolerance);
    }

    @Test
    void pfqn_xzgsblowTest1(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.45);
        double N = 10;
        double Z = 1;
        assertEquals(2.162551897779923, PFQN.pfqn_xzgsblow(L, N, Z), tolerance);
    }

    @Test
    void pfqn_xzgsblowTest2(){
        Matrix L = new Matrix(2, 1);
        L.set(0,0,0.4667);
        L.set(1,0,0.2222);
        double N = 10;
        double Z = 1.5278;
        assertEquals(1.987795408944852, PFQN.pfqn_xzgsblow(L, N, Z), tolerance);
    }

    @Test
    void pfqn_xzgsblowTest3(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,0.2155);
        L.set(1,0,0.2874);
        L.set(2,0,0.2634);
        L.set(3,0,0.2351);
        double N = 10;
        double Z = 2.1880;
        assertEquals(2.268014636747800, PFQN.pfqn_xzgsblow(L,N,Z), tolerance);
    }

    @Test
    void pfqn_qzgblowTest1(){
        Matrix L = new Matrix(2, 1);
        L.set(0,0,0.4667);
        L.set(1,0,0.2222);
        double N = 9;
        double Z = 1.5278;
        int i = 1;
        assertEquals(0.452717683849012, PFQN.pfqn_qzgblow(L,N,Z,i), tolerance);
    }

    @Test
    void pfqn_qzgblowTest2(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,0.2155);
        L.set(1,0,0.2874);
        L.set(2,0,0.2634);
        L.set(3,0,0.2351);
        double N = 9;
        double Z = 2.1880;
        int i = 0;
        assertEquals(0.505511466823690, PFQN.pfqn_qzgblow(L,N,Z,i), tolerance);
    }

    @Test
    void pfqn_qzgblowTest3(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,0.2155);
        L.set(1,0,0.2874);
        L.set(2,0,0.2634);
        L.set(3,0,0.2351);
        double N = 9;
        double Z = 2.1880;
        int i = 2;
        assertEquals(0.695899654061226, PFQN.pfqn_qzgblow(L,N,Z,i), tolerance);
    }

    @Test
    void pfqn_qzgblowTest4(){
        Matrix L = new Matrix(4,1);
        L.set(0,0,0.2155);
        L.set(1,0,0.2874);
        L.set(2,0,0.2634);
        L.set(3,0,0.2351);
        double N = 9;
        double Z = 2.1880;
        int i = 3;
        assertEquals(0.578030281656193, PFQN.pfqn_qzgblow(L,N,Z,i), tolerance);
    }

    @Test
    void pfqn_mvaldmx_ecTest1(){
        Matrix lambda = new Matrix(1,1);
        Matrix D = new Matrix(2,1);
        D.set(0,0,1);
        D.set(1,0,1.5);
        Matrix mu = new Matrix(2,17);
        mu.set(0,0, 1);
        mu.set(0,1, 2);
        mu.set(0,2, 3);
        mu.set(0,3, 4);
        mu.set(0,4, 5);
        mu.set(0,5, 6);
        mu.set(0,6, 7);
        mu.set(0,7, 8);
        mu.set(0,8, 9);
        mu.set(0,9, 10);
        mu.set(0,10, 11);
        mu.set(0,11, 12);
        mu.set(0,12, 13);
        mu.set(0,13, 14);
        mu.set(0,14, 15);
        mu.set(0,15, 16);
        mu.set(0,16, 16);

        mu.set(1,0, 1);
        mu.set(1,1, 2);
        mu.set(1,2, 2);
        mu.set(1,3, 2);
        mu.set(1,4, 2);
        mu.set(1,5, 2);
        mu.set(1,6, 2);
        mu.set(1,7, 2);
        mu.set(1,8, 2);
        mu.set(1,9, 2);
        mu.set(1,10, 2);
        mu.set(1,11, 2);
        mu.set(1,12, 2);
        mu.set(1,13, 2);
        mu.set(1,14, 2);
        mu.set(1,15, 2);
        mu.set(1,16, 2);
        /*
        * Note: the result returned by pfqn_mvaldmx_ec in this case is very large. I checked it manually by using a
        * debugger and verifying the values of the variables against the values obtained in MATLAB
        * */
        pfqnMVALDMXECReturn ret =  pfqn_mvaldmx_ec(lambda, D, mu);
        //System.out.println(ret);
    }

    @Test
    void pfqn_mvaldmx_ecTest2(){
        Matrix lambda = new Matrix(1,2);
        Matrix D = new Matrix(1,2);
        Matrix mu = new Matrix(1, 3);
        lambda.set(0,1,0.1);
        D.set(0,0,0.19);
        D.set(0,1,1);
        mu.set(0,0,1);
        mu.set(0,1,2);
        mu.set(0,2,2);
        pfqnMVALDMXECReturn ret =  pfqn_mvaldmx_ec(lambda, D, mu);
        assertEquals(1.00250626566416, ret.EC.get(0,0), tolerance);
        assertEquals(0.526315789473684, ret.EC.get(0,1), tolerance);
        assertEquals(0.526315789473684, ret.EC.get(0,2), tolerance);

        assertEquals(1.105263157894737, ret.E.get(0,0), tolerance);
        assertEquals(1.108033240997230, ret.E.get(0,1), tolerance);
        assertEquals(1.166350779997084, ret.E.get(0,2), tolerance);
        assertEquals(1.227737663154826, ret.E.get(0,3), tolerance);

        assertEquals(1.052631578947368, ret.Eprime.get(0,0), tolerance);
        assertEquals(0.554016620498615, ret.Eprime.get(0,1), tolerance);
        assertEquals(0.583175389998542, ret.Eprime.get(0,2), tolerance);
        assertEquals(0.613868831577413, ret.Eprime.get(0,3), tolerance);

        assertEquals(0.1, ret.Lo.get(0,0), tolerance);
    }

    @Test
    void pfqn_mvaldmx_ecTest3(){
        Matrix lambda = new Matrix(1,2);
        Matrix D = new Matrix(2,2);
        Matrix mu = new Matrix(2, 4);
        lambda.set(0,1,0.3);
        D.set(0,0,1);
        D.set(0,1,1);
        D.set(1,0,0.5);
        D.set(1,1,0.7);
        mu.set(0,0,1);
        mu.set(0,1,1);
        mu.set(0,2,1);
        mu.set(0,3,1);
        mu.set(1,0,1);
        mu.set(1,1,2);
        mu.set(1,2,2);
        mu.set(1,3,2);

        /*
         * Checking this testcase yields 30 assertions. I will check it manually once again.
         */
        pfqnMVALDMXECReturn ret =  pfqn_mvaldmx_ec(lambda, D, mu);
        //System.out.println(ret);
    }

    @Test
    void pfqn_mvaldmxTest1(){
        Matrix lambda = new Matrix(1,1);
        Matrix D = new Matrix(2,1);
        D.set(0,0,1);
        D.set(1,0,1.5);
        Matrix N = new Matrix(1,1);
        N.set(0,0,16);
        Matrix Z = new Matrix(1,1);
        Matrix mu = new Matrix(2,16);
        mu.set(0,0, 1);
        mu.set(0,1, 2);
        mu.set(0,2, 3);
        mu.set(0,3, 4);
        mu.set(0,4, 5);
        mu.set(0,5, 6);
        mu.set(0,6, 7);
        mu.set(0,7, 8);
        mu.set(0,8, 9);
        mu.set(0,9, 10);
        mu.set(0,10, 11);
        mu.set(0,11, 12);
        mu.set(0,12, 13);
        mu.set(0,13, 14);
        mu.set(0,14, 15);
        mu.set(0,15, 16);

        mu.set(1,0, 1);
        mu.set(1,1, 2);
        mu.set(1,2, 2);
        mu.set(1,3, 2);
        mu.set(1,4, 2);
        mu.set(1,5, 2);
        mu.set(1,6, 2);
        mu.set(1,7, 2);
        mu.set(1,8, 2);
        mu.set(1,9, 2);
        mu.set(1,10, 2);
        mu.set(1,11, 2);
        mu.set(1,12, 2);
        mu.set(1,13, 2);
        mu.set(1,14, 2);
        mu.set(1,15, 2);
        Matrix S = new Matrix(2,1);
        S.set(0,0,Double.POSITIVE_INFINITY);
        S.set(1,0,1);
        pfqnMVALDMXReturn ret = pfqn_mvaldmx(lambda, D, N, Z, mu, S);

        assertEquals(1.333333333322439, ret.XN.get(0,0), tolerance);

        assertEquals(1.333333333322439, ret.QN.get(0,0), tolerance);
        assertEquals(14.666666666677560, ret.QN.get(1,0), tolerance);

        assertEquals(0.736402861884074, ret.UN.get(0,0), tolerance);
        assertEquals(0.999999999999371, ret.UN.get(1,0), tolerance);

        assertEquals(1.000000000000000, ret.CN.get(0,0), tolerance);
        assertEquals(11.000000000098048, ret.CN.get(1,0), tolerance);

        assertEquals(Double.POSITIVE_INFINITY, ret.lGN);

        assertEquals(0.263597138115920, ret.Pc.get(0,0), tolerance);
        assertEquals(0.351462850821228, ret.Pc.get(0,1), tolerance);
        assertEquals(0.234308567214152, ret.Pc.get(0,2), tolerance);
        assertEquals(0.104137140984067, ret.Pc.get(0,3), tolerance);
        assertEquals(0.034712380328022, ret.Pc.get(0,4), tolerance);
        assertEquals(0.009256634754139, ret.Pc.get(0,5), tolerance);
        assertEquals(0.002057029945364, ret.Pc.get(0,6), tolerance);
        assertEquals(3.918152276884368e-04, ret.Pc.get(0,7), tolerance);
        assertEquals(6.530253794807273e-05, ret.Pc.get(0,8), tolerance);
        assertEquals(9.674450066381162e-06, ret.Pc.get(0,9), tolerance);
        assertEquals(1.289926675517486e-06, ret.Pc.get(0,10), tolerance);
        assertEquals(1.563547485475743e-07, ret.Pc.get(0,11), tolerance);
        assertEquals(1.737274983861937e-08, ret.Pc.get(0,12), tolerance);
        assertEquals(1.781820496268652e-09, ret.Pc.get(0,13), tolerance);
        assertEquals(1.696971901208239e-10, ret.Pc.get(0,14), tolerance);
        assertEquals(1.508419467740657e-11, ret.Pc.get(0,15), tolerance);
        assertEquals(6.285081115586075e-13, ret.Pc.get(0,16), tolerance);

        assertEquals(6.288303211476887e-13, ret.Pc.get(1,0), tolerance);
        assertEquals(1.508393410164398e-11, ret.Pc.get(1,1), tolerance);
        assertEquals(1.696975892968902e-10, ret.Pc.get(1,2), tolerance);
        assertEquals(1.781820022942017e-09, ret.Pc.get(1,3), tolerance);
        assertEquals(1.737275016432679e-08, ret.Pc.get(1,4), tolerance);
        assertEquals(1.563547485555210e-07, ret.Pc.get(1,5), tolerance);
        assertEquals(1.289926675371038e-06, ret.Pc.get(1,6), tolerance);
        assertEquals(9.674450066209837e-06, ret.Pc.get(1,7), tolerance);
        assertEquals(6.530253794837930e-05, ret.Pc.get(1,8), tolerance);
        assertEquals(3.918152276883405e-04, ret.Pc.get(1,9), tolerance);
        assertEquals(0.002057029945364, ret.Pc.get(1,10), tolerance);
        assertEquals(0.009256634754139, ret.Pc.get(1,11), tolerance);
        assertEquals(0.034712380328023, ret.Pc.get(1,12), tolerance);
        assertEquals(0.104137140984067, ret.Pc.get(1,13), tolerance);
        assertEquals(0.234308567214151, ret.Pc.get(1,14), tolerance);
        assertEquals(0.351462850821227, ret.Pc.get(1,15), tolerance);
        assertEquals(0.263597138115920, ret.Pc.get(1,16), tolerance);

    }

    @Test
    void pfqn_mvaldmxTest2(){
        Matrix lambda = new Matrix(1,2);
        lambda.set(0,1,0.1000);
        Matrix D = new Matrix(1,2);
        D.set(0,0,0.1900);
        D.set(0,1,1);
        Matrix N = new Matrix(1,2);
        N.set(0,0,2);
        N.set(0,1,Double.POSITIVE_INFINITY);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,0.6667);
        Z.set(0,1,0.2167);
        Matrix mu = new Matrix(1,2);
        mu.set(0,0, 1);
        mu.set(0,1, 2);
        Matrix S = new Matrix(1,1);
        S.set(0,0,2);
        pfqnMVALDMXReturn ret = pfqn_mvaldmx(lambda, D, N, Z, mu, S);

        assertEquals(2.327496139029629, ret.XN.get(0,0), tolerance);
        assertEquals(0.100000000000000, ret.XN.get(0,1), tolerance);

        assertEquals(0.448258324108947, ret.QN.get(0,0), tolerance);
        assertEquals(0.104960398558463, ret.QN.get(0,1), tolerance);

        assertEquals(0.376711289073051, ret.UN.get(0,0), tolerance);
        assertEquals(0.050000000000000, ret.UN.get(0,1), tolerance);

        assertEquals(0.192592510291266, ret.CN.get(0,0), tolerance);
        assertEquals(1.049603985584633, ret.CN.get(0,1), tolerance);

        assertEquals(Double.POSITIVE_INFINITY, ret.lGN);

        assertEquals(0.603461800975736, ret.Pc.get(0,0), tolerance);
        assertEquals(0.344818073939581, ret.Pc.get(0,1), tolerance);
        assertEquals(0.051720125084683, ret.Pc.get(0,2), tolerance);
    }

    @Test
    void pfqn_mvaldmsTest1(){
        Matrix lambda = new Matrix(1,2);
        lambda.set(0,1,0.1000);
        Matrix D = new Matrix(1,2);
        D.set(0,0,0.1900);
        D.set(0,1,1);
        Matrix N = new Matrix(1,2);
        N.set(0,0,2);
        N.set(0,1,Double.POSITIVE_INFINITY);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,0.6667);
        Z.set(0,1,0.2167);
        Matrix S = new Matrix(1,1);
        S.set(0,0,2);

        pfqnMVAReturn ret = pfqn_mvaldms(lambda, D, N, Z, S);

        assertEquals(2.327496139029629, ret.XN.get(0,0), tolerance);
        assertEquals(0.100000000000000, ret.XN.get(0,1), tolerance);

        assertEquals(0.448258324108947, ret.QN.get(0,0), tolerance);
        assertEquals(0.104960398558463, ret.QN.get(0,1), tolerance);

        assertEquals(0.221112133207815, ret.UN.get(0,0), tolerance);
        assertEquals(0.050000000000000, ret.UN.get(0,1), tolerance);

        assertEquals(0.192592510291266, ret.CN.get(0,0), tolerance);
        assertEquals(1.049603985584633, ret.CN.get(0,1), tolerance);

        assertEquals(Double.POSITIVE_INFINITY, ret.lGN);
    }

    @Test
    void pfqn_mvaldmsTest2(){
        Matrix lambda = new Matrix(1,2);
        lambda.set(0,1,0.3000);
        Matrix D = new Matrix(5,2);
        D.set(0,0,1);
        D.set(0,1,1);
        D.set(1,0,0.5);
        D.set(1,1,0.7071);
        D.set(2,0,0.3333);
        D.set(2,1,0.5774);
        D.set(3,0,0.2500);
        D.set(4,1,0.4472);
        Matrix N = new Matrix(1,2);
        N.set(0,0,3);
        N.set(0,1,Double.POSITIVE_INFINITY);
        Matrix Z = new Matrix(1,2);
        Matrix S = new Matrix(5,1);
        S.set(0,0,1);
        S.set(1,0,2);
        S.set(2,0,3);
        S.set(3,0,4);
        S.set(4,0,5);

        pfqnMVAReturn ret = pfqn_mvaldms(lambda, D, N, Z, S);

        assertEquals(0.672904140932310, ret.XN.get(0,0), tolerance);
        assertEquals(0.300000000000000, ret.XN.get(0,1), tolerance);

        assertEquals(2.250509620423859, ret.QN.get(0,0), tolerance);
        assertEquals(1.393075551610225, ret.QN.get(0,1), tolerance);

        assertEquals(0.356561446006589, ret.QN.get(1,0), tolerance);
        assertEquals(0.228239460939755, ret.QN.get(1,1), tolerance);

        assertEquals(0.224702898336474, ret.QN.get(2,0), tolerance);
        assertEquals(0.173696552633359, ret.QN.get(2,1), tolerance);

        assertEquals(0.168226035233078, ret.QN.get(3,0), tolerance);
        assertEquals(0, ret.QN.get(3,1), tolerance);

        assertEquals(0, ret.QN.get(4,0), tolerance);
        assertEquals(0.134177246123306, ret.QN.get(4,1), tolerance);

        assertEquals(0.672904140932310, ret.UN.get(0,0), tolerance);
        assertEquals(0.300000000000000, ret.UN.get(0,1), tolerance);

        assertEquals(0.168226035233078, ret.UN.get(1,0), tolerance);
        assertEquals(0.106065000000000, ret.UN.get(1,1), tolerance);

        assertEquals(0.074759650057580, ret.UN.get(2,0), tolerance);
        assertEquals(0.057740000000000, ret.UN.get(2,1), tolerance);

        assertEquals(0.042056508808269, ret.UN.get(3,0), tolerance);
        assertEquals(0, ret.UN.get(3,1), tolerance);

        assertEquals(0, ret.UN.get(4,0), tolerance);
        assertEquals(0.026832000000000, ret.UN.get(4,1), tolerance);

        assertEquals(3.344472835767919, ret.CN.get(0,0), tolerance);
        assertEquals(4.643585172034085, ret.CN.get(0,1), tolerance);

        assertEquals(0.529884458004037, ret.CN.get(1,0), tolerance);
        assertEquals(0.760798203132517, ret.CN.get(1,1), tolerance);

        assertEquals(0.333930027574429, ret.CN.get(2,0), tolerance);
        assertEquals(0.578988508777863, ret.CN.get(2,1), tolerance);

        assertEquals(0.250000000000000, ret.CN.get(3,0), tolerance);
        assertEquals(0, ret.CN.get(3,1), tolerance);

        assertEquals(0, ret.CN.get(4,0), tolerance);
        assertEquals(0.447257487077687, ret.CN.get(4,1), tolerance);

        assertEquals(Double.POSITIVE_INFINITY, ret.lGN);
    }

    @Test
    void pfqn_mvaldTest1(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.450000000000000);
        Matrix N = new Matrix(1,1);
        N.set(0,0,10);
        Matrix Z = new Matrix(1,1);
        Z.set(0,0,1);
        Matrix mu = new Matrix(1,10);
        mu.set(0,0,1);
        mu.set(0,1,2);
        mu.set(0,2,2);
        mu.set(0,3,2);
        mu.set(0,4,2);
        mu.set(0,5,2);
        mu.set(0,6,2);
        mu.set(0,7,2);
        mu.set(0,8,2);
        mu.set(0,9,2);
        pfqnMVALDReturn ret = pfqn_mvald(L, N, Z, mu);

        assertEquals(4.373375855772732, ret.XN.get(0,0), tolerance);

        assertEquals(5.626624144227267, ret.QN.get(0,0), tolerance);

        assertEquals(0.995079866938112, ret.UN.get(0,0), tolerance);

        assertEquals(1.286563133328750, ret.CN.get(0,0), tolerance);

        assertEquals(-9.789992869331030, ret.lGN.get(ret.lGN.size() - 1), tolerance);

        assertTrue(ret.isNumStable);

        assertEquals(0.004920133061888, ret.pi.get(0,0), tolerance);
        assertEquals(0.022140598778495, ret.pi.get(0,1), tolerance);
        assertEquals(0.044834712526452, ret.pi.get(0,2), tolerance);
        assertEquals(0.080702482547615, ret.pi.get(0,3), tolerance);
        assertEquals(0.127106410012493, ret.pi.get(0,4), tolerance);
        assertEquals(0.171593653516866, ret.pi.get(0,5), tolerance);
        assertEquals(0.193042860206474, ret.pi.get(0,6), tolerance);
        assertEquals(0.173738574185827, ret.pi.get(0,7), tolerance);
        assertEquals(0.117273537575433, ret.pi.get(0,8), tolerance);
        assertEquals(0.052773091908945, ret.pi.get(0,9), tolerance);
        assertEquals(0.011873945679513, ret.pi.get(0,10), tolerance);
    }

    @Test
    void pfqn_mvamxTest1(){
        Matrix lambda = new Matrix(1,2);
        lambda.set(0,1,0.100000000000000);
        Matrix D = new Matrix(1,2);
        D.set(0,0,0.190000000000000);
        D.set(0,1,1);
        Matrix N = new Matrix(1,2);
        N.set(0,0,2);
        N.set(0,1,Double.POSITIVE_INFINITY);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,0.6667);
        Z.set(0,1,0.2167);
        Matrix mi = new Matrix(1,1);
        mi.set(0,0,1);
        pfqnMVAReturn ret = pfqn_mvamx(lambda, D, N, Z, mi);

        assertEquals(2.153819913658974, ret.XN.get(0,0), tolerance);
        assertEquals(0.100000000000000, ret.XN.get(0,1), tolerance);

        assertEquals(0.564048263563562, ret.QN.get(0,0), tolerance);
        assertEquals(0.173783140395951, ret.QN.get(0,1), tolerance);

        assertEquals(0.409225783595205, ret.UN.get(0,0), tolerance);
        assertEquals(0.100000000000000, ret.UN.get(0,1), tolerance);

        assertEquals(0.261882741443011, ret.CN.get(0,0), tolerance);
        assertEquals(1.737831403959514, ret.CN.get(0,1), tolerance);

        assertEquals(-0.897566813596126, ret.lGN, tolerance);
    }

    @Test
    void pfqn_mvamxTest2(){
        Matrix lambda = new Matrix(1,2);
        lambda.set(0,1,0.3000);
        Matrix D = new Matrix(5,2);
        D.set(0,0,1);
        D.set(0,1,1);
        D.set(1,0,0.5);
        D.set(1,1,0.7071);
        D.set(2,0,0.3333);
        D.set(2,1,0.5774);
        D.set(3,0,0.2500);
        D.set(4,1,0.4472);
        Matrix N = new Matrix(1,2);
        N.set(0,0,3);
        N.set(0,1,Double.POSITIVE_INFINITY);
        Matrix Z = new Matrix(1,2);
        Matrix mi = Matrix.ones(5,1);

        pfqnMVAReturn ret = pfqn_mvamx(lambda, D, N, Z, mi);

        assertEquals(0.624131294218873, ret.XN.get(0,0), tolerance);
        assertEquals(0.300000000000000, ret.XN.get(0,1), tolerance);

        assertEquals(1.942579022933237, ret.QN.get(0,0), tolerance);
        assertEquals(1.261105295542816, ret.QN.get(0,1), tolerance);

        assertEquals(0.563762716028655, ret.QN.get(1,0), tolerance);
        assertEquals(0.421035177061138, ret.QN.get(1,1), tolerance);

        assertEquals(0.314591483776149, ret.QN.get(2,0), tolerance);
        assertEquals(0.275422164081986, ret.QN.get(2,1), tolerance);

        assertEquals(0.179066777261960, ret.QN.get(3,0), tolerance);
        assertEquals(0, ret.QN.get(3,1), tolerance);

        assertEquals(0, ret.QN.get(4,0), tolerance);
        assertEquals(0.154947796359605, ret.QN.get(4,1), tolerance);

        assertEquals(0.624131294218873, ret.UN.get(0,0), tolerance);
        assertEquals(0.300000000000000, ret.UN.get(0,1), tolerance);

        assertEquals(0.312065647109436, ret.UN.get(1,0), tolerance);
        assertEquals(0.212130000000000, ret.UN.get(1,1), tolerance);

        assertEquals(0.208022960363150, ret.UN.get(2,0), tolerance);
        assertEquals(0.173220000000000, ret.UN.get(2,1), tolerance);

        assertEquals(0.156032823554718, ret.UN.get(3,0), tolerance);
        assertEquals(0, ret.UN.get(3,1), tolerance);

        assertEquals(0, ret.UN.get(4,0), tolerance);
        assertEquals(0.134160000000000, ret.UN.get(4,1), tolerance);

        assertEquals(3.112452525496992, ret.CN.get(0,0), tolerance);
        assertEquals(4.203684318476053, ret.CN.get(0,1), tolerance);

        assertEquals(0.903275835149763, ret.CN.get(1,0), tolerance);
        assertEquals(1.403450590203792, ret.CN.get(1,1), tolerance);

        assertEquals(0.504046963659904, ret.CN.get(2,0), tolerance);
        assertEquals(0.918073880273287, ret.CN.get(2,1), tolerance);

        assertEquals(0.286905622135275, ret.CN.get(3,0), tolerance);
        assertEquals(0, ret.CN.get(3,1), tolerance);

        assertEquals(0, ret.CN.get(4,0), tolerance);
        assertEquals(0.516492654532015, ret.CN.get(4,1), tolerance);

        assertEquals(2.085520693339201, ret.lGN, tolerance);
    }

    @Test
    void pfqn_mvamsTest1(){
        Matrix lambda = new Matrix(1,2);
        Matrix L = new Matrix(2,2);
        L.set(0,0,10);
        L.set(0,1,5);
        L.set(1,0,5);
        L.set(1,1,9);
        Matrix N = new Matrix(1,2);
        N.set(0,0,1);
        N.set(0,1,2);
        Matrix Z = new Matrix(2,2);
        Z.set(0,0,91);
        Z.set(0,1,92);
        Z.set(1,0,93);
        Z.set(1,1,94);
        Matrix mi = new Matrix(2,1);
        mi.set(0,0,1);
        mi.set(1,0,1);
        Matrix S = Matrix.ones(2,1);
        pfqnMVAReturn ret = pfqn_mvams(lambda, L, N, Z, mi, S);

        assertEquals(0.009270178985455, ret.XN.get(0,0), tolerance);
        assertEquals(0.018351023289770, ret.XN.get(0,1), tolerance);

        assertEquals(0.101687164552697, ret.QN.get(0,0), tolerance);
        assertEquals(0.105393631615683, ret.QN.get(0,1), tolerance);
        assertEquals(0.054726547770905, ret.QN.get(1,0), tolerance);
        assertEquals(0.187961202435678, ret.QN.get(1,1), tolerance);

        assertEquals(0.092701789854549, ret.UN.get(0,0), tolerance);
        assertEquals(0.091755116448852, ret.UN.get(0,1), tolerance);
        assertEquals(0.046350894927275, ret.UN.get(1,0), tolerance);
        assertEquals(0.165159209607933, ret.UN.get(1,1), tolerance);

        assertEquals(10.969277369104285, ret.CN.get(0,0), tolerance);
        assertEquals(5.743201888607152, ret.CN.get(0,1), tolerance);
        assertEquals(5.903504976200779, ret.CN.get(1,0), tolerance);
        assertEquals(10.242546122234852, ret.CN.get(1,1), tolerance);

        assertEquals(13.342678933707491, ret.lGN, tolerance);
    }

    @Test
    void pfqn_mvamsTest2(){
        Matrix lambda = new Matrix(1,1);
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.4161);
        Matrix N = new Matrix(1,1);
        N.set(0,0,4);
        Matrix Z = new Matrix(1,1);
        Z.set(0,0,0.4845);
        Matrix mi = new Matrix(1,1);
        mi.set(0,0,1);
        Matrix S = new Matrix(1,1);
        S.set(0,0,1);
        pfqnMVAReturn ret = pfqn_mvams(lambda, L, N, Z, mi, S);

        assertEquals(2.345421806267623, ret.XN.get(0,0), tolerance);

        assertEquals(2.863643134863337, ret.QN.get(0,0), tolerance);

        assertEquals(0.975930013587958, ret.UN.get(0,0), tolerance);

        assertEquals(1.220950162231323, ret.CN.get(0,0), tolerance);

        assertEquals(-2.349815629247493, ret.lGN, tolerance);
    }

    @Test
    void pfqn_linearizermsTest1(){
        Matrix L = new Matrix(2,2);
        L.set(0,0,1);
        L.set(0,1,1);
        L.set(1,0,1);
        L.set(1,1,0.1);
        Matrix N = new Matrix(1,2);
        N.set(0,0,4);
        N.set(0,1,2);
        Matrix Z = new Matrix(1,2);
        Matrix nservers = new Matrix(2,1);
        nservers.set(0,0,1);
        nservers.set(1,0,3);
        Matrix type = new Matrix(2,1);
        type.set(0,0, SchedStrategy.toID(SchedStrategy.PS));
        type.set(1,0, SchedStrategy.toID(SchedStrategy.FCFS));
        double tol = 1e-04;
        int maxiter = 100;
        pfqnLinearizerMSReturn ret = pfqn_linearizerms(L, N, Z, nservers, type, tol, maxiter);

        assertEquals(2.080508946707209, ret.Q.get(0,0), tolerance);
        assertEquals(1.799386597325419, ret.Q.get(0,1), tolerance);
        assertEquals(1.919491053292791, ret.Q.get(1,0), tolerance);
        assertEquals(0.200613402674581, ret.Q.get(1,1), tolerance);

        assertEquals(0.494554229243189, ret.U.get(0,0), tolerance);
        assertEquals(0.493217155141877, ret.U.get(0,1), tolerance);
        assertEquals(0.164851409747730, ret.U.get(1,0), tolerance);
        assertEquals(0.016440571838063, ret.U.get(1,1), tolerance);

        assertEquals(4.206836831404699, ret.R.get(0,0), tolerance);
        assertEquals(3.648264417744783, ret.R.get(0,1), tolerance);
        assertEquals(3.881254956064513, ret.R.get(1,0), tolerance);
        assertEquals(0.406744576061782, ret.R.get(1,1), tolerance);

        assertEquals(8.088091787469212, ret.C.get(0,0), tolerance);
        assertEquals(4.055008993806564, ret.C.get(0,1), tolerance);

        assertEquals(0.494554229243189, ret.X.get(0,0), tolerance);
        assertEquals(0.493217155141877, ret.X.get(0,1), tolerance);

        assertEquals(39, ret.totiter);
    }

    /*
    @Test
    void pfqn_linearizermsTest2(){
        Matrix L = new Matrix(2,2);
        L.set(0,0,1);
        L.set(0,1,0.3833);
        L.set(1,1,0.1667);
        Matrix N = new Matrix(1,2);
        N.set(0,0,4);
        N.set(0,1,3);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,1);
        Z.set(0,1,0.4);
        Matrix nservers = new Matrix(2,1);
        nservers.set(0,0,3);
        nservers.set(1,0,3);
        Matrix type = new Matrix(2,1);
        type.set(0,0, SchedStrategy.toID(SchedStrategy.FCFS));
        type.set(1,0, SchedStrategy.toID(SchedStrategy.FCFS));
        double tol = 1e-04;
        int maxiter = 100;
        pfqnLinearizerMSReturn ret = pfqn_linearizerms(L, N, Z, nservers, type, tol, maxiter);

        assertEquals(2.755492897438664, ret.Q.get(0,0), tolerance);
        assertEquals(2.202320856357018, ret.Q.get(0,1), tolerance);
        assertEquals(0, ret.Q.get(1,0), tolerance);
        assertEquals(0.385787040986767, ret.Q.get(1,1), tolerance);

        assertEquals(0.414835700853779, ret.U.get(0,0), tolerance);
        assertEquals(0.131565202456773, ret.U.get(0,1), tolerance);
        assertEquals(0, ret.U.get(1,0), tolerance);
        assertEquals(0.057218677927326, ret.U.get(1,1), tolerance);

        assertEquals(2.214123882272387, ret.R.get(0,0), tolerance);
        assertEquals(2.138735695250930, ret.R.get(0,1), tolerance);
        assertEquals(0, ret.R.get(1,0), tolerance);
        assertEquals(0.374648640747321, ret.R.get(1,1), tolerance);

        assertEquals(2.214123882272387, ret.C.get(0,0), tolerance);
        assertEquals(2.513384335998251, ret.C.get(0,1), tolerance);

        assertEquals(1.244507102561336, ret.X.get(0,0), tolerance);
        assertEquals(1.029730256640537, ret.X.get(0,1), tolerance);

        assertEquals(67, ret.totiter);
    }
    */
    @Test
    void pfqn_bsTest1(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.45);
        Matrix N = new Matrix(1,1);
        N.set(0,0,10);
        Matrix Z = new Matrix(1,1);
        Z.set(0,0,1);
        double tol = 1.0e-04;
        int maxiter = 100;
        Matrix QN0 = new Matrix(0,0);
        SchedStrategy[] type = new SchedStrategy[1];
        type[0] = SchedStrategy.FCFS;
        pfqnBSReturn ret = pfqn_bs(L, N, Z, tol, maxiter, QN0, type);
        assertEquals(2.162470968879477, ret.XN.get(0,0), tolerance);
        assertEquals(7.837529031120526, ret.QN.get(0,0), tolerance);
        assertEquals(0.973111935995765, ret.UN.get(0,0), tolerance);
        assertEquals(3.624339537460557, ret.RN.get(0,0), tolerance);
        assertEquals(6, ret.it);
    }

    @Test
    void pfqn_bsTest2(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.4161);
        Matrix N = new Matrix(1,1);
        N.set(0,0,4);
        Matrix Z = new Matrix(1,1);
        Z.set(0,0,0.4845);
        double tol = 1.0e-04;
        int maxiter = 100;
        Matrix QN0 = new Matrix(0,0);
        SchedStrategy[] type = new SchedStrategy[1];
        type[0] = SchedStrategy.PS;
        pfqnBSReturn ret = pfqn_bs(L, N, Z, tol, maxiter, QN0, type);
        assertEquals(2.202783207559328, ret.XN.get(0,0), tolerance);
        assertEquals(2.932751535937506, ret.QN.get(0,0), tolerance);
        assertEquals(0.916578092665436, ret.UN.get(0,0), tolerance);
        assertEquals(1.331384552902498, ret.RN.get(0,0), tolerance);
        assertEquals(6, ret.it);
    }

    @Test
    void pfqn_bsTest3(){
        Matrix L = new Matrix(1,2);
        L.set(0,0,0.1068);
        L.set(0,1,0.1667);
        Matrix N = new Matrix(1,2);
        N.set(0,0,2);
        N.set(0,1,1);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,0.5914);
        Z.set(0,1,1.1667);
        double tol = 1.0e-04;
        int maxiter = 100;
        Matrix QN0 = new Matrix(0,0);
        SchedStrategy[] type = new SchedStrategy[1];
        type[0] = SchedStrategy.PS;
        pfqnBSReturn ret = pfqn_bs(L, N, Z, tol, maxiter, QN0, type);
        assertEquals(2.713478142450679, ret.XN.get(0,0), tolerance);
        assertEquals(0.714646929614981, ret.XN.get(0,1), tolerance);
        assertEquals(0.395249026554668, ret.QN.get(0,0), tolerance);
        assertEquals(0.166221427218201, ret.QN.get(0,1), tolerance);
        assertEquals(0.289799465613733, ret.UN.get(0,0), tolerance);
        assertEquals(0.119131643166817, ret.UN.get(0,1), tolerance);
        assertEquals(0.145661400536545, ret.RN.get(0,0), tolerance);
        assertEquals(0.232592375801227, ret.RN.get(0,1), tolerance);
        assertEquals(8, ret.it);
    }

    @Test
    void pfqn_linearizerTest1(){
        Matrix L = new Matrix(1,1);
        L.set(0,0,0.45);
        Matrix N = new Matrix(1,1);
        N.set(0,0,10);
        Matrix Z = new Matrix(1,1);
        Z.set(0,0,1);
        SchedStrategy[] type = new SchedStrategy[1];
        type[0] = SchedStrategy.FCFS;
        double tol = 1.0e-04;
        int maxiter = 100;
        pfqnLinearizerReturn ret = pfqn_linearizer(L, N, Z, type, tol, maxiter);
        assertEquals(7.782194880076756, ret.Q.get(0,0), tolerance);
        assertEquals(0.998012303965460, ret.U.get(0,0), tolerance);
        assertEquals(3.508962446775345, ret.W.get(0,0), tolerance);
        assertEquals(3.508962446775344, ret.C.get(0,0), tolerance);
        assertEquals(2.217805119923245, ret.X.get(0,0), tolerance);
        assertEquals(32, ret.totiter);
    }

    @Test
    void pfqn_linearizerTest2(){
        Matrix L = new Matrix(1,2);
        L.set(0,0,0.2137);
        L.set(0,1,0.3333);
        Matrix N = new Matrix(1,2);
        N.set(0,0,2);
        N.set(0,1,1);
        Matrix Z = new Matrix(1,2);
        Z.set(0,0,0.4845);
        Z.set(0,1,1);
        SchedStrategy[] type = new SchedStrategy[1];
        type[0] = SchedStrategy.PS;
        double tol = 1.0e-04;
        int maxiter = 100;
        pfqnLinearizerReturn ret = pfqn_linearizer(L, N, Z, type, tol, maxiter);
        assertEquals(0.846623318708752, ret.Q.get(0,0), tolerance);
        assertEquals(0.381067955453777, ret.Q.get(0,1), tolerance);
        assertEquals(0.508723625989556, ret.U.get(0,0), tolerance);
        assertEquals(0.206290050447256, ret.U.get(0,1), tolerance);
        assertEquals(0.355641833728742, ret.W.get(0,0), tolerance);
        assertEquals(0.615686259600861, ret.W.get(0,1), tolerance);
        assertEquals(0.355641833728742, ret.C.get(0,0), tolerance);
        assertEquals(0.615686259600861, ret.C.get(0,1), tolerance);
        assertEquals(2.380550425781730, ret.X.get(0,0), tolerance);
        assertEquals(0.618932044546223, ret.X.get(0,1), tolerance);
        assertEquals(109, ret.totiter);
    }

    @Test
    void pfqn_rdTest1(){
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 16);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 2.0 / 3);
        L.set(1, 0, 1.0);
        N.set(0, 0, 16);
        Z.set(0, 0, 0);
        for (int i=0; i<16; i++) {
            mu.set(0, i, i + 1);
        }
        mu.set(1, 0, 1);
        for (int i=1; i<16; i++) {
            mu.set(1, i, 2);
        }
        pfqnRdReturn ret = pfqn_rd(L, N, Z, mu, options);

        assertEquals(-9.120855161272257, ret.lGN, tolerance);

        assertEquals(6.569832586711733, ret.Cgamma, tolerance);

    }

    @Test
    void pfqn_rdTest2(){
        Matrix L = new Matrix(2, 2);
        Matrix N = new Matrix(1, 2);
        Matrix Z = new Matrix(1, 2);
        Matrix mu = new Matrix(2, 6);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 2.0 / 3);
        L.set(0, 1, 0.8);
        L.set(1, 0, 1);
        L.set(1, 1, 1);
        N.set(0, 0, 4);
        N.set(0, 1, 2);
        Z.set(0, 0, 0);
        Z.set(0, 1, 0);
        for (int i=0; i<6; i++) {
            mu.set(0, i, i + 1);
        }
        mu.set(1, 0, 1);
        for (int i=1; i<6; i++) {
            mu.set(1, i, 2);
        }

        pfqnRdReturn ret = pfqn_rd(L, N, Z, mu, options);

        assertEquals(0.526401647774762, ret.lGN, tolerance);

        assertEquals(5.511730546381030, ret.Cgamma, tolerance);

    }

    @Test
    void pfqn_rdTest3(){
        Matrix L = new Matrix(3, 2);
        Matrix N = new Matrix(1, 2);
        Matrix Z = new Matrix(1, 2);
        Matrix mu = new Matrix(3, 6);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 2.0 / 7);
        L.set(0, 1, 4.0 / 9);
        L.set(1, 0, 3.0 / 7);
        L.set(1, 1, 5.0 / 9);
        L.set(2, 0, 1);
        L.set(2, 1, 1);
        N.set(0, 0, 4);
        N.set(0, 1, 2);
        Z.set(0, 0, 0);
        Z.set(0, 1, 0);
        mu.set(0, 0, 1);
        mu.set(0, 1, 2);
        mu.set(0, 2, 3);
        mu.set(0, 3, 4);
        mu.set(0, 4, 5);
        mu.set(0, 5, 6);
        mu.set(1, 0, 1);
        mu.set(1, 1, 2);
        mu.set(1, 2, 3);
        mu.set(1, 3, 3);
        mu.set(1, 4, 3);
        mu.set(1, 5, 3);
        mu.set(2, 0, 1);
        mu.set(2, 1, 2);
        mu.set(2, 2, 3);
        mu.set(2, 3, 3);
        mu.set(2, 4, 3);
        mu.set(2, 5, 3);

        pfqnRdReturn ret = pfqn_rd(L, N, Z, mu, options);

        assertEquals(-0.264097214680063, ret.lGN, tolerance);

        assertEquals(16.566171738940422, ret.Cgamma, tolerance);
    }

    @Test
    void pfqn_rdTest4(){
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 10);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 1);
        L.set(1, 0, 1);
        N.set(0, 0, 10);
        Z.set(0, 0, 0);
        for (int i=0; i<10; i++) {
            mu.set(0, i, i + 1);
        }
        for (int i=0; i<10; i++) {
            mu.set(1, i, Math.min(5,  i + 1));
        }

        pfqnRdReturn ret = pfqn_rd(L, N, Z, mu, options);

        assertEquals(-8.635479492422949, ret.lGN, tolerance);

        assertEquals(8.680426521339080e+02, ret.Cgamma, tolerance);
    }

    @Test
    void pfqn_rdTest5(){
        Matrix L = new Matrix(2, 2);
        Matrix N = new Matrix(1, 2);
        Matrix Z = new Matrix(1, 2);
        Matrix mu = new Matrix(2, 24);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 2.0 / 3);
        L.set(0, 1, 0.8);
        L.set(1, 0, 1);
        L.set(1, 1, 1);
        N.set(0, 0, 16);
        N.set(0, 1, 8);
        Z.set(0, 0, 0);
        Z.set(0, 1, 0);
        for (int i=0; i<24; i++) {
            mu.set(0, i, i + 1);
        }
        for (int i=0; i<24; i++) {
            mu.set(1, i, 1);
        }

        pfqnRdReturn ret = pfqn_rd(L, N, Z, mu, options);

        assertEquals(14.211561934863019, ret.lGN, tolerance);

        assertEquals(1.960536983667250, ret.Cgamma, tolerance);
    }

    @Test
    void pfqn_nrl_Test1(){
        Matrix L = new Matrix(2, 1);
        Matrix N = new Matrix(1, 1);
        Matrix Z = new Matrix(1, 1);
        Matrix mu = new Matrix(2, 16);
        SolverOptions options = new SolverOptions();
        options.method = "default";
        L.set(0, 0, 2.0 / 3);
        L.set(1, 0, 1.0);
        N.set(0, 0, 16);
        Z.set(0, 0, 0);
        for (int i=0; i<16; i++) {
            mu.set(0, i, i + 1);
        }
        mu.set(1, 0, 1);
        for (int i=1; i<16; i++) {
            mu.set(1, i, 2);
        }
        double ret = pfqn_nrl(L, N, Z, mu, options);

        assertEquals(-8.838083105162221, ret, tolerance);


    }

}
