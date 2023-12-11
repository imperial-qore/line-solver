package jline.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test class for the Hash utilities class. These tests were determined by using release 2.0.27 of the LINE
 * software as a baseline. These tests use as inputs and outputs the values observed through the MATLAB
 * workspace whenever the corresponding functions were called/executed.
 */
class HashTest {

    @Test
    void testHashpop1() {
        Matrix n = new Matrix(1,1);
        Matrix N = new Matrix(1,1);
        int R = 1;
        Matrix prods = new Matrix(1,1);
        N.set(0,0,16);
        prods.set(0,0,1);
        assertEquals(0, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop2() {
        Matrix n = new Matrix(1,1);
        Matrix N = new Matrix(1,1);
        int R = 1;
        Matrix prods = new Matrix(1,1);
        n.set(0,0,1);
        N.set(0,0,16);
        prods.set(0,0,1);
        assertEquals(1, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop3() {
        Matrix n = new Matrix(1,1);
        Matrix N = new Matrix(1,1);
        int R = 1;
        Matrix prods = new Matrix(1,1);
        n.set(0,0,2);
        N.set(0,0,16);
        prods.set(0,0,1);
        assertEquals(2, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop4() {
        Matrix n = new Matrix(1,1);
        Matrix N = new Matrix(1,1);
        int R = 1;
        Matrix prods = new Matrix(1,1);
        N.set(0,0,3);
        prods.set(0,0,1);
        assertEquals(0, PopulationLattice.hashpop(n, N, R, prods));
    }
}