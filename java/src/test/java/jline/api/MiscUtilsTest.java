package jline.api;

import jline.util.Maths;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class MiscUtilsTest {

    @Test
    void shouldChooseNothing() {
        Matrix A = new Matrix(1, 4);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(0,2,3);
        A.set(0,3,4);
        Matrix res = Maths.nchoosek(A, 0);
        assertTrue(res.isEmpty());
    }

    @Test
    void shouldChooseOne() {
        Matrix A = new Matrix(1, 4);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(0,2,3);
        A.set(0,3,4);
        Matrix res = Maths.nchoosek(A, 1);
        assertEquals(4, res.getNumRows());
        assertEquals(1, res.getNumCols());
        assertEquals(1, res.get(0,0));
        assertEquals(2, res.get(1,0));
        assertEquals(3, res.get(2,0));
        assertEquals(4, res.get(3,0));
    }

    @Test
    void shouldChooseTwo(){
        Matrix A = new Matrix(1, 4);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(0,2,3);
        A.set(0,3,4);
        Matrix res = Maths.nchoosek(A, 2);
        assertEquals(6, res.getNumRows());
        assertEquals(2, res.getNumCols());

        assertEquals(1, res.get(0,0));
        assertEquals(2, res.get(0,1));

        assertEquals(1, res.get(1,0));
        assertEquals(3, res.get(1,1));

        assertEquals(1, res.get(2,0));
        assertEquals(4, res.get(2,1));

        assertEquals(2, res.get(3,0));
        assertEquals(3, res.get(3,1));

        assertEquals(2, res.get(4,0));
        assertEquals(4, res.get(4,1));

        assertEquals(3, res.get(5,0));
        assertEquals(4, res.get(5,1));
    }

    @Test
    void shouldChooseAll(){
        Matrix A = new Matrix(1, 4);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(0,2,3);
        A.set(0,3,4);
        Matrix res = Maths.nchoosek(A, 4);
        assertEquals(1, res.getNumRows());
        assertEquals(4, res.getNumCols());
        assertEquals(1, res.get(0,0));
        assertEquals(2, res.get(0,1));
        assertEquals(3, res.get(0,2));
        assertEquals(4, res.get(0,3));

    }
}