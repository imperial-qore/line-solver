package jline.util;

import org.junit.jupiter.api.Test;

import java.util.HashSet;

import static org.junit.jupiter.api.Assertions.*;

class MatrixTest {

    @Test
    void shouldCeilOneValue() {
        Matrix m = new Matrix(1, 2);
        m.set(0, 0, 1.2);
        m = m.ceil();
        assertEquals(2, m.get(0, 0));
        assertEquals(0, m.get(0, 1));
    }

    @Test
    void shouldCeilAllValues() {
        Matrix m = new Matrix(1, 2);
        m.set(0, 0, 1.2);
        m.set(0, 1, -1.2);
        m = m.ceil();
        assertEquals(2, m.get(0, 0));
        assertEquals(-1, m.get(0, 1));
    }

    @Test
    void shouldCeilNothing() {
        Matrix m = new Matrix(1, 2);
        m = m.ceil();
        assertEquals(0, m.get(0, 0));
        assertEquals(0, m.get(0, 1));
        m.set(0, 0, 5);
        m.set(0, 1, 5);
        m = m.ceil();
        assertEquals(5, m.get(0, 0));
        assertEquals(5, m.get(0, 1));
    }

    @Test
    void shouldFlattenMatrix() {
        Matrix m = new Matrix(2, 3);
        m = m.columnMajorOrder();
        assertEquals(6, m.getNumRows());
        assertEquals(1, m.getNumCols());
        for(int i = 0; i < 6; i++){
            assertEquals(0, m.get(i, 0));
        }
    }

    @Test
    void shouldFlattenRowVector() {
        Matrix m = new Matrix(1, 3);
        m.set(0, 0, 5);
        m.set(0, 1, 10);
        m.set(0,2, 15);
        m = m.columnMajorOrder();
        assertEquals(3, m.getNumRows());
        assertEquals(1, m.getNumCols());
        assertEquals(5, m.get(0, 0));
        assertEquals(10, m.get(1, 0));
        assertEquals(15, m.get(2, 0));
    }

    @Test
    void shouldNotFlatten() {
        Matrix m = new Matrix(3, 1);
        assertEquals(0, m.get(0, 0));
        assertEquals(0, m.get(1, 0));
        assertEquals(0, m.get(2, 0));
        m = m.columnMajorOrder();
        assertEquals(3, m.getNumRows());
        assertEquals(1, m.getNumCols());
        assertEquals(0, m.get(0, 0));
        assertEquals(0, m.get(1, 0));
        assertEquals(0, m.get(2, 0));
    }

    @Test
    void shouldBeEqual(){
        Matrix m = new Matrix(1,2);
        Matrix n = new Matrix(1, 2);
        assertTrue(m.isEqualTo(n));
        m.set(0, 0, 2);
        n.set(0,0,2);
        assertTrue(m.isEqualTo(n));
    }

    @Test
    void shouldNotBeEqualDifferentDimensions(){
        Matrix m = new Matrix(1,2);
        Matrix n = new Matrix(5, 2);
        assertFalse(m.isEqualTo(n));
    }

    @Test
    void shouldNotBeEqualDifferentValues(){
        Matrix m = new Matrix(1,2);
        Matrix n = new Matrix(1, 2);
        assertTrue(m.isEqualTo(n));
        m.set(0,0, 5);
        n.set(0,1,2);
        assertFalse(m.isEqualTo(n));
    }

    @Test
    void shouldNotHaveAnyNonZeroElements(){
        Matrix m = new Matrix(1, 5);
        assertFalse(m.any());
    }

    @Test
    void shouldNotHaveAnyElements(){
        Matrix m = new Matrix(0, 0);
        assertFalse(m.any());
    }

    @Test
    void shouldHaveNonZeroElements(){
        Matrix m = new Matrix(5, 1);
        m.set(3,0,5);
        assertTrue(m.any());
    }

    @Test
    void shouldThrowErrorNotVector(){
        Matrix m = new Matrix(5, 2);
        assertThrows(IllegalArgumentException.class, () -> {m.prodVector();});
    }

    @Test
    void shouldComputeZeroProduct(){
        Matrix m = new Matrix(1, 5);
        assertEquals(0, m.prodVector());
        m.set(0, 0, 5);
        assertEquals(0, m.prodVector());
    }

    @Test
    void shouldComputeNonZeroProduct(){
        Matrix m = new Matrix(3, 1);
        assertEquals(0, m.prodVector());
        m.set(0, 0, 5);
        m.set(1,0,2);
        m.set(2,0,6);
        assertEquals(60, m.prodVector());
    }

    @Test
    void shouldComputeZeroNorm(){
        Matrix m = new Matrix(3,3);
        assertEquals(0, m.norm());
    }

    @Test
    void normTestOnlyNegatives(){
        Matrix m = new Matrix(1,2);
        m.set(0,0,-3);
        m.set(0,1,-4);
        assertEquals(5, m.norm());
    }

    @Test
    void normTestOnlyPositives(){
        Matrix m = new Matrix(1,2);
        m.set(0,0,3);
        m.set(0,1,4);
        assertEquals(5, m.norm());
    }

    @Test
    void shouldNotChangeMatrixElementDiv(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        Matrix B = Matrix.ones(2,2);
        Matrix res = A.elementDiv(B);
        assertEquals(1, res.get(0,0));
        assertEquals(2, res.get(0,1));
        assertEquals(3, res.get(1,0));
        assertEquals(4, res.get(1,1));
    }

    @Test
    void shouldDivideMatrixElementDiv(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        Matrix B = Matrix.ones(2,2);
        B.set(0,0,4);
        B.set(1,1,8);
        Matrix res = A.elementDiv(B);
        assertEquals(0.25, res.get(0,0));
        assertEquals(2, res.get(0,1));
        assertEquals(3, res.get(1,0));
        assertEquals(0.5, res.get(1,1));
    }

    @Test
    void shouldProduceInfElementDiv(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        Matrix B = new Matrix(2,2);
        Matrix res = A.elementDiv(B);
        assertEquals(Double.POSITIVE_INFINITY, res.get(0,0));
        assertEquals(Double.POSITIVE_INFINITY, res.get(0,1));
        assertEquals(Double.POSITIVE_INFINITY, res.get(1,0));
        assertEquals(Double.POSITIVE_INFINITY, res.get(1,1));
    }

    @Test
    void shouldProduceNanElementDiv(){
        Matrix A = new Matrix(2,2);
        Matrix B = new Matrix(2,2);
        Matrix res = A.elementDiv(B);
        assertEquals(Double.NaN, res.get(0,0));
        assertEquals(Double.NaN, res.get(0,1));
        assertEquals(Double.NaN, res.get(1,0));
        assertEquals(Double.NaN, res.get(1,1));
    }

    @Test
    void shouldThrowErrorElementDiv(){
        Matrix A = new Matrix(2,2);
        Matrix B = new Matrix(1,3);
        assertThrows(IllegalArgumentException.class, () -> {
            A.elementDiv(B);
        });
    }

    @Test
    void shouldNotRemoveRows(){
        Matrix A = new Matrix(2,2);
        A.removeRows(new HashSet<>());
        assertEquals(2, A.getNumRows());
        assertEquals(2, A.getNumCols());
    }

    @Test
    void shouldNotRemoveRows2(){
        Matrix A = new Matrix(2,2);
        HashSet<Integer> rows = new HashSet<>();
        rows.add(3);
        assertThrows(IllegalArgumentException.class, () -> {
            A.removeRows(rows);
        });
    }

    @Test
    void shouldRemoveOneRow(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        HashSet<Integer> rows = new HashSet<>();
        rows.add(1);
        A.removeRows(rows);
        assertEquals(1, A.getNumRows());
        assertEquals(2, A.getNumCols());
        assertEquals(1, A.get(0,0));
        assertEquals(2, A.get(0,1));
    }

    @Test
    void shouldRemoveAllRows(){
        Matrix A = new Matrix(2,2);
        HashSet<Integer> rows = new HashSet<>();
        rows.add(0);
        rows.add(1);
        A.removeRows(rows);
        assertTrue(A.isEmpty());
    }

    @Test
    void shouldNotRemoveCols(){
        Matrix A = new Matrix(2,2);
        A.removeCols(new HashSet<>());
        assertEquals(2, A.getNumRows());
        assertEquals(2, A.getNumCols());
    }

    @Test
    void shouldNotRemoveCols2(){
        Matrix A = new Matrix(2,2);
        HashSet<Integer> cols = new HashSet<>();
        cols.add(3);
        assertThrows(IllegalArgumentException.class, () -> {
            A.removeCols(cols);
        });
    }

    @Test
    void shouldRemoveOneCol(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        HashSet<Integer> cols = new HashSet<>();
        cols.add(1);
        A.removeCols(cols);
        assertEquals(2, A.getNumRows());
        assertEquals(1, A.getNumCols());
        assertEquals(1, A.get(0,0));
        assertEquals(3, A.get(1,0));
    }

    @Test
    void shouldRemoveAllCols(){
        Matrix A = new Matrix(2,2);
        HashSet<Integer> cols = new HashSet<>();
        cols.add(0);
        cols.add(1);
        A.removeCols(cols);
        assertTrue(A.isEmpty());
    }

    @Test
    void shouldConcatNothing(){
        Matrix A = new Matrix(1, 2);
        Matrix B = new Matrix(1, 0);
        A.set(0,0,1);
        A.set(0,1,2);
        A = A.concatCols(B);
        assertEquals(1, A.getNumRows());
        assertEquals(2, A.getNumCols());
        assertEquals(1, A.get(0,0));
        assertEquals(2, A.get(0,1));
    }

    @Test
    void shouldNotConcatDifferentRows(){
        Matrix A = new Matrix(1, 2);
        Matrix B = new Matrix(2,1);
        assertThrows(IllegalArgumentException.class, () -> {
            A.concatCols(B);
        });
    }

    @Test
    void shouldConcatTwoColumns(){
        Matrix A = new Matrix(2,2);
        A.set(0,0,1);
        A.set(0,1,2);
        A.set(1,0,3);
        A.set(1,1,4);
        Matrix B = new Matrix(2, 2);
        B.set(0,0,5);
        B.set(0,1,6);
        B.set(1,0,7);
        B.set(1,1,8);
        A = A.concatCols(B);
        assertEquals(2, A.getNumRows());
        assertEquals(4, A.getNumCols());
        assertEquals(1, A.get(0,0));
        assertEquals(2, A.get(0,1));
        assertEquals(3, A.get(1,0));
        assertEquals(4, A.get(1,1));
        assertEquals(5, A.get(0,2));
        assertEquals(6, A.get(0,3));
        assertEquals(7, A.get(1,2));
        assertEquals(8, A.get(1,3));
    }
}
