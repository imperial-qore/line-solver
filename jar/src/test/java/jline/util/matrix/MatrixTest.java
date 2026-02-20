package jline.util.matrix;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.ejml.data.DMatrixRMaj;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.ZERO_TOL;

/**
 * Comprehensive test suite for Matrix and related classes, consolidating all matrix tests.
 * This replaces: DenseMatrixImplTest, MatrixArithmeticTest, MatrixCoreTest, 
 * MatrixCorrectnessTestSuite, MatrixViewOperationsTest, MatrixViewsAndManipulationTest, ViewTest.
 */
public class MatrixTest {
    
    private Matrix matrix2x2;
    private Matrix matrix2x3;
    private Matrix matrix3x2;
    private Matrix matrix3x3;
    private Matrix matrix1x1;
    private Matrix vector3x1;
    private Matrix rowVector1x3;
    private Matrix identityMatrix;
    private Matrix zeroMatrix;
    private DenseMatrixImpl denseMatrix;
    
    @BeforeEach
    void setUp() {
        matrix2x2 = new Matrix(2, 2);
        matrix2x2.set(0, 0, 1.0); matrix2x2.set(0, 1, 2.0);
        matrix2x2.set(1, 0, 3.0); matrix2x2.set(1, 1, 4.0);
        
        matrix2x3 = new Matrix(2, 3);
        matrix2x3.set(0, 0, 1.0); matrix2x3.set(0, 1, 2.0); matrix2x3.set(0, 2, 3.0);
        matrix2x3.set(1, 0, 4.0); matrix2x3.set(1, 1, 5.0); matrix2x3.set(1, 2, 6.0);
        
        matrix3x2 = new Matrix(3, 2);
        matrix3x2.set(0, 0, 1.0); matrix3x2.set(0, 1, 2.0);
        matrix3x2.set(1, 0, 3.0); matrix3x2.set(1, 1, 4.0);
        matrix3x2.set(2, 0, 5.0); matrix3x2.set(2, 1, 6.0);
        
        matrix3x3 = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                matrix3x3.set(i, j, i * 3 + j + 1);
            }
        }
        
        matrix1x1 = new Matrix(1, 1);
        matrix1x1.set(0, 0, 42.0);
        
        vector3x1 = new Matrix(3, 1);
        vector3x1.set(0, 0, 1.0); vector3x1.set(1, 0, 2.0); vector3x1.set(2, 0, 3.0);
        
        rowVector1x3 = new Matrix(1, 3);
        rowVector1x3.set(0, 0, 1.0); rowVector1x3.set(0, 1, 2.0); rowVector1x3.set(0, 2, 3.0);
        
        identityMatrix = Matrix.eye(3);
        zeroMatrix = Matrix.zeros(3, 3);
        
        denseMatrix = new DenseMatrixImpl(3, 4);
    }
    
    // ========== Matrix Core Tests ==========
    
    @Test
    void testConstructorDimensions() {
        assertEquals(2, matrix2x2.getNumRows());
        assertEquals(2, matrix2x2.getNumCols());
        assertEquals(2, matrix2x3.getNumRows());
        assertEquals(3, matrix2x3.getNumCols());
        assertEquals(3, matrix3x2.getNumRows());
        assertEquals(2, matrix3x2.getNumCols());
        assertEquals(1, matrix1x1.getNumRows());
        assertEquals(1, matrix1x1.getNumCols());
    }
    
    @Test
    void testConstructorWithCapacity() {
        Matrix matrix = new Matrix(2, 3, 10);
        assertEquals(2, matrix.getNumRows());
        assertEquals(3, matrix.getNumCols());
    }
    
    @Test
    void testConstructorFromDoubleArray() {
        double[][] data = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        Matrix matrix = new Matrix(data);
        
        assertEquals(3, matrix.getNumRows());
        assertEquals(2, matrix.getNumCols());
        assertEquals(1.0, matrix.get(0, 0), ZERO_TOL);
        assertEquals(2.0, matrix.get(0, 1), ZERO_TOL);
        assertEquals(3.0, matrix.get(1, 0), ZERO_TOL);
        assertEquals(4.0, matrix.get(1, 1), ZERO_TOL);
        assertEquals(5.0, matrix.get(2, 0), ZERO_TOL);
        assertEquals(6.0, matrix.get(2, 1), ZERO_TOL);
    }
    
    @Test
    void testConstructorFromIntArray() {
        int[] data = {1, 2, 3};
        Matrix matrix = new Matrix(data);
        
        assertEquals(3, matrix.getNumRows());
        assertEquals(1, matrix.getNumCols());
        assertEquals(1.0, matrix.get(0, 0), ZERO_TOL);
        assertEquals(2.0, matrix.get(1, 0), ZERO_TOL);
        assertEquals(3.0, matrix.get(2, 0), ZERO_TOL);
    }
    
    @Test
    void testConstructorFromDoubleArrayColumn() {
        double[] data = {1.5, 2.5, 3.5};
        Matrix matrix = new Matrix(data);
        
        assertEquals(3, matrix.getNumRows());
        assertEquals(1, matrix.getNumCols());
        assertEquals(1.5, matrix.get(0, 0), ZERO_TOL);
        assertEquals(2.5, matrix.get(1, 0), ZERO_TOL);
        assertEquals(3.5, matrix.get(2, 0), ZERO_TOL);
    }
    
    @Test
    void testCopyConstructor() {
        matrix2x3.set(0, 0, 1.0);
        matrix2x3.set(1, 2, 2.0);
        
        Matrix copy = new Matrix(matrix2x3);
        
        assertEquals(matrix2x3.getNumRows(), copy.getNumRows());
        assertEquals(matrix2x3.getNumCols(), copy.getNumCols());
        assertEquals(1.0, copy.get(0, 0), ZERO_TOL);
        assertEquals(2.0, copy.get(1, 2), ZERO_TOL);
        assertEquals(5.0, copy.get(1, 1), ZERO_TOL);
        
        // Verify independence (deep copy)
        matrix2x3.set(0, 0, 999.0);
        assertEquals(1.0, copy.get(0, 0), ZERO_TOL);
    }
    
    @Test
    void testInitialValuesZero() {
        Matrix m = new Matrix(2, 3);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, m.get(i, j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testSetAndGet() {
        Matrix m = new Matrix(2, 3);
        m.set(0, 0, 1.5);
        assertEquals(1.5, m.get(0, 0), ZERO_TOL);
        
        m.set(1, 2, -3.7);
        assertEquals(-3.7, m.get(1, 2), ZERO_TOL);
        
        // Verify other elements remain zero
        assertEquals(0.0, m.get(0, 1), ZERO_TOL);
        assertEquals(0.0, m.get(1, 0), ZERO_TOL);
    }
    
    @Test
    void testInvalidAccess() {
        // Test out-of-bounds access
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.get(-1, 0));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.get(0, -1));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.get(2, 0));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.get(0, 3));
        
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.set(-1, 0, 1.0));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.set(0, -1, 1.0));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.set(2, 0, 1.0));
        assertThrows(IllegalArgumentException.class, () -> matrix2x3.set(0, 3, 1.0));
    }
    
    @Test
    void testSpecialValues() {
        Matrix m = new Matrix(2, 3);
        // Test infinity and NaN values
        m.set(0, 0, Double.POSITIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, m.get(0, 0), ZERO_TOL);
        
        m.set(0, 1, Double.NEGATIVE_INFINITY);
        assertEquals(Double.NEGATIVE_INFINITY, m.get(0, 1), ZERO_TOL);
        
        m.set(0, 2, Double.NaN);
        assertTrue(Double.isNaN(m.get(0, 2)));
    }
    
    // ========== DenseMatrix Implementation Tests ==========
    
    @Test
    void testDenseMatrixBasicCreation() {
        assertEquals(3, denseMatrix.getNumRows());
        assertEquals(4, denseMatrix.getNumCols());
        
        // All elements should be initialized to zero
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(0.0, denseMatrix.get(i, j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testDenseMatrixSetAndGet() {
        denseMatrix.set(1, 2, 5.5);
        assertEquals(5.5, denseMatrix.get(1, 2), ZERO_TOL);
        
        denseMatrix.set(0, 0, -2.3);
        assertEquals(-2.3, denseMatrix.get(0, 0), ZERO_TOL);
        
        // Other elements should remain zero
        assertEquals(0.0, denseMatrix.get(2, 3), ZERO_TOL);
    }
    
    @Test
    void testDenseMatrixConstructorWithEJMLMatrix() {
        DMatrixRMaj ejmlMatrix = new DMatrixRMaj(2, 2);
        ejmlMatrix.set(0, 0, 1.0);
        ejmlMatrix.set(0, 1, 2.0);
        ejmlMatrix.set(1, 0, 3.0);
        ejmlMatrix.set(1, 1, 4.0);
        
        DenseMatrixImpl matrix3 = new DenseMatrixImpl(ejmlMatrix);
        assertEquals(2, matrix3.getNumRows());
        assertEquals(2, matrix3.getNumCols());
        assertEquals(1.0, matrix3.get(0, 0), ZERO_TOL);
        assertEquals(2.0, matrix3.get(0, 1), ZERO_TOL);
        assertEquals(3.0, matrix3.get(1, 0), ZERO_TOL);
        assertEquals(4.0, matrix3.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testDenseMatrixCopyConstructor() {
        denseMatrix.set(0, 0, 1.5);
        denseMatrix.set(2, 3, -0.7);
        
        DenseMatrixImpl copy = new DenseMatrixImpl(denseMatrix);
        
        assertEquals(denseMatrix.getNumRows(), copy.getNumRows());
        assertEquals(denseMatrix.getNumCols(), copy.getNumCols());
        assertEquals(1.5, copy.get(0, 0), ZERO_TOL);
        assertEquals(-0.7, copy.get(2, 3), ZERO_TOL);
        
        // Verify it's a deep copy
        denseMatrix.set(0, 0, 999.0);
        assertEquals(1.5, copy.get(0, 0), ZERO_TOL); // Should not change
    }
    
    @Test
    void testDenseMatrixClone() {
        denseMatrix.set(1, 1, 42.0);
        
        BaseMatrix cloned = denseMatrix.copy();
        
        assertTrue(cloned instanceof DenseMatrixImpl);
        assertEquals(denseMatrix.getNumRows(), cloned.getNumRows());
        assertEquals(denseMatrix.getNumCols(), cloned.getNumCols());
        assertEquals(42.0, cloned.get(1, 1), ZERO_TOL);
        
        // Verify it's a deep copy
        denseMatrix.set(1, 1, 123.0);
        assertEquals(42.0, cloned.get(1, 1), ZERO_TOL); // Should not change
    }
    
    @Test
    void testDenseMatrixEquals() {
        denseMatrix.set(0, 0, 1.0);
        denseMatrix.set(1, 2, 3.5);
        
        DenseMatrixImpl other = new DenseMatrixImpl(3, 4);
        other.set(0, 0, 1.0);
        other.set(1, 2, 3.5);
        
        assertTrue(denseMatrix.equals(other));
        assertTrue(other.equals(denseMatrix));
        
        // Test inequality
        other.set(2, 1, 0.1);
        assertFalse(denseMatrix.equals(other));
        
        // Test different dimensions
        DenseMatrixImpl different = new DenseMatrixImpl(2, 4);
        assertFalse(denseMatrix.equals(different));
    }
    
    // ========== Matrix Arithmetic Tests ==========
    
    @Test
    void testAddMatrix() {
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 1.0); a.set(0, 1, 2.0);
        a.set(1, 0, 3.0); a.set(1, 1, 4.0);
        
        Matrix b = new Matrix(2, 2);
        b.set(0, 0, 0.5); b.set(0, 1, 1.5);
        b.set(1, 0, 2.5); b.set(1, 1, 3.5);
        
        Matrix result = a.add(b);
        
        assertEquals(1.5, result.get(0, 0), ZERO_TOL);
        assertEquals(3.5, result.get(0, 1), ZERO_TOL);
        assertEquals(5.5, result.get(1, 0), ZERO_TOL);
        assertEquals(7.5, result.get(1, 1), ZERO_TOL);
        
        // Test commutativity
        Matrix result2 = b.add(a);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(result.get(i, j), result2.get(i, j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testAddMatrixInPlace() {
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 1.0); a.set(0, 1, 2.0);
        a.set(1, 0, 3.0); a.set(1, 1, 4.0);
        
        Matrix b = new Matrix(2, 2);
        b.set(0, 0, 0.5); b.set(0, 1, 1.5);
        b.set(1, 0, 2.5); b.set(1, 1, 3.5);
        
        a.addEq(b); // In-place addition
        
        assertEquals(1.5, a.get(0, 0), ZERO_TOL);
        assertEquals(3.5, a.get(0, 1), ZERO_TOL);
        assertEquals(5.5, a.get(1, 0), ZERO_TOL);
        assertEquals(7.5, a.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testAddScalar() {
        Matrix result = matrix2x2.add(10.0);
        
        assertEquals(11.0, result.get(0, 0), ZERO_TOL);
        assertEquals(12.0, result.get(0, 1), ZERO_TOL);
        assertEquals(13.0, result.get(1, 0), ZERO_TOL);
        assertEquals(14.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testSubtractMatrix() {
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 5.0); a.set(0, 1, 6.0);
        a.set(1, 0, 7.0); a.set(1, 1, 8.0);
        
        Matrix b = new Matrix(2, 2);
        b.set(0, 0, 1.0); b.set(0, 1, 2.0);
        b.set(1, 0, 3.0); b.set(1, 1, 4.0);
        
        Matrix result = a.sub(b);
        
        assertEquals(4.0, result.get(0, 0), ZERO_TOL);
        assertEquals(4.0, result.get(0, 1), ZERO_TOL);
        assertEquals(4.0, result.get(1, 0), ZERO_TOL);
        assertEquals(4.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testSubtractScalar() {
        Matrix result = matrix2x2.sub(0.5);
        
        assertEquals(0.5, result.get(0, 0), ZERO_TOL);
        assertEquals(1.5, result.get(0, 1), ZERO_TOL);
        assertEquals(2.5, result.get(1, 0), ZERO_TOL);
        assertEquals(3.5, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplyCompatibleMatrices() {
        // 2x3 * 3x2 = 2x2
        Matrix result = matrix2x3.mult(matrix3x2);
        
        assertEquals(2, result.getNumRows());
        assertEquals(2, result.getNumCols());
        
        // [1 2 3] * [1 2]   = [22 28]
        // [4 5 6]   [3 4]     [49 64]
        //           [5 6]
        assertEquals(22.0, result.get(0, 0), ZERO_TOL);
        assertEquals(28.0, result.get(0, 1), ZERO_TOL);
        assertEquals(49.0, result.get(1, 0), ZERO_TOL);
        assertEquals(64.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplySquareMatrix() {
        Matrix result = matrix2x2.mult(matrix2x2);
        
        // [1 2] * [1 2] = [7  10]
        // [3 4]   [3 4]   [15 22]
        assertEquals(7.0, result.get(0, 0), ZERO_TOL);
        assertEquals(10.0, result.get(0, 1), ZERO_TOL);
        assertEquals(15.0, result.get(1, 0), ZERO_TOL);
        assertEquals(22.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplyByIdentity() {
        Matrix identity = Matrix.eye(2);
        Matrix result = matrix2x2.mult(identity);
        
        assertEquals(matrix2x2.get(0, 0), result.get(0, 0), ZERO_TOL);
        assertEquals(matrix2x2.get(0, 1), result.get(0, 1), ZERO_TOL);
        assertEquals(matrix2x2.get(1, 0), result.get(1, 0), ZERO_TOL);
        assertEquals(matrix2x2.get(1, 1), result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplyByZero() {
        Matrix zero = new Matrix(2, 2);
        Matrix result = matrix2x2.mult(zero);
        
        assertEquals(0.0, result.get(0, 0), ZERO_TOL);
        assertEquals(0.0, result.get(0, 1), ZERO_TOL);
        assertEquals(0.0, result.get(1, 0), ZERO_TOL);
        assertEquals(0.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplyVectors() {
        // Row vector * column vector = scalar (1x1 matrix)
        Matrix result = rowVector1x3.mult(vector3x1);
        assertEquals(1, result.getNumRows());
        assertEquals(1, result.getNumCols());
        assertEquals(14.0, result.get(0, 0), ZERO_TOL); // 1*1 + 2*2 + 3*3 = 14
        
        // Column vector * row vector = outer product (3x3 matrix)
        Matrix outer = vector3x1.mult(rowVector1x3);
        assertEquals(3, outer.getNumRows());
        assertEquals(3, outer.getNumCols());
        assertEquals(1.0, outer.get(0, 0), ZERO_TOL); // 1*1
        assertEquals(2.0, outer.get(0, 1), ZERO_TOL); // 1*2
        assertEquals(6.0, outer.get(2, 1), ZERO_TOL); // 3*2
    }
    
    @Test
    void testMultiplyScalar() {
        Matrix result = matrix2x2.scale(2.5);
        
        assertEquals(2.5, result.get(0, 0), ZERO_TOL);
        assertEquals(5.0, result.get(0, 1), ZERO_TOL);
        assertEquals(7.5, result.get(1, 0), ZERO_TOL);
        assertEquals(10.0, result.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testMultiplyScalarInPlace() {
        matrix2x2.scaleEq(3.0);
        
        assertEquals(3.0, matrix2x2.get(0, 0), ZERO_TOL);
        assertEquals(6.0, matrix2x2.get(0, 1), ZERO_TOL);
        assertEquals(9.0, matrix2x2.get(1, 0), ZERO_TOL);
        assertEquals(12.0, matrix2x2.get(1, 1), ZERO_TOL);
    }
    
    @Test
    void testElementMult() {
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 2.0); a.set(0, 1, 3.0);
        a.set(1, 0, 4.0); a.set(1, 1, 5.0);
        
        Matrix original = new Matrix(2, 2);
        original.set(0, 0, 1.0); original.set(0, 1, 2.0);
        original.set(1, 0, 3.0); original.set(1, 1, 4.0);
        
        Matrix result = original.elementMult(a);
        
        assertEquals(2.0, result.get(0, 0), ZERO_TOL);  // 1*2
        assertEquals(6.0, result.get(0, 1), ZERO_TOL);  // 2*3
        assertEquals(12.0, result.get(1, 0), ZERO_TOL); // 3*4
        assertEquals(20.0, result.get(1, 1), ZERO_TOL); // 4*5
    }
    
    // ========== Matrix Factory Methods Tests ==========
    
    @Test
    void testEyeSquare() {
        Matrix eye3 = Matrix.eye(3);
        
        assertEquals(3, eye3.getNumRows());
        assertEquals(3, eye3.getNumCols());
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == j) {
                    assertEquals(1.0, eye3.get(i, j), ZERO_TOL);
                } else {
                    assertEquals(0.0, eye3.get(i, j), ZERO_TOL);
                }
            }
        }
    }
    
    @Test
    void testOnes() {
        Matrix ones23 = Matrix.ones(2, 3);
        
        assertEquals(2, ones23.getNumRows());
        assertEquals(3, ones23.getNumCols());
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(1.0, ones23.get(i, j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testZeros() {
        Matrix zeros23 = Matrix.zeros(2, 3);
        
        assertEquals(2, zeros23.getNumRows());
        assertEquals(3, zeros23.getNumCols());
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, zeros23.get(i, j), ZERO_TOL);
            }
        }
    }
    
    // ========== Matrix Manipulation Tests ==========
    
    @Test
    void testTranspose() {
        Matrix transposed = matrix2x3.transpose();
        
        assertEquals(3, transposed.getNumRows());
        assertEquals(2, transposed.getNumCols());
        
        assertEquals(1.0, transposed.get(0, 0), ZERO_TOL); // was (0,0)
        assertEquals(4.0, transposed.get(0, 1), ZERO_TOL); // was (1,0)
        assertEquals(2.0, transposed.get(1, 0), ZERO_TOL); // was (0,1)
        assertEquals(5.0, transposed.get(1, 1), ZERO_TOL); // was (1,1)
        assertEquals(3.0, transposed.get(2, 0), ZERO_TOL); // was (0,2)
        assertEquals(6.0, transposed.get(2, 1), ZERO_TOL); // was (1,2)
    }
    
    @Test
    void testDoubleTranspose() {
        Matrix doubleTransposed = matrix2x3.transpose().transpose();
        
        assertEquals(matrix2x3.getNumRows(), doubleTransposed.getNumRows());
        assertEquals(matrix2x3.getNumCols(), doubleTransposed.getNumCols());
        
        for (int i = 0; i < matrix2x3.getNumRows(); i++) {
            for (int j = 0; j < matrix2x3.getNumCols(); j++) {
                assertEquals(matrix2x3.get(i, j), doubleTransposed.get(i, j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testGetColumn() {
        Matrix col0 = matrix2x3.getColumn(0);
        assertEquals(2, col0.getNumRows());
        assertEquals(1, col0.getNumCols());
        assertEquals(1.0, col0.get(0, 0), ZERO_TOL);
        assertEquals(4.0, col0.get(1, 0), ZERO_TOL);
        
        Matrix col2 = matrix2x3.getColumn(2);
        assertEquals(2, col2.getNumRows());
        assertEquals(1, col2.getNumCols());
        assertEquals(3.0, col2.get(0, 0), ZERO_TOL);
        assertEquals(6.0, col2.get(1, 0), ZERO_TOL);
    }
    
    @Test
    void testGetRow() {
        Matrix row0 = matrix2x3.getRow(0);
        assertEquals(1, row0.getNumRows());
        assertEquals(3, row0.getNumCols());
        assertEquals(1.0, row0.get(0, 0), ZERO_TOL);
        assertEquals(2.0, row0.get(0, 1), ZERO_TOL);
        assertEquals(3.0, row0.get(0, 2), ZERO_TOL);
    }
    
    @Test
    void testConcatColumns() {
        Matrix left = new Matrix(2, 2);
        left.set(0, 0, 1.0); left.set(0, 1, 2.0);
        left.set(1, 0, 3.0); left.set(1, 1, 4.0);
        
        Matrix right = new Matrix(2, 1);
        right.set(0, 0, 5.0);
        right.set(1, 0, 6.0);
        
        Matrix result = new Matrix(left.getNumRows(), left.getNumCols() + right.getNumCols());
        Matrix.concatColumns(left, right, result);
        
        assertEquals(2, result.getNumRows());
        assertEquals(3, result.getNumCols());
        assertEquals(1.0, result.get(0, 0), ZERO_TOL);
        assertEquals(2.0, result.get(0, 1), ZERO_TOL);
        assertEquals(5.0, result.get(0, 2), ZERO_TOL);
        assertEquals(3.0, result.get(1, 0), ZERO_TOL);
        assertEquals(4.0, result.get(1, 1), ZERO_TOL);
        assertEquals(6.0, result.get(1, 2), ZERO_TOL);
    }
    
    @Test
    void testConcatRows() {
        Matrix top = new Matrix(1, 3);
        top.set(0, 0, 1.0); top.set(0, 1, 2.0); top.set(0, 2, 3.0);
        
        Matrix bottom = new Matrix(1, 3);
        bottom.set(0, 0, 4.0); bottom.set(0, 1, 5.0); bottom.set(0, 2, 6.0);
        
        Matrix result = new Matrix(top.getNumRows() + bottom.getNumRows(), top.getNumCols());
        Matrix.concatRows(top, bottom, result);
        
        assertEquals(2, result.getNumRows());
        assertEquals(3, result.getNumCols());
        assertEquals(1.0, result.get(0, 0), ZERO_TOL);
        assertEquals(2.0, result.get(0, 1), ZERO_TOL);
        assertEquals(3.0, result.get(0, 2), ZERO_TOL);
        assertEquals(4.0, result.get(1, 0), ZERO_TOL);
        assertEquals(5.0, result.get(1, 1), ZERO_TOL);
        assertEquals(6.0, result.get(1, 2), ZERO_TOL);
    }
    
    // ========== Matrix View Tests ==========
    
    @Test
    void testColumnViewBasicAccess() {
        ColumnView col1 = matrix2x3.getColumnView(1);
        
        assertEquals(2, col1.getNumRows());
        assertEquals(1, col1.getColumnIndex());
        assertEquals(2.0, col1.get(0), ZERO_TOL); // matrix2x3(0,1)
        assertEquals(5.0, col1.get(1), ZERO_TOL); // matrix2x3(1,1)
    }
    
    @Test
    void testRowViewBasicAccess() {
        RowView row1 = matrix2x3.getRowView(1);
        
        assertEquals(3, row1.getNumCols());
        assertEquals(1, row1.getRowIndex());
        assertEquals(4.0, row1.get(0), ZERO_TOL); // matrix2x3(1,0)
        assertEquals(5.0, row1.get(1), ZERO_TOL); // matrix2x3(1,1)
        assertEquals(6.0, row1.get(2), ZERO_TOL); // matrix2x3(1,2)
    }
    
    @Test
    void testColumnViewEquivalence() {
        for (int j = 0; j < matrix2x3.getNumCols(); j++) {
            ColumnView colView = matrix2x3.getColumnView(j);
            Matrix colMatrix = matrix2x3.getColumn(j);
            
            assertEquals(colMatrix.getNumRows(), colView.getNumRows());
            for (int i = 0; i < colView.getNumRows(); i++) {
                assertEquals(colMatrix.get(i, 0), colView.get(i), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testRowViewEquivalence() {
        for (int i = 0; i < matrix2x3.getNumRows(); i++) {
            RowView rowView = matrix2x3.getRowView(i);
            Matrix rowMatrix = matrix2x3.getRow(i);
            
            assertEquals(rowMatrix.getNumCols(), rowView.getNumCols());
            for (int j = 0; j < rowView.getNumCols(); j++) {
                assertEquals(rowMatrix.get(0, j), rowView.get(j), ZERO_TOL);
            }
        }
    }
    
    @Test
    void testMultColumnView() {
        Matrix rowVec = new Matrix(1, 2);
        rowVec.set(0, 0, 0.1); rowVec.set(0, 1, 0.2);
        
        ColumnView colView = matrix2x3.getColumnView(0);
        double result = rowVec.multColumnView(colView);
        
        // Verify against traditional multiplication
        Matrix col = matrix2x3.getColumn(0);
        double expected = rowVec.mult(col).get(0, 0);
        assertEquals(expected, result, ZERO_TOL);
    }
    
    @Test
    void testRowViewDotProduct() {
        RowView row0 = matrix2x3.getRowView(0); // [1, 2, 3]
        
        Matrix colVec = new Matrix(3, 1);
        colVec.set(0, 0, 0.1); colVec.set(1, 0, 0.2); colVec.set(2, 0, 0.3);
        
        double result = row0.dotProduct(colVec);
        // Expected: 1*0.1 + 2*0.2 + 3*0.3 = 0.1 + 0.4 + 0.9 = 1.4
        assertEquals(1.4, result, ZERO_TOL);
    }
    
    // ========== Matrix Properties Tests ==========
    
    @Test
    void testIsFinite() {
        assertTrue(matrix2x2.isFinite());
        assertTrue(Matrix.zeros(3, 3).isFinite());
        
        Matrix infinite = new Matrix(2, 2);
        infinite.set(0, 0, 1.0);
        infinite.set(0, 1, Double.POSITIVE_INFINITY);
        assertFalse(infinite.isFinite());
        
        Matrix nan = new Matrix(2, 2);
        nan.set(0, 0, 1.0);
        nan.set(0, 1, Double.NaN);
        assertFalse(nan.isFinite());
    }
    
    @Test
    void testIsDiagonal() {
        Matrix diagonal = Matrix.eye(3);
        diagonal.set(0, 0, 2.0);
        diagonal.set(1, 1, 3.0);
        diagonal.set(2, 2, 4.0);
        assertTrue(diagonal.isDiag());
        
        Matrix nonDiagonal = new Matrix(diagonal);
        nonDiagonal.set(0, 1, 0.1);
        assertFalse(nonDiagonal.isDiag());
        
        Matrix zerosMatrix = Matrix.zeros(3, 3);
        assertFalse(zerosMatrix.isDiag()); // Zero matrix is not considered diagonal in this implementation
    }
    
    @Test
    void testElementSum() {
        assertEquals(10.0, matrix2x2.elementSum(), ZERO_TOL); // 1+2+3+4 = 10
        assertEquals(6.0, Matrix.ones(2, 3).elementSum(), ZERO_TOL);
        assertEquals(0.0, Matrix.zeros(5, 5).elementSum(), ZERO_TOL);
    }
    
    @Test
    void testNorm() {
        // Frobenius norm should be sqrt of sum of squares
        assertEquals(Math.sqrt(30.0), matrix2x2.norm(), ZERO_TOL); // sqrt(1+4+9+16) = sqrt(30)
        assertEquals(Math.sqrt(6.0), Matrix.ones(2, 3).norm(), ZERO_TOL);
        assertEquals(0.0, Matrix.zeros(3, 3).norm(), ZERO_TOL);
    }
    
    @Test
    void testDeterminant() {
        // 2x2: det([1,2; 3,4]) = 1*4 - 2*3 = -2
        assertEquals(-2.0, matrix2x2.det(), ZERO_TOL);
        
        assertEquals(1.0, Matrix.eye(3).det(), ZERO_TOL);
        assertEquals(0.0, Matrix.zeros(3, 3).det(), ZERO_TOL);
    }
    
    // ========== Linear System Solving Tests ==========
    
    @Test
    void testSolveLinearSystem() {
        // Solve Ax = b where A = matrix2x2, b = [5; 11]
        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 1.0); A.set(0, 1, 2.0);
        A.set(1, 0, 3.0); A.set(1, 1, 4.0);
        
        Matrix b = new Matrix(2, 1);
        b.set(0, 0, 5.0);
        b.set(1, 0, 11.0);
        
        Matrix x = new Matrix(2, 1);
        Matrix.solve(A, b, x);
        
        assertEquals(2, x.getNumRows());
        assertEquals(1, x.getNumCols());
        
        // Verify solution by checking Ax = b
        Matrix product = A.mult(x);
        assertEquals(5.0, product.get(0, 0), 1e-12);
        assertEquals(11.0, product.get(1, 0), 1e-12);
    }
    
    @Test
    void testInverse() {
        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 2.0); A.set(0, 1, 1.0);
        A.set(1, 0, 1.0); A.set(1, 1, 3.0);
        
        Matrix inverse = A.inv();
        
        assertEquals(2, inverse.getNumRows());
        assertEquals(2, inverse.getNumCols());
        
        // Verify A * A^(-1) = I
        Matrix product = A.mult(inverse);
        assertEquals(1.0, product.get(0, 0), 1e-12);
        assertEquals(0.0, product.get(0, 1), 1e-12);
        assertEquals(0.0, product.get(1, 0), 1e-12);
        assertEquals(1.0, product.get(1, 1), 1e-12);
    }
    
    // ========== Matrix Choice Tests ==========
    
    @Test
    void shouldChooseNothing() {
        Matrix A = new Matrix(1, 4);
        A.set(0, 0, 1);
        A.set(0, 1, 2);
        A.set(0, 2, 3);
        A.set(0, 3, 4);
        Matrix res = jline.util.Maths.nCk(A, 0);
        assertTrue(res.isEmpty());
    }

    @Test
    void shouldChooseOne() {
        Matrix A = new Matrix(1, 4);
        A.set(0, 0, 1);
        A.set(0, 1, 2);
        A.set(0, 2, 3);
        A.set(0, 3, 4);
        Matrix res = jline.util.Maths.nCk(A, 1);
        assertEquals(4, res.getNumRows());
        assertEquals(1, res.getNumCols());
        assertEquals(1, res.value());
        assertEquals(2, res.get(1, 0));
        assertEquals(3, res.get(2, 0));
        assertEquals(4, res.get(3, 0));
    }

    @Test
    void shouldChooseTwo() {
        Matrix A = new Matrix(1, 4);
        A.set(0, 0, 1);
        A.set(0, 1, 2);
        A.set(0, 2, 3);
        A.set(0, 3, 4);
        Matrix res = jline.util.Maths.nCk(A, 2);
        assertEquals(6, res.getNumRows());
        assertEquals(2, res.getNumCols());
        assertEquals(1, res.value());
    }
    
}