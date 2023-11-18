package jline.api;

import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CACHETest {

    private final double tolerance = 1.0e-15;

    @Test
    void cache_gamma_lp() {
        Matrix[] lambda = new Matrix[3];
        lambda[0] = new Matrix(5, 2);
        lambda[0].set(0, 0, 0.4);
        lambda[0].set(1, 0, 0.4);
        lambda[0].set(2, 0, 0.4);
        lambda[0].set(3, 0, 0.4);
        lambda[0].set(4, 0, 0.4);
        lambda[0].set(0, 1, 0.4);
        lambda[0].set(1, 1, 0.4);
        lambda[0].set(2, 1, 0.4);
        lambda[0].set(3, 1, 0.4);
        lambda[0].set(4, 1, 0.4);
        lambda[1] = new Matrix(5, 2);
        lambda[2] = new Matrix(5, 2);
        Matrix[][] R = new Matrix[3][5];
        Matrix x = new Matrix(2,2);
        x.set(0,1,1);
        x.set(1,1,1);
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 5; j++){
                R[i][j] = new Matrix(x);
            }
        }
        CACHE.cacheGammaLpReturn ret = CACHE.cache_gamma_lp(lambda, R);
        assertEquals(0.4, ret.gamma.get(0,0), tolerance);
        assertEquals(0.4, ret.gamma.get(1,0), tolerance);
        assertEquals(0.4, ret.gamma.get(2,0), tolerance);
        assertEquals(0.4, ret.gamma.get(3,0), tolerance);
        assertEquals(0.4, ret.gamma.get(4,0), tolerance);
        assertEquals(3, ret.u);
        assertEquals(5, ret.n);
        assertEquals(1, ret.h);
    }

    @Test
    void cache_mva_test1(){
        Matrix gamma = new Matrix(5, 1);
        gamma.fill(0.4);
        Matrix m = new Matrix(1,1);
        m.set(0,0,2);
        CACHE.cacheMVAReturn ret = CACHE.cache_mva(gamma, m);
        assertEquals(0.4, ret.pi.get(0,0), tolerance);
        assertEquals(0.4, ret.pi.get(1,0), tolerance);
        assertEquals(0.4, ret.pi.get(2,0), tolerance);
        assertEquals(0.4, ret.pi.get(3,0), tolerance);
        assertEquals(0.4, ret.pi.get(4,0), tolerance);

        assertEquals(0.6, ret.pi0.get(0,0), tolerance);
        assertEquals(0.6, ret.pi0.get(1,0), tolerance);
        assertEquals(0.6, ret.pi0.get(2,0), tolerance);
        assertEquals(0.6, ret.pi0.get(3,0), tolerance);
        assertEquals(0.6, ret.pi0.get(4,0), tolerance);

        assertEquals(0.4, ret.pij.get(0,0), tolerance);
        assertEquals(0.4, ret.pij.get(1,0), tolerance);
        assertEquals(0.4, ret.pij.get(2,0), tolerance);
        assertEquals(0.4, ret.pij.get(3,0), tolerance);
        assertEquals(0.4, ret.pij.get(4,0), tolerance);

        assertEquals(1.25, ret.x.get(0,0), tolerance);

        assertEquals(0.5, ret.u.get(0,0), tolerance);
        assertEquals(0.5, ret.u.get(1,0), tolerance);
        assertEquals(0.5, ret.u.get(2,0), tolerance);
        assertEquals(0.5, ret.u.get(3,0), tolerance);
        assertEquals(0.5, ret.u.get(4,0), tolerance);

        assertEquals(1, ret.E);
    }

    @Test
    void cache_mva_test2(){
        Matrix gamma = new Matrix(5, 2);
        gamma.set(0,0,0.4);
        gamma.set(1,0,0.4);
        gamma.set(2,0,0.4);
        gamma.set(3,0,0.4);
        gamma.set(4,0,0.4);
        gamma.set(0,1,0.16);
        gamma.set(1,1,0.16);
        gamma.set(2,1,0.16);
        gamma.set(3,1,0.16);
        gamma.set(4,1,0.16);

        Matrix m = new Matrix(1,2);
        m.set(0,0,1);
        m.set(0,1,2);
        CACHE.cacheMVAReturn ret = CACHE.cache_mva(gamma, m);
        assertEquals(0.6, ret.pi.get(0,0), tolerance);
        assertEquals(0.6, ret.pi.get(1,0), tolerance);
        assertEquals(0.6, ret.pi.get(2,0), tolerance);
        assertEquals(0.6, ret.pi.get(3,0), tolerance);
        assertEquals(0.6, ret.pi.get(4,0), tolerance);

        assertEquals(0.4, ret.pi0.get(0,0), tolerance);
        assertEquals(0.4, ret.pi0.get(1,0), tolerance);
        assertEquals(0.4, ret.pi0.get(2,0), tolerance);
        assertEquals(0.4, ret.pi0.get(3,0), tolerance);
        assertEquals(0.4, ret.pi0.get(4,0), tolerance);

        assertEquals(0.2, ret.pij.get(0,0), tolerance);
        assertEquals(0.2, ret.pij.get(1,0), tolerance);
        assertEquals(0.2, ret.pij.get(2,0), tolerance);
        assertEquals(0.2, ret.pij.get(3,0), tolerance);
        assertEquals(0.2, ret.pij.get(4,0), tolerance);
        assertEquals(0.4, ret.pij.get(0,1), tolerance);
        assertEquals(0.4, ret.pij.get(1,1), tolerance);
        assertEquals(0.4, ret.pij.get(2,1), tolerance);
        assertEquals(0.4, ret.pij.get(3,1), tolerance);
        assertEquals(0.4, ret.pij.get(4,1), tolerance);

        assertEquals(0.833333333333333, ret.x.get(0,0), tolerance);
        assertEquals(4.166666666666666, ret.x.get(0,1), tolerance);

        assertEquals(0.333333333333333, ret.u.get(0,0), tolerance);
        assertEquals(0.333333333333333, ret.u.get(1,0), tolerance);
        assertEquals(0.333333333333333, ret.u.get(2,0), tolerance);
        assertEquals(0.333333333333333, ret.u.get(3,0), tolerance);
        assertEquals(0.333333333333333, ret.u.get(4,0), tolerance);
        assertEquals(0.666666666666667, ret.u.get(0,1), tolerance);
        assertEquals(0.666666666666667, ret.u.get(1,1), tolerance);
        assertEquals(0.666666666666667, ret.u.get(2,1), tolerance);
        assertEquals(0.666666666666667, ret.u.get(3,1), tolerance);
        assertEquals(0.666666666666667, ret.u.get(4,1), tolerance);

        assertEquals(1, ret.E);
    }

    @Test
    void cache_prob_asy_test1(){
        Matrix gamma = new Matrix(5, 1);
        gamma.fill(0.4);
        Matrix m = new Matrix(1,1);
        m.set(0,0,2);
        Matrix pij = CACHE.cache_prob_asy(gamma, m);
        assertEquals(0.599999999999999, pij.get(0,0), tolerance);
        assertEquals(0.599999999999999, pij.get(1,0), tolerance);
        assertEquals(0.599999999999999, pij.get(2,0), tolerance);
        assertEquals(0.599999999999999, pij.get(3,0), tolerance);
        assertEquals(0.599999999999999, pij.get(4,0), tolerance);
        assertEquals(0.400000000000001, pij.get(0,1), tolerance);
        assertEquals(0.400000000000001, pij.get(1,1), tolerance);
        assertEquals(0.400000000000001, pij.get(2,1), tolerance);
        assertEquals(0.400000000000001, pij.get(3,1), tolerance);
        assertEquals(0.400000000000001, pij.get(4,1), tolerance);
    }

    @Test
    void cache_prob_asy_test2(){
        Matrix gamma = new Matrix(5, 2);
        gamma.set(0,0,0.4);
        gamma.set(1,0,0.4);
        gamma.set(2,0,0.4);
        gamma.set(3,0,0.4);
        gamma.set(4,0,0.4);
        gamma.set(0,1,0.16);
        gamma.set(1,1,0.16);
        gamma.set(2,1,0.16);
        gamma.set(3,1,0.16);
        gamma.set(4,1,0.16);

        Matrix m = new Matrix(1,2);
        m.set(0,0,1);
        m.set(0,1,2);
        Matrix pij = CACHE.cache_prob_asy(gamma, m);
        assertEquals(0.399999999999999, pij.get(0,0), tolerance);
        assertEquals(0.399999999999999, pij.get(1,0), tolerance);
        assertEquals(0.399999999999999, pij.get(2,0), tolerance);
        assertEquals(0.399999999999999, pij.get(3,0), tolerance);
        assertEquals(0.399999999999999, pij.get(4,0), tolerance);
        assertEquals(0.600000000000001, pij.get(0,1), tolerance);
        assertEquals(0.600000000000001, pij.get(1,1), tolerance);
        assertEquals(0.600000000000001, pij.get(2,1), tolerance);
        assertEquals(0.600000000000001, pij.get(3,1), tolerance);
        assertEquals(0.600000000000001, pij.get(4,1), tolerance);
        assertEquals(0.600000000000001, pij.get(0,2), tolerance);
        assertEquals(0.600000000000001, pij.get(1,2), tolerance);
        assertEquals(0.600000000000001, pij.get(2,2), tolerance);
        assertEquals(0.600000000000001, pij.get(3,2), tolerance);
        assertEquals(0.600000000000001, pij.get(4,2), tolerance);
    }

}