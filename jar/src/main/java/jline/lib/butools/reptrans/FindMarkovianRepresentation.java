/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G Horvath, M Telek, "A minimal representation of Markov arrival processes
 * and a moments matching method," Performance Evaluation 64:(9-12) pp. 1153-1168. (2007)
 */
package jline.lib.butools.reptrans;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import java.util.function.BiFunction;

import jline.util.matrix.Matrix;

public final class FindMarkovianRepresentation {
    private FindMarkovianRepresentation() {}

    /** Transformation function applied to a representation. */
    public interface TransFun {
        List<Matrix> apply(List<Matrix> rep, Matrix B);
    }

    /** Evaluation function returning distance from Markovian representation. */
    public interface EvalFun {
        double apply(List<Matrix> rep, int k);
    }

    /**
     * Overload accepting Java functional-interface flavors of the arguments.
     */
    public static List<Matrix> findMarkovianRepresentation(
            List<Matrix> rep,
            final TransFun transfun,
            final EvalFun evalfun,
            double precision) {
        BiFunction<List<Matrix>, Matrix, List<Matrix>> tf =
                new BiFunction<List<Matrix>, Matrix, List<Matrix>>() {
                    @Override
                    public List<Matrix> apply(List<Matrix> r, Matrix B) {
                        return transfun.apply(r, B);
                    }
                };
        BiFunction<List<Matrix>, Integer, Double> ef =
                new BiFunction<List<Matrix>, Integer, Double>() {
                    @Override
                    public Double apply(List<Matrix> r, Integer k) {
                        return Double.valueOf(evalfun.apply(r, k.intValue()));
                    }
                };
        return findMarkovianRepresentation(rep, tf, ef, precision);
    }

    /**
     * Obtains a Markovian representation from a non-Markovian one while keeping
     * the size the same, by applying a series of elementary transformations.
     */
    public static List<Matrix> findMarkovianRepresentation(
            List<Matrix> rep,
            BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun,
            BiFunction<List<Matrix>, Integer, Double> evalfun,
            double precision) {
        if (evalfun.apply(rep, Integer.valueOf(0)).doubleValue() < precision) {
            return rep;
        }

        List<Matrix> nrep = copyList(rep);
        int M = nrep.get(0).getNumCols();
        double b = 0.5;
        double odist = Double.MAX_VALUE;

        while (b > precision / 2) {
            for (int m = 0; m < M * M; m++) {
                for (int k = 0; k < 4; k++) {
                    Pair<List<Matrix>, Double> result = minimize(nrep, M * M, b, k, evalfun, transfun);
                    nrep = result.getFirst();
                    double ddist = result.getSecond().doubleValue();
                    if (ddist < precision) {
                        return nrep;
                    }
                }
                if (odist <= evalfun.apply(nrep, Integer.valueOf(0)).doubleValue()) {
                    break;
                }
                odist = evalfun.apply(nrep, Integer.valueOf(0)).doubleValue();
            }
            b /= 2.0;
        }

        return nrep;
    }

    public static List<Matrix> findMarkovianRepresentation(
            List<Matrix> rep,
            BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun,
            BiFunction<List<Matrix>, Integer, Double> evalfun) {
        return findMarkovianRepresentation(rep, transfun, evalfun, 1e-7);
    }

    /**
     * Minimization helper function.
     */
    private static Pair<List<Matrix>, Double> minimize(
            List<Matrix> orep,
            int iters,
            double b,
            int k,
            BiFunction<List<Matrix>, Integer, Double> evalfun,
            BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun) {
        double lastdist = evalfun.apply(orep, Integer.valueOf(k)).doubleValue();
        List<Matrix> bestrep = copyList(orep);
        List<Matrix> currentRep = copyList(orep);

        for (int i = 0; i < iters; i++) {
            Pair<List<Matrix>, Double> result = elementary(currentRep, b, k, evalfun, transfun);
            List<Matrix> newRep = result.getFirst();
            double dist = result.getSecond().doubleValue();
            if (dist >= lastdist) {
                break;
            } else {
                lastdist = dist;
                bestrep = copyList(newRep);
                currentRep = newRep;
            }
        }

        return new Pair<List<Matrix>, Double>(bestrep, Double.valueOf(lastdist));
    }

    /**
     * Elementary transformation helper function.
     */
    private static Pair<List<Matrix>, Double> elementary(
            List<Matrix> erep,
            double b,
            int k,
            BiFunction<List<Matrix>, Integer, Double> evalfun,
            BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun) {
        double bestdist = evalfun.apply(erep, Integer.valueOf(k)).doubleValue();
        List<Matrix> bestrep = copyList(erep);
        int repSize = erep.get(0).getNumCols();

        for (int i = 0; i < repSize; i++) {
            for (int j = 0; j < repSize; j++) {
                if (i != j) {
                    Matrix Bpos = Matrix.eye(repSize);
                    Bpos.set(i, j, b);
                    Bpos.set(i, i, 1.0 - b);

                    List<Matrix> newrepPos = transfun.apply(erep, Bpos);
                    double newdistPos = evalfun.apply(newrepPos, Integer.valueOf(k)).doubleValue();

                    if (newdistPos < bestdist) {
                        bestrep = copyList(newrepPos);
                        bestdist = newdistPos;
                    }

                    Matrix Bneg = Matrix.eye(repSize);
                    Bneg.set(i, j, -b);
                    Bneg.set(i, i, 1.0 + b);

                    List<Matrix> newrepNeg = transfun.apply(erep, Bneg);
                    double newdistNeg = evalfun.apply(newrepNeg, Integer.valueOf(k)).doubleValue();

                    if (newdistNeg < bestdist) {
                        bestrep = copyList(newrepNeg);
                        bestdist = newdistNeg;
                    }
                }
            }
        }

        return new Pair<List<Matrix>, Double>(bestrep, Double.valueOf(bestdist));
    }

    private static List<Matrix> copyList(List<Matrix> src) {
        List<Matrix> out = new ArrayList<Matrix>(src.size());
        for (int i = 0; i < src.size(); i++) {
            out.add(src.get(i).copy());
        }
        return out;
    }
}
