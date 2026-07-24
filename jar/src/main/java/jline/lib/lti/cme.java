/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.apfloat.Apcomplex;
import org.apfloat.ApcomplexMath;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.apfloat.Apint;

public final class cme {
    private cme() {}

    public static IterativeLegendreGaussIntegrator ul = new IterativeLegendreGaussIntegrator(50, 1e-9, 1e-9);

    public static Apfloat moments(UnivariateFunction fun, final int k) {
        UnivariateFunction tomultiply = new UnivariateFunction() {
            public double value(double t) {
                return FastMath.pow(t, k);
            }
        };
        UnivariateFunction f = FunctionUtils.multiply(fun, tomultiply);
        return new Apfloat(ul.integrate(900000000, f, 0.0, 99.0), (long) abatewhitt.precision);
    }

    public static Apfloat scv(UnivariateFunction fun) {
        Apfloat res = moments(fun, 2).multiply(moments(fun, 0))
                .divide(ApfloatMath.pow(moments(fun, 1), 2));
        return res.subtract(Apfloat.ONE);
    }

    public static <T> ArrayList<T> deepcopy(ArrayList<T> tocopy) {
        ArrayList<T> result = new ArrayList<T>();
        for (T d : tocopy) result.add(d);
        return result;
    }

    public static ArrayList<Double> getnormalrandom(int n) {
        Random rand = new Random();
        ArrayList<Double> res = new ArrayList<Double>();
        for (int i = 0; i < n; i++) res.add(rand.nextGaussian());
        return res;
    }

    public static UnivariateFunction getmefunction(final ArrayList<Double> params) {
        UnivariateFunction temp = new UnivariateFunction() {
            public double value(double t) {
                return FastMath.exp(-t) * FastMath.pow(Math.cos(params.get(0) * t - params.get(1)), 2);
            }
        };
        for (int i = 0; i < params.size() - 2; i++) {
            final int finalI = i;
            UnivariateFunction temp2 = new UnivariateFunction() {
                public double value(double t) {
                    return FastMath.pow(Math.cos(params.get(0) * t - params.get(finalI + 2)), 2);
                }
            };
            temp = FunctionUtils.multiply(temp, temp2);
        }
        return temp;
    }

    public static ArrayList<Double> rechenberg(ArrayList<Double> x) {
        ArrayList<Double> opt = deepcopy(x);
        double sigma = 0.1;
        double c = 0.9;
        int siker = 0;
        int n = x.size();
        ArrayList<Double> best = deepcopy(opt);
        for (int i = 1; i <= 1000; i++) {
            double saveme = opt.get(0);
            ArrayList<Double> rands = getnormalrandom(n);
            ArrayList<Double> y = new ArrayList<Double>();
            for (int j = 0; j < n; j++) y.add(opt.get(j) + sigma * rands.get(j));
            if (y.get(0) <= 0) y.set(0, saveme);
            UnivariateFunction mef1 = getmefunction(y);
            UnivariateFunction mef2 = getmefunction(opt);
            if (scv(mef1).compareTo(scv(mef2)) < 0) {
                siker++;
                opt = deepcopy(y);
            }
            if (i % 20 == 0) {
                if (siker < 4) sigma *= c;
                else if (siker > 4) sigma /= c;
                siker = 0;
            }
            mef2 = getmefunction(opt);
            UnivariateFunction mef3 = getmefunction(best);
            if (scv(mef2).compareTo(scv(mef3)) < 0) {
                best = deepcopy(opt);
            }
        }
        return opt;
    }

    public static boolean[] binaryCheck(long num, int polysize) {
        boolean[] res = new boolean[polysize];
        res[0] = true;
        StringBuilder bin = new StringBuilder(Long.toBinaryString(num));
        while (bin.length() < polysize - 1) {
            bin.insert(0, "0");
        }
        bin.insert(0, "1");
        for (int i = 0; i < res.length; i++) {
            if (bin.charAt(i) == '1') res[i] = true;
        }
        return res;
    }

    public static <T> ArrayList<T> convert_to_arraylist(T[] input) {
        return new ArrayList<T>(Arrays.asList(input));
    }

    public static Pair<ArrayList<Apcomplex>, ArrayList<Apcomplex>> convert_to_laplace(ArrayList<Double> res) {
        Apfloat omega = new Apfloat(res.get(0), (long) abatewhitt.precision);
        ArrayList<Apfloat[]> polycoeff = new ArrayList<Apfloat[]>();
        for (int i = 1; i < res.size(); i++) {
            polycoeff.add(new Apfloat[] {
                    new Apfloat(-1 * res.get(i), (long) abatewhitt.precision),
                    omega
            });
        }
        ArrayList<Apfloat[]> polypart2 = new ArrayList<Apfloat[]>();
        Apfloat to_divide = ApfloatMath.pow(
                new Apfloat(4, (long) abatewhitt.precision),
                new Apfloat((long) (polycoeff.size() - 1)));
        Apfloat constant_term = Apfloat.ONE.divide(to_divide);
        Apfloat int_value = moments(getmefunction(res), 0);
        System.out.println(int_value);
        constant_term = constant_term.divide(int_value);
        for (long l = 0; l < (long) Math.pow(2.0, (double) (polycoeff.size() - 1)); l++) {
            boolean[] b = binaryCheck(l, polycoeff.size());
            Apfloat[] newpoly = new Apfloat[2];
            newpoly[0] = new Apfloat(0, (long) abatewhitt.precision);
            newpoly[1] = new Apfloat(0, (long) abatewhitt.precision);
            for (int i = 0; i < b.length; i++) {
                if (b[i]) {
                    newpoly[0] = newpoly[0].add(polycoeff.get(i)[0]);
                    newpoly[1] = newpoly[1].add(polycoeff.get(i)[1]);
                } else {
                    newpoly[0] = newpoly[0].subtract(polycoeff.get(i)[0]);
                    newpoly[1] = newpoly[1].subtract(polycoeff.get(i)[1]);
                }
            }
            polypart2.add(newpoly);
        }
        ArrayList<Apfloat[]> polypart3 = new ArrayList<Apfloat[]>(polypart2.size() * polypart2.size() + polypart2.size());
        HashSet<ArrayList<Apfloat>> hs = new HashSet<ArrayList<Apfloat>>();
        for (int i = 0; i < polypart2.size(); i++) {
            for (int j = i; j < polypart2.size(); j++) {
                if (i != j) {
                    Apfloat[] res2 = new Apfloat[2];
                    res2[0] = polypart2.get(i)[0].add(polypart2.get(j)[0]);
                    res2[1] = polypart2.get(i)[1].add(polypart2.get(j)[1]);
                    polypart3.add(res2);
                    res2 = new Apfloat[2];
                    res2[0] = polypart2.get(i)[0].subtract(polypart2.get(j)[0]);
                    res2[1] = polypart2.get(i)[1].subtract(polypart2.get(j)[1]);
                    polypart3.add(res2);
                } else {
                    Apfloat[] res2 = new Apfloat[2];
                    res2[0] = new Apfloat(2).multiply(polypart2.get(i)[0]);
                    res2[1] = new Apfloat(2).multiply(polypart2.get(i)[1]);
                    polypart3.add(res2);
                    hs.add(convert_to_arraylist(res2));
                    res2 = new Apfloat[2];
                    res2[0] = Apfloat.ZERO;
                    res2[1] = Apfloat.ZERO;
                    polypart3.add(res2);
                    hs.add(convert_to_arraylist(res2));
                }
            }
        }
        int ptr = 0;
        HashMap<Apcomplex, Apcomplex> eta_alphacombo = new HashMap<Apcomplex, Apcomplex>();
        while (ptr < polypart3.size()) {
            Apfloat[] doubles = polypart3.get(ptr++);
            if (doubles[1].equals(Apfloat.ZERO)) {
                Apcomplex C = new Apcomplex(ApfloatMath.cos(doubles[0]), Apfloat.ZERO);
                if (hs.contains(convert_to_arraylist(doubles))) {
                    C = C.divide(new Apint(2));
                }
                C = C.multiply(constant_term);
                if (eta_alphacombo.containsKey(Apcomplex.ONE)) {
                    eta_alphacombo.put(Apcomplex.ONE, eta_alphacombo.get(Apcomplex.ONE).add(C));
                } else {
                    eta_alphacombo.put(Apcomplex.ONE, C);
                }
            } else {
                Apcomplex C1 = new Apcomplex(new Apfloat(0, (long) abatewhitt.precision), doubles[0]);
                C1 = ApcomplexMath.exp(C1);
                C1 = C1.divide(new Apint(2));
                Apcomplex A1 = new Apcomplex(Apfloat.ONE, new Apfloat(-1).multiply(doubles[1]));
                Apcomplex C2 = new Apcomplex(new Apfloat(0, (long) abatewhitt.precision),
                        new Apfloat(-1).multiply(doubles[0]));
                C2 = ApcomplexMath.exp(C2);
                C2 = C2.divide(new Apint(2));
                Apcomplex A2 = new Apcomplex(Apfloat.ONE, doubles[1]);
                if (hs.contains(convert_to_arraylist(doubles))) {
                    C1 = C1.divide(new Apint(2));
                    C2 = C2.divide(new Apint(2));
                }
                C1 = C1.multiply(constant_term);
                C2 = C2.multiply(constant_term);
                if (eta_alphacombo.containsKey(A1)) {
                    eta_alphacombo.put(A1, eta_alphacombo.get(A1).add(C1));
                } else {
                    eta_alphacombo.put(A1, C1);
                }
                if (eta_alphacombo.containsKey(A2)) {
                    eta_alphacombo.put(A2, eta_alphacombo.get(A2).add(C2));
                } else {
                    eta_alphacombo.put(A2, C2);
                }
            }
        }
        ArrayList<Apcomplex> alpha = new ArrayList<Apcomplex>();
        ArrayList<Apcomplex> eta = new ArrayList<Apcomplex>();
        for (Apcomplex e : eta_alphacombo.keySet()) {
            alpha.add(e);
            eta.add(eta_alphacombo.get(e));
        }
        return new Pair<ArrayList<Apcomplex>, ArrayList<Apcomplex>>(eta, alpha);
    }
}
