/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.util.FastMath;
import org.apfloat.Apcomplex;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;
import org.apfloat.Apint;

public final class abatewhitt {
    private abatewhitt() {}

    public static final int precision = 35;

    public static UnaryOperator<Apcomplex> to_int = null;

    public static HashMap<Integer, BigInteger> factorial_memo = new HashMap<Integer, BigInteger>();

    public static Apcomplex laplace_result(Apcomplex d) {
        // to_int usage: see _kb/03-api-layer.md ("jline.lib.lti")
        return to_int.apply(d);
    }

    public static Apcomplex optimiser_laplace(Apcomplex d) {
        return Apcomplex.ONE.divide(Apcomplex.ONE.add(d.multiply(d)));
    }

    public static Apfloat optimise_lambda(ArrayList<Apcomplex> alpha, ArrayList<Apcomplex> omega) {
        Apfloat lambda = Apfloat.ONE;
        Apfloat bound_l = new Apfloat(0.78, (long) precision);
        Apfloat bound_r = new Apfloat(1, (long) precision);
        Apfloat bound = (bound_l.add(bound_r)).divide(new Apfloat(2));
        Apfloat max_lambda = Apfloat.ONE;
        HashMap<Apfloat, Apfloat> lambdavals = new HashMap<Apfloat, Apfloat>();
        while (true) {
            bound = (bound_l.add(bound_r)).divide(new Apint(2));
            Apfloat d = ApfloatMath.pi((long) precision).divide(new Apint(2));
            Apfloat ans = new Apfloat(1, (long) precision).divide(d);
            Apcomplex partialsum = new Apfloat(0, (long) precision);
            for (int i = 0; i < alpha.size(); i++) {
                Apcomplex scaled_alpha = alpha.get(i).multiply(lambda);
                partialsum = partialsum.add(omega.get(i).multiply(optimiser_laplace(scaled_alpha.divide(d))));
            }
            ans = ans.multiply(partialsum.real());
            if ((ans.multiply(lambda)).compareTo(bound) >= 0) {
                bound_l = bound;
                max_lambda = lambda;
                lambda = Apfloat.ONE;
                if ((bound_r.subtract(bound_l)).compareTo(new Apfloat(0.0001, (long) precision)) < 0) {
                    break;
                }
            } else if (lambda.compareTo(new Apfloat(9)) > 0) {
                lambda = Apfloat.ONE;
                bound_r = bound;
                if ((bound_r.subtract(bound_l)).compareTo(new Apfloat(0.0001, (long) precision)) < 0) {
                    // LATENT DEFECT, left as-is deliberately: the sibling branch above
                    // breaks on the same convergence test, so a right-bound convergence
                    // never terminates this loop. Behaviour is unchanged pending a real
                    // caller -- optimise_lambda is reached only through getResult, whose
                    // sole reference is the commented-out euler.java:52. See
                    // git show 449847e7b:_kb/log.md.
                }
            } else {
                lambda = lambda.add(new Apfloat(0.01, (long) precision));
            }
            if (!lambdavals.containsKey(lambda)) {
                lambdavals.put(lambda, ans.multiply(lambda));
            }
        }
        return max_lambda;
    }

    public static void getResult(ArrayList<Apcomplex> alpha, ArrayList<Apcomplex> omega) {
        TreeMap<Apfloat, Apfloat> res = new TreeMap<Apfloat, Apfloat>();
        int to_shift = 1;
        double max_val = 0.0;
        double unshifted_max = 0.0;
        System.out.println("Entered getResult");
        System.out.println("Size of alpha " + alpha.size());
        Apfloat lambda = optimise_lambda(alpha, omega);
        System.out.println("optimised lambda " + lambda);
        Apcomplex shifter = new Apcomplex(new Apfloat(-1, (long) precision).divide(new Apfloat(1, (long) precision)), Apfloat.ZERO);
        Apfloat d = new Apfloat(0.01, (long) precision);
        while (d.compareTo(new Apfloat(70.01, (long) precision)) <= 0) {
            Apfloat ans = Apfloat.ONE.divide(d);
            Apfloat ans_unshifted = Apfloat.ONE.divide(d);
            Apcomplex partialsum = Apcomplex.ZERO;
            Apcomplex partialsum_unshifted = Apcomplex.ZERO;
            for (int i = 0; i < alpha.size(); i++) {
                if (to_shift == 0) {
                    Apcomplex scaled_alpha = alpha.get(i).multiply(lambda);
                    partialsum = partialsum.add(omega.get(i).multiply(laplace_result(scaled_alpha.divide(d))));
                } else {
                    partialsum_unshifted = partialsum_unshifted.add(omega.get(i).multiply(laplace_result(alpha.get(i).divide(d))));
                    partialsum = partialsum.add(omega.get(i).multiply(laplace_result((alpha.get(i).divide(d)).add(shifter))));
                }
            }
            if (to_shift == 0) {
                ans = ans.multiply(partialsum.real()).multiply(lambda);
                res.put(d, ans);
                System.out.println(d.toString() + " " + res.get(d).doubleValue());
            } else {
                ans = ans.multiply(partialsum.real());
                ans_unshifted = ans_unshifted.multiply(partialsum_unshifted.real());
                ans = ans.multiply(ApfloatMath.exp(shifter.real().multiply(d)));
                max_val = FastMath.max(ans.doubleValue(), max_val);
                unshifted_max = FastMath.max(ans_unshifted.doubleValue(), unshifted_max);
                res.put(d, ans);
            }
            d = d.add(new Apfloat(0.05, (long) precision));
        }
        if (to_shift == 1) {
            System.out.println("Shifted max " + max_val);
            System.out.println("Unshifted max " + unshifted_max);
            for (Apfloat k : res.keySet()) {
                System.out.println(k.toString() + " " + res.get(k).doubleValue() / (max_val / unshifted_max));
            }
        }
    }

    public static BigInteger factorial(int n) {
        if (factorial_memo.containsKey(Integer.valueOf(n))) {
            return factorial_memo.get(Integer.valueOf(n));
        }
        BigInteger val = new BigInteger(BigInteger.ONE.toString());
        if (n == 1) {
            return val;
        }
        for (int i = 1; i <= n; i++) {
            val = val.multiply(new BigInteger(Integer.toString(i)));
        }
        factorial_memo.put(Integer.valueOf(n), val);
        return val;
    }

    public static BigInteger binom(int n, int k) {
        return factorial(n).divide(factorial(n - k)).divide(factorial(k));
    }
}
