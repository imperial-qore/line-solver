package jline.lib.ltinversion;

import org.apfloat.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.function.UnaryOperator;

public class abatewhitt {
    static UnaryOperator<Apcomplex> to_int = null;
    static final int precision = 35;

    public static Apcomplex laplace_result(Apcomplex d) {
        /**
         * If you want to evaluate a specific function, override to_int. For example,
         * to_int = (t -> (Apcomplex.ONE.divide(t)).divide(ApcomplexMath.exp(t).add(Apcomplex.ONE))); is (1/t)/(exp(t) + 1)
         * to_int = (t -> Apcomplex.ONE.divide(Apcomplex.ONE.add(d.multiply(d)))); is 1/(1 + s^2) (sine function)
         */
        return to_int.apply(d);
    }

    public static Apcomplex optimiser_laplace(Apcomplex d) {
        /**
         * This is used for optimising lambda
         */
        return Apcomplex.ONE.divide(Apcomplex.ONE.add(d.multiply(d)));
    }

    public static Apfloat optimise_lambda(ArrayList<Apcomplex> alpha, ArrayList<Apcomplex> omega) {
        Apfloat lambda = Apfloat.ONE;
        Apfloat bound_l = new Apfloat(0.78, precision);
        Apfloat bound_r = new Apfloat(1, precision);
        Apfloat bound = (bound_l.add(bound_r)).divide(new Apfloat(2));
        Apfloat max_lambda = Apfloat.ONE;
        HashMap<Apfloat, Apfloat> lambdavals = new HashMap<>();
        while (true) {
            bound = (bound_l.add(bound_r)).divide(new Apint(2));
            Apfloat d = ApfloatMath.pi(precision).divide(new Apint(2)); // pi/2
            Apfloat ans = new Apfloat(1, precision).divide(d);
            Apcomplex partialsum = new Apfloat(0, precision);
            for (int i = 0; i < alpha.size(); i++) {
                Apcomplex scaled_alpha = alpha.get(i).multiply(lambda);
                partialsum = partialsum.add(omega.get(i).multiply(abatewhitt.optimiser_laplace(scaled_alpha.divide(d))));
            }
            ans = ans.multiply(partialsum.real());
            if ((ans.multiply(lambda)).compareTo(bound) >= 0) {
                //       System.out.println(bound + " " + bound_l + " " + bound_r + " ");
                bound_l = bound;
                max_lambda = lambda;
                lambda = Apfloat.ONE;
                if ((bound_r.subtract(bound_l)).compareTo(new Apfloat(0.0001, precision)) < 0) {
                    break;
                }
            } else if (lambda.compareTo(new Apfloat(9)) > 0) {
                lambda = Apfloat.ONE;
                //     System.out.println(bound + " " + bound_l + " " + bound_r + " ");
                bound_r = bound;
                if ((bound_r.subtract(bound_l)).compareTo(new Apfloat(0.0001, precision)) < 0) {
                    System.out.println(lambdavals);
                }
            } else
                lambda = lambda.add(new Apfloat(0.01, precision));
            if (!lambdavals.containsKey(lambda))
                lambdavals.put(lambda, ans.multiply(lambda));
        }
        //   System.out.println(lambdavals);
        return max_lambda;
    }

    public static void getResult(ArrayList<Apcomplex> alpha, ArrayList<Apcomplex> omega) {
        TreeMap<Apfloat, Apfloat> res = new TreeMap<>();
        int to_shift = 1; // determine if shifting is necessary
        double max_val = 0;
        double unshifted_max = 0;
        System.out.println("Entered getResult");
        System.out.println("Size of alpha " + alpha.size());
        Apfloat lambda = optimise_lambda(alpha, omega);
        System.out.println("optimised lambda " + lambda);
        Apcomplex shifter = new Apcomplex(new Apfloat(-1, precision).divide(new Apfloat(1, precision)), Apfloat.ZERO);
        //       System.out.println("Printing results");
        for (Apfloat d = new Apfloat(0.01, precision); d.compareTo(new Apfloat(70.01, precision)) <= 0; d = d.add(new Apfloat(0.05, precision))) {
            //    BigDecimal ans = BigMath.ONE.divide(d, precision, RoundingMode.HALF_EVEN);
            Apfloat ans = Apfloat.ONE.divide(d);
            Apfloat ans_unshifted = Apfloat.ONE.divide(d);
            Apcomplex partialsum = Apcomplex.ZERO;
            Apcomplex partialsum_unshifted = Apcomplex.ZERO;
            for (int i = 0; i < alpha.size(); i++) {
                if (to_shift == 0) {
                    Apcomplex scaled_alpha = alpha.get(i).multiply(lambda);
                    partialsum = partialsum.add(omega.get(i).multiply(abatewhitt.laplace_result(scaled_alpha.divide(d))));
                } else {
                    partialsum_unshifted = partialsum_unshifted.add(omega.get(i).multiply(abatewhitt.laplace_result(alpha.get(i).divide(d))));
                    partialsum = partialsum.add(omega.get(i).multiply(abatewhitt.laplace_result((alpha.get(i).divide(d)).add(shifter))));
                }
            }
            if (to_shift == 0) {
                ans = ans.multiply((partialsum.real())).multiply(lambda);
                res.put(d, ans);
                System.out.println(d + " " + res.get(d).doubleValue());
            } else {
                ans = ans.multiply((partialsum.real()));
                ans_unshifted = ans_unshifted.multiply((partialsum_unshifted.real()));
                //     ans = ApfloatMath.exp(ans.multiply(shifter).real());
                ans = ans.multiply(ApfloatMath.exp(shifter.real().multiply(d)));
                max_val = Math.max(ans.doubleValue(), max_val);
                unshifted_max = Math.max(ans_unshifted.doubleValue(), unshifted_max);
                res.put(d, ans);
            }
        }
        if (to_shift == 1) {
            System.out.println("Shifted max " + max_val);
            System.out.println("Unshifted max " + unshifted_max);
            for (Apfloat d : res.keySet()) {
                System.out.println(d + " " + res.get(d).doubleValue() / (max_val / unshifted_max));
            }
        }
    }

    static HashMap<Integer, BigInteger> factorial_memo = new HashMap<>();

    public static BigInteger factorial(int n) {
        if (factorial_memo.containsKey(n))
            return factorial_memo.get(n);
        BigInteger val = new BigInteger(String.valueOf(BigInteger.ONE));
        if (n == 1)
            return val;
        for (int i = 1; i <= n; i++)
            val = val.multiply(new BigInteger(String.valueOf(i)));
        factorial_memo.put(n, val);
        return val;
    }

    public static BigInteger binom(int n, int k) {
        // binomial formula: n!/(n - k!)k!
        return factorial(n).divide(factorial(n - k)).divide(factorial(k));
    }
}
