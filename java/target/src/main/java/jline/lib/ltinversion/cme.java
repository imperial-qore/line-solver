package jline.lib.ltinversion;
import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.util.Pair;
import org.apfloat.*;

import java.math.BigDecimal;
import java.util.*;

public class cme {
    static IterativeLegendreGaussIntegrator ul = new IterativeLegendreGaussIntegrator(50, 1e-9, 1e-9);

    public static Apfloat moments(UnivariateFunction fun, int k) {
        //     UnivariateIntegrator integrator = new RombergIntegrator(1e-9, 1e-9, 3, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction tomultiply = (t -> Math.pow(t, k));
        fun = FunctionUtils.multiply(fun, tomultiply);
        //    System.out.println(fun);
        //  return (laguerre.Laguerre(fun, 5000));
        //    return (customromberg.starter(fun, 0, 99));

        // return BigDecimal.valueOf(integrator.integrate(900000000, fun, 0, 99));
        // Legendre-Gaussian only
        // return BigDecimal.valueOf(ul.integrate(90000000, gaussiancombine(fun, 0, 99), -1, 1));
        return new Apfloat((ul.integrate(900000000, fun, 0, 99)), abatewhitt.precision);
    }

//    public static UnivariateFunction gaussiancombine(UnivariateFunction parent, double lb, double ub) {
//        // parent(child))
//        double eta = (ub - lb) / 2;
//        double mu = (ub + lb) / 2;
//        UnivariateFunction child = (t -> eta * (eta * t + mu));
//        return FunctionUtils.compose(parent, child);
//    }

    public static Apfloat scv(UnivariateFunction fun) {
        Apfloat res = moments(fun, 2).multiply(moments(fun, 0)).divide(ApfloatMath.pow(moments(fun, 1), 2));
        return res.subtract(Apfloat.ONE);
    }

    public static <T> ArrayList<T> deepcopy(ArrayList<T> tocopy) {
        ArrayList<T> result = new ArrayList<>();
        for (T d : tocopy)
            result.add(d);
        return result;
    }

    public static ArrayList<Double> getnormalrandom(int n) {
        // returns n randomly distributed numbers
        Random rand = new Random();
        ArrayList<Double> res = new ArrayList<>();
        for (int i = 0; i < n; i++)
            res.add(rand.nextGaussian());
        return res;
    }

    public static UnivariateFunction getmefunction(ArrayList<Double> params) {
        UnivariateFunction temp = t -> Math.exp(-t) * Math.pow(Math.cos(params.get(0) * t - params.get(1)), 2);
        for (int i = 0; i < params.size() - 2; i++) {
            int finalI = i;
            UnivariateFunction temp2 = t -> Math.pow(Math.cos(params.get(0) * t - params.get(finalI + 2)), 2);
            temp = FunctionUtils.multiply(temp, temp2);
        }
        return temp;
    }

    public static ArrayList<Double> linspace(double min, double max, int points) {
        // https://stackoverflow.com/a/17817524/19370273
        //    double[] d = new double[points];
        ArrayList<Double> d = new ArrayList<>();
        for (int i = 0; i < points; i++) {
            d.add(min + i * (max - min) / (points - 1));
        }
        return d;
    }

    // only portion not under high precision
    // code directly adapted from http://webspn.hit.bme.hu/~illes/mincvnum.zip (as part of the "Concentrated matrix exponential distributions" paper)
    public static ArrayList<Double> rechenberg(ArrayList<Double> x) {
        ArrayList<Double> scv_opts = new ArrayList<>();
        ArrayList<Double> opt = deepcopy(x);
        double sigma = 0.1;
        double c = 0.9;
        int siker = 0;
        int n = x.size();
        ArrayList<Double> best = deepcopy(opt);
        for (int i = 1; i <= 1000; i++) {
            double saveme = opt.get(0);
            ArrayList<Double> rands = getnormalrandom(n);
            ArrayList<Double> y = new ArrayList<>();
            for (int j = 0; j < n; j++)
                y.add(opt.get(j) + sigma * rands.get(j));
            if (y.get(0) <= 0)
                y.set(0, saveme);
            UnivariateFunction mef1 = getmefunction(y); // in preparation for f(y) in rechenberg.m
            UnivariateFunction mef2 = getmefunction(opt); // as above
            if (scv(mef1).compareTo(scv(mef2)) < 0) { // to prevent cases of negative SCV escaping
                siker++;
                opt = deepcopy(y);
            }
            if (i % 20 == 0) {
                if (siker < 4)
                    sigma *= c;
                else if (siker > 4)
                    sigma /= c;
                siker = 0;
            }
            mef2 = getmefunction(opt);
            UnivariateFunction mef3 = getmefunction(best);
            if (scv(mef2).compareTo(scv(mef3)) < 0) {
                best = deepcopy(opt);
            }
            //   scv_opts.add(scv(getmefunction(opt)));
        }
        //     System.out.println("SCV opts");
        //    System.out.println(scv_opts);
        //    double arr = scv(getmefunction(opt));
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
            if (bin.charAt(i) == '1')
                res[i] = true;
        }
        return res;
    }

    public static <T> ArrayList<T> convert_to_arraylist(T[] input) {
        return new ArrayList<>(Arrays.asList(input));
    }

    public static Pair<ArrayList<Apcomplex>, ArrayList<Apcomplex>> convert_to_laplace(ArrayList<Double> res) {
        // res[0] = omega, the rest are phi variables
        Apfloat omega = new Apfloat(res.get(0), abatewhitt.precision);
        /*
        Given problem is
        (pi_{i = 1}^n cos (omega - phi_i))^2
        Then, use fact that (cos A cos B) = (1/2 (cos (A + B) + cos (A - B))) so that we get something of the form
        cos (A +- B +- C +- ...) + summed all over
        After that, the goal is to square this result. This will get a O((2^(n - 1))^2) sum.
        Then, use the fact that 2 cos(P) cos(Q) = (cos(P + Q) + cos(P - Q)). Where P = Q, this will naturally resolve to 1/2 (cos 2P + 1)
        At this point, if the result is real, we can mark it as such and sum them all together; if not, that's Apcomplex and should be resolved together.
        For the Apcomplex variables, we now have the form W = e^-t cos R, where R has both a real and an Apcomplex aspect. Note that e^(x + iy) = e^x (cos y + i sin y).
        Then, W can be split as W = 1/2 (U + V), where U = e^-t (cos R + i sin R) and V = e^-t (cos R - i sin R). Note that this is the same as -R in U
        After that, U and V can be written as e^(-t + R) and e^(-t - R), respectively. This is the same as e^-t(e^R + e^(-R)). Finally, take the t part outside.
         */
        // STEP 1: Get the polynomial cos(P) where P = omega*t - phi_i
        ArrayList<Apfloat[]> polycoeff = new ArrayList<>();
        for (int i = 1; i < res.size(); i++) {
            polycoeff.add(new Apfloat[]{new Apfloat(-1 * res.get(i), abatewhitt.precision), omega}); // -phi_i + omega_t
        }
        // cos (B + AX) where A = omega and B = phi
        // get part 2: this would be 2^{n - 1} polynomials of the form A +- B +- ...
        ArrayList<Apfloat[]> polypart2 = new ArrayList<>();
        Apfloat to_divide = ApfloatMath.pow(new Apfloat(4, abatewhitt.precision), new Apfloat(polycoeff.size() - 1));
        Apfloat constant_term = Apfloat.ONE.divide(to_divide);
        // remember that the PDFs are unnormalised, so we will need to normalise
        // c int(f(t)) = 1
        Apfloat int_value = moments(getmefunction(res), 0);
        System.out.println(int_value);
        constant_term = constant_term.divide(int_value);
        //  constant_term = constant_term * (1 / int_value);
        for (long l = 0; l < (long) (Math.pow(2, polycoeff.size() - 1)); l++) {
            boolean[] b = binaryCheck(l, polycoeff.size());
            // loop and add as appropriate
            Apfloat[] newpoly = new Apfloat[2];
            newpoly[0] = new Apfloat(0, abatewhitt.precision);
            newpoly[1] = new Apfloat(0, abatewhitt.precision);
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
        // part 3: do the squaring
        // also part 4: cos^2 x = (1/2 (cos 2x + 1)) -> constant term is 1/2
        ArrayList<Apfloat[]> polypart3 = new ArrayList<>(polypart2.size() * polypart2.size() + polypart2.size());
        HashSet<ArrayList<Apfloat>> hs = new HashSet<>();
        for (int i = 0; i < polypart2.size(); i++) {
            for (int j = i; j < polypart2.size(); j++) {
                // if (i != j), then 2 cos P cos Q = cos(P + Q) + cos (P - Q)
                if (i != j) {
                    Apfloat[] res2 = new Apfloat[2];
                    res2[0] = polypart2.get(i)[0].add(polypart2.get(j)[0]);
                    res2[1] = polypart2.get(i)[1].add(polypart2.get(j)[1]);
                    polypart3.add(res2);
                    // (cos (P - Q))
                    res2 = new Apfloat[2];
                    res2[0] = polypart2.get(i)[0].subtract(polypart2.get(j)[0]);
                    res2[1] = polypart2.get(i)[1].subtract(polypart2.get(j)[1]);
                    polypart3.add(res2);
                } else {
                    // cos^2 P = 1/2 (cos 2P + 1).
                    // handle 1/2 term, looks like we will need to do it via a hashset
                    // this is because the entire structure is based on polynomials till date
                    Apfloat[] res2 = new Apfloat[2];
                    res2[0] = new Apfloat(2).multiply(polypart2.get(i)[0]);
                    res2[1] = new Apfloat(2).multiply(polypart2.get(i)[1]);
                    polypart3.add(res2);
                    // add to hashset so that we can handle this specifically in the end
                    hs.add(convert_to_arraylist(res2));
                    // also handle the cos (0) part
                    res2 = new Apfloat[2];
                    res2[0] = Apfloat.ZERO;
                    res2[1] = Apfloat.ZERO;
                    polypart3.add(res2);
                    hs.add(convert_to_arraylist(res2));
                }
            }
        }
        //   constant_term /= 2;
        /* now, for all the polypart3 cases, check whether it is dependent on $t$.
        if it is, get a conjugate and finalise in that way [real, Apcomplex] + t*[real, Apcomplex]
        hint: e^-t (cos P) = 1/2 e^-t (cos P + i sin P) + 1/2 e^-t (cos P - i sin P)
        = 1/2 e^(-t + iP) + 1/2 e^(-t - iP)
        */
        int ptr = 0;
        HashMap<Apcomplex, Apcomplex> eta_alphacombo = new HashMap<>();
        while (ptr < polypart3.size()) {
            // check if it even depends on t
            Apfloat[] doubles = polypart3.get(ptr++);
            if (doubles[1].equals(Apfloat.ZERO)) {
                // only real! alpha_i is easily 1
                // these terms can be directly evaluated
                //    alpha.add(Apcomplex.ONE);
                Apcomplex C = new Apcomplex((ApfloatMath.cos(doubles[0])), Apfloat.ZERO);
                if (hs.contains(convert_to_arraylist(doubles))) {
                    // needs to be divided by 2 as before
                    C = C.divide(new Apint(2));
                }
                C = C.multiply(constant_term);
                //       eta.add(C);
                if (eta_alphacombo.containsKey(Apcomplex.ONE)) {
                    eta_alphacombo.put(Apcomplex.ONE, eta_alphacombo.get(Apcomplex.ONE).add(C));
                } else {
                    eta_alphacombo.put(Apcomplex.ONE, C);
                }
            } else {
                // computing eta1
                Apcomplex C1 = new Apcomplex(new Apfloat(0, abatewhitt.precision), doubles[0]); //
                C1 = ApcomplexMath.exp(C1);
                C1 = C1.divide(new Apint(2));
                // computing alpha1
                Apcomplex A1 = new Apcomplex(Apfloat.ONE, new Apfloat(-1).multiply(doubles[1])); // {1, -i omega}
                // eta2
                Apcomplex C2 = new Apcomplex(new Apfloat(0, abatewhitt.precision), new Apfloat(-1).multiply(doubles[0])); //
                C2 = ApcomplexMath.exp(C2);
                C2 = C2.divide(new Apint(2));
                // alpha2
                Apcomplex A2 = new Apcomplex(Apfloat.ONE, doubles[1]); // {1, i omega}
                if (hs.contains(convert_to_arraylist(doubles))) {
                    // the result MUST BE DIVIDED BY 2
                    // (only for etas though!)
                    C1 = C1.divide(new Apint(2));
                    C2 = C2.divide(new Apint(2));
                }
                // divide C1 and C2 with old normalising constant
                C1 = C1.multiply(constant_term);
                C2 = C2.multiply(constant_term);
                // and add it to alpha and eta
                //   eta.add(C1);
                //  alpha.add(A1);
                if (eta_alphacombo.containsKey(A1)) {
                    eta_alphacombo.put(A1, eta_alphacombo.get(A1).add(C1));
                } else {
                    eta_alphacombo.put(A1, C1);
                }
                //  eta.add(C2);
                //   alpha.add(A2);
                if (eta_alphacombo.containsKey(A2)) {
                    eta_alphacombo.put(A2, eta_alphacombo.get(A2).add(C2));
                } else {
                    eta_alphacombo.put(A2, C2);
                }
            }
        }
        ArrayList<Apcomplex> alpha = new ArrayList<>();
        ArrayList<Apcomplex> eta = new ArrayList<>();
        for (Apcomplex e : eta_alphacombo.keySet()) {
            alpha.add(e);
            eta.add(eta_alphacombo.get(e));
        }
        return new Pair<>(eta, alpha);
    }

    public static void main(String[] args) {
        double best_fin = Double.MAX_VALUE;
        int number_of_iterations = 0; // how many iterations to run
        cyclic_queues.test_function(); // get the function that needs to be integrated - result should be PDF
        for (int i = 0; i < number_of_iterations; i++) {
            System.out.println("Iteration number " + (i + 1));
            // this code is directly adapted from http://webspn.hit.bme.hu/~illes/mincvnum.zip
            ArrayList<Double> input = new ArrayList<>();
            input.add(0.3);
            input.add(0.1);
            ArrayList<Double> lin = linspace(Math.PI / 2, Math.PI, 12);
            input.addAll(lin);
            ArrayList<Double> res = rechenberg(input);
            System.out.println(res);
            Apfloat arr = scv(getmefunction(res));
            // END of adapted code

            Pair<ArrayList<Apcomplex>, ArrayList<Apcomplex>> laplacetest = convert_to_laplace(res);
            abatewhitt.getResult(laplacetest.getValue(), laplacetest.getKey());
            best_fin = Math.min(best_fin, arr.doubleValue());
            System.out.println("SCV " + arr);
        }
        ArrayList<Double> test_optimal = new ArrayList<>();
        System.out.println("custom_test");
        Collections.addAll(test_optimal, 0.3730, 0.0794, 1.6030, 1.7370, 1.9543, 2.2255, 2.5323, 2.8652);
        //  Collections.addAll(test_optimal, 0.3387, 0.1990, 0.4595, 1.5901, 1.6718, 1.8077, 1.9803, 2.1765, 2.3889, 2.6132, 2.8471, 3.0896);
        System.out.println("test " + test_optimal);
        Pair<ArrayList<Apcomplex>, ArrayList<Apcomplex>> laplacetest = convert_to_laplace(test_optimal);
        abatewhitt.getResult(laplacetest.getValue(), laplacetest.getKey());
        //     System.out.println("Mean SCV " + fin / number_of_iterations);
        //    System.out.println("Best SCV " + best_fin);
    }
}