package jline.lib.ltinversion;

import org.apfloat.Apcomplex;
import org.apfloat.ApcomplexMath;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

import java.util.ArrayList;
import java.util.function.UnaryOperator;

public class cyclic_queues {
    static ArrayList<ArrayList<Integer>> Z = new ArrayList<>();
    public static void getcombinations(ArrayList<Integer> current_list, int M, int N, int current_sum) {
        if (current_list.size() == M) {
            if (current_sum == N - 1) {
                Z.add(current_list);
            }
            return;
        }
        for (int i = 0; i <= N - 1 - current_sum; i++) {
            // deep copy
            ArrayList<Integer> dc = cme.deepcopy(current_list);
            dc.add(i);
            getcombinations(dc, M, N, current_sum + i);
        }
    }

    public static Apfloat getP(ArrayList<Integer> J, ArrayList<Double> alpha) {
        // normalising constant
        Apfloat denom = Apfloat.ZERO;
        for (ArrayList<Integer> ki : Z) {
            Apfloat prod = Apfloat.ONE;
            for (int i = 0; i < J.size(); i++) {
                prod = prod.multiply(ApfloatMath.pow(new Apfloat(alpha.get(i)), ki.get(i)));
            }
            denom = denom.add(prod);
        }
        // handle the DENOMINATOR
        Apfloat res = Apfloat.ONE.divide(denom);
        for (int j = 0; j < J.size(); j++) {
            res = res.multiply(ApfloatMath.pow(new Apfloat(alpha.get(j)), J.get(j)));
        }
        return res;
    }

    public static UnaryOperator<Apcomplex> to_invert(int M, int N, ArrayList<Double> alpha) {
        Z = new ArrayList<>();
        getcombinations(new ArrayList<>(), M, N, 0);
        UnaryOperator<Apcomplex> f = (t -> Apcomplex.ZERO);
        for (ArrayList<Integer> J : Z) {
            Apfloat p = getP(J, alpha);
            UnaryOperator<Apcomplex> toadd = (t -> new Apcomplex(p, Apfloat.ZERO));
            // product term
            for (int i = 0; i < M; i++) {
                int finalI = i;
                UnaryOperator<Apcomplex> tomul = t -> ApcomplexMath.pow(Apcomplex.ONE.divide((Apcomplex.ONE.add(new Apcomplex(new Apfloat(alpha.get(finalI)), Apfloat.ZERO).multiply(t)))), J.get(finalI) + 1);
                //   toadd = toadd.multiply(tomul);
                toadd = function_wrapper.multiply(toadd, tomul);
                //      toadd = FunctionUtils.multiply(toadd, tomul);
            }
            f = function_wrapper.add(f, toadd);
            // f = FunctionUtils.add(f, toadd);
        }
        return f; // inverting this should return the PDF of the queue
    }

    public static void test_function() {
        int m = 2;
        int n = 3;
        getcombinations(new ArrayList<>(), m, n, 0);
        System.out.println(Z);
        System.out.println("m = " + m + " n = " + n);
        ArrayList<Double> alpha = new ArrayList<>();
        // related to arrival rate parameters
//        alpha.add(2d);
//        alpha.add(3d);
//        alpha.add(2d);
//        alpha.add(1d);
//        alpha.add(4d);
        alpha.add(1d);
        alpha.add(1d);
        UnaryOperator<Apcomplex> uf = to_invert(m, n, alpha);
        abatewhitt.to_int = uf;
        System.out.println("alpha = " + alpha);
        // this is to be integrated to find the transform
    }
    public static void main(String[] args) {
        /*
         */
        //  System.out.println(uf.toString());
        System.out.println("done");
    }
}