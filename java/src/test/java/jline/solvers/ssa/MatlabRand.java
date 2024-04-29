package jline.solvers.ssa;

import jline.util.Maths;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class MatlabRand {

    // percentage difference: allow for rounding inconsistency between matlab and jline
    protected static final double allowedDeviation = 0.1;

    // Mersenne-twister generated numbers with starting seed 1 from matlab
//    final static String random_numbers_file = "RNG";

//    public static void readInRandNums() {
//        ClassLoader classloader = Thread.currentThread().getContextClassLoader();
//        InputStream is = classloader.getResourceAsStream(random_numbers_file);
//        List<Double> numbers = new ArrayList<Double>();
//
//        Scanner myReader = new Scanner(is);
//        while (myReader.hasNextLine()) {
//            String data = myReader.nextLine();
//            numbers.add(Double.parseDouble(data));
//        }
//        myReader.close();
//        Maths.setRandomNumbers(numbers);
//    }

    protected static boolean closeEnough(List<Double> correct, List<Double> test, double deviation) {
        for (int i = 0; i < correct.size(); i++) {
            if (!closeEnough(correct.get(i), test.get(i), deviation)) {
                System.out.println("BROKEN");
                System.out.println(correct);
                System.out.println(test);
                System.out.println(i);
                return false;
            }
        }
        return true;
    }

    // deviation is the percentage difference we allow between test and correct
    private static boolean closeEnough(double correct, double test, double deviation) {
        return Math.abs(correct - test) <= deviation * correct / 100;
    }
}
