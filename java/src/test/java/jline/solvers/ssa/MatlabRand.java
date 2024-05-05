package jline.solvers.ssa;

import jline.util.Maths;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class MatlabRand {

    // percentage difference: allow for rounding inconsistency between matlab and jline
    protected static final double allowedDeviation = 0.1;

    protected static boolean closeEnough(List<Double> correct, List<Double> test, double deviation) {
        for (int i = 0; i < correct.size(); i++) {
            if (!closeEnough(correct.get(i), test.get(i), deviation)) {
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
