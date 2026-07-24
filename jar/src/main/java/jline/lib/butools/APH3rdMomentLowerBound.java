package jline.lib.butools;

import org.apache.commons.math3.complex.Complex;

import jline.GlobalConstants;

public final class APH3rdMomentLowerBound {
    private APH3rdMomentLowerBound() {}

    public static double APH3rdMomentLowerBound(double m1, double m2, int n) {
        double ni2 = m2 / m1 / m1;
        if (ni2 < (n + 1.0) / n) {
            return GlobalConstants.Inf;
        } else if (ni2 < (n + 4.0) / (n + 1.0)) {
            Complex n2 = new Complex(ni2, 0.0);
            Complex three = new Complex(3.0, 0.0);
            Complex negTwo = new Complex(-2.0, 0.0);

            Complex numeratorP = (new Complex(n + 1.0, 0.0)).multiply(n2.subtract(2.0));
            Complex denominatorP = three.multiply(n2).multiply(new Complex(n - 1.0, 0.0));
            Complex sqrtPart = (negTwo.multiply(new Complex(Math.sqrt(n + 1.0), 0.0)))
                    .divide(new Complex(Math.sqrt(-3.0 * n * ni2 + 4.0 * n + 4.0), 0.0))
                    .subtract(1.0);
            Complex p = numeratorP.divide(denominatorP).multiply(sqrtPart);

            Complex numeratorA = n2.subtract(2.0);
            Complex sqrtPartA = p.multiply(p)
                    .add(p.multiply(new Complex(n, 0.0)).multiply(n2.subtract(2.0)).divide(new Complex(n - 1.0, 0.0)))
                    .sqrt();
            Complex denominatorA = p.multiply(new Complex(1.0 - n2.getReal(), 0.0)).add(sqrtPartA);
            Complex a = numeratorA.divide(denominatorA);

            Complex numeratorL1 = (new Complex(3.0 + a.getReal(), 0.0)).multiply(new Complex(n - 1.0, 0.0)).add(a.multiply(2.0));
            Complex denominatorL1 = new Complex((n - 1.0) * (1.0 + a.multiply(p).getReal()), 0.0);
            Complex numeratorL2 = a.multiply(new Complex(2.0 * (n + 1.0), 0.0));
            Complex denominatorL2 = new Complex(
                    2.0 * (n - 1.0)
                            + a.multiply(p).multiply(new Complex(n * a.getReal() + 2.0 * n - 2.0, 0.0)).getReal(),
                    0.0);
            Complex l = numeratorL1.divide(denominatorL1).subtract(numeratorL2.divide(denominatorL2));

            return l.getReal() * m1 * m2;
        } else {
            return (n + 1.0) / n * ni2 * m1 * m2;
        }
    }
}
