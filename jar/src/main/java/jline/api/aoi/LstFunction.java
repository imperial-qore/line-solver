/**
 * @file Laplace-Stieltjes Transform functional interface
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

/**
 * Functional interface for Laplace-Stieltjes Transform evaluation.
 *
 * Evaluates the LST at a given point s in the complex plane.
 */
public interface LstFunction {
    /**
     * Evaluate the LST at the given value of s.
     *
     * @param s The point at which to evaluate the LST (real, non-negative)
     * @return The value of the LST at s
     */
    double evaluate(double s);
}
