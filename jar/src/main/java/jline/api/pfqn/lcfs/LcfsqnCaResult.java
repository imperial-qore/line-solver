/**
 * @file Result class for LCFS convolution algorithm
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.lcfs;

public final class LcfsqnCaResult {
    private final double G;
    private final double V;

    public LcfsqnCaResult(double G, double V) {
        this.G = G;
        this.V = V;
    }

    public double getG() { return G; }
    public double getV() { return V; }
}
