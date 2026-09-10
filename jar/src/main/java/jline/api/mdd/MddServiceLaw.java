package jline.api.mdd;

/**
 * Phase-type service law of one station, as a Markovian (D0,D1) pair.
 *
 * <p>D0 holds the phase transitions that do not complete a service and D1 those
 * that do, so the exit-rate vector is t0 = D1*1. For a renewal law D1 = t0*pie,
 * and the entry law pie is then DERIVED from D1 rather than assumed: see
 * {@link Mdd_descriptor#entryLaw}. Set {@link #pie} only to override that.</p>
 */
public class MddServiceLaw {

    /** Phase transitions without a completion. */
    public final double[][] D0;
    /** Phase transitions with a completion, D1 = t0*pie for a renewal law. */
    public final double[][] D1;
    /** Optional entry law; null to derive it from D1. */
    public final double[] pie;

    public MddServiceLaw(double[][] D0, double[][] D1) {
        this(D0, D1, null);
    }

    public MddServiceLaw(double[][] D0, double[][] D1, double[] pie) {
        this.D0 = D0;
        this.D1 = D1;
        this.pie = pie;
    }

    /** Number of phases. */
    public int phases() {
        return D0.length;
    }
}
