package jline.api.mam;

/**
 * Quadratic Programming result containing solution and status.
 */
public final class QPResult {
    public final double[] solution;
    public final double objective;
    public final boolean success;

    public QPResult(double[] solution, double objective, boolean success) {
        this.solution = solution;
        this.objective = objective;
        this.success = success;
    }

    public double[] getSolution() { return solution; }
    public double getObjective() { return objective; }
    public boolean isSuccess() { return success; }
}
