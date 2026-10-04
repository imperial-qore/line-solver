package jline.opt.solver.de;

import java.util.ArrayList;
import java.util.List;

/**
 * Self-contained port of scipy's {@code differential_evolution} that reproduces
 * its trajectory bit-for-bit for a given integer seed, using
 * {@link NumpyRandomState} for all draws.
 *
 * <p>Covers the configuration used by line-opt: the binomial strategies
 * (default {@code best1bin}), {@code updating='immediate'}, dithered mutation,
 * latinhypercube initialization, {@code polish=false}, {@code workers=1}, and
 * penalty-based constraints handled inside the objective (no scipy-level
 * constraints). The exact numpy draw order is preserved: LHS init
 * ({@code uniform} grid then per-column {@code permutation}), a per-generation
 * dither {@code uniform}, and per candidate {@code randint} (fill point),
 * {@code shuffle} (sample selection) and {@code uniform} (crossover), plus a
 * data-dependent {@code uniform} in the bound-repair step.</p>
 *
 * <p>Java 8 compatible.</p>
 */
public class DifferentialEvolution {

    /** Objective returning the (penalized) energy; may be {@code +inf}. */
    public interface Objective {
        double evaluate(double[] x);
    }

    /** Per-generation callback; return true to stop early. */
    public interface Callback {
        boolean call(double[] bestX, int nit);
    }

    public static class Result {
        public double[] x;
        public double fun;
        public int nit;
        public int nfev;
        public boolean success;
        public List<double[]> bestPerGen = new ArrayList<double[]>();
    }

    private final Objective objective;
    private final double[] low;
    private final double[] high;
    private final int paramCount;

    private final String strategy;
    private final int popsizeMult;
    private final int maxiter;
    private final double ditherLow;
    private final double ditherHigh;
    private final double recombination;
    private final double tol;
    private final double atol;
    private Callback callback;

    private final NumpyRandomState rng;

    private int numMembers;
    private double[][] population;      // [member][param], scaled to [0,1]
    private double[] energies;
    private int[] randomIndex;          // persistent shuffle buffer
    private double scale;
    private int nfev;

    private static final double MACHEPS = Math.ulp(1.0);

    public DifferentialEvolution(Objective objective, double[] low, double[] high,
                                 String strategy, int popsizeMult, int maxiter,
                                 double ditherLow, double ditherHigh,
                                 double recombination, double tol, long seed) {
        this.objective = objective;
        this.low = low.clone();
        this.high = high.clone();
        this.paramCount = low.length;
        this.strategy = strategy;
        this.popsizeMult = popsizeMult;
        this.maxiter = maxiter;
        this.ditherLow = Math.min(ditherLow, ditherHigh);
        this.ditherHigh = Math.max(ditherLow, ditherHigh);
        this.recombination = recombination;
        this.tol = tol;
        this.atol = 0.0;
        this.scale = ditherLow;
        this.rng = new NumpyRandomState(seed);
        initPopulationLhs();
    }

    public void setCallback(Callback cb) {
        this.callback = cb;
    }

    public NumpyRandomState getRng() {
        return rng;
    }

    // ---- scaling between population [0,1] and parameter space --------------

    private double scaleArg1(int j) {
        return 0.5 * (low[j] + high[j]);
    }

    private double scaleArg2(int j) {
        return Math.abs(low[j] - high[j]);
    }

    private double[] scaleParameters(double[] trial) {
        double[] out = new double[paramCount];
        for (int j = 0; j < paramCount; j++) {
            out[j] = scaleArg1(j) + (trial[j] - 0.5) * scaleArg2(j);
        }
        return out;
    }

    // ---- initialization ----------------------------------------------------

    private void initPopulationLhs() {
        int ebCount = 0;
        for (int j = 0; j < paramCount; j++) {
            if (low[j] == high[j]) {
                ebCount++;
            }
        }
        numMembers = Math.max(5, popsizeMult * Math.max(1, paramCount - ebCount));

        double segsize = 1.0 / numMembers;
        // samples = segsize * uniform(size=(numMembers, paramCount))
        //           + linspace(0,1,numMembers,endpoint=False)[:,None]
        double[][] samples = new double[numMembers][paramCount];
        double[] flat = rng.uniform(numMembers * paramCount);
        int p = 0;
        for (int i = 0; i < numMembers; i++) {
            double offset = ((double) i) / numMembers;
            for (int j = 0; j < paramCount; j++) {
                samples[i][j] = segsize * flat[p++] + offset;
            }
        }
        population = new double[numMembers][paramCount];
        for (int j = 0; j < paramCount; j++) {
            int[] order = rng.permutation(numMembers);
            for (int i = 0; i < numMembers; i++) {
                population[i][j] = samples[order[i]][j];
            }
        }
        energies = new double[numMembers];
        for (int i = 0; i < numMembers; i++) {
            energies[i] = Double.POSITIVE_INFINITY;
        }
        randomIndex = new int[numMembers];
        for (int i = 0; i < numMembers; i++) {
            randomIndex[i] = i;
        }
        nfev = 0;
    }

    // ---- helpers -----------------------------------------------------------

    private double[] best() {
        return scaleParameters(population[0]);
    }

    private void calculateInitialEnergies() {
        for (int i = 0; i < numMembers; i++) {
            energies[i] = objective.evaluate(scaleParameters(population[i]));
            nfev++;
        }
    }

    private void promoteLowestEnergy() {
        int l = 0;
        double best = energies[0];
        for (int i = 1; i < numMembers; i++) {
            if (energies[i] < best) {
                best = energies[i];
                l = i;
            }
        }
        if (l != 0) {
            double te = energies[0];
            energies[0] = energies[l];
            energies[l] = te;
            double[] tp = population[0];
            population[0] = population[l];
            population[l] = tp;
        }
    }

    private double std(double[] a) {
        double mean = 0.0;
        for (double v : a) {
            mean += v;
        }
        mean /= a.length;
        double var = 0.0;
        for (double v : a) {
            var += (v - mean) * (v - mean);
        }
        return Math.sqrt(var / a.length);
    }

    private double mean(double[] a) {
        double m = 0.0;
        for (double v : a) {
            m += v;
        }
        return m / a.length;
    }

    private boolean anyInf(double[] a) {
        for (double v : a) {
            if (Double.isInfinite(v)) {
                return true;
            }
        }
        return false;
    }

    private boolean converged() {
        if (anyInf(energies)) {
            return false;
        }
        return std(energies) <= atol + tol * Math.abs(mean(energies));
    }

    // ---- sample selection and mutation ------------------------------------

    private int[] selectSamples(int candidate, int numberSamples) {
        rng.shuffle(randomIndex);
        int[] pick = new int[numberSamples + 1];
        System.arraycopy(randomIndex, 0, pick, 0, numberSamples + 1);
        int[] out = new int[numberSamples];
        int k = 0;
        for (int i = 0; i < pick.length && k < numberSamples; i++) {
            if (pick[i] != candidate) {
                out[k++] = pick[i];
            }
        }
        return out;
    }

    private double[] bprime(int candidate, int[] s) {
        double[] b = new double[paramCount];
        if ("rand1bin".equals(strategy) || "rand1exp".equals(strategy)) {
            for (int j = 0; j < paramCount; j++) {
                b[j] = population[s[0]][j] + scale
                        * (population[s[1]][j] - population[s[2]][j]);
            }
        } else if ("randtobest1bin".equals(strategy) || "randtobest1exp".equals(strategy)) {
            for (int j = 0; j < paramCount; j++) {
                double v = population[s[0]][j];
                v += scale * (population[0][j] - v);
                v += scale * (population[s[1]][j] - population[s[2]][j]);
                b[j] = v;
            }
        } else if ("currenttobest1bin".equals(strategy) || "currenttobest1exp".equals(strategy)) {
            for (int j = 0; j < paramCount; j++) {
                b[j] = population[candidate][j] + scale
                        * (population[0][j] - population[candidate][j]
                           + population[s[0]][j] - population[s[1]][j]);
            }
        } else if ("best2bin".equals(strategy) || "best2exp".equals(strategy)) {
            for (int j = 0; j < paramCount; j++) {
                b[j] = population[0][j] + scale
                        * (population[s[0]][j] + population[s[1]][j]
                           - population[s[2]][j] - population[s[3]][j]);
            }
        } else if ("rand2bin".equals(strategy) || "rand2exp".equals(strategy)) {
            for (int j = 0; j < paramCount; j++) {
                b[j] = population[s[0]][j] + scale
                        * (population[s[1]][j] + population[s[2]][j]
                           - population[s[3]][j] - population[s[4]][j]);
            }
        } else {
            // default best1bin / best1exp
            for (int j = 0; j < paramCount; j++) {
                b[j] = population[0][j] + scale
                        * (population[s[0]][j] - population[s[1]][j]);
            }
        }
        return b;
    }

    private double[] mutate(int candidate) {
        int fillPoint = (int) rng.randint(0, paramCount);
        int[] samples = selectSamples(candidate, 5);
        double[] bp = bprime(candidate, samples);
        double[] trial = population[candidate].clone();
        double[] cross = rng.uniform(paramCount);
        boolean[] doCross = new boolean[paramCount];
        for (int j = 0; j < paramCount; j++) {
            doCross[j] = cross[j] < recombination;
        }
        // binomial: force the fill point
        doCross[fillPoint] = true;
        for (int j = 0; j < paramCount; j++) {
            if (doCross[j]) {
                trial[j] = bp[j];
            }
        }
        return trial;
    }

    private void ensureConstraint(double[] trial) {
        int oob = 0;
        for (int j = 0; j < paramCount; j++) {
            if (trial[j] > 1 || trial[j] < 0) {
                oob++;
            }
        }
        if (oob > 0) {
            double[] repl = rng.uniform(oob);
            int k = 0;
            for (int j = 0; j < paramCount; j++) {
                if (trial[j] > 1 || trial[j] < 0) {
                    trial[j] = repl[k++];
                }
            }
        }
    }

    // ---- main loop ---------------------------------------------------------

    public Result solve() {
        Result result = new Result();
        boolean warningFlag = false;

        if (anyInf(energies)) {
            calculateInitialEnergies();
            promoteLowestEnergy();
        }

        int nit = 0;
        for (nit = 1; nit <= maxiter; nit++) {
            next();

            if (callback != null) {
                double[] bx = best();
                boolean stop = callback.call(bx, nit);
                if (stop) {
                    warningFlag = true;
                }
            }
            result.bestPerGen.add(best());

            if (warningFlag || converged()) {
                break;
            }
        }
        if (nit > maxiter) {
            nit = maxiter;
            warningFlag = true;
        }

        result.x = best();
        result.fun = energies[0];
        result.nit = nit;
        result.nfev = nfev;
        result.success = !warningFlag;
        return result;
    }

    private void next() {
        // dither each generation (mutation is a tuple -> dither active)
        scale = rng.uniform(ditherLow, ditherHigh);

        for (int candidate = 0; candidate < numMembers; candidate++) {
            double[] trial = mutate(candidate);
            ensureConstraint(trial);
            double[] parameters = scaleParameters(trial);
            double energy = objective.evaluate(parameters);
            nfev++;
            if (energy <= energies[candidate]) {
                population[candidate] = trial;
                energies[candidate] = energy;
                if (energy <= energies[0]) {
                    promoteLowestEnergy();
                }
            }
        }
    }
}
