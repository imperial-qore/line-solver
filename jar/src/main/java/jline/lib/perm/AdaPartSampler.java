package jline.lib.perm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import jline.util.matrix.Matrix;

/**
 * Sampling method to approximate the permanent using Adaptive Partitioning (AdaPart).
 */
public class AdaPartSampler extends PermSolver {
    private final int maximumAcceptedSamples;
    private final long maximumTime; // milliseconds
    private final int maximumSamples;
    private final String mode;

    /** Draw budget of the classic mode; exceeding it throws rather than spins. */
    private static final int MAX_DRAWS = 1000000;

    private final Random random = new Random();

    public List<Integer> sampleAccepted = Collections.emptyList();
    public List<Long> sampleTime = Collections.emptyList();
    public List<Double> permStep = Collections.emptyList();

    public AdaPartSampler(Matrix matrix) {
        this(matrix, 100, 30000L, 450, "classic", false);
    }

    public AdaPartSampler(Matrix matrix, int maximumAcceptedSamples, long maximumTime,
                          int maximumSamples, String mode, boolean solve) {
        super(matrix);
        if (matrix.getNumRows() > 0) {
            PermSupport.requireFullSupport(matrix, "adapart");
        }
        this.maximumAcceptedSamples = maximumAcceptedSamples;
        this.maximumTime = maximumTime;
        this.maximumSamples = maximumSamples;
        this.mode = mode;
        if (solve) solve();
    }

    @Override
    public void compute() {
        if ("classic".equals(mode)) value = samplingClassic();
        else if ("time".equals(mode)) value = samplingTime();
        else if ("sample".equals(mode)) value = samplingSample();
        else value = samplingClassic();
    }

    private double samplingClassic() {
        double zUB = soulesBound(matrix);
        int accepted = 0;
        int total = 0;
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        // Bounded independently of the scaling: with perm(A) = 0 the acceptance
        // probability is 0 and this loop would never terminate. A cap that
        // RETURNS a number would be a workaround, so it throws.
        while (accepted < maximumAcceptedSamples) {
            if (total >= MAX_DRAWS) {
                throw new IllegalArgumentException("Only " + accepted + " of the "
                        + maximumAcceptedSamples + " required acceptances were obtained in "
                        + total + " draws. Use the exact engine.");
            }
            int sample = sample();
            accepted += sample;
            total++;
            acceptedList.add(sample);
            timeList.add(System.currentTimeMillis() - startTime);
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep(zUB);
        return zUB * accepted / (double) total;
    }

    private double samplingTime() {
        double zUB = soulesBound(matrix);
        int accepted = 0;
        int total = 0;
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        while (System.currentTimeMillis() - startTime < maximumTime) {
            int sample = sample();
            accepted += sample;
            total++;
            acceptedList.add(sample);
            timeList.add(System.currentTimeMillis() - startTime);
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep(zUB);
        return total > 0 ? zUB * accepted / (double) total : 0.0;
    }

    private double samplingSample() {
        double zUB = soulesBound(matrix);
        int accepted = 0;
        int total = 0;
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        while (acceptedList.size() < maximumSamples) {
            int sample = sample();
            accepted += sample;
            total++;
            acceptedList.add(sample);
            timeList.add(System.currentTimeMillis() - startTime);
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep(zUB);
        return zUB * accepted / (double) total;
    }

    private int sample() {
        Set<List<Integer>> S = new HashSet<List<Integer>>();
        List<Integer> initial = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) initial.add(n);
        S.add(initial);
        long startTime = System.currentTimeMillis();
        while (anyContainsN(S) && (!"time".equals(mode) || System.currentTimeMillis() - startTime < maximumTime)) {
            List<Integer> sInit = S.iterator().next();
            Matrix sMatrix = modifyMatrix(matrix, sInit);
            double ub = soulesBound(sMatrix);
            double zubS = ub;
            boolean init = true;
            while ((ub >= zubS || init) && (!"time".equals(mode) || System.currentTimeMillis() - startTime < maximumTime)) {
                init = false;
                // Only elements with a free column can be refined; expanding a
                // complete assignment yields no children and stalls the sampler.
                List<List<Integer>> sList = new ArrayList<List<Integer>>();
                for (List<Integer> cand : S) {
                    if (cand.contains(n)) sList.add(cand);
                }
                if (sList.isEmpty()) break;
                List<Integer> sSub = sList.get(random.nextInt(sList.size()));
                S.remove(sSub);
                Matrix subMatrix = modifyMatrix(matrix, sSub);
                // Discount the bound of the element actually removed, not the
                // root bound; see _kb/03-api-layer.md.
                double subUb = soulesBound(subMatrix);
                double[] selRes = selectColumn(subMatrix, subUb, ub, sSub);
                double newUb = selRes[0];
                int j = (int) selRes[1];
                for (int i = 0; i < n; i++) {
                    List<Integer> sAdd = new ArrayList<Integer>(sSub);
                    if (!sSub.contains(i)) {
                        sAdd.set(j, i);
                        S.add(sAdd);
                    }
                }
                boolean noProgress = newUb >= ub;
                ub = newUb;
                // The Soules bound is tight on matrices with equal entries, so
                // refinement cannot improve it and "refine until improved" would
                // never exit. Stop on the first non-improving expansion instead;
                // a tight bound means the draw is accepted with probability 1.
                if (noProgress) break;
            }
            int c = computeProbabilities(S, zubS);
            if (c == S.size()) return 0;
            S = subset(S, c);
        }
        return 1;
    }

    private boolean anyContainsN(Set<List<Integer>> S) {
        for (List<Integer> l : S) if (l.contains(n)) return true;
        return false;
    }

    /**
     * Picks the column whose expansion minimizes the summed Soules bound.
     *
     * Only columns still unassigned in sSub are candidates. Scoring an
     * already-assigned column just re-derives its own constraint, which always
     * looks cheapest, so the sampler would re-split the same column forever and
     * never complete an assignment.
     */
    private double[] selectColumn(Matrix sMatrix, double removedUb, double ub,
                                  List<Integer> sSub) {
        double[] ubi = new double[n];
        for (int i = 0; i < n; i++) {
            if (sSub.get(i) != n) {
                ubi[i] = Double.POSITIVE_INFINITY;
                continue;
            }
            for (int j = 0; j < n; j++) {
                List<Integer> L = new ArrayList<Integer>();
                for (int k = 0; k < n; k++) L.add(n);
                L.set(i, j);
                Matrix fB = modifyMatrix(sMatrix, L);
                ubi[i] += soulesBound(fB);
            }
        }
        int j = 0;
        double minVal = ubi[0];
        for (int i = 1; i < n; i++) {
            if (ubi[i] < minVal) { minVal = ubi[i]; j = i; }
        }
        double newUb = ub - removedUb + ubi[j];
        return new double[]{newUb, (double) j};
    }

    private int computeProbabilities(Set<List<Integer>> S, double zubS) {
        List<List<Integer>> sList = new ArrayList<List<Integer>>(S);
        List<Double> p = new ArrayList<Double>();
        for (List<Integer> s : sList) p.add(soulesBound(modifyMatrix(matrix, s)));
        double sum = 0.0;
        for (double v : p) sum += v;
        if (sum > 0) {
            for (int i = 0; i < p.size(); i++) p.set(i, p.get(i) / zubS);
        }
        double slack = 1.0;
        for (double v : p) slack -= v;
        p.add(slack);
        double totalProb = 0.0;
        for (double v : p) totalProb += v;
        if (totalProb > 0) {
            for (int i = 0; i < p.size(); i++) p.set(i, Math.abs(p.get(i)) / totalProb);
        }
        double rand = random.nextDouble();
        double cumSum = 0.0;
        for (int i = 0; i < p.size(); i++) {
            cumSum += p.get(i);
            if (rand <= cumSum) return i;
        }
        return p.size() - 1;
    }

    private Set<List<Integer>> subset(Set<List<Integer>> S, int c) {
        List<List<Integer>> sList = new ArrayList<List<Integer>>(S);
        List<Integer> sInter = new ArrayList<Integer>(sList.get(c));
        int missingCount = 0;
        for (int v : sInter) if (v == n) missingCount++;
        if (missingCount == 1) {
            Set<Integer> existingValues = new HashSet<Integer>();
            for (int v : sInter) if (v != n) existingValues.add(v);
            int missingValue = -1;
            for (int v = 0; v < n; v++) if (!existingValues.contains(v)) { missingValue = v; break; }
            int missingIndex = sInter.indexOf(n);
            sInter.set(missingIndex, missingValue);
        }
        Set<List<Integer>> result = new HashSet<List<Integer>>();
        result.add(sInter);
        return result;
    }

    private Matrix modifyMatrix(Matrix M, List<Integer> t) {
        Matrix result = Matrix.zeros(n, n);
        boolean[][] mask = new boolean[n][n];
        for (int j = 0; j < t.size(); j++) {
            if (t.get(j) != n) mask[t.get(j)][j] = true;
        }
        for (int i = 0; i < n; i++) {
            if (!t.contains(i)) {
                for (int j = 0; j < t.size(); j++) {
                    if (t.get(j) == n) mask[i][j] = true;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (mask[i][j]) result.set(i, j, M.get(i, j));
            }
        }
        return result;
    }

    private double soulesBound(Matrix M) {
        double[] gamma = new double[n + 1];
        gamma[0] = 0.0;
        for (int k = 1; k <= n; k++) {
            double factorial = 1.0;
            for (int i = 1; i <= k; i++) factorial *= i;
            gamma[k] = Math.pow(factorial, 1.0 / k);
        }
        double[] delta = new double[n];
        for (int i = 0; i < n; i++) delta[i] = gamma[n - i] - gamma[n - i - 1];
        double[] mSum = new double[n];
        for (int j = 0; j < n; j++) {
            double[] column = new double[n];
            for (int i = 0; i < n; i++) column[i] = M.get(i, j);
            Arrays.sort(column);
            for (int i = 0; i < n; i++) mSum[j] += column[i] * delta[i];
        }
        double product = 1.0;
        for (double sum : mSum) product *= sum;
        return product;
    }

    private void computePermStep(double zUB) {
        List<Double> steps = new ArrayList<Double>();
        int cumAccepted = 0;
        for (int i = 0; i < sampleAccepted.size(); i++) {
            cumAccepted += sampleAccepted.get(i);
            steps.add(zUB * cumAccepted / (i + 1));
        }
        permStep = steps;
    }
}
