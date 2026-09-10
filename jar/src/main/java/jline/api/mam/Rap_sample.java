package jline.api.mam;

import java.util.Random;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Rational Arrival Process (RAP) sampling by conditional inversion.
 *
 * <p>A RAP (H0, H1) generalises a MAP by dropping the sign constraints on the
 * two matrices, so the CTMC walk of {@link Map_sample} does not apply. Sampling
 * is instead performed on the conditional (vector) representation: with a row
 * vector v satisfying {@code v*e = 1}, the next inter-event time has survival
 * function {@code S_v(x) = v*expm(H0*x)*e} and density
 * {@code f_v(x) = -v*expm(H0*x)*H0*e = v*expm(H0*x)*H1*e}, and after an event at
 * x the vector is updated to</p>
 *
 * <pre>{@code v <- v*expm(H0*x)*H1 / (v*expm(H0*x)*H1*e)}</pre>
 *
 * <p>The generated sequence is therefore <em>correlated</em>: it reproduces the
 * autocorrelation of the process, not merely its ME marginal. Starting v from
 * the embedded stationary vector {@link Map_pie} makes each inter-event time
 * marginally stationary. A MAP is the special case in which H0 has non-negative
 * off-diagonal entries and H1 is non-negative; the algorithm is exact there too,
 * although {@link Map_sample.MapSampler} is cheaper for that case.</p>
 *
 * <p>Because v changes at every draw, no inversion table can be cached across
 * draws. What is cached is the sequence {@code expm(H0*t_i)} on a fixed time
 * grid (built with a single matrix exponential and repeated matrix products),
 * together with the column vectors {@code expm(H0*t_i)*e}. A draw then costs
 * only vector products: the conditional survival at a grid point is a dot
 * product against the cached column, so no matrix exponential is evaluated per
 * draw.</p>
 *
 * @since LINE 3.0
 */
public final class Rap_sample {
    private Rap_sample() {}

    /**
     * Stateful RAP sampler that carries the conditional phase vector across
     * draws, so successive inter-event times are correlated as in the process.
     */
    public static final class RapSampler {
        private final int n;
        private final Matrix d0Matrix;
        private final double[][] d0;
        private final double[][] d1;
        private final double[] rowSumD0;
        private final double normD0;
        private final double[] pie;
        private double[] v;
        private final int gridSize;
        private final double h;
        private final double[] tGrid;
        private final double[][][] eGrid;
        private final double[][] gGrid;
        private final double tEnd;
        private final double eta;
        private final boolean exponential;
        private final double expRate;

        /**
         * Builds the cached grid of matrix exponentials of H0 and initialises the
         * conditional vector to the embedded stationary vector.
         *
         * @param D0 the hidden-transition matrix H0
         * @param D1 the visible-transition matrix H1
         */
        public RapSampler(Matrix D0, Matrix D1) {
            this.n = D0.getNumRows();
            this.d0Matrix = D0;
            this.d0 = Me_sample.toArray2D(D0);
            this.d1 = Me_sample.toArray2D(D1);
            this.rowSumD0 = new double[n];
            double maxRowAbs = 0.0;
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                double abs = 0.0;
                for (int j = 0; j < n; j++) {
                    s += d0[i][j];
                    abs += Math.abs(d0[i][j]);
                }
                rowSumD0[i] = s;
                if (abs > maxRowAbs) maxRowAbs = abs;
            }
            this.normD0 = maxRowAbs;

            boolean isExp = (n == 1) && d0[0][0] < 0.0;
            this.exponential = isExp;
            this.expRate = isExp ? -d0[0][0] : 0.0;

            this.pie = new double[n];
            if (isExp) {
                this.pie[0] = 1.0;
                this.v = this.pie.clone();
                this.gridSize = 0;
                this.h = 0.0;
                this.tGrid = null;
                this.eGrid = null;
                this.gGrid = null;
                this.tEnd = 0.0;
                this.eta = -expRate;
                return;
            }

            Matrix pieMatrix = Map_pie.map_pie(D0, D1);
            for (int j = 0; j < n; j++) {
                this.pie[j] = pieMatrix.get(0, j);
            }
            this.v = this.pie.clone();

            double mean = Map_mean.map_mean(D0, D1);
            double var = Map_var.map_var(D0, D1);
            double sigma = (var > 0.0 && !Double.isNaN(var)) ? Math.sqrt(var) : 0.0;
            double horizon = mean + 10.0 * sigma;
            if (!(horizon > 0.0) || Double.isInfinite(horizon) || Double.isNaN(horizon)) {
                horizon = 1.0;
            }
            for (int k = 0; k < Me_sample.MAX_DOUBLINGS; k++) {
                if (maxRowSurvival(horizon) < Me_sample.TAIL_MASS_TOL) {
                    break;
                }
                horizon *= 2.0;
            }

            this.gridSize = Me_sample.GRID_POINTS;
            this.h = horizon / (gridSize - 1);
            this.tGrid = new double[gridSize];
            this.eGrid = new double[gridSize][][];
            this.gGrid = new double[gridSize][];

            double[][] eh = Me_sample.toArray2D(D0.scale(h).expm_higham());
            double[][] identity = new double[n][n];
            for (int i = 0; i < n; i++) {
                identity[i][i] = 1.0;
            }
            eGrid[0] = identity;
            tGrid[0] = 0.0;
            gGrid[0] = rowSums(identity);
            for (int i = 1; i < gridSize; i++) {
                eGrid[i] = matTimesMat(eGrid[i - 1], eh);
                tGrid[i] = i * h;
                gGrid[i] = rowSums(eGrid[i]);
            }
            this.tEnd = tGrid[gridSize - 1];
            this.eta = Me_sample.dominantRate(D0, mean);
        }

        /** Largest row sum of expm(H0*t), an upper bound on any conditional survival. */
        private double maxRowSurvival(double t) {
            double[][] e = Me_sample.toArray2D(d0Matrix.scale(t).expm_higham());
            double[] rs = rowSums(e);
            double m = 0.0;
            for (int i = 0; i < n; i++) {
                m = Math.max(m, Math.abs(rs[i]));
            }
            return m;
        }

        /**
         * Draws one inter-event time and advances the conditional phase vector.
         *
         * @param random the uniform random source
         * @return the next inter-event time
         */
        public double next(Random random) {
            double u = random.nextDouble();
            if (exponential) {
                return -Math.log(1.0 - u) / expRate;
            }
            double target = 1.0 - u;

            double x;
            int lo;
            double sLast = Me_sample.dot(v, gGrid[gridSize - 1]);
            if (sLast >= target) {
                lo = gridSize - 1;
                if (sLast > 0.0 && target > 0.0 && eta < 0.0) {
                    x = tEnd + Math.log(sLast / target) / (-eta);
                    if (!(x > tEnd)) {
                        x = tEnd;
                    }
                } else {
                    x = tEnd;
                }
            } else {
                lo = 0;
                int hi = gridSize - 1;
                while (hi - lo > 1) {
                    int mid = (lo + hi) >>> 1;
                    if (Me_sample.dot(v, gGrid[mid]) >= target) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                double sLo = Me_sample.dot(v, gGrid[lo]);
                double sHi = Me_sample.dot(v, gGrid[lo + 1]);
                double denom = sLo - sHi;
                x = (denom > 0.0) ? tGrid[lo] + (sLo - target) / denom * h : tGrid[lo];
                double left = tGrid[lo];
                double right = left + h;
                double[] wLo = Me_sample.vecTimesMat(v, eGrid[lo]);
                for (int k = 0; k < Me_sample.NEWTON_STEPS; k++) {
                    double[] w = Me_sample.expmPropagate(wLo, d0, d0Matrix, normD0, x - left);
                    double surv = Me_sample.sum(w);
                    double f = -Me_sample.dot(w, rowSumD0);
                    if (!(f > 0.0)) {
                        break;
                    }
                    // Newton on g(x) = S(x) - target, with g'(x) = S'(x) = -f:
                    // x <- x - g/g' = x + (S(x) - target)/f. When the survival at
                    // x still exceeds the target the root lies to the right, so
                    // the step must be positive.
                    double err = surv - target;
                    if (Math.abs(err) < 1e-14) {
                        break;
                    }
                    double xn = x + err / f;
                    if (!(xn > left) || !(xn < right)) {
                        break;
                    }
                    boolean converged = Math.abs(xn - x) <= 1e-15 * Math.max(1.0, Math.abs(x));
                    x = xn;
                    if (converged) {
                        break;
                    }
                }
            }

            double[] wStart = Me_sample.vecTimesMat(v, eGrid[lo]);
            double[] w = Me_sample.expmPropagate(wStart, d0, d0Matrix, normD0, x - tGrid[lo]);
            double[] vn = Me_sample.vecTimesMat(w, d1);
            double norm = Me_sample.sum(vn);
            if (norm > 1e-300 && !Double.isNaN(norm) && !Double.isInfinite(norm)) {
                for (int j = 0; j < n; j++) {
                    vn[j] /= norm;
                }
                v = vn;
            } else {
                v = pie.clone();
            }
            return x;
        }

        /** Resets the conditional vector to the embedded stationary vector. */
        public void reset() {
            this.v = this.pie.clone();
        }
    }

    /**
     * Generates a correlated sample sequence from a Rational Arrival Process.
     *
     * <p>The returned inter-event times are <em>not</em> i.i.d.: they carry the
     * autocorrelation of the process, obtained by propagating the conditional
     * phase vector across events. A MAP is the special case in which H0 has
     * non-negative off-diagonal entries and H1 is non-negative.</p>
     *
     * @param RAP    the process as a MatrixCell {H0, H1}
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of n correlated inter-event times
     */
    public static double[] rap_sample(MatrixCell RAP, long n, Random random) {
        return rap_sample(RAP.get(0), RAP.get(1), n, random);
    }

    /**
     * Generates a correlated sample sequence from a Rational Arrival Process
     * given its two matrices.
     *
     * @param H0     the hidden-transition matrix
     * @param H1     the visible-transition matrix
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of n correlated inter-event times
     */
    public static double[] rap_sample(Matrix H0, Matrix H1, long n, Random random) {
        RapSampler sampler = new RapSampler(H0, H1);
        double[] samples = new double[(int) n];
        for (int i = 0; i < (int) n; i++) {
            samples[i] = sampler.next(random);
        }
        return samples;
    }

    /** Row sums of a square matrix. */
    private static double[] rowSums(double[][] m) {
        double[] out = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            double s = 0.0;
            for (int j = 0; j < m[i].length; j++) {
                s += m[i][j];
            }
            out[i] = s;
        }
        return out;
    }

    /** Dense matrix product. */
    private static double[][] matTimesMat(double[][] x, double[][] y) {
        int rows = x.length;
        int inner = y.length;
        int cols = y[0].length;
        double[][] out = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            double[] xi = x[i];
            double[] oi = out[i];
            for (int k = 0; k < inner; k++) {
                double xik = xi[k];
                if (xik == 0.0) continue;
                double[] yk = y[k];
                for (int j = 0; j < cols; j++) {
                    oi[j] += xik * yk[j];
                }
            }
        }
        return out;
    }
}
