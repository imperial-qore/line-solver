/**
 * @file Censored GE/GE/c/K;N queue solution
 *
 * Result of the ME censored queue building block.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import java.util.Arrays;

/**
 * Solution of a censored GE/GE/c/K;N queue by entropy maximisation.
 */
public final class MeGegecnResult {
    private final double[] p;
    private final int K;
    private final int N;
    private final double L;
    private final double U;
    private final double PB;
    private final double Lq;

    public MeGegecnResult(double[] p, int K, int N, double L, double U, double PB, double Lq) {
        this.p = p;
        this.K = K;
        this.N = N;
        this.L = L;
        this.U = U;
        this.PB = PB;
        this.Lq = Lq;
    }

    /** Queue length distribution, p[idx] = Pr{n = K+idx}. */
    public double[] getP() { return p; }
    /** Minimum number of jobs in the queue. */
    public int getK() { return K; }
    /** Buffer capacity in jobs. */
    public int getN() { return N; }
    /** Mean number of jobs in the queue. */
    public double getL() { return L; }
    /** Utilization, E[min(n,c)]/c. */
    public double getU() { return U; }
    /** Probability that an arrival of the queue's own stream finds it full. */
    public double getPB() { return PB; }
    /** Mean number of jobs waiting. */
    public double getLq() { return Lq; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MeGegecnResult)) return false;
        MeGegecnResult that = (MeGegecnResult) o;
        return K == that.K && N == that.N
                && Double.compare(L, that.L) == 0
                && Double.compare(U, that.U) == 0
                && Double.compare(PB, that.PB) == 0
                && Double.compare(Lq, that.Lq) == 0
                && Arrays.equals(p, that.p);
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(p);
        result = 31 * result + K;
        result = 31 * result + N;
        long bits = Double.doubleToLongBits(L);
        result = 31 * result + (int) (bits ^ (bits >>> 32));
        return result;
    }
}
