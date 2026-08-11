/**
 * @file Quasi-Birth-Death process setup delays and server switch-off analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;

import jline.lang.processes.Coxian;
import jline.lang.processes.Exp;
import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_pi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_setupdelayoff {
    private Qbd_setupdelayoff() {}

    /**
     * Canonical PH representation of a phase given its rate, entered at phase 1
     * for every SCV. Exponential phases are built from the rate directly, which
     * is exact and avoids the mean-based FineTol cutoff in the fitters: an
     * Immediate setup or delay-off has rate GlobalConstants.Immediate whose mean
     * is exactly FineTol, so the round trip turned a finite 1e8 rate into an
     * infinite one. Mirrors coxian_phase in MATLAB qbd_setupdelayoff.m.
     *
     * <p>APH.fitMeanAndSCV must not be used here: the chain below overloads its
     * phase indices, so a level-up transition cannot also redistribute the
     * phase and an arrival to an off server has to enter the setup at phase 1.
     * APH violates that for SCV &gt; 1, where it returns a hyperexponential
     * entered at phase 2 with probability 3/4, which the chain then silently
     * entered at phase 1 regardless.
     *
     * @param rate the phase rate
     * @param scv the squared coefficient of variation of the phase
     * @return the (D0, D1) representation of the phase
     */
    private static MatrixCell coxianPhase(double rate, double scv) {
        if (scv == 1.0) {
            return new Exp(rate).getProcess();
        }
        return Coxian.fitMeanAndSCV(1.0 / rate, scv).getProcess();
    }

    /** Completion rate of each phase: t(i) = -sum_j D0(i,j). */
    private static Matrix completionRates(Matrix D0) {
        int n = D0.getNumRows();
        Matrix t = Matrix.zeros(n, 1);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += D0.get(i, j);
            }
            t.set(i, 0, -rowSum);
        }
        return t;
    }

    /**
     * Analyze a queue with setup delays and server switch-off using QBD approach.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param alphaRate rate of the setup phase
     * @param alphaScv squared coefficient of variation of the setup phase
     * @param betaRate rate of the delay-off phase
     * @param betaScv squared coefficient of variation of the delay-off phase
     * @return the mean queue length
     */
    public static double qbd_setupdelayoff(double lambda, double mu, double alphaRate, double alphaScv,
                                           double betaRate, double betaScv) {
        MatrixCell alpha = coxianPhase(alphaRate, alphaScv);
        Matrix alphaD0 = alpha.get(0);
        int na = alphaD0.getNumRows();
        // Completion rate of each phase. A Coxian may complete from any phase,
        // not only from the last one, so these are read per phase rather than
        // off the diagonal.
        Matrix ta = completionRates(alphaD0);

        MatrixCell beta = coxianPhase(betaRate, betaScv);
        Matrix betaD0 = beta.get(0);
        int nb = betaD0.getNumRows();
        Matrix tb = completionRates(betaD0);

        int n = na + nb;

        Matrix F = Matrix.zeros(n, n); // forward transitions
        for (int i = 0; i < na; i++) {
            F.set(i, i, lambda);
        }
        for (int i = 0; i < nb; i++) {
            F.set(na + i, na, lambda);
        }
        F.set(na, na, lambda);

        Matrix B = Matrix.zeros(n, n); // backward transitions
        B.set(na, na, mu);

        Matrix L = Matrix.zeros(n, n); // local transitions
        for (int i = 0; i < na; i++) {
            // Whole generator row, so a phase that both advances and completes
            // (any Coxian with SCV > 1) contributes both; reading only the
            // diagonal and the strict upper triangle assumed a pure series,
            // which holds only for SCV <= 1.
            for (int j = 0; j < na; j++) {
                L.set(i, j, alphaD0.get(i, j));
            }
            L.set(i, i, L.get(i, i) - lambda);
            L.set(i, na, ta.get(i, 0)); // setup completes from phase i -> busy server
        }
        L.set(na, na, -mu - lambda);
        for (int i = 1; i < nb; i++) {
            L.set(na + i, na + i, -lambda);
        }

        Matrix L0 = Matrix.zeros(n, n); // local transitions at the boundary level
        for (int i = 0; i < na; i++) {
            L0.set(i, i, -lambda);
        }
        for (int i = 0; i < nb; i++) {
            // As above, the whole generator row: the delay-off may expire from
            // any phase, and its phase-to-phase rate is the generator entry,
            // not the diagonal (they coincide only for a series phase).
            for (int j = 0; j < nb; j++) {
                L0.set(na + i, na + j, betaD0.get(i, j));
            }
            L0.set(na + i, na + i, L0.get(na + i, na + i) - lambda);
            L0.set(na + i, 0, tb.get(i, 0)); // delay-off expires from phase i -> server off
        }

        Map<String, Matrix> qbdResult = QBD_CR.QBD_CR(B, L, F, null, null, null, null);
        Matrix R = qbdResult.get("R");
        if (R == null) {
            throw new RuntimeException("QBD_CR failed to compute R matrix");
        }

        Matrix pn = QBD_pi.QBD_pi(B, L0, R);

        // Mean queue length: QBD_pi returns the level probabilities as a flat
        // vector of n phases per level, so level ni occupies
        // pn(ni*n .. ni*n+n-1), exactly n entries. Summing n+1 of them while
        // advancing by n let each window reach into the next level, so one phase
        // per level was counted twice under two different weights, and the
        // terminating check then dropped the last level: the queue length came
        // out high by up to 15%, worse the slower the setup, which is where the
        // overlapped phases hold most mass. Level 0 is skipped, holding no jobs.
        double QN = 0.0;
        int pnLength = pn.getNumCols();
        int j = n;
        int ni = 0;
        while (j + n <= pnLength) {
            ni++;
            double levelSum = 0.0;
            for (int k = j; k < j + n; k++) {
                levelSum += pn.get(0, k);
            }
            QN += ni * levelSum;
            j += n;
        }

        return QN;
    }
}
