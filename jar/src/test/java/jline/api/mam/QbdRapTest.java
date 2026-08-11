package jline.api.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Validates {@link Qbd_rap}, the general QBD with Rational Arrival Process
 * components, against the published benchmark of N. G. Bean and B. F. Nielsen,
 * "Quasi-Birth-and-Death Processes with Rational Arrival Process Components",
 * Stochastic Models, 26(3), 2010, pp. 309-334: the marginal level distribution
 * of Table 1, the closed form of R, the eigenvalues of R and the stability
 * boundary gamma = 1/2. Also covers the general (non rank-one) path for G by
 * cross-checking a MAP/MAP/1 queue against {@link Qbd_raprap1}, and adds the
 * PH-in-ME-clothing identity for {@link Qbd_raprap1} itself.
 */
public class QbdRapTest {

    /** Marginal level distribution of Table 1 of Bean and Nielsen (2010), gamma = 0.25. */
    static final double[] TABLE1 = {0.6736, 0.2175, 0.0726, 0.0242, 0.0081, 0.0027, 0.0009, 0.0003, 0.0001};

    /** Local block A1 = Ca = Cs of the example process. */
    static Matrix exampleA1() {
        Matrix T = new Matrix(3, 3);
        T.set(0, 0, -1.0);   T.set(0, 1, 0.0);    T.set(0, 2, 0.0);
        T.set(1, 0, -2.0/3); T.set(1, 1, -1.0);   T.set(1, 2, 1.0);
        T.set(2, 0, 2.0/3);  T.set(2, 1, -1.0);   T.set(2, 2, -1.0);
        return T;
    }

    /** Jump matrix Da of the arrival RAP of the example process. */
    static Matrix exampleDa() {
        Matrix Da = new Matrix(3, 3);
        Da.set(0, 0, 14.0/5);  Da.set(0, 1, -9.0/10);  Da.set(0, 2, -9.0/10);
        Da.set(1, 0, 26.0/15); Da.set(1, 1, -8.0/15);  Da.set(1, 2, -8.0/15);
        Da.set(2, 0, 58.0/15); Da.set(2, 1, -19.0/15); Da.set(2, 2, -19.0/15);
        return Da;
    }

    /** Jump matrix Ds of the service RAP, the rank-one product (1,2/3,4/3)'*(3,-1,-1). */
    static Matrix exampleDs() {
        double[] u = {1.0, 2.0/3, 4.0/3};
        double[] v = {3.0, -1.0, -1.0};
        Matrix Ds = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Ds.set(i, j, u[i] * v[j]);
            }
        }
        return Ds;
    }

    /** Closed form of R given in Section 5 of Bean and Nielsen (2010). */
    static Matrix exampleR(double g) {
        Matrix R = new Matrix(3, 3);
        R.set(0, 0, -1.0/15*g*(33+2*g)/(g-1)); R.set(0, 1, 0.0); R.set(0, 2, 1.0/10*g*(g+9)/(g-1));
        R.set(1, 0, -2.0/45*g*(31+4*g)/(g-1)); R.set(1, 1, 0.0); R.set(1, 2, 2.0/15*g*(g+4)/(g-1));
        R.set(2, 0, -4.0/45*g*(34+g)/(g-1));   R.set(2, 1, 0.0); R.set(2, 2, 1.0/15*g*(g+19)/(g-1));
        return R;
    }

    static QbdRapResult solveExample(double g, int numLevels) {
        Matrix A1 = exampleA1();
        Matrix A0 = exampleDa().scale(g);
        Matrix A2 = exampleDs().scale(1.0 - g);
        return Qbd_rap.qbd_rap(A0, A1, A2, A0.copy(), A1.scale(g), numLevels);
    }

    @Test
    public void testTable1LevelDistribution() {
        QbdRapResult res = solveExample(0.25, 8);
        Matrix lp = res.getLevelProb();
        for (int n = 0; n <= 8; n++) {
            assertEquals(TABLE1[n], lp.get(0, n), 5e-5,
                    "Table 1 of Bean and Nielsen (2010), level " + n);
        }
        // The reported levels 0..8 carry essentially all of the mass.
        double tot = 0.0;
        for (int n = 0; n <= 8; n++) {
            tot += lp.get(0, n);
        }
        assertEquals(1.0, tot, 1e-4);
    }

    @Test
    public void testRClosedFormAndEigenvalues() {
        double g = 0.25;
        QbdRapResult res = solveExample(g, 8);
        Matrix R = res.getR();
        Matrix Rp = exampleR(g);
        assertEquals(0.0, R.sub(Rp).normFrobenius(), 1e-12, "R against the closed form of the paper");
        // Eigenvalues of R are (0, -gamma/15, gamma/(1-gamma)); the dominant one
        // is Sp(R) = gamma/(1-gamma).
        assertEquals(g / (1.0 - g), res.getSpr(), 1e-12, "Sp(R)");
    }

    @Test
    public void testGIsRankOneClosedFormAndSolvesTheQuadratic() {
        double g = 0.25;
        QbdRapResult res = solveExample(g, 4);
        Matrix A1 = exampleA1();
        Matrix A0 = exampleDa().scale(g);
        Matrix A2 = exampleDs().scale(1.0 - g);
        Matrix G = res.getG();
        assertEquals(0.0, A0.mult(G).mult(G).add(A1.mult(G)).add(A2).normFrobenius(), 1e-12,
                "residual of A0*G^2 + A1*G + A2");
        // Rank-one A2 gives G = e*v/(v*e) with v = (3,-1,-1) and v*e = 1.
        for (int i = 0; i < 3; i++) {
            assertEquals(3.0, G.get(i, 0), 1e-12);
            assertEquals(-1.0, G.get(i, 1), 1e-12);
            assertEquals(-1.0, G.get(i, 2), 1e-12);
        }
        // U = A1 + A0*G.
        assertEquals(0.0, res.getU().sub(A1.add(A0.mult(G))).normFrobenius(), 1e-12);
    }

    @Test
    public void testStabilityBoundaryAtGammaOneHalf() {
        // The queue is stable exactly when gamma < 1/2 (Corollary 8).
        assertTrue(solveExample(0.49, 4).getQN() > 0.0);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                solveExample(0.5, 4);
            }
        });
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                solveExample(0.55, 4);
            }
        });
    }

    @Test
    public void testMeanQueueLengthIncreasesWithGamma() {
        // Figure 1 of the paper: the mean queue length grows monotonically to
        // infinity as gamma approaches 1/2 from below.
        double prev = -1.0;
        double[] gammas = {0.05, 0.15, 0.25, 0.35, 0.45, 0.49};
        for (int i = 0; i < gammas.length; i++) {
            double qn = solveExample(gammas[i], 4).getQN();
            assertTrue(qn > prev, "mean queue length must increase with gamma");
            prev = qn;
        }
        assertTrue(prev > 20.0, "mean queue length must blow up near gamma = 1/2");
    }

    @Test
    public void testGeneralGPathAgainstMapMap1() {
        // A MAP/MAP/1 queue whose service jump matrix has full rank: A2 is not
        // rank one, so the general functional-iteration plus Newton path for G
        // is exercised. The oracle is Qbd_mapmap1, which reaches the same
        // answer through cyclic reduction and shares no code with Qbd_rap.
        // Qbd_raprap1 would NOT be an independent oracle here: it now delegates
        // to Qbd_rap itself.
        Matrix C0 = new Matrix(1, 1); C0.set(0, 0, -0.5);
        Matrix C1 = new Matrix(1, 1); C1.set(0, 0, 0.5);
        Matrix S0 = new Matrix(2, 2);
        S0.set(0, 0, -2.0); S0.set(0, 1, 0.1); S0.set(1, 0, 0.2); S0.set(1, 1, -3.0);
        Matrix S1 = new Matrix(2, 2);
        S1.set(0, 0, 1.9); S1.set(0, 1, 0.0); S1.set(1, 0, 0.0); S1.set(1, 1, 2.8);
        assertEquals(2, S1.rank(), "the service jump matrix must have full rank for this test");

        Matrix I2 = Matrix.eye(2);
        Matrix A0 = C1.kron(I2);
        Matrix A1 = C0.kron(I2).add(Matrix.eye(1).kron(S0));
        Matrix A2 = Matrix.eye(1).kron(S1);
        Matrix B1 = C0.kron(I2);
        QbdRapResult res = Qbd_rap.qbd_rap(A0, A1, A2, A0.copy(), B1, 400);

        QbdMapMap1Result ref = Qbd_mapmap1.qbd_mapmap1(new MatrixCell(C0, C1), new MatrixCell(S0, S1));
        assertEquals(ref.getQN(), res.getQN(), 1e-8, "mean queue length against qbd_mapmap1");
        double lvl0 = 0.0;
        for (int j = 0; j < 2; j++) {
            lvl0 += ref.getPqueue().get(0, j);
        }
        assertEquals(lvl0, res.getLevelProb().get(0, 0), 1e-9, "empty-system probability");
        // The G found by Newton must solve the quadratic and satisfy G*e = e.
        Matrix G = res.getG();
        assertEquals(0.0, A0.mult(G).mult(G).add(A1.mult(G)).add(A2).normFrobenius(), 1e-10);
        Matrix e = Matrix.ones(2, 1);
        assertEquals(0.0, G.mult(e).sub(e).normFrobenius(), 1e-10);
    }

    @Test
    public void testNonConservativeBlocksAreRejected() {
        Matrix A1 = exampleA1();
        Matrix A0 = exampleDa().scale(0.25);
        Matrix A2 = exampleDs().scale(0.75);
        final Matrix badA2 = A2.add(Matrix.eye(3));
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Qbd_rap.qbd_rap(A0, A1, badA2);
            }
        });
    }

    /** Builds the ME process pair (A, (-A e) alpha) from a representation. */
    static MatrixCell meProcess(double[] alpha, double[][] A) {
        int n = alpha.length;
        Matrix H0 = new Matrix(n, n);
        Matrix H1 = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                H0.set(i, j, A[i][j]);
                rowSum += A[i][j];
            }
            for (int j = 0; j < n; j++) {
                H1.set(i, j, -rowSum * alpha[j]);
            }
        }
        return new MatrixCell(H0, H1);
    }

    static MatrixCell expProcess(double rate) {
        Matrix H0 = new Matrix(1, 1); H0.set(0, 0, -rate);
        Matrix H1 = new Matrix(1, 1); H1.set(0, 0, rate);
        return new MatrixCell(H0, H1);
    }

    /** Genuine (non phase-type) ME: poles -0.5 and -1 +- 2*pi*i. */
    static MatrixCell genuineME() {
        double w = 2 * Math.PI;
        double[] alpha = {0.984694494294579, -0.040430911430916, 0.0557364171363366};
        double[][] A = {{-0.5, 0.0, 0.0}, {0.0, -1.0, w}, {0.0, -w, -1.0}};
        return meProcess(alpha, A);
    }

    /**
     * Production regression for {@link Qbd_raprap1} after its internals were
     * replaced by the Theorem 7 construction of {@link Qbd_rap}. These are the
     * four models the MAM ME/RAP path is validated on; the values are the ones
     * the cyclic-reduction implementation produced and must not move. They are
     * the truncated level series, not the exact closed form, which is why
     * qbd_raprap1 keeps its own truncation rule.
     */
    @Test
    public void testRapRap1ProductionRegression() {
        MatrixCell erlangAsME = meProcess(new double[]{1.0, 0.0}, new double[][]{{-2.0, 2.0}, {0.0, -2.0}});
        MatrixCell hyperAsME = meProcess(new double[]{0.6, 0.4}, new double[][]{{-2.0, 0.0}, {0.0, -0.5}});
        MatrixCell gen = genuineME();

        QbdRapRap1Result a = Qbd_raprap1.qbd_raprap1(expProcess(0.5), erlangAsME);
        assertEquals(0.874999998681, a.getQN(), 1e-9, "(a) M/ME/1 with Erlang-clothed ME");
        assertEquals(27, a.getPqueue().getNumRows(), "(a) truncated level count");

        QbdRapRap1Result b = Qbd_raprap1.qbd_raprap1(expProcess(0.5), hyperAsME);
        assertEquals(1.522222217439, b.getQN(), 1e-9, "(b) M/ME/1 with HyperExp-clothed ME");
        assertEquals(55, b.getPqueue().getNumRows(), "(b) truncated level count");

        QbdRapRap1Result c = Qbd_raprap1.qbd_raprap1(expProcess(0.255775446238906), gen);
        assertEquals(1.015214675095, c.getQN(), 1e-9, "(c) M/ME/1 with a genuine non-PH ME");
        assertEquals(34, c.getPqueue().getNumRows(), "(c) truncated level count");

        QbdRapRap1Result d = Qbd_raprap1.qbd_raprap1(gen, expProcess(1.0 / 0.977420));
        assertEquals(1.037045549083, d.getQN(), 1e-9, "(d) ME/M/1 with a genuine non-PH ME arrival");
        assertEquals(35, d.getPqueue().getNumRows(), "(d) truncated level count");
    }

    /**
     * Utilization of a single-server queue is exactly rho = lambda*E[S],
     * whatever the correlation structure. This is a free exact oracle and it is
     * what exposed the pre-rewrite defect: cyclic reduction silently returned a
     * G with residual 1.19 and ||G*e-e|| = 0.39 for a multi-phase RAP arrival,
     * giving UN = 0.4040 against the true 0.5000000733.
     */
    @Test
    public void testRapRap1UtilizationEqualsRho() {
        MatrixCell gen = genuineME();
        MatrixCell svc = expProcess(1.0 / 0.977420);
        QbdRapRap1Result r = Qbd_raprap1.qbd_raprap1(gen, svc);
        double rho = Map_lambda.map_lambda(gen) * Map_mean.map_mean(svc);
        assertEquals(rho, r.getUN(), 1e-12, "utilization must equal lambda*E[S]");
        // And the G that produced it must genuinely solve the quadratic.
        Matrix G = r.getG();
        assertEquals(0.0, r.getF().mult(G).mult(G).add(r.getL().mult(G)).add(r.getB()).normFrobenius(),
                1e-10, "G must solve F*G^2 + L*G + B = 0");
        Matrix e = Matrix.ones(G.getNumRows(), 1);
        assertEquals(0.0, G.mult(e).sub(e).normFrobenius(), 1e-10, "G*e must equal e");
    }

    /**
     * PH-in-ME-clothing identity for {@link Qbd_raprap1}: an M/ME/1 queue whose
     * matrix-exponential service is a similarity transform of an Erlang-2, and
     * therefore not a nonnegative representation, must give the same answer as
     * the M/Erlang/1 queue solved through the phase-type path of
     * {@link Qbd_mapmap1}.
     */
    @Test
    public void testRapRap1MatchesMapMap1ForPhInMeClothing() {
        double lambda = 0.4;
        double mu = 2.0;
        Matrix C0 = new Matrix(1, 1); C0.set(0, 0, -lambda);
        Matrix C1 = new Matrix(1, 1); C1.set(0, 0, lambda);
        // Erlang-2 as a phase-type MAP.
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -mu); D0.set(0, 1, mu); D0.set(1, 0, 0.0); D0.set(1, 1, -mu);
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 0.0); D1.set(0, 1, 0.0); D1.set(1, 0, mu); D1.set(1, 1, 0.0);
        // Similarity transform with S*e = e: preserves the process but destroys
        // nonnegativity, so (H0,H1) is a genuine ME representation.
        Matrix S = new Matrix(2, 2);
        S.set(0, 0, 1.0); S.set(0, 1, 0.0); S.set(1, 0, -0.5); S.set(1, 1, 1.5);
        Matrix Sinv = S.inv();
        Matrix H0 = Sinv.mult(D0).mult(S);
        Matrix H1 = Sinv.mult(D1).mult(S);
        boolean hasNegative = false;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if (i != j && H0.get(i, j) < 0.0) {
                    hasNegative = true;
                }
                if (H1.get(i, j) < 0.0) {
                    hasNegative = true;
                }
            }
        }
        assertTrue(hasNegative,
                "the transformed representation must have a negative entry, otherwise it is still a MAP");

        QbdMapMap1Result ph = Qbd_mapmap1.qbd_mapmap1(new MatrixCell(C0, C1), new MatrixCell(D0, D1));
        QbdRapRap1Result me = Qbd_raprap1.qbd_raprap1(new MatrixCell(C0, C1), new MatrixCell(H0, H1));

        assertEquals(ph.getXN(), me.getXN(), 1e-6, "throughput");
        assertEquals(ph.getUN(), me.getUN(), 1e-6, "utilization");
        assertEquals(ph.getQN(), me.getQN(), 1e-6, "mean queue length");
    }
}
