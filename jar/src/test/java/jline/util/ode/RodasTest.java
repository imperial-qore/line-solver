package jline.util.ode;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Conformance of the ported RODAS against the unmodified Hairer-Wanner Fortran.
 *
 * <p>THE ORACLE IS THE ORIGINAL, NOT A CLOSED FORM. {@link Rodas} is a
 * transliteration, so the property that matters is not that it solves these
 * problems well but that it solves them IDENTICALLY to the Fortran it came
 * from: a translation defect shows up as a last-digit drift long before it
 * shows up as a wrong answer. The expected values below were produced by
 *
 * <pre>gfortran -O2 -std=legacy rodas.f dc_decsol.f decsol.f</pre>
 *
 * <p>printed at 17 significant digits so the decimal round-trips through a
 * double. They are the SAME constants {@code cpp/tests/test_rodas.cpp} and
 * {@code python/tests/test_rodas.py} carry, so every port is pinned to one
 * oracle rather than to the others.
 *
 * <p>THE STEP COUNTERS ARE CHECKED EXACTLY, and are the more sensitive test:
 * nfcn, njac and nstep diverge on any change to the step-size controller or the
 * error estimate, while the final value can absorb a small perturbation and
 * still look right.
 *
 * <p>THE VALUES ARE CHECKED TO A TOLERANCE RATHER THAN BIT FOR BIT, and the
 * reason is one function call, measured rather than assumed. The step-size
 * controller evaluates {@code err**0.25}. C++ and CPython both reach glibc's
 * {@code pow}, which is correctly rounded, and both reproduce the Fortran
 * BIT FOR BIT on every case below. {@code Math.pow} on HotSpot is an intrinsic
 * accurate to under 1 ulp but NOT correctly rounded: over 200000 arguments
 * drawn from the range {@code err} actually takes, it differs from glibc's
 * {@code pow} in the last bit for 0.107% of them ({@code StrictMath.pow}, the
 * reproducible-by-spec alternative, differs for 9.8% and is therefore worse,
 * not better).
 *
 * <p>What that costs is bounded and was measured too. Sweeping the Van der Pol
 * case over six tolerances, this port reproduces the bit-exact ports' step
 * counters EXACTLY up to about 1100 steps -- see
 * {@link #testStiffOdeStepSequenceIsExactUntilThePowDriftAccumulates} -- and
 * only past that does one differing step size flip an accept/reject decision.
 * So the exactness claim is real and is asserted where it holds; on the two
 * longest runs the claim weakens to accuracy, and the tolerances below are the
 * measured drift rounded up, not a number chosen to make a test pass. A genuine
 * transliteration defect -- an off-by-one in a loop bound, a mistyped tableau
 * entry -- moves these answers in the third digit or destroys them outright,
 * orders of magnitude outside any of these bounds.
 *
 * <p>The fluid {@code dae} transient runs hundreds of steps, not thousands, so
 * it sits in the range where this port and the other three agree step for step.
 *
 * <p>The cases reach three different IJOB paths in DECOMR/SLVROD -- identity,
 * banded and full mass matrix -- crossed with analytic and numerical Jacobians,
 * because those are separate code paths and a transliteration error can hit one
 * and not the others. The singular mass matrix (index-1 DAE) is the case LINE
 * needs; the rest guard it against regressions elsewhere in the file.
 */
public class RodasTest {

    /** DAE cases: the drift has nowhere to accumulate over 47..252 steps. */
    private static final double REL = 1e-13;
    /** Van der Pol, analytic Jacobian: 5409 steps, measured drift 1.3e-13. */
    private static final double REL_LONG = 1e-11;
    /**
     * Van der Pol, numerical Jacobian: 5467 steps and the one case whose step
     * sequence parts company. Measured drift 5.7e-11, and still three orders
     * below the 1e-8 the run was integrated to, so both trajectories are within
     * the accuracy asked of them.
     */
    private static final double REL_LONG_NUMJAC = 1e-9;

    /**
     * Robertson's problem as an index-1 DAE: M y' = f(y), M = diag(1,1,0).
     * The third equation is the algebraic constraint y1+y2+y3 = 1.
     */
    private static final Rodas.Fcn FROB = new Rodas.Fcn() {
        public void eval(double x, double[] y, double[] f) {
            f[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            f[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
            f[2] = y[0] + y[1] + y[2] - 1.0;
        }
    };

    private static final Rodas.Jac JROB = new Rodas.Jac() {
        public void eval(double x, double[] y, double[][] dfy) {
            dfy[0][0] = -0.04;
            dfy[0][1] = 1.0e4 * y[2];
            dfy[0][2] = 1.0e4 * y[1];
            dfy[1][0] = 0.04;
            dfy[1][1] = -1.0e4 * y[2] - 6.0e7 * y[1];
            dfy[1][2] = -1.0e4 * y[1];
            dfy[2][0] = 1.0;
            dfy[2][1] = 1.0;
            dfy[2][2] = 1.0;
        }
    };

    /** Banded storage, MLMAS=MUMAS=0, i.e. the diagonal -- reaches IJOB=3. */
    private static final Rodas.Mas MDIAG3 = new Rodas.Mas() {
        public void eval(double[][] am) {
            am[0][0] = 1.0;
            am[0][1] = 1.0;
            am[0][2] = 0.0;
        }
    };

    /** The same mass matrix stored full -- reaches IJOB=5, and must agree. */
    private static final Rodas.Mas MFULL3 = new Rodas.Mas() {
        public void eval(double[][] am) {
            for (int i = 0; i < am.length; i++) {
                for (int j = 0; j < am[i].length; j++) {
                    am[i][j] = 0.0;
                }
            }
            am[0][0] = 1.0;
            am[1][1] = 1.0;
        }
    };

    /** Van der Pol at eps=1e-6: stiff, M = identity, reaches IJOB=1. */
    private static final Rodas.Fcn FVDP = new Rodas.Fcn() {
        public void eval(double x, double[] y, double[] f) {
            double eps = 1.0e-6;
            f[0] = y[1];
            f[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / eps;
        }
    };

    private static final Rodas.Jac JVDP = new Rodas.Jac() {
        public void eval(double x, double[] y, double[][] dfy) {
            double eps = 1.0e-6;
            dfy[0][0] = 0.0;
            dfy[0][1] = 1.0;
            dfy[1][0] = (-2.0 * y[0] * y[1] - 1.0) / eps;
            dfy[1][1] = (1.0 - y[0] * y[0]) / eps;
        }
    };

    private static Rodas.Result runRob(int ijac, boolean massFull, int mlmas, int mumas,
                                       double rt, double at) {
        Rodas.Options o = new Rodas.Options();
        o.ijac = ijac;
        o.jac = JROB;
        o.mljac = 3;
        o.mujac = 3;
        o.imas = 1;
        o.mas = massFull ? MFULL3 : MDIAG3;
        o.mlmas = mlmas;
        o.mumas = mumas;
        o.iout = 0;
        double[] y = {1.0, 0.0, 0.0};
        return Rodas.integrate(3, FROB, 0.0, y, 0.4, 1.0e-6,
                new double[]{rt}, new double[]{at}, 0, o);
    }

    private static Rodas.Result runVdp(int ijac, double rt, double at) {
        Rodas.Options o = new Rodas.Options();
        o.ijac = ijac;
        o.jac = JVDP;
        o.mljac = 2;
        o.mujac = 2;
        o.imas = 0;
        o.iout = 0;
        double[] y = {2.0, -0.6};
        return Rodas.integrate(2, FVDP, 0.0, y, 2.0, 1.0e-6,
                new double[]{rt}, new double[]{at}, 0, o);
    }

    private static void assertClose(double got, double want) {
        assertClose(got, want, REL);
    }

    private static void assertClose(double got, double want, double rel) {
        assertEquals(want, got, Math.abs(want) * rel);
    }

    @Test
    public void testDaeBandedMassNumericalJacobian() {
        Rodas.Result r = runRob(0, false, 0, 0, 1.0e-8, 1.0e-10);
        assertEquals(1, r.idid);
        assertEquals(281, r.nfcn);
        assertEquals(46, r.njac);
        assertEquals(47, r.nstep);
        assertEquals(46, r.naccpt);
        assertClose(r.y[0], 0.98517211385739623);
        assertClose(r.y[1], 0.33863953610620744e-4);
        assertClose(r.y[2], 0.14794022188993190e-1);
    }

    @Test
    public void testDaeBandedMassAnalyticJacobian() {
        Rodas.Result r = runRob(1, false, 0, 0, 1.0e-8, 1.0e-10);
        assertEquals(1, r.idid);
        assertEquals(281, r.nfcn);
        assertEquals(46, r.njac);
        assertEquals(47, r.nstep);
        assertClose(r.y[0], 0.98517211385689285);
        assertClose(r.y[1], 0.33863953595594975e-4);
        assertClose(r.y[2], 0.14794022189511580e-1);
    }

    @Test
    public void testDaeTightTolerance() {
        Rodas.Result r = runRob(1, false, 0, 0, 1.0e-11, 1.0e-13);
        assertEquals(1, r.idid);
        assertEquals(1510, r.nfcn);
        assertEquals(250, r.njac);
        assertEquals(252, r.nstep);
        assertClose(r.y[0], 0.98517211386097781);
        assertClose(r.y[1], 0.33863953787836347e-4);
        assertClose(r.y[2], 0.14794022185234504e-1);
    }

    @Test
    public void testFullMassMatrixMatchesBanded() {
        // Same problem, different storage, therefore a different IJOB path
        // through DECOMR and SLVROD. Agreement here is what says the two paths
        // were ported correctly, without needing an external reference for
        // either -- and unlike the cross-language comparison it IS bit for bit,
        // because both paths run in this JVM.
        Rodas.Result full = runRob(1, true, 3, 3, 1.0e-8, 1.0e-10);
        Rodas.Result band = runRob(1, false, 0, 0, 1.0e-8, 1.0e-10);
        assertEquals(1, full.idid);
        assertEquals(band.nfcn, full.nfcn);
        assertEquals(band.nstep, full.nstep);
        for (int i = 0; i < 3; i++) {
            assertTrue(full.y[i] == band.y[i],
                    "full and banded storage must agree exactly, component " + i);
        }
        assertClose(full.y[0], 0.98517211385689285);
    }

    @Test
    public void testStiffOdeIdentityMassAnalyticJacobian() {
        Rodas.Result r = runVdp(1, 1.0e-8, 1.0e-10);
        assertEquals(1, r.idid);
        // The step sequence is still exactly the Fortran's here.
        assertEquals(32446, r.nfcn);
        assertEquals(5401, r.njac);
        assertEquals(5409, r.nstep);
        assertClose(r.y[0], 0.17061674639889246e1, REL_LONG);
        assertClose(r.y[1], -0.89280998821562574, REL_LONG);
    }

    @Test
    public void testStiffOdeIdentityMassNumericalJacobian() {
        Rodas.Result r = runVdp(0, 1.0e-8, 1.0e-10);
        assertEquals(1, r.idid);
        // THE COUNTERS ARE NOT ASSERTED HERE, and this is the only case in the
        // file where they are not. A numerical Jacobian is itself a function of
        // y, so a last-ulp difference in one step size perturbs the Jacobian,
        // which perturbs the error estimate -- a faster amplification path than
        // the analytic case has. The Fortran takes 32755/5420/5467 and this
        // port 32685/5420/5453. The exactness claim is carried by
        // testStiffOdeStepSequenceIsExactUntilThePowDriftAccumulates, which
        // pins the same problem and the same code path at a tolerance the drift
        // has not yet reached.
        assertEquals(5420, r.njac);
        assertClose(r.y[0], 0.17061674637556288e1, REL_LONG_NUMJAC);
        assertClose(r.y[1], -0.89280998834260528, REL_LONG_NUMJAC);
    }

    @Test
    public void testStiffOdeStepSequenceIsExactUntilThePowDriftAccumulates() {
        // The same Van der Pol problem and the same IJOB=1 numerical-Jacobian
        // path as the case above, at a tolerance whose run is short enough that
        // Math.pow's last-ulp difference has not yet flipped a decision. The
        // counters here are what the BIT-EXACT ports (cpp/tests/test_rodas.cpp's
        // solver and python/tests/test_rodas.py's) produce on this input, so
        // this asserts that the transliteration itself is step-for-step
        // faithful -- the property the long run can no longer isolate.
        //
        // This is cross-port agreement, not Fortran conformance: the oracle is
        // the port that IS bit-exact against the Fortran, one link further out.
        Rodas.Options o = new Rodas.Options();
        o.ijac = 0;
        o.jac = JVDP;
        o.mljac = 2;
        o.mujac = 2;
        o.imas = 0;
        o.iout = 0;
        double[] y = {2.0, -0.6};
        Rodas.Result r = Rodas.integrate(2, FVDP, 0.0, y, 2.0, 1.0e-6,
                new double[]{1.0e-6}, new double[]{1.0e-8}, 0, o);
        assertEquals(1, r.idid);
        assertEquals(6733, r.nfcn);
        assertEquals(1113, r.njac);
        assertEquals(1124, r.nstep);
        // and the same at the next tolerance up, with the analytic Jacobian
        Rodas.Options o2 = new Rodas.Options();
        o2.ijac = 1;
        o2.jac = JVDP;
        o2.mljac = 2;
        o2.mujac = 2;
        o2.imas = 0;
        o2.iout = 0;
        double[] y2 = {2.0, -0.6};
        Rodas.Result r2 = Rodas.integrate(2, FVDP, 0.0, y2, 2.0, 1.0e-6,
                new double[]{1.0e-7}, new double[]{1.0e-9}, 0, o2);
        assertEquals(1, r2.idid);
        assertEquals(14102, r2.nfcn);
        assertEquals(2342, r2.njac);
        assertEquals(2352, r2.nstep);
    }

    @Test
    public void testDenseOutputSatisfiesTheAlgebraicRow() {
        // CONTRO, the interpolant the fluid DAE reads its output grid off. The
        // algebraic row is an EQUATION of the DAE, not a quantity that drifts,
        // so an interpolant that only got the differential rows right would
        // fail here and nowhere else.
        final double[] grid = new double[8];
        for (int k = 0; k < 8; k++) {
            grid[k] = 0.05 * (k + 1);
        }
        final List<double[]> seen = new ArrayList<double[]>();

        Rodas.Options o = new Rodas.Options();
        o.ijac = 1;
        o.jac = JROB;
        o.mljac = 3;
        o.mujac = 3;
        o.imas = 1;
        o.mas = MDIAG3;
        o.mlmas = 0;
        o.mumas = 0;
        o.iout = 1;
        o.solout = new Rodas.Solout() {
            public int eval(int nr, double xold, double x, double[] y, Rodas.Dense d) {
                while (seen.size() < grid.length && grid[seen.size()] <= x + 1e-13) {
                    double tq = grid[seen.size()];
                    double[] xs = new double[3];
                    for (int i = 0; i < 3; i++) {
                        // nr <= 1 is the call made BEFORE any step, where cont
                        // holds no coefficients yet and the state is y itself.
                        xs[i] = (nr <= 1) ? y[i] : d.value(i, tq);
                    }
                    seen.add(xs);
                }
                return 0;
            }
        };
        double[] y = {1.0, 0.0, 0.0};
        Rodas.Result r = Rodas.integrate(3, FROB, 0.0, y, 0.4, 1.0e-6,
                new double[]{1.0e-8}, new double[]{1.0e-10}, 0, o);
        assertEquals(1, r.idid);
        assertEquals(grid.length, seen.size());
        for (int k = 0; k < seen.size(); k++) {
            double[] xs = seen.get(k);
            assertEquals(1.0, xs[0] + xs[1] + xs[2], 1e-9,
                    "constraint at t=" + grid[k]);
        }
    }

    @Test
    public void testIoutDoesNotChangeTheTrajectory() {
        // Turning dense output on must cost steps, not accuracy: cont is written
        // from stage values that already exist, so the step sequence is the same
        // one. A port that fed CONT back into the step would fail here.
        Rodas.Result quiet = runRob(1, false, 0, 0, 1.0e-8, 1.0e-10);

        Rodas.Options o = new Rodas.Options();
        o.ijac = 1;
        o.jac = JROB;
        o.mljac = 3;
        o.mujac = 3;
        o.imas = 1;
        o.mas = MDIAG3;
        o.iout = 1;
        o.solout = new Rodas.Solout() {
            public int eval(int nr, double xold, double x, double[] y, Rodas.Dense d) {
                return 0;
            }
        };
        double[] y = {1.0, 0.0, 0.0};
        Rodas.Result loud = Rodas.integrate(3, FROB, 0.0, y, 0.4, 1.0e-6,
                new double[]{1.0e-8}, new double[]{1.0e-10}, 0, o);
        assertEquals(quiet.nstep, loud.nstep);
        assertEquals(quiet.nfcn, loud.nfcn);
        for (int i = 0; i < 3; i++) {
            assertTrue(loud.y[i] == quiet.y[i]);
        }
    }

    @Test
    public void testSoloutCanStopTheIntegration() {
        Rodas.Options o = new Rodas.Options();
        o.ijac = 1;
        o.jac = JROB;
        o.mljac = 3;
        o.mujac = 3;
        o.imas = 1;
        o.mas = MDIAG3;
        o.iout = 1;
        o.solout = new Rodas.Solout() {
            public int eval(int nr, double xold, double x, double[] y, Rodas.Dense d) {
                return x > 0.1 ? -1 : 0;
            }
        };
        double[] y = {1.0, 0.0, 0.0};
        Rodas.Result r = Rodas.integrate(3, FROB, 0.0, y, 0.4, 1.0e-6,
                new double[]{1.0e-8}, new double[]{1.0e-10}, 0, o);
        assertEquals(2, r.idid);
        assertTrue(r.x < 0.4);
    }

    @Test
    public void testTolerancesTooSmallAreRefusedByName() {
        final Rodas.Options o = new Rodas.Options();
        o.ijac = 1;
        o.jac = JROB;
        o.mljac = 3;
        o.mujac = 3;
        o.imas = 1;
        o.mas = MDIAG3;
        final double[] y = {1.0, 0.0, 0.0};
        assertThrows(Rodas.RodasInputException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Rodas.integrate(3, FROB, 0.0, y, 0.4, 1.0e-6,
                        new double[]{1.0e-18}, new double[]{1.0e-20}, 0, o);
            }
        });
    }
}
