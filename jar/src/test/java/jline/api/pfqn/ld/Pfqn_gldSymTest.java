package jline.api.pfqn.ld;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.Map;

import org.junit.jupiter.api.Test;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.SymMatrix;
import jline.util.symbolic.SymContext;
import jline.util.symbolic.SymExpr;

/**
 * The SYMBOLIC arm of the gld family, against the numeric one.
 *
 * EVERY ASSERTION PINS THE SUBSTITUTED VALUE against the double routine, as the
 * native python test_pfqn_symbolic.py does and for the same reason: a wrong
 * expression must not be able to pass by being merely well-formed. Checking the
 * printed form instead would pin a normal form Rings is free to change.
 */
public class Pfqn_gldSymTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++)
            for (int j = 0; j < v[i].length; j++) m.set(i, j, v[i][j]);
        return m;
    }

    private static Map<String, Double> at(String[] names, double... vals) {
        Map<String, Double> m = new HashMap<String, Double>();
        for (int i = 0; i < names.length; i++) m.put(names[i], Double.valueOf(vals[i]));
        return m;
    }

    /** A symbolic RATE row needs no threshold declared anywhere. */
    @Test
    public void symbolicRatesSingleClass() {
        SymContext ctx = SymContext.of("m1", "m2", "c");
        SymMatrix L = new SymMatrix(ctx, 2, 1);
        L.set(0, 0, ctx.constant(2.0));
        L.set(1, 0, ctx.constant(3.0));
        SymMatrix mu = new SymMatrix(ctx, 2, 3);
        mu.set(0, 0, ctx.var("m1"));
        mu.set(0, 1, ctx.var("m2"));
        mu.set(0, 2, ctx.var("c"));
        for (int j = 0; j < 3; j++) mu.set(1, j, ctx.one());

        Ret.pfqnNcSym s = Pfqn_gldsingle_sym.pfqn_gldsingle_sym(L, row(3), mu);

        Matrix Ln = mat(new double[][]{{2.0}, {3.0}});
        Matrix mu1 = new Matrix(2, 3);
        for (int i = 0; i < 2; i++) for (int j = 0; j < 3; j++) mu1.set(i, j, 1.0);
        double ref1 = Math.exp(Pfqn_gldsingle.pfqn_gldsingle(Ln, row(3), mu1, new SolverOptions()).lG);
        assertEquals(ref1, s.G.evaluate(at(new String[]{"m1", "m2", "c"}, 1, 1, 1)), TOL);

        // THE POINT OF THE gld FAMILY: a multiserver row given as 1,2,2 is served
        // without anyone identifying where it settles.
        Matrix mu2 = mat(new double[][]{{1.0, 2.0, 2.0}, {1.0, 1.0, 1.0}});
        double ref2 = Math.exp(Pfqn_gldsingle.pfqn_gldsingle(Ln, row(3), mu2, new SolverOptions()).lG);
        assertEquals(ref2, s.G.evaluate(at(new String[]{"m1", "m2", "c"}, 1, 2, 2)), TOL);
    }

    /** Symbolic DEMANDS, multiclass, against the numeric recursion. */
    @Test
    public void symbolicDemandsMulticlass() {
        String[] names = {"L11", "L12", "L21", "L22"};
        SymContext ctx = SymContext.of(names);
        SymMatrix L = new SymMatrix(ctx, 2, 2);
        L.set(0, 0, ctx.var("L11"));
        L.set(0, 1, ctx.var("L12"));
        L.set(1, 0, ctx.var("L21"));
        L.set(1, 1, ctx.var("L22"));
        SymMatrix mu = new SymMatrix(ctx, 2, 4);
        for (int i = 0; i < 2; i++) for (int j = 0; j < 4; j++) mu.set(i, j, ctx.one());

        Ret.pfqnNcSym s = Pfqn_gld_sym.pfqn_gld_sym(L, row(2, 2), mu);

        Matrix Ln = mat(new double[][]{{1.0, 0.6}, {0.5, 1.1}});
        Matrix mun = new Matrix(2, 4);
        for (int i = 0; i < 2; i++) for (int j = 0; j < 4; j++) mun.set(i, j, 1.0);
        double ref = Math.exp(Pfqn_gld.pfqn_gld(Ln, row(2, 2), mun, new SolverOptions()).lG);
        assertEquals(ref, s.G.evaluate(at(names, 1.0, 0.6, 0.5, 1.1)), TOL);
    }

    /**
     * THE DENOMINATOR IS A MONOMIAL, which is the structural fact that makes a
     * rational-function field the right home for this recursion: every division
     * is by a bare mu symbol, so nothing here needs a gcd to stay in normal form.
     */
    @Test
    public void denominatorIsAMonomial() {
        SymContext ctx = SymContext.of("m1", "m2", "c");
        SymMatrix L = new SymMatrix(ctx, 2, 1);
        L.set(0, 0, ctx.constant(2.0));
        L.set(1, 0, ctx.constant(3.0));
        SymMatrix mu = new SymMatrix(ctx, 2, 3);
        mu.set(0, 0, ctx.var("m1"));
        mu.set(0, 1, ctx.var("m2"));
        mu.set(0, 2, ctx.var("c"));
        for (int j = 0; j < 3; j++) mu.set(1, j, ctx.one());

        Ret.pfqnNcSym s = Pfqn_gldsingle_sym.pfqn_gldsingle_sym(L, row(3), mu);
        assertEquals(1, s.G.value().denominator().size(),
                "every divisor is a bare mu symbol, so the denominator is one monomial");
    }

    /** The M==1 closed form: exact factorials, not exp(factln). */
    @Test
    public void singleStationClosedFormIsExact() {
        String[] names = {"a", "b"};
        SymContext ctx = SymContext.of(names);
        SymMatrix L = new SymMatrix(ctx, 1, 2);
        L.set(0, 0, ctx.var("a"));
        L.set(0, 1, ctx.var("b"));
        SymMatrix mu = new SymMatrix(ctx, 1, 5);
        for (int j = 0; j < 5; j++) mu.set(0, j, ctx.one());

        Ret.pfqnNcSym s = Pfqn_gld_sym.pfqn_gld_sym(L, row(2, 3), mu);
        // multinomial(5; 2,3) = 10, so G = 10*a^2*b^3 at mu == 1
        assertEquals(10.0 * 4.0 * 27.0, s.G.evaluate(at(names, 2.0, 3.0)), TOL);

        Matrix Ln = mat(new double[][]{{2.0, 3.0}});
        Matrix mun = new Matrix(1, 5);
        for (int j = 0; j < 5; j++) mun.set(0, j, 1.0);
        double ref = Math.exp(Pfqn_gld.pfqn_gld(Ln, row(2, 3), mun, new SolverOptions()).lG);
        assertEquals(ref, s.G.evaluate(at(names, 2.0, 3.0)), TOL);
    }

    /** A zero population is one, and the constant carries no variable. */
    @Test
    public void emptyPopulationIsOne() {
        SymContext ctx = SymContext.of("x");
        SymMatrix L = new SymMatrix(ctx, 2, 1);
        L.set(0, 0, ctx.var("x"));
        L.set(1, 0, ctx.one());
        SymMatrix mu = new SymMatrix(ctx, 2, 1);
        mu.set(0, 0, ctx.one());
        mu.set(1, 0, ctx.one());
        Ret.pfqnNcSym s = Pfqn_gld_sym.pfqn_gld_sym(L, row(0), mu);
        assertTrue(s.G.isOne());
    }

    /** Two contexts cannot be combined, and the message says why. */
    @Test
    public void contextsDoNotMix() {
        SymContext a = SymContext.of("x");
        SymContext b = SymContext.of("x");
        SymExpr xa = a.var("x");
        SymExpr xb = b.var("x");
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                xa.add(xb);
            }
        });
    }

    /** A partial substitution is still symbolic and has no double to return. */
    @Test
    public void partialSubstitutionIsRefused() {
        SymContext ctx = SymContext.of("p", "q");
        final SymExpr e = ctx.var("p").add(ctx.var("q"));
        final Map<String, Double> partial = new HashMap<String, Double>();
        partial.put("p", Double.valueOf(1.0));
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                e.evaluate(partial);
            }
        });
    }
}
