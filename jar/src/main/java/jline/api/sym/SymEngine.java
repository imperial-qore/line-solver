package jline.api.sym;

import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * Computer algebra operations LINE needs, as seen by the JAR.
 *
 * <p>Java has no computer algebra system, so every operation here is delegated
 * to an external engine; {@link SageRestEngine} is the SageMath implementation
 * and is resolved by {@link SymEngines}. Expressions cross the interface as
 * plain ASCII infix strings, e.g. {@code "2*x1 - 3*x2"}, which is the format
 * {@code SolverCTMC.symbolicGeneratorResult.getSymbolicEntry} already
 * produces.</p>
 *
 * <p>Expression strings are <em>not</em> comparable across codebases: the
 * symbol numbering x1..xE follows event enumeration order, and printed normal
 * forms depend on the engine version. Compare by substituting values with
 * {@link #eval} and comparing numbers.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public interface SymEngine {

    /** Name of the backing engine, e.g. "sage". */
    String name();

    /** True if the engine answers a health probe. */
    boolean isAvailable();

    /**
     * Symbolic stationary distribution of a CTMC, pi Q = 0 with sum(pi) = 1.
     *
     * @param Q       generator entries as expression strings, row major
     * @param symbols the symbols appearing in Q, e.g. x1..xE
     * @return the solution, with pi also split over a common denominator
     * @throws IOException if the engine is unreachable or rejects the request
     */
    CTMCSolution solveCTMC(String[][] Q, List<String> symbols) throws IOException;

    /**
     * Exact parametric sensitivity of a steady-state reward.
     *
     * @param Q       generator entries as expression strings, row major
     * @param symbols the symbols appearing in Q
     * @param theta   the symbol to differentiate with respect to
     * @param reward  reward rate per state, or null for the distribution alone
     * @return the sensitivity of the distribution and, if a reward is given, of its mean
     * @throws IOException if the engine is unreachable or rejects the request
     */
    Sensitivity ctmcSensitivity(String[][] Q, List<String> symbols, String theta,
                                List<String> reward) throws IOException;

    /**
     * Stationary distribution and a set of weighted sums of it.
     *
     * <p>Every mean SolverCTMC reports is a NUMERIC linear functional of pi, so
     * the caller builds one weight vector per measure with no algebra at all
     * and this one round trip returns each {@code w . pi} as a rational
     * function: a throughput is {@code pi . depRates}, a queue length is
     * {@code pi . stateSpaceAggr}, a marginal probability is
     * {@code pi . indicator} and a mean reward is {@code pi . r}.</p>
     *
     * @param Q       generator entries as expression strings, row major
     * @param symbols the symbols appearing in Q
     * @param weights one row per measure, each of length Q.length; entries are
     *                parsed in the same field as Q, so an expression is
     *                admissible and not only a number
     * @param names   a name per weight row, or null for w0, w1, ...
     * @param ratios  ratios of two measures, or null for none
     * @return the distribution, the measures and the ratios
     * @throws IOException if the engine is unreachable or rejects the request
     */
    CTMCMeasures ctmcMeasures(String[][] Q, List<String> symbols, String[][] weights,
                              List<String> names, List<RatioSpec> ratios) throws IOException;

    /**
     * First passage time into a target set: the transform and the moments.
     *
     * <p>With A the complement of the target, S = Q(A,A) the sub-generator,
     * s0 = -S*1 the exit vector and alpha the initial law restricted to A,
     * following Harrison and Knottenbelt (2002),</p>
     *
     * <pre>  L(s) = alpha (sI - S)^-1 s0 + atom        Eqs. 1-2
     *  (-S) M(n) = n M(n-1),  M(0) = 1           Eq. 3</pre>
     *
     * <p>The MOMENTS come from the recursion, not from differentiating L: it
     * carries no transform symbol, stays in the same fraction field and is
     * {@code nmax} right solves against one matrix. Differentiating instead
     * would leave the fraction field, and the evaluation at s = 0 can hit 0/0
     * where a factor of s failed to cancel.</p>
     *
     * @param S       the sub-generator on the non-target states, row major
     * @param s0      the exit vector -S*1
     * @param alpha   the initial law restricted to the non-target states
     * @param atom    the initial mass already inside the target set
     * @param symbols the rate symbols appearing in S
     * @param svar    the transform symbol; adjoined to the field only when a
     *                transform is asked for, and required to differ from every
     *                rate symbol
     * @param want    any of lst, lstall, moments, momall
     * @param nmax    highest moment order, at least 1
     * @return the requested items
     * @throws IOException if the engine is unreachable or rejects the request
     */
    Passage ctmcPassage(String[][] S, String[] s0, String[] alpha, String atom,
                        List<String> symbols, String svar, List<String> want,
                        int nmax) throws IOException;

    /**
     * Rewrites expressions into a normal form.
     *
     * @param exprs the expressions
     * @param form  one of simplify, factor, together, cancel, expand, latex
     * @return the rewritten expressions, in the input order
     * @throws IOException if the engine is unreachable or rejects the request
     */
    List<String> simplify(List<String> exprs, String form) throws IOException;

    /**
     * Differentiates expressions.
     *
     * @param exprs    the expressions
     * @param variable the differentiation variable
     * @param order    the order of the derivative, at least 1
     * @return the derivatives, in the input order
     * @throws IOException if the engine is unreachable or rejects the request
     */
    List<String> diff(List<String> exprs, String variable, int order) throws IOException;

    /**
     * Substitutes values for symbols and evaluates.
     *
     * @param exprs      the expressions
     * @param assignment value of each symbol
     * @return the numeric values, NaN where a free symbol remains
     * @throws IOException if the engine is unreachable or rejects the request
     */
    double[] eval(List<String> exprs, Map<String, Double> assignment) throws IOException;

    /**
     * Jacobian, LaTeX form and equilibria of a fluid vector field.
     *
     * @param rhs  the right hand side of dx/dt, one expression per state variable
     * @param vars the state variable names
     * @param want any of jacobian, latex, equilibria
     * @return the requested items
     * @throws IOException if the engine is unreachable or rejects the request
     */
    FluidODEs fluidODEs(List<String> rhs, List<String> vars, List<String> want) throws IOException;

    /** Symbolic stationary distribution of a CTMC. */
    class CTMCSolution {
        /** Stationary probability of each state, as an expression string */
        public final List<String> pi;
        /** Numerator of each entry over the common denominator {@link #den} */
        public final List<String> num;
        /** Common denominator of the whole vector */
        public final String den;
        /** Number of weakly connected components of the generator */
        public final int nConnComp;
        /** Component index of each state, one based */
        public final int[] connComp;

        public CTMCSolution(List<String> pi, List<String> num, String den,
                            int nConnComp, int[] connComp) {
            this.pi = pi;
            this.num = num;
            this.den = den;
            this.nConnComp = nConnComp;
            this.connComp = connComp;
        }
    }

    /** One weighted sum of the stationary distribution, or a ratio of two. */
    class Measure {
        /** Caller supplied name */
        public final String name;
        /** The value as a rational function, or "undefined" */
        public final String expr;
        /** Numerator of {@link #expr}, null when undefined */
        public final String num;
        /** Denominator of {@link #expr}, null when undefined */
        public final String den;
        /**
         * Why the value is undefined, null otherwise.
         *
         * <p>A ratio whose denominator measure is identically zero over the
         * whole rate space has no value. The numeric arms sweep such a 0/0 to
         * zero with isnan, which is a defensible cleanup of floating point dust
         * and an indefensible answer for an exact one.</p>
         */
        public final String reason;

        public Measure(String name, String expr, String num, String den, String reason) {
            this.name = name;
            this.expr = expr;
            this.num = num;
            this.den = den;
            this.reason = reason;
        }

        /** @return true if this measure has a value */
        public boolean isDefined() {
            return this.reason == null;
        }
    }

    /** A ratio of two measures, by their index into the weight block. */
    class RatioSpec {
        /** Name of the resulting measure */
        public final String name;
        /** Index of the numerator measure */
        public final int num;
        /** Index of the denominator measure */
        public final int den;

        public RatioSpec(String name, int num, int den) {
            this.name = name;
            this.num = num;
            this.den = den;
        }
    }

    /** Stationary distribution together with the measures taken from it. */
    class CTMCMeasures {
        /** Stationary probability of each state */
        public final List<String> pi;
        /** Numerator of each entry over the common denominator {@link #den} */
        public final List<String> num;
        /** Common denominator of the whole vector */
        public final String den;
        /** Number of weakly connected components of the generator */
        public final int nConnComp;
        /** Component index of each state, one based */
        public final int[] connComp;
        /** One entry per weight vector, in the order given */
        public final List<Measure> measures;
        /** One entry per requested ratio, in the order given */
        public final List<Measure> ratios;

        public CTMCMeasures(List<String> pi, List<String> num, String den, int nConnComp,
                            int[] connComp, List<Measure> measures, List<Measure> ratios) {
            this.pi = pi;
            this.num = num;
            this.den = den;
            this.nConnComp = nConnComp;
            this.connComp = connComp;
            this.measures = measures;
            this.ratios = ratios;
        }
    }

    /** First passage time transform and moments. */
    class Passage {
        /** The transform L(s), null if not requested */
        public final String lst;
        /** Numerator of {@link #lst}, null if not requested */
        public final String lstNum;
        /** Denominator of {@link #lst}, null if not requested */
        public final String lstDen;
        /** Per start state transform, null if not requested */
        public final List<String> lstAll;
        /** Moments of order 1..nmax for the initial law, null if not requested */
        public final List<String> moments;
        /** Per start state moments, null if not requested */
        public final List<List<String>> momAll;
        /** States that cannot reach the target at all */
        public final int[] unreachable;

        public Passage(String lst, String lstNum, String lstDen, List<String> lstAll,
                       List<String> moments, List<List<String>> momAll, int[] unreachable) {
            this.lst = lst;
            this.lstNum = lstNum;
            this.lstDen = lstDen;
            this.lstAll = lstAll;
            this.moments = moments;
            this.momAll = momAll;
            this.unreachable = unreachable;
        }
    }

    /**
     * Exact parametric sensitivity, following Trivedi and Bobbio (2017), Sec. 9.7.
     *
     * <p>As in {@code @SolverCTMC/getSensitivity.m}, dr/dtheta is taken to be
     * zero: a reward whose rates themselves depend on theta needs the second
     * term of Eq. (9.83) and is not covered here.</p>
     */
    class Sensitivity {
        /** Stationary distribution */
        public final List<String> pi;
        /** Derivative of the stationary distribution with respect to theta */
        public final List<String> dpi;
        /** Mean reward, null if no reward was given */
        public final String Er;
        /** Unscaled sensitivity d(E[r])/dtheta, Eq. (9.79); null if no reward was given */
        public final String S;
        /** Scaled sensitivity (theta/E[r]) d(E[r])/dtheta, Eq. (9.80); null if no reward was given */
        public final String SS;

        public Sensitivity(List<String> pi, List<String> dpi, String Er, String S, String SS) {
            this.pi = pi;
            this.dpi = dpi;
            this.Er = Er;
            this.S = S;
            this.SS = SS;
        }
    }

    /** Symbolic analysis of a fluid vector field. */
    class FluidODEs {
        /** Jacobian d f_i / d x_j, or null if not requested */
        public final String[][] jacobian;
        /** LaTeX form of each right hand side, or null if not requested */
        public final List<String> latex;
        /** Equilibria as variable to expression maps, or null if not requested */
        public final List<Map<String, String>> equilibria;

        public FluidODEs(String[][] jacobian, List<String> latex,
                         List<Map<String, String>> equilibria) {
            this.jacobian = jacobian;
            this.latex = latex;
            this.equilibria = equilibria;
        }
    }
}
