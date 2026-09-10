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
