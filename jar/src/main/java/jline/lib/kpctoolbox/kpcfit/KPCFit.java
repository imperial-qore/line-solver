package jline.lib.kpctoolbox.kpcfit;

import jline.api.mam.Map_acf;
import jline.api.mam.Map_embedded;
import jline.api.mam.Map_erlang;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_feasblock;
import jline.api.mam.Map_isfeasible;
import jline.api.mam.Map_joint;
import jline.api.mam.Map_kpc;
import jline.api.mam.Map_moment;
import jline.api.mam.Map_normalize;
import jline.api.mam.Map_scale;
import jline.api.mam.Map_scv;
import jline.api.mam.Map2_fit;
import jline.io.Ret;
import jline.io.Ret.mamMAPFitReturn;
import jline.api.mam.Aph_fit;
import jline.lib.kpctoolbox.basic.BasicUtils;
import jline.lib.kpctoolbox.mmpp.MMPP;
import jline.lib.kpctoolbox.trace.TraceAnalysis;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Maths;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * KPC-Toolbox fitting functions.
 * Mechanically translated from KPCFit.kt.
 */
public final class KPCFit {

    public static final double KPCFIT_TOL = 1e-10;

    private KPCFit() {}

    private static MatrixCell arrayToMatrixCell(Matrix[] arr) {
        MatrixCell cell = new MatrixCell(arr.length);
        for (int i = 0; i < arr.length; i++) {
            cell.set(i, arr[i]);
        }
        return cell;
    }

    private static Matrix[] matrixCellToArray(MatrixCell cell) {
        int size = cell.size();
        Matrix[] arr = new Matrix[size];
        for (int i = 0; i < size; i++) {
            arr[i] = cell.get(i);
        }
        return arr;
    }

    public static int[] distinctInts(int[] arr) {
        Set<Integer> seen = new LinkedHashSet<Integer>();
        for (int v : arr) seen.add(v);
        int[] out = new int[seen.size()];
        int i = 0;
        for (Integer v : seen) out[i++] = v;
        return out;
    }

    private static int maxOfInt(int[] arr, int defaultValue) {
        if (arr.length == 0) return defaultValue;
        int m = arr[Integer.MIN_VALUE >>> 31]; // 0
        m = arr[0];
        for (int v : arr) if (v > m) m = v;
        return m;
    }

    public static class TraceData {
        public final double[] S;
        public final double[] E;
        public final double[] AC;
        public final double[] ACFull;
        public final int[] ACLags;
        public final double[] BC;
        public final int[] BCGridLags;
        public final int[][] BCLags;

        public TraceData(double[] S, double[] E, double[] AC, double[] ACFull, int[] ACLags,
                         double[] BC, int[] BCGridLags, int[][] BCLags) {
            this.S = S;
            this.E = E;
            this.AC = AC;
            this.ACFull = ACFull;
            this.ACLags = ACLags;
            this.BC = BC;
            this.BCGridLags = BCGridLags;
            this.BCLags = BCLags;
        }
    }

    public static class KPCFitOptions {
        public boolean onlyAC = false;
        public Integer numMAPs = null;
        public Integer numStates = null;
        public int maxIterAC = 300;
        public int maxIterBC = 10;
        public int maxRunsAC = 50;
        public int maxRunsBC = 30;
        public int maxResAC = 10;
        public int maxRetMAPs = 1;

        public KPCFitOptions() {}
    }

    public static class KPCFitResult {
        public final MatrixCell MAP;
        public final double fac;
        public final double fbc;
        public final List<MatrixCell> subMAPs;
        public final List<MatrixCell> otherMAPs;
        public final double[] otherFACs;
        public final double[] otherFBCs;

        public KPCFitResult(MatrixCell MAP, double fac, double fbc, List<MatrixCell> subMAPs,
                            List<MatrixCell> otherMAPs, double[] otherFACs, double[] otherFBCs) {
            this.MAP = MAP;
            this.fac = fac;
            this.fbc = fbc;
            this.subMAPs = subMAPs;
            this.otherMAPs = otherMAPs;
            this.otherFACs = otherFACs;
            this.otherFBCs = otherFBCs;
        }
    }

    public static TraceData kpcfit_init(double[] S, int[] acLags, int[] bcGridLags) {
        int n = S.length;
        int nMinSupportAC = 10;

        int[] defaultACLags;
        if (acLags != null) {
            defaultACLags = acLags;
        } else {
            double upper = Math.max(1.0, Math.ceil((double) n / nMinSupportAC));
            defaultACLags = distinctInts(BasicUtils.logspacei(1.0, upper, 500));
        }

        int maxACLag = 1;
        for (int v : defaultACLags) if (v > maxACLag) maxACLag = v;

        int[] defaultBCGridLags;
        if (bcGridLags != null) {
            defaultBCGridLags = bcGridLags;
        } else {
            defaultBCGridLags = distinctInts(BasicUtils.logspacei(1.0, (double) maxACLag, 5));
        }

        // Compute moments E[X], E[X^2], E[X^3]
        double[] E = new double[3];
        for (int j = 1; j <= 3; j++) {
            double sum = 0.0;
            for (double v : S) sum += FastMath.pow(v, j);
            E[j - 1] = sum / S.length;
        }

        // Filter ACLags to <= n - 2
        List<Integer> validList = new ArrayList<Integer>();
        for (int v : defaultACLags) if (v <= n - 2) validList.add(v);
        int[] validACLags = new int[validList.size()];
        for (int i = 0; i < validACLags.length; i++) validACLags[i] = validList.get(i);

        double[] AC = (validACLags.length > 0) ? TraceAnalysis.trace_acf(S, validACLags) : new double[0];

        int upperAcFull = Math.max(1, (int) Math.ceil((double) n / nMinSupportAC));
        int[] fullLags = new int[upperAcFull];
        for (int i = 0; i < upperAcFull; i++) fullLags[i] = i + 1;
        double[] ACFull = TraceAnalysis.trace_acf(S, fullLags);

        int cutIdx = validACLags.length;
        for (int i = 0; i < AC.length; i++) {
            if (Math.abs(AC[i]) < 1e-6) {
                cutIdx = i + 1;
                break;
            }
        }

        int[] trimmedACLags = new int[Math.min(cutIdx, validACLags.length)];
        System.arraycopy(validACLags, 0, trimmedACLags, 0, trimmedACLags.length);
        double[] trimmedAC = new double[Math.min(cutIdx, AC.length)];
        System.arraycopy(AC, 0, trimmedAC, 0, trimmedAC.length);

        int trimmedMax = 1;
        for (int v : trimmedACLags) if (v > trimmedMax) trimmedMax = v;

        List<Integer> validBCList = new ArrayList<Integer>();
        for (int v : defaultBCGridLags) if (v <= trimmedMax) validBCList.add(v);
        int[] validBCGridLags = new int[validBCList.size()];
        for (int i = 0; i < validBCGridLags.length; i++) validBCGridLags[i] = validBCList.get(i);

        jline.util.Pair<double[], int[][]> bcResult = TraceAnalysis.trace_bicov(S, validBCGridLags);

        return new TraceData(S, E, trimmedAC, ACFull, trimmedACLags,
                bcResult.getFirst(), validBCGridLags, bcResult.getSecond());
    }

    public static TraceData kpcfit_init(double[] S) {
        return kpcfit_init(S, null, null);
    }

    public static KPCFitResult kpcfit_auto(TraceData trace, KPCFitOptions options) {
        int numMAPs;
        if (options.numMAPs != null) {
            numMAPs = options.numMAPs;
        } else if (options.numStates != null) {
            numMAPs = (int) FastMath.ceil(FastMath.log(2.0, options.numStates.doubleValue()));
        } else {
            numMAPs = kpcfit_sub_bic(trace.ACFull, new int[]{2, 4, 8, 16, 32, 64, 128});
        }

        return kpcfit_manual(numMAPs, trace.E, trace.AC, trace.ACLags, trace.BC, trace.BCLags, options);
    }

    public static KPCFitResult kpcfit_auto(TraceData trace) {
        return kpcfit_auto(trace, new KPCFitOptions());
    }

    public static class AcfitResult {
        public final List<double[]> resSCV;
        public final List<double[]> resG2;
        public final double[] fobjAC;

        public AcfitResult(List<double[]> resSCV, List<double[]> resG2, double[] fobjAC) {
            this.resSCV = resSCV;
            this.resG2 = resG2;
            this.fobjAC = fobjAC;
        }
    }

    public static class BcfitResult {
        public final double[] E1j;
        public final double[] E3j;
        public final double f;

        public BcfitResult(double[] E1j, double[] E3j, double f) {
            this.E1j = E1j;
            this.E3j = E3j;
            this.f = f;
        }
    }

    public static class ComposeResult {
        public final MatrixCell map;
        public final List<MatrixCell> subMAPs;
        public final int errorCode;

        public ComposeResult(MatrixCell map, List<MatrixCell> subMAPs, int errorCode) {
            this.map = map;
            this.subMAPs = subMAPs;
            this.errorCode = errorCode;
        }
    }

    public static class AcfitEvalResult {
        public final double SCVcum;
        public final double[] acfCoeff;

        public AcfitEvalResult(double SCVcum, double[] acfCoeff) {
            this.SCVcum = SCVcum;
            this.acfCoeff = acfCoeff;
        }
    }

    public static KPCFitResult kpcfit_manual(int numMAPs, double[] E, double[] AC, int[] ACLags,
                                             double[] BC, int[][] BCLags, KPCFitOptions options) {
        AcfitResult acRes = kpcfit_sub_acfit(E, AC, ACLags, numMAPs,
                options.maxIterAC, options.maxRunsAC, options.maxResAC);
        List<double[]> resSCV = acRes.resSCV;
        List<double[]> resG2 = acRes.resG2;

        ArrayList<double[]> resE1 = new ArrayList<double[]>();
        ArrayList<double[]> resE3 = new ArrayList<double[]>();
        double[] fobjBC = new double[resSCV.size()];

        if (!options.onlyAC) {
            for (int i = 0; i < resSCV.size(); i++) {
                BcfitResult bc = kpcfit_sub_bcfit(E, resSCV.get(i), resG2.get(i), BC, BCLags,
                        options.maxIterBC, options.maxRunsBC);
                resE1.add(bc.E1j);
                resE3.add(bc.E3j);
                fobjBC[i] = bc.f;
            }
        } else {
            for (int i = 0; i < resSCV.size(); i++) {
                double[] E1j = new double[numMAPs];
                Arrays.fill(E1j, 1.0);
                double[] E3j = new double[numMAPs];
                for (int j = 0; j < numMAPs; j++) {
                    E3j[j] = (1.5 + 0.01) * FastMath.pow(1 + resSCV.get(i)[j], 2.0);
                }
                resE1.add(E1j);
                resE3.add(E3j);
                fobjBC[i] = -1.0;
            }
        }

        // Sort indices by fobjBC
        Integer[] sortedIndices = new Integer[fobjBC.length];
        for (int i = 0; i < fobjBC.length; i++) sortedIndices[i] = i;
        final double[] fobjBCFinal = fobjBC;
        Arrays.sort(sortedIndices, (a, b) -> Double.compare(fobjBCFinal[a], fobjBCFinal[b]));

        ArrayList<MatrixCell> MAPs = new ArrayList<MatrixCell>();
        ArrayList<List<MatrixCell>> subs = new ArrayList<List<MatrixCell>>();
        ArrayList<Double> FACs = new ArrayList<Double>();
        ArrayList<Double> FBCs = new ArrayList<Double>();

        for (int kIdx : sortedIndices) {
            ComposeResult cr = kpcfit_sub_compose(resE1.get(kIdx), resSCV.get(kIdx),
                    resE3.get(kIdx), resG2.get(kIdx));
            if (cr.errorCode != 0 || cr.map == null) continue;

            MatrixCell scaledMAP = Map_scale.map_scale(cr.map, E[0]);
            double[] objs = evaluateObjFunction(scaledMAP, E, AC, ACLags, BC, BCLags);

            MAPs.add(scaledMAP);
            subs.add(cr.subMAPs);
            FACs.add(objs[0]);
            FBCs.add(objs[1]);

            if (MAPs.size() >= options.maxRetMAPs) break;
        }

        if (MAPs.isEmpty()) {
            throw new IllegalStateException("KPC fitting failed - no valid MAP found");
        }

        List<MatrixCell> otherMAPs = (MAPs.size() > 1) ? new ArrayList<MatrixCell>(MAPs.subList(1, MAPs.size())) : new ArrayList<MatrixCell>();
        double[] otherFACs = new double[Math.max(0, FACs.size() - 1)];
        for (int i = 1; i < FACs.size(); i++) otherFACs[i - 1] = FACs.get(i);
        double[] otherFBCs = new double[Math.max(0, FBCs.size() - 1)];
        for (int i = 1; i < FBCs.size(); i++) otherFBCs[i - 1] = FBCs.get(i);

        return new KPCFitResult(MAPs.get(0), FACs.get(0), FBCs.get(0), subs.get(0),
                otherMAPs, otherFACs, otherFBCs);
    }

    public static KPCFitResult kpcfit_manual(int numMAPs, double[] E, double[] AC, int[] ACLags,
                                             double[] BC, int[][] BCLags) {
        return kpcfit_manual(numMAPs, E, AC, ACLags, BC, BCLags, new KPCFitOptions());
    }

    public static int kpcfit_sub_bic(double[] ACFull, int[] orders) {
        int nlags = ACFull.length;
        int ordermax = 2;
        for (int v : orders) if (v > ordermax) ordermax = v;

        int nlagsend = nlags;
        for (int i = 0; i < ACFull.length; i++) {
            if (ACFull[i] < 1e-6) {
                nlagsend = i + ordermax;
                break;
            }
        }

        int NLAGSMAX = 10000;

        int[] effectiveOrders = orders;
        int[] SAlags;

        if (nlagsend > NLAGSMAX) {
            int[] logLags = BasicUtils.logspacei(1.0, (double) (nlagsend - ordermax), NLAGSMAX);
            SAlags = distinctInts(logLags);
        } else {
            if (nlagsend > ordermax) {
                SAlags = new int[nlagsend - ordermax];
                for (int i = 0; i < SAlags.length; i++) SAlags[i] = i + 1;
            } else {
                int effectiveOrdermax = nlagsend - 2;
                if (effectiveOrdermax < 1) return 1;
                int upper = Math.max(1, nlagsend - effectiveOrdermax);
                SAlags = new int[upper];
                for (int i = 0; i < upper; i++) SAlags[i] = i + 1;
                List<Integer> filtered = new ArrayList<Integer>();
                for (int v : orders) if (v <= effectiveOrdermax) filtered.add(v);
                if (filtered.isEmpty()) return 1;
                effectiveOrders = new int[filtered.size()];
                for (int i = 0; i < effectiveOrders.length; i++) effectiveOrders[i] = filtered.get(i);
            }
        }

        int nSamples = SAlags.length;
        if (nSamples < 2) return 1;

        double[] Y = new double[nSamples];
        for (int i = 0; i < nSamples; i++) {
            int idx = SAlags[i] - 1;
            if (idx >= 0 && idx < nlags) {
                Y[i] = ACFull[idx];
            }
        }

        double[] SBC = new double[effectiveOrders.length];
        Arrays.fill(SBC, Double.MAX_VALUE);

        for (int j = 0; j < effectiveOrders.length; j++) {
            int order = effectiveOrders[j];

            double[][] X = new double[nSamples][order];
            boolean validMatrix = true;
            for (int col = 0; col < order; col++) {
                for (int row = 0; row < nSamples; row++) {
                    int lagIdx = SAlags[row] + (col + 1) - 1;
                    if (lagIdx >= 0 && lagIdx < nlags) {
                        X[row][col] = ACFull[lagIdx];
                    } else {
                        validMatrix = false;
                    }
                }
            }

            if (!validMatrix) continue;

            double[] resid = regressResiduals(Y, X);
            if (resid != null) {
                double sse = 0.0;
                for (double r : resid) sse += r * r;
                if (sse > 0) {
                    SBC[j] = nSamples * Math.log(sse) - nSamples * Math.log(nSamples)
                            + Math.log(nSamples) * order;
                }
            }
        }

        int bestIdx = 0;
        double bestSBC = SBC[0];
        for (int j = 1; j < SBC.length; j++) {
            if (SBC[j] < bestSBC) {
                bestSBC = SBC[j];
                bestIdx = j;
            }
        }

        int bestOrder = effectiveOrders[bestIdx];
        return (int) FastMath.round(FastMath.log(2.0, (double) bestOrder));
    }

    private static double[] regressResiduals(double[] Y, double[][] X) {
        int n = Y.length;
        if (n == 0 || X.length == 0 || X[0].length == 0) return null;
        int p = X[0].length;
        try {
            RealMatrix xMatrix = MatrixUtils.createRealMatrix(n, p);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < p; j++) {
                    xMatrix.setEntry(i, j, X[i][j]);
                }
            }
            RealVector yVector = MatrixUtils.createRealVector(Y);

            QRDecomposition qr = new QRDecomposition(xMatrix);
            RealVector b = qr.getSolver().solve(yVector);

            RealVector fitted = xMatrix.operate(b);
            double[] residuals = new double[n];
            for (int i = 0; i < n; i++) {
                residuals[i] = Y[i] - fitted.getEntry(i);
            }
            return residuals;
        } catch (Exception e) {
            return null;
        }
    }

    public static AcfitResult kpcfit_sub_acfit(final double[] E, final double[] SA, final int[] SALags,
                                               final int J, int maxIterAC, int maxRunsAC, int maxResAC) {
        final double SCV = (E[1] - E[0] * E[0]) / (E[0] * E[0]);
        final double NSA = norm(SA, 2);
        final Random random = new Random();

        ArrayList<Double> fset = new ArrayList<Double>();
        ArrayList<double[]> xparamset = new ArrayList<double[]>();

        MultivariateFunction objective = new MultivariateFunction() {
            @Override
            public double value(double[] x) {
                double[] SCVj = new double[J];
                double[] G2j = new double[J];
                for (int j = 0; j < J; j++) {
                    SCVj[j] = x[j];
                    G2j[j] = x[J + j];
                }
                if (SCVj[0] < 0.5 - KPCFIT_TOL) return 1e10;
                for (int j = 1; j < J; j++) {
                    if (SCVj[j] < 1.0 + KPCFIT_TOL) return 1e10;
                }
                for (int j = 0; j < J; j++) {
                    if (G2j[j] < KPCFIT_TOL || G2j[j] > 1 - KPCFIT_TOL) return 1e10;
                }
                AcfitEvalResult er = kpcfit_sub_eval_acfit(SCVj, G2j, SALags);
                double[] diff = new double[SA.length];
                for (int i = 0; i < SA.length; i++) diff[i] = SA[i] - er.acfCoeff[i];
                return norm(diff, 1) / NSA + FastMath.pow(er.SCVcum - SCV, 2.0) / FastMath.pow(SCV, 2.0);
            }
        };

        for (int run = 0; run < maxRunsAC; run++) {
            double[] x0 = new double[2 * J];
            for (int j = 0; j < J; j++) {
                x0[j] = 1.0 + random.nextDouble();
                x0[J + j] = random.nextDouble();
            }

            try {
                SimplexOptimizer optimizer = new SimplexOptimizer(1e-8, 1e-8);
                NelderMeadSimplex simplex = new NelderMeadSimplex(2 * J, 0.1);

                PointValuePair result = optimizer.optimize(
                        new MaxEval(maxIterAC * 100),
                        new MaxIter(maxIterAC),
                        new ObjectiveFunction(objective),
                        GoalType.MINIMIZE,
                        new InitialGuess(x0),
                        simplex
                );
                xparamset.add(result.getPoint());
                fset.add(result.getValue());
            } catch (Exception e) {
                // Skip failed runs
            }
        }

        if (fset.isEmpty()) {
            double[] defaultSCV = new double[J];
            defaultSCV[0] = SCV;
            for (int i = 1; i < J; i++) defaultSCV[i] = 1.5;
            double[] defaultG2 = new double[J];
            Arrays.fill(defaultG2, 0.5);
            List<double[]> rsc = new ArrayList<double[]>();
            rsc.add(defaultSCV);
            List<double[]> rg = new ArrayList<double[]>();
            rg.add(defaultG2);
            return new AcfitResult(rsc, rg, new double[]{1e10});
        }

        Integer[] sortedIndices = new Integer[fset.size()];
        for (int i = 0; i < fset.size(); i++) sortedIndices[i] = i;
        final List<Double> fsetFinal = fset;
        Arrays.sort(sortedIndices, (a, b) -> Double.compare(fsetFinal.get(a), fsetFinal.get(b)));

        ArrayList<double[]> resSCV = new ArrayList<double[]>();
        ArrayList<double[]> resG2 = new ArrayList<double[]>();
        int total = Math.min(maxResAC, sortedIndices.length);
        double[] fobjAC = new double[total];

        for (int i = 0; i < total; i++) {
            double[] x = xparamset.get(sortedIndices[i]);
            double[] SCVj = new double[J];
            double[] G2j = new double[J];
            for (int j = 0; j < J; j++) {
                SCVj[j] = x[j];
                G2j[j] = Math.min(x[J + j], 1 - KPCFIT_TOL);
            }
            resSCV.add(SCVj);
            resG2.add(G2j);
            fobjAC[i] = fset.get(sortedIndices[i]);
        }

        return new AcfitResult(resSCV, resG2, fobjAC);
    }

    public static AcfitEvalResult kpcfit_sub_eval_acfit(double[] SCVj, double[] G2j, int[] lags) {
        int J = SCVj.length;
        double SCVcum = SCVj[0];
        double[] acfCoeff = new double[lags.length];
        for (int idx = 0; idx < lags.length; idx++) {
            acfCoeff[idx] = 0.5 * (1.0 - 1.0 / SCVj[0]) * FastMath.pow(G2j[0], (double) lags[idx]);
        }
        for (int j = 1; j < J; j++) {
            double SCVj_1 = SCVcum;
            SCVcum = (1.0 + SCVcum) * (1.0 + SCVj[j]) / 2.0 - 1.0;
            double r0j = 0.5 * (1.0 - 1.0 / SCVj[j]);
            for (int idx = 0; idx < lags.length; idx++) {
                double X = SCVj[j] * r0j * FastMath.pow(G2j[j], (double) lags[idx]);
                acfCoeff[idx] = (X + SCVj_1 * acfCoeff[idx] * (1.0 + X)) / SCVcum;
            }
        }
        return new AcfitEvalResult(SCVcum, acfCoeff);
    }

    public static BcfitResult kpcfit_sub_bcfit(final double[] E, final double[] SCVj, final double[] G2j,
                                               final double[] BC, final int[][] BCLags,
                                               int maxIterBC, int maxRunsBC) {
        final int NumMAPs = SCVj.length;
        final double TOL = 1e-9;
        final double EPSTOL = 10 * TOL;

        double NBC = norm(BC, 2);
        if (NBC < 1e-30) {
            double[] E1j = new double[NumMAPs];
            double[] E3j = new double[NumMAPs];
            for (int j = 0; j < NumMAPs; j++) {
                E1j[j] = FastMath.pow(E[0], 1.0 / NumMAPs);
                double E2j = (1 + SCVj[j]) * E1j[j] * E1j[j];
                E3j[j] = 1.501 * E2j * E2j / E1j[j];
            }
            return new BcfitResult(E1j, E3j, 0.0);
        }

        ArrayList<Integer> validIndices = new ArrayList<Integer>();
        for (int index = 0; index < BCLags.length; index++) {
            int[] lags = BCLags[index];
            boolean hasNegDiff = false;
            for (int i = 1; i < lags.length; i++) {
                if (lags[i] - lags[i - 1] < 0) {
                    hasNegDiff = true;
                    break;
                }
            }
            if (!hasNegDiff) validIndices.add(index);
        }
        final double[] filteredBC = new double[validIndices.size()];
        final int[][] filteredBCLags = new int[validIndices.size()][];
        for (int i = 0; i < validIndices.size(); i++) {
            filteredBC[i] = BC[validIndices.get(i)];
            filteredBCLags[i] = BCLags[validIndices.get(i)];
        }
        final double filteredNBC = norm(filteredBC, 2);
        if (filteredNBC < 1e-30) {
            double[] E1j = new double[NumMAPs];
            double[] E3j = new double[NumMAPs];
            for (int j = 0; j < NumMAPs; j++) {
                E1j[j] = FastMath.pow(E[0], 1.0 / NumMAPs);
                double E2j = (1 + SCVj[j]) * E1j[j] * E1j[j];
                E3j[j] = 1.501 * E2j * E2j / E1j[j];
            }
            return new BcfitResult(E1j, E3j, 0.0);
        }

        final double E1jBase = FastMath.pow(E[0], 1.0 / NumMAPs);
        double[] E1j0 = new double[NumMAPs];
        Arrays.fill(E1j0, E1jBase);

        final Random random = new Random();
        double tInit = E[0] * random.nextDouble();

        double[] E2j0 = new double[NumMAPs];
        double[] E3j0 = new double[NumMAPs];
        for (int j = 0; j < NumMAPs; j++) {
            E2j0[j] = (1 + SCVj[j]) * E1j0[j] * E1j0[j];
            E3j0[j] = (1.5 + tInit) * E2j0[j] * E2j0[j] / E1j0[j];
        }

        final double[] x0base = new double[2 * NumMAPs];
        for (int j = 0; j < NumMAPs; j++) {
            x0base[j] = E1j0[j];
            x0base[NumMAPs + j] = E3j0[j];
        }

        final double[] foldRef = new double[]{0.0};
        final double Eprime = E[0];

        // xtopar: extract E1, E3 from x; E1[0] = E[0]/prod(E1[1:])
        // objfun: compose, scale, compute BC, return diff norm
        MultivariateFunction objective = new MultivariateFunction() {
            @Override
            public double value(double[] x) {
                double[] e1 = new double[NumMAPs];
                double[] e3 = new double[NumMAPs];
                for (int idx = 0; idx < NumMAPs; idx++) {
                    e1[idx] = x[idx];
                    e3[idx] = x[NumMAPs + idx];
                }
                double prodE1 = 1.0;
                for (int j = 1; j < NumMAPs; j++) prodE1 *= e1[j];
                e1[0] = (prodE1 > 0) ? Eprime / prodE1 : E1jBase;

                for (int j = 0; j < NumMAPs; j++) {
                    if (e1[j] <= EPSTOL || e3[j] <= EPSTOL) return Math.max(2 * foldRef[0], 1e6);
                }
                double[] e2 = new double[NumMAPs];
                for (int j = 0; j < NumMAPs; j++) e2[j] = (1 + SCVj[j]) * e1[j] * e1[j];
                for (int j = 1; j < NumMAPs; j++) {
                    if ((2 + EPSTOL) * e1[j] * e1[j] > e2[j]) return Math.max(2 * foldRef[0], 1e6);
                    if ((1.5 + EPSTOL) * e2[j] * e2[j] / e1[j] > e3[j]) return Math.max(2 * foldRef[0], 1e6);
                }
                if (SCVj[0] > 1) {
                    if ((2 + EPSTOL) * e1[0] * e1[0] > e2[0]) return Math.max(2 * foldRef[0], 1e6);
                    if ((1.5 + EPSTOL) * e2[0] * e2[0] / e1[0] > e3[0]) return Math.max(2 * foldRef[0], 1e6);
                }

                ComposeResult cr = kpcfit_sub_compose(e1, SCVj, e3, G2j);
                if (cr.errorCode != 0 || cr.map == null) return Math.max(2 * foldRef[0], 1e6);

                MatrixCell scaledMAP = Map_scale.map_scale(cr.map, E[0]);
                double[] BCj = new double[filteredBCLags.length];
                for (int indexL = 0; indexL < filteredBCLags.length; indexL++) {
                    try {
                        BCj[indexL] = Map_joint.map_joint(scaledMAP, filteredBCLags[indexL], new int[]{1, 1, 1});
                    } catch (Exception ex) {
                        return Math.max(2 * foldRef[0], 1e6);
                    }
                }
                double[] diff = new double[filteredBC.length];
                for (int i = 0; i < filteredBC.length; i++) diff[i] = filteredBC[i] - BCj[i];
                double f = norm(diff, 2) / filteredNBC;
                if (Double.isNaN(f)) {
                    return 2 * foldRef[0];
                } else {
                    foldRef[0] = f;
                    return f;
                }
            }
        };

        double fBest = Double.MAX_VALUE;
        double[] xBest = null;

        for (int ind = 0; ind < maxRunsBC; ind++) {
            double[] x0;
            if (ind == 0) {
                x0 = x0base.clone();
            } else {
                x0 = new double[2 * NumMAPs];
                for (int j = 0; j < NumMAPs; j++) {
                    x0[j] = x0base[j] * (0.25 + 1.75 * random.nextDouble());
                }
                double tLoc = random.nextDouble() * E[0];
                for (int j = 0; j < NumMAPs; j++) {
                    double e2jLocal = (1 + SCVj[j]) * x0[j] * x0[j];
                    x0[NumMAPs + j] = (1.5 + tLoc) * e2jLocal * e2jLocal / x0[j];
                }
            }
            for (int j = 0; j < 2 * NumMAPs; j++) {
                if (x0[j] <= EPSTOL) x0[j] = EPSTOL;
            }

            try {
                SimplexOptimizer optimizer = new SimplexOptimizer(TOL, TOL);
                NelderMeadSimplex simplex = new NelderMeadSimplex(2 * NumMAPs, 0.01);
                PointValuePair result = optimizer.optimize(
                        new MaxEval(maxIterBC * 1000),
                        new MaxIter(maxIterBC),
                        new ObjectiveFunction(objective),
                        GoalType.MINIMIZE,
                        new InitialGuess(x0),
                        simplex
                );
                if (result.getValue() < fBest) {
                    fBest = result.getValue();
                    xBest = result.getPoint();
                }
            } catch (Exception e) {
                // Skip failed runs
            }
        }

        double[] xExtract = (xBest != null) ? xBest : x0base;
        double[] e1Out = new double[NumMAPs];
        double[] e3Out = new double[NumMAPs];
        for (int idx = 0; idx < NumMAPs; idx++) {
            e1Out[idx] = xExtract[idx];
            e3Out[idx] = xExtract[NumMAPs + idx];
        }
        double prodE1Out = 1.0;
        for (int j = 1; j < NumMAPs; j++) prodE1Out *= e1Out[j];
        e1Out[0] = (prodE1Out > 0) ? E[0] / prodE1Out : E1jBase;

        if (xBest != null) {
            return new BcfitResult(e1Out, e3Out, fBest);
        }
        return new BcfitResult(e1Out, e3Out, Double.MAX_VALUE);
    }

    public static ComposeResult kpcfit_sub_compose(double[] E1j, double[] SCVj, double[] E3j, double[] G2j) {
        int J = G2j.length;
        ArrayList<MatrixCell> subMAPs = new ArrayList<MatrixCell>();

        MatrixCell kpcMAP;
        try {
            double E2_1 = (1 + SCVj[0]) * E1j[0] * E1j[0];
            kpcMAP = MMPP.mmpp2_fit3(E1j[0], E2_1, E3j[0], G2j[0]);

            boolean hasImagOrNaN = hasNaNOrInfinite(kpcMAP);
            if (hasImagOrNaN || !isMapFeasible(kpcMAP)) {
                if (SCVj[0] < 0.5) {
                    kpcMAP = Map_erlang.map_erlang(E1j[0], 2);
                } else {
                    Ret.mamMAPFitReturn fitResult1 = Map2_fit.map2_fit(E1j[0], E2_1, -1.0, G2j[0]);
                    if (fitResult1.MAP != null && fitResult1.MAP.size() > 0 && (int) fitResult1.error == 0) {
                        kpcMAP = fitResult1.MAP;
                    } else {
                        Ret.mamMAPFitReturn fitResult2 = Map2_fit.map2_fit(E1j[0], E2_1, -1.0, 0.0);
                        if (fitResult2.MAP != null && fitResult2.MAP.size() > 0 && (int) fitResult2.error == 0) {
                            kpcMAP = fitResult2.MAP;
                        } else {
                            return new ComposeResult(null, subMAPs, 1);
                        }
                    }
                }
            }
        } catch (Exception e) {
            try {
                double E2_1 = (1 + SCVj[0]) * E1j[0] * E1j[0];
                if (SCVj[0] < 0.5) {
                    kpcMAP = Map_erlang.map_erlang(E1j[0], 2);
                } else {
                    Ret.mamMAPFitReturn fitResult1 = Map2_fit.map2_fit(E1j[0], E2_1, -1.0, G2j[0]);
                    if (fitResult1.MAP != null && fitResult1.MAP.size() > 0 && (int) fitResult1.error == 0) {
                        kpcMAP = fitResult1.MAP;
                    } else {
                        Ret.mamMAPFitReturn fitResult2 = Map2_fit.map2_fit(E1j[0], E2_1, -1.0, 0.0);
                        if (fitResult2.MAP != null && fitResult2.MAP.size() > 0 && (int) fitResult2.error == 0) {
                            kpcMAP = fitResult2.MAP;
                        } else {
                            return new ComposeResult(null, new ArrayList<MatrixCell>(), 1);
                        }
                    }
                }
            } catch (Exception e2) {
                return new ComposeResult(null, new ArrayList<MatrixCell>(), 1);
            }
        }

        subMAPs.add(kpcMAP);

        for (int j = 1; j < J; j++) {
            MatrixCell MAPj;
            try {
                double E2_j = (1 + SCVj[j]) * E1j[j] * E1j[j];
                Matrix[] mapArr = Map_feasblock.map_feasblock(E1j[j], E2_j, E3j[j], G2j[j]);
                MAPj = arrayToMatrixCell(mapArr);

                boolean feasible;
                try {
                    feasible = isMapFeasible(MAPj);
                } catch (Exception ex) {
                    feasible = false;
                }

                if (!feasible) {
                    if (SCVj[j] < 1.0) {
                        MAPj = Map_exponential.map_exponential(E1j[j]);
                    } else {
                        Ret.mamMAPFitReturn fitResult1 = Map2_fit.map2_fit(E1j[j], E2_j, -1.0, G2j[j]);
                        if ((int) fitResult1.error != 0 || fitResult1.MAP == null || fitResult1.MAP.size() == 0) {
                            Ret.mamMAPFitReturn fitResult2 = Map2_fit.map2_fit(E1j[j], E2_j, -1.0, 0.0);
                            if ((int) fitResult2.error != 0 || fitResult2.MAP == null || fitResult2.MAP.size() == 0) {
                                return new ComposeResult(null, subMAPs, 5);
                            }
                            MAPj = fitResult2.MAP;
                        } else {
                            MAPj = fitResult1.MAP;
                        }
                        try {
                            Matrix D0orig = MAPj.get(0);
                            int n = D0orig.getNumRows();
                            if (n == 2) {
                                Matrix negInvD0 = D0orig.inv();
                                for (int ii = 0; ii < n; ii++) {
                                    for (int jj = 0; jj < n; jj++) {
                                        negInvD0.set(ii, jj, -negInvD0.get(ii, jj));
                                    }
                                }
                                java.util.List<org.apache.commons.math3.complex.Complex> eigenValuesList = negInvD0.eig();
                                org.apache.commons.math3.complex.Complex[] eigenValues =
                                        (eigenValuesList == null) ? null
                                                : eigenValuesList.toArray(new org.apache.commons.math3.complex.Complex[0]);
                                if (eigenValues != null && eigenValues.length >= n) {
                                    Matrix D0new = new Matrix(n, n);
                                    for (int ii = 0; ii < n; ii++) {
                                        D0new.set(ii, ii, -1.0 / eigenValues[ii].getReal());
                                    }
                                    Matrix P = Map_embedded.map_embedded(MAPj);
                                    java.util.List<org.apache.commons.math3.complex.Complex> pEigList = P.eig();
                                    org.apache.commons.math3.complex.Complex[] pEig =
                                            (pEigList == null) ? null
                                                    : pEigList.toArray(new org.apache.commons.math3.complex.Complex[0]);
                                    if (pEig != null && pEig.length >= n) {
                                        double pMinReal = pEig[0].getReal();
                                        for (int ii = 1; ii < pEig.length; ii++) {
                                            if (pEig[ii].getReal() < pMinReal) pMinReal = pEig[ii].getReal();
                                        }
                                        Matrix D1new = new Matrix(n, n);
                                        for (int ii = 0; ii < n; ii++) {
                                            D1new.set(ii, 0, -D0new.get(ii, ii) * pMinReal);
                                            D1new.set(ii, 1, -D0new.get(ii, ii) * (1 - pMinReal));
                                        }
                                        MatrixCell newCell = new MatrixCell(2);
                                        newCell.set(0, D0new);
                                        newCell.set(1, D1new);
                                        MAPj = Map_normalize.map_normalize(newCell);
                                    }
                                }
                            }
                        } catch (Exception ex) {
                            // ignore
                        }
                    }
                }
            } catch (Exception e) {
                MAPj = Map_exponential.map_exponential(E1j[j]);
            }

            if (MAPj == null || MAPj.size() == 0) {
                MAPj = Map_exponential.map_exponential(E1j[j]);
            }
            subMAPs.add(MAPj);

            Matrix[] kpcArr = Map_kpc.map_kpc(matrixCellToArray(kpcMAP), matrixCellToArray(MAPj));
            kpcMAP = arrayToMatrixCell(kpcArr);
        }

        kpcMAP = Map_normalize.map_normalize(kpcMAP);

        int error = 0;
        for (MatrixCell subMAP : subMAPs) {
            if (!isMapFeasible(subMAP)) {
                error = 10;
            }
        }

        return new ComposeResult(kpcMAP, subMAPs, error);
    }

    private static boolean hasNaNOrInfinite(MatrixCell MAP) {
        try {
            Matrix D0 = MAP.get(0);
            Matrix D1 = MAP.get(1);
            for (int i = 0; i < D0.getNumRows(); i++) {
                for (int j = 0; j < D0.getNumCols(); j++) {
                    double v = D0.get(i, j);
                    if (Double.isNaN(v) || Double.isInfinite(v)) return true;
                }
            }
            for (int i = 0; i < D1.getNumRows(); i++) {
                for (int j = 0; j < D1.getNumCols(); j++) {
                    double v = D1.get(i, j);
                    if (Double.isNaN(v) || Double.isInfinite(v)) return true;
                }
            }
        } catch (Exception e) {
            return true;
        }
        return false;
    }

    private static boolean isMapFeasible(MatrixCell MAP) {
        try {
            return Map_isfeasible.map_isfeasible(MAP);
        } catch (Exception e) {
            return false;
        }
    }

    private static double[] evaluateObjFunction(MatrixCell map, double[] E, double[] AC, int[] ACLags,
                                                double[] BC, int[][] BCLags) {
        double tSCV = (E[1] - E[0] * E[0]) / (E[0] * E[0]);

        double[] mapACF = new double[ACLags.length];
        for (int i = 0; i < ACLags.length; i++) {
            try {
                Matrix lagMatrix = Matrix.singleton((double) ACLags[i]);
                Matrix acfMatrix = Map_acf.map_acf(map, lagMatrix);
                mapACF[i] = (acfMatrix.getNumElements() > 0) ? acfMatrix.get(0) : 0.0;
            } catch (Exception e) {
                mapACF[i] = 0.0;
            }
        }
        double[] diffAC = new double[AC.length];
        for (int i = 0; i < AC.length; i++) diffAC[i] = AC[i] - mapACF[i];
        double objAC = norm(diffAC, 1) / norm(AC, 2)
                + FastMath.pow(Map_scv.map_scv(map) - tSCV, 2.0) / FastMath.pow(tSCV, 2.0);

        double[] mapBC = new double[BCLags.length];
        for (int i = 0; i < BCLags.length; i++) {
            try {
                mapBC[i] = Map_joint.map_joint(map, BCLags[i], new int[]{1, 1, 1});
            } catch (Exception e) {
                mapBC[i] = 1.0;
            }
        }
        double[] diffBC = new double[BC.length];
        for (int i = 0; i < BC.length; i++) diffBC[i] = BC[i] - mapBC[i];
        double objBC = norm(diffBC, 2) / norm(BC, 2);

        return new double[]{objAC, objBC};
    }

    private static double norm(double[] v, int p) {
        if (p == 1) {
            double s = 0.0;
            for (double x : v) s += Math.abs(x);
            return s;
        } else if (p == 2) {
            double s = 0.0;
            for (double x : v) s += x * x;
            return Math.sqrt(s);
        } else {
            double s = 0.0;
            for (double x : v) s += FastMath.pow(Math.abs(x), p);
            return FastMath.pow(s, 1.0 / p);
        }
    }

    public static double[] kpcfit_hyper_charpoly(double[] E, int n) {
        double[] Ep = new double[1 + E.length];
        Ep[0] = 1.0;
        for (int i = 0; i < E.length; i++) Ep[i + 1] = E[i];

        double[] f = new double[2 * n];
        for (int i = 0; i < 2 * n; i++) f[i] = factorial(i);

        double[][] A = new double[n + 1][n + 1];
        for (int i = 1; i <= n; i++) {
            for (int col = 0; col <= n; col++) {
                int epIdx = (n + i) - col;
                int fIdx = (n + i) - col;
                int epIndex = epIdx - 1;
                int fIndex = fIdx - 1;
                if (epIndex >= 0 && epIndex < Ep.length && fIndex >= 0 && fIndex < f.length) {
                    A[i - 1][col] = Ep[epIndex] / f[fIndex];
                }
            }
        }
        for (int col = 0; col <= n; col++) A[n][col] = 0.0;
        A[n][0] = 1.0;

        double[] b = new double[n + 1];
        b[n] = 1.0;

        RealMatrix aMatrix = MatrixUtils.createRealMatrix(n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                aMatrix.setEntry(i, j, A[i][j]);
            }
        }
        RealVector bVector = MatrixUtils.createRealVector(b);
        QRDecomposition qr = new QRDecomposition(aMatrix);
        RealVector mVector = qr.getSolver().solve(bVector);

        double[] m = new double[n + 1];
        for (int i = 0; i <= n; i++) m[i] = mVector.getEntry(i);
        return m;
    }

    public static MatrixCell kpcfit_ph_prony(double[] E, int n) {
        double[] f = new double[2 * n];
        for (int i = 0; i < 2 * n; i++) f[i] = factorial(i);

        double[] m = kpcfit_hyper_charpoly(E, n);

        double[] mReversed = new double[m.length];
        for (int i = 0; i < m.length; i++) mReversed[i] = m[m.length - 1 - i];
        org.apache.commons.math3.complex.Complex[] thetaComplex = Maths.roots(mReversed);

        double[] theta = new double[n];
        for (int i = 0; i < n; i++) theta[i] = thetaComplex[i].getReal();

        double[][] C = new double[n][n];
        for (int i = 1; i <= n; i++) {
            for (int j = 0; j < n; j++) {
                C[i - 1][j] = f[i] * FastMath.pow(theta[j], (double) i);
            }
        }

        RealMatrix cMatrix = MatrixUtils.createRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cMatrix.setEntry(i, j, C[i][j]);
            }
        }
        double[] eArr = new double[n];
        for (int i = 0; i < n; i++) eArr[i] = E[i];
        RealVector eVector = MatrixUtils.createRealVector(eArr);
        QRDecomposition qr = new QRDecomposition(cMatrix);
        RealVector mSolve = qr.getSolver().solve(eVector);

        double[] M = new double[n];
        for (int i = 0; i < n; i++) M[i] = mSolve.getEntry(i);

        Matrix D0 = new Matrix(n, n);
        for (int i = 0; i < n; i++) D0.set(i, i, -1.0 / theta[i]);

        Matrix D1 = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D1.set(i, j, (1.0 / theta[i]) * M[j]);
            }
        }

        MatrixCell PH = new MatrixCell(2);
        PH.set(0, D0);
        PH.set(1, D1);
        return PH;
    }

    public static class KPCFitPhOptions {
        public boolean verbose = true;
        public int runs = 5;
        public int minNumStates = 2;
        public int maxNumStates = 32;
        public int minExactMom = 3;

        public KPCFitPhOptions() {}

        public KPCFitPhOptions(boolean verbose, int runs, int minNumStates, int maxNumStates, int minExactMom) {
            this.verbose = verbose;
            this.runs = runs;
            this.minNumStates = minNumStates;
            this.maxNumStates = maxNumStates;
            this.minExactMom = minExactMom;
        }
    }

    public static KPCFitPhOptions kpcfit_ph_options(double[] E, boolean verbose, int runs, int minNumStates,
                                                    int maxNumStates, int minExactMom) {
        KPCFitPhOptions options = new KPCFitPhOptions(verbose, runs, minNumStates, maxNumStates, minExactMom);

        if (options.runs < 1) options.runs = 1;
        if (options.minNumStates < 2) options.minNumStates = 2;
        if (options.maxNumStates < 2) options.maxNumStates = 2;
        if (options.minExactMom < 1) options.minExactMom = 1;
        if (options.minExactMom > E.length) options.minExactMom = E.length;

        if ((options.minNumStates & (options.minNumStates - 1)) != 0) {
            options.minNumStates = nextPowerOf2(options.minNumStates);
            if (options.verbose) {
                System.err.println("Warning: MinNumStates not a power of 2, fixed to " + options.minNumStates + ".");
            }
        }

        if ((options.maxNumStates & (options.maxNumStates - 1)) != 0) {
            options.maxNumStates = nextPowerOf2(options.maxNumStates);
            if (options.verbose) {
                System.err.println("Warning: MaxNumStates not a power of 2, fixed to " + options.maxNumStates + ".");
            }
        }

        if (2 * options.maxNumStates - 1 > E.length) {
            // largest power-of-2 number of states fittable from the supplied
            // moments (matching MATLAB kpcfit_ph_options 0.4.0) instead of throwing
            int maxFeasible = 1 << (int) Math.floor(Math.log(Math.max(1, (E.length + 1) / 2)) / Math.log(2));
            if (options.verbose) {
                System.err.println("Warning: MaxNumStates of " + options.maxNumStates + " requires "
                        + (2 * options.maxNumStates - 1) + " moments but only " + E.length
                        + " supplied; reduced to " + maxFeasible + ".");
            }
            options.maxNumStates = maxFeasible;
            if (options.minNumStates > maxFeasible) {
                options.minNumStates = maxFeasible;
            }
        }

        if (options.minNumStates > options.maxNumStates) {
            if (options.verbose) {
                System.err.println("Warning: MaxNumStates < MinNumStates, fixed.");
            }
            int tmp = options.minNumStates;
            options.minNumStates = options.maxNumStates;
            options.maxNumStates = tmp;
        }

        return options;
    }

    public static KPCFitPhOptions kpcfit_ph_options(double[] E) {
        return kpcfit_ph_options(E, true, 5, 2, 32, 3);
    }

    private static int nextPowerOf2(int n) {
        int v = n - 1;
        v = v | (v >> 1);
        v = v | (v >> 2);
        v = v | (v >> 4);
        v = v | (v >> 8);
        v = v | (v >> 16);
        return v + 1;
    }

    private static double factorial(int n) {
        if (n < 0) return 1.0;
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) result *= (double) i;
        return result;
    }

    public static List<MatrixCell> kpcfit_ph_exact(double[] E, KPCFitPhOptions options) {
        ArrayList<MatrixCell> phExact = new ArrayList<MatrixCell>();
        double SCV = (E[1] - E[0] * E[0]) / (E[0] * E[0]);

        if (SCV > 1.0) {
            if (options.verbose) {
                System.out.println("kpcfit_ph: HIGHER variability than an exponential (var/mean^2 = " + SCV + ")");
                System.out.println();
                System.out.println("kpcfit_ph: starting exact hyper-exponential fitting method (Prony's method)");
            }
            for (int n = 2; n <= options.maxNumStates; n++) {
                if (E.length < 2 * n - 1) {
                    if (options.verbose) {
                        System.out.println("kpcfit_ph: not enough moments given in input to fit hyper-exp(" + n + ")");
                    }
                    break;
                }
                try {
                    MatrixCell PH = kpcfit_ph_prony(E, n);
                    if (Map_isfeasible.map_isfeasible(PH)) {
                        phExact.add(PH);
                        if (options.verbose) {
                            System.out.println("\t\t\thyper-exp(" + n + "): feasible, matched exactly " + (2 * n - 1) + " moments. result saved.");
                        }
                    } else {
                        if (options.verbose) {
                            System.out.println("\t\t\thyper-exp(" + n + "): infeasible to fit exactly.");
                        }
                        break;
                    }
                } catch (Exception e) {
                    if (options.verbose) {
                        System.out.println("\t\t\thyper-exp(" + n + "): infeasible to fit exactly.");
                    }
                    break;
                }
            }
        } else if (SCV < 1.0) {
            if (options.verbose) {
                System.out.println("kpcfit_ph: LOWER variability than an exponential (var/mean^2 = " + SCV + ")");
                System.out.println();
            }
            int n = 1;
            while (1.0 / n > SCV) n++;
            if (options.verbose) {
                System.out.println("kpcfit_ph: exact fitting of E[X^2] requires at least " + n + " states");
            }
            if (options.minExactMom >= 2 && options.maxNumStates < n) {
                if (options.verbose) {
                    System.out.println("kpcfit_ph: impossible to fit exactly E[X^2] with MaxNumStates = "
                            + options.maxNumStates + ", increasing to MaxNumStates = " + n + ".");
                }
                return phExact;
            }

            if (n == 2) {
                if (options.verbose) {
                    System.out.println("kpcfit_ph: attempting PH(2) fitting method");
                }
                Ret.mamMAPFitReturn fitResult = Map2_fit.map2_fit(E[0], E[1], E[2], 0.0);
                if (fitResult.MAP == null || fitResult.MAP.size() == 0) {
                    Ret.mamMAPFitReturn fitResult2 = Map2_fit.map2_fit(E[0], E[1], -1.0, 0.0);
                    if (fitResult2.MAP != null && fitResult2.MAP.size() > 0
                            && Map_isfeasible.map_isfeasible(fitResult2.MAP)) {
                        phExact.add(fitResult2.MAP);
                        if (options.verbose) {
                            System.out.println("\t\t\tph(2): feasible, matched exactly 2 moments. result saved.");
                        }
                    } else {
                        throw new IllegalStateException("anomalous set of moments, please check.");
                    }
                } else {
                    phExact.add(fitResult.MAP);
                    if (options.verbose) {
                        System.out.println("\t\t\tph(2): feasible, matched exactly 3 moments. result saved.");
                    }
                }
            } else if (Math.abs(SCV - 1.0 / n) < KPCFIT_TOL) {
                MatrixCell ERL = Map_erlang.map_erlang(E[0], n);
                double momDist = 0.0;
                for (int k = 1; k <= E.length; k++) {
                    double diff = E[k - 1] - Map_moment.map_moment(ERL, k);
                    momDist += diff * diff;
                }
                if (Math.sqrt(momDist) < KPCFIT_TOL) {
                    if (options.verbose) {
                        System.out.println("kpcfit_ph: erlang moment set. fitted erlang-" + n + ". result saved.");
                    }
                    phExact.add(ERL);
                }
            } else {
                int maxorder = options.maxNumStates;
                if (options.verbose) {
                    System.out.println("kpcfit_ph: fitting APH distribution (best effort, max order = " + maxorder + ").");
                }
                try {
                    MatrixCell PH = Aph_fit.aph_fit(E[0], E[1], E[2], maxorder);
                    if (Map_isfeasible.map_isfeasible(PH)) {
                        int aphMatched = 0;
                        int phSize = PH.get(0).getNumRows();
                        for (int k = 1; k <= 2 * phSize - 1; k++) {
                            if (k <= E.length) {
                                double momK = Map_moment.map_moment(PH, k);
                                if (Math.abs(E[k - 1] - momK) < KPCFIT_TOL * momK) {
                                    aphMatched++;
                                }
                            }
                        }
                        if (options.verbose) {
                            System.out.println("\t\t\t      aph(" + phSize + "): feasible, matched exactly " + aphMatched + " moments. result saved.");
                        }
                        phExact.add(PH);
                    } else {
                        if (options.verbose) {
                            System.out.println("kpcfit_ph: cannot fit APH distribution.");
                        }
                    }
                } catch (Exception e) {
                    if (options.verbose) {
                        System.out.println("kpcfit_ph: cannot fit APH distribution.");
                    }
                }
            }
        } else {
            if (options.verbose) {
                System.out.println("kpcfit_ph: SAME variability as an exponential (var/mean^2 = " + SCV + ")");
                System.out.println();
            }
            MatrixCell EXP = new MatrixCell(2);
            Matrix d0 = new Matrix(1, 1);
            d0.set(0, 0, -1.0 / E[0]);
            EXP.set(0, d0);
            Matrix d1 = new Matrix(1, 1);
            d1.set(0, 0, 1.0 / E[0]);
            EXP.set(1, d1);
            double momDist = 0.0;
            for (int k = 1; k <= E.length; k++) {
                double diff = E[k - 1] - Map_moment.map_moment(EXP, k);
                momDist += diff * diff;
            }
            if (Math.sqrt(momDist) < KPCFIT_TOL) {
                if (options.verbose) {
                    System.out.println("kpcfit_ph: exponential moment set. fitted exponential. result saved.");
                }
            }
            phExact.add(EXP);
        }

        return phExact;
    }

    // __TRANSLATION_CONTINUES_PART_4__
}
