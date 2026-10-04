package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.m3pp.M3pp22_fitc_approx_cov;
import jline.api.mam.m3pp.M3pp2m_fitc;
import jline.api.mam.m3pp.M3pp2m_fitc_approx;
import jline.api.mam.m3pp.M3pp2m_fitc_approx_ag_multiclass;
import jline.api.mam.Mmpp2_fitc;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Compresses an MMAP using various approximation methods.
 */
public final class Mmap_compress {
    private Mmap_compress() {}

    public static Matrix[] mmap_compress(Matrix[] mmap) {
        return mmap_compress(mmap, "default");
    }

    public static Matrix[] mmap_compress(Matrix[] mmap, String method) {
        if (mmap.length < 2) throw new IllegalArgumentException("MMAP must have at least D0 and D1 matrices");
        int m = mmap.length - 2;
        if (m == 0) return new Matrix[]{mmap[0], mmap[1], mmap[1]};
        String lc = method.toLowerCase();
        Matrix[] result;
        if ("default".equals(lc) || "mixture".equals(lc) || "mixture.order1".equals(lc)) result = compressMixtureOrder1(mmap);
        else if ("mixture.order2".equals(lc)) result = compressMixtureOrder2(mmap);
        else if ("mamap2".equals(lc)) result = compressMAMAP2(mmap);
        else if ("mamap2.fb".equals(lc)) result = compressMAMAP2FB(mmap);
        else if ("m3pp.approx_cov".equals(lc)) result = compressM3PPApproxCov(mmap);
        else if ("m3pp.approx_ag".equals(lc)) result = compressM3PPApproxAG(mmap);
        else if ("m3pp.exact_delta".equals(lc)) result = compressM3PPExactDelta(mmap);
        else if ("m3pp.approx_delta".equals(lc)) result = compressM3PPApproxDelta(mmap);
        else throw new IllegalArgumentException("Unknown compression method: " + method);
        // MATLAB mmap_compress.m:74 normalizes on every return path; the JAR used to skip this.
        MatrixCell norm = Mmap_normalize.mmap_normalize(new MatrixCell(result));
        Matrix[] out = new Matrix[norm.size()];
        for (int i = 0; i < norm.size(); i++) out[i] = norm.get(i);
        return out;
    }

    /**
     * Order-1 mixture (M3A). Builds K components, one per class, and recombines them with
     * {@link Mmap_mixture} using the class probabilities p_c as mixing weights.
     *
     * <p>Component c must carry the law of the inter-arrival time <em>conditioned on the arrival
     * that ends it being of class c</em>: mmap_mixture marks the arrival leaving component c with
     * class c, so the class of an arrival and the interval preceding it are both governed by the
     * component active during that interval. That conditional law is the class-c BACKWARD moment
     * set B(c,1:3) = E[T^k | class of the ending arrival = c]. It is NOT the forward moment
     * (which conditions on the class of the STARTING arrival) and it is NOT the class-c marginal
     * MAP from mmap_maps (whose mean is 1/lambda_c). This mirrors 'mixture.order2', which
     * conditions on the (last,next) class pair and fits the cross moments; order 1 drops the
     * "last" index.
     *
     * <p>Preserved exactly: aggregate moments 1..3 via the M3A mixture law M_k = sum_c B(k,c)*p_c
     * (M1 always, since aph2_adjust never alters M1; M2/M3 when APH(2)-feasible); the class
     * probabilities p_c and hence the per-class rates lambda_c = p_c/M1; marking consistency
     * D1 = sum_c D1^(c); and MAP feasibility.
     *
     * <p>Lost by construction: every autocorrelation. Each component is re-entered at its map_pie
     * on every arrival, so the intervals are i.i.d. and the result is a RENEWAL process
     * (acf -> 0, IDC -> the SCV-determined renewal value, class sequence i.i.d.).
     */
    private static Matrix[] compressMixtureOrder1(Matrix[] mmap) {
        int m = mmap.length - 2;
        MatrixCell cell = new MatrixCell(mmap);
        Matrix pc = Mmap_pc.mmap_pc(cell);
        Matrix orders = new Matrix(1, 3, 3);
        orders.set(0, 0, 1.0);
        orders.set(0, 1, 2.0);
        orders.set(0, 2, 3.0);
        Matrix B = Mmap_backward_moment.mmap_backward_moment(cell, orders, 1);
        Map<Integer, MatrixCell> maps = new HashMap<Integer, MatrixCell>();
        for (int c = 0; c < m; c++) {
            if (pc.get(c) <= GlobalConstants.Zero) {
                // Class c never arrives, so B(c,:) is an 0/0 normalization. The component
                // carries zero mixture weight: any proper MAP leaves the result unchanged.
                maps.put(c, Map_exponential.map_exponential(1.0));
            } else {
                maps.put(c, Aph2_fit.aph2_fit(B.get(c, 0), B.get(c, 1), B.get(c, 2)).APH);
            }
        }
        MatrixCell res = Mmap_mixture.mmap_mixture(pc, maps);
        Matrix[] result = new Matrix[res.size()];
        for (int i = 0; i < res.size(); i++) result[i] = res.get(i);
        return result;
    }

    private static Matrix[] compressMixtureOrder2(Matrix[] mmap) { return mmap; }

    private static Matrix[] compressMAMAP2(Matrix[] mmap) {
        Matrix[] aggregateMap = new Matrix[]{mmap[0], mmap[1]};
        double M1 = Map_moment.map_moment(new MatrixCell(aggregateMap), 1);
        double M2 = Map_moment.map_moment(new MatrixCell(aggregateMap), 2);
        double M3 = Map_moment.map_moment(new MatrixCell(aggregateMap), 3);
        Matrix acfLags = new Matrix(1, 1);
        acfLags.set(0, 0, 1.0);
        Matrix acf = Map_acf.map_acf(new MatrixCell(aggregateMap), acfLags);
        double GAMMA = acf.get(0, 0);
        Matrix pie = Mmap_pie.mmap_pie(new MatrixCell(mmap));
        double[] P = pie.toArray1D();
        Matrix momMat = new Matrix(1, 1);
        momMat.set(0, 0, 1.0);
        double[] F = Mmap_forward_moment.mmap_forward_moment(new MatrixCell(mmap), momMat).toArray1D();
        double[] B = Mmap_backward_moment.mmap_backward_moment(new MatrixCell(mmap), momMat).toArray1D();
        Matrix S = Mmap_sigma.mmap_sigma(new MatrixCell(mmap));
        return Mamap2m_fit.mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S);
    }

    private static Matrix[] compressMAMAP2FB(Matrix[] mmap) {
        MatrixCell result = Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb_mmap(new MatrixCell(mmap));
        Matrix[] arr = new Matrix[result.size()];
        for (int i = 0; i < result.size(); i++) arr[i] = result.get(i);
        return arr;
    }

    private static Matrix[] compressM3PPApproxCov(Matrix[] mmap) {
        int m = mmap.length - 2;
        if (m != 2) throw new IllegalArgumentException("m3pp.approx_cov method only supports 2 classes, found " + m);
        Matrix[] aggregateMap = new Matrix[]{mmap[0], mmap[1]};
        double moments1 = Map_moment.map_moment(new MatrixCell(aggregateMap), 1);
        double a = 1.0 / moments1;
        double t1 = 1.0, t2 = 10.0, t3 = 100.0;
        double v1 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t1);
        double v2 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t2);
        double vinf = Map_count_var.map_count_var(new MatrixCell(aggregateMap), GlobalConstants.Inf);
        double bt1 = v1 / (a * t1);
        double bt2 = v2 / (a * t2);
        double binf = Double.isFinite(vinf) ? vinf / (a * 1000.0) : 1.0;
        double m3t2 = Map_count_moment.map_count_moment(new MatrixCell(aggregateMap), t2, 3);
        Matrix pie = Mmap_pie.mmap_pie(new MatrixCell(mmap));
        double[] ai = new double[]{pie.get(0, 0) * a, pie.get(0, 1) * a};
        Matrix cov = Mmap_count_mcov.mmap_count_mcov(mmap, t3);
        double ct3 = cov.get(0, 1);
        return M3pp22_fitc_approx_cov.m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, ct3, t3);
    }

    private static Matrix[] compressM3PPApproxAG(Matrix[] mmap) {
        int m = mmap.length - 2;
        Matrix[] aggregateMap = new Matrix[]{mmap[0], mmap[1]};
        double moments1 = Map_moment.map_moment(new MatrixCell(aggregateMap), 1);
        double a = 1.0 / moments1;
        double t1 = 1.0, t2 = 10.0, t3 = 100.0;
        double v1 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t1);
        double v2 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t2);
        double vinf = Map_count_var.map_count_var(new MatrixCell(aggregateMap), GlobalConstants.Inf);
        double bt1 = v1 / (a * t1);
        double bt2 = v2 / (a * t2);
        double binf = Double.isFinite(vinf) ? vinf / (a * 1000.0) : 1.0;
        double m3t2 = Map_count_moment.map_count_moment(new MatrixCell(aggregateMap), t2, 3);
        Matrix pie = Mmap_pie.mmap_pie(new MatrixCell(mmap));
        double[] ai = new double[m];
        for (int i = 0; i < m; i++) ai[i] = pie.get(0, i) * a;
        Matrix variances = Mmap_count_var.mmap_count_var(new MatrixCell(mmap), t3);
        Matrix covariances = Mmap_count_mcov.mmap_count_mcov(mmap, t3);
        double[] gt3 = new double[m];
        for (int i = 0; i < m; i++) {
            double sum = variances.get(i, 0);
            for (int j = 0; j < m; j++) {
                if (i != j) sum += covariances.get(i, j);
            }
            gt3[i] = sum;
        }
        Matrix[] mmpp = Mmpp2_fitc.mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2);
        return M3pp2m_fitc_approx_ag_multiclass.m3pp2m_fitc_approx_ag_multiclass(mmpp, ai, gt3, t3);
    }

    private static Matrix[] compressM3PPExactDelta(Matrix[] mmap) {
        int m = mmap.length - 2;
        Matrix[] aggregateMap = new Matrix[]{mmap[0], mmap[1]};
        double moments1 = Map_moment.map_moment(new MatrixCell(aggregateMap), 1);
        double a = 1.0 / moments1;
        double t1 = 1.0, t2 = 10.0, t3 = 100.0;
        double v1 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t1);
        double v2 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t2);
        double vinf = Map_count_var.map_count_var(new MatrixCell(aggregateMap), GlobalConstants.Inf);
        double bt1 = v1 / (a * t1);
        double bt2 = v2 / (a * t2);
        double binf = Double.isFinite(vinf) ? vinf / (a * 1000.0) : 1.0;
        double m3t2 = Map_count_moment.map_count_moment(new MatrixCell(aggregateMap), t2, 3);
        Matrix pie = Mmap_pie.mmap_pie(new MatrixCell(mmap));
        double[] ai = new double[m];
        for (int i = 0; i < m; i++) ai[i] = pie.get(0, i) * a;
        Matrix variances = Mmap_count_var.mmap_count_var(new MatrixCell(mmap), t3);
        double totalVar = 0.0;
        for (int i = 0; i < m; i++) totalVar += variances.get(i, 0);
        double[] dvt3 = new double[m];
        for (int i = 0; i < m; i++) dvt3[i] = variances.get(i, 0) - (totalVar - variances.get(i, 0));
        return M3pp2m_fitc.m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
    }

    private static Matrix[] compressM3PPApproxDelta(Matrix[] mmap) {
        int m = mmap.length - 2;
        Matrix[] aggregateMap = new Matrix[]{mmap[0], mmap[1]};
        double moments1 = Map_moment.map_moment(new MatrixCell(aggregateMap), 1);
        double a = 1.0 / moments1;
        double t1 = 1.0, t2 = 10.0, t3 = 100.0;
        double v1 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t1);
        double v2 = Map_count_var.map_count_var(new MatrixCell(aggregateMap), t2);
        double vinf = Map_count_var.map_count_var(new MatrixCell(aggregateMap), GlobalConstants.Inf);
        double bt1 = v1 / (a * t1);
        double bt2 = v2 / (a * t2);
        double binf = Double.isFinite(vinf) ? vinf / (a * 1000.0) : 1.0;
        double m3t2 = Map_count_moment.map_count_moment(new MatrixCell(aggregateMap), t2, 3);
        Matrix pie = Mmap_pie.mmap_pie(new MatrixCell(mmap));
        double[] ai = new double[m];
        for (int i = 0; i < m; i++) ai[i] = pie.get(0, i) * a;
        Matrix variances = Mmap_count_var.mmap_count_var(new MatrixCell(mmap), t3);
        double totalVar = 0.0;
        for (int i = 0; i < m; i++) totalVar += variances.get(i, 0);
        double[] dvt3 = new double[m];
        for (int i = 0; i < m; i++) dvt3[i] = variances.get(i, 0) - (totalVar - variances.get(i, 0));
        return M3pp2m_fitc_approx.m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
    }

    /** MMAP compress algorithms. */
    public static final class MmapCompressAlgo {}
}
