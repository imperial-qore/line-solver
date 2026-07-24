/**
 * @file Gauss-Legendre integration for multi-class repairman model normalizing constants
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mmint2_gausslegendre {
    private Pfqn_mmint2_gausslegendre() {}

    /**
     * Compute the normalizing constant of a repairmen model using Gauss-Legendre integration.
     */
    public static Ret.pfqnNc pfqn_mmint2_gausslegendre(Matrix L, Matrix N, Matrix Z, Integer m) {
        int mVal = (m == null) ? 1 : m.intValue();

        List<Double> gausslegendreNodes = new ArrayList<Double>();
        List<Double> gausslegendreWeights = new ArrayList<Double>();

        try {
            InputStream nodeStream = Pfqn_mmint2_gausslegendre.class.getResourceAsStream("/gausslegendre-nodes.txt");
            if (nodeStream == null) {
                throw new FileNotFoundException("Resource gausslegendre-nodes.txt not found.");
            }
            BufferedReader nodeReader = new BufferedReader(new InputStreamReader(nodeStream));
            String line;
            while ((line = nodeReader.readLine()) != null) {
                gausslegendreNodes.add(Double.parseDouble(line));
            }

            InputStream weightStream = Pfqn_mmint2_gausslegendre.class.getResourceAsStream("/gausslegendre-weights.txt");
            BufferedReader weightReader = new BufferedReader(new InputStreamReader(weightStream));
            while ((line = weightReader.readLine()) != null) {
                gausslegendreWeights.add(Double.parseDouble(line));
            }
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        int n = (int) FastMath.max(300.0,
                FastMath.min((double) gausslegendreNodes.size(),
                        2 * (N.sumRows().sumCols().get(0) + mVal - 1) - 1));
        Matrix y = new Matrix(1, n);
        y.fill(0.0);

        if (!(Z.getNumRows() == L.getNumRows() && Z.getNumCols() == L.getNumCols())) {
            throw new RuntimeException("The dimensions of Z and L are not the same.");
        }
        for (int i = 0; i < n; i++) {
            Matrix tmp = L.copy();
            for (int j = 0; j < tmp.getNumRows(); j++) {
                for (int k = 0; k < tmp.getNumCols(); k++) {
                    tmp.set(j, k, FastMath.log(Z.get(j, k) + gausslegendreNodes.get(i) * tmp.get(j, k)));
                }
            }
            y.set(i, (N.mult(tmp.transpose())).value());
        }

        Matrix g = y.copy();
        Matrix nodes = new Matrix(1, n);
        Matrix logNodes = new Matrix(1, n);
        Matrix logWeights = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            nodes.set(i, gausslegendreNodes.get(i));
            logNodes.set(i, FastMath.log(gausslegendreNodes.get(i)));
            logWeights.set(i, FastMath.log(gausslegendreWeights.get(i)));
        }

        for (int i = 0; i < n; i++) {
            g.set(i, g.get(i) + logWeights.get(i) - nodes.get(i));
        }

        double coeff = 0.0;
        for (int i = 0; i < N.length(); i++) {
            coeff -= Maths.factln(N.get(i));
        }
        coeff -= Maths.factln(mVal - 1);
        coeff += (mVal - 1) * logNodes.elementSum();

        double lG = 0.0;
        for (int i = 0; i < g.length(); i++) {
            lG += FastMath.exp(g.get(i));
        }
        lG = FastMath.log(lG) + coeff;
        if (!Double.isFinite(lG)) {
            lG = Matrix.logsumexp(g) + coeff;
        }
        double G = FastMath.exp(lG);
        return new Ret.pfqnNc(G, lG);
    }
}
