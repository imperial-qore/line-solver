/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.FactorialMomsFromMoms;
import jline.lib.butools.JFactorialMomsFromJMoms;
import jline.lib.butools.dph.MGFromMoments;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class DMRAPFromMoments {
    private DMRAPFromMoments() {}

    /**
     * Creates a discrete marked rational arrival process that has the same
     * marginal and lag-1 joint moments as given.
     */
    public static MatrixCell dmrapFromMoments(double[] moms, MatrixCell Nm) {
        MGRepresentation mgResult = MGFromMoments.mgFromMoments(moms);
        Matrix v = mgResult.alpha;
        Matrix H0 = mgResult.A;

        int N = H0.getNumRows();
        Matrix I = Matrix.eye(N);

        Matrix H0i = I.sub(H0).inv();
        Matrix Ge = new Matrix(N, N);
        Matrix G1 = new Matrix(N, N);

        Matrix H0ip = Matrix.eye(N);
        for (int i = 0; i < N; i++) {
            Matrix row = v.mult(H0ip);
            for (int j = 0; j < N; j++) {
                Ge.set(i, j, row.get(0, j));
            }
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int k = 0; k < N; k++) {
                    sum += H0ip.get(j, k);
                }
                G1.set(j, i, sum);
            }
            H0ip = H0ip.scale((double) (i + 1)).mult(H0i);
            if (i > 0) {
                H0ip = H0ip.mult(H0);
            }
        }

        Matrix Gei = Ge.inv();
        Matrix G1i = G1.inv();

        int numTypes = Nm.size();
        MatrixCell H = new MatrixCell(numTypes + 1);
        H.set(0, H0);

        for (int i = 0; i < numTypes; i++) {
            Matrix Nmi = Nm.get(i);

            Matrix row1Input = new Matrix(1, N - 1);
            for (int j = 0; j < N - 1; j++) {
                row1Input.set(0, j, Nmi.get(0, j + 1));
            }
            Matrix row1 = FactorialMomsFromMoms.factorialMomsFromMoms(row1Input);

            Matrix col1Input = new Matrix(N - 1, 1);
            for (int j = 0; j < N - 1; j++) {
                col1Input.set(j, 0, Nmi.get(j + 1, 0));
            }
            Matrix col1 = FactorialMomsFromMoms.factorialMomsFromMoms(col1Input);

            Matrix midInput = new Matrix(N - 1, N - 1);
            for (int r = 0; r < N - 1; r++) {
                for (int c = 0; c < N - 1; c++) {
                    midInput.set(r, c, Nmi.get(r + 1, c + 1));
                }
            }
            Matrix mid = JFactorialMomsFromJMoms.jFactorialMomsFromJMoms(midInput);

            Matrix NmiTransformed = new Matrix(N, N);
            NmiTransformed.set(0, 0, Nmi.get(0, 0));
            for (int j = 0; j < N - 1; j++) {
                NmiTransformed.set(0, j + 1, row1.get(0, j));
                NmiTransformed.set(j + 1, 0, col1.get(j, 0));
            }
            for (int r = 0; r < N - 1; r++) {
                for (int c = 0; c < N - 1; c++) {
                    NmiTransformed.set(r + 1, c + 1, mid.get(r, c));
                }
            }

            H.set(i + 1, I.sub(H0).mult(Gei).mult(NmiTransformed).mult(G1i));
        }

        return H;
    }

    public static MatrixCell dmrapFromMoments(double[] moms, Matrix[] Nm) {
        MatrixCell cell = new MatrixCell(Nm.length);
        for (int i = 0; i < Nm.length; i++) {
            cell.set(i, Nm[i]);
        }
        return dmrapFromMoments(moms, cell);
    }
}
