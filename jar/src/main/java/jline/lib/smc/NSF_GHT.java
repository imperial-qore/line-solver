/**
 * @file NSF_GHT - Non-Skip-Free Gail-Hantler-Taylor Algorithm
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class NSF_GHT {
    private NSF_GHT() {}

    /**
     * Gail-Hantler-Taylor algorithm for Non-Skip-Free Markov chains.
     */
    public static Matrix nsfGht(Matrix A, int N, NSFGHTOptions options) {
        int m = A.getNumRows();
        int K = A.getNumCols() / m - 1;

        if (N < 1) throw new IllegalArgumentException("N must be at least 1");
        if (K < 0) throw new IllegalArgumentException("A must have at least one block");

        int numit = 0;
        double check = 1.0;

        int initialCols = Math.min(m * N, A.getNumCols());
        Matrix G;
        if (initialCols > 0) {
            G = A.extractCols(0, initialCols);
        } else {
            G = Matrix.zeros(m, m * N);
        }
        if (G.getNumCols() < m * N) {
            Matrix paddedG = Matrix.zeros(m, m * N);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < G.getNumCols(); j++) {
                    paddedG.set(i, j, G.get(i, j));
                }
            }
            G = paddedG;
        }

        while (check > 1e-14 && numit < options.getMaxNumIt()) {
            numit++;
            Matrix Gold = G.copy();
            Matrix temp = G.copy();

            Matrix AN;
            if ((N + 1) * m <= A.getNumCols()) {
                AN = A.extractCols(N * m, (N + 1) * m);
            } else {
                AN = Matrix.zeros(m, m);
            }
            G = A.extractCols(0, m * N).add(AN.mult(G));

            for (int j = N + 1; j <= K; j++) {
                Matrix newTemp = Matrix.zeros(m, m * N);
                for (int i = 0; i < m; i++) {
                    for (int col = 0; col < (N - 1) * m; col++) {
                        newTemp.set(i, col + m, temp.get(i, col));
                    }
                }
                Matrix tempRight = temp.extractCols((N - 1) * m, N * m);
                Matrix product = tempRight.mult(Gold);
                for (int i = 0; i < m; i++) {
                    for (int col = 0; col < m * N; col++) {
                        newTemp.set(i, col, newTemp.get(i, col) + product.get(i, col));
                    }
                }
                temp = newTemp;

                if ((j + 1) * m <= A.getNumCols()) {
                    Matrix Aj = A.extractCols(j * m, (j + 1) * m);
                    G = G.add(Aj.mult(temp));
                }
            }

            check = Matrix.infNorm(G.sub(Gold));
        }

        if (numit == options.getMaxNumIt() && check > 1e-14) {
            InputOutput.line_warning("NSF_GHT", "Maximum Number of Iterations %d reached", numit);
        }

        if (!options.isFirstBlockRow() && N > 1) {
            Matrix fullG = Matrix.zeros(m * N, m * N);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m * N; j++) {
                    fullG.set(i, j, G.get(i, j));
                }
            }

            Matrix G1 = G.copy();
            for (int j = 2; j <= N; j++) {
                Matrix prevRow = new Matrix(m, m * N);
                for (int i = 0; i < m; i++) {
                    for (int c = 0; c < m * N; c++) {
                        prevRow.set(i, c, fullG.get((j - 2) * m + i, c));
                    }
                }
                Matrix shiftedPart = Matrix.zeros(m, m * N);
                for (int i = 0; i < m; i++) {
                    for (int c = 0; c < (N - 1) * m; c++) {
                        shiftedPart.set(i, c + m, prevRow.get(i, c));
                    }
                }
                Matrix rightPart = prevRow.extractCols((N - 1) * m, N * m);
                Matrix multPart = rightPart.mult(G1);
                for (int i = 0; i < m; i++) {
                    for (int c = 0; c < m * N; c++) {
                        fullG.set((j - 1) * m + i, c, shiftedPart.get(i, c) + multPart.get(i, c));
                    }
                }
            }
            G = fullG;
        }

        return G;
    }

    public static Matrix nsfGht(Matrix A, int N) {
        return nsfGht(A, N, new NSFGHTOptions());
    }

    /**
     * NSF_pi - Stationary vector of a Non-Skip-Free Markov chain.
     */
    public static Matrix nsfPi(Matrix B, Matrix A, Matrix G, NSFPiOptions options) {
        int m = A.getNumRows();
        int K = A.getNumCols() / m - 1;
        int N = G.getNumCols() / m;

        Matrix fullG;
        if (options.isFirstBlockRow() && G.getNumRows() == m) {
            fullG = expandGBlockRows(G, N, m);
        } else {
            fullG = G;
        }

        int extraBlocksC = N - 1 - ((K - 1) % N);
        int CWidth = m * (N - 1) + A.getNumCols() + m * extraBlocksC;
        Matrix C = Matrix.zeros(N * m, CWidth);

        for (int row = 0; row < m; row++) {
            for (int col = 0; col < (K + 1) * m; col++) {
                C.set(row, col, A.get(row, col));
            }
        }
        for (int i = 1; i < N; i++) {
            for (int row = 0; row < m; row++) {
                for (int col = 0; col < (K + 1) * m; col++) {
                    C.set(i * m + row, i * m + col, A.get(row, col));
                }
            }
        }

        if (B != null) {
            int extraBlocksB = N - 1 - (K % N);
            Matrix paddedB;
            if (extraBlocksB > 0) {
                Matrix newB = Matrix.zeros(N * m, B.getNumCols() + m * extraBlocksB);
                for (int i = 0; i < B.getNumRows(); i++) {
                    for (int j = 0; j < B.getNumCols(); j++) {
                        newB.set(i, j, B.get(i, j));
                    }
                }
                paddedB = newB;
            } else {
                paddedB = B;
            }
            return MG1_pi.mg1_pi(paddedB, C, new MG1PiOptions(options.getMaxNumComp() / N));
        } else {
            Matrix firstRowSumsG = new Matrix(m, N * m);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    firstRowSumsG.set(i, j, fullG.get(i, j));
                }
            }
            for (int block = 1; block < N; block++) {
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < m; j++) {
                        firstRowSumsG.set(i, block * m + j,
                                firstRowSumsG.get(i, (block - 1) * m + j) + fullG.get(i, block * m + j));
                    }
                }
            }

            Matrix lastBlock = firstRowSumsG.extractCols((N - 1) * m, N * m);
            Matrix ghat = Stat.stat(lastBlock);

            Matrix G1Row = fullG.extractRows(0, m);
            Matrix pi0 = ghat.mult(G1Row);

            Matrix g = ghat.mult(firstRowSumsG);

            Matrix beta = Matrix.zeros(m * N, 1);
            Matrix CSum = Matrix.zeros(m * N, m * N);
            int numCBlocks = CWidth / (m * N);

            for (int i = numCBlocks - 1; i >= 1; i--) {
                for (int r = 0; r < m * N; r++) {
                    for (int c = 0; c < m * N; c++) {
                        CSum.set(r, c, CSum.get(r, c) + C.get(r, i * m * N + c));
                    }
                }
                Matrix ones = Matrix.ones(m * N, 1);
                Matrix CSumOnes = CSum.mult(ones);
                for (int r = 0; r < m * N; r++) {
                    beta.set(r, 0, beta.get(r, 0) + CSumOnes.get(r, 0));
                }
            }

            for (int r = 0; r < m * N; r++) {
                for (int c = 0; c < m * N; c++) {
                    CSum.set(r, c, CSum.get(r, c) + C.get(r, c));
                }
            }

            Matrix eye = Matrix.eye(m * N);
            Matrix ones = Matrix.ones(m * N, 1);
            Matrix onesMinusBeta = ones.sub(beta);
            Matrix outerProduct = onesMinusBeta.mult(g);
            Matrix matrixToInvert = eye.sub(CSum).add(outerProduct);
            Matrix invMatrix = matrixToInvert.inv();
            Matrix temp = invMatrix.mult(ones);

            Matrix eyeZeros = new Matrix(m, N * m);
            for (int i = 0; i < m; i++) {
                eyeZeros.set(i, i, 1.0);
            }
            Matrix onesM = Matrix.ones(m, 1);
            Matrix finalMatrix = eyeZeros.sub(G1Row).add(onesM.mult(g));
            temp = finalMatrix.mult(temp);

            double normFactor = ghat.mult(temp).get(0, 0);
            pi0 = pi0.scale(1.0 / normFactor);

            Matrix Bmatrix = new Matrix(N * m, A.getNumCols());
            for (int i = 0; i < N; i++) {
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < A.getNumCols(); c++) {
                        Bmatrix.set(i * m + r, c, A.get(r, c));
                    }
                }
            }

            int extraBlocksB = N - 1 - (K % N);
            Matrix paddedB;
            if (extraBlocksB > 0) {
                Matrix newB = Matrix.zeros(N * m, Bmatrix.getNumCols() + m * extraBlocksB);
                for (int i = 0; i < Bmatrix.getNumRows(); i++) {
                    for (int j = 0; j < Bmatrix.getNumCols(); j++) {
                        newB.set(i, j, Bmatrix.get(i, j));
                    }
                }
                paddedB = newB;
            } else {
                paddedB = Bmatrix;
            }

            return MG1_pi.mg1_pi(paddedB, C,
                    new MG1PiOptions(options.getMaxNumComp() / N, options.isVerbose()));
        }
    }

    public static Matrix nsfPi(Matrix B, Matrix A, Matrix G) {
        return nsfPi(B, A, G, new NSFPiOptions());
    }

    private static Matrix expandGBlockRows(Matrix G, int N, int m) {
        if (G.getNumRows() >= N * m) {
            return G;
        }
        Matrix fullG = Matrix.zeros(m * N, m * N);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m * N; j++) {
                fullG.set(i, j, G.get(i, j));
            }
        }
        for (int j = 2; j <= N; j++) {
            Matrix prevRow = new Matrix(m, m * N);
            for (int i = 0; i < m; i++) {
                for (int c = 0; c < m * N; c++) {
                    prevRow.set(i, c, fullG.get((j - 2) * m + i, c));
                }
            }
            Matrix shiftedPart = Matrix.zeros(m, m * N);
            for (int i = 0; i < m; i++) {
                for (int c = 0; c < (N - 1) * m; c++) {
                    shiftedPart.set(i, c + m, prevRow.get(i, c));
                }
            }
            Matrix rightPart = prevRow.extractCols((N - 1) * m, N * m);
            Matrix G1 = fullG.extractRows(0, m);
            Matrix multPart = rightPart.mult(G1);
            for (int i = 0; i < m; i++) {
                for (int c = 0; c < m * N; c++) {
                    fullG.set((j - 1) * m + i, c, shiftedPart.get(i, c) + multPart.get(i, c));
                }
            }
        }
        return fullG;
    }
}
