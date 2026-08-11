/**
 * @file Spectral decomposition result for the KPC-Toolbox.
 *
 * @since LINE 3.0
 */
package jline.lib.kpctoolbox.basic;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of spectral decomposition.
 */
public final class SpectralDecomposition {
    private final double[] spectrum;
    private final List<Matrix> projectors;
    private final Matrix eigenvectors;
    private final Matrix eigenvalueMatrix;

    public SpectralDecomposition(double[] spectrum, List<Matrix> projectors,
                                 Matrix eigenvectors, Matrix eigenvalueMatrix) {
        this.spectrum = spectrum;
        this.projectors = projectors;
        this.eigenvectors = eigenvectors;
        this.eigenvalueMatrix = eigenvalueMatrix;
    }

    public double[] getSpectrum() { return spectrum; }
    public List<Matrix> getProjectors() { return projectors; }
    public Matrix getEigenvectors() { return eigenvectors; }
    public Matrix getEigenvalueMatrix() { return eigenvalueMatrix; }

    public double[] component1() { return spectrum; }
    public List<Matrix> component2() { return projectors; }
    public Matrix component3() { return eigenvectors; }
    public Matrix component4() { return eigenvalueMatrix; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SpectralDecomposition)) return false;
        SpectralDecomposition that = (SpectralDecomposition) o;
        return Arrays.equals(spectrum, that.spectrum)
                && Objects.equals(projectors, that.projectors);
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(spectrum);
        result = 31 * result + Objects.hashCode(projectors);
        return result;
    }
}
