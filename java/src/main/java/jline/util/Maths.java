package jline.util;

import java.util.*;

import jline.api.PFQN;
import jline.lang.state.ThreadLocalRandom;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.special.Gamma;

/**
 * Mathematical functions and utilities.
 */
public class Maths {

  // these to be used only in testing to replicate matlab random numbers
  // keep false in standard operation as Mersenne Twister slower than Math.random()
    private static Boolean matlabRandom = false;


    /**
   * Softmin function.
   *
   * @param x     first term to compare
   * @param y     second term to compare
   * @param alpha softmin smoothing parameter
   * @return Softmin function value
   */
  public static double softmin(double x, double y, double alpha) {
    return -((-x) * Math.exp(-alpha * x) - y * Math.exp(-alpha * y))
        / (Math.exp(-alpha * x) + Math.exp(-alpha * y));
  }

  /**
   * Returns the max of two numbers. If one of the two is NaN, it returns the other.
   *
   * @param x - the first number to be compared
   * @param y - the second number to be compared
   * @return - the max between the two
   */
  public static double max(double x, double y) {
    if (!Double.isNaN(x) && !Double.isNaN(y))
      return Math.max(x, y);
    else if (Double.isNaN(x))
      return y;
    return x;
  }

  /**
   * Returns the min of two numbers. If one of the two is NaN, it returns the other.
   *
   * @param x - the first number to be compared
   * @param y - the second number to be compared
   * @return - the min between the two
   */
  public static double min(double x, double y) {
    if (!Double.isNaN(x) && !Double.isNaN(y))
      return Math.min(x, y);
    else if (Double.isNaN(x))
      return y;
    return x;
  }

  // Error function
  public static double erf(double x) {
    // Constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    // Save the sign of x
    int sign = (x < 0) ? -1 : 1;
    x = Math.abs(x);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t;

    return sign * (1 - y * Math.exp(-x * x));
  }

  public static void setRandomNumbersMatlab(boolean matlab_style) {
      matlabRandom = matlab_style;
  }

  public static boolean randomMatlabStyle() {
    return matlabRandom;
  }

  public static void setMatlabRandomSeed(final long seed) {
      ThreadLocalRandom.setSeed((int)seed);
  }

  // returns a random double in interal (0,0)
  public static double random() {
    if (matlabRandom) {
      return ThreadLocalRandom.random();
    } else {
      return Math.random();
    }
  }

  public static double factln(double n) {
    return Math.log(Gamma.gamma(1 + n));
  }

  // Gamma function via Lanczos approximation formula
  public static double gammaFunction(double x) {
    double[] p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
        1.5056327351493116e-7};
    int g = 7;
    if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gammaFunction(1 - x));
    x -= 1;
    double a = p[0];
    double t = x + g + 0.5;
    for (int i = 1; i < p.length; i++) {
      a += p[i] / (x + i);
    }
    return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
  }

  // Returns a circulant matrix of order c
  public static Matrix circul(int c) {
    if (c == 1) {
      Matrix C = new Matrix(1, 1);
      C.set(0, 0, 1);
      return C;
    }

    Matrix R = new Matrix(c, c);
    R.set(0, c - 1, 1);
    for (int i = 1; i < c; i++) {
      R.set(i, i - 1, 1);
    }
    R = R.transpose();
    return R;
  }

  // Returns a circulant matrix of order c
  public static Matrix circul(Matrix c) {
    int n = c.length();
    Matrix R = new Matrix(n, n);
    R.set(0, n - 1, 1);
    for (int i = 1; i < n; i++) {
      R.set(i, i - 1, 1);
    }
    Matrix C = new Matrix(n, n);
    for (int t = 0; t < n; t++) {
      if (t > 1) {
        R = new Matrix(R.mult(R, null));
      }
      Matrix tmpC = new Matrix(0, 0);
      R.scale(c.get(0, t), tmpC);
      C = C.add(1, tmpC);
    }
    return C;
  }

  /**
   * Implementation of MATLAB "hist": puts elements of v into k bins
   *
   * @param v      - vector of elements
   * @param minVal - lowest number in v
   * @param maxVal - highest number in v
   * @return - vector containing number of elements in each bin
   */
  public static Matrix hist(Matrix v, int minVal, int maxVal) {
    Matrix result = new Matrix(1, maxVal - minVal + 1);
    for (int j = 0; j < v.getNumCols(); j++) {
      int value = (int) v.get(j);
      if (value >= minVal && value <= maxVal) {
        int index = value - minVal;
        result.set(0, index, result.get(0, index) + 1);
      }
    }
    return result;
  }

  public static Matrix perms(Matrix vec) {
    // copy to avoid issues with needing vec in outer code as swap modifies inputted vector
    Matrix vecCopy = vec.clone();
    // use length to accept row or column vector
    Matrix indexes = new Matrix(1, vecCopy.length());
    indexes.zero();
    List<Matrix> lPermutations = new ArrayList<>();
    lPermutations.add(vecCopy.clone());

    int ind = 0;
    while (ind < vecCopy.length()) {
      if (indexes.get(ind) < ind) {
        swap(vecCopy, ind % 2 == 0 ? 0 : (int) indexes.get(ind), ind);
        lPermutations.add(vecCopy.clone());
        indexes.set(ind, indexes.get(ind) + 1);
        ind = 0;
      } else {
        indexes.set(ind, 0);
        ind++;
      }
    }
    lPermutations.sort((a, b) -> {
      for (int i = 0; i < a.length(); i++) {
        if (a.get(i) != b.get(i)) {
          return Double.compare(b.get(i), a.get(i));
        }
      }
      return 0;
    });


    // Construct permuations matrix from perm vectors
    Matrix permutations = new Matrix(lPermutations.size(), vec.length());
    for (int i = 0; i < lPermutations.size(); i++) {
      Matrix permutation = lPermutations.get(i);
      for (int j = 0; j < permutations.getNumCols(); j++) {
        permutations.set(i, j, permutation.get(j));
      }
    }
    return permutations;
  }

  public static Matrix uniquePerms(Matrix vec) {

    // Vector is empty
    if (vec.isEmpty()) {
      return new Matrix(0, 0);
    }

    // Vector is not a single row
    if (vec.getNumRows() != 1) {
      throw new RuntimeException("Matrix passed to uniquePerms has more than one row. Unsupported.");
    }

    // Number of elements in the vector
    int n = vec.length();

    // Number of unique elements in the vector
    int nu = 0;
    LinkedList<Double> uniqueElements = new LinkedList<>();
    for (int i = 0; i < vec.length(); i++) {
      double element = vec.get(0, i);
      boolean unique = true;
      for (Double uniqueElement : uniqueElements) {
        if (element == uniqueElement) {
          unique = false;
          break;
        }
      }
      if (unique) {
        uniqueElements.add(element);
      }
    }
    nu = uniqueElements.size();

    // Only one unique element
    if (nu == 1) {
      return vec;
    }

    // Every element is unique
    if (n == nu) {
      return perms(vec);
    }
    Matrix[] output = new Matrix[nu];
    for (int i = 0; i < nu; i++) {
      Matrix v = vec.clone();

      int ind = -1;
      double target = uniqueElements.get(i);
      for (int j = 0; j < v.length(); j++) {
        if (target == v.get(j)) {
          ind = j;
          break;
        }
      }
      Matrix newV = new Matrix(1, v.getNumCols() - 1);
      for (int j = 0, destCol = 0; j < v.getNumCols(); j++) {
        if (j != ind) {
          newV.set(0, destCol++, v.get(j));
        }
      }
      v = newV.clone();
      Matrix temp = uniquePerms(v);
      Matrix repeated = new Matrix(temp.getNumRows(), 1);
      repeated.fill(uniqueElements.get(i));
      output[i] = repeated.concatCols(temp);
    }

    // TODO: change this to avoid concat rows and fill in a new matrix
    Matrix result = output[0];
    for (int i = 1; i < output.length; i++) {
      result = Matrix.concatRows(result, output[i], null).clone();
    }
    return result;
  }


  /**
   * Computes the combinations of the elements in v taken k at a time
   *
   * @param v - vector of elements
   * @param k - how many elements to pick at a time
   * @return - the combinations of the elements in v taken k at a time
   */
  public static Matrix nchoosek(Matrix v, int k) {
    int n = v.length();
    if (k < 0 || k > n) {
      return null;
    }
    int a = 1, b = 1, c = 1; // a == n!, b == k!, c == (n-k)!
    for (int i = 2; i <= n; i++) {
      a *= i;
      if (i <= k) {
        b *= i;
      }
      if (i <= n - k) {
        c *= i;
      }
    }
    Matrix res = new Matrix(a / (b * c), k);
    int row = 0;
    int[] indexes = new int[k];
    for (int i = 0; i < indexes.length; i++) {
      indexes[i] = i;
    }
    while (true) {
      for (int i = 0; i < k; i++) {
        res.set(row, i, v.get(indexes[i]));
      }
      row++;
      int last = k - 1;
      while (last >= 0 && indexes[last] == n - k + last) {
        last--;
      }
      if (last == -1) {
        break;
      }
      indexes[last]++;
      for (int i = last + 1; i < k; i++) {
        indexes[i] = indexes[i - 1] + 1;
      }
    }
    return res;
  }

  public static double multinomialln(Matrix n) {
    return factln(n.elementSum()) - Matrix.factln(n).elementSum();
  }

  public static Matrix multiChoose(double n, double k) {

    Matrix v = new Matrix(1, (int) n);
    v.zero();

    if (n == 1) {
      v = new Matrix(1, 1);
      v.set(0, 0, k);
    } else if (k != 0) {


      List<Matrix> tmpSSRows = new ArrayList<>();
      for (int i = 0; i <= k; i++) {
        Matrix w = multiChoose(n - 1, k - i);

        Matrix tmpSSRow = new Matrix(w.getNumRows(), w.getNumCols() + 1);
        for (int j = 0; j < w.getNumRows(); j++) {
          tmpSSRow.set(j, 0, i);
          for (int l = 1; l < w.getNumCols() + 1; l++) {
            tmpSSRow.set(j, l, w.get(j, l - 1));
          }
        }
        tmpSSRows.add(tmpSSRow);
      }
      int totalRows = 0;
      for (Matrix m : tmpSSRows) {
        totalRows += m.getNumRows();
      }
      v = new Matrix(totalRows, (int) n);
      int rowForV = 0;
      for (int i = 0; i < tmpSSRows.size(); i++) {
        int rowForTmpSSRows = 0;
        int initialRowForV = rowForV;
        for (int j = rowForV; j < initialRowForV + tmpSSRows.get(i).getNumRows(); j++) {
          for (int l = 0; l < tmpSSRows.get(i).getNumCols(); l++) {
            v.set(j, l, tmpSSRows.get(i).get(rowForTmpSSRows, l));
          }
          rowForV++;
          rowForTmpSSRows++;
        }
      }
    }

    return v;
  }


  public static Matrix multiChooseCon(Matrix n, double S) {
    // n row vector
    // n[i] is number of elements at ith position we have to choose from
    Matrix v = new Matrix(0, 0);
    List<Integer> indices = new ArrayList<>();
    for (int i = 0; i < n.getNumElements(); i++) {
      if (n.get(i) != 0) {
        indices.add(i);
      }
    }
    if (S == 1) {
      for (int i : indices) {
        v.expandMatrix(v.getNumRows() + 1, n.getNumCols(),
            (v.getNumRows() + 1) * n.getNumCols());
        int lastRow = v.getNumRows() - 1;
        for (int col = 0; col < n.getNumCols(); col++) {
          v.set(lastRow, col, 0);
        }
        v.set(lastRow, i, 1);
      }
      return v;
    } else {
      for (int i : indices) {
        Matrix n1 = n.clone();
        n1.set(0, i, n1.get(i) - 1);
        Matrix T = multiChooseCon(n1, S - 1);
        Matrix y = new Matrix(T.getNumRows(), n.getNumCols());
        y.zero();
        for (int row = 0; row < y.getNumRows(); row++) {
          y.set(row, i, 1);
        }
        if (v.isEmpty()) {
          v = y.add(1, T);
        } else {
          v = Matrix.concatRows(v, y.add(1, T), null);
        }
      }
      return v;
    }
  }

  /**
   * Cumulative sum of an array, where the value at each index of the result is the
   * sum of all the previous values of the input.
   *
   * @param ar Array we want the cumulative sum of.
   * @return New array that is the cumulative sum of the input.
   */
  public static int[] cumSum(int[] ar) {
    int[] result = new int[ar.length];
    result[0] = ar[0];
    for (int i = 1; i < ar.length; i++) {
      result[i] = result[i - 1] + ar[i];
    }
    return result;
  }

  /**
   * Cumulative sum of an array, where the value at each index of the result is the
   * sum of all the previous values of the input.
   *
   * @param ar Array we want the cumulative sum of.
   * @return New array that is the cumulative sum of the input.
   */
  public static double[] cumSum(double[] ar) {
    double[] result = new double[ar.length];
    result[0] = ar[0];
    for (int i = 1; i < ar.length; i++) {
      result[i] = result[i - 1] + ar[i];
    }
    return result;
  }

  public static double[] logspace(int start, int stop, int n) {
    double base = 10;
    double logMax = Math.log10(stop);
    double logMin = Math.log10(start);
    double delta = (logMax - logMin) / (n - 1);
    double exp = logMin;
    double[] vector = new double[n];
    for (int i = 0; i < n; i++) {
      vector[i] = Math.pow(base, exp);
      exp += delta;
    }
    return vector;
  }

  /**
   * Generates an integer list of n logarithmically spaced values.
   *
   * @param start First element of array.
   * @param stop  Last element of array.
   * @param n     Number of values.
   * @return Array of logarithmically spaced values.
   */
  public static int[] logspacei(int start, int stop, int n) {
    // I know this sucks but im making it java 7 compatible on the last possible day, was using lambdas
    ArrayList<Integer> vector = new ArrayList<>(n);
    int i = 0;
    for (double log : logspace(start, stop, n)) {
      vector.add((int) Math.round(log));
    }
    vector = new ArrayList<Integer>(new LinkedHashSet<Integer>(vector));
    int[] result = new int[n];
    int j = 0;
    for (Integer integer : vector) {
      result[j] = integer;
      j++;
    }
    return Arrays.copyOf(result, j);
  }

  /**
   * Transposes a 2D array of values.
   *
   * @param data 2D array to be transposed.
   * @return New 2D array of transposed data.
   */
  public static double[][] transpose(double[][] data) {
    int n = data.length;
    int m = data[0].length;
    double[][] transpose = new double[m][n];
    for (int r = 0; r < m; r++) {
      for (int c = 0; c < n; c++) {
        transpose[r][c] = data[c][r];
      }
    }
    return transpose;
  }

  /**
   * Helper method that adds an integer to every element of an
   * integer array.
   *
   * @param ar Array to be added to.
   * @param i  Integer to add.
   */
  public static void elementAdd(int[] ar, int i) {
    for (int j = 0; j < ar.length; j++) {
      ar[j] += i;
    }
  }

  /**
   * Given an integer x input, find the next integer y greater than x
   * such that y is a power of 2.
   *
   * @param x Input to find next power
   * @return Closest integer greater than input that is a power of two.
   */
  protected static int nextPowerOfTwo(int x) {
    int highestOneBit = Integer.highestOneBit(x);
    return (x == highestOneBit) ? x : highestOneBit << 1;
  }

  /**
   * helper method that swaps elements in a matrix vector around
   *
   * @param vec - the row/column we are swapping elements within
   * @param a   - index into vec
   * @param b   - index into vec
   */
  private static void swap(Matrix vec, int a, int b) {
    double tmp = vec.get(a);
    vec.set(a, vec.get(b));
    vec.set(b, tmp);
  }

  public static class nrl_h implements jline.util.SerializableFunction<Matrix, Complex> {
    public Matrix L;
    public Matrix tsubtb;
    double Nt;
    Matrix alpha;

    public nrl_h(Matrix L, Matrix tsubtb, double Nt, Matrix alpha) {
      this.L = L;
      this.tsubtb = tsubtb;
      this.Nt = Nt;
      this.alpha = alpha;
    }

    @Override
    public Complex apply(Matrix x) {
      int M = L.getNumRows();
      ComplexMatrix c = new ComplexMatrix(tsubtb.getNumRows(), tsubtb.getNumCols());
      for (int i = 0; i < c.getNumElements(); i++) {
        c.set(i, new Complex(0, 2 * Math.PI * tsubtb.get(i)).exp());
      }
      ComplexMatrix L_complex = new ComplexMatrix(L.elementMult(c.real.repmat(M, 1), null), L.elementMult(c.im.repmat(M, 1), null));
      Matrix exp_x = x.clone();
      for (int i = 0; i < x.getNumElements(); i++) {
        exp_x.set(i, Math.exp(x.get(i)) / Math.pow(1 + Math.exp(x.get(i)), 2));
      }
      return PFQN.pfqn_gld_complex(L_complex.sumRows(), new Matrix(Nt), alpha, null).G.multiply(exp_x.prodVector());
    }
  }

  public static class nrp_h implements jline.util.SerializableFunction<Matrix, Complex> {
    public Matrix L;
    public Matrix tsubtb;
    double Nt;
    Matrix alpha;

    public nrp_h(Matrix L, Matrix tsubtb, double Nt, Matrix alpha) {
      this.L = L;
      this.tsubtb = tsubtb;
      this.Nt = Nt;
      this.alpha = alpha;
    }

    @Override
    public Complex apply(Matrix x) {
      int M = L.getNumRows();
      ComplexMatrix c = new ComplexMatrix(tsubtb.getNumRows(), tsubtb.getNumCols());
      for (int i = 0; i < c.getNumElements(); i++) {
        c.set(i, new Complex(0, 2 * Math.PI * tsubtb.get(i)).exp());
      }
      ComplexMatrix L_complex = new ComplexMatrix(L.elementMult(c.real.repmat(M, 1), null), L.elementMult(c.im.repmat(M, 1), null));
      NormalDistribution Z = new NormalDistribution();
      Matrix normpdf_x = x.clone();
      for (int i = 0; i < x.getNumElements(); i++) {
        normpdf_x.set(i, Z.density(x.get(i)));
      }
      return PFQN.pfqn_gld_complex(L_complex.sumRows(), new Matrix(Nt), alpha, null).G.multiply(normpdf_x.prodVector());
    }
  }

  public static class infradius_h implements jline.util.SerializableFunction<Matrix, ComplexMatrix> {
    public Matrix L;
    public Matrix N;
    public Matrix alpha;

    public infradius_h(Matrix L, Matrix N, Matrix alpha) {
      this.L = L;
      this.N = N;
      this.alpha = alpha;
    }

    @Override
    public ComplexMatrix apply(Matrix x) {
      int M = L.getNumRows();
      double Nt = N.elementSum();
      Matrix beta = new Matrix(N.getNumRows(), N.getNumCols());
      N.divide(Nt, beta, true);
      Matrix t = new Matrix(x.getNumRows(), x.getNumCols());
      Matrix tb_mat = new Matrix(x.getNumRows(), 1);
      for (int i = 0; i < x.getNumRows(); i++) {
        double tb_row = 0.0;
        for (int j = 0; j < x.getNumCols(); j++) {
          t.set(i, j, Math.exp(x.get(i, j)) / (1 + Math.exp(x.get(i, j))));
          tb_row += beta.get(i, j) * Math.exp(x.get(i, j)) / (1 + Math.exp(x.get(i, j)));
        }
        tb_mat.set(i, tb_row);
      }
      double tb = tb_mat.elementSum();
      Matrix tsubtb = t.clone();
      for (int i = 0; i < t.getNumElements(); i++) {
        tsubtb.set(i, tsubtb.get(i) - tb);
      }
      nrl_h h = new nrl_h(L, tsubtb, Nt, alpha);
      ComplexMatrix y = new ComplexMatrix(new Matrix(x.getNumRows(), 1), new Matrix(x.getNumRows(), 1));
      for (int i = 0; i < x.getNumRows(); i++) {
        y.set(i, h.apply(Matrix.extractRows(x, i, x.getNumRows(), null)));
      }
      return y;
    }
  }

  public static class infradius_hnorm implements jline.util.SerializableFunction<Matrix, ComplexMatrix> {
    public Matrix L;
    public Matrix N;
    public Matrix alpha;

    public infradius_hnorm(Matrix L, Matrix N, Matrix alpha) {
      this.L = L;
      this.N = N;
      this.alpha = alpha;
    }

    @Override
    public ComplexMatrix apply(Matrix x) {
      int M = L.getNumRows();
      double MU = 0;
      double SIGMA = 1;
      double Nt = N.elementSum();
      Matrix beta = new Matrix(N.getNumRows(), N.getNumCols());
      N.divide(Nt, beta, true);
      NormalDistribution Z = new NormalDistribution(MU, SIGMA);
      Matrix t = x.clone();
      for (int i = 0; i < t.getNumRows(); i++) {
        for (int j = 0; j < t.getNumCols(); j++) {
          t.set(i, j, Z.cumulativeProbability(x.get(i, j)));
        }
      }
      double tb = 0;
      for (int i = 0; i < t.getNumRows(); i++) {
        for (int j = 0; j < t.getNumCols(); j++) {
          tb += beta.get(i, j) * t.get(i, j);
        }
      }
      Matrix tsubtb = t.clone();
      for (int i = 0; i < t.getNumElements(); i++) {
        tsubtb.set(i, tsubtb.get(i) - tb);
      }
      ComplexMatrix y = new ComplexMatrix(new Matrix(x.getNumRows(), 1), new Matrix(x.getNumRows(), 1));
      nrp_h h = new nrp_h(L, tsubtb, Nt, alpha);
      for (int i = 0; i < x.getNumRows(); i++) {
        y.set(i, h.apply(Matrix.extractRows(x, i, x.getNumRows(), null)));
      }
      return new ComplexMatrix(y.real);
    }
  }

  public static ComplexMatrix num_grad_h_complex(Matrix x0, double h, SerializableFunction<Matrix, ComplexMatrix> hfun) {
    ComplexMatrix df = new ComplexMatrix(x0.clone(), x0.clone());
    df.zero();
    for (int i = 0; i < x0.getNumCols(); i++) {
      Matrix x1 = x0.clone();
      Matrix x2 = x0.clone();
      x1.set(i, x0.get(i) - h);
      x2.set(i, x0.get(i) + h);
      Complex y1;
      Complex y2;
      y1 = hfun.apply(x1).get(0).log();
      y2 = hfun.apply(x2).get(0).log();
      df.set(i, y2.subtract(y1).divide(2 * h));
    }
    return df;
  }

  public static Matrix num_grad_h(Matrix x0, double h, SerializableFunction<Matrix, Matrix> hfun) {
    Matrix df = new Matrix(x0.clone());
    df.zero();
    for (int i = 0; i < x0.getNumCols(); i++) {
      Matrix x1 = x0.clone();
      Matrix x2 = x0.clone();
      x1.set(i, x0.get(i) - h);
      x2.set(i, x0.get(i) + h);
      double y1;
      double y2;
      y1 = Math.log(hfun.apply(x1).get(0));
      y2 = Math.log(hfun.apply(x2).get(0));
      df.set(i, y2 - y1 / (2 * h));
    }
    return df;
  }

    public static List<Integer> sub2ind(Matrix dims, Matrix row, Matrix col) {
        List<Integer> result = new ArrayList<>();
        int numRows = (int) dims.get(0);
        int numCols = (int) dims.get(1);
        for (int i = 0; i < row.getNumElements(); i++) {

            int index = (int) (col.get(i) * numRows + row.get(i));
            result.add(index);
        }
        return result;
    }

  public static ComplexMatrix num_hess_h_complex(Matrix x0, double h, SerializableFunction<Matrix, ComplexMatrix> hfun) {
    ComplexMatrix H = new ComplexMatrix(0, x0.getNumElements());
    H.zero();
    for (int i = 0; i < x0.getNumElements(); i++) {
      Matrix x1 = x0.clone();
      x1.set(i, x1.get(i) - h);
      ComplexMatrix df1 = num_grad_h_complex(x1, h, hfun);

      Matrix x2 = x0.clone();
      x2.set(i, x2.get(i) + h);
      ComplexMatrix df2 = num_grad_h_complex(x2, h, hfun);
      ComplexMatrix d2f = new ComplexMatrix(1, x0.getNumElements());
      for (int j = 0; j < d2f.getNumElements(); j++) {
        d2f.set(j, df2.get(j).subtract(df1.get(j)).divide(2 * h));
      }
      H = ComplexMatrix.concatRows(H, d2f, null);
    }
    return H;
  }

  public static Matrix num_hess_h(Matrix x0, double h, SerializableFunction<Matrix, Matrix> hfun) {
    Matrix H = new Matrix(0, x0.getNumElements());
    H.zero();
    for (int i = 0; i < x0.getNumElements(); i++) {
      Matrix x1 = x0.clone();
      x1.set(i, x1.get(i) - h);
      Matrix df1 = num_grad_h(x1, h, hfun);

      Matrix x2 = x0.clone();
      x2.set(i, x2.get(i) + h);
      Matrix df2 = num_grad_h(x2, h, hfun);
      Matrix d2f = new Matrix(1, x0.getNumElements());
      for (int j = 0; j < d2f.getNumElements(); j++) {
        d2f.set(j, (df2.get(j) - df1.get(j)) / (2 * h));
      }
      H = Matrix.concatRows(H, d2f, null);
    }
    return H;
  }

  public static laplaceApproxComplexReturn laplaceapprox_h_complex(Matrix x0, SerializableFunction<Matrix, ComplexMatrix> h) {
    int d = x0.getNumCols();
    double tol = 1e-5;
    Complex detnH = new Complex(-1);
    ComplexMatrix H = new ComplexMatrix(1, 1);
    H.set(0, 0, -1);
    while (detnH.getReal() < 0 && tol <= 1e-3) {
      H = num_hess_h_complex(x0, tol, h);
      ComplexMatrix nH = H.clone();
      nH.scale(-1);
      detnH = nH.det();
      tol *= 10;
    }
    if (detnH.getReal() < 0) {
      System.out.println("Warning: laplaceapprox_h_complex: det(-H)<0.");
    }
    ComplexMatrix infrad = h.apply(x0);
    Complex I = infrad.get(0).multiply(new Complex(Math.pow(2 * Math.PI, d)).divide(detnH).sqrt());
    Complex logI = infrad.get(0).log().add((new Complex((double) d / 2 * Math.log(2 * Math.PI)))).subtract(detnH.log());
    return new laplaceApproxComplexReturn(H, I, logI);
  }

  public static laplaceApproxReturn laplaceapprox_h(Matrix x0, SerializableFunction<Matrix, Matrix> h) {
    int d = x0.getNumCols();
    double tol = 1e-5;
    double detnH = -1;
    Matrix H = new Matrix(1, 1);
    H.set(0, 0, -1);
    while (detnH < 0 && tol <= 1e-3) {
      H = num_hess_h(x0, tol, h);
      Matrix nH = H.clone();
      nH.scale(-1);
      detnH = nH.det();
      tol *= 10;
    }
    if (detnH < 0) {
      System.out.println("Warning: laplaceapprox_h: det(-H)<0.");
    }
    Matrix infrad = h.apply(x0);
    double I = infrad.get(0) * (Math.pow(2 * Math.PI, d) / Math.sqrt(detnH));
    double logI = Math.log(infrad.get(0)) + ((double) d / 2 * Math.log(2 * Math.PI)) - Math.log(detnH);
    return new laplaceApproxReturn(H, I, logI);
  }


  public static class laplaceApproxComplexReturn {
    public ComplexMatrix H;
    public Complex I;
    public Complex logI;

    public laplaceApproxComplexReturn(ComplexMatrix H, Complex I, Complex logI) {
      this.H = H;
      this.I = I;
      this.logI = logI;
    }
  }

  public static class laplaceApproxReturn {
    public Matrix H;
    public double I;
    public double logI;

    public laplaceApproxReturn(Matrix H, double I, double logI) {
      this.H = H;
      this.I = I;
      this.logI = logI;
    }
  }

  public static void main(String[] args) {
     Maths.circul((int) 3).print();
 }
}
