package jline.api;

import jline.util.Maths;
import jline.util.Matrix;
import org.apache.commons.lang3.NotImplementedException;

import java.util.ArrayList;
import java.util.Collections;

/**
 * APIs for stochastic models of caches
 */
public class CACHE {

    /**
     * Computes access factors for the cache
     * @param lambda - request arrival rates from users to items of individual lists
     * @param R - reachability graph of a list
     * @return - access factors, number of users, items, and lists
     */
    public static cacheGammaLpReturn cache_gamma_lp(Matrix[] lambda, Matrix[][] R){
        int u = lambda.length; // Number of users
        int n = lambda[0].getNumRows(); // Number of items
        int h = lambda[0].getNumCols() - 1; // Number of lists
        Matrix gamma = new Matrix(n, h);

        for(int i = 0; i < n; i++){ // for all items
            for(int j = 0; j < h; j++){ // for all levels
                // Compute gamma(i,j)
                Matrix Rvi = new Matrix(R[0][i].getNumRows(), R[0][i].getNumCols());
                for(int v = 0; v < u; v++){
                    Rvi = Rvi.add(1, R[v][i]);
                }
                ArrayList<Integer> Pij = new ArrayList<>();
                Pij.add(1 + j);
                ArrayList<Integer> pr_j = par(Rvi, 1 + j);
                while(pr_j.size() > 0){
                    Pij.add(0, pr_j.get(0));
                    pr_j = par(Rvi, pr_j.get(0));
                }
                if(Pij.size() == 0){
                    gamma.set(i, j, 0);
                } else {
                    gamma.set(i, j, 1);
                    for(int li = 1; li < Pij.size(); li++){ // For all levels up to the current one
                        double y = 0;
                        int l_1 = Pij.get(li - 1);
                        int l = Pij.get(li);
                        for(int v = 0; v < u; v++){ // For all streams
                            for(int t = 0; t <= l_1; t++){
                                y += lambda[v].get(i, t) * R[v][i].get(t, l);
                            }
                        }
                        gamma.set(i, j, gamma.get(i, j) * y);
                    }
                }
            }
        }
        return new cacheGammaLpReturn(gamma, u, n, h);
    }

    /**
     * Finds the parent of j according to the access probabilities in R
     * @param R - the access probabilities
     * @param j - the list for which the parent will be found
     * @return - the parent of j
     */
    public static ArrayList<Integer> par(Matrix R, int j){
        ArrayList<Integer> parent = new ArrayList<>();
        for(int i = 0; i <= j - 1; i++){
            if(R.get(i, j) != 0){
                parent.add(i);
            }
        }
        if(parent.size() > 1){
            throw new RuntimeException("A cache has a list with more than one parent, but the structure must be a tree.");
        }
        return parent;
    }

    /**
     * Exact recursive solution of the caching model
     * @param gamma - access factors for item i to access list j
     * @param m - cache capacity vector
     * @return - cache state probabilities at equilibrium
     *
     * Reference: Giuliano Casale, Nicolas Gast. Performance analysis methods for list-based caches
     * with non-uniform access. IEEE/ACM Transactions on Networking, 2020, pp.1-18, for details.
     */
    public static cacheMVAReturn cache_mva(Matrix gamma, Matrix m){
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        Matrix SS = new Matrix(0,0);
        for(int l = 0; l < h; l++){
            Matrix arg = new Matrix((int) m.get(l) + 1, 1);
            for(int i = 1; i <= m.get(l) + 1; i++){
                arg.set(i - 1, i);
            }
            SS = ssg_decorate(SS, arg);
        }
        for(int i = 0; i < SS.getNumRows(); i++){
            for(int j = 0; j < SS.getNumCols(); j++){
                SS.set(i, j, SS.get(i, j) - 1);
            }
        }
        Matrix pi = new Matrix(SS.getNumRows(), n);
        Matrix[] pij = new Matrix[SS.getNumRows()];
        for(int i = 0; i < SS.getNumRows(); i++){
            pij[i] = new Matrix(n, h);
        }
        Matrix x = new Matrix(1, h);
        int E = 1;
        for(int s = 0; s < SS.getNumRows(); s++){
            Matrix mcur = Matrix.extractRows(SS, s, s+1, null);
            for(int l = 0; l < h; l++){
                Matrix mcur_l = Matrix.oner(mcur, new ArrayList<>(Collections.singletonList(l)));
                int s_l = Matrix.matchrow(SS, mcur_l);
                if(s_l >= 0){
                    Matrix one_pi = Matrix.extractRows(pi, s_l, s_l + 1, null);
                    for(int i = 0; i < one_pi.getNumCols(); i++){
                        one_pi.set(i, 1 - one_pi.get(i));
                    }
                    x.set(l, mcur.get(l)/Matrix.extractColumn(gamma, l, null).transpose().mult(one_pi.transpose()).get(0));
                    Matrix elMult = Matrix.extractColumn(gamma, l, null).transpose().elementMult(one_pi, null);
                    for(int i = 0; i < n; i++){
                        pij[s].set(i, l, elMult.get(i) * x.get(l));
                        pi.set(s, i, pi.get(s, i) + pij[s].get(i, l));
                    }

                }
            }
        }
        int s = Matrix.matchrow(SS, m);
        pi = Matrix.extractRows(pi, s, s + 1, null).transpose();
        Matrix newPij = new Matrix(n, h);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < h; j++){
                newPij.set(i, j, pij[s].get(i, j));
            }
        }
        Matrix pi0 = new Matrix(pi.getNumRows(), pi.getNumCols());
        for(int i = 0; i < pi0.getNumRows(); i++){
            for(int j = 0; j < pi0.getNumCols(); j++){
                pi0.set(i, j, 1 - pi.get(i, j));
            }
        }
        Matrix u = new Matrix(n, h);
        for(int l = 0; l < h; l++){
            for(int k = 0; k < n; k++){
                u.set(k, l, x.get(l) * gamma.get(k, l));
            }
        }
        return new cacheMVAReturn(pi, pi0, newPij, x, u, E);
    }

    /**
     * Appends states as columns to an existing state space
     * @param SS - existing state space
     * @param SS2 - new state columns
     * @return - new state space
     */
    public static Matrix ssg_decorate(Matrix SS, Matrix SS2){
        if(SS.isEmpty())
            return SS2;
        if(SS2.isEmpty())
            return SS;
        int n1 = SS.getNumRows();
        int m1 = SS.getNumCols();

        int n2 = SS2.getNumRows();
        int m2 = SS2.getNumCols();
        SS = SS.repmat(n2, 1);

        Matrix newSS = new Matrix(SS.getNumRows(), SS.getNumCols() + (m1 + m2) - (m1 + 1) + 1);
        for(int i = 0; i < SS.getNumRows(); i++){
            for(int j = 0; j < SS.getNumCols(); j++){
                newSS.set(i, j, SS.get(i, j));
            }
        }
        SS = newSS;
        int curStatesLeft = 0, curStatesRight = n1;
        for(int s = 0; s < n2; s++){
            Matrix rep = Matrix.extractRows(SS2, s, s+1, null).repmat(curStatesRight - curStatesLeft, 1);
            int repi = 0;
            for(int i = curStatesLeft; i < curStatesRight; i++){
                int repj = 0;
                for(int j = m1; j < m1 + m2; j++){
                    SS.set(i, j, rep.get(repi, repj));
                    repj++;
                }
                repi++;
            }
            curStatesLeft += n1;
            curStatesRight += n1;
        }
        return SS;
    }

    // TODO
    public static Matrix cache_ttl_lrum(Matrix[] gamma, Matrix m){
        throw new NotImplementedException("cache_ttl_lrum not implemented!");
    }

    // TODO
    public static Matrix cache_ttl_hlru(Matrix[] gamma, Matrix m){
        throw new NotImplementedException("cache_ttl_hlru not implemented!");
    }

    // TODO
    public static Matrix cache_ttl_tree(Matrix[] lambda, Matrix[][] R, Matrix m){
        throw new NotImplementedException("cache_ttl_tree not implemented!");
    }

    /**
     * Estimate asymptotic values of the cache state probabilities at equilibrium
     * @param gamma - cache access factors
     * @param m -  cache capacity vector
     * @return - cache state probabilities at equilibrium
     */
    public static Matrix cache_prob_asy(Matrix gamma, Matrix m){
        // FPI method
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        Matrix xi = cache_xi_fp(gamma, m, null).xi;
        Matrix prob = new Matrix(n, h + 1);
        for(int i = 0; i < n; i++){
            Matrix mul = Matrix.extractRows(gamma, i, i + 1, null).mult(xi.columnMajorOrder(), null);
            prob.set(i, 0, 1/(1 + mul.get(0)));
            for(int j = 1; j < h + 1; j++){
                prob.set(i, j, mul.get(0) / (1 + mul.get(0)));
            }
        }
        return prob;
    }

    /**
     * Estimate cache xi terms using a fixed-point algorithm. The xi terms are duals of the throughput of a queueing
     * network but in the caching setting.
     * @param gamma - cache access factors
     * @param m -  cache capacity vector
     * @param xi - optional initial value of the fixed point iteration
     * @return - cache xi terms
     */
    public static cacheXiFpReturn cache_xi_fp(Matrix gamma, Matrix m, Matrix xi){
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        double tol = 1.0e-14;
        Matrix pi0 = new Matrix(1, n);
        pi0.fill(1.0/(h+1));
        Matrix pij = new Matrix(n, h);
        if(xi == null || xi.isEmpty()){
            xi = new Matrix(1, h);
            for(int l = 0; l < h; l++){
                double mean = Matrix.extractColumn(gamma, l, null).elementSum() / gamma.getNumRows();
                xi.set(l, m.get(l) / mean / (n + m.elementSum() - 1));
            }
        }
        int it = 1;
        for(; it < 1.0e4; it++){
            Matrix pi0_1 = new Matrix(pi0);
            Matrix mul = pi0_1.mult(gamma, null);
            xi = new Matrix(m.getNumRows(), m.getNumCols());
            for(int i = 0; i < xi.getNumRows(); i++){
                for(int j = 0; j < xi.getNumCols(); j++){
                    xi.set(i, j, m.get(i, j) / mul.get(i, j));
                }
            }
            Matrix intermediate = gamma.elementMult(xi.repmat(n, 1), null);
            mul = gamma.mult(xi.repmat(n, 1).transpose());
            for(int i = 0; i < pij.getNumRows(); i++){
                for(int j = 0; j < pij.getNumCols(); j++){
                    pij.set(i, j, Math.abs(intermediate.get(i, j)) / Math.abs(1 + mul.get(i, j)));
                }
            }
            for(int i = 0; i < n; i++){
                pi0.set(0, i, Maths.max(tol, 1 - pij.sumRows(i)));
            }
            double DELTA = 0;
            for(int i = 0; i < n; i++){
                DELTA += Math.abs(1 - pi0.get(i) / pi0_1.get(i));
            }
            if(DELTA < tol){
                break;
            }
        }
        for(int i = 0; i < xi.getNumRows(); i++){
            for(int j = 0; j < xi.getNumCols(); j++){
                if(xi.get(i, j) < 0){
                    xi.set(i, j, tol);
                }
            }
        }
        return new cacheXiFpReturn(xi, pi0, pij, it);
    }

    public static class cacheXiFpReturn {
        public Matrix xi;
        public Matrix pi0;
        public Matrix pij;
        public int it;

        public cacheXiFpReturn(Matrix xi, Matrix pi0, Matrix pij, int it) {
            this.xi = xi;
            this.pi0 = pi0;
            this.pij = pij;
            this.it = it;
        }
    }

    public static class cacheGammaLpReturn {
        public Matrix gamma;
        public int u;
        public int n;
        public int h;

        public cacheGammaLpReturn(Matrix gamma, int u, int n, int h) {
            this.gamma = gamma;
            this.u = u;
            this.n = n;
            this.h = h;
        }
    }

    public static class cacheMVAReturn {
        public Matrix pi;
        public Matrix pi0;
        public Matrix pij;
        public Matrix x;
        public Matrix u;
        public int E;

        public cacheMVAReturn(Matrix pi, Matrix pi0, Matrix pij, Matrix x, Matrix u, int E) {
            this.pi = pi;
            this.pi0 = pi0;
            this.pij = pij;
            this.x = x;
            this.u = u;
            this.E = E;
        }
    }
}
