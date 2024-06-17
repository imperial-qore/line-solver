package jline.solvers.mva.analyzers;

import jline.api.CACHE;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;
import jline.lang.NodeParam;

/**
 * MVA Analyzer class for non-rentrant caches
 */
public class SolverMVACacheAnalyzer implements MVAAnalyzer{
    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        long startTime = System.currentTimeMillis();
        Matrix QN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix UN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix RN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix TN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix CN = new Matrix(1, sn.nclasses);
        Matrix AN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix WN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix XN = new Matrix(1, sn.nclasses);
        double lG = Double.NaN;
        int iter = 1;

        int source_ist = -1;
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Source){
                source_ist = (int) sn.nodeToStation.get(i);
                break;
            }
        }
        Matrix sourceRate = new Matrix(1, sn.rates.getNumCols());
        for(int i = 0; i < sn.rates.getNumCols(); i++){
            if(!Double.isNaN(sn.rates.get(source_ist, i))){
                sourceRate.set(i, sn.rates.get(source_ist, i));
            }
        }
        TN = new Matrix(sn.nodetypes.size(), sn.nclasses);
        for(int i = 0; i < sourceRate.getNumCols(); i++){
            TN.set(source_ist, i, sourceRate.get(i));
        }
        Cache cache = null;
        for(int i = 0; i < sn.nodetypes.size(); i++){
            if(sn.nodetypes.get(i) == NodeType.Cache){
                cache = (Cache) sn.nodes.get(i);
                break;
            }
        }
        NodeParam ch = sn.nodeparam.get(cache);
        Matrix m = ch.itemcap;
        int n = ch.nitems; // Number of items
        int h = m.length(); // Number of lists
        int u = sn.nclasses; // Number of users
        Matrix[] lambda = new Matrix[u];
        for(int i = 0; i < u; i++){
            lambda[i] = new Matrix(n, h + 1);
        }
        for(int v = 0; v < u; v++){
            for(int k = 0; k < n; k++){
                for(int l = 0; l < h + 1; l++){
                    if(ch.pread.getOrDefault(v, null) != null){
                        lambda[v].set(k, l, sourceRate.get(v) * ch.pread.get(v).get(k));
                    }
                }
            }
        }

        Matrix[][] R = ch.accost;
        Matrix gamma = CACHE.cache_gamma_lp(lambda, R).gamma;
        Matrix pij = null, missRate;
        switch(options.method){
            case "exact":
                pij = CACHE.cache_mva(gamma, m).pij;
                Matrix newPij = new Matrix(pij.getNumRows(), 1 + pij.getNumCols());
                for(int i = 0; i < pij.getNumRows(); i++){
                    for(int j = 0; j < newPij.getNumCols(); j++){
                        if(j < 1){
                            newPij.set(i, j, Math.abs(1 - pij.sumRows(i)));
                        } else {
                            newPij.set(i, j, pij.get(i, j - 1));
                        }
                    }
                }
                pij = newPij;
                missRate = new Matrix(1, u);
                for(int v = 0; v < u; v++){
                    missRate.set(v, Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null)).get(0));
                }
                break;
            case "ttl.lrum":
                pij = CACHE.cache_ttl_lrum(lambda, m); // without considering different graph of different items linear
                missRate = new Matrix(1, u);
                for(int v = 0; v < u; v++){
                    missRate.set(v, Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null)).get(0));
                }
                break;
            case "ttl.hlru":
                pij = CACHE.cache_ttl_hlru(lambda, m); // without considering different graph of different items linear
                missRate = new Matrix(1, u);
                for(int v = 0; v < u; v++){
                    missRate.set(v, Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null)).get(0));
                }
                break;
            case "ttl.tree":
                pij = CACHE.cache_ttl_tree(lambda, R, m); // considering different graphs of different items
                missRate = new Matrix(1, u);
                for(int v = 0; v < u; v++){
                    missRate.set(v, Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null)).get(0));
                }
                break;
            default:
                pij = CACHE.cache_prob_asy(gamma, m); // FPI method
                missRate = new Matrix(1, u);
                for(int v = 0; v < u; v++){
                    missRate.set(v, Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null)).get(0));
                }
                break;
        }
        for(int r = 0; r < sn.nclasses; r++){
            if(ch.hitclass.length() > r && ch.missclass.get(r) > -1 && ch.hitclass.get(r) > -1){
                XN.set((int) ch.missclass.get(r), XN.get((int) ch.missclass.get(r)) + missRate.get(r));
                XN.set((int) ch.hitclass.get(r), XN.get((int) ch.hitclass.get(r)) + sourceRate.get(r) - missRate.get(r));
            }
        }
        long endTime = System.currentTimeMillis();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000.0;
        res.iter = iter;
    }
}
