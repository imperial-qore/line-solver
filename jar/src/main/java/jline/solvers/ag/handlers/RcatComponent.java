package jline.solvers.ag.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.lang.processes.DiscreteDistribution;
import jline.util.matrix.Matrix;

/**
 * One RCAT component: the (station, class) pair and the Markovian processes its
 * QBD is assembled from.
 *
 * The phase of the QBD is the pair (arrival phase, service phase) in the
 * Kronecker order of qbd_mapmap1, so the level blocks are kron(D1a, I) up,
 * krons(D0a, D0s) locally at a busy level, kron(D0a, I) at level 0 where no
 * server runs, and kron(I, D1s) down.
 */
public final class RcatComponent {
    public int ist;
    public int r;
    public Matrix Da0;   // arrival MAP of the external streams reaching (ist,r)
    public Matrix Da1;
    public Matrix Ds0;   // service MAP of (ist,r)
    public Matrix Ds1;
    public Matrix Dsvc;  // kron(I_na, D1s): a completion, level down
    public int na;
    public int ns;
    public int mph;
    public int nlev;
    public int N;
    public double lamNeg;    // single-removal negative signal rate
    public double lamCat;    // catastrophe rate
    public final List<Double> batchRates = new ArrayList<Double>();
    public final List<DiscreteDistribution> batchDists = new ArrayList<DiscreteDistribution>();
}
