package jline.lang.processes;


import jline.util.Matrix;
import jline.lang.distributions.MarkovianDistribution;
import jline.lang.distributions.CumulativeDistribution;

import java.io.Serializable;
import java.util.*;

import static jline.lib.KPCToolbox.*;

/**
 * A Markovian Arrival Process
 */
public class MAP extends MarkovianDistribution implements Serializable {
    //TODO: several methods missing
    List<Double> totalDepartureRate;
    List<Double> totalPhaseRate;
    private final int nPhases;

    public MAP(Matrix D0, Matrix D1) {
        super("MAP", 2);
        int nPhases = D0.getNumCols();
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);

        this.totalDepartureRate = new ArrayList<Double>(nPhases);
        this.totalPhaseRate = new ArrayList<Double>(nPhases);

        for (int i = 0; i < nPhases; i++) {
            double tpr = 0.0;
            double tdr = 0.0;
            for (int j = 0; j < nPhases; j++) {
                tdr += D1.get(i, j);
                if (i == j) {
                    continue;
                }
                tpr += D0.get(i, j);
            }
            this.totalPhaseRate.add(tpr);
            this.totalDepartureRate.add(tdr);
        }
        this.nPhases = nPhases;
    }

    public Matrix getD0() {
        return (Matrix) this.getParam(1).getValue();
    }

    public Matrix getD1() {
        return (Matrix) this.getParam(2).getValue();
    }

    public void normalize() {
        // TODO: incomplete, this is not normalizing so that D0+D1 rows are an infinitesimal generator
        Matrix D0 = getD0();
        Matrix D1 = getD1();
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (D0.get(i, j) < 0) {
                    D0.set(i, j, 0.0);
                }
                if (D1.get(i, j) < 0) {
                    D1.set(i, j, 0.0);
                }
            }
        }
    }

    public long getNumberOfPhases() {
        return ((Matrix) this.getParam(1).getValue()).getNumCols();
    }

    public double getMean() {
        double E1 = map_moment(getD0(), getD1(), 1);
        return E1;
    }

    @Override
    public Matrix sample(long n) {
        return this.sample(n, new Random());
    }

    @Override
    public Matrix sample(long n, Random random) {
        throw new RuntimeException("Not implemented");
    }

    public List<Double> getMoments() {
        List<Double> moments = new ArrayList<>();
        for (int i = 1; i <= 3; i++) {
            moments.add(map_moment(getD0(), getD1(), i));
        }
        return moments;
    }

    public double getVar() {
        double E1 = map_moment(getD0(), getD1(), 1);
        double E2 = map_moment(getD0(), getD1(), 2);
        return E2 - E1 * E1;
    }

    public double getSkew() {
        double E1 = map_moment(getD0(), getD1(), 1);
        double E2 = map_moment(getD0(), getD1(), 2);
        double E3 = map_moment(getD0(), getD1(), 3);
        double skew = E3 - 3 * E2 * E1 + 2 * E1 * E1 * E1;
        double scv = (E2 - E1 * E1) / E1 / E1;
        skew = skew / Math.pow(Math.sqrt(scv) * E1, 3);
        return skew;
    }

    public double getSCV() {
        double mean = this.getMean();
        return this.getVar() / mean / mean;
    }

    public double getRate() {
        return 1.0 / this.getMean();
    }

    public double getDepartureRate(int phase) {
        return this.totalDepartureRate.get(phase);
    }

    public double getTotalPhaseRate(int phase) {
        return this.totalPhaseRate.get(phase);
    }

    public int getNextPhaseAfterDeparture(int curPhase, Random random) {
        List<List<Double>> phaseRates = ((Matrix) this.getParam(2).getValue()).toDoubleList();
        List<Double> phaseTransitions = phaseRates.get(curPhase);
        double tdr = this.totalDepartureRate.get(curPhase);

        CumulativeDistribution<Integer> phaseCumulativeDistribution = new CumulativeDistribution<Integer>(random);

        for (int i = 0; i < this.nPhases; i++) {
            phaseCumulativeDistribution.addElement(i, phaseTransitions.get(i) / tdr);
        }

        return phaseCumulativeDistribution.sample(random);
    }

    public int getNextPhase(int curPhase, Random random) {
        List<List<Double>> phaseRates = ((Matrix) this.getParam(2).getValue()).toDoubleList();
        List<Double> phaseTransitions = phaseRates.get(curPhase);

        CumulativeDistribution<Integer> phaseCumulativeDistribution = new CumulativeDistribution<Integer>(random);

        double tpr = this.getTotalPhaseRate(curPhase);

        for (int i = 0; i < this.nPhases; i++) {
            if (i == curPhase) {
                continue;
            }
            phaseCumulativeDistribution.addElement(i, phaseTransitions.get(i) / tpr);
        }

        return phaseCumulativeDistribution.sample(random);
    }

    public double evalCDF(double t) {
        throw new RuntimeException("Not Implemented!");
    }

    @Override
    public Map<Integer, Matrix> getPH() {
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
        Matrix D0 = this.getD0();
        Matrix D1 = this.getD1();
        res.put(0,D0);
        res.put(1,D1);
        return res;
    }

    public double evalLST(double s) {
        throw new RuntimeException("Not Implemented!");
    }
}
