package jline.lang.distributions;

import java.util.*;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;
import static jline.lib.KPCToolbox.*;

@SuppressWarnings("unchecked")
public class PH extends MarkovianDistribution {

    List<Double> totalPhaseRate;

    public PH(List<Double> alpha, Matrix T) {
        super("PH", 3);

        int nPhases = alpha.size();

        this.setParam(1, "n", nPhases);
        this.setParam(2, "alpha", alpha);
        this.setParam(3, "T", T);
        this.totalPhaseRate = new ArrayList<Double>(nPhases);

        for (int i = 0; i < nPhases; i++) {
            double tpr = 0.0;
            tpr -= T.get(i,i);
            this.totalPhaseRate.add(tpr);
        }
    }

    public Matrix getInitProb() {
        List<Double> param1 = (List<Double>) this.getParam(2).getValue();
        Matrix alpha = new Matrix(1, param1.size(), param1.size());
        for(int i = 0; i < param1.size(); i++)
            alpha.set(0, i, param1.get(i));
        return alpha;
    }

    @Override
    public List<Double> sample(long n) {
        return this.sample(n,null);
    }
    @Override
    public List<Double> sample(long n, Random random) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double getMean() {
        return super.getMean();
    }

    @Override
    public double getRate() {
        return 1/getMean();
    }

    @Override
    public double getSCV() {
        return super.getSCV();
    }

    @Override
    public double getVar() {
        return this.getSCV()*Math.pow(this.getMean(), 2);
    }

    @Override
    public double getSkew() {
        return super.getSkew();
    }

    @Override
    public double evalCDF(double t) {
        //Since currently no function to support calculating eigen values, thus calculating expm. This method is now not implemented.
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double evalLST(double s) {
        return super.evalLST(s);
    }

    @Override
    public Map<Integer, Matrix> getPH() {
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
        Matrix T = getSubgenerator();

        Matrix ones = new Matrix(T.numCols,1,T.numCols);
        Matrix Te = new Matrix(0,0,0);
        Matrix Tepie = new Matrix(0,0,0);
        ones.fill(1.0);
        T.mult(ones, Te);
        Te.mult(this.getInitProb(), Tepie);
        //Tepie.removeZeros(0);
        Tepie.changeSign();

        res = map_normalize(T,Tepie);
        return res;
    }

    @Override
    public long getNumberOfPhases() {
        return (long) this.getParam(1).getValue();
    }

    public Matrix getSubgenerator() {
        return (Matrix) this.getParam(3).getValue();
    }

    public double getTotalPhaseRate(int i) {
        return totalPhaseRate.get(i);
    }
}
