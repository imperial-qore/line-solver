package jline.lang.nodes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.constant.GlobalConstants;
import jline.lang.distributions.*;
import jline.lang.processes.Replayer;
import jline.util.SerializableFunction;
import jline.util.Matrix;

import jline.lang.*;
import jline.lang.constant.DropStrategy;
import jline.solvers.ssa.events.DepartureEvent;

/**
 * A node where jobs can spend time stationing there
 */
public class Station extends StatefulNode implements Serializable {
    protected int numberOfServers;
    protected double cap; // double to allow infinite values.
    protected Map<JobClass, Double> classCap;
    protected Map<JobClass, DropStrategy> dropRule;
    protected Matrix lldScaling;
    protected SerializableFunction<Matrix, Double> lcdScaling;

    protected Map<JobClass, DepartureEvent> departureEvents;

    public Station(String name) {
        super(name);

        this.classCap = new HashMap<JobClass, Double>();
        this.departureEvents = new HashMap<JobClass, DepartureEvent>();

        this.cap = Double.POSITIVE_INFINITY;
        this.dropRule = new HashMap<JobClass, DropStrategy>();
        
        this.lldScaling = new Matrix(0,0);
        this.lcdScaling = null;
    }

    public void setNumberOfServers(int numberOfServers) {
        if (numberOfServers == -1) { // may result of a double Inf is casted as input
            numberOfServers = Integer.MAX_VALUE;
        }
        this.numberOfServers = numberOfServers;
    }

    @Override
    public int getNumberOfServers() {
        return this.numberOfServers;
    }

    public void setCap(double cap) {
        this.cap = cap;
    }

    public void setCap(int cap) {
        this.cap = cap;
    }

    public void setClassCap(JobClass jobClass, double cap) {
        this.classCap.put(jobClass, cap);
    }
    public void setClassCap(JobClass jobClass, int cap) {
        this.classCap.put(jobClass, (double)cap);
    }

    public void setChainCapacity() {

    }

    @Override
    public double getClassCap(JobClass jobClass) {
        if (classCap.containsKey(jobClass)) {
            return Math.min(this.classCap.get(jobClass), cap);
        }

        return cap;
    }

    @Override
    public double getCap() {
        return cap;
    }

    public boolean[] isServiceDefined() {
        throw new RuntimeException("Not Implemented!"); // TODO: not implemented
    }

    public boolean isServiceDefined(JobClass j_class)  {
        throw new RuntimeException("Not Implemented!"); // TODO: not implemented
    }

    public boolean[] isServiceDisabled()  {
        throw new RuntimeException("Not Implemented!"); // TODO: not implemented
    }

    public boolean isServiceDisabled(JobClass j_class)  {
        throw new RuntimeException("Not Implemented!"); // TODO: not implemented
    }

    public List<Object> getSourceRates()  {
    	int nClasses = this.model.getNumberOfClasses();
    	Map<JobClass, Map<Integer, Matrix>> map = new HashMap<JobClass, Map<Integer, Matrix>>();
    	Map<JobClass, Matrix> mu = new HashMap<JobClass, Matrix>();
    	Map<JobClass, Matrix> phi = new HashMap<JobClass, Matrix>();
    	for(int r = 0; r < nClasses; r++) {
    		Source source = (Source) this;
    		JobClass jobclass = this.model.getClassByIndex(r);
			if (!source.containsJobClass(jobclass)) {
				source.setArrival(jobclass, new Disabled());
    			Matrix nan_matrix = new Matrix(1,1,1);
    			nan_matrix.fill(Double.NaN);
    			Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
    			tmp.put(0, nan_matrix);
    			tmp.put(1, nan_matrix.clone());
    			map.put(jobclass, tmp);
    			mu.put(jobclass, nan_matrix.clone());
    			phi.put(jobclass, nan_matrix.clone());
            } else if ((source.getArrivalDistribution(jobclass) instanceof Det)) {
                Distribution distr = source.getArrivalDistribution(jobclass);
                Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
                tmp.put(0, new Matrix(-distr.getRate()));
                tmp.put(1, new Matrix(distr.getRate()));
                map.put(jobclass, tmp);
                mu.put(jobclass, new Matrix(distr.getRate()));
                phi.put(jobclass, new Matrix(1.0));
            } else if (!(source.getArrivalDistribution(jobclass) instanceof Disabled)) {
                Distribution distr = source.getArrivalDistribution(jobclass);
                map.put(jobclass, ((MarkovianDistribution)distr).getRepres());
                mu.put(jobclass, ((MarkovianDistribution)distr).getMu());
                phi.put(jobclass, ((MarkovianDistribution)distr).getPhi());
    		} else {
    			Matrix nan_matrix = new Matrix(1,1,1);
    			nan_matrix.fill(Double.NaN);
    			Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
    			tmp.put(0, nan_matrix);
    			tmp.put(1, nan_matrix.clone());
    			map.put(jobclass, tmp);
    			mu.put(jobclass, nan_matrix.clone());
    			phi.put(jobclass, nan_matrix.clone());
    		}
    		
    	}
    	return new ArrayList<Object>(Arrays.asList(map, mu, phi));
    }

    public List<Object> getServiceRates()  {
    	int nClasses = this.model.getNumberOfClasses();
    	Map<JobClass, Map<Integer, Matrix>> map = new HashMap<JobClass, Map<Integer, Matrix>>();
    	Map<JobClass, Matrix> mu = new HashMap<JobClass, Matrix>();
    	Map<JobClass, Matrix> phi = new HashMap<JobClass, Matrix>();
    	for(int i = 0; i < nClasses; i++) {
    		JobClass jobclass = this.model.getClassByIndex(i);
    		Queue queue = (Queue) this;
    		// Since Delay, Join are all sub-class of Queue, and only queue and source
            // are the sub-class of station, we could cast this to queue to call setService method
    		if(!queue.containsJobClass(jobclass)) {
    			queue.setService(jobclass, new Disabled());
    			Matrix nan_matrix = new Matrix(1,1,1);
    			nan_matrix.fill(Double.NaN);
    			Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
    			tmp.put(0, nan_matrix);
    			tmp.put(1, nan_matrix.clone());
    			map.put(jobclass, tmp);
    			mu.put(jobclass, nan_matrix.clone());
    			phi.put(jobclass, nan_matrix.clone());
    		} else if (queue.getServiceProcess(jobclass) instanceof Immediate) {
    			Distribution distr = this.server.getServiceDistribution(jobclass);
    			Matrix map_matrix_1 = new Matrix(1,1,1);
    			map_matrix_1.set(0, 0, -GlobalConstants.Immediate);
    			Matrix map_matrix_2 = new Matrix(1,1,1);
    			map_matrix_2.set(0, 0, GlobalConstants.Immediate);
    			Matrix mu_matrix = new Matrix(1,1,1);
    			mu_matrix.set(0, 0, GlobalConstants.Immediate);
    			Matrix phi_matrix = new Matrix(1,1,1);
    			phi_matrix.set(0, 0, 1);
    			Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
    			tmp.put(0, map_matrix_1);
    			tmp.put(1, map_matrix_2);
    			map.put(jobclass, tmp);
    			mu.put(jobclass, mu_matrix);
    			phi.put(jobclass, phi_matrix);
            } else if ((queue.getServiceProcess(jobclass) instanceof Det)) {
                Distribution distr = this.server.getServiceDistribution(jobclass);
                Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
                tmp.put(0, new Matrix(-distr.getRate()));
                tmp.put(1, new Matrix(distr.getRate()));
                map.put(jobclass, tmp);
                mu.put(jobclass, new Matrix(distr.getRate()));
                phi.put(jobclass, new Matrix(1.0));
            } else if ((queue.getServiceProcess(jobclass) instanceof Replayer)) {
                Distribution distr = this.server.getServiceDistribution(jobclass);
                Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
                tmp.put(0, new Matrix(-distr.getRate()));
                tmp.put(1, new Matrix(distr.getRate()));
                map.put(jobclass, tmp);
                mu.put(jobclass, new Matrix(distr.getRate()));
                phi.put(jobclass, new Matrix(1.0));
            } else if (!(queue.getServiceProcess(jobclass) instanceof Disabled)) {
                //Current JLine only support Exp, Coxian, Erlang, HyperEXP and APH.
                Distribution distr = this.server.getServiceDistribution(jobclass);
                map.put(jobclass, ((MarkovianDistribution)distr).getRepres());
                mu.put(jobclass, ((MarkovianDistribution)distr).getMu());
                phi.put(jobclass, ((MarkovianDistribution)distr).getPhi());
    		} else { // Disabled
    			Matrix nan_matrix = new Matrix(1,1,1);
    			nan_matrix.fill(Double.NaN);
    			Map<Integer, Matrix> tmp = new HashMap<Integer, Matrix>();
    			tmp.put(0, nan_matrix);
    			tmp.put(1, nan_matrix.clone());
    			map.put(jobclass, tmp);
    			mu.put(jobclass, nan_matrix.clone());
    			phi.put(jobclass, nan_matrix.clone());
    		}
    	}
    	return new ArrayList<Object>(Arrays.asList(map, mu, phi));
    }
    
    public DepartureEvent getDepartureEvent(JobClass jobClass) {
        if (!this.departureEvents.containsKey(jobClass)) {
            this.departureEvents.put(jobClass, new DepartureEvent(this, jobClass));
        }
        return this.departureEvents.get(jobClass);
    }

    @Override
    public boolean isRefstat() {
        for (JobClass jobClass : this.model.getClasses()) {
            if (jobClass instanceof ClosedClass) {
                if (((ClosedClass)jobClass).getRefstat() == this) {
                    return true;
                }
            }
        }

        return false;
    }
       
    public DropStrategy getDropRule(JobClass jobclass) {
    	return this.dropRule.getOrDefault(jobclass, null);
    }
    
    public void setDropRule(JobClass jobclass, DropStrategy drop) {
    	this.dropRule.put(jobclass, drop);
    }

    public void setLimitedLoadDependence(Matrix alpha) {
    	this.lldScaling = alpha;
    }
    
    public Matrix getLimitedLoadDependence() {
    	return this.lldScaling;
    }
    
    public void setLimitedClassDependence(SerializableFunction<Matrix, Double> gamma) {
    	this.lcdScaling = gamma;
    }
    
    public SerializableFunction<Matrix, Double> getLimitedClassDependence(){
    	return this.lcdScaling;
    }
}
