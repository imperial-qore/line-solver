package jline.lang.nodes;
import java.io.Serializable;

import jline.lang.*;
import jline.solvers.ssa.events.ArrivalEvent;
import jline.util.Matrix;

/**
 * A node that can have a state
 */
public class StatefulNode extends Node implements Serializable {
    private Integer statefulIndex;
    private Matrix state; // TODO: this seems misaligned with the existence of a State class
    private Matrix statePrior;
    
    public StatefulNode(String name) {
        super(name);
        statefulIndex = null;
        state = new Matrix(0,0,0);
        statePrior = new Matrix(0,0,0);
    }

    protected void clearState()  {
        this.state.reshape(0, 0, 0);
    }

    public int getStatefulIndex() {
        if (this.statefulIndex == null) {
            this.statefulIndex = this.model.getStatefulNodeIndex(this);
        }
        return this.statefulIndex;
    }

    public int getNumberOfServers() {
        return 1;
    }
    
    public Matrix getState(){
    	return this.state;
    }
    
    public void setState(Matrix state) {
    	this.state = state;
    	if (state.getNumRows() != statePrior.getNumRows()) { 		
    		Matrix initPrior = new Matrix(state.getNumRows(), 1, state.getNumRows());
        	initPrior.set(0, 0, 1);
        	this.setStatePrior(initPrior);
    	}
    }
    
    public void setStatePrior(Matrix prior) {
    	this.statePrior = prior;
    	try {
    		if(state.getNumRows() != statePrior.getNumRows())
    			throw new Exception("The prior probability vector must have the same rows of the station state vector");
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    }
    
    public Matrix getStatePrior(){
    	return this.statePrior;
    }

    @Override
    public ArrivalEvent getArrivalEvent(JobClass jobClass) {
        if (!this.arrivalEvents.containsKey(jobClass)) {
            this.arrivalEvents.put(jobClass, new ArrivalEvent(this, jobClass));
        }
        return this.arrivalEvents.get(jobClass);
    }
}
