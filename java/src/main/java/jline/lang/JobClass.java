package jline.lang;

import java.io.Serializable;

import jline.lang.constant.JobClassType;
import jline.lang.nodes.Node;

/**
 * Superclass representing a class of jobs
 */
public class JobClass extends NetworkElement implements Serializable {
    protected JobClassType type;
    protected int priority;
    protected boolean completes;
    protected Node refstat;
    protected boolean isrefclass;
    protected int index;

    private Integer[] attribute;

    public JobClass(JobClassType type, String name) {
        super(name);
        this.priority = 0;
        this.refstat = new Node("Unallocated");
        this.type = type;
        this.completes = true;
        this.isrefclass = false;
        this.attribute = new Integer[] {null,null};
    }

    public void setReferenceStation(Node source) throws Exception {
        this.refstat = source;
    }

    public boolean isReferenceStation(Node node) {
        return name.equals(node.getName());
    }

    public void printSummary() {
        System.out.format("Job Class: %s\n",this.getName());
        System.out.println();
    }

    public double getNumberOfJobs() {
        return Double.POSITIVE_INFINITY;
    }

    public int getJobClassIdx() {
        return -1;
    }
    
    public boolean isReferenceClass() {
    	return this.isrefclass;
    }
    
    public void setReferenceClass(boolean isrefclass) {
    	this.isrefclass = isrefclass;
    }
    
    public int getPriority() {
    	return this.priority;
    }
    
    public void setPriority(int p) {
    	this.priority = p;
    }
    public JobClassType getJobClassType() { return this.type; }

    public boolean getCompletes(){ return this.completes; }
    public void  setCompletes(boolean completes){
        this.completes =completes;
    }

    public int getIndex() {
        return index;
    }

    public Integer[] getAttribute() {
        return attribute;
    }

    public void setAttribute(Integer[] attribute) {
        this.attribute = attribute;
    }

}
