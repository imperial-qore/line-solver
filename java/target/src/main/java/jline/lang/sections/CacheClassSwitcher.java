package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.util.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class switcher section based on cache hits and misses
 */
public class CacheClassSwitcher extends StatefulClassSwitcher implements Serializable {
    private final List<JobClass> jobClasses;
    private final int items;
    private final Matrix cap;
    private final int levels;
    public Matrix hitClass;
    public Matrix missClass;
    public Matrix actualHitProb;
    public Matrix actualMissProb;
    public Map<Integer, InputJobClassesObj> inputJobClasses;

    public CacheClassSwitcher(List<JobClass> jobClasses, int items, Matrix capacity) {
        this(jobClasses, items, capacity, 1);
    }

    public CacheClassSwitcher(List<JobClass> jobClasses, int items, Matrix capacity, int levels){
        super(jobClasses, "Cache");
        this.jobClasses = jobClasses;
        this.items = items;
        this.cap = capacity;
        this.levels = levels;
        this.csFun = (input) -> this.simpleHitMiss(input.r, input.s, input.state, input.statedep);
        this.hitClass = new Matrix(1, jobClasses.size());
        this.hitClass.fill(-1);
        this.missClass = new Matrix(1, jobClasses.size());
        this.missClass.fill(-1);
        this.actualHitProb = new Matrix(0,0);
        this.actualMissProb = new Matrix(0,0); // These fields will be filled after obtaining the model solution
        this.inputJobClasses = new HashMap<>();
    }

    public double simpleHitMiss(int r, int s){
        return simpleHitMiss(r, s, null, null);
    }

    public double simpleHitMiss(int r, int s, Matrix state){
        return simpleHitMiss(r, s, null, null);
    }

    public double simpleHitMiss(int r, int s, Matrix state, Matrix statep){
        boolean hitFlag, missFlag;
        hitFlag = missFlag = false;
        for(int i = 0; i < hitClass.getNumCols(); i++){
            if(r == hitClass.get(i))
                hitFlag = true;
        }
        for(int i = 0; i < missClass.getNumCols(); i++){
            if(r == missClass.get(i))
                missFlag = true;
        }
        if(state == null || state.isEmpty()){
            if((r == s // hit and miss in the cache can depart in the same class
                || ((r < this.hitClass.getNumCols() && r < this.missClass.getNumCols()) && (s == this.hitClass.get(r) || s == this.missClass.get(r))))
                && (hitFlag || missFlag)){
                return 1;
            } else {
                return 0;
            }
        }
        return 0;
    }

    public static class InputJobClassesObj implements Serializable{
        public JobClass jobClass;
        public SchedStrategyType schedPolicy;
        public DropStrategy dropStrategy;

        public InputJobClassesObj(JobClass jobClass, SchedStrategyType schedPolicy, DropStrategy dropStrategy) {
            this.jobClass = jobClass;
            this.schedPolicy = schedPolicy;
            this.dropStrategy = dropStrategy;
        }
    }
}
