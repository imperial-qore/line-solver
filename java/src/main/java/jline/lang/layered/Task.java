package jline.lang.layered;


import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;

import java.util.ArrayList;
import java.util.List;

/**
 * A LayeredNetwork entity that can host services specified in the form of Entry objects and that runs on a Host
 */
public class Task extends LayeredNetworkElement{
    protected Processor parent;
    protected int multiplicity;
    protected int replication;
    protected SchedStrategy scheduling;//Enum
    protected Distribution thinkTime;
    protected double thinkTimeMean;
    protected double thinkTimeSCV;
    protected List<Entry> entries;
    protected List<Activity> activities;
    protected List<ActivityPrecedence> precedences;
    private List<Entry> replyEntry;

    public Task(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, Distribution thinkTime) {
        super(name);
        this.parent = null;
        this.setReplication(1);
        this.model = model;
        this.multiplicity = multiplicity;
        this.scheduling = scheduling;
        this.setThinkTime(thinkTime);
        this.entries = new ArrayList<>();
        this.activities =new ArrayList<>();
        this.precedences = new ArrayList<>();
        this.replyEntry = new ArrayList<>();
        // link within model
        model.tasks.put(model.tasks.size(),this);
        model.nodes.put(model.tasks.size(),this);
        if(scheduling == SchedStrategy.REF) {
            model.reftasks.put(model.reftasks.size(),this);
        }
    }

    public Task(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        this(model, name, multiplicity, scheduling, new Immediate());
    }

    public Task(LayeredNetwork model, String name, int multiplicity) {
        this(model, name, multiplicity, SchedStrategy.INF, new Immediate());
    }

    public Task(LayeredNetwork model, String name) {
        this(model, name, 1, SchedStrategy.INF, new Immediate());
    }

    public Task setReplication(int replication) {
        this.replication = replication;
        return this;
    }

    public Task on(Processor parent){
        if (this.parent != null) {
            this.parent.removeTask(this);
        }
        this.parent = parent;
        this.parent.addTask(this);
        return this;
    }

    public Task setAsReferenceTask(){
        this.scheduling = SchedStrategy.REF;
        return this;
    }

    public Task setThinkTime(Distribution thinkTime) {
        this.thinkTime = thinkTime;
        this.thinkTimeMean = thinkTime.getMean();
        this.thinkTimeSCV = thinkTime.getSCV();
        return this;
    }

    public Task setThinkTime(double thinkTime){
        if(thinkTime<=GlobalConstants.Zero){
            this.thinkTime = new Immediate();
            this.thinkTimeMean = GlobalConstants.Zero;
            this.thinkTimeSCV = GlobalConstants.Zero;
        }else {
            this.thinkTime = new Exp(1/thinkTime);
            this.thinkTimeMean = thinkTime;
            this.thinkTimeSCV = 1.0;
        }
        return this;
    }

    public Task addEntry(Entry newEntry){
        this.entries.add(newEntry);
        return this;
    }

    public Task addActivity(Activity newActivity) {
        newActivity.setParent(this);
        this.activities.add(newActivity);
        return this;
    }

    public Task setActivity(Activity newActivity, int index){
        this.activities.set(index,newActivity);
        return this;
    }

    public Task removeActivity(int index){
        this.activities.remove(index);
        return this;
    }

    public Task addPrecedence(ActivityPrecedence newPrec){
        this.precedences.add(newPrec);
        return this;
    }

    public Task addPrecedence(List<ActivityPrecedence> newPrec){
        this.precedences.addAll(newPrec);
        return this;
    }

    public Task setReplyEntry(List<Entry> replyEntry) {
        this.replyEntry = replyEntry;
        return this;
    }

    public double getMeanHostDemand(String entryName){
        double meanHostDemand = -1;
        for(int j=0;j<this.entries.size();j++){
            if(this.entries.get(j).getName().equals(entryName)){
                break;
            }
        }
        return meanHostDemand;
    }

}
