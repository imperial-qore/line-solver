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
        if(scheduling == SchedStrategy.REF) {
            model.reftasks.put(model.reftasks.size(),this);
        }
    }

    public Task(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        super(name);
        this.parent = null;
        this.setReplication(1);
        this.model = model;
        this.multiplicity = multiplicity;
        this.scheduling = scheduling;
        setThinkTime(GlobalConstants.Zero);
        model.tasks.put(model.tasks.size(),this);
        this.entries = new ArrayList<>();
        this.activities =new ArrayList<>();
        this.precedences = new ArrayList<>();
        this.replyEntry = new ArrayList<>();
        if(scheduling == SchedStrategy.REF) model.reftasks.put(model.tasks.size(),this);
    }

    public Task(LayeredNetwork model, String name, int multiplicity) {
        super(name);
        this.parent = null;
        this.setReplication(1);
        this.model = model;
        this.multiplicity = multiplicity;
        this.scheduling = SchedStrategy.INF;
        setThinkTime(GlobalConstants.Zero);
        this.entries = new ArrayList<>();
        this.activities =new ArrayList<>();
        this.precedences = new ArrayList<>();
        this.replyEntry = new ArrayList<>();
        model.tasks.put(model.tasks.size(),this);
    }

    public Task(LayeredNetwork model, String name) {
        super(name);
        this.parent = null;
        this.setReplication(1);
        this.model = model;
        this.multiplicity = 1;
        this.scheduling = SchedStrategy.INF;
        setThinkTime(GlobalConstants.Zero);
        this.entries = new ArrayList<>();
        this.activities =new ArrayList<>();
        this.precedences = new ArrayList<>();
        this.replyEntry = new ArrayList<>();
        model.tasks.put(model.tasks.size(),this);
    }

    public void setReplication(int replication) {
        this.replication = replication;
    }

    public void on(Processor parent){
        if (this.parent != null) {
            this.parent.removeTask(this);
        }
        this.parent = parent;
        this.parent.addTask(this);
    }

    public void setAsReferenceTask(){
        this.scheduling = SchedStrategy.REF;
    }

    public void setThinkTime(Distribution thinkTime) {
        this.thinkTime = thinkTime;
        this.thinkTimeMean = thinkTime.getMean();
        this.thinkTimeSCV = thinkTime.getSCV();
    }

    public void setThinkTime(double thinkTime){
        if(thinkTime<=GlobalConstants.Zero){
            this.thinkTime = new Immediate();
            this.thinkTimeMean = GlobalConstants.Zero;
            this.thinkTimeSCV = GlobalConstants.Zero;
        }else {
            this.thinkTime = new Exp(1/thinkTime);
            this.thinkTimeMean = thinkTime;
            this.thinkTimeSCV = 1.0;
        }
    }

    public void addEntry(Entry newEntry){
        this.entries.add(newEntry);
    }

    public void addActivity(Activity newActivity) {

        newActivity.setParent(this);
        this.activities.add(newActivity);
    }

    public void setActivity(Activity newActivity, int index){
        this.activities.set(index,newActivity);
    }

    public void removeActivity(int index){
        this.activities.remove(index);
    }

    public void addPrecedence(ActivityPrecedence newPrec){
        this.precedences.add(newPrec);
    }

    public void addPrecedence(List<ActivityPrecedence> newPrec){
        this.precedences.addAll(newPrec);
    }

    public void setReplyEntry(List<Entry> replyEntry) {
        this.replyEntry = replyEntry;
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
