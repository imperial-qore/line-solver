package jline.lang.layered;


import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;
import jline.util.Matrix;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

/**
 * An element modeling an individual service activity
 */
public class Activity extends LayeredNetworkElement {
    protected Distribution hostDemand;//TODO [Distribution]
    protected double hostDemandMean;
    protected double hostDemandSCV;
    protected Task parent;//
    protected String parentName;
    protected String boundToEntry;
    protected String callOrder;
    protected Map<Integer,String> syncCallDests = new HashMap<>();
    protected Matrix syncCallMeans = new Matrix(1,1,0);
    protected Map<Integer,String> asyncCallDests = new HashMap<>();
    protected Matrix asyncCallMeans = new Matrix(1,1,0);
    protected Matrix scheduling = new Matrix(0,0,0);

    public Activity(LayeredNetwork model, String name, Distribution hostDemand, String boundToEntry, String callOrder) {
        super(name);

        this.setHostDemand(hostDemand);
        this.boundToEntry = boundToEntry;
        this.setCallOrder(callOrder);
        model.activities.put(model.activities.size(),this);
        this.model = model;
    }

    public Activity(LayeredNetwork model, String name, Distribution hostDemand, String boundToEntry) {
        super(name);

        this.setHostDemand(hostDemand);
        this.boundToEntry = boundToEntry;
        this.setCallOrder("STOCHASTIC");
        model.activities.put(model.activities.size(),this);
        this.model = model;
    }

    public Activity(LayeredNetwork model, String name, Distribution hostDemand) {
        super(name);
        this.setHostDemand(hostDemand);
        this.boundToEntry = "";
        this.setCallOrder("STOCHASTIC");
        model.activities.put(model.activities.size(),this);
        this.model = model;
    }

    public Activity(LayeredNetwork model, String name) {
        super(name);
        this.setHostDemand(1e-8);
        this.boundToEntry = "";
        this.setCallOrder("STOCHASTIC");
        model.activities.put(model.activities.size(),this);
        this.model = model;
    }


    public void setParent(Task parent){
        this.parentName = parent.getName();//TODO

        this.parent = parent;

    }

    public void on(Task parent){
        parent.addActivity(this);
        this.parent = parent;
    }



    public void setHostDemand(double hostDemand) {
        if (hostDemand <= GlobalConstants.Zero){
            this.hostDemand = new Immediate();
            this.hostDemandMean = 1e-8;
            this.hostDemandSCV = 1e-8;
        }else{
            this.hostDemand =  new Exp(1/hostDemand);
            this.hostDemandMean = hostDemand;
            this.hostDemandSCV = 1.0;
        }
    }

    public void setHostDemand(Distribution hostDemand) {
        this.hostDemand = hostDemand;
        this.hostDemandMean = hostDemand.getMean();
        this.hostDemandSCV = hostDemand.getSCV();
    }

    public void repliesTo(Entry entry) throws Exception {
        if(this.parent!=null){
            if (Objects.requireNonNull(this.parent.scheduling) == SchedStrategy.REF) {
                throw new Exception("Activities in reference tasks cannot reply.");
            } else {
                entry.replyActivity.put(entry.replyActivity.size(), this.getName());
            }
        }else{
            entry.replyActivity.put( entry.replyActivity.size(), this.getName());
        }
    }

    public  void boundTo(Entry entry){
            this.boundToEntry = entry.getName();
    }

    public void boundTo(String entry) {
        this.boundToEntry = entry;
    }

    public void setCallOrder(String callOrder) {
        if (callOrder.equals("STOCHASTIC") || callOrder.equals("DETERMINISTIC")){
            this.callOrder = callOrder.toUpperCase(Locale.ROOT);
        }else{
            this.callOrder = "STOCHASTIC";
        }
    }

    public void synchCall(Entry synchCallDest, double synchCallMean){//TODO:condition?
        this.syncCallDests.put(syncCallDests.size(),synchCallDest.getName());
        this.syncCallMeans.growMaxColumns(syncCallMeans.getNumCols()+1,true);
        this.syncCallMeans.set(0,syncCallMeans.length()-1,synchCallMean);
    }

    public void synchCall(String synchCallDest, double synchCallMean){//TODO:condition?
        this.syncCallDests.put(syncCallDests.size(),synchCallDest);
        this.syncCallMeans.growMaxColumns(syncCallMeans.getNumCols()+1,true);
        this.syncCallMeans.set(0,syncCallMeans.length()-1,synchCallMean);
    }

    public void synchCall(Entry synchCallDest){//TODO:condition?
        this.syncCallDests.put(syncCallDests.size(),synchCallDest.getName());
        this.syncCallMeans.growMaxColumns(syncCallMeans.getNumCols()+1,true);
        this.syncCallMeans.growMaxLength(syncCallMeans.getNonZeroLength()+1,true);
        this.syncCallMeans.set(0,syncCallMeans.getNumCols()-1,1);
    }

    public void synchCall(String synchCallDest){//TODO:condition?
        this.syncCallDests.put(syncCallDests.size(),synchCallDest);
        this.syncCallMeans.growMaxColumns(syncCallMeans.getNumCols()+1,true);
        this.syncCallMeans.set(0,syncCallMeans.length()-1,1);
    }

    public void asynchCall(Entry asynchCallDest, double asynchCallMean){//TODO: type checking?
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest.getName());
        this.asyncCallMeans.growMaxColumns(asyncCallMeans.getNumCols()+1,true);
        this.asyncCallMeans.set(0,asyncCallMeans.length()-1,asynchCallMean);
    }

    public void asynchCall(String asynchCallDest, double asynchCallMean){//TODO: type checking?
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest);
        this.asyncCallMeans.growMaxColumns(asyncCallMeans.getNumCols()+1,true);
        this.asyncCallMeans.set(0,asyncCallMeans.length()-1,asynchCallMean);
    }

    public void asynchCall(Entry asynchCallDest){
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest.getName());
        this.asyncCallMeans.growMaxColumns(asyncCallMeans.getNumCols()+1,true);
        this.asyncCallMeans.set(0,asyncCallMeans.length()-1,1);
    }

    public void asynchCall(String asynchCallDest){
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest);
        this.asyncCallMeans.growMaxColumns(asyncCallMeans.getNumCols()+1,true);
        this.asyncCallMeans.set(0,asyncCallMeans.length()-1,1);
    }
}
