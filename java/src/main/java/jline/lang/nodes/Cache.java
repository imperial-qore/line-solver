package jline.lang.nodes;

import jline.util.Maths;
import jline.lang.ItemSet;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.*;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Zipf;
import jline.lang.sections.Buffer;
import jline.lang.sections.CacheClassSwitcher;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.Section;
import jline.util.Matrix;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * A class switch node based on cache hits or misses
 */
public class Cache extends StatefulNode implements Serializable {
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;
    private final ItemSet items;
    private final int nLevels;
    private final double cap;
    private final Matrix itemLevelCap;
    public Matrix[][] accessProb;
    private final ReplacementStrategy replacementPolicy;
    private final Map<PopularityKey, Distribution> popularity;
    private int popularityRows;
    private int popularityColumns;
    private final Matrix[] graph;
    private final CacheClassSwitcher cacheServer;

    public Cache(Network model, String name, int nitems, int itemLevelCap, ReplacementStrategy replPolicy) {
        this(model, name, nitems, new Matrix(itemLevelCap), replPolicy, null);
    }

    public Cache(Network model, String name, int nitems, Matrix itemLevelCap, ReplacementStrategy replPolicy) {
        this(model, name, nitems, itemLevelCap, replPolicy, null);
    }

    public Cache(Network model, String name, int nitems, Matrix itemLevelCap, ReplacementStrategy replPolicy, Matrix[] graph) {
        super(name);
        List<JobClass> classes = model.getClasses();
        this.input = new Buffer(classes);
        this.output = new Dispatcher(classes);
        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.items = new ItemSet(model, name + "_Items", nitems, this);
        this.nLevels = itemLevelCap.getNonZeroLength();
        this.cap = Double.POSITIVE_INFINITY; // job capacity
        this.accessProb = null;
        this.itemLevelCap = itemLevelCap; // item capacity
        if(this.itemLevelCap.elementSum() > nitems){
            throw new RuntimeException("The number of items is smaller than the capacity of " + name);
        }
        this.replacementPolicy = replPolicy;
        this.cacheServer = new CacheClassSwitcher(classes, this.nLevels, itemLevelCap);
        this.server = this.cacheServer;
        this.popularity = new HashMap<>();
        this.popularityRows = 0;
        this.popularityColumns = 0;
        this.setModel(model);
        this.model.addNode(this);
        this.graph = graph;
    }

    /**
     * Reset the internal data structures when the network model is reset
     */
    public void reset(){
        this.cacheServer.actualHitProb = new Matrix(0,0);
        this.cacheServer.actualMissProb = new Matrix(0,0);
    }

    public void setScheduling(int jobClass, SchedStrategy strategy){}

    public void setResultHitProb(Matrix actualHitProb){
        this.cacheServer.actualHitProb = actualHitProb;
    }

    public void setResultMissProb(Matrix actualMissProb){
        this.cacheServer.actualMissProb = actualMissProb;
    }

    public Matrix getHitRatio(){
        return this.cacheServer.actualHitProb;
    }

    public Matrix getMissRatio(){
        return this.cacheServer.actualMissProb;
    }

    public void setHitClass(JobClass jobinclass, JobClass joboutclass){
        int jobinindex = jobinclass.getIndex() - 1;
        int joboutindex = joboutclass.getIndex() - 1;
        if(this.model.getClasses().size() > this.cacheServer.hitClass.getNumCols()){
            Matrix newHitClass = new Matrix(1, this.model.getClasses().size());
            newHitClass.fill(-1);
            for(int i = 0; i < this.cacheServer.hitClass.getNumCols(); i++){
                newHitClass.set(i, this.cacheServer.hitClass.get(i));
            }
            this.cacheServer.hitClass = newHitClass;
        }
        this.cacheServer.hitClass.set(jobinindex, joboutindex);
    }

    public void setMissClass(JobClass jobinclass, JobClass joboutclass){
        int jobinindex = jobinclass.getIndex() - 1;
        int joboutindex = joboutclass.getIndex() - 1;
        if(this.model.getClasses().size() > this.cacheServer.missClass.getNumCols()){
            Matrix newMissClass = new Matrix(1, this.model.getClasses().size());
            newMissClass.fill(-1);
            for(int i = 0; i < this.cacheServer.missClass.getNumCols(); i++){
                newMissClass.set(i, this.cacheServer.missClass.get(i));
            }
            this.cacheServer.missClass = newMissClass;
        }
        this.cacheServer.missClass.set(jobinindex, joboutindex);
    }

    public void setRead(JobClass jobClass, Distribution distribution){
        ItemSet itemclass = this.items;
        if(distribution.isDiscrete()){
            this.cacheServer.inputJobClasses.put(jobClass.getIndex(), new CacheClassSwitcher.InputJobClassesObj(jobClass, this.schedPolicy, DropStrategy.WaitingQueue));
            Distribution distribution_copy = null;
            try{
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);
                out.writeObject(distribution);
                ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                ObjectInputStream in = new ObjectInputStream(bis);
                distribution_copy = (Distribution) in.readObject();
            } catch (IOException | ClassNotFoundException e) {
                System.err.println("Could not copy the distribution in the setRead method of CacheExamples.java");
                e.printStackTrace();
            }
            this.popularitySet(itemclass.getIndex(), jobClass.getIndex() - 1, distribution_copy);
            if(distribution_copy.getSupport().getRight() != itemclass.getNumberOfItems()){
                throw new RuntimeException("The reference model is defined on a number of items different from the ones used to instantiate " + this.name);
            }
            if(distribution instanceof Zipf){
                this.popularityGet(itemclass.getIndex(), jobClass.getIndex() - 1).setParam(2, "n", itemclass.getNumberOfItems());
            }
        } else {
            throw new RuntimeException("A discrete popularity distribution is required.");
        }
    }

    public void setReadItemEntry(JobClass jobClass, Distribution popularity, int cardinality){
        if(popularity.isDiscrete()){
            this.cacheServer.inputJobClasses.put(jobClass.getIndex(), new CacheClassSwitcher.InputJobClassesObj(jobClass, this.schedPolicy, DropStrategy.WaitingQueue));
            Distribution popularity_copy = null;
            try{
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);
                out.writeObject(popularity);
                ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                ObjectInputStream in = new ObjectInputStream(bis);
                popularity_copy = (Distribution) in.readObject();
            } catch (IOException | ClassNotFoundException e) {
                System.err.println("Could not copy the distribution in the setReadItemEntry method of CacheExamples.java");
                e.printStackTrace();
            }
            this.popularitySet(jobClass.getIndex(), popularity_copy);
            if(popularity instanceof Zipf){
                this.popularity.get(jobClass.getIndex()).setParam(2, "n", cardinality);
            }
        } else {
            throw new RuntimeException("A discrete popularity distribution is required.");
        }
    }

    public void setAccessProb(Matrix[][] R){
        this.accessProb = R;
    }

    @Override
    public void setProbRouting(JobClass jobClass, Node destination, double probability){
        setRouting(jobClass, RoutingStrategy.PROB, destination, probability);
    }

    public Section[] getSections(){
        Section[] ret = new Section[3];
        ret[0] = this.input;
        ret[1] = this.cacheServer;
        ret[2] = this.output;
        return ret;
    }

    /**
     * For an incoming job of class r, HITCLASS[r] is the new class of that job after a hit
     * @return - the matrix of hit classes
     */
    public Matrix getHitClass(){
        return this.cacheServer.hitClass;
    }

    /**
     * For an incoming job of class r, MISSCLASS[r] is the new class of that job after a miss
     * @return - the matrix of miss classes
     */
    public Matrix getMissClass(){
        return this.cacheServer.missClass;
    }

//    public CellArray<Distribution> getPopularity() {
//        return popularity;
//    }

    public ItemSet getItems() {
        return items;
    }

    public Matrix[] getGraph() {
        return graph;
    }

    public int getnLevels() {
        return nLevels;
    }

    public CacheClassSwitcher getCacheServer() {
        return cacheServer;
    }

    public Matrix getItemLevelCap() {
        return itemLevelCap;
    }

    public ReplacementStrategy getReplacementPolicy() {
        return replacementPolicy;
    }

    public void popularitySet(int i, Distribution o){
        if(this.popularityRows == 0 && this.popularityColumns == 0){
            // New cell array
            popularitySet(0, i, o);
        } else {
            popularitySet(i % this.popularityRows, i / this.popularityColumns, o);
        }
    }

    public void popularitySet(int i, int j, Distribution o){
        if(i >= this.popularityRows){
            this.popularityRows = i + 1;
        }
        if(j >= this.popularityColumns){
            this.popularityColumns = j + 1;
        }
        this.popularity.put(new PopularityKey(i, j), o);
    }

    public Distribution popularityGet(int i){
        return this.popularityGet(i % this.popularityRows, i / this.popularityColumns);
    }

    public Distribution popularityGet(int i, int j){
        return this.popularity.getOrDefault(new PopularityKey(i, j), null);
    }

    public int popularityLength(){
        return (int) Maths.max(this.popularityRows, this.popularityColumns);
    }

    public static class PopularityKey implements Serializable{

        private final int x;
        private final int y;

        public PopularityKey(int x, int y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof PopularityKey)) return false;
            PopularityKey key = (PopularityKey) o;
            return x == key.x && y == key.y;
        }

        @Override
        public int hashCode() {
            return Objects.hash(x, y);
        }
    }
}
