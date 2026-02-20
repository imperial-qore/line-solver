/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.ItemSet;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.*;
import jline.lang.processes.Distribution;
import jline.lang.processes.Zipf;
import jline.lang.sections.Buffer;
import jline.lang.sections.CacheClassSwitcher;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.Section;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.io.*;
import java.util.*;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A cache node that implements cache replacement policies and class switching based on cache hits and misses.
 * 
 * <p>The Cache node models a caching system where incoming jobs request items from a finite cache.
 * When a requested item is found in the cache (hit), the job may be routed to one class; when the
 * item is not found (miss), it may be routed to a different class. This enables modeling of
 * cache-aware queueing networks where performance depends on cache hit rates.</p>
 * 
 * <p>Key features:
 * <ul>
 *   <li>Multi-level cache support with configurable capacity per level</li>
 *   <li>Various replacement strategies (LRU, FIFO, RANDOM, etc.)</li>
 *   <li>Popularity-based item access patterns (e.g., Zipf distribution)</li>
 *   <li>Class switching based on hit/miss outcomes</li>
 *   <li>Optional graph-based cache structures</li>
 * </ul>
 * </p>
 * 
 * @see ItemSet
 * @see ReplacementStrategy
 * @see CacheClassSwitcher
 * @since 1.0
 */
public class Cache extends StatefulNode implements Serializable {
    private final ItemSet items;
    private final int nLevels;
    private final int cap;
    private final Matrix itemLevelCap;
    private final ReplacementStrategy replcStrategy;
    private final Map<PopularityKey, Distribution> popularity;
    private final Matrix[] graph;
    private final CacheClassSwitcher cacheServer;
    public Matrix[][] accessProb;
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;
    private int popularityRows;
    private int popularityColumns;

    /**
     * Creates a single-level cache with the specified item capacity and replacement policy.
     * 
     * @param model The network model this cache belongs to
     * @param name The name of the cache node
     * @param nitems The total number of items that can be requested from this cache
     * @param itemLevelCap The capacity of the cache (number of items it can hold)
     * @param replPolicy The replacement strategy to use when the cache is full
     */
    public Cache(Network model, String name, int nitems, int itemLevelCap, ReplacementStrategy replPolicy) {
        this(model, name, nitems, Matrix.singleton(itemLevelCap), replPolicy, null);
    }

    /**
     * Creates a single-level cache with the specified item capacity, replacement policy, and graph structure.
     * 
     * @param model The network model this cache belongs to
     * @param name The name of the cache node
     * @param nitems The total number of items that can be requested from this cache
     * @param itemLevelCap The capacity of the cache (number of items it can hold)
     * @param replPolicy The replacement strategy to use when the cache is full
     * @param graph Optional graph structure defining cache organization
     */
    public Cache(Network model, String name, int nitems, int itemLevelCap, ReplacementStrategy replPolicy, Matrix[] graph) {
        this(model, name, nitems, Matrix.singleton(itemLevelCap), replPolicy, graph);
    }

    /**
     * Creates a multi-level cache with different capacities per level and a replacement policy.
     * 
     * @param model The network model this cache belongs to
     * @param name The name of the cache node
     * @param nitems The total number of items that can be requested from this cache
     * @param itemLevelCap A matrix specifying the capacity of each cache level
     * @param replPolicy The replacement strategy to use when the cache is full
     */
    public Cache(Network model, String name, int nitems, Matrix itemLevelCap, ReplacementStrategy replPolicy) {
        this(model, name, nitems, itemLevelCap, replPolicy, null);
    }

    /**
     * Creates a multi-level cache with different capacities per level, replacement policy, and graph structure.
     * 
     * <p>This is the main constructor that all other constructors delegate to. It initializes
     * the cache with all necessary components including input/output buffers, the cache server,
     * and item management.</p>
     * 
     * @param model The network model this cache belongs to
     * @param name The name of the cache node
     * @param nitems The total number of items that can be requested from this cache
     * @param itemLevelCap A matrix specifying the capacity of each cache level
     * @param replPolicy The replacement strategy to use when the cache is full
     * @param graph Optional graph structure defining cache organization
     * @throws RuntimeException if the total item capacity exceeds the number of items
     */
    public Cache(Network model, String name, int nitems, Matrix itemLevelCap, ReplacementStrategy replPolicy, Matrix[] graph) {
        super(name);
        List<JobClass> classes = model.getClasses();
        this.input = new Buffer(classes);
        this.output = new Dispatcher(classes);
        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.items = new ItemSet(model, name + "_Items", nitems, this);
        this.nLevels = itemLevelCap.getNonZeroLength();
        this.cap = Integer.MAX_VALUE; // job capacity
        this.accessProb = null;
        this.itemLevelCap = itemLevelCap; // item capacity
        if (this.itemLevelCap.elementSum() > nitems) {
            throw new RuntimeException("The number of items is smaller than the capacity of " + name);
        }
        this.replcStrategy = replPolicy;
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
     * Gets the access probability matrix for a specific cache level and job class.
     * 
     * @param i The cache level index
     * @param j The job class index
     * @return The access probability matrix for the specified indices
     */
    public Matrix getAccessProb(int i, int j) {
        return this.accessProb[i][j];
    }

    /**
     * Gets the internal cache server that handles class switching logic.
     * 
     * @return The CacheClassSwitcher instance managing hit/miss class transitions
     */
    public CacheClassSwitcher getCacheServer() {
        return cacheServer;
    }

    /**
     * Gets the graph structure defining the cache organization.
     * 
     * @return Array of matrices representing the cache graph structure, or null if not defined
     */
    public Matrix[] getGraph() {
        return graph;
    }

    /**
     * For an incoming job of class r, HITCLASS[r] is the new class of that job after a hit
     *
     * @return - the matrix of hit classes
     */
    public Matrix getHitClass() {
        return this.cacheServer.hitClass;
    }

    /**
     * Gets the actual hit probability/ratio for each job class.
     * 
     * <p>This returns the observed hit rates from simulation or analysis,
     * not the theoretical expected values.</p>
     * 
     * @return A matrix containing the hit ratio for each job class
     */
    public Matrix getHitRatio() {
        return this.cacheServer.actualHitProb;
    }

    /**
     * Gets the capacity configuration for each cache level.
     * 
     * @return A matrix where each element specifies the item capacity of a cache level
     */
    public Matrix getItemLevelCap() {
        return itemLevelCap;
    }

    /**
     * Gets the set of items that can be stored in this cache.
     * 
     * @return The ItemSet object managing the cache items
     */
    public ItemSet getItems() {
        return items;
    }

    /**
     * For an incoming job of class r, MISSCLASS[r] is the new class of that job after a miss
     *
     * @return - the matrix of miss classes
     */
    public Matrix getMissClass() {
        return this.cacheServer.missClass;
    }

    /**
     * Gets the actual miss probability/ratio for each job class.
     * 
     * <p>This returns the observed miss rates from simulation or analysis,
     * not the theoretical expected values.</p>
     * 
     * @return A matrix containing the miss ratio for each job class
     */
    public Matrix getMissRatio() {
        return this.cacheServer.actualMissProb;
    }

    /**
     * Gets the total number of items that can be requested from this cache.
     * 
     * @return The number of distinct items in the item set
     */
    public int getNumberOfItems() {
        return this.items.getNumberOfItems();
    }

    /**
     * Gets the replacement strategy used when the cache is full.
     * 
     * @return The cache replacement strategy (e.g., LRU, FIFO, RANDOM)
     */
    public ReplacementStrategy getReplacementStrategy() {
        return replcStrategy;
    }

    /**
     * Gets the internal sections of this cache node.
     * 
     * <p>Returns the three main sections: input buffer, cache server, and output dispatcher.</p>
     * 
     * @return A list containing the input, server, and output sections
     */
    public List<Section> getSections() {
        List<Section> ret = new ArrayList<>();
        ret.add(this.input);
        ret.add(this.cacheServer);
        ret.add(this.output);
        return ret;
    }

    /**
     * Gets the number of cache levels.
     * 
     * @return The number of levels in this multi-level cache
     */
    public int getnLevels() {
        return nLevels;
    }

    /**
     * Gets the popularity distribution for a linear index.
     * 
     * <p>Converts a linear index to 2D coordinates and retrieves the distribution.</p>
     * 
     * @param i The linear index
     * @return The popularity distribution at the specified index
     */
    public Distribution popularityGet(int i) {
        return this.popularityGet(i % this.popularityRows, i / this.popularityRows);
    }

    /**
     * Gets the popularity distribution for a specific item class and job class.
     * 
     * @param i The item class index
     * @param j The job class index
     * @return The popularity distribution, or null if not set
     */
    public Distribution popularityGet(int i, int j) {
        return this.popularity.getOrDefault(new PopularityKey(i, j), null);
    }

    /**
     * Gets the maximum dimension of the popularity matrix.
     * 
     * @return The maximum of rows and columns in the popularity matrix
     */
    public int popularityLength() {
        return (int) Maths.max(this.popularityRows, this.popularityColumns);
    }

//    public CellArray<Distribution> getPopularity() {
//        return popularity;
//    }

    /**
     * Sets the popularity distribution for a linear index.
     * 
     * <p>Converts a linear index to 2D coordinates and sets the distribution.</p>
     * 
     * @param i The linear index
     * @param o The popularity distribution to set
     */
    public void popularitySet(int i, Distribution o) {
        if (this.popularityRows == 0 && this.popularityColumns == 0) {
            // New cell array
            popularitySet(0, i, o);
        } else {
            popularitySet(i % this.popularityRows, i / this.popularityRows, o);
        }
    }

    /**
     * Sets the popularity distribution for a specific item class and job class.
     * 
     * <p>Automatically expands the popularity matrix dimensions if necessary.</p>
     * 
     * @param i The item class index
     * @param j The job class index
     * @param o The popularity distribution to set
     */
    public void popularitySet(int i, int j, Distribution o) {
        if (i >= this.popularityRows) {
            this.popularityRows = i + 1;
        }
        if (j >= this.popularityColumns) {
            this.popularityColumns = j + 1;
        }
        this.popularity.put(new PopularityKey(i, j), o);
    }

    /**
     * Resets the internal data structures when the network model is reset.
     * 
     * <p>Clears the actual hit and miss probability matrices to prepare for
     * a new simulation or analysis run.</p>
     */
    public void reset() {
        this.cacheServer.actualHitProb = new Matrix(0, 0);
        this.cacheServer.actualMissProb = new Matrix(0, 0);
    }

    /**
     * Sets the access probability matrices for all cache levels and job classes.
     * 
     * @param R A 2D array of matrices containing access probabilities
     */
    public void setAccessProb(Matrix[][] R) {
        this.accessProb = R;
    }

    /**
     * Sets the output class for jobs that experience a cache hit.
     * 
     * <p>When a job of class jobinclass hits in the cache, it will be
     * transformed to class joboutclass.</p>
     * 
     * @param jobinclass The incoming job class
     * @param joboutclass The job class after a cache hit
     */
    public void setHitClass(JobClass jobinclass, JobClass joboutclass) {
        int jobinindex = jobinclass.getIndex() - 1;
        int joboutindex = joboutclass.getIndex() - 1;
        if (this.model.getClasses().size() > this.cacheServer.hitClass.getNumCols()) {
            Matrix newHitClass = new Matrix(1, this.model.getClasses().size());
            newHitClass.fill(-1);
            for (int i = 0; i < this.cacheServer.hitClass.getNumCols(); i++) {
                newHitClass.set(i, this.cacheServer.hitClass.get(i));
            }
            this.cacheServer.hitClass = newHitClass;
        }
        this.cacheServer.hitClass.set(jobinindex, joboutindex);
    }

    /**
     * Sets the output class for jobs that experience a cache miss.
     * 
     * <p>When a job of class jobinclass misses in the cache, it will be
     * transformed to class joboutclass.</p>
     * 
     * @param jobinclass The incoming job class
     * @param joboutclass The job class after a cache miss
     */
    public void setMissClass(JobClass jobinclass, JobClass joboutclass) {
        int jobinindex = jobinclass.getIndex() - 1;
        int joboutindex = joboutclass.getIndex() - 1;
        if (this.model.getClasses().size() > this.cacheServer.missClass.getNumCols()) {
            Matrix newMissClass = new Matrix(1, this.model.getClasses().size());
            newMissClass.fill(-1);
            for (int i = 0; i < this.cacheServer.missClass.getNumCols(); i++) {
                newMissClass.set(i, this.cacheServer.missClass.get(i));
            }
            this.cacheServer.missClass = newMissClass;
        }
        this.cacheServer.missClass.set(jobinindex, joboutindex);
    }

    /**
     * Sets probabilistic routing for a job class to a destination node.
     * 
     * @param jobClass The job class to configure routing for
     * @param destination The destination node
     * @param probability The routing probability
     */
    @Override
    public void setProbRouting(JobClass jobClass, Node destination, double probability) {
        setRouting(jobClass, RoutingStrategy.PROB, destination, probability);
    }

    /**
     * Sets the read policy for a job class using a popularity distribution.
     * 
     * <p>The distribution determines which items are requested by jobs of this class.
     * Common distributions include Zipf for modeling popularity skew.</p>
     *
     * @param jobClass The job class to configure
     * @param distribution The discrete popularity distribution over items
     * @throws RuntimeException if the distribution is not discrete or has wrong support
     */
    public void setRead(JobClass jobClass, Distribution distribution) {
        ItemSet itemclass = this.items;
        if (distribution.isDiscrete()) {
            this.cacheServer.inputJobClasses.put(jobClass.getIndex(), new CacheClassSwitcher.InputJobClassesObj(jobClass, this.schedPolicy, DropStrategy.WaitingQueue));
            Distribution distribution_copy = null;
            try {
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);
                out.writeObject(distribution);
                ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                ObjectInputStream in = new ObjectInputStream(bis);
                distribution_copy = (Distribution) in.readObject();
            } catch (IOException | ClassNotFoundException e) {
                e.printStackTrace();
                line_error(mfilename(new Object() {
                }), "Could not copy the distribution in the setRead method of CacheExamples.java");
            }
            this.popularitySet(itemclass.getIndex(), jobClass.getIndex() - 1, distribution_copy);
            if (distribution_copy.getSupport().getRight() != itemclass.getNumberOfItems()) {
                throw new RuntimeException("The reference model is defined on a number of items different from the ones used to instantiate " + this.name);
            }
            if (distribution instanceof Zipf) {
                this.popularityGet(itemclass.getIndex(), jobClass.getIndex() - 1).setParam(2, "n", itemclass.getNumberOfItems());
            }
        } else {
            throw new RuntimeException("A discrete popularity distribution is required.");
        }
    }

    /**
     * Sets the read policy for a job class with explicit item cardinality.
     * 
     * <p>Similar to setRead but allows specifying the number of items explicitly,
     * useful when the distribution needs to be configured with a specific cardinality.</p>
     * 
     * @param jobClass The job class to configure
     * @param popularity The discrete popularity distribution
     * @param cardinality The number of items in the distribution support
     * @throws RuntimeException if the distribution is not discrete
     */
    public void setReadItemEntry(JobClass jobClass, Distribution popularity, int cardinality) {
        if (popularity.isDiscrete()) {
            this.cacheServer.inputJobClasses.put(jobClass.getIndex(), new CacheClassSwitcher.InputJobClassesObj(jobClass, this.schedPolicy, DropStrategy.WaitingQueue));
            Distribution popularity_copy = null;
            try {
                ByteArrayOutputStream bos = new ByteArrayOutputStream();
                ObjectOutputStream out = new ObjectOutputStream(bos);
                out.writeObject(popularity);
                ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
                ObjectInputStream in = new ObjectInputStream(bis);
                popularity_copy = (Distribution) in.readObject();
            } catch (IOException | ClassNotFoundException e) {
                line_error(mfilename(new Object() {
                }), "Could not copy the distribution in the setReadItemEntry method of CacheExamples.java");
                e.printStackTrace();
            }
            this.popularitySet(jobClass.getIndex() - 1, popularity_copy);
            if (popularity instanceof Zipf) {
                this.popularityGet(jobClass.getIndex() - 1).setParam(2, "n", cardinality);
            }
        } else {
            throw new RuntimeException("A discrete popularity distribution is required.");
        }
    }

    /**
     * Sets the actual hit probabilities from simulation or analysis results.
     * 
     * @param actualHitProb Matrix containing the observed hit probabilities
     */
    public void setResultHitProb(Matrix actualHitProb) {
        this.cacheServer.actualHitProb = actualHitProb;
    }

    /**
     * Sets the actual miss probabilities from simulation or analysis results.
     * 
     * @param actualMissProb Matrix containing the observed miss probabilities
     */
    public void setResultMissProb(Matrix actualMissProb) {
        this.cacheServer.actualMissProb = actualMissProb;
    }

    /**
     * Sets the scheduling strategy for a job class.
     * 
     * <p>Note: Currently this method has no implementation as caches use FCFS scheduling.</p>
     * 
     * @param jobClass The job class index
     * @param strategy The scheduling strategy (unused)
     */
    public void setScheduling(int jobClass, SchedStrategy strategy) {
    }

    /**
     * A key class for storing popularity distributions in a 2D coordinate system.
     * 
     * <p>Used internally to map (item class, job class) pairs to their popularity distributions.</p>
     */
    public static class PopularityKey implements Serializable {

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
