/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Zipf;
import jline.lang.sections.*;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.io.*;
import java.util.*;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

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
    private final int totalCacheCapacity;
    private int retrievalSystemCapacity;
    private final Map<Integer, List<Integer>> retrievalSystemQueueIndices;
    private final ReplacementStrategy replcStrategy;
    private double admissionProb = 1.0; // q-LRU admission probability on a miss
    private Matrix itemSize;            // per-item storage cost (size); null when unset
    private Matrix costCap;             // per-list storage cost cap; null when unset
    private boolean costCapGlobal;      // true when costCap came from a single cache-wide cap
    private final Map<PopularityKey, Distribution> popularity;
    private final Matrix[] graph;
    private final CacheClassSwitcher cacheServer;
    private final Set<Integer> retrievalClassIndices;
    // Deferred retrieval-system routing edges (mirrors MATLAB retrievalRoutingEntries).
    // The auto-generated retrieval classes are not part of the user-supplied routing
    // matrix P; their edges (cache->queue, queue->queue, queue->cache) are recorded
    // here and injected into P by Network.link before P.setRouting is applied. Later
    // entries override earlier ones for the same (from,to,src,dst).
    private final List<RetrievalRoutingEntry> retrievalRoutingEntries;
    /** Per-item classes of a cache network, keyed by the chain's first read class. */
    private final Map<Integer, int[]> itemClasses;
    /** Item each per-item class reads (1-based), keyed by 0-based class index. */
    private final Map<Integer, Integer> itemOfClass;
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
        this.totalCacheCapacity = (int) itemLevelCap.elementSum();
        this.retrievalSystemCapacity = 0; // By default there is no retrieval system
        if (this.totalCacheCapacity > nitems) {
            throw new RuntimeException("The number of items is smaller than the capacity of " + name);
        }
        this.retrievalSystemQueueIndices = new HashMap<>();
        this.replcStrategy = replPolicy;
        this.cacheServer = new CacheClassSwitcher(classes, nitems, itemLevelCap);
        this.retrievalClassIndices = new HashSet<>();
        this.retrievalRoutingEntries = new ArrayList<>();
        this.itemClasses = new HashMap<Integer, int[]>();
        this.itemOfClass = new HashMap<Integer, Integer>();
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
     * Gets the delayed-hit fraction per class (retrieval system); empty when the
     * cache has no retrieval system.
     *
     * @return A matrix containing the delayed-hit fraction for each job class
     */
    public Matrix getDelayedHitRatio() {
        return this.cacheServer.actualDelayedHitProb;
    }

    /**
     * Gets the per-class, per-list (per-level) hit fraction matrix
     * [classes x lists]; empty when not computed by the solver.
     *
     * @return A matrix of per-list hit fractions
     */
    public Matrix getHitRatioByList() {
        return this.cacheServer.actualHitProbList;
    }

    /**
     * Gets the per-item occupancy matrix [items x (lists+1)]: column 0 is the
     * miss probability, columns 1.. the per-list probabilities; empty when not
     * computed by the solver.
     *
     * @return A matrix of per-item, per-list occupancy probabilities
     */
    public Matrix getItemProb() {
        return this.cacheServer.actualItemProb;
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
     * Gets the total capacity of the cache.
     *
     * @return The total capacity of the cache (sum of capacities of all levels)
     */
    public int getTotalCacheCapacity() {
        return this.totalCacheCapacity;
    }

    /**
     * Gets the retrieval system capacity of the cache.
     * <p>
     *     By default the cache has no retrieval system, so it will be 0. However, if setRetrievalSystem is called with
     *     any arrival class, then the cache state must contain a list of items at the retrieval system.
     * </p>
     *
     * @return The total number of items that can be in the retrieval system at any one time.
     */
    public int getRetrievalSystemCapacity() {
        return this.retrievalSystemCapacity;
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
     * For an incoming job of class r, RETRIEVALCLASSES[i, r] is the new class of the job after a retrieval
     * begins
     *
     * @return - the matrix of retrieval classes
     */
    public Matrix getRetrievalClasses() {
        return this.cacheServer.retrievalClasses;
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
     * Gets the actual expected latency of an item request for each job class.
     *
     * @return A matrix containing the expected latency for each job class
     */
    public Matrix getResidT() {
        return this.cacheServer.actualResidT;
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
     * Sets the storage cost (size) of each item, a positive integer vector with
     * one entry per item. Used together with setCostCaps to bound the storage
     * held by each cache list.
     *
     * @param sizes the per-item storage costs, or a single value for every item
     */
    public void setItemSizes(Matrix sizes) {
        int nitems = this.items.getNumberOfItems();
        Matrix sz;
        if (sizes.length() == 1) {
            sz = new Matrix(1, nitems);
            for (int i = 0; i < nitems; i++) {
                sz.set(0, i, sizes.get(0));
            }
        } else {
            sz = new Matrix(1, sizes.length());
            for (int i = 0; i < sizes.length(); i++) {
                sz.set(0, i, sizes.get(i));
            }
        }
        if (sz.length() != nitems) {
            throw new IllegalArgumentException("The item size vector of " + this.getName()
                    + " must have one entry per item (" + nitems + ").");
        }
        for (int i = 0; i < nitems; i++) {
            double v = sz.get(i);
            if (v <= 0 || v != Math.rint(v)) {
                throw new IllegalArgumentException("Item sizes must be positive integers.");
            }
        }
        this.itemSize = sz;
    }

    /**
     * Gets the per-item storage costs (sizes).
     *
     * @return the per-item storage costs, null when unset
     */
    public Matrix getItemSizes() {
        return this.itemSize;
    }

    /**
     * Sets the per-list cap on the total storage cost of the resident items. A
     * single value declares one cap for the whole cache, modelled as the same
     * cap on every list.
     *
     * @param caps the per-list storage cost caps
     */
    public void setCostCaps(Matrix caps) {
        int h = this.itemLevelCap.length();
        Matrix cp;
        if (caps.length() == 1) {
            this.costCapGlobal = true;
            cp = new Matrix(1, h);
            for (int j = 0; j < h; j++) {
                cp.set(0, j, caps.get(0));
            }
        } else {
            this.costCapGlobal = false;
            cp = new Matrix(1, caps.length());
            for (int j = 0; j < caps.length(); j++) {
                cp.set(0, j, caps.get(j));
            }
        }
        if (cp.length() != h) {
            throw new IllegalArgumentException("The cost cap vector of " + this.getName()
                    + " must have one entry per cache list (" + h + ").");
        }
        for (int j = 0; j < h; j++) {
            double v = cp.get(j);
            if (v < 0 || v != Math.rint(v)) {
                throw new IllegalArgumentException("Storage cost caps must be non-negative integers.");
            }
        }
        if (this.replcStrategy == ReplacementStrategy.CLIMB) {
            throw new IllegalArgumentException("Storage cost caps are not supported with the CLIMB replacement strategy.");
        }
        this.costCap = cp;
    }

    /**
     * Gets the per-list storage cost caps.
     *
     * @return the per-list storage cost caps, null when unset
     */
    public Matrix getCostCaps() {
        return this.costCap;
    }

    /**
     * Reports whether the cost caps were declared as a single cache-wide cap.
     *
     * @return true when a single cap was declared
     */
    public boolean isCostCapGlobal() {
        return this.costCapGlobal;
    }

    /**
     * Sets the mean storage cost held by each list, as computed by a solver.
     *
     * @param listCost the mean per-list storage cost
     */
    public void setResultListCost(Matrix listCost) {
        this.cacheServer.actualListCost = listCost;
    }

    /**
     * Gets the mean storage cost held by each list.
     *
     * @return the mean per-list storage cost, null when not computed
     */
    public Matrix getListCost() {
        return this.cacheServer.actualListCost;
    }

    /**
     * Probability q in [0,1] of admitting a missed item into the cache (q-LRU).
     * Only used when the replacement strategy is QLRU.
     */
    public void setAdmissionProb(double q) {
        if (q < 0 || q > 1) {
            throw new IllegalArgumentException("The admission probability q must lie in [0,1].");
        }
        this.admissionProb = q;
    }

    public double getAdmissionProb() {
        return this.admissionProb;
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
     * Get the node indices of the queues that comprise the retrieval system for a specific job class
     *
     * @param jobClass: The job class which read requests arrive in before routing to the retrieval system
     * @return List containing all node indices for the retrieval system.
     */
    public List<Integer> getRetrievalSystemQueueIndices(int jobClass) {
        return this.retrievalSystemQueueIndices.get(jobClass);
    }

    /**
     * Get the map of arrival class to node indices of the queues that comprise the retrieval system
     *
     * @return Map containing all node indices for all retrieval systems
     */
    public Map<Integer, List<Integer>> getRetrievalSystemQueueIndices() {
        return this.retrievalSystemQueueIndices;
    }

    /**
     * Gets the set of indices of job classes that have completed retrieval
     *
     * @return The set of indices of job classes that have completed retrieval
     */
    public Set<Integer> getRetrievalClassIndices() {
        return this.retrievalClassIndices;
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
        this.cacheServer.actualDelayedHitProb = new Matrix(0, 0);
        this.cacheServer.actualHitProbList = new Matrix(0, 0);
        this.cacheServer.actualItemProb = new Matrix(0, 0);
        this.cacheServer.actualResidT = new Matrix(0, 0);
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
     * <p>When a job of class jobinClass hits in the cache, it will be
     * transformed to class joboutclass.</p>
     * 
     * @param jobinClass The incoming job class
     * @param joboutclass The job class after a cache hit
     */
    public void setHitClass(JobClass jobinClass, JobClass joboutclass) {
        int jobinindex = jobinClass.getIndex() - 1;
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
     * <p>When a job of class jobinClass misses in the cache, it will be
     * transformed to class joboutclass.</p>
     * 
     * @param jobinClass The incoming job class
     * @param joboutclass The job class after a cache miss
     */
    public void setMissClass(JobClass jobinClass, JobClass joboutclass) {
        int jobinindex = jobinClass.getIndex() - 1;
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
     * Sets the output class for jobs that experience a cache miss and need to be retrieved
     *
     * <p>When a job of class jobinClass misses in the cache, it is routed to the retrieval system
     * with class joboutclass. This method sets the output class for a specific item.</p>
     *
     * @param jobinClass The incoming job class
     * @param joboutclass The outgoing job class
     * @param item The item to be retrieved
     */
    public void setRetrievalClass(JobClass jobinClass, JobClass joboutclass, int item) {
        int jobinindex = jobinClass.getIndex() - 1;
        int joboutindex = joboutclass.getIndex() - 1;

        if (this.model.getClasses().size() > this.cacheServer.retrievalClasses.getNumCols()) {
            int nItems = this.cacheServer.retrievalClasses.getNumRows();
            int newCols = this.model.getClasses().size();
            Matrix newRetrievalClasses = new Matrix(nItems, newCols);
            newRetrievalClasses.fill(-1);
            for (int r = 0; r < this.cacheServer.retrievalClasses.getNumRows(); r++) {
                for (int c = 0; c < this.cacheServer.retrievalClasses.getNumCols(); c++) {
                    newRetrievalClasses.set(r, c, this.cacheServer.retrievalClasses.get(r, c));
                }
            }
            this.cacheServer.retrievalClasses = newRetrievalClasses;
        }
        this.cacheServer.retrievalClasses.set(item, jobinindex, joboutindex);
    }

    /**
     * Reconstruct the retrieval-system bookkeeping on a cache whose retrieval
     * classes, routing and queue service already exist in the model (e.g. after
     * deserialization). Unlike {@link #setRetrievalSystem}, this does NOT create
     * new job classes or routing; it only restores the cache-internal state that
     * the simulator reads: the retrieval-system capacity, the per-arrival-class
     * retrieval queue indices, the retrieval-class index set, and the per-item
     * retrieval-class mapping.
     *
     * @param jobinClass     arrival job class routed through the retrieval system
     * @param queueIndices   node indices of the retrieval-system queues
     * @param retrievalClassByItem retrieval class for each item (length = nItems)
     */
    public void attachRetrievalSystem(JobClass jobinClass, List<Integer> queueIndices,
                                      JobClass[] retrievalClassByItem) {
        int nItems = this.items.getNumberOfItems();
        this.retrievalSystemCapacity = nItems - this.totalCacheCapacity;
        this.retrievalSystemQueueIndices.put(jobinClass.getIndex() - 1, queueIndices);
        for (int i = 0; i < retrievalClassByItem.length; i++) {
            JobClass rc = retrievalClassByItem[i];
            if (rc != null) {
                setRetrievalClass(jobinClass, rc, i);
                this.retrievalClassIndices.add(rc.getIndex() - 1);
            }
        }
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
    @Override
    public void removeJobClass(JobClass jobClass) {
        throw new UnsupportedOperationException(
            "Cannot dynamically remove classes in models with caches. You need to re-instantiate the model.");
    }

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

    private boolean isCacheOrRetrievalNode(JobClass jobinClass, Node node) {
        if (node == this) {
            return true;
        }
        List<Integer> queueIndices = retrievalSystemQueueIndices.get(jobinClass.getIndex() - 1);
        return queueIndices == null || queueIndices.contains(node.getNodeIndex());
    }

    /**
     * A deferred retrieval-system routing edge, injected into the routing matrix by
     * Network.link. Mirrors a MATLAB retrievalRoutingEntries tuple.
     */
    public static class RetrievalRoutingEntry implements Serializable {
        public final JobClass fromClass;
        public final JobClass toClass;
        public final Node srcNode;
        public final Node destNode;
        public final double prob;
        public RetrievalRoutingEntry(JobClass fromClass, JobClass toClass, Node srcNode, Node destNode, double prob) {
            this.fromClass = fromClass;
            this.toClass = toClass;
            this.srcNode = srcNode;
            this.destNode = destNode;
            this.prob = prob;
        }
    }

    public List<RetrievalRoutingEntry> getRetrievalRoutingEntries() {
        return this.retrievalRoutingEntries;
    }

    /**
     * Register a deferred routing edge for an auto-generated retrieval class. With
     * allowZero=true an explicit prob==0 entry is recorded so it can override (delete)
     * a default edge inherited from the read class; negative probabilities are always
     * dropped and the internal broadcast path keeps allowZero=false (zeros = no edge).
     */
    public void addRetrievalRoutingEntry(JobClass fromClass, JobClass toClass, Node srcNode, Node destNode,
                                         double prob, boolean allowZero) {
        if (prob < 0 || (prob == 0 && !allowZero)) {
            return;
        }
        this.retrievalRoutingEntries.add(new RetrievalRoutingEntry(fromClass, toClass, srcNode, destNode, prob));
    }

    /**
     * Sets the routing probability for an item retrieval transition between two nodes of the
     * retrieval system. A node is either one of the retrieval queues or the cache itself: pass
     * the cache as source for a cache-&gt;queue entry, or as destination for a queue-&gt;cache exit.
     *
     * @param jobinClass The job class for a request that can be routed through this system.
     * @param item The item whose retrieval routing is being updated.
     * @param source The node the retrieval departs from (a retrieval queue or the cache).
     * @param dest The node the retrieval is routed to (a retrieval queue or the cache).
     * @param probability The probability of routing the retrieval from source to dest.
     */
    public void setItemRoutingProbability(JobClass jobinClass, int item, Node source, Node dest,
                                          double probability) {
        int jobinClassIdx = jobinClass.getIndex() - 1;

        if (!isCacheOrRetrievalNode(jobinClass, source) || !isCacheOrRetrievalNode(jobinClass, dest)) {
            throw new RuntimeException("Node is not the cache or a queue in the retrieval system for " +
                    jobinClass.getName() + " at " + name);
        }

        int retrievalClassIdx = (int) this.getRetrievalClasses().get(item, jobinClassIdx);
        JobClass retrievalClass = this.model.getJobClassFromIndex(retrievalClassIdx);

        this.addRetrievalRoutingEntry(retrievalClass, retrievalClass, source, dest, probability, true);
    }

    /** Short alias for {@link #setItemRoutingProbability}. */
    public void setItemRoutingProb(JobClass jobinClass, int item, Node source, Node dest, double probability) {
        setItemRoutingProbability(jobinClass, item, source, dest, probability);
    }

    /**
     * Initializes the retrieval system when a job of class jobinClass arrives at the cache node.
     * Neither service rates nor routing matrices are initialized and both must be set later.
     *
     * <p>
     *     To set service rates for an item at a queue call queue.setItemServiceRate(Cache cache,
     *     JobClass jobinClass, int item, double rate).
     *     To set the retrieval routing for an item call cacheNode.setItemRoutingProb(JobClass jobinClass,
     *     int item, Node source, Node dest, double probability). source/dest are a retrieval queue or the
     *     cache itself: pass the cache as source for a cache->queue entry, or as dest for a queue->cache exit.
     * </p>
     *
     * @param jobinClass The job class for a request that can be routed through this system.
     * @param missClass The job class for retrieval misses.
     * @param queues Array of queues in the retrieval system.
     */
    /**
     * Normalise a hit/miss class argument: one class shared by every item, or one per item.
     */
    private JobClass[] perItemClasses(JobClass[] spec, int nItems, String what) {
        if (spec.length == nItems) {
            return spec;
        }
        if (spec.length == 1) {
            JobClass[] out = new JobClass[nItems];
            for (int i = 0; i < nItems; i++) {
                out[i] = spec[0];
            }
            return out;
        }
        throw new RuntimeException(spec.length + " " + what + " classes were given for " + nItems
                + " items; pass one class per item or a single class shared by all.");
    }

    /**
     * Record that a class reads a given item (1-based), as read back from JSON.
     */
    public void setItemOfClass(JobClass jobClass, int item) {
        this.itemOfClass.put(jobClass.getIndex() - 1, item);
    }

    /**
     * Item each per-item class reads (1-based), 0 where the class is not one.
     */
    public int getItemOfClass(int classIdx0) {
        Integer v = this.itemOfClass.get(classIdx0);
        return v == null ? 0 : v.intValue();
    }

    /**
     * Declare that readClasses[i] is the request stream for item i at this cache. Use at
     * the cache the exogenous requests enter, where the per-item classes are the user's
     * own; popularity is then carried by the per-class request rates. Keeps a cache
     * network free of arc-level class switching, so no class acquires a default route
     * into the cache the model never intended.
     */
    public void setItemReadClasses(JobClass[] readClasses, JobClass[] hitClasses) {
        int nItems = this.items.getNumberOfItems();
        if (readClasses.length != nItems) {
            throw new RuntimeException(this.name + " holds " + nItems + " items but "
                    + readClasses.length + " read classes were given; pass exactly one class per item.");
        }
        JobClass[] hitArr = perItemClasses(hitClasses, nItems, "hit");
        int[] idx = new int[nItems];
        for (int i = 0; i < nItems; i++) {
            JobClass cls = readClasses[i];
            double[] pop = new double[nItems];
            pop[i] = 1.0;
            this.popularitySet(this.items.getIndex(), cls.getIndex() - 1,
                    new DiscreteSampler(new Matrix(pop)));
            setHitClass(cls, hitArr[i]);
            idx[i] = cls.getIndex();
            this.itemOfClass.put(cls.getIndex() - 1, i + 1);
        }
        this.itemClasses.put(readClasses[0].getIndex(), idx);
    }

    /**
     * Mint one class per item at a cache FED BY ANOTHER CACHE, so item identity survives
     * the miss hop. Idempotent.
     */
    public int[] setItemClasses(JobClass jobinClass, JobClass[] hitClasses) {
        int[] existing = this.itemClasses.get(jobinClass.getIndex());
        if (existing != null) {
            return existing;
        }
        int nItems = this.items.getNumberOfItems();
        JobClass[] hitArr = perItemClasses(hitClasses, nItems, "hit");
        int[] idx = new int[nItems];
        for (int i = 0; i < nItems; i++) {
            JobClass cls;
            if (jobinClass instanceof ClosedClass) {
                cls = new ClosedClass(this.model, this.name + "_item" + (i + 1), 0,
                        jobinClass.getReferenceStation(), 0);
            } else {
                cls = new OpenClass(this.model, this.name + "_item" + (i + 1));
            }
            double[] pop = new double[nItems];
            pop[i] = 1.0;
            this.popularitySet(this.items.getIndex(), cls.getIndex() - 1,
                    new DiscreteSampler(new Matrix(pop)));
            setHitClass(cls, hitArr[i]);
            idx[i] = cls.getIndex();
            this.itemOfClass.put(cls.getIndex() - 1, i + 1);
        }
        this.itemClasses.put(jobinClass.getIndex(), idx);
        return idx;
    }

    /**
     * Send this cache's misses to nextCache preserving item identity: the miss class of
     * this cache for item i IS the read class of nextCache for item i.
     */
    public JobClass[] setMissCache(JobClass jobinClass, Cache nextCache, JobClass[] hitClassAtNext) {
        if (nextCache.getItems().getNumberOfItems() != this.items.getNumberOfItems()) {
            throw new RuntimeException("Cache " + this.name + " holds "
                    + this.items.getNumberOfItems() + " items but " + nextCache.getName() + " holds "
                    + nextCache.getItems().getNumberOfItems()
                    + "; a cache network requires one common item set.");
        }
        int[] selfIdx = this.itemClasses.get(jobinClass.getIndex());
        if (selfIdx == null) {
            throw new RuntimeException("No per-item classes at " + this.name + " for class "
                    + jobinClass.getName() + "; call setItemReadClasses before setMissCache.");
        }
        int[] nextIdx = nextCache.setItemClasses(jobinClass, hitClassAtNext);
        List<JobClass> all = this.model.getClasses();
        JobClass[] nextClasses = new JobClass[selfIdx.length];
        for (int i = 0; i < selfIdx.length; i++) {
            JobClass selfCls = all.get(selfIdx[i] - 1);
            JobClass nextCls = all.get(nextIdx[i] - 1);
            setMissClass(selfCls, nextCls);
            addRetrievalRoutingEntry(nextCls, nextCls, this, nextCache, 1.0, false);
            nextClasses[i] = nextCls;
        }
        return nextClasses;
    }

    /**
     * Terminate a cache network: every per-item class of this cache reports a miss as the
     * matching entry of missClasses, which the user routes onward.
     */
    public void setItemMissClass(JobClass jobinClass, JobClass[] missClasses) {
        int[] selfIdx = this.itemClasses.get(jobinClass.getIndex());
        if (selfIdx == null) {
            throw new RuntimeException("No per-item classes at " + this.name + " for class "
                    + jobinClass.getName() + "; call setItemClasses before setItemMissClass.");
        }
        JobClass[] missArr = perItemClasses(missClasses, selfIdx.length, "miss");
        List<JobClass> all = this.model.getClasses();
        for (int i = 0; i < selfIdx.length; i++) {
            setMissClass(all.get(selfIdx[i] - 1), missArr[i]);
        }
    }

    public void setRetrievalSystem(JobClass jobinClass, JobClass missClass, Queue[] queues) {
        int nQueues = queues.length;
        int nItems = this.items.getNumberOfItems();
        this.retrievalSystemCapacity = nItems - this.totalCacheCapacity;

        if (nQueues == 0) {
            throw new RuntimeException("Retrieval system cannot be initialised with no stations " + name);
        }

        // Inherit the read class's service distribution at each queue as the per-item default.
        Distribution[] readService = new Distribution[nQueues];
        for (int q = 0; q < nQueues; q++) {
            Distribution d = queues[q].getService(jobinClass);
            if (d == null || d instanceof jline.lang.processes.Disabled) {
                throw new RuntimeException("No service distribution for the read class at queue " +
                        queues[q].getName() + "; call queue.setService(readClass, ...) before setRetrievalSystem.");
            }
            readService[q] = d;
        }

        List<Integer> queueIndices = new ArrayList<Integer>(nQueues);
        for (Queue queue : queues) {
            queueIndices.add(queue.getNodeIndex());
        }
        retrievalSystemQueueIndices.put(jobinClass.getIndex() - 1, queueIndices);

        List<JobClass> retrievalClassesList = new ArrayList<JobClass>();
        for (int i = 0; i < nItems; i++) {
            JobClass retrievalClass;
            if (jobinClass instanceof ClosedClass) {
                Station refStation = jobinClass.getReferenceStation();
                retrievalClass = new ClosedClass(this.model, jobinClass.getName() + "_retrievalClass_" + (i + 1), 0, refStation, 0);
            } else {
                retrievalClass = new OpenClass(this.model, jobinClass.getName() + "_retrievalClass_" + (i + 1));
            }
            retrievalClassesList.add(retrievalClass);
            this.retrievalClassIndices.add(retrievalClass.getIndex() - 1);
        }

        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        for (int i = 0; i < nItems; i++) {
            JobClass retrievalClass = retrievalClassesList.get(i);

            // When an item needs to be retrieved, a job of jobinClass switches to retrievalClass
            setRetrievalClass(jobinClass, retrievalClass, i);

            // On arrival of the retrieval, it will always trigger a read of item i
            double[] itemPopularity = new double[this.getNumberOfItems()];
            itemPopularity[i] = 1.0;
            this.popularitySet(this.items.getIndex(), retrievalClass.getIndex() - 1, new DiscreteSampler(new Matrix(itemPopularity)));

            // After the arrival of the retrieval, the job will switch to a miss
            setMissClass(retrievalClass, missClass);

            for (int q = 0; q < nQueues; q++) {
                queues[q].setService(retrievalClass, readService[q]);
                queues[q].setClassCap(retrievalClass, 1);
            }
        }
    }

    /**
     * Sets the actual hit probabilities from simulation or analysis results.
     *
     * @param actualHitProb Matrix containing the observed hit probabilities
     */
    /**
     * Sets the per-item delayed-hit queue length: the mean number of secondary
     * requests waiting on the in-flight fetch of each item (d1), and the same count
     * including the request that triggered the fetch (dfull).
     *
     * @param d1 mean merged secondary requests per item
     * @param dfull d1 plus the triggering request
     */
    public void setResultDelayedHitQLen(Matrix d1, Matrix dfull) {
        this.cacheServer.actualDelayedHitQLen = d1;
        this.cacheServer.actualDelayedHitQLenFull = dfull;
    }

    /**
     * Per-item delayed-hit queue length, empty when the solver does not compute it.
     *
     * @return mean merged secondary requests per item
     */
    public Matrix getDelayedHitQLen() {
        return this.cacheServer.actualDelayedHitQLen;
    }

    /**
     * Per-item delayed-hit queue length including the triggering request.
     *
     * @return d1 plus the triggering request, per item
     */
    public Matrix getDelayedHitQLenFull() {
        return this.cacheServer.actualDelayedHitQLenFull;
    }

    public void setResultHitProb(Matrix actualHitProb) {
        this.cacheServer.actualHitProb = actualHitProb;
    }

    /**
     * Sets the delayed-hit fraction per class (retrieval system).
     *
     * @param actualDelayedHitProb Matrix containing the delayed-hit fractions
     */
    public void setResultDelayedHitProb(Matrix actualDelayedHitProb) {
        this.cacheServer.actualDelayedHitProb = actualDelayedHitProb;
    }

    /**
     * Sets the per-class, per-list (per-level) hit fraction matrix
     * [classes x lists]; rows sum to the aggregate hit fraction.
     *
     * @param actualHitProbList Matrix of per-list hit fractions
     */
    public void setResultHitProbList(Matrix actualHitProbList) {
        this.cacheServer.actualHitProbList = actualHitProbList;
    }

    /**
     * Sets the per-item occupancy matrix [items x (lists+1)]; column 0 = miss,
     * columns 1.. = per-list.
     *
     * @param actualItemProb Matrix of per-item, per-list occupancy probabilities
     */
    public void setResultItemProb(Matrix actualItemProb) {
        this.cacheServer.actualItemProb = actualItemProb;
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
     * Sets the actual expected latency from simulation or analysis results.
     *
     * @param actualResidT Matrix containing the observed expected latencies
     */
    public void setResultResidT(Matrix actualResidT) {
        this.cacheServer.actualResidT = actualResidT;
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
