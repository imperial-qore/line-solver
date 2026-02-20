/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.*;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.processes.DiscreteDistribution;
import jline.util.Pair;
import jline.lang.reward.RewardFunction;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.api.sn.SnPrintKt.snPrint;

/**
 * Class summarizing the characteristics of a Network object
 */
public class NetworkStruct implements Copyable, Cloneable {
    //For data structure, {} is represented by HashMap, [] is represented by ArrayList, [][] is represented by matrix;
    //For the matrix that stores Constant. Use double list instead.
    public int nstations;
    public int nstateful;
    public int nnodes;
    public int nclasses;
    public int nclosedjobs;
    public int nchains;

    public Map<JobClass, Map<JobClass, Matrix>> rtorig;
    public Map<Station, Map<JobClass, SerializableFunction<Double, Double>>> lst;
    public Map<StatefulNode, Matrix> state;
    public Map<StatefulNode, Matrix> stateprior;
    public Map<StatefulNode, Matrix> space;
    public Map<StatefulNode, Map<String, Integer>> spaceHash;
    public Map<Node, Map<JobClass, RoutingStrategy>> routing;
    public Map<Station, Map<JobClass, ProcessType>> procid;
    public Map<Station, Map<JobClass, Matrix>> mu;
    public Map<Station, Map<JobClass, Matrix>> phi;
    public Map<Station, Map<JobClass, MatrixCell>> proc;
    public Map<Station, Map<JobClass, Matrix>> pie;
    public Map<Station, SchedStrategy> sched;
    public Map<Integer, Matrix> inchain;
    public Map<Integer, Matrix> visits;    //The integer represents the chain's ID (inchain)
    public Map<Integer, Matrix> nodevisits; //The integer represents the chain's ID (inchain)
    public Map<Station, Map<JobClass, DropStrategy>> droprule;    //This represents dropid in LINE
    public Map<Node, NodeParam> nodeparam;
    public Map<Integer, Sync> sync;
    public Map<Integer, GlobalSync> gsync;
    public Map<Station, SerializableFunction<Matrix, Double>> cdscaling;

    // Impatience (customer abandonment) parameters - Reneging
    public Map<Station, Map<JobClass, ProcessType>> impatienceType;
    public Map<Station, Map<JobClass, Matrix>> impatienceMu;
    public Map<Station, Map<JobClass, Matrix>> impatiencePhi;
    public Map<Station, Map<JobClass, MatrixCell>> impatienceProc;
    public Map<Station, Map<JobClass, Matrix>> impatiencePie;
    public Map<Station, Map<JobClass, Integer>> impatiencePhases;

    // Balking parameters (queue-length or wait-time based refusal to join)
    public Map<Station, Map<JobClass, BalkingStrategy>> balkingStrategy;
    public Map<Station, Map<JobClass, List<BalkingThreshold>>> balkingThresholds;

    // Retrial parameters (orbit and retry behavior)
    public Map<Station, Map<JobClass, ProcessType>> retrialType;
    public Map<Station, Map<JobClass, Matrix>> retrialMu;
    public Map<Station, Map<JobClass, Matrix>> retrialPhi;
    public Map<Station, Map<JobClass, MatrixCell>> retrialProc;
    public Map<Station, Map<JobClass, Integer>> retrialMaxAttempts;

    public Matrix refstat;
    public Matrix njobs;
    public Matrix nservers;
    public Matrix connmatrix;
    public Matrix scv;
    public Matrix isstation;
    public Matrix isstateful;
    public Matrix isstatedep;
    public Matrix nodeToStateful;
    public Matrix nodeToStation;
    public Matrix stationToNode;
    public Matrix stationToStateful;
    public Matrix statefulToStation;
    public Matrix statefulToNode;
    public Matrix rates;
    public Matrix classprio;
    public Matrix classdeadline;
    public Matrix phases;
    public Matrix phasessz;
    public Matrix phaseshift;
    public Matrix schedparam;
    public Matrix chains;
    public Matrix rt;
    public Matrix nvars;
    public Matrix rtnodes;
    public Matrix csmask;
    public Matrix isslc;
    public Matrix immfeed;  // (M x K) matrix where 1.0 indicates immediate feedback enabled for class at station
    public Matrix issignal;  // (nclasses x 1) matrix where 1.0 indicates signal class
    public List<SignalType> signaltype;  // SignalType for each class, null for non-signal classes
    public Matrix syncreply;  // (nclasses x 1) matrix where entry is reply signal class index, -1.0 if no reply expected
    public List<DiscreteDistribution> signalRemovalDist;  // Removal distribution for each signal class (null for single removal)
    public List<RemovalPolicy> signalRemovalPolicy;  // Removal policy for each signal class (null for non-signals)
    public Matrix isCatastrophe;  // (nclasses x 1) matrix where 1.0 indicates catastrophe signal
    public Matrix cap;
    public Matrix classcap;
    public int nregions;  // Number of finite capacity regions (F)
    public MatrixCell region;  // CellMatrix of size F; region.get(f) is Matrix(M, R+1) where entry (i,r) is max jobs of class r at station i in region f; (i,R) is global max at station i; -1 = infinite
    public Matrix regionrule;  // Matrix(F, R) where entry (f,r) is DropStrategy id for class r in region f
    public Matrix regionweight;  // Matrix(F, R) where entry (f,r) is class weight for class r in region f (default 1.0)
    public Matrix regionsz;  // Matrix(F, R) where entry (f,r) is class size/memory for class r in region f (default 1)
    public Matrix refclass;
    public Matrix lldscaling;
    public Map<Station, Matrix> ljdscaling;  // Linearized joint-dependent scaling per station
    public Map<Station, Matrix> ljdcutoffs;  // Per-class cutoffs per station
    public Map<Station, Map<JobClass, Matrix>> ljcdscaling;  // Per-class joint-dependent scaling
    public Map<Station, Matrix> ljcdcutoffs;  // Per-class cutoffs for LJCD
    public Matrix fj;
    public Matrix varsparam;
    public List<NodeType> nodetype;
    public List<String> classnames;
    public List<String> nodenames;
    SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun;

    // Objects - to be removed
    public List<Station> stations;
    public List<StatefulNode> stateful;
    public List<JobClass> jobclasses;
    public List<Node> nodes;

    // Reward definitions for CTMC reward computation
    // Map from reward name to reward function
    public Map<String, RewardFunction> reward;

    // ==================== Heterogeneous Server Fields ====================

    /**
     * Number of server types per station.
     * Matrix (M x 1) where entry (i) is the count of server types at station i.
     * A value of 0 indicates a homogeneous (traditional) queue.
     */
    public Matrix nservertypes;

    /**
     * Server type names per station.
     * servertypenames.get(station).get(t) returns the name of server type t at that station.
     */
    public Map<Station, List<String>> servertypenames;

    /**
     * Number of servers per server type at each station.
     * serverspertype.get(station) returns a Matrix (nTypes x 1) where entry (t)
     * is the count of servers of type t.
     */
    public Map<Station, Matrix> serverspertype;

    /**
     * Server-class compatibility matrix per station.
     * servercompat.get(station) returns a Matrix (nTypes x K) where entry (t, r)
     * is 1.0 if server type t can serve job class r, 0.0 otherwise.
     */
    public Map<Station, Matrix> servercompat;

    /**
     * Heterogeneous service rates per station, server type, and class.
     * heterorates.get(station).get(serverTypeId).get(classId) returns the service rate.
     */
    public Map<Station, Map<Integer, Map<Integer, Double>>> heterorates;

    /**
     * Heterogeneous service process parameters per station, server type, and class.
     * heteroproc.get(station).get(serverTypeId).get(classId) returns the PH matrices.
     */
    public Map<Station, Map<Integer, Map<Integer, MatrixCell>>> heteroproc;

    /**
     * Heterogeneous service process type per station, server type, and class.
     * heteroprocid.get(station).get(serverTypeId).get(classId) returns the ProcessType.
     */
    public Map<Station, Map<Integer, Map<Integer, ProcessType>>> heteroprocid;

    /**
     * Heterogeneous scheduling policy per station.
     * heteroschedpolicy.get(station) returns the HeteroSchedPolicy for that station.
     */
    public Map<Station, HeteroSchedPolicy> heteroschedpolicy;


    /**
     * Returns the list of stations in the network.
     * 
     * @return List of Station objects
     */
    public List<Station> getStations() {
        return this.stations;
    }

    /**
     * Validates the structural consistency of stations, stateful nodes, and nodes
     * according to MATLAB implementation requirements.
     * 
     * This method ensures:
     * 1. All stations are stateful nodes (stations ⊆ stateful)
     * 2. All stateful nodes are nodes (stateful ⊆ nodes)
     * 3. Hash mapping consistency between node types
     * 4. Count consistency for nstations, nstateful, nnodes
     * 
     * @throws IllegalStateException if structural consistency is violated
     */
    public void validateStructuralConsistency() {
        // Check basic counts
        if (nstations < 0 || nstateful < 0 || nnodes < 0) {
            throw new IllegalStateException("Node counts must be non-negative");
        }
        
        if (nstations > nstateful) {
            throw new IllegalStateException("Number of stations cannot exceed number of stateful nodes");
        }
        
        if (nstateful > nnodes) {
            throw new IllegalStateException("Number of stateful nodes cannot exceed total number of nodes");
        }
        
        // Validate matrix dimensions if matrices exist
        if (isstation != null) {
            if (isstation.getNumRows() != nnodes || isstation.getNumCols() != 1) {
                throw new IllegalStateException("isstation matrix must be nnodes x 1");
            }
            
            // Check that station count matches isstation sum
            double stationSum = isstation.elementSum();
            if (Math.abs(stationSum - nstations) > 1e-10) {
                throw new IllegalStateException("nstations must equal sum of isstation matrix");
            }
        }
        
        if (isstateful != null) {
            if (isstateful.getNumRows() != nnodes || isstateful.getNumCols() != 1) {
                throw new IllegalStateException("isstateful matrix must be nnodes x 1");
            }
            
            // Check that stateful count matches isstateful sum
            double statefulSum = isstateful.elementSum();
            if (Math.abs(statefulSum - nstateful) > 1e-10) {
                throw new IllegalStateException("nstateful must equal sum of isstateful matrix");
            }
        }
        
        // Validate hierarchy: all stations must be stateful
        if (isstation != null && isstateful != null) {
            for (int i = 0; i < nnodes; i++) {
                if (isstation.get(i, 0) > 0 && isstateful.get(i, 0) == 0) {
                    throw new IllegalStateException("All stations must be stateful nodes (violation at node " + i + ")");
                }
            }
        }
        
        // Validate node type lists consistency
        if (stations != null && stations.size() != nstations) {
            throw new IllegalStateException("stations list size must match nstations");
        }
        
        if (stateful != null && stateful.size() != nstateful) {
            throw new IllegalStateException("stateful list size must match nstateful");
        }
        
        if (nodes != null && nodes.size() != nnodes) {
            throw new IllegalStateException("nodes list size must match nnodes");
        }
        
        // Validate hash mappings dimensions
        validateHashMappings();
    }
    
    /**
     * Validates the hash mapping matrices used for node type conversions.
     * These mappings must be consistent with MATLAB implementation.
     */
    private void validateHashMappings() {
        if (nodeToStateful != null) {
            if (nodeToStateful.getNumRows() != nnodes || nodeToStateful.getNumCols() != 1) {
                throw new IllegalStateException("nodeToStateful must be nnodes x 1");
            }
        }
        
        if (nodeToStation != null) {
            if (nodeToStation.getNumRows() != nnodes || nodeToStation.getNumCols() != 1) {
                throw new IllegalStateException("nodeToStation must be nnodes x 1");
            }
        }
        
        if (stationToNode != null) {
            if (stationToNode.getNumRows() != nstations || stationToNode.getNumCols() != 1) {
                throw new IllegalStateException("stationToNode must be nstations x 1");
            }
        }
        
        if (stationToStateful != null) {
            if (stationToStateful.getNumRows() != nstations || stationToStateful.getNumCols() != 1) {
                throw new IllegalStateException("stationToStateful must be nstations x 1");
            }
        }
        
        if (statefulToStation != null) {
            if (statefulToStation.getNumRows() != nstateful || statefulToStation.getNumCols() != 1) {
                throw new IllegalStateException("statefulToStation must be nstateful x 1");
            }
        }
        
        if (statefulToNode != null) {
            if (statefulToNode.getNumRows() != nstateful || statefulToNode.getNumCols() != 1) {
                throw new IllegalStateException("statefulToNode must be nstateful x 1");
            }
        }
        
        // Validate hash mapping consistency
        if (nodeToStation != null && stationToNode != null && isstation != null) {
            for (int i = 0; i < nnodes; i++) {
                if (isstation.get(i, 0) > 0) {
                    int stationIdx = (int) nodeToStation.get(i, 0);
                    if (stationIdx < 0 || stationIdx >= nstations) {
                        throw new IllegalStateException("Invalid station index in nodeToStation mapping");
                    }
                    
                    int nodeIdx = (int) stationToNode.get(stationIdx, 0);
                    if (nodeIdx != i) {
                        throw new IllegalStateException("Inconsistent station-node mapping");
                    }
                }
            }
        }
    }
    
    /**
     * Checks if this NetworkStruct is consistent with MATLAB implementation.
     * This is a convenience method that calls validateStructuralConsistency
     * and returns true if no exceptions are thrown.
     * 
     * @return true if structure is consistent, false otherwise
     */
    public boolean isConsistentWithMatlab() {
        try {
            validateStructuralConsistency();
            return true;
        } catch (IllegalStateException e) {
            return false;
        }
    }
    
    /**
     * Print comprehensive information about this NetworkStruct.
     * This method displays all fields, matrices, lists, and maps in a formatted manner
     * useful for debugging and inspection.
     */
    public void print() {
        snPrint(this);
    }
    

}
