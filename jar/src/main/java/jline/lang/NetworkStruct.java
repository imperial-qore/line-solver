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

import static jline.api.sn.SnPrint.snPrint;

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
    public Map<Integer, FJSync> fjsync; // fork firing synchronizations (FJ tag-augmented structs only, see ModelAdapter.fjtag)
    public Matrix fjclassmap; // (1,nclasses) original class of each FJ auxiliary class, -1 for originals (FJ tag-augmented structs only)
    public boolean isfjaugmented = false; // true on FJ tag-augmented structs (Join/Fork carry count-vector states)
    public Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling;
    // Declared peak (max) class-dependent rate scaling per class: a 1xR row
    // vector per class-dependent station, used to normalize Util = T*S/peak.
    public Map<Station, Matrix> cdscalingpeak;
    // Joint-dependence (non-product-form) scalings eta_i(n): a handle of the
    // joint per-class population vector, scalar-shared or per-class. Twin of
    // cdscaling, kept separate to preserve the product-form vs joint distinction
    // (see setJointDependence). Evaluated numerically identically to cdscaling.
    public Map<Station, SerializableFunction<Matrix, Matrix>> jdscaling;
    // Declared peak joint-dependent rate scaling per class, twin of cdscalingpeak.
    public Map<Station, Matrix> jdscalingpeak;

    // Impatience (customer abandonment) parameters - Reneging
    public Map<Station, Map<JobClass, ProcessType>> impatienceType;
    public Map<Station, Map<JobClass, Matrix>> impatienceMu;
    public Map<Station, Map<JobClass, Matrix>> impatiencePhi;
    public Map<Station, Map<JobClass, MatrixCell>> impatienceProc;
    public Map<Station, Map<JobClass, Matrix>> impatiencePie;
    public Map<Station, Map<JobClass, Integer>> impatiencePhases;
    // ImpatienceType (RENEGING, BALKING) per station-class; absent entry = none.
    // Mirrors MATLAB sn.impatienceClass(i,r) and Python sn.impatienceClass[i,r],
    // which are 0 where no impatience type is set.
    public Map<Station, Map<JobClass, ImpatienceType>> impatienceClass;

    // Balking parameters (queue-length or wait-time based refusal to join)
    public Map<Station, Map<JobClass, BalkingStrategy>> balkingStrategy;
    public Map<Station, Map<JobClass, List<BalkingThreshold>>> balkingThresholds;

    // Retrial parameters (orbit and retry behavior)
    public Map<Station, Map<JobClass, ProcessType>> retrialType;
    public Map<Station, Map<JobClass, Matrix>> retrialMu;
    public Map<Station, Map<JobClass, Matrix>> retrialPhi;
    public Map<Station, Map<JobClass, MatrixCell>> retrialProc;
    public Map<Station, Map<JobClass, Integer>> retrialMaxAttempts;

    /**
     * Retrial policy per station-class: {@link jline.lang.constant.RetrialPolicy#LINEAR}
     * (per-customer timers, aggregate rate n*nu) or
     * {@link jline.lang.constant.RetrialPolicy#CONSTANT} (one controller retries for
     * the whole orbit at rate nu). An absent entry means LINEAR. Mirrors the
     * (nstations x nclasses) sn.retrialPolicy of MATLAB.
     */
    public Map<Station, Map<JobClass, Integer>> retrialPolicy;

    /**
     * Orbit capacity per station-class; -1 (or an absent entry) means unbounded.
     * A job that finds the orbit full is lost. Mirrors sn.orbitMaxJobs of MATLAB.
     */
    public Map<Station, Map<JobClass, Integer>> orbitMaxJobs;

    // Server breakdown / repair. Failure and repair are properties of the SERVER,
    // so they are per station; the optional degraded service used while the server
    // is down is a service distribution and is therefore per class.

    /** (nnodes x 1) 1 iff the node's server is subject to breakdowns. */
    public Matrix hasbreakdown;

    /** (nstations x 1) failure rate of an up server (0 = never fails). */
    public Matrix breakdownMu;

    /** (nstations x 1) repair rate of a down server (0 = never repaired). */
    public Matrix repairMu;

    /** Failure-time (D0,D1) representation per station. */
    public Map<Station, MatrixCell> breakdownProc;

    /** Repair-time (D0,D1) representation per station. */
    public Map<Station, MatrixCell> repairProc;

    /** (nstations x nclasses) service rate while the server is down (0 = no service). */
    public Matrix downServiceRates;

    // Orbit impatience parameters (abandonment from the retrial orbit); (D0,D1) process per station-class
    public Map<Station, Map<JobClass, MatrixCell>> orbitImpatience;

    // Batch rejection probability per station-class (retrial queues); 0 = partial admission allowed
    public Map<Station, Map<JobClass, Double>> batchRejectProb;

    public Matrix refstat;
    public Matrix njobs;
    public Matrix nservers;
    public Matrix connmatrix;
    public Matrix scv;
    public Matrix isstation;
    public Matrix isstateful;
    public Matrix isstatedep;
    /**
     * Station-indexed mask of queue stations carrying setup/delay-off times
     * (function stations), as in MATLAB refreshStruct.m.
     */
    public Matrix isfunction;
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
    /** (nnodes,1) 1 iff the node is the upstream/blocking side of a true-BAS relation (BUG-83). */
    public Matrix isbasblocking;
    /**
     * (nstations,nclasses) 1 iff refusing an arrival of that class here must BLOCK an
     * upstream BAS station rather than lose the job. Complement of {@link #isbasblocking}:
     * a refusal is resolved by {@code State.arrivalIsLost} at the DESTINATION, which can
     * only see that station's own drop rule, but under the upstream declaration form the
     * BAS rule sits on the blocking station instead. Without this field an open-class
     * arrival refused at the destination is declared LOST, the become-blocked edge never
     * fires, and the blocking station degenerates into an isolated M/M/1/K.
     */
    public Matrix isbasdestination;
    public Matrix rtnodes;
    public Matrix csmask;
    public Matrix isslc;
    public Matrix immfeed;  // (M x K) matrix where 1.0 indicates immediate feedback enabled for class at station
    public Matrix issignal;  // (nclasses x 1) matrix where 1.0 indicates signal class
    public Matrix signaltarget;  // (nclasses x 1) 0-based index of the class a signal removes; -1 if none
    public List<SignalType> signaltype;  // SignalType for each class, null for non-signal classes
    public Matrix syncreply;  // (nclasses x 1) matrix where entry is reply signal class index, -1.0 if no reply expected
    public Matrix classspawn; // (nclasses x 1) matrix where entry is the class injected at the same station on each completion, -1.0 if none
    public List<DiscreteDistribution> signalremdist;  // Removal distribution for each signal class (null for single removal)
    public List<DiscreteDistribution> arrivalbatch;  // (nclasses) batch-size law released at each Source arrival epoch; null = single arrival
    public List<RemovalPolicy> signalrempolicy;  // Removal policy for each signal class (null for non-signals)
    public Matrix iscatastrophe;  // (nclasses x 1) matrix where 1.0 indicates catastrophe signal
    public Matrix cap;
    public Matrix classcap;
    public int nregions;  // Number of finite capacity regions (F)
    public MatrixCell region;  // CellMatrix of size F; region.get(f) is Matrix(M, R+1) where entry (i,r) is max jobs of class r at station i in region f; (i,R) is global max at station i; -1 = infinite
    public Matrix regionrule;  // Matrix(F, R) where entry (f,r) is DropStrategy id for class r in region f
    public Matrix regionweight;  // Matrix(F, R) where entry (f,r) is class weight for class r in region f (default 1.0)
    public Matrix regionsz;  // Matrix(F, R) where entry (f,r) is class size/memory for class r in region f (default 1)
    public MatrixCell regionmaxmem;  // CellMatrix(F); get(f) is Matrix(M, 1) where entry (i,0) is the region global memory budget replicated on each member station i, -1 = unbounded
    // CellMatrix(F); get(f) is Matrix(M, 1) with entry (i,0) = 1 iff station i belongs to region f.
    // Membership must be recorded explicitly because it cannot be recovered from region.get(f):
    // -1 there means "unbounded", which is indistinguishable from "not a member", so a region
    // constrained only by regionlincon would read as empty and be silently ignored by every engine.
    public MatrixCell regionmembers;
    public Map<Integer, MatrixCell> regionlincon;  // get(f) is a MatrixCell holding the (A,b) pair for region f: get(0) is Matrix(C_f, K) constraint matrix, get(1) is Matrix(C_f, 1) capacity vector
    public Matrix refclass;
    public Matrix lldscaling;
    public Matrix fj;
    public Matrix varsparam;
    public Matrix markidx; // (nstations x nclasses) mark index (1-based) of class r at source station i, -1 = not marked
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
    // Per-station heterogeneous-server parameters (server-type lists, compatibility
    // matrices, per-type/per-class rates and processes) are ragged, node-type-conditional
    // structures and are therefore stored in the nodeparam container rather than as flat
    // root fields. They live on the station's ServiceNodeParam (see getServiceParam) and
    // are populated only for Queue stations declaring server types.

    /**
     * Returns the {@link ServiceNodeParam} attached to the given station, or null if
     * the station has no nodeparam entry or its entry is not a ServiceNodeParam.
     *
     * <p>This is the single access point for the per-station heterogeneous-server
     * fields (nservertypes, servertypenames, serverspertype, servercompat,
     * heterorates, heteroproc, heteroprocid, heteroschedpolicy).</p>
     *
     * @param st the station
     * @return the station's ServiceNodeParam, or null
     */
    public jline.lang.nodeparam.ServiceNodeParam getServiceParam(Station st) {
        if (nodeparam == null || st == null) {
            return null;
        }
        NodeParam p = nodeparam.get(st);
        if (p instanceof jline.lang.nodeparam.ServiceNodeParam) {
            return (jline.lang.nodeparam.ServiceNodeParam) p;
        }
        return null;
    }


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
