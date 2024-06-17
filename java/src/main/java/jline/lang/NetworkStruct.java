package jline.lang;

import java.io.Serializable;
import java.util.List;
import java.util.Map;

import jline.lang.constant.DropStrategy;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.state.EventCacheKey;
import jline.lang.state.EventResult;
import jline.util.*;

/**
 * Class summarizing the characteristics of a Network object
 */
public class NetworkStruct implements Serializable, Cloneable {
    //For data structure, {} is represented by HashMap, [] is represented by ArrayList, [][] is represented by matrix;
    //For the matrix that stores Constant. Use double list instead.
    public int nstations;
    public int nstateful;
    public int nnodes;
    public int nclasses;
    public double nclosedjobs;
    public int nchains;

    //TODO: check if use of station, stateful and node are consistent with MATLAB
    public Map<JobClass, Map<JobClass, Matrix>> rtorig;
    public Map<Station, Map<JobClass, SerializableFunction<Double, Double>>> lst;
    public Map<StatefulNode, Matrix> state;
    public Map<StatefulNode, Matrix> stateprior;
    public Map<Station, Matrix> space;
    public Map<Node, Map<JobClass, RoutingStrategy>> routing;
    public Map<Station, Map<JobClass, ProcessType>> proctype;
    public Map<Station, Map<JobClass, Matrix>> mu;
    public Map<Station, Map<JobClass, Matrix>> phi;
    public Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc;
    public Map<Station, Map<JobClass, Matrix>> pie;
    public Map<Station, SchedStrategy> sched;
    public Map<Integer, Matrix> inchain;
    public Map<Integer, Matrix> visits;	//The integer represents the chain's ID (inchain)
    public Map<Integer, Matrix> nodevisits; //The integer represents the chain's ID (inchain)
    public Map<Station, Map<JobClass, DropStrategy>> droprule;	//This represents dropid in LINE
	public Map<Node, NodeParam> nodeparam;
	public Map<Integer, Sync> sync;
	public Map<Station, SerializableFunction<Matrix, Double>> cdscaling;
    
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
    public Matrix statefulToNode;
    public Matrix rates;
    public Matrix classprio;
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
    public Matrix cap;
    public Matrix classcap;
    public Matrix refclass;
    public Matrix lldscaling;
    public Matrix fj;
    
    SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun;
    
    public List<NodeType> nodetypes;
    public List<String> classnames;
    public List<String> nodenames;
    public List<Station> stations;
    public List<StatefulNode> stateful;
    public List<JobClass> jobclasses;
    public List<Node> nodes;
    public Map<EventCacheKey, EventResult> eventCache;
    public boolean eventCacheEnabled = true;

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

}
