package jline.lang;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.nodes.Station;
import jline.util.Matrix;

/**
 * Set of reachable classes for a given job (a chain)
 */
public class Chain extends NetworkElement implements Serializable{

	protected List<Station> stations;
	protected List<JobClass> classes;
	protected List<String> classnames;
	protected Matrix visits; //visist.get(i,r) means the number of visist that a job in chain c pays to station i in class r
	protected Map<JobClass, Integer> classIndexMap;
	protected Map<Station, Integer> stationIndexMap;
	protected Matrix completes;
	protected Matrix njobs;
	
	public Chain(String neName) {
		super(neName);
		
		stations = new ArrayList<Station>();
		classes = new ArrayList<JobClass>();
		classnames = new ArrayList<String>();
		visits = new Matrix(0,0,0);
		classIndexMap = new HashMap<JobClass, Integer>();
		stationIndexMap = new HashMap<Station, Integer>();
		completes = new Matrix(0,0,0);
		njobs = new Matrix(0,0,0);
	}
	
	public Chain(String neName, List<JobClass> classes, List<Station> stations) {
		this(neName);
		
		int nClasses = classes.size();
		int nStations = stations.size();
		visits.reshape(nStations, nClasses, nStations * nClasses);
		completes.reshape(nStations, nClasses, nStations * nClasses);
		njobs.reshape(nStations, nClasses, nStations * nClasses);
		for(int i = 0; i < nClasses; i++) {
			JobClass jobclass = classes.get(i);
			this.classes.add(jobclass);
			this.classnames.add(jobclass.getName());
			this.classIndexMap.put(jobclass, i);
			this.njobs.set(0, i, jobclass.getNumberOfJobs());
			if (jobclass.completes)
				this.completes.set(0, i, 1);
		}
		for(int i = 0; i < nStations; i++) {
			Station station = stations.get(i);
			this.stations.add(station);
			this.stationIndexMap.put(station, i);
		}
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public void setVisits(JobClass jobclass, Station station, double val) {
		int jobClassIdx = this.classIndexMap.getOrDefault(jobclass, -1);
		int stationIdx = this.stationIndexMap.getOrDefault(station, -1);
		
		if (jobClassIdx == -1 || stationIdx == -1)
			return;
		
		this.visits.set(stationIdx, jobClassIdx, val);
	}
	
	public void addClass(JobClass jobclass) {
		int jobClassIdx = this.classIndexMap.getOrDefault(jobclass, -1);
		
		if(jobClassIdx == -1) {
			int idx = this.classes.size();
			this.visits.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
			this.completes.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
			this.njobs.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
			
			this.classes.add(jobclass);
			this.classnames.add(jobclass.getName());
			this.classIndexMap.put(jobclass, idx);
			this.njobs.set(0, idx, jobclass.getNumberOfJobs());
			if (jobclass.completes)
				this.completes.set(0, idx, 1);
		}
	}
	
	public void addStation(Station station) {
		int stationIdx = this.stationIndexMap.getOrDefault(station, -1);
		
		if(stationIdx == -1) {
			int idx = this.stations.size();
			this.stations.add(station);
			this.stationIndexMap.put(station, idx);
			this.visits.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
			this.completes.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
			this.njobs.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
		}
	}
}
