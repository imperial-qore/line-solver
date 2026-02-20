/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class modelling a set of reachable classes for a given job (a chain)
 */
public class Chain extends NetworkElement implements Serializable {

    protected List<Station> stations;
    protected List<JobClass> classes;
    protected List<String> classnames;
    protected Matrix visits; //visist.get(i,r) means the number of visist that a job in chain c pays to station i in class r
    protected Map<JobClass, Integer> classIndexMap;
    protected Map<Station, Integer> stationIndexMap;
    protected Matrix completes;
    protected Matrix njobs;

    /**
     * Creates a new chain with the specified name.
     * Initializes empty collections for stations, classes, and visit matrices.
     * 
     * @param neName the name for this chain
     */
    public Chain(String neName) {
        super(neName);

        stations = new ArrayList<Station>();
        classes = new ArrayList<JobClass>();
        classnames = new ArrayList<String>();
        visits = new Matrix(0, 0, 0);
        classIndexMap = new HashMap<JobClass, Integer>();
        stationIndexMap = new HashMap<Station, Integer>();
        completes = new Matrix(0, 0, 0);
        njobs = new Matrix(0, 0, 0);
    }

    /**
     * Creates a new chain with the specified name, job classes, and stations.
     * Initializes visit, completion, and job count matrices based on the provided classes and stations.
     * 
     * @param neName the name for this chain
     * @param classes the list of job classes in this chain
     * @param stations the list of stations visited by this chain
     */
    public Chain(String neName, List<JobClass> classes, List<Station> stations) {
        this(neName);

        int nClasses = classes.size();
        int nStations = stations.size();
        visits.reshape(nStations, nClasses, nStations * nClasses);
        completes.reshape(nStations, nClasses, nStations * nClasses);
        njobs.reshape(nStations, nClasses, nStations * nClasses);
        for (int i = 0; i < nClasses; i++) {
            JobClass jobclass = classes.get(i);
            this.classes.add(jobclass);
            this.classnames.add(jobclass.getName());
            this.classIndexMap.put(jobclass, i);
            this.njobs.set(0, i, jobclass.getNumberOfJobs());
            if (jobclass.completes)
                this.completes.set(0, i, 1);
        }
        for (int i = 0; i < nStations; i++) {
            Station station = stations.get(i);
            this.stations.add(station);
            this.stationIndexMap.put(station, i);
        }
    }

    /**
     * Adds a job class to this chain.
     * Expands the visit, completion, and job count matrices if the class is new.
     * 
     * @param jobclass the job class to add to this chain
     */
    public void addClass(JobClass jobclass) {
        int jobClassIdx = this.classIndexMap.getOrDefault(jobclass, -1);

        if (jobClassIdx == -1) {
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

    /**
     * Adds a station to this chain.
     * Expands the visit, completion, and job count matrices if the station is new.
     * 
     * @param station the station to add to this chain
     */
    public void addStation(Station station) {
        int stationIdx = this.stationIndexMap.getOrDefault(station, -1);

        if (stationIdx == -1) {
            int idx = this.stations.size();
            this.stations.add(station);
            this.stationIndexMap.put(station, idx);
            this.visits.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
            this.completes.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
            this.njobs.expandMatrix(this.stations.size(), this.classes.size(), this.stations.size() * this.classes.size());
        }
    }

    /**
     * Sets the name of this chain.
     * 
     * @param name the new name for this chain
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Sets the number of visits that a job in this chain makes to a station in a specific class.
     * 
     * @param jobclass the job class making the visits
     * @param station the station being visited
     * @param val the number of visits
     */
    public void setVisits(JobClass jobclass, Station station, double val) {
        int jobClassIdx = this.classIndexMap.getOrDefault(jobclass, -1);
        int stationIdx = this.stationIndexMap.getOrDefault(station, -1);

        if (jobClassIdx == -1 || stationIdx == -1)
            return;

        this.visits.set(stationIdx, jobClassIdx, val);
    }

    /**
     * Gets the list of job classes in this chain.
     * 
     * @return the list of job classes
     */
    public List<JobClass> getClasses() {
        return this.classes;
    }
}
