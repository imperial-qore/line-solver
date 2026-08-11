/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Data structure representing a finite capacity region with all its constraints and properties.
 * This structure is extracted from Region objects and stored in NetworkStruct.
 */
public class RegionStruct {
    public List<Integer> stationIndices;  // Station indices in this region
    public int globalMaxJobs;              // Global capacity limit for region (-1 = unbounded)
    public int globalMaxMemory;            // Global memory limit for region (-1 = unbounded)
    public Map<Integer, Integer> classMaxJobs;   // [classIdx] = per-class max jobs (-1 = unbounded)
    public Map<Integer, Integer> classMaxMemory; // [classIdx] = per-class max memory (-1 = unbounded)
    public Map<Integer, Boolean> dropRule;      // [classIdx] = drop rule (true=drop, false=buffer)
    public Map<Integer, Integer> classSize;     // [classIdx] = memory size per job

    public RegionStruct() {
        this.stationIndices = new ArrayList<>();
        this.globalMaxJobs = Region.UNBOUNDED;
        this.globalMaxMemory = Region.UNBOUNDED;
        this.classMaxJobs = new HashMap<>();
        this.classMaxMemory = new HashMap<>();
        this.dropRule = new HashMap<>();
        this.classSize = new HashMap<>();
    }
}
