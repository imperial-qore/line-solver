package jline.examples.java.advanced;

import jline.lang.layered.*;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.DiscreteSampler;
import java.util.ArrayList;
import java.util.List;

/**
 * Model definitions for Layered CQ (Contention Queue) examples.
 * Companion class to LayeredCQExamples.
 */
public class LayeredCQModel {

    /**
     * Model for lcq_singlehost.ipynb - matches MATLAB lcq_singlehost.m
     * This model includes a cache task with hit/miss activities
     */
    public static LayeredNetwork lcq_singlehost() {
        LayeredNetwork model = new LayeredNetwork("cacheInLayeredNetwork");

        // Client processor and task (matching MATLAB)
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        // Cache task (matching MATLAB)
        int totalitems = 4;
        int cachecapacity = 2;
        jline.util.matrix.Matrix probabilities = new jline.util.matrix.Matrix(1, totalitems);
        for (int i = 0; i < totalitems; i++) {
            probabilities.set(0, i, 1.0 / totalitems);
        }
        DiscreteSampler pAccess = new DiscreteSampler(probabilities);

        Processor PC = new Processor(model, "PC", 1, SchedStrategy.PS);
        CacheTask C2 = new CacheTask(model, "C2", totalitems, cachecapacity, ReplacementStrategy.RR, 1);
        C2.on(PC);
        ItemEntry I2 = new ItemEntry(model, "I2", totalitems, pAccess).on(C2);

        // Activities (matching MATLAB)
        Activity A1 = new Activity(model, "A1", new Immediate()).on(T1).boundTo(E1).synchCall(I2, 1.0);
        Activity AC2 = new Activity(model, "AC2", new Immediate()).on(C2).boundTo(I2);
        Activity AC2h = new Activity(model, "AC2h", new Exp(1.0)).on(C2).repliesTo(I2);
        Activity AC2m = new Activity(model, "AC2m", new Exp(0.5)).on(C2).repliesTo(I2);

        // Set up cache access precedence (hit/miss branching)
        List<Activity> postActs = new ArrayList<>();
        postActs.add(AC2h);
        postActs.add(AC2m);
        C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, postActs));

        return model;
    }

    /**
     * Model for lcq_threehosts.ipynb - matches MATLAB lcq_threehosts.m
     * This model includes a cache task with three hosts (processors)
     */
    public static LayeredNetwork lcq_threehosts() {
        LayeredNetwork model = new LayeredNetwork("LQNwithCaching");
        int nusers = 1;
        int ntokens = 1;
        
        // Client processor and task (matching MATLAB)
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", nusers, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);
        
        // Cache task with two-level cache (matching MATLAB)
        int totalitems = 4;
        int[] cachecapacity = new int[]{1, 1}; // Two-level cache matching MATLAB's [1,1]
        jline.util.matrix.Matrix probabilities = new jline.util.matrix.Matrix(1, totalitems);
        for (int i = 0; i < totalitems; i++) {
            probabilities.set(0, i, 1.0 / totalitems);
        }
        DiscreteSampler pAccess = new DiscreteSampler(probabilities);
        
        Processor PC = new Processor(model, "Pc", 1, SchedStrategy.PS);
        CacheTask C2 = new CacheTask(model, "CT", totalitems, cachecapacity, ReplacementStrategy.RR, ntokens);
        C2.on(PC);
        ItemEntry I2 = new ItemEntry(model, "IE", totalitems, pAccess).on(C2);
        
        // Third processor and task (matching MATLAB)
        Processor P3 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T3 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P3);
        Entry E3 = new Entry(model, "E2").on(T3);
        Activity A3 = new Activity(model, "A2", new Exp(5.0)).on(T3).boundTo(E3).repliesTo(E3);
        
        // Activities (matching MATLAB)
        Activity A1 = new Activity(model, "A1", new Immediate()).on(T1).boundTo(E1).synchCall(I2, 1.0);
        
        Activity AC2 = new Activity(model, "Ac", new Immediate()).on(C2).boundTo(I2);
        Activity AC2h = new Activity(model, "Ac_hit", new Exp(1.0)).on(C2).repliesTo(I2);
        Activity AC2m = new Activity(model, "Ac_miss", new Exp(0.5)).on(C2).synchCall(E3, 1.0).repliesTo(I2);
        
        // Set up cache access precedence (hit/miss branching)
        List<Activity> postActs = new ArrayList<>();
        postActs.add(AC2h);
        postActs.add(AC2m);
        C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, postActs));
        
        return model;
    }
}