package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.constant.ReplacementStrategy;
import jline.solvers.ctmc.CTMC;
import jline.solvers.ssa.SSA;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.fluid.FLD;
import jline.solvers.NetworkAvgNodeTable;

public class CacheExamples {
    
    public static void cache_replc_rr() {
        //System.out.println("=== Cache RR Example ===");
        Network model = CacheModel.cache_replc_rr();
        
        try {
            CTMC solver1 = new CTMC(model, "keep", false, "cutoff", 1, "seed", 1);
            NetworkAvgNodeTable avgTable1 = solver1.getAvgNodeTable();
            //System.out.println("--- CTMC Solver ---");
            avgTable1.print();
            
            model.reset();
            SSA solver2 = new SSA(model, "samples", 10000, "verbose", true, "method", "serial", "seed", 1);
            NetworkAvgNodeTable avgTable2 = solver2.getAvgNodeTable();
            //System.out.println("--- SSA Solver ---");
            avgTable2.print();
            
            model.reset();
            MVA solver3 = new MVA(model, "seed", 1);
            NetworkAvgNodeTable avgTable3 = solver3.getAvgNodeTable();
            //System.out.println("--- MVA Solver ---");
            avgTable3.print();
            
            model.reset();
            NC solver4 = new NC(model, "seed", 1); 
            NetworkAvgNodeTable avgTable4 = solver4.getAvgNodeTable();
            //System.out.println("--- NC Solver ---");
            avgTable4.print();

            model.reset();
            // The reference runs the RMF fluid limit here, which is what its golden
            // holds; the default fluid closure answers a different question.
            FLD solverFLD = new FLD(model, "method", "rmf");
            solverFLD.getAvgNodeTable().print();
            
            jline.lang.nodes.Cache cacheNode = (jline.lang.nodes.Cache) model.getNodeByName("Cache");
            double hitRatio = cacheNode.getHitRatio().get(0);
            double missRatio = cacheNode.getMissRatio().get(0);
            //System.out.println("Hit Ratio: " + hitRatio);
            //System.out.println("Miss Ratio: " + missRatio);
            
        } catch (Exception e) {
            //System.out.println("Error in cache_replc_rr: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    /**
     * LRU cache replacement (cache_replc_lru.ipynb).
     *
     * <p>The reference's three solvers, at its own run length and seed: the SSA
     * row of the golden is that measurement, not a longer or shorter one.
     */
    public static void cache_replc_lru() {
        Network model = CacheModel.cache_replc_lru();

        try {
            new CTMC(model, "keep", false).getAvgNodeTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("CTMC failed: " + e.getMessage());
        }
        try {
            new SSA(model, "samples", 100000, "verbose", true, "method", "serial", "seed", 23000)
                    .getAvgNodeTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("SSA failed: " + e.getMessage());
        }
        try {
            new MVA(model).getAvgNodeTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("MVA failed: " + e.getMessage());
        }
    }

    public static void cache_replc_fifo() {
        //System.out.println("=== Cache FIFO Example ===");
        Network model = CacheModel.cache_replc_fifo();
        
        try {
            CTMC solver1 = new CTMC(model, "keep", false, "seed", 1);
            NetworkAvgNodeTable avgTable1 = solver1.getAvgNodeTable();
            //System.out.println("--- CTMC Solver ---");
            avgTable1.print();
            
            model.reset();
            SSA solver2 = new SSA(model, "samples", 100000, "verbose", true, "method", "serial", "seed", 23000);
            NetworkAvgNodeTable avgTable2 = solver2.getAvgNodeTable();
            //System.out.println("--- SSA Solver ---");
            avgTable2.print();
            
            model.reset();
            MVA solver3 = new MVA(model, "seed", 1);
            NetworkAvgNodeTable avgTable3 = solver3.getAvgNodeTable();
            //System.out.println("--- MVA Solver ---");
            avgTable3.print();

            model.reset();
            // The reference runs the RMF fluid limit here, which is what its golden
            // holds; the default fluid closure answers a different question.
            FLD solverFLD = new FLD(model, "method", "rmf");
            solverFLD.getAvgNodeTable().print();
            
        } catch (Exception e) {
            //System.out.println("Error in cache_replc_fifo: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    public static void cache_replc_routing() {
        //System.out.println("=== Cache Routing Example ===");
        Network model = CacheModel.cache_replc_routing();
        
        try {
            CTMC solver1 = new CTMC(model, "keep", false, "cutoff", 1, "seed", 1);
            NetworkAvgNodeTable avgTable1 = solver1.getAvgNodeTable();
            //System.out.println("--- CTMC Solver ---");
            avgTable1.print();
            
            model.reset();
            SSA solver2 = new SSA(model, "samples", 10000, "verbose", true, "method", "serial", "seed", 1);
            NetworkAvgNodeTable avgTable2 = solver2.getAvgNodeTable();
            //System.out.println("--- SSA Solver ---");
            avgTable2.print();
            
            model.reset();
            MVA solver3 = new MVA(model, "seed", 1);
            NetworkAvgNodeTable avgTable3 = solver3.getAvgNodeTable();
            //System.out.println("--- MVA Solver ---");
            avgTable3.print();
            
            model.reset();
            NC solver4 = new NC(model, "seed", 1);
            NetworkAvgNodeTable avgTable4 = solver4.getAvgNodeTable();
            //System.out.println("--- NC Solver ---");
            avgTable4.print();

            model.reset();
            // The reference runs the RMF fluid limit here, which is what its golden
            // holds; the default fluid closure answers a different question.
            FLD solverFLD = new FLD(model, "method", "rmf");
            solverFLD.getAvgNodeTable().print();
            
            jline.lang.nodes.Cache cacheNode = (jline.lang.nodes.Cache) model.getNodeByName("Cache");
            double hitRatio = cacheNode.getHitRatio().get(0);
            double missRatio = cacheNode.getMissRatio().get(0);
            //System.out.println("Hit Ratio: " + hitRatio);
            //System.out.println("Miss Ratio: " + missRatio);
            
        } catch (Exception e) {
            //System.out.println("Error in cache_replc_routing: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    public static void cache_compare_replc() {
        System.out.println("=== Cache Compare Replacement Strategies ===");

        ReplacementStrategy[] replStrat = {ReplacementStrategy.RR, ReplacementStrategy.FIFO, ReplacementStrategy.LRU};

        for (int s = 0; s < replStrat.length; s++) {

            try {
                // Create a fresh model for each strategy (replacement strategy is set in constructor)
                Network model;
                if (replStrat[s] == ReplacementStrategy.RR) {
                    model = CacheModel.cache_replc_rr();
                } else if (replStrat[s] == ReplacementStrategy.FIFO) {
                    model = CacheModel.cache_replc_fifo();
                } else if (replStrat[s] == ReplacementStrategy.LRU) {
                    model = CacheModel.cache_replc_lru();
                } else {
                    model = CacheModel.cache_compare_replc();
                }
                jline.lang.nodes.Cache cacheNode = (jline.lang.nodes.Cache) model.getNodeByName("Cache");
                
                // Skip CTMC solver as commented in the source notebook
                
                model.reset();
                MVA solver2 = new MVA(model, "seed", 1, "verbose", false);
                NetworkAvgNodeTable avgTable2 = solver2.getAvgNodeTable();
                double mvaHitRatio = cacheNode.getHitRatio().get(0);
                
                double ncHitRatio = Double.NaN;
                try {
                    model.reset();
                    NC solver3 = new NC(model, "seed", 1, "verbose", false);
                    NetworkAvgNodeTable avgTable3 = solver3.getAvgNodeTable();
                    ncHitRatio = cacheNode.getHitRatio().get(0);
                } catch (Exception ncException) {
                    ncHitRatio = Double.NaN;
                }
                
                System.out.printf("%s: MVA=%.8f, NC=%.8f%n", replStrat[s], mvaHitRatio, ncHitRatio);

            } catch (Exception e) {
                System.out.println("Error testing " + replStrat[s] + ": " + e.getMessage());
            }
        }
    }
    
    public static void main(String[] args) {
        //System.out.println("Running all cache examples...\n");
        
        //cache_replc_rr();
        //System.out.println();
        
        //cache_replc_fifo();
        //System.out.println();
        
        cache_replc_routing();
        //System.out.println();
        
        cache_compare_replc();
        //System.out.println();
        
        //System.out.println("All cache examples completed.");
    }
}