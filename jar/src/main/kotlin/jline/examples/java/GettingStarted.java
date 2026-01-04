/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java;

import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;
import jline.VerboseLevel;
import jline.examples.java.advanced.RandomEnvExamples;
import jline.examples.java.basic.OpenExamples;
import jline.examples.java.basic.OpenModel;
import jline.lang.*;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.io.Ret.DistributionResult;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.qns.QNS;
import jline.solvers.ssa.SSA;
import jline.solvers.env.ENV;
import jline.lang.constant.SolverType;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import java.util.function.Function;
import java.util.List;

/**
 * Getting started examples
 */
public class GettingStarted {

    /**
     * Getting started example 1: Basic M/M/1 queue.
     * <p>
     * Features:
     * - Simple M/M/1 queueing system
     * - Exponential arrivals (rate 1.0) and service (rate 2.0)
     * - FCFS scheduling with single server
     * - JMT solver demonstration with result printing
     * - Shows basic network construction and solution
     *
     * @return configured M/M/1 network model
     */
    public static Network tut01_mm1_basics() {
        Network model = new Network("M/M/1");
        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        // Block 2: classes
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(1.0)); // (source, jobclass)
        queue.setService(jobclass, new Exp(2.0)); // (queue, jobclass)
        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));
        // Block 4: solution
        NetworkAvgTable avgTable = new JMT(model, "seed", 23000).getAvgTable();
        avgTable.print();
        //model.jsimgView();
        avgTable.tget(queue, jobclass).print();
        avgTable.tget("Queue", "Class1").print();

        return model;
    }


    /**
     * Getting started example 2: M/G/1 queue with two classes.
     * <p>
     * Features:
     * - Two open classes with different service distributions
     * - Class1: Erlang service distribution (SCV = 1/3)
     * - Class2: Trace-driven service using APH fitting
     * - Demonstrates empirical data integration via Replayer
     * - Manual routing matrix construction
     *
     * @return configured M/G/1 multi-class network model
     */
    public static Network tut02_mg1_multiclass_solvers() {
        Network model = new Network("M/G/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass1 = new OpenClass(model, "Class1");
        OpenClass jobclass2 = new OpenClass(model, "Class2");
        source.setArrival(jobclass1, new Exp(0.5));
        source.setArrival(jobclass2, new Exp(0.5));
        queue.setService(jobclass1, Erlang.fitMeanAndSCV(1, (double) 1 / 3));
        try {
            URI fileURI = GettingStarted.class.getResource("/example_trace.txt").toURI();
            String fileName = Paths.get(fileURI).toString();
            queue.setService(jobclass2, new Replayer(fileName).fitAPH());
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass1, Network.serialRouting(source, queue, sink));
        P.set(jobclass2, Network.serialRouting(source, queue, sink));
        model.link(P);
        NetworkAvgTable avgTable = new JMT(model, "seed", 23000).getAvgTable();
        avgTable.print();

        return model;
    }

    /**
     * Getting started example 3: Machine Repair Problem (MRP).
     * <p>
     * Features:
     * - Closed network with 3 machines
     * - Working state (delay) and repair queue with 2 servers
     * - CTMC solver for exact state space analysis
     * - Steady state probability computation
     * - Demonstrates closed system modeling
     *
     * @return configured machine repair model
     */
    public static Network tut03_repairmen() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "WorkingState");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Exp(4.0));
        model.link(Network.serialRouting(delay, queue));
        CTMC solver = new CTMC(model);
        NetworkAvgTable ctmcAvgTable = solver.getAvgTable();
        ctmcAvgTable.print();
        Matrix stateSpace = solver.getStateSpace().stateSpace;
        stateSpace.print();

        Matrix infGen = solver.getGenerator().infGen;
        infGen.print();
        CTMC.printInfGen(infGen, stateSpace);
        return model;
    }

    /**
     * Getting started example 4: Load balancing with routing strategies.
     * <p>
     * Features:
     * - Load balancer (router) distributing jobs to two PS queues
     * - Demonstrates RAND (random) and RROBIN (round-robin) routing
     * - Identical service rates for fair comparison
     * - Shows how to change routing strategies dynamically
     * - Model reset and re-solve for strategy comparison
     *
     * @return configured load balancing network model
     */
    public static Network tut04_lb_routing() {
        Network model = new Network("RRLB");

        Source source = new Source(model, "Source");
        Router lb = new Router(model, "LB");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1.0));
        queue1.setService(oclass, new Exp(2.0));
        queue2.setService(oclass, new Exp(2.0));

        model.addLink(source, lb);
        model.addLink(lb, queue1);
        model.addLink(lb, queue2);
        model.addLink(queue1, sink);
        model.addLink(queue2, sink);

        lb.setRouting(oclass, RoutingStrategy.RAND);
        new JMT(model, "seed", 23000).getAvgTable().print();

        lb.setRouting(oclass, RoutingStrategy.RROBIN);
        model.reset();
        new JMT(model, "seed", 23000).getAvgTable().print();

        return model;
    }

    /**
     * Getting started example 5: Class switching in single-node network.
     * <p>
     * Features:
     * - Three closed classes cycling through Class1→Class2→Class3→Class1
     * - Single FCFS queue with different Erlang service for each class
     * - Demonstrates class switching routing matrix
     * - Shows impact of setCompletes(false) on performance metrics
     * - NC solver for exact normalizing constant computation
     *
     * @return configured class switching model
     */
    public static Network tut05_completes_flag() {
        Network model = new Network("RL");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass jobClass1 = new ClosedClass(model, "Class1", 1, queue);
        ClosedClass jobClass2 = new ClosedClass(model, "Class2", 0, queue);
        ClosedClass jobClass3 = new ClosedClass(model, "Class3", 0, queue);
        queue.setService(jobClass1, Erlang.fitMeanAndOrder(1, 2));
        queue.setService(jobClass2, Erlang.fitMeanAndOrder(2, 2));
        queue.setService(jobClass3, Erlang.fitMeanAndOrder(3, 2));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass1, jobClass2, queue, queue, 1.0);
        P.set(jobClass2, jobClass3, queue, queue, 1.0);
        P.set(jobClass3, jobClass1, queue, queue, 1.0);
        model.link(P);
        new NC(model).getAvgTable().print();
        new NC(model).getAvgSysTable().print();
        jobClass1.setCompletes(false);
        jobClass2.setCompletes(false);
        new NC(model).getAvgSysTable().print();
        return model;
    }

    /**
     * Getting started example 6: Cache modeling with Zipf access pattern.
     * <p>
     * Features:
     * - LRU cache with 1000 items and capacity 50
     * - Zipf distribution (α=1.4) for realistic access patterns
     * - Three classes: ClientClass, HitClass, MissClass
     * - Different service times for cache hits vs misses
     * - Demonstrates cache performance analysis
     *
     * @return configured cache network model
     */
    public static Network tut06_cache_lru_zipf() {
        Network model = new Network("Model");

        // Block 1: nodes
        Delay clientDelay = new Delay(model, "Client");
        Cache cacheNode = new Cache(model, "Cache", 10, 5, ReplacementStrategy.LRU);
        Delay cacheDelay = new Delay(model, "CacheDelay");

        // Block 2: classes
        ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

        clientDelay.setService(clientClass, Immediate.getInstance()); // (Client,ClientClass)
        cacheDelay.setService(hitClass, Exp.fitMean(0.200)); // (CacheDelay,HitClass)
        cacheDelay.setService(missClass, Exp.fitMean(1.00)); // (CacheDelay,MissClass)

        cacheNode.setRead(clientClass, new Zipf(1.4, 10));
        cacheNode.setHitClass(clientClass, hitClass);
        cacheNode.setMissClass(clientClass, missClass);

        // Block 3: topology
        RoutingMatrix P = model.initRoutingMatrix();

        P.set(clientClass, clientClass, clientDelay, cacheNode, 1.00); // (Client,ClientClass) -> (Cache,ClientClass)

        P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.00); // (Client,HitClass) -> (Cache,HitClass)
        P.set(missClass, missClass, cacheNode, cacheDelay, 1.00); // (Cache,MissClass) -> (CacheDelay,MissClass)

        P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.00); // (Cache,HitClass) -> (CacheDelay,HitClass)
        P.set(missClass, clientClass, cacheDelay, clientDelay, 1.00); // (Client,MissClass) -> (Cache,MissClass)

        model.link(P);
        //new SSA(model,"samples",20000,"seed",1,"verbose", VerboseLevel.STD).getAvgTable().print();
        return model;
    }

    /**
     * Getting started example 7: Simple closed network for optimization.
     * <p>
     * Features:
     * - Basic closed network with delay and PS queue
     * - Single closed class with 5 jobs
     * - Demonstrates foundation for optimization studies
     * - Shows basic structure for parameter analysis
     *
     * @return configured simple closed network model
     */
    public static Network tut07_respt_cdf() {
        Network model = new Network("Model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 5, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.00)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(2.00)); // (Queue1,Class1)

        // Block 3: topology
        model.link(Network.serialRouting(node1, node2));

        // Block 4: Response time distribution analysis
        // Demonstrate fluid and simulation-based CDF computation
        try {
            FLD fluidSolver = new FLD(model);
            DistributionResult RDfluid = fluidSolver.getCdfRespT();
            
            // Store fluid CDF results if available
            if (RDfluid != null) {
                System.out.println("Fluid CDF analysis completed successfully");
            }
        } catch (Exception e) {
            // CDF analysis is optional for this tutorial
            System.out.println("Note: Fluid CDF analysis skipped - " + e.getMessage());
        }
        
        try {
            SolverOptions jmtOptions = new SolverOptions();
            jmtOptions.seed(23000);
            jmtOptions.samples(10000); // Use 1e4 samples
            JMT jmtSolver = new JMT(model, jmtOptions);
            DistributionResult RDsim = jmtSolver.getCdfRespT();
            
            // Store simulation CDF results if available
            if (RDsim != null) {
                System.out.println("JMT CDF analysis completed successfully");
            }
        } catch (Exception e) {
            // CDF analysis is optional for this tutorial
            System.out.println("Note: JMT CDF analysis skipped - " + e.getMessage());
        }

        return model;
    }

    /**
     * Getting started example 8: Optimization with COBYLA algorithm.
     * <p>
     * Features:
     * - Closed load balancing network with optimization
     * - Two PS queues with different service rates
     * - COBYLA optimizer to find optimal routing probability
     * - Minimizes average system response time
     * - Demonstrates integration of queueing models with optimization
     *
     * @return optimized load balancing network model
     */
    public static Network tut08_opt_load_balancing() {
        Network model = new Network("LoadBalCQN");

        // Block 1: nodes
        Delay delay = new Delay(model, "Think");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass cclass = new ClosedClass(model, "Job1", 16, delay);
        delay.setService(cclass, new Exp(1));
        queue1.setService(cclass, new Exp(0.75));
        queue2.setService(cclass, new Exp(0.50));

        // Block 3: topology
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cclass, queue1, delay, 1.0);
        P.set(cclass, queue2, delay, 1.0);
        model.link(P);

        Function<Double, Double> routingModel = (p) -> {
            P.set(cclass, delay, queue1, p);
            P.set(cclass, delay, queue2, 1 - p);
            model.relink(P);
            // Block 4: solution
            MVA solver = new MVA(model, "exact");
            return (Double) solver.getAvgSysRespT().get(0);
        };

        Calcfc objFun = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                return routingModel.apply(x[0]);
            }
        };

        double[] p = {0.5};
        Cobyla.findMinimum(objFun, 1, 0, p, 0.5, 1.0e-8, 0, 10000);

        System.out.println("Optimal p: " + p[0]);
        return model;
    }

    /**
     * Helper method to create M/E/1 model (equivalent to gallery_merl1).
     * <p>
     * Features:
     * - M/E/1 queueing system with exponential arrivals and Erlang service
     * - Arrival rate: 1.0, Service: Erlang with mean 0.5 and order 2
     * - Basic serial routing topology
     * 
     * @return array containing [model, source, queue, sink, oclass]
     */
    public static Object[] gallery_merl1() {
        Network model = new Network("M/E/1");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1.0));
        queue.setService(oclass, Erlang.fitMeanAndOrder(0.5, 2));
        
        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));
        
        return new Object[]{model, source, queue, sink, oclass};
    }

    /**
     * Getting started example 9: Studying a departure process.
     * <p>
     * Features:
     * - Analysis of departure process from M/E/1 queue
     * - CTMC solver with system sampling (sampleSysAggr)
     * - Departure event filtering and inter-departure time analysis
     * - Comparison of empirical vs theoretical SCV (Marshall's formula)
     * - Demonstrates advanced sampling and stochastic analysis
     * 
     * The Java implementation fully supports departure process analysis including:
     * - Filtering departure events from sampled trajectories  
     * - Computing inter-departure times and empirical SCV
     * - Calculating theoretical SCV using Marshall's formula
     * - Comparing empirical vs theoretical results
     */
    public static void tut09_dep_process_analysis() {
        // Create M/E/1 model using gallery_merl1 equivalent
        Object[] modelComponents = gallery_merl1();
        Network model = (Network) modelComponents[0];
        Source source = (Source) modelComponents[1];
        Queue queue = (Queue) modelComponents[2];
        Sink sink = (Sink) modelComponents[3];
        OpenClass oclass = (OpenClass) modelComponents[4];
        
        // Block 4: solution
        CTMC solver = new CTMC(model, "cutoff", 150, "seed", 23000);
        
        // Note: Using sampleSys instead of sampleSysAggr as the latter is not fully implemented
        // In MATLAB: sa = solver.sampleSysAggr(1e5);
        jline.io.Ret.SampleResult sa = solver.sampleSys((int) 5e3);
        
        int queueIndex = model.getNodeIndex(queue);
        
        System.out.println("M/E/1 departure process analysis:");
        System.out.println("Queue index: " + queueIndex);
        System.out.println("Number of samples: " + sa.numSamples);
        System.out.println("Sampling completed using: " + sa.handle);
        
        // Get basic performance metrics for context
        NetworkAvgTable avgTable = solver.getAvgTable();
        avgTable.print();
        
        // Calculate utilization and other metrics for Marshall's formula
        Matrix util = solver.getAvgUtil();
        double queueUtil = util.get(queueIndex, 0);
        System.out.println("Queue utilization: " + queueUtil);
        
        // Get service parameters for theoretical SCV calculation
        double svcRate = queue.getServiceProcess(oclass).getRate();
        double svcSCV = queue.getServiceProcess(oclass).getSCV();
        double arrSCV = source.getArrivalProcess(oclass).getSCV();
        
        System.out.println("Service rate: " + svcRate);
        System.out.println("Service SCV: " + svcSCV);
        System.out.println("Arrival SCV: " + arrSCV);
        
        // Departure event filtering and analysis
        System.out.println("\n=== Departure Process Analysis ===");
        
        
        // 1. Filter departure events from the sample result
        Matrix departureEvents = filterDepartureEvents(sa, queueIndex);
        System.out.println("Number of departure events: " + departureEvents.getNumRows());
        
        if (departureEvents.getNumRows() > 1) {
            // 2. Calculate inter-departure times
            double[] interDepartureTimes = calculateInterDepartureTimes(departureEvents);
            
            // 3. Compute empirical SCV of departures
            double empiricalSCV = calculateEmpiricalSCV(interDepartureTimes);
            System.out.println("Empirical departure SCV: " + empiricalSCV);
            
            // 4. Compare with Marshall's theoretical formula
            double theoreticalSCV = calculateMarshallSCV(queueUtil, svcSCV, arrSCV, solver, queue, oclass);
            System.out.println("Theoretical departure SCV (Marshall): " + theoreticalSCV);
            System.out.println("SCV difference (empirical - theoretical): " + (empiricalSCV - theoreticalSCV));
            
            // Statistical summary
            double meanInterDepartureTime = calculateMean(interDepartureTimes);
            System.out.println("Mean inter-departure time: " + meanInterDepartureTime);
            System.out.println("Expected inter-departure time (1/throughput): " + (1.0 / svcRate * (1 - queueUtil)));
        } else {
            System.out.println("Insufficient departure events for analysis");
        }
    }

    /**
     * Filters departure events from a sample result for a specific node.
     * 
     * @param sampleResult the sample result containing event data
     * @param nodeIndex the index of the node to filter departure events for
     * @return matrix containing [time, nodeIndex] pairs for departure events
     */
    private static Matrix filterDepartureEvents(jline.io.Ret.SampleResult sampleResult, int nodeIndex) {
        Matrix events = sampleResult.event;
        Matrix times = sampleResult.t;
        
        if (events == null || times == null || events.isEmpty() || times.isEmpty()) {
            return new Matrix(0, 2);
        }
        
        // Find departure events for the specified node
        // Event matrix format: [eventType, nodeIndex, classIndex, ...]
        java.util.List<Double> depTimes = new java.util.ArrayList<>();
        java.util.List<Integer> depNodeIndices = new java.util.ArrayList<>();
        
        // Debug: Check event type enum values
        // EventType enum: INIT=0, LOCAL=1, ARV=2, DEP=3, PHASE=4, etc.
        int depEventType = jline.lang.constant.EventType.DEP.ordinal();
        
        for (int i = 0; i < events.getNumRows(); i++) {
            // Event matrix format: [eventType, nodeIndex, classIndex, ...]
            int eventType = (int) events.get(i, 0);
            int eventNodeIndex = (int) events.get(i, 1);
            
            if (eventType == depEventType && eventNodeIndex == nodeIndex) { // DEP event for our queue
                depTimes.add(times.get(i, 0));
                depNodeIndices.add(nodeIndex);
            }
        }
        
        // Convert to matrix
        Matrix departureEvents = new Matrix(depTimes.size(), 2);
        for (int i = 0; i < depTimes.size(); i++) {
            departureEvents.set(i, 0, depTimes.get(i));
            departureEvents.set(i, 1, depNodeIndices.get(i));
        }
        
        return departureEvents;
    }
    
    /**
     * Calculates inter-departure times from a matrix of departure events.
     * 
     * @param departureEvents matrix with departure times in first column
     * @return array of inter-departure times
     */
    private static double[] calculateInterDepartureTimes(Matrix departureEvents) {
        if (departureEvents.getNumRows() <= 1) {
            return new double[0];
        }
        
        double[] interTimes = new double[departureEvents.getNumRows() - 1];
        for (int i = 1; i < departureEvents.getNumRows(); i++) {
            interTimes[i - 1] = departureEvents.get(i, 0) - departureEvents.get(i - 1, 0);
        }
        
        return interTimes;
    }
    
    /**
     * Calculates the empirical squared coefficient of variation (SCV) from a sample.
     * 
     * @param samples array of sample values
     * @return empirical SCV (variance / mean^2)
     */
    private static double calculateEmpiricalSCV(double[] samples) {
        if (samples.length <= 1) {
            return 0.0;
        }
        
        double mean = calculateMean(samples);
        double variance = calculateVariance(samples, mean);
        
        return variance / (mean * mean);
    }
    
    /**
     * Calculates the mean of a sample array.
     * 
     * @param samples array of sample values
     * @return sample mean
     */
    private static double calculateMean(double[] samples) {
        if (samples.length == 0) {
            return 0.0;
        }
        
        double sum = 0.0;
        for (double sample : samples) {
            sum += sample;
        }
        
        return sum / samples.length;
    }
    
    /**
     * Calculates the variance of a sample array.
     * 
     * @param samples array of sample values
     * @param mean the sample mean
     * @return sample variance
     */
    private static double calculateVariance(double[] samples, double mean) {
        if (samples.length <= 1) {
            return 0.0;
        }
        
        double sumSquaredDiffs = 0.0;
        for (double sample : samples) {
            double diff = sample - mean;
            sumSquaredDiffs += diff * diff;
        }
        
        return sumSquaredDiffs / (samples.length - 1);
    }
    
    /**
     * Calculates the theoretical departure SCV using Marshall's formula.
     * Marshall's formula: SCVd = SCVa + 2*ρ²*SCVs - 2*ρ*(1-ρ)*μ*W
     * 
     * @param utilization queue utilization (ρ)
     * @param serviceSCV service process SCV
     * @param arrivalSCV arrival process SCV  
     * @param solver the CTMC solver for getting wait times
     * @param queue the queue node
     * @param jobClass the job class
     * @return theoretical departure SCV
     */
    private static double calculateMarshallSCV(double utilization, double serviceSCV, double arrivalSCV,
                                             CTMC solver, Queue queue, OpenClass jobClass) {
        try {
            // Get average wait time from solver using tget method
            NetworkAvgTable avgTable = solver.getAvgTable();
            NetworkAvgTable queueTable = avgTable.tget(queue, jobClass);
            
            // Get wait time from response time and service time
            List<Double> respTimes = queueTable.getRespT();
            double avgRespTime = respTimes.get(0);
            double avgServiceTime = 1.0 / queue.getServiceProcess(jobClass).getRate();
            double avgWaitTime = avgRespTime - avgServiceTime;
            
            // Get service rate
            double serviceRate = queue.getServiceProcess(jobClass).getRate();
            
            // Marshall's formula
            double rho = utilization;
            double scvDeparture = arrivalSCV + 2 * rho * rho * serviceSCV - 2 * rho * (1 - rho) * serviceRate * avgWaitTime;
            
            return scvDeparture;
        } catch (Exception e) {
            System.out.println("Warning: Could not calculate Marshall's formula: " + e.getMessage());
            // Fallback approximation when wait time is unavailable
            double rho = utilization;
            return arrivalSCV + 2 * rho * rho * serviceSCV;
        }
    }

    /**
     * Getting started example 10: Basic layered queueing network.
     * <p>
     * Features:
     * - Layered queueing network (LQN) modeling
     * - Two-tier client-server application
     * - Client processor and database processor
     * - Reference task with think time
     * - Synchronous calls between tasks
     * - Multiple database calls per client request
     * - Demonstrates hierarchical performance analysis
     *
     * @return configured layered network model
     */
    public static LayeredNetwork tut10_lqn_basics() {
        // Create the layered network model
        LayeredNetwork model = new LayeredNetwork("ClientDBSystem");

        // Create processors
        Processor P1 = new Processor(model, "ClientProcessor", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "DBProcessor", 1, SchedStrategy.PS);

        // Create tasks
        Task T1 = new Task(model, "ClientTask", 10, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Exp.fitMean(5.0)); // 5-second think time
        Task T2 = new Task(model, "DBTask", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);

        // Create entries that represent service interfaces
        Entry E1 = new Entry(model, "ClientEntry").on(T1);
        Entry E2 = new Entry(model, "DBEntry").on(T2);

        // Define activities that specify the work performed and synchronous calls
        // Client activity: processes request and calls DB
        Activity A1 = new Activity(model, "ClientActivity", Exp.fitMean(1.0)).on(T1);
        A1.boundTo(E1).synchCall(E2, 2.5); // 2.5 DB calls on average

        // DB activity: processes database request
        Activity A2 = new Activity(model, "DBActivity", Exp.fitMean(0.8)).on(T2);
        A2.boundTo(E2).repliesTo(E2);

        return model;
    }

    /**
     * Getting started example 11: Random environments and SolverENV.
     * <p>
     * This tutorial has been moved to RandomEnvExamples.example_random_env_basics().
     * This method delegates to that implementation for backward compatibility.
     * <p>
     * This tutorial illustrates how to model a queueing system operating in
     * a random environment, where system parameters (e.g., service rates)
     * change according to an underlying environmental process.
     * <p>
     * Scenario: A server that alternates between "Fast" and "Slow" modes.
     * In Fast mode, service rate is 4.0. In Slow mode, service rate is 1.0.
     * The environment switches from Fast->Slow at rate 0.5 and Slow->Fast at rate 1.0.
     * <p>
     * Features:
     * - Random environment with two operational stages
     * - Stage-dependent service rates
     * - Exponential transitions between environment stages
     * - ENV solver with Fluid (FLD) sub-solver
     * - Environment-averaged performance metrics
     *
     * @return configured environment model
     * @throws Exception if solver encounters an error
     * @see jline.examples.java.advanced.RandomEnvExamples#example_random_env_basics()
     */
    public static Environment tut11_random_env() throws Exception {
        return RandomEnvExamples.example_random_env_basics();
    }

    /**
     * Main method for testing and demonstrating getting started examples.
     *
     * <p>Currently configured to run tut09_dep_process_analysis() and solve it
     * using the SSA solver with specified sampling parameters.
     *
     * @param args command line arguments (not used)
     * @throws IllegalAccessException       if access restrictions prevent execution
     * @throws ParserConfigurationException if XML parsing configuration fails
     */
    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        Network model = GettingStarted.tut04_lb_routing();
    }
}
