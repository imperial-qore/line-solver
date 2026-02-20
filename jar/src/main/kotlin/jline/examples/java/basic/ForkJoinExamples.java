/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.solvers.jmt.JMT;
import jline.solvers.mva.MVA;
import java.util.Scanner;

/**
 * Fork-join network examples mirroring the Kotlin notebooks in forkJoin.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/forkJoin/. Each method 
 * demonstrates a specific fork-join concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic fork-join synchronization patterns
 * - Asymmetric and nested fork-join structures
 * - Class switching within fork-join topologies
 * - Complex serial and parallel compositions
 * - Deep nesting and routing overlaps
 */
public class ForkJoinExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        // Skip pause if running in non-interactive mode (e.g., Maven exec)
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors in case of pipe or redirection
        }
    }

    /**
     * Asymmetric fork-join network (fj_asymm.ipynb).
     * <p>
     * Demonstrates an asymmetric fork-join where one branch has serial queues
     * (Queue2→Queue3) while the other has a single queue (Queue1).
     * <p>
     * Features:
     * - Closed network with 10 jobs circulating through delay→fork-join
     * - Branch 1: Single Queue1, Branch 2: Serial Queue2→Queue3
     * - Different service rates creating bottleneck effects
     * - Multiple solver comparison (JMT and MVA)
     */
    public static void fj_asymm() throws Exception {
        Network model = ForkJoinModel.fj_asymm();
        
        // Solve with multiple methods as in notebook
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (int i = 0; i < solvers.length; i++) {
            try {
                if (solvers[i] instanceof JMT) {
                    JMT solver = (JMT) solvers[i];
                    solver.getAvgTable().print();
                } else if (solvers[i] instanceof MVA) {
                    MVA solver = (MVA) solvers[i];
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error: " + e.getMessage());
            }
        }
        pauseForUser();
    }
    
    /**
     * Basic closed fork-join network (fj_basic_closed.ipynb).
     * <p>
     * Demonstrates symmetric fork-join with PS scheduling in closed system.
     * <p>
     * Features:
     * - 5 jobs circulating through Delay→Fork→{Queue1,Queue2}→Join→Delay
     * - Symmetric parallel branches with identical service rates (1.0)
     * - PS scheduling for fair sharing within each queue
     * - Multiple solver comparison (JMT and MVA with fork-join extensions)
     */
    public static void fj_basic_closed() throws Exception {
        Network model = ForkJoinModel.fj_basic_closed();
        
        // Solve with multiple methods as in notebook
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model, "fork_join", "ht", "method", "amva");
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (int i = 0; i < solvers.length; i++) {
            try {
                if (solvers[i] instanceof JMT) {
                    JMT solver = (JMT) solvers[i];
                    solver.getAvgTable().print();
                } else if (solvers[i] instanceof MVA) {
                    MVA solver = (MVA) solvers[i];
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error: " + e.getMessage());
            }
        }
        pauseForUser();
    }
    
    /**
     * Basic nested fork-join network (fj_basic_nesting.ipynb).
     * <p>
     * Shows hierarchical fork-join structures with multiple levels.
     */
    public static void fj_basic_nesting() throws Exception {
        Network model = ForkJoinModel.fj_basic_nesting();
        JMT solver = new JMT(model, "seed", 23000);
        solver.getAvgTable().print();
        MVA solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Basic open fork-join network (fj_basic_open.ipynb).
     * <p>
     * Demonstrates fundamental fork-join synchronization with open arrivals,
     * showing parallel processing and synchronization effects.
     * <p>
     * Features:
     * - Fork splits jobs into parallel streams for processing
     * - Two FCFS queues with different service rates (1.0 and 2.0)
     * - Join synchronizes completion from both parallel branches
     * - Queue2 is faster, so Queue1 becomes the bottleneck
     * - Multiple solver comparison (JMT and MVA)
     */
    public static void fj_basic_open() throws Exception {
        Network model = ForkJoinModel.fj_basic_open();
        
        // Solve with multiple solvers exactly as in the notebook
        JMT solverJMT = new JMT(model, "seed", 23000, "verbose", false, "keep", false);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Complex serial fork-join network (fj_complex_serial.ipynb).
     * <p>
     * Demonstrates complex fork-join with serial processing within branches.
     * <p>
     * Features:
     * - Fork to two branches with series queues within each branch
     * - Branch 1: Queue1 → Queue4 → Queue5 (3 queues in series)
     * - Branch 2: Queue2 → Queue3 (2 queues in series)
     * - Different service rates creating complex performance interactions
     * - Shows how serial chains within parallel branches affect overall performance
     */
    public static void fj_complex_serial() throws Exception {
        Network model = ForkJoinModel.fj_complex_serial();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork-join with class switching and multiple visits (fj_cs_multi_visits.ipynb).
     * <p>
     * Demonstrates feedback loop with class switching in fork-join.
     * <p>
     * Features:
     * - Two open classes: Class1 with arrivals, Class2 generated internally
     * - Router node creates feedback from Join back to Fork
     * - Class switching capabilities showing job transformation
     * - PS scheduling for fair resource sharing
     * - Feedback loop creates multiple visits to fork-join structure
     */
    public static void fj_cs_multi_visits() throws Exception {
        Network model = ForkJoinModel.fj_cs_multi_visits();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork-join with post-fork class switching (fj_cs_postfork.ipynb).
     * <p>
     * Shows class switching after fork splitting.
     * <p>
     * Features:
     * - Two closed classes with 1 job each
     * - Class switching occurs after fork via router nodes
     * - Router nodes between Fork and Queue nodes enable class transformation
     * - PS queues with identical service rates for fairness
     * - Demonstrates how jobs can change class during fork-join processing
     */
    public static void fj_cs_postfork() throws Exception {
        Network model = ForkJoinModel.fj_cs_postfork();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork-join with pre-fork class switching (fj_cs_prefork.ipynb).
     * <p>
     * Demonstrates class switching before fork operation.
     * <p>
     * Features:
     * - Two closed classes with 10 jobs each
     * - Class switching happens before entering fork-join structure
     * - Delay1 → Delay2 path with class transformation
     * - Only transformed jobs enter the fork-join section
     * - Shows selective fork-join participation based on job class
     */
    public static void fj_cs_prefork() throws Exception {
        Network model = ForkJoinModel.fj_cs_prefork();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Deep nested fork-join network (fj_deep_nesting.ipynb).
     * <p>
     * Shows complex hierarchical fork-join structures.
     * <p>
     * Features:
     * - Nested fork-join: Fork1 splits to Queue1 and Queue2
     * - Queue1 leads to Fork2 which splits to Queue3 and Queue4
     * - Creates asymmetric tree structure with deeper nesting on one side
     * - Single job to clearly show synchronization behavior
     * - Demonstrates how nested fork-joins affect response times
     */
    public static void fj_deep_nesting() throws Exception {
        Network model = ForkJoinModel.fj_deep_nesting();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork-join with multiple delay stages (fj_delays.ipynb).
     * <p>
     * Demonstrates impact of sequential delays before fork-join.
     * <p>
     * Features:
     * - Two delay nodes in series before fork-join structure
     * - Delay1 (rate 0.5) → Delay2 (rate 2.0) → Fork
     * - Symmetric parallel queues after fork (both rate 1.0)
     * - Shows how pre-processing delays affect fork-join performance
     * - 10 jobs circulating in closed system
     */
    public static void fj_delays() throws Exception {
        Network model = ForkJoinModel.fj_delays();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork without join synchronization (fj_nojoin.ipynb).
     * <p>
     * Shows fork splitting without synchronization requirement.
     * <p>
     * Features:
     * - Fork splits arrivals to three PS queues
     * - No Join node - jobs exit independently to sink
     * - Probabilistic routing: 1/3 to each queue
     * - Different service rates: Queue1(1.0), Queue2(0.5), Queue3(0.333)
     * - Demonstrates load balancing without synchronization overhead
     */
    public static void fj_nojoin() throws Exception {
        Network model = ForkJoinModel.fj_nojoin();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Fork-join with overlapping routes (fj_route_overlap.ipynb).
     * <p>
     * Three-way fork-join demonstrating route overlap effects.
     * <p>
     * Features:
     * - Fork splits to three FCFS queues simultaneously
     * - Different service rates: Queue1(1.0), Queue2(2.0), Queue3(3.0)
     * - Join waits for all three branches to complete
     * - Queue1 is fastest, Queue3 is slowest (bottleneck)
     * - Shows impact of slowest branch on overall response time
     */
    public static void fj_route_overlap() throws Exception {
        Network model = ForkJoinModel.fj_route_overlap();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Closed serial fork-join stages (fj_serialfjs_closed.ipynb).
     * <p>
     * Two fork-join stages in series within closed network.
     * <p>
     * Features:
     * - First stage: Fork1 → {Queue1, Queue2} → Join1
     * - Second stage: Fork2 → {Queue3, Queue4} → Join2
     * - 10 jobs circulating through Delay → Stage1 → Stage2 → Delay
     * - All queues PS with identical service rates (1.0)
     * - Shows cumulative effect of cascaded fork-join stages
     */
    public static void fj_serialfjs_closed() throws Exception {
        Network model = ForkJoinModel.fj_serialfjs_closed();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Open serial fork-join stages (fj_serialfjs_open.ipynb).
     * <p>
     * Cascaded fork-join stages with open arrivals.
     * <p>
     * Features:
     * - External arrivals at rate 2.5
     * - Two fork-join stages in series
     * - Stage 1: Fork1 → {Queue1, Queue2} → Join1
     * - Stage 2: Fork2 → {Queue3, Queue4} → Join2
     * - All queues FCFS with identical service rates
     * - Demonstrates pipeline parallelism with synchronization points
     */
    public static void fj_serialfjs_open() throws Exception {
        Network model = ForkJoinModel.fj_serialfjs_open();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Multi-class three-branch fork-join (fj_threebranches.ipynb).
     * <p>
     * Asymmetric fork-join with three branches and two classes.
     * <p>
     * Features:
     * - Two closed classes with 10 jobs each
     * - Fork splits to three PS queues
     * - Queue2 → Queue3 creates series within middle branch
     * - Different service rates for each class at each queue
     * - Shows multi-class interaction in asymmetric fork-join topology
     */
    public static void fj_threebranches() throws Exception {
        Network model = ForkJoinModel.fj_threebranches();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Two-class fork-join with task multiplication (fj_twoclasses_forked.ipynb).
     * <p>
     * Multi-class fork-join with different task multiplicity.
     * <p>
     * Features:
     * - Two open classes with arrival rate 0.25 each
     * - Fork creates 2 tasks per link (task multiplication)
     * - PS scheduling for resource sharing between classes
     * - Class2 has immediate service at Queue1, exponential at Queue2
     * - Shows effect of task multiplication on system performance
     */
    public static void fj_twoclasses_forked() throws Exception {
        Network model = ForkJoinModel.fj_twoclasses_forked();
        JMT solverJMT = new JMT(model, "seed", 23000);
        MVA solverMVA = new MVA(model);
        
        Object[] solvers = {solverJMT, solverMVA};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }

    /**
     * Main method demonstrating all fork-join examples.
     */
    public static void main(String[] args) throws Exception {
        fj_nojoin();
        System.out.println("\n=== Running example: fj_basic_open ===");
        try {
            fj_basic_open();
        } catch (Exception e) {
            System.err.println("fj_basic_open failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_basic_closed ===");
        try {
            fj_basic_closed();
        } catch (Exception e) {
            System.err.println("fj_basic_closed failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_asymm ===");
        try {
            fj_asymm();
        } catch (Exception e) {
            System.err.println("fj_asymm failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_basic_nesting ===");
        try {
            fj_basic_nesting();
        } catch (Exception e) {
            System.err.println("fj_basic_nesting failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_complex_serial ===");
        try {
            fj_complex_serial();
        } catch (Exception e) {
            System.err.println("fj_complex_serial failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_cs_multi_visits ===");
        try {
            fj_cs_multi_visits();
        } catch (Exception e) {
            System.err.println("fj_cs_multi_visits failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_cs_postfork ===");
        try {
            fj_cs_postfork();
        } catch (Exception e) {
            System.err.println("fj_cs_postfork failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_cs_prefork ===");
        try {
            fj_cs_prefork();
        } catch (Exception e) {
            System.err.println("fj_cs_prefork failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_deep_nesting ===");
        try {
            fj_deep_nesting();
        } catch (Exception e) {
            System.err.println("fj_deep_nesting failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_delays ===");
        try {
            fj_delays();
        } catch (Exception e) {
            System.err.println("fj_delays failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_nojoin ===");
        try {
            fj_nojoin();
        } catch (Exception e) {
            System.err.println("fj_nojoin failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_route_overlap ===");
        try {
            fj_route_overlap();
        } catch (Exception e) {
            System.err.println("fj_route_overlap failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_serialfjs_closed ===");
        try {
            fj_serialfjs_closed();
        } catch (Exception e) {
            System.err.println("fj_serialfjs_closed failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_serialfjs_open ===");
        try {
            fj_serialfjs_open();
        } catch (Exception e) {
            System.err.println("fj_serialfjs_open failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_threebranches ===");
        try {
            fj_threebranches();
        } catch (Exception e) {
            System.err.println("fj_threebranches failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fj_twoclasses_forked ===");
        try {
            fj_twoclasses_forked();
        } catch (Exception e) {
            System.err.println("fj_twoclasses_forked failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}