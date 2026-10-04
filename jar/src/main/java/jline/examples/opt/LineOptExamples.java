package jline.examples.opt;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.opt.OptimizationProblem;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.MinimizeSystemResponseTime;
import jline.opt.objectives.ResponseTimeConstraint;
import jline.opt.objectives.UtilizationConstraint;
import jline.opt.pareto.ParetoPoint;
import jline.opt.pareto.ParetoSweep;
import jline.opt.results.OptimizationResult;
import jline.opt.results.WorkflowResult;
import jline.opt.decomposition.DecompositionWorkflow;
import jline.opt.solver.BisectionSolver;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.JobPopulation;
import jline.opt.variables.RoutingProbabilities;
import jline.opt.variables.ServerAllocation;
import jline.opt.variables.ServiceRate;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Runnable demonstrations of the line-opt optimization framework, mirroring the
 * scripts in {@code python/examples/opt/} and {@code matlab/examples/opt/}.
 * Each method builds a LINE model, defines an optimization problem, solves it,
 * and prints the result.
 */
public class LineOptExamples {

    public static void main(String[] args) {
        serverSizing();
        bisectionSizing();
        serviceRate();
        loadBalancing();
        populationSizing();
        robustSizing();
        paretoFrontier();
        decomposition();
    }

    private static Map<Station, Double> cost(Station s, double c) {
        Map<Station, Double> m = new LinkedHashMap<Station, Double>();
        m.put(s, c);
        return m;
    }

    private static Network mmc(double arrivalRate) {
        Network model = new Network("MMc");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(arrivalRate));
        queue.setService(jobs, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Minimum servers under a 50% utilization SLA (differential evolution). */
    public static void serverSizing() {
        Network model = mmc(3.0);
        Queue queue = (Queue) model.getNodes().get(1);
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new ServerAllocation(queue, 1, 10));
        p.setObjective(new MinimizeCost(cost(queue, 10.0), null, null, null));
        p.addConstraint(new UtilizationConstraint(queue, 0.5));
        OptimizationResult r = p.solve(new LineOptSolverOptions().setSeed(42));
        System.out.println("[server_sizing] servers=" + r.getVariableValue("Server_servers")
                + " cost=" + r.objectiveValue + " evals=" + r.modelEvaluations);
    }

    /** Same sizing solved exactly with BisectionSolver. */
    public static void bisectionSizing() {
        Network model = mmc(3.0);
        Queue queue = (Queue) model.getNodes().get(1);
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new ServerAllocation(queue, 1, 10));
        p.setObjective(new MinimizeCost(cost(queue, 10.0), null, null, null));
        p.addConstraint(new UtilizationConstraint(queue, 0.5));
        OptimizationResult r = new BisectionSolver(p).solve();
        System.out.println("[bisection_sizing] servers=" + r.getVariableValue("Server_servers")
                + " cost=" + r.objectiveValue + " probes=" + r.modelEvaluations);
    }

    /** Cheapest continuous service rate meeting an RT SLA (theory 5.0). */
    public static void serviceRate() {
        Network model = new Network("MM1");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(3.0));
        queue.setService(jobs, new Exp(4.0));
        model.link(Network.serialRouting(source, queue, sink));
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new ServiceRate(queue, jobs, 3.5, 8.0));
        p.setObjective(new MinimizeCost(null, cost(queue, 20.0), null, null));
        p.addConstraint(new ResponseTimeConstraint(queue, jobs, 0.5));
        OptimizationResult r = p.solve(new LineOptSolverOptions().setSeed(42).setMaxIterations(60));
        System.out.println("[service_rate] rate=" + r.getVariableValue("Server_Jobs_rate")
                + " cost=" + r.objectiveValue);
    }

    /** Routing split minimizing end-to-end response time (theory fast 0.743). */
    public static void loadBalancing() {
        Network model = new Network("LoadBalance");
        Source source = new Source(model, "S");
        Queue fast = new Queue(model, "Fast", SchedStrategy.PS);
        Queue slow = new Queue(model, "Slow", SchedStrategy.PS);
        Sink sink = new Sink(model, "K");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(2.0));
        fast.setService(jobs, new Exp(4.0));
        slow.setService(jobs, new Exp(2.5));
        model.addLink(source, fast);
        model.addLink(source, slow);
        model.addLink(fast, sink);
        model.addLink(slow, sink);
        OptimizationProblem p = new OptimizationProblem(model);
        List<Node> targets = new ArrayList<Node>();
        targets.add(fast);
        targets.add(slow);
        p.addVariable(new RoutingProbabilities(jobs, source, targets));
        p.setObjective(new MinimizeSystemResponseTime(jobs));
        OptimizationResult r = p.solve(new LineOptSolverOptions().setSeed(42).setMaxIterations(60));
        double[] probs = (double[]) r.getVariableValue("Jobs_routing_from_S");
        System.out.println("[load_balancing] fast=" + probs[0] + " slow=" + probs[1]);
    }

    /** Largest population sustaining an interactive SLA (bisection max_feasible). */
    public static void populationSizing() {
        Network model = new Network("Interactive");
        Delay think = new Delay(model, "Think");
        Queue server = new Queue(model, "AppServer", SchedStrategy.PS);
        ClosedClass users = new ClosedClass(model, "Users", 1, think);
        think.setService(users, new Exp(1.0));
        server.setService(users, new Exp(2.0));
        model.link(Network.serialRouting(think, server));
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new JobPopulation(users, 1, 50));
        p.setObjective(new MinimizeCost(null, null, null, null));
        p.addConstraint(new ResponseTimeConstraint(server, users, 2.0));
        OptimizationResult r = new BisectionSolver(p, "max_feasible").solve();
        System.out.println("[population_sizing] users=" + r.getVariableValue("Users_population")
                + " evals=" + r.modelEvaluations);
    }

    /** Worst-case sizing across an average and a peak workload scenario. */
    public static void robustSizing() {
        Network base = mmc(3.0);
        Network peak = mmc(4.5);
        Queue queue = (Queue) base.getNodes().get(1);
        OptimizationProblem p = new OptimizationProblem(base);
        p.addVariable(new ServerAllocation(queue, 1, 12));
        p.setObjective(new MinimizeCost(cost(queue, 10.0), null, null, null));
        p.addConstraint(new UtilizationConstraint(queue, 0.5));
        p.addScenario(peak, 1.0);
        OptimizationResult r = new BisectionSolver(p).solve();
        System.out.println("[robust_sizing] servers=" + r.getVariableValue("Server_servers")
                + " cost=" + r.objectiveValue);
    }

    /** Cost vs utilization-bound frontier via ParetoSweep. */
    public static void paretoFrontier() {
        Network model = mmc(3.0);
        final Queue queue = (Queue) model.getNodes().get(1);
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new ServerAllocation(queue, 1, 10));
        p.setObjective(new MinimizeCost(cost(queue, 10.0), null, null, null));
        double[] eps = {0.3, 0.4, 0.5, 0.6, 0.75, 0.9};
        ParetoSweep sweep = new ParetoSweep(p, new ParetoSweep.ConstraintFactory() {
            public Constraint make(double epsilon) {
                return new UtilizationConstraint(queue, epsilon);
            }
        }, eps, "bisection");
        sweep.solve();
        System.out.println("[pareto_frontier]");
        for (ParetoPoint pt : sweep.getFrontier()) {
            System.out.println("  util <= " + pt.epsilon + " -> " + pt.objectiveValue
                    + " (" + pt.result.getVariableValue("Server_servers") + ")");
        }
    }

    /** Two-variable problem via DecompositionWorkflow (Gauss-Seidel). */
    public static void decomposition() {
        Network model = mmc(3.0);
        Queue queue = (Queue) model.getNodes().get(1);
        OpenClass jobs = (OpenClass) model.getClasses().get(0);
        OptimizationProblem p = new OptimizationProblem(model);
        p.addVariable(new ServerAllocation(queue, 1, 10));
        p.addVariable(new ServiceRate(queue, jobs, 1.0, 4.0));
        p.setObjective(new MinimizeCost(cost(queue, 10.0), cost(queue, 20.0), null, null));
        p.addConstraint(new ResponseTimeConstraint(queue, jobs, 0.5));
        DecompositionWorkflow workflow = p.decompose();
        workflow.autoDecompose();
        workflow.setSolverOptions(new LineOptSolverOptions().setSeed(42));
        WorkflowResult r = workflow.solveSequential(4, 1e-3);
        System.out.println("[decomposition] converged=" + r.converged
                + " servers=" + r.getFinalVariableValue("Server_servers")
                + " rate=" + r.getFinalVariableValue("Server_Jobs_rate")
                + " obj=" + r.finalObjective);
    }
}
