/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Disabled;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVA;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

/**
 * Examples of models with load-dependent stations
 */
public class LoadDependentModel {
    /**
     * Basic load-dependent queue with FCFS scheduling.
     * <p>
     * Features:
     * - Closed network with 16 jobs and delay node
     * - FCFS queue with load-dependent service capacity
     * - Service capacity increases linearly up to 2 servers
     * - Alpha function: min(jobs+1, 2) servers available
     * - Demonstrates load-dependent server allocation
     *
     * @return configured load-dependent network model
     */
    public static Network ld_multiserver_fcfs() {
        int N = 16; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node2.setService(jobclass1, Exp.fitMean(1.500));
        Matrix alpha = new Matrix(1, N);
        for (int i = 0; i < N; i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Multi-class load-dependent network with PS scheduling.
     * <p>
     * Features:
     * - Two closed classes: Class1 (4 jobs), Class2 (2 jobs)
     * - PS queue with total load-dependent capacity
     * - Service capacity based on total population across both classes
     * - Different service rates for each class
     * - Demonstrates multi-class load dependence
     *
     * @return configured multi-class load-dependent model
     */
    public static Network ld_multiserver_ps_twoclasses() {
        int N = 4; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N / 2.0, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node1.setService(jobclass2, Exp.fitMean(2.00));
        node2.setService(jobclass1, Exp.fitMean(1.500));
        node2.setService(jobclass2, Exp.fitMean(2.500));
        Matrix alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Three-station network with multiple load-dependent queues.
     * <p>
     * Features:
     * - Two closed classes with different populations
     * - Two PS queues (Queue1, Queue2) both with load dependence
     * - Serial routing: Delay → Queue1 → Queue2 → Delay
     * - Each queue has capacity up to 3 servers
     * - Different service rates at each station per class
     *
     * @return configured multi-station load-dependent model
     */
    public static Network ld_multiserver_ps() {
        int N = 4; // number of jobs
        int c = 3; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N / 2.0, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node1.setService(jobclass2, Exp.fitMean(2.00));

        node2.setService(jobclass1, Exp.fitMean(1.500));
        node2.setService(jobclass2, Exp.fitMean(2.500));
        Matrix alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        node3.setService(jobclass1, Exp.fitMean(3.500));
        node3.setService(jobclass2, Exp.fitMean(4.500));
        alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node3.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }


    /**
     * Class-dependent service capacity model.
     * <p>
     * Features:
     * - Two closed classes with different populations
     * - PS queue with class-dependent capacity function
     * - Service capacity depends only on Class1 population
     * - Demonstrates selective load dependence by class
     * - Custom SerializableFunction for capacity calculation
     *
     * @return configured class-dependent network model
     */
    public static Network ld_class_dependence() {
        int N = 16;
        int c = 2;

        Network cdmodel = new Network("model");

        Delay node1 = new Delay(cdmodel, "Delay");
        Queue node2 = new Queue(cdmodel, "Queue1", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(cdmodel, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(cdmodel, "Class2", N / 2.0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.0));
        node1.setService(jobclass2, Exp.fitMean(2.0));

        node2.setService(jobclass1, Exp.fitMean(1.5));
        node2.setService(jobclass2, Exp.fitMean(2.5));

        // beta_{i,r}(n_{i,r}): class 1 scales up to c servers with its OWN count,
        // class 2 is always single-server. Both entries read only their own
        // marginal, so the demands satisfy the product-form recurrence and this
        // is setClassDependence, not setJointDependence. Peak rate scaling
        // [c 1] per class normalizes Util = T*S/peak. See ld_joint_dependence()
        // for the non-product-form twin, whose scalar eta reads a foreign
        // marginal.
        SerializableFunction<Matrix, Matrix> beta =
                ni -> {
                    double class1Jobs = ni.get(0, 0); // ni is a 1xR row of per-class counts
                    Matrix out = new Matrix(1, 2);
                    out.set(0, 0, Math.min(class1Jobs, c));
                    out.set(0, 1, 1.0);
                    return out;
                };
        Matrix cdPeak = new Matrix(1, 2);
        cdPeak.set(0, 0, c);
        cdPeak.set(0, 1, 1.0);
        node2.setClassDependence(beta, cdPeak);

        RoutingMatrix routingMatrix = cdmodel.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);

        cdmodel.link(routingMatrix);
        return cdmodel;
    }

    /**
     * Joint-dependent (non-product-form) model, the twin of
     * {@link #ld_class_dependence()} and the port of MATLAB
     * ld_joint_dependence.m.
     *
     * <p>eta_i(n) reads the class-1 marginal only and returns a scalar shared by
     * every class. Because it depends on a foreign class marginal rather than
     * the own-class one, it is NON-product-form and must be declared through
     * setJointDependence; setClassDependence is reserved for the product-form
     * beta_{i,r}(n_{i,r}).</p>
     */
    public static Network ld_joint_dependence() {
        int N = 16;
        int c = 2;

        Network jdmodel = new Network("model");

        Delay node1 = new Delay(jdmodel, "Delay");
        Queue node2 = new Queue(jdmodel, "Queue1", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(jdmodel, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(jdmodel, "Class2", N / 2.0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.0));
        node1.setService(jobclass2, Exp.fitMean(2.0));

        node2.setService(jobclass1, Exp.fitMean(1.5));
        node2.setService(jobclass2, Exp.fitMean(2.5));

        SerializableFunction<Matrix, Matrix> eta =
                ni -> {
                    double class1Jobs = ni.get(0, 0); // ni is a 1xR row of per-class counts
                    Matrix out = new Matrix(1, 1);
                    out.set(0, 0, Math.min(class1Jobs, c));
                    return out;
                };
        Matrix jdPeak = new Matrix(1, 1);
        jdPeak.set(0, 0, c); // peak rate scaling = c (Util = T*S/c)
        node2.setJointDependence(eta, jdPeak);

        RoutingMatrix routingMatrix = jdmodel.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);

        jdmodel.link(routingMatrix);
        return jdmodel;
    }

    /**
     * Globally state-dependent (Whittle) model: setGlobalDependence declares a
     * rate scaling phi(n) over the FULL (nstations x nclasses) population matrix,
     * not just the population local to one station. Here two PS stations share one
     * unit of capacity, phi_s(n) = n_s/|n|, the single-link allocation every
     * alpha-fair rule collapses to.
     *
     * <p>phi satisfies the Whittle balance property
     * phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t), so the chain is reversible,
     * has the product form pi(n) ~ Phi(n) prod rho^n and is INSENSITIVE: the means
     * do not change when the exponential service is replaced by an Erlang or a
     * hyperexponential of the same mean.
     *
     * <p>Only SolverCTMC and SolverSSA declare the GlobalDependence feature; every
     * other solver rejects the model rather than solving it unscaled. SolverSSA
     * carries the same factorization on the sample path, on its serial engine:
     * the NRM's propensity closures see one station's population slice and never
     * the whole population matrix phi reads.
     *
     * @return configured globally state-dependent network model
     */
    public static Network ld_global_dependence() {
        int N = 3;

        Network gdmodel = new Network("model");

        Queue node1 = new Queue(gdmodel, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(gdmodel, "Queue2", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(gdmodel, "Class1", N, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.0));
        node2.setService(jobclass1, Exp.fitMean(0.5));

        RoutingMatrix routingMatrix = gdmodel.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        gdmodel.link(routingMatrix);

        final int M = gdmodel.getNumberOfStations();
        final int K = gdmodel.getNumberOfClasses();
        // phi as SolverCTMC sees it: an (nstations x nclasses) matrix of scalings
        SerializableFunction<Matrix, Matrix> phi =
                n -> {
                    Matrix v = Matrix.ones(M, K);
                    double tot = 0;
                    for (int i = 0; i < M; i++) {
                        for (int r = 0; r < K; r++) tot += n.get(i, r);
                    }
                    if (tot > 0) {
                        for (int i = 0; i < M; i++) {
                            double ni = 0;
                            for (int r = 0; r < K; r++) ni += n.get(i, r);
                            for (int r = 0; r < K; r++) v.set(i, r, ni / tot);
                        }
                    }
                    return v;
                };
        // peak scaling is 1: no station ever receives more than the whole link
        Matrix gdPeak = new Matrix(1, 1);
        gdPeak.set(0, 0, 1.0);
        gdmodel.setGlobalDependence(phi, gdPeak);

        return gdmodel;
    }

    /**
     * Balance function of balanced fairness over the bounded lattice:
     * Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s), with Phi(0) = 1.
     *
     * @param A      link/route incidence, one row per link
     * @param C      link capacities
     * @param cutoff per-route truncation
     * @return the balance function indexed by the mixed-radix state code
     */
    private static double[] bfBalanceFunction(int[][] A, double[] C, int cutoff) {
        int S = A[0].length;
        int base = cutoff + 1;
        int ns = 1;
        for (int s = 0; s < S; s++) ns *= base;
        double[] Phi = new double[ns];
        Phi[0] = 1.0;
        // walk the states in increasing total population, so that every
        // Phi(n-e_s) the recursion reads has already been assigned
        java.util.List<int[]> order = new java.util.ArrayList<int[]>();
        for (int k = 0; k < ns; k++) {
            int[] n = new int[S];
            int rem = k;
            for (int s = 0; s < S; s++) {
                n[s] = rem % base;
                rem /= base;
            }
            order.add(n);
        }
        java.util.Collections.sort(order, new java.util.Comparator<int[]>() {
            public int compare(int[] a, int[] b) {
                int sa = 0, sb = 0;
                for (int i = 0; i < a.length; i++) { sa += a[i]; sb += b[i]; }
                return Integer.compare(sa, sb);
            }
        });
        for (int oi = 0; oi < order.size(); oi++) {
            int[] n = order.get(oi);
            int tot = 0;
            for (int s = 0; s < S; s++) tot += n[s];
            if (tot == 0) continue;
            double best = 0;
            for (int l = 0; l < A.length; l++) {
                double acc = 0;
                for (int s = 0; s < S; s++) {
                    if (A[l][s] > 0 && n[s] > 0) {
                        n[s] -= 1;
                        acc += Phi[bfIndex(n, base)];
                        n[s] += 1;
                    }
                }
                best = Math.max(best, acc / C[l]);
            }
            Phi[bfIndex(n, base)] = best;
        }
        return Phi;
    }

    private static int bfIndex(int[] n, int base) {
        int idx = 0, mult = 1;
        for (int s = 0; s < n.length; s++) {
            idx += n[s] * mult;
            mult *= base;
        }
        return idx;
    }

    /**
     * Open Whittle network: a bandwidth-sharing model in which one route holds
     * SEVERAL links at once, which no per-station rate scaling can express. This
     * is the 2-link linear network -- route 1 crosses both links, routes 2 and 3
     * use one link each -- shared by BALANCED FAIRNESS, whose rates
     * x_s(n) = Phi(n-e_s)/Phi(n) satisfy the Whittle balance property by
     * construction. The stationary law is therefore pi(n) ~ Phi(n) prod rho_s^n_s
     * and is insensitive.
     *
     * <p>Each route is modelled as its own PS queue fed by its own open class, so
     * the queue populations ARE the coordinates of the Whittle state n. Solve with
     * SolverCTMC and an options.cutoff matching CUTOFF below.
     *
     * <p>CUTOFF is 3 here rather than the 6 the MATLAB twin uses: the JAR sizes the
     * open state space of this model more conservatively and its memory guard
     * refuses 6. At the same cutoff all four codebases agree to the digit, and to
     * the closed-form product form.
     *
     * @return configured open bandwidth-sharing model
     */
    public static Network ld_whittle_bandwidth() {
        final int[][] A = {{1, 1, 0}, {1, 0, 1}};
        final double[] C = {1.0, 1.0};
        final double[] nu = {0.20, 0.30, 0.30};
        final double[] mu = {1.00, 1.00, 1.00};
        final int CUTOFF = 3;
        final int S = 3;
        final int base = CUTOFF + 1;
        final double[] Phi = bfBalanceFunction(A, C, CUTOFF);

        Network model = new Network("model");
        Source source = new Source(model, "Source");
        Queue[] routes = new Queue[S];
        for (int s = 0; s < S; s++) {
            routes[s] = new Queue(model, "Route" + (s + 1), SchedStrategy.PS);
        }
        Sink sink = new Sink(model, "Sink");
        OpenClass[] classes = new OpenClass[S];
        for (int s = 0; s < S; s++) {
            classes[s] = new OpenClass(model, "Route" + (s + 1) + "Flows");
            source.setArrival(classes[s], new Exp(nu[s]));
        }
        for (int s = 0; s < S; s++) {
            for (int t = 0; t < S; t++) {
                if (s == t) {
                    routes[t].setService(classes[s], new Exp(mu[s]));
                } else {
                    routes[t].setService(classes[s], new Disabled());
                }
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int s = 0; s < S; s++) {
            P.set(classes[s], classes[s], source, routes[s], 1.0);
            P.set(classes[s], classes[s], routes[s], sink, 1.0);
        }
        model.link(P);

        NetworkStruct sn = model.getStruct();
        final int[] idx = new int[S];
        for (int s = 0; s < S; s++) {
            idx[s] = (int) sn.nodeToStation.get(model.getNodeIndex(routes[s]));
        }
        final int M = model.getNumberOfStations();
        final int K = model.getNumberOfClasses();
        model.setGlobalDependence(n -> {
            int[] npop = new int[S];
            for (int s = 0; s < S; s++) npop[s] = (int) Math.round(n.get(idx[s], s));
            Matrix v = Matrix.ones(M, K);
            int tot = 0;
            for (int s = 0; s < S; s++) tot += npop[s];
            if (tot == 0) return v;
            double den = Phi[bfIndex(npop, base)];
            for (int s = 0; s < S; s++) {
                if (npop[s] > 0) {
                    npop[s] -= 1;
                    v.set(idx[s], s, Phi[bfIndex(npop, base)] / den);
                    npop[s] += 1;
                } else {
                    v.set(idx[s], s, 0.0);
                }
            }
            return v;
        }, Matrix.ones(1, 1));

        return model;
    }

    /**
     * The 4-station tandem the `fes_*` and `ld_fes_*` scripts aggregate.
     *
     * <p>Delay -> Q1 -> Q2 -> Q3 -> Delay, every queue processor-sharing. The
     * second class is built only when a second mean is supplied, which is what
     * separates the single-class scripts from the two-class ones.
     *
     * @param name        model name
     * @param delayName   name of the think-time station, which the reference varies
     * @param n1          population of Class1
     * @param n2          population of Class2, ignored when the model is single-class
     * @param delayMeans  per-class think time; length 1 makes the model single-class
     * @param q1Means     per-class mean service at the first queue
     * @param q2Means     per-class mean service at the second queue
     * @param q3Means     per-class mean service at the third queue
     * @param queueNames  the three queue names
     * @return the tandem, linked
     */
    public static Network fes_tandem(String name, String delayName, double n1, double n2,
                                     double[] delayMeans, double[] q1Means, double[] q2Means,
                                     double[] q3Means, String[] queueNames) {
        Network model = new Network(name);
        Delay delay = new Delay(model, delayName);
        Queue q1 = new Queue(model, queueNames[0], SchedStrategy.PS);
        Queue q2 = new Queue(model, queueNames[1], SchedStrategy.PS);
        Queue q3 = new Queue(model, queueNames[2], SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", n1, delay, 0);
        delay.setService(class1, Exp.fitMean(delayMeans[0]));
        q1.setService(class1, Exp.fitMean(q1Means[0]));
        q2.setService(class1, Exp.fitMean(q2Means[0]));
        q3.setService(class1, Exp.fitMean(q3Means[0]));

        ClosedClass class2 = null;
        if (delayMeans.length > 1) {
            class2 = new ClosedClass(model, "Class2", n2, delay, 0);
            delay.setService(class2, Exp.fitMean(delayMeans[1]));
            q1.setService(class2, Exp.fitMean(q1Means[1]));
            q2.setService(class2, Exp.fitMean(q2Means[1]));
            q3.setService(class2, Exp.fitMean(q3Means[1]));
        }

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, q1, 1.0);
        P.set(class1, class1, q1, q2, 1.0);
        P.set(class1, class1, q2, q3, 1.0);
        P.set(class1, class1, q3, delay, 1.0);
        if (class2 != null) {
            P.set(class2, class2, delay, q1, 1.0);
            P.set(class2, class2, q1, q2, 1.0);
            P.set(class2, class2, q2, q3, 1.0);
            P.set(class2, class2, q3, delay, 1.0);
        }
        model.link(P);
        return model;
    }

    /** The two-class tandem of `fes_aggregation`, N = (3, 2). */
    public static Network fes_aggregation() {
        return fes_tandem("OriginalModel", "ThinkTime", 3, 2,
                          new double[] {5.0, 4.0}, new double[] {1.5, 2.0},
                          new double[] {1.0, 1.2}, new double[] {0.8, 1.0},
                          new String[] {"Queue1", "Queue2", "Queue3"});
    }

    /** The single-class tandem of `fes_single_class`, N = 5. */
    public static Network fes_single_class() {
        return fes_tandem("OriginalModel", "ThinkTime", 5, 0,
                          new double[] {5.0}, new double[] {1.5},
                          new double[] {1.0}, new double[] {0.8},
                          new String[] {"Queue1", "Queue2", "Queue3"});
    }

    /** The single-class tandem of `ld_fes_singleclass`, N = 5, queues named Q1..Q3. */
    public static Network ld_fes_singleclass() {
        return fes_tandem("OriginalModel", "Delay", 5, 0,
                          new double[] {5.0}, new double[] {1.5},
                          new double[] {1.0}, new double[] {0.8},
                          new String[] {"Q1", "Q2", "Q3"});
    }

    /** The two-class tandem of `ld_fes_multiclass`, N = (3, 2). */
    public static Network ld_fes_multiclass() {
        return fes_tandem("OriginalModel", "Delay", 3, 2,
                          new double[] {1.0, 1.5}, new double[] {0.5, 0.8},
                          new double[] {0.3, 0.6}, new double[] {0.4, 0.7},
                          new String[] {"Q1", "Q2", "Q3"});
    }

    /**
     * Main method for testing and demonstrating load-dependent examples.
     *
     * <p>Currently configured to run ld_multiserver_ps_twoclasses() and solve it
     * using the MVA solver with default method settings.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Network model = ld_multiserver_ps_twoclasses();
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.method = "default";
        MVA solver = new MVA(model, options);
        NetworkAvgTable t = solver.getAvgTable();
        t.print(options);
    }
}
