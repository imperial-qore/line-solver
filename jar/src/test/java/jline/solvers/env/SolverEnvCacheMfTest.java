/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.ENV;
import jline.solvers.fluid.FLD;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Verifies the ENV mean-field FLD cache-hit aggregation on the mmap_rr_env
 * stage models. Target (MATLAB mf/FLD row): Read1 ~ 0.5740, Read2 ~ 0.4830.
 */
public class SolverEnvCacheMfTest {

    private static Network stage(double l1, double l2) {
        Network model = new Network("MMAPCache");
        Source source = new Source(model, "Source");
        Cache cache = new Cache(model, "Cache", 4, 2, ReplacementStrategy.RR);
        Sink sink = new Sink(model, "Sink");

        OpenClass rd1 = new OpenClass(model, "Read1", 0);
        OpenClass rd2 = new OpenClass(model, "Read2", 0);
        OpenClass h1 = new OpenClass(model, "Hit1", 0);
        OpenClass s1 = new OpenClass(model, "Miss1", 0);
        OpenClass h2 = new OpenClass(model, "Hit2", 0);
        OpenClass s2 = new OpenClass(model, "Miss2", 0);

        cache.setRead(rd1, new DiscreteSampler(new Matrix(new double[]{8.0 / 15, 4.0 / 15, 2.0 / 15, 1.0 / 15})));
        cache.setRead(rd2, new DiscreteSampler(new Matrix(new double[]{1.0 / 15, 2.0 / 15, 4.0 / 15, 8.0 / 15})));
        cache.setHitClass(rd1, h1);
        cache.setMissClass(rd1, s1);
        cache.setHitClass(rd2, h2);
        cache.setMissClass(rd2, s2);

        source.setArrival(rd1, new Exp(l1));
        source.setArrival(rd2, new Exp(l2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(rd1, rd1, source, cache, 1.0);
        P.set(rd2, rd2, source, cache, 1.0);
        P.set(h1, h1, cache, sink, 1.0);
        P.set(s1, s1, cache, sink, 1.0);
        P.set(h2, h2, cache, sink, 1.0);
        P.set(s2, s2, cache, sink, 1.0);
        model.link(P);
        return model;
    }

    @Test
    public void meanfieldFldCacheHitRatio() {
        int E = 2;
        Environment env = new Environment("MMPPphase", E);
        env.addStage(0, "Phase1", "bursty", stage(3.6, 0.4));
        env.addStage(1, "Phase2", "calm", stage(0.2, 0.8));
        env.addTransition(0, 1, new Exp(0.5));
        env.addTransition(1, 0, new Exp(0.5));
        env.init();

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            SolverOptions o = new SolverOptions(SolverType.FLUID);
            o.method = "rmf";
            o.timespan = new double[]{0, 50};
            o.verbose = VerboseLevel.SILENT;
            solvers[e] = new FLD(env.getModel(e), o);
        }
        SolverOptions envOpt = new SolverOptions(SolverType.ENV);
        envOpt.verbose = VerboseLevel.SILENT;
        envOpt.iter_max = 100;
        envOpt.iter_tol = 1e-4;

        ENV solver = new ENV(env, solvers, envOpt);
        solver.getAvg();

        Cache refCache = (Cache) env.getModel(0).getNodeByName("Cache");
        Matrix hr = refCache.getHitRatio();
        double read1 = hr.get(0);
        double read2 = hr.get(1);

        assertEquals(0.5740, read1, 0.02, "Read1 hit ratio");
        assertEquals(0.4830, read2, 0.02, "Read2 hit ratio");
    }
}
