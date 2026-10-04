/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.lang.ModelAdapter;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Chain-aggregation strategy for {@link TransformSolve}.
 *
 * <p>Collapses every chain onto a single class, class switching disappearing
 * with it, solves that model with the CALLER'S OWN solver, and maps the chain
 * metrics back onto the classes through alpha, the per-station share of the
 * chain's visits each class carries.
 *
 * <p>WHAT IS TRADED. The aggregation is EXACT on a product-form model: the chain
 * is the unit MVA and convolution already solve in, and the deaggregation is the
 * same alpha-weighted split those solvers apply. It is an APPROXIMATION
 * otherwise, because one aggregate service law, fitted to the alpha-weighted
 * first two moments, replaces the per-class ones. A caller who needs the exact
 * multiclass answer leaves the transform off and pays the state space.
 *
 * <p>SINGLE PASS: one solve of the aggregate determines the answer, so the
 * context is not iterated and couple/converged are never called.
 *
 * <p>Mirrors MATLAB {@code solver_tr_chains_analyzer.m} and python
 * {@code transform_driver._chains_strategy}.
 */
public class ChainsStrategy implements TransformSolve.Strategy {

    /** Carries what the deaggregation needs from the aggregation. */
    public static final class ChainsContext extends TransformSolve.Expanded {
        final NetworkStruct snOrig;
        final Matrix alpha;
        final ModelAdapter.DeaggInfo deagg;

        ChainsContext(List<Network> submodels, NetworkStruct snOrig, Matrix alpha,
                      ModelAdapter.DeaggInfo deagg) {
            super(submodels, false);
            this.snOrig = snOrig;
            this.alpha = alpha;
            this.deagg = deagg;
        }
    }

    @Override
    public TransformSolve.Expanded expand(Network model, NetworkStruct sn, SolverOptions options) {
        // With one class per chain the aggregation is the IDENTITY, and
        // aggregateChains then returns no deaggregation tables at all. Refusing
        // by name beats failing later on a missing field; SolverCTMC's own
        // chain_aggregation entry applies the same guard before reaching here.
        if (sn.nchains >= sn.nclasses) {
            throw new RuntimeException("chain aggregation needs more classes than chains: this "
                    + "model has " + sn.nclasses + " classes in " + sn.nchains + " chains, so the "
                    + "transform is the identity.");
        }
        ModelAdapter.AggregateChainResult agg = ModelAdapter.aggregateChains(model, "");
        List<Network> submodels = new ArrayList<Network>();
        submodels.add(agg.getChainModel());
        return new ChainsContext(submodels, sn, agg.getAlpha(), agg.getDeaggInfo());
    }

    @Override
    public TransformSolve.Lifted lift(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res) {
        ChainsContext c = (ChainsContext) ctx;
        TransformSolve.Inner r = res.get(0);
        // ST is left empty so the per-class service times are recovered from
        // sn.rates. X is the per-chain SYSTEM throughput, which is why the
        // driver collects it on its own channel rather than from getAvg.
        Ret.snDeaggregateChainResults d = jline.api.sn.SnDeaggregateChainResults.snDeaggregateChainResults(
                c.snOrig, c.deagg.Lchain, new Matrix(0, 0), c.deagg.STchain, c.deagg.Vchain,
                c.alpha, r.Q, r.U, r.R, r.T, new Matrix(0, 0), r.X);
        return new TransformSolve.Lifted(d.Q, d.U, d.R, d.T, d.C, d.X);
    }

    @Override
    public TransformSolve.Expanded couple(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res, int e) {
        throw new UnsupportedOperationException("chain aggregation is single pass and does not couple.");
    }

    @Override
    public boolean converged(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res, int it) {
        return true;
    }
}
