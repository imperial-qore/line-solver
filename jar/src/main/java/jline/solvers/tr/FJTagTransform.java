/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import jline.lang.ModelAdapter;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.VerboseLevel;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.io.InputOutput.line_debug;

/**
 * Fork-join TAG AUGMENTATION, the transform/lift pair shared by CTMC and SSA.
 *
 * <p>{@code ModelAdapter.fjtag} is the OTHER fork-join route: where {@code mmt}
 * and {@code ht} drive an outer fixed point for MVA, NC and FLD (see
 * {@code jline.solvers.fj.FJFixedPoint}), the tag augmentation is EXACT and
 * single pass. It rewrites the model so each sibling branch carries its own
 * auxiliary class, the solver runs unchanged on that struct, and the auxiliary
 * columns are folded back at the end.
 *
 * <p>This is deliberately NOT built on the {@code FJFixedPoint} shape. That
 * driver owns a loop and reaches the inner solve through
 * {@code FJFixedPoint.InnerSolve} because MMT re-solves a transformed model
 * repeatedly. {@code fjtag} substitutes the struct and then the CALLER'S OWN
 * analyzer runs on it to completion: there is no callback seam and no second
 * pass, so the reusable unit is the transform/lift PAIR, not a driver.
 *
 * <p>Mirrors MATLAB {@code matlab/src/solvers/TR/solver_tr_fjtag_analyzer.m}
 * and python {@code line_solver/solvers/fjtag_transform.py}.
 */
public final class FJTagTransform {

    private FJTagTransform() {
    }

    /** What the expand phase hands the lift phase. */
    public static final class Context {
        /** The tag-augmented struct the solver runs on. */
        public final NetworkStruct fjsn;
        /** The original struct, whose class indices the lifted matrices use. */
        public final NetworkStruct snOrig;
        /** Auxiliary-to-original class map produced by the augmentation. */
        public final Matrix fjclassmap;
        /** Number of classes before the augmentation. */
        public final int korig;

        Context(NetworkStruct fjsn, NetworkStruct snOrig, Matrix fjclassmap, int korig) {
            this.fjsn = fjsn;
            this.snOrig = snOrig;
            this.fjclassmap = fjclassmap;
            this.korig = korig;
        }
    }

    /** The lifted metrics, in ORIGINAL class coordinates. */
    public static final class Lifted {
        public final Matrix QN;
        public final Matrix UN;
        public final Matrix RN;
        public final Matrix TN;
        public final Matrix CN;
        public final Matrix XN;
        public final Matrix AN;

        Lifted(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN, Matrix AN) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.CN = CN;
            this.XN = XN;
            this.AN = AN;
        }
    }

    /**
     * Builds the tag-augmented struct the solver is to run on.
     *
     * <p>Solver-specific refusals (a CTMC transient, an SSA parallel method)
     * stay at the call site, ahead of this, exactly as they did before.
     *
     * @param model       the fork-join model
     * @param sn          its struct, before augmentation
     * @param verbose     the caller's verbosity, for the debug line
     * @param solverLabel the solver name to print, e.g. "CTMC"
     * @return the augmentation context
     */
    public static Context expand(Network model, NetworkStruct sn, VerboseLevel verbose, String solverLabel) {
        int korig = sn.nclasses;
        ModelAdapter.FJTagResult fjRet = ModelAdapter.fjtag(model);
        Context ctx = new Context(fjRet.fjsn, sn, fjRet.fjclassmap, korig);
        line_debug(verbose, String.format(
                "%s: fork-join tag augmentation, %d classes (%d auxiliary), %d fork firings",
                solverLabel, ctx.fjsn.nclasses, ctx.fjsn.nclasses - korig, ctx.fjsn.fjsync.size()));
        return ctx;
    }

    /**
     * Folds the auxiliary sibling classes back onto the original ones.
     *
     * <p>The lift must see the ORIGINAL struct, not the augmented one: a Place
     * counts tokens rather than firings, so {@code snPnAvgRates} has to rescale
     * before the arrival rates are derived, and the class indices it uses are
     * the original ones. A Join then reports the per-sibling waiting time
     * {@code QN/AN} (the JMT convention), since it sees one arrival per sibling
     * for every job it releases.
     *
     * <p>The caller is responsible for restoring its own {@code sn} reference to
     * {@link Context#snOrig} afterwards: a Java {@code NetworkStruct} is a
     * reference where the MATLAB struct is a value copy.
     */
    public static Lifted lift(Context ctx, Matrix QN, Matrix UN, Matrix RN, Matrix TN,
                              Matrix CN, Matrix XN, AvgHandle T) {
        NetworkStruct snOrig = ctx.snOrig;
        int korig = ctx.korig;
        ModelAdapter.fjFoldback(QN, UN, RN, TN, ctx.fjclassmap, korig);
        QN = Matrix.extract(QN, 0, QN.getNumRows(), 0, korig);
        UN = Matrix.extract(UN, 0, UN.getNumRows(), 0, korig);
        RN = Matrix.extract(RN, 0, RN.getNumRows(), 0, korig);
        TN = Matrix.extract(TN, 0, TN.getNumRows(), 0, korig);
        CN = Matrix.extract(CN, 0, CN.getNumRows(), 0, Math.min(korig, CN.getNumCols()));
        XN = Matrix.extract(XN, 0, XN.getNumRows(), 0, Math.min(korig, XN.getNumCols()));
        jline.api.sn.SnPnAvgRates.snPnAvgRates(snOrig, QN, TN, null, RN);
        Matrix AN = snGetArvRFromTput(snOrig, TN, T);
        for (int ind = 0; ind < snOrig.nnodes; ind++) {
            if (snOrig.nodetype.get(ind) == NodeType.Join) {
                int ist = (int) snOrig.nodeToStation.get(ind);
                for (int r = 0; r < korig; r++) {
                    if (AN.get(ist, r) > 0) {
                        RN.set(ist, r, QN.get(ist, r) / AN.get(ist, r));
                    }
                }
            }
        }
        return new Lifted(QN, UN, RN, TN, CN, XN, AN);
    }
}
