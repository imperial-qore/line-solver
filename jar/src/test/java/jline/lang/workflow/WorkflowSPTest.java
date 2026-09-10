/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.workflow;

import jline.lang.processes.APH;
import jline.lang.processes.Exp;
import jline.lang.processes.Markovian;
import jline.lang.processes.PH;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Closed-form checks of the series-parallel composition of Workflow and of the
 * geometric loop.
 */
public class WorkflowSPTest {

    private static final double TOL = 1e-9;

    private static List<WorkflowActivity> list(WorkflowActivity... acts) {
        return new ArrayList<WorkflowActivity>(Arrays.asList(acts));
    }

    @Test
    public void serialOfTwoExponentialsIsErlang2() {
        Workflow wf = new Workflow("serial2");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(1.0));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        wf.addPrecedence(Workflow.Serial(a, b));

        Markovian ph = wf.toPH();
        assertEquals(2.0, ph.getMean(), TOL);
        assertEquals(0.5, ph.getSCV(), TOL);
        assertNotNull(wf.getSPTree());
    }

    @Test
    public void parallelIsTheMaximumOfTheBranches() {
        double ra = 0.5;
        double rb = 1.0 / 1.5;
        Pair<Matrix, Matrix> par = Workflow.composeParallel(Matrix.singleton(1.0), Matrix.singleton(-ra),
                Matrix.singleton(1.0), Matrix.singleton(-rb));
        PH ph = new PH(par.getLeft(), par.getRight());
        assertEquals(1 / ra + 1 / rb - 1 / (ra + rb), ph.getMean(), TOL);
    }

    @Test
    public void geometricLoopOfAnExponentialIsExponential() {
        double mean = 2.0;
        double count = 3.0;
        Pair<Matrix, Matrix> law = Workflow.composeLoopGeometric(Matrix.singleton(1.0),
                Matrix.singleton(-1 / mean), count);
        PH ph = new PH(law.getLeft(), law.getRight());
        assertEquals(count * mean, ph.getMean(), TOL);
        assertEquals(1.0, ph.getSCV(), TOL);
    }

    @Test
    public void geometricLoopMatchesTheCompoundMoments() {
        // Erlang-2 body of mean 2 and SCV 1/2
        Matrix alpha = new Matrix(1, 2, 2);
        alpha.set(0, 0, 1.0);
        Matrix T = new Matrix(2, 2, 4);
        T.set(0, 0, -1.0);
        T.set(0, 1, 1.0);
        T.set(1, 1, -1.0);

        double count = 3.0;
        double scvBody = 0.5;
        Pair<Matrix, Matrix> law = Workflow.composeLoopGeometric(alpha, T, count);
        PH ph = new PH(law.getLeft(), law.getRight());

        assertEquals(count * 2.0, ph.getMean(), TOL);
        assertEquals(scvBody / count + 1 - 1 / count, ph.getSCV(), TOL);
        // the order is that of the body, and the generator is cyclic
        assertEquals(2, law.getRight().getNumRows());
        assertFalse(Workflow.isAcyclicGenerator(law.getRight()));
    }

    @Test
    public void fractionalLoopCountRunsTheBodyWithThatProbability() {
        Pair<Matrix, Matrix> law = Workflow.composeLoopGeometric(Matrix.singleton(1.0),
                Matrix.singleton(-0.5), 0.25);
        PH ph = new PH(law.getLeft(), law.getRight());
        assertEquals(0.25 * 2.0, ph.getMean(), 1e-6);
    }

    @Test
    public void loopWorkflowKeepsTheMeanAndTheBodyOrder() {
        Workflow wf = new Workflow("LoopWorkflow");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(1.0));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(2.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(0.5));
        wf.addPrecedence(Workflow.Loop(a, list(b, c), 3.0));

        Markovian ph = wf.toPH();
        assertEquals(1.0 + 3 * 2.0 + 0.5, ph.getMean(), TOL);
        // A, the geometric loop over B, and C
        assertEquals(3, ph.getNumberOfPhases());
        double varTot = 1.0 + 36.0 + 0.25;
        assertEquals(varTot / (7.5 * 7.5), ph.getSCV(), TOL);
    }

    @Test
    public void orForkMixesTheBranches() {
        Workflow wf = new Workflow("BranchingWorkflow");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(1.0));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(2.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(5.0));
        WorkflowActivity d = wf.addActivity("D", Exp.fitMean(0.5));
        wf.addPrecedence(Workflow.OrFork(a, list(b, c), new double[]{0.6, 0.4}));
        wf.addPrecedence(Workflow.OrJoin(list(b, c), d));

        Markovian ph = wf.toPH();
        assertEquals(1.0 + 0.6 * 2.0 + 0.4 * 5.0 + 0.5, ph.getMean(), TOL);
        assertNotNull(wf.getSPTree());
    }

    @Test
    public void forkNestedInsideALoopIsReducedExactly() {
        Workflow wf = new Workflow("ForkInLoop");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(0.5));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(2.0));
        WorkflowActivity d = wf.addActivity("D", Exp.fitMean(1.5));
        WorkflowActivity e = wf.addActivity("E", Exp.fitMean(0.25));
        WorkflowActivity f = wf.addActivity("F", Exp.fitMean(0.75));
        wf.addPrecedence(Workflow.Loop(a, list(b, f), 2.0));
        wf.addPrecedence(Workflow.AndFork(b, list(c, d)));
        wf.addPrecedence(Workflow.AndJoin(list(c, d), e));

        Markovian ph = wf.toPH();
        double emax = 2.0 + 1.5 - 1 / (0.5 + 1 / 1.5);
        double body = 1.0 + emax + 0.25;
        assertEquals(0.5 + 2 * body + 0.75, ph.getMean(), TOL);
        assertNotNull(wf.getSPTree());
    }

    @Test
    public void incrementalRefreshMatchesARebuiltWorkflow() {
        Workflow wf = new Workflow("ForkInLoop");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(0.5));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(2.0));
        WorkflowActivity d = wf.addActivity("D", Exp.fitMean(1.5));
        WorkflowActivity e = wf.addActivity("E", Exp.fitMean(0.25));
        WorkflowActivity f = wf.addActivity("F", Exp.fitMean(0.75));
        wf.addPrecedence(Workflow.Loop(a, list(b, f), 2.0));
        wf.addPrecedence(Workflow.AndFork(b, list(c, d)));
        wf.addPrecedence(Workflow.AndJoin(list(c, d), e));
        wf.toPH();

        Workflow ref = new Workflow("ForkInLoopRef");
        WorkflowActivity a2 = ref.addActivity("A", Exp.fitMean(0.5));
        WorkflowActivity b2 = ref.addActivity("B", Exp.fitMean(1.0));
        WorkflowActivity c2 = ref.addActivity("C", Exp.fitMean(3.0));
        WorkflowActivity d2 = ref.addActivity("D", Exp.fitMean(1.5));
        WorkflowActivity e2 = ref.addActivity("E", Exp.fitMean(0.25));
        WorkflowActivity f2 = ref.addActivity("F", Exp.fitMean(0.75));
        ref.addPrecedence(Workflow.Loop(a2, list(b2, f2), 2.0));
        ref.addPrecedence(Workflow.AndFork(b2, list(c2, d2)));
        ref.addPrecedence(Workflow.AndJoin(list(c2, d2), e2));
        Markovian phRef = ref.toPH();

        wf.setActivityDemand("C", Exp.fitMean(3.0));
        Markovian phInc = wf.refreshPH();

        assertEquals(phRef.getMean(), phInc.getMean(), TOL);
        assertEquals(phRef.getSCV(), phInc.getSCV(), TOL);
        assertEquals(phRef.getNumberOfPhases(), phInc.getNumberOfPhases());
    }

    @Test
    public void meanOnlyRescaleKeepsTheShape() {
        Workflow wf = new Workflow("Rescale");
        WorkflowActivity a = wf.addActivity("A", APH.fitMeanAndSCV(2.0, 0.3));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        wf.addPrecedence(Workflow.Serial(a, b));
        wf.toPH();

        wf.setActivityDemandMean("A", 5.0);
        Markovian scaled = wf.refreshPH();

        Workflow ref = new Workflow("RescaleRef");
        WorkflowActivity a2 = ref.addActivity("A", APH.fitMeanAndSCV(5.0, 0.3));
        WorkflowActivity b2 = ref.addActivity("B", Exp.fitMean(1.0));
        ref.addPrecedence(Workflow.Serial(a2, b2));
        Markovian refPh = ref.toPH();

        assertEquals(6.0, scaled.getMean(), TOL);
        assertEquals(refPh.getSCV(), scaled.getSCV(), 1e-6);
        assertEquals(refPh.getNumberOfPhases(), scaled.getNumberOfPhases());
    }

    @Test
    public void quorumJoinIsRefusedByName() {
        Workflow wf = new Workflow("Quorum");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(1.0));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(1.0));
        WorkflowActivity d = wf.addActivity("D", Exp.fitMean(1.0));
        wf.addPrecedence(Workflow.AndFork(a, list(b, c)));
        wf.addPrecedence(Workflow.AndJoin(list(b, c), d, 1));

        IllegalStateException ex = assertThrows(IllegalStateException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                wf.toPH();
            }
        });
        assertTrue(ex.getMessage().contains("quorum"));
    }

    @Test
    public void aGraphThatIsNotSeriesParallelFallsBack() {
        Workflow wf = new Workflow("NotSP");
        WorkflowActivity a = wf.addActivity("A", Exp.fitMean(1.0));
        WorkflowActivity b = wf.addActivity("B", Exp.fitMean(1.0));
        WorkflowActivity c = wf.addActivity("C", Exp.fitMean(1.0));
        wf.addPrecedence(Workflow.Serial(a, b));
        wf.addPrecedence(Workflow.Serial(a, c));

        Markovian ph = wf.toPH();
        assertTrue(ph.getMean() > 0);
        assertNull(wf.getSPTree());
    }
}
