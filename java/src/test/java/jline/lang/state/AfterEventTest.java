package jline.lang.state;

import jline.examples.ClosedModel;
import jline.examples.OpenModel;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


public class AfterEventTest {


    @Test
    public void returnsCorrectSpaceRateProb1() {
        NetworkStruct sn = ClosedModel.ex4_line().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 0, 1, 1, 0, 1, 1}});

        Matrix outspace = new Matrix(1, 7);
        outspace.fromArray2D(new int[][]{{0, 4, 1, 1, 0, 1, 1}});
        Matrix outrate = new Matrix(1, 1);
        outrate.set(0, 0, -1);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1);

        EventResult result = State.afterEvent(sn, 2, inspace, EventType.ARV, 3, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);

    }

    @Test
    public void returnsCorrectSpaceRateProb2() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 4);
        inspace.fromArray2D(new int[][]{{4, 8, 3, 2}});

        Matrix outspace = new Matrix(1, 4);
        outspace.fromArray2D(new int[][]{{4, 8, 3, 1}});
        Matrix outrate = new Matrix(1, 1);
        outrate.set(0, 0, 2);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1);

        EventResult result = State.afterEvent(sn, 0, inspace, EventType.DEP, 3, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);

    }

    @Test
    public void returnsCorrectSpaceRateProb3() {
        NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
        Matrix inspace = new Matrix(1, 3);
        inspace.fromArray2D(new int[][]{{1, 0, 0}});

        Matrix outspace = new Matrix(1, 3);
        outspace.fromArray2D(new int[][]{{0, 0, 0}});
        Matrix outrate = new Matrix(1, 1);
        outrate.set(0, 0, 1);
        Matrix outprob = new Matrix(1, 1);
        outprob.set(0, 0, 1);
        EventResult result = State.afterEvent(sn, 1, inspace, EventType.DEP, 0, false);


        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb4() {
        NetworkStruct sn = ClosedModel.ex4_line().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{4, 2, 1, 1, 1, 1, 1}});

        Matrix outspace = new Matrix(4, 7);
        outspace.fromArray2D(new int[][]{{0, 4, 1, 2, 0, 1, 1}, {0, 4, 1, 2, 0, 1, 1}, {0, 4, 1, 2, 1, 0, 1}, {0, 4, 1, 2, 1, 0, 1}});
        Matrix outrate = new Matrix(4, 1);
        outrate.fill(Double.NaN);
        Matrix outprob = new Matrix(4, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 2, inspace, EventType.DEP, 2, false);
        assertEquals(outspace, result.outspace);
        for (int i = 0; i < 4; i++) {
            assertTrue(Double.isNaN(result.outrate.get(i)));
        }
        assertEquals(outprob, result.outprob);

    }

    @Test
    public void returnsCorrectSpaceRateProb5() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{4, 2, 1, 1, 1, 1, 1}});

        Matrix outspace = new Matrix(2, 7);
        outspace.fromArray2D(new int[][]{{0, 4, 0, 2, 1, 1, 1}, {0, 4, 0, 1, 2, 1, 1}});
        Matrix outrate = new Matrix(2, 1);
        outrate.fromArray2D(new int[][]{{1}, {0}});
        Matrix outprob = new Matrix(3, 1);
        outprob.fromArray2D(new int[][]{{1}, {1}, {0}});


        EventResult result = State.afterEvent(sn, 1, inspace, EventType.DEP, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb6() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 2, 1, 1, 0, 1, 1}});


        Matrix outspace = new Matrix(4, 7);
        outspace.fromArray2D(new int[][]{{0, 0, 1, 1, 0, 1, 1}, {0, 0, 1, 0, 1, 1, 1}, {0, 0, 1, 1, 0, 1, 1}, {0, 0, 1, 0, 1, 1, 1}});
        Matrix outrate = new Matrix(4, 1);
        outrate.fill(0.0);
        Matrix outprob = new Matrix(10, 1);
        outprob.fill(0.0);


        EventResult result = State.afterEvent(sn, 1, inspace, EventType.DEP, 1, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb7() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 0, 0, 1, 0, 0, 1}});


        Matrix outspace = new Matrix(2, 7);
        outspace.fromArray2D(new int[][]{{0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 1}});
        Matrix outrate = new Matrix(2, 1);
        outrate.fill(0.0);
        Matrix outprob = new Matrix(2, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 1, inspace, EventType.DEP, 1, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb8() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 0, 0, 1, 0, 0, 1}});

        Matrix outspace = new Matrix(1, 7);
        outspace.fromArray2D(new int[][]{{0, 0, 0, 1, 0, 0, 0}});
        Matrix outrate = new Matrix(1, 1);
        outrate.fill(1.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 2, inspace, EventType.DEP, 3, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb9() {
        NetworkStruct sn = ClosedModel.ex4().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 0, 1, 1, 0, 1, 1}});

        Matrix outspace = new Matrix(1, 7);
        outspace.fromArray2D(new int[][]{{0, 0, 1, 0, 1, 1, 1}});
        Matrix outrate = new Matrix(1, 1);
        outrate.fill(1.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 1, inspace, EventType.PHASE, 1, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb10() {
        NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
        Matrix inspace = new Matrix(1, 4);
        inspace.fromArray2D(new int[][]{{1, 0, 1, 0}});

        Matrix outspace = new Matrix(1, 4);
        outspace.fromArray2D(new int[][]{{0,1,1,0}});
        Matrix outrate = new Matrix(1, 1);
        outrate.fill(3.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 0, inspace, EventType.PHASE, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb11() {
        NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
        Matrix inspace = new Matrix(1, 3);
        inspace.fromArray2D(new int[][]{{1, 1, 1}});

        Matrix outspace = new Matrix(2, 3);
        outspace.fromArray2D(new int[][]{{0,2,1}, {2,0,1}});
        Matrix outrate = new Matrix(2, 1);
        outrate.fromArray2D(new int[][]{{0}, {0}});
        Matrix outprob = new Matrix(2, 1);
        outprob.fromArray2D(new int[][]{{1}, {1}});

        EventResult result = State.afterEvent(sn, 1, inspace, EventType.PHASE, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb12() {
        NetworkStruct sn = ClosedModel.ex4_line().getStruct(false);
        Matrix inspace = new Matrix(1, 7);
        inspace.fromArray2D(new int[][]{{0, 2, 0, 0, 0, 0, 0}});

        Matrix outspace = new Matrix(1, 7);
        outspace.fromArray2D(new int[][]{{0,2,1,0,0,0,0}});
        Matrix outrate = new Matrix(1, 1);
        outrate.fill(-1.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 1, inspace, EventType.ARV, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb13() {
        NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
        Matrix inspace = new Matrix(1, 4);
        inspace.fromArray2D(new int[][]{{0, 0, 2, 0}});

        Matrix outspace = new Matrix(2, 4);
        outspace.fromArray2D(new int[][]{{1,0,2,0}, {0,1,2,0,}});
        Matrix outrate = new Matrix(2, 1);
        outrate.fill(-1.0);
        Matrix outprob = new Matrix(2, 1);
        outprob.fromArray2D(new int[][]{{1}, {0}});

        EventResult result = State.afterEvent(sn, 0, inspace, EventType.ARV, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }

    @Test
    public void returnsCorrectSpaceRateProb14() {
        NetworkStruct sn = OpenModel.ex1_line().getStruct(false);
        Matrix inspace = new Matrix(1, 2);
        inspace.fromArray2D(new double[][]{{Double.POSITIVE_INFINITY, 1.0}});

        Matrix outspace = new Matrix(1, 2);
        outspace.fromArray2D(new double[][]{{Double.POSITIVE_INFINITY, 1.0}});
        Matrix outrate = new Matrix(1, 1);
        outrate.fill(0.0);
        Matrix outprob = new Matrix(1, 1);
        outprob.fill(1.0);

        EventResult result = State.afterEvent(sn, 2, inspace, EventType.ARV, 0, false);
        assertEquals(outspace, result.outspace);
        assertEquals(outrate, result.outrate);
        assertEquals(outprob, result.outprob);
    }


}
