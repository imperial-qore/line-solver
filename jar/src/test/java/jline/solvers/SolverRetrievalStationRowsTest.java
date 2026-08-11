package jline.solvers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The analytical retrieval analyzers must report the RETRIEVAL STATION, not just
 * the Source.
 *
 * <p>{@code QN/UN/RN/TN} are indexed by station. Both analyzers wrote them at the
 * NODE index of the retrieval queue instead, which on any model where a
 * non-station node (the Cache itself) sits between the Source and the queue put
 * the metrics outside the station range: the row vanished and the whole model
 * reported a single Source row with throughput 1. MATLAB converts first
 * ({@code sst = sn.nodeToStation(queueNodes(s))} in
 * {@code solver_nc_retrieval_analyzer.m}) and is the reference here.
 *
 * <p>A second defect on the same row: the FETCH throughput. A retrieval station
 * is entered only on a MISS, so its arrival rate is the request rate weighted by
 * the access-weighted miss probability. Both analyzers dropped that weight and
 * reported the request rate itself, which on a model whose access probabilities
 * sum to one is exactly the source rate.
 *
 * <p>Expected values are MATLAB, run on
 * {@code examples/basic/cacheModel/retrieval_simple} (3 items, capacity 1, FIFO,
 * arrival rate 1, retrieval service rate 2, INF station):
 * NC QLen = Util = 0.231399 and Tput = 0.462797; MVA 0.243091 and 0.486182.
 * python-native agrees to the last digit shown.
 */
public class SolverRetrievalStationRowsTest {

    private static final double NC_QLEN = 0.231399;
    private static final double NC_TPUT = 0.462797;
    private static final double MVA_QLEN = 0.243091;
    private static final double MVA_TPUT = 0.486182;
    private static final double TOL = 1e-4;

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** examples/basic/cacheModel/retrieval_simple, built inline. */
    private static Network retrievalSimple() {
        double[] accessProb = {0.6, 0.3, 0.1};
        int n = accessProb.length;

        Network model = new Network("DelayedHits");
        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, new Matrix("[1]"), ReplacementStrategy.FIFO);
        Queue queue = new Queue(model, "Queue", SchedStrategy.INF);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));
        queue.setService(jobClass, new Exp(2.0));

        cacheNode.setRead(jobClass, new DiscreteSampler(new Matrix(accessProb)));
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);
        cacheNode.setRetrievalSystem(jobClass, missClass, new Queue[]{queue});

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, cacheNode, 1.0);
        P.set(jobClass, jobClass, cacheNode, queue, 1.0);
        P.set(jobClass, jobClass, queue, cacheNode, 1.0);
        P.set(hitClass, hitClass, cacheNode, sink, 1.0);
        P.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(P);
        return model;
    }

    private static int rowOf(NetworkAvgTable table, String station) {
        List<String> names = table.getStationNames();
        for (int i = 0; i < names.size(); i++) {
            if (station.equals(names.get(i))) {
                return i;
            }
        }
        return -1;
    }

    private static void assertRetrievalStationReported(NetworkAvgTable table, double expectedQLen,
                                                       double expectedTput, String solver) {
        int row = rowOf(table, "Queue");
        assertTrue(row >= 0, solver + ": the retrieval station is missing from the table entirely, "
                + "stations reported were " + table.getStationNames());
        double qlen = table.getQLen().get(row);
        double util = table.getUtil().get(row);
        assertEquals(expectedQLen, qlen, TOL,
                solver + ": mean queue length at the retrieval station");
        // The analyzer defines station occupancy phi_s and reports it as both
        // QLen and Util, so a mismatch means only one of the two was written.
        assertEquals(qlen, util, TOL,
                solver + ": QLen and Util disagree at the retrieval station");
        // The FETCH throughput: a retrieval station is entered only on a MISS,
        // so this is the request rate weighted by the access-weighted miss
        // probability, NOT the request rate. Reporting the latter (1.0 here)
        // was the second defect on this row.
        double tput = table.getTput().get(row);
        assertEquals(expectedTput, tput, TOL,
                solver + ": fetch throughput at the retrieval station");
        // Little's law closes the loop: the station is INF with mean service
        // 0.5, so QLen/Tput must be that service time whatever the numbers are.
        assertEquals(0.5, qlen / tput, TOL,
                solver + ": QLen/Tput is not the retrieval service time");
    }

    @Test
    public void ncReportsTheRetrievalStationRow() {
        SolverNC solver = new SolverNC(retrievalSimple());
        assertRetrievalStationReported(solver.getAvgTable(), NC_QLEN, NC_TPUT, "NC");
    }

    @Test
    public void mvaReportsTheRetrievalStationRow() {
        SolverMVA solver = new SolverMVA(retrievalSimple());
        assertRetrievalStationReported(solver.getAvgTable(), MVA_QLEN, MVA_TPUT, "MVA");
    }
}
