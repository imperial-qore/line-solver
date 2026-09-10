package jline.api.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.FieldMatrix;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.solvers.mam.handlers.Solver_mam_transient_qbd;
import jline.util.matrix.Matrix;

/**
 * Validates the Java transient-QBD engine (TransientQbd.transient2Open /
 * transient2) at the matrix level against reference values produced by the
 * independently validated native-Python solver, at a fixed complex argument
 * s = 0.5 + 0.7i. This checks transcription fidelity cheaply (single transform
 * evaluations) without the full inverse-Laplace sweep. Also checks the
 * getTranAvg auto-select predicate.
 */
public class TransientQbdTest {

    // Reference V(s,0,1) matrices from native-Python at s = 0.5 + 0.7i.
    static final double[][] OPEN_V01_RE = {
        {0.0486233179484,0.0120369930579,0.0275828363011,0.00848691346067,0.0594600870893,0.0175210998847},
        {0.0323160475015,0.0306497704178,0.0196925936078,0.0172001236519,0.0515573573412,0.0294105413076},
        {0.0208462192626,0.00588941328497,0.0390009857181,0.00856463685308,0.0569725506561,0.0140032292305},
        {0.0157506569403,0.0121590649,0.0248784177795,0.0247449806786,0.0469069551294,0.028867936684},
        {0.0302067400922,0.00751158453641,0.0141524366213,0.00440361939784,0.0644407982493,0.0138514759482},
        {0.0211776970287,0.0182496343918,0.0100696809317,0.00885341967688,0.0491097889479,0.0346078985934},
    };
    static final double[][] OPEN_V01_IM = {
        {-0.0441747179455,-0.0225012658448,-0.0240839049337,-0.0133328874031,-0.0979065387245,-0.0433254531235},
        {-0.0411174001079,-0.0260329703284,-0.0223714742826,-0.0150734774173,-0.0953490801685,-0.0473998940308},
        {-0.0408854182152,-0.0210120423217,-0.0265296858754,-0.0133507854597,-0.0984119530745,-0.0419467164585},
        {-0.038324268397,-0.0240732380358,-0.0237088960893,-0.0165108065481,-0.0947542146505,-0.0477952103726},
        {-0.0432595086094,-0.0214294355201,-0.0222909292873,-0.0122904118785,-0.0993108496747,-0.0416193542081},
        {-0.0398207869969,-0.0255082006573,-0.0205254343939,-0.0140817727236,-0.0947909123559,-0.0484993794336},
    };
    static final double[][] FIN_V01_RE = {
        {0.0486069719849,0.0124481913137,0.0278662119692,0.0090601896255,0.0593702739918,0.0182325878594},
        {0.0324558617162,0.0311293188559,0.0200978163473,0.0178253951923,0.051802973826,0.0302935685998},
        {0.0198697084489,0.00548671863564,0.0386065838837,0.00854564644628,0.0549546466413,0.0135353273747},
        {0.0148626703909,0.0118072878119,0.0245663799484,0.0247626748754,0.0451058200615,0.0285179445754},
        {0.0293496738186,0.00724411965777,0.0137056177488,0.00431715590664,0.0624633539973,0.0134091332257},
        {0.0204200545271,0.0180312361643,0.00969645153182,0.0088027051446,0.0473485162621,0.034280770487},
    };
    static final double[][] FIN_V01_IM = {
        {-0.0497906146806,-0.0265677964374,-0.0282231382768,-0.0166663696401,-0.110652688535,-0.0504924854403},
        {-0.0469324654609,-0.030217920601,-0.0266473994406,-0.0184933445422,-0.108546843535,-0.0547888838529},
        {-0.0455329407125,-0.0244289812224,-0.0300177393559,-0.0162128296015,-0.109047595992,-0.0480611989418},
        {-0.0431363447074,-0.0275902818236,-0.0273141811212,-0.0194468797108,-0.105770869627,-0.0541002175461},
        {-0.0480180892525,-0.0249378704172,-0.0257997089339,-0.0151610829695,-0.110096812632,-0.0478102769512},
        {-0.0447552619505,-0.0291216128453,-0.0241532470918,-0.0170277010764,-0.105974011375,-0.0548880408565},
    };

    private static Matrix mat(double[][] a) { return new Matrix(a); }
    private static Matrix eye(int n) { Matrix I = new Matrix(n, n); for (int i = 0; i < n; i++) I.set(i, i, 1.0); return I; }

    private void assertMatches(FieldMatrix<Complex> V, double[][] re, double[][] im) {
        double maxErr = 0.0;
        for (int i = 0; i < re.length; i++) {
            for (int j = 0; j < re[0].length; j++) {
                Complex z = V.getEntry(i, j);
                maxErr = Math.max(maxErr, Math.hypot(z.getReal() - re[i][j], z.getImaginary() - im[i][j]));
            }
        }
        assertEquals(0.0, maxErr, 1e-6);
    }

    // Blocks for the correlated MAP/MAP/1 (rho = 0.6), shared by the checks.
    private Matrix[] blocks() {
        Matrix D0 = mat(new double[][]{{-8, 1, 3}, {0, -6, 4}, {2, 0, -3}});
        Matrix D1 = mat(new double[][]{{3, 1, 0}, {0, 2, 0}, {0, 0, 1}});
        Matrix S0 = mat(new double[][]{{-3, 1}, {6, -7}});
        Matrix S1 = mat(new double[][]{{0, 2}, {1, 0}});
        double fac = (Map_lambda.map_lambda(D0, D1) / 0.6) / Map_lambda.map_lambda(S0, S1);
        Matrix Ds0 = S0.scale(fac), Ds1 = S1.scale(fac);
        Matrix Ins = eye(2), Ina = eye(3);
        Matrix Lrep = D0.krons(Ds0), Frep = D1.kron(Ins), Brep = Ina.kron(Ds1);
        Matrix Lv0 = D0.kron(Ins), F0 = D1.kron(Ins), B0 = Ina.kron(Ds1);
        Matrix LvTop = D0.add(1.0, D1).kron(Ins).add(1.0, Ina.kron(Ds0));
        return new Matrix[]{Lrep, Frep, Brep, Lv0, F0, B0, LvTop};
    }

    @SafeVarargs
    private static FieldMatrix<Complex>[] fm(FieldMatrix<Complex>... a) { return a; }

    @Test
    public void openEngineMatchesPython() {
        Matrix[] b = blocks();
        FieldMatrix<Complex>[] B = fm(null, TransientQbd.cm(b[5]), TransientQbd.cm(b[2]));
        FieldMatrix<Complex>[] L = fm(null, null, TransientQbd.cm(b[0]));
        FieldMatrix<Complex>[] F = fm(null, TransientQbd.cm(b[4]), TransientQbd.cm(b[1]));
        FieldMatrix<Complex>[] Lv = fm(null, TransientQbd.cm(b[3]), TransientQbd.cm(b[0]));
        int[] T = {0, 0, 1};
        FieldMatrix<Complex> V = TransientQbd.transient2Open(B, L, F, Lv, T, 0, 1, new Complex(0.5, 0.7));
        assertMatches(V, OPEN_V01_RE, OPEN_V01_IM);
    }

    @Test
    public void finiteEngineMatchesPython() {
        Matrix[] b = blocks();
        FieldMatrix<Complex>[] B = fm(null, TransientQbd.cm(b[2]));
        FieldMatrix<Complex>[] L = fm(null, TransientQbd.cm(b[0]));
        FieldMatrix<Complex>[] F = fm(null, TransientQbd.cm(b[1]));
        FieldMatrix<Complex>[] Lv = fm(null, TransientQbd.cm(b[3]), TransientQbd.cm(b[6]));
        int[] T = {0, 0, 2};
        FieldMatrix<Complex> V = TransientQbd.transient2(B, L, F, Lv, T, 0, 1, new Complex(0.5, 0.7));
        assertMatches(V, FIN_V01_RE, FIN_V01_IM);
    }

    @Test
    public void autoSelectPredicate() {
        // Correlated MAP arrival + MAP service -> Laplace transient QBD applies.
        Network m1 = new Network("MAPMAP1");
        Source s1 = new Source(m1, "S");
        Queue q1 = new Queue(m1, "Q", SchedStrategy.FCFS);
        Sink k1 = new Sink(m1, "K");
        OpenClass c1 = new OpenClass(m1, "C");
        s1.setArrival(c1, new MAP(mat(new double[][]{{-8, 1, 3}, {0, -6, 4}, {2, 0, -3}}),
                mat(new double[][]{{3, 1, 0}, {0, 2, 0}, {0, 0, 1}})));
        q1.setService(c1, new MAP(mat(new double[][]{{-3, 1}, {6, -7}}), mat(new double[][]{{0, 2}, {1, 0}})));
        m1.link(Network.serialRouting(c1, s1, q1, k1));
        assertTrue(Solver_mam_transient_qbd.applicable(m1.getStruct()));

        // M/M/1 (Poisson arrival + exp service) -> stays on the libQBD fast path.
        Network m2 = new Network("MM1");
        Source s2 = new Source(m2, "S");
        Queue q2 = new Queue(m2, "Q", SchedStrategy.FCFS);
        Sink k2 = new Sink(m2, "K");
        OpenClass c2 = new OpenClass(m2, "C");
        s2.setArrival(c2, new Exp(1.0));
        q2.setService(c2, new Exp(2.0));
        m2.link(Network.serialRouting(c2, s2, q2, k2));
        assertFalse(Solver_mam_transient_qbd.applicable(m2.getStruct()));
    }
}
