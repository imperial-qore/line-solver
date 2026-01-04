package jline.io;

import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.VerboseLevel;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.solvers.mam.SolverMAM;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ssa.SolverSSA;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;

public class M2MTest {

    @TempDir
    Path tempDir;
  
    @BeforeEach
  public void matlabRandomSeedSetUp() {
    Maths.setRandomNumbersMatlab(true);
  }

  @AfterEach
  public void matlabRandomSeedClear() {
    Maths.setRandomNumbersMatlab(false);
  }

  @Test
  public void test_oqn_cs_routing_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/routing/oqn_cs_routing.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 100000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    
    // Verify that ArvR values are computed for Queue stations
    List<Double> arvRValues = avgTable.getArvR();
    List<String> stationNames = avgTable.getStationNames();
    
    // Check that Queue stations have non-zero arrival rates where expected
    for (int i = 0; i < stationNames.size(); i++) {
      if (stationNames.get(i).startsWith("Queue")) {
        // Queue stations should have non-zero arrival rates (except for classes with no arrivals)
        assertTrue(arvRValues.get(i) >= 0, "ArvR for " + stationNames.get(i) + " should be non-negative");
      } else if (stationNames.get(i).startsWith("Source")) {
        // Source stations should have zero arrival rates
        assertEquals(0.0, arvRValues.get(i), MID_TOL, "ArvR for Source should be zero");
      }
    }
    
    // Additionally verify that at least some Queue stations have positive arrival rates
    boolean hasPositiveQueueArvR = false;
    for (int i = 0; i < stationNames.size(); i++) {
      if (stationNames.get(i).startsWith("Queue") && arvRValues.get(i) > 0) {
        hasPositiveQueueArvR = true;
        break;
      }
    }
    assertTrue(hasPositiveQueueArvR, "At least one Queue station should have positive arrival rate");
  }

  @Test
  public void test_open_1class_1stat_mg1fcfs_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE(
                    "src/test/resources/jsimg/open/open_1class_1stat_mg1fcfs.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {1.60063498577707, 0};
    double[] expectedUtil = {0.480218372101896, 0};
    double[] expectedRespT = {1.62937393977393, 0};
    double[] expectedResidT = {1.62937393977393, 0};
    double[] expectedArvR = {1.02047837301846, 0};
    double[] expectedTput = {1.01939454021389, 1.02047837301846};
    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
                assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
                assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
                assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
                assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
                assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
                assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
              }    expectedQLen = new double[]{1.60063498577707, 0};
    expectedUtil = new double[]{0.480218372101896, 0};
    expectedRespT = new double[]{1.62937393977393, 0};
    expectedResidT = new double[]{1.62937393977393, 0};
    expectedArvR = new double[]{1.02047837301846, 0};
    expectedTput = new double[]{1.01939454021389, 1.02047837301846};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_open_1class_1stat_mg1ps_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/open/open_1class_1stat_mg1ps.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {1.02016952268552, 0};
    double[] expectedUtil = {0.503397183461647, 0};
    double[] expectedRespT = {0.950490831944015, 0};
    double[] expectedResidT = {0.950490831944015, 0};
    double[] expectedArvR = {1.03384956429008, 0};
    double[] expectedTput = {1.02026494745937, 1.03384956429008};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_open_1class_1stat_mm1k_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/open/open_1class_1stat_mm1k.jsimg");
    
    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.seed = 23000;
    options.samples = 50000;
    SolverJMT solver = new SolverJMT(model, options);
    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = new double[]{4.28940632367192, 0.0};
    double[] expectedUtil = new double[]{0.856218025631651, 0.0};
    double[] expectedRespT = new double[]{0.491719973644026, 0.0};
    double[] expectedResidT = new double[]{0.491719973644026, 0.0};
    double[] expectedArvR = new double[]{9.03190797919717, 0.0};
    double[] expectedTput = new double[]{8.66123703688577, 9.03190797919717};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_open_1class_1stat_mm1k_fcr_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE(
                    "src/test/resources/jsimg/open/open_1class_1stat_mm1k_fcr.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {5.6360, 0.0};
    double[] expectedUtil = {0.8910, 0.0};
    double[] expectedRespT = {0.6676, 0.0};
    double[] expectedResidT = {0.6676, 0.0};
    double[] expectedArvR = {8.9983, 0.0};
    double[] expectedTput = {8.9794, 8.9983};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }


//  @Test
// This JMT model is unstable in the first place, https://sourceforge.net/p/jmt/bugs/146/
//  public void test_open_1class_3stat_advfork_jsimg() {
//    ModelTransformer modelTransformer = new ModelTransformer();
//    Network model =
//            modelTransformer.JSIM2LINE(
//                    "src/test/java/jline/io/examples/open_1class_3stat_advfork.jsimg");
//
//    SolverOptions options = Solver.defaultOptions();
//    options.verbose = VerboseLevel.SILENT;
//    options.cutoff = Matrix.singleton(7);
//    options.seed = 23000;
//    options.samples = 10000;
//    SolverJMT solver = new SolverJMT(model, options);
//
//    NetworkAvgTable avgTable = solver.getAvgTable();
//    // Expected values
//    double[] expectedQLen = {0.0, 719.8878, 842.1135, 761.4572, 179.6867};
//    double[] expectedUtil = {0.0, 0.9444, 0.9368, 0.9374, 0.0};
//    double[] expectedRespT = {0.0, 4.9013e+03, 5.4148e+03, 5.0091e+03, 120.9071};
//    double[] expectedResidT = {0.0, 4.9013e+03, 5.4148e+03, 5.0091e+03, 362.7214};
//    double[] expectedArvR = {0.0, 0.2509, 0.2509, 0.2509, 0.5061};
//    double[] expectedTput = {0.0835, 0.1897, 0.1872, 0.1891, 0.0650};
//
//    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
//      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
//      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
//      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
//      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
//      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
//      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
//    }
//  }

//  @Test
//  public void test_open_1class_3stat_fork_jsimg() {
//    ModelTransformer modelTransformer = new ModelTransformer();
//    Network model =
//            modelTransformer.JSIM2LINE("src/test/java/jline/io/examples/open_1class_3stat_fork.jsimg");
//
//    SolverOptions options = Solver.defaultOptions();
//    options.verbose = VerboseLevel.SILENT;
//    options.cutoff = Matrix.singleton(7);
//    options.seed = 23000;
//    options.samples = 10000;
//    SolverJMT solver = new SolverJMT(model, options);
//
//    NetworkAvgTable avgTable = solver.getAvgTable();
//    // Expected values
//    double[] expectedQLen = {
//            0, 719.887792794558, 842.113532153842, 761.45716706121, 179.686721528542
//    };
//    double[] expectedUtil = {0, 0.944407981986922, 0.936810739057026, 0.93740099753255, 0};
//    double[] expectedRespT = {
//            0, 4901.28423530701, 5414.75218291496, 5009.06259329689, 120.907130545827
//    };
//    double[] expectedResidT = {
//            0, 4901.28423530701, 5414.75218291496, 5009.06259329689, 362.721391637481
//    };
//    double[] expectedArvR = {
//            0, 0.250889465466965, 0.250889465466965, 0.250889465466965, 0.506104987467522
//    };
//    double[] expectedTput = {
//            0.0834565298138268,
//            0.189699205440928,
//            0.187237431868481,
//            0.189111962264473,
//            0.0649751842983006
//    };
//
//    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
//      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
//      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
//      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
//      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
//      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
//      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
//    }
//  }

//  @Test
//  public void test_open_2class_3stat_fcr_jsimg() {
//    ModelTransformer modelTransformer = new ModelTransformer();
//    Network model =
//            modelTransformer.JSIM2LINE("src/test/java/jline/io/examples/open_2class_3stat_fcr.jsimg");
//
//    SolverOptions options = Solver.defaultOptions();
//    options.verbose = VerboseLevel.SILENT;
//    options.cutoff = Matrix.singleton(7);
//    options.seed = 23000;
//    options.samples = 10000;
//    SolverJMT solver = new SolverJMT(model, options);
//
//    NetworkAvgTable avgTable = solver.getAvgTable();
//    // Expected values
//    double[] expectedQLen = {
//            0,
//            0,
//            0,
//            81.3318138697168,
//            82.3144727644165,
//            0,
//            0,
//            966.492659701031,
//            0,
//            0.294293251196322,
//            0.366423134659199,
//            0
//    };
//    double[] expectedUtil = {
//            0,
//            0,
//            0,
//            0.587497219849474,
//            0.406917682098421,
//            0,
//            0,
//            0.995447620212352,
//            0,
//            0.294293251196322,
//            0.366423134659199,
//            0
//    };
//    double[] expectedRespT = {
//            0,
//            0,
//            0,
//            270.902131198466,
//            395.167728757376,
//            0,
//            0,
//            5173.43761572356,
//            0,
//            0.991777979874804,
//            0.996132902137542,
//            0
//    };
//    double[] expectedResidT = {
//            0,
//            0,
//            0,
//            270.902131198466,
//            395.167728757376,
//            0,
//            0,
//            5173.43761572356,
//            0,
//            1.98355595974961,
//            1.99226580427508,
//            0
//    };
//    double[] expectedArvR = {
//            0,
//            0,
//            0,
//            0.624252924244009,
//            0.434707677062558,
//            0,
//            0,
//            0.432847699132658,
//            0,
//            0.303106480455566,
//            0.360366671869717,
//            0
//    };
//    double[] expectedTput = {
//            0.313693936229892,
//            0.502128973505335,
//            0,
//            0.597173034178503,
//            0.416671626243924,
//            0.0004638516289053,
//            0.0196318121245519,
//            0.307609774660577,
//            0.000235591703365263,
//            0.296459185370978,
//            0.361400052368512,
//            0.00049555401380304
//    };
//
//    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
//      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
//      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
//      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
//      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
//      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
//      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
//    }
//  }

  @Test
  public void test_open_3class_3stat_loadbal_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE(
                    "src/test/resources/jsimg/open/open_3class_3stat_loadbal.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {
            0,
            0,
            0,
            0.451124217972402,
            0.774720481704585,
            0.977885095101402,
            0.917107746550702,
            0.401418330826574,
            0.994030022805306,
            1.81575802275549,
            0.51085237381917,
            0.910461504340658
    };
    double[] expectedUtil = {
            0,
            0,
            0,
            0.0510514224721988,
            0.369441511736153,
            0.345579223505567,
            0.303481201583139,
            0.10481066368623,
            0.33941945619938,
            0.504034464953145,
            0.248481538623085,
            0.0875655854398051
    };
    double[] expectedRespT = {
            0,
            0,
            0,
            2.12859186147664,
            2.27109540412058,
            3.01545748741143,
            3.08446251214976,
            1.0718222217919,
            3.09578468632334,
            3.61309546956878,
            2.05221145831205,
            2.72978185856044
    };
    double[] expectedResidT = {
            0,
            0,
            0,
            0.425718372295328,
            0.757031801373527,
            1.00515249580381,
            0.925338753644928,
            0.357274073930633,
            1.03192822877445,
            1.80654773478439,
            0.684070486104016,
            0.909927286186812
    };
    double[] expectedArvR = {
            0,
            0,
            0,
            0.201519443729802,
            0.373918009699397,
            0.340093555388765,
            0.302130310618254,
            0.397189106841162,
            0.335799701920208,
            0.507340429489571,
            0.24658822543575,
            0.327832719438403
    };
    double[] expectedTput = {
            1.01656949090911,
            1.01810312143833,
            1.004334710012,
            0.20168438462512,
            0.374487907223116,
            0.33995206957991,
            0.305782391035163,
            0.397174752530571,
            0.337922698743859,
            0.518950293718848,
            0.246661695241061,
            0.334586707291289
    };

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_open_3class_4stat_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/open/open_3class_4stat.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {
            0,
            0,
            0,
            2.16224826238001,
            1.42154029678101,
            1.88668377609624,
            0.552402817393378,
            0.374707193143169,
            0.482584320393151,
            4.23838713694232,
            2.8828551508072,
            3.30261164881393,
            1.36445139884974,
            0.718857029041157,
            1.08639872845597
    };
    double[] expectedUtil = {
            0,
            0,
            0,
            0.237921915719099,
            0.263544037276593,
            0.363958618192865,
            0.204936394024811,
            0.158117027342323,
            0.204799255682884,
            0.384308861454798,
            0.257661189449258,
            0.275299898240821,
            0.293936200624196,
            0.107410584425337,
            0.321224298626909
    };
    double[] expectedRespT = {
            0,
            0,
            0,
            2.62624689069238,
            2.74653029960433,
            2.83616690768689,
            3.0678122205049,
            3.16846710091424,
            3.37919071104924,
            24.5748936361156,
            25.1729653352128,
            24.7978253336472,
            6.74242666363641,
            5.39020709364954,
            7.76455448182952
    };
    double[] expectedResidT = {
            0,
            0,
            0,
            10.5049875627695,
            10.9861211984173,
            11.3446676307476,
            3.0678122205049,
            3.16846710091424,
            3.37919071104924,
            24.5748936361156,
            25.1729653352128,
            24.7978253336472,
            6.74242666363641,
            5.39020709364954,
            7.76455448182952
    };
    double[] expectedArvR = {
            0,
            0,
            0,
            0.802196122549262,
            0.500707752190375,
            0.605725881094821,
            0.198200825321537,
            0.123812069657006,
            0.138497438471086,
            0.201928466422825,
            0.12419574493127,
            0.135648296683338,
            0.197220294997684,
            0.122519556633639,
            0.143484493230681
    };
    double[] expectedTput = {
            0.197838562315363,
            0.123801666447271,
            0.143359367626595,
            0.812682562538064,
            0.499950569066624,
            0.605684493639118,
            0.193552704358862,
            0.121055850422637,
            0.138510915699188,
            0.199883161031852,
            0.124152422831554,
            0.135641854565324,
            0.196105702422049,
            0.122556094513403,
            0.141840236077065
    };

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_closed_1class_2stat_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/closed/closed_1class_2stat.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {0.993993723776265, 4.00600627622373};
    double[] expectedUtil = {0.993993723776265, 0.996541426583246};
    double[] expectedRespT = {0.992635618705519, 4.01029163257552};
    double[] expectedResidT = {0.992635618705519, 4.01029163257552};
    double[] expectedArvR = {1.00000220723877, 0.998623468705891};
    double[] expectedTput = {0.998623468705891, 1.00000220723877};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_closed_2class_2stat_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/closed/closed_2class_2stat.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {0.319638143118524, 0.341054155335189, 4.68036185688148, 9.65894584466481};
    double[] expectedUtil = {0.319638143118524, 0.341054155335189, 0.325863207100071, 0.674137222517124};
    double[] expectedRespT = {0.984096298570179, 0.986850885743574, 14.0319418298598, 28.927362172926};
    double[] expectedResidT = {0.984096298570179, 0.986850885743574, 14.0319418298598, 28.927362172926};
    double[] expectedArvR = {0.326523610092565, 0.336314914196545, 0.325240001586613, 0.336532509089059};
    double[] expectedTput = {0.325240001586613, 0.336312131209576, 0.326523610092565, 0.336314914196545};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  @Test
  public void test_closed_2class_2stat_switch_jsimg() {
    M2M modelTransformer = new M2M();
    Network model =
            modelTransformer.JSIM2LINE("src/test/resources/jsimg/closed/closed_2class_2stat_switch.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {0.374121470512981, 0.303007036866454, 5.36783693056382, 8.9310053645355};
    double[] expectedUtil = {0.374121470512981, 0.303007036866454, 0.366129250256917, 0.633873151888559};
    double[] expectedRespT = {0.991385535955688, 0.977514237953504, 14.2403197201617, 29.1885484720112};
    double[] expectedResidT = {0.540755746884921, 0.444324653615229, 7.76744712008821, 13.2675220327324};
    double[] expectedArvR = {0.373844096195864, 0.308132060808013, 0.376338133306122, 0.309474956694674};
    double[] expectedTput = {0.3771803390718, 0.30811801032872, 0.376721457846725, 0.308132060808013};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }

  /*
  // TODO: Has the same error in MATLAB
  @Test
  public void test_closed_2class_3stat_switch_jsimg() {
    ModelTransformer modelTransformer = new ModelTransformer();
    Network model =
            modelTransformer.JSIM2LINE("src/test/java/jline/io/examples/closed_2class_3stat_switch.jsimg");

    SolverOptions options = Solver.defaultOptions();
    options.verbose = VerboseLevel.SILENT;
    options.cutoff = Matrix.singleton(7);
    options.seed = 23000;
    options.samples = 10000;
    SolverJMT solver = new SolverJMT(model, options);

    NetworkAvgTable avgTable = solver.getAvgTable();
    // Expected values
    double[] expectedQLen = {0.374121470512981,0.303007036866454,5.36783693056382,8.9310053645355};
    double[] expectedUtil = {0.374121470512981,0.303007036866454,0.366129250256917,0.633873151888559};
    double[] expectedRespT = {0.991385535955688,0.977514237953504,14.2403197201617,29.1885484720112};
    double[] expectedResidT = {0.540755746884921,0.444324653615229,7.76744712008821,13.2675220327324};
    double[] expectedArvR = {0.373844096195864,0.308132060808013,0.376338133306122,0.309474956694674};
    double[] expectedTput = {0.3771803390718,0.30811801032872,0.376721457846725,0.308132060808013};

    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
    }
  }
   */

  //  @Test
  //  public void test_open_1class_1stat_mg1fcfs_jsimg() {
  //    ModelTransformer modelTransformer = new ModelTransformer();
  //    Network model =
  //
  // modelTransformer.JSIM2LINE("src/test/java/jline/io/examples/open_1class_1stat_mg1fcfs.jsimg");
  //
  //    SolverOptions options = Solver.defaultOptions();
  //    options.verbose = VerboseLevel.SILENT;
  //    options.cutoff = Matrix.singleton(7);
  //    options.seed = 23000;
  //    options.samples = 10000;
  //    SolverJMT solver = new SolverJMT(model, options);
  //
  //    NetworkAvgTable avgTable = solver.getAvgTable();
  //    // Expected values
  //    double[] expectedQLen = {1.6006, 0.0};
  //    double[] expectedUtil = {0.4802, 0.0};
  //    double[] expectedRespT = {1.6294, 0.0};
  //    double[] expectedResidT = {1.6294, 0.0};
  //    double[] expectedArvR = {1.0205, 0.0};
  //    double[] expectedTput = {1.0194, 1.0205};
  //
  //    for (int i = 0; i < avgTable.getClassNames().size(); i++) {
  //      assertEquals(expectedQLen[i], avgTable.get(0).get(i), MID_TOL);
  //      assertEquals(expectedUtil[i], avgTable.get(1).get(i), MID_TOL);
  //      assertEquals(expectedRespT[i], avgTable.get(2).get(i), MID_TOL);
  //      assertEquals(expectedResidT[i], avgTable.get(3).get(i), MID_TOL);
  //      assertEquals(expectedArvR[i], avgTable.get(4).get(i), MID_TOL);
  //      assertEquals(expectedTput[i], avgTable.get(5).get(i), MID_TOL);
  //    }
  //  }

    // ==============================================
    // LINE2JSIMG Tests (from LINE2JSIMTest.java)
    // ==============================================



    /**
     * Test round-trip conversion: LINE model -> JSIMG -> LINE model
     * Compares SolverMVA results before and after conversion to ensure they match
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_ex3() {
        // Step 1: Create the original tut03 model
        Network originalModel = M2MTestFixtures.closed_ex3();

        // Step 2: Solve original model with SolverMVA
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        options.cutoff = Matrix.singleton(7);
        options.seed = 23000;

        SolverMVA originalSolver = new SolverMVA(originalModel, options);
        NetworkAvgTable originalAvgTable = originalSolver.getAvgTable();
        originalAvgTable.print();

        // Step 3: Convert original model to JSIMG using LINE2JSIMG
        M2M m2m = new M2M();
        File jsimgFile = tempDir.resolve("tut03_test.jsimg").toFile();

        boolean saveSuccess = m2m.LINE2JSIMG(originalModel, jsimgFile.getAbsolutePath());
        assertTrue(saveSuccess, "Failed to save model to JSIMG format");
        assertTrue(jsimgFile.exists(), "JSIMG file was not created");
        
        // Debug: Print file size and content preview for debugging
        if (jsimgFile.length() == 0) {
            fail("Generated JSIMG file is empty");
        }
        
        // Print the JSIMG file contents to screen
        // try {
        //     System.out.println("=== JSIMG File Contents ===");
        //     java.nio.file.Files.lines(jsimgFile.toPath()).forEach(System.out::println);
        //     System.out.println("=== End JSIMG File Contents ===");
        // } catch (java.io.IOException e) {
        //     System.out.println("Error reading JSIMG file: " + e.getMessage());
        // }

        // Step 4: Load the JSIMG file back using JSIM2LINE
        Network loadedModel = null;
        try {
            loadedModel = m2m.JSIM2LINE(jsimgFile.getAbsolutePath());
        } catch (Exception e) {
            fail("Exception during JSIM2LINE: " + e.getMessage());
        }
        assertNotNull(loadedModel, "Failed to load model from JSIMG file");

        // Step 5: Solve loaded model with SolverMVA
        SolverMVA loadedSolver = new SolverMVA(loadedModel, options);
        NetworkAvgTable loadedAvgTable = loadedSolver.getAvgTable();

        // Step 6: Compare the results - they should match within tolerance
        assertNotNull(originalAvgTable, "Original avg table should not be null");
        assertNotNull(loadedAvgTable, "Loaded avg table should not be null");

        // Compare dimensions
        assertEquals(originalAvgTable.getStationNames().size(), loadedAvgTable.getStationNames().size(),
                    "Number of stations should match");
        assertEquals(originalAvgTable.getClassNames().size(), loadedAvgTable.getClassNames().size(),
                    "Number of classes should match");

        // Compare all metrics: QLen, Util, RespT, ResidT, ArvR, Tput
        // The metrics are indexed as: 0=QLen, 1=Util, 2=RespT, 3=ResidT, 4=ArvR, 5=Tput
        for (int metricIndex = 0; metricIndex < 6; metricIndex++) {
            List<Double> originalMetric = originalAvgTable.get(metricIndex);
            List<Double> loadedMetric = loadedAvgTable.get(metricIndex);

            assertEquals(originalMetric.size(), loadedMetric.size(),
                        String.format("Size mismatch for metric %s", getMetricName(metricIndex)));

            for (int valueIndex = 0; valueIndex < originalMetric.size(); valueIndex++) {
                double originalValue = originalMetric.get(valueIndex);
                double loadedValue = loadedMetric.get(valueIndex);

                String metricName = getMetricName(metricIndex);
                int stationIndex = valueIndex / originalAvgTable.getClassNames().size();
                int classIndex = valueIndex % originalAvgTable.getClassNames().size();
                String stationName = stationIndex < originalAvgTable.getStationNames().size() ?
                                   originalAvgTable.getStationNames().get(stationIndex) : "Unknown";
                String className = classIndex < originalAvgTable.getClassNames().size() ?
                                 originalAvgTable.getClassNames().get(classIndex) : "Unknown";

                assertEquals(originalValue, loadedValue, MID_TOL,
                    String.format("Mismatch in %s for station %s, class %s: original=%.6f, loaded=%.6f",
                                metricName, stationName, className, originalValue, loadedValue));
            }
        }

        // Clean up
        if (jsimgFile.exists()) {
            jsimgFile.delete();
        }
    }

    /**
     * Helper method to get metric name by index
     */
    private String getMetricName(int index) {
        switch (index) {
            case 0: return "QLen";
            case 1: return "Util";
            case 2: return "RespT";
            case 3: return "ResidT";
            case 4: return "ArvR";
            case 5: return "Tput";
            default: return "Unknown";
        }
    }

    @Test
    public void test_LINE2JSIMG_model_structure_preservation() {
        // Create original model
        Network originalModel = M2MTestFixtures.closed_ex3();

        // Convert to JSIMG and back
        M2M m2m = new M2M();
        File jsimgFile = tempDir.resolve("tut03_structure_test.jsimg").toFile();
        
        assertTrue(m2m.LINE2JSIMG(originalModel, jsimgFile.getAbsolutePath()), 
                  "Failed to save model to JSIMG");
        
        Network loadedModel = m2m.JSIM2LINE(jsimgFile.getAbsolutePath());
        assertNotNull(loadedModel, "Failed to load model from JSIMG");

        // Check basic structure
        assertEquals(originalModel.getNodes().size(), loadedModel.getNodes().size(), 
                    "Number of nodes should match");
        assertEquals(originalModel.getClasses().size(), loadedModel.getClasses().size(), 
                    "Number of classes should match");

        // Check node names are preserved (closed_ex3 has 8 nodes: Delay1, Queue1, Delay2, Queue2, Delay3, Queue3, Delay4, Queue4)
        assertEquals("Delay1", loadedModel.getNodes().get(0).getName(), 
                    "First node name should be preserved");
        assertEquals("Queue1", loadedModel.getNodes().get(1).getName(), 
                    "Second node name should be preserved");
        assertEquals("Delay2", loadedModel.getNodes().get(2).getName(), 
                    "Third node name should be preserved");
        assertEquals("Queue2", loadedModel.getNodes().get(3).getName(), 
                    "Fourth node name should be preserved");
        assertEquals("Delay3", loadedModel.getNodes().get(4).getName(), 
                    "Fifth node name should be preserved");
        assertEquals("Queue3", loadedModel.getNodes().get(5).getName(), 
                    "Sixth node name should be preserved");
        assertEquals("Delay4", loadedModel.getNodes().get(6).getName(), 
                    "Seventh node name should be preserved");
        assertEquals("Queue4", loadedModel.getNodes().get(7).getName(), 
                    "Eighth node name should be preserved");

        // Check class properties (closed_ex3 has Class1 with population 32)
        assertEquals("Class1", loadedModel.getClasses().get(0).getName(), 
                    "Class name should be preserved");
        assertTrue(loadedModel.getClasses().get(0) instanceof ClosedClass, 
                  "Class type should be preserved as ClosedClass");

        ClosedClass loadedClass = (ClosedClass) loadedModel.getClasses().get(0);
        assertEquals(32, loadedClass.getPopulation(), 
                    "Class population should be preserved");

        // Clean up
        if (jsimgFile.exists()) {
            jsimgFile.delete();
        }
    }

    // ==============================================
    // Comprehensive Tests for all M2MTestFixtures models
    // ==============================================

    /**
     * Test roundtrip conversion for closed_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex1() {
        testRoundtripConversion("closed_ex1", M2MTestFixtures.closed_ex1());
    }

    /**
     * Test roundtrip conversion for closed_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex2() {
        testRoundtripConversion("closed_ex2", M2MTestFixtures.closed_ex2());
    }

    /**
     * Test roundtrip conversion for closed_ex4 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex4() {
        testRoundtripConversion("closed_ex4", M2MTestFixtures.closed_ex4());
    }

    /**
     * Test roundtrip conversion for closed_ex5 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex5() {
        testRoundtripConversion("closed_ex5", M2MTestFixtures.closed_ex5());
    }

    /**
     * Test roundtrip conversion for closed_ex6 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex6() {
        testRoundtripConversion("closed_ex6", M2MTestFixtures.closed_ex6());
    }

    /**
     * Test roundtrip conversion for closed_ex7 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex7() {
        testRoundtripConversion("closed_ex7", M2MTestFixtures.closed_ex7());
    }

    /**
     * Test roundtrip conversion for closed_ex8 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex8() {
        testRoundtripConversion("closed_ex8", M2MTestFixtures.closed_ex8());
    }

    /**
     * Test roundtrip conversion for closed_ex9 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_closed_ex9() {
        testRoundtripConversion("closed_ex9", M2MTestFixtures.closed_ex9());
    }

    /**
     * Test roundtrip conversion for open_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex1() {
        testRoundtripConversion("open_ex1", M2MTestFixtures.open_ex1());
    }

    /**
     * Test roundtrip conversion for open_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex2() {
        testRoundtripConversion("open_ex2", M2MTestFixtures.open_ex2());
    }

    /**
     * Test roundtrip conversion for open_ex3 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex3() {
        testRoundtripConversion("open_ex3", M2MTestFixtures.open_ex3());
    }

    /**
     * Test roundtrip conversion for open_ex4 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex4() {
        testRoundtripConversion("open_ex4", M2MTestFixtures.open_ex4());
    }

    /**
     * Test roundtrip conversion for open_ex5 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex5() {
        testRoundtripConversion("open_ex5", M2MTestFixtures.open_ex5());
    }

    /**
     * Test roundtrip conversion for open_ex6 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex6() {
        testRoundtripConversion("open_ex6", M2MTestFixtures.open_ex6());
    }

    /**
     * Test roundtrip conversion for open_ex7 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex7() {
        testRoundtripConversion("open_ex7", M2MTestFixtures.open_ex7());
    }

    /**
     * Test roundtrip conversion for open_ex8 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex8() {
        testRoundtripConversion("open_ex8", M2MTestFixtures.open_ex8());
    }

    /**
     * Test roundtrip conversion for open_ex9 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex9() {
        testRoundtripConversion("open_ex9", M2MTestFixtures.open_ex9());
    }

    /**
     * Test roundtrip conversion for open_ex10 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_open_ex10() {
        testRoundtripConversion("open_ex10", M2MTestFixtures.open_ex10());
    }

    /**
     * Test roundtrip conversion for other_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_other_ex1() {
        testRoundtripConversion("other_ex1", M2MTestFixtures.other_ex1());
    }

    /**
     * Test roundtrip conversion for other_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_other_ex2() {
        testRoundtripConversion("other_ex2", M2MTestFixtures.other_ex2());
    }

    /**
     * Test roundtrip conversion for ruuskanenExample2 model
     * Note: Disabled due to JSIM2LINE parsing issues ("Outside of matrix bounds")
     */
    @Test
    public void test_LINE2JSIMG_roundtrip_ruuskanenExample2() {
        testRoundtripConversion("ruuskanenExample2", M2MTestFixtures.ruuskanenExample2());
    }

    /**
     * Generic helper method to test roundtrip conversion for any model
     */
    private void testRoundtripConversion(String modelName, Network originalModel) {
        // Step 1: Solve original model with SolverMVA
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = 23000;

        SolverMVA originalSolver = new SolverMVA(originalModel, options);
        NetworkAvgTable originalAvgTable = originalSolver.getAvgTable();

        // Step 2: Convert original model to JSIMG using LINE2JSIMG
        M2M m2m = new M2M();
        File jsimgFile = tempDir.resolve(modelName + "_test.jsimg").toFile();
        boolean saveSuccess = m2m.LINE2JSIMG(originalModel, jsimgFile.getAbsolutePath());
        assertTrue(saveSuccess, "Failed to save " + modelName + " model to JSIMG format");
        assertTrue(jsimgFile.exists(), "JSIMG file was not created for " + modelName);

        // Step 3: Load the JSIMG file back using JSIM2LINE
        Network loadedModel = null;
        try {
            loadedModel = m2m.JSIM2LINE(jsimgFile.getAbsolutePath());
        } catch (Exception e) {
            fail("Exception during JSIM2LINE for " + modelName + ": " + e.getMessage());
        }
        assertNotNull(loadedModel, "Failed to load " + modelName + " model from JSIMG file");

        // Step 4: Solve loaded model with SolverMVA
        SolverMVA loadedSolver = new SolverMVA(loadedModel, options);
        NetworkAvgTable loadedAvgTable = loadedSolver.getAvgTable();

        // Step 5: Compare the results - they should match within tolerance
        assertNotNull(originalAvgTable, modelName + " original avg table should not be null");
        assertNotNull(loadedAvgTable, modelName + " loaded avg table should not be null");

        // Compare dimensions
        assertEquals(originalAvgTable.getStationNames().size(), loadedAvgTable.getStationNames().size(),
                modelName + ": Number of stations should match");
        assertEquals(originalAvgTable.getClassNames().size(), loadedAvgTable.getClassNames().size(),
                modelName + ": Number of classes should match");

        // Compare all metrics: QLen, Util, RespT, ResidT, ArvR, Tput
        for (int metricIndex = 0; metricIndex < 6; metricIndex++) {
            List<Double> originalMetric = originalAvgTable.get(metricIndex);
            List<Double> loadedMetric = loadedAvgTable.get(metricIndex);

            assertEquals(originalMetric.size(), loadedMetric.size(),
                    String.format("%s: Size mismatch for metric %s", modelName, getMetricName(metricIndex)));

            for (int valueIndex = 0; valueIndex < originalMetric.size(); valueIndex++) {
                double originalValue = originalMetric.get(valueIndex);
                double loadedValue = loadedMetric.get(valueIndex);

                String metricName = getMetricName(metricIndex);
                int stationIndex = valueIndex / originalAvgTable.getClassNames().size();
                int classIndex = valueIndex % originalAvgTable.getClassNames().size();
                String stationName = stationIndex < originalAvgTable.getStationNames().size() ?
                        originalAvgTable.getStationNames().get(stationIndex) : "Unknown";
                String className = classIndex < originalAvgTable.getClassNames().size() ?
                        originalAvgTable.getClassNames().get(classIndex) : "Unknown";

                assertEquals(originalValue, loadedValue, MID_TOL,
                        String.format("%s: Mismatch in %s for station %s, class %s: original=%.6f, loaded=%.6f",
                                modelName, metricName, stationName, className, originalValue, loadedValue));
            }
        }

        // Clean up
        if (jsimgFile.exists()) {
            jsimgFile.delete();
        }
    }

    // ==============================================
    // Structure Preservation Tests for all M2MTestFixtures models
    // ==============================================

    /**
     * Test structure preservation for closed_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex1() {
        testStructurePreservation("closed_ex1", M2MTestFixtures.closed_ex1());
    }

    /**
     * Test structure preservation for closed_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex2() {
        testStructurePreservation("closed_ex2", M2MTestFixtures.closed_ex2());
    }

    /**
     * Test structure preservation for closed_ex4 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex4() {
        testStructurePreservation("closed_ex4", M2MTestFixtures.closed_ex4());
    }

    /**
     * Test structure preservation for closed_ex5 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex5() {
        testStructurePreservation("closed_ex5", M2MTestFixtures.closed_ex5());
    }

    /**
     * Test structure preservation for closed_ex6 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex6() {
        testStructurePreservation("closed_ex6", M2MTestFixtures.closed_ex6());
    }

    /**
     * Test structure preservation for closed_ex7 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex7() {
        testStructurePreservation("closed_ex7", M2MTestFixtures.closed_ex7());
    }

    /**
     * Test structure preservation for closed_ex8 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex8() {
        testStructurePreservation("closed_ex8", M2MTestFixtures.closed_ex8());
    }

    /**
     * Test structure preservation for closed_ex9 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_closed_ex9() {
        testStructurePreservation("closed_ex9", M2MTestFixtures.closed_ex9());
    }

    /**
     * Test structure preservation for open_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex1() {
        testStructurePreservation("open_ex1", M2MTestFixtures.open_ex1());
    }

    /**
     * Test structure preservation for open_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex2() {
        testStructurePreservation("open_ex2", M2MTestFixtures.open_ex2());
    }

    /**
     * Test structure preservation for open_ex3 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex3() {
        testStructurePreservation("open_ex3", M2MTestFixtures.open_ex3());
    }

    /**
     * Test structure preservation for open_ex4 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex4() {
        testStructurePreservation("open_ex4", M2MTestFixtures.open_ex4());
    }

    /**
     * Test structure preservation for open_ex5 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex5() {
        testStructurePreservation("open_ex5", M2MTestFixtures.open_ex5());
    }

    /**
     * Test structure preservation for open_ex6 model
     * Note: Disabled due to solver inability to compute results
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex6() {
        testStructurePreservation("open_ex6", M2MTestFixtures.open_ex6());
    }

    /**
     * Test structure preservation for open_ex7 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex7() {
        testStructurePreservation("open_ex7", M2MTestFixtures.open_ex7());
    }

    /**
     * Test structure preservation for open_ex8 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex8() {
        testStructurePreservation("open_ex8", M2MTestFixtures.open_ex8());
    }

    /**
     * Test structure preservation for open_ex9 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex9() {
        testStructurePreservation("open_ex9", M2MTestFixtures.open_ex9());
    }

    /**
     * Test structure preservation for open_ex10 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_open_ex10() {
        testStructurePreservation("open_ex10", M2MTestFixtures.open_ex10());
    }

    /**
     * Test structure preservation for other_ex1 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_other_ex1() {
        testStructurePreservation("other_ex1", M2MTestFixtures.other_ex1());
    }

    /**
     * Test structure preservation for other_ex2 model
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_other_ex2() {
        testStructurePreservation("other_ex2", M2MTestFixtures.other_ex2());
    }

    /**
     * Test structure preservation for ruuskanenExample2 model
     * Note: Disabled due to JSIM2LINE parsing issues ("Outside of matrix bounds")
     */
    @Test
    public void test_LINE2JSIMG_structure_preservation_ruuskanenExample2() {
        testStructurePreservation("ruuskanenExample2", M2MTestFixtures.ruuskanenExample2());
    }

    /**
     * Generic helper method to test structure preservation for any model
     */
    private void testStructurePreservation(String modelName, Network originalModel) {
        // Convert to JSIMG and back
        M2M m2m = new M2M();
        File jsimgFile = tempDir.resolve(modelName + "_structure_test.jsimg").toFile();

        assertTrue(m2m.LINE2JSIMG(originalModel, jsimgFile.getAbsolutePath()),
                "Failed to save " + modelName + " model to JSIMG");

        Network loadedModel = m2m.JSIM2LINE(jsimgFile.getAbsolutePath());
        assertNotNull(loadedModel, "Failed to load " + modelName + " model from JSIMG");

        // Check basic structure
        assertEquals(originalModel.getNodes().size(), loadedModel.getNodes().size(),
                modelName + ": Number of nodes should match");
        assertEquals(originalModel.getClasses().size(), loadedModel.getClasses().size(),
                modelName + ": Number of classes should match");

        // Check node names are preserved
        for (int i = 0; i < originalModel.getNodes().size(); i++) {
            assertEquals(originalModel.getNodes().get(i).getName(), loadedModel.getNodes().get(i).getName(),
                    String.format("%s: Node %d name should be preserved", modelName, i));
        }

        // Check class names and types are preserved
        for (int i = 0; i < originalModel.getClasses().size(); i++) {
            assertEquals(originalModel.getClasses().get(i).getName(), loadedModel.getClasses().get(i).getName(),
                    String.format("%s: Class %d name should be preserved", modelName, i));

            assertEquals(originalModel.getClasses().get(i).getClass(), loadedModel.getClasses().get(i).getClass(),
                    String.format("%s: Class %d type should be preserved", modelName, i));

            // Check population for closed classes
            if (originalModel.getClasses().get(i) instanceof ClosedClass) {
                ClosedClass originalClass = (ClosedClass) originalModel.getClasses().get(i);
                ClosedClass loadedClass = (ClosedClass) loadedModel.getClasses().get(i);
                assertEquals(originalClass.getPopulation(), loadedClass.getPopulation(),
                        String.format("%s: Class %d population should be preserved", modelName, i));
            }
        }

        // Clean up
        if (jsimgFile.exists()) {
            jsimgFile.delete();
        }
    }
}
