package jline.util;

import jline.lang.distributions.CumulativeDistribution;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

class CumulativeDistributionTest {
    CumulativeDistribution<Integer> testCumulativeDistribution;

    @BeforeEach
    void setUp() {
        testCumulativeDistribution = new CumulativeDistribution<Integer>(new Random());
    }

    @AfterEach
    void tearDown() {
    }

    @Test
    void addElement() {
        testCumulativeDistribution.addElement(61, 1);
        assertEquals((int) testCumulativeDistribution.sample(new Random()), 61);
    }

    @Test
    void generateEven() {
        for (int i = 0; i < 10; i++) {
            testCumulativeDistribution.addElement(i, 0.1);
        }

        int times_five = 0;
        for (int i = 0; i < 500; i++) {
            if (testCumulativeDistribution.sample(new Random()) == 5) {
                times_five++;
            }
        }

        assertTrue(times_five > 30);
        assertTrue(times_five < 70);
    }

    @Test
    void generateUnEven() {
        for (int i = 0; i <= 9; i++) {
            testCumulativeDistribution.addElement(i, ((double)i)/(45.0));
        }

        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
        for (int i = 0; i < 500; i++) {
            int n = testCumulativeDistribution.sample(new Random());
            if (counts.containsKey(n)) {
                counts.put(n, counts.get(n) + 1);
            } else {
                counts.put(n, 1);
            }
        }

        int testNum = (new Random()).nextInt(9)+1;
        int testCt = counts.get(testNum);
        int estCt = (int)((((double)testNum)/(45.0))*500.0);

        assertTrue(testCt > estCt-20);
        assertTrue(testCt < estCt+20);
    }

    @Test
    void generateRandom() {
        for (int i = 0; i < 100; i++) {
            Random r=  new Random();
            this.testCumulativeDistribution = new CumulativeDistribution<Integer>(r);

            int nMembers = r.nextInt(29)+1;
            List<Double> pdf = new ArrayList<Double>();
            double totalProb = 0;
            for (int j = 0; j < nMembers; j++) {
                double p = r.nextDouble();
                totalProb += p;
                pdf.add(p);
            }
            for (int j = 0; j < nMembers; j++) {
                pdf.set(j, pdf.get(j)/totalProb);
                this.testCumulativeDistribution.addElement(j, pdf.get(j));
            }

            List<Integer> counts = new ArrayList<Integer>();
        }
    }
}