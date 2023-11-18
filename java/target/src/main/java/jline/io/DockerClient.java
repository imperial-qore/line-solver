package jline.io;

import jline.lang.layerednetworks.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Erlang;
import jline.lang.distributions.Exp;

import java.io.File;

/**
 * Examples of client for the Dockerized MATLAB service
 */
public class DockerClient {
    public static void main(String[] args) throws Exception {
        long s1 = System.nanoTime();
        LayeredNetwork model = new LayeredNetwork("LQN1");
        long t1 = System.nanoTime()-s1;

        // definition of processors, tasks and entries
        long s2 = System.nanoTime();
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.INF);
        long t2 = System.nanoTime()-s2;

        long s3 = System.nanoTime();
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        long t3 = System.nanoTime()-s3;


        long s4 = System.nanoTime();
        Entry E1 = new Entry(model, "E1");
        E1.on(T1);
        long t4 = System.nanoTime()-s4;

        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF);
        T2.on(P2);
        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        //definition of activities
        long s5 = System.nanoTime();
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001,2));
        long t5 = System.nanoTime()-s5;

        long s6 = System.nanoTime();
        Activity A1 = new Activity(model, "A1", new Exp(1.0));
        A1.on(T1);A1.boundTo(E1);A1.synchCall(E2,3);
        long t6 = System.nanoTime() - s6;
        Activity A2 = new Activity(model, "A2", new Exp(2.0));
        A2.on(T2);A2.boundTo(E2);
        A2.repliesTo(E2);


        File file = new File("");
        String filePath = file.getCanonicalPath();
        long s7 = System.nanoTime();
        model.writeXML(filePath+"/ex2.xml",true);
        long t7 = System.nanoTime() - s7;
        System.out.println(t1+" "+t2+" "+t3+" "+t4+" "+t5+" "+t6+" "+t7+" ");

        //model.writeXML("test.xml",false);
        model.sendModel("test.xml","5462");
    }
}
