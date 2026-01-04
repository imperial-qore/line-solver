package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Zipf;

public class SolverCTMCTestFixtures {
  public static Network test_tut03_repairmen() {
    Network model = new Network("MRP");
    Delay delay = new Delay(model, "WorkingState");
    Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
    queue.setNumberOfServers(2);
    ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
    delay.setService(closedClass, new Exp(0.5));
    queue.setService(closedClass, new Exp(4.0));
    model.link(model.serialRouting(delay, queue));
    return model;
  }

  public static Network test_tut05_completes_flag() {
      Network model = new Network("RL");
      Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
      ClosedClass jobClass1 = new ClosedClass(model, "Class1", 1, queue);
      ClosedClass jobClass2 = new ClosedClass(model, "Class2", 0, queue);
      ClosedClass jobClass3 = new ClosedClass(model, "Class3", 0, queue);
      queue.setService(jobClass1, Erlang.fitMeanAndOrder(1, 2));
      queue.setService(jobClass2, Erlang.fitMeanAndOrder(2, 2));
      queue.setService(jobClass3, Erlang.fitMeanAndOrder(3, 2));
      RoutingMatrix P = model.initRoutingMatrix();
      P.set(jobClass1, jobClass2, queue, queue, 1.0);
      P.set(jobClass2, jobClass3, queue, queue, 1.0);
      P.set(jobClass3, jobClass1, queue, queue, 1.0);
      model.link(P);
      //new SolverNC(model).getAvgTable().print();
      //new SolverNC(model).getAvgSysTable().print();
      jobClass1.setCompletes(false);
      jobClass2.setCompletes(false);
      //new SolverNC(model).getAvgSysTable().print();
      return model;
  }

  public static Network test_tut06_cache_lru_zipf() {
    Network model = new Network("Model");

    // Block 1: nodes
    Delay clientDelay = new Delay(model, "Client");
    Cache cacheNode = new Cache(model, "Cache", 5, 3, ReplacementStrategy.LRU);
    Delay cacheDelay = new Delay(model, "CacheDelay");

    // Block 2: classes
    ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
    ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
    ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

    clientDelay.setService(clientClass, Immediate.getInstance()); // (Client,ClientClass)
    cacheDelay.setService(hitClass, Exp.fitMean(0.2)); // (CacheDelay,HitClass)
    cacheDelay.setService(missClass, Exp.fitMean(1.0)); // (CacheDelay,MissClass)

    cacheNode.setRead(clientClass, new Zipf(1.4, 5));
    cacheNode.setHitClass(clientClass, hitClass);
    cacheNode.setMissClass(clientClass, missClass);

    // Block 3: topology
    RoutingMatrix P = model.initRoutingMatrix();

    P.set(clientClass, clientClass, clientDelay, cacheNode, 1.00); // (Client,ClientClass) -> (Cache,ClientClass)

    P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.00); // (Client,HitClass) -> (Cache,HitClass)
    P.set(missClass, missClass, cacheNode, cacheDelay, 1.00); // (Cache,MissClass) -> (CacheDelay,MissClass)

    P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.00); // (Cache,HitClass) -> (CacheDelay,HitClass)
    P.set(missClass, clientClass, cacheDelay, clientDelay, 1.00); // (Client,MissClass) -> (Cache,MissClass)

    model.link(P);
    //TODO: cache SSA
    //new SolverSSA(model,"samples",2e4,"seed",1,"verbose",true).getAvgTable().print();
    return model;
  }
}
