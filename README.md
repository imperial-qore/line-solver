## LINE Solver: Queueing Theory Algorithms 
Website: http://line-solver.sourceforge.net/

Latest stable release: https://sourceforge.net/projects/line-solver/files/latest/download

[![License](https://img.shields.io/badge/License-BSD%203--Clause-red.svg)](https://github.com/imperial-qore/line-solver/blob/master/LICENSE)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fimperial-qore%2Fline-solver&count_bg=%23FFC401&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

Main distribution of the LINE solver for [MATLAB](https://github.com/imperial-qore/line-solver/tree/main/matlab) (stable version), [Java](https://github.com/imperial-qore/line-solver/tree/main/java) (alpha version), and [Python](https://github.com/imperial-qore/line-solver/tree/main/python) (alpha version). Past releases can be found on [Sourceforge](https://sourceforge.net/projects/line-solver/files/).

### What is LINE?
LINE is an open source package to analyze queueing models via analytical methods and simulation. The tool features algorithms for the solution of open queueing systems (e.g., M/M/1, M/M/k, M/G/1, ...), open and closed queueing networks, and layered queueing networks. 

### Documentation
Check out the [LINE manual](https://github.com/imperial-qore/line-solver/tree/main/doc) and the README files in the [java/](https://github.com/imperial-qore/line-solver/tree/main/java), [matlab/](https://github.com/imperial-qore/line-solver/tree/main/matlab), and [python/](https://github.com/imperial-qore/line-solver/tree/main/python) folders for getting started information. 

### License
LINE is released as open source under the [BSD-3 license](https://raw.githubusercontent.com/imperial-qore/line-solver/main/matlab/LICENSE).

### Acknowledgement
The development of LINE has been partially funded by the European Commission grants FP7-318484 (MODAClouds), H2020-644869 (DICE), H2020-825040 (RADON), and by the EPSRC grant EP/M009211/1 (OptiMAM).

### Example: Solving a basic queueing system
We illustrate how to simulate an M/M/1 queue:

**MATLAB**: 
```
lineStart;
model = Network('M/M/1');

source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

jobclass = OpenClass(model, 'Class1');
source.setArrival(jobclass, Exp(1.0));
queue.setService(jobclass, Exp(2.0));

model.link(Network.serialRouting(source,queue,sink));

AvgTable = SolverJMT(model,'seed',23000).getAvgTable
```
**Java**: 
```
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.solvers.jmt.JMTOptions;
import jline.solvers.jmt.SolverJMT;

public class MM1 {
    public static void main(String[] args){
        Network model = new Network("M/M/1");
        
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(1.0)); 
        queue.setService(jobclass, new Exp(2.0)); 
        
        model.link(Network.serialRouting(source, queue, sink));
        
        new SolverJMT(model, "seed", 23000).getAvgTable().print();
    }
}
```
**Python**: 
```
from line_solver import *

if __name__ == "__main__":
    model = Network("M/M/1")

    source = Source(model, "Source")
    queue = Queue(model, "Queue", SchedStrategy.FCFS)
    sink = Sink(model, "Sink")
    
    jobclass = OpenClass(model, "Class1")
    source.setArrival(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))
   
    model.link(Network.serialRouting(source, queue, sink))
    
    avgTable = SolverJMT(model,'seed',23000).getAvgTable()
```
To extra a particular value you may use LINE's table get function (tget), for example: 
```
    print(tget(table,"Queue"))
    print(tget(table,"RespT","Queue"))
```
