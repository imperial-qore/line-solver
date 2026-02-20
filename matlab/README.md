## LINE Solver for MATLAB

This repository includes the MATLAB version of the [LINE solver](https://line-solver.sourceforge.net/). A MATLAB version of the [manual](https://line-solver.sourceforge.net/doc/LINE-matlab.pdf) includes getting started information.

### Getting started

To install LINE, first expand the archive (or clone the repository) in the chosen installation folder. Then start MATLAB, change the active directory to the installation folder and run:
```
lineInstall
```
To begin using LINE, add all LINE folders to the path using the following command:
```
lineStart
```
The lineStart command will be required at the beginning of every MATLAB session. You can now use LINE. 

For example, to solve a basic M/M/1 queue, type:
```
model = Network('M/M/1');
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');

oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, Exp(1.25));

model.link(Network.serialRouting(source,queue,sink));

AvgTable = SolverMVA(model).getAvgTable
```
This will provide the following output, showing for example an 80% utilization:
```
MVA analysis (method: default) completed. Runtime: 0.002546 seconds. Iterations: 1.
AvgTable =
  2×8 table
    Station     JobClass    QLen    Util    RespT    ResidT    ArvR    Tput
    ________    ________    ____    ____    _____    ______    ____    ____
    mySource    myClass      0        0       0        0        0       1  
    myQueue     myClass      4      0.8       4        4        1       1  
```    
Additional examples are available in the manual.

### License
LINE is released as open source under the [BSD-3 license](http://opensource.org/licenses/BSD-3-Clause).

### Acknowledgement
The development of LINE has been partially funded by the European Commission grants FP7-318484 ([MODAClouds](http://multiclouddevops.com/)), H2020-644869 ([DICE](http://www.dice-h2020.eu/)), H2020-825040 ([RADON](http://radon-h2020.eu)), and by the EPSRC grant EP/M009211/1 ([OptiMAM](https://wp.doc.ic.ac.uk/optimam/)).

