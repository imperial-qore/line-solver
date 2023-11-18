## LINE: Queueing Theory Algorithms for MATLAB

Website: http://line-solver.sourceforge.net/

Latest stable release: https://sourceforge.net/projects/line-solver/files/latest/download

[![License](https://img.shields.io/badge/License-BSD%203--Clause-red.svg)](https://github.com/imperial-qore/line-solver/blob/master/LICENSE)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fimperial-qore%2Fline-solver&count_bg=%23FFC401&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
[![View LINE on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/71486-line)

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

A [getting started](https://github.com/imperial-qore/line-solver/wiki/Getting-started) walkthrough is available.


### License
LINE is released as open source under the [BSD-3 license](https://raw.githubusercontent.com/line-solver/line/master/matlab/LICENSE).

### Acknowledgement
The development of LINE has been partially funded by the European Commission grants FP7-318484 ([MODAClouds](http://multiclouddevops.com/)), H2020-644869 ([DICE](http://www.dice-h2020.eu/)), H2020-825040 ([RADON](http://radon-h2020.eu)), and by the EPSRC grant EP/M009211/1 ([OptiMAM](https://wp.doc.ic.ac.uk/optimam/)).

