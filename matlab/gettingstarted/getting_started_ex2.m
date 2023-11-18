% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-2-a-multiclass-mg1-queue
cwd = fileparts(mfilename('fullpath'));

model = Network('M/G/1');
source = Source(model,'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model,'Sink');

jobclass1 = OpenClass(model, 'Class1');
jobclass2 = OpenClass(model, 'Class2');

source.setArrival(jobclass1, Exp(0.5));
source.setArrival(jobclass2, Exp(0.5));

queue.setService(jobclass1, Erlang.fitMeanAndSCV(1,1/3));
queue.setService(jobclass2, Replayer([cwd,filesep,'example_trace.txt']));

P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue,sink);
P{2} = Network.serialRouting(source,queue,sink);
model.link(P);

jmtAvgTable = SolverJMT(model,'seed',23000).getAvgTable

queue.setService(jobclass2, Replayer([cwd,filesep,'example_trace.txt']).fitAPH());
ctmcAvgTable2 = SolverCTMC(model,'cutoff',2,'verbose',true).getAvgTable
ctmcAvgTable4 = SolverCTMC(model,'cutoff',4,'verbose',true).getAvgTable
mamAvgTable = SolverMAM(model).getAvgTable
