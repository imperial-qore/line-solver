% Example 2: A multiclass M/G/1 queue
GlobalConstants.setVerbose(VerboseLevel.STD);
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
tracePath = fullfile(cwd,'example_trace.txt');
queue.setService(jobclass2, Replayer(tracePath));

% Use direct serial routing for multiclass model
% model.link(Network.serialRouting(source,queue,sink));
P = model.initRoutingMatrix();
P.set(jobclass1, Network.serialRouting(source,queue,sink));
P.set(jobclass2, Network.serialRouting(source,queue,sink));
model.link(P);

jmtAvgTable = JMT(model,'seed',23000,'samples',10000).avgTable()

queue.setService(jobclass2, Replayer(tracePath).fitAPH());

% Use options struct to set nested config fields
ctmcOptions = CTMC.defaultOptions;
ctmcOptions.cutoff = 2;
ctmcOptions.verbose = true;
ctmcOptions.config.nonmkv = 'none';  % Disable automatic non-Markovian conversion
ctmcAvgTable2 = CTMC(model, ctmcOptions).avgTable()

ctmcOptions.cutoff = 4;
ctmcAvgTable4 = CTMC(model, ctmcOptions).avgTable()

mamOptions = MAM.defaultOptions;
mamOptions.config.nonmkv = 'none';
mamAvgTable = MAM(model, mamOptions).avgTable()
