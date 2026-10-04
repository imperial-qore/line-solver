function tranSysState = sampleSys(self, numEvents)
% TRANSYSSTATE = SAMPLESYS(NUMSAMPLES)


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    % A chain given directly to SolverCTMC is not a Network, so linemodel_save
    % has nothing to write for it; only the network case routes.
    if self.isChainSolver()
        CPPLINE.cppUnsupported(self.name, 'sampleSys', ...
            ['a chain given directly to SolverCTMC is not a Network, so there is no ' ...
            'model.json for line-cli to read']);
    end
    tranSysState = CPPLINE.sysSamplePath(self.name, self.model, self.options, numEvents, false);
    return
end

if self.isChainSolver()
    % Chain mode: a sample path of the user-supplied chain, started from
    % options.init_sol when given. A DTMC advances one unit of time per step.
    options = self.getOptions;
    Solver.resetRandomGeneratorSeed(options.seed);
    stateSpace = self.getStateSpace();
    n = size(stateSpace,1);
    pi0 = options.init_sol;
    if isempty(pi0)
        % Uniform start, as in the transient analysis, so that a sample path is
        % reproducible under options.seed.
        pi0 = ones(1,n)/n;
    else
        pi0 = reshape(pi0,1,[]);
    end
    tranSysState = struct();
    tranSysState.handle = {};
    if self.isDiscreteChain()
        sts = dtmc_simulate(self.chainModel.getTransMat(), pi0, numEvents);
        tranSysState.t = (0:(numEvents-1))';
    else
        [sjt, sts] = ctmc_simulate(self.chainModel.getGenerator(), pi0, numEvents);
        tranSysState.t = cumsum([0; reshape(sjt(1:end-1),[],1)]);
    end
    tranSysState.state = {stateSpace(sts(:),:)};
    tranSysState.event = {};
    tranSysState.isaggregate = false;
    return
end

self.assertPhaseTypeStates('sampleSys');

options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    runAnalyzer(self);
end
[infGen, eventFilt] = getGenerator(self);
[stateSpace, localStateSpace] = getStateSpace(self);
stateSpaceAggr = getStateSpaceAggr(self);

sn = self.getStruct;
initState = sn.state;
% see _kb/11-conventions-and-gotchas.md (CTMC state-vector padding) for rationale
spaceWidths = zeros(1,length(initState));
s0 = [];
for isf=1:length(initState)
    w = size(localStateSpace{isf},2);
    spaceWidths(isf) = w;
    s0 = [s0, zeros(1,w-length(initState{isf})), initState{isf}];
end
nst = cumsum([1,spaceWidths]);

% set initial state
pi0 = zeros(1,size(stateSpace,1));
pi0(matchrow(stateSpace,s0))=1;

% filter all CTMC events as a marked Markovian arrival process
D1 = cellsum(eventFilt);
D0 = infGen-D1;
MMAP = mmap_normalize([{D0},{D1},eventFilt(:)']);

% now sampel the MMAP
[sjt,event,~,~,sts] = mmap_sample(MMAP,numEvents, pi0);

sn = self.getStruct;
tranSysState = struct();
tranSysState.handle = self.model.getStatefulNodes';
tranSysState.t = cumsum([0,sjt(1:end-1)']');
for isf=1:length(initState)
    tranSysState.state{isf} = stateSpace(sts,(nst(isf):nst(isf+1)-1));
end

tranSysState.event = {};
for e = 1:length(event)    
    for a=1:length(sn.sync{event(e)}.active)
        tranSysState.event{end+1} = sn.sync{event(e)}.active{a};
        tranSysState.event{end}.t = tranSysState.t(e);
    end
    for p=1:length(sn.sync{event(e)}.passive)
        tranSysState.event{end+1} = sn.sync{event(e)}.passive{p};
        tranSysState.event{end}.t = tranSysState.t(e);
    end
end
tranSysState.isaggregate = false;

end