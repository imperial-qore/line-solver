% MMAP-fed small RR cache with two correlated classes.
%
% A marked MMPP2 arrival stream feeds a small Round-Robin (RR) cache. Its two
% marks are bound to two open read classes that share the modulating chain, so
% the classes are cross-correlated and autocorrelated in time. Each class reads
% the cache with a DIFFERENT item-popularity distribution.
%
% Phase 1 (bursty) emits mostly class-1 references at a high rate; phase 2
% (calm) emits mostly class-2 references at a low rate. The shared modulating
% chain therefore couples "which class arrives" with "how fast requests arrive".
%
% We solve the same system three ways:
%   (1) LDES  - discrete-event simulation of the true MMAP-fed cache;
%   (2) CTMC  - exact continuous-time Markov chain of the true system;
%   (3) ENV   - three random-environment methods that view the MMPP2 phase as
%               an environment modulating phase-conditional Poisson arrivals
%               (D1 diagonal):
%                 'avg'   - fast-environment limit: replace the modulation by its
%                           time-average per-class Poisson rates and solve ONE
%                           rate-averaged cache;
%                 'dec'   - slow-environment limit: quasi-stationary
%                           decomposition, solve each phase independently and
%                           average per-phase hit ratios weighted by phase prob;
%                 'blend' - state-vector coupling: carry the cache-state
%                           distribution across phase switches and average each
%                           phase's sojourn-weighted distribution. For a
%                           Markovian environment this recovers the exact joint
%                           (cache x phase) solution, i.e. it matches CTMC.
% The 'avg' and 'dec' limits discard the within-phase temporal correlation of
% references, so neither is guaranteed to bracket the exact value.

clear solver AvgTable;

n = 4; % number of items
m = 2; % cache capacity

% Per-class item-popularity distributions (deliberately different)
pAccess1 = DiscreteSampler([8 4 2 1]/15); % class 1 favors low-index items
pAccess2 = DiscreteSampler([1 2 4 8]/15); % class 2 favors high-index items

% Marked MMPP2 (M3A layout D = {D0, D11, D12}); D1 = D11 + D12 diagonal.
% Phase 1 bursty (rate 4, 90% class1); phase 2 calm (rate 1, 80% class2).
% Off-diagonal of D0 are the phase-switch rates (both 0.5).
D0  = [-4.5, 0.5; 0.5, -1.5];
D11 = [3.6, 0; 0, 0.2];   % class-1 arrivals per phase
D12 = [0.4, 0; 0, 0.8];   % class-2 arrivals per phase
mmap = MarkedMAP({D0, D11, D12}, 2);

%% (1) LDES - simulation of the true MMAP-fed cache
trueModel = buildCacheModel(n, m, pAccess1, pAccess2);
src = trueModel.getNodeByName('Source');
src.setMarkedArrival(mmap, {trueModel.classes{1}, trueModel.classes{2}});
cacheTrue = trueModel.getNodeByName('Cache');

solver{1} = LDES(trueModel, 'samples', 2e5, 'seed', 23000, 'verbose', true);
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}
hitLDES = cacheTrue.getHitRatio;

%% (2) CTMC - exact solution of the true system
trueModel.reset;
solver{2} = CTMC(trueModel, 'keep', false, 'cutoff', 1);
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}
hitCTMC = cacheTrue.getHitRatio;

%% Random environment: the MMPP2 phase modulates phase-conditional Poisson
% arrivals (D1 diagonal); the environment switches at the MMPP2 phase-transition
% rates (-D0 diagonal minus the total arrival rate).
envBase = buildCacheModel(n, m, pAccess1, pAccess2);
env = Environment('MMPPphase');
env.addStage('Phase1', 'bursty', local_setRates(envBase, D11(1,1), D12(1,1))); % 3.6/0.4
env.addStage('Phase2', 'calm',   local_setRates(envBase, D11(2,2), D12(2,2))); % 0.2/0.8
env.addTransition('Phase1', 'Phase2', Exp(-D0(1,1) - (D11(1,1)+D12(1,1)))); % 0.5
env.addTransition('Phase2', 'Phase1', Exp(-D0(2,2) - (D11(2,2)+D12(2,2)))); % 0.5
env.init();
env.getStageTable()

solverFactory = @(mdl) CTMC(mdl, 'keep', false, 'cutoff', 1);

% (3a) ENV method 'avg' - fast-environment limit (rate-averaged single model)
optAvg = Solver.defaultOptions; optAvg.method = 'avg'; optAvg.verbose = false;
solver{3} = ENV(env, solverFactory, optAvg);
solver{3}.getAvg();
hitAVG = solver{3}.ensemble{1}.getNodeByName('Cache').getHitRatio;

% (3b) ENV method 'dec' - slow-environment quasi-stationary decomposition
optDec = Solver.defaultOptions; optDec.method = 'dec'; optDec.verbose = false;
solver{4} = ENV(env, solverFactory, optDec);
solver{4}.getAvg();
hitDEC = solver{4}.ensemble{1}.getNodeByName('Cache').getHitRatio;

% (3c) ENV method 'blend' - state-vector coupling: carries the cache-state
% distribution across phase switches and averages each phase's sojourn-weighted
% distribution. Needs a finite-timespan CTMC inner solver.
blendFactory = @(mdl) CTMC(mdl, 'keep', false, 'cutoff', 1, 'timespan', [0, 1e3]);
optBlend = Solver.defaultOptions; optBlend.method = 'blend'; optBlend.verbose = false;
optBlend.iter_max = 100; optBlend.iter_tol = 1e-4;
solver{5} = ENV(env, blendFactory, optBlend);
solver{5}.getAvg();
hitBLEND = solver{5}.ensemble{1}.getNodeByName('Cache').getHitRatio;

% (3d) ENV default mean-field with an FLD (refined mean-field) inner solver.
% The cache is analyzed by the RMF drift; the mean occupancy is carried across
% phase switches (the cache analog of the queue-length handoff) and the hit
% ratio is the probEnv-weighted, sojourn-averaged (arrival x hit-prob). Needs a
% finite-timespan fluid inner solver. This is the fully mean-field counterpart
% of 'blend', trading the exact joint distribution for a fluid approximation.
fldFactory = @(mdl) FLD(mdl, 'method', 'rmf', 'timespan', [0, 50]);
optFLD = Solver.defaultOptions; optFLD.verbose = false;
optFLD.iter_max = 100; optFLD.iter_tol = 1e-4;
solver{6} = ENV(env, fldFactory, optFLD);
solver{6}.getAvg();
hitFLD = solver{6}.ensemble{1}.getNodeByName('Cache').getHitRatio;

%% Summary: per-read-class actual hit ratio (classes 1=Read1, 2=Read2)
fprintf('\n--- Actual cache hit ratio per read class ---\n');
fprintf('                 Read1      Read2\n');
fprintf('LDES (sim)   : %8.4f  %8.4f\n', hitLDES(1),  hitLDES(2));
fprintf('CTMC (true)  : %8.4f  %8.4f\n', hitCTMC(1),  hitCTMC(2));
fprintf('ENV (avg)    : %8.4f  %8.4f\n', hitAVG(1),   hitAVG(2));
fprintf('ENV (dec)    : %8.4f  %8.4f\n', hitDEC(1),   hitDEC(2));
fprintf('ENV (blend)  : %8.4f  %8.4f\n', hitBLEND(1), hitBLEND(2));
fprintf('ENV (mf/FLD) : %8.4f  %8.4f\n', hitFLD(1),   hitFLD(2));

%% Local functions
function model = buildCacheModel(n, m, pAccess1, pAccess2)
    model = Network('MMAPCache');
    source    = Source(model, 'Source');
    cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.RR);
    sink      = Sink(model, 'Sink');

    rd1  = OpenClass(model, 'Read1', 0);
    rd2  = OpenClass(model, 'Read2', 0);
    hit1 = OpenClass(model, 'Hit1', 0);
    mis1 = OpenClass(model, 'Miss1', 0);
    hit2 = OpenClass(model, 'Hit2', 0);
    mis2 = OpenClass(model, 'Miss2', 0);

    cacheNode.setRead(rd1, pAccess1);
    cacheNode.setRead(rd2, pAccess2);
    cacheNode.setHitClass(rd1, hit1);  cacheNode.setMissClass(rd1, mis1);
    cacheNode.setHitClass(rd2, hit2);  cacheNode.setMissClass(rd2, mis2);

    P = model.initRoutingMatrix;
    P{rd1,rd1}(source, cacheNode)  = 1.0;
    P{rd2,rd2}(source, cacheNode)  = 1.0;
    P{hit1,hit1}(cacheNode, sink)  = 1.0;
    P{mis1,mis1}(cacheNode, sink)  = 1.0;
    P{hit2,hit2}(cacheNode, sink)  = 1.0;
    P{mis2,mis2}(cacheNode, sink)  = 1.0;
    model.link(P);
end

function model = local_setRates(baseModel, lambda1, lambda2)
    model  = baseModel.copy();
    source = model.getNodeByName('Source');
    source.setArrival(model.classes{1}, Exp(lambda1)); % Read1 Poisson
    source.setArrival(model.classes{2}, Exp(lambda2)); % Read2 Poisson
end
