% Layered Queueing Network (LQN) - Production 4-Tier J2EE Architecture with Cache
%
% This example extends the 4-tier J2EE model (example_layered_production_4tier.m)
% with a dedicated CACHE LAYER (Layer 4), demonstrating:
%
%   Layer 1: CLIENT LAYER
%   - Reference task with 1 user (closed class)
%   - Represents browser clients with sequential user operations
%
%   Layer 2: APPLICATION LAYER
%   - Service task handling web requests
%   - Multiple entries for different operations
%   - Complex activity precedence with branching
%
%   Layer 3: DATABASE LAYER
%   - Primary data storage
%   - Cache lookup operations (via calls to cache layer)
%   - Write and log operations
%
%   Layer 4: CACHE LAYER (NEW)
%   - Dedicated cache with LRU replacement strategy
%   - Cache capacity: 2 items out of 10 total items
%   - Uniform access distribution
%   - Cache hit/miss modeling with different service times
%
% Model demonstrates:
% - Cache hit/miss precedence patterns (CacheAccess)
% - ItemEntry for cache-aware workloads
% - CacheTask with replacement strategy
% - Integration of caching into multi-layer system
%
% KNOWN LIMITATION:
% SolverLN has issues with routing synchronous cache calls in complex multi-layer
% networks with precedence constraints. This model will fail during layer decomposition.
% For cache support in production systems, use:
%   - Simpler 2-3 layer models (example_layeredcache_async.m)
%   - Or implement caching at application layer (no cross-layer cache calls)
%
% Compare with: example_layered_production_4tier.m (non-cached variant - WORKS)

clear all;

%% Create the model
model = LayeredNetwork('testLQN3_Cache');

%% Layer 1: CLIENT LAYER
P0 = Processor(model, 'P0', 1, SchedStrategy.PS);
T0 = Task(model, 'T0', 1, SchedStrategy.REF).on(P0);  % Reference task with 1 user
E0 = Entry(model, 'E0').on(T0);  % Single entry point

%% Layer 2: APPLICATION SERVER LAYER
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1);
E10 = Entry(model, 'E10').on(T1);  % renderHomePage
E11 = Entry(model, 'E11').on(T1);  % processProductPage
E12 = Entry(model, 'E12').on(T1);  % handleLogin
E13 = Entry(model, 'E13').on(T1);  % checkoutWorkflow

%% Layer 3: DATABASE LAYER
P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2);
E20 = Entry(model, 'E20').on(T2);  % readUserProfile
E21 = Entry(model, 'E21').on(T2);  % getProductInfo
E22 = Entry(model, 'E22').on(T2);  % writeOrderData
E23 = Entry(model, 'E23').on(T2);  % logSessionEvent

%% Layer 4: CACHE LAYER (NEW)
% Cache configuration: 10 items total, capacity of 2, LRU replacement
totalitems = 10;
cachecapacity = 2;
pAccess = DiscreteSampler((1/totalitems)*ones(1, totalitems));  % Uniform access
P3 = Processor(model, 'P3', 1, SchedStrategy.PS);
T3 = CacheTask(model, 'T3', totalitems, cachecapacity, ReplacementStrategy.LRU, 1).on(P3);
E3 = ItemEntry(model, 'E3', totalitems, pAccess).on(T3);  % accessCache

%% CLIENT LAYER ACTIVITIES
% Sequential processing flow: login → browse → shop → checkout
A0 = Activity(model, 'A0', Exp(1.0)).on(T0).boundTo(E0).synchCall(E12, 1.0);  % login -> handleLogin
A1 = Activity(model, 'A1', Exp(1.0)).on(T0).synchCall(E10, 1.0);              % browseSession -> renderHomePage
A2 = Activity(model, 'A2', Exp(1.0)).on(T0).synchCall(E11, 1.0);              % viewProductPage -> processProductPage
A3 = Activity(model, 'A3', Exp(1.0)).on(T0).synchCall(E13, 1.0);              % submitCheckout -> checkoutWorkflow

%% APPLICATION LAYER ACTIVITIES
% Entry E10: renderHomePage
B0 = Activity(model, 'B0', Exp(1.0)).on(T1).boundTo(E10);      % renderStaticContent
B1 = Activity(model, 'B1', Exp(1.0)).on(T1).repliesTo(E10);    % fetchRecommendations

% Entry E11: processProductPage
B2 = Activity(model, 'B2', Exp(1.0)).on(T1).boundTo(E11);                      % parseRequest
B3 = Activity(model, 'B3', Exp(1.0)).on(T1).synchCall(E21, 1.0).repliesTo(E11);  % getProductDetails -> getProductInfo

% Entry E12: handleLogin
B4 = Activity(model, 'B4', Exp(1.0)).on(T1).boundTo(E12).synchCall(E20, 1.0).repliesTo(E12);  % validateCredentials -> readUserProfile

% Entry E13: checkoutWorkflow
B5 = Activity(model, 'B5', Exp(1.0)).on(T1).boundTo(E13);           % verifyCart
B6 = Activity(model, 'B6', Exp(1.0)).on(T1);                        % computeTotalPrice
B7 = Activity(model, 'B7', Exp(1.0)).on(T1).synchCall(E22, 1.0);   % WriteOrder (no reply - or-fork will reply)
B7a = Activity(model, 'B7a', Exp(1.0)).on(T1).repliesTo(E13);       % skipLog (70% probability)
B7b = Activity(model, 'B7b', Exp(1.0)).on(T1).synchCall(E23, 1.0).repliesTo(E13);  % logOrder (30% probability)

%% DATABASE LAYER ACTIVITIES
% Entry E20: readUserProfile
C0 = Activity(model, 'C0', Exp(1.0)).on(T2).boundTo(E20);                  % parseReadQuery
C1 = Activity(model, 'C1', Exp(1.0)).on(T2).synchCall(E3, 1.0).repliesTo(E20);  % fetchFromDisk -> accessCache

% Entry E21: getProductInfo (cache-aware)
C2 = Activity(model, 'C2', Exp(1.0)).on(T2).boundTo(E21).synchCall(E3, 1.0).repliesTo(E21);  % readCache -> accessCache

% Entry E22: writeOrderData
C3 = Activity(model, 'C3', Exp(1.0)).on(T2).boundTo(E22);       % parseWriteQuery
C4 = Activity(model, 'C4', Exp(1.0)).on(T2);                    % writeToDisk
C5 = Activity(model, 'C5', Exp(1.0)).on(T2).repliesTo(E22);     % updateIndex

% Entry E23: logSessionEvent
C6 = Activity(model, 'C6', Exp(1.0)).on(T2).boundTo(E23).repliesTo(E23);  % appendToLog

%% CACHE LAYER ACTIVITIES
% Immediate arrival processing
D0 = Activity(model, 'D0', Immediate()).on(T3).boundTo(E3);      % arrival
D1a = Activity(model, 'D1a', Exp(1.0)).on(T3).repliesTo(E3);     % cache hit (fast)
D1b = Activity(model, 'D1b', Exp(0.5)).on(T3).repliesTo(E3);     % cache miss (slow)

%% ACTIVITY PRECEDENCES - CLIENT LAYER
% Sequential flow through all client operations
T0.addPrecedence(ActivityPrecedence.Serial(A0, A1, A2, A3));

%% ACTIVITY PRECEDENCES - APPLICATION LAYER
% Path 1: renderHomePage flow
T1.addPrecedence(ActivityPrecedence.Serial(B0, B1));

% Path 2: processProductPage flow
T1.addPrecedence(ActivityPrecedence.Serial(B2, B3));

% Path 3: checkoutWorkflow flow with branching
T1.addPrecedence(ActivityPrecedence.Serial(B5, B6, B7));              % Sequential: verify -> compute -> write
T1.addPrecedence(ActivityPrecedence.OrFork(B7, {B7a, B7b}, [0.7, 0.3]));  % 70% skip logging, 30% log
% Note: Both B7a and B7b reply to E13, providing implicit join

%% ACTIVITY PRECEDENCES - DATABASE LAYER
% Query path 1: readUserProfile
T2.addPrecedence(ActivityPrecedence.Serial(C0, C1));

% Query path 2: writeOrderData with index update
T2.addPrecedence(ActivityPrecedence.Serial(C3, C4, C5));

%% ACTIVITY PRECEDENCES - CACHE LAYER
% Cache access precedence: arrival determines hit/miss
T3.addPrecedence(ActivityPrecedence.CacheAccess(D0, {D1a, D1b}));

%% SOLVE WITH SOLVERMVA (Layer solver) inside SolverLN
fprintf('\n=== Solving 4-Tier LQN with Cache using SolverLN + MVA ===\n');
lnoptions = LN.defaultOptions;
lnoptions.verbose = VerboseLevel.STD;
mvaopt = MVA.defaultOptions;
mvaopt.verbose = VerboseLevel.SILENT;

solver = LN(model, @(layer) MVA(layer, mvaopt), lnoptions);
AvgTable = solver.getAvgTable;

fprintf('\n=== Results ===\n');
disp(AvgTable);

%% Display cache statistics if available
fprintf('\n=== Cache Performance (if available) ===\n');
if isprop(solver, 'cache_stats')
    disp(solver.cache_stats);
else
    fprintf('Cache statistics not available from this solver.\n');
end

%% Comparison notes
fprintf('\n=== Comparison with Non-Cached Variant ===\n');
fprintf('Expected impacts of cache layer:\n');
fprintf('- Lower product info query response times (cache hits)\n');
fprintf('- Higher overall throughput (cached lookups are fast)\n');
fprintf('- Cache replacement activity (LRU eviction)\n');
fprintf('See: example_layered_production_4tier.m for non-cached baseline\n');
