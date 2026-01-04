% Layered Queueing Network (LQN) - Production 4-Tier J2EE Architecture Model
%
% This example demonstrates a comprehensive 4-layer J2EE system architecture:
%
%   Layer 1: CLIENT LAYER
%   - Reference task with 1 user (closed class)
%   - Represents browser clients with sequential user operations
%
%   Layer 2: APPLICATION LAYER
%   - Service task handling 4 concurrent web requests
%   - Multiple entries (renderHomePage, processProductPage, handleLogin, checkoutWorkflow)
%   - Complex activity precedence with or-fork/or-join branching
%
%   Layer 3: DATABASE LAYER
%   - Service task managing database operations
%   - Multiple entries for different query types (read, write, query, log)
%   - Sequential activity chains
%
% The model demonstrates:
% - Multiple entries per task
% - Complex activity precedence patterns (serial, or-fork, or-join)
% - Synchronous call routing between layers
% - Load-dependent processing characteristics
% - Production-grade LQN features
%
% This is a reference model for validating SolverLN against complex layered systems.
% For cache support variant, see: example_layered_production_4tier_cache.m

clear all;

%% Create the model
model = LayeredNetwork('testLQN3');

%% Layer 1: CLIENT LAYER
P0 = Processor(model, 'P0', 1, SchedStrategy.PS);
T0 = Task(model, 'T0', 1, SchedStrategy.REF).on(P0);  % Reference task with 1 user
E0 = Entry(model, 'E0').on(T0);  % Single entry point for all client requests

%% Layer 2: APPLICATION SERVER LAYER
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1);  % FCFS for serialization
E10 = Entry(model, 'E10').on(T1);  % renderHomePage
E11 = Entry(model, 'E11').on(T1);  % processProductPage
E12 = Entry(model, 'E12').on(T1);  % handleLogin
E13 = Entry(model, 'E13').on(T1);  % checkoutWorkflow

%% Layer 3: DATABASE LAYER
P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2);  % Database serialization
E20 = Entry(model, 'E20').on(T2);  % readUserProfile
E21 = Entry(model, 'E21').on(T2);  % getProductInfo
E22 = Entry(model, 'E22').on(T2);  % writeOrderData
E23 = Entry(model, 'E23').on(T2);  % logSessionEvent

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
B6 = Activity(model, 'B6', Exp(1.0)).on(T1);                         % computeTotalPrice
B7 = Activity(model, 'B7', Exp(1.0)).on(T1).synchCall(E22, 1.0);    % WriteOrder
B7a = Activity(model, 'B7a', Exp(1.0)).on(T1);                       % skipLog (30% probability)
B7b = Activity(model, 'B7b', Exp(1.0)).on(T1).synchCall(E23, 1.0);  % logOrder (70% probability)
B8 = Activity(model, 'B8', Exp(1.0)).on(T1).repliesTo(E13);          % success

%% DATABASE LAYER ACTIVITIES
% Entry E20: readUserProfile
C0 = Activity(model, 'C0', Exp(1.0)).on(T2).boundTo(E20);      % parseReadQuery
C1 = Activity(model, 'C1', Exp(1.0)).on(T2).repliesTo(E20);    % fetchFromDisk

% Entry E21: getProductInfo
C2 = Activity(model, 'C2', Exp(1.0)).on(T2).boundTo(E21).repliesTo(E21);  % readCache (in-DB cache)

% Entry E22: writeOrderData
C3 = Activity(model, 'C3', Exp(1.0)).on(T2).boundTo(E22);       % parseWriteQuery
C4 = Activity(model, 'C4', Exp(1.0)).on(T2);                    % writeToDisk
C5 = Activity(model, 'C5', Exp(1.0)).on(T2).repliesTo(E22);     % updateIndex

% Entry E23: logSessionEvent
C6 = Activity(model, 'C6', Exp(1.0)).on(T2).boundTo(E23).repliesTo(E23);  % appendToLog

%% ACTIVITY PRECEDENCES - CLIENT LAYER
% Sequential flow through all client operations
T0.addPrecedence(ActivityPrecedence.Serial(A0, A1, A2, A3));

%% ACTIVITY PRECEDENCES - APPLICATION LAYER
% Path 1: renderHomePage flow
T1.addPrecedence(ActivityPrecedence.Serial(B0, B1));

% Path 2: processProductPage flow
T1.addPrecedence(ActivityPrecedence.Serial(B2, B3));

% Path 3: checkoutWorkflow flow with branching
T1.addPrecedence(ActivityPrecedence.Serial(B5, B6, B7));      % Sequential: verify -> compute -> write
T1.addPrecedence(ActivityPrecedence.OrFork(B7, {B7a, B7b}, [0.7, 0.3]));  % 70% skip, 30% log
T1.addPrecedence(ActivityPrecedence.OrJoin({B7a, B7b}, B8));  % Converge branches to success

%% ACTIVITY PRECEDENCES - DATABASE LAYER
% Query path 1: readUserProfile
T2.addPrecedence(ActivityPrecedence.Serial(C0, C1));

% Query path 2: writeOrderData with index update
T2.addPrecedence(ActivityPrecedence.Serial(C3, C4, C5));

%% SOLVE WITH SOLVERMVA (Layer solver) inside SolverLN
fprintf('\n=== Solving 4-Tier LQN with SolverLN + MVA ===\n');
lnoptions = LN.defaultOptions;
lnoptions.verbose = VerboseLevel.STD;
mvaopt = MVA.defaultOptions;
mvaopt.verbose = VerboseLevel.SILENT;

solver = LN(model, @(layer) MVA(layer, mvaopt), lnoptions);
AvgTable = solver.getAvgTable;

fprintf('\n=== Results ===\n');
disp(AvgTable);

%% Alternative: Solve with SolverFluid (approximation)
% Uncomment to compare with fluid approximation
% fprintf('\n=== Solving 4-Tier LQN with SolverLN + SolverFluid (Approximation) ===\n');
% lnoptions = LN.defaultOptions;
% lnoptions.verbose = VerboseLevel.STD;
% fluidopt = SolverFluid.defaultOptions;
% fluidopt.verbose = VerboseLevel.SILENT;
%
% solver_fluid = LN(model, @(layer) SolverFluid(layer, fluidopt), lnoptions);
% AvgTable_fluid = solver_fluid.getAvgTable;
% disp('Fluid Approximation Results:');
% disp(AvgTable_fluid);
