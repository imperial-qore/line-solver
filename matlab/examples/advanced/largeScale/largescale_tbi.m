%% Large-scale transient fluid analysis with trajectory-based iteration (TBI)
% This example builds a vehicle-sharing-style closed network with M
% queueing stations and one Erlang transit delay per ordered station pair,
% giving O(M^2) stations and, with 16 Erlang phases per delay, a fluid ODE
% with about 1500 state variables. At this size the monolithic stiff fluid
% solution (method 'closing') takes several minutes, dominated by the
% Jacobian factorizations, while trajectory-based iteration (TBI) solves
% each station cell separately against frozen inbound trajectories and
% completes in seconds. See:
%
%   M. Sheldon, D. Tuncer, G. Casale, "TBI: Transient Hierarchical
%   Modeling of Large-Scale Vehicle Sharing Systems", IEEE Transactions
%   on Intelligent Transportation Systems.
%
% Each cell holds one queueing station together with its outbound transit
% delays, mirroring the spatial submodels of the paper.

clear model;
M = 10;    % queueing stations; the model has M + M*(M-1) stations in total
N = 250;   % closed population (vehicles)
kph = 16;  % Erlang phases per transit delay
Tend = 20; % transient horizon

model = Network('tbi_largescale');
rng(1,'twister');
Q = cell(1,M);
for i = 1:M
    Q{i} = Queue(model, sprintf('Q%d',i), SchedStrategy.PS);
end
D = cell(M,M);
for i = 1:M
    for j = 1:M
        if i ~= j
            D{i,j} = Delay(model, sprintf('D%d_%d',i,j));
        end
    end
end
job = ClosedClass(model, 'C1', N, Q{1});
lambda = zeros(M,M);
for i = 1:M
    Q{i}.setService(job, Erlang.fitMeanAndOrder(1/(1 + 3*rand()), 4));
    for j = 1:M
        if i ~= j
            D{i,j}.setService(job, Erlang.fitMeanAndOrder(0.2 + 2*rand(), kph));
            lambda(i,j) = rand();
        end
    end
end
P = model.initRoutingMatrix;
for i = 1:M
    pr = lambda(i,:) / sum(lambda(i,:));
    for j = 1:M
        if i ~= j
            P{1}(Q{i}, D{i,j}) = pr(j);
            P{1}(D{i,j}, Q{j}) = 1;
        end
    end
end
model.link(P);
model.initDefault;

% one cell per queueing station plus its outbound transit delays
names = cellfun(@(s) s.name, model.stations, 'UniformOutput', false);
cells = cell(1, M);
for i = 1:M
    idx = find(strcmp(names, sprintf('Q%d',i)));
    for j = 1:M
        if i ~= j
            idx(end+1) = find(strcmp(names, sprintf('D%d_%d',i,j))); %#ok<AGROW>
        end
    end
    cells{i} = idx;
end

solver = SolverFLD(model, 'method', 'tbi', 'timespan', [0 Tend], 'stiff', true);
solver.options.config.tbi_cells = cells;
tic;
AvgTable = solver.getAvgTable();
fprintf('TBI solved %d stations (%d ODE variables) in %.1f seconds.\n', ...
    model.getNumberOfStations(), M*(M-1)*kph + M*4, toc);
AvgTable(1:M,:) % queue-length summary at the queueing stations

% For comparison, the undecomposed solution of the same model:
%   solver = SolverFLD(model, 'method', 'closing', 'timespan', [0 Tend], 'stiff', true);
% takes several minutes on the same machine (about 80 seconds already at
% M=8, and beyond 10 minutes at M=12).
