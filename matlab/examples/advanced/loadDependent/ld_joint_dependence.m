clear node jobclass

% Joint-dependent (LJD) model example
% Service rate scaling depends on per-class population vector (n1, n2)
% Scaling table is stored as linearized lookup table

N1 = 3; % number of class-1 jobs
N2 = 2; % number of class-2 jobs

%% Create model
model = Network('JointDependenceModel');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', N1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(2.0));
node{2}.setService(jobclass{1}, Exp.fitMean(1.5));
node{2}.setService(jobclass{2}, Exp.fitMean(2.5));

% Define cutoffs for per-class population
cutoffs = [N1, N2];  % max population per class

% Create scaling table for joint dependence
% Table size = (N1+1) * (N2+1) = 4 * 3 = 12
% Linearized index: idx = 1 + n1 + n2*(N1+1) (MATLAB 1-indexed)

% Example scaling: more jobs of either class -> lower scaling (higher congestion)
% Row-major order would be: (0,0), (1,0), (2,0), (3,0), (0,1), ...
% But we use column-major linearization for the lookup
scalingTable = zeros(1, prod(cutoffs + 1));
for n1 = 0:N1
    for n2 = 0:N2
        idx = ljd_linearize([n1, n2], cutoffs);
        % Example: scaling decreases as total population increases
        totalN = n1 + n2;
        if totalN == 0
            scalingTable(idx) = 1.0;
        elseif totalN <= 2
            scalingTable(idx) = 1.0 / totalN;  % slower service with more jobs
        else
            scalingTable(idx) = 0.5 / totalN;  % even slower at high load
        end
    end
end

% Set joint dependence on Queue1
node{2}.setJointDependence(scalingTable, cutoffs);

% Routing
P = model.initRoutingMatrix();
P{1,1} = model.serialRouting(node);
P{2,2} = model.serialRouting(node);
model.link(P);

%% Solve with MVA (qd method supports LJD)
fprintf('Solving joint-dependent model with MVA (qd):\n');
AvgTableMVA = MVA(model, 'method', 'qd').getAvgTable;
disp(AvgTableMVA);

%% For comparison, solve an equivalent load-dependent model
% This is just for validation - not exact equivalence
fprintf('\nNote: Joint dependence allows class-specific scaling.\n');
fprintf('This example uses total population scaling as a simpler illustration.\n');
