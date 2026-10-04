% Product-form state-dependent routing, Krzesinski (1987), Perform. Eval. 7:125-143.
%
% Section 5.1.1, adaptive routing with balanced loading: a central server and
% four peripheral centers forming four single-center branches under one level
% of subnetwork nesting. With C_1 = -1 and d_{1,i} = d the entry center routes
%
%   P_{1,i}(N) = (d - n_i) / (4d - V),      V = n_2 + n_3 + n_4 + n_5,
%
% which prefers the least congested peripheral center and caps each of them at
% d customers. The state-independent counterpart of eq. (15) is P_{1,i} = 1/4.
% The throughput gain of SDR over SIR grows with the CPU/IO speed ratio, which
% is the paper's Fig. 4.

clear node jobclass;

d = 1;              % routing coefficient d_{1,i}
N = 4;              % network population
mu1 = 100;          % central server rate
mup = 20;           % peripheral rate

model = Network('sdroute_krzesinski');
node{1} = Queue(model, 'CPU', SchedStrategy.FCFS);
for i = 2:5
    node{i} = Queue(model, sprintf('IO%d',i-1), SchedStrategy.FCFS);
end
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp(mu1));
for i = 2:5
    node{i}.setService(jobclass{1}, Exp(mup));
end

model.addLink(node{1}, node{1});    % denied entry: the busy form of waiting
for i = 2:5
    model.addLink(node{1}, node{i});
    model.addLink(node{i}, node{1});
    node{i}.setProbRouting(jobclass{1}, node{1}, 1.0);
end

% One level of nesting, V = V_1 = {2,3,4,5}, each peripheral center a branch
dcoeff = zeros(1,5);
dcoeff(2:5) = d;
node{1}.setStateDepRouting(jobclass{1}, node{1}, ...
    {[], {node{2}}, {node{3}}, {node{4}}, {node{5}}}, [0 1 1 1 1], -1, dcoeff);

AvgTable = SolverNC(model).getAvgTable

% The same network under the state-independent routing of eq. (15)
sirmodel = Network('sdroute_krzesinski_sir');
snode{1} = Queue(sirmodel, 'CPU', SchedStrategy.FCFS);
for i = 2:5
    snode{i} = Queue(sirmodel, sprintf('IO%d',i-1), SchedStrategy.FCFS);
end
sclass{1} = ClosedClass(sirmodel, 'Class1', N, snode{1}, 0);
snode{1}.setService(sclass{1}, Exp(mu1));
for i = 2:5
    snode{i}.setService(sclass{1}, Exp(mup));
    sirmodel.addLink(snode{1}, snode{i});
    sirmodel.addLink(snode{i}, snode{1});
    snode{1}.setProbRouting(sclass{1}, snode{i}, 0.25);
    snode{i}.setProbRouting(sclass{1}, snode{1}, 1.0);
end
SIRTable = SolverNC(sirmodel).getAvgTable

Tsdr = sum(AvgTable.Tput(2:5));
Tsir = sum(SIRTable.Tput(2:5));
line_printf('\nPeripheral throughput: SDR %.4f, SIR %.4f, gain %.1f%%\n', ...
    Tsdr, Tsir, 100*(Tsdr-Tsir)/Tsir);
