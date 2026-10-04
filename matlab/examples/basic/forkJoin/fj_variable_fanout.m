function model = fj_variable_fanout(mode)
% FJ_VARIABLE_FANOUT Fork-join with a variable forking level
%
% MODEL = FJ_VARIABLE_FANOUT(MODE) builds a closed fork-join network whose Fork
% emits a number of tasks per link that is NOT one fixed number. MODE selects
% which way the degree varies:
%
%   'fixed'   the classic fork, one task on each of two links (the baseline)
%   'vector'  three tasks towards Queue2 and one towards Queue1
%   'random'  one or three tasks per link, each with probability 1/2
%   'prob'    the branch towards Queue2 fires only half the time
%
% Exact under SolverJMT and SolverLDES, which draw the degree at the fork epoch.
% SolverMVA's MMT method sees the EXPECTED degree. The exact CTMC/SSA path
% accepts 'fixed' and 'vector' -- on 'vector' its answer matches JMT to
% simulation noise -- and refuses 'random' and 'prob' by name, because its tag
% construction fixes the sibling count, and the sibling SET, when the state
% space is built.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1
    mode = 'fixed';
end

model = Network('model');

delay = Delay(model,'Delay');
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);
fork = Fork(model,'Fork');
join = Join(model,'Join',fork);

jobclass1 = ClosedClass(model,'class1',4,delay);

delay.setService(jobclass1,Exp(1.0));
queue1.setService(jobclass1,Exp(2.0));
queue2.setService(jobclass1,Exp(2.0));

switch mode
    case 'fixed'
        % no override: the baseline every other mode is compared against
    case 'vector'
        fork.setTasksPerLink(jobclass1, 3, queue2);
    case 'random'
        fork.setTasksPerLinkDistribution(jobclass1, DiscreteSampler([0.5 0.5],[1 3]));
    case 'prob'
        % an uncertain branch needs a Join that does not wait for it
        join.setStrategy(jobclass1, JoinStrategy.PARTIAL);
        join.setRequired(jobclass1, 1);
        fork.setBranchProbability(jobclass1, queue2, 0.5);
    otherwise
        line_error(mfilename,'Unknown mode; use fixed, vector, random or prob.');
end

P = model.initRoutingMatrix();
P{jobclass1,jobclass1}(delay,fork) = 1.0;
P{jobclass1,jobclass1}(fork,queue1) = 1.0;
P{jobclass1,jobclass1}(fork,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join) = 1.0;
P{jobclass1,jobclass1}(queue2,join) = 1.0;
P{jobclass1,jobclass1}(join,delay) = 1.0;
model.link(P);

if nargout == 0
    for m = {'fixed','vector','random','prob'}
        mm = fj_variable_fanout(m{1});
        line_printf('\n<strong>mode = %s</strong>', m{1});
        AvgTable = SolverJMT(mm,'seed',23000).getAvgTable();
        disp(AvgTable);
    end
end
end
