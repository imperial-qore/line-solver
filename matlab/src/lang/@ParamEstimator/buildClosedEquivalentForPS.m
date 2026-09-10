function [eqModel, eqNode] = buildClosedEquivalentForPS(self, node)
% BUILDCLOSEDEQUIVALENTFORPS Build a closed equivalent model for open/mixed networks.
% For open classes, uses the equivalence Z_r = N_r / lambda_r where
% N_r is the effective population (from options.openPopulation) and
% lambda_r is the arrival rate from the Source node.
% Closed classes retain their original population and Delay rates.
%
% Returns a 2-node closed queueing network (Delay + PS Queue) suitable
% for use with MLPS and FMLPS estimators.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;
R = sn.nclasses;

% Get effective population for open classes
if isfield(self.options, 'openPopulation') && ~isempty(self.options.openPopulation)
    Nopen = self.options.openPopulation;
else
    Nopen = 100;
end

% Determine populations and think times per class
N = zeros(1, R);
Z = zeros(1, R);

allNodes = self.model.getNodes;
for n = 1:length(allNodes)
    if isa(allNodes{n}, 'Delay')
        svcProc = allNodes{n}.getService;
        for r = 1:R
            if sn.njobs(r) < Inf
                Z(r) = Z(r) + svcProc{r}.getMean();
            end
        end
    elseif isa(allNodes{n}, 'Source')
        svcProc = allNodes{n}.getService;
        for r = 1:R
            if sn.njobs(r) == Inf
                lambda_r = 1 / svcProc{r}.getMean();
                Z(r) = Nopen / lambda_r;
            end
        end
    end
end

for r = 1:R
    if sn.njobs(r) < Inf
        N(r) = sn.njobs(r);
    else
        N(r) = Nopen;
    end
end

% Convert think times to rates
delayRate = 1 ./ Z;

% Get queue service processes from original node
queueSvcProc = node.getService;

% Build closed equivalent model
eqModel = Network('closed_equiv');
eqDelay = Delay(eqModel, 'Think');
eqQueue = Queue(eqModel, 'Queue1', SchedStrategy.PS);
eqQueue.setNumberOfServers(node.getNumberOfServers());

eqClass = cell(1, R);
for r = 1:R
    eqClass{r} = ClosedClass(eqModel, sprintf('Class%d', r), N(r), eqDelay, 0);
    eqDelay.setService(eqClass{r}, Exp(delayRate(r)));
    eqQueue.setService(eqClass{r}, queueSvcProc{r});
end

P = eqModel.initRoutingMatrix;
for r = 1:R
    P{r} = Network.serialRouting({eqDelay, eqQueue});
end
eqModel.link(P);

eqNode = eqQueue;

end
