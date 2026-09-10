function [simDoc, section] = saveRetrialDistributions(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVERETRIALDISTRIBUTIONS(SIMDOC, SECTION, IND)
%
% Generates XML for retrial delay distributions for JMT Queue sections.
% The retrial distribution tells JMT how long a blocked customer waits
% in the orbit before retrying.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

retrialNode = simDoc.createElement('parameter');
retrialNode.setAttribute('array', 'true');
retrialNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategy');
retrialNode.setAttribute('name', 'retrialDistributions');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();

% Get the actual Queue node to access retrialDelays
nodes = self.model.getNodes();
currentNode = nodes{ind};

for r=1:numOfClasses
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    retrialNode.appendChild(refClassNode);

    serviceTimeStrategyNode = simDoc.createElement('subParameter');
    serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
    serviceTimeStrategyNode.setAttribute('name', 'ServiceTimeStrategy');

    % Check if retrial delay is defined for this class
    hasRetrial = false;
    retrialDist = [];
    if isa(currentNode, 'Queue') && ~isempty(currentNode.retrialDelays)
        if r <= size(currentNode.retrialDelays, 2) && ~isempty(currentNode.retrialDelays{1, r})
            hasRetrial = true;
            retrialDist = currentNode.retrialDelays{1, r};
        end
    end

    if ~hasRetrial
        % Default: Exp(1) placeholder
        distributionNode = simDoc.createElement('subParameter');
        distributionNode.setAttribute('classPath', 'jmt.engine.random.Exponential');
        distributionNode.setAttribute('name', 'Exponential');
        serviceTimeStrategyNode.appendChild(distributionNode);

        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', 'jmt.engine.random.ExponentialPar');
        distrParNode.setAttribute('name', 'distrPar');
        subParNodeLambda = simDoc.createElement('subParameter');
        subParNodeLambda.setAttribute('classPath', 'java.lang.Double');
        subParNodeLambda.setAttribute('name', 'lambda');
        subParValue = simDoc.createElement('value');
        subParValue.appendChild(simDoc.createTextNode('1.000000000000'));
        subParNodeLambda.appendChild(subParValue);
        distrParNode.appendChild(subParNodeLambda);
        serviceTimeStrategyNode.appendChild(distrParNode);
    else
        % see _kb/06-solver-catalog.md (Wrappers: JMT phase-type distribution export is shared)
        emitPhaseType = isa(retrialDist, 'PH') || isa(retrialDist, 'APH') || ...
            isa(retrialDist, 'Coxian') || ...
            (isa(retrialDist, 'HyperExp') && retrialDist.getNumberOfPhases() > 2);

        if emitPhaseType
            javaClass = 'jmt.engine.random.PhaseTypeDistr';
            javaParClass = 'jmt.engine.random.PhaseTypePar';
            distName = 'Phase-Type';
        elseif isa(retrialDist, 'Exp')
            javaClass = 'jmt.engine.random.Exponential';
            javaParClass = 'jmt.engine.random.ExponentialPar';
            distName = 'Exponential';
        elseif isa(retrialDist, 'Erlang')
            javaClass = 'jmt.engine.random.Erlang';
            javaParClass = 'jmt.engine.random.ErlangPar';
            distName = 'Erlang';
        elseif isa(retrialDist, 'HyperExp')
            javaClass = 'jmt.engine.random.HyperExp';
            javaParClass = 'jmt.engine.random.HyperExpPar';
            distName = 'Hyperexponential';
        elseif isa(retrialDist, 'Det')
            javaClass = 'jmt.engine.random.DeterministicDistr';
            javaParClass = 'jmt.engine.random.DeterministicDistrPar';
            distName = 'Deterministic';
        elseif isa(retrialDist, 'Gamma')
            javaClass = 'jmt.engine.random.GammaDistr';
            javaParClass = 'jmt.engine.random.GammaDistrPar';
            distName = 'Gamma';
        elseif isa(retrialDist, 'Uniform')
            javaClass = 'jmt.engine.random.Uniform';
            javaParClass = 'jmt.engine.random.UniformPar';
            distName = 'Uniform';
        else
            % Fallback: treat as exponential with the distribution's rate
            javaClass = 'jmt.engine.random.Exponential';
            javaParClass = 'jmt.engine.random.ExponentialPar';
            distName = 'Exponential';
        end

        distributionNode = simDoc.createElement('subParameter');
        distributionNode.setAttribute('classPath', javaClass);
        distributionNode.setAttribute('name', distName);
        serviceTimeStrategyNode.appendChild(distributionNode);

        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', javaParClass);
        distrParNode.setAttribute('name', 'distrPar');

        if emitPhaseType
            [alpha, T, phases] = phaseTypeRepres(retrialDist);

            % Alpha vector
            subParNodeAlpha = simDoc.createElement('subParameter');
            subParNodeAlpha.setAttribute('array', 'true');
            subParNodeAlpha.setAttribute('classPath', 'java.lang.Object');
            subParNodeAlpha.setAttribute('name', 'alpha');
            subParNodeAlphaVec = simDoc.createElement('subParameter');
            subParNodeAlphaVec.setAttribute('array', 'true');
            subParNodeAlphaVec.setAttribute('classPath', 'java.lang.Object');
            subParNodeAlphaVec.setAttribute('name', 'vector');
            for k=1:phases
                subParNodeAlphaElem = simDoc.createElement('subParameter');
                subParNodeAlphaElem.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlphaElem.setAttribute('name', 'entry');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', alpha(k))));
                subParNodeAlphaElem.appendChild(subParValue);
                subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
            end

            % T matrix
            subParNodeT = simDoc.createElement('subParameter');
            subParNodeT.setAttribute('array', 'true');
            subParNodeT.setAttribute('classPath', 'java.lang.Object');
            subParNodeT.setAttribute('name', 'T');
            for k=1:phases
                subParNodeTvec = simDoc.createElement('subParameter');
                subParNodeTvec.setAttribute('array', 'true');
                subParNodeTvec.setAttribute('classPath', 'java.lang.Object');
                subParNodeTvec.setAttribute('name', 'vector');
                for j=1:phases
                    subParNodeTElem = simDoc.createElement('subParameter');
                    subParNodeTElem.setAttribute('classPath', 'java.lang.Double');
                    subParNodeTElem.setAttribute('name', 'entry');
                    subParValue = simDoc.createElement('value');
                    if k==j
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -abs(T(k,j)))));
                    else
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', abs(T(k,j)))));
                    end
                    subParNodeTElem.appendChild(subParValue);
                    subParNodeTvec.appendChild(subParNodeTElem);
                end
                subParNodeT.appendChild(subParNodeTvec);
            end

            subParNodeAlpha.appendChild(subParNodeAlphaVec);
            distrParNode.appendChild(subParNodeAlpha);
            distrParNode.appendChild(subParNodeT);

        elseif isa(retrialDist, 'HyperExp')
            % 2-phase only: HyperExpPar carries exactly (p, lambda1, lambda2).
            % Emitting a bare `lambda` here, as the generic fallback below did,
            % silently degraded the retrial delay to an exponential.
            [alpha, T] = phaseTypeRepres(retrialDist);
            subParNodeP = simDoc.createElement('subParameter');
            subParNodeP.setAttribute('classPath', 'java.lang.Double');
            subParNodeP.setAttribute('name', 'p');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', alpha(1))));
            subParNodeP.appendChild(subParValue);
            distrParNode.appendChild(subParNodeP);

            subParNodeLambda1 = simDoc.createElement('subParameter');
            subParNodeLambda1.setAttribute('classPath', 'java.lang.Double');
            subParNodeLambda1.setAttribute('name', 'lambda1');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -T(1,1))));
            subParNodeLambda1.appendChild(subParValue);
            distrParNode.appendChild(subParNodeLambda1);

            subParNodeLambda2 = simDoc.createElement('subParameter');
            subParNodeLambda2.setAttribute('classPath', 'java.lang.Double');
            subParNodeLambda2.setAttribute('name', 'lambda2');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -T(2,2))));
            subParNodeLambda2.appendChild(subParValue);
            distrParNode.appendChild(subParNodeLambda2);

        elseif isa(retrialDist, 'Exp')
            subParNodeLambda = simDoc.createElement('subParameter');
            subParNodeLambda.setAttribute('classPath', 'java.lang.Double');
            subParNodeLambda.setAttribute('name', 'lambda');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', retrialDist.getRate())));
            subParNodeLambda.appendChild(subParValue);
            distrParNode.appendChild(subParNodeLambda);

        elseif isa(retrialDist, 'Det')
            subParNodeT = simDoc.createElement('subParameter');
            subParNodeT.setAttribute('classPath', 'java.lang.Double');
            subParNodeT.setAttribute('name', 't');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', retrialDist.getMean())));
            subParNodeT.appendChild(subParValue);
            distrParNode.appendChild(subParNodeT);

        elseif isa(retrialDist, 'Erlang')
            subParNodeAlpha = simDoc.createElement('subParameter');
            subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
            subParNodeAlpha.setAttribute('name', 'alpha');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', retrialDist.getRate())));
            subParNodeAlpha.appendChild(subParValue);
            distrParNode.appendChild(subParNodeAlpha);

            subParNodeR = simDoc.createElement('subParameter');
            subParNodeR.setAttribute('classPath', 'java.lang.Long');
            subParNodeR.setAttribute('name', 'r');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%d', retrialDist.getNumberOfPhases())));
            subParNodeR.appendChild(subParValue);
            distrParNode.appendChild(subParNodeR);

        else
            % Generic fallback: use rate as exponential
            subParNodeLambda = simDoc.createElement('subParameter');
            subParNodeLambda.setAttribute('classPath', 'java.lang.Double');
            subParNodeLambda.setAttribute('name', 'lambda');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', 1/retrialDist.getMean())));
            subParNodeLambda.appendChild(subParValue);
            distrParNode.appendChild(subParNodeLambda);
        end

        serviceTimeStrategyNode.appendChild(distrParNode);
    end

    retrialNode.appendChild(serviceTimeStrategyNode);
end

section.appendChild(retrialNode);
end

function [alpha, T, phases] = phaseTypeRepres(dist)
% [ALPHA, T, PHASES] = PHASETYPEREPRES(DIST)
%
% Extracts the (alpha, T) acyclic phase-type representation of a Markovian
% distribution DIST from its own (D0,D1) process, so that the exported JMT
% distribution is the SAME distribution rather than a moment-matched refit.
% T = D0 and, since D1 = (-D0*e)*alpha for a PH renewal process, alpha is
% recovered as the row of D1 rescaled by that phase's exit rate. This is the
% same derivation used by MNetwork/refreshStruct.m for sn.impatiencePie.
proc = dist.getProcess();
if ~iscell(proc) || numel(proc) < 2
    line_error(mfilename, sprintf(['Retrial distribution %s has no (D0,D1) ', ...
        'representation, so it cannot be exported as a phase-type.'], class(dist)));
end
D0 = proc{1};
D1 = proc{2};
T = D0;
phases = size(D0, 1);
exitRates = -D0 * ones(phases, 1);
idx = find(exitRates > 1e-10, 1);
if ~isempty(idx)
    alpha = abs(D1(idx, :) / exitRates(idx));
else
    alpha = ones(1, phases) / phases;
end
end
