function [simDoc, section] = saveImpatience(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEIMPATIENCE(SIMDOC, SECTION, IND)
%
% Generates XML for impatience (reneging) distributions for JMT Queue sections
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

impatienceNode = simDoc.createElement('parameter');
impatienceNode.setAttribute('array', 'true');
impatienceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Impatience');
impatienceNode.setAttribute('name', 'Impatience');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
i = sn.nodeToStation(ind);

for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    impatienceNode.appendChild(refClassNode);

    % see _kb/06-solver-catalog.md (Wrappers: JMT impatience export)
    hasBalking = false;
    if i > 0 && isfield(sn, 'balkingStrategy') && ~isempty(sn.balkingStrategy)
        if sn.balkingStrategy(i, r) == BalkingStrategy.QUEUE_LENGTH
            hasBalking = true;
        end
    end
    % see _kb/06-solver-catalog.md (Wrappers: JMT impatience export, RENEGING detection)
    hasImpatience = false;
    if i > 0 && isfield(sn, 'impatienceClass') && ~isempty(sn.impatienceClass)
        if sn.impatienceClass(i, r) == ImpatienceType.RENEGING
            hasImpatience = true;
        end
    end

    impatienceStrategyNode = simDoc.createElement('subParameter');

    if hasBalking
        % Balking: emit a Balking strategy wrapping a load-dependent range
        % table {from -> probability} plus the priorityActivated flag.
        impatienceStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Balking');
        impatienceStrategyNode.setAttribute('name', 'Balking');
        [simDoc, impatienceStrategyNode] = saveBalkingStrategy(self, simDoc, impatienceStrategyNode, sn.balkingThresholds{i, r}, sn.nservers(i));
        impatienceNode.appendChild(impatienceStrategyNode);
        continue;
    end

    impatienceStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Reneging');
    impatienceStrategyNode.setAttribute('name', 'Reneging');

    if ~hasImpatience
        % No impatience defined - use null
        subParValue = simDoc.createElement('value');
        subParValue.appendChild(simDoc.createTextNode('null'));
        impatienceStrategyNode.appendChild(subParValue);
    else
        % Impatience is defined - generate distribution XML
        procType = sn.impatienceType(i, r);
        impatiencePhases = sn.impatiencePhases(i, r);

        % see _kb/06-solver-catalog.md (Wrappers: JMT impatience export, HyperExp arity)
        emitPhaseType = (procType == ProcessType.PH) || ...
            (procType == ProcessType.APH) || ...
            (procType == ProcessType.COXIAN) || ...
            (impatiencePhases > 2 && procType == ProcessType.HYPEREXP);

        distributionNode = simDoc.createElement('subParameter');

        if emitPhaseType
            javaClass = 'jmt.engine.random.PhaseTypeDistr';
            javaParClass = 'jmt.engine.random.PhaseTypePar';
        else
        switch procType
            case ProcessType.DET
                javaClass = 'jmt.engine.random.DeterministicDistr';
                javaParClass = 'jmt.engine.random.DeterministicDistrPar';
            case ProcessType.ERLANG
                javaClass = 'jmt.engine.random.Erlang';
                javaParClass = 'jmt.engine.random.ErlangPar';
            case ProcessType.EXP
                javaClass = 'jmt.engine.random.Exponential';
                javaParClass = 'jmt.engine.random.ExponentialPar';
            case ProcessType.GAMMA
                javaClass = 'jmt.engine.random.GammaDistr';
                javaParClass = 'jmt.engine.random.GammaDistrPar';
            case ProcessType.HYPEREXP
                javaClass = 'jmt.engine.random.HyperExp';
                javaParClass = 'jmt.engine.random.HyperExpPar';
            case ProcessType.PARETO
                javaClass = 'jmt.engine.random.Pareto';
                javaParClass = 'jmt.engine.random.ParetoPar';
            case ProcessType.WEIBULL
                javaClass = 'jmt.engine.random.Weibull';
                javaParClass = 'jmt.engine.random.WeibullPar';
            case ProcessType.LOGNORMAL
                javaClass = 'jmt.engine.random.Lognormal';
                javaParClass = 'jmt.engine.random.LognormalPar';
            case ProcessType.UNIFORM
                javaClass = 'jmt.engine.random.Uniform';
                javaParClass = 'jmt.engine.random.UniformPar';
            otherwise
                line_error(mfilename, sprintf('Unsupported impatience distribution type: %s', ProcessType.toText(procType)));
        end
        end

        distributionNode.setAttribute('classPath', javaClass);
        if emitPhaseType
            distributionNode.setAttribute('name', 'Phase-Type');
        else
        switch procType
            case ProcessType.EXP
                distributionNode.setAttribute('name', 'Exponential');
            case ProcessType.HYPEREXP
                distributionNode.setAttribute('name', 'Hyperexponential');
            otherwise
                distributionNode.setAttribute('name', ProcessType.toText(procType));
        end
        end
        impatienceStrategyNode.appendChild(distributionNode);

        % Create distribution parameters
        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', javaParClass);
        distrParNode.setAttribute('name', 'distrPar');

        % Get impatience parameters
        impatienceMu = sn.impatienceMu(i, r);
        impatiencePhi = sn.impatiencePhi(i, r);

        if emitPhaseType
            % see _kb/06-solver-catalog.md (Wrappers: JMT impatience export, cell-array deref)
            impatienceProc = sn.impatienceProc{i, r};
            impatiencePie = sn.impatiencePie{i, r};
            phases = impatiencePhases;
            PH = impatienceProc{1};
            alpha = abs(impatiencePie);

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
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -abs(PH(k,j)))));
                    else
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', abs(PH(k,j)))));
                    end
                    subParNodeTElem.appendChild(subParValue);
                    subParNodeTvec.appendChild(subParNodeTElem);
                end
                subParNodeT.appendChild(subParNodeTvec);
            end

            subParNodeAlpha.appendChild(subParNodeAlphaVec);
            distrParNode.appendChild(subParNodeAlpha);
            distrParNode.appendChild(subParNodeT);
        else
        switch procType
            case ProcessType.DET
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 't');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', 1/impatienceMu(1))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);

            case ProcessType.EXP
                subParNodeLambda = simDoc.createElement('subParameter');
                subParNodeLambda.setAttribute('classPath', 'java.lang.Double');
                subParNodeLambda.setAttribute('name', 'lambda');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', impatienceMu(1))));
                subParNodeLambda.appendChild(subParValue);
                distrParNode.appendChild(subParNodeLambda);

            case ProcessType.ERLANG
                phases = sn.impatiencePhases(i, r);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', impatienceMu(1) * phases)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeR = simDoc.createElement('subParameter');
                subParNodeR.setAttribute('classPath', 'java.lang.Long');
                subParNodeR.setAttribute('name', 'r');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%d', phases)));
                subParNodeR.appendChild(subParValue);
                distrParNode.appendChild(subParNodeR);

            case ProcessType.HYPEREXP
                % 2-phase only: HyperExpPar carries exactly (p,lambda1,lambda2).
                % Braces, not (): these sn fields are cell arrays.
                impatienceProc = sn.impatienceProc{i, r};
                impatiencePie = sn.impatiencePie{i, r};
                subParNodeP = simDoc.createElement('subParameter');
                subParNodeP.setAttribute('classPath', 'java.lang.Double');
                subParNodeP.setAttribute('name', 'p');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', impatiencePie(1))));
                subParNodeP.appendChild(subParValue);
                distrParNode.appendChild(subParNodeP);
                subParNodeLambda1 = simDoc.createElement('subParameter');
                subParNodeLambda1.setAttribute('classPath', 'java.lang.Double');
                subParNodeLambda1.setAttribute('name', 'lambda1');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -impatienceProc{1}(1,1))));
                subParNodeLambda1.appendChild(subParValue);
                distrParNode.appendChild(subParNodeLambda1);
                subParNodeLambda2 = simDoc.createElement('subParameter');
                subParNodeLambda2.setAttribute('classPath', 'java.lang.Double');
                subParNodeLambda2.setAttribute('name', 'lambda2');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', -impatienceProc{1}(2,2))));
                subParNodeLambda2.appendChild(subParValue);
                distrParNode.appendChild(subParNodeLambda2);

            case ProcessType.GAMMA
                scv = impatiencePhi(1);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', 1/scv)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeBeta = simDoc.createElement('subParameter');
                subParNodeBeta.setAttribute('classPath', 'java.lang.Double');
                subParNodeBeta.setAttribute('name', 'beta');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', scv/impatienceMu(1))));
                subParNodeBeta.appendChild(subParValue);
                distrParNode.appendChild(subParNodeBeta);

            case ProcessType.PARETO
                scv = impatiencePhi(1);
                shape = sqrt(1 + 1/scv) + 1;
                scale = 1/impatienceMu(1) * (shape - 1) / shape;
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', shape)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeK = simDoc.createElement('subParameter');
                subParNodeK.setAttribute('classPath', 'java.lang.Double');
                subParNodeK.setAttribute('name', 'k');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', scale)));
                subParNodeK.appendChild(subParValue);
                distrParNode.appendChild(subParNodeK);

            case ProcessType.WEIBULL
                scv = impatiencePhi(1);
                c = sqrt(scv);
                rval = c^(-1.086); % Justus approximation (1976)
                alpha = 1/impatienceMu(1) / gamma(1+1/rval);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', alpha)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeR = simDoc.createElement('subParameter');
                subParNodeR.setAttribute('classPath', 'java.lang.Double');
                subParNodeR.setAttribute('name', 'r');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', rval)));
                subParNodeR.appendChild(subParValue);
                distrParNode.appendChild(subParNodeR);

            case ProcessType.LOGNORMAL
                scv = impatiencePhi(1);
                c = sqrt(scv);
                mu = log(1/impatienceMu(1) / sqrt(c*c + 1));
                sigma = sqrt(log(c*c + 1));
                subParNodeMu = simDoc.createElement('subParameter');
                subParNodeMu.setAttribute('classPath', 'java.lang.Double');
                subParNodeMu.setAttribute('name', 'mu');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', mu)));
                subParNodeMu.appendChild(subParValue);
                distrParNode.appendChild(subParNodeMu);
                subParNodeSigma = simDoc.createElement('subParameter');
                subParNodeSigma.setAttribute('classPath', 'java.lang.Double');
                subParNodeSigma.setAttribute('name', 'sigma');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', sigma)));
                subParNodeSigma.appendChild(subParValue);
                distrParNode.appendChild(subParNodeSigma);

            case ProcessType.UNIFORM
                mean = 1/impatienceMu(1);
                % For uniform [a,b]: mean = (a+b)/2, use a=0 for simplicity
                b = 2 * mean;
                subParNodeMin = simDoc.createElement('subParameter');
                subParNodeMin.setAttribute('classPath', 'java.lang.Double');
                subParNodeMin.setAttribute('name', 'min');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode('0.0'));
                subParNodeMin.appendChild(subParValue);
                distrParNode.appendChild(subParNodeMin);
                subParNodeMax = simDoc.createElement('subParameter');
                subParNodeMax.setAttribute('classPath', 'java.lang.Double');
                subParNodeMax.setAttribute('name', 'max');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f', b)));
                subParNodeMax.appendChild(subParValue);
                distrParNode.appendChild(subParNodeMax);

        end
        end

        impatienceStrategyNode.appendChild(distrParNode);
    end

    impatienceNode.appendChild(impatienceStrategyNode);
end

section.appendChild(impatienceNode);
end

function [simDoc, balkingNode] = saveBalkingStrategy(self, simDoc, balkingNode, thresholds, nservers) %#ok<INUSL>
% Populate a JMT <Balking> impatience element from LINE queue-length balking
% thresholds. LINE stores a list of {minJobs, maxJobs, probability} closed
% intervals; JMT's Balking reads a LoadDependentStrategy whose LDParameter
% ranges are {from -> probability}, selecting the last range with from <=
% queueLength (default 0). We translate the closed intervals into from-based
% breakpoints, inserting explicit 0-probability breakpoints at gaps so that
% queue lengths outside any interval do not inherit a neighbour's probability.
%
% Queue-length convention: LINE evaluates balking against the TOTAL station
% population (in-service + waiting), as do the CTMC/SSA/LDES solvers. JMT's
% Balking evaluates against the number WAITING only (the buffer occupancy).
% When the servers are busy (the regime where balking is meaningful) the two
% differ by exactly the server count S, so we shift every `from` down by S:
% JMT waiting w maps to LINE total n = w + S. This makes the exported JMT
% model reproduce LINE's balking numerically.

% Build sorted (from, prob) breakpoints
lo = zeros(numel(thresholds),1);
hi = zeros(numel(thresholds),1);
pr = zeros(numel(thresholds),1);
for ti = 1:numel(thresholds)
    th = thresholds{ti};
    lo(ti) = th{1};
    hi(ti) = th{2};
    pr(ti) = th{3};
end
[lo, order] = sort(lo);
hi = hi(order);
pr = pr(order);

if isempty(nservers) || ~isfinite(nservers) || nservers < 1
    S = 1;
else
    S = nservers;
end
froms = [];
probs = [];
for ti = 1:numel(lo)
    froms(end+1) = max(0, lo(ti) - S); %#ok<AGROW>
    probs(end+1) = pr(ti);             %#ok<AGROW>
    if ~isinf(hi(ti))
        if ti < numel(lo)
            nextLo = lo(ti+1);
        else
            nextLo = Inf;
        end
        if hi(ti)+1 < nextLo
            froms(end+1) = max(0, hi(ti)+1 - S); %#ok<AGROW>
            probs(end+1) = 0.0;                  %#ok<AGROW>
        end
    end
end

% LoadDependentStrategy wrapper (Parameter #1 of Balking)
ldStrategyNode = simDoc.createElement('subParameter');
ldStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LoadDependentStrategy');
ldStrategyNode.setAttribute('name', 'LoadDependentStrategy');

ldArrayNode = simDoc.createElement('subParameter');
ldArrayNode.setAttribute('array', 'true');
ldArrayNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LDParameter');
ldArrayNode.setAttribute('name', 'LDParameter');

for ti = 1:numel(froms)
    rangeNode = simDoc.createElement('subParameter');
    rangeNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LDParameter');
    rangeNode.setAttribute('name', 'LDParameter');

    % from (Integer)
    fromNode = simDoc.createElement('subParameter');
    fromNode.setAttribute('classPath', 'java.lang.Integer');
    fromNode.setAttribute('name', 'from');
    v = simDoc.createElement('value');
    v.appendChild(simDoc.createTextNode(sprintf('%d', froms(ti))));
    fromNode.appendChild(v);
    rangeNode.appendChild(fromNode);

    % dummy distribution (only the function/probability is read for balking)
    distrNode = simDoc.createElement('subParameter');
    distrNode.setAttribute('classPath', 'jmt.engine.random.Exponential');
    distrNode.setAttribute('name', 'Exponential');
    rangeNode.appendChild(distrNode);

    distrParNode = simDoc.createElement('subParameter');
    distrParNode.setAttribute('classPath', 'jmt.engine.random.ExponentialPar');
    distrParNode.setAttribute('name', 'distrPar');
    lambdaNode = simDoc.createElement('subParameter');
    lambdaNode.setAttribute('classPath', 'java.lang.Double');
    lambdaNode.setAttribute('name', 'lambda');
    v = simDoc.createElement('value');
    v.appendChild(simDoc.createTextNode('1.0'));
    lambdaNode.appendChild(v);
    distrParNode.appendChild(lambdaNode);
    rangeNode.appendChild(distrParNode);

    % function (String) = balking probability for this range
    funcNode = simDoc.createElement('subParameter');
    funcNode.setAttribute('classPath', 'java.lang.String');
    funcNode.setAttribute('name', 'function');
    v = simDoc.createElement('value');
    v.appendChild(simDoc.createTextNode(sprintf('%.12f', probs(ti))));
    funcNode.appendChild(v);
    rangeNode.appendChild(funcNode);

    ldArrayNode.appendChild(rangeNode);
end

ldStrategyNode.appendChild(ldArrayNode);
balkingNode.appendChild(ldStrategyNode);

% priorityActivated flag (Parameter #2 of Balking)
prioNode = simDoc.createElement('subParameter');
prioNode.setAttribute('classPath', 'java.lang.Boolean');
prioNode.setAttribute('name', 'priorityActivated');
v = simDoc.createElement('value');
v.appendChild(simDoc.createTextNode('false'));
prioNode.appendChild(v);
balkingNode.appendChild(prioNode);
end
