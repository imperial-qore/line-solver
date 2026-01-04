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

    impatienceStrategyNode = simDoc.createElement('subParameter');
    impatienceStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Impatience');
    impatienceStrategyNode.setAttribute('name', 'Impatience');

    % Check if impatience is defined for this station-class pair
    % impatienceType is a matrix indexed by (station, class)
    % i > 0 check: nodeToStation returns 0 for non-station nodes (Fork, Join, etc.)
    hasImpatience = false;
    if i > 0 && isfield(sn, 'impatienceType') && ~isempty(sn.impatienceType)
        if sn.impatienceType(i, r) ~= 0
            hasImpatience = true;
        end
    end

    if ~hasImpatience
        % No impatience defined - use null
        subParValue = simDoc.createElement('value');
        subParValue.appendChild(simDoc.createTextNode('null'));
        impatienceStrategyNode.appendChild(subParValue);
    else
        % Impatience is defined - generate distribution XML
        procType = sn.impatienceType(i, r);

        distributionNode = simDoc.createElement('subParameter');

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
            case {ProcessType.PH, ProcessType.APH, ProcessType.COXIAN}
                javaClass = 'jmt.engine.random.PhaseTypeDistr';
                javaParClass = 'jmt.engine.random.PhaseTypePar';
            otherwise
                line_error(mfilename, sprintf('Unsupported impatience distribution type: %s', ProcessType.toText(procType)));
        end

        distributionNode.setAttribute('classPath', javaClass);
        switch procType
            case ProcessType.EXP
                distributionNode.setAttribute('name', 'Exponential');
            case ProcessType.HYPEREXP
                distributionNode.setAttribute('name', 'Hyperexponential');
            case {ProcessType.PH, ProcessType.APH, ProcessType.COXIAN}
                distributionNode.setAttribute('name', 'Phase-Type');
            otherwise
                distributionNode.setAttribute('name', ProcessType.toText(procType));
        end
        impatienceStrategyNode.appendChild(distributionNode);

        % Create distribution parameters
        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', javaParClass);
        distrParNode.setAttribute('name', 'distrPar');

        % Get impatience parameters
        impatienceMu = sn.impatienceMu(i, r);
        impatiencePhi = sn.impatiencePhi(i, r);

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
                impatienceProc = sn.impatienceProc(i, r);
                impatiencePie = sn.impatiencePie(i, r);
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

            case {ProcessType.PH, ProcessType.APH, ProcessType.COXIAN}
                impatienceProc = sn.impatienceProc(i, r);
                impatiencePie = sn.impatiencePie(i, r);
                phases = sn.impatiencePhases(i, r);
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
        end

        impatienceStrategyNode.appendChild(distrParNode);
    end

    impatienceNode.appendChild(impatienceStrategyNode);
end

section.appendChild(impatienceNode);
end
