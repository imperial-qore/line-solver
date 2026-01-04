function [simDoc, section] = saveServiceStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVESERVICESTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategy');
strategyNode.setAttribute('name', 'ServiceStrategy');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
i = sn.nodeToStation(ind);

for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode2 = simDoc.createElement('refClass');
    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
    strategyNode.appendChild(refClassNode2);    
    serviceTimeStrategyNode = simDoc.createElement('subParameter');
    if sn.procid(i,r) == ProcessType.DISABLED
        serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy');
        serviceTimeStrategyNode.setAttribute('name', 'DisabledServiceTimeStrategy');
    elseif sn.procid(i,r) == ProcessType.IMMEDIATE
        serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy');
        serviceTimeStrategyNode.setAttribute('name', 'ZeroServiceTimeStrategy');
    elseif sn.procid(i,r) == ProcessType.PH || sn.procid(i,r) == ProcessType.APH || sn.procid(i,r) == ProcessType.COXIAN || (sn.phases(i,r)>2 && sn.procid(i,r) == ProcessType.HYPEREXP) %|| (sn.phases(i,r)>2 && sn.procid(i,r) == ProcessType.COXIAN) || (sn.phases(i,r)>2 && sn.procid(i,r) == ProcessType.HYPEREXP)
        serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
        serviceTimeStrategyNode.setAttribute('name', 'ServiceTimeStrategy');
        distributionNode = simDoc.createElement('subParameter');
        distributionNode.setAttribute('classPath', 'jmt.engine.random.PhaseTypeDistr');
        distributionNode.setAttribute('name', 'Phase-Type');
        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', 'jmt.engine.random.PhaseTypePar');
        distrParNode.setAttribute('name', 'distrPar');
        
        subParNodeAlpha = simDoc.createElement('subParameter');
        subParNodeAlpha.setAttribute('array', 'true');
        subParNodeAlpha.setAttribute('classPath', 'java.lang.Object');
        subParNodeAlpha.setAttribute('name', 'alpha');
        subParNodeAlphaVec = simDoc.createElement('subParameter');
        subParNodeAlphaVec.setAttribute('array', 'true');
        subParNodeAlphaVec.setAttribute('classPath', 'java.lang.Object');
        subParNodeAlphaVec.setAttribute('name', 'vector');
        PH = sn.proc{i}{r};
        alpha = abs(sn.pie{i}{r});
        for k=1:sn.phases(i,r)
            subParNodeAlphaElem = simDoc.createElement('subParameter');
            subParNodeAlphaElem.setAttribute('classPath', 'java.lang.Double');
            subParNodeAlphaElem.setAttribute('name', 'entry');
            subParValue = simDoc.createElement('value');
            subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',alpha(k))));
            subParNodeAlphaElem.appendChild(subParValue);
            subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
        end
        
        subParNodeT = simDoc.createElement('subParameter');
        subParNodeT.setAttribute('array', 'true');
        subParNodeT.setAttribute('classPath', 'java.lang.Object');
        subParNodeT.setAttribute('name', 'T');
        T = PH{1};
        for k=1:sn.phases(i,r)
            subParNodeTvec = simDoc.createElement('subParameter');
            subParNodeTvec.setAttribute('array', 'true');
            subParNodeTvec.setAttribute('classPath', 'java.lang.Object');
            subParNodeTvec.setAttribute('name', 'vector');
            for j=1:sn.phases(i,r)
                subParNodeTElem = simDoc.createElement('subParameter');
                subParNodeTElem.setAttribute('classPath', 'java.lang.Double');
                subParNodeTElem.setAttribute('name', 'entry');
                subParValue = simDoc.createElement('value');
                if k==j
                    subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',-abs(T(k,j)))));
                else
                    subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',abs(T(k,j)))));
                end
                subParNodeTElem.appendChild(subParValue);
                subParNodeTvec.appendChild(subParNodeTElem);
            end
            subParNodeT.appendChild(subParNodeTvec);
        end
        
        subParNodeAlpha.appendChild(subParNodeAlphaVec);
        distrParNode.appendChild(subParNodeAlpha);
        distrParNode.appendChild(subParNodeT);
        serviceTimeStrategyNode.appendChild(distributionNode);
        serviceTimeStrategyNode.appendChild(distrParNode);
        
    elseif sn.procid(i,r) == ProcessType.MAP || sn.procid(i,r) == ProcessType.MMPP2
        serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
        serviceTimeStrategyNode.setAttribute('name', 'ServiceTimeStrategy');
        distributionNode = simDoc.createElement('subParameter');
        distributionNode.setAttribute('classPath', 'jmt.engine.random.MAPDistr');
        distributionNode.setAttribute('name', 'Burst (MAP)');
        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', 'jmt.engine.random.MAPPar');
        distrParNode.setAttribute('name', 'distrPar');
        
        MAP = sn.proc{i}{r};
        
        subParNodeD0 = simDoc.createElement('subParameter');
        subParNodeD0.setAttribute('array', 'true');
        subParNodeD0.setAttribute('classPath', 'java.lang.Object');
        subParNodeD0.setAttribute('name', 'D0');
        D0 = MAP{1};
        for k=1:sn.phases(i,r)
            subParNodeD0vec = simDoc.createElement('subParameter');
            subParNodeD0vec.setAttribute('array', 'true');
            subParNodeD0vec.setAttribute('classPath', 'java.lang.Object');
            subParNodeD0vec.setAttribute('name', 'vector');
            for j=1:sn.phases(i,r)
                subParNodeD0Elem = simDoc.createElement('subParameter');
                subParNodeD0Elem.setAttribute('classPath', 'java.lang.Double');
                subParNodeD0Elem.setAttribute('name', 'entry');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',D0(k,j))));
                subParNodeD0Elem.appendChild(subParValue);
                subParNodeD0vec.appendChild(subParNodeD0Elem);
            end
            subParNodeD0.appendChild(subParNodeD0vec);
        end
        distrParNode.appendChild(subParNodeD0);
        
        subParNodeD1 = simDoc.createElement('subParameter');
        subParNodeD1.setAttribute('array', 'true');
        subParNodeD1.setAttribute('classPath', 'java.lang.Object');
        subParNodeD1.setAttribute('name', 'D1');
        D1 = MAP{2};
        for k=1:sn.phases(i,r)
            subParNodeD1vec = simDoc.createElement('subParameter');
            subParNodeD1vec.setAttribute('array', 'true');
            subParNodeD1vec.setAttribute('classPath', 'java.lang.Object');
            subParNodeD1vec.setAttribute('name', 'vector');
            for j=1:sn.phases(i,r)
                subParNodeD1Elem = simDoc.createElement('subParameter');
                subParNodeD1Elem.setAttribute('classPath', 'java.lang.Double');
                subParNodeD1Elem.setAttribute('name', 'entry');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',D1(k,j))));
                subParNodeD1Elem.appendChild(subParValue);
                subParNodeD1vec.appendChild(subParNodeD1Elem);
            end
            subParNodeD1.appendChild(subParNodeD1vec);
        end
        distrParNode.appendChild(subParNodeD1);
        serviceTimeStrategyNode.appendChild(distributionNode);
        serviceTimeStrategyNode.appendChild(distrParNode);
    else
        serviceTimeStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
        serviceTimeStrategyNode.setAttribute('name', 'ServiceTimeStrategy');
        
        distributionNode = simDoc.createElement('subParameter');
        switch sn.procid(i,r)
            case ProcessType.DET
                javaClass = 'jmt.engine.random.DeterministicDistr';
                javaParClass = 'jmt.engine.random.DeterministicDistrPar';
            case ProcessType.COXIAN
                javaClass = 'jmt.engine.random.CoxianDistr';
                javaParClass = 'jmt.engine.random.CoxianPar';
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
            case ProcessType.MMPP2
                javaClass = 'jmt.engine.random.MMPP2Distr';
                javaParClass = 'jmt.engine.random.MMPP2Par';
            case {ProcessType.REPLAYER, ProcessType.TRACE}
                javaClass = 'jmt.engine.random.Replayer';
                javaParClass = 'jmt.engine.random.ReplayerPar';
        end
        distributionNode.setAttribute('classPath', javaClass);
        switch sn.procid(i,r)
            case {ProcessType.REPLAYER, ProcessType.TRACE}
                distributionNode.setAttribute('name', 'Replayer');
            case ProcessType.EXP
                distributionNode.setAttribute('name', 'Exponential');
            case ProcessType.HYPEREXP
                distributionNode.setAttribute('name', 'Hyperexponential');
            otherwise
                distributionNode.setAttribute('name', ProcessType.toText(ProcessType.fromId(sn.procid(i,r))));
        end
        serviceTimeStrategyNode.appendChild(distributionNode);
        
        distrParNode = simDoc.createElement('subParameter');
        distrParNode.setAttribute('classPath', javaParClass);
        distrParNode.setAttribute('name', 'distrPar');
                
        switch sn.procid(i,r)
            case ProcessType.DET
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 't');
                subParValue = simDoc.createElement('value');                
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',1/sn.rates(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.EXP
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'lambda');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.rates(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.HYPEREXP
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'p');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.pie{i}{r}(1))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'lambda1');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',-sn.proc{i}{r}{1}(1,1))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'lambda2');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',-sn.proc{i}{r}{1}(2,2))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.ERLANG
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.rates(i,r)*sn.phases(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Long');
                subParNodeAlpha.setAttribute('name', 'r');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%d',sn.phases(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.MMPP2
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'lambda0');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.proc{i}{r}{2}(1,1))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'lambda1');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.proc{i}{r}{2}(2,2))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'sigma0');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.proc{i}{r}{1}(1,2))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'sigma1');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.proc{i}{r}{1}(2,1))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.GAMMA
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',1/sn.scv(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'beta');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sn.scv(i,r)/sn.rates(i,r))));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.PARETO
                shape = sqrt(1+1/sn.scv(i,r))+1;
                scale = 1/sn.rates(i,r) *  (shape - 1) / shape;
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha'); % shape
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',shape)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'k'); % scale
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',scale)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.WEIBULL
                c = sqrt(sn.scv(i,r));
                rval = c^(-1.086); % Justus approximation (1976)
                alpha =  1/sn.rates(i,r) / gamma(1+1/rval);                
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'alpha'); % shape
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',alpha)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'r'); % scale
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',rval)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.LOGNORMAL
                c = sqrt(sn.scv(i,r));
                mu = log(1/sn.rates(i,r)  / sqrt(c*c + 1));
                sigma = sqrt(log(c*c + 1));                            
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'mu'); % shape
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',mu)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'sigma'); % scale
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',sigma)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            case ProcessType.UNIFORM
                maxVal = ((sqrt(12*sn.scv(i,r)/sn.rates(i,r)^2))+2/sn.rates(i,r))/2;
                minVal = 2/sn.rates(i,r)-maxVal;
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'min'); % shape
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',minVal)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlpha.setAttribute('name', 'max'); % scale
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',maxVal)));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);                
            case {ProcessType.REPLAYER, ProcessType.TRACE}
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', 'java.lang.String');
                subParNodeAlpha.setAttribute('name', 'fileName');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sn.nodeparam{ind}{r}.fileName));
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
        end        
        serviceTimeStrategyNode.appendChild(distrParNode);
    end
    strategyNode.appendChild(serviceTimeStrategyNode);
    section.appendChild(strategyNode);
end
end
