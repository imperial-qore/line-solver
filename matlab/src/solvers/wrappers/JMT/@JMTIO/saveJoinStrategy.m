function [simDoc, section] = saveJoinStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEJOINSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategy');
strategyNode.setAttribute('name', 'JoinStrategy');

sn = self.getStruct;
numOfClasses = sn.nclasses;
for r=1:numOfClasses    
    switch sn.nodeparam{ind}.joinStrategy{r}
        case JoinStrategy.STD
            refClassNode2 = simDoc.createElement('refClass');
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
            strategyNode.appendChild(refClassNode2);
            
            joinStrategyNode = simDoc.createElement('subParameter');
            joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.NormalJoin');
            joinStrategyNode.setAttribute('name', 'Standard Join');
            reqNode = simDoc.createElement('subParameter');
            reqNode.setAttribute('classPath', 'java.lang.Integer');
            reqNode.setAttribute('name', 'numRequired');
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.fanIn{r})));
            reqNode.appendChild(valueNode);
            joinStrategyNode.appendChild(reqNode);
            strategyNode.appendChild(joinStrategyNode);
            section.appendChild(strategyNode);
        case JoinStrategy.PARTIAL
            refClassNode2 = simDoc.createElement('refClass');
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
            strategyNode.appendChild(refClassNode2);
            
            joinStrategyNode = simDoc.createElement('subParameter');
            joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.PartialJoin');
            joinStrategyNode.setAttribute('name', 'Quorum');
            reqNode = simDoc.createElement('subParameter');
            reqNode.setAttribute('classPath', 'java.lang.Integer');
            reqNode.setAttribute('name', 'numRequired');
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.joinRequired{r})));
            reqNode.appendChild(valueNode);
            joinStrategyNode.appendChild(reqNode);
            strategyNode.appendChild(joinStrategyNode);
            section.appendChild(strategyNode);
        % case JoinStrategy.Guard
        %     refClassNode2 = simDoc.createElement('refClass');
        %     refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
        %     strategyNode.appendChild(refClassNode2);
        % 
        %     joinStrategyNode = simDoc.createElement('subParameter');
        %     joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.PartialJoin');
        %     joinStrategyNode.setAttribute('name', 'Quorum');
        %     reqNode = simDoc.createElement('subParameter');
        %     reqNode.setAttribute('classPath', 'java.lang.Integer');
        %     reqNode.setAttribute('name', 'numRequired');
        %     valueNode = simDoc.createElement('value');
        %     valueNode.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.joinRequired{r})));
        %     reqNode.appendChild(valueNode);
        %     joinStrategyNode.appendChild(reqNode);
        %     strategyNode.appendChild(joinStrategyNode);
        %     section.appendChild(strategyNode);
    end
end
end
