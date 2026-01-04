function [simDoc, section] = saveCacheStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVECACHESTRATEGY(SIMDOC, SECTION, NODEIDX)
%
% Saves cache strategy parameters for JMT simulation
% This includes: maxItems, cacheCapacity, transition matrix,
% job classes, hit/miss classes, replacement policy, and popularity distributions

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
K = sn.nclasses;
nodeparam = sn.nodeparam{ind};

% Parameter 1: maxItems (Integer)
maxItemsNode = simDoc.createElement('parameter');
maxItemsNode.setAttribute('classPath', 'java.lang.Integer');
maxItemsNode.setAttribute('name', 'maxItems');
valNode = simDoc.createElement('value');
valNode.appendChild(simDoc.createTextNode(sprintf('%d', nodeparam.nitems)));
maxItemsNode.appendChild(valNode);
section.appendChild(maxItemsNode);

% Parameter 2: cacheCapacity (Integer[])
capNode = simDoc.createElement('parameter');
capNode.setAttribute('array', 'true');
capNode.setAttribute('classPath', 'java.lang.Integer');
capNode.setAttribute('name', 'cacheCapacity');
itemcap = nodeparam.itemcap;
for lvl = 1:length(itemcap)
    subParNode = simDoc.createElement('subParameter');
    subParNode.setAttribute('classPath', 'java.lang.Integer');
    subParNode.setAttribute('name', 'capacity');
    valNode = simDoc.createElement('value');
    valNode.appendChild(simDoc.createTextNode(sprintf('%d', itemcap(lvl))));
    subParNode.appendChild(valNode);
    capNode.appendChild(subParNode);
end
section.appendChild(capNode);

% Parameter 3: matrix (Object[] of Float[]) - Cache transition matrix
% For LRU: items always move to level 0 (front) on hit
% For FIFO/RR: items stay at their current level (identity matrix)
matrixNode = simDoc.createElement('parameter');
matrixNode.setAttribute('array', 'true');
matrixNode.setAttribute('classPath', 'java.lang.Object');
matrixNode.setAttribute('name', 'matrix');
nLevels = length(itemcap);
for lvl1 = 1:nLevels
    rowNode = simDoc.createElement('subParameter');
    rowNode.setAttribute('array', 'true');
    rowNode.setAttribute('classPath', 'java.lang.Float');
    rowNode.setAttribute('name', 'row');
    for lvl2 = 1:nLevels
        cellNode = simDoc.createElement('subParameter');
        cellNode.setAttribute('classPath', 'java.lang.Float');
        cellNode.setAttribute('name', 'cell');
        valNode = simDoc.createElement('value');
        % For LRU, items move to next level on hit (climb up)
        % Hit at list i → move to list i+1; hit at last list → stay
        % For other policies, use identity matrix (stay in place)
        if nodeparam.replacestrat == ReplacementStrategy.LRU
            if lvl1 < nLevels  % not the last list
                % Move to next level (lvl1+1)
                if lvl2 == lvl1 + 1
                    valNode.appendChild(simDoc.createTextNode('1.0'));
                else
                    valNode.appendChild(simDoc.createTextNode('0.0'));
                end
            else  % last list: stay in place
                if lvl2 == nLevels
                    valNode.appendChild(simDoc.createTextNode('1.0'));
                else
                    valNode.appendChild(simDoc.createTextNode('0.0'));
                end
            end
        else
            % FIFO, RR: identity matrix
            if lvl1 == lvl2
                valNode.appendChild(simDoc.createTextNode('1.0'));
            else
                valNode.appendChild(simDoc.createTextNode('0.0'));
            end
        end
        cellNode.appendChild(valNode);
        rowNode.appendChild(cellNode);
    end
    matrixNode.appendChild(rowNode);
end
section.appendChild(matrixNode);

% Parameter 4: jobClasses (JobClass[])
% Classes that access the cache (have hit/miss behavior)
jobClassesNode = simDoc.createElement('parameter');
jobClassesNode.setAttribute('array', 'true');
jobClassesNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
jobClassesNode.setAttribute('name', 'jobClasses');
for r = 1:K
    if nodeparam.hitclass(r) > 0 || nodeparam.missclass(r) > 0
        subParNode = simDoc.createElement('subParameter');
        subParNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
        subParNode.setAttribute('name', 'jobClass');
        valNode = simDoc.createElement('value');
        valNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
        subParNode.appendChild(valNode);
        jobClassesNode.appendChild(subParNode);
    end
end
section.appendChild(jobClassesNode);

% Parameter 5: hitClasses (JobClass[])
hitClassesNode = simDoc.createElement('parameter');
hitClassesNode.setAttribute('array', 'true');
hitClassesNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
hitClassesNode.setAttribute('name', 'hitClasses');
for r = 1:K
    if nodeparam.hitclass(r) > 0
        subParNode = simDoc.createElement('subParameter');
        subParNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
        subParNode.setAttribute('name', 'hitClass');
        valNode = simDoc.createElement('value');
        hitIdx = nodeparam.hitclass(r);
        valNode.appendChild(simDoc.createTextNode(sn.classnames{hitIdx}));
        subParNode.appendChild(valNode);
        hitClassesNode.appendChild(subParNode);
    end
end
section.appendChild(hitClassesNode);

% Parameter 6: missClasses (JobClass[])
missClassesNode = simDoc.createElement('parameter');
missClassesNode.setAttribute('array', 'true');
missClassesNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
missClassesNode.setAttribute('name', 'missClasses');
for r = 1:K
    if nodeparam.missclass(r) > 0
        subParNode = simDoc.createElement('subParameter');
        subParNode.setAttribute('classPath', 'jmt.engine.QueueNet.JobClass');
        subParNode.setAttribute('name', 'missClass');
        valNode = simDoc.createElement('value');
        missIdx = nodeparam.missclass(r);
        valNode.appendChild(simDoc.createTextNode(sn.classnames{missIdx}));
        subParNode.appendChild(valNode);
        missClassesNode.appendChild(subParNode);
    end
end
section.appendChild(missClassesNode);

% Parameter 7: replacePolicy (CacheStrategy subclass)
replaceNode = simDoc.createElement('parameter');
replaceNode.setAttribute('name', 'replacePolicy');
% Set classPath directly to the specific cache strategy class
switch nodeparam.replacestrat
    case ReplacementStrategy.LRU
        replaceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.CacheStrategies.LRUCache');
    case ReplacementStrategy.FIFO
        replaceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.CacheStrategies.FIFOCache');
    case ReplacementStrategy.SFIFO
        % Strict FIFO also uses FIFOCache in JMT (see FIFOCache.java comment)
        replaceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.CacheStrategies.FIFOCache');
    case ReplacementStrategy.RR
        % Random Replacement maps to RandomCache in JMT
        replaceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.CacheStrategies.RandomCache');
    otherwise
        % Default to LRU
        replaceNode.setAttribute('classPath', 'jmt.engine.NetStrategies.CacheStrategies.LRUCache');
end
section.appendChild(replaceNode);

% Parameter 8: popularity (DiscreteDistribution[]) - Popularity distributions
% Access the original Cache node to get distribution objects
cacheNode = self.model.nodes{ind};
popNode = simDoc.createElement('parameter');
popNode.setAttribute('array', 'true');
popNode.setAttribute('classPath', 'jmt.engine.random.discrete.DiscreteDistribution');
popNode.setAttribute('name', 'popularity');
for r = 1:K
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    popNode.appendChild(refClassNode);

    subParNode = simDoc.createElement('subParameter');
    if ~iscell(nodeparam.pread) || length(nodeparam.pread) < r || isnan(nodeparam.pread{r}(1))
        % No popularity distribution for this class
        subParNode.setAttribute('classPath', 'jmt.engine.random.discrete.DiscreteDistribution');
        subParNode.setAttribute('name', 'null');
        valNode = simDoc.createElement('value');
        valNode.appendChild(simDoc.createTextNode('null'));
        subParNode.appendChild(valNode);
    else
        % Get the original distribution from the Cache node
        % popularity is indexed by {itemclass, jobclass}
        distrib = [];
        if isprop(cacheNode, 'popularity') && ~isempty(cacheNode.popularity)
            % Find the distribution for this class
            for iItem = 1:size(cacheNode.popularity, 1)
                if size(cacheNode.popularity, 2) >= r && ~isempty(cacheNode.popularity{iItem, r})
                    distrib = cacheNode.popularity{iItem, r};
                    break;
                end
            end
        end

        if ~isempty(distrib) && isa(distrib, 'Zipf')
            % Use JMT's native Zipf distribution
            % Structure matches JLINE Java: single subParameter with classPath=Zipf, name=popularity
            % containing alpha and numberOfElements as child subParameters
            subParNode.setAttribute('classPath', 'jmt.engine.random.discrete.Zipf');
            subParNode.setAttribute('name', 'popularity');

            % Get Zipf parameters: alpha (shape) and n (number of elements)
            alpha = distrib.getParam(3).paramValue;  % 's' parameter
            n = distrib.getParam(4).paramValue;      % 'n' parameter

            % alpha parameter (comes first in Java implementation)
            alphaNode = simDoc.createElement('subParameter');
            alphaNode.setAttribute('classPath', 'java.lang.Double');
            alphaNode.setAttribute('name', 'alpha');
            valNode = simDoc.createElement('value');
            valNode.appendChild(simDoc.createTextNode(sprintf('%.12f', alpha)));
            alphaNode.appendChild(valNode);
            subParNode.appendChild(alphaNode);

            % numberOfElements parameter
            numElemNode = simDoc.createElement('subParameter');
            numElemNode.setAttribute('classPath', 'java.lang.Integer');
            numElemNode.setAttribute('name', 'numberOfElements');
            valNode = simDoc.createElement('value');
            valNode.appendChild(simDoc.createTextNode(sprintf('%d', n)));
            numElemNode.appendChild(valNode);
            subParNode.appendChild(numElemNode);
        elseif ~isempty(distrib) && isa(distrib, 'DiscreteSampler')
            % Use JMT's Uniform distribution for DiscreteSampler
            % This maps the discrete sampler to a uniform distribution over its support
            subParNode.setAttribute('classPath', 'jmt.engine.random.discrete.Uniform');
            subParNode.setAttribute('name', 'popularity');

            % Get min/max from the distribution's support
            minVal = distrib.support(1);
            maxVal = distrib.support(2);

            % min parameter
            minNode = simDoc.createElement('subParameter');
            minNode.setAttribute('classPath', 'java.lang.Integer');
            minNode.setAttribute('name', 'min');
            valNode = simDoc.createElement('value');
            valNode.appendChild(simDoc.createTextNode(sprintf('%d', minVal)));
            minNode.appendChild(valNode);
            subParNode.appendChild(minNode);

            % max parameter
            maxNode = simDoc.createElement('subParameter');
            maxNode.setAttribute('classPath', 'java.lang.Integer');
            maxNode.setAttribute('name', 'max');
            valNode = simDoc.createElement('value');
            valNode.appendChild(simDoc.createTextNode(sprintf('%d', maxVal)));
            maxNode.appendChild(valNode);
            subParNode.appendChild(maxNode);
        else
            % Fallback: use null for unsupported distributions
            % JMT doesn't have a generic empirical discrete distribution
            subParNode.setAttribute('classPath', 'jmt.engine.random.discrete.DiscreteDistribution');
            subParNode.setAttribute('name', 'null');
            valNode = simDoc.createElement('value');
            valNode.appendChild(simDoc.createTextNode('null'));
            subParNode.appendChild(valNode);
        end
    end
    popNode.appendChild(subParNode);
end
section.appendChild(popNode);

end
