% this macro will need refactoring to decouple the observation from the Model class
function [loggerBefore,loggerAfter] = linkAndLog(self, P, isNodeLogged, logPath)
% [LOGGERBEFORE,LOGGERAFTER] = LINKANDLOG(P, ISNODELOGGED, LOGPATH)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.resetStruct; % important to regenerate the sn with the loggers

if self.hasState
    line_warning(mfilename,'The network state cannot be initialized, or must be reset, before calling linkAndLog.');
end

if ~isempty(self.connections)
    line_warning(mfilename,'Network topology already instantiated. Calling resetNetwork automatically before adding loggers.\n');
    self.resetNetwork;
end
R = self.getNumberOfClasses;
Mnodes = self.getNumberOfNodes;

if Mnodes ~= numel(isNodeLogged)
    line_error(mfilename,'The size of the isNodeLogged vector does not match the number of nodes.');
end

isNodeLogged = [isNodeLogged(:)'];
if ~isempty(self.getSource)
    sinkIndex = self.getIndexSinkNode;
    if isNodeLogged(sinkIndex)
        line_warning(mfilename,'Sink station cannot be logged, ignoring.\n');
        isNodeLogged(sinkIndex) = false;
    end
end
if ~isempty(self.getSource)
    sourceIndex = self.getIndexSourceNode;
    if isNodeLogged(sourceIndex)
        line_warning(mfilename,'Source station cannot be logged, ignoring.\n');
        isNodeLogged(sourceIndex) = false;
    end
end

if nargin>=4 %exist('logPath','var')
    self.setLogPath(logPath);
else
    logPath = getLogPath(self);
end

loggerBefore = cell(1,0);
loggerAfter = cell(1,0);
for ind=1:Mnodes
    if isNodeLogged(ind)
        if ispc
            loggerBefore{end+1} = Logger(self,sprintf('Arv_%s',self.getNodeNames{ind}),[logPath,filesep,sprintf('%s-Arv.csv',self.getNodeNames{ind})]);
            for r=1:R
                loggerBefore{end}.setRouting(self.classes{r}, RoutingStrategy.RAND);
            end
        elseif isunix
            loggerBefore{end+1} = Logger(self,sprintf('Arv_%s',self.getNodeNames{ind}),[logPath,filesep,sprintf('%s-Arv.csv',self.getNodeNames{ind})]);
            for r=1:R
                loggerBefore{end}.setRouting(self.classes{r}, RoutingStrategy.RAND);
            end
        end
    end
end
for ind=1:Mnodes
    if isNodeLogged(ind)
        if ispc
            loggerAfter{end+1} = Logger(self,sprintf('Dep_%s',self.getNodeNames{ind}),[logPath,filesep,sprintf('%s-Dep.csv',self.getNodeNames{ind})]);
            for r=1:R
                loggerAfter{end}.setRouting(self.classes{r}, RoutingStrategy.RAND);
            end
        elseif isunix
            loggerAfter{end+1} = Logger(self,sprintf('Dep_%s',self.getNodeNames{ind}),[logPath,filesep,sprintf('%s-Dep.csv',self.getNodeNames{ind})]);
            for r=1:R
                loggerAfter{end}.setRouting(self.classes{r}, RoutingStrategy.RAND);
            end
        end
    end
end

Mnodesnew = 3*Mnodes;
newP = cellzeros(R,R,Mnodesnew,Mnodesnew);
for r=1:R
    for s=1:R
        for ind=1:Mnodes
            for jnd=1:Mnodes
                if P{r,s}(ind,jnd)>0
                    if isNodeLogged(ind) && isNodeLogged(jnd)
                        % link loggerArvi to loggerDepj
                        newP{r,s}(2*Mnodes+ind,Mnodes+jnd) = P{r,s}(ind,jnd);
                    elseif isNodeLogged(ind) && ~isNodeLogged(jnd)
                        % link logAi to j
                        newP{r,s}(2*Mnodes+ind,jnd) = P{r,s}(ind,jnd);
                    elseif ~isNodeLogged(ind) && isNodeLogged(jnd)
                        % link i to logBj
                        newP{r,s}(ind,Mnodes+jnd) = P{r,s}(ind,jnd);
                    else
                        % link i to j
                        newP{r,s}(ind,jnd) = P{r,s}(ind,jnd);
                    end
                end
            end
        end
        for ind=1:Mnodes
            if isNodeLogged(ind)
                newP{r,r}(Mnodes+ind,ind) = 1.0; % logBi -> i
                newP{r,r}(ind,2*Mnodes+ind) = 1.0; % i -> logAi
            end
        end
    end
end
for r=1:R
    for s=1:R
        idx = find(isNodeLogged);
        newP{r,s} = newP{r,s}([1:Mnodes,Mnodes+idx,2*Mnodes+idx],[1:Mnodes,Mnodes+idx,2*Mnodes+idx]);
    end
end
self.link(newP);
end
