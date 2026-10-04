function [RD,logData] = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

sn = self.getStruct;
if GlobalConstants.DummyMode
    RD = cell(sn.nstations,sn.nclasses);
    logData = NaN;
    return
end

if nargin<2 %~exist('R','var')
    R = getAvgRespTHandles(self);
end
RD = cell(sn.nstations, sn.nclasses);
% The steady-state marginal only seeds initFromMarginal for the logging run
% below, so it belongs to the same native-JMT computation and takes the same
% lang. Left on the caller's lang it is the ONE solve here that is not forced.
% Pinned to 'matlab' because the logging run that follows is this wrapper's
% own (it rewrites the JSIM document with per-node loggers), so the seed must
% come from the same implementation whatever the caller's lang says.
origLang = self.options.lang;
self.options.lang = 'matlab';
QN = getAvgQLen(self); % steady-state qlen
self.options.lang = origLang;
n = QN;
for r=1:sn.nclasses
    if isinf(sn.njobs(r))
        n(:,r) = floor(QN(:,r));
    else
        n(:,r) = floor(QN(:,r));
        if sum(n(:,r)) < sn.njobs(r)
            imax = maxpos(n(:,r)); % put jobs on the bottleneck
            n(imax,r) = n(imax,r) + sn.njobs(r) - sum(n(:,r));
        end
    end
end
cdfmodel = self.model.copy;
cdfmodel.resetNetwork;
cdfmodel.reset;
isNodeClassLogged = false(cdfmodel.getNumberOfNodes, cdfmodel.getNumberOfClasses);
for i= 1:cdfmodel.getNumberOfStations
    for r=1:cdfmodel.getNumberOfClasses
        if ~R{i,r}.disabled
            ni = self.model.getNodeIndex(cdfmodel.getStationNames{i});
            isNodeClassLogged(ni,r) = true;
        end
    end
end
Plinked = sn.rtorig;
isNodeLogged = max(isNodeClassLogged,[],2);
logpath = lineTempDir;
cdfmodel.linkAndLog(Plinked, isNodeLogged, logpath);
cdfmodel.initFromMarginal(n);
cdfOptions = self.getOptions;
cdfOptions.lang = 'matlab'; % CDF computation requires log files from native JMT
SolverJMT(cdfmodel, cdfOptions).getAvg(); % log data
logData = SolverJMT.parseLogs(cdfmodel, isNodeLogged, MetricType.toText(MetricType.RespT));
% from here convert from nodes in logData to stations
for i= 1:cdfmodel.getNumberOfStations
    ni = cdfmodel.getNodeIndex(cdfmodel.getStationNames{i});
    for r=1:cdfmodel.getNumberOfClasses
        if isNodeClassLogged(ni,r)
            if ~isempty(logData{ni,r}) && ~isempty(logData{ni,r}.RespT)
                [F,X] = ecdf(logData{ni,r}.RespT);
                RD{i,r} = [F,X];
            end
        end
    end
end
end