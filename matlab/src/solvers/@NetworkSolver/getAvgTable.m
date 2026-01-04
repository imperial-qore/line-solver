function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,TT,AT] = GETAVGTABLE(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Return table of average station metrics
%
% When confint option is enabled, values are displayed as 'mean ± ci_halfwidth'
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if GlobalConstants.DummyMode
    [AvgTable, QT, UT, RT, TT, WT, AT] = deal(Table());
    AvgTable = IndexedTable(AvgTable);
    return
end

% Check if confidence intervals should be displayed
[confintEnabled, ~] = Solver.parseConfInt(self.options.confint);
hasCI = confintEnabled && self.hasResults && isfield(self.result, 'Avg') && isfield(self.result.Avg, 'QCI');

if strcmp(self.options.lang,'java') && ~strcmp(self.name,'SolverDES')
    sn = self.model.getStruct;
    M = sn.nstations;
    R = sn.nclasses;
    [QN,UN,RN,TN,AN,WN] = self.getAvg([],[],[],[],[],[]);

    % Update hasCI after getAvg call
    hasCI = confintEnabled && self.hasResults && isfield(self.result, 'Avg') && isfield(self.result.Avg, 'QCI');

    [Qval, Uval, Rval, Tval, Wval, Aval] = deal([]);
    [QCIval, UCIval, RCIval, TCIval, WCIval, ACIval] = deal([]);
    JobClass = {};
    Station = {};
    for ist=1:M
        for k=1:R
            if any(sum([QN(ist,k),UN(ist,k),RN(ist,k),TN(ist,k),AN(ist,k),WN(ist,k)])>0)
                c = sn.chains(:,k)>0;
                inchain = sn.inchain{c};
                JobClass{end+1,1} = sn.classnames{k};
                Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)};
                Qval(end+1) = QN(ist,k);
                Uval(end+1) = UN(ist,k);
                Rval(end+1) = RN(ist,k);
                Tval(end+1) = TN(ist,k);
                Wval(end+1) = WN(ist,k);
                Aval(end+1) = AN(ist,k);
                if hasCI
                    QCIval(end+1) = getCIValue(self.result.Avg, 'QCI', ist, k);
                    UCIval(end+1) = getCIValue(self.result.Avg, 'UCI', ist, k);
                    RCIval(end+1) = getCIValue(self.result.Avg, 'RCI', ist, k);
                    TCIval(end+1) = getCIValue(self.result.Avg, 'TCI', ist, k);
                    WCIval(end+1) = getCIValue(self.result.Avg, 'WCI', ist, k);
                    ACIval(end+1) = getCIValue(self.result.Avg, 'ACI', ist, k);
                end
            end
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);

    if hasCI
        % Format values with CI as strings
        QLen = formatWithCI(Qval(:), QCIval(:));
        Util = formatWithCI(Uval(:), UCIval(:));
        RespT = formatWithCI(Rval(:), RCIval(:));
        ResidT = formatWithCI(Wval(:), WCIval(:));
        Tput = formatWithCI(Tval(:), TCIval(:));
        ArvR = formatWithCI(Aval(:), ACIval(:));
    else
        QLen = Qval(:);
        Util = Uval(:);
        RespT = Rval(:);
        ResidT = Wval(:);
        Tput = Tval(:);
        ArvR = Aval(:);
    end
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput);
    AvgTable = IndexedTable(AvgTable);
    return
end

sn = getStruct(self);

if nargin<8
    keepDisabled = false;
elseif isempty(Q) && isempty(U) && isempty(R) && isempty(T) && isempty(A) && isempty(W)
    [Q,U,R,T,A,W] = getAvgHandles(self);
end

M = sn.nstations;
K = sn.nclasses;
if nargin == 2
    if islogical(Q)
        keepDisabled = Q;
    elseif iscell(Q) && ~isempty(Q)
        param = Q;
        keepDisabled = param{5};
        % case where varargin is passed as input
    end
    [Q,U,R,T,A,W] = getAvgHandles(self);
elseif nargin == 1
    [Q,U,R,T,A,W] = getAvgHandles(self);
end
if isfinite(self.getOptions.timespan(2))
    [Qt,Ut,Tt] = getTranHandles(self);
    [QNt,UNt,TNt] = self.getTranAvg(Qt,Ut,Tt);
    QN = cellfun(@(c) c.metric(end),QNt);
    UN = cellfun(@(c) c.metric(end),UNt);
    TN = cellfun(@(c) c.metric(end),TNt);
    [RN, WN, AN] = deal(zeros(size(QN)));
else
    [QN,UN,RN,TN,AN,WN] = self.getAvg(Q,U,R,T,A,W);
end

% Update hasCI after getAvg call
hasCI = confintEnabled && self.hasResults && isfield(self.result, 'Avg') && isfield(self.result.Avg, 'QCI');

% this is required because getAvg can alter the chain structure in the
% presence of caches
sn = self.model.getStruct;

if isempty(QN)
    [AvgTable, QT, UT, RT, TT, AT, WT] = deal(Table());
elseif ~keepDisabled
    [Qval, Uval, Rval, Tval, Wval, Aval] = deal([]);
    [QCIval, UCIval, RCIval, TCIval, WCIval, ACIval] = deal([]);
    JobClass = {};
    Station = {};
    for ist=1:M
        for k=1:K
            if any(sum([QN(ist,k),UN(ist,k),RN(ist,k),TN(ist,k),AN(ist,k),WN(ist,k)])>0)
                c = sn.chains(:,k)>0;
                inchain = sn.inchain{c};
                JobClass{end+1,1} = char(Q{ist,k}.class.name);
                Station{end+1,1} = char(Q{ist,k}.station.name);
                Qval(end+1) = QN(ist,k);
                Uval(end+1) = UN(ist,k);
                Rval(end+1) = RN(ist,k);
                Wval(end+1) = WN(ist,k);
                Tval(end+1) = TN(ist,k);
                Aval(end+1) = AN(ist,k);
                if hasCI
                    QCIval(end+1) = getCIValue(self.result.Avg, 'QCI', ist, k);
                    UCIval(end+1) = getCIValue(self.result.Avg, 'UCI', ist, k);
                    RCIval(end+1) = getCIValue(self.result.Avg, 'RCI', ist, k);
                    WCIval(end+1) = getCIValue(self.result.Avg, 'WCI', ist, k);
                    TCIval(end+1) = getCIValue(self.result.Avg, 'TCI', ist, k);
                    ACIval(end+1) = getCIValue(self.result.Avg, 'ACI', ist, k);
                end
            end
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);

    if hasCI
        % Format values with CI as strings
        QLen = formatWithCI(Qval(:), QCIval(:));
        Util = formatWithCI(Uval(:), UCIval(:));
        RespT = formatWithCI(Rval(:), RCIval(:));
        ResidT = formatWithCI(Wval(:), WCIval(:));
        Tput = formatWithCI(Tval(:), TCIval(:));
        ArvR = formatWithCI(Aval(:), ACIval(:));
    else
        QLen = Qval(:);
        Util = Uval(:);
        RespT = Rval(:);
        ResidT = Wval(:);
        Tput = Tval(:);
        ArvR = Aval(:);
    end
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput);
else
    [Qval, Uval, Rval, Wval, Tval, Aval] = deal(zeros(M, K));
    [QCIval, UCIval, RCIval, WCIval, TCIval, ACIval] = deal(zeros(M, K));
    JobClass = cell(K*M,1);
    Station = cell(K*M,1);
    for ist=1:M
        for k=1:K
            JobClass{(ist-1)*K+k} = Q{ist,k}.class.name;
            Station{(ist-1)*K+k} = Q{ist,k}.station.name;
            Qval((ist-1)*K+k) = QN(ist,k);
            Uval((ist-1)*K+k) = UN(ist,k);
            Rval((ist-1)*K+k) = RN(ist,k);
            Wval((ist-1)*K+k) = WN(ist,k);
            Tval((ist-1)*K+k) = TN(ist,k);
            Aval((ist-1)*K+k) = AN(ist,k);
            if hasCI
                QCIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'QCI', ist, k);
                UCIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'UCI', ist, k);
                RCIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'RCI', ist, k);
                WCIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'WCI', ist, k);
                TCIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'TCI', ist, k);
                ACIval((ist-1)*K+k) = getCIValue(self.result.Avg, 'ACI', ist, k);
            end
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);

    if hasCI
        % Format values with CI as strings
        QLen = formatWithCI(Qval(:), QCIval(:));
        Util = formatWithCI(Uval(:), UCIval(:));
        RespT = formatWithCI(Rval(:), RCIval(:));
        ResidT = formatWithCI(Wval(:), WCIval(:));
        Tput = formatWithCI(Tval(:), TCIval(:));
        ArvR = formatWithCI(Aval(:), ACIval(:));
    else
        QLen = Qval(:);
        Util = Uval(:);
        RespT = Rval(:);
        ResidT = Wval(:);
        Tput = Tval(:);
        ArvR = Aval(:);
    end
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    AvgTable = Table(Station, JobClass,  QLen, Util, RespT, ResidT, ArvR, Tput);
end
AvgTable = IndexedTable(AvgTable);
end

function strVals = formatWithCI(vals, ciVals)
% FORMATWITHCI Format values with confidence intervals as strings
% Returns a cell array of strings in format 'mean ± ci_halfwidth'
    n = length(vals);
    strVals = cell(n, 1);
    for i = 1:n
        if isnan(vals(i)) || vals(i) == 0
            strVals{i} = sprintf('%.4g', vals(i));
        elseif isnan(ciVals(i)) || ciVals(i) == 0
            strVals{i} = sprintf('%.4g', vals(i));
        else
            strVals{i} = sprintf('%.4g ± %.4g', vals(i), ciVals(i));
        end
    end
end

function val = getCIValue(avgStruct, fieldName, ist, k)
% GETCIVALUE Safely retrieve a CI value from the Avg structure
% Returns 0 if the field doesn't exist or is empty
    if isfield(avgStruct, fieldName) && ~isempty(avgStruct.(fieldName))
        val = avgStruct.(fieldName)(ist, k);
    else
        val = 0;
    end
end
