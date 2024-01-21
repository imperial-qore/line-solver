function [AvgTable,QT,UT,RT,WT,TT,AT] = getAvgTable(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,TT,AT] = GETAVGTABLE(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Return table of average station metrics
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
if GlobalConstants.DummyMode
    AvgTable = Table();
    QT = Table(); UT = Table();
    RT = Table(); TT = Table();
    WT = Table(); AT = Table();
    return
end
if ~isempty(self.model.obj)
    sn = self.model.getStruct;
    M = sn.nstations;
    R = sn.nclasses;
    V = cellsum(sn.visits);
    [QN,UN,RN,TN,AN,WN] = self.getAvg([],[],[],[],[],[]);
    Qval = []; Uval = [];
    Rval = []; Tval = [];
    Wval = []; Aval = [];
    Residval =[];
    JobClass = {};
    Station = {};
    for i=1:M
        for k=1:R
            if any(sum([QN(i,k),UN(i,k),RN(i,k),TN(i,k),AN(i,k),WN(i,k)])>0)
                c = find(sn.chains(:,k));
                inchain = sn.inchain{c};
                JobClass{end+1,1} = sn.classnames{k};
                Station{end+1,1} = sn.nodenames{sn.stationToNode(i)};
                Qval(end+1) = QN(i,k);
                Uval(end+1) = UN(i,k);
                Rval(end+1) = RN(i,k);
                if RN(i,k)<GlobalConstants.FineTol
                    Residval(end+1) = RN(i,k);
                else
                    if sn.refclass(c)>0
                        Residval(end+1) = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.refclass(c)));
                    else
                        Residval(end+1) = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.chains(c,:)));
                    end
                end
                Tval(end+1) = TN(i,k);
                Wval(end+1) = WN(i,k);
                Aval(end+1) = AN(i,k);
            end
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    Util = Uval(:); % we need to save first in a variable named like the column
    RespT = Rval(:); % we need to save first in a variable named like the column
    ResidT = Residval(:); % we need to save first in a variable named like the column
    Tput = Tval(:); % we need to save first in a variable named like the column
    ArvR = Aval(:); % we need to save first in a variable named like the column
    % QLen = num2str(QLen, ['%' sprintf('.%df', 5)]);
    % Util = num2str(Util, ['%' sprintf('.%df', 5)]);
    % RespT = num2str(RespT, ['%' sprintf('.%df', 5)]);
    % ResidT = num2str(ResidT, ['%' sprintf('.%df', 5)]);
    % ArvR = num2str(ArvR, ['%' sprintf('.%df', 5)]);
    % Tput = num2str(Tput, ['%' sprintf('.%df', 5)]);
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    %AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput);
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
        Q = param{1};
        U = param{2};
        R = param{3};
        T = param{4};
        A = param{5};
        W = param{6};
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
    RN = zeros(size(QN));
    WN = zeros(size(QN));
    AN = zeros(size(QN));
else
    [QN,UN,RN,TN,AN,WN] = self.getAvg(Q,U,R,T,A,W);
end

% this is required because getAvg can alter the chain structure in the
% presence of caches
sn = self.model.getStruct;

if isempty(QN)
    AvgTable = Table();
    QT = Table();
    UT = Table();
    RT = Table();
    TT = Table();
    AT = Table();
    WT = Table();
elseif ~keepDisabled
    V = cellsum(sn.visits);
    if isempty(V) % SSA
        for i=1:M
            for c=1:sn.nchains
                chain_classes = find(sn.chains(c,:));
                k = chain_classes(1);
                Tchain=sum(TN(sn.refstat(k),chain_classes));
                for k=chain_classes
                    V(i,k)=TN(i,k)/Tchain;
                end
            end
        end
    end

    Qval = []; Uval = [];
    Rval = []; Tval = [];
    Wval = []; Aval = [];
    Residval =[];
    JobClass = {};
    Station = {};
    for i=1:M
        for k=1:K
            if any(sum([QN(i,k),UN(i,k),RN(i,k),TN(i,k),AN(i,k),WN(i,k)])>0)
                c = find(sn.chains(:,k));
                inchain = sn.inchain{c};
                JobClass{end+1,1} = Q{i,k}.class.name;
                Station{end+1,1} = Q{i,k}.station.name;
                Qval(end+1) = QN(i,k);
                Uval(end+1) = UN(i,k);
                Rval(end+1) = RN(i,k);
                if RN(i,k)<GlobalConstants.FineTol
                    Residval(end+1) = RN(i,k);
                else
                    if sn.refclass(c)>0
                        Residval(end+1) = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.refclass(c)));
                    else
                        Residval(end+1) = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.chains(c,:)));
                    end
                end
                Tval(end+1) = TN(i,k);
                Wval(end+1) = WN(i,k);
                Aval(end+1) = AN(i,k);
            end
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);

    QLen = Qval(:); % we need to save first in a variable named like the column
    Util = Uval(:); % we need to save first in a variable named like the column
    RespT = Rval(:); % we need to save first in a variable named like the column
    ResidT = Residval(:); % we need to save first in a variable named like the column
    Tput = Tval(:); % we need to save first in a variable named like the column
    ArvR = Aval(:); % we need to save first in a variable named like the column
    % QLen = num2str(QLen, ['%' sprintf('.%df', 5)]);
    % Util = num2str(Util, ['%' sprintf('.%df', 5)]);
    % RespT = num2str(RespT, ['%' sprintf('.%df', 5)]);
    % ResidT = num2str(ResidT, ['%' sprintf('.%df', 5)]);
    % ArvR = num2str(ArvR, ['%' sprintf('.%df', 5)]);
    % Tput = num2str(Tput, ['%' sprintf('.%df', 5)]);
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    %AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput);
else
    V = cellsum(sn.visits);
    if isempty(V) % SSA
        for i=1:M
            for c=1:sn.nchains
                chain_classes = find(sn.chains(c,:));
                k = chain_classes(1);
                Tchain=sum(TN(sn.refstat(k),chain_classes));
                for k=chain_classes
                    V(i,k)=TN(i,k)/Tchain;
                end
            end
        end
    end
    Qval = zeros(M,K); Uval = zeros(M,K);
    Rval = zeros(M,K); Tval = zeros(M,K);
    Residval = zeros(M,K);
    JobClass = cell(K*M,1);
    Station = cell(K*M,1);
    for i=1:M
        for k=1:K
            c = find(sn.chains(:,k));
            inchain = sn.inchain{c};
            JobClass{(i-1)*K+k} = Q{i,k}.class.name;
            Station{(i-1)*K+k} = Q{i,k}.station.name;
            Qval((i-1)*K+k) = QN(i,k);
            Uval((i-1)*K+k) = UN(i,k);
            Rval((i-1)*K+k) = RN(i,k);
            if RN(i,k)<GlobalConstants.FineTol
                Residval((i-1)*K+k) = RN(i,k);
            else
                if sn.refclass(c)>0
                    Residval((i-1)*K+k)  = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.refclass(c)));
                else
                    Residval((i-1)*K+k)  = RN(i,k)*V(i,k)/sum(V(sn.refstat(k),sn.chains(c,:)));
                end
            end
            Tval((i-1)*K+k) = TN(i,k);
            Aval((i-1)*K+k) = AN(i,k);
        end
    end
    Station = label(Station);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    Util = Uval(:); % we need to save first in a variable named like the column
    RespT = Rval(:); % we need to save first in a variable named like the column
    ResidT = Residval(:); % we need to save first in a variable named like the column
    Tput = Tval(:); % we need to save first in a variable named like the column
    ArvR = Aval(:); % we need to save first in a variable named like the column
    % QLen = num2str(QLen, ['%' sprintf('.%df', 5)]);
    % Util = num2str(Util, ['%' sprintf('.%df', 5)]);
    % RespT = num2str(RespT, ['%' sprintf('.%df', 5)]);
    % ResidT = num2str(ResidT, ['%' sprintf('.%df', 5)]);
    % ArvR = num2str(ArvR, ['%' sprintf('.%df', 5)]);
    % Tput = num2str(Tput, ['%' sprintf('.%df', 5)]);
    QT = Table(Station,JobClass,QLen);
    UT = Table(Station,JobClass,Util);
    RT = Table(Station,JobClass,RespT);
    WT = Table(Station,JobClass,ResidT);
    TT = Table(Station,JobClass,Tput);
    AT = Table(Station,JobClass,ArvR);
    AvgTable = Table(Station, JobClass,  QLen, Util, RespT, ResidT, ArvR, Tput);
end
end
