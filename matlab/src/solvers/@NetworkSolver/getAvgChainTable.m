function [AvgChain,QTc,UTc,RTc,WTc,TTc] = getAvgChainTable(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,WTc,TTC] = GETAVGCHAINTABLE(SELF,Q,U,R,T)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    AvgChain = Table();
    QTc = Table(); UTc = Table();
    RTc = Table(); TTc = Table();
    WTc = Table();
    return
end

sn = self.model.getStruct;
M = sn.nstations;
C = sn.nchains;
if nargin == 1
    [Q,U,R,T] = getAvgHandles(self);
end
if nargin == 2
    if iscell(Q) && ~isempty(Q)
        param = Q;
        Q = param{1};
        U = param{2};
        R = param{3};
        T = param{4};
        % case where varargin is passed as input
    elseif iscell(Q) && isempty(Q)
        [Q,U,R,T] = getAvgHandles(self);
    end
end
[QNc,UNc,RNc,TNc] = self.getAvgChain(Q,U,R,T);

ChainObj = self.model.getChains();
ChainName = cellfun(@(c) c.name,ChainObj,'UniformOutput',false);
ChainClasses = cell(1,length(ChainName));
for c=1:length(ChainName)
    ChainClasses{c} = ChainObj{c}.classnames;
end
if isempty(QNc)
    AvgChain = Table();
    QTc = Table();
    UTc = Table();
    RTc = Table();
    TTc = Table();
else
    Qval = zeros(M,C); Uval = zeros(M,C);
    Rval = zeros(M,C); Tval = zeros(M,C);
    Resval = zeros(M,C);
    Chain = cell(C*M,1);
    JobClasses = cell(C*M,1);
    Station = cell(C*M,1);
    for c=1:sn.nchains
        for i=1:M
            Chain{(i-1)*C+c} = ChainName{c};
            JobClasses((i-1)*C+c,1) = {label(ChainClasses{c}(:))};
            Station{(i-1)*C+c} = Q{i,c}.station.name;
            Qval((i-1)*C+c) = QNc(i,c);
            Uval((i-1)*C+c) = UNc(i,c);
            Rval((i-1)*C+c) = RNc(i,c);
            Resval((i-1)*C+c) = RNc(i,c)/sum(sn.visits{c}(sn.stationToStateful(i),:));
            Tval((i-1)*C+c) = TNc(i,c);
        end
    end
    Chain = label(Chain);
    Station = label(Station);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QTc = Table(Station,Chain,JobClasses,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UTc = Table(Station,Chain,JobClasses,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RTc = Table(Station,Chain,JobClasses,RespT);
    ResidT = Resval(:); % we need to save first in a variable named like the column
    WTc = Table(Station,Chain,JobClasses,ResidT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TTc = Table(Station,Chain,JobClasses,Tput);
    AvgChain = Table(Station,Chain,JobClasses,QLen,Util,RespT,ResidT,Tput);
end
end
