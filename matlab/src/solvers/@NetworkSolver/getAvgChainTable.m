function [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = getAvgChainTable(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,WTc,TTC] = GETAVGCHAINTABLE(SELF,Q,U,R,T)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    [AvgChain, QTc, UTc, RTc, TTc, WTc, ATc] = deal(Table());
    AvgChain = IndexedTable(AvgChain);
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
        [Q, U, R, T] = deal(param{1:4});
        % case where varargin is passed as input
    elseif iscell(Q) && isempty(Q)
        [Q,U,R,T] = getAvgHandles(self);
    end
end
[QNc,UNc,RNc,WNc,ANc,TNc] = self.getAvgChain(Q,U,R,T);

ChainObj = self.model.getChains();
ChainName = cellfun(@(c) c.name,ChainObj,'UniformOutput',false);
ChainClasses = cell(1,length(ChainName));
for c=1:length(ChainName)
    ChainClasses{c} = ChainObj{c}.classnames;
end
if isempty(QNc)
    [AvgChain, QTc, UTc, RTc, TTc, ATc] = deal(Table());
    AvgChain = IndexedTable(AvgChain);
else
    [Qval, Uval, Rval, Tval, Resval, Aval] = deal(zeros(M, C));
    Chain = cell(C*M,1);
    JobClasses = cell(C*M,1);
    Station = cell(C*M,1);
    for c=1:sn.nchains
        for ist=1:M
            Chain{(ist-1)*C+c} = ChainName{c};
            JobClasses((ist-1)*C+c,1) = {label(ChainClasses{c}(:))};
            Station{(ist-1)*C+c} = Q{ist,c}.station.name;
            Qval((ist-1)*C+c) = QNc(ist,c);
            Uval((ist-1)*C+c) = UNc(ist,c);
            Rval((ist-1)*C+c) = RNc(ist,c);
            Resval((ist-1)*C+c) = WNc(ist,c);
            Aval((ist-1)*C+c) = ANc(ist,c);
            Tval((ist-1)*C+c) = TNc(ist,c);
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
    ArvR = Aval(:); % we need to save first in a variable named like the column
    ATc = Table(Station,Chain,JobClasses,ArvR);
    AvgChain = Table(Station,Chain,JobClasses,QLen,Util,RespT,ResidT,ArvR,Tput);
    AvgChain = IndexedTable(AvgChain);
end
end
