function [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = getAvgNodeChainTable(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,WTC,ATC,TTC] = GETAVGCHAINNODETABLE(SELF,Q,U,R,T)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    [AvgChain, QTc, UTc, RTc, TTc, WTc] = deal(Table());
    AvgChain = IndexedTable(AvgChain);
    return;
end

sn = self.model.getStruct;
I = sn.nnodes;
C = sn.nchains;
if nargin == 1
    [Q,U,R,T] = getAvgHandles(self);
end
if nargin == 2 && iscell(Q)
    if ~isempty(Q)
        [Q, U, R, T] = deal(Q{:});
    else
        [Q, U, R, T] = getAvgHandles(self);
    end
end

[QNc,UNc,RNc,WNc,ANc,TNc] = self.getAvgNodeChain(Q,U,R,T);

ChainObj = self.model.getChains();
ChainName = cellfun(@(c) c.name,ChainObj,'UniformOutput',false);
ChainClasses = cell(1,length(ChainName));
for c=1:length(ChainName)
    ChainClasses{c} = ChainObj{c}.classnames;
end
if isempty(QNc)
    [AvgChain, QTc, UTc, RTc, TTc] = deal(Table());
    AvgChain = IndexedTable(AvgChain);
else
    [Qval, Uval, Rval, Tval, Resval, Aval] = deal(zeros(I, C));
    [Chain, JobClasses, Node] = deal(cell(C*I, 1));
    for c = 1:sn.nchains
        for ind = 1:I
            idx = (ind-1)*C + c;
            Chain{idx} = ChainName{c};
            JobClasses{idx} = label(ChainClasses{c}(:));
            Node{idx} = sn.nodenames{ind};
            Qval(idx) = QNc(ind, c);
            Uval(idx) = UNc(ind, c);
            Rval(idx) = RNc(ind, c);
            Aval(idx) = ANc(ind, c);
            Tval(idx) = TNc(ind, c);
            Resval(idx) = WNc(ind,c);
        end
    end
    [Chain, Node] = deal(label(Chain), label(Node));
    QTc = table(Node, Chain, JobClasses, Qval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'QLen'});
    UTc = table(Node, Chain, JobClasses, Uval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'Util'});
    RTc = table(Node, Chain, JobClasses, Rval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'RespT'});
    WTc = table(Node, Chain, JobClasses, Resval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'ResidT'});
    TTc = table(Node, Chain, JobClasses, Tval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'Tput'});
    ATc = table(Node, Chain, JobClasses, Aval(:), 'VariableNames', {'Node', 'Chain', 'JobClasses', 'ArvR'});

    AvgChain = table(Node, Chain, JobClasses, Qval(:), Uval(:), Rval(:), Resval(:), Aval(:), Tval(:), ...
        'VariableNames', {'Node', 'Chain', 'JobClasses', 'QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput'});
    AvgChain = IndexedTable(AvgChain);
end
end
