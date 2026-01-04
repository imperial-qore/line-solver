function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgNodeTable(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,AT,TT] = GETNODEAVGTABLE(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Return table of average node metrics
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<8 %~exist('keepDisabled','var')
    keepDisabled = false;
end

if GlobalConstants.DummyMode
    [AvgTable, QT, UT, RT, TT, WT, AT] = deal(Table());
    AvgTable = IndexedTable(AvgTable);
    return
end

sn = self.model.getStruct;
I = sn.nnodes + sn.nregions;  % Include FCR virtual nodes
K = sn.nclasses;

if nargin == 1
    [Q,U,R,T,A,W] = getAvgHandles(self);
elseif isempty(Q) && isempty(U) && isempty(R) && isempty(T) && isempty(A) && isempty(W)
    [Q,U,R,T,A,W] = getAvgHandles(self);
end

[QN,UN,RN,TN,AN,WN] = self.getAvgNode(Q,U,R,T,A,W);

if isempty(QN)
    [AvgTable, QT, UT, RT, WT, TT, AT] = deal(Table());
    AvgTable = IndexedTable(AvgTable);
elseif ~keepDisabled
    [Qval, Uval, Rval, Tval, Aval, Wval] = deal([]);
    JobClass = {};
    Node = {};
    for ind=1:I
        for k=1:K
            % Check if any metric is non-zero (use nansum to handle NaN values for FCR nodes)
            filterVal = nansum([QN(ind,k),UN(ind,k),RN(ind,k),TN(ind,k),AN(ind,k)]);
            if filterVal>0
                c = sn.chains(:,k)>0;
                inchain = sn.inchain{c};
                JobClass{end+1,1} = sn.classnames{k};
                Node{end+1,1} = sn.nodenames{ind};
                Qval(end+1) = QN(ind,k);
                Uval(end+1) = UN(ind,k);
                Rval(end+1) = RN(ind,k);
                Wval(end+1) = WN(ind,k);
                Tval(end+1) = TN(ind,k);
                Aval(end+1) = AN(ind,k);
            end
        end
    end
    Node = label(Node);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    ResidT = Wval(:); % we need to save first in a variable named like the column
    WT = Table(Node,JobClass,ResidT);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ResidT,ArvR,Tput);
    AvgTable = IndexedTable(AvgTable);
else
    [Qval, Uval, Rval, Tval, Aval, Wval] = deal(zeros(I, K));
    JobClass = {};
    Node = {};
    for ind=1:I
        for k=1:K
            c = sn.chains(:,k)>0;
            inchain = sn.inchain{c};
            JobClass{end+1,1} = sn.classnames{k};
            Node{end+1,1} = sn.nodenames{ind};
            Qval((ind-1)*K+k) = QN(ind,k);
            Uval((ind-1)*K+k) = UN(ind,k);
            Rval((ind-1)*K+k) = RN(ind,k);
            Wval((ind-1)*K+k) = WN(ind,k);
            Tval((ind-1)*K+k) = TN(ind,k);
            Aval((ind-1)*K+k) = AN(ind,k);
        end
    end
    Node = label(Node);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    ResidT = Wval(:); % we need to save first in a variable named like the column
    WT = Table(Node,JobClass,ResidT);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ArvR,Tput);
    AvgTable = IndexedTable(AvgTable);
end
end
