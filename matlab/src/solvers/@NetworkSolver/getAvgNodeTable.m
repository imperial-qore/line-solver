function [AvgTable,QT,UT,RT,WT,TT,AT] = getAvgNodeTable(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,TT] = GETNODEAVGTABLE(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Return table of average node metrics
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
if nargin<8 %~exist('keepDisabled','var')
    keepDisabled = false;
end
sn = self.model.getStruct;
I = sn.nnodes;
K = sn.nclasses;
if nargin == 1
    [Q,U,R,T,A,W] = getAvgHandles(self);
elseif isempty(Q) && isempty(U) && isempty(R) && isempty(T) && isempty(A) && isempty(W)
    [Q,U,R,T,A,W] = getAvgHandles(self);
end

[QN,UN,RN,TN,AN,WN] = self.getAvgNode(Q,U,R,T,A,W);

if isempty(QN)
    AvgTable = Table();
    QT = Table();
    UT = Table();
    RT = Table();
    WT = Table();
    TT = Table();
    AT = Table();
elseif ~keepDisabled
    V = cellsum(sn.nodevisits);
    if isempty(V) % SSA
        for i=1:I
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
    Rval = []; Tval = []; Aval=[];
    Residval = [];
    JobClass = {};
    Node = {};
    for i=1:I
        for k=1:K
            if any(sum([QN(i,k),UN(i,k),RN(i,k),TN(i,k),AN(i,k)])>0)
                c = find(sn.chains(:,k));
                inchain = sn.inchain{c};
                JobClass{end+1,1} = sn.classnames{k};
                Node{end+1,1} = sn.nodenames{i};
                Qval(end+1) = QN(i,k);
                Uval(end+1) = UN(i,k);
                Rval(end+1) = RN(i,k);
                if RN(i,k)<GlobalConstants.FineTol
                    Residval(end+1) = WN(i,k);
                else
                    if sn.refclass(c)>0
                        Residval(end+1) = WN(i,k)/sum(V(sn.refstat(k),sn.refclass(c)));
                    else
                        Residval(end+1) = WN(i,k)/sum(V(sn.refstat(k),sn.chains(c,:)));
                    end
                end
                Tval(end+1) = TN(i,k);
                Aval(end+1) = AN(i,k);
            end
        end
    end
    Node = label(Node);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    ResidT = Residval(:); % we need to save first in a variable named like the column
    WT = Table(Node,JobClass,ResidT);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ResidT,ArvR,Tput);
else
    V = cellsum(sn.visits);
    if isempty(V) % SSA
        for i=1:I
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

    Qval = zeros(I,K); Uval = zeros(I,K);
    Rval = zeros(I,K); Tval = zeros(I,K); Aval = zeros(I,K);
    Residval = zeros(I,K);
    JobClass = {};
    Node = {};
    for i=1:I
        for k=1:K
            c = find(sn.chains(:,k));
            inchain = sn.inchain{c};
            JobClass{end+1,1} = sn.classnames{k};
            Node{end+1,1} = sn.nodenames{i};
            Qval((i-1)*K+k) = QN(i,k);
            Uval((i-1)*K+k) = UN(i,k);
            Rval((i-1)*K+k) = RN(i,k);
            if RN(i,k)<GlobalConstants.FineTol
                Residval((i-1)*K+k) = RN(i,k);
            else
                if sn.refclass(c)>0
                    Residval((i-1)*K+k)  = WN(i,k)/sum(V(sn.refstat(k),sn.refclass(c)));
                else
                    Residval((i-1)*K+k)  = WN(i,k)/sum(V(sn.refstat(k),sn.chains(c,:)));
                end
            end

            Tval((i-1)*K+k) = TN(i,k);
            Aval((i-1)*K+k) = AN(i,k);
        end
    end
    Node = label(Node);
    JobClass = label(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    ResidT = Residval(:); % we need to save first in a variable named like the column
    WT = Table(Node,JobClass,ResidT);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ArvR,Tput);
end
end
