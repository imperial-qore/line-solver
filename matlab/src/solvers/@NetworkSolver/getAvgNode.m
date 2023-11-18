function [QNn,UNn,RNn,TNn,ANn,WNn] = getAvgNode(self, Q, U, R, T, A, W)
% [QNN,UNN,RNN,TNN,ANn,WNn] = GETNODEAVG(Q, U, R, T, A, W)

% Compute average utilizations at steady-state for all nodes
if nargin == 1 % no parameter
    if isempty(self.model.handles) || ~isfield(self.model.handles,'Q') || ~isfield(self.model.handles,'U') || ~isfield(self.model.handles,'R') || ~isfield(self.model.handles,'T') || ~isfield(self.model.handles,'A')
        resetResults(self); % reset in case there are partial results saved
    end
    [Q,U,R,T,A,W] = self.getAvgHandles;
elseif nargin == 2
    handlers = Q;
    Q=handlers{1};
    U=handlers{2};
    R=handlers{3};
    T=handlers{4};
    A=handlers{5};
    W=handlers{6};
end
[QN,UN,RN,TN,AN,WN] = self.getAvg(Q,U,R,T,A,W);
if isempty(QN)
    QNn=[];
    UNn=[];
    RNn=[];
    TNn=[];
    ANn=[];
    WNn=[];
    return
end

sn = self.model.getStruct; % must be called after getAvg

I = sn.nnodes;
M = sn.nstations;
R = sn.nclasses;
C = sn.nchains;
QNn = zeros(I,R);
UNn = zeros(I,R);
RNn = zeros(I,R);
TNn = zeros(I,R);
ANn = zeros(I,R);
WNn = zeros(I,R);

for ist=1:M
    ind = sn.stationToNode(ist);
    QNn(ind,:) = QN(ist,:);
    UNn(ind,:) = UN(ist,:);
    RNn(ind,:) = RN(ist,:);
    WNn(ind,:) = WN(ist,:);
    %TNn(ind,:) = TN(ist,:); % this is built later from ANn
    ANn(ind,:) = AN(ist,:);
end

%if any(sn.isstatedep(:,3)) || any(sn.nodetype == NodeType.Cache)
%    line_warning(mfilename,'Node-level metrics not available in models with state-dependent routing. Returning station-level metrics only.');
%    return
%end

% update tputs for all nodes but the sink and the joins
for ind=1:I
    for c = 1:C
        inchain = sn.inchain{c};
        for r = inchain
            %anystat = find(sn.visits{c}(:,r));
            refstat = sn.refstat(c);
            %if ~isempty(anystat)
            if sn.nodetype(ind) ~= NodeType.Source
                %if sn.isstation(ind)
                %    ist = sn.nodeToStation(ind);
                %    ANn(ind, r) =  TN(ist,r);
                %else
                switch sn.nodetype(ind)
                    case NodeType.Cache
                        if any(find(r==sn.nodeparam{ind}.hitclass))
                        TNn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                        elseif any(find(r==sn.nodeparam{ind}.missclass))
                        TNn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                        else
                        ANn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                        end
                    otherwise
                        ANn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                end
                %end
            end
            %end
        end
    end
end

for ind=1:I
    for c = 1:C
        inchain = sn.inchain{c};
        for r = inchain
            anystateful = find(sn.visits{c}(:,r));
            if ~isempty(anystateful)
                if sn.nodetype(ind) ~= NodeType.Sink && sn.nodetype(ind) ~= NodeType.Join
                    for s = inchain
                        for jnd=1:I
                            switch sn.nodetype(ind)
                                case NodeType.Source
                                    ist = sn.nodeToStation(ind);
                                    TNn(ind, s) = TN(ist,s);
                                case NodeType.Cache
                                     if ind~=jnd
                                         TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                                     end
                                otherwise
                                    TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                            end
                        end
                    end
                elseif sn.nodetype(ind) == NodeType.Join
                    for s = inchain
                        for jnd=1:I
                            if sn.nodetype(ind) ~= NodeType.Source
                                TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                            else
                                ist = sn.nodeToStation(ind);
                                TNn(ind, s) = TN(ist,s);
                            end
                        end
                    end
                end
            end
        end
    end
end

ANn(isnan(ANn)) = 0;
TNn(isnan(TNn)) = 0;
end
