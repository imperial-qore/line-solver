function [QNn,UNn,RNn,TNn,ANn,WNn] = getAvgNode(self, Q, U, R, T, A, W)
% [QNN,UNN,RNN,TNN,ANn,WNn] = GETNODEAVG(Q, U, R, T, A, W)
%
% Compute average utilizations at steady-state for all nodes

if nargin == 1 % no parameter
    if   isempty(self.model.handles) || ~isfield(self.model.handles,'Q') || ...
        ~isfield(self.model.handles,'U') || ~isfield(self.model.handles,'R') || ...
        ~isfield(self.model.handles,'T') || ~isfield(self.model.handles,'A') || ...
        ~isfield(self.model.handles,'W')
        reset(self); % reset in case there are partial results saved
    end
    [Q,U,R,T,A,W] = self.getAvgHandles;
elseif nargin == 2
    handlers = Q;
    [Q,U,R,T,A,W] = deal(handlers{:}); % set Q=handlers{1}, U=handlers{2}, ...
end

[QN,UN,RN,TN,AN,WN] = self.getAvg(Q,U,R,T,A,W);
if isempty(QN)
    [QNn, UNn, RNn, TNn, ANn, WNn] = deal([]); % set each to []
    return
end

sn = self.model.getStruct; % must be called after getAvg eg for caches

I = sn.nnodes;  % Physical nodes only
M = sn.nstations;
R = sn.nclasses;
F = sn.nregions;  % Number of FCR virtual nodes

% Total nodes includes physical nodes + FCR virtual nodes
totalNodes = I + F;

% get average metrics that are zero at non-station nodes
QNn = zeros(totalNodes, R);
UNn = zeros(totalNodes, R);
RNn = zeros(totalNodes, R);
WNn = zeros(totalNodes, R);

% Map station metrics to node metrics
for ist=1:M
    ind = sn.stationToNode(ist);
    QNn(ind,:) = QN(ist,:);
    UNn(ind,:) = UN(ist,:);
    RNn(ind,:) = RN(ist,:);
    WNn(ind,:) = WN(ist,:);
end

% get the remaining node average metrics
ANn = sn_get_node_arvr_from_tput(sn, TN(1:M,:), T, AN);
TNn = sn_get_node_tput_from_tput(sn, TN(1:M,:), T, ANn);

% Fix arrival rates for ClassSwitch and Sink nodes for cache hit/miss classes
% The arrival rate at these nodes for hit/miss classes equals the Cache throughput
for cacheInd = 1:I
    if sn.nodetype(cacheInd) == NodeType.Cache
        hitclass = sn.nodeparam{cacheInd}.hitclass;
        missclass = sn.nodeparam{cacheInd}.missclass;
        % Update arrival rates at ClassSwitch and Sink nodes for hit/miss classes
        for ind = 1:I
            if sn.nodetype(ind) == NodeType.ClassSwitch || sn.nodetype(ind) == NodeType.Sink
                for classIdx = 1:R
                    % Check if this class is a hit class from the cache
                    if any(classIdx == hitclass(hitclass > 0))
                        ANn(ind, classIdx) = TNn(cacheInd, classIdx);
                    end
                    % Check if this class is a miss class from the cache
                    if any(classIdx == missclass(missclass > 0))
                        ANn(ind, classIdx) = TNn(cacheInd, classIdx);
                    end
                end
            end
        end
    end
end

% Extend ANn and TNn to include FCR rows
if F > 0
    ANn = [ANn; NaN(F, R)];  % FCR A is NaN (JMT doesn't provide)
    TNn = [TNn; zeros(F, R)];
end

% Map FCR pseudo-station metrics to FCR node indices
% FCR metrics are stored at station indices M+1 to M+F in self.result.Avg
% FCR node indices are I+1 to I+F in node-level matrices
% Note: getAvg() only returns station-level metrics (1:M), so we access FCR
% metrics directly from self.result.Avg
if F > 0 && isfield(self.result, 'Avg') && size(self.result.Avg.Q, 1) > M
    for f = 1:F
        fcrStationIdx = M + f;  % FCR pseudo-station index in result matrices
        fcrNodeIdx = I + f;     % FCR node index in nodenames
        QNn(fcrNodeIdx,:) = self.result.Avg.Q(fcrStationIdx,:);
        UNn(fcrNodeIdx,:) = self.result.Avg.U(fcrStationIdx,:);
        RNn(fcrNodeIdx,:) = self.result.Avg.R(fcrStationIdx,:);
        WNn(fcrNodeIdx,:) = self.result.Avg.W(fcrStationIdx,:);
        TNn(fcrNodeIdx,:) = self.result.Avg.T(fcrStationIdx,:);
        % ANn for FCR is NaN (already set above)
    end
end

end
