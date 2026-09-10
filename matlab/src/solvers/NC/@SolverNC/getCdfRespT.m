function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

% lang='cpp' takes the sojourn law from line-cli (-s nc -a cdf), which
% evaluates it on the same logarithmic time grid this function uses. It is the
% same quantity, so recomputing it here would report a MATLAB evaluation as a
% C++ one.
%
% A model with NO FCFS station carries no tagged-job law at all, and the branch
% below then returns {} rather than a grid of empties. line-cli warns and sends
% no curve, which the shared reader turns into an (nstations x nclasses) cell of
% empties, so the shape is restated here: a caller testing isempty(RD) would
% otherwise read "no law" as "a law everywhere, all of it missing". The two are
% the same case here because the branch below fills EVERY station row whenever
% it fills any.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    RD = CPPLINE.cdfRespT(self.name, self.model, self.options);
    if all(cellfun(@isempty, RD(:)))
        RD = {};
    end
    return
end

T0 = tic;
if nargin<2 || isempty(R) %~exist('R','var')
    R = self.getAvgRespTHandles;
end

config = self.getOptions.config;
if ~isfield(config,'algorithm')
	config.algorithm = 'exact';
end
RD = {};
sn = self.getStruct;
[~,D,N,Z,~,S]= sn_get_product_form_params(sn);
fcfsNodes = find(sn.sched(sn.sched ~= SchedStrategy.INF) == SchedStrategy.FCFS);
fcfsNodeIds = find(sn.sched == SchedStrategy.FCFS);
delayNodeIds = find(sn.sched == SchedStrategy.INF);
if ~isempty(fcfsNodes)
    % The tagged-job passage time is defined on a CLOSED network: pfqn_stdf
    % conditions on a job circulating among a fixed population. Without this
    % gate an open class carried its infinite population into the quadrature
    % and MATLAB failed with the unspecific 'MATLAB:nonaninf' ("NaN and Inf not
    % allowed"), naming neither the model nor the requirement.
    if any(isinf(N))
        line_error(mfilename, 'The tagged-job sojourn-time law is defined on a CLOSED network; this model has an open class. Use SolverFluid or SolverMAM for a response-time distribution with open classes.');
    end
    % THREE index spaces meet here and mixing them is a crash, not a wrong
    % number. sn.rates is indexed by STATION, so the horizon uses fcfsNodeIds.
    % fcfsNodes is relative to the NON-DELAY stations and is what pfqn_stdf
    % indexes L(k,:), S(k) and rates(k,:) with, D having one row per Queue
    % node. The rates passed below must therefore be in the NON-DELAY space
    % too: slicing sn.sched==FCFS gives one row per FCFS station instead, and
    % the two coincide only when every non-delay station is FCFS. On
    % Delay+PS+FCFS that made fcfsNodes=2 index a 1-row matrix and pfqn_stdf
    % died with "Index in position 1 exceeds array bounds".
    T = max(sum(N) * mean(1./sn.rates(fcfsNodeIds,:)));
    tset = logspace(0,2*log10(T),100);
    rates = sn.rates(sn.sched ~= SchedStrategy.INF,:);
    switch config.algorithm
        case 'exact'
            RDout = pfqn_stdf(D,N,Z,S,fcfsNodes,rates,tset);
        case 'rd'
            RDout = pfqn_stdf_heur(D,N,Z,S,fcfsNodes,rates,tset);
    end
    % RDout rows are in the NON-DELAY space (pfqn_stdf fills row fcfsNodes(x)),
    % while RD rows are in the STATION space. Walk the two lists together:
    % iterating over size(RDout,1) and indexing fcfsNodeIds with that counter
    % mixes the spaces again and overruns as soon as some non-delay station is
    % not FCFS, which is the same Delay+PS+FCFS model as above.
    for idx=1:length(fcfsNodes)
        for j=1:size(RDout,2)
            RD{fcfsNodeIds(idx),j} = real(RDout{fcfsNodes(idx),j}); % drop complex round-offs
        end
    end
    for i=1:length(delayNodeIds)
        for j=1:size(RDout,2)
            RD{delayNodeIds(i),j} = [map_cdf(sn.proc{delayNodeIds(i)}{j}, tset(:))' tset(:)];
        end
    end
    runtime = toc(T0);
    self.setDistribResults(RD, runtime);
else
    line_warning(mfilename, 'getCdfRespT applies only to FCFS nodes.\n');
    return
end
end