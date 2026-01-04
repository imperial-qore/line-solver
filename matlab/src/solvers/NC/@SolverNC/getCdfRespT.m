function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)
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
    T = max(sum(N) * mean(1./sn.rates(fcfsNodes,:)));
    tset = logspace(0,2*log10(T),100);
    rates = sn.rates(sn.sched == SchedStrategy.FCFS,:);
    switch config.algorithm
        case 'exact'
            RDout = pfqn_stdf(D,N,Z,S,fcfsNodes,rates,tset);
        case 'rd'
            RDout = pfqn_stdf_heur(D,N,Z,S,fcfsNodes,rates,tset);
    end
    for i=1:size(RDout,1)
        for j=1:size(RDout,2)
            RD{fcfsNodeIds(i),j} = real(RDout{i,j}); % remove complex number round-offs
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