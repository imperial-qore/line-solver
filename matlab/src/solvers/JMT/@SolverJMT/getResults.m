function [result, parsed] = getResults(self)
% [RESULT, PARSED] = GETRESULTS()

options = self.getOptions;
switch options.method
    case {'jsim','default'}
        [result, parsed] = self.getResultsJSIM;
    otherwise
        [result, parsed] = self.getResultsJMVA;
end

sn = self.model.getStruct;

result.Avg.Q = zeros(sn.nstations, sn.nclasses);
result.Avg.U = zeros(sn.nstations, sn.nclasses);
result.Avg.R = zeros(sn.nstations, sn.nclasses);
result.Avg.T = zeros(sn.nstations, sn.nclasses);
result.Avg.A = zeros(sn.nstations, sn.nclasses);
result.Avg.W = zeros(sn.nstations, sn.nclasses);

for m=1:length(result.metric)
    metric = result.metric{m};
    switch metric.measureType
        case MetricType.QLen
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.Q(i,r) = metric.meanValue;
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.Q(i,r) = metric.meanValue;
                else
                    result.Avg.Q(i,r) = 0;
                end
            end
        case MetricType.Util
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.U(i,r) = metric.meanValue;
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.U(i,r) = metric.meanValue;
                else
                    result.Avg.U(i,r) = 0;
                end
            end
        case MetricType.RespT
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.R(i,r) = metric.meanValue;
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.R(i,r) = metric.meanValue;
                else
                    result.Avg.R(i,r) = 0;
                end
            end
        case MetricType.ResidT
            % JMT ResidT is inconsistently defined with LINE's on some
            % difficult class switching cases, hence we recompute it at the
            % level of the NetworkSolver class to preserve consistency

%             i = sn.nodeToStation(find(sn.nodenames == metric.station));
%             r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
%             if isinf(sn.njobs(r))
%                 result.Avg.W(i,r) = metric.meanValue;
%             else % 'closed'
%                 N = sn.njobs;
%                 chainIdx = find(sn.classnames == metric.class);
%                 if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
%                     result.Avg.W(i,r) = metric.meanValue;
%                 else
%                     result.Avg.W(i,r) = 0;
%                 end
%             end
        case MetricType.ArvR
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.A(i,r) = metric.meanValue;
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.A(i,r) = metric.meanValue;
                else
                    result.Avg.A(i,r) = 0;
                end
            end            
        case MetricType.Tput
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.T(i,r) = metric.meanValue;
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.T(i,r) = metric.meanValue;
                else
                    result.Avg.T(i,r) = 0;
                end
            end
            
    end
self.result = result;    
end
end