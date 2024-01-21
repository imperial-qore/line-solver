%% getTranHandles: add all transient mean performance indexes
% Q{i,r}: timeseries of mean queue-length of class r at node i
% U{i,r}: timeseries of mean utilization of class r at node i
% R{i,r}: timeseries of mean response time of class r at node i (summed across visits)
% T{i,r}: timeseries of mean throughput of class r at node i
function [Qt,Ut,Tt] = getTranHandles(self)
% [QT,UT,TT] = GETTRANHANDLES()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
M = getNumberOfStations(self);
K = getNumberOfClasses(self);

Tt = cell(1,K); % throughputs
Qt = cell(M,K); % queue-length
%Rt = cell(M,K); % response times
Ut = cell(M,1); % utilizations
for i=1:M
    for r=1:K
        Tt{i,r} = Metric(MetricType.TranTput, self.classes{r}, self.stations{i});
        Qt{i,r} = Metric(MetricType.TranQLen, self.classes{r}, self.stations{i});
        %        Rt{i,r} = Metric(MetricType.TranRespT, self.classes{r}, self.stations{i});
        Ut{i,r} = Metric(MetricType.TranUtil, self.classes{r}, self.stations{i});
        if isa(self.stations{i},'Source')
            Qt{i,r}.disabled = true;
            %            Rt{i,r}.disabled = true;
            Ut{i,r}.disabled = true;
        end
        if isa(self.stations{i},'Sink')
            Qt{i,r}.disabled = true;
            %            Rt{i,r}.disabled = true;
            Ut{i,r}.disabled = true;
        end
        if isa(self.stations{i},'Join') || isa(self.stations{i},'Fork')
            Ut{i,r}.disabled = true;
        end
        if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
            if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                Tt{i,r}.disabled = true;
                Qt{i,r}.disabled = true;
                %                Rt{i,r}.disabled = true;
                Ut{i,r}.disabled = true;
            end
        end
    end
end
end
