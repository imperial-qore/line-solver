% T(i,r): mean throughput of class r at station i
function [T] = getAvgTputHandles(self)
% [T] = GETAVGTPUTHANDLES()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
if isempty(self.handles) || ~isfield(self.handles,'T')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    T = cell(1,K); % throughputs
    for i=1:M
        for r=1:K
            T{i,r} = Metric(MetricType.Tput, self.classes{r}, self.stations{i});
            if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
                if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                    T{i,r}.disable=true;
                end
            end
        end
    end
    self.handles.T = T;
else
    T = self.handles.T;
end
end
