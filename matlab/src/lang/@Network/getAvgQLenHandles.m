function [Q] = getAvgQLenHandles(self)
% [Q] = GETAVGQLENHANDLES()

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
if isempty(self.handles) || ~isfield(self.handles,'Q')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    Q = cell(M,K); % queue-length
    for i=1:M
        for r=1:K
            Q{i,r} = Metric(MetricType.QLen, self.classes{r}, self.stations{i});                        
            if isa(self.stations{i},'Source')
                Q{i,r}.disable=true;
            end            
            if isa(self.stations{i},'Sink')
                Q{i,r}.disable=true;
            end            
            if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
                if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                    Q{i,r}.disable=true;
                end
            end
        end
    end
    self.handles.Q = Q;
else
    Q = self.handles.Q;
end
end
