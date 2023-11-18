% W(i,r): mean residence time of class r at node i (summed across visits)
function W = getAvgResidTHandles(self)
% W = GETAVGRESIDTHANDLES()

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
if isempty(self.handles) || ~isfield(self.handles,'W')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    W = cell(M,K); % response times
    for i=1:M
        for r=1:K
            W{i,r} = Metric(MetricType.ResidT, self.classes{r}, self.stations{i});
            if isa(self.stations{i},'Source')
                W{i,r}.disable=true;
            end
            if isa(self.stations{i},'Sink')
                W{i,r}.disable=true;
            end
            if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
                if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                    W{i,r}.disable=true;
                end
            end
        end
    end
    self.handles.W = W;
else
    W = self.handles.W;
end
end
