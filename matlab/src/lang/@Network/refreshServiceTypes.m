function [procid] = refreshServiceTypes(self, statSet, classSet)
% [PROCTYPE] = REFRESHSERVICETYPES()

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
elseif nargin==2
    line_error(mfilename,'refreshRates requires either 0 or 2 parameters.')
elseif nargin==3 && isfield(self.sn,'rates') && isfield(self.sn,'scv')
    % we are only updating selected stations and classes so use the
    % existing ones for the others
    procid = self.sn.procid;
else
    procid = nan(M,K); % disabled
end
hasOpenClasses = self.hasOpenClasses;

stations = self.stations;
% determine rates
for i=statSet
    
    for r=classSet
        switch stations{i}.server.className
            case 'ServiceTunnel'
                % do nothing
                switch class(stations{i})
                    case 'Source'
                        if isempty(stations{i}.input.sourceClasses{r}) || stations{i}.input.sourceClasses{r}{end}.isDisabled
                            procid(i,r) = ProcessType.ID_DISABLED;
                        elseif stations{i}.input.sourceClasses{r}{end}.isImmediate
                            procid(i,r) = ProcessType.ID_IMMEDIATE;
                        else
                            procid(i,r) = ProcessType.toId(ProcessType.fromText(class(stations{i}.input.sourceClasses{r}{end})));
                        end
                    case 'Join'
                        procid(i,r) = ProcessType.ID_IMMEDIATE;
                end
            otherwise
                if ~hasOpenClasses || i ~= self.getIndexSourceStation
                    if isempty(stations{i}.server.serviceProcess{r}) || stations{i}.server.serviceProcess{r}{end}.isDisabled
                        procid(i,r) = ProcessType.ID_DISABLED;
                    elseif stations{i}.server.serviceProcess{r}{end}.isImmediate
                        procid(i,r) = ProcessType.ID_IMMEDIATE;
                    else
                        processtype = class(stations{i}.server.serviceProcess{r}{end});
                        processtypeid = ProcessType.fromText(processtype);
                        procid(i,r) = ProcessType.toId(processtypeid);
                    end
                end
        end
    end
end

if ~isempty(self.sn)
    self.sn.procid = procid;
end
end