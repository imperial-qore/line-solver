function [ph, phases] = refreshProcessRepresentations(self)
% [PH, PHASES] = REFRESHPROCESSREPRESENTATIONS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


M = getNumberOfStations(self);
K = getNumberOfClasses(self);
ph = cell(M,1);
for ist=1:M
    ph{ist,1} = cell(1,K);
end
phases = zeros(M,K);
stations = self.stations;
for ist=1:M
    if ist == self.getIndexSourceStation
        ph_i = stations{ist}.getSourceRates();
    else
        switch class(stations{ist})
            case 'Fork'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
                ph_i = Coxian(mu_i,phi_i).getProcess;
            case 'Join'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                    ph_i{r} = Coxian(mu_i{r},phi_i{r}).getProcess;
                end
            otherwise
                ph_i = stations{ist}.getServiceRates();
        end
    end
    for r=1:K
        ph{ist}{r} = ph_i{r};
        if isempty(ph{ist}{r}) % fluid fails otherwise
            phases(ist,r) = 1;
        elseif ~isMAP(ph{ist}{r})
            % Non-Markovian distribution: convert to MAP representation
            ph{ist}{r} = convertToMAP(stations{ist}, r, ph{ist}{r});
            phases(ist,r) = length(ph{ist}{r}{1});
        elseif any(isnan(ph{ist}{r}{1}(:))) || any(isnan(ph{ist}{r}{2}(:))) % disabled
            phases(ist,r) = 0;
        else
            phases(ist,r) = length(ph{ist}{r}{1});
        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'mu')
    proc = ph;
    pie = cell(size(ph));
    for ist=1:M
        for r=1:K
            map_ir = ph{ist}{r};
            if ~isempty(map_ir)
                %.proc{i}{r} = map_normalize(map_ir);
                proc{ist}{r} = map_ir;
                pie{ist}{r} = map_pie(map_ir);
            else
                pie{ist}{r} = NaN;
            end
        end
    end
    self.sn.proc = proc;
    self.sn.pie = pie;
    self.sn.phases = phases;
    self.sn.phasessz = max(self.sn.phases,ones(size(self.sn.phases)));
    self.sn.phasessz(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:)=phases(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:);
    self.sn.phaseshift = [zeros(size(phases,1),1),cumsum(self.sn.phasessz,2)];
end
end

function result = isMAP(proc)
% ISMAP Check if a process representation is a valid MAP {D0, D1}
%
% A valid MAP has exactly 2 cell elements, both of which are square matrices
% of the same size.
result = false;
if ~iscell(proc) || length(proc) ~= 2
    return;
end
D0 = proc{1};
D1 = proc{2};
if ~isnumeric(D0) || ~isnumeric(D1)
    return;
end
if ~ismatrix(D0) || ~ismatrix(D1)
    return;
end
[m0, n0] = size(D0);
[m1, n1] = size(D1);
if m0 ~= n0 || m1 ~= n1 || m0 ~= m1
    return;
end
result = true;
end

function MAP = convertToMAP(station, classIdx, proc)
% CONVERTTOMAP Convert non-Markovian distribution parameters to MAP
%
% For non-Markovian distributions, the process representation contains
% distribution parameters rather than {D0, D1} matrices. This function
% converts them to an Erlang approximation.

% Get the distribution object from the station
if isa(station, 'Source')
    dist = station.input.sourceClasses{classIdx}{end};
else
    dist = station.server.serviceProcess{classIdx}{end};
end

% Get mean for Erlang approximation
targetMean = dist.getMean();

% Determine number of phases based on SCV
% For Det (SCV=0), use high number of phases; for others, match SCV
scv = dist.getSCV();
if scv < GlobalConstants.CoarseTol
    % Deterministic or near-deterministic: use 20 phases
    nPhases = 20;
else
    % Match SCV: for Erlang, SCV = 1/n, so n = 1/SCV
    nPhases = max(1, ceil(1/scv));
    nPhases = min(nPhases, 100); % Cap at 100 phases
end

% Create Erlang MAP approximation
MAP = map_erlang(targetMean, nPhases);
end
