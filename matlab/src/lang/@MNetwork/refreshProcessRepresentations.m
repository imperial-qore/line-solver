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
        % NHPP carries a rate schedule, not D0/D1: always 1 active phase
        isSchedule_ir = false;
        dist_ir = [];
        if isa(stations{ist}, 'Source') && any(ist == self.getIndexSourceStation) ...
                && length(stations{ist}.input.sourceClasses) >= r ...
                && ~isempty(stations{ist}.input.sourceClasses{r})
            dist_ir = stations{ist}.input.sourceClasses{r}{end};
            isSchedule_ir = ismethod(dist_ir, 'getRateSchedule');
        elseif isa(stations{ist}, 'ServiceStation') ...
                && length(stations{ist}.server.serviceProcess) >= r ...
                && ~isempty(stations{ist}.server.serviceProcess{r})
            dist_ir = stations{ist}.server.serviceProcess{r}{end};
            isSchedule_ir = ismethod(dist_ir, 'getRateSchedule');
        end
        % Gamma/Weibull/Lognormal/Pareto/Uniform return the raw distribution
        % PARAMETERS from getProcess, not a (D0,D1) pair. Two 1x1 parameters
        % satisfy isMAP's "square matrices of equal size" shape test, so they
        % were classified as a valid single-phase MAP and convertToMAP was never
        % reached: sn.proc then held {D0=shape, D1=scale}, phases was 1, and
        % map_pie ran on parameters. The shape test cannot tell a 1x1 (D0,D1)
        % from two scalars, so the process type has to be asked directly. Det
        % returns a single parameter and already fails the length>=2 test.
        isRawParam_ir = ~isempty(dist_ir) && (isa(dist_ir,'Gamma') ...
            || isa(dist_ir,'Weibull') || isa(dist_ir,'Lognormal') ...
            || isa(dist_ir,'Pareto') || isa(dist_ir,'Uniform'));
        if isempty(ph{ist}{r}) % fluid fails otherwise
            phases(ist,r) = 1;
        elseif isSchedule_ir
            phases(ist,r) = 1;
        elseif isRawParam_ir || ~isMAP(ph{ist}{r})
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
            % NHPP carries a rate schedule, not D0/D1: skip map_pie
            isSchedule_ir = false;
            if isa(stations{ist}, 'Source') && any(ist == self.getIndexSourceStation) ...
                    && length(stations{ist}.input.sourceClasses) >= r ...
                    && ~isempty(stations{ist}.input.sourceClasses{r})
                isSchedule_ir = ismethod(stations{ist}.input.sourceClasses{r}{end}, 'getRateSchedule');
            elseif isa(stations{ist}, 'ServiceStation') ...
                    && length(stations{ist}.server.serviceProcess) >= r ...
                    && ~isempty(stations{ist}.server.serviceProcess{r})
                isSchedule_ir = ismethod(stations{ist}.server.serviceProcess{r}{end}, 'getRateSchedule');
            end
            if ~isempty(map_ir)
                proc{ist}{r} = map_ir;
                if isSchedule_ir
                    pie{ist}{r} = NaN;
                else
                    pie{ist}{r} = map_pie(map_ir);
                end
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
    % Marked (MMAP) source classes share the carrier's modulating chain: the
    % non-carrier classes (mark index > 1) contribute a single always-zero
    % state column rather than their own phase block.
    if isfield(self.sn,'markidx') && ~isempty(self.sn.markidx)
        self.sn.phasessz(self.sn.markidx > 1) = 1;
    end
    self.sn.phaseshift = [zeros(size(phases,1),1),cumsum(self.sn.phasessz,2)];
end
end

function result = isMAP(proc)
% ISMAP Check if a process representation is a valid MAP {D0, D1}
%
% A valid representation has at least 2 cell elements (D0, D1, plus
% optional marked/batch matrices D_k as in BMAP/MarkedMAP), all square
% matrices of the same size.
result = false;
if ~iscell(proc) || length(proc) < 2
    return;
end
n0 = size(proc{1}, 1);
for e = 1:length(proc)
    De = proc{e};
    if ~isnumeric(De) || ~ismatrix(De) || size(De,1) ~= size(De,2) || size(De,1) ~= n0
        return;
    end
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
