function [envModel, info] = map2renv(model, options)
% [ENVMODEL, INFO] = MAP2RENV(MODEL, OPTIONS)
%
% Markov-modulated image of a network with MAP/MMPP/MMAP arrival or service
% processes as a queueing network in a random environment.
%
% Every non-renewal process is a point process modulated by the CTMC with
% generator Q = D0 + D1, whose conditional intensity in phase k is
% lambda(k) = sum_j D1(k,j). The transformation freezes each phase into an
% environment stage in which the process is the Poisson process of that
% intensity, i.e. an exponential arrival or service time, and lets the
% environment switch stages at the rates of Q. With P modulated processes the
% stage set is the Cartesian product of their phase spaces and the environment
% generator is the Kronecker sum of the individual Q's, so only one process
% changes phase at a time, as in the original model.
%
% The image is exact in structure for an MMPP (diagonal D1): the modulating
% chain, its stationary distribution and the phase-conditional intensities are
% all preserved. For a general MAP the phase jumps that occur AT an event
% epoch (off-diagonal D1) are aggregated into Q and their correlation with the
% event stream is lost, so the image matches the modulating chain and the
% conditional intensities but not the full inter-event autocorrelation.
%
% Populations are carried across stage switches unchanged (identity reset), as
% a phase switch moves no job.
%
% OPTIONS is a solver options structure; OPTIONS.config.map_env_maxstages caps
% the number of environment stages (default 64).
%
% INFO reports nstages, the per-process phase orders, whether every process was
% an MMPP, and the modulation records themselves.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = Solver.defaultOptions;
end

if ~isa(model, 'Network')
    line_error(mfilename, 'map2renv requires a Network model.');
end

maxStages = 64;
if isfield(options,'config') && isfield(options.config,'map_env_maxstages') ...
        && ~isempty(options.config.map_env_maxstages)
    maxStages = options.config.map_env_maxstages;
end

sn = model.getStruct();
mods = sn_map_modulation(sn);
if isempty(mods)
    line_error(mfilename, 'The model declares no MAP, MMPP2 or MMAP process, so it has no random-environment image.');
end

P = numel(mods);
orders = [mods.order];
nstages = prod(orders);
if nstages > maxStages
    line_error(mfilename, sprintf(['The random-environment image of this model has %d stages (phase orders %s), ' ...
        'above the options.config.map_env_maxstages cap of %d. Reduce the order of the modulating processes ' ...
        'or raise the cap.'], nstages, mat2str(orders), maxStages));
end

% Stage s enumerates the phase tuples in column-major order, phaseOf(s,p) is
% the phase of process p in stage s.
phaseOf = zeros(nstages, P);
for s = 1:nstages
    rem_s = s - 1;
    for p = 1:P
        phaseOf(s,p) = mod(rem_s, orders(p)) + 1;
        rem_s = floor(rem_s / orders(p));
    end
end

envModel = Environment([model.getName, '_renv']);
stageNames = cell(1,nstages);
for s = 1:nstages
    stageNames{s} = stageName(phaseOf(s,:));
    envModel.addStage(stageNames{s}, 'item', buildStage(model, mods, phaseOf(s,:), stageNames{s}));
end

% Kronecker sum of the phase generators: a transition changes the phase of one
% process only, at the rate that process assigns to it.
for s = 1:nstages
    for p = 1:P
        Qp = mods(p).D0 + sumCells(mods(p).D1);
        k = phaseOf(s,p);
        for l = setdiff(1:orders(p), k)
            if Qp(k,l) > GlobalConstants.Zero
                t = s + (l - k) * prod(orders(1:p-1));
                envModel.addTransition(stageNames{s}, stageNames{t}, Exp(Qp(k,l)));
            end
        end
    end
end

envModel.init();

info = struct('nstages',nstages,'orders',orders,'isMMPP',all([mods.isMMPP]),'mods',mods);
end

function name = stageName(phases)
name = 'Phase';
for p = 1:numel(phases)
    name = [name, '_', num2str(phases(p))]; %#ok<AGROW>
end
end

function S = sumCells(C)
S = C{1};
for k = 2:numel(C)
    S = S + C{k};
end
end

function stageNet = buildStage(model, mods, phases, stageName)
% Copy of the base model in which every modulated process is the exponential
% process of its phase-conditional intensity.
stageNet = model.copy();
stageNet.setName([model.getName, '_', stageName]);
classes = stageNet.getClasses();
for p = 1:numel(mods)
    station = stageNet.stations{mods(p).ist};
    for c = 1:numel(mods(p).classes)
        r = mods(p).classes(c);
        rate = sum(mods(p).D1{c}(phases(p), :));
        switch mods(p).kind
            case 'arrival'
                % A silent phase (zero intensity) is an ON/OFF source: keep it
                % as a rate rather than a Disabled class, so that the class
                % still exists in every stage and the rate-averaged limit
                % averages a zero instead of skipping the station.
                station.setArrival(classes{r}, Exp(max(rate, GlobalConstants.Zero)));
            case 'service'
                if rate <= GlobalConstants.Zero
                    line_error(mfilename, sprintf(['Phase %d of the service process of class %d at station %d ' ...
                        'has zero completion rate: the station never empties while the environment sits in that ' ...
                        'stage, so the stage has no steady state and the random-environment image is not ' ...
                        'defined. Model the stalled server as a breakdown stage instead.'], phases(p), r, mods(p).ist));
                end
                station.setService(classes{r}, Exp(rate));
        end
    end
end
% The copy inherited the base model's cached NetworkStruct, so the edits above
% are invisible until the struct is rebuilt. see _kb/06-solver-catalog.md
stageNet.refreshStruct(true);
end
