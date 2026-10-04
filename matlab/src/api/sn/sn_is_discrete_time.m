%{ @file sn_is_discrete_time.m
 %  @brief Decides whether a model lives on a discrete (slotted) time scale
 %
 %  @author LINE Development Team
%}

%{
 % @brief Decides whether a model lives on a discrete (slotted) time scale
 %
 % @details
 % A model is discrete-time when every enabled interarrival and service law is
 % lattice-valued on a common slot length d, and at least one of them is
 % intrinsically discrete. The lattice families are Geometric (support
 % {1,2,...} slots), DMAP, DiscreteUniform with integral bounds, and Det whose
 % value is a positive integral number of slots. Immediate is deliberately not
 % a lattice law: a zero interval is not a point of {d,2d,...}, which is the
 % same refusal the LDES slotted engine makes in slotSnap.
 %
 % The test runs on sn.procid, sn.rates and sn.scv rather than on sn.proc,
 % because refreshProcessRepresentations already replaced a Geometric by a
 % continuous MAP fit through convertToMAP. procid keeps the requested family
 % and (mean, SCV) identify the member of it exactly for every family above,
 % so the discrete law can be rebuilt losslessly. DMAP is the exception: its
 % (D0,D1) pair survives verbatim in sn.proc because it already has MAP shape.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_discrete_time(sn)
 % [bool, slotLength, info] = sn_is_discrete_time(sn, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>options<td>(Optional) solver options; reads config.timescale
 %                    ('auto' default, 'discrete', 'continuous') and
 %                    config.slotlength (default 1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True when the model is discrete-time on slotLength
 % <tr><td>slotLength<td>Slot length in model time units
 % <tr><td>info<td>Struct with fields hasLattice, hasContinuous, mixed,
 %                 hasDmap and reason (empty when bool is true)
 % </table>
%}
function [bool, slotLength, info] = sn_is_discrete_time(sn, options)

if nargin < 2
    options = struct();
end

timescale = 'auto';
slotLength = 1;
if isfield(options, 'config') && ~isempty(options.config)
    if isfield(options.config, 'timescale') && ~isempty(options.config.timescale)
        timescale = lower(char(options.config.timescale));
    end
    if isfield(options.config, 'slotlength') && ~isempty(options.config.slotlength)
        slotLength = options.config.slotlength;
    end
end

if ~any(strcmp(timescale, {'auto','discrete','continuous'}))
    line_error(mfilename, 'options.config.timescale must be ''auto'', ''discrete'' or ''continuous''.');
end
if ~isnumeric(slotLength) || ~isscalar(slotLength) || ~isfinite(slotLength) || slotLength <= 0
    line_error(mfilename, 'options.config.slotlength must be a positive finite scalar.');
end

info = struct('hasLattice', false, 'hasContinuous', false, 'mixed', false, ...
    'hasDmap', false, 'reason', '');

if strcmp(timescale, 'continuous')
    bool = false;
    return;
end

M = sn.nstations;
K = sn.nclasses;
tol = GlobalConstants.FineTol;

for ist = 1:M
    for r = 1:K
        procType = sn.procid(ist, r);
        if isnan(procType) || procType == ProcessType.DISABLED
            continue;
        end
        rate = sn.rates(ist, r);
        if isnan(rate) || rate <= 0
            % a disabled or absent class carries no interval law
            continue;
        end
        meanSlots = 1 / (rate * slotLength);
        switch procType
            case ProcessType.GEOMETRIC
                % mean 1/p slots; p in (0,1] is recovered exactly from the mean
                info.hasLattice = true;
                if meanSlots < 1 - tol
                    info.hasContinuous = true;
                    info.reason = sprintf(['Geometric at station %d class %d has mean %g slots, ' ...
                        'below the one-slot minimum of its support {1,2,...}.'], ist, r, meanSlots);
                end
            case ProcessType.DMAP
                info.hasLattice = true;
                info.hasDmap = true;
            case ProcessType.DUNIFORM
                info.hasLattice = true;
                [lo, hi] = duniform_bounds(meanSlots, sn.scv(ist, r));
                if isnan(lo) || lo < 1 - tol
                    info.hasContinuous = true;
                    info.reason = sprintf(['DiscreteUniform at station %d class %d spans [%g,%g] slots, ' ...
                        'which is not contained in {1,2,...}.'], ist, r, lo, hi);
                end
            case ProcessType.DET
                if abs(meanSlots - round(meanSlots)) <= tol * max(1, meanSlots) && round(meanSlots) >= 1
                    info.hasLattice = true;
                else
                    % a Det off the lattice is what makes the model continuous
                    info.hasContinuous = true;
                end
            otherwise
                info.hasContinuous = true;
        end
    end
end

info.mixed = info.hasLattice && info.hasContinuous;

if strcmp(timescale, 'discrete')
    if info.mixed
        line_error(mfilename, ['options.config.timescale=''discrete'' was requested but the model ' ...
            'mixes lattice and non-lattice laws. %s'], info.reason);
    end
    if ~info.hasLattice
        line_error(mfilename, ['options.config.timescale=''discrete'' was requested but no ' ...
            'interarrival or service law is lattice-valued on a slot of %g.'], slotLength);
    end
    bool = true;
    return;
end

bool = info.hasLattice && ~info.hasContinuous;

if ~bool && info.hasDmap
    % A DMAP has no continuous-time reading: its (D0,D1) are probability
    % matrices, so the CTMC machinery would compute inv(-D0) where the law
    % needs inv(I-D0) and return a wrong number in silence.
    line_error(mfilename, ['The model mixes a DMAP with continuous-time laws. A DMAP is only ' ...
        'defined on a slotted time scale, so no solver can interpret this model. %s'], info.reason);
end

if ~bool && isempty(info.reason) && info.mixed
    info.reason = 'the model mixes lattice-valued and continuous laws';
end

end

%% Recovers the integral bounds of a DiscreteUniform from its mean and SCV
function [lo, hi] = duniform_bounds(meanSlots, scv)
if isnan(scv) || scv < 0
    lo = NaN; hi = NaN;
    return;
end
varSlots = scv * meanSlots^2;
% var = ((hi-lo+1)^2-1)/12 for integral bounds
width = sqrt(max(0, 12 * varSlots + 1)) - 1;
lo = round(meanSlots - width / 2);
hi = round(meanSlots + width / 2);
end
