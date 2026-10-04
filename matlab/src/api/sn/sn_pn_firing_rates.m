%{ @file sn_pn_firing_rates.m
 %  @brief Recovers per-mode transition firing rates of a Petri net from the
 %         Place throughputs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Recovers per-mode transition firing rates from the Place throughputs
 %
 % @details
 % The firing rates of a Petri net are not carried by the network structure,
 % but they are determined by the Place throughputs together with the net
 % structure. Writing x for the vector of per-mode firing rates, two families
 % of equations hold at steady state, for every Place p and class k:
 %
 %   departure    sum over the modes consuming (p,k) of x, weighted by the
 %                input arc multiplicity when TPUTISTOKENS is true and
 %                unweighted when it is false, equals TN(p,k)
 %   balance      sum over all modes of x times (produced minus consumed)
 %                equals zero
 %
 % The system is solved in least squares. That is deliberate: an exact solver
 % supplies throughputs that satisfy it exactly and the fit is then the exact
 % answer, whereas a simulator supplies estimates that satisfy it only up to
 % sampling error and the least-squares fit is the right estimator there. A
 % residual test would reject every simulated run.
 %
 % @par Syntax:
 % @code
 % [x, consumed, produced, placeNodes] = sn_pn_firing_rates(sn, TN, tputIsTokens)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>TN<td>Average throughputs at stations
 % <tr><td>tputIsTokens<td>True when TN counts tokens, false when it counts firing events
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>x<td>Firing rate per (transition, mode) pair, empty when undetermined
 % <tr><td>consumed<td>Tokens consumed, indexed (mode, place, class)
 % <tr><td>produced<td>Tokens produced, indexed (mode, place, class)
 % <tr><td>placeNodes<td>Node indices of the Places, in the order used above
 % </table>
%}
function [x, consumed, produced, placeNodes] = sn_pn_firing_rates(sn, TN, tputIsTokens)

x = [];
consumed = [];
produced = [];

R = sn.nclasses;
placeNodes = find(sn.nodetype == NodeType.Place);
transNodes = find(sn.nodetype == NodeType.Transition);
if isempty(placeNodes) || isempty(transNodes) || isempty(TN)
    return
end

% see _kb/04-networkstruct.md (api/sn derived-field helpers) for rationale
if any(sn.nodetype == NodeType.Source) || any(sn.nodetype == NodeType.Sink)
    return
end
statefulNodes = find(sn.isstateful);
for pp = 1:length(placeNodes)
    sfp = find(statefulNodes == placeNodes(pp), 1);
    if isempty(sfp)
        return
    end
    for sfj = 1:length(statefulNodes)
        if sfj == sfp
            continue
        end
        blockOut = sn.rt((sfp-1)*R+(1:R), (sfj-1)*R+(1:R));
        blockIn = sn.rt((sfj-1)*R+(1:R), (sfp-1)*R+(1:R));
        if (any(blockOut(:) > 0) || any(blockIn(:) > 0)) ...
                && sn.nodetype(statefulNodes(sfj)) ~= NodeType.Transition
            return
        end
    end
end

% Enumerate the (transition, mode) pairs: a mode is what carries a firing
% rate, and a transition may hold several.
modeTrans = [];
modeIdx = [];
modeTimed = [];
for tt = 1:length(transNodes)
    ind = transNodes(tt);
    param = sn.nodeparam{ind};
    if isempty(param) || ~isfield(param, 'nmodes')
        return
    end
    for m = 1:param.nmodes
        modeTrans(end+1) = ind; %#ok<AGROW>
        modeIdx(end+1) = m; %#ok<AGROW>
        % An immediate firing takes zero time and is not a timed event, so
        % the analyzers never count it in TN. Its rate is an unknown to be
        % recovered from the balance equations, not a measured quantity.
        timed = true;
        if isfield(param, 'timing') && numel(param.timing) >= m
            timed = param.timing(m) ~= TimingStrategy.IMMEDIATE;
        end
        modeTimed(end+1) = timed; %#ok<AGROW>
    end
end
modeTimed = logical(modeTimed);
nModes = length(modeTrans);
if nModes == 0
    return
end

consumed = zeros(nModes, length(placeNodes), R);
produced = zeros(nModes, length(placeNodes), R);
for mm = 1:nModes
    param = sn.nodeparam{modeTrans(mm)};
    enab = reshape(param.enabling{modeIdx(mm)}, [], R);
    fire = reshape(param.firing{modeIdx(mm)}, [], R);
    for pp = 1:length(placeNodes)
        pind = placeNodes(pp);
        for k = 1:R
            consumed(mm, pp, k) = max(0, enab(pind, k));
            produced(mm, pp, k) = max(0, fire(pind, k));
        end
    end
end

% see _kb/04-networkstruct.md (api/sn derived-field helpers) for rationale
nEq = 2 * length(placeNodes) * R;
A = zeros(nEq, nModes);
b = zeros(nEq, 1);
row = 0;
nMeasured = 0;
for pp = 1:length(placeNodes)
    ist = sn.nodeToStation(placeNodes(pp));
    for k = 1:R
        if tputIsTokens
            arow = consumed(:, pp, k)';
        else
            arow = double(consumed(:, pp, k)' > 0);
        end
        arow(~modeTimed) = 0;
        if any(arow ~= 0)
            row = row + 1;
            A(row, :) = arow;
            b(row) = TN(ist, k);
            nMeasured = nMeasured + 1;
        end

        row = row + 1;
        A(row, :) = produced(:, pp, k)' - consumed(:, pp, k)';
        b(row) = 0;
    end
end
A = A(1:row, :);
b = b(1:row);

% With no measured row the system is homogeneous and pinv returns the zero
% vector, which would report every Place as idle. Keep what the caller had.
if nMeasured == 0
    return
end

xfit = pinv(A) * b;

% A negative firing rate means the net structure was not read as intended;
% reporting a rate that cannot occur would be worse than reporting nothing.
if any(xfit < -1e-6 * max(1, max(abs(xfit))))
    return
end

x = xfit;
end
