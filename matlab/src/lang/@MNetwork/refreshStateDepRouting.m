function sn = refreshStateDepRouting(self)
% SN = REFRESHSTATEDEPROUTING()
%
% Builds sn.sdr, the station-indexed twin of the Krzesinski state-dependent
% routing structure declared on the entry node with Node.setStateDepRouting.
% The node-indexed copy stays in sn.nodeparam{ind}{r}.sdr, where the routing
% functions of refreshRoutingMatrix read it; the station-indexed copy is what
% the product-form API (pfqn_sdr, pfqn_sdrvisits, pfqn_sdrmva) consumes.
%
% sn.sdr is empty when no node declares state-dependent routing. A network may
% declare at most one SDR subnetwork Q(V,V), and every class routed by it must
% declare the same one: the routing probabilities of Krzesinski (1987) are
% chain independent, so a per-class topology has no product form.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.sn;
sn.sdr = [];
if isempty(sn.nodeparam)
    self.sn = sn;
    return
end

decl = [];
declnode = 0;
for ind = 1:sn.nnodes
    if ind > numel(sn.nodeparam) || isempty(sn.nodeparam{ind}) || ~iscell(sn.nodeparam{ind})
        continue
    end
    for r = 1:min(sn.nclasses, numel(sn.nodeparam{ind}))
        if isempty(sn.nodeparam{ind}{r}) || ~isstruct(sn.nodeparam{ind}{r}) || ~isfield(sn.nodeparam{ind}{r},'sdr')
            continue
        end
        cand = sn.nodeparam{ind}{r}.sdr;
        if isempty(decl)
            decl = cand;
            declnode = ind;
        elseif ~isequaln(decl, cand)
            line_error(mfilename, sprintf(['Two different state-dependent routing structures are declared (nodes %s and %s). ', ...
                'The routing probabilities of Krzesinski (1987) are chain independent, so a network admits one subnetwork Q(V,V) ', ...
                'and every class routed by it must declare the same branches, nesting and coefficients.'], ...
                sn.nodenames{declnode}, sn.nodenames{ind}));
        end
    end
end
if isempty(decl)
    self.sn = sn;
    return
end

% Translate node indices into station indices. Every center of an SDR network
% must be a station: the product form is over queue lengths, and a stateless
% node holds none.
    function ist = nd2station(ind)
        ist = sn.nodeToStation(ind);
        if isnan(ist) || ist < 1
            line_error(mfilename, sprintf('Node %s takes part in state-dependent routing but is not a station.', sn.nodenames{ind}));
        end
    end

B = numel(decl.branch);
sdr = struct();
sdr.entry = nd2station(decl.entry);
sdr.departure = nd2station(decl.departure);
sdr.branch = cell(1,B);
sdr.entryOf = zeros(1,B);
sdr.departureOf = zeros(1,B);
for b = 2:B
    sdr.branch{b} = arrayfun(@nd2station, decl.branch{b});
    sdr.entryOf(b) = nd2station(decl.entryOf(b));
    sdr.departureOf(b) = nd2station(decl.departureOf(b));
end
sdr.level = decl.level;
sdr.C = decl.C;
sdr.d = decl.d;
sdr.entrynode = decl.entry;
sdr.departurenode = decl.departure;
sdr.branchnode = decl.branch;

pfqn_sdrcoeff(sdr); % validates the declaration and its population bounds
sn.sdr = sdr;
self.sn = sn;
end
