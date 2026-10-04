function [binds, node, cap, r, isOpen] = findBindingCapacity(self)
% [BINDS, NODE, CAP, R, ISOPEN] = FINDBINDINGCAPACITY()
%
% The first station whose finite capacity can actually BIND, or BINDS false
% when no buffer in the model can refuse a job.
%
% ONE PREDICATE, TWO CALLERS. NetworkSolver.checkBindingCapacity turns the
% answer into the refusal the product-form solvers raise, and
% getUsedLangFeatures marks the registry name 'FiniteCapacity' on it, so a
% solver that does not declare the name refuses exactly the models the
% structural gate refuses. The test reads the node-level cap / classCap the
% user set and the class populations from the CLASS OBJECTS (getNumberOfJobs),
% never sn.cap / sn.classcap / sn.njobs: refreshCapacity derives a FINITE
% classcap (the chain population) for every closed model, so an sn-level test
% would call every closed model capped, and reading the struct from the
% recorder would trigger a refresh on every feature query.
%
% Only a capacity that can bind counts. A closed model whose station capacity
% is at least the total population can never block a job, so the declaration
% is a no-op (setCap(N) on an order-independent station of an N-job model is a
% common idiom). The population of an open class is Inf, so any finite
% capacity an open class can reach binds.
%
% A Cache model is exempt: Cache.m sets classCap = 1 on the retrieval queues
% it builds, and the cache analyzers solve those rather than treating them as
% a buffer constraint.
%
% NODE is the binding station, CAP the value that binds, R the class index for
% a per-class buffer and 0 for a station-level one, and ISOPEN says whether an
% open class reaches the buffer (station level: any open class in the model;
% class level: class R is open), which decides the fallback advice.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

binds = false;
node = [];
cap = Inf;
r = 0;
isOpen = false;
nodes = self.nodes;
for i = 1:numel(nodes)
    if isa(nodes{i}, 'Cache')
        return
    end
end
njobs = getNumberOfJobs(self);
njobs = njobs(:)';
totalJobs = sum(njobs); % Inf as soon as one class is open
for i = 1:numel(nodes)
    nd = nodes{i};
    if ~isa(nd, 'Station') || isa(nd, 'Source') || isa(nd, 'Sink')
        continue
    end
    if ~isempty(nd.cap) && ~isinf(nd.cap) && nd.cap >= 0 && nd.cap < totalJobs
        binds = true;
        node = nd;
        cap = nd.cap;
        r = 0;
        isOpen = any(isinf(njobs));
        return
    end
    ccap = nd.classCap;
    for k = 1:min(numel(ccap), numel(njobs))
        if ~isinf(ccap(k)) && ccap(k) > 0 && ccap(k) < njobs(k)
            binds = true;
            node = nd;
            cap = ccap(k);
            r = k;
            isOpen = isinf(njobs(k));
            return
        end
    end
end
end
