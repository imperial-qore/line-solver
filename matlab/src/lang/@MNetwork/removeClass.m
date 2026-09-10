function self = removeClass(self, jobclass)
% SELF = REMOVECLASS(SELF, CLASS)
%
% Remove the specified CLASS from the model. Each node drops its own
% per-class configuration through Node.removeJobClass and its overrides
% (service, capacity, arrival, class switching), mirroring
% Network.removeClass in the JAR and Network.remove_class in python.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Resolve the class on THIS model: a caller holding the class object of
% another model (ModelAdapter.removeClass works on a copy, whose class
% objects are distinct) must still resolve, hence the lookup by name.
target = self.getClassByName(jobclass.name);
if isempty(target) || (isobject(target) && isnan(target.index))
    return
end

if hasSingleClass(self)
    line_error(mfilename,'The network has a single class, it cannot be removed from the model.');
end

r = target.index;
remaining = setdiff(1:length(self.classes), r);

for i=1:length(self.nodes)
    self.nodes{i}.removeJobClass(target);
end

self.classes = self.classes(remaining);
% RENUMBER the survivors. addJobClass sets index = nClasses+1, so the invariant
% every consumer relies on is classes{k}.index == k; dropping a class without
% restoring it leaves a surviving class claiming a position that no longer
% exists, and Queue.getService/setService then index the one-column service
% table by the OLD position. The JAR and python twins renumber, so this also
% removes a silent three-way divergence.
for k = 1:length(self.classes)
    self.classes{k}.index = k;
end
% The class-switching mask cached by link() is (K x K): leaving it at the
% old size makes getRoutingMatrix return chains over classes that no longer
% exist, and refreshRoutingMatrix then indexes sn.refstat out of range
% instead of reporting a model error.
if ~isempty(self.csMatrix)
    self.csMatrix = self.csMatrix(remaining, remaining);
end
self.reset(true); % require a complete re-initialization including state
end
