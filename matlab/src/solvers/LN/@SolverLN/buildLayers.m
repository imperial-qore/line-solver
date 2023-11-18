function buildLayers(self)
lqn = self.lqn;
self.ensemble = cell(lqn.ntasks,1);

self.servt_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.call_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.arvproc_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.thinkt_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.route_prob_updmap = cell(lqn.nhosts+lqn.ntasks,1);

%% build one subnetwork for every processor
for hidx = 1:lqn.nhosts
    callers = lqn.tasksof{hidx};
    self.buildLayersRecursive(hidx, callers, true);
end

%% build one subnetwork for every task
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if ~lqn.isref(tidx) & ~(isempty(find(self.lqn.iscaller(tidx,:), 1)) & isempty(find(self.lqn.iscaller(:,tidx), 1)))  %#ok<OR2,AND2> % ignore isolated tasks and ref tasks
        % obtain the activity graph of each task that calls some entry in t
        [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
        callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
        if ~isempty(callers) % true if the server is a software task
            self.buildLayersRecursive(tidx, callers, false);
        else
            self.ensemble{tidx} = [];
        end
    else
        self.ensemble{tidx} = [];
    end
end

self.thinkt_classes_updmap = cell2mat(self.thinkt_classes_updmap);
self.call_classes_updmap = cell2mat(self.call_classes_updmap);
self.servt_classes_updmap = cell2mat(self.servt_classes_updmap);
self.arvproc_classes_updmap = cell2mat(self.arvproc_classes_updmap);
self.route_prob_updmap = cell2mat(self.route_prob_updmap);

% we now calculate the new index of the models after removing the empty
% models associated to 'ref' tasks
emptymodels = cellfun(@isempty,self.ensemble);
self.ensemble(emptymodels) = [];
self.idxhash = [1:length(emptymodels)]' - cumsum(emptymodels);
self.idxhash(emptymodels) = NaN;

self.model.ensemble = self.ensemble;
end
