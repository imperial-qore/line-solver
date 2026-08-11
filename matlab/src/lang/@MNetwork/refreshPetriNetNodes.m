function refreshPetriNetNodes(self)
% REFRESHPETRINETNODES()

for ind=1:self.getNumberOfNodes
    node = self.getNodeByIndex(ind);
    switch class(node)
        case 'Place'
            % noop
        case 'Transition'
            self.sn.nodeparam{ind}.nmodes = length(node.modeNames);
            self.sn.nodeparam{ind}.modenames = node.modeNames;
            self.sn.nodeparam{ind}.enabling = {};
            self.sn.nodeparam{ind}.inhibiting = {};
            self.sn.nodeparam{ind}.firing = {};
            nnodes = self.getNumberOfNodes();
            nclasses = self.getNumberOfClasses();
            for m = 1:self.sn.nodeparam{ind}.nmodes
                % Pad enabling/inhibiting/firing to (nnodes x nclasses) since
                % addMode initializes with node count at creation time
                en = node.enablingConditions{m};
                inh = node.inhibitingConditions{m};
                fir = node.firingOutcomes{m};
                if size(en,1) < nnodes || size(en,2) < nclasses
                    en_full = zeros(nnodes, nclasses);
                    en_full(1:size(en,1), 1:size(en,2)) = en;
                    en = en_full;
                end
                if size(inh,1) < nnodes || size(inh,2) < nclasses
                    inh_full = Inf*ones(nnodes, nclasses);
                    inh_full(1:size(inh,1), 1:size(inh,2)) = inh;
                    inh = inh_full;
                end
                if size(fir,1) < nnodes || size(fir,2) < nclasses
                    fir_full = zeros(nnodes, nclasses);
                    fir_full(1:size(fir,1), 1:size(fir,2)) = fir;
                    fir = fir_full;
                end
                self.sn.nodeparam{ind}.enabling{m} = en;
                self.sn.nodeparam{ind}.inhibiting{m} = inh;
                self.sn.nodeparam{ind}.firing{m} = fir;
            end
            self.sn.nodeparam{ind}.nmodeservers = node.numberOfServers;
            self.sn.nodeparam{ind}.firingprio = node.firingPriorities;
            self.sn.nodeparam{ind}.fireweight = node.firingWeights;
            self.sn.nodeparam{ind}.timing = node.timingStrategies;
            for m = 1:self.sn.nodeparam{ind}.nmodes
                if isa(node.distributions{m},'Markovian')
                    self.sn.nodeparam{ind}.firingproc{m} = node.distributions{m}.getProcess;
                    self.sn.nodeparam{ind}.firingpie{m} = node.distributions{m}.getInitProb;
                    self.sn.nodeparam{ind}.firingphases(m) = node.distributions{m}.getNumberOfPhases;
                elseif isa(node.distributions{m},'ContinuousDistribution')
                    % For non-Markovian distributions, store actual parameters
                    self.sn.nodeparam{ind}.firingproc{m} = node.distributions{m}.getProcess();
                    self.sn.nodeparam{ind}.firingpie{m} = {};
                    self.sn.nodeparam{ind}.firingphases(m) = NaN;
                else
                    % Fallback for other distribution types
                    self.sn.nodeparam{ind}.firingproc{m} = {};
                    self.sn.nodeparam{ind}.firingpie{m} = {};
                    self.sn.nodeparam{ind}.firingphases(m) = NaN;
                end
                self.sn.nodeparam{ind}.firingprocid(m) = ProcessType.toId(ProcessType.fromText(class(node.distributions{m})));
            end
            % Marking-dependent firing-rate multiplier g_m(marking) per mode;
            % empty entry == unit multiplier. See Transition.setFiringRateDependence.
            self.sn.nodeparam{ind}.firingdep = cell(1, self.sn.nodeparam{ind}.nmodes);
            for m = 1:self.sn.nodeparam{ind}.nmodes
                if numel(node.firingRateDependence) >= m
                    self.sn.nodeparam{ind}.firingdep{m} = node.firingRateDependence{m};
                end
            end
    end
end
end