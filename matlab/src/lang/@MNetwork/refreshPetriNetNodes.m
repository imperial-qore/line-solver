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
            for m = 1:self.sn.nodeparam{ind}.nmodes
                self.sn.nodeparam{ind}.enabling{m} = node.enablingConditions{m};
                self.sn.nodeparam{ind}.inhibiting{m} = node.inhibitingConditions{m};
                self.sn.nodeparam{ind}.firing{m} = node.firingOutcomes{m};
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
    end
end
end