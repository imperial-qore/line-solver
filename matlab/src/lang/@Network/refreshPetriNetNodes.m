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
            self.sn.nodeparam{ind}.nmodeservers = node.numbersOfServers;
            self.sn.nodeparam{ind}.firingprio = node.firingPriorities;
            self.sn.nodeparam{ind}.fireweight = node.firingWeights;
            self.sn.nodeparam{ind}.firingid = node.timingStrategies;
            for m = 1:self.sn.nodeparam{ind}.nmodes
                self.sn.nodeparam{ind}.firingproc{m} = node.distributions{m}.getRepres;
                self.sn.nodeparam{ind}.firingprocid(m) = ProcessType.toId(ProcessType.fromText(class(node.distributions{m})));
                self.sn.nodeparam{ind}.firingphases(m) = node.distributions{m}.getNumberOfPhases;
            end
    end
end
end