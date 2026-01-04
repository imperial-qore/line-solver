function bool = hasInitState(self)
% BOOL = HASINITSTATE()

bool = true;
if ~self.hasState % check if all stations are initialized
    for ind=1:self.getNumberOfNodes
        if isa(self.nodes{ind},'StatefulNode') && isempty(self.nodes{ind}.state)
            bool = false;
            break;
        end
    end
end
end