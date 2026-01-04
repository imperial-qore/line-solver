function self = removeClass(self, jobclass)
% SELF = REMOVECLASS(SELF, CLASS)
%
% Remove the specified CLASS from the model

if hasSingleClass(self)
    if self.classes{1}.name == jobclass.name
        line_error(mfilename,'The network has a single class, it cannot be removed from the model.');
    else
        % no changes
    end
else
    nClasses = length(self.classes);
    r = self.getClassByName(jobclass.name).index; % class to remove
    remaining = setdiff(1:nClasses, r);
    if ~isnan(r)
        % check with SEPT/LEPT
        for i=1:length(self.nodes)
            switch class(self.nodes{i})
                case {'Delay','DelayStation','Queue'}
                    self.nodes{i}.schedStrategyPar = self.nodes{i}.schedStrategyPar(remaining);
                    self.nodes{i}.serviceProcess = self.nodes{i}.serviceProcess(remaining);
                    self.nodes{i}.classCap = self.nodes{i}.classCap(remaining);
                    self.nodes{i}.server.serviceProcess = self.nodes{i}.server.serviceProcess(remaining);
                    self.nodes{i}.output.outputStrategy = self.nodes{i}.output.outputStrategy(remaining);
                case 'ClassSwitch'
                    self.nodes{i}.server.updateClassSwitch(self.nodes{i}.server.csFun(remaining,remaining));
                    self.nodes{i}.output.outputStrategy = self.nodes{i}.output.outputStrategy(remaining);
                case 'Cache'
                    line_error(mfilename,'Cannot dynamically remove classes in models with caches. You need to re-instantiate the model.');
                case 'Source'
                    self.nodes{i}.arrivalProcess = self.nodes{i}.arrivalProcess(remaining);
                    self.nodes{i}.classCap = self.nodes{i}.classCap(remaining);
                    self.nodes{i}.input.sourceClasses = self.nodes{i}.input.sourceClasses(remaining);
                    %self.nodes{i}.server.serviceProcess = self.nodes{i}.server.serviceProcess(remaining);
                    self.nodes{i}.output.outputStrategy = self.nodes{i}.output.outputStrategy(remaining);
                case 'Sink'
                    self.nodes{i}.output.outputStrategy = self.nodes{i}.output.outputStrategy(remaining);
            end
        end
        self.classes = self.classes(remaining);
        self.reset(true); % require a complete re-initialization including state
    end
end
end