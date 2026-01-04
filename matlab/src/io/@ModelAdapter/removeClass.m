function newmodel = removeClass(model, jobclass)
% SELF = REMOVECLASS(SELF, CLASS)
%
% Remove the specified CLASS from the model
newmodel = model.copy();

if hasSingleClass(newmodel)
    if newmodel.classes{1}.name == jobclass.name
        line_error(mfilename,'The network has a single class, it cannot be removed from the model.');
    else
        % no changes
    end
else
    nClasses = length(newmodel.classes);
    r = newmodel.getClassByName(jobclass.name).index; % class to remove
    remaining = setdiff(1:nClasses, r);
    if ~isnan(r)
        % check with SEPT/LEPT
        for i=1:length(newmodel.nodes)
            switch class(newmodel.nodes{i})
                case {'Delay','DelayStation','Queue'}
                    newmodel.nodes{i}.schedStrategyPar = newmodel.nodes{i}.schedStrategyPar(remaining);
                    newmodel.nodes{i}.serviceProcess = newmodel.nodes{i}.serviceProcess(remaining);
                    newmodel.nodes{i}.classCap = newmodel.nodes{i}.classCap(remaining);
                    newmodel.nodes{i}.server.serviceProcess = newmodel.nodes{i}.server.serviceProcess(remaining);
                    newmodel.nodes{i}.output.outputStrategy = newmodel.nodes{i}.output.outputStrategy(remaining);
                case 'ClassSwitch'
                    newmodel.nodes{i}.server.updateClassSwitch(newmodel.nodes{i}.server.csFun(remaining,remaining));
                    newmodel.nodes{i}.output.outputStrategy = newmodel.nodes{i}.output.outputStrategy(remaining);
                case 'Cache'
                    line_error(mfilename,'Cannot dynamically remove classes in models with caches. You need to re-instantiate the model.');
                case 'Source'
                    newmodel.nodes{i}.arrivalProcess = newmodel.nodes{i}.arrivalProcess(remaining);
                    newmodel.nodes{i}.classCap = newmodel.nodes{i}.classCap(remaining);
                    newmodel.nodes{i}.input.sourceClasses = newmodel.nodes{i}.input.sourceClasses(remaining);
                    %self.nodes{i}.server.serviceProcess = self.nodes{i}.server.serviceProcess(remaining);
                    newmodel.nodes{i}.output.outputStrategy = newmodel.nodes{i}.output.outputStrategy(remaining);
                case 'Sink'
                    newmodel.nodes{i}.output.outputStrategy = newmodel.nodes{i}.output.outputStrategy(remaining);
            end
        end
        newmodel.classes = newmodel.classes(remaining);
        newmodel.reset(true); % require a complete re-initialization including state
    end
end
end