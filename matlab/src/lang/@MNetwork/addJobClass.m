function addJobClass(self, customerClass)
% ADDJOBCLASS(CUSTOMERCLASS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if self.enableChecks
    if sum(cellfun(@(x) strcmp(x.name,customerClass.name), {self.classes{1:end}}))>0
        line_error(mfilename,sprintf('A class with name %s already exists.\n', customerClass.name));
    end
end
nClasses = length(self.classes);
customerClass.index = nClasses+1;
self.classes{end+1,1} = customerClass;
self.setUsedLangFeature(class(customerClass)); % open or closed
end
