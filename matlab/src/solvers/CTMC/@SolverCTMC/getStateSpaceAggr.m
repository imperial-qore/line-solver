function stateSpaceAggr = getStateSpaceAggr(self)
% STATESPACEAGGR = GETSTATESPACEAGGR()

options = self.getOptions;
if options.force
    self.run;
end
if isempty(self.result) || ~isfield(self.result,'spaceAggr')
    line_warning(mfilename,'The model has not been cached. Either solve it or use the ''force'' option to require this is done automatically, e.g., SolverCTMC(model,''force'',true).getStateSpaceAggr()\n');
    stateSpaceAggr = [];
else
    stateSpaceAggr = self.result.spaceAggr;
end
end