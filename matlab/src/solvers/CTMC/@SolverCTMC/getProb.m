function Pnir = getProb(self, node, state)
% PNIR = GETPROB(NODE, STATE)

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

if nargin<2 %~exist('node','var')
    line_error(mfilename,'getProb requires to pass a parameter the station of interest.');
end
if ~isfield(self.options,'keep')
    self.options.keep = false;
end
T0 = tic;
sn = self.getStruct;
sn.state = sn.state;
if nargin>=3 %exist('state','var')
    sn.state{node} = state;
end
ind = node.index;
for isf=1:length(sn.state)
    isf_param = sn.nodeToStateful(ind);
    if isf ~= isf_param
        sn.state{isf} = sn.state{isf}*0 -1;
    end
end
Pnir = solver_ctmc_marg(sn, self.options);
self.result.('solver') = getName(self);
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(node);
end