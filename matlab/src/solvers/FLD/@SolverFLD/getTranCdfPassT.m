function RD = getTranCdfPassT(self, R)
% RD = GETTRANCDFPASST(R)


% lang='cpp' cannot serve this getter; the reason is named, not blanket.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.cppUnsupported(self.name, 'getTranCdfPassT', ...
        ['the C++ fluid port carries no transient passage-time law']);
end

T0 = tic;
if nargin<2 %~exist('R','var')
    R = self.getAvgRespTHandles;
end
sn = self.getStruct;
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        if nnz(sn.stateprior{isf})>1
            line_error(mfilename,'getTranCdfPassT: multiple initial states have non-zero prior - unsupported.');
        end
        sn.state{isf} = sn.state{isf}(1,:); % assign initial state to network
    end
end
options = self.getOptions;
[odeStateVec] = solver_fluid_initsol(sn, options);
options.init_sol = odeStateVec;
RD = solver_fluid_passage_time(sn, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end