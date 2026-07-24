function Pnir = getProb(self, node, state)
% PNIR = GETPROB(NODE, STATE) Returns state probabilities at equilibrium
%
% @brief Returns state probabilities for a given node, including phase information
%
% This method computes the steady-state marginal state probability distribution
% for jobs at a specified station. Unlike getProb_aggr(), this method returns
% full state-space probabilities including phase information from service distributions.
%
% @param self SolverCTMC instance
% @param node Queue or node object (or node index as integer)
% @param state (optional) State specification as vector. If provided, returns
%             probability of this specific state. If omitted, returns probabilities
%             for all states at the node.
%
% @return Pnir Matrix of state probabilities. If state parameter was provided,
%         returns a scalar probability value. Otherwise returns a matrix where
%         rows represent different states and columns represent job classes.
%         For phase-type distributions, includes phase dimension in state space.
%
% @note Only supported by CTMC and SSA solvers due to their state-space formulation.
%       For aggregated probabilities (without phase info), use getProb_aggr() instead.
%
% @see getProb_aggr - Returns probabilities aggregated over phases
% @see getProbSys - Returns joint system state probabilities
%
% Example:
% @code
% solver = SolverCTMC(model);
% queue1 = model.nodes{1};
%
% % Get all state probabilities for queue1
% prob_matrix = solver.getProb(queue1);
%
% % Get probability for specific state [2 jobs of class 0, 1 job of class 1]
% queue1.setState([2, 1]);
% prob_specific = solver.getProb(queue1);
% @endcode

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