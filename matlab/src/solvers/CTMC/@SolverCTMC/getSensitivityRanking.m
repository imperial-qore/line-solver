function RankTable = getSensitivityRanking(self, params, reward)
% RANKTABLE = GETSENSITIVITYRANKING(PARAMS, REWARD)
%
% Rank model parameters by their influence on a steady-state reward, as in
% Trivedi and Bobbio (2017), Table 9.3.
%
% Both the unscaled sensitivity of Eq. (9.79) and the scaled sensitivity of
% Eq. (9.80) are reported. The ranking is by descending absolute scaled
% sensitivity, since that is the comparison the book makes: the scaled form
% is dimensionless, so it is the one that can be compared across parameters
% measured in different units. The sign is retained in the table because it
% says whether increasing a parameter helps or hurts.
%
% @param params Cell array of parameter structs, see SolverCTMC.getSensitivity
% @param reward Reward rate vector or handle over the state space
% @return RankTable Table with columns Parameter, Value, Sens, ScaledSens, sorted
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~iscell(params)
    line_error(mfilename, 'params must be a cell array of parameter structs');
end
if nargin < 3 || isempty(reward)
    line_error(mfilename, 'a reward is required to rank parameters');
end

L = length(params);
Parameter = cell(L, 1);
Value = zeros(L, 1);
Sens = zeros(L, 1);
ScaledSens = zeros(L, 1);

for l = 1:L
    p = params{l};
    if isfield(p, 'name') && ~isempty(p.name)
        Parameter{l} = p.name;
    else
        Parameter{l} = sprintf('theta%d', l);
    end
    [S, SS] = self.getSensitivity(p, reward);
    Value(l) = p.value;
    Sens(l) = S;
    ScaledSens(l) = SS;
end

[~, ord] = sort(abs(ScaledSens), 'descend');
Parameter = categorical(Parameter(ord));
Value = Value(ord);
Sens = Sens(ord);
ScaledSens = ScaledSens(ord);

RankTable = table(Parameter, Value, Sens, ScaledSens);
end
