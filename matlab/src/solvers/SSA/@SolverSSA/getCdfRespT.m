function RD = getCdfRespT(self, R) %#ok<INUSD>
% RD = GETCDFRESPT(R)
%
% Not available: SolverSSA does not record per-job response times.
%
% A simulator must report what it measured. The inherited NetworkSolver
% implementation fabricates an exponential law with the right mean, which
% carries no information about the tail and would be indistinguishable, to the
% caller, from a measured distribution. SSA samples state trajectories, not
% per-job sojourn times, so there is nothing to build an empirical CDF from.
%
% Use SolverJMT (whose getCdfRespT is the ecdf of the logged per-job response
% times), or getPerctRespT(...,'forktail') for the analytical fork-join tail.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


% lang='cpp' cannot serve this getter; the reason is named, not blanket.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.cppUnsupported(self.name, 'getCdfRespT', ...
        ['SolverSSA records no per-job response times on either side, so an empirical CDF ' ...
        'would be an inherited exponential fit carrying no tail information while looking ' ...
        'measured; line-cli refuses -s ssa -a cdf for the same reason. Use SolverJMT for a ' ...
        'measured CDF, SolverFLD or SolverMAM for an analytical one']);
end

line_error(mfilename, ['SolverSSA does not record per-job response times, so it cannot return an ' ...
    'empirical response time CDF. Use SolverJMT for a measured CDF, or ' ...
    'getPerctRespT(...,''forktail'') for the analytical fork-join tail.']);
end
