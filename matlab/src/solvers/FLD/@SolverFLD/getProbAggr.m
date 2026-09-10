function varargout = getProbAggr(self,varargin)
% [PNIR,LOGPNIR] = GETPROBAGGR(IST)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for current state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Two evaluations are available and the analysis that ran decides which:
%
%   moment closure ('minnormal', 'refined') -- the solved state carries a
%       covariance, so the JOINT law of the per-class populations at the
%       station is the multivariate normal of the linear noise approximation
%       and the answer is the probability it assigns to the unit cell around
%       n. Correlation between the classes is accounted for, which is the
%       whole point of asking for a joint probability rather than a product
%       of marginals.
%
%   first-order methods -- no second moment exists, so the classes can only
%       be treated as independent: Schmidt's binomial per closed class, and
%       the product form (Poisson at a Delay, multinomial-geometric at a
%       queue) per open class.
%
% Input:
%   ist - Station index
%
% Output:
%   Pnir    - Scalar probability in [0,1]
%   logPnir - Log probability for numerical stability
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probAggr', varargout{1});
end

function [Pnir,logPnir] = getProbAggr_impl(self, ist)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob, the probability that a station holds the
% marginal population of the state the model carries; a prior over several rows
% is refused by name. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbAggr', self.model);
    if nargin < 2
        line_error(mfilename,'getProbAggr requires to pass a parameter the station of interest.');
    end
    if ~isnumeric(ist), istc = ist.index; else, istc = ist; end
    Pnir = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'ProbAggr', istc, self.name, 'getProbAggr');
    if nargout > 1
        logPnir = log(Pnir);
    end
    return
end

if nargin<2 %~exist('ist','var')
    line_error(mfilename,'getProbAggr requires to pass a parameter the station of interest.');
end
if isempty(self.result)
    self.runAnalyzer;
end
% Read the struct AFTER the analysis: sn.state is empty until the model is
% initialized, and State.toMarginal indexes the phase block off its width.
sn = self.getStruct;
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
Q = self.result.Avg.Q;
N = sn.njobs;

state = sn.state{sn.stationToStateful(ist)};
[~, nir, ~, ~] = State.toMarginal(sn, ist, state);

% The moment closure supplies the joint law; a Source is excluded because its
% coordinate is a normalisation constant rather than a population and carries
% no covariance (see FLUID_MOMENT_TERMS).
moments = [];
if isfield(self.result,'solverSpecific') && isstruct(self.result.solverSpecific) ...
        && isfield(self.result.solverSpecific,'moments')
    moments = self.result.solverSpecific.moments;
end
if ~isempty(moments) && isfield(moments,'Sigma') && ~isempty(moments.Sigma) ...
        && isfield(moments,'classBlock') && sn.sched(ist) ~= SchedStrategy.EXT ...
        && ~local_has_open_class(sn, ist, moments)
    [Pnir, logPnir] = local_gaussian_cell(sn, ist, nir(:), Q, moments);
    return
end

openClasses = find(isinf(N));
closedClasses = find(isfinite(N));
logPnir = 0;

% Product-form probability for open classes
if ~isempty(openClasses)
    U = self.result.Avg.U;
    if sn.sched(ist) == SchedStrategy.INF
        % Delay (infinite server): independent Poisson per class
        for r = openClasses
            if Q(ist,r) > 0
                logPnir = logPnir + nir(r)*log(Q(ist,r)) - Q(ist,r) - gammaln(nir(r)+1);
            elseif nir(r) > 0
                logPnir = -Inf;
            end
        end
    elseif sn.sched(ist) ~= SchedStrategy.EXT
        % Queue station: multinomial-geometric product form
        % P(n_1,...,n_R) = (1-rho) * n!/prod(n_r!) * prod(rho_r^n_r)
        % where rho_r = U(ist,r) is per-server utilization
        rho_total = sum(U(ist, openClasses));
        n_total = sum(nir(openClasses));
        if rho_total < 1
            logPnir = logPnir + log(1 - rho_total) + gammaln(n_total + 1);
            for r = openClasses
                rho_r = U(ist, r);
                if nir(r) > 0
                    if rho_r > 0
                        logPnir = logPnir + nir(r)*log(rho_r) - gammaln(nir(r)+1);
                    else
                        logPnir = -Inf;
                    end
                end
            end
        else
            logPnir = -Inf;
        end
    end
end

% Binomial approximation for closed classes
% Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
for r = closedClasses
    logPnir = logPnir + nchoosekln(N(r),nir(r));
    logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
    logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
end

Pnir = real(exp(logPnir));
end

function tf = local_has_open_class(sn, ist, moments)
% Whether an OPEN class is served at station IST.
%
% The Gaussian cell is used only where it beats the alternative. For an open
% class the first-order path is not an independence heuristic but the exact
% product form of the underlying queue -- geometric at a queue, Poisson at a
% Delay -- so replacing it by a normal approximation of the same law would be a
% loss: on M/M/1 at rho = 0.5 the product form is exact where the cell of the
% linear noise approximation returns 0.39 for the empty queue against 0.50.
% The closure earns its place on the CLOSED populations, where the alternative
% is Schmidt's binomial, itself an approximation, and where correlation between
% the classes is real.
tf = false;
for r = 1:sn.nclasses
    if ~isempty(moments.classBlock{ist,r}) && isinf(sn.njobs(r))
        tf = true;
        return
    end
end
end

function [Pnir, logPnir] = local_gaussian_cell(sn, ist, nir, Q, moments)
% Joint probability of the per-class populations at station IST under the
% linear noise approximation solved by SOLVER_FLUID_MOMENTS.
%
% The state coordinates of class r at the station are MOMENTS.CLASSBLOCK{ist,r}
% (one per service phase), so the class population is their sum: its mean is
% the reported Q(ist,r) and the class-to-class covariance is the sum of the
% corresponding block of MOMENTS.SIGMA. The integer count n is then read off
% the continuous law as the unit cell [n-1/2, n+1/2], with the two ends
% extended to infinity at the boundaries of the state space, so that the mass
% the normal puts on negative populations lands on the empty station and the
% mass above a closed population lands on the full one.
K = sn.nclasses;
N = sn.njobs;
Sigma = moments.Sigma;
classBlock = moments.classBlock;

idx = zeros(1,0);
m = zeros(0,1);
a = zeros(0,1);
b = zeros(0,1);
for r = 1:K
    blk = classBlock{ist,r};
    if isempty(blk)
        % the class has no service process here, so it has no coordinate: any
        % positive count is impossible rather than improbable
        if nir(r) > 0
            Pnir = 0; logPnir = -Inf;
            return
        end
        continue
    end
    idx(end+1) = r; %#ok<AGROW>
    m(end+1,1) = Q(ist,r); %#ok<AGROW>
    if nir(r) <= 0
        a(end+1,1) = -Inf; %#ok<AGROW>
    else
        a(end+1,1) = nir(r) - 0.5; %#ok<AGROW>
    end
    if isfinite(N(r)) && nir(r) >= N(r)
        b(end+1,1) = Inf; %#ok<AGROW>
    else
        b(end+1,1) = nir(r) + 0.5; %#ok<AGROW>
    end
end

if isempty(idx)
    Pnir = 1; logPnir = 0;
    return
end

nr = numel(idx);
C = zeros(nr);
for u = 1:nr
    for v = u:nr
        C(u,v) = sum(sum(Sigma(classBlock{ist,idx(u)}, classBlock{ist,idx(v)})));
        C(v,u) = C(u,v);
    end
end

[Pnir, logPnir] = fluid_mvn_rectangle(m, C, a, b);
end
