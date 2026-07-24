% Example 15: Bound analysis with SolverBA
%
% A bounding solver answers a different question from SolverMVA. Instead of a
% single point estimate it returns one guaranteed side of an interval, and the
% .lower/.upper pair of a family brackets the exact solution. Bounds need only
% the service demands and the population, never the service distributions, so
% they are cheap enough to sit inside an optimization loop where a full solve
% would be too slow.
%
% The model is a closed network with a think Delay and two Queues of unequal
% speed, so that the bottleneck is well defined.

%% Block 1: model
N = 5;                                  % number of jobs in the closed chain
model = Network('BoundsDemo');
delay = Delay(model,'Think');
q1 = Queue(model,'Q1', SchedStrategy.PS);
q2 = Queue(model,'Q2', SchedStrategy.PS);
jobs = ClosedClass(model,'C', N, delay);
delay.setService(jobs, Exp(1/2));       % think time  Z = 2
q1.setService(jobs, Exp(1/1.0));        % demand D1 = 1.0
q2.setService(jobs, Exp(1/1.5));        % demand D2 = 1.5  (bottleneck)
model.link(Network.serialRouting(delay,q1,q2));

%% Block 2: the exact reference
% The throughput at the reference station is the system throughput of the
% closed chain, and is what every bound below brackets.
Xexact = MVA(model,'method','exact').getAvgTable().Tput(1);

%% Block 3: the bounds table
% getBoundsTable is the counterpart of getAvgTable: it reports the bracket per
% station and class, with columns Qlower/Qupper and Tlower/Tupper. A single
% call evaluates both sides of the family of the selected method.
gbTable = SolverBA(model,'method','gb.upper').getBoundsTable()

%% Block 4: comparing families
% getBounds is the programmatic accessor, returning the raw bracket rather
% than a formatted table. aba uses only the bottleneck demand and the total
% demand and is the crudest; bjb and gb exploit more structure. Which family
% is sharpest is model dependent, so it is worth comparing several.
fprintf('\n%-6s %12s %12s %12s\n','family','Tlower','Texact','Tupper');
for family = {'aba','bjb','gb'}
    b = SolverBA(model,'method',[family{1} '.upper']).getBounds();
    fprintf('%-6s %12.6f %12.6f %12.6f\n', family{1}, ...
        b.Tlower(1), Xexact, b.Tupper(1));
end

%% Block 5: a bound hierarchy tightening with the level option
% Hierarchical families are parameterized by options.level: raising it spends
% more work and returns a tighter pair. The Eager-Sevcik hierarchy pbh becomes
% exact once the level reaches the population, so the bracket width falls to
% zero at level N.
fprintf('\n%-6s %12s %12s %12s\n','level','Tlower','Tupper','width');
for level = 1:N
    b = SolverBA(model,'method','pbh.upper','level',level).getBounds();
    fprintf('%-6d %12.6f %12.6f %12.6f\n', level, ...
        b.Tlower(1), b.Tupper(1), b.Tupper(1)-b.Tlower(1));
end

%% Block 6: one-sided families
% Not every family is two-sided. cub is upper-only and mbjb and ldbcmp are
% lower-only, so the missing side is reported as NaN rather than as zero,
% which keeps "no bound" distinguishable from "the bound is zero".
cubTable = SolverBA(model,'method','cub.upper').getBoundsTable()

% listValidMethods reports every method name the solver accepts. Some carry
% structural restrictions beyond the feature set: sb and sib are delay-free,
% and lr additionally requires a single-server single-class closed model, so
% on the model above they raise an error rather than return a wrong answer.
fprintf('SolverBA advertises %d methods.\n', ...
    numel(SolverBA(model).listValidMethods()));
