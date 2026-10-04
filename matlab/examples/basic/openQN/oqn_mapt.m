% Open queueing network with a MAP_t arrival process and the Ko-Pender limits.
%
% MAPt is a time-inhomogeneous Markovian arrival process: segment k covers
% [breakpoints(k), breakpoints(k+1)) and carries the pair (D0{k}, D1{k}), so the
% stream is both non-renewal, through the modulating phase, and non-stationary,
% through the schedule. Setting one phase recovers an NHPP; setting one segment
% recovers an ordinary MAP.
%
% SolverFLD's 'kp' method integrates the fluid and diffusion limits of Ko and
% Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper.
% Res. Lett. 45 (2017) 248-253: the mean and the covariance of the queue length
% are integrated jointly, so it is the only fluid method that returns a second
% moment. For infinite-server stations the rate functions are affine in the
% state, so both are exact rather than asymptotic -- the variance below is the
% exact variance of the queue length, not an approximation.

model = Network('model');

source = Source(model,'Source');
delay = Delay(model,'Delay');
sink = Sink(model,'Sink');

jobclass = OpenClass(model, 'OpenClass', 0);

% Two segments of a 2-phase MAP, held for 1 and 1.5 time units and repeating.
% The second segment runs the same phase graph at roughly twice the rate.
breakpoints = [0 1 2.5];
D0 = {[-5 1; 2 -4], [-12 3; 5 -9]};
D1 = {[3 1; 1 1],   [7 2; 2 2]};
arrival = MAPt(breakpoints, D0, D1, true);

source.setArrival(jobclass, arrival);
delay.setService(jobclass, Exp(2));

model.link(Network.serialRouting(source,delay,sink));

% Steady state of a cyclic schedule is the average over one period.
% tol below the default 1e-4: the JAR's DormandPrince and MATLAB's ode15s differ
% by ~0.2% at the default, which is integrator tolerance rather than a modelling
% difference and would show up as a cross-codebase disagreement.
AvgTable = SolverFLD(model,'method','kp','tol',1e-9).getAvgTable();
AvgTable

% Transient mean AND variance over two periods.
options = SolverFLD.defaultOptions;
options.method = 'kp';
options.timespan = [0 5];
options.tol = 1e-9;
[~,~,~,~,~,QNt,~,~,~,t,~,~,QVart] = solver_fluid_kp(model.getStruct(), options);
fprintf('\n   t     lambda(t)    QLen mean     QLen var\n');
for probe = [0.5 1.0 1.5 2.5 3.0 4.0 5.0]
    fprintf('%6.2f  %10.4f  %11.6f  %11.6f\n', probe, arrival.getRateAt(probe), ...
        interp1(t, QNt{2,1}, probe), interp1(t, QVart{2,1}, probe));
end
