% LQN_TRANSIENT  Transient (time-dependent) analysis of a layered network.
%
% Demonstrates LN.getTranAvg, which returns the transient mean queue
% length, utilization and throughput of every ensemble layer over time. The
% traces are assembled block-diagonally: layer e occupies a disjoint block of
% rows (its stations) and columns (its classes). Transient output is only
% produced by transient-capable per-layer solvers (Fluid, CTMC, SSA); here the
% layers are solved with the fluid ODE solver. Do not set a 'timespan' on the
% per-layer factory: the steady-state fixed point rejects it, and getTranAvg
% auto-selects the timespan per layer.
%
% getTranAvg first converges the LN fixed point (getAvg), pinning the
% inter-layer demands to equilibrium, then runs each layer's transient with
% those demands frozen. The initial point is the layer's default state (all
% closed jobs at the reference station, via State.initDefault), NOT the
% converged occupancy, so each curve relaxes from all-at-reference to the
% layer steady state. Use setState to seed a different start.
%
% The model is a simple 3-layer LQN: a reference task T1 on processor P1 whose
% activity synchronously calls entry E2 of task T2 on processor P2. The three
% ensemble layers are the two processor (host) layers and the T2 task layer.

model = LayeredNetwork('lqn_transient');

P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
P2 = Processor(model, 'P2', 1, SchedStrategy.PS);

T1 = Task(model, 'T1', 5, SchedStrategy.REF).on(P1).setThinkTime(Exp(1.0));
T2 = Task(model, 'T2', 5, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1.0));

E1 = Entry(model, 'E1').on(T1);
E2 = Entry(model, 'E2').on(T2);

A1 = Activity(model, 'A1', Exp(2.0)).on(T1).boundTo(E1).synchCall(E2, 1);
A2 = Activity(model, 'A2', Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2);

% Solve the ensemble with a fluid solver on each layer, then obtain the
% per-layer transient averages in a single call.
solver = LN(model, @(m) FLD(m, 'verbose', false), 'verbose', false);
[QNt, UNt, TNt] = solver.getTranAvg();

E = solver.nlayers;
fprintf('LN.getTranAvg returned transient traces for %d layers.\n', E);

% Plot the transient mean queue length E[N](t) of each station, one panel per
% layer. The block-diagonal offsets (r0,c0) advance by each layer's station and
% class counts, mirroring how getTranAvg stacks the per-layer blocks.
figure('Name', 'LN.getTranAvg: per-layer transient E[N](t)');
r0 = 0;
c0 = 0;
for e = 1:E
    sn = solver.ensemble{e}.getStruct();
    M = sn.nstations;
    K = sn.nclasses;
    subplot(E, 1, e);
    hold on;
    box on;
    grid on;
    tsettle = 0;
    for i = 1:M
        for r = 1:K
            trace = QNt{r0 + i, c0 + r};
            if isstruct(trace) && isfield(trace, 'metric') && any(trace.metric > 1e-6)
                stname = trace.handle{1}.name;
                clname = trace.handle{2}.name;
                plot(trace.t, trace.metric, 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('%s [%s]', stname, clname));
                % Track when this trace stops changing so the axis can be
                % zoomed onto the transient (getTranAvg auto-selects a wide
                % timespan and the curves are flat once settled).
                m = trace.metric;
                if numel(m) > 1
                    tol = 0.01 * (max(m) - min(m)) + 1e-9;
                    klast = find(abs(m - m(end)) > tol, 1, 'last');
                    if ~isempty(klast)
                        tsettle = max(tsettle, trace.t(min(klast + 1, numel(m))));
                    end
                end
            end
        end
    end
    r0 = r0 + M;
    c0 = c0 + K;
    if tsettle > 0
        xlim([0, 1.3 * tsettle]);
    end
    title(sprintf('Layer %d', e));
    xlabel('time t');
    ylabel('E[N](t)');
    legend('Location', 'best');
end
