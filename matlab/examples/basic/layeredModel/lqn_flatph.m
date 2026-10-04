% Method 'flat.ph': the squashed layering with the composed phase-type encoding.
%
% A method name of SolverLN carries TWO decisions. The part before the dot is
% the LAYERING, which fixes what a submodel is; the part after it is the
% ENCODING, which fixes how an activity graph is written into that submodel. The
% four combinations are therefore 'srvn.cs', 'srvn.ph', 'flat.cs' and, as of
% this example, 'flat.ph'.
%
% Under 'flat.ph' ONE submodel holds a station for every processor and every
% called task, and a caller task is ONE closed class that visits each server it
% uses once per invocation, served there by the composed law of the demand it
% places on that server. Squashing does not conflict with the composition: a
% task station's service law is ALREADY the inflated entry law -- host demand
% plus call counts times callee service and waiting -- which the outer fixed
% point updates, exactly as lqns does under --squashed-layering.
%
% The model below is a three-deep call chain, where the middle task is at once a
% server to T1 and a caller of T3, so the two levels of the LQN accounting are
% both exercised. Against lqsim on this model 'flat.ph' reads +0.87 per cent
% where 'flat.cs' reads +2.95, 'srvn.cs' +3.59 and lqns +4.28.

clear solver AvgTable

model = LayeredNetwork('flatphExample');

P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
P3 = Processor(model, 'P3', 1, SchedStrategy.PS);

T1 = Task(model, 'T1', 4, SchedStrategy.REF).on(P1).setThinkTime(Exp(1));
T2 = Task(model, 'T2', 2, SchedStrategy.FCFS).on(P2);
T3 = Task(model, 'T3', 1, SchedStrategy.FCFS).on(P3);

E1 = Entry(model, 'E1').on(T1);
E2 = Entry(model, 'E2').on(T2);
E3 = Entry(model, 'E3').on(T3);

Activity(model, 'A1', Exp(5)).on(T1).boundTo(E1).synchCall(E2, 1);
Activity(model, 'A2', Exp(5)).on(T2).boundTo(E2).synchCall(E3, 1).repliesTo(E2);
Activity(model, 'A3', Exp(5)).on(T3).boundTo(E3).repliesTo(E3);

% The four encodings side by side. 'flat' is the ALIAS of 'flat.cs' and resolves
% unconditionally rather than probing 'flat.ph', because a model is squashed in
% order to express the routed call groups that only the routing encoding
% dispatches; naming 'flat.ph' is therefore not the same as naming 'flat'.
for m = {'srvn.cs', 'srvn.ph', 'flat.cs', 'flat.ph'}
    solver = SolverLN(model, 'method', m{1}, 'verbose', false);
    AvgTable = solver.getAvgTable();
    line_printf('\n--- method=%s (built as %s, %d submodel(s))\n', ...
        m{1}, solver.lnmethod, length(solver.ensemble));
    disp(AvgTable);
end
