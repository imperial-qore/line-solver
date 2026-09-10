function lqn_server_pools()
% LQN_SERVER_POOLS  A processor whose servers are not interchangeable.
%
% P1 declares three servers, but they are not a homogeneous pool: S1 is
% dedicated to task T2, S3 to task T3, and only S2 can take either. Neither task
% can therefore reach more than two of the three servers, and the model is a
% different system from a plain multiplicity-3 processor even though it holds
% the same number of servers.
%
%   P1 (3 servers)          S1 --- T2
%                           S2 --< T2, T3
%                           S3 --- T3
%
% The pools are declared with ServerType, the same class the heterogeneous
% queueing station uses, with the compatible entities being the OPERANDS of the
% layered server: the tasks of a processor, or the entries of a task.
%
% SolverLN lowers the declaration to the activated-server rate of
% SN_COMPAT_RATE, carried onto the layer station as a joint dependence, so a
% compatibility declaration is an APPROXIMATION inside a layer and is admitted
% only under the class-switching layerings ('srvn.cs', 'flat.cs').
%
% Reference:
%   Dorsman, Gardner (2024). New directions in pass-and-swap queues. Queueing
%   Systems 107(3):205-256, Fig. 1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

opt = SolverLN.defaultOptions;
opt.method = 'srvn.cs';

for variant = {'homogeneous','compatibility'}
    kind = variant{1};

    model = LayeredNetwork('LQNserverPools');
    P0 = Processor(model, 'P0', 1, SchedStrategy.PS);
    P1 = Processor(model, 'P1', 3, SchedStrategy.PS);
    T1 = Task(model, 'T1', 3, SchedStrategy.REF).on(P0).setThinkTime(Exp(1.0));
    T2 = Task(model, 'T2', 3, SchedStrategy.FCFS).on(P1);
    T3 = Task(model, 'T3', 3, SchedStrategy.FCFS).on(P1);
    E1 = Entry(model, 'E1').on(T1);
    E2 = Entry(model, 'E2').on(T2);
    E3 = Entry(model, 'E3').on(T3);
    Activity(model, 'A1', Exp(2.0)).on(T1).boundTo(E1).synchCall(E2,1).synchCall(E3,1);
    Activity(model, 'A2', Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2);
    Activity(model, 'A3', Exp(2.0)).on(T3).boundTo(E3).repliesTo(E3);

    switch kind
        case 'homogeneous'
            % one pool of three, every task eligible on every server: this is
            % the neutral declaration and reproduces the plain multiserver
            P1.addServerType(ServerType('All', 3, [T2 T3]));
        case 'compatibility'
            P1.addServerType(ServerType('S1', 1, [T2]));      % dedicated to T2
            P1.addServerType(ServerType('S2', 1, [T2 T3]));   % shared
            P1.addServerType(ServerType('S3', 1, [T3]));      % dedicated to T3
    end

    fprintf('\n--- %s pool on P1 ---\n', kind);
    LN(model, @(m) MVA(m), opt).getAvgTable
end
end
