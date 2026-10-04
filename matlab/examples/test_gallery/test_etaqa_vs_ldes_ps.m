% Test ETAQA PS departure process: MAM dec.mmap vs LDES
%
% Erl(3)/M/1-PS -> M/1-PS tandem
% Compares MAM (dec.mmap with ETAQA PS departures) against LDES simulation

model = Network('ETAQA-PS-Tandem');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
sink   = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Erlang.fitMeanAndOrder(2, 3));  % mean=2, SCV=1/3
queue1.setService(oclass, Exp(2));    % rho1=0.5
queue2.setService(oclass, Exp(1.5));  % rho2=0.667

model.link(Network.serialRouting(source, queue1, queue2, sink));

%% Solve
solverMMAP = MAM(model, 'method', 'dec.mmap');
T_mmap = solverMMAP.getAvgTable();

solverMVA = MVA(model);
T_mva = solverMVA.getAvgTable();

solverLDES = LDES(model, 'seed', 23000, 'samples', 5e5);
T_ldes = solverLDES.getAvgTable();

%% Display
fprintf('\n=== ETAQA PS: MAM vs LDES ===\n');
fprintf('%-10s %10s %10s %10s\n', 'Station', 'LDES', 'MVA', 'dec.mmap');

fprintf('--- Queue Lengths ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f %10.4f\n', stn, ...
            T_ldes.QLen(i), T_mva.QLen(i), T_mmap.QLen(i));
    end
end

fprintf('--- Response Times ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f %10.4f\n', stn, ...
            T_ldes.RespT(i), T_mva.RespT(i), T_mmap.RespT(i));
    end
end

fprintf('--- Utilizations ---\n');
for i = 1:height(T_mva)
    stn = string(T_mva.Station(i));
    if contains(stn, 'Queue')
        fprintf('%-10s %10.4f %10.4f %10.4f\n', stn, ...
            T_ldes.Util(i), T_mva.Util(i), T_mmap.Util(i));
    end
end
