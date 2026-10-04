function qnet = renv_genqn(rate, N)
%% sn1
qnet = Network('qn1');

node{1} = Delay(qnet, 'Queue1');
node{2} = Queue(qnet, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(qnet, 'Class1', N, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(rate(1)));
node{2}.setService(jobclass{1}, Exp(rate(2)));

K = 1;
P = cell(K,K);
P{1} = circul(length(node));

qnet.link(P);
end
