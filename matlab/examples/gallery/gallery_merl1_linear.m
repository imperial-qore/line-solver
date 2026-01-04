function model = gallery_merl1_linear(n,Umax)
if ~exist('Umax','var')
    Umax = 0.9;
end
if ~exist('n','var')
    n = 2;
end
model = Network('M/Erl/1-Linear');
%% Block 1: nodes
line{1} = Source(model, 'mySource');
for i=1:n
    line{1+i} = Queue(model, ['Queue',num2str(i)], SchedStrategy.FCFS);
end
line{end+1} = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
line{1}.setArrival(oclass, Exp(1));
means = linspace(0.1,Umax,n/2); 
if mod(n,2)==0
    means=[means,means(end:-1:1)];
else
    means=[means,Umax,means(end:-1:1)];
end
for i=1:n            
    line{1+i}.setService(oclass, Erlang.fitMeanAndOrder(means(i),i));
end
%% Block 3: topology
model.link(Network.serialRouting(line));
end
