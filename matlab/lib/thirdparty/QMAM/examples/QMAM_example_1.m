
% inter-arrival times (PH distribution)
alpha = [1/2 1/2];
T = [-1 0;
     0  -2];

% service times (PH distribution)
beta = [1 0];
S = [-10 10; 
     0  -10];

[ql,wait_alpha,wait_T]=Q_CT_PH_PH_1(alpha,T,beta,S);

subplot(1,2,1);
bar(ql);
xlabel('Queue length')
ylabel('PMF')

subplot(1,2,2);
x = 0.01:0.01:10;
y = zeros(size(x));
exitT = -sum(wait_T,2); 
for i = 1:length(x);
    y(i) = wait_alpha*expm(wait_T*x(i))*exitT; 
end
plot(x,y);
xlabel('Waiting time')
ylabel('PDF')