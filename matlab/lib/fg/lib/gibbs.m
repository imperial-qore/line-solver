function demand = gibbs(data,nbCores,tol)

if exist('tol','var') == 0
    tol = 10^-3;
end

alg = 'TE'; % or MCI

data_needed = 200000;
likelihood_sample = 5000;
nbSamples = 2000;

nbClasses = size(data,2)-1;
nbNodes = 2;
nbJobs = zeros(1,nbClasses);
[prob, nbJobs, N0] = analyseData(data, nbJobs, nbClasses, nbNodes, data_needed);

usedCores = 0;
for k = 1:size(prob,1)
    if (sum(prob(k,nbClasses+1:nbClasses*2)) > nbCores)
        usedCores = usedCores + nbCores*prob(k,end);
    else
        usedCores = usedCores + sum(prob(k,nbClasses+1:nbClasses*2))*prob(k,end);
    end
end
usedCores = usedCores/(1-prob(end,end));

think_time = zeros(1,nbClasses);
for k = 1:nbClasses
    think_time(k) = (nbJobs(k)-N0(k))/mean(data{6,k});
end

range_size = ones(1,nbClasses*(nbNodes-1));

%         testset = [];
%         for k = 1:size(prob,1)
%             testset = [testset; repmat(prob(k,1:(nbClasses*nbNodes)),round(prob(k,end)*likelihood_sample),1)];
%         endinterval

cum_prob = cumsum(prob(:,nbClasses*nbNodes+1));
testset = zeros(likelihood_sample,nbClasses*nbNodes);
for k = 1:likelihood_sample
    uni_value = rand(1);
    index = find(uni_value<cum_prob);
    testset(k,:) = prob(index(1),1:(nbClasses*nbNodes));
end

LV(1) = 0; %log(0!)
LV(2) = 0; %log(1!)
for k = 3:sum(nbJobs)+1
    LV(k) = LV(k-1)+log(k-1);
end

A=feval(@(x) LV(x+1), testset);
sumA = sum(A(:));

initial = zeros(1,nbClasses*nbNodes-nbClasses);

logG_initial = sum(nbJobs.*log(think_time));
for k = 1:nbClasses
    logG_initial = logG_initial - sum(log(1:nbJobs(k)));
end

smpl = zeros(nbSamples,nbClasses*(nbNodes-1));
sample_index = 0;

for k = 1:round(nbSamples/50)
    for s = 1:50
        sample_index = sample_index + 1;
        for h = 1:nbClasses*(nbNodes-1)
            if sample_index==1
                theta = [smpl(sample_index,1:h-1),initial(h:end)];
            else
                theta = [smpl(sample_index,1:h-1),smpl(sample_index-1,h:end)];
            end
            
            [smpl(sample_index,h), logG_initial, range_size_dim]= gibbsSamplerSimple(alg,think_time,theta,testset,h,nbNodes,nbClasses,nbJobs,logG_initial,tol,range_size(h),LV,sumA);
            
            range_size(h) = range_size_dim*2;
            
        end
    end
    
    if  k == 2
        demand_old = mean(smpl(51:sample_index,:));
    elseif k > 2
        demand_now = mean(smpl((k-1)*50+1:sample_index,:));
        demand_now = demand_now/(k+1)+demand_old/(k+1)*k;
        if mean(abs((demand_now-demand_old)./demand_old)) < tol
            nbSample =  sample_index-1;
            N = round(nbSample/2)+1;
            demand = mean(smpl(N:nbSample,:)*usedCores);
            return
        else
            demand_old = demand_now;
        end
    end
    
end
nbSample =  sample_index-1;
N = round(nbSample/2)+1;
demand = mean(smpl(N:nbSample,:)*usedCores);

end

function [prob_nbCustomer, N, N0] = analyseData( data, nbJobs, nbClasses, nbNodes, data_needed)

%number of customer classes, start from 1.
K = nbClasses;

%total number of jobs in the system
N = nbJobs;
N0 = zeros(1,K);

tempTS=[];
tempClass=[];
tempLogger=[];
for i = 1:K
    temp_length = size(data{3,i},1);
    tempTS = [tempTS;data{3,i};data{3,i}+data{4,i}*1000];
    tempClass = [tempClass;ones(temp_length*2,1)*i];
    tempLogger = [tempLogger;ones(temp_length,1);ones(temp_length,1)*2];
end

[ts index] = sort(tempTS);
class_id = tempClass(index);
logger_id = tempLogger(index);

burnin = length(ts)-data_needed;

if burnin < 0 || data_needed == 0
    burnin = 1;
end

%Initialise
total_length = length(ts);
count = zeros(total_length,K,nbNodes); %number of customers in the queue, start from time 0
count(1,:,1) = N; %initialise delay center with N jobs

% serial
for i = 1:total_length-1
    count(i+1,:,:) = count(i,:,:);
    count(i+1,:,:) = count(i,:,:);
    
    count(i+1,class_id(i),logger_id(i)) = count(i,class_id(i),logger_id(i))-1;
    
    if logger_id(i) == nbNodes
        count(i+1,class_id(i),1) = count(i,class_id(i),1)+1;
    else
        count(i+1,class_id(i),logger_id(i)+1) = count(i,class_id(i),logger_id(i)+1)+1;
    end
end

if sum(N) == 0
    for i = 1:K
        N(i) = max(max(count(:,i,:)));
    end
end

for i = 1:total_length
    for j = 1:K
        count(i,j,1) = count(i,j,1) + N(j);
    end
end


% parallel
% for i = 1:total_length-1
%     count(i+1,:,:) = count(i,:,:);
%
%     if logger_id(i) < 5
%         count(i+1,class_id(i),1) = count(i,class_id(i),1)-1;
%         count(i+1,class_id(i),logger_id(i)+1) = count(i,class_id(i),logger_id(i)+1)+1;
%     end
%
%     if logger_id(i) > 5
%         count(i+1,class_id(i),1) = count(i,class_id(i),1)+1;
%         count(i+1,class_id(i),logger_id(i)-9) = count(i,class_id(i),logger_id(i)-9)-1;
%     end
%
% end

count = reshape(count,total_length,nbClasses*nbNodes);

%calculate the interval between each timestamp
%time_interval(1) = ts(1);
time_interval(1) = 0;
time_interval(2:total_length) = diff(ts);

count(:,end+1) = time_interval';
count = count(burnin:end,:);

count = sortrows(count,[1:size(count,2)-1]);

[C ia ic] = unique(count(:,1:end-1),'rows','legacy');

%the first one
Time = C;
Time(1,end+1) = sum(count(1:ia(1),end));
for i =2:size(C,1)
    Time(i,end) = sum(count(ia(i-1)+1:ia(i),end));
end

%observed time period
obs_length = ts(end)-ts(burnin);
%calculate the probability
prob_nbCustomer = Time;
prob_nbCustomer(:,end) = prob_nbCustomer(:,end)/obs_length;

for i = 1:K
    N0(i) = sum(prob_nbCustomer(:,end).*prob_nbCustomer(:,K+i));
end
end

function [value, logG_current, range_size_dim] = gibbsSamplerSimple(alg,think_time,theta,testset,index,nbNodes,nbClasses,nbJobs,logG_initial,interval,range_size,LV,sumA)

range = (0:interval:range_size);
N = length(range);

if strcmp(alg,'TE')
    x = zeros(nbNodes,nbClasses);
    x(1,:) = think_time;
    for i = 1:nbNodes-1
        x(i+1,:) = theta(i*nbClasses+1-nbClasses:i*nbClasses);
    end
    
    index_i = floor((index-1)/nbClasses)+2;
    index_j = index-(index_i-2)*nbClasses;
    
    logG = zeros(1,N);
    
    [~,QN]=pfqn_bs(x(2:end,:),nbJobs,x(1,:));
    
    index_previous = find(range==theta(index));
    logG(index_previous) = logG_initial;
    
    for i = index_previous-1:-1:1
        x(index_i,index_j) = range(i+1);
        %[~,QN]=aql(x(2:end,:),nbJobs,x(1,:),interval);
        [~,QN]=pfqn_bs(x(2:end,:),nbJobs,x(1,:),interval,1000,QN);
        if 1+QN(index_i-1,index_j)/(range(i+1)+eps)*-interval < 0
            logG(i) = logG(i+1);
        else
            logG(i) = logG(i+1) + log(1+QN(index_i-1,index_j)/(range(i+1)+eps)*-interval);
        end
    end
    
    x(index_i,index_j) = theta(index);
    [~,QN]=pfqn_bs(x(2:end,:),nbJobs,x(1,:));
    
    for i = index_previous+1:N
        x(index_i,index_j) = range(i-1);
        %[~,QN]=aql(x(2:end,:),nbJobs,x(1,:),interval);
        [~,QN]=pfqn_bs(x(2:end,:),nbJobs,x(1,:),interval,1000,QN);
        if 1+QN(index_i-1,index_j)/(range(i-1)+eps)*interval < 0
            logG(i) = logG(i-1);
        else
            logG(i) = logG(i-1) + log(1+QN(index_i-1,index_j)/(range(i-1)+eps)*interval);
        end
    end
end

log_prob = zeros(1,N);
for i = 1:N
    theta(index) = range(i);
    if strcmp(alg,'TE')
        log_prob(i) = sum(testset(:,index+nbClasses))*log(range(i))-logG(i)*size(testset,1);
        %log_prob(i) = pdf_slice(alg,method,interval,tol,theta,nbJobs,think_time,testset,nbClasses,nbNodes,LV,sumA,logG(i));
    end
    if strcmp(alg,'MCI')
        log_prob(i) = pdf_slice(alg,interval,theta,nbJobs,think_time,testset,nbClasses,nbNodes,LV,sumA);
    end
end
log_prob = log_prob-max(log_prob);

prob = exp(log_prob);
prob = prob/sum(prob);

cum_prob = cumsum(prob);
rand_variable = rand(1);
index_prob = find(rand_variable<cum_prob);

%    if strcmp(alg,'TE')
range_size_dim = find(cum_prob > 1-1e-10, 1, 'first');
range_size_dim = range(range_size_dim)*2;
%     end
%     if strcmp(alg,'MCI')
%         range_size_dim = 0;
%     end

if isempty(index_prob)
    value = theta(index);
    if strcmp(alg,'TE')
        logG_current = logG_initial;
    end
else
    value = range(index_prob(1));
    if strcmp(alg,'TE')
        logG_current = logG(index_prob(1));
    end
end

if strcmp(alg,'MCI')
    logG_current = 0;
end
end

function [result] = pdf_slice(alg,interval,theta,nbJobs,think_time,testset,nbClasses,nbNodes,LV,sumA,logG)

if sum(theta < 0) > 0
    result = -inf;
    return;
end

theta = reshape(theta,nbClasses,nbNodes-1)';

if strcmp(alg,'TE') && exist('logG','var') == 0
    logG = approLogG([think_time;theta],nbJobs,interval);
end
%logG = log(gmva(theta,nbJobs,think_time));

result = 0;
n_node = zeros(size(testset,1),nbNodes);
for j = 2:nbNodes
    n_node(:,j) = sum(testset(:, nbClasses*j-nbClasses+1:nbClasses*j),2);
end
B=feval(@(x) LV(x+1), n_node(:,2:nbNodes));
result = result + sum(B(:));

y = [log(think_time+eps);log(theta+eps)];
temp = reshape(y',nbClasses*nbNodes,1);

result = result + sum(testset*temp);

result = result - logG*size(testset,1);

result = result - sumA;
end
