function ql = getQLArrival(data)
% GETQLARRIVAL determines the queue lengths observed at arrival times
% when the arrival and response times are available.  
% Assumes data is available in standard format. 

precision = 1E-6; %0.1 ms
K = size(data,2) - 1;
at = data{3,1}/1000; % to secs 
rt = data{4,1};
classes = ones(size(data{4,1},1) ,1);
I = eye(K);
delta = ones(size(data{4,1},1),1)*I(1,:);
for j = 2:K
    at = [at; data{3,j}/1000];  % to secs 
    rt = [rt; data{4,j}];
    classes = [classes; j*ones( size(data{4,j},1),1 ) ]; 
    delta =   [delta;     ones( size(data{4,j},1),1 )*I(j,:) ];
end
alltimes = [at rt classes delta]; %arrival time - response time - classes - delta
alltimes = sortrows(alltimes,1);
exittimes = sum(alltimes(:,1:2),2); 
n = size(alltimes,1);
jobIDs = [1:n]'; 


arrexit = [jobIDs alltimes(:,[1 3]) ones(n,1) alltimes(:,4:end)];         % jobID - arrival time - class - +1 - +delta
% add epsilon to differentiate from departure times (first departures, then arrivals)
arrexit(:,2) = arrexit(:,2) + precision;
arrexit = [arrexit;   
          [jobIDs exittimes alltimes(:,3) -ones(n,1) -alltimes(:,4:end)] ];% jobID - exit time    - class - -1 - -delta


arrexit= sortrows(arrexit,2); 

state = zeros(1,K);
stateObs = zeros(n,K); %state observed upon arrival
for i = 1:n*2
    state = state+arrexit(i,5:end);
    if arrexit(i,4) == 1
        stateObs(arrexit(i,1),:) = state; % - I(arrexit(i,3),:);
    end
end


origOrder = cell(1,K);
numObs= zeros(1,K);
for j = 1:K
    numObs(1,j) = size(data{3,j},1);
    jobOrder = [1:numObs(1,j)]'; 
    arr = [jobOrder data{3,j}];
    arr = sortrows(arr,2);
    origOrder{j} = arr(:,1);
end
ql = cell(1,K);
arrexit(:,4) = arrexit(:,4)*(-1); %invert sign for later ordering
arrexit = sortrows(arrexit, [4 3 2]);
counter = 0;
for j = 1:K
    jobIDsClass = arrexit(counter+1:counter+numObs(j), 1);
    ql{1,j} = stateObs(jobIDsClass,:);
    ql{1,j}(origOrder{j},:) = ql{1,j};
    counter = counter + numObs(j); 
end
end