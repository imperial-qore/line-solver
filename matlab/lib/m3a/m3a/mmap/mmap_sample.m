function [T,A,LAST,FIRST,STS] = mmap_sample(MMAP, nSamples, pi, seed)
% Generates samples from a Marked MAP.
% Input:
%  - MMAP: representation of the Marked MAP
%  - nSamples: number of random samples to collect
%  - pi: (optional) initial probability. If not provided, steady state
%        probabilities are used.
%  - seed: (optional) seed for the random number generator

% default: stationary inialization (as in single-class MAP)
if (nargin < 3)
   pi = map_pie(MMAP); 
end
LAST=[];
FIRST=[];

% check if single-class
% MAP has two matrices
% MMAP with K=1 has three matrices
if (length(MMAP) == 2 || length(MMAP) == 3)
    if (nargin == 4)
        [T,LAST,FIRST]=map_sample({MMAP{1}, MMAP{2}}, nSamples, pi, seed);
    else
        [T,LAST,FIRST]=map_sample({MMAP{1}, MMAP{2}}, nSamples, pi);
    end
    A = ones(nSamples,1);
    return
end

% if provided, initialize RNG with the given seed
if (nargin == 4)
    s = RandStream('mt19937ar', 'seed', seed);
    RandStream.setGlobalStream(s);
end

% number of states
M = size(MMAP{1},1);

% number of classes
C = length(MMAP)-2;

% single state -> exponential
if M==1
    T=exprnd(map_mean(MMAP),nSamples,1);
    crate = zeros(1,C);
    for c = 1:C
       crate(c) = MMAP{c+2}; 
    end
    A=randp(crate,nSamples,1);
    LAST=1;
    FIRST=1;
    return
end

% rates of the exponential distributions in each state
D0 = MMAP{1};
rates=-diag(D0);
% row s contains the discrete probability distribution for the transitions
% from state s (the first M are the embedded transitions, the following 
% M*C are the transitions that mark an arrival)
targets=zeros(M,M+M*C);
for i = 1:M
   for j = 1:M
       if i ~= j
           targets(i,j) = D0(i,j) / rates(i);
       end
   end
   for c = 1:C
       for j = 1:M
           targets(i,c*M+j) = MMAP{2+c}(i,j) / rates(i);
       end
   end
end

% pick a random initial state according to pi
initState=find(rand()<cumsum(pi),1);

T = zeros(nSamples,1);
A = zeros(nSamples,1);

% simulate
s=initState;
tTime=0;
h = 1;
batchSize=zeros(1,M);
batchCounter=zeros(1,M);
timeBatch=cell(1,M);
tranBatch=cell(1,M);
for p = 1:M
    batchSize(p)=0;
    batchCounter(p)=1;
    timeBatch{p}=[];
    tranBatch{p}=[];
end
STS = zeros(nSamples,1);
while h <= nSamples
    STS(h) = s;
    % regenerate exponential batch if needed
    if (h == 1 || batchCounter(s) > batchSize(s))
        batchSize(s)=min(100,nSamples-h+1);
        timeBatch{s}=exprnd(1/rates(s),batchSize(s),1);
        tranBatch{s}=randp(targets(s,:),batchSize(s),1);
        batchCounter(s)=1;
    end
    % pick from batch random values
    time = timeBatch{s}(batchCounter(s));
    transition = tranBatch{s}(batchCounter(s));
    % increase batch counter usage
    batchCounter(s)=batchCounter(s)+1;
    % choose transition
    % should be equivalent to but much faster than idivide
    % type=idivide(int32(transition-1),int32(M),'floor');
    type=int32(floor((transition-1)/M));
    if (type) == 0
        % no event, internal transition
        tTime=tTime + time;
        sNext=mod(transition,M);
        if (sNext == 0)
            sNext = M;
        end
        %fprintf('No event at time %0.2f, transition %d->%d\n', ...
        %    tTime, s, sNext);
        % change state
        s=sNext;
    else
        % event
        c=type;
        T(h)=tTime + time;
        A(h)=c;
        sNext = mod(transition,M);
        if (sNext == 0)
            sNext = M;
        end
        %fprintf('Event %2d at time %.2f, transition %d->%d, class %d\n', ...
        %    h, T(h), s, sNext, A(h));
        % reset time accumulator
        tTime=0;
        % change state
        s=sNext;
        % next sample
        h=h+1;
    end    
end

end