function [Se, Sestar, R0, Ke, Kc, S_long, S_last] = constructSRK( C, services, service_h, S)

dim = length(service_h.beta);
m = length(services.tau_st);

indexes_notbusy = build_index(m,1);
dim_NB = size(indexes_notbusy,1);

dim_C = C+1;

newdim = dim_C*dim;
dim_notbusy = (dim_C-1)*dim_NB+1;

Se = zeros(newdim+dim_notbusy,newdim+dim_notbusy);
Sestar = Se;

Se(1:newdim,1:newdim) = S;

%% from busy to not-busy
% when the job in the shorter queue complete service
S_long = zeros(dim,dim_NB);

for row = 1 : dim
    
    countvect = service_h.service_phases(row,:);
    
    for i = m+1:2*m % the starting phase
        if countvect(i) > 0
            
            tovect = countvect(1:m);
            
            col = vectmatch(tovect,indexes_notbusy);
            
            % the job in the shorter queue completes service
            S_long(row,col) = S_long(row,col)+countvect(i)*services.St(i-m);
            
        end
    end
end

Se(1:dim,newdim+1:newdim+dim_NB) = S_long;
for row = 2 : dim_C-1
    Se((row-1)*dim+1:row*dim,newdim+(row-2)*dim_NB+1:newdim+(row-1)*dim_NB) = S_long;
end

% when the 2 queue have the same length
S_last = zeros(dim,dim_NB);

for row = 1 : dim
    
    countvect = service_h.service_phases(row,:);
    
    for k = 1 : 2
        for i = (k-1)*m+1:k*m
            if countvect(i) > 0
                
                tovect = countvect((2-k)*m+1:(2-k+1)*m);
                
                col = vectmatch(tovect,indexes_notbusy);
                
                S_last(row,col) = S_last(row,col)+countvect(i)*services.St(i-(k-1)*m);
                
            end
        end
    end
end

Se((dim_C-1)*dim+1:dim_C*dim,newdim+(dim_C-2)*dim_NB+1:newdim+(dim_C-1)*dim_NB) = S_last;

Sestar(1:newdim,newdim+1:end) = Se(1:newdim,newdim+1:end);

%% from not-busy to not-busy

for row = 1 : dim_C-1
    Se(newdim+(row-1)*dim_NB+1:newdim+row*dim_NB,newdim+(row-1)*dim_NB+1:newdim+row*dim_NB) = services.ST;
end

A = -sum(services.ST,2)*services.tau_st;

for row = 1 : dim_C-2
    Se(newdim+(row-1)*dim_NB+1:newdim+row*dim_NB,newdim+row*dim_NB+1:newdim+(row+1)*dim_NB) = A;
end

Se(newdim+(dim_C-2)*dim_NB+1:newdim+(dim_C-1)*dim_NB,newdim+(dim_C-1)*dim_NB+1:end) = -sum(services.ST,2);

for row = newdim+1 : newdim+dim_notbusy
    Se(row,row) = 0;
    Se(row,row) = -sum(Se(row,:));
end


%% Construct R0 matrix: from the not-busy period to a busy period
% focus on the shortest queue, the system has only one idle server, thus an
% arrival leads the system back to a busy-period
% first dim phases: busy; second dim phases: not-busy
% the queue jumps to the same phase in a busy period after an arrival
R0 = zeros(newdim+dim_notbusy,newdim+dim_notbusy);
R_NB = zeros(dim_NB,dim);

for row = 1 : dim_NB
    countvect = indexes_notbusy(row,:);
    for i = 1 : m
        tovect_short = zeros(1,m);
        tovect_short(i) = 1;
        tovect = [countvect,tovect_short];
        col = vectmatch(tovect,service_h.service_phases);
        R_NB(row,col) = R_NB(row,col)+services.tau_st(i);
    end
end


for row = 1 : dim_C-1
    R0(newdim+(row-1)*dim_NB+1:newdim+row*dim_NB,(row-1)*dim+1:row*dim) = R_NB;
end

R0(end,newdim-dim+1:newdim) = service_h.beta;

Ke = cat(2,eye(newdim),zeros(newdim,dim_notbusy));
Kc = cat(1,eye(newdim),zeros(dim_notbusy,newdim));

Se = sparse(Se);
Sestar = sparse(Sestar);
R0 = sparse(R0);
Ke = sparse(Ke);
Kc = sparse(Kc);

end
