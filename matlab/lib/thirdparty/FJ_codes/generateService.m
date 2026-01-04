function [ T, newdim, dim_notbusy ] = generateService( services, service_h, C, S )

% C: limit in the difference between the shortest and longest queue lengths;
%
% For example: for the case K = 2 servers, C = 2, phases are ordered as
% (2,(1,0,1)), (2,(1,1,0)), (2,(2,0,0)), (1,(0,0,1)) (1,(0,1,0))
% (1,(1,0,0)), S

dim = length(service_h.beta);
m = length(services.tau_st);

indexes_notbusy = build_index(m,1);
dim_NB = size(indexes_notbusy,1);

dim_C = C+1;

newdim = dim_C*dim;
dim_notbusy = dim_C*dim_NB;

T = zeros(newdim+dim_notbusy,newdim+dim_notbusy);
t = zeros(newdim+dim_notbusy,1);

t(end-dim_NB+1:end) = services.St;

T(1:newdim,1:newdim) = S;

%% S_long: the job in the shorter queue completes service
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

for row = 1 : dim_C-1
    T((row-1)*dim+1:row*dim,newdim+(row-1)*dim_NB+1:newdim+row*dim_NB) = S_long;
end

%% S_last: all queues have the same length
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

T((dim_C-1)*dim+1:dim_C*dim,newdim+(dim_C-1)*dim_NB+1:newdim+dim_C*dim_NB) = S_last;

%% only one subtask in service

for row = 1 : dim_C
    T(newdim+(row-1)*dim_NB+1:newdim+row*dim_NB,newdim+(row-1)*dim_NB+1:newdim+row*dim_NB) = services.ST;
end

A = -sum(services.ST,2)*services.tau_st;

for row = 1 : dim_C-1
    T(newdim+(row-1)*dim_NB+1:newdim+row*dim_NB,newdim+row*dim_NB+1:newdim+(row+1)*dim_NB) = A;
end

for row = newdim+1 : newdim+dim_notbusy
    T(row,row) = 0;
    T(row,row) = -sum(T(row,:))-t(row);
end

