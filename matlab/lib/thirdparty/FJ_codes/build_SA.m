function [S,A_jump] = build_SA(services, service_h, C)

dim = length(service_h.beta);
m = length(services.tau_st);

dim_C = C+1;

newdim = dim_C*dim;
S = zeros(newdim,newdim);
A_jump = zeros(newdim,newdim);

% for example: 2 subtasks of different or same job, each can choose
% from 2 phases, then states are:
% (1,0;1,0); (1,0;0,1); (0,1;1,0); (0,1;0,1)
% since the two subtasks are not symmetric, the first is the for the longer
% one, the latter is for the shorter one
% the phase changes within each c in C is service_h.S
% no job completes service
for row = 1 : dim_C
    S((row-1)*dim+1:row*dim,(row-1)*dim+1:row*dim) = service_h.S;
end

% when a job completes service, and a new job starts service
A = -sum(services.ST,2)*services.tau_st;

% when the job in the longer queue complete service
S_Cminus1 = zeros(dim,dim);
for row = 1 : dim
    
    countvect = service_h.service_phases(row,:);
    
    for i = 1:m % the starting phase
        if countvect(i) > 0
            for j = 1:m % the resulting phase
                
                tovect=countvect;
                tovect(i) = tovect(i)-1; % jump from phase i to phase j
                tovect(j) = tovect(j)+1;
                
                col = vectmatch(tovect,service_h.service_phases);
                
                % the job in the longer queue completes service
                S_Cminus1(row,col) = S_Cminus1(row,col)+countvect(i)*A(i,j);
                
            end
        end
    end
end

% form c to c-1
for c = C : -1 : 1
    S((C-c)*dim+1:(C-c+1)*dim,(C-c+1)*dim+1:(C-c+2)*dim) = S_Cminus1;
end

% when the job in the shorter queue complete service
A_Cplus1 = zeros(dim,dim);
for row = 1 : dim
    
    countvect = service_h.service_phases(row,:);
    
    for i = m+1:2*m % the starting phase
        if countvect(i) > 0
            for j = m+1:2*m % the resulting phase
                
                tovect=countvect;
                tovect(i) = tovect(i)-1; % jump from phase i to phase j
                tovect(j) = tovect(j)+1;
                
                col = vectmatch(tovect,service_h.service_phases);
                
                % the job in the shorter queue completes service
                A_Cplus1(row,col) = A_Cplus1(row,col)+A(i-m,j-m);
                
            end
        end
    end
end

for c = C-1 : -1 : 1
    A_jump((C-c)*dim+1:(C-c+1)*dim,(C-c-1)*dim+1:(C-c)*dim) = A_Cplus1;
end

A_jump(1:dim,1:dim) = A_Cplus1;

% both queue have a same length
% the other one that didn't finish becomes the longer queue
A_last = zeros(dim,dim);
for row = 1 : dim
    
    countvect = service_h.service_phases(row,:);
    
    for k = 1 : 2
        for i = 1 : m % the starting phase
            if countvect((k-1)*m+i) > 0
                
                tovect = countvect((2-k)*m+1:(2-k+1)*m);
                for j = 1 : m
                    temp_tovector = zeros(1,m);
                    temp_tovector(j) = 1;
                    temp_tovector = [tovect,temp_tovector];
                    col = vectmatch(temp_tovector,service_h.service_phases);
                    
                    A_last(row,col) = A_last(row,col)+countvect((k-1)*m+i)*A(i,j);
                    
                end
            end
        end
    end
end

A_jump(C*dim+1:end,(C-1)*dim+1:C*dim) = A_last;
