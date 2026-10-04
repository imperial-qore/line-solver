function [S_notallbusy] = constructNotAllBusy( C, services, service_h )

dim = length(service_h.beta);
m = length(services.tau_st);

indexes_notbusy = build_index(m,1);
dim_NB = size(indexes_notbusy,1);
dim_C = C+1;
dim_notbusy = (dim_C-1)*dim_NB+1;

S_notallbusy = zeros(dim_notbusy,dim_notbusy);

%% from not-busy to not-busy

for row = 1 : dim_C-1
    S_notallbusy((row-1)*dim_NB+1:row*dim_NB,(row-1)*dim_NB+1:row*dim_NB) = services.ST;
end

A = -sum(services.ST,2)*services.tau_st;

for row = 1 : dim_C-2
    S_notallbusy((row-1)*dim_NB+1:row*dim_NB,row*dim_NB+1:(row+1)*dim_NB) = A;
end

S_notallbusy((dim_C-2)*dim_NB+1:(dim_C-1)*dim_NB,(dim_C-1)*dim_NB+1:end) = -sum(services.ST,2);

for row = 1 : dim_notbusy
    S_notallbusy(row,row) = 0;
    S_notallbusy(row,row) = -sum(S_notallbusy(row,:));
end

end