function [ service_h ] = build_Service_h( service )

% This function builds the PH representations for a 2-node FJ job
% phase for one single subtask

dim_single = length(service.tau_st);
phases_single = build_index(dim_single,1);

% Possible service phases for the 2-node FJ job
service_h.service_phases = zeros(dim_single^2,2*size(phases_single,2));
k = 1;
for i = 1 : dim_single
    for j = 1 : dim_single
        % [longest, shortest]
        service_h.service_phases(k,:) = [phases_single(i,:),phases_single(j,:)];
        k = k+1;
    end
end
% PH representation for the service time of a 2-node FJ job 
service_h.beta = kron(service.tau_st,service.tau_st);
service_h.S = kronsum(service.ST,service.ST);

end

