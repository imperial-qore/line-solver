function [ percentileRTs ] = returnRT2( arrival, service, pers, C, T_mode )

% the RTs for K = 2

service_h = build_Service_h(service);

[ T, S, A_jump, S_Arr, sum_Ajump ] = computeT( arrival, service, service_h, C, T_mode);
mWait = size(A_jump,1);
phi = sum(T-S_Arr,2);
[ pi0, En1 ] = computePi(T, arrival, service, service_h, C, S, A_jump);
[ wait_alpha, wait_Smat, prob_wait, alfa ] = returnWait( En1, pi0, T, phi, sum_Ajump );

% Waiting time percentiles
percentileWait = returnPer(wait_alpha, wait_Smat, pers);

%% Response time

[ ST, dim, dim_notbusy ] = generateService(service, service_h, C, S);
ST = kron(ST,arrival.Ia);

% starting state of the service process
pi0 = pi0/sum(pi0);
dim = arrival.ma*dim;

dim_notbusy = arrival.ma*dim_notbusy;
% the starting state according to arrival a the not-all-busy period
dim_service = dim + dim_notbusy; % except the empty state
% the starting state of service for a job in a not-all-busy period
notbusy_start = zeros(1,dim_service);
notbusy_start(1:dim) = notbusy_start(1:dim)+(1-prob_wait)*pi0;
% the starting state of service for a job in all-busy
busy_start = zeros(1,dim_service);

TS = T-S_Arr;
busy_start(1:dim) = busy_start(1:dim)+prob_wait*alfa*TS/sum(alfa*TS);

%% PH of response time
Tr = size(ST,1);
Sc = size(wait_Smat,2);
TS_sum = sum(TS,2)*ones(1,arrival.ma*mWait);
TS = TS./TS_sum;

stat_service_phase = busy_start/(-ST);
busy_nz = stat_service_phase>0;
tr_start_state = -sum(ST,2)'*diag(stat_service_phase);
tr_ST = zeros(size(ST));
tr_ST(busy_nz,busy_nz) = diag(1./stat_service_phase(busy_nz))...
    *ST(busy_nz,busy_nz)'*diag(stat_service_phase(busy_nz));

tr_ST_exit = sum(-tr_ST,2);
tr_ST_exit_mat = tr_ST_exit(1:Sc)*ones(1,Sc);


TS2 = (TS')*diag(alfa);
tr_TS_ind2 = sum(TS2,2); 
tr_TS_ind2(abs(tr_TS_ind2) < 10E-12) = ones(sum(abs(tr_TS_ind2) < 10E-12),1);
tr_TS_ind2 = TS2./( tr_TS_ind2*ones(1,arrival.ma*mWait) );

tildeP = [tr_ST_exit_mat.*tr_TS_ind2;zeros(Tr-Sc,Sc)];

% removal of unreachable phases in tr_ST
m_tr_ST = sum(busy_nz);
tr_ST = tr_ST(busy_nz,busy_nz);
tr_start_state = tr_start_state(busy_nz);
tildeP = tildeP(busy_nz,:);

gamma_res = [notbusy_start, tr_start_state zeros(1,Sc)];
C_res = [ST,            zeros(Tr,m_tr_ST),  zeros(Tr,Sc);
    zeros(m_tr_ST,Tr),  full(tr_ST),        tildeP;
    zeros(Sc,Tr),       zeros(Sc,m_tr_ST),  wait_Smat];

% Response time percentiles
percentileRTs = returnPer( gamma_res, C_res, pers );

