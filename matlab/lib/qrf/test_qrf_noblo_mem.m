%% QRF BAS test from code release
%function [UUB,ULB,ID]=runquadraticfcrbas(MAP,N,P,F,MM)
addpath C:\Users\csg\Dropbox\experiments\09_mapqnfc\
rt=[0,1; 1,0];
MAPs={map_exponential(0.1); map_exponential(10)};
M = length(MAPs);
N=3; %(F1+1)
%F = [N-1;N];
%[XN,QN,UN,EN,p,p1,p1c,p2,MAPQN,MM]=mapqnfc_ezsolve_bas(MAPs,N,P,ones(size(MAPs,1),N),F)
[XNctmc,QNctmc,UNctmc]=mapqn_ezsolve(MAPs,N,P,ones(size(MAPs,1),N))

%[UNqrf,QNqrf] = ctmc_mmi_noblo(M,MR,K,N,mu,v,rt);
[UNqrf,QNqrf] = qrf_noblo_mem(M,MR,K,N,mu,v,rt);
UNex=sum(UNctmc,2)'
UNap=UNqrf
QNex=sum(QNctmc,2)'
QNap=QNqrf/sum(QNqrf)*N