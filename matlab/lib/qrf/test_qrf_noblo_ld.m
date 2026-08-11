%% QRF BAS test from code release
%function [UUB,ULB,ID]=runquadraticfcrbas(MAP,N,P,F,MM)
%addpath C:\Users\csg\Dropbox\experiments\09_mapqnfc\
P=[0,1;1,0];
MAPs={map_hyperexp(2,20); map_exponential(1)};
N = 5; 
M = 2;
alpha = ones(M,N);
[XN,QN,UN]=mapqn_ezsolve(MAPs,N,P,alpha);
UN, QN
%%
MR=1;
[UNqrf,QNqrf] = qrf_noblo_mmi_ld(MAPs,N,P,alpha);
[sum(UN,2)';UNqrf]
QNqrf=QNqrf/sum(QNqrf)*N;
[sum(QN,2)';QNqrf]