function [ percentileRTs ] = returnRT1( arrival, service, pers )

% The following function needs to use the QMAM toolbox, which can be 
% downloaded from http://win.ua.ac.be/~vanhoudt/.

% the RTs for K = 1

[ ql, res_alpha, wait_alpha, Smat ] = Q_CT_MAP_MAP_1( arrival.lambda0, arrival.lambda1,service.ST, service.St*service.tau_st );

percentileRTs = returnPer( res_alpha, Smat, pers );
