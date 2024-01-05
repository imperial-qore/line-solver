from numpy import *
from ..util import *

def pfqn_stdf(L = None,N = None,Z = None,S = None,fcfsNodes = None,rates = None,tset = None):
    # sojourn time distribution function at multiserver FCFS nodes
    # J. McKenna 1987 JACM
    stabilityWarnIssued = False
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    T = len(tset)
    mu = zeros((M,sum(N)))
    for k in arange(0,M+1):
        for n in arange(0,sum(N)+1):
            mu[k,n] = amin(S[k],n)

    # when t==0 the hkc function is not well defined
    tset[tset == 0] = Distrib.Zero
    for k in transpose(fcfsNodes):
        if range_(rates[k,:]) > Distrib.Zero:
        line_error(mfilename,'The FCFS stations has distinct service rates, the model is invalid.')
    hMAP = cell(1,1 + sum(N))
    for n in arange(0,sum(N)+1):
        if n < S[k]:
            hMAP[1 + n] = map_exponential(1 / rates(k,1))
        else:
            hMAP[1 + n] = map_sumind(array([map_exponential(1 / rates(k,1)),map_erlang((n - S[k] + 1) / (S[k] * rates(k,1)),n - S[k] + 1)]))
        hkc[arange[1,T+1],1 + n] = transpose(map_cdf(hMAP[1 + n],tset))
    for r in arange(0,R+1):
        if L[k,r] > Distrib.Zero:
            Nr = oner(N,r)
            __,__,__,__,lGr,isNumStable = pfqn_mvald(L,Nr,Z,mu)
            lGr = lGr(end())
            if not isNumStable  and not stabilityWarnIssued :
                stabilityWarnIssued = True
                line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable')
            RD[k,r] = zeros((len(tset),2))
            RD[k,r][:,2] = tset
            ## this is the original code in the paper
            #             Gkrt = zeros(T,1);
            #             nvec = pprod(Nr);
            #             while nvec >= 0
            #                 [~,~,~,~,lGk,isNumStable]  = pfqn_mvald(L(setdiff1d(1:M,k),:),Nr-nvec,Z,mu(setdiff1d(1:M,k),:));
            #                 lGk = lGk(end);
            #                 if ~isNumStable && ~stabilityWarnIssued
            #                     stabilityWarnIssued = true;
            #                     line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable');
            #                 end
            #                 for t=1:T
            #                     if sum(nvec)==0
            #                         Gkrt(t) = Gkrt(t) + hkc(t,1+0) * exp(lGk);
            #                     else
            #                         lFk = nvec*log(L[k,:])' +factln(sum(nvec)) - sum(factln(nvec)) - sum(log(mu(k,1:sum(nvec))));
            #                         Gkrt(t) = Gkrt(t) + hkc(t,1+sum(nvec)) * exp(lFk + lGk);
            #                     end
            #                 end # nvec
            #                 nvec = pprod(nvec, Nr);
            #             end
            #             lGkrt = log(Gkrt);
            #             RD{k,r}(1:T,1) = exp(lGkrt - lGr);
            ## this is faster as it uses the recursive form for LD models
            Hkrt = []
            __,__,__,__,lGk = pfqn_mvald(L(setdiff1d(arange(0,M+1),k),:),Nr,Z,mu(setdiff1d(arange(0,M+1),k),:))
            lGk = lGk(end())
            for t in arange(0,T+1):
                #[t,T]
                gammat = mu
                for m in arange(0,sum(Nr)+1):
                    gammat[k,m] = mu(k,m) * hkc(t,1 + m - 1) / hkc(t,1 + m)
                gammak = mushift(gammat,k)
                gammak = gammak(:,arange(0,(sum(Nr) - 1)+1))
                Hkrt[t] = hkc(t,1 + 0) * exp(lGk)
                for s in arange(0,R+1):
                    if Nr[s] > 0:
                        #gammak
                        #lYks_t  = pfqn_rd(L,oner(Nr,s),Z,gammak);
                        #RD = lYks_t
                        #NRP
                        #lYks_t = pfqn_nrp(L,oner(Nr,s),Z,gammak);
                        __,__,__,__,lYks_t = pfqn_mvald(L,oner(Nr,s),Z,gammak)
                        lYks_t = lYks_t(end())
                        #if t==10
                        #    keyboard
                        #end
                        Hkrt[t] = Hkrt(t) + (L(k,s) * hkc(t,1 + 0) / gammat(k,1)) * exp(lYks_t)
            Hkrt[isnan[Hkrt]] = Distrib.Zero
            lHkrt = log(Hkrt)
            RD[k,r][arange[1,T+1],1] = amin(1,exp(lHkrt - lGr))

return RD


def Fm(m = None,x = None):
    if m == 1:
        Fmx = 1 - exp(- x)
    else:
        A = 0
        for j in arange(0,(m - 1)+1):
            A = A + x ** j / factorial[j]
        Fmx = 1 - exp(- x) * A

    return Fmx