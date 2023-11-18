function [infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(self,invertSymbol,primeNumbers)
if nargin<2
    invertSymbol = false;
end
if nargin<3
    primeNumbers = false;
end
if ~exist('sym')
    line_error(mfilename,'This method requires MATLAB''s Symbolic Toolbox.');
end
if ~isdeployed
    [~, F] = getGenerator(self);
    [ stateSpace, nodeStateSpace] = getStateSpace(self);
    infGen = sym(zeros(size(F{1})));
    eventFilt = cell(1, length(F));
    if primeNumbers % use prime numbers instead than symbolic variables
        N=10;
        xprime = primes(N);
        while length(xprime) < length(F)
            N = N*2;
            xprime = primes(N);
        end
    end
    for e = 1:length(F)
        F{e} = full(F{e});
        minF = min(min(F{e}(F{e}>0)));
        if ~isempty(minF)
            F{e} = F{e} / minF;
            if invertSymbol
                if ~primeNumbers
                    eventFilt{e} = F{e} / sym(['x',num2str(e)],'real');
                else
                    eventFilt{e} = F{e} / xprime(e);
                end
            else
                if ~primeNumbers
                    eventFilt{e} = F{e} * sym(['x',num2str(e)],'real');
                else
                    eventFilt{e} = F{e} * xprime(e);
                end
            end
            infGen = infGen + eventFilt{e};
        end
    end
    infGen = ctmc_makeinfgen(infGen);
    syncInfo = self.getStruct.sync;
else
    infGen = [];
    eventFilt = [];
    syncInfo = [];
end
end