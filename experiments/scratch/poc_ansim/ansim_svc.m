function d = ansim_svc(meanval, law)
% ANSIM_SVC  Service law at the priority queue, at a fixed mean MEANVAL.
%   'exp'      SCV = 1     Exp
%   'hyperexp' SCV = 4     two-phase HyperExp
%   'hypoexp'  SCV = 0.5   Erlang-2
switch law
    case 'exp'
        d = Exp(1/meanval);
    case 'hyperexp'
        d = HyperExp.fitMeanAndSCV(meanval, 4.0);
    case 'hypoexp'
        d = Erlang.fitMeanAndSCV(meanval, 0.5);
    otherwise
        error('unknown service law %s', law);
end
end
