function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,METHOD] = SOLVER_NC_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;

line_debug('NC analyzer starting: method=%s, nstations=%d, nclasses=%d, njobs=%s', options.method, sn.nstations, sn.nclasses, mat2str(sn.njobs));

nservers = sn.nservers;
if max(nservers(nservers<Inf))>1 & any(isinf(sn.njobs)) & strcmpi(options.method,'exact') %#ok<AND2>
    line_error(mfilename,'NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the ''exact'' option.');
end

% interpolate for non-integer closed populations
eta = abs(sn.njobs - floor(sn.njobs));
if any(eta>GlobalConstants.FineTol)
    line_debug('Fractional populations detected, using interpolation');
    sn_floor = sn; sn_floor.njobs = floor(sn.njobs);
    [Qf,Uf,Rf,Tf,Cf,Xf,lGf,~,iterf] = solver_nc(sn_floor, options);
    sn_ceil = sn; sn_ceil.njobs = ceil(sn.njobs);
    [Qc,Uc,Rc,Tc,Cc,Xc,lGc,~,iterc,method] = solver_nc(sn_ceil, options);
    Q = Qf +  eta .* (Qc-Qf);
    U = Uf +  eta .* (Uc-Uf);
    R = Rf +  eta .* (Rc-Rf);
    T = Tf +  eta .* (Tc-Tf);
    C = Cf +  eta .* (Cc-Cf);
    X = Xf +  eta .* (Xc-Xf);
    lG = lGf +  eta .* (lGf-lGc);
    iter = iterc + iterf;
    line_debug('NC interpolation complete: used %d + %d iterations', iterf, iterc);
else % if integers or open model
    if any(isinf(sn.njobs))
        line_debug('Open/mixed model detected, calling solver_nc');
    else
        line_debug('Using exact integer populations, calling solver_nc');
    end
    [Q,U,R,T,C,X,lG,~,iter,method] = solver_nc(sn, options);
end
runtime = toc(Tstart);
end