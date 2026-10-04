function ILchain = sn_interlock_chain(sn, ILclass)
% ILCHAIN = SN_INTERLOCK_CHAIN(SN, ILCLASS)
%
% Aggregate a class-indexed interlock matrix to the chain basis the MVA
% solvers work in. ILCLASS(r,s) is the share of the class-s queue that a
% class-r arrival must not see, the interlocked flow of Franks (1999),
% Eq. (4.7). Two classes of the same chain belong to the same client, so the
% diagonal blocks carry no information and the chain diagonal is zero: an
% arrival always sees its own chain in full.
%
% Reference: G. Franks, "Performance Analysis of Distributed Server Systems",
% PhD thesis, Carleton University, 1999, Ch. 4.

ILchain = [];
if isempty(ILclass)
    return
end
R = sn.nclasses;
if size(ILclass,1) ~= R || size(ILclass,2) ~= R
    line_error(mfilename, sprintf('the interlock matrix is %dx%d but the model has %d classes.', size(ILclass,1), size(ILclass,2), R));
end

K = sn.nchains;
ILchain = zeros(K,K);
for cr = 1:K
    memr = find(sn.chains(cr,:));
    if isempty(memr)
        continue
    end
    for cs = 1:K
        if cs == cr
            continue
        end
        mems = find(sn.chains(cs,:));
        if isempty(mems)
            continue
        end
        ILchain(cr,cs) = max(max(ILclass(memr,mems)));
    end
end

if ~any(ILchain(:) > 0)
    ILchain = [];
end
end
