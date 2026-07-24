%{ @file cache_build_item_graphs.m
 %  @brief Per-item access graph aggregated over users; {} when linear/absent.
 %
 %  @author LINE Development Team
%}

%{
 % @brief Build per-item (h+1)x(h+1) access graphs from ACCOST.
 %
 % @details
 % ACCOST is the per-(user,item) access graph cell{v,k} (each an (h+1)x(h+1)
 % matrix). Returns a 1xN cell of per-item graphs aggregated over users by
 % request rate, or {} when ACCOST is absent or the standard linear chain (so
 % the caller keeps the refined/linear path). Row 1 is miss admission, row 1+i
 % a hit in list i.
%}
function G = cache_build_item_graphs(accost, lambda, n, h)
G = {};
if isempty(accost)
    return;
end
lin = zeros(h+1, h+1); lin(1, 2) = 1;
for a = 1:(h-1), lin(a+1, a+2) = 1; end
lin(h+1, h+1) = 1;
u = size(accost, 1);
Gc = cell(1, n);
isLinear = true;
for k = 1:n
    num = zeros(h+1, h+1); den = 0;
    for v = 1:u
        wv = sum(lambda(v, k, 1));
        if ~isfinite(wv), wv = 0; end
        gvk = accost{v, k};
        if isempty(gvk), continue; end
        num = num + wv * gvk;
        den = den + wv;
    end
    if den > 0
        gk = num / den;
    else
        gk = accost{1, k};
        if isempty(gk), gk = lin; end
    end
    for a = 1:(h+1)
        srow = sum(gk(a, :));
        if srow > 0, gk(a, :) = gk(a, :) / srow; end
    end
    Gc{k} = gk;
    if ~all(all(abs(gk - lin) < 1e-9))
        isLinear = false;
    end
end
if ~isLinear
    G = Gc;
end
end
