%{ @file cache_pos_drift_graph.m
 %  @brief General position-resolved mean-field drift honouring a per-item
 %         cache access graph, for FIFO(m) and strict FIFO(m).
 %
 %  @author LINE Development Team
%}

%{
 % @brief Position-resolved DDPP drift with a per-item access graph.
 %
 % @details
 % G{k} is the (h+1)x(h+1) access graph of item k: row 1 is miss admission
 % (col 1 = reject, col 1+l = admit to list l), row 1+i is a hit in list i
 % (col 1+b = promote to list b>=i; b==i means STAY in place, the FIFO/SFIFO
 % convention). A miss admits at the head of the target list (its tail evicted);
 % a hit at position j of list i promotes to the head of target b>i (the tail of
 % b demoted to list i -- to the vacated position j for FIFO, to the head with a
 % 1..j-1 shift for SFIFO). REINSERT is 'head' (SFIFO) or 'pos' (FIFO). Reduces
 % exactly to the linear drift when G is the standard chain; requires a cold
 % (empty) initial state so non-admissible items drain.
%}
function dX = cache_pos_drift_graph(x, p, G, m, n, h, slots, sidx, S, reinsert)
x = max(0, min(1, x));
isHead = strcmp(reinsert, 'head');
Kidx = @(k, i, j) (k-1)*S + sidx(i, j);

MI = zeros(1, h);        % MI(l): miss admission into list l
HP = zeros(h, h);        % HP(i,b): promotion i->b (b>i)
for k = 1:n
    ok = 1 - pos_sumocc(x, k, S, size(slots,1));
    gk = G{k};
    for l = 1:h
        MI(l) = MI(l) + p(k) * ok * gk(1, l+1);
    end
    for i = 1:h
        oc = 0;
        for jj = 1:m(i), oc = oc + x(Kidx(k, i, jj)); end
        for b = (i+1):h
            HP(i, b) = HP(i, b) + p(k) * oc * gk(i+1, b+1);
        end
    end
end
Sin = zeros(1, h);
for l = 1:h
    Sin(l) = MI(l);
    for s = 1:(l-1), Sin(l) = Sin(l) + HP(s, l); end
end
POp = zeros(h, max(m));
for s = 1:S
    i = slots(s,1); j = slots(s,2);
    acc = 0;
    for k = 1:n
        acc = acc + p(k) * x(Kidx(k, i, j)) * (1 - G{k}(i+1, i+1));
    end
    POp(i, j) = acc;
end

dX = zeros(n*S, 1);
for k = 1:n
    gk = G{k};
    ok = 1 - pos_sumocc(x, k, S, size(slots,1));
    for s = 1:S
        i = slots(s,1); j = slots(s,2);
        xk = x(Kidx(k, i, j));
        o = p(k) * xk * (1 - gk(i+1, i+1));
        if isHead
            o = o + (Sin(i) + pos_gg(POp, i, j, m)) * xk;
        else
            o = o + Sin(i) * xk;
        end
        dX(Kidx(k, i, j)) = dX(Kidx(k, i, j)) - o;
        if j >= 2
            if isHead
                dX(Kidx(k,i,j)) = dX(Kidx(k,i,j)) + (Sin(i) + pos_gg(POp,i,j-1,m)) * x(Kidx(k,i,j-1));
            else
                dX(Kidx(k,i,j)) = dX(Kidx(k,i,j)) + Sin(i) * x(Kidx(k,i,j-1));
            end
        else
            dX(Kidx(k,i,1)) = dX(Kidx(k,i,1)) + p(k) * ok * gk(1, i+1);
            for ss = 1:(i-1)
                occ_s = 0;
                for jj = 1:m(ss), occ_s = occ_s + x(Kidx(k, ss, jj)); end
                dX(Kidx(k,i,1)) = dX(Kidx(k,i,1)) + p(k) * occ_s * gk(ss+1, i+1);
            end
        end
        for b = (i+1):h
            if isHead
                if j == 1
                    dX(Kidx(k,i,1)) = dX(Kidx(k,i,1)) + HP(i,b) * x(Kidx(k, b, m(b)));
                end
            else
                poj = 0;
                for kk = 1:n
                    poj = poj + p(kk) * x(Kidx(kk, i, j)) * G{kk}(i+1, b+1);
                end
                dX(Kidx(k,i,j)) = dX(Kidx(k,i,j)) + poj * x(Kidx(k, b, m(b)));
            end
        end
    end
end
end

function s = pos_sumocc(x, k, S, ns)
s = 0;
for t = 1:ns
    s = s + x((k-1)*S + t);
end
end

function g = pos_gg(POp, i, jp, m)
g = 0;
for jj = (jp+1):m(i)
    g = g + POp(i, jj);
end
end
