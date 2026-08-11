function MMAP=mmap_mixture(alpha,MAPs)
% probabilistic composition of MAPs - MAP=map_pcompose(P,MAPs)
Dk = {};
I = length(MAPs);
for j=1:I+2
    Dk{j}= [];
end
for i = 1:I
    if isempty(MAPs{i})
        MAPs{i} = map_exponential(1e6);
    end
end
for i = 1:I
    
    if i==1
        Dk{1} = MAPs{i}{1};
    else
        Dk{1} = blkdiag(Dk{1},MAPs{i}{1});
    end
    D1i = alpha(1)*MAPs{i}{2}*e(length(MAPs{i}{2}))*map_pie(MAPs{1});
    for j=2:I
        D1i = horzcat(D1i,MAPs{i}{2}*alpha(j)*e(length(MAPs{i}{2}))*map_pie(MAPs{j}));
    end
    Dk{2} = vertcat(Dk{2},D1i);
    for j=1:I
        if i==j
            Dk{2+j} = vertcat(Dk{2+j},   D1i);
        else
            Dk{2+j} = vertcat(Dk{2+j}, 0*D1i);
        end
    end
end

MMAP = mmap_normalize(Dk);
end