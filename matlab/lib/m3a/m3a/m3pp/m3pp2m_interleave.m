function s = m3pp2m_interleave(m3pps)
% Computes the interleaved MMAP obtained by multiple M3PP(2,m).

L = length(m3pps);

symbolic = 0;
for i = 1:L
    if map_issym(m3pps{i})
        symbolic = 1;
        break;
    end
end

if symbolic
    r = sym(zeros(2,L));
else
    r = zeros(2,L);
end
r(1,L) = m3pps{L}{1}(1,2);
for i = (L-1):-1:1
    r(1,i) = m3pps{i}{1}(1,2) - sum(r(1,(i+1):end));
end
r(2,1) = m3pps{1}{1}(2,1);
for i = 2:L
    r(2,i) = m3pps{i}{1}(2,1) - sum(r(2,1:(i-1)));
end

M = 0;
for i = 1:L
    M = M + size(m3pps{i},2) - 2;
end

n = 2+(L-1);
s = cell(1,2+M);
for i = 1:n
    for j = 1:n
        if j > i
            s{1}(i,j) = r(1,j-1);
        elseif j < i
            s{1}(i,j) = r(2,j);
        end
    end
end
% compute D1c, c = (i-1)*2+j
c = 1;
for i = 1:L % for each M3PP[2]
    m = size(m3pps{i},2)-2;
    for j = 1:m % for each class in the i-th M3PP[m]
        if symbolic
            s{2+c} = sym(zeros(n,n));
        else
            s{2+c} = zeros(n,n);
        end
        for h = 1:n
            % set transition rates
            if h <= i
                s{2+c}(h,h) = m3pps{i}{2+j}(1,1);
            else
                s{2+c}(h,h) = m3pps{i}{2+j}(2,2);
            end
        end
        c = c + 1;
    end
end
% compute D1
if symbolic
    s{2} = sym(zeros(n,n));
else
    s{2} = zeros(n,n);
end
for i = 1:M
    s{2} = s{2} + s{2+i};
end
% compute diagonal of D0
for h = 1:n
    s{1}(h,h) = -sum(s{1}(h,:)+s{2}(h,:));
end

end