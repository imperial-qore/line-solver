function qbd_bmapbmap1(MAPa, pbatcha, MAPs)

na = length(MAPa{1});
ns = length(MAPs{1});
maxbatch = length(pbatcha);

for b=1:maxbatch
    A1{b} = kron(MAPa{2}*pbatcha(b),eye(ns));
end
A0 = krons(MAPa{1},MAPs{1});
A_1 = kron(eye(na),MAPs{2});
A0bar = kron(MAPa{1},eye(ns));
B0 = krons(MAPa{1},eye(ns));
for b=1:maxbatch
    B1{b} = kron(MAPa{2}*pbatcha(b),eye(ns));
end

% first few levels of Q matrix
%Z = 0*B0;
%Q=[B0,cell2mat(B1),Z,Z,Z; 
% A_1,A0,cell2mat(A1),Z,Z;
% Z,A_1,A0,cell2mat(A1),Z;
% Z,Z,A_1,A0,cell2mat(A1)];

end