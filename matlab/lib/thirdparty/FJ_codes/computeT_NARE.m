function T2 = computeT_NARE(D0, D1, S, A_jump)
tic
m = size(S,1);
ma = size(D0,1);
ms = m/ma;

H = [kron(eye(ms),D0)  kron(eye(ms), D1);
    -kron(A_jump,eye(ma)) -S];

% H = [kron(eye(ms),D0)  kron(eye(ms), D1);
%     -kron(A_jump,eye(ma)) -S];

%n = size(A,1);
%m = size(D,1);
%H = [A,-B;-C,-D];
% H = full(H);

[U,Q] = schur(H,'real');
time_Schur = toc;
e = ordeig(Q); % Eigenvalues of quasitriangular matrices
%[es,is] = sort(real(e), 'descend');
[es,is] = sort(real(e), 'ascend');
sel = zeros(2*m,1);
sel(is(1:m)) = 1;
Q1 = ordschur(U,Q,sel); % Sorting the Schur form of H
timeOrdSchur = toc;
timeOrdSchur = timeOrdSchur-time_Schur;
X = Q1(m+1:2*m,1:m)/Q1(1:m,1:m); % \bar{L}

T2 = S + X*kron(eye(ms), D1);
e2 = toc;


res_norm=norm(T2*X+X*kron(eye(ms),D0)+kron(A_jump,eye(ma)),inf)
fprintf('Final Residual Error for T matrix: %d\n',res_norm);