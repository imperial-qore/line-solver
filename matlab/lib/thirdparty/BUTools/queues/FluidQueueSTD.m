function [alpha, A] = FluidQueueSTD (Q, Rin, Rout, Q0, transToPH)

    if ~exist('Q0','var')
        Q0 = [];
    end

    if ~exist('transToPH','var')
        transToPH = false;
    end
    
    % obtain solution
    [mass0, ini, K, clo] = GeneralFluidSolve (Q, Rin-Rout, Q0);
   
    N = size(Q,1);
    iniKi = linsolve(K',-ini')'; % iniki = ini*inv(-K);
    lambda = sum(mass0*Rin + iniKi*clo*Rin);
    if transToPH
        % transform it to PH
        Delta = diag(iniKi/lambda);
        alpha = reshape(clo*Rin,1,N*length(ini))*kron(eye(N),Delta);
        A = kron(Rout, inv(Delta)*K'*Delta) + kron(Q, eye(size(K)));
    else
        B = TransformToOnes(reshape(inv(-K)*clo*Rin,N*length(ini),1));
        Bi = inv(B);
        alpha = kron(ones(1,N), ini/lambda)*Bi;
        A = B*(kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K)))*Bi;        
    end
end

