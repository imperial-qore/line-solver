function [alpha, A] = FluFluSTD (Qin, Rin, Qout, Rout, srv0stop, transToPH)

    if ~exist('transToPH','var')
        transToPH = false;
    end
    
    % solve special fluid queue
    Iin = eye(size(Qin));
    Iout = eye(size(Qout));

    Rh = kron(Rin,Iout) - kron(Iin,Rout);
    Qh = kron(Qin, Rout) + kron(Rin, Qout);       
    [massh, inih, Kh, cloh] = GeneralFluidSolve (Qh, Rh);

    % sojourn time density in case of 
    % srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda
    % srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu    
    
    lambda = sum(CTMCSolve(Qin)*Rin);
    mu = sum(CTMCSolve(Qout)*Rout);
    
    if transToPH   
        % convert result to PH representation
        Delta = diag(linsolve(Kh',-inih')); % Delta = diag (inih*inv(-Kh));
        A = inv(Delta)*Kh'*Delta;       
        if ~srv0stop        
            alpha = sum(Delta*cloh*kron(Rin,Iout)/lambda,2)';
        else
            alpha = sum(Delta*cloh*kron(Rin,Rout)/lambda/mu,2)';
        end        
    else
        % convert result to ME representation
        if ~srv0stop
            B = TransformToOnes(sum(cloh*kron(Rin,Iout)/lambda,2));
        else
            B = TransformToOnes(sum(cloh*kron(Rin,Rout)/lambda/mu,2));
        end
        iB = inv(B);
        A = B*Kh*iB;
        alpha = inih*inv(-Kh)*iB;
    end
end

