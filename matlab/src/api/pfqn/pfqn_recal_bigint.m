function [G,lG]=pfqn_recal_bigint(Din,Nin,Zin)
try
    import DataStructures.*; %#ok<SIMPT>
    import QueueingNet.*; %#ok<SIMPT>
    import DataStructures.*; %#ok<SIMPT>
    import Utilities.*; %#ok<SIMPT>
catch
    javaaddpath(which('pfqn_nclib.jar'));
    import DataStructures.*; %#ok<SIMPT>
    import QueueingNet.*; %#ok<SIMPT>
    import DataStructures.*; %#ok<SIMPT>
    import Utilities.*; %#ok<SIMPT>
end

Din=Din(sum(Din,2)>GlobalConstants.FineTol,:);
% rescale
numdigits = max(arrayfun(@(e) numel(num2str(e)), [Din(:);Zin(:)]));
scaleexp = min(numdigits,8);  % java.lang.Integer takes max 10 digits
scale = 10^(scaleexp);
Din = round(Din*scale); 
Zin = round(sum(Zin*scale,1));

[M,R]=size(Din);
[~,I]=sort(sum(Din,1),'descend');
Din=Din(:,I);
Nin=Nin(1,I);
Zin=Zin(:,I);

N = javaArray('java.lang.Integer',R);
for r=1:R
    N(r) = java.lang.Integer(Nin(r));
end

mult = javaArray('java.lang.Integer',M);
for i=1:M
    mult(i) = java.lang.Integer(1);
end

Z= javaArray('java.lang.Integer',R);
for r=1:R
    Z(r) = java.lang.Integer(Zin(r));
end

D = javaArray('java.lang.Integer',M,R);
for i=1:M
    for r=1:R
        D(i,r) = java.lang.Integer(Din(i,r));
    end
end

qnm = QNModel(M,R);
qnm.N = PopulationVector(N);
qnm.Z = EnhancedVector(Z);
qnm.multiplicities = MultiplicitiesVector(mult);
qnm.D = D;
comom = RECALSolver(qnm);
comom.computeNormalisingConstant();
G = qnm.getNormalisingConstant();
lG = G.log();
lG = lG - sum(Nin)*log(scale);
G = exp(lG);
end
