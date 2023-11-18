function [G,lG,XN,QN]=pfqn_comombtf_java(Din,Nin,Zin)
% global isNCLibLoaded;
% if isempty(isNCLibLoaded)
%     javaaddpath(which('pfqn_nclib.jar')); % now in lineStart
%     isNCLibLoaded = true;
% end
import DataStructures.*; %#ok<SIMPT>
import QueueingNet.*; %#ok<SIMPT>
import Utilities.*; %#ok<SIMPT>

%Din=Din(sum(Din,2)>GlobalConstants.FineTol,:);

% rescale
numdigits = ceil(max(abs(log10(abs([Din(:);Zin(:)])))));
scaleexp = min(numdigits,10);  % java.lang.Integer takes max 10 digits
scale = 10^(scaleexp);
Din = round(Din*scale);
Zin = round(sum(Zin*scale,1));

[M,R]=size(Din);
[~,I]=sort(sum(Din>0,1),'descend');
Din = Din(:,I);
Nin = Nin(I);
Zin = Zin(:,I);

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

qnm = DataStructures.QNModel(M,R);
qnm.N = PopulationVector(N);
qnm.Z = EnhancedVector(Z);
qnm.multiplicities = MultiplicitiesVector(mult);
qnm.D = D;
comom = CoMoMBTFSolver(qnm);
comom.computeNormalisingConstant();
G = qnm.getNormalisingConstant();
lG = G.log();
lG = lG - sum(Nin)*log(scale);
G = exp(lG);

if nargout > 2
    comom.computePerformanceMeasures();
    XNbig = qnm.getMeanThroughputs();
    XN = zeros(1,R);
    for r=1:length(XNbig)
        XN(r) = XNbig(r).approximateAsDouble;
    end
    XN = XN * scale;
    XN(:,I) = XN;
    QNbig = qnm.getMeanQueueLengths();
    QN = zeros(M,R);
    for i=1:M
        for r=1:R
            QN(i,r) = QNbig(i,r).approximateAsDouble;
        end
    end
    QN(:,I)=QN;
end
end