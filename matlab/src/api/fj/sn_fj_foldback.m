function [QN,UN,RN,TN,CN,XN] = sn_fj_foldback(QN,UN,RN,TN,CN,XN,fjclassmap,Korig)
% [QN,UN,RN,TN,CN,XN] = SN_FJ_FOLDBACK(QN,UN,RN,TN,CN,XN,FJCLASSMAP,KORIG)
%
% Fold the auxiliary-class columns of the average metrics computed on an
% FJ tag-augmented struct (ModelAdapter.fjtag) back into the original
% classes: queue lengths, utilizations and throughputs of the sibling
% classes are exact aggregates of the original class they were forked
% from; response times are recomputed by Little's law after folding.
%
% QN,UN,RN,TN are (nstations x Kaug); CN,XN are (1 x Kaug); the outputs
% retain only the first Korig columns.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Kaug = length(fjclassmap);
for a=1:Kaug
    r = fjclassmap(a);
    if r > 0
        QN(:,r) = QN(:,r) + QN(:,a);
        UN(:,r) = UN(:,r) + UN(:,a);
        TN(:,r) = TN(:,r) + TN(:,a);
    end
end
QN = QN(:,1:Korig);
UN = UN(:,1:Korig);
TN = TN(:,1:Korig);
RN = zeros(size(QN));
for r=1:Korig
    for ist=1:size(QN,1)
        if TN(ist,r) > 0
            RN(ist,r) = QN(ist,r) / TN(ist,r);
        end
    end
end
% system metrics are measured on the original class columns only
CN = CN(:,1:Korig);
XN = XN(:,1:Korig);

end
