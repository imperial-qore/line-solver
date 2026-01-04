function [ lowerboundsolutions, upperboundsolutions, integralsolutions,Lzeromulti,Lnegmulti,Lposmulti,Anegmulti,Aposmulti] = AdditiveDecomposition(Qmulti,driftregimesmulti,Bmulti)

for multi=1:size(Qmulti,3)
    %ORDERED SCHUR 
    zerodrift=find(driftregimesmulti(multi,:)==0);
    drifts=driftregimesmulti(multi,:);
    Q=Qmulti(:,:,multi);
    
    if ~isempty(zerodrift)
        for i=1:length(zerodrift)
            last=length(Q)+1-i;
            nextzero=zerodrift(i);
            tempQrow=Q(last,:);
            Q(last,:)=Q(nextzero,:);
            Q(nextzero,:)=tempQrow;
            tempQcolumn=Q(:,last);
            Q(:,last)=Q(:,nextzero);
            Q(:,nextzero)=tempQcolumn;
            tempdriftsrow=drifts(last);
            drifts(last)=drifts(nextzero);
            drifts(nextzero)=tempdriftsrow;
        end

        R=diag(drifts);
        nindex=length(drifts)-length(zerodrift);
        driftsn=drifts(1:nindex);
        driftsz=drifts(nindex+1:length(drifts));
        Qnn=Q(1:nindex,1:nindex);
        Qzz=Q(nindex+1:length(drifts),nindex+1:length(drifts));
        Qnz=Q(1:nindex,nindex+1:length(drifts));
        Qzn=Q(nindex+1:length(drifts),1:nindex);

        Qin=Qnn-(Qnz/Qzz)*Qzn;
        Rn=diag(driftsn);
        zeroconverter=-Qnz/Qzz;
    else
        Qin=Q;
        Rn=diag(drifts);
        R=diag(drifts);
    end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55   
    %%Q=Qmulti(:,:,multi);
    %%R=diag(driftregimesmulti(multi,:));
    
    upperB=Bmulti(multi);
    if multi==1
        lowerB=0;
    else
        lowerB=Bmulti(multi-1);
    end
    A=Qin/Rn;
%     right=Rn*ones(size(Rn,1),1); 
%     e1=[1;zeros(size(Rn,1)-1,1)];
%     u=right-norm(right,2)*e1;
%     Q1=eye(size(Rn))-2*(u*u')/(u'*u);
%     temp=Q1*A*Q1;
%     Qbar=temp(2:size(Rn,1),2:size(Rn,1));
%     [Z,M]=schur(Qbar);
%     [Z,M]=ordschur(Z,M,'lhp');
%     Z1=Q1*blkdiag(eye(1),Z);
%     D1=Z1'*A*Z1;
     [Z,D1]=schur(A); 

    lengthD=length(D1);
     
    c=zeros(1, lengthD);
    positivecount=0; 
    negativecount=0;
    zerocount=0;
    for i=1: lengthD
    if D1(i,i)>0.0000001
        positivecount=positivecount+1;
        c(i)=1;
    elseif D1(i,i)<-0.0000001
        negativecount=negativecount+1;
        c(i)=2;
    else
        zerocount=zerocount+1;
        c(i)=3;
    end
    end
          [Z1,D1]=ordschur(Z,D1,c);
    negindex=zerocount+1; 
    posindex=negindex+negativecount;
    
    %compute X1
    A0=D1(1:zerocount,1:zerocount);
    k1=D1(1:zerocount,negindex :  lengthD);
    k2=D1(negindex:  lengthD,negindex : lengthD);
    X1=sylvester(A0,-k2,-k1);
   
    %compute X2
    lengthk2=length(k2);
    posk2index=1+negativecount;
    Aneg= k2(1:negativecount,1:negativecount);
    Apos= k2(posk2index:lengthk2,posk2index:lengthk2);
    Aposneg= k2(1:negativecount,posk2index:lengthk2);
    X2=sylvester(Aneg,-Apos,-Aposneg);

    %Compute Y and Corresponding A(T)
    e1=eye(lengthD);
    e1(1:size(X1,1),(size(e1)-size(X1,2)+1):lengthD)=X1;

    e2=eye(size(X2,1)+size(X2,2));
    e2(1:size(X2,1),(size(e2)-size(X2,2)+1):size(e2))=X2;

    e3=eye(lengthD);
    e3((lengthD-size(e2)+1):lengthD,(lengthD-size(e2)+1):lengthD)=e2;

    Y=Z1*e1*e3;
    Yinv=inv(Y);

    T=Y\A*Y;

%     if zerocount ~= 0
%         X1_temp = sylvester(-D1(1:zerocount,1:zerocount),D1((zerocount+1):length(D1),(zerocount+1):length(D1)),-D1(1:zerocount,(zerocount+1):length(D1)));
% 
% %         X1_temp = -D1(1:zerocount,(zerocount+1):length(D1))/D1((zerocount+1):length(D1),(zerocount+1):length(D1));
%         X2_temp = sylvester(-D1((zerocount+1):(zerocount+negativecount),(zerocount+1):(zerocount+negativecount)),D1((zerocount+negativecount+1):length(D1),(zerocount+negativecount+1):length(D1)),-D1((zerocount+1):(zerocount+negativecount),(zerocount+negativecount+1):length(D1)));
%         Y_temp = Z1*[eye(size(X1_temp,1)) -X1_temp;zeros(size(X1_temp,2),length(D1)-size(X1_temp,2)) eye(size(X1_temp,2))]*[eye(zerocount) zeros(zerocount,negativecount) zeros(zerocount,positivecount);zeros(negativecount,zerocount) eye(negativecount) -X2_temp;zeros(positivecount,zerocount) zeros(positivecount,negativecount) eye(positivecount)];
%     else
%         X2_temp = sylvester(-D1(1:negativecount,1:negativecount),D1((negativecount+1):length(D1),(negativecount+1):length(D1)),-D1(1:negativecount,(negativecount+1):length(D1)));
%         Y_temp = Z1*[eye(negativecount) -X2_temp;zeros(positivecount,negativecount) eye(positivecount)];
%     end
%     
%     
%     % diagonal_mat=(inv(Y_temp))*A*Y_temp;
%     diagonalized_A=double((Y_temp\A)*Y_temp);
%     diagonalized_A_num(1:length(Qin),length(Qin)*(i-1)+1:length(Qin)*(i-1)+length(Qin)) = diagonalized_A;
%     
%     diag_A_zero = diagonalized_A(1:zerocount,1:zerocount);
%     diag_A_negative = diagonalized_A(zerocount+1:zerocount+negativecount,zerocount+1:zerocount+negativecount);
%     diag_A_positive = diagonalized_A(zerocount+negativecount+1:length(diagonalized_A),zerocount+negativecount+1:length(diagonalized_A));
%     
%     Y_temp_inv_trans = (Y_temp\eye(length(Y_temp)))';

    %RETURN Lneg Lpos Lzero 
    Lzero=Yinv(1:zerocount,:); 
    Lneg=Yinv(negindex:zerocount+negativecount,:);
    Lpos=Yinv(posindex:lengthD,:); 
    
    if ~isempty(zerodrift)
        Lzero=[Lzero,Lzero*zeroconverter];
        Lpos=[Lpos,Lpos*zeroconverter];
        Lneg=[Lneg,Lneg*zeroconverter];
        for i=1:length(zerodrift)
            last=length(Q)+1-i;
            nextzero=zerodrift(i);
            templzerorow=Lzero(:,last);
            templposrow=Lpos(:,last);
            templnegrow=Lneg(:,last);
            Lzero(:,last)=Lzero(:,nextzero);
            Lneg(:,last)=Lneg(:,nextzero);
            Lpos(:,last)=Lpos(:,nextzero);
            Lpos(:,nextzero)=templposrow;
            Lneg(:,nextzero)=templnegrow;
            Lzero(:,nextzero)=templzerorow;

        end
    end
    
    %Return values f(0),f(b) and F(B) when coefficients a=1 
    syms  x
    integralsolution=[Lzero*(upperB-lowerB); (Aneg\(expm(Aneg*(upperB-lowerB))-eye(size(Aneg,1))))*Lneg;((-Apos)\(expm(-Apos*(upperB-lowerB))-eye(size(Apos,1))))*Lpos];
    lowsolution=[Lzero; expm(Aneg*0)*Lneg;expm((-Apos)*(upperB-lowerB))*Lpos];
    upsolution=[Lzero; expm(Aneg*(upperB-lowerB))*Lneg;expm((-Apos)*(0))*Lpos];
   % f_bar_result = [Lzero; expm(Aneg*x)*Lneg;expm(-Apos*(B-x))*Lpos]
    
    lowerboundsolutions{multi}=lowsolution;
    upperboundsolutions{multi}=upsolution;
    integralsolutions{multi}=integralsolution;
    Lzeromulti{multi}=Lzero; 
    Lposmulti{multi}=Lpos;
    Lnegmulti{multi}=Lneg; 
    Anegmulti{multi}=Aneg;
    Aposmulti{multi}=Apos;
end
end