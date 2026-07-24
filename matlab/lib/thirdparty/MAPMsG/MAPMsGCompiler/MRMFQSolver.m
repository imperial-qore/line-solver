function [coefficients,boundaries,Lzeromulti,Lnegmulti,Lposmulti,Anegmulti,Aposmulti] = MRMFQSolver(Qregimes,Qboundaries,driftregimes,driftboundaries,B)

    syms x
    nonzeroboundaries=ones(size(driftboundaries));
    zerolowerpdf=zeros(size(driftregimes));
    zeroupperpdf=zeros(size(driftregimes));
    for level=1:size(nonzeroboundaries,1)
        for state=1:size(nonzeroboundaries,2)
            if level==1 
                if driftregimes(level,state)>0
                    nonzeroboundaries(level,state)=0;
                end
            elseif level==size(nonzeroboundaries,1)
                if driftregimes(level-1,state)<0
                    nonzeroboundaries(level,state)=0;
                end
            else
                if (driftregimes(level,state)>0 && driftregimes(level-1,state)>0) || (driftregimes(level,state)<0 && driftregimes(level-1,state)<0)||(driftregimes(level,state)>0 && driftregimes(level-1,state)<0 && (driftboundaries(level,state)~=0))
                    nonzeroboundaries(level,state)=0;
                end
                if (driftregimes(level,state)>0 && (driftboundaries(level,state)<=0))
                    zerolowerpdf(level,state)=1;
                end
                if (driftregimes(level-1,state)<0 && (driftboundaries(level,state)>=0))
                    zeroupperpdf(level-1,state)=1;
                end
            end
            
        end
    end
       Rregimes=Qregimes;
    for k=1:size(driftregimes,1)
        Rregimes(:,:,k)=diag(driftregimes(k,:));
    end

%obtainin L, A matrices and f(0),f(b) and F(B) when coefficients a=1
    [lowerboundsolutions, upperboundsolutions, integralsolutions,Lzeromulti,Lnegmulti,Lposmulti,Anegmulti,Aposmulti]=AdditiveDecomposition(Qregimes,driftregimes,B);
    rowindex=1;
    columnindex=1;
    for funcindex=1:size(driftboundaries,1)
        zeroboundariesindexes=find(nonzeroboundaries(funcindex,:)==0);
        Qeliminated=Qboundaries(:,:,funcindex);
        Qeliminated(zeroboundariesindexes,:)=[];
        if funcindex==1
            H(1:size(Qeliminated,1),1:size(Rregimes(:,:,funcindex),2))=-Qeliminated;
            rowindex=rowindex+size(Qeliminated,1);
            H(rowindex:rowindex+size(lowerboundsolutions{funcindex},1)-1,1:size(Rregimes(:,:,funcindex),2))=lowerboundsolutions{funcindex}*Rregimes(:,:,funcindex);
            columnindex=columnindex+size(Rregimes(:,:,funcindex),2);
        elseif funcindex==size(driftboundaries,1)
            H(rowindex:rowindex+size(upperboundsolutions{funcindex-1},1)-1,columnindex:columnindex+size(Rregimes(:,:,funcindex-1),2)-1)=upperboundsolutions{funcindex-1}*Rregimes(:,:,funcindex-1);
            rowindex=rowindex+size(upperboundsolutions{funcindex-1},1);

            H(rowindex:rowindex+size(Qeliminated,1)-1,columnindex:columnindex+size(Rregimes(:,:,funcindex-1),2)-1)=Qeliminated;
        else
            for state=1:size(driftboundaries,2)
                if zeroupperpdf(funcindex-1,state)==1;
                    temp=upperboundsolutions{funcindex-1};
                    H(rowindex:rowindex+size(upperboundsolutions{funcindex-1},1)-1,columnindex:columnindex)=temp(:,state);
                     columnindex=columnindex+1;
                end
            end
            H(rowindex:rowindex+size(upperboundsolutions{funcindex-1},1)-1,columnindex:columnindex+size(Rregimes(:,:,funcindex-1),2)-1)=upperboundsolutions{funcindex-1}*Rregimes(:,:,funcindex-1);
            rowindex=rowindex+size(upperboundsolutions{funcindex-1},1);
            H(rowindex:rowindex+size(Qeliminated,1)-1,columnindex:columnindex+size(Rregimes(:,:,funcindex),2)-1)=Qeliminated;
            rowindex=rowindex+size(Qeliminated,1);
            H(rowindex:rowindex+size(lowerboundsolutions{funcindex},1)-1,columnindex:columnindex+size(Rregimes(:,:,funcindex),2)-1)=-lowerboundsolutions{funcindex}*Rregimes(:,:,funcindex);
            columnindex=columnindex+size(Rregimes(:,:,funcindex),2);
            for state=1:size(driftboundaries,2)
                if zerolowerpdf(funcindex,state)==1;
                    temp2=lowerboundsolutions{funcindex};
                    H(rowindex:rowindex+size(lowerboundsolutions{funcindex},1)-1,columnindex:columnindex)=temp2(:,state);
                     columnindex=columnindex+1;
                end
            end
        end
    end

        Hbar=H;
%         Hbar(:,1)=0;
%         Hbar(1,1)=1;
%         b=Hbar(:,1);
%         ters=inv(Hbar);
% %         ters2=inv(ters);
%             z=Hbar(:,1)'/Hbar;
%          z=null(H');
%          z=z';
Hbar(:,1) = ones(size(H,1),1); 

z = eye(1,size(H,2))/Hbar;

        rowindex2=1; 
        columnindex2=1; 
%         for i=1:size(driftboundaries,1)
%             if i<size(driftboundaries,1)
%                 if i==1
%                     ithblocksize=sum(nonzeroboundaries(i,:) == 1);
%                 else
%                     ithblocksize=nextblocksize;
%                 end
%                     nextblocksize=sum(nonzeroboundaries(i+1,:) == 1)+size(lowerboundsolutions{i},1);
%                     
%                     Hl{i}=Hbar(rowindex2+ithblocksize:rowindex2+ithblocksize+nextblocksize-1,columnindex2:columnindex2+ithblocksize-1);
%                     Hu{i}=Hbar(rowindex2:rowindex2+ithblocksize-1,columnindex2+ithblocksize:columnindex2+ithblocksize+nextblocksize-1);
%             else
%                  ithblocksize=nextblocksize;   
%             end
%             Hm{i}=Hbar(rowindex2:rowindex2+ithblocksize-1,columnindex2:columnindex2+ithblocksize-1);
%             rowindex2=rowindex2+ithblocksize;
%             columnindex2=columnindex2+ithblocksize;
%         end
        
        index=1;
        normalizationcoefmasses=0;
        normalizationcoefintegrals=0;
        for level=1:size(driftboundaries,1)

            countnonzeroboundaries=sum(nonzeroboundaries(level,:) == 1);
            if countnonzeroboundaries>0
                normalizationcoefmasses=normalizationcoefmasses+sum(z(index:index+countnonzeroboundaries-1));
                %boundaries{level}=z(index:index+countnonzeroboundaries-1);
                index=index+countnonzeroboundaries;
            end
           
            if level<size(driftboundaries,1)
                normalizationcoefintegrals=normalizationcoefintegrals+ sum(z(index:index+size(integralsolutions{level},1)-1)*integralsolutions{level});
                %coefficients{level}=z(index:index+size(integralsolutions(:,:,level),1)-1);
                index=index+size(integralsolutions{level},1);
              
            end
        end
        normalizationcoef=normalizationcoefintegrals+normalizationcoefmasses;
        finalZ=z/normalizationcoef;
        index2=1;
         for level=1:size(driftboundaries,1)
             countnonzeroboundaries=sum(nonzeroboundaries(level,:) == 1);
             boundaries{level}=finalZ(index2:index2+countnonzeroboundaries-1);
            index2=index2+countnonzeroboundaries;
            if level<size(driftboundaries,1)
               
                coefficients{level}=finalZ(index2:index2+size(integralsolutions{level},1)-1);
                index2=index2+size(integralsolutions{level},1);
              
            end
         end
        for level=1:size(driftboundaries,1)
            l=1;
            for k=1:size(nonzeroboundaries,2)
                if nonzeroboundaries(level,k) == 1
                    bound=boundaries{level};
                    nonzeroboundaries(level,k) = bound(l);
                    l=l+1;
                end
            end
            boundaries{level}=nonzeroboundaries(level,:);
        end

   

end
