function [G,lG]=pfqn_rgf(L,N,mi)
if nargin==2
    [L,~,J]=unique(L,'rows');
    mi=[];
    for i=1:size(L,1)
        mi(1,i)=sum(J==i);
    end
end
[M,R]=size(L);
Nck=zeros(sum(N)+sum(mi),sum(N)+sum(mi));
Dnk=cell(1,sum(N));
Gsingle=zeros(1,2*M+1);
iset=pprod(ones(1,R-1)*(M-1));
G=0;
Ir=cell(R,1);
Jr=zeros(M,1);
ctr=0;gctr=0;
while iset~=-1
    g=rgfaux(L,N,mi,iset);
    G=G+g;
    ctr=ctr+1;
    iset=pprod(iset,ones(1,R-1)*(M-1));
end
lG=log(G);

    function G=rgfaux(L,N,mi,iset)
        [M,R]=size(L);
        X=L;
        G=0;
        %% determine the loadings at the different steps
        for r=R:-1:2
            for l=1:M
                if l~=(iset(r-1)+1)
                    for p=1:(r-1)
                        X(l,p)=(X(iset(r-1)+1,r)*X(l,p)-X(l,r)*X(iset(r-1)+1,p))/(X(iset(r-1)+1,r)-X(l,r));
                    end
                end
            end
        end
        %% determine the multiplicities at the different steps
        % for each class we keep track of the couple (j,k), where k is the
        % index in the summation of the combinations that sum to j-1
        jstate=zeros(1,R);
        kstate=zeros(1,R);
        jmax=zeros(1,R); % maximum range of jstate
        kmax=zeros(1,R); % maximum range of kstate
        mistate=zeros(M,R); % multiplicities in the step of class r

        % initialize the elements of the state vectors in class R
        mistate(:,R)=mi(:);
        r=R;
        g=zeros(1,R);
        c=zeros(1,R);
        Jr=zeros(M,1);

        % initialize DFS recursion on jstate and kstate
        while r>1
            jmax(r)=mistate(iset(r-1)+1,r);
            jstate(r)=1;
            I{r}=mchoose(M-1,jstate(R)-1);
            kmax(r)=size(I{r},1);
            kstate(r)=1;
            if iset(r-1)+1>1
                Jr(1:(iset(r-1)+1-1),1)=I{r}(kstate(r),1:(iset(r-1)+1-1));
            end
            Jr(iset(r-1)+1,1)=N(r)+1-jstate(r);
            if iset(r-1)+1<M
                Jr((iset(r-1)+1+1):M,1)=I{r}(kstate(r),iset(r-1)+1:(M-1));
            end
            mistate(:,r-1)=mistate(:,r)+Jr;
            % update g
            g(r)=(-1)^(jstate(r)+1)*nck(mistate(iset(r-1)+1,r)+N(r)-jstate(r),N(r))*X(iset(r-1)+1,r)^N(r);
            for l=1:M
                if l~=iset(r-1)+1
                    g(r) = g(r)*(X(iset(r-1)+1,r)/(X(iset(r-1)+1,r)-X(l,r)))^mistate(l,r);
                end
            end
            
            c(r)=1;
            for l=1:M
                if l~=iset(r-1)+1
                    c(r) = c(r) * nck(mistate(l,r)+Jr(l)-1,Jr(l)) * (X(l,r)/(X(iset(r-1)+1,r)-X(l,r)))^Jr(l);
                end
            end
            r=r-1;
        end

        r=2;
        while r<R+1
            while jstate(r)<=jmax(r)
                % update g
                g(r)=(-1)^(jstate(r)+1)*nck(mistate(iset(r-1)+1,r)+N(r)-jstate(r),N(r))*X(iset(r-1)+1,r)^N(r);
                for l=1:M
                    if l~=iset(r-1)+1
                        g(r) = g(r) * (X(iset(r-1)+1,r)/(X(iset(r-1)+1,r)-X(l,r)))^mistate(l,r);
                    end
                end
                while kstate(r)<=kmax(r)
                    % update multiplicities of class r
                    if iset(r-1)+1>1
                        Jr(1:(iset(r-1)+1-1),1)=I{r}(kstate(r),1:(iset(r-1)+1-1));
                    end
                    Jr(iset(r-1)+1,1)=N(r)+1-jstate(r);
                    if iset(r-1)+1<M
                        Jr((iset(r-1)+1+1):M,1)=I{r}(kstate(r),iset(r-1)+1:(M-1));
                    end
                    mistate(:,r-1)=mistate(:,r)+Jr;
                    % update c
                    c(r)=1;
                    for l=1:M
                        if l~=iset(r-1)+1
                            c(r) = c(r) * nck(mistate(l,r)+Jr(l)-1,Jr(l)) * (X(l,r)/(X(iset(r-1)+1,r)-X(l,r)))^Jr(l);
                        end
                    end
                    if r==2 % if processing the last class accumulate G
                        G=G+prod(g(2:R))*prod(c(2:R))*gsingle(X(:,1),N(1),mistate(:,1)');
                        kstate(r)=kstate(r)+1; % go on
                    else
                        break
                    end
                end
                if r==2
                    jstate(r)=jstate(r)+1;
                    I{r}=mchoose(M-1,jstate(r)-1);
                    kmax(r)=size(I{r},1);
                    kstate(r)=1;
                elseif r>2 % if just updating jstate and kstate for r>2
                    r=r-1;
                    jmax(r)=mistate(iset(r-1)+1,r);
                    kstate(r)=1;
                    jstate(r)=1;
                    I{r}=mchoose(M-1,jstate(r)-1);
                    kmax(r)=size(I{r},1);
                end
            end
            % find the largest class that needs updating
            r=r+1;
            if r>R
                break;
            end
            while (kmax(r)<=kstate(r) && jmax(r)<=jstate(r))
                r=r+1;
                if r>R
                    break;
                end
            end
            % if completed go to the next iset element
            if r>R
                break;
            end
            % update class r state
            if kmax(r)>kstate(r)
                kstate(r)=kstate(r)+1;
            else
                jstate(r)=jstate(r)+1;
                kstate(r)=1;
                I{r}=mchoose(M-1,jstate(r)-1);
                kmax(r)=size(I{r},1);
            end
            jmax(1:(r-1))=0;
            jstate(1:(r-1))=0;
            kmax(1:(r-1))=0;
            kstate(1:(r-1))=0;
            mistate(:,1:(r-1))=0;
        end

        function I=mchoose(n,k)
            I=Dnk{1,k+1};
            if isempty(I)
                I=multichoose(n,k);
                Dnk{1,k+1}=I;
            end
        end

        function b=nck(c,k)
            b=Nck(c+1,k+1);
            if b==0
                b=nchoosek(c,k);
                Nck(c+1,k+1)=b;
            end
        end

        function G=gsingle(L,n,mult)
            L=round(L*1e6)/1e6;
            row=matchrow(Gsingle(:,1:2*M),[L',mult]);
            if row==-1
                gctr=gctr+1;
                [~,~,~,~,lG]=pfqn_mva(L,n,0*n,mult);
                G=exp(lG);
                Gsingle(end+1,:)=[L',mult,G];
                return;
            end
            G=Gsingle(row,2*M+1);
        end

    end
end

