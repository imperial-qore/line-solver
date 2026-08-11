% O. Gursoy, K. A. Mehr, N. Akar:
% The MAP/M/s + G Call Center Model with Generally Distributed Patience Times: Steady-state Solution and First Passage Time Distribution
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Steady-state and first passage time solver for MAP/M/s + G queues
%
% INPUTS
%-----------------------------
% -SOLUTION: Solution type(steady state or first passage actual/virtual)
% -SERVERSIZE: server size of the call center problem
% -mu
% -C: C of MAP(C,D) 
% -D: D of MAP(C,D) 
% -MAPSIZE: Order of MAP(C,D) distribution, Service Rate(mu)
% -QUANTIZATION: number of regimes that will quantize abandonment function ga
% -BoundaryLevels: Boundary levels of Regimes
% -ga: Abandonment Probability in Each Regime (starts with a 0 value)
% -gabanfunction: for continious abandonment functions enter ccdf of abandonement as a function of x
%
%
% INPUTS for first passage times(these are not necessary for steady state
% calculations)
%-----------------------------
% -b: threshold level 
% -tau: time horizon
% -CMEorErlang: indicator of which distribution is chosen for approximation
% -OrderofPHCME: order of the chosen approximator (25, 51 or 101)
% -pi0: initial server occupancy distribution
% -theta0: initial MAP state distribution
%
% OUTPUTS
%---------------------------
% -steadyStateResult: 1x7 steady state result vector including PR{W=0}
% PR{W=0|S}   PR{A}    E{W|S}   Var{W|S} Fw|s,w>0(.1) Fw|s,w>0(.2)
% respectively(when SOLUTION=1).
%
% -FirstPassageTimeofVirtual: First passage time probability for the virtual waiting time(when SOLUTION=2).
%
% -FirstPassageTimeofActual: First passage time probability for the actual
% waiting time(when SOLUTION=3).
%
%
% IMPORTANT:
%-----------
% For first passage time cases the initial wait time, a=0.
%
% 
%_________________________________________________________________
%/////////////////////////////////////////////////////////////////

clc; 
close all;
clear all;
 
tol=10^-5;

%INPUTS___________________________________________________________
%/////////////////////////////////////////////////////////////////

%solution=1 if steady state,  solution=2 if first passage distribution of virtual waiting time,
%solution=3 if first passage distribution of actual waiting time,
SOLUTION=3;

%server size, state size of MAP(C,D), service rate and QUANTIZATION (regime count)
SERVERSIZE=10;
MAPSIZE=2;
mu=1;
QUANTIZATION=11;

% C,D matrices (2 state MAP(C,D) is default)
%-------------------------------------------
decay=0.95;
    rho=0.99;
    S=SERVERSIZE;
    m=MAPSIZE;
    C_onsquared = 16;
    mean_arrival = 1/(S*mu*rho); 
    e = ones(m,1);

    p_1 = 0.5*(1 + sqrt((C_onsquared-1)/(C_onsquared+1)));
    p_2 = 1-p_1;
    mu_1 = 2*p_1 / mean_arrival;
    mu_2 = 2*p_2 / mean_arrival;

    v = [p_1 p_2];
    T = [-mu_1 0;0 -mu_2];
    T1 = T;
    T0 = [mu_1; mu_2];
    s = T0;
    lambda = 1/mean_arrival;
    D0 = T1;
    D1 = (1-decay)*T0*v - decay*T1 ; 
C=D0;
D=D1; 
em=ones(MAPSIZE,1);
lmap=D*em;



%abandonment prob for regimes (ga) and boundary levels (BoundaryLevels)
%---------------------------------------------------------------------
%(piecewise constant abandonment function is default) 
BoundaryLevels = linspace(0,10,QUANTIZATION);
ga=0.10*floor(BoundaryLevels);
BoundaryLevels(1)=[];
BoundaryLevels(end+1)=10000000;
ga=[0 ga];


%Uncomment for continious abandonment function**
% interval=1/(QUANTIZATION);
% ga(1)=0;
% for k=1:QUANTIZATION
%     ga(k+1)=interval*k+ga(1);
% end
% regime=1;
% for i=1:100000000
%     x=i/1000000;
%     %enter ccdf of abandonement as a function of x**
%     gabanfunction=1-exp(-x)-x*exp(-x)-x^2*exp(-x)/2;
%     if (ga(regime+1)+tol>gabanfunction)&&(ga(regime+1)-tol<gabanfunction)
%         BoundaryLevels(regime)=x;
%         if regime<QUANTIZATION
%             regime=regime+1;
%         end
%     end
% end
% gatemp=ga;
% for i=2: length(ga)
%     gatemp(i)=(ga(i-1)+ga(i))/2;
% end
% ga=gatemp;



%if first passage distribution is chosen (SOLUTION=3) enter the parameters**
%------------------------------------------------------------
b=0.25;
tau=1;
%initial server occupancy distribution 
pi0=[1 0 0 0 0 0 0 0 0 0 0];
%initial MAP state distribution
theta0=[0 1];
%CMEorErlang=1 for CME approximation and CMEorErlang=0 for Erlangization
CMEorErlang=0;
%order of PH or CME, order=25 if OrderofPHCME=1, order=51 if OrderofPHCME=2, order=101 if OrderofPHCME=3,
OrderofPHCME=1;

%END OF INPUTS _______________________________________________________
%/////////////////////////////////////////////////////////////////


%CALCULATIONS_____________________________________________________
%/////////////////////////////////////////////////////////////////

%construction of infinitesial generators and drift matrices for both regimes (Qy,ydriftregimes)and
%boundaries (Qybounds,ydriftbounds)
I=eye(MAPSIZE);
Qy0=zeros((SERVERSIZE+1)*size(I,1),(SERVERSIZE+1)*size(I,1));
Qy=zeros((SERVERSIZE+1)*size(I,1),(SERVERSIZE+1)*size(I,1),QUANTIZATION);
Rydiag=-ones(1,size(Qy0,1));
Rydiag(size(Qy0,1)-size(I,1)+1:size(Qy0,1))=-Rydiag(size(Qy0,1)-size(I,1)+1:size(Qy0,1)); 
Ry=diag(Rydiag);
for row=1:SERVERSIZE+1    
    if row ==1
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE)=C;
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row)*MAPSIZE+1: (row)*MAPSIZE+MAPSIZE)=D;
    elseif row==SERVERSIZE+1
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE)=-(row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1: (row-2)*MAPSIZE+MAPSIZE)=(row-1)*mu*I;
    else
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE)=C-(row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1: (row-2)*MAPSIZE+MAPSIZE)=(row-1)*mu*I;
        Qy0((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row)*MAPSIZE+1: (row)*MAPSIZE+MAPSIZE)=D;
    end
end

for regimecount=1:QUANTIZATION
    for row=SERVERSIZE:SERVERSIZE+1 
        if row ==SERVERSIZE
            Qy((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE,regimecount)=(ga(regimecount+1))*D+C;
            Qy((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row)*MAPSIZE+1: (row)*MAPSIZE+MAPSIZE,regimecount)=(1-(ga(regimecount+1)))*D;
        elseif row==SERVERSIZE+1
            Qy((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE,regimecount)=-(row-1)*mu*I;
            Qy((row-1)*MAPSIZE+1: (row-1)*MAPSIZE+MAPSIZE, (row-2)*MAPSIZE+1: (row-2)*MAPSIZE+MAPSIZE,regimecount)=(row-1)*mu*I;
        end
    end
    ydriftregimes(regimecount,:)=Rydiag;
    Ryregimes(:,:,regimecount)=Ry;
end
Qybounds=cat(3,Qy0,Qy);
ydriftbounds=cat(1,Rydiag,ydriftregimes);
Rybounds=cat(3,Ry,Ryregimes);

if OrderofPHCME==1
    OrderofPH=25;
elseif OrderofPHCME==2
    OrderofPH=51;
elseif OrderofPHCME==3
    OrderofPH=101;
end

%Based on the chosen SOLUTION make adjustments for the input of MRMFQ solver
if SOLUTION==1 
    driftbound=ydriftbounds;
    driftregimes=ydriftregimes;
    Qregimes=Qy;
    Qbounds=Qybounds;
    BoundaryLevelsLast=BoundaryLevels;
elseif SOLUTION==2
     if CMEorErlang==1 
            n=OrderofPH;
            [MESystem,cv]= CMEParameterCalculator(n,tau);
            S=MESystem.A;
            S0=MESystem.b;
            alpha=MESystem.c;
     else
            alpha=zeros(1,OrderofPH);
            alpha(1)=1;
            ss=-ones(1,OrderofPH);
            S=diag(ss);
            for y=1:OrderofPH-1
               S(y,y+1)=1;
            end
            S=OrderofPH*S/tau;
            e=ones(length(S),1);
            S0=-S*e;
     end
     count=1;
     while BoundaryLevels(count)<b
        count=count+1;
     end
     BoundaryLevelsLast=BoundaryLevels(1:count);
     BoundaryLevelsLast(count)=b;
   
     for reg=1:length(BoundaryLevelsLast)
            eTilde= ones(1,length(Qybounds(:,:,reg+1)));
            eTilde(end-MAPSIZE+1:end)=0;
            Itilde=diag(eTilde);

            Qz(:,:,reg)=[0  zeros(1,length(kron(alpha,kron(pi0,theta0))));kron(S0,diag(Itilde)) (kron(eye(length(S)),Qy(:,:,reg))+kron(S,Itilde))]; 
            Qzbounds(:,:,reg+1)=[0  zeros(1,length(kron(alpha,kron(pi0,theta0))));kron(S0,diag(Itilde)) (kron(eye(length(S)),Qybounds(:,:,reg+1))+kron(S,Itilde))]; 
            
            Rz(:,:,reg)=diag([-1 ones(1,size(kron(eye(length(S)),Ryregimes(:,:,reg)),1))* kron(eye(length(S)),Ryregimes(:,:,reg))]);
            Rzbounds(:,:,reg+1)=diag([-1 ones(1,size(kron(eye(length(S)),Rybounds(:,:,reg+1)),1))* kron(eye(length(S)),Rybounds(:,:,reg+1))]);
     end
        Qzbounds(:,:,1)=[-1 kron(alpha,kron(pi0,theta0));kron(S0,diag(Itilde)) (kron(eye(length(S)),Qybounds(:,:,1))+kron(S,Itilde))];
        Rzbounds(:,:,1)=diag([-1 ones(1,size(kron(eye(length(S)),Rybounds(:,:,1)),1))* kron(eye(length(S)),Rybounds(:,:,1))]);
        Qzbounds(:,:,end)=[0  zeros(1,length(kron(alpha,kron(pi0,theta0)))); ones(size(Qzbounds(:,:,1),1)-1,1) -eye(size(Qzbounds(:,:,1),1)-1)];
        Rzbounds(:,:,end)=Rzbounds(:,:,1).*0;
        Rzbounds(1,1,end)=-1;
        
        Qregimes=Qz;
        Qbounds=Qzbounds; 
        driftbound(1,:)=[diag(Rzbounds(:,:,1))'];
       
        for reg=1:length(BoundaryLevelsLast)
            driftregimes(reg,:)=[diag(Rz(:,:,reg))'];
            driftbound(reg+1,:)=[diag(Rzbounds(:,:,reg+1))'];
        end 
        
elseif SOLUTION==3
         if CMEorErlang==1 
            n=OrderofPH;
            [MESystem,cv]= CMEParameterCalculator(n,tau);
            S=MESystem.A;
            S0=MESystem.b;
            alpha=MESystem.c;
         else
            alpha=zeros(1,OrderofPH);
            alpha(1)=1;
            ss=-ones(1,OrderofPH);
            S=diag(ss);
            for y=1:OrderofPH-1
               S(y,y+1)=1;
            end
            S=OrderofPH*S/tau;
            e=ones(length(S),1);
            S0=-S*e;
         end
         count=1;
        while BoundaryLevels(count)<b
            count=count+1;
        end        
        BoundaryLevelsLast=BoundaryLevels(1:count);
        BoundaryLevelsLast(count)=b;
        BoundaryLevelsLast(count+1: length(BoundaryLevels)+1)=BoundaryLevels(count: length(BoundaryLevels));
        Rytemp=Ryregimes;
        Ryboundstemp=Rybounds;
        Qytemp=Qy;
        Qyboundstemp=Qybounds;
        Rytemp(:,:,end+1)=Ry;
        Ryboundstemp(:,:,end+1)=Ry;
        Qytemp(:,:,count+1:end+1)=Qy(:,:,count:end);
        Qyboundstemp(:,:,count+2:end+1)=Qybounds(:,:,count+1:end);
        Qyboundstemp(:,:,count+1)=Qybounds(:,:,count+1);
        
         for reg=1:length(BoundaryLevelsLast)
            eTilde= ones(1,length(Qybounds(:,:,reg)));
            eTilde(end-MAPSIZE+1:end)=0;
            ITilde=diag(eTilde);

            QzTilde(:,:,reg)=[0 0  zeros(1,length(kron(alpha,kron(pi0,theta0))));0 0  zeros(1,length(kron(alpha,kron(pi0,theta0))));zeros(length(kron(S0,diag(ITilde))),1) kron(S0,diag(ITilde)) (kron(eye(length(S)),Qytemp(:,:,reg))+kron(S,ITilde))]; 
            QzboundsTilde(:,:,reg+1)=[0 0  zeros(1,length(kron(alpha,kron(pi0,theta0))));0 0  zeros(1,length(kron(alpha,kron(pi0,theta0))));zeros(length(kron(S0,diag(ITilde))),1) kron(S0,diag(ITilde)) (kron(eye(length(S)),Qyboundstemp(:,:,reg+1))+kron(S,ITilde))]; 
            
            RzTilde(:,:,reg)=diag([-1 -1 ones(1,size(kron(eye(length(S)),Rytemp(:,:,reg)),1))* kron(eye(length(S)),Rytemp(:,:,reg))]);
            RzboundsTilde(:,:,reg+1)=diag([-1 -1 ones(1,size(kron(eye(length(S)),Ryboundstemp(:,:,reg+1)),1))* kron(eye(length(S)),Ryboundstemp(:,:,reg+1))]);
        end
        QzboundsTilde(:,:,1)=[-1 1  zeros(1,length(kron(alpha,kron(pi0,theta0))));0 -1 kron(alpha,kron(pi0,theta0));zeros(length(kron(S0,diag(ITilde))),1) kron(S0,diag(ITilde)) (kron(eye(length(S)),Qybounds(:,:,1))+kron(S,ITilde))];
        RzboundsTilde(:,:,1)=diag([-1 -1 ones(1,size(kron(eye(length(S)),Ryboundstemp(:,:,1)),1))* kron(eye(length(S)),Ryboundstemp(:,:,1))]);

        Qregimes=QzTilde;
        Qbounds=QzboundsTilde; 
         
        driftbound(1,:)=[diag(RzboundsTilde(:,:,1))'];
       
        for reg=1:length(BoundaryLevelsLast)
            driftregimes(reg,:)=[diag(RzTilde(:,:,reg))'];
            driftbound(reg+1,:)=[diag(RzboundsTilde(:,:,reg+1))'];
        end 
        for reg=count+1:length(BoundaryLevelsLast)
            for d=1:OrderofPH
                index=2+(d-1)*size(Qy,2)+size(Qy,2)-2*MAPSIZE+1;
                Qbounds(index:index+MAPSIZE-1,1,reg)=Qbounds(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,reg)*ones(MAPSIZE,1);
                Qbounds(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,reg)=zeros(MAPSIZE,MAPSIZE);
                Qregimes(index:index+MAPSIZE-1,1,reg)=Qregimes(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,reg)*ones(MAPSIZE,1);
                Qregimes(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,reg)=zeros(MAPSIZE,MAPSIZE);
            end
        end
       Qbounds(index:index+MAPSIZE-1,1,end)=Qbounds(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,end)*ones(MAPSIZE,1);
       Qbounds(index:index+MAPSIZE-1,index+MAPSIZE:index+2*MAPSIZE-1,end)=zeros(MAPSIZE,MAPSIZE);
       
end


%MRMFQ solver
[coefficients,boundaries,Lzeromulti,Lnegmulti,Lposmulti,Anegmulti,Aposmulti]=MRMFQSolver(Qregimes,Qbounds,driftregimes,driftbound,BoundaryLevelsLast);
Qregimes=[]; 
Qbounds=[];  
driftregimes =[]; 
driftbound=[];
QzTilde=[];
RzTilde=[];
QzboundsTilde=[];
RzboundsTilde=[];
Qz=[];
Rz=[];
Qzbounds=[];
Rzbounds=[];
B=BoundaryLevelsLast;

%OUTPUT CALCULATION_______________________________________________________
%/////////////////////////////////////////////////////////////////
%Based on the chosen SOLUTION get the results from the outputs of MRMFQ solver
if SOLUTION==1 
    zeromass=boundaries{1};
    integral(1,:) =coefficients{1}*[Lzeromulti{1}*(B(1)-0); Anegmulti{1}\(expm(Anegmulti{1}*((B(1)-0)))-eye(size(Anegmulti{1},1)))*Lnegmulti{1};(Aposmulti{1})\(eye(size(Aposmulti{1},1))-expm((-Aposmulti{1})*(B(1)-0)))*Lposmulti{1}];
    waitintegral(1,:) =(B(1)/2)*(1-(ga(2)))*coefficients{1}*[Lzeromulti{1}*(B(1)-0); Anegmulti{1}\(expm(Anegmulti{1}*((B(1)-0)))-eye(size(Anegmulti{1},1)))*Lnegmulti{1};(Aposmulti{1})\(eye(size(Aposmulti{1},1))-expm((-Aposmulti{1})*(B(1)-0)))*Lposmulti{1}];
    abandonintegral(1,:) =(ga(2))*coefficients{1}*[Lzeromulti{1}*(B(1)-0); Anegmulti{1}\(expm(Anegmulti{1}*((B(1)-0)))-eye(size(Anegmulti{1},1)))*Lnegmulti{1};(Aposmulti{1})\(eye(size(Aposmulti{1},1))-expm((-Aposmulti{1})*(B(1)-0)))*Lposmulti{1}];
    successfulintegral(1,:) =(1-(ga(2)))*coefficients{1}*[Lzeromulti{1}*(B(1)-0); Anegmulti{1}\(expm(Anegmulti{1}*((B(1)-0)))-eye(size(Anegmulti{1},1)))*Lnegmulti{1};(Aposmulti{1})\(eye(size(Aposmulti{1},1))-expm((-Aposmulti{1})*(B(1)-0)))*Lposmulti{1}];
    for d=2:QUANTIZATION
          integral(d,:) =coefficients{d}*[Lzeromulti{d}*(B(d)-B(d-1)); Anegmulti{d}\(expm(Anegmulti{d}*((B(d)-B(d-1))))-eye(size(Anegmulti{d},1)))*Lnegmulti{d};(Aposmulti{d})\(eye(size(Aposmulti{d},1))-expm((-Aposmulti{d})*(B(d)-B(d-1))))*Lposmulti{d}];
          waitintegral(d,:) =((B(d)+B(d-1))/2)*(1-(ga(d+1)))*coefficients{d}*[Lzeromulti{d}*(B(d)-B(d-1)); Anegmulti{d}\(expm(Anegmulti{d}*((B(d)-B(d-1))))-eye(size(Anegmulti{d},1)))*Lnegmulti{d};(Aposmulti{d})\(eye(size(Aposmulti{d},1))-expm((-Aposmulti{d})*(B(d)-B(d-1))))*Lposmulti{d}];
          abandonintegral(d,:) =(ga(d+1))*coefficients{d}*[Lzeromulti{d}*(B(d)-B(d-1)); Anegmulti{d}\(expm(Anegmulti{d}*((B(d)-B(d-1))))-eye(size(Anegmulti{d},1)))*Lnegmulti{d};(Aposmulti{d})\(eye(size(Aposmulti{d},1))-expm((-Aposmulti{d})*(B(d)-B(d-1))))*Lposmulti{d}];
          successfulintegral(d,:) =(1-(ga(d+1)))*coefficients{d}*[Lzeromulti{d}*(B(d)-B(d-1)); Anegmulti{d}\(expm(Anegmulti{d}*((B(d)-B(d-1))))-eye(size(Anegmulti{d},1)))*Lnegmulti{d};(Aposmulti{d})\(eye(size(Aposmulti{d},1))-expm((-Aposmulti{d})*(B(d)-B(d-1))))*Lposmulti{d}];
    end 
    
    x=[0.1 0.2];
    ArrivalsWaitLessThanX=zeros(length(x),size(integral,2));
       for d=1:length(x)
            count=1;
            while B(count)<x(d)
                    count=count+1;
            end
            if count>1
                int=(1-(ga(count+1)))*coefficients{count}*[Lzeromulti{count}*(x(d)-B(count-1)); Anegmulti{count}\(expm(Anegmulti{count}*(x(d)-B(count-1)))-eye(size(Anegmulti{count},1)))*Lnegmulti{count};(Aposmulti{count})\(expm((-Aposmulti{count})*(B(count)-x(d)))-expm((-Aposmulti{count})*(B(count)-B(count-1))))*Lposmulti{count}];
            else
                int=(1-(ga(count+1)))*coefficients{count}*[Lzeromulti{count}*(x(d)-0); Anegmulti{count}\(expm(Anegmulti{count}*(x(d)-0))-eye(size(Anegmulti{count},1)))*Lnegmulti{count};(Aposmulti{count})\(expm((-Aposmulti{count})*(B(count)-x(d)))-expm((-Aposmulti{count})*(B(count)-0)))*Lposmulti{count}];
            end
            for l=1:count
                if l>1
                    ArrivalsWaitLessThanX(d,:)=ArrivalsWaitLessThanX(d,:)+successfulintegral(l-1,:);
                end
            end
            ArrivalsWaitLessThanX(d,:)=ArrivalsWaitLessThanX(d,:)+int;
       end

        for r=1:length(zeromass)
         ArrivalsWaitLessThanXmapped(:,r)=ArrivalsWaitLessThanX(:,r)*lmap(rem(r-1,2)+1);
         IntegralMapped(:,r)=integral(:,r)*lmap(rem(r-1,2)+1); 
         AbandonIntegralMapped(:,r)=abandonintegral(:,r)*lmap(rem(r-1,2)+1);
         WaitIntegralMapped(:,r)=waitintegral(:,r)*lmap(rem(r-1,2)+1);
         ZeroMassMapped(r)=zeromass(r)*lmap(rem(r-1,2)+1); 
        end
         normalization=SERVERSIZE*MAPSIZE;
        
        w0=sum(ZeroMassMapped)/(sum(sum(IntegralMapped(:,1:normalization)))+sum(ZeroMassMapped));       
        AbandonProb=sum(sum(AbandonIntegralMapped(:,1:normalization)))/(sum(ZeroMassMapped)+sum(sum(IntegralMapped(:,1:normalization))));
        w0s=sum(ZeroMassMapped)/((sum(sum(IntegralMapped(:,1:normalization)))+sum(ZeroMassMapped))*(1-AbandonProb));
        ProbArrivalsWaitLessThan01=sum(ArrivalsWaitLessThanXmapped(1,1:normalization))/((sum(sum(IntegralMapped(:,1:normalization)))+sum(ZeroMassMapped))*(1-AbandonProb)*(1-w0s));
        ProbArrivalsWaitLessThan02=sum(ArrivalsWaitLessThanXmapped(2,1:normalization))/((sum(sum(IntegralMapped(:,1:normalization)))+sum(ZeroMassMapped))*(1-AbandonProb)*(1-w0s));
        ExpectedWait=sum(sum(WaitIntegralMapped(:,1:normalization)))/((sum(ZeroMassMapped)+sum(sum(IntegralMapped(:,1:normalization))))*(1-AbandonProb));

       VarianceIntegral(1,:) =((B(1)/2)-ExpectedWait)^2*(1-(ga(2)))*coefficients{1}*[Lzeromulti{1}*(B(1)-0); Anegmulti{1}\(expm(Anegmulti{1}*((B(1)-0)))-eye(size(Anegmulti{1},1)))*Lnegmulti{1};(Aposmulti{1})\(eye(size(Aposmulti{1},1))-expm((-Aposmulti{1})*(B(1)-0)))*Lposmulti{1}];
       for d=2:QUANTIZATION
            VarianceIntegral(d,:) =(((B(d)+B(d-1))/2)-ExpectedWait)^2*(1-(ga(d+1)))*coefficients{d}*[Lzeromulti{d}*(B(d)-B(d-1)); Anegmulti{d}\(expm(Anegmulti{d}*((B(d)-B(d-1))))-eye(size(Anegmulti{d},1)))*Lnegmulti{d};(Aposmulti{d})\(eye(size(Aposmulti{d},1))-expm((-Aposmulti{d})*(B(d)-B(d-1))))*Lposmulti{d}];
       end 
       for r=1:length(zeromass)
           VarianceIntegralMapped(:,r)=VarianceIntegral(:,r)*lmap(rem(r-1,2)+1);
       end
       Variance=(sum(sum(VarianceIntegralMapped(:,1:normalization)))+sum(ZeroMassMapped)*ExpectedWait^2)/((sum(ZeroMassMapped)+sum(sum(IntegralMapped(:,1:normalization))))*(1-AbandonProb));
       
       steadyStateResult=[w0 w0s AbandonProb ExpectedWait Variance ProbArrivalsWaitLessThan01  ProbArrivalsWaitLessThan02];
       
       disp('    PR{W=0}  PR{W=0|S}   PR{A}    E{W|S}   Var{W|S} Fw|s,w>0(.1)Fw|s,w>0(.2)')
      disp([w0 w0s AbandonProb ExpectedWait Variance ProbArrivalsWaitLessThan01  ProbArrivalsWaitLessThan02])

        integral=[]; 
        integrr=[]; 
        integraa=[];
        integrall=[];
        integrass=[];
        bound11=[];
        denom=[];
        nomin=[];
        sumnomin=[];
        
elseif SOLUTION==2
        ca=boundaries{1};
        cb=boundaries{end};
        FirstPassageTimeofVirtual=(sum(cb))/ca(1) 
        
elseif SOLUTION==3
        ca=boundaries{1};
        FirstPassageTimeofActual=ca(1)/ca(2) 
end
