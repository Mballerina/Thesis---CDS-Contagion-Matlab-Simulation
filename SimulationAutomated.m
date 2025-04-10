%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    4. AUTOMATED SIMULATION FILE:                        %
%                  MAIN FILE FOR SIMULATION OF MODEL                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutine Explanation
% Program 3 calls this subroutine as a function giving all inputs
% automatically and calling as many runs as needed. The program gives all
% calculated model outputs for one run and one set of inputs.

function Values=SimulationAutomated(I,K,MT,ML,D_k,CCPFail,Pred,Dist,Collusion,CCPTrade,DISTTrade)

Space=' ';
Intro='WELCOME TO THE CCP LIQUIDATION AND RECOVERY SIMULATION';
Underbar='---------------------------------------------------------';
disp(Space)
disp(Underbar)
disp(Intro)
disp(Underbar)
disp(Space)

%____________________________________________________________
%%%%%%%%%%%%%%%%%%%MODEL SETUP%%%%%%%%%%%%%%%%%%%
%____________________________________________________________
%Set random seed for repeatability
rng(3) %goes back to initial value, reseeds generator
J=I; %Counterparties - must be square matrix
%I=14; %Number of Banks
%K=100; %Span of assets

%SET UP THE LOOPS FOR THE PERIODS
%setup periods 0,1,2,3
T=3;
%the timestep is 1 day
tao=1;
%liquidation period timesteps are 0,1,2,3,4,5
L=6;

%% Time 0 - Setup
%START WITH HOMOGENOUS CASE W/ HETEROGENEOUS POSITIONS
%t=0
%____________________________________________________________
%%%%%%%%%%%%%%%%%%%TIME 0 - NO LIQUIDATION%%%%%%%%%%%%%%%%%%%
%____________________________________________________________
%Initialize matrix of endowments for banks C_i(0)
endow=10e10;
%C_i0=round(endow*rand(i))
C_i0=endow.*ones(I,1); %kind of inconsequential

%The value and matrix of outside assets Q
Q=1.1e10;
%Q_i=round(Q*rand(i))
Q_i=Q.*ones(I,1);
%Need a hybrid and pure version for looping
HQ_i=Q_i;
PQ_i=Q_i;

%The cash endowment of banks (starts with cash=endowment)
Gamma=endow;
Gamma_i=Gamma.*ones(I,1);

%The cash endowment of the CCP
Gamma_0=5e9;

%The defaulted banks' liability and surplus held by CCP
%Hybrid:
EDHL_0=0;
EDHA_0=0;

%Pure:
EDPL_0=0;
EDPA_0=0;

HCCPLoss=0;
PCCPLoss=0;
THC_0=0;
TPC_0=0;
THC_i=0;
TPC_i=0;
flag=0;
lell=0;
hell=0;
pell=0;
bell=0;
bhell=0;
bpell=0;
sell=0;
LSurplus=0;
BSurplus=0;

%CALIBRATION:
%Nominal worth from 2010 BIS data
%---------------------------------
N=19e12; %Need full market value to get sufficient a and b values *(test)
ND=19e12.*.85;  %Nominal value of top 14 dealers
N_i=ND/I %Nominal value of each bank
N_ij=N_i/(I-1); %Nominal value across counterparties
N_ijk=N_ij/(K); %Nominal value distr for each cp and each asset (but this is yearly value)

%Calibration from Duffie and from Oehmke data
%-----------------------------------------------
%Value of nominal holdings in USD (Should be in billions)
%Get an integer matrix based on normal distribution
%matrix must be negative reflections of triangle with zero diagonal
grosscds=13.3e10/I;
m=grosscds/K;
volgrosscds=12.19e10/(I)
s=volgrosscds/(K);
pd = makedist('Normal','mu', m/(I-1),'sigma',s/(I-1));
for k=1:K
    X=round(random(pd,I,J,K));
    X1(:,:,k)=triu(X(:,:,k));
    X2(:,:,k)=permute(-X1(:,:,k),[2,1]);
    X_ij(:,:,k)=(X1(:,:,k)+X2(:,:,k));
end

%Some relationships must also randomly be zero
%Set the saturation of of the matrix randomly
sat=rand(K,1); %random level of saturation for each asset
saturation=round(sat*I); %# pairs not holding

%Saturation function randomly knocks out bank pairs
if saturation>0
    %Pair (b1,b2) and (b2,b1) must be zero
    for k=1:K
        for s=1:saturation(k)
            b1=1;
            b2=1;
            %Choose bank to get a zero in the holdings randomly
            while b1==b2 | b1==0 | b2==0
                b1=round(I*rand(1));
                b2=round(J*rand(1));
            end
            X_ij(b1,b2,k)=0;
            X_ij(b2,b1,k)=0;
        end
    end
end

%set X_ji to be the opposite side of positions
X_ji=permute(X_ij,[2,1,3]);

%Get an initial value for size of holdings;
X_0=X_ij;

%Test of sufficient holding
%---------------------------
Xdealer_i=sum(abs(X_ij),2);
Xdealer_k=sum(abs(Xdealer_i),3);
Xmarket_k=sum(Xdealer_i); % gives max availability of each asset k in market
Xmarket=sum(Xdealer_k);

disp('The Notional Holding of the Market in the model is:')
disp(Xmarket)
disp('The Notional Holding of the Market from data is:')
disp(N)


%CALIBRATION
%--------------
%Need an original value of positions at time zero
%Parties enter into contract at zero and immediately price moves (bps)
%In normal times the bps can swing from 0 to 50 (20/1000) bps and in turmoil from 50
%500 bps= 0.5% of notional
%Set the Fundamental price randomly

%Inialize a random price matrix for each asset
rng(5)
%We set our values going into the crisis range which has avg of 250 bps.
%and a max of about 550 bps with a min of about 50 bps.
pd = makedist('Normal','mu',249.2/1000,'sigma',269/1000);
for l=1:L
    S=random(pd,1,1,K,L);
    S_l=repmat(S,I,J);
end
S0=S_l(I,1,:,1);

%Set the guarantee fund contribution according to holding
%Initial margin is currently 50% of holdings by Treasury
%10-30 % in bilateral CDS
g=11.2e9;
n=N/I;
margin=g/n*100; %for I=14 banks this is a margin of approx. 0.0954 % of holding.
g_i=initmargin(margin,l,S_l,X_ij); %total guarantee for each i
%Note: this should be just taken from original endowment.
%Needed to do loops over the periods
hg_i=g_i;
pg_i=g_i;

%The default fund should be about 1/10th=10% of the default fund.
d_i=0.10.*g_i;
hd_i=d_i;
pd_i=d_i;

%Total cash holdings of banks
HGamma_i=Gamma_i; %Since accounted for in Z
PGamma_i=Gamma_i;

%Initial liquidation matrix must be zero
a_ij0=zeros(I,J,K);
a_ji0=-a_ij0;

%Time 1 matrix must have one exogenous default with four possible
%Bank 4 defaults first
dnum=1;
Dflt=[0];

%Set exogenous default
Dflt(dnum)=4;
D0=Dflt;

%Set depletion for defaulted party
hg_i(Dflt(1))=0;
HGamma_i(Dflt(1))=0;
HQ_i(Dflt(1))=0;
pg_i(Dflt(1))=0;
PGamma_i(Dflt(1))=0;
PQ_i(Dflt(1))=0;

%Set future vector of defaults.
Hnum=1;
Pnum=1; %index number of default vector
HDefault=[Dflt(1)];
HDefaultdummy=zeros(I,1,L);
HDefaultdummy(Dflt(1),1,1)=Dflt(1);
PDefault=[Dflt(1)];
PDefaultdummy=zeros(I,1,L);
PDefaultdummy(Dflt(1),1,1)=Dflt(1);

%--------------------------------------------------------------
%% Time 1 - Liquidation Period
%t=1
%___________________________________________________________
%%%%%%%%%%%%TIME 1 - REALIZATION OF LIABILITIES%%%%%%%%%%%%%%
%                        LIQUIDATION PERIOD
%____________________________________________________________
%get an integer matrix based on normal distribution
rng(10);
%a=(-1).*sqrt(N_ijk); %might need more room
%a=(-1).*(N_ijk./10)
%b=sqrt(N_ijk);
%b=(N_ijk./10) % calibrated to market trading level (see word doc)
m=0.516e9/I;
m=m/K;
s=0.646e9/(I);
s=s/K;
pd = makedist('Normal','mu', m/(I-1),'sigma',s/(I-1));
a_ij=round(random(pd,I,J,K)); %next round it should be the one from previous round.
for k=1:K
    A(:,:,k)=eye(I).*a_ij(:,:,k);
    a_ij(:,:,k)=-A(:,:,k) + a_ij(:,:,k);
end
%Make sure that cannot liquidation more than have in posession
Dummy1=(abs(a_ij)>abs(X_ij));
%Dummy2=(-1).*(X_ij<0)+(X_ij>0); %Gets the direction of holdings
%Dummy=(-1).*(Dummy2).*(Dummy1); %Gives you the opposite direction for thos holdings
V1=(-1).*X_ij.*Dummy1; %gives the opposite direction for liquidation
V2=a_ij-a_ij.*Dummy1; %removes the current entries
a_ij=V2+V1;


%LIQUIDATION DECISION FOR PREDATORS AND DISTRESSED.
%----------------------------------------------------------------
%We will have distressed always liquidate their side for now.
%We will have predators always predate (for now)
S_l1=0;
S_l2=0;

%Start the loop for periods
%If you want loop for more periods you can do L=20 then see break condition
for l=2:L
    
    Space=' ';
    disp(Space)
    Period='THE LIQUIDATION TIMESTEP IS ';
    period=l-1;
    Underbar='---------------------------------------------------------';
    disp(Period), disp(period)
    disp(Underbar)
    disp(Space)
    
    %CALCULATE THE PREDATORS
    %------------------------
    %Any entity related directly to default is distressed
    Distress_k=distress(Dist,Dflt,X_ij,I,J,K);
    %Any entity not connected to default can be predator
    Predator_k=predator(Dist,Pred,Distress_k,Dflt,X_ij,I,J,K);
    %Note:entities can be in both groups for an asset holding.
    
    %choose largest predator with highest holding and most income
    LPred=lpred(Predator_k,HGamma_i,X_ij,I,J,K);
    %--------------------------------------------------------------------
    
    %LIQUIDATION MATRIX ASSEMBLED
    %-----------------------------
    %can change defaults by changing random seed
    a_ij=liquidmatrix(N_ijk,Distress_k,Dist,Predator_k,Pred,LPred,Collusion,a_ij,Dflt,X_ij,I,J,K); %calls function
    a_ji=permute(a_ij,[2,1,3]);
    
    %TO GET THE PRICE
    %-----------------
    if l>2
        S_l2=S_l1;
    end
    
    %To get the movement in the spread alone
    Spread=spread(l,D_k,Dflt,S_l,X_ij,X_ji,a_ij,a_ji,I,J,K);
    S_ij=Spread(1,1,:);
    S_ji=Spread(2,1,:);
    S_l1=repmat(S_ij,I,J);
    
    if l>1
        S_l1Plot(1,1,:,l)=S_l1(1,1,:);
    end
    
    %THE PRICING FUNCTIONAL
    %---------------------------------------------------------------------
    %Collect of each i, to create either a liability or recievable
    %By construction:liability is negative and a recievable is positive
    Pricetwo=pricetwo(l,D_k,Dflt,S_l,S_l1,S_l2,X_ij,X_ji,a_ij,a_ji,I,J,K);
    P_IJ=Pricetwo(1,:);
    P_JI=Pricetwo(2,:);
    P=cat(1,P_JI,-P_IJ);
    
    %Market Liquidity needs to be adapted (decrease) each liquidation time-step
    D_k=D_k-MT;
    
    %We get the Net Exposure
    Lambda_i=sum(P); %over banks i for all exposures
    dummyLambda=1-Lambda_i.*0;
    NetLambdaPlus_i=max(Lambda_i,zeros(size(dummyLambda))); %Surplus
    NetLambdaMinus_i=abs(min(Lambda_i,zeros(size(dummyLambda)))); %Liability
    %____________________________________________________________
    
    %Time 1 BUDGET OF BANKS:
    %-------------------------
    %Different ones for each l - timestep
    %Hybrid Fund
    hA_i=max(HGamma_i + transpose(P_JI) + HQ_i, zeros(I,1));
    hL_i=max(transpose(P_IJ) + abs(min(HGamma_i + transpose(Lambda_i), zeros(I,1))), zeros(I,1));
    %The nominal wealthis (Gamma_i + transpose(Lambda_i))
    %Pure Fund
    pA_i=max(PGamma_i + transpose(P_JI) + PQ_i, zeros(I,1));
    pL_i=max(transpose(P_IJ) + abs(min(PGamma_i + transpose(Lambda_i), zeros(I,1))), zeros(I,1)); %can't have negative liability
    
    %BANK LIABILITY TO CCP
    %------------------------
    %Set the Liability to the CCP
    HL_i0=max(transpose(NetLambdaMinus_i)-hg_i ,zeros(I,1)); %Hybrid Fund
    PL_i0=transpose(NetLambdaMinus_i); %Pure Fund
    
    %LIQUIDATION FRACTION OF BANKS
    %------------------------------
    %Assume Fraction recovered in liquidation is 10%
    R=.10*Q;
    R_i=0.10*Q.*ones(I,1);
    R_i(Dflt(:),1)=0; % A defaulted bank with no assets has no recovery value.
    
    %LIQUIDATION FRACTION FOR BANKS
    %-------------------------------
    %For Hybrid fund:
    Dummy=(HQ_i~=0); %want to pick out which banks have no external asset
    HZ_i=min(abs(min(HGamma_i-hg_i-HL_i0,zeros(I,1))./R), ones(I,1));
    %Note margin is used and available and so not subtracted from cash
    
    %For Pure fund: assume the guarantee comes from cash assets
    Dummy=(PQ_i~=0); %want to pick out which banks have no external asset
    PZ_i=min(abs(min(PGamma_i-pg_i-transpose(NetLambdaMinus_i),zeros(I,1))./R),ones(I,1));
    %Note: margin is subtracted from the available cash/cash is earmarked for margin
    
    %GUARANTEE FUND CONTRIBUTION OF BANKS:
    %-------------------------------------
    %In Hybrid Fund:
    HG_i=max(transpose(Lambda_i)+hg_i, zeros(I,1))-transpose(NetLambdaPlus_i);
    HG_tot=sum(HG_i);
    
    %In Pure Fund:
    PG_i=max(transpose(Lambda_i)+PGamma_i+R_i, zeros(I,1))-max(transpose(Lambda_i)+PGamma_i+R_i-pg_i,zeros(I,1));
    PG_tot=sum(PG_i);
    
    %DEFAULT FUND CONTRIBUTION OF BANKS:
    %-------------------------------------
    %In Hybrid Fund:
    HD_i=max(transpose(Lambda_i)+hg_i+hd_i, zeros(I,1))-max(transpose(Lambda_i)+hg_i,zeros(I,1));
    HD_tot=sum(HD_i);
    
    %In Pure Fund:
    PD_i=max(transpose(Lambda_i)+PGamma_i+R_i-pg_i, zeros(I,1))-max(transpose(Lambda_i)+PGamma_i+R_i-pg_i-pd_i,zeros(I,1));
    PD_tot=sum(PD_i);
    %_________________________________________________________________________
    
    %Time 1 BUDGET OF CCP
    %----------------------
    %Set Volume fee level for service
    f=0.05; % 5% max. fee for clearing
    
    %Set Liability of CCP
    L_0i=(1-f).*transpose(NetLambdaPlus_i);
    L_0=abs(sum(L_0i));
    
    %Hybrid Case:
    HA_0=Gamma_0 + sum(hg_i) + sum(HL_i0);
    HL_0=abs(min(sum(L_0i) + sum(HG_i)+(Gamma_0 + f.*sum(NetLambdaPlus_i)), 0));
    
    %Pure Case:
    PA_0= Gamma_0 + sum(pg_i) + sum(PL_i0);
    PL_0=abs(min(sum(L_0i) + sum(PG_i)+(Gamma_0 + f.*sum(NetLambdaPlus_i)), 0));
    %The nominal wealth is (Gamma_0 + f.*sum(NetLambdaPlus_i))
    %_____________________________________________________________________
    
    %EQUILIBRIUM PAYMENT OF BANK TO CCP
    %----------------------------------
    %Hybrid Case:
    EHL_i0=abs(min(HL_i0,(HGamma_i-hg_i+R_i)));
    %Pure Case:
    EPL_i0=abs(min(PL_i0,(PGamma_i-pg_i+R_i)+pg_i));
    
    %Shorfall to CCP
    HLMinus_i0=abs(HL_i0-EHL_i0);
    PLMinus_i0=abs(PL_i0-EPL_i0);
    %______________________________________________________________________
    
    %------------------------------------------------------
    %                TOTAL EQUIB HOLDINGS
    %                   END OF PERIOD
    %------------------------------------------------------
    %Time 1 EQUILIBRIUM BUDGET OF CCP
    %--------------------------------
    %Hybrid Case:
    EHA_0= Gamma_0 + sum(hg_i) + sum(EHL_i0);
    EHL_0=abs(min(EHA_0,HL_0));
    
    %Pure Case:
    EPA_0= Gamma_0 + sum(pg_i) + sum(EPL_i0);
    EPL_0=abs(min(EPA_0,PL_0));
    
    %Proportionality Rule for Disbursements with Limited Liability:
    if sum(L_0i)>0
        HPI_0i=max(L_0i/sum(L_0i),zeros(I,1));
        PPI_0i=max(L_0i/sum(L_0i),zeros(I,1));
    else
        HPI_0i=zeros(I,1);
        PPI_0i=zeros(I,1);
    end
    
    %THE TOTAL GUARANTEE/DEFAULT FUND HELD BY THE CCP
    %---------------------------------------
    %Hybrid Fund
    EHG_tot=min(HG_tot,max(HA_0+EHL_0-Gamma_0-f.*sum(NetLambdaPlus_i),0));
    
    
    %Pure Fund
    %EPG_tot=sum(EPG_i)
    EPG_tot=min(PG_tot,max(PA_0+EPL_0-Gamma_0-f.*sum(NetLambdaPlus_i),0));
    
    %THE PORTION OF GUARANTEE/DEFAULT FUND RETURN TO BANKS
    %-----------------------------------------------
    %Hybrid Fund
    if HG_tot~=0
        EHG_i=(HG_i./HG_tot).*EHG_tot;
    else
        EHG_i=0.*HG_i;
    end
    EHG_tot2=sum(EHG_i);
    
    %Pure Funds
    if PG_tot~=0
        EPG_i=(PG_i./PG_tot).*EPG_tot;
    else
        EPG_i=0.*PG_i;
    end
    EPG_tot2=sum(EPG_i);
    
    
    %THE ACCOUNT AND SHORTFALL OF CCP
    %----------------------------------
    %Need to determine shortfall of CCP into the next round.
    %From default condition and Seizure we need a CCP liquidation decision.
    %Terminal Net Worth
    HC_0=(EHA_0)-HL_0-EHG_tot;
    PC_0=(EPA_0)-PL_0-EPG_tot;
    
    %Terminal Shortfall
    HCMinus_0=HL_0-(EHL_0);  %won't add loss from previous period until the end.
    PCMinus_0=PL_0-(EPL_0);
    
    %THE ACCOUNT AND SHORTFALL OF BANKS
    %-----------------------------------
    %Hybrid Terminal Assets
    HA_i=HGamma_i+HZ_i.*R_i+(1-HZ_i).*HQ_i+EHL_0.*HPI_0i+EHG_i-hg_i;
    HC_i=HA_i-HL_i0;
    HCMinus_i=HL_i0-EHL_i0;
    HCPlus_i=HGamma_i+HZ_i.*R_i+(1-HZ_i).*HQ_i;
    
    %Pure Terminal Assets
    PA_i=PGamma_i+PZ_i.*R_i+(1-PZ_i).*PQ_i+EPL_0.*PPI_0i+EPG_i-pg_i; %problem
    PC_i=PA_i-PL_i0;
    PCMinus_i=PL_i0-EPL_i0;
    PCPlus_i=PGamma_i+PZ_i.*R_i+(1-PZ_i).*PQ_i;
    
    %In the aggregate
    HCPlus_tot=sum(HCPlus_i);
    PCPlus_tot=sum(PCPlus_i);
    Surplus=HCPlus_tot-PCPlus_tot;
    LSurplus=LSurplus+Surplus;
    
    Text3='The difference in the surplus of banks due to liquidations in Hybrid vs. Pure Fund:';
    disp(Text3)
    disp(Surplus)
    %_________________________________________________________________________
    
    
    %--------------------------------------------------------------------
    %   !!!!!!!!!!!!!!!!!!! DEFAULT AND SEIZURE !!!!!!!!!!!!!!!!!!!
    %--------------------------------------------------------------------
    %DETERMINE THE SHORTFALL THE CCP TAKES ON FOR NEXT ROUND
    
    %THE BANKS DEFAULT CONDITION
    %-----------------------------------------------
    %In Hybrid Case
    for i=1:I
        if  HL_i0(i)>=(HGamma_i(i)+R_i(i))
            Hnum=Hnum+1;
            HDefaultdummy(i,1,l)=i; %fills by default per time step
            HDefault=sum(HDefaultdummy,3); %sums defaults
            numdummy=[1:I]; %sets up a index matrix of banks
            numdummy=transpose(numdummy); %makes it into column
            for j=1:length(numdummy)
                if HDefault(j)~=numdummy(j) %if there are multiples
                    if HDefault(j)~=0 %checks that nonzero
                        HDefault(j)=numdummy(j); %makes sure that entry=index
                    end
                end
            end
        end
        
        %In Pure Case
        if PL_i0(i)>=PGamma_i(i)+R_i(i)+pg_i(i)
            Pnum=Pnum+1;
            PDefaultdummy(i,1,l)=i;
            PDefault=sum(PDefaultdummy,3);
            numdummy=[1:I];
            numdummy=transpose(numdummy);
            for j=1:length(numdummy)
                if PDefault(j)~=numdummy(j)
                    if PDefault(j)~=0
                        PDefault(j)=numdummy(j);
                    end
                end
            end
        end
    end
    
    disp('Is the Hybrid Default the same as the Pure Default')
    disp(HDefault')
    disp(PDefault')
    
    %% THE CCP SEIZURE RULE
    %-----------------------
    %allows default fund updating unambiguously
    HSeizure=(HL_i0-HGamma_i-R_i>=hg_i);
    PSeizure=(PL_i0-PGamma_i-R_i>=pg_i);
    dnum=1;
    Default=[D0];
    for i=1:I
        if HSeizure(i,1)==1 & PSeizure(i,1)==1;
            dnum=dnum+1;
            Default(dnum,1)=i;
            %Dflt(dnum,1)=i;
        end
    end
    
    LDefault=Default;
    %---------------------------------------------------------------
    
    
    %THE LIQUIDATION DECISION OF BANKS FOR NEXT ROUND
    %------------------------------------------------
    %Can determine the new X matrix with liquidation
    X_ij=X_ij+a_ij; %each round we replace the matrix
    X_ji=X_ji+a_ji;
    
    %END OF DAY MARGIN REFILL
    %-------------------------------------------------
    %HYBRID:
    %------------
    %%Reduction of amount from g_i which was taken
    %hg_i=HG_i Check on amount
    
    %%Need to reduce the amount initial margin fund (if used) in Hybrid Case
    hgMin=HL_i0; %remaining liability not covered by initial margin
    hg_i=max(hg_i-transpose(NetLambdaMinus_i),zeros(I,1)); %amount of margin used and left
    
    %%Need to reduce cash by amount that still owed to CCP
    HGammaMin=abs(min(HGamma_i-hgMin, zeros(I,1)));
    HGamma_i=max(HGamma_i-hgMin,zeros(I,1)); %amount of cash reduced by that not covered by guarantee
    %HGammaMin=min(HGamma_i-HL_i0, zeros(I,1)); %amount of liability cash cannot cover
    
    %%Need to reduce the external asset by the amount that must be liquidated
    HQMin_i=abs(min(HQ_i-HGammaMin,zeros(I,1)));%amount the external fund cannot cover should equal CMin_i
    HQ_i=HQ_i.*(1-HZ_i); %fraction needed needed to be liquidated
    %HQ_i=max(HQ_i-HGammaMin, zeros(I,1)) should match above
    
    %%Any leftover liability goes to the CCP in an DL_i0=DL_0 (line #352)
    HDL_i0=abs(HQMin_i); %amt transferred to CCP from liabilities
    %%Those that don't have enough to cover margin call are in Default
    for i=1:length(HDefault)
        if transpose(HDL_i0(i))~=0
            %if transpose(HCMinus_i(i))~=0
            HDefault(i)=i;
        end
    end
    
    %Left over default liability not covered by bank must go to the CCP
    EDHL_0=EDHL_0+sum(HDL_i0);
    
    %Any surplus of Defaulted banks must go toward min their liability
    dummy=(HDefault~=0); %non-zero defaulted banks
    EDHA_0=EDHA_0+sum(HCPlus_i.*dummy);
    
    %PURE:
    %-----------
    %%Need to reduce amount in initial margin if used in Pure Case:
    %pg_i=PG_i Check on amount
    
    %%Need to reduce cash by amount that still owed to CCP
    PGammaMin=abs(min(PGamma_i-PL_i0,zeros(I,1))); %Amount lacking to pay CCP in cash
    PGamma_i=max(PGamma_i-PL_i0,zeros(I,1)); %Amount left in account
    
    %%Need to reduce the external asset by the amount liquidated
    PQMin_i=abs(min(PQ_i-PGammaMin,zeros(I,1))); %amount the external fund cannot cover
    PQ_i=PQ_i.*(1-PZ_i); %fraction needed needed to be liquidated
    
    %%If cash&external asset not enough for liability then use initial margin
    pg_iMin=abs(min(pg_i-PQMin_i,zeros(I,1))); %amt margin can't cover
    pg_i=max(pg_i-PQMin_i,zeros(I,1)); %amt removed from initial margin
    
    %%Any leftover liability goes to the CCP in an DL_i0=DL_0 (line #352)
    PDL_i0=abs(pg_iMin); %amt transferred to CCP from liabilities
    %%Those that don't have enough to cover margin call are in Default
    for i=1:length(PDefault)
        if transpose(PDL_i0(i))~=0
            % if transpose(PCMinus_i(i))~=0
            PDefault(i)=i;
        end
    end
    
    %Left over default liability not covered by bank must go to the CCP
    EDPL_0=EDPL_0+sum(PDL_i0); %should equal PCMinus_i
    
    %Any surplus of Defaulted banks must go toward min their liability
    Dummy=(PDefault~=0); %non-zero defaulted banks
    EDPA_0=EDPA_0+sum(PCPlus_i.*Dummy);
    
    %Liquidation loss incurred vs. assets made on liquidation only up to second
    HLiquidationGain=max(EDHA_0-EDHL_0,zeros(1,1)) %these are cumulative
    HLiquidationLoss=abs(min(EDHA_0-EDHL_0,zeros(1,1)))
    PLiquidationGain=max(EDPA_0-EDPL_0,zeros(1,1))
    PLiquidationLoss=abs(min(EDPA_0-EDPL_0,zeros(1,1)))
    HCCPLoss=min(HC_0+HLiquidationGain-HLiquidationLoss,0) %can't add previous CCP Loss because Liquidation is cumulative
    PCCPLoss=min(PC_0+PLiquidationGain-PLiquidationLoss,0)
    
    if HCCPLoss<0
        Flag='The Hybrid CCP has FAILED! At liquidation time-step: '
        hell=l-1;
        disp(hell)
        flag=L; % will end the looping after this round
    end
    if PCCPLoss<0
        Flag='The Pure CCP has FAILED! At liquidation time-step:'
        pell=l-1;
        disp(pell)
        flag=L; % will end the looping after this round
    end
    
    %for plotting we need separate values that don't cumulate
    LHLiquidationGain= HLiquidationGain;
    LHLiquidationLoss=HLiquidationLoss;
    LPLiquidationGain=PLiquidationGain;
    LPLiquidationLoss=PLiquidationLoss;
    LHCCPLoss=HCCPLoss;
    LPCCPLoss=PCCPLoss;
    %___________________________________________________________________
    
    %% Time 1 - Break Conditions
    if flag==L & CCPFail=='Y';
        disp('CCP Failed due to Insolvency in l step')
        lell=l-1;
        break
    end
end

%% Time 1 - Liquidation Defaults
ell=l+1;
disp('Defaulted this liquidation round ')
disp(ell)
disp(transpose(Default))

for i=1:length(HDefault)
    if HDefault(i)~=0
        Dflt(dnum)=HDefault(i);
        dnum=dnum+1;
    end
    
    disp('The liquidation evolution of default is: ')
    DfltEvo=Dflt;
    disp(DfltEvo)
    Dflt=sort(Dflt);
    Dflt=Dflt(cumsum(Dflt(:),2)>0);
    Dflt = unique(Dflt,'sorted')
    disp('The Dflt is: ')
    disp(Dflt)
    LDflt=Dflt; %for plotting
    LDflt=sum(Dflt>0);
end

%--------------------------------------------------------------------------
%% Time 1 - Terminating Values
%----------------End of Period 1/Liquidation Period time steps--------------

%THE TERMINATING VALUE OF THE BANK AND CCP FOR LIQUIDATION DECISIONS
%--------------------------------------------------------------------
THC_0=THC_0+Gamma_0+f*sum(transpose(NetLambdaPlus_i))-abs(min(HG_tot-sum(HCMinus_i),0))+HLiquidationGain-HLiquidationLoss; %Check this addition
TPC_0=TPC_0+Gamma_0+f*sum(transpose(NetLambdaPlus_i))-abs(min(PG_tot-sum(PCMinus_i),0))+PLiquidationGain-PLiquidationLoss;

DHG_tot=HG_tot;
if HG_tot==0
    DHG_tot=1;
end
THC_i=THC_i+HGamma_i+HQ_i+transpose(Lambda_i)-HPI_0i.*HCMinus_0-HZ_i.*(HQ_i-R_i)-f.*transpose(NetLambdaPlus_i)-(HG_i./DHG_tot).*(HG_tot-EHG_tot);

DPG_tot=PG_tot;
if PG_tot==0
    DPG_tot=1;
end
TPC_i=TPC_i+PGamma_i+PQ_i+transpose(Lambda_i)-PPI_0i.*PCMinus_0-PZ_i.*(PQ_i-R_i)-f.*transpose(NetLambdaPlus_i)-(PG_i./DPG_tot).*(PG_tot-EPG_tot);
%These values should be carried over to new period.

%For Plotting
LTHC_0=THC_0;
LTPC_0=TPC_0;
LTHC_i=sum(THC_i);
LTPC_i=sum(TPC_i);
LEHG_tot=EHG_tot2;
LEPG_tot=EPG_tot2;

%***********************************************************
%***********************************************************
%% Time=2 - BuyBack Period

%--------------------Start of new period--------------------
%t=2
%___________________________________________________________
%%%%%%%%%%%%%%%%%TIME 2 - BUY-BACK PERIOD%%%%%%%%%%%%%%%%%%%
%____________________________________________________________
if Pred~=0
    %Need to keep S value from previous case to serve as the zero value.
    L1=L; %number of last round in previous period
    SF1=S_l(1,1,:,L1); %gets last value prices
    
    %[Calibration]
    %--------------
    X00=sum(X_0,2); %sum of org. holding for each bank
    A_k=sum(abs(sum(a_ij,2)),1); %sum of all trading times
    tao=abs(X00.*I./repmat(A_k,[I,1])); %max time each bank can buyback
    Xmax=sum(sum(abs(X00),2),1)./I; %highest possible holding
    Xmin=min(abs(X00));
    if Pred==1 | Collusion=='Y'
        L=round(min(Xmax./A_k)); %gets largest liquidation time needed
    elseif Pred>1 & Collusion=='N'
        L=min(round(mean(tao))); %better approach based on the mean rate
    end
    
    if L<=2
        L=3;
    end
    
    disp('The number of buyback rounds needed is ')
    disp(L)
    
    %CALIBRATION
    %Need an original value of positions at time zero
    %Parties enter into contract at zero and immediately price moves (bps)
    %Set the Fundamental price randomly
    %Inialize a random price matrix for each asset
    rng(7)
    %Try with a distribution fit to data
    pd = makedist('Normal','mu',249.2/1000,'sigma',569/1000);
    for l=1:L
        S=random(pd,1,1,K,L);
        S(1,1,:,1)=SF1;
    end
    S_l=repmat(S,I,J);
end %only for pred=0

if Pred==0
    L=3
end

%Start the loop for periods
for l=2:L
    if Pred~=0
        
        Space=' ';
        disp(Space)
        Period='THE BUY-BACK TIME-STEP IS ';
        period=l-1;
        Underbar='---------------------------------------------------------';
        disp(Period), disp(period)
        disp(Underbar)
        disp(Space)
        
        %choose largest predator with highest holding and most income
        LPred=lpred(Predator_k,HGamma_i,X_ij,I,J,K);
        
        %LIQUIDATION MATRIX ASSEMBLED
        %---------------------------
        %can change defaults by changing random seed
        a_ij=buybackmatrix(CCPTrade,DISTTrade,N_ijk,X_0,Xmax,Distress_k,Dist,Predator_k,Pred,LPred,Collusion,a_ij,Dflt,X_ij,I,J,K); %calls function
        a_ji=permute(a_ij,[2,1,3]);
        
        %THE PRICING FUNCTIONAL
        %---------------------------------------------------------------------
        if l>2
            S_l2=S_l1;
        end
        
        %To get the movement in the spread alone
        Spread=spread(l,D_k,Dflt,S_l,X_ij,X_ji,a_ij,a_ji,I,J,K);
        S_ij=Spread(1,1,:);
        S_ji=Spread(2,1,:);
        S_l1=repmat(S_ij,I,J);
        
        if l>1
            l2=5+l;
            S_l1Plot(1,1,:,l2)=S_l1(1,1,:);
        end
        
        %Collect of each i, to create either a liability or recievable
        %By construction:liability is negative and a recievable is positive
        Pricetwo=pricetwo(l,D_k,Dflt,S_l,S_l1,S_l2,X_ij,X_ji,a_ij,a_ji,I,J,K);
        P_IJ=Pricetwo(1,:);
        P_JI=Pricetwo(2,:);
        P=cat(1,P_JI,-P_IJ);
        
        %Market Liquidity needs to be adapted (increase) each buyback time-step
        D_k=D_k+MT;
        
        %!Need separate for H and P but first figure out default entries!
        %We get the Net Exposure
        Lambda_i=sum(P); %over banks i for all exposures
        dummyLambda=1-Lambda_i.*0;
        NetLambdaPlus_i=max(Lambda_i,zeros(size(dummyLambda))); %Surplus
        NetLambdaMinus_i=abs(min(Lambda_i,zeros(size(dummyLambda)))); %Liability
        %_____________________________________________________________
        
        %Time 1 BUDGET OF BANKS:
        %-------------------------
        %Different ones for each l - timestep
        %Hybrid Fund
        hA_i=max(HGamma_i + transpose(P_JI) + HQ_i, zeros(I,1));
        hL_i=max(transpose(P_IJ) + abs(min(HGamma_i + transpose(Lambda_i), zeros(I,1))), zeros(I,1));
        %The nominal wealth is (Gamma_i + transpose(Lambda_i))
        %Pure Fund
        pA_i=max(PGamma_i + transpose(P_JI) + PQ_i, zeros(I,1));
        pL_i=max(transpose(P_IJ) + abs(min(PGamma_i + transpose(Lambda_i), zeros(I,1))), zeros(I,1)); %can't have negative liability
        
        %BANK LIABILITY TO CCP
        %------------------------
        %Set the Liability to the CCP
        HL_i0=max(transpose(NetLambdaMinus_i)-hg_i ,zeros(I,1)); %Hybrid Fund
        PL_i0=transpose(NetLambdaMinus_i); %Pure Fund
        
        %LIQUIDATION FRACTION OF BANKS
        %------------------------------
        %Assume Fraction recovered in liquidation is 10%
        R=.10*Q;
        R_i=0.10*Q.*ones(I,1);
        R_i(Dflt(:),1)=0; % A defaulted bank with no assets has no recovery value.
        
        %LIQUIDATION FRACTION FOR BANKS
        %-------------------------------
        %For Hybrid fund:
        Dummy=(HQ_i~=0); %want to pick out which banks have no external asset
        HZ_i=min(abs(min(HGamma_i-hg_i-HL_i0,zeros(I,1))./R), ones(I,1) );
        %Note margin is used and available and so not subtracted from cash
        
        %For Pure fund: assume the guarantee comes from cash assets
        Dummy=(PQ_i~=0); %want to pick out which banks have no external asset
        PZ_i=min(abs(min(PGamma_i-pg_i-transpose(NetLambdaMinus_i),zeros(I,1))./R),ones(I,1));
        %Note: margin is subtracted from the available cash/cash is earmarked for margin
        
        %GUARANTEE FUND CONTRIBUTION OF BANKS:
        %-------------------------------------
        %In Hybrid Fund:
        HG_i=max(transpose(Lambda_i)+hg_i, zeros(I,1))-transpose(NetLambdaPlus_i);
        HG_tot=sum(HG_i);
        
        %In Pure Fund:
        PG_i=max(transpose(Lambda_i)+PGamma_i+R_i, zeros(I,1))-max(transpose(Lambda_i)+PGamma_i+R_i-pg_i,zeros(I,1));
        PG_tot=sum(PG_i);
    end
    
    %DEFAULT FUND CONTRIBUTION OF BANKS:
    %-------------------------------------
    %In Hybrid Fund:
    HD_i=max(transpose(Lambda_i)+hg_i+hd_i, zeros(I,1))-max(transpose(Lambda_i)+hg_i,zeros(I,1));
    HD_tot=sum(HD_i);
    
    %In Pure Fund:
    PD_i=max(transpose(Lambda_i)+PGamma_i+R_i-pg_i, zeros(I,1))-max(transpose(Lambda_i)+PGamma_i+R_i-pg_i-pd_i,zeros(I,1));
    PD_tot=sum(PD_i);
    %_________________________________________________________________________
    
    if Pred~=0
        %Time 1 BUDGET OF CCP
        %----------------------
        %Set Volume fee level for service
        f=0.10; %10% fee for clearing
        
        %Set Liability of CCP
        L_0i=(1-f).*transpose(NetLambdaPlus_i);
        L_0=abs(sum(L_0i));
        
        %Hybrid Case:
        HA_0=Gamma_0 + sum(hg_i) + sum(HL_i0);
        HL_0=abs(min(sum(L_0i) + sum(HG_i)+(Gamma_0 + f.*sum(NetLambdaPlus_i)), 0));
        
        %Pure Case:
        PA_0= Gamma_0 + sum(pg_i) + sum(PL_i0);
        PL_0=abs(min(sum(L_0i) + sum(PG_i)+(Gamma_0 + f.*sum(NetLambdaPlus_i)), 0));
        %The nominal wealth is (Gamma_0 + f.*sum(NetLambdaPlus_i))
        %_____________________________________________________________________
        
        %EQUILIBRIUM PAYMENT OF BANK TO CCP
        %----------------------------------
        %Hybrid Case:
        EHL_i0=abs(min(HL_i0,(HGamma_i-hg_i+R_i)));
        %Pure Case:
        EPL_i0=abs(min(PL_i0,(PGamma_i-pg_i+R_i)+pg_i));
        
        %Shorfall to CCP
        HLMinus_i0=abs(HL_i0-EHL_i0);
        PLMinus_i0=abs(PL_i0-EPL_i0);
        %______________________________________________________________________
        
        %------------------------------------------------------
        %                TOTAL EQUIB HOLDINGS
        %                   END OF PERIOD
        %------------------------------------------------------
        %Time 1 EQUILIBRIUM BUDGET OF CCP
        %--------------------------------
        %Hybrid Case:
        EHA_0= Gamma_0 + sum(hg_i) + sum(EHL_i0);
        EHL_0=abs(min(EHA_0,HL_0));
        
        %Pure Case:
        EPA_0= Gamma_0 + sum(pg_i) + sum(EPL_i0);
        EPL_0=abs(min(EPA_0,PL_0));
        
        %Proportionality Rule for Disbursements with Limited Liability:
        if sum(L_0i)>0
            HPI_0i=max(L_0i/sum(L_0i),zeros(I,1));
            PPI_0i=max(L_0i/sum(L_0i),zeros(I,1));
        else
            HPI_0i=zeros(I,1);
            PPI_0i=zeros(I,1);
        end
    end
    
    %THE TOTAL GUARANTEE & DEFAULT FUND HELD BY THE CCP
    %---------------------------------------
    %Hybrid Fund
    EHG_tot=sum(min(HG_tot,max(HA_0+EHL_0-Gamma_0-f.*sum(NetLambdaPlus_i),0)));
    EHD_tot=sum(min(HD_tot,max(HD_i-sum(max(HG_i-HCMinus_i,zeros(I,1))),zeros(I,1))));
    
    %Pure Fund
    EPG_tot=sum(min(PG_tot,max(PA_0+EPL_0-Gamma_0-f.*sum(NetLambdaPlus_i),0)));
    EPD_tot=sum(min(PD_tot,sum(max(PD_i-min(PG_i-PCMinus_i,zeros(I,1)),zeros(I,1)))));
    
    
    %THE PORTION OF GUARANTEE/DEFAULT FUND RETURN TO BANKS
    %-----------------------------------------------
    %Hybrid Fund
    %---------------
    %(Guarantee Fund)
    if HG_tot~=0
        EHG_i=(HG_i./HG_tot).*EHG_tot;
    else
        EHG_i=0.*HG_i;
    end
    EHG_tot2=sum(EHG_i);
    %(Default Fund)
    if HD_tot~=0
        EHD_i=(HD_i./HD_tot).*EHD_tot;
    else
        EHD_i=0.*HD_i;
    end
    EHD_tot2=sum(EHD_i);
    
    %Pure Funds
    %----------------
    %(Guarantee Fund)
    if PG_tot~=0
        EPG_i=(PG_i./PG_tot).*EPG_tot;
    else
        EPG_i=0.*PG_i;
    end
    EPG_tot2=sum(EPG_i);
    %(Default Fund)
    if PD_tot~=0
        EPD_i=(PD_i./PD_tot).*EPD_tot;
    else
        EPD_i=0.*PD_i;
    end
    EPD_tot2=sum(EPD_i);
    
    if Pred~=0
        %THE ACCOUNT AND SHORTFALL OF CCP
        %----------------------------------
        %Need to determine shortfall of CCP into the next round.
        %From default condition and Seizure we need a CCP liquidation decision.
        %From default condition and Seizure we need a CCP liquidation decision.
        %Terminal Net Worth
        HC_0=(EHA_0)-HL_0-EHG_tot;
        PC_0=(EPA_0)-PL_0-EPG_tot;
        
        %Terminal Shortfall %might need to do something here with the tranch.
        HCMinus_0=abs(min(HC_0,0));
        HCPlus_0=max(HC_0,0);%won't add loss from previous period until the end.
        PCMinus_0=abs(min(PC_0,0));
        PCPlus_0=max(PC_0,0);
        
        %THE ACCOUNT AND SHORTFALL OF BANKS
        %-----------------------------------
        %Hybrid Terminal Assets
        HA_i=HGamma_i+HZ_i.*R_i+(1-HZ_i).*HQ_i+EHL_0.*HPI_0i+EHG_i-hg_i;
        HC_i=HA_i-HL_i0;
        HCMinus_i=HL_i0-EHL_i0;
        HCPlus_i=HGamma_i+HZ_i.*R_i+(1-HZ_i).*HQ_i;
        
        %Pure Terminal Assets
        PA_i=PGamma_i+PZ_i.*R_i+(1-PZ_i).*PQ_i+EPL_0.*PPI_0i+EPG_i-pg_i;
        PC_i=PA_i-PL_i0;
        PCMinus_i=PL_i0-EPL_i0;
        PCPlus_i=PGamma_i+PZ_i.*R_i+(1-PZ_i).*PQ_i;
        
        %In the aggregate
        HCPlus_tot=sum(HCPlus_i);
        PCPlus_tot=sum(PCPlus_i);
        Surplus=HCPlus_tot-PCPlus_tot;
        BSurplus=BSurplus+Surplus; %for plotting
        
        Text3='The difference in the surplus of banks due to liquidations in Hybrid vs. Pure Fund:';
        disp(Text3)
        disp(Surplus)
        
        %_________________________________________________________________________
        %% Time 1 - Liquidation Defaults
        
        %DETERMINE THE SHORTFALL THE CCP TAKES ON FOR NEXT ROUND
        %---------------------------------------------------------
        %   !!!!!!!!!!!!!!!!!!! DEFAULT AND SEIZURE !!!!!!!!!!!!!!!!!!!
        %--------------------------------------------------------------------
        
        %THE BANKS DEFAULT CONDITION
        %-----------------------------------------------
        %In Hybrid Case
        for i=1:I
            if  HL_i0(i)>=(HGamma_i(i)+R_i(i))
                Hnum=Hnum+1;
                HDefaultdummy(i,1,l)=i; %fills by default per time step
                HDefault=sum(HDefaultdummy,3); %sums defaults
                numdummy=[1:I]; %sets up a index matrix of banks
                numdummy=transpose(numdummy); %makes it into column
                for j=1:length(numdummy)
                    if HDefault(j)~=numdummy(j) %if there are multiples
                        if HDefault(j)~=0 %checks that nonzero
                            HDefault(j)=numdummy(j); %makes sure that entry=index
                        end
                    end
                end
            end
            
            %In Pure Case
            if PL_i0(i)>=PGamma_i(i)+R_i(i)+pg_i(i)
                Pnum=Pnum+1;
                PDefaultdummy(i,1,l)=i;
                PDefault=sum(PDefaultdummy,3);
                numdummy=[1:I];
                numdummy=transpose(numdummy);
                for j=1:length(numdummy)
                    if PDefault(j)~=numdummy(j)
                        if PDefault(j)~=0
                            PDefault(j)=numdummy(j);
                        end
                    end
                end
            end
        end
        
        EqualCheck='Is the Hybrid Default the same as the Pure Default';
        disp(EqualCheck)
        disp(HDefault')
        disp(PDefault')
        
        
        ell=l+1;
        disp('Defaulted this buyback round ')
        disp(ell)
        disp(transpose(Default))
        
        for i=1:length(HDefault)
            if HDefault(i)~=0
                Dflt(dnum)=HDefault(i);
                dnum=dnum+1;
            end
        end
        
        disp('The buyback evolution of default is: ')
        DfltEvo=Dflt;
        disp(DfltEvo)
        Dflt=sort(Dflt);
        Dflt=Dflt(cumsum(Dflt(:),2)>0);
        Dflt = unique(Dflt,'sorted')
        disp('The Dflt is: ')
        disp(Dflt)
        
        Dfltcount=length(Dflt)
        
        %% THE CCP SEIZURE RULE
        %-----------------------
        %allows default fund updating unambiguously
        HSeizure=(HL_i0-HGamma_i-R_i>=hg_i);
        PSeizure=(PL_i0-PGamma_i-R_i>=pg_i);
        dnum=1;
        Default=[D0];
        for i=1:I
            if HSeizure(i,1)==1 & PSeizure(i,1)==1;
                dnum=dnum+1;
                Default(dnum,1)=i;
            end
        end
        %---------------------------------------------------------------
        
        %for plot
        BDefault=Default;
        BDefault=sum(BDefault>0);
        
        %THE LIQUIDATION DECISION OF BANKS FOR NEXT ROUND
        %------------------------------------------------
        %Can determine the new X matrix with liquidation
        X_ij=X_ij+a_ij; %each round we replace the matrix
        X_ji=X_ji+a_ji;
        
        %END OF DAY MARGIN REFILL
        %-------------------------------------------------
        %HYBRID:
        %------------
        %%Need to reduce the amount initial margin fund (if used) in Hybrid Case
        hgMin=HL_i0; %remaining liability not covered by initial margin
        hg_i=max(hg_i-transpose(NetLambdaMinus_i),zeros(I,1)); %amount of margin used and left
        
        %%Need to reduce cash by amount that still owed to CCP
        HGammaMin=abs(min(HGamma_i-hgMin, zeros(I,1)));
        HGamma_i=max(HGamma_i-hgMin,zeros(I,1)); %amount of cash reduced by that not covered by guarantee
        
        %%Need to reduce the external asset by the amount that must be liquidated
        HQMin_i=abs(min(HQ_i-HGammaMin,zeros(I,1)));%amount the external fund cannot cover should equal CMin_i
        HQ_i=HQ_i.*(1-HZ_i); %fraction needed needed to be liquidated
        
        %%Need to reduce own default fund contribution
        hd_iMin=abs(min(hd_i-HQMin_i,zeros(I,1)));
        Hd_i=max(hd_i-HQMin_i,zeros(I,1));
        
        %%Any leftover liability goes to the CCP in an DL_i0=DL_0 (line #352)
        HDL_i0=abs(hd_iMin); %amt transferred to CCP from liabilities
        %%Those that don't have enough to cover margin call are in Default
        for i=1:length(HDefault)
            if transpose(HDL_i0(i))~=0
                HDefault(i)=i;
            end
        end
        
        %Left over default liability not covered by bank must go to the CCP
        EDHL_0=EDHL_0+sum(HDL_i0);
        
        %Any surplus of Defaulted banks must go toward min their liability
        dummy=(HDefault~=0); %non-zero defaulted banks
        EDHA_0=EDHA_0+sum(HCPlus_i.*dummy);
        
        %PURE:
        %-----------
        %%Need to reduce cash by amount that still owed to CCP
        PGammaMin=abs(min(PGamma_i-PL_i0,zeros(I,1))); %Amount lacking to pay CCP in cash
        PGamma_i=max(PGamma_i-PL_i0,zeros(I,1)); %Amount left in account
        
        %%Need to reduce the external asset by the amount liquidated
        PQMin_i=abs(min(PQ_i-PGammaMin,zeros(I,1))); %amount the external fund cannot cover
        PQ_i=PQ_i.*(1-PZ_i); %fraction needed needed to be liquidated
        
        %%If cash&external asset not enough for liability then use initial margin
        pg_iMin=abs(min(pg_i-PQMin_i,zeros(I,1))); %amt margin can't cover
        pg_i=max(pg_i-PQMin_i,zeros(I,1)); %amt removed from initial margin
    end
    
    %%Need to reduce own default fund contribution
    pd_iMin=abs(min(pd_i-PQMin_i,zeros(I,1)));
    pd_i=max(pd_i-PQMin_i,zeros(I,1));
    
    %%Any leftover liability goes to the CCP in an DL_i0=DL_0 (line #352)
    PDL_i0=abs(pd_iMin); %amt transferred to CCP from liabilities
    %%Those that don't have enough to cover margin call are in Default
    for i=1:length(PDefault)
        if transpose(PDL_i0(i))~=0
            % if transpose(PCMinus_i(i))~=0
            PDefault(i)=i;
        end
    end
    
    %Left over default liability not covered by bank must go to the CCP
    EDPL_0=EDPL_0+sum(PDL_i0); %should equal PCMinus_i
    
    %Any surplus of Defaulted banks must go toward min their liability
    Dummy=(PDefault~=0); %non-zero defaulted banks
    EDPA_0=EDPA_0+sum(PCPlus_i.*Dummy);
    
    %Liquidation loss incurred vs. assets made on liquidation only up to second
    HLiquidationGain=max(EDHA_0-EDHL_0,zeros(1,1))
    HLiquidationLoss=abs(min(EDHA_0-EDHL_0,zeros(1,1)))
    PLiquidationGain=max(EDPA_0-EDPL_0,zeros(1,1))
    PLiquidationLoss=abs(min(EDPA_0-EDPL_0,zeros(1,1)))
    HCCPLoss=min(HC_0+HLiquidationGain-HLiquidationLoss,0)
    PCCPLoss=min(PC_0+PLiquidationGain-PLiquidationLoss,0)
    
    if HCCPLoss<0
        Flag='The Hybrid CCP has FAILED! At buyback time-step: '
        bhell=l-1;
        disp(bhell)
        flag=L; % will end the looping after this round
    end
    if PCCPLoss<0
        Flag='The Pure CCP has FAILED! At buyback time-step:'
        bpell=l-1;
        disp(bpell)
        flag=L; % will end the looping after this round
    end
    
    %for plotting we need separate values that don't cumulate
    BHLiquidationGain=HLiquidationGain;
    BHLiquidationLoss=HLiquidationLoss;
    BPLiquidationGain=PLiquidationGain;
    BPLiquidationLoss=PLiquidationLoss;
    BHCCPLoss=HCCPLoss;
    BPCCPLoss=PCCPLoss;
    
    if Pred~=0
        %A Method of Checking Predation Profits.
        XF=sum(X_ij,2);%gets nominal value of positions.
        VF=XF.*repmat(S_l(I,1,:,l),[I,1]).*Predator_k; %gets final value of interest paid
        VF=XF.*repmat(S_l(I,1,:,l),[I,1]).*Predator_k; %gets final value of interest paid
        V00=X00.*repmat(S0,[I,1]).*Predator_k; %value of starting position size at start
        VF0=X00.*repmat(S_l(I,1,:,l),[I,1]).*Predator_k;%value of starting position now.
        BuybackProfitsVF=VF-VF0; %Buyback profits on using end price
        BuybackProfitsV0=VF-V00; %buyback profits on using starting price and end price
        
        %For Plotting
        BBuybackProfitsVF=sum(sum(BuybackProfitsVF,3))/I; %Buyback profits on using end price
        BBuybackProfitsV0=sum(sum(BuybackProfitsV0,3))/I;
        
        %BREAK CONDITION: To break for more than five buyback periods.
        %Only allow for break on 5 days of exogenous default Dflt(1)
        last=length(Dflt);
        for k=1:K
            if  abs(sum(sum(X_ij(:,:,k),2)))==Xmarket_k(:,:,k)
                Flag='The holdings of predators have reached the maxx allowable holding'; %since distressed can't do anything here
                display(Flag) %sets max holding to be the maximum of asset in market, slightly larger than start value.
                bell=l-1;
                display(bell)
                flag='P'; % will end the looping after this round
            end
            if sum(sum(a_ij(:,:,k),2))==0
                Flag='The predators have reached capacity';
                display(Flag)
                sell=l-1;
                display(sell)
                flag='P'; % will end the looping after this round
            end
        end
        %___________________________________________________________________
        
        %% Time 2 - Break Conditions
        if flag==L & CCPFail=='Y';
            disp('CCP Failed due to Insolvency in buyback step')
            break
        end
        
        if flag=='P';
            disp('The Predation Threshold was reached')
            break
        end
        
    end
end

%% Time 2- Terminating Values
%THE TERMINATING VALUE OF THE BANK AND CCP FOR DISTRIBUTION OF G & D FUND DECISIONS
%--------------------------------------------------------------------
%These values should be carried over to new period.
HEpsilon=0.10;  %needs to be adjust for amount needed by CCP
EHTranche=(HEpsilon).*(Gamma_0+f.*sum(NetLambdaPlus_i));
NetEHTranche=max(EHTranche-abs(sum(min(EHG_i+EHD_i+HCMinus_i,zeros(I,1)))),0);
THC_0=THC_0+(Gamma_0+f.*sum(NetLambdaPlus_i))-abs(sum(min(EHG_i+EHD_i+HCMinus_i,zeros(I,1))))+HLiquidationGain-HLiquidationLoss; %Check this addition

PEpsilon=0.10;
EPTranche=(PEpsilon).*(Gamma_0+f.*sum(NetLambdaPlus_i));
NetEPTranche=max(EPTranche-abs(sum(min(EPG_i+EPD_i+PCMinus_i,zeros(I,1)))),0);
TPC_0=TPC_0+(Gamma_0+f.*sum(NetLambdaPlus_i))-abs(sum(min(EPG_i+EPD_i+PCMinus_i,zeros(I,1))))+PLiquidationGain-PLiquidationLoss;

DHG_tot=HG_tot;
DHD_tot=HD_tot;
if HG_tot==0
    DHG_tot=1;
end
if HD_tot==0
    DHD_tot=1;
end
THC_i=THC_i+(HGamma_i+HQ_i+transpose(Lambda_i))-((HPI_0i.*HCMinus_0+HZ_i.*(HQ_i-R_i)+f.*transpose(NetLambdaPlus_i))+(HG_i/DHG_tot.*(HG_tot-EHG_tot)-HD_i/DHD_tot.*(HD_tot-EHD_tot)));

DPG_tot=PG_tot;
DPD_tot=PD_tot;
if PG_tot==0
    DPG_tot=1;
end
if PD_tot==0
    DPD_tot=1;
end
TPC_i=TPC_i+(PGamma_i+PQ_i+transpose(Lambda_i))-((PPI_0i.*PCMinus_0+PZ_i.*(PQ_i-R_i)+f.*transpose(NetLambdaPlus_i))+(PG_i/PG_tot.*(DPG_tot-EPG_tot)-PD_i/DPD_tot.*(PD_tot-EPD_tot)));


%For Plotting
BTHC_0=THC_0;
BTPC_0=TPC_0;
BTHC_i=sum(THC_i);
BTPC_i=sum(TPC_i);
BEHG_tot=EHG_tot2;
BEPG_tot=EPG_tot2;
BEHD_tot=EHD_tot2;
BEPD_tot=EPD_tot2;

%***********************************************************************
%***********************************************************************
%% Time=3 - Recovery Period
%--------------------Start of last period--------------------
%t=3
%___________________________________________________________
%%%%%%%%%%%%%%%%%TIME 3 - RECOVERY PERIOD%%%%%%%%%%%%%%%%%%%
%____________________________________________________________
%*************************************************************************
%EXPLANATION: THE CCP MAKES A MARGIN CALL ON SURVIVING BANKS NEW POSITION,
%THEY MUST ADD TO THEIR INITIAL MARGIN NOT ONLY WHAT THEY USED, BUT WHAT
%WAS DEPLETED BY OTHERS, BUT ONLY IN THE HYBRID CASE.
%*************************************************************************
%%Adjust for new evaluation of initial margin -- "Margin Call"
margin=0.60;
Newhg_i=initmargin(margin,l,S_l,X_ij); %new margin level
Newhg_i(Dflt(:))=0; %sets all defaulted banks to zero.
hgMissing=min(hg_i-Newhg_i,zeros(I,1)); %amount needed to replenish fun
hgReturn=max(hg_i-Newhg_i,zeros(I,1)); %amount to return to bank

%%Need to reduce cash by amount taken to replenish fund or to return
hgHGammaMin=min(HGamma_i+hgMissing+hgReturn, zeros(I,1)); %amount bank
%cannot meet and transferred to external fund.
HQ_i=max(HQ_i+hgHGammaMin,zeros(I,1)); %amt removed from external fund
HQMin_i=min(HQ_i+hgHGammaMin,zeros(I,1)); %amt external fund cannot cover
HGamma_i=max(HGamma_i+hgMissing+hgReturn,zeros(I,1)); %amt to
%increase/decrease initial margin

%Amount needed to replenish margin fund
HGR_i=abs(min(EHG_i-Newhg_i,zeros(I,1))); %only own but need to see how depleted
NewHCPlus_i=max(HCPlus_i-HGR_i,zeros(I,1));
HCPlusLoss_i=max(HCPlus_i-NewHCPlus_i,zeros(I,1))

%For Plotting
AvgHGR_i=max(HGR_i,zeros(I,1))
AvgHGR_i=sort(AvgHGR_i)
AvgHGR_i=AvgHGR_i(cumsum(AvgHGR_i(:),2)>0)
AvgHGR=sum(AvgHGR_i)/sum(AvgHGR_i>0)
AvgHCPlusLoss_i=max(HCPlusLoss_i,zeros(I,1))
AvgHCPlusLoss_i=sort(AvgHCPlusLoss_i)
AvgHCPlusLoss_i=AvgHGR_i(cumsum(AvgHGR_i(:),2)>0)
AvgHCPlusLoss=sum(AvgHCPlusLoss_i)/sum(AvgHCPlusLoss_i>0)


%PURE:
%-----------
%%Need to reduce amount in initial margin if used in Pure Case:
%%Adjust for new evaluation of initial margin -- "Margin Call"
margin=.60; %set through input
Newpg_i=initmargin(margin,l,S_l,X_ij); %new margin level
Newpg_i(Dflt(:))=0; %sets all defaulted banks to zero.
pgMissing=min(pg_i-Newpg_i,zeros(I,1)); %amount needed to replenish fun
pgReturn=max(pg_i-Newpg_i,zeros(I,1)); %amount to return to bank

%%Need to reduce cash by amount taken to replenish fund or to return
pgPGammaMin=min(PGamma_i+pgMissing+pgReturn, zeros(I,1)); %amount bank
%cannot meet and transferred to external fund.
PQ_i=max(PQ_i+pgPGammaMin,zeros(I,1)); %amt removed from external fund
PQMin_i=min(PQ_i+pgPGammaMin,zeros(I,1)); %amt external fund cannot cover
PGamma_i=max(PGamma_i+pgMissing+pgReturn,zeros(I,1)); %amt to
%increase/decrease initial margin

%Amount needed to replenish margin fund with what has been returned to you
PGR_i=abs(min(EPG_i-Newpg_i,zeros(I,1)));
NewPCPlus_i=max(PCPlus_i-PGR_i,zeros(I,1));


%% Time 3 - Terminating Values
%----------------End of Period 3/Recovery Period time -------------------

%THE TERMINATING VALUE OF THE BANK AND CCP FOR DISTRIBUTION OF G & D FUND DECISIONS
%--------------------------------------------------------------------
%These values should be carried over to new period.
THC_0=THC_0+(Gamma_0+f.*sum(NetLambdaPlus_i))-abs(min(EHG_tot+EHD_tot+sum(HCMinus_i),0))+HLiquidationGain-HLiquidationLoss; %Check this addition
TPC_0=TPC_0+(Gamma_0+f.*sum(NetLambdaPlus_i))-abs(min(EPG_tot+EPD_tot+sum(PCMinus_i),0))+PLiquidationGain-PLiquidationLoss;

DHG_tot=HG_tot;
DHD_tot=HD_tot;
if HG_tot==0
    DHG_tot=1;
end
if HD_tot==0
    DHD_tot=1;
end
NewTHC_i=THC_i-HGR_i;
THC_i=THC_i+(HGamma_i+HQ_i+transpose(Lambda_i))-((HPI_0i.*HCMinus_0+HZ_i.*(HQ_i-R_i)+f.*transpose(NetLambdaPlus_i))+(HG_i/DHG_tot.*(HG_tot-EHG_tot)-HD_i/DHD_tot.*(HD_tot-EHD_tot))+HGR_i);
HProfitLoss=(NewTHC_i./THC_i);

DPG_tot=PG_tot;
DPD_tot=PD_tot;
if PG_tot==0
    DPG_tot=1;
end
if PD_tot==0
    DPD_tot=1;
end
NewTPC_i=TPC_i-PGR_i;
TPC_i=TPC_i+(PGamma_i+PQ_i+transpose(Lambda_i))-((PPI_0i.*PCMinus_0+PZ_i.*(PQ_i-R_i)+f.*transpose(NetLambdaPlus_i))+(PG_i/PG_tot.*(DPG_tot-EPG_tot)-PD_i/DPD_tot.*(PD_tot-EPD_tot))+PGR_i);
PProfitLoss=((TPC_i-NewTPC_i)./TPC_i);

%For Plotting
RTHC_0=THC_0;
RTPC_0=TPC_0;
RTHC_i=sum(THC_i);
RTPC_i=sum(TPC_i);
AvgRTHC=sum(THC_i)/I;
AvgRTPC=sum(TPC_i)/I;
REHG_tot=EHG_tot2;
REPG_tot=EPG_tot2;
REHD_tot=EHD_tot2;
REPD_tot=EPD_tot2;
RNewTHC_i=sum(NewTHC_i);
RNewTPC_i=sum(NewTPC_i);
AvgRNewTHC_i=sum(NewTHC_i)/I;
AvgRNewTPC_i=sum(NewTPC_i)/I;

%Predator ProfitLoss only
dummy1=(THC_i>0);
PredHPL=dummy1.*((THC_i-NewTHC_i)./THC_i);
HSurvNum=sum(PredHPL>0);
if HSurvNum==0
    PredHPL=0;
else
    PredHPL=sum(PredHPL)/I;
end
PredHProfitLoss=PredHPL(1,1);


dummy2=(TPC_i>0);
PredPPL=dummy2.*((TPC_i-NewTPC_i)./TPC_i);
PSurvNum=sum(PredPPL>0);
if PSurvNum==0
    PredPPL=0;
else
    PredPPL=sum(PredPPL)/I;
end
PredPProfitLoss=PredPPL(1,1);

%% Time 3 - Return Function
Predator=Pred;
Distress=Dist;

Values0=[Xmarket;N];
Values1=[Predator;Distress;Dfltcount]; %setup values
ValuesL1=[LHCCPLoss;LPCCPLoss;lell;pell;hell;LDflt;LTHC_0;LTPC_0;LTHC_i;LTPC_i];
ValuesL2=[LHLiquidationLoss;LPLiquidationLoss;LHLiquidationGain;LPLiquidationGain;LEHG_tot;LEPG_tot];
ValuesB1=[BHCCPLoss;BPCCPLoss;bhell;bpell;BDefault;BBuybackProfitsVF;BBuybackProfitsV0;bell;sell];
ValuesB2=[BTHC_0;BTPC_0;BTHC_i;BTPC_i;BEHG_tot;BEPG_tot;BEHD_tot;BEPD_tot];
ValuesB3=[BHLiquidationLoss;BPLiquidationLoss;BHLiquidationGain;BPLiquidationGain;AvgHGR;AvgHCPlusLoss];
ValuesR1=[RTHC_0;RTPC_0;RTHC_i;RTPC_i;AvgRTHC; AvgRTPC;REHG_tot;REPG_tot;REHD_tot;REPD_tot];
ValuesR2=[RNewTHC_i;RNewTPC_i;AvgRNewTHC_i;AvgRNewTPC_i;LSurplus;BSurplus];
ValuesR3=[PredHProfitLoss;PredPProfitLoss]
Values=[Values0;Values1;ValuesL1;ValuesL2;ValuesB1;ValuesB2;ValuesB3;ValuesR1;ValuesR2;ValuesR3];
end