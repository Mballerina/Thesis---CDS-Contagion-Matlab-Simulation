%*******************************************************************************************************%
%*                                       BUY-BACK MATRIX MAKER                                         *%                                            
%*                                                                                                     *% 
%*******************************************************************************************************%


function a_ij=buybackmatrix(CCPTrade,DISTTrade,N_ijk,X_0,Xmax,Distress_k,Dist,Predator_k,Pred,LPred,Collusion,a_ij,Dflt,X_ij,I,J,K)

%Ideally the CCP should be done liquidating defaulted assets, but may need
%to keep going in the event of failures. In that case the buyback period
%has more periods. The defaulted assets from the liquidation period should
%be zero at this point.

%There is the option that distressed banks may be able to liquidate (not entertained here)

%CCP DETERMINES ITS TRADING RATE AND STARTING MARKET TRADING RATE
A_ki=abs(sum(a_ij,2)); %get the aggregate rate of trading in market, for each bank, over each asset
A_k=sum(A_ki,1); %aggregate trading rate for the period, for each asset, over all banks

%trading rate for CCP:at each period the CCP will adjust its trading rate until it converges 
a_0=A_k/I; %the average market rate (agg. rate divided by number of banks)

s=randi(100,1,1); %sets a random but reproducable seed for each call value trading rate
rng(s); %calls seed
%a=N_ijk/10; %lower bound rate taken from the empirical market size
%b=N_ijk/10; %upper bound for rate taken from above
%a_ij=round((b-a).*randn(I,J,K)+a); %generates random trading rate for predators based on bounds.

%CALIBRATION:Distribution from empirical literature (Oehmke)
m=0.516e9/I; % mean divided by number of bank number
m=m/K; %mean further divided by the number of asset number
s=0.646e9/(I); %standard deviation divided by bank number
s=s/K; %standard deviation further divided by asset number
pd = makedist('Normal','mu', m/(I-1),'sigma',s/(I-1)); %creates distribution from above values
a_ij=round(random(pd,I,J,K)); %generates random trading rate for preds based on distribution

%For each period we have a trading rate matrix:
%To set the diagonal zero running over k:
for k=1:K
    A(:,:,k)=eye(I).*a_ij(:,:,k);
    a_ij(:,:,k)=-A(:,:,k) + a_ij(:,:,k);
end

%Make sure that cannot have a trade size of more than have in posession
%(Ensures that CCP cannot liquidate beyond zero.)
Dummy1=(abs(a_ij)>abs(X_ij));
V1=Dummy1.*X_ij;
V2=a_ij-a_ij.*Dummy1;
a_ij=V2+V1;

%Set a_ji to be the opposite side of liquidations (counterparty liquidation)
a_ji=permute(a_ij,[2,1,3]);

%As a check: Each bank needs time X_i for full liquidation
X_i=sum(X_ij,2);
%tao=abs(X_i./repmat(a_0,[I,1]));
   

%CCP liquidates for the defaulted banks through transposed default column 
Dflt=sort(Dflt); %gets zeros first
Dflt=Dflt(cumsum(Dflt(:),2)>0); %resets default vector to have only non-zero entries
Dflt = unique(Dflt,'sorted'); %makes sure there are no doubles
D=length(Dflt); %gets the number of defaulters

if CCPTrade=='Y' %if CCP is allowed to continue trading
a_dflt=zeros(I,J,K); %creates 4 rows, 4 columns of zeros
a_dflt(Dflt(1:D),:,:)=repmat(a_0,D,I); %inserts one column liquidations for each defaulter
DummyNorm=(-1)*(X_ij(:,:,:)<0)+(X_ij(:,:,:)>0); %gets the direction of the holdings
Dummy=(-1).*DummyNorm; %Need the opposite direction for liquidation
a_dflt(Dflt(1:D),:,:)=Dummy(Dflt(1:D),:,:).*a_dflt(Dflt(1:D),:,:); %demands liquidation in opposite direction
Dummy1=(X_ij(:,:,:)~=0); %Want the opposite, but this gives zero for case TRUE
a_dflt(Dflt(1:D),:,:)=a_dflt(Dflt(1:D),:,:).*Dummy1(Dflt(1:D),:,:); %Inputs back into matrix
%CCP initiates all defaults for bank 4 with ji liquidating completely.
a_ij(Dflt(1:D),:,:)=0; %take out the line for the CCP liquidation so there is no extra addition
a_ij=a_ij+a_dflt; %combines the two matrices (liquidation row plus full matrix w/o liquidation row)
end


%BUYBACK BEGINS IN THIS CASE:
D=length(Dflt); %gets the number of defaulters
for k=1:K
    if DISTTrade=='N' %otherwise random assignment of trading direction above as all banks in market
        %All Distressed banks are not in a position to continue sell or buyback
        distress=sort(Distress_k); %get zeros first
        distress=distress(cumsum(distress(:,:,k),2)>0,:,k); %create smaller matrix with nonzero entries
        a_ij(distress(:),Dflt(1:D),k)=0; %set the distressed banks with defaulted counterparties column to zero
    end
    
    %If there are no predators than everyone, but CCP has random selling
    X0=sum(X_0,2); %sums the original holdings for each bank with all counterparties
    count=sum(Predator_k~=0); %counts the number of non-zero predators
    dummya=(X_0>1); %obtains any positive (buy position) 
    dummyb=(-1).*(X_0<1 & X_0~=0); %gives original holdings that are non-zero and negative (sell position)
    dummy=dummya+dummyb; %combines so that all buy positions and all sell positions have correct sign
    check=(X_ij>=repmat(Xmax,[I,J])); %gets entries where any entry is larger or size of max threshold for predators
    DP=length(Predator_k(:,:,k)); %gets the length of the Predator vector (will contain zeros also so always I=banks)
    
    if count(:,:,k)==0 | LPred(:,:,k)==0  %if number of predators is zero or there is no largest predator
         a_ij(:,:,k)=0; %there is no buyback
    end
    %IF ONE PREDATOR THEN ACTS ALONG WITH CCP 
    if count(:,:,k)==1 %if there is one predator for the asset
       a_ij(LPred(:,:,k),Dflt(1:D),k)=dummy(LPred(:,:,k),Dflt(1:D),k).*A_k(:,:,k); %the largest predator will buyback the defaulted asset at avg market rate
       if LPred~=0
            if check(LPred(:,:,k),Dflt(1:D),k)==1 %if the holding is larger than the allowable max
                X_ij(LPred(:,:,k),Dflt(1:D),k)=Xmax(:,:,k); %then the holding is limited to Xmax = max threshold for predators
            end
       end
    %IF MULTIPLE PREDATORS
    elseif count(:,:,k)>1 & Collusion=='Y' %if there are more predators and they act as one,
        %If collusion then all predators liquidate at the same rate as CCP
        for d=1:count(:,:,k) %cycles over non-zero count of predators
            if Predator_k(d)~=0 %if the predator entry is not zero 
                a_ij(Predator_k(d),Dflt(1:D),k)=dummy(Predator_k(d),Dflt(1:D),k).*a_0(:,:,k); %trades at avg trading rate
                if check(Predator_k(d),Dflt(1:D),k)==1 %if the holding is larger than the allowable max threshold
                    X_ij(Predator_k(d),Dflt(1:D),k)=Xmax(:,:,k); %then the holding turns to Xmax
                end
            end
        end
    elseif count(:,:,k)>1 & Collusion=='N' %if multiple predators competing
        %If collusion doesn't exist they all liquidate at different rate 
        countd=sum(Distress_k~=0); %determine number of non-zero distressed banks
        a_k=(A_k.*countd)/I.*(count-1); %determines the rate for predators liquidation (Brunnermeier)
        for d=1:count(:,:,k) %counts over non-zero predators
            if Predator_k(d)~=0 %picks a non-zero predator for each k
                a_ij(Predator_k(d),Dflt(1:D),k)=dummy(Predator_k(d),Dflt(1:D),k).*abs(a_k(:,:,k)); %buyback in direction of org. holding (dummy determines the holding)
                if check(Predator_k(d),Dflt(1:D),k)==1 %if the holding aquired is larger than the allowable max
                    X_ij(Predator_k(d),Dflt(1:D),k)=Xmax(:,:,k); %then the holding turns to Xmax
                end
            end
        end
    end
end

%[MAY NEED TO REMOVE THIS ENTRY FOR THE BUYBACK PERIOD, AS IT IS A LOWER BOUND]
%Check again that no one is liquidating or buying more than they have in possession
clear Dummy1
Dummy1=abs(a_ij)>abs(X_ij); %if buyback matrix is larger than holdings matrix
V1=Dummy1.*X_ij; %isolates the holdings which are larger and maxes them at the holding level
V2=a_ij-a_ij.*Dummy1; %removes the higher entries from the matrix
a_ij=V2+V1; %pastes the two matrices together

%[ENSURE THAT NO ONE LIQUIDATES MORE THAN CCP AND CAUSES IMPACT] 
Dummy=(a_ij~=0).*repmat(a_0,I,J); %extends CCP liquid matrix to check against
Dummy1=(abs(a_ij)>abs(Dummy));
V1=Dummy1.*X_ij;
V2=a_ij-a_ij.*Dummy1;
a_ij=V2+V1;


%set a_ji to be the opposite side of liquidations
a_ji=permute(a_ij,[2,1,3]);

end

%%EXPLANATION: Since the predators can at most liquidate the same amount as
%%the CCP each day, if they have a smaller or equal holding, they are
%%liquidating as fast as possible without making an impact.

%%BASE CALIBRATION: For a run of 2 predators and random distressed in a
%%stable and constant financial market the with options:
%1.HCCPLoss=-5.3693e14
%2.HCCPLoss=-2.1424e12 add in limit buyback to less than size of holding
%3.HCCPLoss=-1.6466e12 add in limit to buyback rate below CCP
%4.HCCPLoss=-7.1767e+51  add in the CCP not being able to liquidate further
%5.HCCPLoss=-1.9198e+10  if effect 4 is not added but in the distressed trading randomly
%6.HCCPLoss=-7.1767e+51 all effects together