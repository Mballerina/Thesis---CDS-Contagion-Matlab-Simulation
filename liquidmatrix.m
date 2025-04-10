%*******************************************************************************************************%
%*                                    LIQUIDATION MATRIX MAKER                                         *%                                            
%*                                                                                                     *% 
%*******************************************************************************************************%

function a_ij=liquidmatrix(N_ijk,Distress_k,Dist,Predator_k,Pred,LPred,Collusion,a_ij,Dflt,X_ij,I,J,K)

%CCP DETERMINES ITS TRADING RATE AND STARTING MARKET TRADING RATE
A_ki=abs(sum(a_ij,2));
A_k=sum(A_ki,1); %trading rate for the period for each asset

%trading rate for CCP
a_0=A_k/I; %at each period the CCP will adjust its trading rate until it converges 

s=randi(100,1,1);
rng(s);
%a=(-1).*N_ijk./10; %might need more room
%b=N_ijk/10;
%a_ij=round((b-a).*randn(I,J,K)+a); %need this to generate random for non predators.
m=0.516e9/I;
m=m/K;
s=0.646e9/(I);
s=s/K;
pd = makedist('Normal','mu', m/(I-1),'sigma',s/(I-1));
a_ij=round(random(pd,I,J,K)); %next round it should be the one from previous round.

%For each period we have a liquidation matrix:
%have to set the diagonal zero running over k
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

%set a_ji to be the opposite side of liquidations
a_ji=permute(a_ij,[2,1,3]);

%As a check each bank would need to liquidate for the following liquid time
X_i=sum(X_ij,2);
%tao=abs(X_i./repmat(a_0,[I,1]));

%the CCP then has for each i's assets a column vector
D=length(Dflt);
a_dflt=zeros(I,J,K); %creates 4 rows, 4 columns of zeros
a_dflt(Dflt(1:D),:,:)=repmat(a_0,I,1); %inserts one column liquidations
Dummy=(-1)*(X_ij(:,:,:)<0)+(X_ij(:,:,:)>0); %gets the direction of the holdings

Dummy=(-1).*Dummy; %Need the opposite direction for liquidation

a_dflt(Dflt(1:D),:,:)=Dummy(Dflt(1:D),:,:).*a_dflt(Dflt(1:D),:,:); %demands liquidation in opposite direction
Dummy1=(X_ij(:,:,:)~=0); %Want the opposite, but this gives zero for case TRUE
a_dflt(Dflt(1:D),:,:)=a_dflt(Dflt(1:D),:,:).*Dummy1(Dflt(1:D),:,:); %Inputs back into matrix

%CCP initiates all defaults for bank 4 with ji liquidating completely.
a_ij(Dflt(1:D),:,:)=0; %take out the line for the CCP so no addition
a_ij=a_ij+a_dflt; %combines the two matrices


for k=1:K
%all distressed must have random selling
distress=sort(Distress_k); %get zeros first
distress=distress(cumsum(distress(:,:,k),2)>0,:,k); %create smaller matrix with less
a_ij(distress(:,:),Dflt(1:D),k)=a_0(:,:,k); %set the distressed to zero
%If no predators than everyone, but CCP has random selling
count=sum(Predator_k~=0);
DP=length(Predator_k(:,:,k));
    if count(:,:,k)==0 | LPred(:,:,k)~=0
        a_ij=a_ij;
    %If predator is 1 then predator acts alone with CCP 
    elseif count(:,:,k)==1 & LPred(:,:,k)~=0
       a_ij(LPred(:,:,k),Dflt(1:D),k)=Dummy(LPred(:,:,k),Dflt(1:D),k).*repmat(a_0(:,:,k),[D,1]);
    %choose predator with largets holding and most money
%If multiple predators exist
    elseif count>1 & Collusion=='Y'
        %If collusion then all predators liquidate at the same rate as CCP
        for d=1:DP
            if Predator_k(d)~=0
                a_ij(Predator_k(d),Dflt(1:D),k)=Dummy(Predator_k(d),Dflt(1:D),k).*repmat(a_0(:,:,k),[D,1]);
            end
        end
    elseif count>1 & Collusion=='N'
        %If collusion doesn't exist they all liquidate at different rate 
        for d=1:DP
            if Predator_k(d)~=0
                a_ij(Predator_k(d),Dflt(1:D),k)=a_0(:,:,k); %this rate will be different for buyback
            end
        end
    end
end


%Check again that no one is liquidating more than they have in possession
%Make sure that cannot liquidation more than have in posession
Dummy1=(abs(a_ij)>abs(X_ij));
V1=(-1).*Dummy1.*X_ij;
V2=a_ij-a_ij.*Dummy1;
a_ij=V2+V1;

%Check that no one liquidates more than CCP (May not need this)
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