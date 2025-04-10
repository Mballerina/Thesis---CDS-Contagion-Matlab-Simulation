function C=pricetwo(l,D_k,Dflt,S_l,S_l1,S_l2,X_ij,X_ji,a_ij,a_ji,I,J,K)


%Pricing Functionals 
z=zeros(I,J,K);
%--------------

%P_0: Get the fundamental value of postions
Dummyi=(X_ij~=0);
Dummyj=(X_ji~=0);
P0_i=X_ij.*Dummyi.*S_l(:,:,:,l-1);
P0_i=max(P0_i,zeros(I,J,K));
P0_j=X_ji.*Dummyj.*S_l(:,:,:,l-1);
P0_j=max(P0_j,zeros(I,J,K));
%----------------------------------------------------------------------

D=length(Dflt); %determines the current length of default vector
D_k=D_k(1,1,1); %picks out one entry of market liquidity
D_k=D_k*ones(I,D,K); %makes a liquidity/market depth matrix of the size of number of defaulters

%P_1: Get the Primary Price Impact [KEY: Separate if statements]
epsilon=.98; %let this be user input
%sum over defaulted banks
p11=z; %initiates a zero matrix
Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
p11(:,Dflt(1:D),:)=Factor.*X_ij(:,Dflt(1:D),:); %gets the results for defaulted banks' counterparties
p11(Dflt(1:D),:,:)=permute(Factor,[2,1,3]).*X_ij(Dflt(1:D),:,:); %gets a result for defaulted banks

%sum over undefaulted banks
p12=z; %initiates zero matrix
p12=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
p12(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
p12(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
p12=p12.*epsilon.*Factor.*X_ij; %gets result for undefaulted banks

%Multiplier from CCP liquidation over defaulted banks
p1common=z; %initiates zero matrix
Factor_j=(X_ji(:,:,:)~=0);
if l<2
    p1common(:,Dflt(:),:)=abs(Dummyj(:,Dflt(:),:).*S_l(:,Dflt(:),:,l-1)).*(a_ji(:,Dflt(:),:)./(D_k(:,:,:))).*Factor_j(:,Dflt(:),:);
    p1common(Dflt(:),:,:)=abs(Dummyj(Dflt(:),:,:).*S_l(Dflt(:),:,:,l-1)).*(a_ji(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3])).*Factor_j(Dflt(:),:,:);                                    
else
    p1common(:,Dflt(:),:)=abs(Dummyj(:,Dflt(:),:).*S_l1(:,Dflt(:),:)).*(a_ji(:,Dflt(:),:)./(D_k(:,:,:))).*Factor_j(:,Dflt(:),:);
    p1common(Dflt(:),:,:)=abs(Dummyj(Dflt(:),:,:).*S_l1(Dflt(:),:,:)).*(a_ji(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3])).*Factor_j(Dflt(:),:,:);                                     
end
p1common=sum(p1common(Dflt(:),:,:),2);
if D>=2
    p1common=sum(p1common);
end 
p1common=repmat(p1common,I);


%Multiplier from Bank i in distress
p1common2a=z; %initiates zero matrix
Factor_i=(X_ij(:,:,:)~=0);
if l<2
    p1common2a(:,Dflt(:),:)=abs(Dummyi(:,Dflt(:),:).*S_l(:,Dflt(:),:,l-1)).*(a_ij(:,Dflt(:),:)./(D_k(:,:,:))).*Factor_i(:,Dflt(:),:);
    p1common2a(Dflt(:),:,:)=abs(Dummyi(Dflt(:),:,:).*S_l(Dflt(:),:,:,l-1)).*(a_ij(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3])).*Factor_i(Dflt(:),:,:);
else
    p1common2a(:,Dflt(:),:)=abs(Dummyi(:,Dflt(:),:).*S_l1(:,Dflt(:),:)).*(a_ij(:,Dflt(:),:)./(D_k(:,:,:))).*Factor_i(:,Dflt(:),:);
    p1common2a(Dflt(:),:,:)=abs(Dummyi(Dflt(:),:,:).*S_l1(Dflt(:),:,:)).*(a_ij(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3])).*Factor_i(Dflt(:),:,:);  
end
p1common2a=sum(p1common2a(Dflt(:),:,:),2); %sums across rows which are j's hear //holding default with each cp
if D>=2
    p1common2a=sum(p1common2a); %sums across colums //holdings with all defaulters and cps
end
p1common2a=repmat(p1common2a,I); %replicates to be size of other matrices

%From Bank i in predation
p1common2b=z; %initiates zero matrix
p1common2b=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
p1common2b(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
p1common2b(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
dummy=repmat(D_k(:,1,:),1,I);
if l<2
    p1common2b=p1common2b.*abs(Dummyj(:,:,:).*S_l(:,:,:,l-1)).*(a_ji(:,:,:)./dummy).*Factor_j(:,:,:);
else
    p1common2b=p1common2b.*abs(Dummyj(:,:,:).*S_l1(:,:,:)).*(a_ji(:,:,:)./dummy).*Factor_j(:,:,:);
end
p1common2b=sum(p1common2b(:,:,:),2);
p1common2b=repmat(p1common2b,1,I);
if D>=2
    p1common2b=sum(p1common2b);
    p1common2b=repmat(p1common2b,I,1);
end
                
%Put together for the final term
            P1_i=(p11+p12).*p1common;
            P1_j=permute(P1_i,[2,1,3]);
            P12_i=p11.*p1common2a+p12.*p1common2b;
            P12_j=permute(P12_i,[2,1,3]);
        
%---------------------------------------------------------------------------

%PP: Get the Predatory Price Impact
pp1=z; %initiates zero matrix
pp1=Dummyi.*Factor_i.*X_ij;

%Multiplier from others predation liquidation 
%over non-defaulted banks holding defaulted asset
ppcommon=z;
dummy=repmat(D_k(:,1,:),1,I);
if l<2
    ppcommon(:,:,:)=abs(Dummyj(:,:,:).*S_l(:,:,:,l-1)).*(a_ji(:,:,:)./dummy).*Factor_j(:,:,:);
else
    ppcommon(:,:,:)=abs(Dummyj(:,:,:).*S_l1(:,:,:)).*(a_ji(:,:,:)./dummy).*Factor_j(:,:,:);
end
ppcommon(Dflt(:),:,:)=0;
ppcommon(:,Dflt(:),:)=0;
ppcommon=sum(ppcommon(:,:,:),2);
%if D>=2
    ppcommon=sum(ppcommon);
%end 
ppcommon=repmat(ppcommon,I);
                    
%Put together for the final term                
            PP_i=(pp1).*ppcommon;
            PP_j=permute(PP_i,[2,1,3]);
    
%---------------------------------------------------------------------------
%P2: Get the Secondary Price Impact
P2_i=z;
P2_j=z;
P3_i=z;
P3_j=z;

if l>2
    
%Over the defaulted banks j
p21=z;
Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
p21(:,Dflt(:),:)=(3/factorial(2)).*Factor.*X_ij(:,Dflt(:),:); %gets the results for defaulted banks' counterparties
p21(Dflt(:),:,:)=(3/factorial(2)).*permute(Factor,[2,1,3]).*X_ij(Dflt(:),:,:); %gets a result for defaulted banks

%Over the undefaulted banks j
p22=z;
p22=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
p22(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
p22(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
p22=p22.*Factor.*X_ij; %gets result for undefaulted banks
            
%Multiplier over the other assets with all banks
p2common=z;
dummy=repmat(D_k(:,1,:),1,I);
Factor_i=(X_ij(:,:,:)~=0);
%Factor_i2=circshift(Factor_i,1,3)
Factor_i2=sum(Factor_i,2);
Factor_i2=sum(Factor_i2,3);
p2common=abs(Dummyj(:,:,:).*S_l2).*(a_ji(:,:,:)./dummy).*Factor_j(:,:,:);
p2common=sum(p2common(:,:,:),2);
p2common=sum(p2common(:,:,:),1);
p2common=sum(p2common(:,:,:),3);
p2common=repmat(p2common,J,1);
p2common=Factor_i2.*p2common;
%if D>=2
    p2common=sum(p2common);
%end 
p2common=repmat(p2common,I,J,K);

                P2_i=(1/factorial(2)).*(p21+p22).*p2common;
                P2_j=permute(P2_i,[2,1,3]);  
end
   
if l>2
%P3: Get the Tertiary Price Impact
%Over the defaulted banks j
p31=z;
Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
p31(:,Dflt(:),:)=(9/factorial(3)).*Factor.*X_ij(:,Dflt(:),:); %gets the results for defaulted banks' counterparties
p31(Dflt(:),:,:)=(9/factorial(3)).*permute(Factor,[2,1,3]).*X_ij(Dflt(:),:,:); %gets a result for defaulted banks
Factor31=(p31~=0);
Nfactor31=(p31==0);
Nfactor31=sum(Nfactor31,3);
Nfactor31=repmat(Nfactor31,1,1,K);

%Over defaulted multiplier
p31common=abs(Factor31(:,:,:).*S_l2).*(a_ji(:,:,:)./dummy).*Factor31(:,:,:);
p31common=sum(p31common(:,:,:),2);
p31common=sum(p31common(:,:,:),3);
p31common=repmat(p31common,1,J,K);               

%Over the undefaulted banks j
p32=z;  
p32=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
p32(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
p32(:,Dflt(:),:)=0; %defaulted banks

%Over defaulted multiplier
p31common=abs(Factor31(:,:,:).*S_l2).*(a_ji(:,:,:)./dummy).*Factor31(:,:,:);
p31common=sum(p31common(:,:,:),2);
p31common=sum(p31common(:,:,:),3);
p31common=repmat(p31common,1,J,K);               

%Over the undefaulted banks j
p32=z;  
p32=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
p32(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
p32(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
p32=p32.*Factor.*X_ij; %gets result for undefaulted banks
Factor32=(p32~=0);
Nfactor32=(p32==0);
Nfactor32=sum(Nfactor32,3);
Nfactor32=repmat(Nfactor32,1,1,K);

%Over undefaulted multiplier
p32common=abs(Factor32(:,:,:).*S_l2).*(a_ji(:,:,:)./dummy).*Factor32(:,:,:);
p32common=sum(p32common(:,:,:),2);
p32common=sum(p32common(:,:,:),3);
p32common=repmat(p32common,1,J,K);  

%Summation 
P3_i=(1/factorial(3)).*(p31.*Nfactor31.*p31common+p32.*Nfactor32.*p32common);
P3_j=permute(P3_i,[2,1,3]);
end


%OBTAIN LIABILITY OF THE BANK AND CCP
%--------------------------------------
%Gather all P_ij over all i (the sum function goes down the column)
 P0_ji=P0_j(1,:,:);  %I think i need to do it by ji and ij rather than k
 P0_ij=P0_i(1,:,:);
 %---------------
 P1_ji=P1_j(1,:,:);
 P1_ij=P1_i(1,:,:);
  %---------------
 P12_ji=P12_j(1,:,:);
 P12_ij=P12_i(1,:,:);
  %---------------
 PP_ji=PP_j(1,:,:);
 PP_ij=PP_i(1,:,:);
  %---------------
 P2_ji=P2_j(1,:,:);
 P2_ij=P2_i(1,:,:);
  %---------------
 P3_ji=P3_j(1,:,:);
 P3_ij=P3_i(1,:,:);
 
 %First sum over all assets for one k the assets
 P_ij=[P0_ij;P1_ij;P12_ij;PP_ij;P2_ij;P3_ij];
 P_ji=[P0_ji;P1_ji;P12_ji;PP_ji;P2_ji;P3_ij];
 P_ij=sum(P_ij);
 P_ji=sum(P_ji);
 
 P_ij=P_ij(1,:,:);
 P_ji=P_ji(1,:,:);
 
%Sum for each i over the k
 P_IJ=sum(P_ij,3);
 P_JI=sum(P_ji,3);
 
 C=[P_IJ;P_JI];
end
