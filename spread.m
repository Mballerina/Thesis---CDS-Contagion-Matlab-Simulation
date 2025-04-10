function B=spread(l,D_k,Dflt,S_l,X_ij,X_ji,a_ij,a_ji,I,J,K)

%Pricing Functionals 
z=zeros(I,J,K);
%--------------

%SP_0: Get the fundamental value of postions
Dummyi=(X_ij~=0);
Dummyj=(X_ji~=0);
%SP0_i=Dummyi.*S_l(:,:,:,l-1);
SP0_i=S_l(:,:,:,l-1);
%SP0_i=max(SP0_i,zeros(I,J,K));
%SP0_j=Dummyj.*S_l(:,:,:,l-1);
SP0_j=S_l(:,:,:,l-1);
%SP0_j=max(SP0_j,zeros(I,J,K));
%----------------------------------------------------------------------

D=length(Dflt); %determines the current length of default vector
D_k=D_k(1,1,1); %picks out one entry of market liquidity
D_k=D_k*ones(I,D,K); %makes a liquidity/market depth matrix of the size of number of defaulters
%P_1: Get the Primary Price Impact [KEY: Separate if statements]
epsilon=.98; %let this be user input
%sum over defaulted banks
%sp11=z; %initiates a zero matrix
%Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
%sp11(:,Dflt(1:D),:)=Factor.*X_ij(:,Dflt(1:D),:); %gets the results for defaulted banks' counterparties
%sp11(Dflt(1:D),:,:)=permute(Factor,[2,1,3]).*X_ij(Dflt(1:D),:,:); %gets a result for defaulted banks

%sum over undefaulted banks
%sp12=z; %initiates zero matrix
%sp12=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
%sp12(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
%sp12(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
%Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
%sp12=sp12.*epsilon.*Factor.*X_ij; %gets result for undefaulted banks

%Multiplier from CCP liquidation over defaulted banks
sp1common=z; %initiates zero matrix
Factor_j=(X_ji(:,:,:)~=0);
Inverse_j=Factor_j.*(1/X_ji);
Inverse_j(isnan(Inverse_j))=0;
sp1common(:,Dflt(:),:)=abs(Factor_j(:,Dflt(:),:).*S_l(:,Dflt(:),:,l-1)).*(a_ji(:,Dflt(:),:)./(D_k(:,:,:)));%.*Inverse_j(:,Dflt(:),:);
sp1common(Dflt(:),:,:)=abs(Factor_j(Dflt(:),:,:).*S_l(Dflt(:),:,:,l-1)).*(a_ji(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3]));%.*Inverse_j(Dflt(:),:,:);                                    
sp1common=sum(sp1common(Dflt(:),:,:),2);
if D>=2
    sp1common=sum(sp1common);
end 
sp1common=repmat(sp1common,I);


%Multiplier from Bank i in distress
sp1common2a=z; %initiates zero matrix
Factor_i=(X_ij(:,:,:)~=0);
Inverse_i=Factor_i.*(1/X_ij);
Inverse_i(isnan(Inverse_i))=0;
sp1common2a(:,Dflt(:),:)=abs(Factor_i(:,Dflt(:),:).*S_l(:,Dflt(:),:,l-1)).*(a_ij(:,Dflt(:),:)./(D_k(:,:,:)));%.*Inverse_i(:,Dflt(:),:);
sp1common2a(Dflt(:),:,:)=abs(Factor_i(Dflt(:),:,:).*S_l(Dflt(:),:,:,l-1)).*(a_ij(Dflt(:),:,:)./permute(D_k(:,:,:),[2,1,3]));%.*Inverse_i(Dflt(:),:,:);
sp1common2a=sum(sp1common2a(Dflt(:),:,:),2); %sums across rows which are j's hear //holding default with each cp
if D>=2
    sp1common2a=sum(sp1common2a); %sums across colums //holdings with all defaulters and cps
end
sp1common2a=repmat(sp1common2a,I); %replicates to be size of other matrices


%From Bank i in predation
sp1common2b=z; %initiates zero matrix
sp1common2b=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
sp1common2b(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
sp1common2b(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
dummy=repmat(D_k(:,1,:),1,I);
sp1common2b=sp1common2b.*abs(Factor_j(:,:,:).*S_l(:,:,:,l-1)).*(a_ji(:,:,:)./dummy);%.*Inverse_j(:,:,:);
sp1common2b=sum(sp1common2b(:,:,:),2);
sp1common2b=repmat(sp1common2b,1,I);
if D>=2
    sp1common2b=sum(sp1common2b);
     sp1common2b=repmat(sp1common2b,I,1);
end


                
%Put together for the final term
            SP1_i=sp1common;
            SP1_j=permute(SP1_i,[2,1,3]);
            SP12_i=sp1common2a+sp1common2b;
            SP12_j=permute(SP12_i,[2,1,3]);
        
%---------------------------------------------------------------------------

%PP: Get the Predatory Price Impact
%spp1=z; %initiates zero matrix
%spp1=Dummyi.*Factor_i.*X_ij;

%Multiplier from others predation liquidation 
%over non-defaulted banks holding defaulted asset
sppcommon=z;
dummy=repmat(D_k(:,1,:),1,I);
sppcommon(:,:,:)=abs(Factor_j(:,:,:).*S_l(:,:,:,l-1)).*(a_ji(:,:,:)./dummy);%.*Inverse_j(:,:,:);
sppcommon(Dflt(:),:,:)=0;
sppcommon(:,Dflt(:),:)=0;
sppcommon=sum(sppcommon(:,:,:),2);
%if D>=2
    sppcommon=sum(sppcommon);
%end 
sppcommon=repmat(sppcommon,I);
                    
%Put together for the final term                
            SPP_i=sppcommon;
            SPP_j=permute(SPP_i,[2,1,3]);
    
%---------------------------------------------------------------------------
%P2: Get the Secondary Price Impact
SP2_i=z;
SP2_j=z;
SP3_i=z;
SP3_j=z;

if l>2
    
%Over the defaulted banks j
%p21=z;
%Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
%p21(:,Dflt(:),:)=(3/factorial(2)).*Factor.*X_ij(:,Dflt(:),:); %gets the results for defaulted banks' counterparties
%p21(Dflt(:),:,:)=(3/factorial(2)).*permute(Factor,[2,1,3]).*X_ij(Dflt(:),:,:); %gets a result for defaulted banks

%Over the undefaulted banks j
%p22=z;
%p22=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
%p22(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
%p22(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
%Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
%p22=p22.*Factor.*X_ij; %gets result for undefaulted banks
            
%Multiplier over the other assets with all banks
sp2common=z;
dummy=repmat(D_k(:,1,:),1,I);
Factor_i=(X_ij(:,:,:)~=0);
%Factor_i2=circshift(Factor_i,1,3)
%Factor_i2=sum(Factor_i,2);
%Factor_i2=sum(Factor_i2,3);
sp2common=abs(Factor_j(:,:,:).*S_l(:,:,:,l-2)).*(a_ji(:,:,:)./dummy);%.*Inverse_j(:,:,:);
sp2common=sum(sp2common(:,:,:),2);
sp2common=sum(sp2common(:,:,:),1);
sp2common=sum(sp2common(:,:,:),3);
sp2common=repmat(sp2common,J,1);
%sp2common=Factor_i2.*sp2common;
%if D>=2
    sp2common=sum(sp2common);
%end 
sp2common=repmat(sp2common,I,J,K);

                SP2_i=(1/factorial(2)).*sp2common;
                SP2_j=permute(SP2_i,[2,1,3]);  
end
   
if l>2
%P3: Get the Tertiary Price Impact
%Over the defaulted banks j
%sp31=z;
%Factor=(X_ij(:,Dflt(1:D),:)~=0); %gets a factor for banks holding the asset
%sp31(:,Dflt(:),:)=(9/factorial(3)).*Factor.*X_ij(:,Dflt(:),:); %gets the results for defaulted banks' counterparties
%sp31(Dflt(:),:,:)=(9/factorial(3)).*permute(Factor,[2,1,3]).*X_ij(Dflt(:),:,:); gets a result for defaulted banks
%Factor31=(sp31~=0);
%Nfactor31=(sp31==0);
%Nfactor31=sum(Nfactor31,3);
%Nfactor31=repmat(Nfactor31,1,1,K);

%Over defaulted multiplier
Factor_j=z; %initiates zero matrix
Factor=(X_ji(Dflt(:),:,:)~=0);
Factor_j(Dflt(:),:,:)=Factor;
Inverse_j=Factor_j.*(1/X_ji);
Inverse_j(isnan(Inverse_j))=0;

sp31common=abs(Factor_j(:,:,:).*S_l(:,:,:,l-2)).*(a_ji(:,:,:)./dummy);%.*Inverse_j(:,:,:);
sp31common=sum(sp31common(:,:,:),2);
sp31common=sum(sp31common(:,:,:),3);
sp31common=repmat(sp31common,1,J,K);               

%Over the undefaulted banks j
%sp32=z;  
%sp32=(X_ij(:,:,:)~=0); %gets matrix of 1's for held assets
%sp32(Dflt(:),:,:)=0; %defaulted banks' holdings are 0
%sp32(:,Dflt(:),:)=0; %defaulted banks' counterparties holdings are zero
%Factor=(X_ij(:,:,:)~=0);  %gets factor of ones for held assets
%sp32=sp32.*Factor.*X_ij; %gets result for undefaulted banks
%Factor32=(sp32~=0);
%Nfactor32=(sp32==0);
%Nfactor32=sum(Nfactor32,3);
%Nfactor32=repmat(Nfactor32,1,1,K);

%Over undefaulted multiplier
Factor_j=z; %initiates zero matrix
Factor_j=(X_ji~=0);
Factor_j(Dflt(:),:,:)=0;
Inverse_j=Factor_j.*(1/X_ji);
Inverse_j(isnan(Inverse_j))=0;

sp32common=abs(Factor_j(:,:,:).*S_l(:,:,:,l-1)).*(a_ji(:,:,:)./dummy);%.*Inverse_j(:,:,:);
sp32common=sum(sp32common(:,:,:),2);
sp32common=sum(sp32common(:,:,:),3);
sp32common=repmat(sp32common,1,J,K);  

%Summation 
SP3_i=(1/factorial(3)).*(sp31common+sp32common);
SP3_j=permute(SP3_i,[2,1,3]);
end


%OBTAIN LIABILITY OF THE BANK AND CCP
%--------------------------------------
%Gather all P_ij over all i (the sum function goes down the column)
 SP0_ji=SP0_j(1,:,:);  %I think i need to do it by ji and ij rather than k
 SP0_ij=SP0_i(1,:,:);
 %---------------
 SP1_ji=SP1_j(1,:,:);
 SP1_ij=SP1_i(1,:,:);
  %---------------
 SP12_ji=SP12_j(1,:,:);
 SP12_ij=SP12_i(1,:,:);
  %---------------
 SPP_ji=SPP_j(1,:,:);
 SPP_ij=SPP_i(1,:,:);
  %---------------
 SP2_ji=SP2_j(1,:,:);
 SP2_ij=SP2_i(1,:,:);
  %---------------
 SP3_ji=SP3_j(1,:,:);
 SP3_ij=SP3_i(1,:,:);
 
 %First sum over all assets for one k the assets
 SP_ij=[SP0_ij;SP1_ij;SP12_ij;SPP_ij;SP2_ij;SP3_ij];
 SP_ji=[SP0_ji;SP1_ji;SP12_ji;SPP_ji;SP2_ji;SP3_ij];
 SP_ij=sum(SP_ij);
 SP_ji=sum(SP_ji);
 
 SP_ij=SP_ij(1,:,:);
 SP_ji=SP_ji(1,:,:);
 
%Sum for each i over the k
 %SP_IJ=sum(SP_ij,3);
 %SP_JI=sum(SP_ji,3);
 
 B=[SP_ij;SP_ji];
end
