function Predator_k=predator(Dist,Pred,Distress_k,Dflt,X_ij,I,J,K)
    D=length(Dflt);
    Dummy=(Distress_k==0); %gets set of distressed guys
    Dummy1=[1:I]'; % gets set of banks
    Dummy1=repmat(Dummy1, [1,1,K]);
    Dummy2=Dummy1.*Dummy; %puts them together
    Dummy2(Dflt(1:D),:,:)=0; %takes out the defaulters
    Predator_k=Dummy2; %sets this as predator
    
    dummy=transpose([1:I]); %dummy of banks
    msize = numel(dummy); %get size
    idx = randperm(msize); %permute the size randomly
    for k=1:K
        if Pred<sum(Predator_k(:,:,k)~=0) %if distress is smaller
           Diff=sum(Predator_k(:,:,k)~=0)-Dist; %get the difference
            while Pred~=sum(Predator_k(:,:,k)~=0) %while the diff is non zero
                idx = randperm(msize);
                idx(1);      
                Predator_k(idx(1),:,k)=0;  %keep knocking out pieces
            end
        end
        % Pred=sum(Predator_k(:,:,k)~=0) %no. distress is equal to full
        % matrix but not needed
    end
end