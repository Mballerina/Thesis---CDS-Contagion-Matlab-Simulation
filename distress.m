function Distress_k=distress(Dist,Dflt,X_ij,I,J,K)
    D=length(Dflt); %determines the current length of default vector
    Dummy1=(X_ij(:,Dflt(1:D),:)~=0); %binary matrix; gives nonempty holdings with defaulter (distressed banks)
    size0=transpose(repmat([1:I],D,1)); %build index matrix of bank numbers
    size=cat(3,size0,size0); %concatenate so that array is over k assets
    for k=1:K-2
        size=cat(3,size,size0); %repeat this process K-1 times in order to get dimension 3 to equal K
    end
    
    %REMOVING THE MULTIPLES IN THE DISTRESSED GROUP
    size1=size; %set this matrix in another variable to manipulate it
    size1(Dflt(1:D),:,:)=0; %knock out all defaulted banks which cannot be distressed
    Dummy1=(Dummy1.*size1); %turn binary matrix into an indexed matrix and collapse columns
    half=(sum(Dummy1,2)); %sums the columns 
    s=transpose(1:I);
    increment=half./repmat(s,[1,1,K]); %determines how much bigger than index with multiplier
    minus=(increment==0); %binary matrix to turn non-zero entries to one for division
    increment=(increment+minus); %adjust increment for division
    half=half./increment; %divides to get proper index
    Distress_k=half; % gives one list of all possible distressed banks
    %Distress_k=reshape(Distress_k, [], K); %gives a compression
    
    %MAKE SURE THE ARRAY OF DISTRESSED BANKS NOT LARGER THAN CHOSEN BY USER
    %If a number is chosen for Dist (otherwise given empty case given in main):
    dummy=transpose([1:I]); %dummy of bank numbers
    msize = numel(dummy); %get size of the dummy
    idx = randperm(msize); %permute the numbers of the banks randomly
    for k=1:K %for each asset, cycles over dimension 3
        if Dist<sum(Distress_k(:,:,k)~=0) %if no. of distressed banks is smaller than those array of distressed banks
           Diff=sum(Distress_k(:,:,k)~=0)-Dist; %get the difference between the current set and the user chosen set
            while Dist<sum(Distress_k(:,:,k)~=0) %while the diff is non-zero
                idx = randperm(msize); %shuffles the bank number
                %idx(1); %first bank is picked out from random shuffle
                Distress_k(idx(1),:,k)=0;  %sets entry to zero to until narrows done number of banks user chose.
            end
        end
        % Dist=sum(Distress_k(:,:,k)~=0) %no. distress is equal to full
        % matrix not needed
    end
    
end