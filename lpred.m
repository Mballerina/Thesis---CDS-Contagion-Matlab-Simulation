function LPred=lpred(Predator_k,HGamma_i,X_ij,I,J,K)
    %choose largest predator with highest holding and most income
for i=1:I
    X_i=sum(X_ij,2); %sums rows for i to get largest holding
    dummy=(Predator_k~=0); %picks out those that are predators
    hdummy1=X_i.*dummy; %picks out relevant holdings
    hdummy2=repmat(HGamma_i,[1,1,K]).*dummy; %does the same for cash
    fdummy=hdummy1+hdummy2; %adds them together
        [Max,LPred]=max(abs(fdummy)); %takes on first indices which satisfies max
    for k=1:K
        if Max(1,1,k)==0
            LPred(1,1,k)=0; %restores zero to the party that had zero holdings ie. no preds
        end
    end
end

end