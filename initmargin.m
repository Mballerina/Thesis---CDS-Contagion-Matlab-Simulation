function g_i=initmargin(margin,l,S_l,X_ij)
         g_i=sum((margin.*sum(abs(S_l(:,:,:,l).*X_ij),3)),2);
end