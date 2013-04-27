function class = marg_dwt_concat_voies(class,nbf,h,deb,chan)
% nbf=nbf+1;
nbclasses = size(class,2);
nbchan = length(chan);
for i=1:nbclasses
    clear xK;
    for K=1:nbchan
        xK(:,:) = class(i).signal_tronc(:,:,chan(K));
        hK = h(K,:); 
        longueur_h = length(hK);
        class(i).des(:,(K-1)*nbf+1:nbf*K) = marginal_dwt_popu(xK, hK, longueur_h, deb);
    end
end