function e = marginal_dwt_popu (x, h, lh_filtre, deb)
%
% function e = marginal_dwt_popu (x, h, lh_filtre, deb);
%
% marginales des DWT d'une population x de plusieurs classes de signaux
%
% ENTREE
% x (long_sig, nbsigclas, nbclas) : ensemble des signaux de la population d'apprentissage
% h : filtre pour le calcul de la DWT
% lh : longueurs de h
% deb : nombre de niveaux profonds negliges pour le calcul de la dwt

% SORTIE
% e (nbniv, nbclas, nbsigclas) : marginales des DWT des signaux de x

nbsig = size(x,1);

h=h(1:lh_filtre);
for k = 1:nbsig,
     z = x (k,:);
     e_z = marginal_dwt_signal (z, h, deb);
     e (k,:) = e_z;
end
