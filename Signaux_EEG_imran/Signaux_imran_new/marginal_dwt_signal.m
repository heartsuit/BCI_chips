function e = marginal_dwt_signal (x, h, deb)
%
% function e = marginal_dwt_signal (x, h, deb);
%
% marginales de abs(DWT) d'un signal x
% la DWT est calculée par le filtre h
%
% ENTREE
% x  : signal
% h : filtre
% deb : nombre de niveaux profonds negliges pour le calcul de l'entropie
% SORTIE
% e : marginales de la DWT du signal x
N = length (x);
NL = fix(log2(N)); 

% à la mani?re de WAVELET toolbox
mode = 'per';
g = h(end:-1:1);
for ii=1:length(g)
    g(ii) = (-1)^ii * g(ii); 
end
[y,l] = wavedec_simplif(x,mode,NL-deb,h,g,0,0);
% 
% % RICE toolbox
% y=mdwt (x, h, NL - deb); 

% % wavelab toolbox
% L = deb;
% y = FWT_PO(x,L,h);  
% 
e = marginal(abs(y),deb);
