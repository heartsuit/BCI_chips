function [h,g]=cons_wav_filter(longueur,par)
%
%h=construction_filtreh(longueur,par);
%
%ENTREE
%%%%%%%
% longueur	: taille du filtre (4 ou 6)
% par			: parametre libre, si vide tiré au hasard
%
%SORTIE
%%%%%%%
% h			: filtre h
%

% Tirage aleatoire des parametres libres
if nargin<2,
   par =rand(1,3)*2*pi;
end

if longueur==4,
   a=par(1);  
   h(1)=(1-cos(a)+sin(a))/(2*sqrt(2));
   h(2)=(1+cos(a)+sin(a))/(2*sqrt(2));
   h(3)=(1+cos(a)-sin(a))/(2*sqrt(2));
   h(4)=(1-cos(a)-sin(a))/(2*sqrt(2));
elseif longueur==6,
   a=par(1);
   b=par(2);
   h(1)=( (1+cos(a)+sin(a))*(1-cos(b)-sin(b))+2*sin(b)*cos(a) )/(4*sqrt(2));
   h(2)=( (1-cos(a)+sin(a))*(1+cos(b)-sin(b))-2*sin(b)*cos(a) )/(4*sqrt(2));
   h(3)=(1+cos(a-b)+sin(a-b))/(2*sqrt(2));
   h(4)=(1+cos(a-b)-sin(a-b))/(2*sqrt(2));
   h(5)=(1/sqrt(2))-h(1)-h(3);
   h(6)=(1/sqrt(2))-h(2)-h(4);
elseif longueur==8, 
   t0 = par(1);
   t1 = par(2);
   t2 = par(3);
   t3 = 1/4*pi-t0-t1-t2;
   h(1) = cos(t3)*cos(t2)*cos(t1)*cos(t0);
   h(2) = cos(t3)*cos(t2)*cos(t1)*sin(t0);
   h(3) = -cos(t3)*cos(t2)*sin(t1)*sin(t0)-cos(t3)*sin(t2)*sin(t1)*cos(t0)-...
       sin(t3)*sin(t2)*cos(t1)*cos(t0);
   h(4) = cos(t3)*cos(t2)*sin(t1)*cos(t0)-cos(t3)*sin(t2)*sin(t1)*sin(t0)-...
       sin(t3)*sin(t2)*cos(t1)*sin(t0);
   h(5) = -cos(t3)*sin(t2)*cos(t1)*sin(t0)+sin(t3)*sin(t2)*sin(t1)*sin(t0)-...
       sin(t3)*cos(t2)*sin(t1)*cos(t0);
   h(6) = cos(t3)*sin(t2)*cos(t1)*cos(t0)-sin(t3)*sin(t2)*sin(t1)*cos(t0)-...
       sin(t3)*cos(t2)*sin(t1)*sin(t0);
   h(7) = -sin(t3)*cos(t2)*cos(t1)*sin(t0);
   h(8) = sin(t3)*cos(t2)*cos(t1)*cos(t0);
end


g = h(end:-1:1);
for i=1:length(g)
    g(i) = (-1)^i * g(i); 
end


