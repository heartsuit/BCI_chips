function [xcentrenorme, lambda,coscarre,U,Y,qualvar,avg,std_x]=facp09(x,xnom,classe,nomvar,op)
% facp(x,xnom)
% 
% ENTREES
% x         tableau des donnee (individus,caracteres)
% xnom		vecteur colonne des noms des individus
% classe	vecteur (1xn)representant le nombre d'individu des n classes (couleur)
%           les individus ayant ete classe dans l'ordre des classes dans
%           x et xnom
%           ex : classe=(1,2,5)	il y a 3 classes avec 1 individu dans la classe 
%				1, 2 dans la classe 2 et 5 dans la 3
% op		option connection des individus d'un meme couleur entre eux
%           par defaut ils ne le sont pas
%
% SORTIES
% xcentrenorme : donnees centrees normees
% lambda		valeurs propres
% coscarre		matrice des cosinus2 (qualité de projection des individus sur les axes
% U             matrice des vecteurs propres
% Y             Y = xcentrenorme * U (donnees projetees
% qualvar       correlation des variables initiales avec les composantes principales
%

[m,n] = size(x);
if nargin==5, op=1;end
if nargin<5, op=2;end
if nargin<3, classe=m;end
if nargin<2, classe=ones(1,m);end

% centrage et normalisation du nuage de points
avg = mean(x);
xcentrenorme = zeros(size(x));
for j=1:n
  xcentrenorme(:,j) = (x(:,j)-avg(j));
  std_x = std(xcentrenorme(:,j),1);
  xcentrenorme(:,j) = xcentrenorme(:,j)/std(xcentrenorme(:,j),1);
end
%calcul dela matrice de variance covariance
xcov=cov(xcentrenorme,1);

%calcul des valeurs propres et vecteurs propres
[v,e]=eig(xcov);
%valeures propres mises dans l'ordre decroissant
ed=diag(e);
j=1;
while j<=n,
	ma=max(ed);
	i=1;
	while ed(i)~=ma
		i=i+1;		
	end
	rg(j)=i;
	ed(i)=0;
	j=j+1;
end	
ed=diag(e);
U=zeros(size(v));
for i=1:n,
	lambda(i)=ed(rg(i));
	U(:,i)= v(:,rg(i));
end
% lambda
Inertie = lambda/sum(lambda)*100;
disp('les 7 premieres valeurs propres en pourcent')
Inertie(1:min([7 n]))
% figure; plot([1:size(lambda,2)],lambda,'+-','LineWidth',3); title('valeurs propres')

%calcul du nuage de points dans la base des vecteurs propres
Y=xcentrenorme*U;

%calcul de la qualite ponctuelle des individus axe par axe
coscarre=zeros(size(x));
for i=1:m,
	coscarre(i,:)=Y(i,:).^2/norm(xcentrenorme(i,:))^2;
end
%disp('matrice des cosinus limitee a 7 colonnes')
coscarre(:,1:min([7 n]));
%
% qualité de représentation des variables
% utiliseŽ pour la reprŽsentation des variables dans le cercle des
% corrŽlations.
for k=1:n
   for j=1:n
      qualvar(j,k)=sqrt(lambda(k))*U(j,k);
   end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   representations graphiques                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=2;
choixaxe=1;
axe1=1;
axe2=2;

while choixaxe==1,
% 	figure(f),close(f),figure(f)
    figure;
	axis([min(Y(:,axe1)) max(Y(:,axe1)) min(Y(:,axe2)),...
   max(Y(:,axe2))]) 
	hold on
	if classe~=0,
		b=0;
      fx=['ko';'ro';'b+';'g*';'ms';'yo';'wo'];
      fx2=['r+';'g+';'b+';'r*';'g*';'b*';'m*'];

%		fx=['go';'c+';'yx';'yo';'y*';'go';'bo'];
		for i=1:length(classe),
			a=b+1;
			j=classe(i);
			b=b+j;
			plot([Y(a:b,axe1)' Y(a,axe1)]',...
			[Y(a:b,axe2)' Y(a,axe2)]',fx(rem(i,7)+1,1:2),'LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% BOUCLE A NE PAS EXECUTER SI ON NE VEUT PAS LES NOMS DES INDIVIDUS
%         for k=a:b,
% 			       text(Y(k,axe1),Y(k,axe2),['  ' int2str(xnom(k))])
% 			end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
	else
		plot(Y(:,axe1),Y(:,axe2),'o')
   end
% 	line([min(Y(:,axe1)) max(Y(:,axe1))],[0 0])
% 	line([0 0],[min(Y(:,axe2)) max(Y(:,axe2))])
%    title(['axe ' int2str(axe1) ' / axe ' int2str(axe2)]),grid on
   xlabel([num2str(round(Inertie(axe1))) ' %']);
   ylabel([num2str(round(Inertie(axe2))) ' %']);
% 
%    % cercle des corrélations
%    figure;
%    teta=[0:2*pi/200:2*pi];
%    xx=cos(teta);yy=sin(teta);
%    hold on
%    plot(xx,yy);
%    for j=1:n
%             compass(qualvar(j,axe1),qualvar(j,axe2),'-');
%             h1=text(qualvar(j,axe1),qualvar(j,axe2),['  ' nomvar(j,:)]);		
%             set(h1,'fontsize',10)%,'color','j');
% 
%    end;
%    title(['axe ' int2str(axe1) ' / axe ' int2str(axe2)]),grid on
%       xlabel([num2str(round(Inertie(axe1))) ' %']);
%    ylabel([num2str(round(Inertie(axe2))) ' %']);
% 
%    hold off
% 	choixaxe=menu('visualisation d''un autre plan','oui','non');
    choixaxe = 0;
% 	if choixaxe==1
% 		axe1=input('axe 1: ');
% 		axe2=input('axe 2: ');
% 		f=f+2;
% 	end
end

