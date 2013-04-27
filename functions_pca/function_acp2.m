function [Y,U,x_mean,x_std]=function_acp2(x,y)
% 
% ENTREES
% x         tableau des donnee (individus,caracteres)
% y         classe des individus


marker_color = {'bx','r.','g*','ksquare'};

[m,n] = size(x);

% centrage et normalisation du nuage de points
x_mean = mean(x,1);
x_std = std(x,1,1);
x_centre = (x - x_mean(ones(m,1),:))./x_std(ones(m,1),:);

% %calcul dela matrice de variance covariance
xcov=cov(x_centre,1);

% %calcul des valeurs propres et vecteurs propres
[v,e]=eig(xcov);
% %valeures propres mises dans l'ordre decroissant
ed=diag(e);
[lambda, rg] = sort(ed,1,'descend');
U = v(:,rg);
Y=x_centre*U;
Inertie = lambda/sum(lambda)*100;
% disp('les 7 premieres valeurs propres en pourcent')
% disp(Inertie(1:min([7 n]))')

axe1 = 1;
axe2 = 2;
% figure;
% axis([min(Y(:,axe1)) max(Y(:,axe1)) min(Y(:,axe2)), max(Y(:,axe2))]); 
% hold on;
% for i=1:max(y)
%     plot(Y(y==i,axe1),Y(y==i,axe2),marker_color{i},'LineWidth',1.5)
% end
% hold off;


