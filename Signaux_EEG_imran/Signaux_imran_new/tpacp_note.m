% tpacp_note
% exemple programme d'appel de la fonction ACP
% ici la matrice donnees est de taille 52 (individus) x 6 (variables)
% penser à fermer la fenêtre MENU (en répondant NON) quand on veut terminer le programme

clear
close all

load donnees; % matrice individus x variables
[n,m] = size(donnees);
disp('taille des donnees')
[n,m]
nomindividu = (1:1:n)'; %nom des individus 
% on visualise la population avec 4 couleurs (classes) ou 1 couleur
classe=[13 13 13 13]; % nombre d'individus par couleur à visualiser
% classe=[52]; % nombre d'individus par couleur à visualiser

nomvariable = [
'  Note   '
'  Taille '
'  Poids  '
'  Alea   '
'  Sport1 '
'  Sport2 '
]; 
% appel de la fonction acp
[X, lambda, coscarre, U, Y, qualvar] = facp09(donnees,nomindividu,classe,nomvariable,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable supplementaire a traiter
variablesup = ['  Genre   '];
Varsup=-ones(n,1); Varsup(1:n/2,1) = 1; % genre feminin = 1 ; masculin = -1;

