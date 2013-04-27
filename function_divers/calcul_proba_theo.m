function [Err_glob, P_repeat, BpT_final] = calcul_proba_theo(E_SVM,P_C_C,P_W_W,P_rep,n)
% calcul_proba_theo(E_SVM,P_C_C,P_W_W,P_rep,n)
% E_SVM : Probabilité d'erreur du SVM (Valeur entre [0 1])
% P_C_C : Probabilité que le détecteur dise juste sachant que la réponse
% est juste (Valeur entre [0 1])
% P_W_W : Probabilité que le détecteur dise faux sachant que la réponse
% est faux (Valeur entre [0 1])
% P_rep : Probabilité de réponse du détecteur. Pour l'approche classif 
% P_rep = 1 (Valeur entre [0 1])
% n : nombre de classes

P_W_C = 1-P_C_C;
P_C_W = 1-P_W_W;

P_W = E_SVM;
P_C = 1-E_SVM;

P = 1/n;

P_post_W_W = (P_W_W*P_W*P)/(P_C*P_W_C*P + P_W*P_W_W*P);
P_post_C_C = (P_C_C*P_C*P)/(P_C*P_C_C*P + P_W*P_C_W*P);
P_post_C_W = (P_W_C*P_C*P)/(P_C*P_W_C*P + P_W*P_W_W*P);
P_post_W_C = (P_C_W*P_W*P)/(P_C*P_C_C*P + P_W*P_C_W*P);

disp(['P_post_wrong_estim_wrong :' num2str(P_post_W_W)]);
disp(['P_post_correct_estim_correct :' num2str(P_post_C_C)]);
disp(['P_post_correct_estim_wrong :' num2str(P_post_C_W)]);
disp(['P_post_wrong_estim_correct :' num2str(P_post_W_C)]);
disp('--------------------------------------------------');
P_estim_C = P_C_W*P_W + P_C_C*P_C;
P_estim_W = P_W_W*P_W + P_W_C*P_C;

Err_glob = ((1-P_rep)*P_W+P_rep*P_C_W)/((1-P_rep)+P_rep*P_estim_C)*P_W;
P_repeat = P_rep*P_estim_W;

disp(['proba erreur total : ' num2str(Err_glob)]);
disp(['proba de répéter : ' num2str(P_repeat)]);

BpT_ini = log2(n)+(1-P_W)*log2(1-P_W)+P_W*log2(P_W/(n-1));
BpT_final = (log2(n)+(1-Err_glob)*log2(1-Err_glob)+Err_glob*log2(Err_glob/(n-1)))*(1-P_repeat);

disp(['BpT initial : ' num2str(BpT_ini)]);
disp(['BpT final : ' num2str(BpT_final)]);
