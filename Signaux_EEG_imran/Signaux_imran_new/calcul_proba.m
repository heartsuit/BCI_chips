clc;

detect.P_W_W = 0.63;
detect.P_C_C = 0.80;
detect.P_W_C = 1-detect.P_C_C;
detect.P_C_W = 1-detect.P_W_W;
detect.P_rep = 1;
svm.P_C = 0.78;
svm.P_W = 1-svm.P_C;
nbclass = 2;

P_W_W = detect.P_W_W;
P_C_C = detect.P_C_C;
P_W_C = detect.P_W_C;
P_C_W = detect.P_C_W;
P_rep = detect.P_rep;

P_C = svm.P_C;
P_W = svm.P_W;

P = 1/nbclass;

PWestimEestim_C = P*P_C*P_C_C+P*P_W*P_C_W;
PWestimEestim_W = P*P_C*P_W_C+P*P_W*P_W_W;

P_post_W_W = (P_W_W*P_W*P)/PWestimEestim_W;
P_post_C_C = (P_C_C*P_C*P)/PWestimEestim_C;
P_post_C_W = (P_W_C*P_C*P)/PWestimEestim_W;
P_post_W_C = (P_C_W*P_W*P)/PWestimEestim_C;

% disp(['P_post_wrong_estim_wrong :' num2str(P_post_W_W)]);
disp(['P_post_correct_estim_correct :' num2str(P_post_C_C)]);
disp(['P_post_correct_estim_wrong :' num2str(P_post_C_W)]);
% disp(['P_post_wrong_estim_correct :' num2str(P_post_W_C)]);

P_estim_C = P_C_W*P_W + P_C_C*P_C;
P_estim_W = P_W_W*P_W + P_W_C*P_C;

P_err = (P_C_W*P_W)/P_estim_C;

Err_glob = ((1-P_rep)*P_W+P_rep*P_estim_C*P_err)/((1-P_rep)+P_rep*P_estim_C);
P_repeat = P_rep*P_estim_W;

disp(['proba erreur total : ' num2str(Err_glob)]);
disp(['proba de répéter : ' num2str(P_repeat)]);

p = svm.P_C;
c = detect.P_C_C;
e = detect.P_W_W;
pt = 1-P_repeat;
Nc = nbclass;
p_prime = 1-Err_glob;

B_r_ini = log2(Nc)+p*log2(p)+(1-p)*log2((1-p)/(Nc-1));
B_r_new = pt*(log2(Nc)+p_prime*log2(p_prime)+(1-p_prime)*log2((1-p_prime)/(Nc-1)));

disp(['Bit rate sans Errp : ' num2str(B_r_ini)])
disp(['Bit rate avec Errp : ' num2str(B_r_new)])

