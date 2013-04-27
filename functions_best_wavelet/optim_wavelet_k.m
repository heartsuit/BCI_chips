function [features_opt,h_opt, g_opt, c_opt] =...
    optim_wavelet_k(y_learn, features_dic, h_dic, g_dic, measure)
% optim_wavelet_k(x_learn(:,:,k), y_learn,...
%             features_opt, h_opt(k,:), g_opt(k,:),...
%             features_dic(:,:,k,:), h_dic(:,:,k), g_dic(:,:,k),...
%             k, ind_des, N_dec, measure);

%measure = 'pce' ou 'fisher'

c_opt=inf;
nb_waves = size(h_dic,1);
i_opt = 1;

for ii=1:nb_waves
    features = features_dic(:,:,ii);
    if strcmp(measure,'pce')
        %calcul du crière du pce
        [y_estim, y_real] = xval_procedure(features,y_learn,10);
        c_q = sum(y_estim ~= y_real)/length(y_estim);
    elseif strcmp(measure,'fisher')
        %clacul du critère de fisher
        c_q = fct_calc_fisher(features,y_learn);
    end
    if (c_q < c_opt)
        c_opt=c_q;
        i_opt = ii;
    end
end

% "optimal" filter
%choisi le set de features qui correspond au meilleur theta ( troisième
%coomposante du feature_dic dans ce programme (4ème dans le programme du
%dessus)
features_opt = features_dic(:,:,i_opt);
h_opt = h_dic(i_opt,:);
g_opt = g_dic(i_opt,:);

    