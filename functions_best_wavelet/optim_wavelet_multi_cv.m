function [features_opt, h_opt, g_opt, c_opt] =optim_wavelet_multi_cv(...
    y_learn,features_ini, h_ini, g_ini,features_dic, h_dic, g_dic, measure)
          

nb_trials = size(features_ini,1);
nb_chans = size(features_ini,3);
nb_sets = 10;
%h_opt : nb_chans x nb_coeff_filter
%Initialisation du filtre
h_opt = h_ini;
g_opt = g_ini;
features_opt = features_ini;
features_reshape = reshape(features_opt,nb_trials,[]);
if strcmp(measure,'fisher')
    %Critères fisher
    c_ini = fct_calc_fisher(features_reshape,y_learn);
elseif strcmp(measure,'pce')
    %Critère pce
    [y_estim, y_real] = xval_procedure(features_reshape,y_learn,nb_sets);
    c_ini = sum(y_estim ~= y_real)/length(y_estim);
end

K = 0;
c_opt = c_ini;
%Tant que le critère n'est pas amélioré K fois on optimise selon une voie
while K < nb_chans
    for k=1:nb_chans
        features_new = features_opt;
        %calcul le theta pselon le critère de fisher en renvoyant h,g et
        [features_new_k,h_new,g_new] = ...
            optim_wavelet_k(y_learn, features_dic(:,:,k,:),...
            h_dic(:,:,k), g_dic(:,:,k), measure);
        %         h_opt_temp(k,:) = h_new;
        %         g_opt_temp(k,:) = g_new;
        features_new(:,:,k) = features_new_k;
        features_reshape = reshape(features_new,nb_trials,[]);
        if strcmp(measure,'fisher')
            %         Critères fisher
            c_new = fct_calc_fisher(features_reshape,y_learn);
        elseif strcmp(measure,'pce')
            %         Critère pce
            [y_estim, y_real] = xval_procedure(features_reshape,...
                y_learn,nb_sets);
            c_new = sum(y_estim ~= y_real)/length(y_estim);
        end
        if c_new < c_opt
            features_opt = features_new;
            c_opt = c_new;
            h_opt(k,:) = h_new;
            g_opt(k,:) = g_new;
            K=1;
        else
            K=K+1;
            if K == nb_chans
                break;
            end
        end
    end
end

