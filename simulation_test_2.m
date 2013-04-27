close all
clear all
clc

addpath(genpath(cd));


load toy2_3v_2c.mat
% 
% [x_app, y_app, x_test, y_test] = fct_signal_simulation(N, h1, h2, nApp, nTest, proba, RSB);
kernelparam.ktype = 1;
kernelparam.kscale = 0.2;

x_app = permute(reshape(permute(x,[1 3 2 4]),size(x,1),size(x,3),size(x,2)*size(x,4)),[3 1 2]);
y_app = [ones(size(x,2),1); 2*ones(size(x,2),1)];
x_test = permute(reshape(permute(y,[1 3 2 4]),size(y,1),size(y,3),size(y,2)*size(y,4)),[3 1 2]);
y_test = [ones(size(y,2),1); 2*ones(size(y,2),1)];

n_perm_app = randperm(length(y_app));
x_app = x_app(n_perm_app,:,:);
y_app = y_app(n_perm_app);

n_perm_test = randperm(length(y_test));
x_test = x_test(n_perm_test,:,:);
y_test = y_test(n_perm_test);

%Paramètres de l'ensemble d'apprentissage
[nb_trials_app,nb_samples,nb_chans] = size(x_app);
nb_trials_test = size(x_test,1);
nb_class = max(y_app);


N_dec_max = floor(log2(nb_samples));
N_dec = N_dec_max-2;
disp(['niveau de décomposition :' num2str(N_dec)]);
ind_des = [2:N_dec+1]; %indice des marginales ou coeff dwt utilisés vide pour toutes
if isempty(ind_des)
    ind_des = 1:N_dec+1;
end


wname_ini = 'db2';

optim_fisher = 1;
optim_pce = 0;
if optim_fisher || optim_pce
    %on calcul les descripteurs de tous les individus pour
    %tous les paramètres d'ondelettes possibles
    [features_dic,h_dic,g_dic] = create_dic(x_app,N_dec,ind_des);
end

[h_ini, g_ini] = wfilters(wname_ini);
h_ini = ones(nb_chans,1)*h_ini;
g_ini = ones(nb_chans,1)*g_ini;
% features_ini_app = zeros(nb_trials_app,length(ind_des),nb_chans);
% features_ini = zeros(nb_trials,nb_samples,nb_chans);
for i=1:nb_trials_app
    for k=1:nb_chans
        features_ini_app(i,:,k) = calc_features(x_app(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
    end
end

% features_ini_test = zeros(nb_trials_test,length(ind_des),nb_chans);
for i=1:nb_trials_test
    for k=1:nb_chans
        features_ini_test(i,:,k) = calc_features(x_test(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
    end
end

features_reshape_ini_app = reshape(features_ini_app,nb_trials_app,[]);
features_reshape_ini_test = reshape(features_ini_test,nb_trials_test,[]);

sep_ini = svm_learning(features_reshape_ini_app,y_app,[],nb_class,kernelparam);
y_estim_ini = test_class_ovr(features_reshape_ini_test, sep_ini, nb_class, features_reshape_ini_app, y_app, kernelparam);

pce_db2 = sum(y_estim_ini~=y_test)/length(y_test);

% wname_ini = 'db2';
% [h_ini, g_ini] = wfilters(wname_ini);
% h_ini = ones(nb_chans,1)*h_ini;
% g_ini = ones(nb_chans,1)*g_ini;

y_pce_app = y_app;
[features_pce_app,h_pce,g_pce,critere_pce] = optim_wavelet_multi_cv(y_pce_app, features_ini_app, h_ini, g_ini, features_dic, h_dic, g_dic, 'pce');
for i=1:nb_trials_test
    for k=1:nb_chans
        features_pce_test(i,:,k) = calc_features(x_test(i,:,k),N_dec,ind_des,h_pce(k,:),g_pce(k,:));
    end
end
features_reshape_pce_app = reshape(features_pce_app,nb_trials_app,[]);
features_reshape_pce_test = reshape(features_pce_test,nb_trials_test,[]);
sep_pce = svm_learning(features_reshape_pce_app,y_pce_app,[],nb_class,kernelparam);
y_estim_pce = test_class_ovr(features_reshape_pce_test, sep_pce, nb_class, features_reshape_pce_app, y_pce_app, kernelparam);

pce_pce = sum(y_estim_pce~=y_test)/length(y_test);

y_fisher_app = y_app;
[features_fisher_app,h_fisher,g_fisher,critere_fisher] = optim_wavelet_multi_cv(y_fisher_app, features_ini_app, h_ini, g_ini, features_dic, h_dic, g_dic, 'fisher');
for i=1:nb_trials_test
    for k=1:nb_chans
        features_fisher_test(i,:,k) = calc_features(x_test(i,:,k),N_dec,ind_des,h_fisher(k,:),g_fisher(k,:));
    end
end

features_reshape_fisher_app = reshape(features_fisher_app,nb_trials_app,[]);
features_reshape_fisher_test = reshape(features_fisher_test,nb_trials_test,[]);
sep_fisher = svm_learning(features_reshape_fisher_app,y_fisher_app,[],nb_class,kernelparam);
y_estim_fisher = test_class_ovr(features_reshape_fisher_test, sep_fisher, nb_class, features_reshape_fisher_app, y_fisher_app, kernelparam);

pce_fisher = sum(y_estim_fisher~=y_test)/length(y_test);



%     N_dec = N_dec-2;
    N_dec = 5;
features_bb_app = set_wp_marg(x_app,wname_ini,N_dec);
features_bb_test = set_wp_marg(x_test,wname_ini,N_dec);
K = 0.2;

y_bb_app = y_app;
basis_bb = best_basis_search(features_bb_app,y_bb_app,N_dec,K);
%     disp(sum(basis_bb == 1));
features_reshape_bb_app = reshape(features_bb_app(:,boolean(basis_bb),:),nb_trials_app,[]);
features_reshape_bb_test = reshape(features_bb_test(:,boolean(basis_bb),:),nb_trials_test,[]);

%initialisation de la séparatrice avec tous les élément de l'apprentissage
sep_bb = svm_learning(features_reshape_bb_app,y_bb_app,[],nb_class,kernelparam);
y_estim_bb = test_class_ovr(features_reshape_bb_test, sep_bb, nb_class, features_reshape_bb_app, y_bb_app, kernelparam);

pce_bb = sum(y_estim_bb~=y_test)/length(y_test);



% 
% 
% disp(['db2        pce         fisher       best basis'])
% disp([num2str(mean(pce_db2)) '   ' num2str(mean(pce_pce)) '   ' num2str(mean(pce_fisher)) '   ' num2str(mean(pce_bb))])


disp(['db2         fisher       best basis'])
disp([num2str(mean(pce_db2)) '   ' num2str(mean(pce_fisher)) '   ' num2str(mean(pce_bb))])


rmpath(genpath(cd));