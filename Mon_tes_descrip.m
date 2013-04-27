%programme_test_descripteur_xval
% appelé directement et modifié pour signaux CHU

clear all
close all
%clc
dwt_without_optim = 1;
optim_fisher = 0;
optim_pce = 0;
optim_basis = 0;

nb_subsets = 10;

rand('state',0);
% 
% disp(sub_name);
% figure('Name',sub_name);

%choix des canaux de mesure
canal = [1];
%ratio de chevauchement de la fenêtre coulisaante
overl = 0.5;%7/8;
%longueur (s) de la fenêtre coulissante
leng = 1;
%[filename, pathname, filterindex] = uigetfile('*.gdf', 'Pick a GDF-file');
filename = 'record-[2011.07.08-10.03.51].gdf'; %main gauche 52
%filename = 'record-[2011.07.07-11.21.26].gdf'; % pied gauche 26
%filename = 'record-[2011.07.06-16.36.28].gdf';% main gauche 14
filename = 'record-[2011.07.07-15.25.52].gdf'; % reel pied
%file = '[2011.07.07-15.25.52]'; % reel pied

[matrice_signal, activite, fe] = loadgdf2matrice (filename, canal, leng, overl);
size(matrice_signal)
size(activite)

% figure;hold on
% for i=1:size(matrice_signal,1)
%     plot(matrice_signal(i,:,1))
% end;
% 
% signal = [];
% for i=1:size(matrice_signal,1)
%     signal = [signal matrice_signal(i,:)];
% end;
% figure; plot(signal)
% pause

y_learn = activite;
signal_learn = matrice_signal;

wname_ini = 'db2';
Fe = 512;
chans = 1;
if ~exist('chans','var')
    nb_chans = size(matrice_signal,3);
    chans = 1:nb_chans;
end
nb_class = 2;
%--------------------------------------------------------------------------

%On mélange les essais pour ne pas avoir tous les signaux d'une même classe
%ensemble
ind_perm = randperm(length(y_learn)); %Vecteur de permutation
signal_learn = signal_learn(ind_perm,:,:);
y_learn = y_learn(ind_perm);

%Paramètres de l'ensemble d'apprentissage
[nb_trials,nb_samples,nb_chans] = size(signal_learn);

N_dec_max = floor(log2(nb_samples));
N_dec = N_dec_max;
disp(['Niveau de décomposition :' num2str(N_dec)]);
ind_des = [1:4]; %indice des marginales ou coeff dwt utilisés vide pour toutes
if isempty(ind_des)
    ind_des = 1:N_dec+1;
end



tic
if optim_fisher || optim_pce
    %on calcul les descripteurs de tous les individus pour
    %tous les paramètres d'ondelettes possibles
    [features_dic,h_dic,g_dic] = create_dic(signal_learn,N_dec,ind_des);
end
toc

[h_ini, g_ini] = wfilters(wname_ini);
h_ini = ones(nb_chans,1)*h_ini;
g_ini = ones(nb_chans,1)*g_ini;
% features_ini = zeros(nb_trials,length(ind_des),nb_chans);
% features_ini = zeros(nb_trials,nb_samples,nb_chans);
for i=1:nb_trials
    for k=1:nb_chans
        features_ini(i,:,k) = calc_features(signal_learn(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
    end
end

% measure = 'fisher';
% [features_opt,h_opt,g_opt,c_opt] = optim_wavelet_multi_cv(y_learn, features_ini, h_ini, g_ini, features_dic, h_dic, g_dic, measure);


y_estim_ini = [];
y_real_ini = [];
y_estim_pce = [];
y_real_pce = [];
y_estim_fisher = [];
y_real_fisher = [];
y_estim_bb = [];
y_real_bb = [];
subset_length = floor(length(y_learn)/nb_subsets);
critere_ini = zeros(nb_subsets,1);
critere_pce = zeros(nb_subsets,1);
critere_fisher = zeros(nb_subsets,1);

if dwt_without_optim
    wbar = waitbar(0,'Computation...');
    tic
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar)
        %     signal_test = signal_learn(i,:,k);
        %     y_test = y_learn(i);
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        %     y_app = y_learn(ind_app);
        features_ini_app = features_ini(ind_app,:,:);
        features_ini_test = features_ini(ind_test,:,:);
        y_ini_app = y_learn(ind_app);
        y_real_ini = [y_real_ini; y_learn(ind_test)];
        
        features_reshape_ini_app = reshape(features_ini_app,length(ind_app),[]);
        features_reshape_ini_test = reshape(features_ini_test,length(ind_test),[]);
        
        
        critere_ini(ii) = fct_calc_fisher(features_reshape_ini_app,y_ini_app);
        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);
        %initialisation de la séparatrice avec tous les élément de
        %l'apprentissage
        sep_ini = svm_learning(features_reshape_ini_app,y_ini_app,[],nb_class,kernelparam);
        y_estim_ini = [y_estim_ini; test_class_ovr(features_reshape_ini_test, sep_ini, nb_class, features_reshape_ini_app, y_ini_app, kernelparam)];
    end
    toc
    close(wbar)
 
    
    subplot(2,2,1)
    function_acp(features_reshape_ini_app,y_ini_app);
    title('Dwt sans optim')
    drawnow
    P_1_1_ini = 100*sum(y_estim_ini == 1 & y_real_ini == 1)/sum(y_real_ini == 1);
    P_2_1_ini = 100*sum(y_estim_ini == 2 & y_real_ini == 1)/sum(y_real_ini == 1);
    P_2_2_ini = 100*sum(y_estim_ini == 2 & y_real_ini == 2)/sum(y_real_ini == 2);
    P_1_2_ini = 100*sum(y_estim_ini == 1 & y_real_ini == 2)/sum(y_real_ini == 2);
    
    fprintf('P_C_estim_C_real\t | C_estim = 1\t | C_estim = 2\n')
    fprintf('C_real = 1\t | %f\t\t\t | %f\n', P_1_1_ini, P_2_1_ini)
    fprintf('C_real = 2\t | %f\t\t\t | %f\n', P_1_2_ini, P_2_2_ini)
    disp(['mean of optim citerion : ' num2str(mean(critere_ini))]);
    c_q_ini = sum(y_estim_ini ~= y_real_ini)/length(y_estim_ini);
    disp(['PCE_db2 : ' num2str(c_q_ini)]);
    disp('-----------------------------------------------------------------');
end

if optim_pce
    
    [h_ini, g_ini] = wfilters(wname_ini);
    h_ini = ones(nb_chans,1)*h_ini;
    g_ini = ones(nb_chans,1)*g_ini;
    
    %--------------------------------------------------------------------------
    wbar = waitbar(0,'computation...');
    tic
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar)
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        features_ini_app = features_ini(ind_app,:,:);
        y_pce_app = y_learn(ind_app);
        [features_pce_app,h_pce,g_pce,critere_pce(ii)] = optim_wavelet_multi_cv(y_pce_app, features_ini_app, h_ini, g_ini, features_dic(ind_app,:,:,:), h_dic, g_dic, 'pce');
        features_reshape_pce_app = reshape(features_pce_app,length(ind_app),[]);
        %     features_pce_test = zeros(length(ind_test),length(ind_des),nb_chans);
        % features_pce_test = zeros(length(ind_test),nb_samples,nb_chans);
        for i=1:length(ind_test)
            for k=1:nb_chans
                features_pce_test(i,:,k) = calc_features(signal_learn(ind_test(i),:,k),N_dec,ind_des,h_pce(k,:),g_pce(k,:));
            end
        end
        features_reshape_pce_test = reshape(features_pce_test,length(ind_test),[]);
        
        y_real_pce = [y_real_pce; y_learn(ind_test)];
        
        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);
        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_pce = svm_learning(features_reshape_pce_app,y_pce_app,[],nb_class,kernelparam);
        y_estim_pce = [y_estim_pce; test_class_ovr(features_reshape_pce_test, sep_pce, nb_class, features_reshape_pce_app, y_pce_app, kernelparam)];
        
    end
    toc
    close(wbar)
    
    p_1_1_pce = 100*sum(y_estim_pce == 1 & y_real_pce == 1)/sum(y_real_pce == 1);
    p_2_1_pce = 100*sum(y_estim_pce == 2 & y_real_pce == 1)/sum(y_real_pce == 1);
    p_2_2_pce = 100*sum(y_estim_pce == 2 & y_real_pce == 2)/sum(y_real_pce == 2);
    p_1_2_pce = 100*sum(y_estim_pce == 1 & y_real_pce == 2)/sum(y_real_pce == 2);
    
    subplot(2,2,2)
    function_acp(features_reshape_pce_app,y_pce_app);
    title('Dwt optim Fisher')
    drawnow
    fprintf('p_c_estim_c_real\t | c_estim = 1\t | c_estim = 2\n')
    fprintf('c_real = 1\t | %f\t\t\t | %f\n', p_1_1_pce, p_2_1_pce)
    fprintf('c_real = 2\t | %f\t\t\t | %f\n', p_1_2_pce, p_2_2_pce)
    disp(['mean of optim citerion : ' num2str(mean(critere_pce))]);
    c_q_pce = sum(y_estim_pce ~= y_real_pce)/length(y_estim_pce);
    disp(['pce_pce : ' num2str(c_q_pce)]);
    disp('-----------------------------------------------------------------');
end

%--------------------------------------------------------------------------
if optim_fisher
    wbar = waitbar(0,'Computation...');
    tic
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar);
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        features_ini_app = features_ini(ind_app,:,:);
        y_fisher_app = y_learn(ind_app);
        [features_fisher_app,h_fisher,g_fisher,critere_fisher(ii)] = optim_wavelet_multi_cv(y_fisher_app, features_ini_app, h_ini, g_ini, features_dic(ind_app,:,:,:), h_dic, g_dic, 'fisher');
        features_reshape_fisher_app = reshape(features_fisher_app,length(ind_app),[]);
        %     features_fisher_test = zeros(length(ind_test),length(ind_des),nb_chans);
        % features_fisher_test = zeros(length(ind_test),nb_samples,nb_chans);
        for i=1:length(ind_test)
            for k=1:nb_chans
                features_fisher_test(i,:,k) = calc_features(signal_learn(ind_test(i),:,k),N_dec,ind_des,h_fisher(k,:),g_fisher(k,:));
            end
        end
        features_reshape_fisher_test = reshape(features_fisher_test,length(ind_test),[]);
        
        y_real_fisher = [y_real_fisher; y_learn(ind_test)];
        
        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);
        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_fisher = svm_learning(features_reshape_fisher_app,y_fisher_app,[],nb_class,kernelparam);
        y_estim_fisher = [y_estim_fisher; test_class_ovr(features_reshape_fisher_test, sep_fisher, nb_class, features_reshape_fisher_app, y_fisher_app, kernelparam)];
        
        
    end
    toc
    close(wbar)
    
    subplot(2,2,3)
    function_acp(features_reshape_fisher_app,y_fisher_app);
    title('Dwt optim Fisher')
    drawnow
    P_1_1_fisher = 100*sum(y_estim_fisher == 1 & y_real_fisher == 1)/sum(y_real_fisher == 1);
    P_2_1_fisher = 100*sum(y_estim_fisher == 2 & y_real_fisher == 1)/sum(y_real_fisher == 1);
    P_2_2_fisher = 100*sum(y_estim_fisher == 2 & y_real_fisher == 2)/sum(y_real_fisher == 2);
    P_1_2_fisher = 100*sum(y_estim_fisher == 1 & y_real_fisher == 2)/sum(y_real_fisher == 2);
    
    fprintf('P_C_estim_C_real\t | C_estim = 1\t | C_estim = 2\n')
    fprintf('C_real = 1\t | %f\t\t\t | %f\n', P_1_1_fisher, P_2_1_fisher)
    fprintf('C_real = 2\t | %f\t\t\t | %f\n', P_1_2_fisher, P_2_2_fisher)
    disp(['mean of optim citerion : ' num2str(mean(critere_fisher))]);
    c_q_fisher = sum(y_estim_fisher ~= y_real_fisher)/length(y_estim_fisher);
    disp(['PCE_fisher : ' num2str(c_q_fisher)]);
    disp('-----------------------------------------------------------------');
end


if optim_basis
    N_dec = N_dec-2;
%     N_dec = 3;
    dic_features_bb = set_wp_marg(signal_learn,wname_ini,N_dec);
    k = 0.2;
    wbar = waitbar(0,'computation...');
    tic
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar);
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        features_bb_app = dic_features_bb(ind_app,:,:);
        features_bb_test = dic_features_bb(ind_test,:,:);
        y_bb_app = y_learn(ind_app);
        basis_bb = best_basis_search(features_bb_app,y_bb_app,N_dec,k);
        %     disp(sum(basis_bb == 1));
        features_reshape_bb_app = reshape(features_bb_app(:,boolean(basis_bb),:),length(ind_app),[]);
        features_reshape_bb_test = reshape(features_bb_test(:,boolean(basis_bb),:),length(ind_test),[]);
        y_real_bb = [y_real_bb; y_learn(ind_test)];
        
        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);
        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_bb = svm_learning(features_reshape_bb_app,y_bb_app,[],nb_class,kernelparam);
        y_estim_bb = [y_estim_bb; test_class_ovr(features_reshape_bb_test, sep_bb, nb_class, features_reshape_bb_app, y_bb_app, kernelparam)];
        
        
    end
    % keyboard
    toc
    close(wbar)
    
    subplot(2,2,4)
    function_acp(features_reshape_bb_app,y_bb_app);
    title('Dwt optim Fisher')
    drawnow
    p_1_1_bb = 100*sum(y_estim_bb == 1 & y_real_bb == 1)/sum(y_real_bb == 1);
    p_2_1_bb = 100*sum(y_estim_bb == 2 & y_real_bb == 1)/sum(y_real_bb == 1);
    p_2_2_bb = 100*sum(y_estim_bb == 2 & y_real_bb == 2)/sum(y_real_bb == 2);
    p_1_2_bb = 100*sum(y_estim_bb == 1 & y_real_bb == 2)/sum(y_real_bb == 2);
    
    fprintf('p_c_estim_c_real\t | c_estim = 1\t | c_estim = 2\n')
    fprintf('c_real = 1\t | %f\t\t\t | %f\n', p_1_1_bb, p_2_1_bb)
    fprintf('c_real = 2\t | %f\t\t\t | %f\n', p_1_2_bb, p_2_2_bb)
    % disp(['mean of optim citerion : ' num2str(mean(critere_bb))]);
    c_q_bb = sum(y_estim_bb ~= y_real_bb)/length(y_estim_bb);
    disp(['pce_bb : ' num2str(c_q_bb)]);
    disp('-----------------------------------------------------------------');
end