%programme_test_bb_dwt_chu2
% appelé directement et modifié pour signaux CHU

clear all
close all
warning off
%clc

tabvarpca = [2 3 4 5 6 7 8 10 15 20];
nb_subsets = 10;
wname_ini = 'db8';

rand('state',0);

%choix des canaux de mesure
canal = [1 2 3];
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

y_learn = activite;
signal_learn = matrice_signal;

Fe = 512;
chans = canal;
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
%disp(['Niveau de décomposition :' num2str(N_dec)]);
disp([filename]);
disp([wname_ini]);
disp(['canaux utilisées : ' num2str(canal)]);

%disp(['Nombre de variables pca : ' num2str(nbvarpca)]);

ind_des = [1:4]; %indice des marginales ou coeff dwt utilisés vide pour toutes
if isempty(ind_des)
    ind_des = 1:N_dec+1;
end


y_estim_ini = [];
y_real_ini = [];
y_estim_bb = [];
y_real_bb = [];
c_q_bb = [];
subset_length = floor(length(y_learn)/nb_subsets);
N_dec = N_dec-2;

%     N_dec = 3;
    dic_features_bb = set_wp_marg(signal_learn,wname_ini,N_dec);
    k = 0.2;
    wbar = waitbar(0,'computation...');
    %tic
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
        
        [features_acp_app,U,x_mean,x_std] = function_acp2(features_reshape_bb_app,y_bb_app);
        m=size(features_reshape_bb_test,1);
        xtest_centre = (features_reshape_bb_test - x_mean(ones(m,1),:))./x_std(ones(m,1),:);
        features_acp_test = xtest_centre*U;
        
        y_estim = [];
        for kpca = 1:length(tabvarpca)
            xapp = features_acp_app(:,1:tabvarpca(kpca));
            xtest = features_acp_test(:,1:tabvarpca(kpca));
            kernelparam.ktype = 1;
            kernelparam.kscale = 0.2;
            nb_class = max(y_learn);
            %initialisation de la séparatrice avec tous les élément de l'apprentissage
            sep_bb = svm_learning(xapp,y_bb_app,[],nb_class,kernelparam);
            y_estim = [y_estim  test_class_ovr(xtest, sep_bb, nb_class, xapp, y_bb_app, kernelparam)];
        end;
        y_estim_bb = [y_estim_bb; y_estim];
        
    end
    % keyboard
    %toc
    close(wbar)
   for kpca = 1:length(tabvarpca) disp(['Best Basis']);
        disp(['nbvar pca : ' num2str(tabvarpca(kpca))]);
        p_1_1_bb = 100*sum(y_estim_bb (:,kpca) == 1 & y_real_bb == 1)/sum(y_real_bb == 1);
        p_2_1_bb = 100*sum(y_estim_bb (:,kpca) == 2 & y_real_bb == 1)/sum(y_real_bb == 1);
        p_2_2_bb = 100*sum(y_estim_bb (:,kpca) == 2 & y_real_bb == 2)/sum(y_real_bb == 2);
        p_1_2_bb = 100*sum(y_estim_bb (:,kpca) == 1 & y_real_bb == 2)/sum(y_real_bb == 2);
        fprintf('p_c_estim_c_real\t | c_estim = 1\t | c_estim = 2\n')
        fprintf('c_real = 1\t | %f\t\t\t | %f\n', p_1_1_bb, p_2_1_bb)
        fprintf('c_real = 2\t | %f\t\t\t | %f\n', p_1_2_bb, p_2_2_bb)
        % disp(['mean of optim citerion : ' num2str(mean(critere_bb))]);
        c_q_bb = [c_q_bb;sum(y_estim_bb(:,kpca) ~= y_real_bb)/length(y_estim_bb)];
        disp(['pce_bb : ' num2str(c_q_bb(kpca))]);
        disp('-----------------------------------------------------------------');
   end;

    function_acp(features_reshape_bb_app,y_bb_app);
    title('Best Basis')
    drawnow
[tabvarpca' c_q_bb]
figure; plot (tabvarpca,c_q_bb)