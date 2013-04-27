%programme_test_descripteur_xval
% appelé directement et modifié pour signaux CHU

clear all
close all
clc
dwt_without_optim = 1;
optim_fisher = 0;
optim_pce = 0;
optim_basis = 0;

rand('state',0);
% 
% disp(sub_name);
% figure('Name',sub_name);

param = 1;
%choix des canaux de mesure
canal = [1 2];
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

[matrice_signal, activite, Fe] = loadgdf2matrice (filename, canal, leng, overl);
nb_subsets = 10; 
%nb_subsets = length(activite); si leave one out

y_learn = activite;
signal_learn = matrice_signal;

Nfft=5000;
Fe = 512;
tabf_uniforme=[0 5 10 15 20 25 30 40 50 60];
tabf_dyadique=[0 2.5 5 10 20 40 80 160];
%chans = 1;
if ~exist('chans','var')
    nb_chans = size(matrice_signal,3);
    chans = 1:nb_chans;
end
nb_class = 2;
%--------------------------------------------------------------------------

%On mélange les essais pour ne pas avoir tous les signaux d'une même classe
%
ind_perm = randperm(length(y_learn)); %Vecteur de permutation
signal_learn = signal_learn(ind_perm,:,:);
y_learn = y_learn(ind_perm);

%Paramètres de l'ensemble d'apprentissage
[nb_trials,nb_samples,nb_chans] = size(signal_learn);

for i=1:nb_trials
    for k=1:nb_chans
        features_ini(i,:,k) = calcfeatures_fourier(signal_learn(i,:,k),Fe,Nfft,tabf_uniforme,tabf_dyadique,param);%origine
        
    end
end

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
%         initialisation de la séparatrice avec tous les élément de
%         l'apprentissage
        sep_ini = svm_learning(features_reshape_ini_app,y_ini_app,[],nb_class,kernelparam);
        y_estim_ini = [y_estim_ini; test_class_ovr(features_reshape_ini_test, sep_ini, nb_class, features_reshape_ini_app, y_ini_app, kernelparam)];
    end
    toc
    close(wbar)
 
    
    function_acp(features_reshape_ini_app,y_ini_app);
    title('Puissances par bandes')
    drawnow
    P_1_1_ini = 100*sum(y_estim_ini == 1 & y_real_ini == 1)/sum(y_real_ini == 1);
    P_2_1_ini = 100*sum(y_estim_ini == 2 & y_real_ini == 1)/sum(y_real_ini == 1);
    P_2_2_ini = 100*sum(y_estim_ini == 2 & y_real_ini == 2)/sum(y_real_ini == 2);
    P_1_2_ini = 100*sum(y_estim_ini == 1 & y_real_ini == 2)/sum(y_real_ini == 2);
    
    fprintf('P_C_estim_C_real\t | C_estim = 1\t | C_estim = 2\n')
    fprintf('C_real = 1\t | %f\t\t\t | %f\n', P_1_1_ini, P_2_1_ini)
    fprintf('C_real = 2\t | %f\t\t\t | %f\n', P_1_2_ini, P_2_2_ini)
    c_q_ini = sum(y_estim_ini ~= y_real_ini)/length(y_estim_ini);
    if param == 0
        disp(['PCE_puissance_bande uniforme : ' num2str(c_q_ini)]);
    else
        disp(['PCE_puissance_bande dyadique : ' num2str(c_q_ini)]);
    end;

    disp('-----------------------------------------------------------------');
end
