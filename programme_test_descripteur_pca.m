function programme_test_descripteur_pca(file,path)
%Function pour charger les signaux
%Signaux sous la forme de structure
%Parametres de la structure class
%Doit contenir le champ signal de taille: nb_trials x nb_samples x nb_chans
if nargin<2
    class = load_class;
else
    class = load_class(file,path);
end

Fe = 1024;
t_start = 3.5;
t_end = 5.5;
N_start = round(t_start*Fe);
N_end = round(t_end*Fe)-1;
N_selec = N_start:N_end;
chans = 4:6;
if ~exist('chans','var')
    nb_chans = size(class(1).signal,3);
    chans = 1:nb_chans;
end
nb_class = size(class,2);
for i=1:nb_class
    class(i).signal = class(i).signal(:,N_selec,:);
end
%--------------------------------------------------------------------------
%Création de la matrice contenant les éléments de chaque classe et le
%vecteur des labels associé
field  = 'signal';
signal_learn = [];  %Matrice concaténant tous les essais de chaque classe
%Dim : nb_trials_all x nb_samples x nb_chans
y_learn = [];
for i=1:nb_class
    signal_learn = cat(1,signal_learn,class(i).(field));
    y_learn = cat(1, y_learn, i*ones(size(class(i).(field),1),1));
end
signal_learn = signal_learn(:,:,chans);

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

wname_ini = 'db2';

optim_fisher = 1;
optim_pce = 0;
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
for i=1:nb_trials
    for k=1:nb_chans
        features_ini(i,:,k) = calc_features(signal_learn(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
    end
end



features_ini_app = features_ini;
y_ini_app = y_learn;
features_reshape_ini_app = reshape(features_ini_app,length(y_ini_app),[]);
function_acp(features_reshape_ini_app,y_ini_app);



features_ini_app = features_ini;
y_pce_app = y_learn;
[features_pce_app] = optim_wavelet_multi_cv(y_pce_app, features_ini_app, h_ini, g_ini, features_dic, h_dic, g_dic, 'pce');
features_reshape_pce_app = reshape(features_pce_app,length(y_pce_app),[]);
function_acp(features_reshape_pce_app,y_pce_app)


features_ini_app = features_ini;
y_fisher_app = y_learn;
[features_fisher_app] = optim_wavelet_multi_cv(y_fisher_app, features_ini_app, h_ini, g_ini, features_dic, h_dic, g_dic, 'fisher');
features_reshape_fisher_app = reshape(features_fisher_app,length(y_fisher_app),[]);
function_acp(features_reshape_fisher_app,y_fisher_app)


N_dec = N_dec-2;
% N_dec = 3;
dic_features_bb = set_wp_marg(signal_learn,wname_ini,N_dec);
K = 0.2;

features_bb_app = dic_features_bb;
y_bb_app = y_learn;
basis_bb = best_basis_search(features_bb_app,y_bb_app,N_dec,K);
features_reshape_bb_app = reshape(features_bb_app(:,boolean(basis_bb),:),length(y_bb_app),[]);
function_acp(features_reshape_bb_app,y_bb_app)
