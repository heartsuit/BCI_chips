function programme_test_descripteur_xval_ep(file,path)
%Function pour charger les signaux
%Signaux sous la forme de structure
%Parametres de la structure class
%Doit contenir le champ signal de taille: nb_trials x nb_samples x nb_chans
if nargin<2
    [class_mvt, class_ep] = load_class;
else
    [class_mvt, class_ep] = load_class(file,path);
end

class = class_ep;
Fe = 1024;
t_start = 0.15;
t_end = 0.65;
N_start = round(t_start*Fe);
N_end = round(t_end*Fe)-1;
N_sample = 16;
N_selec = N_start:N_sample:N_end;
chans = 4:6;
if ~exist('chans','var')
    nb_chans = size(class(1).signal,3);
    chans = 1:nb_chans;
end
nb_class = size(class,2);

Wb = 0.5*2/1024;
Wh = 10*2/1024;
N_ord = 2;
[B,A] = butter(N_ord,Wh);

for i=1:nb_class
    for j=1:size(class(i).signal,1)
        for k=1:size(class(i).signal,3)
            class(i).signal(j,:,k) = filter(B,A,class(i).signal(j,:,k));
        end
    end
end


for i=1:nb_class
    class(i).signal = class(i).signal(:,N_selec,:);
end
%--------------------------------------------------------------------------
%Cr�ation de la matrice contenant les �l�ments de chaque classe et le
%vecteur des labels associ�
field  = 'signal';
signal_learn = [];  %Matrice concat�nant tous les essais de chaque classe
%Dim : nb_trials_all x nb_samples x nb_chans
y_learn = [];
for i=1:nb_class
    signal_learn = cat(1,signal_learn,class(i).(field));
    y_learn = cat(1, y_learn, i*ones(size(class(i).(field),1),1));
end
signal_learn = signal_learn(:,:,chans);

%On m�lange les essais pour ne pas avoir tous les signaux d'une m�me classe
%ensemble
ind_perm = randperm(length(y_learn)); %Vecteur de permutation
signal_learn = signal_learn(ind_perm,:,:);
y_learn = y_learn(ind_perm);

%Param�tres de l'ensemble d'apprentissage
[nb_trials,nb_samples,nb_chans] = size(signal_learn);

N_dec_max = floor(log2(nb_samples));
N_dec = N_dec_max;
disp(['Niveau de d�composition :' num2str(N_dec)]);
ind_des = [1:4]; %indice des marginales ou coeff dwt utilis�s vide pour toutes
if isempty(ind_des)
    ind_des = 1:N_dec+1;
end

% for i=1:size(signal_learn,1)
%     for j=1:nb_chans
%         signal_learn(i,:,j) = filter(B,A,signal_learn(i,:,j));
%     end
% end

features_ini = signal_learn;


% figure
% plot(1:nb_samples,mean(signal_learn(y_learn==1,:,1),1),'b',1:nb_samples,mean(signal_learn(y_learn==2,:,1),1),'r')
nb_subsets = 120;
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
    %initialisation de la s�paratrice avec tous les �l�ment de l'apprentissage
    sep_ini = svm_learning(features_reshape_ini_app,y_ini_app,[],nb_class,kernelparam);
    y_estim_ini = [y_estim_ini; test_class_ovr(features_reshape_ini_test, sep_ini, nb_class, features_reshape_ini_app, y_ini_app, kernelparam)];
end
toc
close(wbar)

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