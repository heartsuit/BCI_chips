function programme_test_descripteur_xval(subject,optim,plot_tree,kept_axes,...
    nb_subsets,simu_para,restric,filtre,wname_ini,plot_raw,plot_real_sig)
%Function pour charger les signaux
%Signaux sous la forme de structure
%Parametres de la structure class
%Doit contenir le champ signal de taille: nb_trials x nb_samples x nb_chans

%lol

dwt_without_optim=optim.dwt_without_optim;
optim_fisher=optim.optim_fisher;
optim_pce=optim.optim_pce;
optim_basis=optim.optim_basis;
optim_basis_ACP=optim.optim_basis_ACP;
optim_algo=optim.algo;

file=subject.filename;
path=subject.pathname;
sub_name=subject.sub_name;

sim_sig=simu_para.on;

% load('../fastcheby.mat');

subset_=1;
Fe = 1024;
t_start = 3.5;
t_end = 5.5;
N_start = round(t_start*Fe);
N_end = round(t_end*Fe)-1;
N_selec = N_start:N_end;

if sim_sig==0
    if ~exist('signal_learn','var')
        if nargin<2
          class = load_class;

        else
          class = load_class(file,path);
        end

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
        plot_tr.on=plot_tree;
        plot_tr.fe=Fe;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if filtre
        % standard wp and ws values: 35 , 45
        % plotèraw values 10 , 20
        [n_cheb,W_cheb]=cheb2ord(2*35/1024,2*45/1024,3,60);
        
        [b_filt,a_filt]=cheby2(n_cheb,60,W_cheb);
        
        coef=floor(1/W_cheb);
        
%         signal_learn_2=filter(b_filt,a_filt,signal_learn(:,:,3)');
        
        for i=1:length(chans)
            signal_learn_2=filter(b_filt,a_filt,signal_learn(:,:,i)');
%             signal_learn_22(:,:,i)=signal_learn_2';
            inter=downsample(signal_learn_2,coef);

%             inter=downsample(signal_learn(:,:,i)',coef);
            signal_learn_inter(:,:,i)=inter';
        end
        Fe=Fe/coef;
        signal_learn=signal_learn_inter;

%%%%%%%%%%%%%%%%%%% RAW SIGNAL VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if plot_raw
            
            mean1=mean(signal_learn(y_learn==1,:,1),1);
            std1=std(signal_learn(y_learn==1,:,1),1);
            mean2=mean(signal_learn(y_learn==2,:,1),1);
            std2=std(signal_learn(y_learn==2,:,1),1);      
            for i=1:3
                mean1(i)=mean1(4);
                std1(i)=std1(4);
                mean2(i)=mean2(4);
                std2(i)=std2(4);    
            end
            
            
            figure;
            hold on;
            title('raw data representation');
            plot([-0.5:1/Fe:1.5],mean1,'-b','LineWidth',2);
            plot([-0.5:1/Fe:1.5],mean1+std1,':b');
            plot([-0.5:1/Fe:1.5],mean1-std1,':b');
            plot([-0.5:1/Fe:1.5],mean2,'-r','LineWidth',2);
            plot([-0.5:1/Fe:1.5],mean2+std2,':r');
            plot([-0.5:1/Fe:1.5],mean2-std2,':r');
            plot([0,0],[min([mean1-std1,mean2-std2]),max([mean1+std1,mean2+std2])],'-k');
            xlabel('Time (s)');
            ylabel('U (mV)');
           
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %représentation fréquentielle de la classe 1 addition des FFT pour
    %chaque élément de la classe. vecteur mis dans plot_tr.f1, vecteur des
    %fréquences dans plot_tr.f_scal
    %PS: on ne se concentre que sur  la voie 1
    

    
    signal_1=signal_learn(y_learn==1,:,1);
    [ plot_tr.f1,plot_tr.f_scal]=fftsum(signal_1,Fe);
    

    %représentation classe 2 dans plot_tr.f2

    signal_2=signal_learn(y_learn==2,:,1);
        plot_tr.f2=fftsum(signal_2,Fe);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SIGNAL VISUALISATION  %%%%%%%%%%%%%%%%%%%%
    if plot_real_sig
        tri=10;
        figure('Name',['temporal representation of the signal number ' num2str(tri) ' with its mean']);

    %         [a,b,c]=analyseFFT(signal_learn(tri,:,1)',Fe);
        [a,b,c]=analyseFFT((signal_learn(tri,:,1)-mean(signal_learn(tri,:,1)))',Fe);

        plot(signal_learn(tri,:,1),'color','b');
        hold on;
        disp(mean(signal_learn(tri,:,1)));
        plot(mean(signal_learn(tri,:,1))*ones(length(signal_learn(tri,:,1)),1),'color','r');
        figure('Name', 'fourrier representation of this signal');
        plot(a(1:floor(c/2)),abs(b(1:floor(c/2))));
        keyboard;
    end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




else

    [y_learn,signal_learn,f_simu,f_class,StN]= simu_creation...
        (simu_para.nb_trials,simu_para.N_dec_max,simu_para.time,simu_para.nb_class,...
        Fe,simu_para.sigma,simu_para.posi,simu_para.perfect);
    
    N_dec_max=simu_para.N_dec_max;
    nb_class=simu_para.nb_class;
    subset_=simu_para.subset;
    
    disp(['Signal-to-Noise ratio mean of simulated signals = ' num2str(mean(StN))]);
    disp(['Signal-to-Noise ratio standard deviation of simulated signals = ' num2str(std(StN))]);
    
    plot_tr.f=f_simu;
    plot_tr.on=plot_tree;
    plot_tr.fe=Fe;
    plot_tr.f_cl=f_class;
end







if nargin == 2
    dwt_without_optim = 1;
    optim_fisher = 0;
    optim_pce = 0;
    optim_basis = 0;
end


%On mélange les essais pour ne pas avoir tous les signaux d'une même classe
%ensemble
ind_perm = randperm(length(y_learn)); %Vecteur de permutation
signal_learn = signal_learn(ind_perm,:,:);
y_learn = y_learn(ind_perm);

%Paramètres de l'ensemble d'apprentissage
[nb_trials,nb_samples,nb_chans] = size(signal_learn);
chans=1:nb_chans;

if sim_sig==0
    N_dec_max = floor(log2(nb_samples));
end
N_dec = N_dec_max;
disp(['depth of the DWT tree :' num2str(N_dec)]);
ind_des = [1:4]; %indice des marginales ou coeff dwt utilisés. Si vide => toutes
if isempty(ind_des)
    ind_des = 1:N_dec+1;
end



tic
if optim_fisher || optim_pce || optim_algo
    %on calcul les descripteurs de tous les individus pour
    %tous les paramètres d'ondelettes possibles
    
%------------------------------------------------------------------
%THE FUNCTION    

if optim_algo
    ind_des = 1:N_dec+1;
end

    [features_dic,h_dic,g_dic] = create_dic(signal_learn,N_dec,ind_des);

%------------------------------------------------------------------

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

if ~exist('nb_subsets','var')
    nb_subsets = 120;
end
    
y_estim_ini = [];
y_real_ini = [];
y_estim_pce = [];
y_real_pce = [];
y_estim_fisher = [];
y_real_fisher = [];
y_estim_bb = [];
y_real_bb = [];
y_estim_bb_acp = [];
y_real_bb_acp = [];
y_estim_algo = [];
y_real_algo = [];
subset_length = floor(length(y_learn)/nb_subsets);
critere_ini = zeros(nb_subsets,1);
critere_pce = zeros(nb_subsets,1);
critere_fisher = zeros(nb_subsets,1);

if dwt_without_optim
    disp('svm without any optimization');
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
    
    for j=N_dec+2-length(ind_des):N_dec
        ind_M(1,j-(N_dec+1-length(ind_des)))=j;
        ind_M(2,j-(N_dec+1-length(ind_des)))=1;
        if j==N_dec
            ind_M(1,j-(N_dec-length(ind_des)))=j;
            ind_M(2,j-(N_dec-length(ind_des)))=0;
        end
    end
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=0;
    function_acp(features_reshape_ini_app(:,1:size(features_reshape_ini_app,2)/length(chans)),y_ini_app,plot_,M_phi,sub_name);
    title('Dwt sans optim')
    drawnow
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_ini{i,j}='X';
                else
                    p_ini{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_ini{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_ini{i,j} = num2str(100*sum(y_estim_ini == j-1 & y_real_ini == i-1)/sum(y_real_ini == i-1));
                end
            end
        end
    end
    
    disp(p_ini);
    disp(['mean of optim citerion : ' num2str(mean(critere_ini))]);
    c_q_ini = sum(y_estim_ini ~= y_real_ini)/length(y_estim_ini);
    disp(['PCE_db2 : ' num2str(c_q_ini)]);
    disp('-----------------------------------------------------------------');
    
    clk=clock;
    if plot_
        saveas(gcf-1,['rec_graph/distrib_ini_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf,['rec_graph/correl_circle_ini_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
    end;
    
end

if optim_pce
    disp('svm with optimization by pce of the mother wavelet');
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

    for j=N_dec+2-length(ind_des):N_dec
        ind_M(1,j-(N_dec+1-length(ind_des)))=j;
        ind_M(2,j-(N_dec+1-length(ind_des)))=1;
        if j==N_dec
            ind_M(1,j-(N_dec-length(ind_des)))=j;
            ind_M(2,j-(N_dec-length(ind_des)))=0;
        end
    end
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=1;
    function_acp(features_reshape_pce_app(:,1:size(features_reshape_pce_app,2)/length(chans)),y_pce_app,plot_,M_phi,sub_name);
    title('Dwt optim Fisher')
    drawnow   
    
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_pce{i,j}='X';
                else
                    p_pce{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_pce{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_pce{i,j} = num2str(100*sum(y_estim_pce == j-1 & y_real_pce == i-1)/sum(y_real_pce == i-1));
                end
            end
        end
    end
    
    disp(p_pce);
    disp(['mean of optim citerion : ' num2str(mean(critere_pce))]);
    c_q_pce = sum(y_estim_pce ~= y_real_pce)/length(y_estim_pce);
    disp(['pce_pce : ' num2str(c_q_pce)]);
    disp('-----------------------------------------------------------------');
    
    clk=clock;
    if plot_
        saveas(gcf-1,['rec_graph/distrib_pce_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf,['rec_graph/correl_circle_pce_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
    end
end

%--------------------------------------------------------------------------
if optim_fisher
    disp('svm with modification of the form of mother wavelet by fisher criterion');
    wbar = waitbar(0,'Computation...');
    tic
    
    [h_ini, g_ini] = wfilters('db2');
    h_ini = ones(nb_chans,1)*h_ini;
    g_ini = ones(nb_chans,1)*g_ini;
    % features_ini = zeros(nb_trials,length(ind_des),nb_chans);
    % features_ini = zeros(nb_trials,nb_samples,nb_chans);
    for i=1:nb_trials
        for k=1:nb_chans
            features_ini(i,:,k) = calc_features(signal_learn(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
        end
    end
    
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar);
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        features_ini_app = features_ini(ind_app,:,:);
        y_fisher_app = y_learn(ind_app);
        [features_fisher_app,h_fisher,g_fisher,critere_fisher(ii)] = ...
            optim_wavelet_multi_cv(y_fisher_app, features_ini_app, h_ini, g_ini, features_dic(ind_app,:,:,:), h_dic, g_dic, 'fisher');
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

    for j=N_dec+2-length(ind_des):N_dec
        ind_M(1,j-(N_dec+1-length(ind_des)))=j;
        ind_M(2,j-(N_dec+1-length(ind_des)))=1;
        if j==N_dec
            ind_M(1,j-(N_dec-length(ind_des)))=j;
            ind_M(2,j-(N_dec-length(ind_des)))=0;
        end
    end
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=0;
    function_acp(features_reshape_fisher_app(:,1:size(features_reshape_fisher_app,2)/length(chans)),y_fisher_app,plot_,M_phi,sub_name);
    title('Dwt optim Fisher')
    drawnow
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_fisher{i,j}='X';
                else
                    p_fisher{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_fisher{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_fisher{i,j} = num2str(100*sum(y_estim_fisher == j-1 & y_real_fisher == i-1)/sum(y_real_fisher == i-1));
                end
            end
        end
    end
    
    disp(p_fisher);
    disp(['mean of optim citerion : ' num2str(mean(critere_fisher))]);
    c_q_fisher = sum(y_estim_fisher ~= y_real_fisher)/length(y_estim_fisher);
    disp(['PCE_fisher : ' num2str(c_q_fisher)]);
    disp('-----------------------------------------------------------------');
    
    clk=clock;
    if plot_
        saveas(gcf-1,['rec_graph/distrib_fish_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf,['rec_graph/correl_circle_fish_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
    end

end


if optim_basis
    disp('svm on features obtained with the optimazed wavelet basis');
    disp(['wavelet used for basis calculation: ' wname_ini ]);
    %N_dec = N_dec-2;
%     N_dec = 3;
%     dic_features_bb = set_wp_marg(signal_learn(:,:,1),wname_ini,N_dec);% vecteur contenant pour chaque essai de chaque voie les coefficient M
    dic_features_bb = set_wp_marg(signal_learn,wname_ini,N_dec);
    k = 0.2;

    if subset_
        wbar = waitbar(0,'computation...');
        tic
        for ii=1:nb_subsets
            waitbar(ii/nb_subsets,wbar);
            ind_test = (ii-1)*subset_length+1:ii*subset_length;
            ind_app = setdiff(1:length(y_learn),ind_test);
            features_bb_app = dic_features_bb(ind_app,:,:);
            features_bb_test = dic_features_bb(ind_test,:,:);
            y_bb_app = y_learn(ind_app);
            [basis_bb, fisher_tree] = best_basis_search(features_bb_app,...
                y_bb_app,N_dec,k,plot_tr,sim_sig,restric);     
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
    else
        a=length(dic_features_bb);
        features_bb_app = dic_features_bb(1:floor(size(dic_features_bb,1)/10),:,:);
        features_bb_test = dic_features_bb(floor(size(dic_features_bb,1)/10)+1:end,:,:);
        y_bb_app = y_learn(1:floor(size(dic_features_bb,1)/10));
        [basis_bb, fisher_tree] = best_basis_search(features_bb_app,...
            y_bb_app,N_dec,k,plot_tr,sim_sig,restric);     
        %     disp(sum(basis_bb == 1));
        features_reshape_bb_app = reshape(features_bb_app(:,boolean(basis_bb),:),floor(size(dic_features_bb,1)/10),[]);
        features_reshape_bb_test = reshape(features_bb_test(:,boolean(basis_bb),:),size(dic_features_bb,1)-floor(size(dic_features_bb,1)/10),[]);
        y_real_bb = [y_real_bb; y_learn(floor(size(dic_features_bb,1)/10)+1:end)];

        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);

        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_bb = svm_learning(features_reshape_bb_app,y_bb_app,[],nb_class,kernelparam);
        y_estim_bb = [y_estim_bb; test_class_ovr(features_reshape_bb_test, sep_bb, nb_class, features_reshape_bb_app, y_bb_app, kernelparam)];

    end 
    toc
    close(wbar)   
    
    %%%%%%%%%%%% DRAWING OF TRIALS INTO MARGINAL SPACES %%%%%%%%%%%%%%%%%%%
    
    plot_trial=0;
    
    if plot_trial
        
        inds_draw=findmarg(fisher_tree,basis_bb,N_dec);

        if size(inds_draw,2)<5
            marg=cell(1,size(inds_draw,2)-1);
            for i=1:size(inds_draw,2)-1
                marg{i}=[node(inds_draw(1,i),inds_draw(2,i)) ,node(inds_draw(1,i+1),inds_draw(2,i+1))];
            %     marg={[node(9,0) ,node(8,1)];[node(9,0),node(6,1)];};

            end
        else 

            marg=cell(1,4);
            for i=1:4
                marg{i}=[node(inds_draw(1,i),inds_draw(2,i)) ,node(inds_draw(1,i+1),inds_draw(2,i+1))];
            %     marg={[node(9,0) ,node(8,1)];[node(9,0),node(6,1)];};

            end

        end
        trial_drawing(dic_features_bb(:,:,1),y_learn,marg,N_dec);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k=find(basis_bb);
    ind_M=nds2ind(k,N_dec);
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=1;
    function_acp(features_reshape_bb_app(:,1:size(features_reshape_bb_app,2)/length(chans)),y_bb_app,plot_,M_phi,sub_name);
    drawnow
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_bb{i,j}='X';
                else
                    p_bb{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_bb{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_bb{i,j} = num2str(100*sum(y_estim_bb == j-1 & y_real_bb == i-1)/sum(y_real_bb == i-1));
                end
            end
        end
    end
    
    disp(p_bb);
    
    c_q_bb = sum(y_estim_bb ~= y_real_bb)/length(y_estim_bb);
    disp(['pce_bb : ' num2str(c_q_bb)]);
    disp('-----------------------------------------------------------------');
    
    clk=clock;
    if plot_
        if exist('gcf-6','var')
            ssaveas(gcf-6,['rec_graph/tree_bb_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-5','var')
            saveas(gcf-5,['rec_graph/1st_space_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-4','var')
            saveas(gcf-4,['rec_graph/2nd_space_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-3','var')
            saveas(gcf-3,['rec_graph/3rd_space_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-2','var')
            saveas(gcf-2,['rec_graph/4th_space_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-1','var')
            saveas(gcf-1,['rec_graph/distrib_bb_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf','var')
            saveas(gcf,['rec_graph/correl_circle_bb_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        end
    end
end

if optim_basis_ACP
    disp('svm with decomposition basis + ACP');
    disp(['wavelet used for basis calculation: ' wname_ini ]);
    disp([num2str(kept_axes) ' axes of the ACP kept for the svm']);
    %N_dec = N_dec-2;
%     N_dec = 3;
    dic_features_bb = set_wp_marg(signal_learn,wname_ini,N_dec);% vecteur contenant pour chaque essai de chaque voie les coefficient M
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
        [basis_bb, fisher_tree] = best_basis_search(features_bb_app,...
                y_bb_app,N_dec,k,plot_tr,sim_sig,restric);
        %     disp(sum(basis_bb == 1));
        features_reshape_bb_app = reshape(features_bb_app(:,boolean(basis_bb),:),length(ind_app),[]);
        features_reshape_bb_test = reshape(features_bb_test(:,boolean(basis_bb),:),length(ind_test),[]);
        y_real_bb_acp = [y_real_bb_acp; y_learn(ind_test)];
        
        [features_acp_app,U,x_mean,x_std]=function_acp(features_reshape_bb_app,y_bb_app,0);
        m=size(features_reshape_bb_test,1);
        x_test_centre=(features_reshape_bb_test - x_mean(ones(m,1),:))./x_std(ones(m,1),:);
        features_acp_test=x_test_centre*U;
        
        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);
        
        if kept_axes>size(features_acp_app,2)
            xapp=features_acp_app;
            xtest=features_acp_test;
        else
            xapp=features_acp_app(:,1:kept_axes);
            xtest=features_acp_test(:,1:kept_axes);
        end
        
        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_bb = svm_learning(xapp,y_bb_app,[],nb_class,kernelparam);
        y_estim_bb_acp = [y_estim_bb_acp; test_class_ovr(xtest, sep_bb, nb_class, xapp, y_bb_app, kernelparam)];
        
    end
    % keyboard
    toc
    close(wbar)
    
        %%%%%%%%%%%% DRAWING OF TRIALS INTO MARGINAL SPACES %%%%%%%%%%%%%%%%%%%
    
    plot_trial=0;
    
    if plot_trial 
        inds_draw=findmarg(fisher_tree,basis_bb,N_dec);

        if size(inds_draw,2)<5
            marg=cell(1,size(inds_draw,2)-1);
            for i=1:size(inds_draw,2)-1
                marg{i}=[node(inds_draw(1,i),inds_draw(2,i)) ,node(inds_draw(1,i+1),inds_draw(2,i+1))];
            %     marg={[node(9,0) ,node(8,1)];[node(9,0),node(6,1)];};

            end
        else 

            marg=cell(1,4);
            for i=1:4
                marg{i}=[node(inds_draw(1,i),inds_draw(2,i)) ,node(inds_draw(1,i+1),inds_draw(2,i+1))];
            %     marg={[node(9,0) ,node(8,1)];[node(9,0),node(6,1)];};

            end

        end
        trial_drawing(dic_features_bb(:,:,1),y_learn,marg,N_dec);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    k=find(basis_bb);
    ind_M=nds2ind(k,N_dec);
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=0;
    function_acp(features_reshape_bb_app(:,1:size(features_reshape_bb_app,2)/length(chans)),y_bb_app,plot_,M_phi,sub_name);
    title('Best Basis + ACP')
    drawnow
    
    %création du tableau des erreurs
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_bb_acp{i,j}='X';
                else
                    p_bb_acp{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_bb_acp{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_bb_acp{i,j} = num2str(100*sum(y_estim_bb_acp == j-1 & y_real_bb_acp == i-1)/sum(y_real_bb_acp == i-1));
                end
            end
        end
    end
    
    disp(p_bb_acp);
    
    % disp(['mean of optim citerion : ' num2str(mean(critere_bb))]);
    c_q_bb_acp = sum(y_estim_bb_acp ~= y_real_bb_acp)/length(y_estim_bb_acp);
    disp(['pce_bb_acp : ' num2str(c_q_bb_acp)]);
    disp('-----------------------------------------------------------------');
    clk=clock;
    if plot_
        if exist('gcf-6','var')
            ssaveas(gcf-6,['rec_graph/tree_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-5','var')
            saveas(gcf-5,['rec_graph/1st_space_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-4','var')
            saveas(gcf-4,['rec_graph/2nd_space_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-3','var')
            saveas(gcf-3,['rec_graph/3rd_space_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-2','var')
            saveas(gcf-2,['rec_graph/4th_space_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf-1','var')
            saveas(gcf-1,['rec_graph/distrib_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        elseif exist('gcf','var')
            saveas(gcf,['rec_graph/correl_circle_bb_acp_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        end
    end
end

if optim_algo

    disp('svm using algorythme based on modification of the form of mother wavelet by fisher criterion and then best basis');
    wbar = waitbar(0,'Computation...');
    tic
    
    [h_ini, g_ini] = wfilters('db2');
    h_ini = ones(nb_chans,1)*h_ini;
    g_ini = ones(nb_chans,1)*g_ini;
    % features_ini = zeros(nb_trials,length(ind_des),nb_chans);
    % features_ini = zeros(nb_trials,nb_samples,nb_chans);
    for i=1:nb_trials
        for k=1:nb_chans
            features_ini(i,:,k) = calc_features(signal_learn(i,:,k),N_dec,ind_des,h_ini(k,:),g_ini(k,:));
        end
    end
    
    for ii=1:nb_subsets
        waitbar(ii/nb_subsets,wbar);
        ind_test = (ii-1)*subset_length+1:ii*subset_length;
        ind_app = setdiff(1:length(y_learn),ind_test);
        features_ini_app = features_ini(ind_app,:,:);
        y_algo_app = y_learn(ind_app);
        
        %%%%%%%%%%%%%%%%%%%%%%% SEARCH OF MOTHER WAVELET %%%%%%%%%%%%%%%%%%
        
        [useless,h_mother,g_mother] = ...
            optim_wavelet_multi_cv(y_algo_app, features_ini_app, h_ini, g_ini, features_dic(ind_app,:,:,:), h_dic, g_dic, 'fisher');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        disp(['wavelet used for basis calculation: mother wavelet' ]);

    %N_dec = N_dec-2;
%     N_dec = 3;
        dic_features_algo = set_wp_marg(signal_learn,wname_ini,N_dec,h_mother,g_mother);% vecteur contenant pour chaque essai de chaque voie les coefficient M
        k = 0.2;

        features_algo_app = dic_features_algo(ind_app,:,:);
        features_algo_test = dic_features_algo(ind_test,:,:);
        [basis_algo, fisher_tree_algo] = best_basis_search(features_algo_app,...
            y_algo_app,N_dec,k,plot_tr,sim_sig,restric);     
        %     disp(sum(basis_algo == 1));
        features_reshape_algo_app = reshape(features_algo_app(:,boolean(basis_algo),:),length(ind_app),[]);
        features_reshape_algo_test = reshape(features_algo_test(:,boolean(basis_algo),:),length(ind_test),[]);
        y_real_algo = [y_real_algo; y_learn(ind_test)];

        kernelparam.ktype = 1;
        kernelparam.kscale = 0.2;
        nb_class = max(y_learn);

        %initialisation de la séparatrice avec tous les élément de l'apprentissage
        sep_algo = svm_learning(features_reshape_algo_app,y_algo_app,[],nb_class,kernelparam);
        y_estim_algo = [y_estim_algo; test_class_ovr(features_reshape_algo_test, sep_algo, nb_class, features_reshape_algo_app, y_algo_app, kernelparam)];



    end


    %%%%%%%%%%%% DRAWING OF TRIALS INTO MARGINAL SPACES %%%%%%%%%%%%%%%%%%%
    
    inds_draw=findmarg(fisher_tree_algo,basis_algo,N_dec);
    
    marg=cell(1,4);
    for i=1:4
        marg{i}=[node(inds_draw(1,i),inds_draw(2,i)) ,node(inds_draw(1,i+1),inds_draw(2,i+1))];
    %     marg={[node(9,0) ,node(8,1)];[node(9,0),node(6,1)];};
        
    end
    trial_drawing(dic_features_algo(:,:,1),y_learn,marg,N_dec);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k=find(basis_algo);
    ind_M=nds2ind(k,N_dec);
    M_phi=cell(length(ind_M),1);
    for i=1:length(M_phi)  
        %création du texte associé à chaque point du cercle de corrélation
        M_phi{i}=['  M\phi_{' num2str(ind_M(1,i)) '}^{' num2str(ind_M(2,i)) '}' ];
    end
    plot_=1;
    function_acp(features_reshape_algo_app(:,1:size(features_reshape_algo_app,2)/length(chans)),y_algo_app,plot_,M_phi,sub_name);
    drawnow
    for i=1:nb_class+1
        for j=1:nb_class+1
            if i==1
                if j==1
                    p_algo{i,j}='X';
                else
                    p_algo{i,j}=['c_estim = ' num2str(j-1)];
                end
            else
                if j==1
                    p_algo{i,j}=['c_real = ' num2str(i-1)];
                else
                    p_algo{i,j} = num2str(100*sum(y_estim_algo == j-1 & y_real_algo == i-1)/sum(y_real_algo == i-1));
                end
            end
        end
    end
    
    disp(p_algo);
    
    c_q_algo = sum(y_estim_algo ~= y_real_algo)/length(y_estim_algo);
    disp(['pce_algo : ' num2str(c_q_algo)]);
    disp('-----------------------------------------------------------------');
    
    clk=clock;
    
    if plot_
        saveas(gcf-6,['rec_graph/tree_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf-5,['rec_graph/1st_space_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf-4,['rec_graph/2nd_space_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf-3,['rec_graph/3rd_space_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf-2,['rec_graph/4th_space_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf-1,['rec_graph/distrib_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
        saveas(gcf,['rec_graph/correl_circle_algo_sim_' num2str(sim_sig) '_if_nosim_name_' sub_name '_date_' num2str(clk(3)) '-' num2str(clk(2)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
    end

end
end