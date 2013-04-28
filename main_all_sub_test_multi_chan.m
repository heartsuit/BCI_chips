clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of loop you want to do for the training set definition (to do the 
%leave-one-out method, ut this param to 120)
nb_subsets = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% OPTIMIZATION METHOD CHOICE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% NO OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optim.dwt_without_optim = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%% PCE OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optim.optim_pce = 0;

%%%%%%%%%%%%%%%%%%%%%%% FISHER OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optim.optim_fisher = 0;

%%%%%%%%%%%%%%%%%%%%%%%% DWPT OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optim.optim_basis = 1;

%because our informations are in the low frequency domain we can choose to
%keep only the first fourth of the tree decomposition (restric=1)
restric=0;

%to be used if you want to plot the tree of marginals (WARNING: also used
%in the best basis + PCA method)
plot_tree=1;

%%%%%%%%%%%%%%%%%% DWPT OPTIMIZATION + PCA REFINEMENT %%%%%%%%%%%%%%%%%%%%%%

optim.optim_basis_ACP = 0;

%number of axes from the acp kept for the SVM input
kept_axes=5;

%%%%%%%%%%%%%%%%%%%%%%%%% HYBRID METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optim.algo= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL PARAMETERIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if you want simulated signals put the following parameter on 1
simu_para.on=1;

%total number of signals that we will train and test
simu_para.nb_trials=120;

%number of classes we want to implement
simu_para.nb_class=5;

%maximum scale of the wavelet decomposition
simu_para.N_dec_max=6;

%time length of the signals
simu_para.time=2;

%paremeters on the shape of the simulated signals:
%the global shape: 
    %do we want to have the same form than the wavelet and
    %doing correlation (using db4 for the signal and the wavelet)
    simu_para.perfect=1;

    %do we want some moise in our simulated signals
    simu_para.sigma=0;
    
    %the different classes will everytime be chosen as close as possible
    %but we still can choose where we want the group of classes on the
    %frequency axis (offset from 0 Hz)
    simu_para.posi=0;
    
%attempt to reduce calculation time removing the leave-one-out option. To
%run the program without bugs, let this parameter to 1
simu_para.subset=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REAL SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filter the raw data with a Cheby2 with cut off freq at 45 Hz.
filtre=1; 

%plot raw data 
plot_raw=0; 

%plot fourrier transfor of raw data for one sample
plot_real_sig=0; 

%projection wavelet for DWPT marginal computation
wname_ini='db2'; 

%subjects you want to deal with (can be a vector following the row numbers of the vector filename_)
subjects=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




addpath(genpath(cd));

rand('state',0);

subject.filename = [];
subject.pathname = [];
subject.sub_name = [];


if simu_para.on
    disp('use of simulated signals');
    sub_name= 'Signaux simulés';
    
    disp(['subset number used for cross validation: ' num2str(nb_subsets)]);
    
    
    programme_test_descripteur_xval(subject,optim,plot_tree,...
        kept_axes,nb_subsets,simu_para,restric,filtre,wname_ini)
    disp('-----------------------------------------------------------------');

    
else
            
    filename_={'gediminas_unvoluntary.mat',...
        'gediminas_unvoluntary.mat',... 
        'marija_unvoluntary.mat',...
        'marija_unvoluntary.mat',...
        'mustafa_unvoluntary.mat',...
        'unvoluntary.mat',...
        'unvoluntary.mat',...
        'unvoluntary.mat'};
        
    pathname_={'Signaux_EEG_imran/subject1_new/Unvoluntary/',...
        'Signaux_EEG_imran/subject1_old/Unvoluntary/',...
        'Signaux_EEG_imran/subject2_new/unvoluntary/',...
        'Signaux_EEG_imran/subject2_old/Unvoluntary/',...
        'Signaux_EEG_imran/subject3/unvoluntary/',...
        'Signaux_EEG_imran/DanTAn/',...
        'Signaux_EEG_imran/duy/',...
        'Signaux_EEG_imran/jakob/'};
        
    sub_name_={'Subject 1 New unvolontary',...
        'Subject 1 Old unvolontary',...
        'Subject 2 New unvolontary',...
        'Subject 2 Old unvolontary',...
        'Subject 3 unvolontary',...
        'Subject Dantan unvolontary',...
        'Subject Duy unvolontary',...
        'Subject Jakob unvolontary'};
        
       
    for i=subjects
        
        subject.filename = filename_{i};
        subject.pathname = pathname_{i};
        subject.sub_name = sub_name_{i};
        disp(subject.sub_name);
        if filtre
            disp('filtre on');
        else
            disp('filtre off');
        end
        disp(['subset number used for cross validation: ' num2str(nb_subsets)]);

        programme_test_descripteur_xval(subject,optim,plot_tree,...
            kept_axes,nb_subsets,simu_para,restric,filtre,wname_ini,plot_raw,plot_real_sig)
        disp('-----------------------------------------------------------------');
    end 
end