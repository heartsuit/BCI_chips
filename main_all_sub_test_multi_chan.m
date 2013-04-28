
% clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%  COMPULSORY PARAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simu_para.on=0;

simu_para.nb_trials=120;
simu_para.nb_class=5;
%     simu_para.N_dec_max = floor(log2(length(N_selec)));
simu_para.N_dec_max=6;
simu_para.time=2;
simu_para.perfect=1; %db2
simu_para.sigma=0; %noise
simu_para.posi=0; %offset from 0Hz
simu_para.subset=1; %always one (works only for best basis)

filtre=1; %filter the raw data with a Cheby2 with cut off freq at 45 Hz.
plot_raw=0; %plot raw data synthesis
plot_real_sig=0; %plot fourrier transfor of raw data for one exp
wname_ini='db2'; %projection wavelet

optim.dwt_without_optim = 0; 
optim.optim_fisher = 0;
optim.optim_pce = 0;

%------------------------------------------------------------------------

optim.optim_basis = 1;
restric=0; %conservation of a tiny part on the left side ...

optim.optim_basis_ACP = 0;
kept_axes=5;    %nombre d'axe de l'acp retenu


optim.algo= 0;

plot_tree=1; %variable dessinant ou non l'arbre de la base retenue pour la DWPT

%%%%%%%%%%%%%%%%%%%%%%%%  ADVISED PARAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_subsets = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%     filename_={'gediminas_unvoluntary.mat',... 
%         'gediminas_voluntary.mat',...
%         'marija_unvoluntary.mat',...
%         'marija_voluntary.mat',...
%         'mustafa_unvoluntary.mat',...
%         'mustafa_voluntary.mat',...
%         'andersen_voluntary.mat'};
%         
%     pathname_={'Signaux_EEG_imran/subject1_new/Unvoluntary/',...
%         'Signaux_EEG_imran/subject1_new/voluntary/',...
%         'Signaux_EEG_imran/subject2_new/unvoluntary/',...
%         'Signaux_EEG_imran/subject2_new/voluntary/',...
%         'Signaux_EEG_imran/subject3/unvoluntary/',...
%         'Signaux_EEG_imran/subject3/voluntary/',...
%         'Signaux_EEG_imran/subject4/voluntary/'};
%         
%     sub_name_={'Subject 1 New unvolontary',...
%         'Subject 1 New volontary',...
%         'Subject 2 New unvolontary',...
%         'Subject 2 New volontary',...
%         'Subject 3 unvolontary',...
%         'Subject 3 volontary',...
%         'Subject 4 volontary'};
        
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
        
    
%     for i=1:length(filename_)   
    for i=4 %CHOICE!!!
        
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


% sub_name = 'Subject 1 Old';
% filename = 'gediminas_unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/subject1_old/Unvoluntary/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');

% sub_name = 'Subject 2 New';
% filename = 'marija_unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/subject2_new/unvoluntary/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% 
% sub_name = 'Subject 2 Old';
% filename = 'marija_unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/subject2_old/Unvoluntary/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% 
% sub_name = 'Subject 3';
% filename = 'mustafa_unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/subject3/unvoluntary/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% 
% sub_name = 'Dantan';
% filename = 'unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/DanTAn/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% 
% sub_name = 'Duy';
% filename = 'unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/duy/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% 
% sub_name = 'Jackob';
% filename = 'unvoluntary.mat';
% pathname = 'Signaux_EEG_imran/jakob/';
% disp(sub_name);
% figure('Name',sub_name);
% programme_test_descripteur_xval(filename,pathname,dwt_without_optim,optim_fisher,optim_pce,optim_basis)
% disp('-----------------------------------------------------------------');
% rmpath(genpath(cd));