clear all
close all
clc

% filtering = 0;
global file_name
file_name = 'resultat_classif.tex';
C = 10;
filtering = 0;
chan = 6;


%SUBJECT 1 OLD
nom_sujet = 'SUBJECT 1 OLD';
disp(['sujet : ', nom_sujet]);
File = {'gediminas_unvoluntary'};
rep ='subject1_old\unvoluntary';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%SUBJECT 1 NEW
nom_sujet = 'SUBJECT 1 NEW';
disp(['sujet : ', nom_sujet]);
File = {'gediminas_unvoluntary'};
rep ='subject1_new\unvoluntary';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%SUBJECT 2 OLD
nom_sujet = 'SUBJECT 2 OLD';
disp(['sujet : ', nom_sujet]);
File = {'marija_unvoluntary'};
rep = 'subject2_old\unvoluntary';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%SUBJECT 2 NEW
nom_sujet = 'SUBJECT 2 NEW';
disp(['sujet : ', nom_sujet]);
File = {'marija_unvoluntary'};
rep = 'subject2_new\unvoluntary';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%SUBJECT 3
nom_sujet = 'SUBJECT 3';
disp(['sujet : ', nom_sujet]);
File = {'mustafa_unvoluntary'};
rep = 'subject3\unvoluntary';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%Duy
nom_sujet = 'DUY';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'duy';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%Jakob
nom_sujet = 'JAKOB';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'jakob';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

%Dantan
nom_sujet = 'DANTAN';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'DanTan';
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');

% for i=1:6
%     disp(['chan : ' num2str(i)])
%     chan = i;
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'rakoto');
% % approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'lda');
% approche_classif_mvt(class_mvt,chan,C,filtering,Fe,'bayes');
% end

