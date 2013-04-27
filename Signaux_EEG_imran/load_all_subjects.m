clear all
close all
clc

%SUBJECT 1 OLD
nom_sujet = 'SUBJECT 1 OLD';
disp(['sujet : ', nom_sujet]);
File = {'gediminas_unvoluntary'};
rep ='subject1_old\unvoluntary';
[class_mvt, class_ep_1, chan_label, Fe] = read_imran_file(File,rep);


%SUBJECT 1 NEW
nom_sujet = 'SUBJECT 1 NEW';
disp(['sujet : ', nom_sujet]);
File = {'gediminas_unvoluntary'};
rep ='subject1_new\unvoluntary';
[class_mvt, class_ep_2, chan_label, Fe] = read_imran_file(File,rep);


%SUBJECT 2 OLD
nom_sujet = 'SUBJECT 2 OLD';
disp(['sujet : ', nom_sujet]);
File = {'marija_unvoluntary'};
rep = 'subject2_old\unvoluntary';
[class_mvt, class_ep_3, chan_label, Fe] = read_imran_file(File,rep);


%SUBJECT 2 NEW
nom_sujet = 'SUBJECT 2 NEW';
disp(['sujet : ', nom_sujet]);
File = {'marija_unvoluntary'};
rep = 'subject2_new\unvoluntary';
[class_mvt, class_ep_4, chan_label, Fe] = read_imran_file(File,rep);


%SUBJECT 3
nom_sujet = 'SUBJECT 3';
disp(['sujet : ', nom_sujet]);
File = {'mustafa_unvoluntary'};
rep = 'subject3\unvoluntary';
[class_mvt, class_ep_5, chan_label, Fe] = read_imran_file(File,rep);

%Duy
nom_sujet = 'DUY';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'duy';
[class_mvt, class_ep_6, chan_label, Fe] = read_imran_file(File,rep);

%Jakob
nom_sujet = 'JAKOB';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'jakob';
[class_mvt, class_ep_7, chan_label, Fe] = read_imran_file(File,rep);


%Dantan
nom_sujet = 'DANTAN';
disp(['sujet : ', nom_sujet]);
File = {'unvoluntary'};
rep = 'DanTan';
[class_mvt, class_ep_8, chan_label, Fe] = read_imran_file(File,rep);

