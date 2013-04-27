clear all
close all
clc

N_start = 2200;
N_end = 3000;
filtering = 0;
global file_name
file_name = 'resultat_classif.tex';
C = 1;
%SUBJECT 1 OLD
% %volontary
% File = {'gediminas_voluntary'};
% rep = 'subject1_old\voluntary';
% chan = 1:6;
%unvolontary
% nom_sujet = 'SUBJECT 1 OLD';
% disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'w');
% fprintf(fid,'\\newpage\n');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
% File = {'gediminas_unvoluntary'};
% rep ='subject1_old\unvoluntary';
% chan = 1:6;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% approche_classif(class_ep,chan,C,'concat',[],[],1);
% approche_classif(class_ep,chan,C,'dwt',[],[],1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,1,1);

%SUBJECT 1 NEW
% %volontary
% File = {'gediminas_voluntary'};
% rep = 'subject1_new\voluntary';
% % unvolontary
nom_sujet = 'SUBJECT 1 NEW';
disp(['sujet : ', nom_sujet]);
fid = fopen(file_name,'a');
fprintf(fid,'\\newpage\n');
fprintf(fid,'\\section*{%s}\n',nom_sujet);
fclose(fid);
File = {'gediminas_unvoluntary'};
rep ='subject1_new\unvoluntary';
chan = [1 6];
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif(class_ep,chan,C,'concat',[],[],1);
% approche_classif(class_ep,chan,C,'dwt',[],[],1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,1,1);

%SUBJECT 2 OLD
% %volontary
% File = {'marija_voluntary'};
% rep = 'subject2_old\voluntary';
% %unvolontary
% nom_sujet = 'SUBJECT 2 OLD';
% disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'a');
% fprintf(fid,'\\newpage\n');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
% File = {'marija_unvoluntary'};
% rep = 'subject2_old\unvoluntary';
% chan = 3;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% approche_classif(class_ep,chan,C,'concat',[],[],1);
% approche_classif(class_ep,chan,C,'dwt',[],[],1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,1,1);

%SUBJECT 2 NEW
% %volontary
% File = {'marija_voluntary'};
% rep = 'subject2_new\voluntary';
% %unvolontary
% nom_sujet = 'SUBJECT 2 NEW';
% disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'a');
% fprintf(fid,'\\newpage\n');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
% File = {'marija_unvoluntary'};
% rep = 'subject2_new\unvoluntary';
% chan = 3;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% approche_classif(class_ep,chan,C,'concat',[],[],1);
% approche_classif(class_ep,chan,C,'dwt',[],[],1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,1,1);


%SUBJECT 3
% %volontary
% File = {'mustafa_voluntary'};
% rep = 'subject3\voluntary';
%unvolontary
% nom_sujet = 'SUBJECT 3';
% disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'a');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
% File = {'mustafa_unvoluntary'};
% rep = 'subject3\unvoluntary';
% chan = 3;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% approche_classif(class_ep,chan,C,'concat',[],[],1);
% approche_classif(class_ep,chan,C,'dwt',[],[],1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,0,1);
% approche_classif(class_ep,chan,C,'spatial_filter',0,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',1,1,1);
% approche_classif(class_ep,chan,C,'spatial_filter',2,1,1);

%SUBJET 4
% %volontary
% File = {'andersen_voluntary'};
% rep = 'subject4\voluntary';
% chan = 1:6;

% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% algo_test_EP_loo(class_ep,chan,0,Fe,filtering,N_start,N_end);
% algo_test_EP_loo(class_ep,chan,1,Fe,filtering,N_start,N_end);
% algo_test_EP_loo(class_ep,chan,2,Fe,filtering,N_start,N_end);

%Duy
nom_sujet = 'DUY';
disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'a');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
File = {'unvoluntary'};
rep = 'duy';
chan = 1:6;
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
approche_classif(class_ep,chan,C,'concat',[],[],0);
approche_classif(class_ep,chan,C,'dwt',[],[],0);
approche_classif(class_ep,chan,C,'spatial_filter',0,0,0);
approche_classif(class_ep,chan,C,'spatial_filter',1,0,0);
approche_classif(class_ep,chan,C,'spatial_filter',2,0,0);
approche_classif(class_ep,chan,C,'spatial_filter',0,1,0);
approche_classif(class_ep,chan,C,'spatial_filter',1,1,0);
approche_classif(class_ep,chan,C,'spatial_filter',2,1,0);

%Jakob
nom_sujet = 'JAKOB';
disp(['sujet : ', nom_sujet]);
fid = fopen(file_name,'a');
fprintf(fid,'\\section*{%s}\n',nom_sujet);
fclose(fid);
File = {'unvoluntary'};
rep = 'jakob';
chan = 1:6;
[class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
