clear all
close all
clc

t_start = 0.15;
t_end = 0.70;
N_start = round((2+t_start)*1024);
N_end = round((2+t_end)*1024);
% N_start = 2200;
% N_end = 2800;
% possible_chan = {1,6,[1 6]};
filtering = 1;
global file_name
file_name = 'resultat.tex';
%SUBJECT 1 OLD
% %volontary
% File = {'gediminas_voluntary'};
% rep = 'subject1_old\voluntary';
% chan = 1:6;
% %unvolontary
nom_sujet = 'SUBJECT 1 OLD';
disp(['sujet : ', nom_sujet]);
% fid = fopen(file_name,'w');
% fprintf(fid,'\\newpage\n');
% fprintf(fid,'\\section*{%s}\n',nom_sujet);
% fclose(fid);
File = {'gediminas_unvoluntary'};
rep ='subject1_old\unvoluntary';
chan = [1 2 3 4 5 6];
[class_mvt , class_ep, chan_label, Fe] = read_imran_file(File,rep);
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
    algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% end
% algo_test_EP_loo_recal(class_ep,chan,1,0,1)
% algo_test_EP_loo_recal(class_ep,chan,2,0,1)
% algo_test_EP_loo_recal(class_ep,chan,0,1,1)
% algo_test_EP_loo_recal(class_ep,chan,1,1,1)
% algo_test_EP_loo_recal(class_ep,chan,2,1,1)
%SUBJECT 1 NEW
% %volontary
% File = {'gediminas_voluntary'};
% rep = 'subject1_new\voluntary';
% % unvolontary
% nom_sujet = 'SUBJECT 1 NEW';
% disp(['sujet : ', nom_sujet]);
% % fid = fopen(file_name,'a');
% % fprintf(fid,'\\newpage\n');
% % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % fclose(fid);
% File = {'gediminas_unvoluntary'};
% rep ='subject1_new\unvoluntary';
% % chan = 1:6;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
%     algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% % end
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1)
% 
% %SUBJECT 2 OLD
% % %volontary
% % File = {'marija_voluntary'};
% % rep = 'subject2_old\voluntary';
% % %unvolontary
% nom_sujet = 'SUBJECT 2 OLD';
% disp(['sujet : ', nom_sujet]);
% % fid = fopen(file_name,'a');
% % fprintf(fid,'\\newpage\n');
% % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % fclose(fid);
% File = {'marija_unvoluntary'};
% rep = 'subject2_old\unvoluntary';
% % chan = 1:6;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
%     algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% % end
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1)
% 
% %SUBJECT 2 NEW
% % %volontary
% % File = {'marija_voluntary'};
% % rep = 'subject2_new\voluntary';
% % %unvolontary
% nom_sujet = 'SUBJECT 2 NEW';
% disp(['sujet : ', nom_sujet]);
% % fid = fopen(file_name,'a');
% % fprintf(fid,'\\newpage\n');
% % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % fclose(fid);
% File = {'marija_unvoluntary'};
% rep = 'subject2_new\unvoluntary';
% % chan = [1 2 3 4 5 6];
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
%     algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% % end
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1)
% 
% %SUBJECT 3
% % %volontary
% % File = {'mustafa_voluntary'};
% % rep = 'subject3\voluntary';
% %unvolontary
% nom_sujet = 'SUBJECT 3';
% disp(['sujet : ', nom_sujet]);
% % fid = fopen(file_name,'a');
% % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % fclose(fid);
% File = {'mustafa_unvoluntary'};
% rep = 'subject3\unvoluntary';
% % chan = 1:6;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
%     algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% % end
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1)
% 
% %SUBJET 4
% % %volontary
% % File = {'andersen_voluntary'};
% % rep = 'subject4\voluntary';
% % chan = 1:6;
% 
% % [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % algo_test_EP_loo(class_ep,chan,0,Fe,filtering,N_start,N_end);
% % algo_test_EP_loo(class_ep,chan,1,Fe,filtering,N_start,N_end);
% % algo_test_EP_loo(class_ep,chan,2,Fe,filtering,N_start,N_end);
% 
% %DUY
% % nom_sujet = 'DUY';
% % disp(['sujet : ', nom_sujet]);
% % % fid = fopen(file_name,'a');
% % % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % % fclose(fid);
% % File = {'unvoluntary'};
% % rep = 'duy';
% % chan = 1:6;
% % [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1,Fe,1)
% 
% % Jakob
% % nom_sujet = 'JAKOB';
% % disp(['sujet : ', nom_sujet]);
% % % fid = fopen(file_name,'a');
% % % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % % fclose(fid);
% % File = {'unvoluntary'};
% % rep = 'jakob';
% % chan = 1:6;
% % [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
% % % algo_test_EP_loo_recal(class_ep,chan,1,0,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1,Fe,1)
% % % algo_test_EP_loo_recal(class_ep,chan,1,1,1,Fe,1)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1,Fe,1)
% 
% nom_sujet = 'DANTAN';
% disp(['sujet : ', nom_sujet]);
% % fid = fopen(file_name,'a');
% % fprintf(fid,'\\section*{%s}\n',nom_sujet);
% % fclose(fid);
% File = {'unvoluntary'};
% rep = 'DanTan';
% % chan = 1:6;
% [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
%     algo_test_EP_loo_recal(class_ep,chan,0,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,1,0,0,Fe,1,N_start,N_end)
%     algo_test_EP_loo_recal(class_ep,chan,2,0,0,Fe,1,N_start,N_end)
% % end
% % algo_test_EP_loo_recal(class_ep,chan,1,0,1,Fe,0)
% % algo_test_EP_loo_recal(class_ep,chan,2,0,1,Fe,0)
% % algo_test_EP_loo_recal(class_ep,chan,0,1,1,Fe,0)
% % algo_test_EP_loo_recal(class_ep,chan,1,1,1,Fe,0)
% % algo_test_EP_loo_recal(class_ep,chan,2,1,1,Fe,0)
