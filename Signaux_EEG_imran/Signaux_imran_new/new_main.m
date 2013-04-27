% clear all
close all
clc

% filtering = 0;
global file_name
file_name = 'resultat_classif.tex';
C = 10;
filtering = 1;
% chan = [1];
% possible_chan = {1,6,[1 6]};
 possible_chan = {[1 6]};
t_start = 0.15;
t_end = 0.70;

% %SUBJECT 1 OLD
% nom_sujet = 'SUBJECT 1 OLD';
% disp(['sujet : ', nom_sujet]);
% File = {'gediminas_unvoluntary'};
% rep ='subject1_old\unvoluntary';
% [class_mvt_1, class_ep_1, chan_label, Fe] = read_imran_file(File,rep);

% w=0;
for i=1:length(possible_chan)
    disp(['chan : ' num2str(possible_chan{i})])
    chan = possible_chan{i};
    disp(['intervalle de temps : ' num2str([t_start t_end])])
    approche_classif_2(class_ep_1,chan,C,filtering,Fe,'rakoto',t_start,t_end);
%     filtering = 1;
%     approche_classif_2(class_ep_1,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% %     for t_start = 0:0.05:0.5
% %         for t_end = t_start+0.3:0.05:1
% %             w = w+1;
% %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% %             sub(1).chan{w} = chan;
% %             sub(1).t_start(w) = t_start;
% %             sub(1).t_end(w) = t_end;
% %             sub(1).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% %         end
% %     end
end

% %SUBJECT 1 NEW
% nom_sujet = 'SUBJECT 1 NEW';
% disp(['sujet : ', nom_sujet]);
% File = {'gediminas_unvoluntary'};
% rep ='subject1_new\unvoluntary';
% [class_mvt, class_ep_2, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
%     disp(['intervalle de temps : ' num2str([t_start t_end])])
%     approche_classif_2(class_ep_2,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(2).chan{w} = chan;
% % %             sub(2).t_start(w) = t_start;
% % %             sub(2).t_end(w) = t_end;
% % %             sub(2).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% end

% %SUBJECT 2 OLD
% nom_sujet = 'SUBJECT 2 OLD';
% disp(['sujet : ', nom_sujet]);
% File = {'marija_unvoluntary'};
% rep = 'subject2_old\unvoluntary';
% [class_mvt, class_ep_3, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
%     disp(['intervalle de temps : ' num2str([t_start t_end])])
%     approche_classif_2(class_ep_3,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(3).chan{w} = chan;
% % %             sub(3).t_start(w) = t_start;
% % %             sub(3).t_end(w) = t_end;
% % %             sub(3).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% end
% 
% %SUBJECT 2 NEW
% nom_sujet = 'SUBJECT 2 NEW';
% disp(['sujet : ', nom_sujet]);
% File = {'marija_unvoluntary'};
% rep = 'subject2_new\unvoluntary';
% [class_mvt, class_ep_4, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
%     disp(['intervalle de temps : ' num2str([t_start t_end])])
%     approche_classif_2(class_ep_4,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(4).chan{w} = chan;
% % %             sub(4).t_start(w) = t_start;
% % %             sub(4).t_end(w) = t_end;
% % %             sub(4).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% end
% 
% %SUBJECT 3
% nom_sujet = 'SUBJECT 3';
% disp(['sujet : ', nom_sujet]);
% File = {'mustafa_unvoluntary'};
% rep = 'subject3\unvoluntary';
% [class_mvt, class_ep_5, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
%     disp(['intervalle de temps : ' num2str([t_start t_end])])
%     approche_classif_2(class_ep_5,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(5).chan{w} = chan;
% % %             sub(5).t_start(w) = t_start;
% % %             sub(5).t_end(w) = t_end;
% % %             sub(5).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% end
% 
% %Duy
% % nom_sujet = 'DUY';
% % disp(['sujet : ', nom_sujet]);
% % File = {'unvoluntary'};
% % rep = 'duy';
% % [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
% %     disp(['intervalle de temps : ' num2str([t_start t_end])])
% %     approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(6).chan{w} = chan;
% % %             sub(6).t_start(w) = t_start;
% % %             sub(6).t_end(w) = t_end;
% % %             sub(6).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% % end
% % 
% %Jakob
% % nom_sujet = 'JAKOB';
% % disp(['sujet : ', nom_sujet]);
% % File = {'unvoluntary'};
% % rep = 'jakob';
% % [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% % for i=1:length(possible_chan)
% %     disp(['chan : ' num2str(possible_chan{i})])
% %     chan = possible_chan{i};
% %     disp(['intervalle de temps : ' num2str([t_start t_end])])
% %     approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(7).chan{w} = chan;
% % %             sub(7).t_start(w) = t_start;
% % %             sub(7).t_end(w) = t_end;
% % %             sub(7).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% % end
% % 
% %Dantan
% nom_sujet = 'DANTAN';
% disp(['sujet : ', nom_sujet]);
% File = {'unvoluntary'};
% rep = 'DanTan';
% [class_mvt, class_ep_6, chan_label, Fe] = read_imran_file(File,rep);
% % w=0;
% for i=1:length(possible_chan)
%     disp(['chan : ' num2str(possible_chan{i})])
%     chan = possible_chan{i};
%     disp(['intervalle de temps : ' num2str([t_start t_end])])
%     approche_classif_2(class_ep_6,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %     for t_start = 0:0.05:0.5
% % %         for t_end = t_start+0.3:0.05:1
% % %             w = w+1;
% % %             disp(['intervalle de temps : ' num2str([t_start t_end])])
% % %             sub(8).chan{w} = chan;
% % %             sub(8).t_start(w) = t_start;
% % %             sub(8).t_end(w) = t_end;
% % %             sub(8).crit_cal(w)= approche_classif_2(class_ep,chan,C,filtering,Fe,'rakoto',t_start,t_end);
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'lda');
% % %         % approche_classif_2(class_ep,chan,C,filtering,Fe,'bayes');
% % %         end
% % %     end
% end

