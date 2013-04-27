function [class_mvt, class_ep, chan_label, Fe] = read_imran_file(File,rep)

% 1st class Flex_fast 2nd class Ext_fast
% File = {'gediminas_unvoluntary'}; %real movement training set
%File = {'wahid_2704_real_NEW_'};% imaginary training set

% [b,a] = butter(4,5/(1024/2),'low');

addpath(path,rep);
w=1;
x=1;
y=1;
z=1;
F=size(File,2);
for f=1:F
    load(File{f});
    Fe = 1024;
    trialNb=length(trialsInformation);
    for i=1:trialNb
        if trialsInformation{1,i}.classNb==1
%             fast_signals_tr{x}=daqread( [File{f} '_'  num2str(i) '.daq']);
%             fast_signals_ep{x}=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
%             fast_signals_tr{x}=fast_signals_tr{x}(2048:7168,:);
            
            class_mvt(1).signal(x,:,:)=daqread( [File{f} '_'  num2str(i) '.daq']);
            x=x+1;
        else
%             slow_signals_tr{y}=daqread( [File{f} '_' num2str(i) '.daq']);
%             slow_signals_ep{y}=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
%             slow_signals_tr{y}=slow_signals_tr{y}(2048:7168,:);
            
            class_mvt(2).signal(y,:,:)=daqread( [File{f} '_' num2str(i) '.daq']);
            y=y+1;
        end
        if trialsInformation{1,i}.rand_result ~= -1
%             flag_time = find(strcmp(trialsInformation{1,i}.time_frames_label,'ep_focus')...
%                 | strcmp(trialsInformation{1,1}.time_frames_label,'result_display'));
%             N_tronc = (trialsInformation{1,i}.time_frames{flag_time} - trialsInformation{1,i}.time_frames{flag_time(1)} + 1) * Fe
%             size(daqread( [File{f} '_'  num2str(i) '_EP.daq']))
            class_ep(1).signal(w,:,:)=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
            w=w+1;
        else
            class_ep(2).signal(z,:,:)=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
            z=z+1;
        end
    end
end
chan_label = {'Fcz','Fz','Cpz','C2','C1','Cz'};
rmpath(rep);
