clear variables
clc
% 1st class Flex_fast 2nd class Ext_fast
File = {'unvoluntary_'}; %real movement training set
%File = {'wahid_2704_real_NEW_'};% imaginary training set

[b,a] = butter(4,5/(1024/2),'low');

x=1;
y=1;
F=size(File,2);
for f=1:F
    currentFile=File{f};
    load(currentFile(1:end-1));
    trialNb=length(trialsInformation);
    for i=1:trialNb
        if trialsInformation{1,i}.classNb==1
            fast_signals_tr{x}=daqread( [File{f}  num2str(i) '.daq']);
            fast_signals_ep{x}=daqread( [File{f}  num2str(i) '_EP.daq']);
           
            x=x+1;
        else
            slow_signals_tr{y}=daqread( [File{f}  num2str(i) '.daq']);
            slow_signals_ep{y}=daqread( [File{f}  num2str(i) '_EP.daq']);
            %                     slow_signals{y}=slow_signals{y}(1:14541,:);
            y=y+1;
        end
    end
    
    %
end
