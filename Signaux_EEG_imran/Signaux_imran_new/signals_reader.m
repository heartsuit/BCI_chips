function [class] = signals_reader(File,Rep)

% 1st class Flex_fast 2nd class Ext_fast
% File = {'gediminas_unvoluntary'}; %real movement training set
% Rep = 
%File = {'wahid_2704_real_NEW_'};% imaginary training set

% [b,a] = butter(4,5/(1024/2),'low');

path(path,Rep);
x=1;
y=1;
F=size(File,2);
for f=1:F
    currentFile=File{f};
    load(currentFile);
    trialNb=length(trialsInformation);
    for i=1:trialNb
        if trialsInformation{1,i}.classNb==1
            class(1).signal{x}=daqread( [File{f} '_'  num2str(i) '.daq']);
            class(1).EP{x}=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
            class(1).true{x} = trialsInformation{1,i}.rand_result ~= -1;
            %             fast_signals_tr{x}=fast_signals_tr{x}(2048:7168,:);
            
            
            x=x+1;
        else
            class(2).signal{y}=daqread( [File{f} '_'  num2str(i) '.daq']);
            class(2).EP{y}=daqread( [File{f} '_'  num2str(i) '_EP.daq']);
            class(2).true{x} = trialsInformation{1,i}.rand_result ~= -1;
            %             slow_signals_tr{y}=slow_signals_tr{y}(2048:7168,:);
            y=y+1;
        end
    end
    
    %
end
