function [class_mvt, class_ep] = load_class(filename,pathname)

if nargin ~= 2
    [filename, pathname] = uigetfile('*.mat', 'Pick an M-file');
end
if filename == 0
    return;
end
if ~strcmp(filename(end-3:end),'.mat')
    disp('choose a .mat')
    return
end
load(fullfile(pathname,filename));
% Fe = trialsInformation{1}.sampfreq;

w=1;
x=1;
y=1;
z=1;

trialNb=length(trialsInformation);
for i=1:trialNb
    if trialsInformation{1,i}.classNb==1
        class_mvt(1).signal(x,:,:)=daqread( [fullfile(pathname,filename(1:end-4)) '_'  num2str(i) '.daq']);
        x=x+1;
    else
        class_mvt(2).signal(y,:,:)=daqread( [fullfile(pathname,filename(1:end-4)) '_' num2str(i) '.daq']);
        y=y+1;
    end
    if trialsInformation{1,i}.rand_result ~= -1
        class_ep(1).signal(w,:,:)=daqread( [fullfile(pathname,filename(1:end-4)) '_'  num2str(i) '_EP.daq']);
        w=w+1;
    else
        class_ep(2).signal(z,:,:)=daqread( [fullfile(pathname,filename(1:end-4)) '_'  num2str(i) '_EP.daq']);
        z=z+1;
    end
end

% chan_label = {'Fcz','Fz','Cpz','C2','C1','Cz'};
