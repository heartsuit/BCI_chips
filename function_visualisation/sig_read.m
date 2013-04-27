function [x,y] = sig_read(filename,pathname,type)

if nargin < 2
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

trialNb=length(trialsInformation);
x = zeros(trialNB,size(daqread([fullfile(pathname,filename(1:end-4)) '_'  num2str(1) '.daq'])));
y = zeros(trialNB,1);
if strcmp(type,'mvt')
    for i=1:trialNb
        x(i,:,:) = daqread( [fullfile(pathname,filename(1:end-4)) '_'  num2str(i) '.daq']);
        if trialsInformation{1,i}.classNb==1
            y(i) = 1;
        else
            y(i) = 2;
        end
    end
elseif strcmp(type,'ep')
    for i=1:trialNb
        x(i,:,:) = daqread( [fullfile(pathname,filename(1:end-4)) '_'  num2str(i) '_EP.daq']);
        if trialsInformation{1,i}.rand_result ~= -1
            y(i) = 1;
        else
            y(i) = 2;
        end
    end
end
x = x(y ~=0,:,:);
y = y(y~=0);

% chan_label = {'Fcz','Fz','Cpz','C2','C1','Cz'};
