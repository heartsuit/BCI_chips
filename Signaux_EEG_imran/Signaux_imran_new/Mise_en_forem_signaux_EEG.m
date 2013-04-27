% makeEEG2_5voies
%
% 4-decimated samples (500 Hz --> 250 Hz)
% Creating the data :
% array (signal length, nb signals per class, nb channels, nb classes)

nbchannel = 1;
chan = 5;
individu = 6;

% load sf500_subD.mat
% bh = size(EEGData.BH);
% bl = size(EEGData.BL);
% mh = size(EEGData.MH);
% ml = size(EEGData.ML);

load imaginary.mat
size_ball = size(imaginary(individu).ballistic);
size_slow = size(imaginary(individu).slow);

% % If the signal length is not the same for each signal,
% % the number of points is ajusted to the minimum
% EEGData.BH = double(EEGData.BH(:,1:min([bh(2) bl(2) mh(2) ml(2)]),:));
% EEGData.BL = double(EEGData.BL(:,1:min([bh(2) bl(2) mh(2) ml(2)]),:));
% EEGData.MH = double(EEGData.MH(:,1:min([bh(2) bl(2) mh(2) ml(2)]),:));
% EEGData.ML = double(EEGData.ML(:,1:min([bh(2) bl(2) mh(2) ml(2)]),:));

% on ne prends que les 1024 premiers points des signaux après 3s
x=512;
imaginary(individu).ballistic = double(imaginary(individu).ballistic(:,x+1:x+1024,:));
imaginary(individu).slow = double(imaginary(individu).slow(:,x+1:x+1024,:));


npt = 256; % signal length


nbclass = 2;
nbsig = min(size_slow(3),size_ball(3)); % number of signals per class
xy = zeros(npt,nbsig,nbchannel,nbclass);
for i=1:nbsig,
        vect = decimate(imaginary(individu).ballistic(chan,:,i),4);
        xy(1:end,i,1,1) = vect;
        vect = decimate(imaginary(individu).slow(chan,:,i),4);
        xy(1:end,i,1,2) = vect;
end;
size(xy)
%save eeg2subD_bhbl_5voies xy
save imaginary_indiv_6 xy