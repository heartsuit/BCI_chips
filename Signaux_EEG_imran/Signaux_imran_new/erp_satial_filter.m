function [f,class] = erp_satial_filter(class,chan)

f=0;

nb_classes = size(class,2);
if nargin < 2
    chan = 1:size(class(1).signal,3);
else
    for i=1:nb_classes
        class(i).signal = class(i).signal(:,:,chan);
    end
end
if (nb_classes ~= 2)
    disp('Must have 2 classes for CSP!')
end

p = [3/4, 1/4];
nb_chan = length(chan);

for i=1:nb_classes
    class(i).signal_mean = squeeze(mean(class(i).signal,1));
end

total_mean = (class(1).signal_mean + class(2).signal_mean)/2;
Sb = zeros(nb_chan,nb_chan);
Sw = zeros(nb_chan,nb_chan);
nb_trials_total = 0;
for i=1:nb_classes
    Sb = Sb + p(i)*(class(i).signal_mean-total_mean)'*(class(i).signal_mean-total_mean);
    nb_trials = size(class(i).signal,1);
    nb_trials_total = nb_trials_total + nb_trials;
    for j=1:nb_trials
        Sw = Sw + (squeeze(class(i).signal(j,:,:))-class(i).signal_mean)'*(squeeze(class(i).signal(j,:,:))-class(i).signal_mean);
    end
    
end
% Sb = Sb*Sb';
% Sw = Sw*Sw';
Sw = Sw/nb_trials_total+1e-3*eye(size(Sw));
L = chol(Sw);
% inv_L = inv(L);

[B,D] = eig(L'\Sb/L);
% [B,D] = eig(Sb,Sw);
[max_D,ind]=max(diag(D));
f=(L\B(:,ind))';
f=f/norm(f);
for i=1:nb_classes
    nb_trials = size(class(i).signal,1);
    for j = 1:nb_trials
        class(i).signal_SpatFilter(j,:) = f*squeeze(class(i).signal(j,:,:))';
    end
end