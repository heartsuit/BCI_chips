function algo_test_EP(class,chan,csp,Fe,filtering,N_start,N_end)
%=====================================================
%|algo_test_EP(class,chan,Fe,filtering,N_start,N_end)|
%=====================================================
%Strucure class : size >= 2
%field : class(i).signal
%chan : channels used for the detection : chan = [1 3 6] => the sum of chan
%1 3 and 6 is used for the intercorrelation.
%Fe : sampling frequency
%Filtering : 1 => preprocessing of EEG signals by filtering | 0 => Nothing
%N_start : first sample of the EEG signal
%N_end : last sample of the EEG signal


if nargin == 3
    filtering = 0;
end
if nargin > 2 && nargin < 6
    N_start = 2200;
    N_end = 2800;
elseif nargin < 3
    error('not enough argument');
end

%Dans le cas non volontaire
% t_debut = 2.8;
% t_fin = 3.5;
% t_debut = 2.150;
% t_fin = 2.650;

%dans le cas volontaire on augmentera afin d'avoir les MRCP
nbclasses = size(class,2);
% N_tronc = round((t_debut:(1/srate):t_fin)*srate);
N_tronc = N_start:N_end;
%Filtrage des signaux
if filtering
    srate = Fe;
    N = 4;  %Ordre du filtre
    R = 20;    %Ripple du filtre
    W1 = (2*0.5)/srate;   %Fr�quence de coupure basse
    W2 = (2*35)/srate;  %Fr�quence de coupure haute
    
    Wp = [W1 W2];
    [B,A] = cheby2(N,R,Wp); %Filtre de chebychev
    for i=1:nbclasses
        class(i).signal = filter(B,A,class(i).signal,[],2);
    end
end

for i=1:nbclasses
    class(i).signal = class(i).signal(:,N_tronc,:);
    class(i).signal_init_mean = squeeze(mean(squeeze(sum(class(i).signal(:,:,chan),3))/length(chan),1));
    class(i).signal_var = squeeze(var(squeeze(sum(class(i).signal(:,:,chan),3))/length(chan),0,1));
end
% figure
% hold on;
% plot(class(1).signal_init_mean)
% plot(class(1).signal_init_mean+2*sqrt(class(1).signal_var),':')
% plot(class(1).signal_init_mean-2*sqrt(class(1).signal_var),':')
% plot(class(2).signal_init_mean,'r')
% plot(class(2).signal_init_mean+2*sqrt(class(2).signal_var),'r:')
% plot(class(2).signal_init_mean-2*sqrt(class(2).signal_var),'r:')
% hold off

if csp == 1
    [f,class] = CSP(class,chan);
elseif csp == 2
    [f,class] = erp_satial_filter(class,chan);
else
    for i=1:nbclasses
        class(i).signal_SpatFilter = squeeze(sum(class(i).signal(:,:,chan),3))/length(chan);
    end
end

for i=1:nbclasses
    class(i).signal_mean = squeeze(mean(class(i).signal_SpatFilter,1));
    class(i).signal_SpatFilter_var = squeeze(var(class(i).signal_SpatFilter,0,1));
end

% figure
% hold on;
% plot(class(1).signal_mean)
% plot(class(1).signal_mean+2*sqrt(class(1).signal_SpatFilter_var),':')
% plot(class(1).signal_mean-2*sqrt(class(1).signal_SpatFilter_var),':')
% plot(class(2).signal_mean,'r')
% plot(class(2).signal_mean+2*sqrt(class(2).signal_SpatFilter_var),'r:')
% plot(class(2).signal_mean-2*sqrt(class(2).signal_SpatFilter_var),'r:')
% hold off

for i=1:nbclasses
    for j=1:size(class(i).signal,1)
        class(i).signal_corr(j,:) = xcorr(class(i).signal_SpatFilter(j,:),class(2).signal_mean);
    end
    class(i).signal_corr_mean = squeeze(mean(class(i).signal_corr,3));
end

% for i=chan
%     if filtering
%         figure
%         hold on;
%         plot(class(1).signal_filter_corr(:,:,i)','b')
%         plot(class(2).signal_filter_corr(:,:,i)','r')
%         hold off;
%     else
%         figure
%         hold on;
%         plot(class(1).signal_corr(:,:,i)','b')
%         plot(class(2).signal_corr(:,:,i)','r')
%         hold off;
%     end
% end

% figure
% hold on;
% plot(class(1).signal_corr','b')
% plot(class(2).signal_corr','r')
% hold off;


taux_max = max(max(class(2).signal_corr_mean));
taux_min = min(min(class(2).signal_corr_mean));
ori = round(size(class(1).signal_corr,2)/2+1)+(-100:100);

nb_taux=100;
i=0;
detect = zeros(2,nb_taux);
for taux=taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max
    i=i+1;
    detect(1,i) = 100*sum(max(class(1).signal_corr_mean(:,ori),[],2) >= taux)/size(class(1).signal_corr,1);
    detect(2,i) = 100*sum(max(class(2).signal_corr_mean(:,ori),[],2) >= taux)/size(class(2).signal_corr,1);
end

detect2 = zeros(2,nb_taux);
for taux=taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max
    i=i+1;
    detect2(2,i) = 100*sum(max(class(1).signal_corr_mean(:,ori),[],2) <= taux)/size(class(1).signal_corr,1);
    detect2(1,i) = 100*sum(max(class(2).signal_corr_mean(:,ori),[],2) <= taux)/size(class(2).signal_corr,1);
end


figure
hold on
plot(detect(1,:), detect(2,:), 'x')
plot(detect2(1,:), detect2(2,:), 'rx')
grid
hold off
xlabel('Taux de fausse alarme')
ylabel('Taux de bonnes d�tections')
legend('classe fausse', 'classe juste')
% title(['Courbe Roc pour l''individu 1 new']);
