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
    W1 = (2*0.5)/srate;   %Fréquence de coupure basse
    W2 = (2*35)/srate;  %Fréquence de coupure haute
    
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
figure
hold on;
plot(class(1).signal_init_mean)
plot(class(1).signal_init_mean+2*sqrt(class(1).signal_var),':')
plot(class(1).signal_init_mean-2*sqrt(class(1).signal_var),':')
plot(class(2).signal_init_mean,'r')
plot(class(2).signal_init_mean+2*sqrt(class(2).signal_var),'r:')
plot(class(2).signal_init_mean-2*sqrt(class(2).signal_var),'r:')
hold off

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

figure
hold on;
plot(class(1).signal_mean)
plot(class(1).signal_mean+2*sqrt(class(1).signal_SpatFilter_var),':')
plot(class(1).signal_mean-2*sqrt(class(1).signal_SpatFilter_var),':')
plot(class(2).signal_mean,'r')
plot(class(2).signal_mean+2*sqrt(class(2).signal_SpatFilter_var),'r:')
plot(class(2).signal_mean-2*sqrt(class(2).signal_SpatFilter_var),'r:')
hold off

%% Loo
for i=1:nbclasses
    nb_trials = size(class(i).signal,1);
    if i~=2
        for j=1:nb_trials
            class(i).signal_corr(j,:) = xcorr(class(i).signal_SpatFilter(j,:),class(2).signal_mean);
        end
    else
        for j=1:nb_trials
            ind_moy = [1:j-1 j+1:nb_trials];
            class(i).signal_corr(j,:) = xcorr(class(i).signal_SpatFilter(j,:),squeeze(mean(class(i).signal_SpatFilter(ind_moy,:),1)));
        end
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
detect_haut = zeros(2,nb_taux);
taux_h=taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max;
for taux=taux_h
    i=i+1;
    detect_haut(1,i) = 100*sum(max(class(1).signal_corr_mean(:,ori),[],2) >= taux)/size(class(1).signal_corr_mean,1);
    detect_haut(2,i) = 100*sum(max(class(2).signal_corr_mean(:,ori),[],2) >= taux)/size(class(2).signal_corr_mean,1);
end
i=0;
detect_bas = zeros(2,nb_taux);
taux_b = taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max;
for taux=taux_b
    i=i+1;
    detect_bas(2,i) = 100*sum(max(class(1).signal_corr_mean(:,ori),[],2) <= taux)/size(class(1).signal_corr_mean,1);
    detect_bas(1,i) = 100*sum(max(class(2).signal_corr_mean(:,ori),[],2) <= taux)/size(class(2).signal_corr_mean,1);
end


figure
hold on
plot(detect_haut(1,:), detect_haut(2,:), 'x')
plot(detect_bas(1,:), detect_bas(2,:), 'rx')
grid
hold off
xlabel('Taux de fausse alarme')
ylabel('Taux de bonnes détections')
legend('classe fausse', 'classe juste')
% title(['Courbe Roc pour l''individu 1 new']);


ind_haut = find(detect_haut(1,:) < 20);
[max_value,ind_max] = max(detect_haut(2,ind_haut));
seuil_h = taux_h(ind_haut(ind_max));
disp(['seuil haut : ' num2str(seuil_h)]);

ind_bas = find(detect_bas(1,:) < 20);
[max_value,ind_max] = max(detect_bas(2,ind_bas));
seuil_b = taux_b(ind_bas(ind_max));
disp(['seuil bas : ' num2str(seuil_b)]);


nb_rep_correct = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h);
nb_rep_wrong = sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h);

P_C_C = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b)/nb_rep_correct;
P_W_W = sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h)/nb_rep_wrong;
P_W_C = sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h)/nb_rep_correct;
P_C_W = sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b)/nb_rep_wrong;
fprintf('P_E_E\t\t\t | E=correct | E=wrong\n')
fprintf('E_estim=correct\t | %f\t | %f\n',P_C_C,P_C_W)
fprintf('E_estim=wrong\t | %f\t | %f\n',P_W_C,P_W_W)
% %Nombre de signaux au dessus du seuil_h
% nb_seuil_h = sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h) + sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h);
% %Nombre de signaux en dessous du seuil_b
% nb_seuil_b = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b);
% % disp(['nb signaux au dessus du seuil_h : ' num2str(nb_seuil_h)]);
% % disp(['nb signaux en dessous du seuil_b : ' num2str(nb_seuil_b)]);
% %Pourcentage de détection :
p_detect_global = (nb_rep_correct+nb_rep_wrong)/(size(class(1).signal_corr,1)+size(class(2).signal_corr,1));
%Pourcentage réponse classe wrong
p_rep_h = nb_rep_wrong/size(class(2).signal_corr,1);
%Pourcentage réponse classe correct
p_rep_b = nb_rep_correct/size(class(1).signal_corr,1);
disp(['pourcentage de réponse global : ' num2str(p_detect_global)])
disp(['pourcentage de réponse wrong : ' num2str(p_rep_h)])
disp(['pourcentage de réponse correct : ' num2str(p_rep_b)])

% %Pourcentage de signaux avec pic réellement détectés
% p_detect_h = sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h)/nb_seuil_h;
% disp(['wrong/wrong : ' num2str(p_detect_h)]);
% %Porcentage de signaux sans pic réeelement détectés
% p_detect_b = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b)/nb_seuil_b;
% disp(['correct/correct : ' num2str(p_detect_b)]);
