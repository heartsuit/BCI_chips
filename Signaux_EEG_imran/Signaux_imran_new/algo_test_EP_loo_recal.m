function algo_test_EP_loo_recal(class,chan,csp,recal,print_result,Fe,filtering,N_start,N_end)
global file_name
close all
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

if nargin == 5
    filtering = 0;
end
if nargin > 2 && nargin < 8
    N_start = 2200;
    N_end = 2800;
elseif nargin < 5
    error('not enough argument');
end
N_sample = 16;
if filtering == 1
    srate = Fe;
    N = 4;  %Ordre du filtre
%     R = 20;    %Ripple du filtre
    W1 = (2*1)/srate;   %Fr�quence de coupure basse (1Hz)
    W2 = (2*10)/srate;  %Fr�quence de coupure haute (10Hz)
    
    Wp = [W1 W2];
    [B,A] = butter(N,W2);
    % [B,A] = cheby2(N,R,Wp,'stop'); %Filtre de chebychev
%     [B,A] = cheby2(N,R,W1,'low');
    for i=chan
        for j=1:size(class(1).signal,1)
            class(1).signal(j,:,i) = filter(B,A,class(1).signal(j,:,i));
        end
        for j=1:size(class(2).signal,1)
            class(2).signal(j,:,i) = filter(B,A,class(2).signal(j,:,i));
        end
    end
end
nbclasses = size(class,2);
for i=1:nbclasses
    class(i).signal = class(i).signal(:,:,chan);
end
N_tronc = N_start:N_end;

% class_spa_filter(1).signal = class(1).signal(:,N_tronc,:);
% class_spa_filter(2).signal = class(2).signal(:,N_tronc,:);

if csp == 1
    spatial_filtering = 'CSP';
elseif csp == 2
    spatial_filtering = 'Fisher';
else
    spatial_filtering = 'Mean of channels';
end

for i=1:nbclasses
    if i==1
        ind_rec = recalage(class(2).signal,N_tronc);
        signal_app = zeros(size(class(2).signal,1),length(N_tronc),length(chan));
        if recal == 1
            for k=1:size(class(2).signal,1)
                signal_app(k,:,:) = class(2).signal(k,N_tronc-ind_rec(k),:);
            end
        else
            sig_recal = 'No';
            signal_app = class(2).signal(:,N_tronc,:);
        end
        nb_trials = size(class(i).signal,1);
        for j=1:nb_trials
            signal_test = squeeze(class(i).signal(j,N_tronc,:));
            ind_app = [1:j-1 j+1:nb_trials];
            class_spa_filter(1).signal = class(1).signal(ind_app,N_tronc,:);
            class_spa_filter(2).signal = signal_app;
            if csp == 1
                f = CSP(class_spa_filter);
            elseif csp == 2
                f = erp_satial_filter(class_spa_filter);
            else
                f = ones(1,length(chan))/length(chan);
            end
            class(i).signal_corr(j,:) = xcorr(f*signal_test',f*squeeze(mean(signal_app,1))');
        end       
    else
        nb_trials = size(class(i).signal,1);
        for j=1:nb_trials
            signal_test = squeeze(class(i).signal(j,N_tronc,:));
            ind_app = [1:j-1 j+1:nb_trials];
            signal_app_inter = class(2).signal(ind_app,:,:);
            if recal == 1
                sig_recal = 'Yes';
                ind_rec = recalage(signal_app_inter,N_tronc);
                signal_app = zeros(length(ind_app),length(N_tronc),length(chan));
                for k=1:length(ind_app)
                    signal_app(k,:,:) = signal_app_inter(k,N_tronc-ind_rec(k),:);
                end
            else
                sig_recal = 'No';
                signal_app = signal_app_inter(:,N_tronc,:);
            end
            class_spa_filter(1).signal = class(1).signal(:,N_tronc,:);
            class_spa_filter(2).signal = signal_app;
            if csp == 1
                f = CSP(class_spa_filter);
            elseif csp == 2
                f = erp_satial_filter(class_spa_filter);
            else
                f = ones(1,length(chan))/length(chan);
            end
%             disp(num2str(f))
            class(i).signal_corr(j,:) = xcorr(f*signal_test',f*squeeze(mean(signal_app,1))');
        end 
    end    
end
% 
% for i=chan
%     figure
%     hold on
%     plot(mean(class(1).signal(:,:,i),1)')
%     plot(mean(class(2).signal(:,:,i),1)','r')
%     hold off
% end


% 
% figure
% hold on
% for i=1:size(class(1).signal_corr,1)
%     plot(class(1).signal_corr(i,:))
% end
% for i=1:size(class(2).signal_corr,1)
%     plot(class(2).signal_corr(i,:),'r')
% end
% hold off

figure
hold on
plot(mean(class(1).signal_corr,1))
plot(mean(class(2).signal_corr,1),'r')
hold off

taux_max = max(max(class(2).signal_corr));
taux_min = min(min(class(2).signal_corr));
ori = round(size(class(1).signal_corr,2)/2+1)+(-100:100);

nb_taux=100;
i=0;
detect_haut = zeros(2,nb_taux);
taux_h=taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max;


for taux = taux_h
    i=i+1;
    detect_haut(1,i) = 100*sum(max(class(1).signal_corr(:,ori),[],2) >= taux)/size(class(1).signal_corr,1);%(sum(max(class(1).signal_corr(:,ori),[],2) >= taux) + sum(max(class(2).signal_corr(:,ori),[],2) >= taux));
    detect_haut(2,i) = 100*sum(max(class(2).signal_corr(:,ori),[],2) >= taux)/size(class(2).signal_corr,1);%(sum(max(class(1).signal_corr(:,ori),[],2) >= taux) + sum(max(class(2).signal_corr(:,ori),[],2) >= taux));
end
% size(class(2).signal_corr,1);%

i=0;
detect_bas = zeros(2,nb_taux);
taux_b = taux_min:(taux_max-taux_min)/(nb_taux-1):taux_max;
for taux = taux_b
    i=i+1;
    detect_bas(2,i) = 100*sum(max(class(1).signal_corr(:,ori),[],2) <= taux)/size(class(1).signal_corr,1);%(sum(max(class(1).signal_corr(:,ori),[],2) <= taux) + sum(max(class(2).signal_corr(:,ori),[],2) <= taux));
    detect_bas(1,i) = 100*sum(max(class(2).signal_corr(:,ori),[],2) <= taux)/size(class(2).signal_corr,1);%(sum(max(class(1).signal_corr(:,ori),[],2) <= taux) + sum(max(class(2).signal_corr(:,ori),[],2) <= taux));
end
% size(class(1).signal_corr,1);%

figure
hold on
plot(detect_haut(1,:), detect_haut(2,:), 'x')
plot(detect_bas(1,:), detect_bas(2,:), 'rx')
grid
hold off
xlabel('Taux de fausse alarme')
ylabel('Taux de bonnes d�tections')
legend('classe fausse', 'classe juste')

P_sig1_seuil_h = 20;
ind_haut = find(detect_haut(1,:) < P_sig1_seuil_h);
[max_value,ind_max] = max(detect_haut(2,ind_haut));
seuil_h = taux_h(ind_haut(ind_max));
disp(['seuil haut : ' num2str(seuil_h)]);

P_sig2_seuil_b = 20;
ind_bas = find(detect_bas(1,:) < P_sig2_seuil_b);
[max_value,ind_max] = max(detect_bas(2,ind_bas));
seuil_b = taux_b(ind_bas(ind_max));
disp(['seuil bas : ' num2str(seuil_b)]);

nb_detect_1 = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h);
nb_detect_2 = sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h);


P_W_W = 100*sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h)/nb_detect_2;
P_C_W = 100*sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b)/nb_detect_2;
P_C_C = 100*sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b)/nb_detect_1;
P_W_C = 100*sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h)/nb_detect_1;

p_rep_global = 100*(nb_detect_1+nb_detect_2)/(size(class(1).signal_corr,1)+size(class(2).signal_corr,1));
p_rep_wrong = 100*nb_detect_2/size(class(2).signal_corr,1);
p_rep_correct = 100*nb_detect_1/size(class(1).signal_corr,1);


fprintf('P_E_estim_E\t | E_estim = Wrong\t | E_estim = Correct\n')
fprintf('E = Wrong\t | %f\t\t\t | %f\n', P_W_W, P_C_W)
fprintf('E = Correct\t | %f\t\t\t | %f\n', P_W_C, P_C_C)

disp(['pourcentage de r�ponse global : ' num2str(p_rep_global)])
disp(['pourcentage de r�ponse correct : ' num2str(p_rep_correct)])
disp(['pourcentage de r�ponse wrong : ' num2str(p_rep_wrong)])

c_q = (P_C_C + P_W_W + p_rep_global)/3;
disp(['Crit�re qualit� : ' num2str(c_q)]);
if print_result
    fid = fopen(file_name,'a');
    fprintf(fid,'Spatial filter : %s \\\\\nTemporal readjustment : %s \\\\\n',spatial_filtering,sig_recal);
    fprintf(fid,'Percentage of signals correct detected as wrong : %6.1f \\\\\n',P_sig1_seuil_h);
    fprintf(fid,'Percentage of signals wrong detected as correct : %6.1f \\\\\n',P_sig2_seuil_b);
    fprintf(fid,'Table of results : \\\\\n');
    fprintf(fid,'\\begin{tabular}{|c|c|c|}\n\\hline\t\t\t\t& $\\hat{e}$ = juste & $\\hat{e}$ = faux \\\\\n');
    fprintf(fid,'\\hline  $e$ = juste\t& %6.0f\\%%\t\t\t& %6.0f\\%%\t\t\\\\\n', P_C_C,P_W_C);
    fprintf(fid,'\\hline  $e$ = faux\t& %6.0f\\%%\t\t\t& %6.0f\\%%\t\t\\\\\n', P_C_W,P_W_W);
    fprintf(fid,'\\hline\n\\end{tabular}\\\\\n');
    fprintf(fid,'Global percentage of responses : %6.1f \\\\\n',p_rep_global);
    fprintf(fid,'Percentage of correct responses : %6.1f \\\\\n',p_rep_correct);
    fprintf(fid,'Percentage of wrong responses : %6.1f \\\\\n',p_rep_wrong);
    fprintf(fid,'Quality criterion : %6.1f \\\\\n\n',c_q);
    fclose(fid);
end
% keyboard
% %Nombre de signaux au dessus du seuil_h
% nb_seuil_h = sum(max(class(1).signal_corr(:,ori),[],2) >= seuil_h) + sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h);
% %Nombre de signaux en dessous du seuil_b
% nb_seuil_b = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b) + sum(max(class(2).signal_corr(:,ori),[],2) <= seuil_b);
% disp(['nb signaux au dessus du seuil_h : ' num2str(nb_seuil_h)]);
% disp(['nb signaux en dessous du seuil_b : ' num2str(nb_seuil_b)]);
% %Pourcentage de d�tection :
% p_detect = (nb_seuil_h+nb_seuil_b)/(size(class(1).signal_corr,1)+size(class(2).signal_corr,1));
% disp(['pourcentage de d�tection : ' num2str(p_detect)])
% 
% %Pourcentage de signaux avec pic r�ellement d�tect�s
% p_detect_h = sum(max(class(2).signal_corr(:,ori),[],2) >= seuil_h)/nb_seuil_h;
% disp(['wrong/wrong : ' num2str(p_detect_h)]);
% %Porcentage de signaux sans pic r�eelement d�tect�s
% p_detect_b = sum(max(class(1).signal_corr(:,ori),[],2) <= seuil_b)/nb_seuil_b;
% disp(['correct/correct : ' num2str(p_detect_b)]);
