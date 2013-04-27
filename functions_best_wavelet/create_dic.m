function [features_dic,h_dic,g_dic] = create_dic(signal,N_dec,ind_des)

pas = 0.2; % discretisation step
a = -pi+pas:pas:pi;
nb_wave = length(a);
lh_filtre = 4;
[nb_trials,nb_samples,nb_chans] = size(signal);
h_dic = zeros(nb_wave,lh_filtre,nb_chans);
g_dic = zeros(nb_wave,lh_filtre,nb_chans);
i_wave = 0;
for ii=1:nb_wave
    hh = construction_filtreh(lh_filtre,a(ii));
    gh = hh(end:-1:1);
    for jj=1:length(gh)
        gh(jj) = (-1)^jj * gh(jj);
    end
    nbz = nbturn(hh); % turncount of mother wavelet (hh)
    if nbz < 30
        i_wave = i_wave+1;
        % Computing of the marginales and the criterion
        %on met a jour la marginale k avec les nouveaux paramètres de
        %l'ondelette du chan k.
        h_dic(i_wave,:,:) = (ones(nb_chans,1)*hh)';
        g_dic(i_wave,:,:) = (ones(nb_chans,1)*gh)';
    end
end

nb_wave = i_wave;
h_dic = h_dic(1:nb_wave,:,:);
g_dic = g_dic(1:nb_wave,:,:);
% dwt_coeff = zeros(nb_trials_all,nb_samples,nb_chans,nb_wave);
% l = zeros(nb_trials_all,N_dec+2,nb_chans,nb_wave);
% features_dic = zeros(nb_trials,length(ind_des),nb_chans,nb_wave);
% features_dic = zeros(nb_trials,nb_samples,nb_chans,nb_wave);
for ii=1:nb_wave
    for i=1:nb_trials
        for k=1:nb_chans
            features_dic(i,:,k,ii) = calc_features(signal(i,:,k),N_dec,ind_des,h_dic(ii,:,k),g_dic(ii,:,k));
        end
    end
end