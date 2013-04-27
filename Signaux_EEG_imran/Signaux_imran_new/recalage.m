function ind_dec = recalage(signal,N_tronc)


signal_tronc_filter = squeeze(mean(signal(:,N_tronc,:),3));
signal_mean = squeeze(mean(signal_tronc_filter,1));
nb_trials = size(signal,1);
ind_dec = ones(length(nb_trials),1)/length(nb_trials);
signal_corr = zeros(nb_trials,2*length(N_tronc)-1);
for i=1:nb_trials
    signal_corr(i,:) = xcorr(signal_tronc_filter(i,:),signal_mean);
    ori = (-50:50)+round(size(signal_corr,2)/2);
    [max_corr, ind_dec(i)] = max(signal_corr(i,ori));
    ind_dec(i) = ori(ind_dec(i))-length(N_tronc)-1;
end