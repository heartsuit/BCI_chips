function fisher_tree = calc_fisher_wp_set (set,y,K)

%set : nb_trials x nb_nodes x nb_chans

[nb_trials,nb_nodes,nb_chans] = size(set);

fisher_tree = zeros(nb_nodes,1);
nb_classes = max(y);
MA_c = zeros(nb_classes,nb_chans);
dist_c = zeros(nb_classes,1);
intra_dist = zeros(nb_trials,1);
for i=1:nb_nodes
    %Ecart entre les centres
%     keyboard
    MA =  squeeze(mean(set(:,i,:),1))';
    for j=1:nb_classes
        MA_c(j,:) = mean(set(y==j,i,:),1);
        dist_c(j) = sum(y==j)*(MA_c(j,:)-MA)*(MA_c(j,:)-MA)';
    end
    distance = sum(dist_c)/nb_trials;
    %Inertie
    for j=1:nb_trials
            intra_dist(j) = (squeeze(set(j,i,:))'-MA_c(y(j),:))*(squeeze(set(j,i,:))'-MA_c(y(j),:))';
    end
    inertia = sum(intra_dist)/nb_trials;
    fisher_tree(i) = K*distance - (1-K)*inertia;
end