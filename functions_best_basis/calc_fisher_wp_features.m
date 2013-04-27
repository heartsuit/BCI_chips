function fisher_tree = calc_fisher_wp_features(features,y,K,restric,N_dec)

%set : nb_trials x nb_nodes x nb_chans

[nb_trials,nb_nodes,nb_chans] = size(features);

fisher_tree = zeros(nb_nodes,1);
nb_classes = max(y);
MA_c = zeros(nb_classes,nb_chans);
dist_c = zeros(nb_classes,1);
% intra_dist = zeros(nb_trials,1);

if restric
    gd_nodes=[];
    for i=N_dec-2:N_dec
        gd_nodes=[gd_nodes [2^i:2^i+2^(i-(N_dec-3))-1]];
    end
    for i=gd_nodes
        %Ecart entre les centres
        MA =  squeeze(mean(features(:,i,:),1))';
        for j=1:nb_classes
            MA_c(j,:) = squeeze(mean(features(y==j,i,:),1));
            dist_c(j) = sum(y==j)*sum((MA_c(j,:)-MA).^2);
        end
        distance = sum(dist_c)/nb_trials;
        %Inertie
    %     for j=1:nb_trials
        intra_dist = sum(sum((squeeze(features(:,i,:))-MA_c(y,:)).^2));
    %     end
        inertia = intra_dist/nb_trials;

        %------------------------------------------------------------
        % formule de Fisher pour d�compo en arbres


        fisher_tree(i) = K*distance - (1-K)*inertia;
        %---------------------------------------------------------
    end
else
    for i=1:nb_nodes
        %Ecart entre les centres
        MA =  squeeze(mean(features(:,i,:),1))';
        for j=1:nb_classes
            MA_c(j,:) = squeeze(mean(features(y==j,i,:),1));
            dist_c(j) = sum(y==j)*sum((MA_c(j,:)-MA).^2);
        end
        distance = sum(dist_c)/nb_trials;
        %Inertie
    %     for j=1:nb_trials
        intra_dist = sum(sum((squeeze(features(:,i,:))-MA_c(y,:)).^2));
    %     end
        inertia = intra_dist/nb_trials;

        %------------------------------------------------------------
        % formule de Fisher pour d�compo en arbres


        fisher_tree(i) = K*distance - (1-K)*inertia;
        %---------------------------------------------------------
    end

end
%%%%%%%%%%%%%%%%%%%%%%%% TRIGERING OF VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fisher_tree(abs(fisher_tree)<10^(-4))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end