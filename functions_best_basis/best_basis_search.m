function [basis, fisher_tree, dyad_basis] = best_basis_search(features,y,...
    N_dec,K,plot_tree,sim_sig,restric)

fisher_tree = calc_fisher_wp_features(features,y,K,restric,N_dec); %calcul les coeff C(A_node) pour tous les noeuds
basis = find_best_basis_fisher(fisher_tree,N_dec,plot_tree,sim_sig,restric); %fait le calcul additif de la remont�e des noeuds pour touver la meilleure base
basis = boolean(basis);
nb_nodes = length(fisher_tree);
dyad_basis = zeros(1,nb_nodes);
for i=1:N_dec
    dyad_basis(node(i,1)) = 1;
end
dyad_basis(node(N_dec,0)) = 1;

%%%%%%%%%%%%%%%%%%%% MARGINAL VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% level=4;
% rank=0;
% class_1=1;
% class_2=2;
% class_3=3;
% 
% 
% marg_chan_1=features(y==class_1,node(4,1),:);
% marg_mean_1=mean(marg_chan_1,3);
% 
% marg_chan_2=features(y==class_2,node(4,1),:);
% marg_mean_2=mean(marg_chan_2,3);
% 
% marg_chan_3=features(y==class_3,node(4,1),:);
% marg_mean_3=mean(marg_chan_3,3);
% 
% figure('Name',['Values of the marginal j=' num2str(level) ' p=' num2str(rank) ' among the learning set'])
% 
% plot(marg_mean_1,ones(1,length(marg_mean_1)),'s','MarkerSize',5,'color','b');
% hold on;
% plot(marg_mean_2,0.9*ones(1,length(marg_mean_2)),'s','MarkerSize',5,'color','r');
% plot(marg_mean_3,0.8*ones(1,length(marg_mean_3)),'s','MarkerSize',5,'color','g');
% ylim([0.4 1.4]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%