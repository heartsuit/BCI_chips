function set = set_wp_marg(x,wname,N_dec,h,g)

%x ensemble des exemples (nb_trials x nb_samples x nb_chans)
[nb_trials,nb_samples,nb_chans] = size(x);

if nargin<2
    wname = 'db2';
end


if nargin<3
    N_dec = floor(log2(nb_samples));
end

[h_ini, g_ini] = wfilters(wname);

nb_nodes = (2^(N_dec+1)-1);
set = zeros(nb_trials,nb_nodes,nb_chans);
for i=1:nb_trials
    for k=1:nb_chans
        if nargin<4
            wp_tree_i = wpdec_simplif(x(i,:,k),'per',N_dec,h_ini,g_ini,0,0); %création de l'arbre des coefficient [a,d]            
        else
            wp_tree_i = wpdec_simplif(x(i,:,k),'per',N_dec,h(k,:),g(k,:),0,0); %création de l'arbre des coefficient [a,d]
        end
        set(i,:,k) = wp_marg(wp_tree_i,'abs');  %conversion de [a,d] à M
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       set(i,:,k)(node(k,d))= M_phi_ind_k_exp_d ; (cf wavelet lesson)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end
