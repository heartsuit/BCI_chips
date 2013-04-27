close all
clear all
clc
x = rand(12,16,3);
y = 1+(rand(12,1)>0.5);


[nb_trials,nb_samples,nb_chans] = size(x);
N_dec = floor(log2(nb_trials));
N_dec = 3;
wname = 'db2';
K = 0.5;
set = set_wp_marg(x,wname,N_dec);
fisher_tree = calc_fisher_wp_set (set,y,K);
basis = find_best_basis_fisher(fisher_tree,N_dec);
