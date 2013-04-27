function features = calc_features(signal,N_dec,ind_des,h_ini,g_ini)

[dwt_coeff_ini, l_ini] = wavedec_simplif(signal,'per',N_dec,h_ini,g_ini,0,0);
features = calc_marg(dwt_coeff_ini,l_ini,ind_des);
% features = calc_dwt_coeff(dwt_coeff_ini,l_ini,ind_des);
