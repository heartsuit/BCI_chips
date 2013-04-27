path_wavelab = 'C:\Users\Xavier\Desktop\mon_algo_best_basis\WAVELAB850';
addpath(genpath(path_wavelab));

name_wav = 'db2';
if strcmp(name_wav,'db2')
    h_1 = MakeONFilter('Daubechies',4);
end


signal = rand(25,128);
[n,D] = dyadlength(signal);
% qmf = MakeONFilter('Coiflet',3);
wp = WPAnalysis(signal,D,h_1);
stree = CalcStatTree_modif(wp,'Marginals');
[btree,vtree] = BestBasis(stree,D);
% 
% h_2 = wfilters(name_wav);
% disp(h_1)
% disp(h_2)

rmpath(genpath(path_wavelab));