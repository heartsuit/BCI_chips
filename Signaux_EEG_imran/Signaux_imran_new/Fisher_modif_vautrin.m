% Classification of EEG signals with a wavelet packet decomposition
% The basis is determined 
% with the strategy of Wickerhauser and Coifman
% and a Fisher type discriminant measure
%
% The deepest decomposition level is limited to 'levelMax'



clc; clear all; warning off;

name_sig = 'imaginary_indiv_6';

name_wav = 'db2';   % Decomposition wavelet
K = 0.2;            % Parameter for the Fisher type measure
C = 100;            % SVM regularization parameter

eval(['load ' name_sig]);

disp(['Signal : ' name_sig]);
disp(['Wavelet : ' name_wav]);

N = size(xy,1);
nbsig = size(xy,2);
nbchann = size(xy,3);

if name_wav =='db2', h = MakeONFilter('Daubechies',4); end;
if name_wav =='db3', h = MakeONFilter('Daubechies',6); end;
if name_wav =='db4', h = MakeONFilter('Daubechies',8); end;

% Deepest decomposition level
deb = floor(log2(length(h))); J = floor(log2(N))-deb;
levelMax = J;

nbsubset = nbsig;   % Crossvalidation parameter (Leave One Out)
nbsigclas_test = floor(nbsig/nbsubset);
effsize = nbsubset*nbsigclas_test;

mcrate_test_cat = zeros(1,nbsubset);

% Calculation of the marginals
marg = fct_calc_marginals(xy,J,h);

for i = 1:nbsubset,
    disp(['Subset ' num2str(i) '/' num2str(nbsubset)]);

    % ind1 : start and ind2 : end index of current test set
    ind1 = 1+(i-1)*nbsigclas_test;
    ind2 = i*nbsigclas_test;
    currentInd = [1:(ind1-1), (ind2+1):effsize];

    %%% Training %%%

    trainingSet = marg(:,:,currentInd,:,:);
    crit_add = fct_calc_FisherTypeMeasure(trainingSet,K);
    mB = fct_search_bestBasis_max(crit_add,levelMax);
    length(find(mB==1))
    vect_trainingSet = fct_feature_extraction(trainingSet,mB);
    net = apprenti_SVM_OVRs(vect_trainingSet, C);

    %%% Test %%%

    testSet = marg(:,:,ind1:ind2,:,:);
    vect_testSet = fct_feature_extraction(testSet,mB);
    mc_test_cat = classer_ovr(vect_testSet,net);
    mcrate_test_cat(i) = mc_test_cat;

end

disp(' ')
disp(['   PCE Test (mean) : ' num2str(mean(mcrate_test_cat))]);


