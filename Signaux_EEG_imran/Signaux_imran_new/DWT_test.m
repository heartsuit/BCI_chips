% main_cat_multimonocvEEG

% Multiclass Classification : SVM (linear) + OVR
% Multichannel
% Catalogue wavelet
% Crossvalidation to evaluate the performances (nbsubset)

clear all %clc;
warning off;

% loading the population EEG
name_sig = 'imaginary_indiv_3';
%name_sig = 'databhmh2_F402T6';
eval(['load ' name_sig]);
x=xy;y=xy;
N = size(xy,1);
nbvoies = size(xy,3);

% ------------------------------ Enter the specific data ----------------------------
% Crossvalidation parameter 
nbsubset = 3;

% representation parameters
name_wav = 'db4'; % if we use a standard wavelet

% svm parameters
C = 100; % SVM regularization parameter
% -------------------------------------------------------------------------
% % ondelette ˆ tester
% h_catK = [0.0602    0.7585    0.6469   -0.0513]
% % toolbox wavelet
% h_catK = wfilters(name_wav)
if name_wav =='db2' h_catK = [-0.1294    0.2241    0.8365    0.4830]; end;
if name_wav =='db3' h_catK = [0.0352   -0.0854   -0.1350    0.4599    0.8069    0.3327]; end;
if name_wav =='db4' h_catK = [-0.0106    0.0329    0.0308   -0.1870   -0.0280    0.6309    0.7148    0.2304
]; end;
lh = length(h_catK);
deb = floor(-log2(1/lh)); % number of deep levels not used for the calculation of dwt
nbfeatures = floor(log2(N))-deb; % number of features for each channel
nbsig = size(xy,2); if nbsubset >= (nbsig/2) disp('leave one out'); end;
nbsigclas_test = floor (nbsig/nbsubset);
% the training set is not exactly the whole data minus the one test set,
% but it is the union of the other test sets. So only effsize examples
% of the data set will ever be used
effsize = nbsubset*nbsigclas_test;
nbsigclas_app = effsize - nbsigclas_test;
if (nbsubset == 1),
  nbsigclas_app = size(x,2);
  nbsigclas_test = size(y,2);
end;
% ----------------------------------------------------------------------------------------------
Nmax = size(x,1); if N > Nmax N = Nmax; end;
disp([' ']); disp(['signal name :     ' name_sig]);
disp(['size signal : ', num2str(size(xy))]);
disp([' ']); disp(['standard wavelet :     ' name_wav]);
disp(['C (svm parameter) = ',num2str(C)]);
disp(['N points = ',num2str(N)]);
disp(['nbfeatures / channel = ',num2str(nbfeatures)]);
disp(['nbsig / clas_app = ',num2str(nbsigclas_app)]);
disp(['nbsig / clas_test = ',num2str(nbsigclas_test)]);
disp(['nbsubset (external CV) = ',num2str(nbsubset)]);
% disp(['nbloc (learning CV) = ',num2str(nbloc)]);
disp([' ']);
% -------------------------------------------------------------------------
svmparam(1) = C;
% ----------------------------------------
t=cputime; 
mcrate_app = []; mcrate_test = [];classif_cat = []; classifmean = 0;

% External Cross-Validation (to better evaluate the the performances of the method)
for i = 1:nbsubset,
  % currentX is the current training set
    % start and end index of current test set
    ind1 = 1+(i-1)*nbsigclas_test;
    ind2 = i*nbsigclas_test;
    currentInd = [1:(ind1-1), (ind2+1):effsize];
    currentX = xy (:,currentInd,:,:);
    currentY = xy (:,ind1:ind2,:,:);
    if (nbsubset == 1),
      currentX = x;
      currentY = y;
    end;
    disp(['size learning_signal : ', num2str(size(currentX))])
    
    for K = 1:nbvoies,    
        h_cat(K,:) = h_catK;
    end;
   
    % representation = marginals of all the channels
    Rx_cat = marg_dwt_concat_voies (currentX,nbfeatures,h_cat,deb); % with catalog wavelet
    % representation = power per band of all the channels
    Sx = specsubband_concat_voies (currentX, nbfeatures); % power per band

% learning step
    [net_cat] = apprenti_SVM_OVRs (Rx_cat, svmparam);
    [net_spec] = apprenti_SVM_OVRs (Sx, svmparam);

% Test on the learning set (x)  
    [mc_app_cat, classif] = classer_ovr (Rx_cat,net_cat);
    [mc_app_spec, classif] = classer_ovr (Sx,net_spec);
    
% Test on the test set (y)  
    Ry_cat = marg_dwt_concat_voies (currentY,nbfeatures,h_cat,deb);
    Sy = specsubband_concat_voies (currentY, nbfeatures);
    
    [mc_test_spec, classif] = classer_ovr (Sy,net_spec);
    [mc_test_cat, classif] = classer_ovr (Ry_cat,net_cat);
    classif_cat = [classif_cat; classif];
    classifmean = classifmean + classif;
    mcrate_app_cat(i) = mc_app_cat;
    mcrate_test_cat(i) = mc_test_cat;
    mcrate_app_spec(i) = mc_app_spec;
    mcrate_test_spec(i) = mc_test_spec;
end;
% End External Cross-Validation (to evaluate the the performances of the method)
% 
classifmean = classifmean/nbsubset;
computation_time=cputime-t;
disp([' '])
disp(['%tage mal classes apprentissage']);
disp(['cat            spec']);
disp([num2str([mcrate_app_cat' mcrate_app_spec'])]);
disp([' ']);
disp(['average']);
disp([num2str([mean(mcrate_app_cat) mean(mcrate_app_spec)])])
disp([' '])
disp(['%tage mal classes test']);
disp(['cat            spec']);
disp([num2str([mcrate_test_cat' mcrate_test_spec'])]);
disp([' ']);
disp(['average']);
disp([num2str([mean(mcrate_test_cat) mean(mcrate_test_spec)])])
disp([' '])

disp(['computation_time : ',num2str(computation_time)]); disp([' ']);

classif_cat
classifmean
probamean = classifmean/nbsigclas_test
