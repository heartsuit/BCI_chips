clc
clear all
close all
clc
addpath(genpath(cd));
rand('state',0);
for i=1:6
    filename = ['imaginary_' num2str(i)];
    programme_test_descripteur_xval_imaginary(filename);
    disp('-----------------------------------------------------------------');
end

rmpath(genpath(cd));