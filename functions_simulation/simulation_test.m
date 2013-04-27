clear all
close all
clc
addpath(genpath(cd));

N = 128;                          % Signal length
h1 = MakeONFilter('Coiflet',4);   % Wavelet for the class 1
h2 = MakeONFilter('Coiflet',5);   % Wavelet for the class 2
proba = 0.05;                     % Parameter of the Bernoulli distribution
nApp = 60;                       % Nb of signals in the training set
nTest = 1000;                     % Nb of signals in the test set
RSB = 10;                         % Signal to noise ratio

[x_app, y_app, x_test, y_test] = fct_signal_simulation(N, h1, h2, nApp, nTest, proba, RSB);





rmpath(genpath(cd));