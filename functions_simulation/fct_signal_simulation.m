function [trainingSet, y_train, testSet, y_test] = fct_signal_simulation(N, h1, h2, nTrain, nTest, proba, SNR)

% [trainingSet, testSet] = fct_signal_simulation(N, h1, h2, nTrain, nTest, proba, SNR)
%
% Simulation of signals for the wavelet packet decomposition analysis
%
% Two classes of mono-channel signals are generated from their wavelet packet decomposition.
% We use the probability of Bernoulli to generate the wavelet packet decompositions.
%
% N = length of each signal
% h1 = wavelet used to generate the signal of the first class
% h2 = wavelet used to generate the signal of the second class
% nTrain = number of signals in the training set 
% nTest = number of signals in the test set 
% proba = probability of success for the Bernoulli process (between 0 and 1)
% SNR = signal to noise ratio (dB)
%
% size(trainingSet) = N x nTrain x nbchannels (=1) x nbclasses (=2)
% size(testSet) = N x nTest x nbchannels (=1) x nbclasses (=2)
%


basis1 = [               1                 ...
                 0               1         ...
             0       0       1       1     ...
           0   0   0   0   0   0   0   0   ...
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

basis2 = [               1                 ...
                 0               1         ...
             0       0       1       1     ...
           0   0   0   0   0   0   0   0   ...
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];



trainingSet = zeros(nTrain,N);
y_train = zeros(nTrain,1);
h = 0;
for i = 1:nTrain
    h = h+1;
    wpdClass1 = zeros(N,5);
    wpdClass1(N/2+1:3*N/4, 4) = (rand(N/4,1)<proba);
    signal = WPSynthesis(basis1,wpdClass1,h1);
    noise = std(signal)*10^(-SNR/20)*randn(1,length(signal));
    trainingSet(h,:) = signal+noise;
    y_train(h) = 1;
    
    h = h+1;
    wpdClass2 = zeros(N,5);
    wpdClass2(5*N/8+1:7*N/8, 4) = (rand(N/4,1)<proba);
    signal = WPSynthesis(basis2,wpdClass2,h2);
    noise = std(signal)*10^(-SNR/20)*randn(1,length(signal));
    trainingSet(h,:) = signal+noise;
    y_train(h) = 2;
end


testSet = zeros(nTest,N);
y_test = zeros(nTest,1);

h=0;
for i = 1:nTest
    h = h+1;
    wpdClass1 = zeros(N,5);
    wpdClass1(N/2+1:3*N/4, 4) = (rand(N/4,1)<proba);
    signal = WPSynthesis(basis1,wpdClass1,h1);
    noise = std(signal)*10^(-SNR/20)*randn(1,length(signal));
    testSet(h,:) = signal+noise;
    y_test(h) = 1;
    
    h = h+1;
    wpdClass2 = zeros(N,5);
    wpdClass2(5*N/8+1:7*N/8 ,4) = (rand(N/4,1)<proba);
    signal = WPSynthesis(basis2,wpdClass2,h2);
    noise = std(signal)*10^(-SNR/20)*randn(1,length(signal));
    testSet(h,:) = signal+noise;    
    y_test(h) = 2;
end


n_perm_train = randperm(length(y_train));
trainingSet = trainingSet(n_perm_train,:);
y_train = y_train(n_perm_train);

n_perm_test = randperm(length(y_test));
testSet = testSet(n_perm_test,:);
y_test = y_test(n_perm_test);
