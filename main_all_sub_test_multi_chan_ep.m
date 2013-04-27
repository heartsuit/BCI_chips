clc
clear all
close all
clc
addpath(genpath(cd));
rand('state',0);
sub_name = 'Subject 1 New';
filename = 'gediminas_unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\subject1_new\Unvoluntary\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Subject 1 Old';
filename = 'gediminas_unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\subject1_old\Unvoluntary\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Subject 2 New';
filename = 'marija_unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\subject2_new\unvoluntary\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Subject 2 Old';
filename = 'marija_unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\subject2_old\Unvoluntary\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Subject 3';
filename = 'mustafa_unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\subject3\unvoluntary\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Dantan';
filename = 'unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\DanTAn\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Duy';
filename = 'unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\duy\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');

sub_name = 'Jackob';
filename = 'unvoluntary.mat';
pathname = 'C:\Users\Xavier\Desktop\optimisation_algo_dictionnaire\Signaux_EEG_imran\jakob\';
disp(sub_name);
programme_test_descripteur_xval_ep(filename,pathname)
disp('-----------------------------------------------------------------');
rmpath(genpath(cd));