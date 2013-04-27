clear all
clc

load imaginary

rand('state',0);


for i = 1:6
    [nb_chan_fast nb_samples_fast,nb_fast] = size(imaginary(i).ballistic);
    [nb_chan_slow nb_samples_slow,nb_slow] = size(imaginary(i).slow);
    nb_samples = min([nb_samples_fast nb_samples_slow]);
    x_mvt = cat(1,permute(imaginary(i).ballistic(:,1:nb_samples,:),[3 2 1]),permute(imaginary(i).slow(:,1:nb_samples,:),[3 2 1]));
    y_mvt = [ones(nb_fast,1); 2*ones(nb_slow,1)];
    n_perm = randperm(length(y_mvt));
    x_mvt = x_mvt(n_perm,:,:);
    y_mvt = y_mvt(n_perm,:,:);
    file_save = ['imaginary_' num2str(i)];
    save(file_save,'x_mvt','y_mvt');
end