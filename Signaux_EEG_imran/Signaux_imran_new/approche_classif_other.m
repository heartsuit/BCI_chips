function approche_classif_other(class,chan,C,filtering,Fe,nb_sub)


N_start = round(2.15*1024);
N_end = round(2.65*1024);
N_sample = 16;

if filtering
    srate = Fe;
    N = 4;  %Ordre du filtre
    R = 20;    %Ripple du filtre
    W1 = (2*1)/srate;   %Fréquence de coupure basse (1Hz)
    W2 = (2*10)/srate;  %Fréquence de coupure haute (10Hz)
    
    Wp = [W1 W2];
    
    [B,A] = butter(N,R,Wp); %Filtre de chebychev
%     [B,A] = cheby2(N,R,W2);
    for i=chan
        for j=1:size(class(1).signal,1)
            class(1).signal(j,:,i) = filter(B,A,class(1).signal(j,:,i));
        end
        for j=1:size(class(2).signal,1)
            class(2).signal(j,:,i) = filter(B,A,class(2).signal(j,:,i));
        end
    end
end

kernelparam.ktype = 1;
kernelparam.kscale = 0.2;
N_tronc = N_start:N_sample:N_end;
%descripteur
nbclasses = size(class,2);

for i=1:nbclasses
    nb_trials = size(class(i).signal,1);
    for j=chan
        for k=1:nb_trials
            class(i).signal_decim(k,:,j) = class(i).signal(k,N_tronc,j);
        end
    end
    class(i).des = [];
    for j=chan
        class(i).des = [class(i).des, class(i).signal_decim(:,:,j)];
    end
end

% class(1).signal_tronc       = class(1).signal(:,N_tronc,:);
% class(2).signal_tronc       = class(2).signal(:,N_tronc,:);
% N                           = length(N_tronc);
% name_wav                    = 'db4';
% h_catK                      = wfilters(name_wav);
% nbchan                      = length(chan);
% h_cat                       = ones(nbchan,1)*h_catK;
% lh                          = length(h_catK);
% deb                         = floor(-log2(1/lh));
% nbfeatures                  = floor(log2(N))-deb;
% class = marg_dwt_concat_voies(class,nbfeatures,h_cat,deb,chan);



for i=1:nb_sub
    ind_1 = [zeros((i-1)*90,1); ones(90,1); zeros((nb_sub-i)*90,1)];
    ind_2 = [zeros((i-1)*30,1); ones(30,1); zeros((nb_sub-i)*30,1)];
    test = cat(1,class(1).des(ind_1 == 1,:),class(2).des(ind_2 == 1,:));
    real_label = [-ones(sum(ind_1 == 1),1);ones(sum(ind_2 == 1),1)];
    app = cat(1,class(1).des(ind_1 ~= 1,:),class(2).des(ind_2 ~= 1,:));
    labels_app = [-ones(sum(ind_1 ~= 1),1);ones(sum(ind_2 ~= 1),1)];
    
    ind_permute = randperm(length(labels_app));
    app = app(ind_permute,:);
    labels_app = labels_app(ind_permute);
    O1 = NaiveBayes.fit(app,labels_app);
    class_estim = O1.predict(test);
%     [a,b] = SVM_rakoto(app,labels_app,C,kernelparam);
%     class_estim = test_class(test,app,labels_app,a,b,kernelparam);
    
    P_W_W = 100*sum(class_estim == 1 & real_label == 1)/sum(real_label == 1);
    P_C_W = 100*sum(class_estim == -1 & real_label == 1)/sum(real_label == 1);
    P_C_C = 100*sum(class_estim == -1 & real_label == -1)/sum(real_label == -1);
    P_W_C = 100*sum(class_estim == 1 & real_label == -1)/sum(real_label == -1);
    
    disp(num2str(i))
    fprintf('P_E_estim_E\t | E_estim = Wrong\t | E_estim = Correct\n')
    fprintf('E = Wrong\t | %f\t\t\t | %f\n', P_W_W, P_C_W)
    fprintf('E = Correct\t | %f\t\t\t | %f\n', P_W_C, P_C_C)
    disp(' ')

end 

