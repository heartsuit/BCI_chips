function [class_estim,real_label] = approche_classif_mvt_experiment(class,chan,C,filtering,Fe,classifier)

print_result = 0;

global file_name;
N_start = round(4*1024);
N_end = round(6*1024)-1;
% if N_end>size(class(1).signal,2)
%     N_start = N_start - 1536;
%     N_end = N_end - 1536;
% end

if filtering
    srate = Fe;
    N = 6;  %Ordre du filtre
    R = 20;    %Ripple du filtre
    W1 = (2*1)/srate;   %Fréquence de coupure basse (1Hz)
    W2 = (2*4)/srate;  %Fréquence de coupure haute (10Hz)
    
    Wp = [W1 W2];
    [B,A] = butter(N,W2);
%     [B,A] = cheby2(N,R,Wp); %Filtre de chebychev
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
N_tronc = (N_start:1:N_end);
% N_tronc_2 = (N_start:16:N_end);
%descripteur
nbclasses = size(class,2);
% %fft

%dwt sans marginales
[lo_d,hi_d] = wfilters('db8');
for i=1:nbclasses
    nb_trials = size(class(i).signal,1);
    for j=chan
        for k=1:nb_trials
            [c,l] = wavedec_simplif(class(i).signal(k,N_tronc,j),'per',11,lo_d,hi_d,0,0);      
%             class(i).dwt_coeff(k,:,j) = marg(c,l,4);
            nb_des = sum(l(1:3));
            class(i).dwt_coeff(k,:,j) = c(1:nb_des);
        end
    end
    class(i).des = [];
    for j=chan
        class(i).des = [class(i).des, class(i).dwt_coeff(:,:,j)];
    end
end

donnees = [class(1).des; class(2).des];
nomindividu = (1:size(donnees,1))';
classe = [60, 60];

facp09(donnees,nomindividu,classe);
drawnow;
w = 0;
for i = 1:nbclasses
    nb_trials = size(class(i).signal,1);
    for j = 1:nb_trials
        w=w+1;
%         disp(['iter : ' num2str(w)])
        if i == 1
            ind_app = [1:j-1 j+1:nb_trials];
            test = class(i).des(j,:);
            real_label(w) = 1;
            app = [class(1).des(ind_app,:); class(2).des];
            labels = [ones(length(ind_app),1); -ones(size(class(2).des,1),1)];
            ind_permute = randperm(length(labels));
            app = app(ind_permute,:);
            labels = labels(ind_permute);
            if strcmp(classifier,'rakoto')
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
            elseif strcmp(classifier,'lda')
                class_estim(w) = classify(test,app,labels,'linear',[3/4 1/4]);
            elseif strcmp(classifier,'bayes')
                O1 = NaiveBayes.fit(app,labels,'Prior',[3/4 1/4]);
                post(w,:) = O1.posterior(test);
                class_estim(w) = O1.predict(test);
            end
            clear app test labels
        elseif i==2
            ind_app = [1:j-1 j+1:nb_trials];
            test = class(2).des(j,:);
            real_label(w) = -1;
            app = [class(1).des; class(2).des(ind_app,:)];
            labels = [ones(size(class(1).des,1),1); -ones(length(ind_app),1)];
            ind_permute = randperm(length(labels));
            app = app(ind_permute,:);
            labels = labels(ind_permute);
            if strcmp(classifier,'rakoto')
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
            elseif strcmp(classifier,'lda')
                class_estim(w) = classify(test,app,labels,'linear',[3/4 1/4]);
            elseif strcmp(classifier,'bayes')
                O1 = NaiveBayes.fit(app,labels,'Prior',[3/4 1/4]);
                post(w,:) = O1.posterior(test);
                class_estim(w) = O1.predict(test);
            end
            clear net app test labels
        end
    end
end


P_W_W = 100*sum(class_estim == -1 & real_label == -1)/sum(real_label == -1);
P_C_W = 100*sum(class_estim == 1 & real_label == -1)/sum(real_label == -1);
P_C_C = 100*sum(class_estim == 1 & real_label == 1)/sum(real_label == 1);
P_W_C = 100*sum(class_estim == -1 & real_label == 1)/sum(real_label == 1);



fprintf('P_E_estim_E\t | E_estim = Wrong\t | E_estim = Correct\n')
fprintf('E = Wrong\t | %f\t\t\t | %f\n', P_W_W, P_C_W)
fprintf('E = Correct\t | %f\t\t\t | %f\n', P_W_C, P_C_C)


c_q = (P_C_C + P_W_W)/2;
disp(['Critère qualité : ' num2str(c_q)]);
if print_result
    fid = fopen(file_name,'a');
    fprintf(fid,'descriptors : %s \\\\\n',des);
    if strcmp(des,'spatial_filter')
         fprintf(fid,'temporal reajustement : %s \\\\\n',sig_recal);
         fprintf(fid,'Type of filter : %s \\\\\n',spatial_filtering);
    end         
    fprintf(fid,'Table of results : \\\\\n');
    fprintf(fid,'\\begin{tabular}{|c|c|c|}\n\\hline\t\t\t\t& $\\hat{e}$ = juste & $\\hat{e}$ = faux \\\\\n');
    fprintf(fid,'\\hline  $e$ = juste\t& %6.0f\\%%\t\t\t& %6.0f\\%%\t\t\\\\\n', P_C_C,P_W_C);
    fprintf(fid,'\\hline  $e$ = faux\t& %6.0f\\%%\t\t\t& %6.0f\\%%\t\t\\\\\n', P_C_W,P_W_W);
    fprintf(fid,'\\hline\n\\end{tabular}\\\\\n');
    fprintf(fid,'Quality criterion : %6.1f \\\\\n\n',c_q);
    fclose(fid);
end

