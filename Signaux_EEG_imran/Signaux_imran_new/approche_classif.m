function [class_estim,real_label] = approche_classif(class,chan,C,des,csp,recal,print_result)


global file_name;
N_start = 2200;
N_end = 3000;

srate = 1024;
N = 3;  %Ordre du filtre
R = 20;    %Ripple du filtre
W1 = (2*1)/srate;   %Fr�quence de coupure basse (1Hz)
W2 = (2*10)/srate;  %Fr�quence de coupure haute (10Hz)

Wp = [W1 W2];

[B,A] = cheby2(N,R,Wp); %Filtre de chebychev
% [B,A] = cheby2(N,R,W1,'low');
for i=chan
    for j=1:size(class(1).signal,1)
        class(1).signal(j,:,i) = filter(B,A,class(1).signal(j,:,i));
    end
    for j=1:size(class(2).signal,1)
        class(2).signal(j,:,i) = filter(B,A,class(2).signal(j,:,i));
    end
end



kernelparam.ktype = 1;
kernelparam.kscale = 0.2;
N_tronc = N_start:16:N_end;
%descripteur
nbclasses = size(class,2);

if strcmp(des,'spatial_filter')
    if csp == 1
        spatial_filtering = 'CSP';
    elseif csp == 2
        spatial_filtering = 'Fisher';
    else
        spatial_filtering = 'Mean of channels';
    end

    w=0;
    for i=1:nbclasses
        if i==1
            if recal == 1
                sig_recal = 'Yes';
                ind_rec = recalage(class(2).signal,N_tronc);
                for k=1:size(class(2).signal,1)
                    class(2).signal_tronc(k,:,:) = class(2).signal(k,N_tronc-ind_rec(k),chan);
                end
                class(1).signal_tronc = class(1).signal(:,N_tronc,chan);
            else
                sig_recal = 'No';
                class(1).signal_tronc = class(1).signal(:,N_tronc,chan);
                class(2).signal_tronc = class(2).signal(:,N_tronc,chan);
            end
            nb_trials = size(class(i).signal,1);
            for j=1:nb_trials
                w=w+1;
                disp(['iter : ' num2str(w)])
                %             signal_test = squeeze(class(i).signal(j,N_tronc,:));
                ind_app = [1:j-1 j+1:nb_trials];
                class_spa_filter(1).signal = class(1).signal_tronc(ind_app,:,:);
                class_spa_filter(2).signal = class(2).signal_tronc;
                if csp == 1
                    f = CSP(class_spa_filter);
                elseif csp == 2
                    f = erp_satial_filter(class_spa_filter);
                else
                    f = ones(1,length(chan))/length(chan);
                end
                for h=ind_app
                    class(1).des(h,:) = decimate(f*(squeeze(class(1).signal_tronc(h,:,:)))',50);
                end
                class(1).des(j,:) = decimate(f*(squeeze(class(1).signal_tronc(j,:,:)))',50);
                for h=1:size(class(2).signal_tronc,1)
                    class(2).des(h,:) = decimate(f*(squeeze(class(2).signal_tronc(h,:,:)))',50);
                end
                test = class(1).des(j,:);
                real_label(w) = 1;
                app = [class(1).des(ind_app,:); class(2).des];
                labels = [ones(length(ind_app),1); -ones(size(class(2).des,1),1)];
                ind_permute = randperm(length(labels));
                app = app(ind_permute,:);
                labels = labels(ind_permute);
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
                clear net app test labels
            end
        else
            nb_trials = size(class(i).signal,1);
            for j=1:nb_trials
                w=w+1;
                disp(['iter : ' num2str(w)])
                ind_app = [1:j-1 j+1:nb_trials];
                signal_app_inter = class(2).signal(ind_app,:,:);
                if recal == 1
                    sig_recal = 'Yes';
                    ind_rec = recalage(signal_app_inter,N_tronc);
                    for k=1:length(ind_app)
                        class(2).signal_tronc(ind_app(k),:,1:length(chan)) = signal_app_inter(k,N_tronc-ind_rec(k),chan);
                    end
                    class(2).signal_tronc(j,:,:) = class(2).signal(j,N_tronc,chan);
                    class(1).signal_tronc = class(1).signal(:,N_tronc,chan);
                else
                    sig_recal = 'No';
                    class(1).signal_tronc = class(1).signal(:,N_tronc,chan);
                    class(2).signal_tronc = class(2).signal(:,N_tronc,chan);
                end
                class_spa_filter(1).signal = class(1).signal_tronc;
                class_spa_filter(2).signal = class(2).signal_tronc;
                if csp == 1
                    f = CSP(class_spa_filter);
                elseif csp == 2
                    f = erp_satial_filter(class_spa_filter);
                else
                    f = ones(1,length(chan))/length(chan);
                end
                for h=1:size(class(1).signal_tronc,1)
                    class(1).des(h,:) = decimate(f*(squeeze(class(1).signal_tronc(h,:,:)))',50);
                end
                for h=ind_app
                    class(2).des(h,:) = decimate(f*(squeeze(class(2).signal_tronc(h,:,:)))',50);
                end
                class(2).des(j,:) = decimate(f*(squeeze(class(2).signal_tronc(j,:,:)))',50);
                test = class(1).des(j,:);
                real_label(w) = -1;
                app = [class(1).des; class(2).des];
                labels = [ones(size(class(1).des,1),1); -ones(length(ind_app),1)];
                ind_permute = randperm(length(labels));
                app = app(ind_permute,:);
                labels = labels(ind_permute);
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
                clear app test labels
            end
        end
    end
else
    if strcmp(des,'concat')
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
    elseif strcmp(des,'pca')
%         for i=1:nbclasses
%             nb_trials = size(class(i).signal,1);
%             for j=chan
%                 for k=1:nb_trials
%                     class(i).signal_decim(k,:,j) = decimate(class(i).signal(k,N_tronc,j),50);
%                 end
%             end
%             class(i).des = [];
%             for j=chan
%                 class(i).des = [class(i).des, class(i).signal_decim(:,:,j)];
%             end
%         end
    elseif strcmp(des,'dwt')
        N_start = 2244;
        N_end = 2755;
        N_tronc = N_start:N_end;
        class(1).signal_tronc       = class(1).signal(:,N_tronc,:);
        class(2).signal_tronc       = class(2).signal(:,N_tronc,:);
        N                           = length(N_tronc);
        name_wav                    = 'db4';
        h_catK                      = wfilters(name_wav);
        nbchan                      = length(chan);
        h_cat                       = ones(nbchan,1)*h_catK;
        lh                          = length(h_catK);
        deb                         = floor(-log2(1/lh));                
        nbfeatures                  = floor(log2(N))-deb;
        class = marg_dwt_concat_voies(class,nbfeatures,h_cat,deb,chan);
    end
    
    w = 0;
    for i = 1:nbclasses
        nb_trials = size(class(i).signal,1);
        for j = 1:nb_trials
            w=w+1;
            disp(['iter : ' num2str(w)])
            if i == 1
                ind_app = [1:j-1 j+1:nb_trials];
                test = class(i).des(j,:);
                real_label(w) = 1;
                app = [class(1).des(ind_app,:); class(2).des];
                labels = [ones(length(ind_app),1); -ones(size(class(2).des,1),1)];
                ind_permute = randperm(length(labels));
                app = app(ind_permute,:);
                labels = labels(ind_permute);
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
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
                [a,b] = SVM_rakoto(app,labels,C,kernelparam);
                class_estim(w) = test_class(test,app,labels,a,b,kernelparam);
%                 net = svm(size(app,2), 'linear', [], C); % init
%                 warning off
%                 net = svmtrainiter100(net, app, labels); % learning
%                 warning on
%                 [class_estim(w),dist(w)]= svmfwd(net, test);
                clear net app test labels
            end
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
disp(['Crit�re qualit� : ' num2str(c_q)]);
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

