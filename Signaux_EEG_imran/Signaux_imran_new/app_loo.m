function app_loo(class)

w = 0;
for i = 1:nbclasses
    nb_trials = size(class(i).signal,1);
    for j = 1:nb_trials
        w=w+1;
        if i == 1
            ind_app = [1:j-1 j+1:nb_trials];
            test = class(i).des(j,:);
            real_label(w) = 1;
            app = [class(1).des(ind_app,:); class(2).des];
            labels = [ones(length(ind_app),1); -ones(size(class(2).des,1),1)];
            ind_permute = randperm(length(labels));
            app = app(ind_permute,:);
            labels = labels(ind_permute);
            net = svm(size(app,2), 'linear', [], C); % init
            warning off
            net = svmtrain(net, app, labels); % learning
            warning on
            [class_estim(w),dist(w)]= svmfwd(net, test);
            clear net app test labels
        elseif i==2
            ind_app = [1:j-1 j+1:nb_trials];
            test = class(i).des(j,:);
            real_label(w) = -1;
            app = [class(1).des; class(2).des(ind_app,:)];
            labels = [-ones(size(class(2).des,1),1); ones(length(ind_app),1)];
            ind_permute = randperm(length(labels));
            app = app(ind_permute,:);
            labels = labels(ind_permute);
            net = svm(size(app,2), 'linear', [], C); % init
            warning off
            net = svmtrain(net, app, labels); % learning
            warning on
            [class_estim(w),dist(w)]= svmfwd(net, test);
            clear net app test labels
        end
    end
end