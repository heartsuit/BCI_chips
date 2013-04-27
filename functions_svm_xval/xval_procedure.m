function [y_estim, y_real] = xval_procedure(x_learn,y_learn,nb_subsets)
kernelparam.ktype = 1;
kernelparam.kscale = 0.2;

nb_class = max(y_learn);
%initialisation de la séparatrice avec tous les élément de l'apprentissage
sep = svm_learning(x_learn,y_learn,[],nb_class,kernelparam);

subset_length = floor(length(y_learn)/nb_subsets);

%test loo
y_estim = [];
y_real = [];
%Pour tous les éléments de l'apprentissage on en retire un on met à jour la
%séparatrice et on retest l'élément retiré.
for ii = 1:nb_subsets
    ind_test = (ii-1)*subset_length+1:ii*subset_length;
%     ind_app = setdiff(1:length(y_learn),ind_test);
    x_test = x_learn(ind_test,:);
    y_real = [y_real; y_learn(ind_test)];
    [new_sep,x_loo,y_loo] = svm_unlearning_set(x_learn,y_learn,ind_test,sep,nb_class,kernelparam); %retire l'élément i calcul la nouvelle séparatrice et retourne l'apprentissage sans l'élément i
    y_estim = [y_estim; test_class_ovr(x_test, new_sep, nb_class, x_loo, y_loo, kernelparam)]; %estime la classe de l'élément retiré avec la séparatrice apprise précédemment
end