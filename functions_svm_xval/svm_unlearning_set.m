function [class,x_learn,y_learn] = svm_unlearning_set(x_learn,y_learn,indrem,class,nbclasses,kernelparam)

%SVM parameters
kernelparam.ktype = 1;      % Choix du type de noyau %voir la fonction kernel pour les différents types
kernelparam.kscale = 0.1;   % 1/(2*sigma)^2 pour le noyau gaussien
C = 10;
x_learn = x_learn(~ismember(1:end,indrem),:);
y_learn = y_learn(~ismember(1:end,indrem));
for i=1:nbclasses
    a_init = class(i).a(~ismember(1:end,indrem));
    y_learn_ovr(y_learn == i,1) = 1;
    y_learn_ovr(y_learn ~= i,1) = -1;
    [class(i).a,class(i).b] = SVM_rakoto(x_learn, y_learn_ovr, C, kernelparam, a_init);
end
