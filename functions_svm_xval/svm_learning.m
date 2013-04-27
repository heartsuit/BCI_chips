function class = svm_learning(x_learn,y_learn,class,nbclasses,kernelparam)

%SVM parameters
kernelparam.ktype = 1;      % Choix du type de noyau %voir la fonction kernel pour les différents types
kernelparam.kscale = 0.1;   % 1/(2*sigma)^2 pour le noyau gaussien
C = 100;


if isempty(class)
    for i=1:nbclasses
        class(i).a = [];                         % coefficients alpha
        class(i).b = 0;                          % Le biais
    end
end


for i=1:nbclasses
    y_learn_ovr(y_learn == i,1) = 1;
    y_learn_ovr(y_learn ~= i,1) = -1;
    [class(i).a,class(i).b] = SVM_rakoto(x_learn,y_learn_ovr,C,kernelparam,[class(i).a; 0]);
end