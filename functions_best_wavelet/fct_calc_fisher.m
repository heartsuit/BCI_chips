function measure = fct_calc_fisher(x,y)
% function measure = fct_calc_FisherMeasureMono(marginals)
%
% Calculation of the discriminant measure
% Measure = Fisher discriminant measure
%         = distance between means / inertias of classes
%
% size(marginals) = nb des x nb classes x nb trials per class
% 
%

nb_des = size(x,2);
nb_class = max(y);
nb_sig = size(x,1);
x_mean = mean(x,1);
x_mean_class = zeros(nb_class,nb_des);
dis_class = zeros(nb_class,1);
% 1. Distance between means
for i=1:nb_class
    x_mean_class(i,:) = mean(x(y==i,:),1);
    dis_class(i) = (sum(y==i)/length(y))*(x_mean_class(i,:)-x_mean)*(x_mean_class(i,:)-x_mean)';
end
distance = sum(dis_class,1);

inertia = 0;
for i = 1:nb_sig
    inertia = inertia + (x(i,:)-x_mean_class(y(i),:))*(x(i,:)-x_mean_class(y(i),:))';
end
    inertia = inertia/length(y);


measure = distance/inertia;
measure = 1/measure;