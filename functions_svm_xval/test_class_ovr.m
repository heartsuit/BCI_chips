function classe = test_class_ovr(xtest, class, nbclasses, xapp, yapp, kernelparam)

Ntest = size(xtest,1);
f = zeros(Ntest,nbclasses);
Kindc = kernel(xtest,xapp,kernelparam);
Y = zeros(size(yapp));
for i = 1:nbclasses
    Y(yapp==i,1) = 1;
    Y(yapp~=i,1) = -1;
    a = class(i).a;
    b = class(i).b;
    f(:,i) = Kindc*(Y.*a)+b;
end

[m,classe] = max(f,[],2);