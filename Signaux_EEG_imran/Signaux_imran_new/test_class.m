function class = test_class(xtest,xapp,yapp,a,b,kernelparam)

Kindc = kernel(xtest,xapp,kernelparam);
f = Kindc*(yapp.*a)+b;
class = sign(f);