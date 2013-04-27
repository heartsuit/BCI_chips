function [a,b] = SVM_rakoto(x_current,y_current,C,kernelparam,a_init)

y_current = y_current(:);
vec_ones = ones(length(y_current),1);
Q = (y_current*y_current').*kernel(x_current,x_current,kernelparam);
[a_new, b_new, pos] = monqp(Q,vec_ones,y_current,0,C*vec_ones,1e-6,0,a_init);
a = zeros(length(y_current),1);
a(pos) = a_new;
b = b_new;