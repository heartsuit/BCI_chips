function y = marg(c,l,ind)


y = [];
% for i=ind
%      y = [y c((sum(l(1:i-1))+1):sum(l(1:i)))];
% end
w = 0;
for i=ind
    w = w+1;
     y(w,1) = sum(abs(c((sum(l(1:i-1))+1):sum(l(1:i)))));
end


y = y(:);