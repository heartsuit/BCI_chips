function y = dwt_coeff(c,l,ind)


y = [];
for i=ind
     y = [y c((sum(l(1:i-1))+1):sum(l(1:i)))];
end

y = y(:);