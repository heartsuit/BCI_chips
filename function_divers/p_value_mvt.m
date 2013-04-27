function p_value = p_value_mvt(p,q,P_err,n)

k=round((1-P_err)*n);
p_value = (factorial(n)/(factorial(k)*factorial(n-k)))*p^k*q^(n-k);