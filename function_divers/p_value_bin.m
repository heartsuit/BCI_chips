function p_value = p_value_bin(p,q,P_1_1,P_2_2,n)

k=round((p*P_1_1+q*P_2_2)*n);
p_value = (factorial(n)/(factorial(k)*factorial(n-k)))*p^k*q^(n-k);