p=0.5;
q=0.5;
n=120;
P_1_1= [92 91 93 87 83 80 88]/100;
P_2_2= [73 43 60 60 37 63 56]/100;

for i=1:length(P_1_1)
    p_ep(i) =  p_value_ep(p,q,P_1_1(i),P_2_2(i),n);
end


p=0.5;
q=0.5;
n=120;
P_err = [28 24 32 30 22 22 26]/100;
for i=1:length(P_err)
    p_mvt(i) = p_value_mvt(p,q,P_err(i),n);
end


x = [0.14 0.20 0.10 0.12 0.24 0.24 0.17;
	0.39 0.29 0.26 0.26 0.26 0.34 0.30]';