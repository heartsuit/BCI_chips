function marginals = calc_marg(c,l,ind)

% c = abs(c);
if nargin < 3 || isempty(ind)
    ind = 1:length(l)-1;
end
marginals = zeros(length(ind),1);

for i=1:length(ind)
    marginals(i) = mean(c((1:l(ind(i)))+sum(l(1:ind(i)-1)))); %Ajout des différents niveaux de détail
end
