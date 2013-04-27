function dwt_coeff = calc_dwt_coeff(c,l,ind)

if nargin < 3 || isempty(ind)
    ind = 1:length(l)-1;
end
% dwt_coeff = zeros(sum(l(ind)),1);
dwt_coeff = [];
for i=1:length(ind)
    dwt_coeff = [dwt_coeff c((1:l(ind(i)))+sum(l(1:ind(i)-1)))]; %Ajout des différents niveaux de détail
end
