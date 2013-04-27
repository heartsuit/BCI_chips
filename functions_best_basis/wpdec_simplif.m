function tree_wp = wpdec_simplif(x,mode,n,Lo_D,Hi_D,shifta,shiftd)
%le vecteur tree_wp est construit sur une structure comportant deux chalps:
%le champs 'data' recueillant les coef de la d�composition par base index�
%par leur profondeur sur l'arbre et rang sur la ligne et par le champ 'L'
%donnant le nombre de coefficients par noeud.


% Initialization.
x = x(:)'; % row vector
L = length(x);

tree_wp(1:2^(n)-1) = struct('data',0,'L',0);

%%%%%%%%%%%%%%%%%%%%%%%%% PRECONCEPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%le coef AO est �gal aux approximation de x(t) par �chantillonage
tree_wp(1).data = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tree_wp(1).L = L;
for k = 0:n-1
    for d = 0:2^k-1
        [cA,cD] = dwt_simplif(tree_wp(node(k,d)).data,mode,Lo_D,Hi_D,shifta,shiftd); % decomposition
        tree_wp(node(k+1,2*d)).data = cA;
        tree_wp(node(k+1,2*d)).L = length(cA);
        tree_wp(node(k+1,2*d+1)).data = cD;
        tree_wp(node(k+1,2*d+1)).L = length(cD);
    end
end
% h = Lo_D(end:-1:1);
% dpo = WPAnalysis(x,n,h);
% for j=0:n
%     for b=0:2^j-1
%         tree_wp(node(j,b)).data = dpo(packet(j,b,L),j+1);
%         tree_wp(node(j,b)).data = length(packet(j,b,L));
%     end
% end


