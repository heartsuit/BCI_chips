function tree_marg = wp_marg(tree_wp,type)

K = size(tree_wp,2);
tree_marg = zeros(size(tree_wp));

if strcmp(type,'abs')
    for i=1:K
        wp_i_coeff = tree_wp(i).data;
        tree_marg(i) = sum(abs(wp_i_coeff));
%         tree_marg(i) = mean(wp_i_coeff);
    end
end