function [inds,marg]=findmarg(tree,basis,N_dec)    

    tree_sorted=sort(tree(find(basis)),'descend');
 
    j=1;
    
    if length(tree_sorted)<5
         marg=zeros(length(tree_sorted),1);
        for i=1:length(tree_sorted)
            marg(i)=find(tree==tree_sorted(i),1,'first');
        end
    else
        marg=zeros(5,1);
        for i=1:5


                if tree_sorted(1:5)==zeros(5,1)
                    k=basis(tree==0);
                    zero_tree=find(tree==0);
                    good_zeros=zero_tree(find(k));
                    inds_zeros=nds2ind(good_zeros,N_dec);

                    while marg(i)==0
                        if j==size(inds_zeros,2)+1
                            marg=good_zeros(1:5);
                        else
                            if inds_zeros(2,j)<=2^(inds_zeros(1,j)-2)
                                disp(['hihi ' num2str(inds_zeros(2,j)) ' haha ' num2str(inds_zeros(1,j))]);
                                marg(i)=good_zeros(j);
                            end
                            j=j+1;
                        end
                    end
                else
                    marg(i)=find(tree==tree_sorted(i),1,'first');
                end

        end
    end  
    inds=nds2ind(marg,N_dec);
    
end