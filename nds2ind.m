function inds=nds2ind(ind,J)
%gives the depth (inds(1,1)) and the index (inds(2,1)) from left to right
%of a binar tree with an depth of J.
    inds=zeros(2,length(ind));
    for i=1:length(ind)
        div=2^J;    
        while(ind(i)/div<1)
             div=div/2;
        end
    inds(1,i)=log2(div);
    inds(2,i)=mod(ind(i),div);
    end
end