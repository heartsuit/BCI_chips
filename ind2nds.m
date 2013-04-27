function nodes=nds2ind(ind,J)
%gives the depth (nodes(1,1)) and the index (nodes(2,1)) from left to right
%of a binar tree
    nodes=zeros(2,length(ind));
    for i=1:length(ind)
        div=2^J;    
        while(ind(i)/div<1)
             div=div/2;
        end
    nodes(1,i)=log2(div);
    nodes(2,i)=mod(ind(i),div);
    end
end