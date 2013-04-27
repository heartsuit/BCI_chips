function trial_drawing(feat,label,marg,N_dec)
%marg doit être un cell de matrice contenant un ou deux éléments

marker_color = {'bx','r.','g*','ksquare','s'};

for j=1:length(marg)
    if length(marg{j})==1
        inds=nds2ind(marg{j},N_dec);
        figure;
        title(['trials represented on  M_{' num2str(inds(1,1)) '}^{' num2str(inds(2,1)) '}']);
        hold on
        for p=1:max(label)
            

            plot(feat(label==p,marg{j},:),ones(length(find(label==p)),1),marker_color{p},'LineWidth',1.5);
        
        end
    else
        inds=nds2ind(marg{j},N_dec);
        figure;
        hold on
        title(['trials represented on  M_{' num2str(inds(1,1)) '}^{' num2str(inds(2,1)) '} and M_{' num2str(inds(1,2)) '}^{' num2str(inds(2,2)) '}']);
        for p=1:max(label)
            
            plot(feat(label==p,marg{j}(1),:),feat(label==p,marg{j}(2),:),marker_color{p},'LineWidth',1.5);
        %comment git
        end
    end
end