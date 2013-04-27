function basis = find_best_basis_fisher(tree,N_dec,plot_tr,sim_sig,restric)
%l'output est un tableau contenant des 0 et des 1 pour chaque noeud, les un sont les noeud retenu pour la base 
%attention lorsque plot_ est à 1, le résultat affichera autant de figure
%que de subset!! (subset=120)...


%tree : fisher_coeff
nb_nodes = length(tree);
basis = ones(1,nb_nodes);
value = tree;

if restric
    for d=N_dec-1:-1:2,
        for b=0:(2^d-3),
            vparent = tree(node(d,b));
            vchild  = value(node(d+1,2*b)) + value(node(d+1,2*b+1));
            if(vparent < vchild),
                basis(node(d,b)) = 0;
                value(node(d,b)) = vchild;
            else
                for d_child = d+1:N_dec
                    for b_child=b*2^(d_child-d):(b+1)*2^(d_child-d)-1
                        basis(node(d_child,b_child)) = 0;
                    end
                end
            end
        end
    end
    bad_nodes=[1,2,3];
    for i=2:N_dec
        bad_nodes=[bad_nodes [2^i+2^(i-2):2^(i+1)-1]];
    end   
    basis(bad_nodes)=0;
else
    for d=N_dec-1:-1:0,
        for b=0:(2^d-1),
            vparent = tree(node(d,b));
            vchild  = value(node(d+1,2*b)) + value(node(d+1,2*b+1));
            if(vparent < vchild),
                basis(node(d,b)) = 0;
                value(node(d,b)) = vchild;
            else
                for d_child = d+1:N_dec
                    for b_child=b*2^(d_child-d):(b+1)*2^(d_child-d)-1
                        basis(node(d_child,b_child)) = 0;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%% FORCED VALUES FOR BEST BASIS %%%%%%%%%%%%%%%%%%%%%%%%%
%         basis(node(9,1))=0;        
%         basis(node(11,4))=1;
%         basis(node(11,5))=1;   
%         basis(node(10,3))=1;   
      
%         basis(node(10,1))=0;   
%         basis(node(11,2))=1;
%         basis(node(11,3))=1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
if plot_tr.on
   
   color_={'-r','-b','-g','-c','-m','-y','-black','--r','--b','--g','--c','--m','--y','--k',':r',':b',':g',':c',':m',':y',':k','-.r','-.b','-.g','-.c','-.m','-.y','-.k'};
   weft=1:(nb_nodes-1)/2;
   weft=[weft;weft];weft=weft(:).';
   [x_weft,y_weft]=treelayout(cat(2,0,weft));
   
   %%%%%%%%%%%%%%%%%%%%%%  WARNING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  il faut modifier treelayout:
   %ligne 81: x = deltax * (xmin+xmax)/2-deltax/2;
   %ligne 79: deltax = 1/(nleaves);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if restric
        gd_nodes=[];
        for i=N_dec-2:N_dec
            gd_nodes=[gd_nodes [2^i:2^i+2^(i-(N_dec-3))-1]];
        end
        x_pre_fail=x_weft.*(-(basis-1));
        x_pre_fail=x_pre_fail(gd_nodes);
        ind_fail_x=x_pre_fail~=0; 
        x_fail=x_pre_fail(ind_fail_x);

        y_pre_fail=y_weft.*(-(basis-1));
        y_pre_fail=y_pre_fail(gd_nodes);
        ind_fail_y=y_pre_fail~=0; 
        y_fail=y_pre_fail(ind_fail_y);

        x_pre_success=x_weft.*basis;
        x_pre_success=x_pre_success(gd_nodes);
        ind_success_x=x_pre_success~=0; 
        x_success=x_pre_success(ind_success_x);

        y_pre_success=y_weft.*basis;
        y_pre_success=y_pre_success(gd_nodes);
        ind_success_y=y_pre_success~=0; 
        y_success=y_pre_success(ind_success_y);       
    else
   
       x_pre_fail=x_weft.*(-(basis-1));
       ind_fail_x=x_pre_fail~=0; 
       x_fail=x_pre_fail(ind_fail_x);

       y_pre_fail=y_weft.*(-(basis-1));
       ind_fail_y=y_pre_fail~=0; 
       y_fail=y_pre_fail(ind_fail_y);

       x_pre_success=x_weft.*basis;
       ind_success_x=x_pre_success~=0; 
       x_success=x_pre_success(ind_success_x);

       y_pre_success=y_weft.*basis;
       ind_success_y=y_pre_success~=0; 
       y_success=y_pre_success(ind_success_y);
   end
   
   figure;
   hold on
   title('tree decomposition');
   
   plot(x_fail*plot_tr.fe/2-0.5,y_fail,'s','MarkerFaceColor','b');

   plot(x_success*plot_tr.fe/2-0.5,y_success,'s','MarkerFaceColor','g','MarkerEdgeColor','g');

    if sim_sig
        plot(plot_tr.f,(min(y_weft)/2)*ones(length(plot_tr.f),1),...
            's','MarkerFaceColor','c','MarkerSize',5); 
        %comparaison base choisi signaux à discerner représenté par leur paquet de fréquences caractéristiques

        plot(plot_tr.f_cl,(min(y_weft)*3/4)*ones(length(plot_tr.f_cl),1),...
            's','MarkerFaceColor','r','MarkerSize',5); 
        %comparaison base choisi signaux à discerner représenter par leur fréquence centrale
   
        label{1}='useless M';
        label{2}='used M';
        label{3}='trials central freq.';
        label{4}='central class freq.';
        
    else
        
        label{1}='useless M';
        label{2}='used M';
        label{3}='class1 trials';
        label{4}='class2 trials';
        
        plot(plot_tr.f_scal,plot_tr.f1-1,'color','g');
        plot(plot_tr.f_scal,plot_tr.f2-1,'color','m');  
    end
   
%%%%%%%%%%%%%%%%%%%%%%%% WAVELET FFT PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [inds,marg_ind]=findmarg(tree,basis,N_dec); 
        
%         for i=1:size(inds,2)
          for i=1:2
            %je pars du principe que les signaux durent 2 secondes
            for k=0:2*plot_tr.fe/2^(inds(1,i))-1
                input(k+1,:)= MakeWaveletPacket(inds(1,i),inds(2,i),k,'Daubechies',4,2*plot_tr.fe); 
            end
            [fft,f_scal]=fftsum(input,plot_tr.fe);
            plot(f_scal,fft-1,color_{i});
        
            label{i+4}=['M\phi_{' num2str(inds(1,i)) '}^{' num2str(inds(2,i)) '}'];
        end


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    legend(label);   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%% TEXT SHAPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   h=num2str(tree,'% .2g');
   for i=1:size(h,1)
       k=strfind(h(i,:),'e-0');
       if k~=0
           C2(i,:)=[h(i,1:k+1),h(i,k+3:end) ' '];
       else
           k=strfind(h(i,:),'e+0');
           if k~=0
                C2(i,:)=[h(i,1:k+1),h(i,k+3:end) ' '];
           else
                C2(i,:)=h(i,:);
           end
       end
       kk=strfind(C2(i,1),' ');
       if kk~=0
           if (kk==1 || kk==2)
                 C2(i,:)=[C2(i,2:end) ' '];
           end
       end
   end
   for i=1:size(C2,1)
       k=strfind(C2(i,:),'e-0');
       if k~=0
            C2(i,:)=[C2(i,1:k+1),C2(i,k+3:end) ' '];
       else
           k=strfind(C2(i,:),'e+0'); 
           if k~=0
                C2(i,:)=[h(i,1:k+1),h(i,k+3:end) ' '];
           end
       end
       kk=strfind(C2(i,1),' ');
       if kk~=0 
           if (kk==1 || kk==2)
                C2(i,:)=[C2(i,2:end) ' '];
           end
       end
   end
       
   for p=1:length(y_weft)
       if y_weft(p)>0.4
           y_C2(p)=y_weft(p)+(-1)^p*(y_weft(1)-y_weft(2))/4;
       elseif y_weft(p)<=0.4 && y_weft(p)>0.3
           y_C2(p)=y_weft(p)+(-1)^p*(y_weft(1)-y_weft(2))/8;
       elseif y_weft(p)<=0.3 && y_weft(p)>0.2
           y_C2(p)=y_weft(p)+(-1)^p*(y_weft(1)-y_weft(2))/16;
       elseif y_weft(p)<=0.2 && y_weft(p)>0.1
           y_C2(p)=y_weft(p)+(-1)^p*(y_weft(1)-y_weft(2))/62;
       elseif y_weft(p)<=0.1
           y_C2(p)=y_weft(p)+(-1)^p*(y_weft(1)-y_weft(2))/64;
       end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_weft_fe=x_weft*plot_tr.fe/2;
        
    if length(x_weft)<200
        
        x_C2_disp=x_weft_fe;
        y_C2_disp=y_C2;
        C2_disp=C2;
    
    else
        
        x_C2_disp=x_weft_fe(1);            
        y_C2_disp=y_C2(1);            
        C2_disp=C2(1,:);
        
        for p=1:N_dec
            
%             %affichage de la partie haute fréquence des coef de l'arbre
% 
%             x_C2_disp=[x_C2_disp x_weft_fe(2^(p+1)-p:2^(p+1)-1)];
% 
%             y_C2_disp=[y_C2_disp y_C2(2^(p+1)-p:2^(p+1)-1)];
% 
%             C2_disp=[C2_disp ; C2(2^(p+1)-p:2^(p+1)-1,:)];
            
            %affichage de la partie basse fréquence
            
            x_C2_disp=[x_C2_disp x_weft_fe(2^p:2^p+(p-1))];

            y_C2_disp=[y_C2_disp y_C2(2^p:2^p+(p-1))];

            C2_disp=[C2_disp ; C2(2^p:2^p+(p-1),:)];
            
%             %affichage des coef du milieu de l'arbre
%             
%             if p>4
%    
%             x_C2_disp=[x_C2_disp x_weft_fe(node(p,2^p/2-5):node(p,2^p/2+5))];
% 
%             y_C2_disp=[y_C2_disp y_C2(node(p,2^p/2-5):node(p,2^p/2+5))];
% 
%             C2_disp=[C2_disp ; C2(node(p,2^p/2-5):node(p,2^p/2+5),:)];

        end
            
            x_C2_disp=[x_C2_disp x_weft_fe(marg_ind)];

            y_C2_disp=[y_C2_disp y_C2(marg_ind)];

            C2_disp=[C2_disp ; C2(marg_ind,:)];
    end 
    
    
%     text(x_C2_disp(basis==1),y_C2_disp(basis==1),C2_disp(basis==1,:));          
end

xlabel('frequency');
ylabel('Amplitude + tree');
 
   
%%%%%%%%%%%%%%%%%%%%% NUMBER OF COEF PER LEVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
%    lvl=zeros(1,N_dec);
%    ind_lvl=N_dec;
%    for i=length(y_success):-1:2
%             lvl(ind_lvl)=lvl(ind_lvl)+1;
%        if y_success(i)~=y_success(i-1)
%            ind_lvl=ind_lvl-1;
%        end
%    end
%    
%    disp('----------------------------------------');
%    wght_coef=N_dec:-1:1;
%    for i=1:N_dec
%        disp(['le niveau ' num2str(i) ' représente ' num2str(lvl(i)/sum(lvl)) '% des coef choisis']);
% %        disp(['le niveau ' num2str(i) ' représente ' num2str(lvl(i)*(N_dec+1-i)/sum(lvl.*wght_coef)) '% des coef choisis']);
%    end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
   
   