function [Y,U,x_mean,x_std]=function_acp(x,y,plot_,M_phi,sub_name)
% facp(x,xnom)
% 
% ENTREES
% x         tableau des donnee (individus,caracteres)
% y         classe des individus


marker_color = {'bx','r.','g*','ksquare','s'};

[m,n] = size(x);

% centrage et normalisation du nuage de points
x_mean = mean(x,1);
x_std = std(x,1,1);
x_centre = (x - x_mean(ones(m,1),:))./x_std(ones(m,1),:);

% %calcul dela matrice de variance covariance
xcov=cov(x_centre,1);

% %calcul des valeurs propres et vecteurs propres
[v,e]=eig(xcov);
% %valeures propres mises dans l'ordre decroissant
ed=diag(e);
[lambda, rg] = sort(ed,1,'descend');
U = v(:,rg);
Y=x_centre*U;
Inertie = lambda/sum(lambda)*100;
% disp('les 7 premieres valeurs propres en pourcent')
% disp(Inertie(1:min([7 n]))')



if plot_

    axe1 = 1;
    axe2 = 2;
    
    figure;
    hold on;
    title(sub_name);


    axis([min(Y(:,axe1)) max(Y(:,axe1)) min(Y(:,axe2)), max(Y(:,axe2))]); 

    for i=1:max(y)
        plot(Y(y==i,axe1),Y(y==i,axe2),marker_color{i},'LineWidth',1.5);
    end
    
    %cr�ation du cercle de corr�lation sur le plan principale:

    x_dec=zeros(1,size(x,2));
    y_dec=zeros(1,size(x,2));
    
    for i=1:size(x,2)
        x_dec(i)=sqrt(lambda(1))*U(i,1);
        y_dec(i)=sqrt(lambda(2))*U(i,2);
    end
    
    
    figure;
    hold on;
    title('correlation circle');
    
    axis([-1 1 -1 1])
    
    base_f=zeros(4,length(M_phi));
    if length(M_phi)>1
        for i=1:length(M_phi)
            op=strfind(M_phi{i},'{');
            cl=strfind(M_phi{i},'}');
            dpth=str2num(M_phi{i}(op(1)+1:cl(1)-1));
            rang=str2num(M_phi{i}(op(2)+1:cl(2)-1));
            if dpth==1
                if rang
                    base_f(4,i)=1;
                else
                    base_f(1,i)=1;
                end
            else
                if(rang+1<=(2^dpth/4))
                    base_f(1,i)=1;
                elseif (rang+1<=2^dpth/2 && rang+1>2^dpth/4)
                    base_f(2,i)=1;
                elseif (rang+1<=3*2^dpth/4 && rang+1>2^dpth/2)
                    base_f(3,i)=1;
                else
                    base_f(4,i)=1;
                end
            end
        end
    else
        plot(x_dec,y_dec,'s','MarkerFaceColor','c','MarkerSize',4);
    end       

    color_={'r','y','g','b'};
    for i=1:4
        x_dec_decom=x_dec(find(base_f(i,:)));
        y_dec_decom=y_dec(find(base_f(i,:)));

        plot(x_dec_decom,y_dec_decom,'s','MarkerFaceColor',color_{i},'MarkerSize',6);
    end
    
    VTheta = 0:0.01:2*pi;
    XCercle = cos(VTheta);
    YCercle = sin(VTheta);
    plot(XCercle, YCercle,'color','m')
    
    if length(x_dec)<200
        text(x_dec,y_dec,M_phi);
    else
        text(x_dec(find(base_f(1,:))),y_dec(find(base_f(1,:))),M_phi(find(base_f(1,:))));
    end
end


